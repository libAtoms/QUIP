
#include "error.inc"
module Potential_Precon_Minim_module


  use Potential_Module, only : Potential, potential_minimise, energy_func, gradient_func, cutoff, max_rij_change, calc, print_hook, constrain_virial_post, fix_atoms_deform_grad, pack_pos_dg, unpack_pos_dg, prep_atoms_deform_grad
  use error_module
  use system_module, only : dp, inoutput, print, PRINT_ALWAYS, PRINT_NORMAL, PRINT_VERBOSE, PRINT_NERD, initialise, finalise, INPUT, &
   optional_default, current_verbosity, mainlog, round, verbosity_push_decrement, verbosity_push,verbosity_push_increment, verbosity_pop, print_warning, system_timer, system_abort, operator(//)
  use units_module, only : GPA
  use periodictable_module, only :  ElementCovRad, ElementMass
  use extendable_str_module, only : extendable_str, initialise, read, string, finalise
  use linearalgebra_module , only : norm, trace, matrix3x3_det, normsq, least_squares, add_identity, inverse, diagonalise, symmetric_linear_solve, operator(.fne.), operator(.mult.), operator(.feq.), print
  use dictionary_module, only : dictionary, STRING_LENGTH, lookup_entry_i, write_string, get_value, has_key, read_string, set_value, remove_value, initialise, finalise
  use paramreader_module, only : param_register, param_read_line
  use mpi_context_module, only : mpi_context
  use table_module, only : table, find, int_part, wipe, append, finalise
  use minimization_module , only : minim, n_minim, fire_minim, test_gradient, n_test_gradient, precon_data, preconminim, precondimer
  use connection_module, only : connection
  use atoms_types_module, only : atoms, assign_pointer, add_property, assign_property_pointer, add_property_from_pointer, diff_min_image, distance_min_image
  use atoms_module, only : has_property, cell_volume, neighbour, n_neighbours, set_lattice, is_nearest_neighbour, &
   get_param_value, remove_property, calc_connect, set_cutoff, set_param_value, calc_dists, atoms_repoint, finalise, assignment(=)
  use cinoutput_module, only : cinoutput, write
  use dynamicalsystem_module, only : dynamicalsystem, ds_print_status, advance_verlet1, advance_verlet2
  use clusters_module, only : HYBRID_ACTIVE_MARK, HYBRID_NO_MARK, HYBRID_BUFFER_MARK, create_embed_and_fit_lists_from_cluster_mark, create_embed_and_fit_lists, &
   create_hybrid_weights, add_cut_hydrogens, create_cluster_info_from_mark, bfs_grow, carve_cluster

  implicit none
  
  public :: Precon_Minim
  interface Precon_Minim
     module procedure Precon_Potential_Minim
  end interface
 
  public :: Precon_Dimer
  interface Precon_Dimer
     module procedure Precon_Potential_Dimer
  end interface

!  public :: Precon_MEP
!  interface Precon_MEP
!     module procedure Precon_Potential_MEP
!  end interface
 
  contains


  subroutine allocate_precon(this,at,precon_id,nneigh,energy_scale,length_scale,cutoff,res2,max_iter,max_sub)
    type (precon_data), intent(inout) :: this
    type (Atoms) :: at
    character(*) :: precon_id
    integer :: nneigh,max_iter,max_sub
    real(dp) :: energy_scale, length_scale, cutoff,res2

    this%precon_id = precon_id
    this%nneigh = nneigh
    this%energy_scale = energy_scale
    this%length_scale = length_scale
    this%cutoff = cutoff
    this%res2 = res2
    this%mat_mult_max_iter = max_iter
    this%max_sub = max_sub


    if (has_property(at, 'move_mask')) then
       this%has_fixed = .TRUE.
    else
       this%has_fixed = .FALSE.
    end if

    allocate(this%preconrowlengths(at%N))
    allocate(this%preconindices(nneigh+1,at%N))
    
    if (trim(precon_id) == "LJ" .OR. trim(precon_id) == "ID" .OR. trim(precon_id) == "C1") then
      this%multI = .TRUE.
    elseif (trim(precon_id) == "LJdense") then
      this%dense = .true.
    else 
      call print ("Unrecognized Preconditioner, exiting")
      call exit()
    end if
    

    if (this%multI .eqv. .true.) then
      allocate(this%preconcoeffs(nneigh+1,at%N,1))
    elseif (this%dense .eqv. .true.) then
      allocate(this%preconcoeffs(nneigh+1,at%N,6))
    end if
  end subroutine
  
  subroutine build_precon(this,am_data)
  
    implicit none

    type (precon_data),intent(inout) :: this
    character(len=1) :: am_data(:)

    type (potential_minimise) :: am

    real(dp) :: conconstant 
    integer :: I,J,thisneighcount,thisneighcountlocal,thisind,thisind2
    real(dp) :: thisdist, thiscoeff
    real(dp) :: thisdiff(3) 
    integer :: nearneighcount
    logical :: did_rebuild, do_this

    call system_timer('build_precon')
    am = transfer(am_data,am)
        
    call atoms_repoint(am%minim_at)    
    call set_cutoff(am%minim_at, max(am%minim_at%cutoff, this%cutoff)) ! JRK do not decrease Atoms%cutoff
    
    call calc_connect(am%minim_at,did_rebuild=did_rebuild)
    if (did_rebuild .eqv. .true.) then    
      call print("Connectivity rebuilt by preconditioner")
      am%connectivity_rebuilt = .true.
    end if

!:     call print(this%precon_id == "C1")
!:    !Bail out if the preconditioner does not need updating
!:    if (am%connectivity_rebuilt .eqv. .false. .and. this%precon_id == "C1") then
!:      call print("Saved a recompute")
!:      return
!:    end if
    !conconstant = 1.0_dp/am%minim_at%N
    !conconstant = 1.0_dp
    conconstant = 0.1_dp

    this%preconrowlengths = 0 
    this%preconindices = 0
    this%preconcoeffs = 0.0
    !call print('Building precon with cutoff ' // this%cutoff) 
    if(this%precon_id /= "ID") then
    do I = 1,(am%minim_at%N)
      !call print("rah")  
      !if(am%minim_at%move_mask(I) == 1) then
      
      if(this%has_fixed) then
        if (am%minim_at%move_mask(I) == 0) then
          this%preconcoeffs(1,I,1) = 1.0
          this%preconindices(1,I) = I
          this%preconrowlengths(I) = 1
          cycle
        end if
      end if
      
      thisneighcount = n_neighbours(am%minim_at,I)
      thisneighcountlocal = n_neighbours(am%minim_at,I,max_dist=this%cutoff)
      !call print(thisneighcount)
       
     ! call print(thisneighcount // ' ' // thisneighcountlocal)
     ! do J=1,thisneighcount
     !   call print(neighbour(am%minim_at,I,J))
     ! end do
     ! call exit()
      
      if(thisneighcountlocal > this%nneigh) then
        call print("Not enough memory was allocated for preconditioner, increase value of nneigh, undefined behaviour (probably runtime error) follows",PRINT_ALWAYS)
      end if

      this%preconindices(1,I) = I
      this%preconcoeffs(1,I,1) = conconstant

      nearneighcount = 1
      do J = 1,thisneighcount

        thisind = neighbour(am%minim_at,I,J,distance=thisdist,diff=thisdiff,index=thisind2) 
        if (thisind > 0 .and. (thisdist <= this%cutoff)) then 
          !call print(thisdist // ' ' // this%cutoff)
          !call print(  I  // ' '// J // ' '// thisneighcount// ' ' // thisind // ' ' // thisdist)
          !call print(  I  //  ' ' // thisind // ' ' //thisind2 // ' ' // thisdist)
          

          !call writevec(reshape(this%preconcoeffs(1,I,1:6),(/6/)),'coeffs0.dat')
          if (this%precon_id == "LJ" .or. this%precon_id == "C1") then
              
            if (this%precon_id == "LJ") then
              thiscoeff = ( thisdist/this%length_scale)**(-6.0_dp)
              thiscoeff = min(thiscoeff,this%energy_scale)
            else if (this%precon_id == "C1") then
              thiscoeff = 1.0_dp
            end if
          
	    do_this = .true.
            if (this%has_fixed) then
	      do_this = (am%minim_at%move_mask(thisind) == 1 .and. thisind .ne. I)
	    end if
                  
	    if (do_this) then
                nearneighcount = nearneighcount+1
                this%preconcoeffs(1,I,1) = this%preconcoeffs(1,I,1) + thiscoeff 
                this%preconcoeffs(nearneighcount,I,1) = -thiscoeff
                this%preconindices(nearneighcount,I) = thisind

            end if
              
            !this%preconcoeffs(1,I,1) = this%preconcoeffs(1,I,1) + thiscoeff 
            !this%preconcoeffs(nearneighcount,I,1) = -thiscoeff
         
            !nearneighcount = nearneighcount+1
!          elseif (this%precon_id == "LJdense") then
!           
!            ! Coeff 1 
!            thiscoeff = this%energy_scale * (  (thisdiff(1)**2.0)*( -42.0*4.0*this%length_scale**12.0/(thisdist**16.0) &
!                                                                    +24.0*4.0*this%length_scale**6.0/(thisdist**10.0)) & 
!                                              + 12.0*this%length_scale**12.0/(thisdist**14.0) - 12.0*this%length_scale**6.0/(thisdist**8.0)) 
!            this%preconcoeffs(nearneighcount,I,1) = thiscoeff
!            this%preconcoeffs(1,I,1) = this%preconcoeffs(1,I,1) - thiscoeff 
!          
!            ! Coeff 2
!            thiscoeff = this%energy_scale * (  (thisdiff(1)*thisdiff(2))*( -42.0*4.0*this%length_scale**12.0/(thisdist**16.0) &
!                                                                    +24.0*4.0*this%length_scale**6.0/(thisdist**10.0))) 
!            this%preconcoeffs(nearneighcount,I,2) = thiscoeff
!            this%preconcoeffs(1,I,2) = this%preconcoeffs(1,I,2) - thiscoeff 
!            
!            ! Coeff 3
!            thiscoeff = this%energy_scale * (  (thisdiff(1)*thisdiff(3))*( -42.0*4.0*this%length_scale**12.0/(thisdist**16.0) &
!                                                                    +24.0*4.0*this%length_scale**6.0/(thisdist**10.0))) 
!            this%preconcoeffs(nearneighcount,I,3) = thiscoeff
!            this%preconcoeffs(1,I,3) = this%preconcoeffs(1,I,3) - thiscoeff 
!            
!            ! Coeff 4
!            thiscoeff = this%energy_scale * (  (thisdiff(2)**2.0)*( -42.0*4.0*this%length_scale**12.0/(thisdist**16.0) &
!                                                                    +24.0*4.0*this%length_scale**6.0/(thisdist**10.0)) & 
!                                              + 12.0*this%length_scale**12.0/(thisdist**14.0) - 12.0*this%length_scale**6.0/(thisdist**8.0)) 
!            this%preconcoeffs(nearneighcount,I,4) = thiscoeff
!            this%preconcoeffs(1,I,4) = this%preconcoeffs(1,I,4) - thiscoeff 
!            
!            ! Coeff 5
!            thiscoeff = this%energy_scale * (  (thisdiff(2)*thisdiff(3))*( -42.0*4.0*this%length_scale**12.0/(thisdist**16.0) &
!                                                                    +24.0*4.0*this%length_scale**6.0/(thisdist**10.0))) 
!            this%preconcoeffs(nearneighcount,I,5) = thiscoeff
!            this%preconcoeffs(1,I,5) = this%preconcoeffs(1,I,5) - thiscoeff 
! 
!            ! Coeff 6
!            thiscoeff = this%energy_scale * (  (thisdiff(3)**2.0)*( -42.0*4.0*this%length_scale**12.0/(thisdist**16.0) &
!                                                                    +24.0*4.0*this%length_scale**6.0/(thisdist**10.0)) & 
!                                              + 12.0*this%length_scale**12.0/(thisdist**14.0) - 12.0*this%length_scale**6.0/(thisdist**8.0)) 
!            this%preconcoeffs(nearneighcount,I,6) = thiscoeff
!            this%preconcoeffs(1,I,6) = this%preconcoeffs(1,I,6) - thiscoeff 
!
!            !call writevec(reshape(this%preconcoeffs(1,I,1:6),(/6/)),'coeffs1.dat')
!            !call writevec(reshape(this%preconcoeffs(nearneighcount,I,1:6),(/6/)),'coeffs2.dat')
!            !call writevec(thisdiff,'diff.dat')
!            !call exit()
!
          end if


        end if
      end do

      !call print(nearneighcount)
      this%preconrowlengths(I) = nearneighcount 
      !end if
    end do
    else if (this%precon_id == 'ID') then
    do I = 1,(am%minim_at%N)
      this%preconcoeffs(1,I,1) = 1.0
      this%preconrowlengths(I) = 1
      this%preconindices(1,I) = I
    end do 
    end if
   am%connectivity_rebuilt = .false.
    !call exit()

    call system_timer('build_precon')

  end subroutine build_precon

  subroutine getdenseC1precon(prmat,at,cutoff)

    implicit none

    real(dp), intent(inout) :: prmat(:,:) !The dense precon matrix will be stored here, it should be (N*3+9)x(N*3+9)
    type (Atoms), target :: at ! The atoms the precon will be built from
    real(dp) :: cutoff
    type (precon_data) :: pr
    character, allocatable :: am_data(:)
    character, dimension(1) :: am_mold
    type (potential_minimise) :: am
    integer :: I, J, target_elements(3), row_elements(3), thisind, am_data_size, K
    real(dp) :: scoeff

    call allocate_precon(pr,at,'C1',125,1.0_dp,1.0_dp,cutoff,10.0_dp**(-10.0),100,20)
    
    am%minim_at => at
    am_data_size = size(transfer(am, am_mold))
    allocate(am_data(am_data_size))
    am_data = transfer(am, am_data)
    call build_precon(pr,am_data) 
    prmat = 0.0
    do I = 1,9
      prmat(I,I) = 1.0
    end do
    do I = 1,size(pr%preconindices,DIM=2)
      
      !call print(pr%preconindices(1:pr%preconrowlengths(I),I))
      target_elements = (/ I*3-2+9, I*3-1+9, I*3+9 /)
      if (pr%preconrowlengths(I) >= 1) then
      do J = 1,(pr%preconrowlengths(I))
        
        thisind = pr%preconindices(J,I)
        row_elements = (/ thisind*3-2+9, thisind*3-1+9, thisind*3+9/)
       
        scoeff = pr%preconcoeffs(J,I,1) 
        !call print(scoeff)
        do K = 1,3
          prmat(row_elements(K),target_elements(K)) = prmat(row_elements(K),target_elements(K)) + scoeff
          prmat(target_elements(K),row_elements(K)) = prmat(target_elements(K),row_elements(K)) + scoeff
        end do 
      end do
      end if 
    end do  
  end subroutine 
   
  function Precon_Potential_Minim(this, at, method, convergence_tol, max_steps,efuncroutine, linminroutine, do_print, print_inoutput, print_cinoutput, &
       do_pos, do_lat, args_str,external_pressure, &
       hook_print_interval, error,precon_id,length_scale,energy_scale,precon_cutoff,nneigh,res2,mat_mult_max_iter,max_sub,infoverride)
    
    implicit none
    
    type(Potential), intent(inout), target :: this !% potential to evaluate energy/forces with
    type(Atoms), intent(inout), target :: at !% starting configuration
    character(*), intent(in)    :: method !% passed to precon_minim()
! Options for method
! 'preconSD' - preconditioned steepest descent
! 'preconCG' - preconditioned safeguarded Polak-Ribiere conjugate gradient
! 'preconLBFGS - preconditioned LBFGS

    real(dp),     intent(in)    :: convergence_tol !% Minimisation is treated as converged once $|\mathbf{\nabla}f|^2 <$
                                                    !% 'convergence_tol'. 
    integer,      intent(in)    :: max_steps  !% Maximum number of steps
    
    character(*), intent(in), optional    :: efuncroutine !% How to evaluate change in energy in a step
!Options for efuncroutine
! 'basic' - naive summation (default)
! 'kahan' - kahan summation of local energy differences (QUIPs potential must provide local_energy, if it does turn this on) 
    
    character(*), intent(in),optional    :: linminroutine !% Name of the line minisation routine to use
! Options for linminroutine
! 'basic' - simple backtracking, each iteration satisfies armijo
! 'standard' - bracket and zoom by cubic interpolation, with bisection as backup, each iteration satisfies wolfe
! 'none' - no linesearch, relies on estimating alpha from previous gradients etc, very dangerous

    logical, optional :: do_print !% if true, print configurations using minim's hook()
    type(inoutput), intent(inout), optional, target :: print_inoutput !% inoutput object to print configs to, needed if do_print is true
    type(cinoutput), intent(inout), optional, target :: print_cinoutput !% cinoutput object to print configs to, needed if do_print is true
    logical, optional :: do_pos, do_lat !% do relaxation w.r.t. positions and/or lattice (if neither is included, do both)
    character(len=*), intent(in), optional :: args_str !% arguments to pass to calc()
    real(dp), dimension(3,3), optional :: external_pressure
    integer, intent(in), optional :: hook_print_interval !% how often to print xyz from hook function
    integer, intent(out), optional :: error !% set to 1 if an error occurred during minimisation
    
    character(*), intent(in),optional :: precon_id
! Eventually will support multiple preconditioners for now just supports:
! 'ID' - identity preconditioner (default) i.e. no preconditioner
! 'C1' - binary connectivity precontioner
! 'LJ' - connectivity preconditioner based on LJ potential 

    real(dp), intent(in), optional :: length_scale !length scale of potential (reference lattice distance, approximately will be probably good enough), default 1.0
    real(dp), intent(in), optional :: energy_scale !prefactor of the potential energy, default 1.0
    real(dp), intent(in), optional :: precon_cutoff !cutoff radius of the preconditioner, default 1.5, probably set this midway between first and second neighbour distances, assuming first neighbours contribute much more to the energy, if the potential is more or less 'flat' then may need to include second neighbours
    integer, intent(in), optional :: nneigh !maximum number of neighbours expected in precon_cutoff radius, may be removed when this becomes automatic, default is 125, only necessary to edit if you run out of memory
    
    real(dp), intent(in), optional :: res2 !criteria for the residual squared of the approximate preconditioner inverse, probably dont need to change this
    integer, intent(in), optional :: mat_mult_max_iter !max number of iterations of the preconditioner inverter, probably dont need to change this
    integer, intent(in), optional :: max_sub !max number of iterations of the inverter before restarting, probably dont need to change this

    real(dp), optional :: infoverride !optional override to max step in infinity norm
    
    integer:: Precon_Potential_Minim

    character(len=STRING_LENGTH) :: use_method

    integer n_iter, n_iter_tot
    real(dp), allocatable :: x(:)
    real(dp) :: deform_grad(3,3)
    logical my_do_print
    logical done

    type(potential_minimise) am
    integer :: am_data_size
    character, allocatable :: am_data(:)
    character, dimension(1) :: am_mold
    integer :: status

    type (precon_data) :: pr
    real(dp) :: my_length_scale, my_energy_scale, my_precon_cutoff, my_res2
    integer :: my_nneigh, my_mat_mult_max_iter, my_max_sub
    character(10) :: my_precon_id



    INIT_ERROR(error)

    my_precon_id = optional_default('ID',precon_id)
    if ((present(length_scale) .eqv. .false.) .and. (trim(my_precon_id) == 'LJ')) then
      call print("No length_scale specified, defaulting to reference lattice distance = 1.0. Unless this is a good estimate preconditioning may work VERY POORLY, if in doubt estimate high, set optional argument length_scale.")
    end if
  
!    if ((present(energy_scale) .eqv. .false.) .and. (trim(my_precon_id) == 'LJ')) then
!      call print("No energy_scale specified, defaulting to prefactor = 1.0. Unless this is a good estimate preconditioning may work very poorly, set optional argument energy_scale.")
!    end if
  
    if ((present(precon_cutoff) .eqv. .false.) .and. (trim(my_precon_id) /= 'ID')) then
      call print("No precon cutoff specified, using the potential's own cutoff = "//cutoff(this)//". Decreasing this may improve performance, set optional argument precon_cutoff.")
    end if


    am_data_size = size(transfer(am, am_mold))
    allocate(am_data(am_data_size))

    use_method = trim(method)

    call calc_connect(at)
    am%minim_at => at
    am%pos_lat_preconditioner_factor = am%minim_pos_lat_preconditioner*am%minim_at%N

    if (present(args_str)) then
      am%minim_args_str = args_str
    else
      am%minim_args_str = ""
    endif
    am%minim_pot => this

    if (.not.present(do_pos) .and. .not. present(do_lat)) then
      am%minim_do_pos = .true.
      am%minim_do_lat = .true.
    else
      am%minim_do_pos = optional_default(.false., do_pos)
      am%minim_do_lat = optional_default(.false., do_lat)
    endif

    am%external_pressure = 0.0_dp
    if (present(external_pressure)) then
       am%external_pressure = external_pressure
       if( (am%external_pressure(1,1) .fne. am%external_pressure(2,2)) .or. &
       & (am%external_pressure(1,1) .fne. am%external_pressure(3,3)) .or. &
       & (am%external_pressure(1,2) .fne. 0.0_dp) .or. &
       & (am%external_pressure(1,3) .fne. 0.0_dp) .or. &
       & (am%external_pressure(2,1) .fne. 0.0_dp) .or. &
       & (am%external_pressure(2,3) .fne. 0.0_dp) .or. &
       & (am%external_pressure(3,1) .fne. 0.0_dp) .or. &
       & (am%external_pressure(3,2) .fne. 0.0_dp) ) then
          if(trim(use_method) /= 'fire') then
             call print_warning('Anisotrpic pressure is being used. Switching to fire_minim.')
             use_method = 'fire'
          endif
       endif
    endif

    my_do_print = optional_default(.false., do_print)

    if (my_do_print .and. .not. present(print_inoutput) .and. .not. present(print_cinoutput)) &
         call system_abort("potential_minim: do_print is true, but no print_inoutput or print_cinoutput present")

    if (my_do_print) then
       if (present(print_cinoutput)) then
          am%minim_cinoutput_movie => print_cinoutput
          if (this%is_forcemixing) &
               this%forcemixing%minim_cinoutput_movie => print_cinoutput
#ifdef HAVE_LOCAL_E_MIX
          if (this%is_local_e_mix) &
               this%local_e_mix%minim_cinoutput_movie => print_cinoutput
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
          if (this%is_oniom) &
               this%oniom%minim_cinoutput_movie => print_cinoutput
#endif /* HAVE_ONIOM */
       end if
    else
      nullify(am%minim_cinoutput_movie)
    endif

    am%minim_n_eval_e = 0
    am%minim_n_eval_f = 0

    allocate(x(9+am%minim_at%N*3))
    allocate(am%last_connect_x(size(x)))
    am%last_connect_x=1.0e38_dp
    deform_grad = 0.0_dp; call add_identity(deform_grad)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor)

    am_data = transfer(am, am_data)
    my_nneigh = optional_default(125,nneigh)
    my_energy_scale = optional_default(1.0_dp,energy_scale)
    my_length_scale = optional_default(1.0_dp,length_scale)
    my_precon_cutoff = optional_default(cutoff(this),precon_cutoff)
    my_res2 = optional_default(10.0_dp**(-5.0_dp),res2)
    my_mat_mult_max_iter = optional_default(am%minim_at%N*3,mat_mult_max_iter)
    my_max_sub = 30
    
    call allocate_precon(pr,at,my_precon_id,my_nneigh,my_energy_scale,my_length_scale,my_precon_cutoff,my_res2,my_mat_mult_max_iter,my_max_sub)
    !call print(use_method)   
    n_iter = preconminim(x, energy_func_local, gradient_func, build_precon, pr, use_method, convergence_tol, max_steps,efuncroutine=efuncroutine, linminroutine=linminroutine, &
            hook=print_hook, hook_print_interval=hook_print_interval, am_data=am_data, status=status, writehessian=writeapproxhessiangrad,gethessian=getapproxhessian,getfdhconnectivity=getfdhconnectivity,infoverride = infoverride)
!       n_iter = minim(x, energy_func, gradient_func, use_method, convergence_tol, max_steps, linminroutine, &
 !           print_hook, hook_print_interval=hook_print_interval, eps_guess=my_eps_guess, data=am_data, status=status)
 
    
    call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)
    call calc_connect(am%minim_at)
    n_iter_tot = n_iter
    done = .true.
    deallocate(am%last_connect_x)
    deallocate(x)
    call print("MINIM_N_EVAL E " // am%minim_n_eval_e // " F " // am%minim_n_eval_f // &
      " EF " // am%minim_n_eval_ef, PRINT_VERBOSE)

    Precon_Potential_Minim = n_iter_tot

    deallocate(am_data)

  end function 
  
  function Precon_Potential_Dimer(this, at, method, convergence_tol, max_steps,efuncroutine, linminroutine, do_print, print_inoutput, print_cinoutput, &
       do_pos, do_lat, args_str,external_pressure, &
       hook_print_interval, error,precon_id,length_scale,energy_scale,precon_cutoff,nneigh,res2,mat_mult_max_iter,max_sub)
    
    implicit none
    
    type(Potential), intent(inout), target :: this !% potential to evaluate energy/forces with
    type(Atoms), intent(inout), target :: at !% starting configuration
    character(*), intent(in)    :: method !% passed to precon_minim()
! Options for method
! 'preconSD' - preconditioned steepest descent
! 'preconCG' - preconditioned safeguarded Polak-Ribiere conjugate gradient

    real(dp),     intent(in)    :: convergence_tol !% Minimisation is treated as converged once $|\mathbf{\nabla}f|^2 <$
                                                    !% 'convergence_tol'. 
    integer,      intent(in)    :: max_steps  !% Maximum number of steps
    
    character(*), intent(in), optional    :: efuncroutine !% How to evaluate change in energy in a step
!Options for efuncroutine
! 'basic' - naive summation (default)
! 'kahan' - kahan summation of local energy differences (QUIPs potential must provide local_energy) 
    
    character(*), intent(in),optional    :: linminroutine !% Name of the line minisation routine to use
! Options for linminroutine
! 'basic' - simple backtracking, each iteration satisfies armijo
! 'basicpp' - cubic backtracking,  tries to satisfy armijo
! 'standard' - bracket and zoom by cubic interpolation, with bisection as backup, each iteration satisfies wolfe
! 'none' - no linesearch, relies on estimating alpha from previous gradients etc, very dangerous

    logical, optional :: do_print !% if true, print configurations using minim's hook()
    type(inoutput), intent(inout), optional, target :: print_inoutput !% inoutput object to print configs to, needed if do_print is true
    type(cinoutput), intent(inout), optional, target :: print_cinoutput !% cinoutput object to print configs to, needed if do_print is true
    logical, optional :: do_pos, do_lat !% do relaxation w.r.t. positions and/or lattice (if neither is included, do both)
    character(len=*), intent(in), optional :: args_str !% arguments to pass to calc()
    real(dp), dimension(3,3), optional :: external_pressure
    integer, intent(in), optional :: hook_print_interval !% how often to print xyz from hook function
    integer, intent(out), optional :: error !% set to 1 if an error occurred during minimisation
    
    character(*), intent(in),optional :: precon_id
! Eventually will support multiple preconditioners for now just supports:
! 'ID' - identity preconditioner (default) i.e. no preconditioner
! 'C1' - binary connectivity precontioner
! 'LJ' - connectivity preconditioner based on LJ potential 

    real(dp), intent(in), optional :: length_scale !length scale of potential (reference lattice distance, approximately will be probably good enough), default 1.0
    real(dp), intent(in), optional :: energy_scale !prefactor of the potential energy, default 1.0
    real(dp), intent(in), optional :: precon_cutoff !cutoff radius of the preconditioner, default 1.5, probably set this midway between first and second neighbour distances, assuming first neighbours contribute much more to the energy, if the potential is more or less 'flat' then may need to include second neighbours
    integer, intent(in), optional :: nneigh !maximum number of neighbours expected in precon_cutoff radius, may be removed when this becomes automatic, default is 125 
    
    real(dp), intent(in), optional :: res2 !criteria for the residual squared of the approximate preconditioner inverse, probably dont need to change this
    integer, intent(in), optional :: mat_mult_max_iter !max number of iterations of the preconditioner inverter, probably dont need to change this
    integer, intent(in), optional :: max_sub !max number of iterations of the inverter before restarting

    integer:: Precon_Potential_Dimer

    character(len=STRING_LENGTH) :: use_method

    integer n_iter, n_iter_tot
    real(dp), allocatable :: x(:)
    real(dp) :: deform_grad(3,3)
    logical my_do_print
    logical done

    type(potential_minimise) am
    integer :: am_data_size
    character, allocatable :: am_data(:)
    character, dimension(1) :: am_mold
    integer :: status

    type (precon_data) :: pr
    real(dp) :: my_length_scale, my_energy_scale, my_precon_cutoff, my_res2
    integer :: my_nneigh, my_mat_mult_max_iter, my_max_sub
    character(10) :: my_precon_id



    INIT_ERROR(error)

    my_precon_id = optional_default('ID',precon_id)
    if ((present(length_scale) .eqv. .false.) .and. (trim(my_precon_id) == 'LJ')) then
      call print("No length_scale specified, defaulting to reference lattice distance = 1.0. Unless this is a good estimate preconditioning may work VERY POORLY, if in doubt estimate high, set optional argument length_scale.")
    end if
  
    if ((present(energy_scale) .eqv. .false.) .and. (trim(my_precon_id) == 'LJ')) then
      call print("No energy_scale specified, defaulting to prefactor = 1.0. Unless this is a good estimate preconditioning may work very poorly, set optional argument energy_scale.")
    end if
  
    if ((present(precon_cutoff) .eqv. .false.) .and. (trim(my_precon_id) /= 'ID')) then
      call print("No precon cutoff specified, using the potential's own cutoff = "//cutoff(this)//". Decreasing this may improve performance, set optional argument precon_cutoff.")
    end if


    am_data_size = size(transfer(am, am_mold))
    allocate(am_data(am_data_size))

    use_method = trim(method)

    call calc_connect(at)
    am%minim_at => at
    am%pos_lat_preconditioner_factor = am%minim_pos_lat_preconditioner*am%minim_at%N

    if (present(args_str)) then
      am%minim_args_str = args_str
    else
      am%minim_args_str = ""
    endif
    am%minim_pot => this
    
!      am%minim_do_lat = .false.
    am%minim_do_lat = .false.
!  if (.not.present(do_pos) .and. .not. present(do_lat)) then
!      am%minim_do_pos = .true.
!      am%minim_do_lat = .true.
!    else
!      am%minim_do_pos = optional_default(.false., do_pos)
!      am%minim_do_lat = optional_default(.false., do_lat)
!    endif

    am%external_pressure = 0.0_dp
    if (present(external_pressure)) then
       am%external_pressure = external_pressure
       if( (am%external_pressure(1,1) .fne. am%external_pressure(2,2)) .or. &
       & (am%external_pressure(1,1) .fne. am%external_pressure(3,3)) .or. &
       & (am%external_pressure(1,2) .fne. 0.0_dp) .or. &
       & (am%external_pressure(1,3) .fne. 0.0_dp) .or. &
       & (am%external_pressure(2,1) .fne. 0.0_dp) .or. &
       & (am%external_pressure(2,3) .fne. 0.0_dp) .or. &
       & (am%external_pressure(3,1) .fne. 0.0_dp) .or. &
       & (am%external_pressure(3,2) .fne. 0.0_dp) ) then
          if(trim(use_method) /= 'fire') then
             call print_warning('Anisotrpic pressure is being used. Switching to fire_minim.')
             use_method = 'fire'
          endif
       endif
    endif

    my_do_print = optional_default(.false., do_print)

    if (my_do_print .and. .not. present(print_inoutput) .and. .not. present(print_cinoutput)) &
         call system_abort("potential_minim: do_print is true, but no print_inoutput or print_cinoutput present")

    if (my_do_print) then
       if (present(print_cinoutput)) then
          am%minim_cinoutput_movie => print_cinoutput
          if (this%is_forcemixing) &
               this%forcemixing%minim_cinoutput_movie => print_cinoutput
#ifdef HAVE_LOCAL_E_MIX
          if (this%is_local_e_mix) &
               this%local_e_mix%minim_cinoutput_movie => print_cinoutput
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
          if (this%is_oniom) &
               this%oniom%minim_cinoutput_movie => print_cinoutput
#endif /* HAVE_ONIOM */
       end if
    else
      nullify(am%minim_cinoutput_movie)
    endif

    am%minim_n_eval_e = 0
    am%minim_n_eval_f = 0

    allocate(x(9+am%minim_at%N*3))
    allocate(am%last_connect_x(size(x)))
    am%last_connect_x=1.0e38_dp
    deform_grad = 0.0_dp; call add_identity(deform_grad)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor)

    am_data = transfer(am, am_data)
    my_nneigh = optional_default(125,nneigh)
    my_energy_scale = optional_default(1.0_dp,energy_scale)
    my_length_scale = optional_default(1.0_dp,length_scale)
    my_precon_cutoff = optional_default(cutoff(this),precon_cutoff)
    my_res2 = optional_default(10.0_dp**(-5.0_dp),res2)
    my_mat_mult_max_iter = optional_default(am%minim_at%N*3,mat_mult_max_iter)
    my_max_sub = 30
    
    call allocate_precon(pr,at,my_precon_id,my_nneigh,my_energy_scale,my_length_scale,my_precon_cutoff,my_res2,my_mat_mult_max_iter,my_max_sub)
    !call print(use_method)   
    n_iter = precondimer(x, energy_func_local, gradient_func, build_precon, pr, use_method, convergence_tol, max_steps,efuncroutine=efuncroutine, linminroutine=linminroutine, &
            hook=print_hook, hook_print_interval=hook_print_interval, am_data=am_data, status=status, writehessian=writeapproxhessiangrad,gethessian=getapproxhessian)
!       n_iter = minim(x, energy_func, gradient_func, use_method, convergence_tol, max_steps, linminroutine, &
 !           print_hook, hook_print_interval=hook_print_interval, eps_guess=my_eps_guess, data=am_data, status=status)
 
    
    call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)
    call calc_connect(am%minim_at)
    n_iter_tot = n_iter
    done = .true.
    deallocate(am%last_connect_x)
    deallocate(x)

    Precon_Potential_Dimer = n_iter_tot

    deallocate(am_data)

  end function 


!  function Precon_Potential_MEP(this,allat,at1,at2, method, convergence_tol, max_steps, linminroutine,efuncroutine,k, do_print, print_inoutput, print_cinoutput, &
!       do_pos, do_lat, args_str,external_pressure, &
!       hook_print_interval, error,precon_id,length_scale,energy_scale,precon_cutoff,nneigh,NatMax,fix_ends,res2,mat_mult_max_iter,max_sub,deltat)
!    
!    implicit none
!    
!    type(Potential), intent(in), target :: this !% potential to evaluate energy/forces with
!    type(Atoms), intent(inout), target, dimension(:) :: allat
!    type(Atoms), intent(in), target :: at1 !% start of the chain
!    type(Atoms), intent(in), target :: at2 !% end of the chain !Number of states in the chain
!    character(*), intent(in)    :: method !% passed to precon_minim()
!! Options for method
!! 'preconSD' - preconditioned steepest descent
!! 'preconCG' - preconditioned safeguarded Polak-Ribiere conjugate gradient
!
!    real(dp),     intent(in)    :: convergence_tol !% Minimisation is treated as converged once $|\mathbf{\nabla}f|^2 <$
!                                                    !% 'convergence_tol'. 
!    integer,      intent(in)    :: max_steps  !% Maximum number of steps
!    character(*), intent(in),optional    :: linminroutine !% Name of the line minisation routine to use
!! Options for linminroutine
!! 'basic' - simple backtracking, each iteration satisfies armijo
!! 'basicpp' - cubic backtracking,  tries to satisfy armijo
!! 'standard' - bracket and zoom by cubic interpolation, with bisection as backup, each iteration satisfies wolfe
!! 'none' - no linesearch, relies on estimating alpha from previous gradients etc, very dangerous
!
!    character(*), intent(in), optional    :: efuncroutine !% How to evaluate change in energy in a step
!    real(dp), optional ::k
!    logical, optional :: do_print !% if true, print configurations using minim's hook()
!    type(inoutput), intent(inout), optional, target :: print_inoutput !% inoutput object to print configs to, needed if do_print is true
!    type(cinoutput), intent(inout), optional, target :: print_cinoutput !% cinoutput object to print configs to, needed if do_print is true
!    logical, optional :: do_pos, do_lat !% do relaxation w.r.t. positions and/or lattice (if neither is included, do both)
!    character(len=*), intent(in), optional :: args_str !% arguments to pass to calc()
!    real(dp), dimension(3,3), optional :: external_pressure
!    integer, intent(in), optional :: hook_print_interval !% how often to print xyz from hook function
!    integer, intent(out), optional :: error !% set to 1 if an error occurred during minimisation
!    
!    character(*), intent(in),optional :: precon_id
!! Eventually will support multiple preconditioners for now just supports:
!! 'ID' - identity preconditioner (default) i.e. no preconditioner
!! 'C1' - binary connectivity precontioner
!! 'LJ' - connectivity preconditioner based on LJ potential 
!
!    real(dp), intent(in), optional :: length_scale !length scale of potential (reference lattice distance, approximately will be probably good enough), default 1.0
!    real(dp), intent(in), optional :: energy_scale !prefactor of the potential energy, default 1.0
!    real(dp), intent(in), optional :: precon_cutoff !cutoff radius of the preconditioner, default 1.5, probably set this midway between first and second neighbour distances
!    integer, intent(in), optional :: nneigh !maximum number of neighbours expected in precon_cutoff radius, may be removed when this becomes automatic, default is 125 
!    integer :: status
!    integer, intent(in), optional :: NatMax
!    logical, intent(in), optional :: fix_ends
!    real(dp), intent(in), optional :: res2 !criteria for the residual squared of the approximate preconditioner inverse, probably dont need to change this
!    integer, intent(in), optional :: mat_mult_max_iter !max number of iterations of the preconditioner inverter, probably dont need to change this
!    integer, intent(in), optional :: max_sub !max number of iterations of the inverter before restarting
!    real(dp), optional :: deltat
!
!    integer :: Precon_Potential_MEP
!    type (potential_minimise), allocatable, dimension(:) :: am 
!    real (dp), allocatable, dimension(:) :: x1, x2
!    real (dp), allocatable, dimension(:,:) :: allx
!    integer :: I, Nd, Nat
!    integer :: am_data_size
!    character, dimension(1) :: am_mold
!    real(dp) :: deform_grad(3,3)
!    logical my_do_print
!    character(len=100) :: use_method
!    character, allocatable, dimension(:,:) :: am_data
!    real(dp) :: var 
!    real(dp) :: myk
!    
!    type (precon_data), allocatable, dimension (:) :: pr
!    real(dp) :: my_length_scale, my_energy_scale, my_precon_cutoff, my_res2
!    integer :: my_nneigh, my_mat_mult_max_iter,my_max_sub
!    character(10) :: my_precon_id
!    integer :: n_iter
!    logical :: my_fix_ends
!
!    my_fix_ends = .true.
!    if( present(fix_ends) ) my_fix_ends = fix_ends
!    
!    my_precon_id = optional_default('ID',precon_id)
!    if ((present(length_scale) .eqv. .false.) .and. (trim(my_precon_id) == 'LJ')) then
!      call print("No length_scale specified, defaulting to reference lattice distance = 1.0. Unless this is a good estimate preconditioning may work VERY POORLY, if in doubt estimate high, set optional argument length_scale.")
!    end if
!  
!    if ((present(energy_scale) .eqv. .false.) .and. (trim(my_precon_id) == 'LJ')) then
!      call print("No energy_scale specified, defaulting to prefactor = 1.0. Unless this is a good estimate preconditioning may work very poorly, set optional argument energy_scale.")
!    end if
!  
!    if ((present(precon_cutoff) .eqv. .false.) .and. (trim(my_precon_id) /= 'ID')) then
!      call print("No precon cutoff specified, using the potential's own cutoff = "//cutoff(this)//". Decreasing this may improve performance, set optional argument precon_cutoff.")
!    end if
!
!    myk = 1.0
!    if (present(k)) myk=k
!     
!    Nat = size(allat)
!    allocate(am(Nat))
!    allocate(pr(Nat))
! 
!    ! build the Atoms types 
!    do I = 1,Nat
!      call initialise(allat(I),at1%N,at1%lattice)
!    end do
!     
!    ! just for now most Atoms are copies of at1
!    do I = 1,(Nat-1)
!      allat(I) = at1
!    end do
!    allat(Nat) = at2
!
!
!    Nd = 9+at1%N*3
!    call print("Nd " // Nd)
!    
!    ! setup minimization
!    use_method = trim(method)
!    do I = 1,Nat
!    call calc_connect(allat(I))
!    am(I)%minim_at => allat(I)
!    am(I)%pos_lat_preconditioner_factor = am(I)%minim_pos_lat_preconditioner*am(I)%minim_at%N
!
!    if (present(args_str)) then
!      am(I)%minim_args_str = args_str
!    else
!      am(I)%minim_args_str = ""
!    endif
!    am(I)%minim_pot => this
!
!    if (.not.present(do_pos) .and. .not. present(do_lat)) then
!      am(I)%minim_do_pos = .true.
!      am(I)%minim_do_lat = .true.
!    else
!      am(I)%minim_do_pos = optional_default(.false., do_pos)
!      am(I)%minim_do_lat = optional_default(.false., do_lat)
!    endif
!
!    am(I)%external_pressure = 0.0_dp
!    if (present(external_pressure)) then
!       am(I)%external_pressure = external_pressure
!       if( (am(I)%external_pressure(1,1) .fne. am(1)%external_pressure(2,2)) .or. &
!       & (am(I)%external_pressure(1,1) .fne. am(1)%external_pressure(3,3)) .or. &
!       & (am(I)%external_pressure(1,2) .fne. 0.0_dp) .or. &
!       & (am(I)%external_pressure(1,3) .fne. 0.0_dp) .or. &
!       & (am(I)%external_pressure(2,1) .fne. 0.0_dp) .or. &
!       & (am(I)%external_pressure(2,3) .fne. 0.0_dp) .or. &
!       & (am(I)%external_pressure(3,1) .fne. 0.0_dp) .or. &
!       & (am(I)%external_pressure(3,2) .fne. 0.0_dp) ) then
!          if(trim(use_method) /= 'fire') then
!             call print_warning('Anisotrpic pressure is being used. Switching to fire_minim.')
!             use_method = 'fire'
!          endif
!       endif
!    endif
!
!    my_do_print = optional_default(.false., do_print)
!
!    if (my_do_print .and. .not. present(print_inoutput) .and. .not. present(print_cinoutput)) &
!         call system_abort("potential_minim: do_print is true, but no print_inoutput or print_cinoutput present")
!
!    if (my_do_print) then
!       if (present(print_cinoutput)) then
!          am(I)%minim_cinoutput_movie => print_cinoutput
!          if (this%is_forcemixing) &
!               this%forcemixing%minim_cinoutput_movie => print_cinoutput
!#ifdef HAVE_LOCAL_E_MIX
!          if (this%is_local_e_mix) &
!               this%local_e_mix%minim_cinoutput_movie => print_cinoutput
!#endif /* HAVE_LOCAL_E_MIX */
!#ifdef HAVE_ONIOM
!          if (this%is_oniom) &
!               this%oniom%minim_cinoutput_movie => print_cinoutput
!#endif /* HAVE_ONIOM */
!       end if
!    else
!      nullify(am(I)%minim_cinoutput_movie)
!    endif
!
!    am(I)%minim_n_eval_e = 0
!    am(I)%minim_n_eval_f = 0
!
!    allocate(am(I)%last_connect_x(Nd))
!    am(I)%last_connect_x=1.0e38_dp
!    end do
!
!    allocate(x1(Nd))
!    allocate(x2(Nd))
!    allocate(allx(Nd,Nat))
!    
!    deform_grad = 0.0_dp; call add_identity(deform_grad)
!    call pack_pos_dg(am(1)%minim_at%pos, deform_grad, x1, am(1)%pos_lat_preconditioner_factor)
!    call pack_pos_dg(am(Nat)%minim_at%pos, deform_grad, x2, am(Nat)%pos_lat_preconditioner_factor)
!    allx(1:Nd,1) = x1
!    allx(1:Nd,Nat) = x2
!        
!    ! place atoms linearly along the path 
!    do I=2,(Nat-1)
!      var = (I-1.0)/(Nat-1.0)
!      allx(1:Nd,I) = x1 + var*(x2-x1) 
!      call unpack_pos_dg(allx(1:Nd,I) , am(I)%minim_at%N, am(I)%minim_at%pos, deform_grad, 1.0_dp/am(I)%pos_lat_preconditioner_factor)
!    end do
!
!    am_data_size = size(transfer(am, am_mold))
!
!    my_nneigh = optional_default(125,nneigh)
!    my_energy_scale = optional_default(1.0_dp,energy_scale)
!    my_length_scale = optional_default(1.0_dp,length_scale)
!    my_precon_cutoff = optional_default(cutoff(this),precon_cutoff)
!    my_res2 = optional_default(10.0_dp**(-10.0_dp),res2)
!    my_mat_mult_max_iter = optional_default(am(1)%minim_at%N*3,mat_mult_max_iter)
!    my_max_sub = optional_default(200,max_sub)
!
!    allocate(am_data(am_data_size,Nat))     
!    do I=1,Nat
!      am_data(1:am_data_size,I) = transfer(am(I), am_data(1:am_data_size,I))
!      call allocate_precon(pr(I),allat(I),my_precon_id,my_nneigh,my_energy_scale,my_length_scale,my_precon_cutoff,my_res2,my_mat_mult_max_iter,my_max_sub)
!    end do
!   
!    n_iter = preconMEP(allx,energy_func, gradient_func, build_precon, pr, use_method,my_fix_ends, convergence_tol, max_steps, linminroutine=linminroutine,efuncroutine=efuncroutine,k=myk, &
!               hook=print_hook, hook_print_interval=hook_print_interval, am_data=am_data,status=status,deltat=deltat)
!
!    do I=1,Nat
!      call unpack_pos_dg(allx(1:Nd,I) , am(I)%minim_at%N, am(I)%minim_at%pos, deform_grad, 1.0_dp/am(I)%pos_lat_preconditioner_factor)
!    end do
!   
!
!  
!  end function
 
  ! compute energy
  function energy_func_local(x, am_data, local_energy_inout,gradient_inout)
    real(dp) :: x(:)
    character(len=1), optional :: am_data(:)
    real(dp) :: energy_func_local
    real(dp), intent(inout),optional :: local_energy_inout(:), gradient_inout(:)

    real(dp) :: max_atom_rij_change
    real(dp) :: deform_grad(3,3)
    type(potential_minimise)  :: am
    real(dp) , allocatable :: f(:,:)
    real(dp) :: virial(3,3),  deform_grad_inv(3,3)
    integer :: i
    integer, pointer, dimension(:) :: move_mask, fixed_pot
    logical :: did_rebuild
    
    call system_timer("energy_func")

    if (.not. present(am_data)) call system_abort("potential_minimise energy_func must have am_data")
    am = transfer(am_data, am)
    call atoms_repoint(am%minim_at)

    am%minim_n_eval_e = am%minim_n_eval_e + 1

    if (size(x) /= am%minim_at%N*3+9) call system_abort("Called energy_func() with size mismatch " // &
      size(x) // " /= " // am%minim_at%N // "*3+9")

    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    max_atom_rij_change = max_rij_change(am%last_connect_x, x, cutoff(am%minim_pot), &
      1.0_dp/am%pos_lat_preconditioner_factor)
    !call print(max_atom_rij_change// ' '//am%minim_at%cutoff // ' '// cutoff(am%minim_pot))
    !call print(am%minim_at%connect%neighbour1(302)%t%int(1,1:))

    if (current_verbosity() >= PRINT_NERD) call print("energy_func got x " // x, PRINT_NERD)

    call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

    call calc_connect(am%minim_at,did_rebuild=did_rebuild)
    am%last_connect_x = x

    if (did_rebuild .eqv. .true.) then
      call print("Connectivity rebuilt in objective function",PRINT_NERD)
      am%connectivity_rebuilt = .true.
    end if
 
    if (current_verbosity() >= PRINT_NERD) then
       call print("energy_func using am%minim_at", PRINT_NERD)
       call write(am%minim_at, 'stdout')
    end if

    if( .not. present(gradient_inout) ) then
      call calc(am%minim_pot, am%minim_at, energy = energy_func_local, local_energy = local_energy_inout, args_str = am%minim_args_str)
    else
      allocate(f(3,am%minim_at%N))
      f = 0.0_dp
      virial = 0.0_dp
      if (am%minim_do_pos .and. am%minim_do_lat) then
        call calc(am%minim_pot, am%minim_at, energy = energy_func_local, local_energy = local_energy_inout, force = f, virial = virial, args_str = am%minim_args_str)
      else if (am%minim_do_pos) then
        call calc(am%minim_pot, am%minim_at, energy = energy_func_local, local_energy = local_energy_inout, force = f, args_str = am%minim_args_str)
      else
        call calc(am%minim_pot, am%minim_at, energy = energy_func_local, local_energy = local_energy_inout, virial = virial, args_str = am%minim_args_str)
      endif


    ! zero forces if fixed by potential
      if (am%minim_do_pos .and. assign_pointer(am%minim_at, "fixed_pot", fixed_pot)) then
        do i=1, am%minim_at%N
          if (fixed_pot(i) /= 0) f(:,i) = 0.0_dp
        end do
      endif

    ! Zero force on any fixed atoms
      if (assign_pointer(am%minim_at, 'move_mask', move_mask)) then
        do i=1,am%minim_at%N
          if (move_mask(i) == 0) f(:,i) = 0.0_dp
        end do
      end if

      if (current_verbosity() >= PRINT_NERD) then
         call print ("gradient_func got f", PRINT_NERD)
         call print(f, PRINT_NERD)
         call print ("gradient_func got virial", PRINT_NERD)
         call print(virial, PRINT_NERD)
      end if

      virial = virial - am%external_pressure*cell_volume(am%minim_at)

      if (current_verbosity() >= PRINT_NERD) then
         call print ("gradient_func got virial, external pressure subtracted", PRINT_NERD)
         call print(virial, PRINT_NERD)
      end if

      f = transpose(deform_grad) .mult. f

      call inverse(deform_grad, deform_grad_inv)
      virial = virial .mult. transpose(deform_grad_inv)

      call constrain_virial_post(am%minim_at, virial)

      call pack_pos_dg(-f, -virial, gradient_inout, 1.0_dp/am%pos_lat_preconditioner_factor)
    end if
  
    call print ("energy_func got energy " // energy_func_local, PRINT_NERD)
    energy_func_local = energy_func_local + cell_volume(am%minim_at)*trace(am%external_pressure) / 3.0_dp
    call print ("energy_func got enthalpy " // energy_func_local, PRINT_NERD)

    call fix_atoms_deform_grad(deform_grad, am%minim_at, am)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor)

    am_data = transfer(am, am_data)

    call system_timer("energy_func")

  end function energy_func_local

  ! utility function to dump a vector into a file (for checking,debugging)
  subroutine writevec(vec,filename)
    real(dp) :: vec(:)
    character(*) :: filename

    integer :: outid = 10
    open(unit=outid,file=filename,action="write",status="replace")
    write(outid,*) vec
    close(outid)
  
  end subroutine
 
  subroutine writemat(mat,filename)
    real(dp) :: mat(:,:)
    character(*) :: filename

    integer :: outid = 10
    open(unit=outid,file=filename,action="write",status="replace")
    write(outid,*) mat
    close(outid)
  
  end subroutine
 
  subroutine writematint(mat,filename)
    integer :: mat(:,:)
    character(*) :: filename

    integer :: outid = 10
    open(unit=outid,file=filename,action="write",status="replace")
    write(outid,*) mat
    close(outid)
  
  end subroutine

  subroutine writeapproxhessian(x,data,filename)
    
    implicit none

    real(dp) :: x(:)
    character(len=1) :: data(:)
    character(*) :: filename

    type(potential_minimise) :: am

    real(dp),parameter :: eps = 10.0**(-5)
    integer,parameter :: nneigh = 500
    real(dp), allocatable :: hess(:)
    integer, allocatable :: hessinds(:,:)! ,hxcounts(:)
    integer :: I, J, thisneighcount, thisind, thisotherind, II, JJ, hxcount
    real(dp), allocatable :: xpp(:), xpm(:), xmm(:), xmp(:), lepp(:), lepm(:), lemm(:), lemp(:)
    real(dp) :: fpp, fpm, fmm, fmp
    integer :: hx,hy
    real(dp) :: fx
    am = transfer(data,am)
    
    !call atoms_repoint(am%minim_at)    
    !call calc_connect(am%minim_at)
    !call calc_dists(am%minim_at)
 
    !call fix_atoms_deform_grad(deform_grad, am%minim_at, am)
    !call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor)

    allocate(xpp(size(x)),xpm(size(x)),xmm(size(x)),xmp(size(x)))
    allocate(lepp(am%minim_at%N),lepm(am%minim_at%N),lemm(am%minim_at%N),lemp(am%minim_at%N))
    allocate(hess(3*am%minim_at%N*nneigh))
    allocate(hessinds(3*am%minim_at%N*nneigh,2))
    !allocate(hxcounts(3*am%minim_at%N))
    hessinds = 0
    hess = 0.0
    fx = energy_func_local(x,data,lepp)
    fx = KahanSum(lepp)
    hxcount = 1
    do I = 1,(am%minim_at%N) ! Loop over atoms
      call print(I)
      thisneighcount = n_neighbours(am%minim_at,I)
      !call print(x)
      !call exit() 
      if(thisneighcount > 3*nneigh) then
        call print("Not enough memory was allocated for Hessian, increase value of nneigh")
      end if
            
      do II = 1,3 ! Loop over coordinates of atom I
        hx = 3*(I-1) + II
                  
          do J = 1,(thisneighcount+1) ! Loop over neighbours of atom I
          if (J > 1) then
          thisind = neighbour(am%minim_at,I,J-1,index=thisotherind) 
          else
          thisind = I
          end if
          do JJ = 1,3 ! Loop over coordinates of atom J
            hy = 3*(thisind-1) + JJ
            
            if(hy .eq. hx) then
            xpp = x
            xmm = x
            xpp(hx+9) = xpp(hx+9) + eps
            xmm(hx+9) = xmm(hx+9) - eps

            fpp = energy_func_local(xpp,data,lepp)
            fmm = energy_func_local(xmm,data,lemm)
            fpp = KahanSum(lepp)
            fmm = KahanSum(lemm)
            hess(hxcount) = (fpp - 2.0*fx + fmm)/(eps**2.0)
            hessinds(hxcount,1) = hx
            hessinds(hxcount,2) = hx
            !call print(I // ' ' // II  // ' '// thisind// ' '//JJ//' ' // J //  ' '// hx // ' '//hy// ' '//hess(hxcount))
            hxcount = hxcount + 1
            
            elseif(hy .gt. hx) then
            xpp = x
            xpm = x
            xmp = x
            xmm = x
            xpp(hx+9) = xpp(hx+9) + eps
            xpp(hy+9) = xpp(hy+9) + eps
            xpm(hx+9) = xpm(hx+9) + eps
            xpm(hy+9) = xpm(hy+9) - eps
            xmp(hx+9) = xmp(hx+9) - eps
            xmp(hy+9) = xmp(hy+9) + eps
            xmm(hx+9) = xmm(hx+9) - eps
            xmm(hy+9) = xmm(hy+9) - eps
         
           
            fpp = energy_func_local(xpp,data,lepp)
            fpm = energy_func_local(xpm,data,lepm)
            fmp = energy_func_local(xmp,data,lemp)
            fmm = energy_func_local(xmm,data,lemm)
            fpp = KahanSum(lepp)
            fpm = KahanSum(lepm)
            fmp = KahanSum(lemp)
            fmm = KahanSum(lemm)
            hess(hxcount) = (fpp + fmm - fmp - fpm)/(4.0*eps**2.0)
            hessinds(hxcount,1) = hx
            hessinds(hxcount,2) = hy
            !call print(I // ' ' // II  // ' '// thisind// ' '//JJ//' ' // J //  ' '// hx // ' '//hy// ' '//hess(hxcount))
            hxcount = hxcount + 1
            endif
            !call print(fpp// ' '//fpm// ' '//fmp // ' '// fmm// ' '//fpp+fmm-fmp-fpm) 
          end do
        end do
      end do
    end do
    
    call writevec(hess,filename // 'hess.dat')
    call writematint(hessinds, filename // 'hessinds.dat')
     
  end subroutine
  
  function KahanSum(vec)
    
    real(dp) :: vec(:)
    real(dp) :: KahanSum

    integer :: I,N
    real(dp) :: C,T,Y

    N = size(vec)

    KahanSum = 0.0
    C = 0.0
    do I = 1,N
      y = vec(I) - C
      T = KahanSum + Y
      C = (T - KahanSum) - Y
      KahanSum = T
    end do
  
  end function

  subroutine writeapproxhessiangrad(x,data,filename)
    
    implicit none

    real(dp) :: x(:)
    character(len=1) :: data(:)
    character(*) :: filename

    type(potential_minimise) :: am

    real(dp),parameter :: eps = 10.0**(-3)
  
    real(dp), allocatable :: hess(:)
    integer, allocatable :: hessinds(:,:)
    real(dp),allocatable :: gp(:,:), gm(:,:), g(:), xm(:)

    integer :: I,J,II,JJ,N
    integer :: hx, hy, thisind, hxcount
    integer :: nneigh = 500
    integer :: thisneighcount, thisotherind

    N = size(x)
    am = transfer(data,am)
    allocate(xm(N),g(N))
    allocate(gp(N,N-9))
    allocate(gm(N,N-9))
    allocate(hessinds(3*am%minim_at%N*nneigh,2))
    allocate(hess(3*am%minim_at%N*nneigh))
  
    hess = 0.0
    hessinds = 0  
    gp = 0.0
    gm = 0.0
    call verbosity_push_decrement()
    do I = 10,N
      xm = x
      xm(I) = xm(I) + 1.0*eps
      g = gradient_func(xm,data)
      !call print(g)
      !call exit()
      gp(1:N,I-9) = g     
      xm(I) = xm(I) - 2.0*eps
      g = gradient_func(xm,data)
      gm(1:N,I-9) = g     
      !call print(gp(1:N,I-9)) 
    end do
    call verbosity_pop()
    !call exit()
   hxcount = 1
    do I = 1,(am%minim_at%N) ! Loop over atoms
      !call print(I)
      thisneighcount = n_neighbours(am%minim_at,I)
      !call print(x)
      !call exit() 
      if(thisneighcount > 3*nneigh) then
        call print("Not enough memory was allocated for Hessian, increase value of nneigh")
      end if
            
      do II = 1,3 ! Loop over coordinates of atom I
        hx = 3*(I-1) + II
                  
          do J = 1,(thisneighcount+1) ! Loop over neighbours of atom I
          if (J > 1) then
          thisind = neighbour(am%minim_at,I,J-1,index=thisotherind) 
          else
          thisind = I
          end if
          do JJ = 1,3 ! Loop over coordinates of atom J
            hy = 3*(thisind-1) + JJ
            if (hy >= hx) then
            !call print(hx //" "//hy//" "//gp(hy+9,hx) // " "// gm(hy+9,hx) // " "// gp(hx+9,hy)// " "//gm(hx+9,hy))
     !call print(gp(1:N,1))
    !call exit()
            hess(hxcount) = (gp(hy+9,hx) - gm(hy+9,hx))/(4.0*eps) + (gp(hx+9,hy) - gm(hx+9,hy))/(4.0*eps)
            hessinds(hxcount,1) = hx
            hessinds(hxcount,2) = hy
            !call print(I // ' ' // II  // ' '// thisind// ' '//JJ//' ' // J //  ' '// hx // ' '//hy// ' '//hess(hxcount))
            hxcount = hxcount + 1
            end if
            !call print(fpp// ' '//fpm// ' '//fmp // ' '// fmm// ' '//fpp+fmm-fmp-fpm) 
          end do
        end do
      end do
    end do
    !call exit()
    call writevec(hess,filename // 'hess.dat')
    call writematint(hessinds, filename // 'hessinds.dat')
    !call exit() 
  end subroutine 
  
  subroutine getapproxhessian(x,data,hessout)
    
    implicit none

    real(dp),intent(in) :: x(:)
    character(len=1),intent(in) :: data(:)
    real(dp),intent(inout) :: hessout(:,:)

    type(potential_minimise) :: am

    real(dp),parameter :: eps = 10.0**(-3)
  
    real(dp),allocatable :: gp(:,:), gm(:,:), g(:), xm(:)
    real(dp) :: hesscoeff
    integer :: I,J,II,JJ,N
    integer :: hx, hy, thisind, hxcount
    integer :: nneigh = 500
    integer :: thisneighcount, thisotherind

    N = size(x)
    am = transfer(data,am)
    allocate(xm(N),g(N))
    allocate(gp(N,N-9))
    allocate(gm(N,N-9))
  
    hessout = 0.0_dp
    
    gp = 0.0
    gm = 0.0
    call verbosity_push_decrement()
    do I = 10,N
      xm = x
      xm(I) = xm(I) + 1.0*eps
      g = gradient_func(xm,data)
      !call print(g)
      !call exit()
      gp(1:N,I-9) = g     
      xm(I) = xm(I) - 2.0*eps
      g = gradient_func(xm,data)
      gm(1:N,I-9) = g     
      !call print(gp(1:N,I-9)) 
    end do
    call verbosity_pop()
    !call exit()
    hxcount = 1
    do I = 1,9
      hessout(I,I) = 1
    end do
    do I = 1,(am%minim_at%N) ! Loop over atoms
      !call print(I)
      thisneighcount = n_neighbours(am%minim_at,I)
      !call print(x)
      !call exit() 
      if(thisneighcount > 3*nneigh) then
        call print("Not enough memory was allocated for Hessian, increase value of nneigh")
      end if
            
      do II = 1,3 ! Loop over coordinates of atom I
        hx = 3*(I-1) + II
                  
          do J = 1,(thisneighcount+1) ! Loop over neighbours of atom I
          if (J > 1) then
          thisind = neighbour(am%minim_at,I,J-1,index=thisotherind) 
          else
          thisind = I
          end if
          do JJ = 1,3 ! Loop over coordinates of atom J
            hy = 3*(thisind-1) + JJ
            if (hy >= hx) then
            !call print(hx //" "//hy//" "//gp(hy+9,hx) // " "// gm(hy+9,hx) // " "// gp(hx+9,hy)// " "//gm(hx+9,hy))
     !call print(gp(1:N,1))
    !call exit()
            hesscoeff = (gp(hy+9,hx) - gm(hy+9,hx))/(4.0*eps) + (gp(hx+9,hy) - gm(hx+9,hy))/(4.0*eps)
            hessout(hx+9,hy+9) = hessout(hx+9,hy+9) + hesscoeff
            if (hy > hx) then
            hessout(hy+9,hx+9) = hessout(hy+9,hx+9) + hesscoeff
            end if
            !call print(I // ' ' // II  // ' '// thisind// ' '//JJ//' ' // J //  ' '// hx // ' '//hy// ' '//hess(hxcount))
            hxcount = hxcount + 1
            end if
            !call print(fpp// ' '//fpm// ' '//fmp // ' '// fmm// ' '//fpp+fmm-fmp-fpm) 
          end do
        end do
      end do
    end do
    deallocate(xm,g,gp,gm)
    !call exit()
    !call exit() 
  end subroutine 

  subroutine getfdhconnectivity(rows,diag,rn,data)

    implicit none

    integer, intent(inout) :: rows(:), diag(:)
    character(len=1),intent(in) :: data(:)
    integer, intent(out) :: rn

    type(potential_minimise) :: am
    integer :: rowsindex
    integer :: thisneighcount,thisind
    integer :: I, J, II, JJ, hx, hy
    integer :: buffer(300)
    integer :: bufferindex
    integer :: cleanN

    diag(1:9) = (/1, 2, 3, 4, 5, 6, 7, 8, 9/)
    rows(1:9) = (/1, 2, 3, 4, 5, 6, 7, 8, 9/)
    am = transfer(data,am)
    
    rowsindex = 10 
    do I = 1,(am%minim_at%N)
      thisneighcount = n_neighbours(am%minim_at,I)
           
      do II = 1,3 ! Loop over coordinates of atom I
         
        hx = 3*(I-1) + II + 9
        diag(hx) = rowsindex 
        buffer = 3*am%minim_at%N + 10
        bufferindex = 1
        do J = 1,(thisneighcount+1) ! Loop over neighbours of atom I
          
          if (J > 1) then
            thisind = neighbour(am%minim_at,I,J-1) 
          else
            thisind = I
          end if
          !if (thisind >= I) then 
          do JJ = 1,3
            hy = 3*(thisind-1) + JJ + 9
            if (hy>=hx) then
            buffer(bufferindex) = hy
            bufferindex = bufferindex+1
            !rows(rowsindex) = hy
            !rowsindex = rowsindex + 1
            end if
          end do
          !end if
        end do
        !call print(" ")
        !call print(buffer)        
        call fdhcleaner(buffer,3*am%minim_at%N+10,cleanN)
        !call print(buffer(1:cleanN))
        rows(rowsindex:(rowsindex+cleanN-1)) = buffer(1:cleanN)
        rowsindex = rowsindex + cleanN
      end do
    end do 
    rn = rowsindex-1
  end subroutine  

 recursive function qsort( data ) result( sorted ) 
    integer, dimension(:), intent(in) :: data 
    integer, dimension(1:size(data))  :: sorted 
    if ( size(data) > 1 ) then 
      sorted = (/ qsort( pack( data(2:), abs(data(2:)) <= abs(data(1)) )  ), data(1), qsort( pack( data(2:), abs(data(2:)) > abs(data(1)) ) ) /) 
    else 
      sorted = data    
    endif 
  end function 

  subroutine fdhcleaner(data,data_max,cleanN)
  
    implicit none

    integer, intent(inout) :: data(:)
    integer, intent(in) :: data_max
    integer, intent(out) :: cleanN

    integer :: I, N

    N = size(data)
    data = qsort(data)
    
    !call print(data)
    !call exit()
    I = 1
    do
      !call print(I // ' ' // N // ' '// data(I:) //  ' '// data_max)
      if (data(I) == data(I+1)) then
        data(I+1:N-1) = data(I+2:N)
      else 
        I = I+1
      end if
      if(I == N .or. data(I) >= data_max) then
        exit
      end if
    end do
    cleanN = I-1

  end subroutine

end module
