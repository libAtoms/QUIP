!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     QUIP: quantum mechanical and interatomic potential simulation package
!X     
!X     Portions written by Noam Bernstein, while working at the
!X     Naval Research Laboratory, Washington DC. 
!X
!X     Portions written by Gabor Csanyi, Copyright 2006-2007.   
!X
!X     When using this software,  please cite the following reference:
!X
!X     reference
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X Potential module 
!X
!% The Potential module handles the first-level selection of the
!% desidered force field (tight-binding \texttt{TB_type}, empirical potential \texttt{IP_type} or 
!% hybrid description \texttt{FilePot_type}).
!% It contains the interfaces \texttt{Initialise}, \texttt{Finalise}, \texttt{cutoff}, 
!% \texttt{Calc}, \texttt{Print}, which have the only role to 
!% re-addressing the calls to the corresponding
!% modules. 
!% A Potential object simply contains a pointer to the desired force field type
!% (only one of the three can be selected). 
!% It is initialised with
!%>    call Initialise(pot,arg_str,io_obj,[mpi_obj])
!% where, \texttt{arg_str} is a string defining the potential type and
!% possible some additional options.
!% For interatomic potentials \texttt{arg_str} will start with 'IP'; there are several different IPs defined:
!% see documentation for the 'IP_module' module. For tight binding, 'arg_str' should start with 'TB';
!% see docuementation for 'TB_module' for more details.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! $Id: Potential.f95,v 1.60 2008-07-11 09:54:52 gc121 Exp $
! $Log: not supported by cvs2svn $
! Revision 1.59  2008/07/11 09:53:38  gc121
! added hardcoded wrapper potential, fixed bug in finite difference force
!
! Revision 1.58  2008/06/27 15:24:04  nb326
! Fix order of arguments calling TB_calc
!
! Revision 1.57  2008/05/15 19:35:03  jrk33
! Added some calc_dists() to finite difference calculation
!
! Revision 1.56  2008/05/15 18:35:31  gc121
! implemented finite difference forces in Potential. untested.
!
! Revision 1.55  2008/05/01 02:27:39  nb326
! Add setup_parallel
!
! Revision 1.54  2008/04/25 18:42:16  nb326
! Init potential_initialise_filename, pass master_only to initialis(io,...)
!
! Revision 1.53  2008/04/08 12:39:19  jrk33
! Added docs
!
! Revision 1.52  2008/03/07 23:08:22  nb326
! Move hybrid into MetaPotential
!
! Revision 1.51  2008/02/26 14:29:43  jrk33
! FilePot does not require param_str
!
! Revision 1.50  2008/02/14 22:47:52  nb326
! Use FilePot
!
! Revision 1.49  2007/12/14 18:27:58  nb326
! Remove setting of hybrid defaults - now done in definition of structure
!
! Revision 1.48  2007/11/28 15:37:51  nb326
! Zero err in potential_calc
!
! Revision 1.47  2007/11/26 17:20:49  nb326
! Add no_parallel option to potential_initialise() for doing parallel I/O but not setting up parallel calculations
!
! Revision 1.46  2007/11/14 18:40:40  nb326
! Just print this%tb in print - internal handling in type(TB) now sensible (I claim)
!
! Revision 1.45  2007/11/14 14:26:30  nb326
!  Add optional err argument to potential_calc, used only with TB so far
!
! Revision 1.44  2007/11/08 13:00:49  gc121
! cosmetics
!
! Revision 1.43  2007/11/05 21:53:23  nb326
! Trivial changes in notations
!
! Revision 1.42  2007/11/01 20:34:12  nb326
! Move hybrid stuff out, and protect with #ifdef HAVE_HYBRID
!
! Revision 1.41  2007/10/16 17:06:13  nb326
! Pass mpi communicator to read(extendable_str...) for reading parameter file in parallel
!
! Revision 1.40  2007/10/03 15:10:55  nb326
! Pass convert_to_string=.true. to read to do automatic stack size increase
!
! Revision 1.39  2007/09/21 14:30:10  nb326
! Only abort in hybrid_calc with virial if there are some QM atoms
!
! Revision 1.38  2007/09/17 15:32:58  nb326
! Lots more calls to system_timer
!
! Revision 1.37  2007/09/13 13:22:03  nb326
! Print k-points for TB object
!
! Revision 1.36  2007/09/06 15:23:30  gc121
! bug in loop termination in create_weights
!
! Revision 1.35  2007/09/06 14:56:19  gc121
! moved clearing of hybrid flags from Potential.f95 to vacancy_map_hybrid.f95
!
! Revision 1.34  2007/09/06 14:28:22  gc121
! clear out flags at the end of create_local_energy_weights
!
! Revision 1.33  2007/09/04 17:06:40  nb326
! Finalise currentlist and nextline in create_local_energy_weights to prevent mem leaks
!
! Revision 1.32  2007/09/04 14:10:15  nb326
! Finalise cluster in hybrid_calc to prevent mem leaks
!
! Revision 1.31  2007/09/04 10:38:04  gc121
! hardwired global_u into the hybrid for now
!
! Revision 1.30  2007/09/03 16:22:36  gc121
! sync
!

module Potential_module

  use libAtoms_module
  use QUIP_Common_module
  use MPI_context_module
  use IP_module
  use TB_module
  use FilePot_module

#define MAX_CUT_BONDS 4

  implicit none
  private

  public :: potential
  type potential
     type(MPI_context) :: mpi

     character(len=124) :: type_args_str

     ! only one of these can be allocated at any given time
     type(TB_type), pointer     :: tb => null()
     type(IP_type), pointer     :: ip => null()
     type(FilePot_type), pointer     :: filepot => null()
     logical :: is_wrapper
  end type potential

!% Initialise a Potential object (selecting the force field) and, if necessary, the  input file for potential parameters.
  public :: Initialise
  public :: Potential_Initialise_filename
  interface Initialise
     module procedure Potential_Initialise_inoutput, Potential_Initialise_str
  end interface Initialise

!% Finalise the Potential object
  public :: Finalise
  interface Finalise
     module procedure Potential_Finalise
  end interface Finalise

!% Print potential details 
  public :: Print
  interface Print
     module procedure Potential_Print
  end interface Print

!% Set potential cutoff
  public :: cutoff
  interface cutoff
     module procedure potential_cutoff
  end interface cutoff

!% Potential calculator for energy, forces and virial 
  public :: Calc
  interface Calc
     module procedure Potential_Calc
  end interface Calc

!% Set up what you need for parallel calculation
  public :: setup_parallel
  interface setup_parallel
     module procedure Potential_setup_parallel
  end interface setup_parallel


contains

  subroutine Potential_Initialise_filename(this, args_str, filename, mpi_obj, no_parallel)
    type(Potential), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    character(len=*), intent(in) :: filename
    type(MPI_context), intent(in), optional :: mpi_obj
    logical, intent(in), optional :: no_parallel

    type(inoutput) :: io

    ! WARNING: setting master_only=.true. may lead to failure if not all processes are calling the initialise
    call Initialise(io, filename, INPUT, master_only = .true.)
    call Initialise(this, args_str, io, mpi_obj, no_parallel)
    call Finalise(io)

  end subroutine Potential_Initialise_filename

  subroutine Potential_Initialise_inoutput(this, args_str, io_obj, mpi_obj, no_parallel)
    type(Potential), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    type(Inoutput), intent(inout) :: io_obj
    type(MPI_context), intent(in), optional :: mpi_obj
    logical, intent(in), optional :: no_parallel

    type(extendable_str) :: es
    logical my_no_parallel

    my_no_parallel = optional_default(.false., no_parallel)
    call Initialise(es)
    if (present(mpi_obj)) then
      call read(es, io_obj%unit, convert_to_string=.true., mpi_comm = mpi_obj%communicator)
    else
      call read(es, io_obj%unit, convert_to_string=.true.)
    endif
    if (my_no_parallel) then
      call Initialise(this, args_str, string(es))
    else
      call Initialise(this, args_str, string(es), mpi_obj)
    endif

    call Finalise(es)

  end subroutine Potential_Initialise_inoutput

  subroutine Potential_Initialise_str(this, args_str, param_str, mpi_obj)
    type(Potential), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    character(len=*), intent(in), optional :: param_str
    type(MPI_context), intent(in), optional :: mpi_obj

    logical is_TB, is_IP, is_FilePot, is_wrapper
    type(Dictionary) :: params

    call Finalise(this)

    call initialise(params)
    call param_register(params, 'TB', 'false', is_TB)
    call param_register(params, 'IP', 'false', is_IP)
    call param_register(params, 'FilePot', 'false', is_FilePot)
    call param_register(params, 'wrapper', 'false', is_wrapper)
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Initialise_str args_str')) then
      call system_abort("Potential_Initialise_str failed to parse args_str='"//trim(args_str)//"'")
    endif
    call finalise(params)

    if (count( (/is_TB, is_IP, is_FilePot, is_wrapper /) ) /= 1) then
      call system_abort("Potential_Initialise_str found too few or too many Potential types args_str='"//trim(args_str)//"'")
    endif

    if (is_TB) then
       allocate(this%tb)
       if(present(param_str)) then
          call Initialise(this%tb, args_str, param_str, mpi_obj = mpi_obj)
       else
          call system_abort('Potential_initialise: no param_str present during TB init')
       endif
    else if (is_IP) then
       allocate(this%ip)
       if(present(param_str)) then
          call Initialise(this%ip, args_str, param_str, mpi_obj)
       else
          call system_abort('Potential_initialise: no param_str present during IP init')
       endif
    else if (is_FilePot) then
       allocate(this%filepot)
       call Initialise(this%filepot, args_str, mpi_obj)
    else if(is_wrapper) then
       this%is_wrapper = .true.
    endif
  end subroutine Potential_Initialise_str

  subroutine Potential_Finalise(this)
    type(Potential), intent(inout) :: this

    if(associated(this%tb)) then
       call Finalise(this%tb)
       deallocate(this%tb)
    elseif(associated(this%ip)) then
       call Finalise(this%ip)
       deallocate(this%ip)
    elseif(associated(this%filepot)) then
       call Finalise(this%filepot)
       deallocate(this%filepot)
    end if
    this%is_wrapper = .false.
  end subroutine Potential_Finalise

  function Potential_cutoff(this)
    type(Potential), intent(in) :: this
    real(dp) :: Potential_cutoff               

    if(associated(this%tb)) then
       Potential_cutoff = cutoff(this%tb)
    elseif(associated(this%ip)) then
       Potential_cutoff = cutoff(this%ip)
    elseif(associated(this%filepot)) then
       Potential_cutoff = cutoff(this%filepot)
    else
       Potential_cutoff = 0.0_dp
    end if
  end function Potential_cutoff

  recursive subroutine Potential_Calc(this, at, e, local_e, f, df, virial, args_str, err, mpi_obj)
    type(Potential), intent(inout) :: this
    type(Atoms), intent(inout) :: at     !% The atoms structure to compute energy and forces
    real(dp), intent(out), optional, target :: e                   !% Total energy
    real(dp), intent(out), optional, target :: local_e(:)          !% Energy per atom    
    real(dp), intent(out), optional, target :: f(:,:)              !% Forces, dimensioned \texttt{(3,at%N)}
    real(dp), intent(out), optional, target :: df(:,:)             !% Finite difference forces, \texttt{(3,at%N)}
    real(dp), intent(out), optional, target :: virial(3,3)         !% Virial
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: err
    type(MPI_context), intent(in), optional :: mpi_obj

    integer:: i,k,n, zero_loc(1)
    real(dp):: e_plus, e_minus, pos_save, r_scale
    real(dp), parameter::delta = 1.0e-4_dp
    type(Dictionary) :: params
    logical :: single_cluster, little_clusters, dummy, do_rescale_r
    character(len=10240) :: my_args_str, new_args_str
    integer, pointer, dimension(:) :: hybrid_mark, cluster_index, termindex
    real(dp), pointer, dimension(:) :: weight_region1
    integer, allocatable, dimension(:) :: hybrid_mark_saved
    real(dp), allocatable, dimension(:) :: weight_region1_saved
    real(dp), allocatable, dimension(:,:) :: f_cluster
    logical :: do_carve_cluster
    type(Table) :: cluster_info, cut_bonds
    integer, pointer :: cut_bonds_p(:,:)
    integer :: i_inner, i_outer, n_non_term
    type(Atoms) :: cluster
    character(len=256) :: prefix_save
    real(dp), pointer :: force_ptr(:,:), df_ptr(:,:), local_e_ptr(:)
    real(dp), pointer :: e_ptr, virial_ptr(:,:)
    real(dp), target :: my_e, my_virial(3,3)
    logical :: calc_force, calc_energy, calc_local_e, calc_df, calc_virial, do_calc_force, do_calc_energy, do_calc_local_e, do_calc_df, do_calc_virial
    integer, pointer :: cluster_mark_p(:)
    integer, pointer :: old_cluster_mark_p(:)
    character(len=FIELD_LENGTH) :: old_cluster_mark_postfix 

    if (at%N <= 0) &
      call system_abort("Potential_Calc called with at%N <= 0")

    if (present(err)) err = 0

    if (present(args_str)) then
      call print('potential_calc got args_str "'//trim(args_str)//'"', VERBOSE)
      my_args_str = args_str
    else
      call print('potential_calc got not args_str', VERBOSE)
      my_args_str = ""
    endif

    call initialise(params)
    call param_register(params, 'single_cluster', 'F', single_cluster)
    old_cluster_mark_postfix=""
    call param_register(params, 'old_cluster_mark_postfix', '', old_cluster_mark_postfix)
    call param_register(params, 'carve_cluster', 'T', do_carve_cluster)
    call param_register(params, 'little_clusters', 'F', little_clusters)
    call param_register(params, 'do_rescale_r', 'F', do_rescale_r)
    call param_register(params, 'r_scale', '1.0', r_scale)

    call param_register(params, 'calc_force', 'F', calc_force)
    call param_register(params, 'calc_energy', 'F', calc_energy)
    call param_register(params, 'calc_local_e', 'F', calc_local_e)
    call param_register(params, 'calc_df', 'F', calc_df)
    call param_register(params, 'calc_virial', 'F', calc_virial)

    if (.not. param_read_line(params, my_args_str, ignore_unknown=.true.,task='Potential_Calc_str args_str') ) &
      call system_abort("Potential_calc failed to parse args_str='"//trim(my_args_str)//"'")
    call finalise(params)

    if (single_cluster .and. little_clusters) &
         call system_abort('Potential_calc: single_cluster and little_clusters options are mutually exclusive')

    if (little_clusters) then

       if (present(e) .or. present(local_e) .or. present(df) .or. present(virial) .or. .not. present(f)) &
            call system_abort('potential_calc: little_clusters option only supports calcualtion of forces, not energies, local energies or virials')

       ! must remove "little_clusters" from args_str so that recursion terminates
       call initialise(params)
       call read_string(params, my_args_str)
       call remove_value(params, 'little_clusters')
       new_args_str = write_string(params)
       call finalise(params)
       
       if (.not. assign_pointer(at, 'hybrid_mark', hybrid_mark)) &
            call system_abort('potential_calc: cannot assign pointer to hybrid_mark property ')

       if (.not. any(hybrid_mark == HYBRID_ACTIVE_MARK)) then
          f = 0.0_dp
          return
       end if

       ! Save hybrid_mark and weight_region1, because we're going to overwrite them
       allocate(hybrid_mark_saved(at%N))
       hybrid_mark_saved = hybrid_mark

       if (has_property(at, 'weight_region1')) then
          dummy = assign_pointer(at, 'weight_region1', weight_region1)
          allocate(weight_region1_saved(at%N))
          weight_region1_saved = weight_region1
       end if
       
       ! call ourselves once for each active or transition atom 
       f = 0.0_dp
       n = 0
       do i=1,at%N
          if (hybrid_mark_saved(i) /= HYBRID_ACTIVE_MARK .and. hybrid_mark_saved(i) /= HYBRID_TRANS_MARK) cycle

          if (present(mpi_obj)) then
             if (mpi_obj%active) then
                n = n + 1
                call print('potential_calc: cluster '//n//' around atom '//i//'  assigned to proc '//mod(n-1,mpi_obj%n_procs)//' of '//(mpi_obj%n_procs), VERBOSE)
                if (mod(n-1, mpi_obj%n_procs) .ne. mpi_obj%my_proc) cycle
             end if
          end if
          call print('potential_calc: constructing little_cluster around atom '//i, VERBOSE)
          hybrid_mark = HYBRID_NO_MARK
          hybrid_mark(i) = HYBRID_ACTIVE_MARK
          call create_hybrid_weights(at, new_args_str)
          cluster_info = create_cluster_info_from_hybrid_mark(at, new_args_str)
	  cluster = carve_cluster(at, new_args_str, cluster_info)
	  call finalise(cluster_info)
          allocate(f_cluster(3,cluster%N))

	  if (current_verbosity() >= NERD) then
	    prefix_save = mainlog%prefix
	    mainlog%prefix="LITTLE_CLUSTER_"//i
	    call print_xyz(cluster, mainlog, all_properties=.true.)
	    mainlog%prefix=prefix_save
	  endif
call print('ARGS0 | '//new_args_str,VERBOSE)
          call calc(this, cluster, f=f_cluster, args_str=new_args_str)
          if (do_rescale_r)  f_cluster = f_cluster*r_scale
          f(:,i) = f_cluster(:,1)
          deallocate(f_cluster)
	   call finalise(cluster)
       end do

       if (present(mpi_obj)) then
          call sum_in_place(mpi_obj, f)
       end if
       hybrid_mark = hybrid_mark_saved
       deallocate(hybrid_mark_saved)

       if (allocated(weight_region1_saved)) then
          weight_region1 = weight_region1_saved
          deallocate(weight_region1_saved)
       end if

    else if (single_cluster) then

       if (present(e) .or. present(local_e) .or. present(df) .or. present(virial) .or. .not. present(f)) &
            call system_abort('potential_calc: single_cluster option only supports calcualtion of forces, not energies, local energies or virials')

       if (.not. assign_pointer(at, 'hybrid_mark', hybrid_mark)) &
            call system_abort('potential_calc: cannot assign pointer to hybrid_mark property ')

       if (.not. any(hybrid_mark == HYBRID_ACTIVE_MARK)) then
          f = 0.0_dp
          return
       end if

       ! must remove "single_cluster" from args_str so that recursion terminates
       call initialise(params)
       call read_string(params, my_args_str)
       call remove_value(params, 'single_cluster')
       new_args_str = write_string(params)
       call finalise(params)

       ! call ourselves on a cluster formed from marked atoms
       call print('potential_calc: constructing single_cluster', VERBOSE)

       if (do_carve_cluster) then
	 call print('potential_calc: carving cluster', VERBOSE)
	 cluster_info = create_cluster_info_from_hybrid_mark(at, new_args_str)

         ! Check there are no repeated indices among the non-termination atoms in the cluster
         n_non_term = count(cluster_info%int(6,1:cluster_info%n) == 0)
         if (multiple_images(int_subtable(cluster_info,(/ (i,i=1,n_non_term) /),(/1/)))) &
              call system_abort('Potential_calc: single_cluster=T not yet implemented when cluster contains repeated periodic images')

	 cluster = carve_cluster(at, new_args_str, cluster_info)
	 call finalise(cluster_info)
	 if (current_verbosity() >= NERD) then
	   prefix_save = mainlog%prefix
	   mainlog%prefix="CLUSTER"
	   call print_xyz(cluster, mainlog, all_properties=.true.)
	   mainlog%prefix=prefix_save
	 endif
	 if (.not. assign_pointer(cluster, 'index', cluster_index)) &
	      call system_abort('potential_calc: cluster is missing index property')
	 if (.not. assign_pointer(cluster, 'termindex', termindex)) &
	      call system_abort('potential_calc: cluster is missing termindex property')
	 allocate(f_cluster(3,cluster%N))
call print('ARGS1 | '//new_args_str,VERBOSE)
	 call calc(this, cluster, f=f_cluster, args_str=new_args_str)
	 if (do_rescale_r)  f_cluster = f_cluster*r_scale

	 ! copy forces for all active and transition atoms
	 f = 0.0_dp
	 do i=1,cluster%N
	    if (termindex(i) /= 0) cycle ! skip termination atoms
	    if (hybrid_mark(cluster_index(i)) == HYBRID_ACTIVE_MARK .or. &
		hybrid_mark(cluster_index(i)) == HYBRID_TRANS_MARK) &
		  f(:,cluster_index(i)) = f_cluster(:,i)
	 end do
	 deallocate(f_cluster)   
	 call finalise(cluster)
       else ! not do_carve_cluster
	 call print('potential_calc: not carving cluster', VERBOSE)
	 cluster_info = create_cluster_info_from_hybrid_mark(at, trim(new_args_str) // " cluster_same_lattice", cut_bonds)

         !save cluster in cluster_mark property
         call add_property(at,'cluster_mark',HYBRID_NO_MARK)
         call add_property(at,'old_cluster_mark'//trim(old_cluster_mark_postfix),HYBRID_NO_MARK)
	 if (.not. assign_pointer(at, 'cluster_mark', cluster_mark_p)) &
	   call system_abort("potential_calc failed to assing pointer for cluster_mark pointer")
	 if (.not. assign_pointer(at, 'old_cluster_mark'//trim(old_cluster_mark_postfix), old_cluster_mark_p)) &
	   call system_abort("potential_calc failed to assing pointer for old_cluster_mark pointer")
         old_cluster_mark_p = cluster_mark_p
         cluster_mark_p = HYBRID_NO_MARK
!**** old_cluster_mark
         do i=1,cluster_info%N
            select case (cluster_info%str(1,i))
            case ('h_active  ')
              if (cluster_mark_p(cluster_info%int(1,i)).eq.HYBRID_NO_MARK) then
                 cluster_mark_p(cluster_info%int(1,i)) = HYBRID_ACTIVE_MARK
              else
                 call system_abort('atom '//cluster_info%int(1,i)//' has already cluster_mark '//cluster_mark_p(cluster_info%int(1,i)))
              endif
            case ('h_buffer  ')
              if (cluster_mark_p(cluster_info%int(1,i)).eq.HYBRID_NO_MARK) then
                 cluster_mark_p(cluster_info%int(1,i)) = HYBRID_BUFFER_MARK
              else
                 call system_abort('atom '//cluster_info%int(1,i)//' has already cluster_mark '//cluster_mark_p(cluster_info%int(1,i)))
              endif
            case ('h_outer_l ')
              if (cluster_mark_p(cluster_info%int(1,i)).eq.HYBRID_NO_MARK) then
                 cluster_mark_p(cluster_info%int(1,i)) = HYBRID_BUFFER_OUTER_LAYER_MARK
              else
                 call system_abort('atom '//cluster_info%int(1,i)//' has already cluster_mark '//cluster_mark_p(cluster_info%int(1,i)))
              endif
            case ('hollow    ')
              if (cluster_mark_p(cluster_info%int(1,i)).eq.HYBRID_NO_MARK) then
                 cluster_mark_p(cluster_info%int(1,i)) = HYBRID_BUFFER_MARK
              else
                 call system_abort('atom '//cluster_info%int(1,i)//' has already cluster_mark '//cluster_mark_p(cluster_info%int(1,i)))
              endif
            case ('clash     ')
              if (cluster_mark_p(cluster_info%int(1,i)).eq.HYBRID_NO_MARK) then
                 cluster_mark_p(cluster_info%int(1,i)) = HYBRID_BUFFER_MARK
              else
                 call system_abort('atom '//cluster_info%int(1,i)//' has already cluster_mark '//cluster_mark_p(cluster_info%int(1,i)))
              endif
            case ('group     ')
              if (cluster_mark_p(cluster_info%int(1,i)).eq.HYBRID_NO_MARK) then
                 cluster_mark_p(cluster_info%int(1,i)) = HYBRID_BUFFER_MARK
              else
                 call system_abort('atom '//cluster_info%int(1,i)//' has already cluster_mark '//cluster_mark_p(cluster_info%int(1,i)))
              endif
            case ('term      ')
              if (cluster_mark_p(cluster_info%int(1,i)).eq.HYBRID_NO_MARK) then
                 cluster_mark_p(cluster_info%int(1,i)) = HYBRID_TERM_MARK
              else
                 call system_abort('atom '//cluster_info%int(1,i)//' has already cluster_mark '//cluster_mark_p(cluster_info%int(1,i)))
              endif
            case default
              call system_abort('unknown mark: '//cluster_info%str(1,i))
            end select
!            mark_i = hybrid_mark_from_name(cluster_info%str(i))
!            if (mark_i.ne.HYBRID_NO_MARK) cluster_mark_p(i) = mark_i
         enddo

         !save cut bonds in cut_bonds property
	 call add_property(at, 'cut_bonds', 0, n_cols=MAX_CUT_BONDS)
	 if (.not. assign_pointer(at, 'cut_bonds', cut_bonds_p)) &
	   call system_abort("potential_calc failed to assing pointer for cut_bonds pointer")
         !zero it
         cut_bonds_p = 0
	 do i=1, cut_bonds%N
	   i_inner = cut_bonds%int(1,i)
	   i_outer = cut_bonds%int(2,i)
	   zero_loc = minloc(cut_bonds_p(:,i_inner))
	   if (cut_bonds_p(zero_loc(1),i_inner) == 0) then ! free space for a cut bond
	     cut_bonds_p(zero_loc(1),i_inner) = i_outer
	   else
	     call print("cut_bonds table:", VERBOSE)
	     call print(cut_bonds, VERBOSE)
	     call print("ERROR: potential_calc ran out of space to store cut_bonds information", ERROR)
	     call print("ERROR: inner atom " // i_inner // " already has cut_bonds to " // cut_bonds_p(:,i_inner) // &
	      " no space to add cut bond to " // i_outer, ERROR)
	     call system_abort("potential_calc out of space to store cut_bonds information")
	   endif
	 end do
	 call finalise(cut_bonds)
	 if (current_verbosity() >= ANAL) then
	   prefix_save = mainlog%prefix
	   mainlog%prefix="UNCARVED_CLUSTER"
	   call print_xyz(at, mainlog, all_properties=.true.)
	   mainlog%prefix=prefix_save
	 endif
call print('ARGS2 | '//new_args_str,VERBOSE)
	 call calc(this, at, f=f, args_str=new_args_str)
	 if (do_rescale_r)  f = f*r_scale
       endif ! do_carve_cluster
    else ! little_clusters is false..

       ! For IP, call setup_atoms() hook now in case any properties must be added.
       ! This must be done *before* we assign pointers to force, local_e etc.
       if (associated(this%ip)) then
          call setup_atoms(this%ip, at)
       end if

       do_calc_force = calc_force .or. present(f)
       do_calc_energy = calc_energy .or. present(e)
       do_calc_local_e = calc_local_e .or. present(local_e)
       do_calc_df = calc_df .or. present(df)
       do_calc_virial = calc_virial .or. present(virial)

       if(do_calc_virial .or. do_calc_energy .or. do_calc_force .or. do_calc_local_e) then

          if (do_calc_force) then
             if (.not. present(f)) then
                if (.not. has_property(at, 'force')) call add_property(at, 'force', 0.0_dp, n_cols=3)
             else
                force_ptr => f
             end if
          end if

          if (do_calc_local_e) then
             if (.not. present(local_e)) then
                if (.not. has_property(at, 'local_e')) call add_property(at, 'local_e', 0.0_dp)
             else
                local_e_ptr => local_e
             end if
          end if

          if (do_calc_energy) then
             if (.not. present(e)) then
                e_ptr => my_e
             else
                e_ptr => e
             end if
          end if

          if (do_calc_virial) then
             if (.not. present(virial)) then
                virial_ptr => my_virial
             else
                virial_ptr => virial
             end if
          end if

          ! Do pointer assignments after all properties have been added
          if (do_calc_force .and. .not. present(f)) then
             if (.not. assign_pointer(at, 'force', force_ptr)) call system_abort('Potential_calc: cannot assign force_ptr')
          end if
          if (do_calc_local_e .and. .not. present(local_e)) then
             if (.not. assign_pointer(at, 'local_e', local_e_ptr)) call system_abort('Potential_calc: cannot assign local_e_ptr')
          end if
          
          if(associated(this%tb)) then

             if (do_calc_virial) then
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, virial=virial_ptr, args_str=args_str, err=err)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, forces=force_ptr, virial=virial_ptr, args_str=args_str, err=err)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str, err=err)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, local_e=local_e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, virial=virial_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, local_e=local_e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, err=err)
                end if
             else
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, args_str=args_str, err=err)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, forces=force_ptr, args_str=args_str, err=err)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, local_e=local_e_ptr, args_str=args_str, err=err)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, local_e=local_e_ptr, forces=force_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, forces=force_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, local_e=local_e_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, local_e=local_e_ptr, forces=force_ptr, args_str=args_str, err=err)
                end if                
             end if

          elseif(associated(this%ip)) then

             if (do_calc_virial) then
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, virial=virial_ptr, args_str=args_str)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, f=force_ptr, virial=virial_ptr, args_str=args_str)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, local_e=local_e_ptr, f=force_ptr, virial=virial_ptr, args_str=args_str)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, virial=virial_ptr, args_str=args_str)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, f=force_ptr, virial=virial_ptr, args_str=args_str)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, local_e=local_e_ptr, f=force_ptr, virial=virial_ptr, args_str=args_str)
                end if
             else
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, args_str=args_str)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, f=force_ptr, args_str=args_str)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, local_e=local_e_ptr, args_str=args_str)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, local_e=local_e_ptr, f=force_ptr, args_str=args_str)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, args_str=args_str)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, f=force_ptr, args_str=args_str)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, local_e=local_e_ptr, args_str=args_str)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, local_e=local_e_ptr, f=force_ptr, args_str=args_str)
                end if
             end if

          elseif(associated(this%filepot)) then

             if (do_calc_virial) then
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, virial=virial_ptr, args_str=args_str, err=err)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, forces=force_ptr, virial=virial_ptr, args_str=args_str, err=err)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str, err=err)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, local_e=local_e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, virial=virial_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, local_e=local_e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, err=err)
                end if
             else
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, args_str=args_str, err=err)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, forces=force_ptr, args_str=args_str, err=err)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, local_e=local_e_ptr, args_str=args_str, err=err)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, local_e=local_e_ptr, forces=force_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, forces=force_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, local_e=local_e_ptr, args_str=args_str, err=err)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, local_e=local_e_ptr, forces=force_ptr, args_str=args_str, err=err)
                end if
             end if

          elseif(this%is_wrapper) then
             !
             ! put here hardcoded energy and force functions
             !
             !if(present(e)) e = wrapper_energy(at)
             !if(present(f)) call wrapper_force(at, f)
             call system_abort("Potential_Calc: hardcoded wrapper functions are not defined")
          else
             call system_abort ("Potential_Calc: potential is not initialised")
          end if

          if (calc_energy) then
             call set_value(at%params, 'energy', e_ptr)
          end if

          if (calc_virial) then
             call set_value(at%params, 'virial', virial_ptr)
          end if

          ! Copy force and local_e to properties as well if necessary
          if (calc_force .and. present(f)) then
             if (.not. has_property(at, 'force')) call add_property(at, 'force', 0.0_dp, n_cols=3)
             if (.not. assign_pointer(at, 'force', force_ptr)) call system_abort('Potential_calc: cannot assign force_ptr')
             force_ptr(:,:) = f(:,:)
          end if

          if (calc_local_e .and. present(local_e)) then
             if (.not. has_property(at, 'local_e')) call add_property(at, 'local_e', 0.0_dp)
             if (.not. assign_pointer(at, 'local_e', local_e_ptr)) call system_abort('Potential_calc: cannot assign local_e_ptr')
             local_e_ptr(:) = local_e(:)
          end if

       end if

       if (do_calc_df) then ! do forces by finite difference

          ! must remove 'calc_df' from args_str if it's there (and the rest, or properties get overwritten)
          call initialise(params)
          call read_string(params, my_args_str)
          call remove_value(params, 'calc_df')
          call remove_value(params, 'calc_force')
          call remove_value(params, 'calc_energy')
          call remove_value(params, 'calc_local_e')
          call remove_value(params, 'calc_virial')
          new_args_str = write_string(params)
          call finalise(params)

          if (.not. present(df)) then
             if (.not. has_property(at, 'df')) call add_property(at, 'df', 0.0_dp, n_cols=3)
             if (.not. assign_pointer(at, 'df', df_ptr)) call system_abort('Potential_calc: cannot assign df_ptr')
          else
             df_ptr => df
          end if

          do i=1,at%N
             do k=1,3
                pos_save = at%pos(k,i)
                at%pos(k,i) = pos_save + delta
                call calc_dists(at)
                call calc(this, at, e=e_plus, args_str=new_args_str, err=err)
                at%pos(k,i) = pos_save - delta
                call calc_dists(at)
                call calc(this, at, e=e_minus, args_str=new_args_str, err=err)
                at%pos(k,i) = pos_save
                call calc_dists(at)
                df_ptr(k,i) = (e_minus-e_plus)/(2.0_dp*delta) ! force is -ve gradient
             end do
          end do

          ! Copy to property as well if necessary
          if (calc_df .and. present(df)) then
             if (.not. has_property(at, 'df')) call add_property(at, 'df', 0.0_dp, n_cols=3)
             if (.not. assign_pointer(at, 'df', df_ptr)) call system_abort('Potential_calc: cannot assign df_ptr')
             df_ptr(:,:) = df(:,:)
          end if
       end if
    end if

  end subroutine Potential_Calc

!  function Potential_Brief_Description(this) result desc
!    type(Potential), intent(inout) :: this
!    character(30), intent(out) :: desc
!    if(associated(this%tb)) then
!       desc = ''
!    elseif(associated(this%ip)) then
!       desc = Brief_Description(this%ip)
!    elseif(associated(this%filepot)) then
!       desc = ''
!    elseif(this%is_wrapper) then
!       desc = ''
!    else
!       call system_abort ("Potential_Brief_Description: potential is not initialised")
!    end if
!
!  end function Potential_Brief_Description


  subroutine Potential_Print(this, file)
    type(Potential), intent(inout) :: this
    type(Inoutput), intent(inout),optional:: file

    if(associated(this%tb)) then
       call Print(this%tb, file=file)
    elseif(associated(this%ip)) then
       call Print(this%ip, file=file)
    elseif(associated(this%filepot)) then
       call Print(this%filepot, file=file)
    elseif(this%is_wrapper) then
       call print("Potential: wrapper potential")
    else
       call system_abort ("Potential_Print: potential is not initialised")
    end if

  end subroutine Potential_Print

  subroutine Potential_setup_parallel(this, at, e, local_e, f, virial, args_str)
    type(Potential), intent(inout) :: this
    type(Atoms), intent(inout) :: at     !% The atoms structure to compute energy and forces
    real(dp), intent(out), optional :: e                   !% Total energy
    real(dp), intent(out), optional :: local_e(:)          !% Energy per atom    
    real(dp), intent(out), optional :: f(:,:)              !% Forces, dimensioned \texttt{(3,at%N)}
    real(dp), intent(out), optional :: virial(3,3)         !% Virial
    character(len=*), intent(in), optional :: args_str

    if(associated(this%tb)) then
       return
    elseif(associated(this%ip)) then
       call setup_parallel(this%ip, at, e, local_e, f, virial)
    elseif(associated(this%filepot)) then
       return
    elseif(this%is_wrapper) then
       return
    else
       call system_abort ("Potential_Print: potential is not initialised")
    end if

  end subroutine Potential_setup_parallel


end module Potential_module
