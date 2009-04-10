program md
  use libAtoms_module
  use QUIP_module 
  use ts_module
  use cp2k_driver_module

  implicit none

  type(Dynamicalsystem) :: ds_in, ds_fin 
  type(Potential)     :: pot
  type(Atoms)         :: at_in, at_fin, at_image
  type(inoutput)      :: xml, in_image, in_in, in_fin, file_res
  type(MetaPotential) :: metapot 
  type(TS)            :: tts 
  type(Dictionary)    :: params_in
  type(MPI_Context)    :: mpi

  real(dp)              :: energy 
  real(dp), allocatable :: forces(:,:)
  integer               :: niter, im, nfix, nprint, nimages,relax_iter 
  integer               :: freq_rep
  real(dp)              :: relax_tol, force_tol, gfac, Connect_Cutoff    
  logical               :: lrestart, lmobile_first, lmobile_last, lclimbing, loptimise,PSF
  character(len=STRING_LENGTH) :: pot_type, pot_model, param_file
  character(len=STRING_LENGTH) :: method, nomefilepos_first, nomefilepos_last
  character(len=STRING_LENGTH) :: args_str, Run_type,string_arg, args_str_calc
  real(dp), allocatable, dimension(:,:) :: conf
  character(len=FIELD_LENGTH)  :: Residue_Library
  integer                      :: Run_type_1,PSF_Print, Topology_Print
  integer, parameter           :: TOPOLOGY_NO_PSF = 0
  integer, parameter           :: TOPOLOGY_CP2K_ONLY_STEP0 = -1
  integer, parameter           :: TOPOLOGY_DRIVER_ONLY_STEP0 = -2
  integer, parameter           :: TOPOLOGY_USE_EXISTING_PSF = -3

  call system_initialise(VERBOSITY=normal)

  call initialise(mpi)

  call parse_input()
  call print_title('Initialising Potential ')
  if(pot_type.eq.'wrapper') then
    call initialise(pot, 'wrapper=.true.')
    call initialise_cp2k()
  elseif(pot_type.eq.'castep') then
    call initialise(pot, "FilePot command=./castep_driver.py property_list=species:pos")
  else
    call initialise(xml, param_file)
    call initialise(pot, trim(pot_type)//" "//trim(pot_model), xml)
  endif

  call print_title('Initialising First and Last Image')
  call Initialise(in_in, nomefilepos_first)
  call read_xyz(at_in, in_in)
  call Initialise(in_fin, nomefilepos_last)
  call read_xyz(at_fin, in_fin)

  if (trim(pot_type).eq.'wrapper'.and.Run_Type_1.ne.QS_RUN) then
   call set_cutoff(at_in, Connect_Cutoff)
   call set_cutoff(at_fin,Connect_Cutoff) 
  elseif(.not.trim(pot_type).eq.'wrapper') then
   call set_cutoff(at_in, cutoff(pot))
   call set_cutoff(at_fin,cutoff(pot)) 
  endif
 
  call initialise(ds_in,at_in)
  call initialise(ds_fin,at_fin)
!  Fix atoms (if necessary) 
  ds_in%atoms%move_mask(ds_in%atoms%N-nfix+1:ds_in%atoms%N) = 0
  allocate(forces(3,at_in%N))

  tts%cos%N = nimages 
!  if you want to initialise interpolating between the first and last image
  if(.not.lrestart) then
     if(loptimise) then
        call Initialise(metapot, "Simple", pot)  
        call Print_title('First Image Optimisation')
        if(.not.trim(pot_type).eq.'wrapper') then
           call calc_connect(ds_in%Atoms)
        else
           if (Run_Type_1.ne.QS_RUN) then
             call set_value(ds_in%atoms%params,'Library',trim(Residue_Library))
             call calc_topology(ds_in%atoms,do_CHARMM=.true.)
           endif
           write (string_arg,'(a,i0)') 'project=image', 1
           args_str_calc=trim(args_str)//" "//trim(string_arg)
        endif
        call calc(pot, ds_in%Atoms, e = energy, f = forces,args_str=args_str_calc)
        call Print("tot= " // energy )
        niter = minim(metapot, ds_in%Atoms, 'cg', relax_tol, relax_iter, 'FAST_LINMIN', do_pos=.true.,args_str=args_str_calc)
        call calc(pot, ds_in%Atoms, e = energy, f = forces,args_str=args_str)
        call Print("Etot=" // energy )
       
        call Print_title('Last Image Optimisation')
        if(.not.trim(pot_type).eq.'wrapper') then
           call calc_connect(ds_fin%Atoms)
        else 
           if(Run_Type_1.ne.QS_RUN) then
             call set_value(ds_fin%atoms%params,'Library',trim(Residue_Library))
             call calc_topology(ds_fin%atoms,do_CHARMM=.true.)
           endif
           write (string_arg,'(a,i0)') 'project=image', nimages
           args_str_calc=trim(args_str)//" "//trim(string_arg)
        endif
        call calc(pot, ds_fin%Atoms, e = energy, f = forces,args_str=args_str_calc)
        call Print("Etot=" // energy )
        niter = minim(metapot, ds_fin%atoms, 'cg', relax_tol, relax_iter, 'FAST_LINMIN', do_pos=.true.,args_str=args_str_calc)
        call calc(pot, ds_fin%Atoms, e = energy, f = forces,args_str=args_str)
        call Print("Etot=" // energy )
     endif
   
     call print_title('Initialisation of the chain of state interpolating between the first and last image')
     call initialise(tts,method,ds_in%atoms,ds_fin%atoms,gfac=gfac,freq_rep=freq_rep,lmobile_first=lmobile_first,lmobile_last=lmobile_last,lclimbing=lclimbing)

! if you want to start from previous configuration for the path
  else
     allocate(conf(tts%cos%N, 3 * at_in%N) )
     do im=1,tts%cos%N
       call Initialise(in_image, 'conf.'//im//'.xyz')
       call read_xyz(at_image, in_image)
       conf(im,:) = reshape(at_image%pos, (/3*at_image%N/) ) 
       call finalise(at_image)   
       call finalise(in_image)
     enddo
     call print_title('Initialisation of the chain of state using the guessed path')
     call initialise(tts,method,ds_in%atoms,conf, lmobile_first,lmobile_last,freq_rep,gfac,lclimbing)
  endif

  call print(tts)
  do im=1, tts%cos%N 
   call print_xyz(tts%cos%image(im)%at,'image.'//im//'.xyz',real_format='f12.5') 
  enddo 

!  call calc(tts,metapot,niter,'cg', relax_tol, relax_iter, 'FAST_LINMIN', do_pos=.true.)
  call initialise(file_res, "out.dat",OUTPUT)

  call print_title('Transition state calculation')
  call calc(tts,pot,niter,relax_tol,force_tol,relax_iter,nprint,file_res,args_str=args_str)
  call print('Number or Iterations :  ' // niter )

  do im=1, tts%cos%N
   call print_xyz(tts%cos%image(im)%at,'final.'//im//'.xyz')
  enddo

  call finalise(ds_in)
  call finalise(ds_fin)
  call finalise(at_in)
  call finalise(at_fin)
  deallocate(forces)

  call system_finalise()

  contains
    
   subroutine parse_input
      call initialise(params_in)
      call param_register(params_in, 'pot', 'castep', pot_type)
      call param_register(params_in, 'pot_model', 'FS', pot_model)
      call param_register(params_in, 'param_file', 'ip.parms.FS.xml', param_file)
      call param_register(params_in, 'first','first.xyz',nomefilepos_first)
      call param_register(params_in, 'last', 'last.xyz', nomefilepos_last)
      call param_register(params_in, 'nimages','10', nimages)
      call param_register(params_in, 'relax_iter', '30', relax_iter)
      call param_register(params_in, 'relax_tol', '0.001', relax_tol)
      call param_register(params_in, 'force_tol', '0.001', force_tol)
      call param_register(params_in, 'gfac', '0.1', gfac)
      call param_register(params_in, 'freq_rep', '5', freq_rep)
      call param_register(params_in, 'nprint', '10', nprint)
      call param_register(params_in, 'nfix', '9', nfix)
      call param_register(params_in, 'move_first', 'F', lmobile_first)
      call param_register(params_in, 'move_last', 'F', lmobile_last)
      call param_register(params_in, 'climbing', 'F', lclimbing)
      call param_register(params_in, 'optimise', 'F', loptimise)
      call param_register(params_in, 'restart', 'F', lrestart)
      call param_register(params_in, 'method', 'neb', method) !%  neb or sm
!  FOR CP2K
      call param_register(params_in, 'Run_type', 'QS', Run_type)
      call param_register(params_in, 'PSF', '.false.', PSF)
      call param_register(params_in, 'Connect_Cutoff', '1.2', Connect_cutoff)
      call param_register(params_in, 'Residue_Library', 'all_res.CHARMM.lib',Residue_Library)
    
      call print("n_args " // cmd_arg_count())
    
      if (.not. param_read_args(params_in, do_check = .true.)) then
        call print("Usage: ts_calculation [method=sm or neb] [first=file(xyz)] [last=file(xyz)]",ERROR)
        call print("  [relax_iter=i] [relax_tol=r] [force_tol=r]", ERROR)
        call print("  [gfac=r] [freq_rep=i] [nprint=i] [move_first] [move_last]")
        call system_abort("Confused by input arguments")
      end if
      call finalise(params_in) 
   end subroutine parse_input

   subroutine initialise_cp2k
      if (trim(Run_Type).eq.'MM') then
      Run_Type_1 = MM_RUN
      if (PSF) then
        Topology_Print= TOPOLOGY_DRIVER_ONLY_STEP0
        PSF_Print     = DRIVER_PRINT_AND_SAVE
        call print('PSF : DRIVER_PRINT_AND_SAVE')
      else
        Topology_Print= TOPOLOGY_USE_EXISTING_PSF
        PSF_Print     = USE_EXISTING_PSF
        call print('PSF : USE_EXISTING_PSF')
        endif
      elseif(trim(Run_Type).eq.'QS') then
        Run_Type_1     = QS_RUN
        Topology_Print = TOPOLOGY_NO_PSF
        PSF_Print      = NO_PSF
        call print('PSF : NO PSF, Topology_Print : TOPOLOGY_NO_PSF')
      else
        stop 'Error'
      endif
      write (args_str,'(a,i0,a,i0,a)') 'Run_Type=',Run_type_1,' PSF_Print=',PSF_Print,' cp2k_program=cp2k_serial'
      call print('args_str to CP2K : ' // args_str)
!!'
    end subroutine initialise_cp2k

  end program md
