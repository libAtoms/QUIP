!FilePotential for the CP2K driver
!Topology must have already been calculated and PSF file written for MM or QM/MM runs.
!FilePot uses go_cp2k from the cp2k_driver_module
!reads filepot.0.xyz  >>  runs CP2K  >>  writes filepot.0.out with forces

!Initialisation for QM/MM or MM:
!   call initialise(CP2K_pot,'FilePot command=/Users/csilla/QUIP/build.darwin_x86_64_g95/cp2k_filepot property_list=pos:atom_type:atom_res_name:atom_mol_name:atom_res_number:atom_charge min_cutoff=0.0')
!Initialisation for QM:
!   call initialise(CP2K_pot,'FilePot command=/Users/csilla/QUIP/build.darwin_x86_64_g95/cp2k_filepot property_list=pos min_cutoff=0.0')

!Calc (args_str is the same as for go_cp2k):
!   call calc(CP2K_pot,ds%atoms,energy=energy,forces=f1,args_str=trim(args_str))


program cp2k_filepot

  use libatoms_module
  use cp2k_driver_module
  use quip_module


  implicit none

  real(dp)                            :: energy
  real(dp), dimension(:,:), allocatable :: f0
  real(dp), pointer                   :: forces_p(:,:)
  type(Table)                         :: embedlist
  integer                             :: i, n
  type(Inoutput)                      :: xyz
  type(Atoms)                         :: my_atoms
  character(len=STRING_LENGTH)        :: args_str
  real(dp)                            :: temp

  type(Dictionary)            :: params_in
  character(len=FIELD_LENGTH) :: Run_Type            !_MM_, QS, QMMM_EXTENDED or QMMM_CORE
  integer                     :: nproc
  character(len=FIELD_LENGTH) :: Print_PSF           !_NO_PSF_, USE_EXISTING_PSF
  character(len=FIELD_LENGTH) :: coord_file
  character(len=FIELD_LENGTH) :: new_coord_file      !output XYZ file

    call system_initialise(verbosity=normal,enable_timing=.true.)
    call system_timer('program')

      call initialise(params_in)
      call param_register(params_in, 'Run_Type', 'MM', Run_Type)
      call param_register(params_in, 'nproc', '1', nproc)
      call param_register(params_in, 'Print_PSF', 'NO_PSF', Print_PSF)
      call param_register(params_in, 'coord_file', 'filepot.0.xyz',coord_file) 
      call param_register(params_in, 'new_coord_file', 'filepot.0.out',new_coord_file) 

!      if (.not. param_read_args(params_in, do_check = .true.)) then
!        call system_abort('could not parse argument line')
!      end if

      call finalise(params_in)

      call print('Run parameters:')
      call print('  Run_Type '//Run_Type)
      call print('  nproc '//nproc)
      call print('  Print_PSF '//Print_PSF)
      call print('  coord_file '//coord_file) 
      call print('  new_coord_file '//new_coord_file) 
      call print('---------------------------------------')
      call print('')

! starts here

    call system_timer('step')
    call initialise(xyz,new_coord_file,action=OUTPUT)
  
    !read in the coordinates, create atoms object with connectivities
    call read_xyz(my_atoms,coord_file)

! CALC. CONNECTIONS

    allocate(f0(3,my_atoms%N))
    if (.not.(trim(Print_PSF).eq.'NO_PSF'.or.trim(Print_PSF).eq.'USE_EXISTING_PSF')) call system_abort('Print_PSF not known: '//trim(Print_PSF))
    write (args_str,'(a,i0,a)') 'Run_Type='//trim(Run_Type)//' nproc=',nproc,' PSF_Print='//trim(Print_PSF)//' cp2k_program=cp2k_serial'
    call print(args_str)
    call go_cp2k(my_atoms=my_atoms, forces=f0, energy=energy, args_str=args_str)

    call set_value(my_atoms%params,'energy',energy)
    call add_property(my_atoms,'force', 0.0_dp, n_cols=3)
    if (.not. assign_pointer(my_atoms, 'force', forces_p)) then
       call system_abort("filepot_read_output needed forces, but couldn't find force in my_atoms ")
    endif

    forces_p = f0

    !PRINT XYZ
    call print_xyz(my_atoms,xyz,all_properties=.true.,real_format='f17.10')

    call system_timer('step')

  deallocate(f0)

  call finalise(xyz)

  call print_title('THE')
  call print('Finished. CP2K is now having a rest, since deserved it. Bye-Bye!')
  call print_title('END')

    call system_timer('program')
  call system_finalise

contains

! momentum conservation, do not care about masses
  function sum0(force) result(force0)

    real(dp), dimension(:,:), intent(in) :: force
    real(dp), allocatable, dimension(:,:) :: force0
    integer :: i
    real(dp) :: sumF(3)

    allocate(force0(size(force,1),size(force,2)))

    do i = 1, size(force,2)
       sumF(1) = sum(force(1,1:size(force,2)))
       sumF(2) = sum(force(2,1:size(force,2)))
       sumF(3) = sum(force(3,1:size(force,2)))
    enddo
    if ((sumF(1).feq.0.0_dp).and.(sumF(2).feq.0.0_dp).and.(sumF(3).feq.0.0_dp)) then
       call print('Sum of the forces are zero.')
       force0(1:3,1:size(force,2)) = force(1:3,1:size(force,2))
    endif

    call print('Sum of the forces was '//sumF(1:3))
    sumF = sumF / size(force,2)

    do i = 1, size(force,2)
       force0(1:3,i) = force(1:3,i) - sumF(1:3)
    enddo

    do i = 1, size(force0,2)
       sumF(1) = sum(force0(1,1:size(force0,2)))
       sumF(2) = sum(force0(2,1:size(force0,2)))
       sumF(3) = sum(force0(3,1:size(force0,2)))
    enddo
    call print('Sum of the forces after mom.cons.: '//sumF(1:3))

  end function sum0

end program cp2k_filepot
