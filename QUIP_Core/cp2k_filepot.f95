!FilePotential for the CP2K driver
!FilePot uses go_cp2k from the cp2k_driver_module
!reads filepot.0.xyz  >>  runs CP2K  >>  writes filepot.0.out with forces

!Initialisation:
!   call initialise(CP2K_pot,'FilePot command=/Users/csilla/QUIP/build.darwin_x86_64_g95/cp2k_filepot property_list=pos min_cutoff=0.0')

!Calc (args_str is the same as for go_cp2k):
!   call calc(CP2K_pot,ds%atoms,energy=energy,forces=f1,args_str=trim(args_str))


program cp2k_filepot

  use libatoms_module
  use cp2k_driver_module
  use quip_module
!  use atoms_module
!  use topology_module

  implicit none

    type(Atoms)                   :: my_atoms
    real(dp)                      :: energy
    real(dp), allocatable         :: f0(:,:)
    real(dp), pointer             :: forces_p(:,:)
    type(Inoutput)                :: xyz
 
    !Input
    type(Dictionary)              :: params_in
    character(len=STRING_LENGTH)  :: args_str
    character(len=FIELD_LENGTH)   :: Run_Type
    character(len=FIELD_LENGTH)   :: Print_PSF
    character(len=FIELD_LENGTH)   :: coord_file
    character(len=FIELD_LENGTH)   :: new_coord_file
    character(len=FIELD_LENGTH)   :: cp2k_program
    character(len=FIELD_LENGTH)   :: fileroot_str
    character(len=FIELD_LENGTH)   :: basis_set_file, potential_file, dft_file, cell_file
!  logical                     :: Delete_Metal_Connections

    call system_initialise(verbosity=SILENT,enable_timing=.true.)
!    call system_initialise(verbosity=NORMAL,enable_timing=.true.)
    call system_timer('program')

    !read parameters
    call initialise(params_in)
    call param_register(params_in, 'Run_Type', 'MM', Run_Type)
    call param_register(params_in, 'Print_PSF', 'NO_PSF', Print_PSF)
    call param_register(params_in, 'coord_file', 'filepot.0.xyz',coord_file) 
    call param_register(params_in, 'new_coord_file', 'filepot.0.out',new_coord_file) 
    call param_register(params_in, 'cp2k_program', 'cp2k_serial',cp2k_program) 
    call param_register(params_in, 'root', '', fileroot_str)
    call param_register(params_in, 'basis_set_file', '', basis_set_file)
    call param_register(params_in, 'potential_file', '', potential_file)
    call param_register(params_in, 'dft_file', '', dft_file)
    call param_register(params_in, 'cell_file', '', cell_file)
!    call param_register(params_in, 'Delete_Metal_Connections', 'T', Delete_Metal_Connections)

    if (.not.param_read_args(params_in, do_check=.true.,ignore_unknown=.true.)) then
      call system_abort('could not parse argument line')
    end if

    call finalise(params_in)

    !print parameters
    call print('Run parameters:')
    call print('  Run_Type '//trim(Run_Type))
    call print('  Print_PSF '//trim(Print_PSF))
    call print('  coord_file '//trim(coord_file))
    call print('  new_coord_file '//trim(new_coord_file))
    call print('  cp2k_program '//trim(cp2k_program))
    call print('  root '//trim(fileroot_str))
    call print('  basis_set_file '//trim(basis_set_file))
    call print('  potential_file '//trim(potential_file))
    call print('  dft_file '//trim(dft_file))
    call print('  cell_file '//trim(cell_file))
    call print('---------------------------------------')
    call print('')

! starts here

    call read_xyz(my_atoms,coord_file)

    !create connectivity & topology
#ifdef HAVE_DANNY
    call set_cutoff(my_atoms,0._dp)
    call calc_connect_danny(my_atoms,cut_3body=2.8_dp)
#else
    call set_cutoff(my_atoms,0._dp)
    call calc_connect(my_atoms)
#endif
!    if (Delete_Metal_Connections) call delete_metal_connects(my_atoms)
    call map_into_cell(my_atoms)
    call calc_dists(my_atoms)
    if (trim(Run_Type).ne.'QS') then
       call create_CHARMM(my_atoms,do_CHARMM=.true.)
    endif

    !generate args_str for go_cp2k
    write (args_str,'(a)') 'Run_Type='//trim(Run_Type)//' PSF_Print='//trim(Print_PSF)
    if (trim(cp2k_program).ne.'') args_str = trim(args_str)//' cp2k_program='//trim(cp2k_program)
    if (trim(fileroot_str).ne.'') args_str = trim(args_str)//' root='//trim(fileroot_str)
    if (trim(basis_set_file).ne.'') args_str = trim(args_str)//' basis_set_file='//trim(basis_set_file)
    if (trim(potential_file).ne.'') args_str = trim(args_str)//' potential_file='//trim(potential_file)
    if (trim(dft_file).ne.'') args_str = trim(args_str)//' dft_file='//trim(dft_file)
    if (trim(cell_file).ne.'') args_str = trim(args_str)//' cell_file='//trim(cell_file)

    call print('FILEPOT |'//args_str,verbosity=SILENT)

    !call CP2K
    allocate(f0(3,my_atoms%N))
    call go_cp2k(my_atoms=my_atoms, forces=f0, energy=energy, args_str=args_str)
    !momentum conservation
    call sum0(f0)

    !write energy and forces
    call set_value(my_atoms%params,'energy',energy)
    call add_property(my_atoms,'force', 0.0_dp, n_cols=3)
    if (.not. assign_pointer(my_atoms, 'force', forces_p)) then
       call system_abort("filepot_read_output needed forces, but couldn't find force in my_atoms ")
    endif
    forces_p = f0

    call initialise(xyz,new_coord_file,action=OUTPUT)
    call print_xyz(my_atoms,xyz,all_properties=.true.,real_format='f17.10')
    call finalise(xyz)

    deallocate(f0)
    call finalise(my_atoms)
    call verbosity_push(NORMAL)
    call system_timer('FILEPOT |')
    call verbosity_push(SILENT)
    call system_finalise

contains

  !momentum conservation
  !   weighing function: 1 (simply subtract sumF/n)
  !?!   weighing function: m (keeping same acceleration on the atoms)
  subroutine sum0(force)

    real(dp), dimension(:,:), intent(inout) :: force
    integer   :: i
    real(dp)  :: sumF(3)

    do i = 1, size(force,2)
       sumF(1) = sum(force(1,1:size(force,2)))
       sumF(2) = sum(force(2,1:size(force,2)))
       sumF(3) = sum(force(3,1:size(force,2)))
    enddo

    if ((sumF(1).feq.0.0_dp).and.(sumF(2).feq.0.0_dp).and.(sumF(3).feq.0.0_dp)) then
       call print('Sum of the forces are zero.')
       return
    endif

    call print('Sum of the forces was '//sumF(1:3))
    sumF = sumF / size(force,2)

    do i = 1, size(force,2)
       force(1:3,i) = force(1:3,i) - sumF(1:3)
    enddo

    do i = 1, size(force,2)
       sumF(1) = sum(force(1,1:size(force,2)))
       sumF(2) = sum(force(2,1:size(force,2)))
       sumF(3) = sum(force(3,1:size(force,2)))
    enddo
    call print('Sum of the forces after mom.cons.: '//sumF(1:3))

  end subroutine sum0

end program cp2k_filepot
