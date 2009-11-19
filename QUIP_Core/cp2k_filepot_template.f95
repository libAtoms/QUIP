program cp2k_filepot_template

  use libatoms_module
  use cp2k_driver_template_module

  implicit none

    type(Atoms)                   :: my_atoms
    real(dp)                      :: energy
    real(dp), allocatable         :: f0(:,:)
    real(dp), pointer             :: forces_p(:,:)
    type(Inoutput)                :: xyz_io

    character(len=1024)  :: infile, outfile
    character(len=10240) :: arg, args_str
    integer :: i

    call system_initialise(verbosity=SILENT,enable_timing=.true.)
    call verbosity_push(NORMAL)
    call system_timer('cp2k_filepot_template')

    if (cmd_arg_count() < 2) &
      call system_abort("Usage: cp2k_filepot_template infile outfile [ other_cli_arg1=v1 other_cli_arg2=v2 ...]")

    call get_cmd_arg(1, infile)
    call get_cmd_arg(2, outfile)
    args_str = ""
    if (cmd_arg_count() > 2) then
      do i=3, cmd_arg_count()
	call get_cmd_arg(i, arg)
	args_str = trim(args_str) // " " // trim(arg)
      end do
    endif

    call read_xyz(my_atoms,infile)

    call print('cp2k_filepot_template args_str '//trim(args_str), ERROR)

    !call CP2K
    allocate(f0(3,my_atoms%N))
    call go_cp2k_template(at=my_atoms, f=f0, e=energy, args_str=args_str)
    !momentum conservation
    call sum0(f0)

    !write energy and forces
    call set_value(my_atoms%params,'energy',energy)
    call add_property(my_atoms,'force', 0.0_dp, n_cols=3)
    if (.not. assign_pointer(my_atoms, 'force', forces_p)) then
       call system_abort("filepot_read_output needed forces, but couldn't find force in my_atoms ")
    endif
    forces_p = f0

    call initialise(xyz_io,outfile,action=OUTPUT)
    call print_xyz(my_atoms,xyz_io,all_properties=.true.,real_format='f20.13')
    call finalise(xyz_io)

    deallocate(f0)
    call finalise(my_atoms)
    call system_timer('cp2k_filepot_template')
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

end program cp2k_filepot_template
