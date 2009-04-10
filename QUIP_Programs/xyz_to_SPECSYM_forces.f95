program xyz_to_SPECSYM_forces
use libatoms_module
implicit none
  type(Atoms) :: at
  type(Dictionary) :: cli_params
  real(dp), pointer :: forces(:,:)
  integer config_n
  integer i
  real(dp) :: dipole_mom(3)

  call system_initialise()

  call initialise(cli_params)
  call param_register(cli_params, 'config_n', '0', config_n)
  if (.not. param_read_args(cli_params, do_check = .true.)) then
    call system_abort("Usage: xyz_to_SPECSYM_forces config_n")
  endif
  call finalise(cli_params)

  call read_xyz(at, "stdin")
  
  if (.not. assign_pointer(at, 'forces', forces)) &
    call system_abort("Failed to assign pointer for forces property")

  if (.not. get_value(at%params, 'Dipole_Moment', dipole_mom)) dipole_mom = 0.0_dp

  mainlog%prefix = "SPECSYMFORCE"
  call print("0.0 " // config_n)
  do i=1, at%N
    call print (at%pos(:,i)/BOHR // " " // at%Z(i) // " " // at%Z(i))
    call print (forces(:,i)/(HARTREE/BOHR))
  end do

  call print (dipole_mom/BOHR)
  call print ("0.0 0.0 0.0")

  mainlog%prefix = ""

  call system_finalise()

end program
