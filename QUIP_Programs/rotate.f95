program do_rotate
use libatoms_module
use frametools_module
implicit none
  type(atoms) :: at
  integer stat
  real(dp) :: axis(3), theta, origin(3)
  type(dictionary) :: cli_params

  call system_initialise(SILENT)

  call initialise(cli_params)
  call param_register(cli_params, "axis", PARAM_MANDATORY, axis)
  call param_register(cli_params, "angle", PARAM_MANDATORY, theta)
  call param_register(cli_params, "origin", "0 0 0", origin)
  if (.not. param_read_args(cli_params, do_check=.true.)) then
    call print('Usage: rotate axis="x y z" angle=theta [ origin="0 0 0" ]', ERROR)
    call system_abort("Failed to parse cli args")
  endif
  call finalise(cli_params)

  call read_xyz(at, "stdin", status=stat)
  do while (stat == 0) 
    call rotate(at%pos, axis, theta, origin)
    call verbosity_push(NORMAL)
    call print_xyz(at, mainlog, all_properties=.true., real_format="f15.10")
    call verbosity_pop()
    call read_xyz(at, "stdin", status=stat)
  end do

  call system_finalise()
end program
