! FrameTools like wrapper for coalesce_in_one_periodic_image()
! moves all atoms so that no bonds occur across a periodic image
program wrap_prog
use libatoms_module
implicit none
  type(Atoms) :: at
  type(Dictionary) :: cli_params
  integer :: seed
  real(dp) :: cutoff_factor

  call system_initialise()

  call initialise(cli_params)
  call param_register(cli_params, 'seed', '1', seed)
  call param_register(cli_params, 'cutoff_factor', '1.0', cutoff_factor)
  if (.not. param_read_args(cli_params, do_check = .true.)) then
    call print("Usage: wrap [seed=1] [cutoff_factor=1.0]", ERROR)
    call system_abort("Confused by CLI parameters")
  endif
  call finalise(cli_params)

  call read_xyz(at, "stdin")

  call set_cutoff_factor(at, cutoff_factor)
  call calc_connect(at)

  call coalesce_in_one_periodic_image(at, seed)

  mainlog%prefix="WRAPPED"
  call print_xyz(at, mainlog, properties="pos")
  mainlog%prefix=""

  call system_finalise()
end program

