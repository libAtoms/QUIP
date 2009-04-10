program est_n_poles
use libatoms_module
use quip_module
implicit none

  type(ApproxFermi) :: af
  character(len=1024) :: arg
  real(dp) :: fermi_T, band_width

  call system_initialise()

  if (cmd_arg_count() /= 2) then
    call system_abort("Usage: est_n_poles fermi_T band_width")
  endif

  call get_cmd_arg(1, arg)
  fermi_T = string_to_real(arg)
  call get_cmd_arg(2, arg)
  band_width = string_to_real(arg)

  call Initialise(af, 0.0_dp, fermi_T, band_width)

  call verbosity_push(VERBOSE)
  call print(af)
  call verbosity_pop()

  call system_finalise()

end program
