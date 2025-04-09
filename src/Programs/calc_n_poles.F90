! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   Portions of this code were written by
! H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! H0 X
! H0 X   Copyright 2006-2010.
! H0 X
! H0 X   These portions of the source code are released under the GNU General
! H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
! H0 X
! H0 X   If you would like to license the source code under different terms,
! H0 X   please contact Gabor Csanyi, gabor@csanyi.net
! H0 X
! H0 X   Portions of this code were written by Noam Bernstein as part of
! H0 X   his employment for the U.S. Government, and are not subject
! H0 X   to copyright in the USA.
! H0 X
! H0 X
! H0 X   When using this software, please cite the following reference:
! H0 X
! H0 X   http://www.libatoms.org
! H0 X
! H0 X  Additional contributions by
! H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! H0 X
! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

program est_n_poles
use libatoms_module
use approxfermi_module
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

  call verbosity_push(PRINT_VERBOSE)
  call print(af)
  call verbosity_pop()

  call system_finalise()

end program
