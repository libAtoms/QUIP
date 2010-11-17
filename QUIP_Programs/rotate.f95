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

program do_rotate
use libatoms_module
use frametools_module
implicit none
  type(atoms) :: at
  integer stat
  real(dp) :: axis(3), theta, origin(3)
  type(dictionary) :: cli_params

  call system_initialise(PRINT_SILENT)

  call initialise(cli_params)
  call param_register(cli_params, "axis", PARAM_MANDATORY, axis, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "angle", PARAM_MANDATORY, theta, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "origin", "0 0 0", origin, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli_params)) then
    call print('Usage: rotate axis="x y z" angle=theta [ origin="0 0 0" ]', PRINT_ALWAYS)
    call system_abort("Failed to parse cli args")
  endif
  call finalise(cli_params)

  call read_xyz(at, "stdin", status=stat)
  do while (stat == 0) 
    call rotate(at%pos, axis, theta, origin)
    call verbosity_push(PRINT_NORMAL)
    call print_xyz(at, mainlog, all_properties=.true., real_format="f15.10")
    call verbosity_pop()
    call read_xyz(at, "stdin", status=stat)
  end do

  call system_finalise()
end program
