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

! print 1 out of every n frames
#include "error.inc"
program decimate
use libatoms_module
implicit none
  type(Atoms) :: at
  integer :: i, n, error
  type(dictionary) :: cli_params
  character(len=FIELD_LENGTH) :: filename
  type(cinoutput) :: xyzfile

  call system_initialise(verbosity=PRINT_SILENT)
  call verbosity_push(PRINT_NORMAL)

  call initialise(cli_params)
  call param_register(cli_params,"n",PARAM_MANDATORY, n, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"file","stdin", filename, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli_params, do_check = .true.)) then
    call system_abort("Usage: decimate n=(1) [file=(stdin)]")
  endif
  call finalise(cli_params)

  call initialise(xyzfile, filename, INPUT)

  call read(at, xyzfile, error=error)
  HANDLE_ERROR(error)
  call write(at, "stdout")
  i = 1
  do while (error == 0)
    i = i + n
    call read(at, xyzfile, frame=i, error=error)
    if (error /= 0) then
       if (error == ERROR_IO_EOF) exit
       HANDLE_ERROR(error)
    endif
    call write(at, "stdout")
  end do

  call verbosity_pop()
  call system_finalise()
end program decimate
