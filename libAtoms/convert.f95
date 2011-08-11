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
#include "error.inc"

program convert

use libAtoms_module
use libatoms_misc_utils_module

implicit none

  type(Atoms) at
  type(CInOutput) :: infile, outfile
  character(len=FIELD_LENGTH) :: infilename, outfilename
  integer error

  call system_initialise(verbosity=PRINT_SILENT)

  if (NUM_COMMAND_ARGS == 2) then
     infilename = COMMAND_ARG(1)
     outfilename = COMMAND_ARG(2)
  else if(NUM_COMMAND_ARGS == 0) then
     infilename = 'stdin'
     outfilename = 'stdout'
  else
     call system_abort("Usage: convert [ infile.{xyz|nc} outfile.{xyz|nc} ]")
  end if

  call initialise(infile, trim(infilename))
  call initialise(outfile, trim(outfilename), action=OUTPUT)
     
  do
     call read(at, infile, error=error)
     if (error /= 0) then
        if (error == ERROR_IO_EOF) then
	   exit
	else
	   HANDLE_ERROR(error)
	endif
     endif
     call print(at)
     call write(at, outfile)
  end do
  call system_finalise()
end program convert
