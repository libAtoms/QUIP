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

! Reads a sequence of configs and writes them back, mostly to convert from xyz to netcdf and back

#include "error.inc"

program file_rewrite
use libatoms_module
implicit none

  type(Dictionary) :: cli_params
  character(len=FIELD_LENGTH) :: infilename
  character(len=FIELD_LENGTH) :: outfilename
  logical :: netcdf4
  type(CInOutput) :: infile, outfile
  type(Atoms) :: at
  integer i
  integer :: error = ERROR_NONE

  call system_initialise(PRINT_NORMAL)

  call initialise(cli_params)
  call param_register(cli_params, 'infile', 'stdin', infilename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'outfile', 'stdout', outfilename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'netcdf4', 'F', netcdf4, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli_params)) then
      call print_usage()
    call system_abort('could not parse argument line')
  end if
  call finalise(cli_params)

  call print("infile " // trim(infilename) // " outfile " // trim(outfilename) // " netcdf4 " // netcdf4)

  call initialise(infile, infilename, action=INPUT, no_compute_index=.true.)
  call initialise(outfile, outfilename, action=OUTPUT, netcdf4=netcdf4)

  i = 1
  call read(infile, at, error=error)
  do while (error == ERROR_NONE) 
    call write(outfile, at, error=error)
    HANDLE_ERROR(error)
    if (mod(i,100) == 0) write (*,'(I1,$)') mod(i/100,10)
    call read(infile, at, error=error)
    HANDLE_ERROR(error)
    i = i + 1
  end do
  ! We do not want to handle this error
  call clear_error(error)

  call finalise(outfile)
  call finalise(infile)
  call system_finalise()

contains

  subroutine print_usage()

    character(len=1024) my_exec_name

    if (EXEC_NAME == "<UNKNOWN>") then
      my_exec_name=""
    else
      my_exec_name=EXEC_NAME
    endif

    call print("Usage: " // trim(my_exec_name)//" infile=filename(stdin) outfile=filename(stdout)[=file.nc for NETCDF] [netcdf4 (for output)]", PRINT_ALWAYS)
  end subroutine print_usage

end program file_rewrite
