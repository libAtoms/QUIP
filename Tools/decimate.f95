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
program decimate
use libatoms_module
implicit none
  type(Atoms) :: at
  integer :: i, n, stat
  type(dictionary) :: cli_params
  character(len=FIELD_LENGTH) :: filename
  type(inoutput) :: xyzfile

  call system_initialise(verbosity=PRINT_SILENT)
  call verbosity_push(PRINT_NORMAL)

  call initialise(cli_params)
  call param_register(cli_params,"n",PARAM_MANDATORY, n)
  call param_register(cli_params,"file","stdin", filename)
  if (.not. param_read_args(cli_params, do_check = .true.)) then
    call system_abort("Usage: decimate n=(1) [file=(stdin)]")
  endif
  call finalise(cli_params)

  call initialise(xyzfile, filename, INPUT)

  call read_xyz(at, xyzfile, status=stat)
  call print_xyz(at, mainlog, all_properties=.true.)
  do while (stat == 0)
    do i=1, n-1
      call read_xyz(xyzfile, status=stat)
      if (stat /= 0) exit
    end do
    if (stat == 0) then
      call read_xyz(at, xyzfile, status=stat)
      if (stat == 0) call print_xyz(at, mainlog, all_properties=.true.)
    endif
  end do

  call verbosity_pop()
  call system_finalise()
end program decimate
