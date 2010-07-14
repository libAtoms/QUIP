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

program cxyzreader
use libatoms_module
implicit none
  type(Atoms) :: prim, at
  type(CInoutput) :: io
  integer i

  call system_initialise()

  call system_command("rm -f test.xyz test.xyz.idx")
  call print("writing test.xyz", PRINT_ALWAYS)
  call initialise(io,"test.xyz", OUTPUT)
  call diamond(prim, 5.43_dp, (/ 14 /) )
  call supercell(at, prim, 3, 3, 3)
  do i=1, 200
    call write(at, io)
    at%pos(1,1) = at%pos(1,1) + 0.001
  end do
  call finalise(io)
  call finalise(at)

  call print("reading test.xyz", PRINT_ALWAYS)
  call initialise(io,"test.xyz", INPUT)
  call read(io, at)
  call print("first read at%pos(1,1) " // at%pos(1,1))
  call read(io, at)
  call print("second read at%pos(1,1) " // at%pos(1,1))
  call finalise(io)
  call finalise(at)

  call initialise(io,"test.xyz", INPUT)
  call read(io, at, frame=50)
  call print("frame=50 read at%pos(1,1) " // at%pos(1,1))
  call read(io, at, frame=25)
  call print("frame=25 read at%pos(1,1) " // at%pos(1,1))
  call finalise(io)
  call finalise(at)

  call system_finalise()
end program cxyzreader
