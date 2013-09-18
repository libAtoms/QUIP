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

program order

use libAtoms_module
use libatoms_misc_utils_module

implicit none

  type(Atoms) at, at2
  logical, dimension(:), allocatable :: atomdone
  integer error, i, j, k, jj

  
  call system_initialise(verbosity=PRINT_SILENT)

  do
     call read(at, 'stdin', error=error)
     if (error /= 0) then
        if (error == ERROR_IO_EOF) then
	   exit
	else
	   HANDLE_ERROR(error)
	endif
     endif

     call set_cutoff(at, 0.0_dp)
     call calc_connect(at)

     call initialise(at2, at%N, at%lattice)
     allocate(atomdone(at%N))
     atomdone = .false.
     k = 1
     do i=1,at%N
        if(atomdone(i) .eqv. .true.) cycle
        at2%species(:,k) = at%species(:,i)
        at2%pos(:,k) = at%pos(:,i)
        k = k+1
        atomdone(i) = .true.
        do j=1,n_neighbours(at, i)
           jj = neighbour(at, i, j)
           at2%species(:,k) = at%species(:,jj)
           at2%pos(:,k) = at%pos(:,jj)
           k = k+1
           atomdone(jj) = .true.
        end do
     end do

     call verbosity_push(PRINT_NORMAL)
     call write(at2, 'stdout')
     call verbosity_pop()

  enddo
  call system_finalise()

end program order
