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

program test_potential
use libAtoms_module
use Potential_module

implicit none

  type(Potential) pot
  type(inoutput) params
  type(Atoms) at

  real(dp) :: E0
  real(dp), allocatable :: F0(:,:)


  integer i

  call system_initialise()

  call Initialise(params, "quip_params.xml")

  call read(at, 'stdin')

  call Initialise(pot, 'IP LJ', params)

  call set_cutoff(at, cutoff(pot)+0.2)

  call calc_connect(at)

  call print(at)

  allocate(F0(3,at%N))
  call calc(pot, at, e = E0, f = F0)

  call print('E0 '//E0)
 
  do i=1, at%N
     call print('F '//i//at%pos(:,i)//F0(:,i))
  end do

  call system_finalise()

end program
