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

program randomise_calc

use libAtoms_module
use Potential_module
use MetaPotential_module

implicit none

  type(Potential) pot
  type(Atoms):: at, at2
  real(dp)::e
  real(dp):: lat(3,3), v(3,3)
  real(dp), allocatable::f(:,:)
  integer :: i

  call system_initialise()
  

  !call Initialise(pot, 'FilePot command=./castep_driver.sh')
  call Initialise(pot, 'FilePot command=./castep_driver.py')

  call read_xyz(at2, 'FeSb3La.xyz')
  allocate(f(3,at2%N))
  do i=1,20
     at=at2
     call randomise(at%pos, 0.2_dp)
     lat = at%lattice
     call randomise(lat, 0.5_dp)
     call set_lattice(at, lat)
     call calc(pot, at, e=e, f=f, virial=v)
     call add_property(at, 'f', f)
     call print_xyz(at, 'output.xyz', properties='pos:f', comment='Energy='//e//' virial={'//reshape(v, (/9/))//'}', append=.true.)
  end do
  deallocate(f)

  call system_finalise()


end program randomise_calc


