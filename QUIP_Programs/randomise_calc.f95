
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


