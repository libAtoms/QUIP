program test_potential
use libAtoms_module
use Potential_module

implicit none

  type(Potential) pot
  type(inoutput) params, in
  type(Atoms) at

  real(dp) :: E0
  real(dp), allocatable :: F0(:,:)


  integer i

  call system_initialise()

  call Initialise(params, "quip_params.xml")

  call Initialise(in, "stdin")
  call read_xyz(at, in)

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
