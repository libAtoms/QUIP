program test_potential
use libAtoms_module
use Potential_module
use atoms_minimisation_module

implicit none

  type(Potential), target :: pot
  type(Potential)         :: hybridpot
  type(inoutput) :: params, in
  type(Atoms)    :: at, dia
  type(Table)    :: qmlist
  real(dp) :: E0,E1
  real(dp), allocatable :: F0(:,:), F1(:,:)
  integer :: i
  logical :: status
  integer, pointer :: hybrid(:)
  
  call system_initialise()

  call Initialise(params, "quip_params.xml")
  call Initialise(pot, 'IP SW', params)
  call Initialise(hybridpot, 'HYBRID')
  hybridpot%hybrid%pot_region1 => pot
  hybridpot%hybrid%pot_region2 => pot

  call diamond(dia, 5.44_dp)
  call supercell(at, dia, 4,4,4)
  call randomise(at%pos, 0.3_dp)

  at%Z = 14
  call set_cutoff(at, cutoff(pot))

  call calc_connect(at)

  allocate(F0(3,at%N))

  call calc(pot, at, e = E0, f = F0)
  call print('E0 '//E0) 

  call add_property(at, 'hybrid', HYBRID_NO_MARK)
  status = assign_pointer(at, 'hybrid', hybrid)
  call bfs_grow(at, qmlist, 1, 2)
  hybrid(int_part(qmlist,1)) = HYBRID_REGION1_MARK

  call create_local_energy_weights(at, 2, 1)

  allocate(F1(3,at%N))
  call calc(hybridpot, at, e = E1, f = F1)
  call print('E1 '//E1) 

  call print('E1-E0 '//(E1-E0))
  call print (transpose(F1-F0))

  status = test_gradient(at, pot, NERD, do_forces = .true.)

  call system_finalise()

end program
