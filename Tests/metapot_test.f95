program metapot_test

  use libAtoms_module
  use QUIP_Module

  implicit None

  type(Atoms) :: dia, at
  type(Table) :: embedlist
  type(MetaPotential) :: metapot
  type(InOutput) :: xml
  type(Potential) :: pot1, pot2
  real(dp), allocatable :: f1(:,:), f2(:,:)
  integer, pointer :: hybrid(:), hybrid_mark(:)

  call system_initialise(VERBOSE,seed=2)
  call initialise(xml, 'metapot_test_params.xml')
  call initialise(pot1, 'IP SW label="PRB_31_plus_H"', xml)
  call rewind(xml)
  call initialise(pot2, 'IP SW label="eps_2.6"', xml)
  
  call diamond(dia, 5.44_dp, (/14/))
  call supercell(at, dia, 4, 4, 4)
  call randomise(at%pos, 0.1_dp)
  
  allocate(f1(3,at%N),f2(3,at%N))

  call set_cutoff(at, cutoff(pot1)+2.0_dp)
  call calc_connect(at)

  ! Mark some atoms for embedding
  call add_property(at, 'hybrid', 0)
  if (.not. assign_pointer(at, 'hybrid', hybrid)) call system_abort('Cannot assign hybrid pointer')
  call append(embedlist, (/1,0,0,0/))
  call bfs_grow(at, embedlist, 2, nneighb_only=.true.)
  hybrid(int_part(embedlist,1)) = 1

  call print('embed list')
  call print(embedlist)

  ! Test embedding methods in Potential calc()
  call add_property(at, 'hybrid_mark', HYBRID_NO_MARK)
  if (.not. assign_pointer(at, 'hybrid_mark', hybrid_mark)) call system_abort('Cannot assign hybrid_mark pointer')
  where (hybrid /= 0) hybrid_mark = HYBRID_ACTIVE_MARK

  call print('single_cluster')
  call create_hybrid_weights(at, 'buffer_hops=3')
  call calc(pot2, at, f=f1, args_str='single_cluster=T cluster_calc_connect=T')
  call print(f1)

  call print('little_clusters')
  call calc(pot2, at, f=f2, args_str='little_clusters=T cluster_calc_connect=T buffer_hops=3')
  call print(f2)

  ! buffering is perfect, forces should be equal
  call print('embedding test force error '//rms_diff(f1, f2)//' '//maxval(abs(f1-f2)))

  ! Now test all the (non-LOTF) force mixing methods

  call print('force_mixing')
  call initialise(metapot, 'ForceMixing method=force_mixing buffer_hops=3 qm_args_str={single_cluster=T cluster_calc_connect=T}', pot1, pot2)
  call print(metapot)
  call calc(metapot, at, f=f1)
  call print(f1)
  call finalise(metapot)

  call print('force_mixing_smooth')
  call initialise(metapot, 'ForceMixing method=force_mixing_smooth buffer_hops=1 qm_args_str={single_cluster=T  cluster_calc_connect=T}', pot1, pot2)
  call print(metapot)
  call calc(metapot, at, f=f2)
  call print(f2)
  call print('force_mixing_smooth force error '//rms_diff(f1, f2)//' '//maxval(abs(f1-f2)))
  call finalise(metapot)

  call print('force_mixing_super_smooth')
  call initialise(metapot, 'ForceMixing method=force_mixing_super_smooth buffer_hops=1 qm_args_str={single_cluster=T cluster_calc_connect=T}', pot1, pot2)
  call print(metapot)
  call calc(metapot, at, f=f2)
  call print(f2)
  call print('force_mixing_super_smooth force error '//rms_diff(f1, f2)//' '//maxval(abs(f1-f2)))
  call finalise(metapot)

  call print('conserve_momentum')
  call initialise(metapot, 'ForceMixing method=conserve_momentum fit_hops=2 buffer_hops=1 qm_args_str={single_cluster=T cluster_calc_connect=T}', pot1, pot2)
  call print(metapot)
  call calc(metapot, at, f=f2)
  call print(f2)
  call print('conserve_momentum force error '//rms_diff(f1, f2)//' '//maxval(abs(f1-f2)))
  call finalise(metapot)


  deallocate(f1, f2)
  call finalise(xml)
  call finalise(pot1)
  call finalise(pot2)
  call finalise(at)
  call finalise(dia)
  call system_finalise

end program metapot_test
