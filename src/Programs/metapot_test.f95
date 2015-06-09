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

program pot_test

  use libAtoms_module
  use Potential_Module

  implicit None

  type(Atoms) :: dia, at
  type(Table) :: embedlist
  type(Potential) :: pot
  type(InOutput) :: xml
  type(Potential) :: pot1, pot2
  real(dp), allocatable :: f1(:,:), f2(:,:)
  integer, pointer :: hybrid(:), hybrid_mark(:)

  call system_initialise(PRINT_VERBOSE,seed=2)
  call initialise(xml, 'pot_test_params.xml')
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
  call set_value(at%params, 'mode', 'embedding')
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
  if (.not. assign_pointer(at, 'hybrid', hybrid)) call system_abort('Cannot assign hybrid pointer')
  where (hybrid /= 0) hybrid_mark = HYBRID_ACTIVE_MARK

  call print('single_cluster')
  call create_hybrid_weights(at, 'buffer_hops=3')
  call calc(pot2, at, force=f1, args_str='single_cluster=T cluster_calc_connect=T')
  call print(f1)

  call print('little_clusters')
  call calc(pot2, at, force=f2, args_str='little_clusters=T cluster_calc_connect=T buffer_hops=3')
  call print(f2)

  ! buffering is perfect, forces should be equal
  call print('embedding test force error '//rms_diff(f1, f2)//' '//maxval(abs(f1-f2)))

  ! Now test all the (non-LOTF) force mixing methods

  call set_value(at%params, 'mode', 'force_mixing')
  call print('force_mixing')
  call initialise(pot, 'ForceMixing method=force_mixing buffer_hops=3 qm_args_str={single_cluster=T cluster_calc_connect=T}', pot1, pot2)
  call print(pot)
  call calc(pot, at, force=f1)
  call print(f1)
  call finalise(pot)

  call set_value(at%params, 'mode', 'force_mixing_smooth')
  call print('force_mixing_smooth')
  call initialise(pot, 'ForceMixing method=force_mixing_smooth buffer_hops=1 qm_args_str={single_cluster=T  cluster_calc_connect=T}', pot1, pot2)
  call print(pot)
  call calc(pot, at, force=f2)
  call print(f2)
  call print('force_mixing_smooth force error '//rms_diff(f1, f2)//' '//maxval(abs(f1-f2)))
  call finalise(pot)

  call set_value(at%params, 'mode', 'force_mixing_super_smooth')
  call print('force_mixing_super_smooth')
  call initialise(pot, 'ForceMixing method=force_mixing_super_smooth buffer_hops=1 qm_args_str={single_cluster=T cluster_calc_connect=T}', pot1, pot2)
  call print(pot)
  call calc(pot, at, force=f2)
  call print(f2)
  call print('force_mixing_super_smooth force error '//rms_diff(f1, f2)//' '//maxval(abs(f1-f2)))
  call finalise(pot)

  call set_value(at%params, 'mode', 'conserve_momentum')
  call print('conserve_momentum')
  call initialise(pot, 'ForceMixing method=conserve_momentum fit_hops=2 buffer_hops=1 qm_args_str={single_cluster=T cluster_calc_connect=T}', pot1, pot2)
  call print(pot)
  call calc(pot, at, force=f2)
  call print(f2)
  call print('conserve_momentum force error '//rms_diff(f1, f2)//' '//maxval(abs(f1-f2)))
  call finalise(pot)


  deallocate(f1, f2)
  call finalise(xml)
  call finalise(pot1)
  call finalise(pot2)
  call finalise(at)
  call finalise(dia)
  call system_finalise

end program pot_test
