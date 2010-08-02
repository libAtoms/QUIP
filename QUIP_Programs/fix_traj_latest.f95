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

program fix_traj
use libatoms_module
  type(atoms) :: at
  integer, pointer :: cluster_mark(:), t_i(:)
  character(len=1), pointer :: t_s(:,:)
  real(dp), pointer :: t_d(:)

  call system_initialise()

  call read(at, "stdin")

  call add_property(at, 'Z', 0)
  call add_property(at, 'mass', 0.0_dp)
  call add_property(at, 'travel', 0, 3)
  call add_property(at, 'move_mask', 1)
  call add_property(at, 'damp_mask', 0)
  call add_property(at, 'thermostat_region', 0)
  call add_property(at, 'avg_ke', 0.0_dp)
  call add_property(at, 'acc', 0.0_dp, 3)
  call add_property(at, 'avgpos', 0.0_dp, 3)
  call add_property(at, 'oldpos', 0.0_dp, 3)
  call add_property(at, 'pot', 0.0_dp)
  call add_property(at, 'hybrid', 1)
  call add_property(at, 'hybrid_mark', 1)
  call add_property(at, 'old_hybrid_mark', 1)
  call add_property(at, 'old_cluster_mark', 1)
  call add_property(at, 'cut_bonds', 0, 4)
  call add_property(at, 'atom_type', '')
  call add_property(at, 'atom_res_name', '')
  call add_property(at, 'atom_mol_name', '')
  call add_property(at, 'atom_res_number', 0)
  call add_property(at, 'atom_charge', 0.0_dp)
  call add_property(at, 'weight_region1', 0.0_dp)
  call add_property(at, 'qm_force', 0.0_dp, 3)
  call add_property(at, 'mm_force', 0.0_dp, 3)
  call add_property(at, 'force', 0.0_dp, 3)


  if (.not. assign_pointer(at, 'cluster_mark', cluster_mark)) call system_abort("no cluster_mark")

  if (.not. assign_pointer(at, 'Z', t_i)) call system_abort("no Z")
  if (.not. assign_pointer(at, 'species', t_s)) call system_abort("no species")
  do i=1, at%N
    t_i(i) = atomic_number_from_symbol(trim(a2s(t_s(:,i))))
  end do

  if (.not. assign_pointer(at, 'mass', t_d)) call system_abort("no mass")
  do i=1, at%N
    t_d(i) = ElementMass(at%Z(i))
  end do

  if (.not. assign_pointer(at, 'hybrid', t_i)) call system_abort("no hybrid")
  where (cluster_mark == HYBRID_ACTIVE_MARK)
    t_i = 1
  else where
    t_i = 0
  end where
  if (.not. assign_pointer(at, 'hybrid_mark', t_i)) call system_abort("no hybrid_mark")
  t_i = cluster_mark
  if (.not. assign_pointer(at, 'old_hybrid_mark', t_i)) call system_abort("no old_hybrid_mark")
  t_i = cluster_mark
  if (.not. assign_pointer(at, 'old_cluster_mark', t_i)) call system_abort("no cluster_mark")
  t_i = cluster_mark

  call write(at, "stdout", prefix="OUT")

  call system_finalise()
end program
