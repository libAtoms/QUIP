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

program test
use libatoms_module
use potential_module
use vacancy_map_module

implicit none
  type(Atoms) at

  type(Dictionary) :: cli_params
  real(dp) :: sphere_qm_center(3), sphere_qm_rad
  real(dp) :: cyl_qm_center(3), cyl_qm_vec(3), cyl_qm_rad
  integer :: qm_hops, buffer_hops, transition_hops
  type(Table) :: core_list
  real(dp)  :: r2, dr(3)
  integer, pointer :: hybrid_mark(:)
  integer :: i

  call system_initialise(enable_timing=.true.,verbosity=PRINT_SILENT)

  call initialise(cli_params)
  call param_register(cli_params, "buffer_hops", "1", buffer_hops, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "transition_hops", "1", transition_hops, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "qm_hops", "0", qm_hops, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "sphere_qm_center", "0 0 0", sphere_qm_center, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "sphere_qm_rad", "0", sphere_qm_rad, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "cyl_qm_center", "0 0 0", cyl_qm_center, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "cyl_qm_vec", "0 0 1", cyl_qm_vec, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "cyl_qm_rad", "0", cyl_qm_rad, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli_params, do_check = .true.)) then
    call print("Usage: vacancy_map_hybrid_generate [buffer_hops=i(1)] [transition_hops=i(1)]", PRINT_ALWAYS)
    call print("   [qm_hops=i(1)] sphere_qm_center={x y z}(0 0 0) sphere_qm_rad=r(1)", PRINT_ALWAYS)
    call print("   cyl_qm_center={x y z}(0 0 0) cyl_qm_vec={x y z}(0 0 0) cyl_qm_rad=r(1)", PRINT_ALWAYS)
    call system_abort("Confused by CLI arguments")
  end if
  call finalise(cli_params)

  call print("Using buffer_hops " // buffer_hops // " transition_hops " // transition_hops // " qm_hops " // qm_hops)
  call print("Using sphere_qm_center " // sphere_qm_center // " qm_rad " // sphere_qm_rad)
  call print("Using cyl_qm_center " // cyl_qm_center // " cyl_qm_vec " // cyl_qm_vec // " cyl_qm_rad " // cyl_qm_rad)

  ! read initial config
  call read_xyz(at, "stdin")
  call set_cutoff(at, 3.2_dp)
  call calc_connect(at)

  call add_property(at, 'hybrid_mark', 0)
  call add_property(at, 'weight_region1', 0.0_dp)

  if (.not. assign_pointer(at, 'hybrid_mark', hybrid_mark)) &
    call system_abort("Couldn't find hybrid_mark in at")

  ! core qm region
  call wipe(core_list)
  cyl_qm_vec = cyl_qm_vec / norm(cyl_qm_vec)
  do i=1,at%N
    if (sphere_qm_rad > 0) then
      r2 = sum(((at%pos(:,i)-sphere_qm_center)/sphere_qm_rad)**2)
      if (r2 < 1.0_dp) call append(core_list, (/i, 0, 0 ,0/))
    endif
    if (cyl_qm_rad > 0) then
      dr = at%pos(:,i) - cyl_qm_center
      dr = dr - cyl_qm_vec*(cyl_qm_vec .dot. dr)
      r2 = sum((dr/cyl_qm_rad)**2)
      if (r2 < 1.0_dp) call append(core_list, (/i, 0, 0 ,0/))
    endif
  end do

  if (qm_hops > 0) call bfs_grow(at, core_list, qm_hops, nneighb_only = .false.)

  call print(core_list, PRINT_VERBOSE)
   
  hybrid_mark = HYBRID_NO_MARK
  hybrid_mark(int_part(core_list, 1)) = HYBRID_ACTIVE_MARK

  call create_local_energy_weights(at, transition_hops, buffer_hops) ! trans_width, buffer_width

  call verbosity_push(PRINT_NORMAL)
  call print_xyz(at, mainlog, properties="pos:Z:hybrid_mark:weight_region1")
  call verbosity_pop()

  call system_finalise()
end program
