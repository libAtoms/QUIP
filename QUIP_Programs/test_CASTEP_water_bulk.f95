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

program test_CASTEP_water_chain
use libatoms_module
use potential_module
implicit none
  type(Atoms) :: cluster_qm, cluster_mm, bulk
  type(Inoutput) :: io
  character(len=FIELD_LENGTH) :: init_args, calc_args, infile
  real(dp) :: g_width_O, g_width_H, vacuum, charge_scale
  integer :: qm_center_i, qm_center_i_cluster
  real(dp) :: r_qm, r
  type(Potential) :: pot
  type(Dictionary) :: cli
  real(dp), allocatable :: f(:,:)
  integer :: i, ii, j
  real(dp) :: charge(128)
  real(dp) :: lat(3,3)

  call system_initialise(seed=1)

  call initialise(cli)
  init_args = ""
  call param_register(cli, "init_args", param_mandatory, init_args, help_string="No help yet.  This source file was $LastChangedBy$")
  calc_args = ""
  call param_register(cli, "calc_args", "", calc_args, help_string="No help yet.  This source file was $LastChangedBy$")
  infile = ""
  call param_register(cli, "infile", "stdin", infile, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "qm_center_i", param_mandatory, qm_center_i, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "r_qm", param_mandatory, r_qm, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "charge_scale", "1.0", charge_scale, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "g_width_O", ""//(1.2_dp/sqrt(2.0_dp)), g_width_O, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "g_width_H", ""//(0.44_dp/sqrt(2.0_dp)), g_width_H, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "vacuum", "15.0", vacuum, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli, ignore_unknown=.false., do_check=.true.)) &
    call system_abort("Failed to parse CLI parameters")

  call print("qm_center_i " // qm_center_i)
  call print("r_qm " // r_qm)
  call print("charge_scale " // charge_scale)
  call print("g_width_O " // g_width_O)
  call print("g_width_H " // g_width_H)
  call print("vacuum " // vacuum)

  charge = 0.0_dp
  charge(1) = 0.417_dp*charge_scale
  charge(8) = -0.834_dp*charge_scale

  call finalise(cli)

  call read(bulk, trim(infile))

  call initialise(pot, init_args)
  call initialise(cluster_qm,0,bulk%lattice)
  call initialise(cluster_mm,0,bulk%lattice)

  call calc_connect(bulk)

  call print("past calc_connect(bulk)", PRINT_ALWAYS)

  ! create qm and mm lists
  do i=1, bulk%N
    if (mod(i,1000) == 1) call print("creating qm and mm lists for atom " // i, PRINT_ALWAYS)
    if (bulk%Z(i) == 8) then
      r = distance_min_image(bulk, i, qm_center_i)
      if (r <= r_qm) then
        call add_atoms(cluster_qm, bulk%pos(:,i), bulk%Z(i))
	if (i == qm_center_i) qm_center_i_cluster = cluster_qm%N
	call print("atom " // i // " n_neigh " // atoms_n_neighbours(bulk, i))
        do ii=1, atoms_n_neighbours(bulk, i)
          j = atoms_neighbour(bulk, i, ii)
          call add_atoms(cluster_qm, bulk%pos(:,j), bulk%Z(j))
        end do
      else
        call add_atoms(cluster_mm, bulk%pos(:,i), bulk%Z(i))
        do ii=1, atoms_n_neighbours(bulk, i)
          j = atoms_neighbour(bulk, i, ii)
          call add_atoms(cluster_mm, bulk%pos(:,j), bulk%Z(j))
        end do
      endif
    endif
  end do

  call print("past creating cluster_mm and cluster_qm", PRINT_ALWAYS)

  ! move qm_center to origin
  do i=1, cluster_qm%N
    cluster_qm%pos(:,i) = cluster_qm%pos(:,i) - bulk%pos(:,qm_center_i)
  end do
  do i=1, cluster_mm%N
    cluster_mm%pos(:,i) = cluster_mm%pos(:,i) - bulk%pos(:,qm_center_i)
  end do

  call print("past shifting positions cluster_mm and cluster_qm", PRINT_ALWAYS)

  ! map back into cell with original (bulk) lattice
  call map_into_cell(cluster_mm)
  call map_into_cell(cluster_qm)

  call print("past map_into_cell positions cluster_mm and cluster_qm", PRINT_ALWAYS)

  call calc_connect(cluster_mm)
  call print("past calc_connect cluster_mm", PRINT_ALWAYS)
  ! make sure MM molecules aren't broken up by PBCs
  do i=1, cluster_mm%N, 3
    if (mod(i,1000) == 1) call print("coalescing molecule for atom " // i, PRINT_ALWAYS)
    call coalesce_in_one_periodic_image(cluster_mm, i)
  end do

  ! calculate qm lattice, and center around middle of qm cell
  ! calulate pbcs for qm region
  lat = 0.0_dp
  lat(1,1) = (maxval(cluster_qm%pos(1,:))-minval(cluster_qm%pos(1,:)))+vacuum
  lat(2,2) = (maxval(cluster_qm%pos(2,:))-minval(cluster_qm%pos(2,:)))+vacuum
  lat(3,3) = (maxval(cluster_qm%pos(3,:))-minval(cluster_qm%pos(3,:)))+vacuum
  call set_lattice(cluster_qm, lat, scale_positions=.false.)

  call print("lat", PRINT_ALWAYS)
  call print(lat, PRINT_ALWAYS)

  do i=1, cluster_qm%N
    cluster_qm%pos(:,i) = cluster_qm%pos(:,i) + sum(lat,2)/2.0_dp
  end do
  do i=1, cluster_mm%N
    cluster_mm%pos(:,i) = cluster_mm%pos(:,i) + sum(lat,2)/2.0_dp
  end do

  call finalise(bulk)

  call initialise(io, "extcharges", OUTPUT)
  do i=1, cluster_mm%N
    if (cluster_mm%Z(i) == 1) then
      call print(charge(cluster_mm%Z(i)) //" " // cluster_mm%pos(:,i) //" " // g_width_H, file=io)
    else if (cluster_mm%Z(i) == 8) then
      call print(charge(cluster_mm%Z(i)) //" " // cluster_mm%pos(:,i) //" " // g_width_O, file=io)
    else
      call system_abort("No width defined for Z="//cluster_mm%Z(i))
    endif
  end do
  call finalise(io)

  call print("BOB "//(cluster_qm%N+cluster_mm%N))
  call print("BOB ")
  do i=1, cluster_qm%N
    call print("BOB "//cluster_qm%Z(i)//" "//cluster_qm%pos(:,i))
  end do
  do i=1, cluster_mm%N
    call print("BOB " //(cluster_mm%Z(i)+1)//" "//cluster_mm%pos(:,i))
  end do

  allocate(f(3,cluster_qm%N))
  f = 0.0_dp
  call calc(pot, cluster_qm, f=f, args_str=calc_args)
 
  do i=1, cluster_qm%N
    if (i == qm_center_i_cluster) then
      call print(" " // r_qm // " " // i // " " // cluster_qm%pos(:,i) // " " // f(:,i) // " CTR")
    else
      call print(" " // r_qm // " " // i // " " // cluster_qm%pos(:,i) // " " // f(:,i))
    endif
  end do
  
  call finalise(pot)
  call system_finalise()

contains

  subroutine add_water(at, p)
    type(Atoms), intent(inout) :: at
    real(dp), intent(in) :: p(3)

    real(dp) :: po(3), ph1(3), ph2(3)
    real(dp), parameter :: r_oh = 0.95_dp, theta_hoh=105.0_dp*PI/180.0_dp
    real(dp) :: v1(3), v2(3)

    v1 = random_unit_vector()
    v2 = random_unit_vector()
    v2 = v2 - (v1.dot.v2)*v1
    v2 = v2 / norm(v2)

    po = -v1*r_oh*cos(theta_hoh/2.0_dp)/3.0_dp
    ph1 = po + v1*r_oh*cos(theta_hoh/2.0_dp) + v2*r_oh*sin(theta_hoh/2.0_dp)
    ph2 = po + v1*r_oh*cos(theta_hoh/2.0_dp) - v2*r_oh*sin(theta_hoh/2.0_dp)

    call add_atoms(at, p+po, 8)
    call add_atoms(at, p+ph1, 1)
    call add_atoms(at, p+ph2, 1)

  end subroutine add_water

end program test_CASTEP_water_chain
