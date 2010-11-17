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
  type(Atoms) :: cluster_qm, cluster_mm
  type(Inoutput) :: io
  character(len=FIELD_LENGTH) :: init_args, calc_args
  real(dp) :: g_width, vacuum, lx, spacing, px, charge_scale
  integer :: n_qm
  type(Potential) :: pot
  type(Dictionary) :: cli
  real(dp), allocatable :: f(:,:)
  integer :: i, nctr
  real(dp) :: charge(128)
  real(dp) :: lat(3,3)

  call system_initialise(seed=1)

  call initialise(cli)
  init_args = ""
  call param_register(cli, "init_args", param_mandatory, init_args, help_string="No help yet.  This source file was $LastChangedBy$")
  calc_args = ""
  call param_register(cli, "calc_args", "", calc_args, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "n_qm", param_mandatory, n_qm, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "spacing", "2.0", spacing, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "charge_scale", "1.0", charge_scale, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "g_width", "0.25", g_width, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "vacuum", "15.0", vacuum, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli, ignore_unknown=.false.)) &
    call system_abort("Failed to parse CLI parameters")

  call print("n_qm " // n_qm)

  charge = 0.0_dp
  charge(1) = 0.417_dp*charge_scale
  charge(8) = -0.834_dp*charge_scale

  call finalise(cli)

  call initialise(pot, init_args)
  lx = 2.0*spacing*n_qm+vacuum
  lat(:,1) = (/ lx, 0.0_dp, 0.0_dp /)
  lat(:,2) = (/ 0.0_dp, vacuum, 0.0_dp /)
  lat(:,3) = (/ 0.0_dp, 0.0_dp, vacuum /)
  call initialise(cluster_qm,0,lat)
  lat(:,1) = (/ 2000.0_dp, 0.0_dp, 0.0_dp /)
  lat(:,2) = (/ 0.0_dp, 2000.0_dp, 0.0_dp /)
  lat(:,3) = (/ 0.0_dp, 0.0_dp, 2000.0_dp /)
  call initialise(cluster_mm,0,lat)
  do i = -50, 50
    px = i*spacing
    if (abs(i) <= n_qm) then
      call add_water(cluster_qm, (/ px+lx/2.0_dp, vacuum/2.0_dp, vacuum/2.0_dp /))
      if (i == 0) nctr = cluster_qm%N - 2
    else
      call add_water(cluster_mm, (/ px+lx/2.0_dp, vacuum/2.0_dp, vacuum/2.0_dp /))
    endif
  end do

  call initialise(io, "extcharges", OUTPUT)
  do i=1, cluster_mm%N
    call print(charge(cluster_mm%Z(i)) //" " // cluster_mm%pos(:,i) //" " // g_width, file=io)
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
    if (i >= nctr .and. i <= nctr+2) then
      call print(i // " " // cluster_qm%pos(:,i) // " " // f(:,i) // " CTR " // (i-nctr))
    else
      call print(i // " " // cluster_qm%pos(:,i) // " " // f(:,i))
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
