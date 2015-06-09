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

program test_CASTEP_MM_buffer_crack
use libatoms_module
use potential_module
implicit none
  type(Atoms) :: at, cluster_qm, cluster_mm
  type(Inoutput) :: io
  character(len=STRING_LENGTH) :: infile, init_args, calc_args
  real(dp) :: rqm, rmm, g_width
  type(Potential) :: pot
  type(Dictionary) :: cli
  type(Table) :: cluster_info
  integer, pointer :: idx(:)
  logical, pointer :: is_in_qm(:)
  real(dp), allocatable :: f(:,:)
  integer :: i
  real(dp) :: charge(128)

  call system_initialise()

  call initialise(cli)
  infile = ""
  call param_register(cli, "infile", "stdin", infile, help_string="No help yet.  This source file was $LastChangedBy$")
  init_args = ""
  call param_register(cli, "init_args", param_mandatory, init_args, help_string="No help yet.  This source file was $LastChangedBy$")
  calc_args = ""
  call param_register(cli, "calc_args", "", calc_args, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "rqm", param_mandatory, rqm, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "rmm", "0.0", rmm, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "g_width", "0.25", rmm, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli, ignore_unknown=.false.)) &
    call system_abort("Failed to parse CLI parameters")

  charge = 0.0_dp
  charge(6) = -1.25_dp
  charge(14) = 1.25_dp

  call finalise(cli)

  call read(at, infile)

  call initialise(pot, init_args)

  call create_hybrid_weights(at, "hysteretic_buffer hysteretic_buffer_inner_radius="//(rqm+rmm)//" hysteretic_buffer_outer_radius="//(rqm+rmm))
  cluster_info = create_cluster_info_from_hybrid_mark(at, "cluster_periodic_z")
  call carve_cluster(at, "randomise_buffer=F cluster_periodic_z" , cluster_info, cluster_mm)

  call create_hybrid_weights(cluster_mm, "hysteretic_buffer hysteretic_buffer_inner_radius="//rqm//" hysteretic_buffer_outer_radius="//rqm)
  cluster_info = create_cluster_info_from_hybrid_mark(cluster_mm, "cluster_periodic_z")
  call carve_cluster(cluster_mm, "randomise_buffer=F cluster_periodic_z" , cluster_info, cluster_qm)
  call finalise(cluster_info)

  call add_property(cluster_mm, "is_in_qm", .false.)
  if (.not. assign_pointer(cluster_mm, "is_in_qm", is_in_qm)) &
    call system_abort("Failed to assign pointer for is_in_qm in cluster_mm")
  if (.not. assign_pointer(cluster_qm, "index", idx)) &
    call system_abort("Failed to assign pointer for index in cluster_qm")

  call initialise(io, "extcharges", OUTPUT)
  do i=1, cluster_qm%N
    is_in_qm(idx(i)) = .true.
  end do
  do i=1, cluster_mm%N
    if (.not. is_in_qm(i)) then
      call print(charge(cluster_mm%Z(i))//" "//cluster_mm%pos(:,i)//" "//g_width, file=io)
    endif
  end do
  call finalise(io)

  allocate(f(3,cluster_qm%N))
  call calc(pot, cluster_qm, f=f)

  do i=1, cluster_qm%N
    call print(i// " " // cluster_qm%pos(:,i) // " " // f(:,i))
  end do

  call system_finalise()

end program test_CASTEP_MM_buffer_crack
