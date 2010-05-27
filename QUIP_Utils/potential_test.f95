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

! NOTES:
!
! pot1 is the "cheap" potential
! pot2 is the "expensive" potential
! at should be cubic and large enough
! qw stuff commented out

! TO DO:
!
! subroutine create_interstitial
! comments
! improve create_???_surface and shearing_positions (could be more efficient)
! improve converge_kpoints (args_list)

! TEST SUBROUTINES:
!
! test_bulk_hot_md(pot1, pot2, at, log, temp)
! test_vacancy_hot_md(pot1, pot2, at, log, temp)
! test_interstitial_hot_md(pot1, pot2, at, log, temp)
! test_100_surface_hot_md(pot1, pot2, at, log, temp)
! test_110_surface_hot_md(pot1, pot2, at, log, temp)
! test_111_surface_hot_md(pot1, pot2, at, log, temp)
! test_shear_x_along_z_hot_md(pot1, pot2, at, log, temp, shear_disp)
! test_shear_x_along_y_hot_md(pot1, pot2, at, log, temp, shear_disp)
!
! test_bulk_forces(pot1, pot2, at, log, disp)
! test_vacancy_forces(pot1, pot2, at, log, disp)
! test_interstitial_forces(pot1, pot2, at, log, disp)
! test_100_surface_forces(pot1, pot2, at, log, disp)
! test_110_surface_forces(pot1, pot2, at, log, disp)
! test_111_surface_forces(pot1, pot2, at, log, disp)
! test_shear_x_along_z_forces(pot1, pot2, at, log, disp, shear_disp)
! test_shear_x_along_y_forces(pot1, pot2, at, log, disp, shear_disp)
!
! test_bulk_e(pot1, pot2, at, log)
! test_vacancy_e(pot1, pot2, at, log)
! test_interstitial_e(pot1, pot2, at, log)
! test_100_surface_e(pot1, pot2, at, log)
! test_110_surface_e(pot1, pot2, at, log)
! test_111_surface_e(pot1, pot2, at, log)
! test_shear_x_along_z_e(pot1, pot2, at, log, disp)
! test_shear_x_along_y_e(pot1, pot2, at, log, disp)
!
! test_bulk_elastic_consts(pot1, pot2, at, log)
!
! converge_kpoints(pot, at, log, kpoints_n)
! converge_fermi_t(pot, at, log, fermi_t)

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
! MODULE
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module potential_test_module

  use libAtoms_module
  use QUIP_module
  use elasticity_module
!  use qw_module

  implicit none

  character(len=10) :: default_print_real_format = 'e0.16'  ! Printing format of real variables
  real(dp)          :: default_md_dt = 1.0_dp               ! Timestep
  integer           :: default_md_calc_connect_sep_iter = 5 ! No of iteration between calc connect calls
  integer           :: default_md_temperature_eq_iter = 100 ! No of iterations to equilibrate temperature
  integer           :: default_md_trials_no = 100           ! Trials no
  integer           :: default_md_trials_sep_iter = 100     ! No of iterations between trials
  integer           :: default_forces_trials_no = 100       ! Trials no
!  integer           :: default_qw_order(6) = (/1, 2, 3, 4, 5, 6/) ! Order of qw order parameters to calculate

  real(dp)          :: default_e_rel_tol = 0.01_dp          ! relaxation tolerance
  integer           :: default_e_rel_iter = 100             ! relaxation max iterations

contains

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
! HOT MD TEST SUBROUTINES
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine test_bulk_hot_md(pot1, pot2, at, log, temp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: temp

    call test_hot_md(pot1, pot2, at, log, temp)

  end subroutine test_bulk_hot_md

  subroutine test_vacancy_hot_md(pot1, pot2, at, log, temp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: temp
    type(Atoms) :: at_vac

    call create_vacancy(at, at_vac)
    call test_hot_md(pot1, pot2, at_vac, log, temp)
    call finalise(at_vac)

  end subroutine test_vacancy_hot_md

  subroutine test_interstitial_hot_md(pot1, pot2, at, log, temp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: temp
    type(Atoms) :: at_int

    call create_interstitial(at, at_int, pot1)
    call test_hot_md(pot1, pot2, at_int, log, temp)
    call finalise(at_int)

  end subroutine test_interstitial_hot_md

  subroutine test_100_surface_hot_md(pot1, pot2, at, log, temp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: temp
    type(Atoms) :: at_surf

    call create_100_surface(at, at_surf)
    call test_hot_md(pot1, pot2, at_surf, log, temp)
    call finalise(at_surf)

  end subroutine test_100_surface_hot_md

  subroutine test_110_surface_hot_md(pot1, pot2, at, log, temp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: temp
    type(Atoms) :: at_surf

    call create_110_surface(at, at_surf)
    call test_hot_md(pot1, pot2, at_surf, log, temp)
    call finalise(at_surf)

  end subroutine test_110_surface_hot_md

  subroutine test_111_surface_hot_md(pot1, pot2, at, log, temp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: temp
    type(Atoms) :: at_surf

    call create_111_surface(at, at_surf)
    call test_hot_md(pot1, pot2, at_surf, log, temp)
    call finalise(at_surf)

  end subroutine test_111_surface_hot_md

  subroutine test_shear_x_along_z_hot_md(pot1, pot2, at, log, temp, shear_disp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: temp, shear_disp
    type(Atoms) :: at_shear

    call shear_positions(at, at_shear, 1, 3, shear_disp)
    call test_hot_md(pot1, pot2, at_shear, log, temp)
    call finalise(at_shear)

  end subroutine test_shear_x_along_z_hot_md

  subroutine test_shear_x_along_y_hot_md(pot1, pot2, at, log, temp, shear_disp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: temp, shear_disp
    type(Atoms) :: at_shear

    call shear_positions(at, at_shear, 1, 2, shear_disp)
    call test_hot_md(pot1, pot2, at_shear, log, temp)
    call finalise(at_shear)

  end subroutine test_shear_x_along_y_hot_md

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
! FORCES TEST SUBROUTINES
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine test_bulk_forces(pot1, pot2, at, log, disp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: disp

    call test_forces(pot1, pot2, at, log, disp)

  end subroutine test_bulk_forces

  subroutine test_vacancy_forces(pot1, pot2, at, log, disp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: disp
    type(Atoms) :: at_vac

    call create_vacancy(at, at_vac)
    call test_forces(pot1, pot2, at_vac, log, disp)
    call finalise(at_vac)

  end subroutine test_vacancy_forces

  subroutine test_interstitial_forces(pot1, pot2, at, log, disp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: disp
    type(Atoms) :: at_int

    call create_interstitial(at, at_int, pot1)
    call test_forces(pot1, pot2, at_int, log, disp)
    call finalise(at_int)

  end subroutine test_interstitial_forces

  subroutine test_100_surface_forces(pot1, pot2, at, log, disp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: disp
    type(Atoms) :: at_surf

    call create_100_surface(at, at_surf)
    call test_forces(pot1, pot2, at_surf, log, disp)
    call finalise(at_surf)

  end subroutine test_100_surface_forces

  subroutine test_110_surface_forces(pot1, pot2, at, log, disp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: disp
    type(Atoms) :: at_surf

    call create_110_surface(at, at_surf)
    call test_forces(pot1, pot2, at_surf, log, disp)
    call finalise(at_surf)

  end subroutine test_110_surface_forces

  subroutine test_111_surface_forces(pot1, pot2, at, log, disp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: disp
    type(Atoms) :: at_surf

    call create_111_surface(at, at_surf)
    call test_forces(pot1, pot2, at_surf, log, disp)
    call finalise(at_surf)

  end subroutine test_111_surface_forces

  subroutine test_shear_x_along_z_forces(pot1, pot2, at, log, disp, shear_disp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: disp, shear_disp
    type(Atoms) :: at_shear

    call shear_positions(at, at_shear, 1, 3, shear_disp)
    call test_forces(pot1, pot2, at_shear, log, disp)
    call finalise(at_shear)

  end subroutine test_shear_x_along_z_forces

  subroutine test_shear_x_along_y_forces(pot1, pot2, at, log, disp, shear_disp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: disp, shear_disp
    type(Atoms) :: at_shear

    call shear_positions(at, at_shear, 1, 2, shear_disp)
    call test_forces(pot1, pot2, at_shear, log, disp)
    call finalise(at_shear)

  end subroutine test_shear_x_along_y_forces

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
! HOT MD AND FORCES SUBROUTINES
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine test_hot_md(pot1, pot2, at, file, temperature)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: file
    real(dp), intent(in) :: temperature
    type(Atoms) :: at_array(0:default_md_trials_no)
    real(dp) :: force_errors(4,at%N,default_md_trials_no), forces1(4,at%N,default_md_trials_no), forces2(4,at%N,default_md_trials_no)

    call do_hot_md(pot1, at, at_array, temperature)

    call compare_forces(pot1, pot2, at_array, force_errors, forces1, forces2)

    call print_hot_md_log(file, at_array, force_errors, forces1, forces2, temperature)

  end subroutine test_hot_md

  subroutine print_hot_md_log(file, at_array, force_errors, forces1, forces2, temperature)

    type(Inoutput), intent(inout) :: file
    type(Atoms), intent(inout) :: at_array(0:)
    real(dp), intent(in) :: force_errors(:,:,:), forces1(:,:,:), forces2(:,:,:)
    real(dp), intent(in) :: temperature

    call print("#### START HOT MD LOG", file = file)

    call print("#### Temperature: " // temperature, file = file)
    call print("#### Timestep: " // default_md_dt, file = file)
    call print("#### No of iteration between calc connect calls: " // default_md_calc_connect_sep_iter, file = file)
    call print("#### No of iterations to equilibrate temperature: " // default_md_temperature_eq_iter, file = file)
    call print("#### Trials no: " // default_md_trials_no, file = file)
    call print("#### No of iterations between trials: " // default_md_trials_sep_iter, file = file)

    call print_xyz_log(file, at_array)

    call print_forces_log(file, force_errors, forces1, forces2)

    call print("#### END HOT MD LOG", file = file)

  end subroutine print_hot_md_log

  subroutine read_hot_md_log(file, at_array, force_errors, forces1, forces2, temperature)

    type(Inoutput), intent(inout) :: file
    type(Atoms), intent(out), allocatable :: at_array(:)
    real(dp), intent(out), allocatable :: force_errors(:,:,:), forces1(:,:,:), forces2(:,:,:)
    real(dp), intent(out) :: temperature
    character(len=1000) :: line
    character(len=100), dimension(10) :: line_fields
    integer :: line_fields_no

    line = read_line(file)
    if (trim(line) /= "#### START HOT MD LOG") call system_abort('Error reading HOT MD LOG')

    call parse_line(file, ' ', line_fields, line_fields_no)
    temperature = string_to_real(line_fields(line_fields_no))

    line = read_line(file)
    line = read_line(file)
    line = read_line(file)
    line = read_line(file)
    line = read_line(file)

    call read_xyz_log(file, at_array)

    call read_forces_log(file, force_errors, forces1, forces2)

    line = read_line(file)
    if (trim(line) /= "#### END HOT MD LOG") call system_abort('Error reading HOT MD LOG')

  end subroutine read_hot_md_log

  subroutine test_forces(pot1, pot2, at, file, displacement)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: file
    real(dp), intent(in) :: displacement
    type(Atoms) :: at_array(0:default_forces_trials_no)
    real(dp) :: force_errors(4,at%N,default_forces_trials_no), forces1(4,at%N,default_forces_trials_no), forces2(4,at%N,default_forces_trials_no)
    integer :: i

    at_array(0) = at
    do i = 1, default_forces_trials_no
       call randomise_positions(at, at_array(i), displacement)
    end do

    call compare_forces(pot1, pot2, at_array, force_errors, forces1, forces2)

    call print_test_forces_log(file, at_array, force_errors, forces1, forces2, displacement)

  end subroutine test_forces

  subroutine print_test_forces_log(file, at_array, force_errors, forces1, forces2, displacement)

    type(Inoutput), intent(inout) :: file
    type(Atoms), intent(inout) :: at_array(0:)
    real(dp), intent(in) :: force_errors(:,:,:), forces1(:,:,:), forces2(:,:,:)
    real(dp), intent(in) :: displacement

    call print("#### START TEST FORCES LOG", file = file)

    call print("#### Displacement: " // displacement, file = file)
    call print("#### Trials no: " // default_forces_trials_no, file = file)

    call print_xyz_log(file, at_array)

    call print_forces_log(file, force_errors, forces1, forces2)

    call print("#### END TEST FORCES LOG", file = file)

  end subroutine print_test_forces_log

  subroutine read_test_forces_log(file, at_array, force_errors, forces1, forces2, displacement)

    type(Inoutput), intent(inout) :: file
    type(Atoms), intent(out), allocatable :: at_array(:)
    real(dp), intent(out), allocatable :: force_errors(:,:,:), forces1(:,:,:), forces2(:,:,:)
    real(dp), intent(out) :: displacement
    character(len=1000) :: line
    character(len=100), dimension(10) :: line_fields
    integer :: line_fields_no

    line = read_line(file)
    if (trim(line) /= "#### START TEST FORCES LOG") call system_abort('Error reading TEST FORCES LOG')

    call parse_line(file, ' ', line_fields, line_fields_no)
    displacement = string_to_real(line_fields(line_fields_no))

    line = read_line(file)

    call read_xyz_log(file, at_array)

    call read_forces_log(file, force_errors, forces1, forces2)

    line = read_line(file)
    if (trim(line) /= "#### END TEST FORCES LOG") call system_abort('Error reading TEST FORCES LOG')

  end subroutine read_test_forces_log

  subroutine do_hot_md(pot, at, at_array, temperature)

    type(Potential), intent(inout) :: pot
    type(Atoms), intent(inout) :: at
    type(Atoms), intent(out) :: at_array(0:default_md_trials_no)
    real(dp), intent(in) :: temperature
    type(DynamicalSystem) :: ds
    integer :: i

    at_array(0) = at

    if (associated(pot%tb)) then
       call setup_atoms(pot%tb, at)
    else
       call set_cutoff(at, cutoff(pot))
    end if

    call initialise(ds, at)

    if (temperature /= 0.0_dp) call rescale_velo(ds, temperature)

    call do_md(pot, ds, default_md_dt, default_md_temperature_eq_iter)

    at_array(1) = ds%atoms

    do i = 2, default_md_trials_no
       call do_md(pot, ds, default_md_dt, default_md_trials_sep_iter)

       at_array(i) = ds%atoms
    end do

    call finalise(ds)

  end subroutine do_hot_md

  subroutine do_md(pot, ds, dt, iter)

    type(Potential), intent(inout) :: pot
    type(DynamicalSystem), intent(inout) :: ds
    real(dp), intent(in) :: dt
    integer, intent(in) :: iter
    integer :: i
    real(dp) :: force(3,ds%atoms%N)

    do i = 1, iter
       if (mod(i, default_md_calc_connect_sep_iter) == 1) call calc_connect(ds%atoms)

       call calc(pot, ds%atoms, f = force)

       call advance_verlet(ds, dt, force)
    end do

  end subroutine do_md

  subroutine compare_forces(pot1, pot2, at_array, force_errors, forces1, forces2)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at_array(0:)
    real(dp), intent(out) :: force_errors(4,at_array(0)%N,size(at_array) - 1)
    real(dp), intent(out), optional :: forces1(4,at_array(0)%N,size(at_array) - 1), forces2(4,at_array(0)%N,size(at_array) - 1)
    real(dp) :: forces1_temporary(4,at_array(0)%N,size(at_array) - 1), forces2_temporary(4,at_array(0)%N,size(at_array) - 1)

    call calculate_forces(pot1, at_array, forces1_temporary)
    call calculate_forces(pot2, at_array, forces2_temporary)

    if (present(forces1)) forces1 = forces1_temporary
    if (present(forces2)) forces2 = forces2_temporary

    force_errors = forces1_temporary - forces2_temporary

  end subroutine compare_forces

  subroutine calculate_forces(pot, at_array, forces)

    type(Potential), intent(inout) :: pot
    type(Atoms), intent(inout) :: at_array(0:)
    real(dp), intent(out) :: forces(4,at_array(0)%N,size(at_array) - 1)
    integer :: i, j

    do i = 1, (size(at_array) - 1)

       if (associated(pot%tb)) then
          call setup_atoms(pot%tb, at_array(i))
       else
          call set_cutoff(at_array(i), cutoff(pot))
       end if

       call calc_connect(at_array(i))

       call calc(pot, at_array(i), f = forces(1:3,:,i))

       do j = 1, at_array(0)%N
          forces(4,j,i) = norm(forces(1:3,j,i))
       end do
    end do

  end subroutine calculate_forces

  subroutine print_xyz_log(file, at_array)

    type(Inoutput), intent(inout) :: file
    type(Atoms), intent(inout) :: at_array(0:)
    integer :: i

    call print("## START XYZ LOG", file = file)

    call print("## Trials no: " // (size(at_array) - 1), file = file)

    call print("# Initial atomic configuration:", file = file)

    call print_xyz(at_array(0), file, all_properties = .true., real_format = default_print_real_format)

    do i = 1, (size(at_array) - 1)
       call print("# Trial " // i // ":", file = file)

       call print_xyz(at_array(i), file, all_properties = .true., real_format = default_print_real_format)
    end do

    call print("## END XYZ LOG", file = file)

  end subroutine print_xyz_log

  subroutine read_xyz_log(file, at_array)

    type(Inoutput), intent(inout) :: file
    type(Atoms), intent(out), allocatable :: at_array(:)
    character(len=1000) :: line
    character(len=100), dimension(10) :: line_fields
    integer :: line_fields_no
    integer :: trials_no
    integer :: i

    line = read_line(file)
    if (trim(line) /= "## START XYZ LOG") call system_abort('Error reading XYZ LOG')

    call parse_line(file, ' ', line_fields, line_fields_no)
    trials_no = string_to_int(line_fields(line_fields_no))

    allocate(at_array(0:trials_no))

    do i = 0, trials_no
       line = read_line(file)

       call read_xyz(at_array(i), file)
    end do

    line = read_line(file)
    if (trim(line) /= "## END XYZ LOG") call system_abort('Error reading XYZ LOG')

  end subroutine read_xyz_log

  subroutine print_forces_log(file, force_errors, forces1, forces2)

    type(Inoutput), intent(inout) :: file
    real(dp), intent(in) :: force_errors(:,:,:)
    real(dp), intent(in), optional :: forces1(:,:,:), forces2(:,:,:)
    integer :: i

    call print("## START FORCES LOG", file = file)

    call print("## Trials no: " // size(force_errors, 3), file = file)
    call print("## Atoms no: " // size(force_errors, 2), file = file)

    do i = 1, size(force_errors, 3)
       call print("# Trial " // i // ":", file = file)

       if (present(forces1)) then
          call print("# Forces using potential 1:", file = file)

          call print_forces(file, forces1(:,:,i))
       end if

       if (present(forces2)) then
          call print("# Forces using potential 2:", file = file)

          call print_forces(file, forces2(:,:,i))
       end if

       call print("# Force errors:", file = file)

       call print_forces(file, force_errors(:,:,i))
    end do

    call print("## END FORCES LOG", file = file)

  end subroutine print_forces_log

  subroutine print_forces(file, forces)

    type(Inoutput), intent(inout) :: file
    real(dp), intent(in) :: forces(:,:)
    integer :: i, j
    character(len=1000) :: line
    character(len=1000) :: line_field

    call print("# N    f_x    f_y    f_z    |f|", file = file)

    do i = 1, size(forces, 2)
       line = '' // i

       do j = 1, 4
          write (line_field,'(' // trim(default_print_real_format) // ')') forces(j,i)

          line = trim(line) // ' ' // trim(line_field)
       end do

       call print(line, file = file)
    end do

  end subroutine print_forces

  subroutine read_forces_log(file, force_errors, forces1, forces2)

    type(Inoutput), intent(inout) :: file
    real(dp), intent(out), allocatable :: force_errors(:,:,:)
    real(dp), intent(out), optional, allocatable :: forces1(:,:,:), forces2(:,:,:)
    character(len=1000) :: line
    character(len=100), dimension(10) :: line_fields
    integer :: line_fields_no
    integer :: trials_no
    integer :: atoms_no
    integer :: i, j

    line = read_line(file)
    if (trim(line) /= "## START FORCES LOG") call system_abort('Error reading FORCES LOG')

    call parse_line(file, ' ', line_fields, line_fields_no)
    trials_no = string_to_int(line_fields(line_fields_no))

    call parse_line(file, ' ', line_fields, line_fields_no)
    atoms_no = string_to_int(line_fields(line_fields_no))

    allocate(force_errors(4,atoms_no,trials_no))
    if (present(forces1)) allocate(forces1(4,atoms_no,trials_no))
    if (present(forces2)) allocate(forces2(4,atoms_no,trials_no))

    do i = 1, trials_no
       line = read_line(file)

       line = read_line(file)

       if (trim(line) == "# Forces using potential 1:") then
          if (present(forces1) .and. allocated(forces1)) then
             call read_forces(file, forces1(:,:,i))
          else
             do j = 1, (atoms_no + 1)
                line = read_line(file)
             end do
          end if

          line = read_line(file)
       else
          if (present(forces1) .and. allocated(forces1)) deallocate(forces1)
       end if

       if (trim(line) == "# Forces using potential 2:") then
          if (present(forces2) .and. allocated(forces2)) then
             call read_forces(file, forces2(:,:,i))
          else
             do j = 1, (atoms_no + 1)
                line = read_line(file)
             end do
          end if

          line = read_line(file)
       else
          if (present(forces2) .and. allocated(forces2)) deallocate(forces2)
       end if

       call read_forces(file, force_errors(:,:,i))
    end do

    line = read_line(file)
    if (trim(line) /= "## END FORCES LOG") call system_abort('Error reading FORCES LOG')

  end subroutine read_forces_log

  subroutine read_forces(file, forces)

    type(Inoutput), intent(inout) :: file
    real(dp), intent(out) :: forces(:,:)
    character(len=1000) :: line
    integer :: i
    integer :: atom

    line = read_line(file)

    do i = 1, size(forces, 2)
       read (file%unit,*) atom, forces(1,i), forces(2,i), forces(3,i), forces(4,i)
    end do

  end subroutine read_forces

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
! QW SUBROUTINES
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!  subroutine calculate_qw(at_array, qw, cutoff)
!
!    type(Atoms), intent(inout) :: at_array(0:)
!    real(dp), intent(out) :: qw(size(default_qw_order) * 2,at_array(0)%N,size(at_array) - 1)
!    real(dp), intent(in) :: cutoff
!
!    call calculate_atoms_qw(at_array, cutoff)
!
!    call copy_atoms_qw(at_array, qw)
!
!    call remove_atoms_qw(at_array)
!
!  end subroutine calculate_qw

!  subroutine calculate_atoms_qw(at_array, cutoff)
!
!    type(Atoms), intent(inout) :: at_array(0:)
!    real(dp), intent(in) :: cutoff
!    integer :: i, j
!
!    do i = 1, (size(at_array) - 1)
!       call set_cutoff(at_array(i), cutoff)
!
!      Not necessary with new continous qw
!       call calc_connect(at_array(i))
!
!       do j = 1, size(default_qw_order)
!          call calc_qw_at(at_array(i), default_qw_order(j))
!       end do
!    end do
!
!  end subroutine calculate_atoms_qw

!  subroutine copy_atoms_qw(at_array, qw)
!
!    type(Atoms), intent(inout) :: at_array(0:)
!    real(dp), intent(out) :: qw(size(default_qw_order) * 2,at_array(0)%N,size(at_array) - 1)
!    integer :: i, j
!    real(dp), pointer :: q(:), w(:)
!
!    do i = 1, (size(at_array) - 1)
!       do j = 1, size(default_qw_order)
!          if (.not. assign_pointer(at_array(i), 'q' // default_qw_order(j), q)) call system_abort('Could not assign pointer q to atoms object')
!          if (.not. assign_pointer(at_array(i), 'w' // default_qw_order(j), w)) call system_abort('Could not assign pointer w to atoms object')
!
!          qw((j * 2) - 1,:,i) = q
!          qw(j * 2,:,i) = w
!       end do
!    end do
!
!  end subroutine copy_atoms_qw

!  subroutine copy_qw_atoms(qw, at_array)
!
!    real(dp), intent(in) :: qw(:,:,:)
!    type(Atoms), intent(inout) :: at_array(0:size(qw, 3))
!    integer :: i, j
!
!    do i = 1, (size(at_array) - 1)
!       do j = 1, size(default_qw_order)
!          call add_property(at_array(i), 'q' // default_qw_order(j), qw((j * 2) - 1,:,i))
!          call add_property(at_array(i), 'w' // default_qw_order(j), qw(j * 2,:,i))
!       end do
!    end do
!
!  end subroutine copy_qw_atoms

!  subroutine remove_atoms_qw(at_array)
!
!    type(Atoms), intent(inout) :: at_array(0:)
!    integer :: i, j
!
!    do i = 1, (size(at_array) - 1)
!       do j = 1, size(default_qw_order)
!          call remove_property(at_array(i), 'q' // default_qw_order(j))
!          call remove_property(at_array(i), 'w' // default_qw_order(j))
!       end do
!    end do
!
!  end subroutine remove_atoms_qw

!  subroutine print_atoms_qw_log(file, at_array)
!
!    type(Inoutput), intent(inout) :: file
!    type(Atoms), intent(inout) :: at_array(0:)
!    character(len=1000) :: properties
!    integer :: i
!
!    call print("## START QW LOG", file = file)
!
!    call print("## Trials no: " // (size(at_array) - 1), file = file)
!
!    properties = 'pos'
!
!    do i = 1, size(default_qw_order)
!       properties = trim(properties) // ':q' // default_qw_order(i) // ':w' // default_qw_order(i)
!    end do
!
!    do i = 1, (size(at_array) - 1)
!       call print("# Trial " // i // ":", file = file)
!
!       call print_xyz(at_array(i), file, properties = trim(properties), real_format = default_print_real_format)
!    end do
!
!    call print("## END QW LOG", file = file)
!
!  end subroutine print_atoms_qw_log

!  subroutine read_atoms_qw_log(file, at_array)
!
!    type(Inoutput), intent(inout) :: file
!    type(Atoms), intent(out), allocatable :: at_array(:)
!    character(len=1000) :: line
!    character(len=100), dimension(10) :: line_fields
!    integer :: line_fields_no
!    integer :: trials_no
!    integer :: i
!
!    line = read_line(file)
!    if (trim(line) /= "## START QW LOG") call system_abort('Error reading QW LOG')
!
!    call parse_line(file, ' ', line_fields, line_fields_no)
!    trials_no = string_to_int(line_fields(line_fields_no))
!
!    allocate(at_array(0:trials_no))
!
!    do i = 1, trials_no
!       line = read_line(file)
!
!       call read_xyz(at_array(i), file)
!    end do
!
!    line = read_line(file)
!    if (trim(line) /= "## END QW LOG") call system_abort('Error reading QW LOG')
!
!  end subroutine read_atoms_qw_log

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
! ENERGY TEST SUBROUTINES
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine test_bulk_e(pot1, pot2, at, log)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log

    call print("# :::: TEST BULK ENERGY", file = log)

    call test_energy(pot1, pot2, at, log)

    call print("", file = log)

  end subroutine test_bulk_e

  subroutine test_vacancy_e(pot1, pot2, at, log)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    type(Atoms) :: at_vac

    call print("# :::: TEST VACANCY ENERGY", file = log)

    call create_vacancy(at, at_vac)

    call test_energy(pot1, pot2, at, log, at_vac)

    call finalise(at_vac)

    call print("", file = log)

  end subroutine test_vacancy_e

  subroutine test_interstitial_e(pot1, pot2, at, log)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    type(Atoms) :: at_int

    call print("# :::: TEST INTERSTITIAL ENERGY", file = log)

    call create_interstitial(at, at_int, pot1)

    call test_energy(pot1, pot2, at, log, at_int)

    call finalise(at_int)

    call print("", file = log)

  end subroutine test_interstitial_e

  subroutine test_100_surface_e(pot1, pot2, at, log)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    type(Atoms) :: at_surf

    call print("# :::: TEST 100 SURFACE ENERGY", file = log)

    call create_100_surface(at, at_surf)

    call test_energy(pot1, pot2, at, log, at_surf, def_no = 2)

    call finalise(at_surf)

    call print("", file = log)

  end subroutine test_100_surface_e

  subroutine test_110_surface_e(pot1, pot2, at, log)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    type(Atoms) :: at_surf

    call print("# :::: TEST 110 SURFACE ENERGY", file = log)

    call create_110_surface(at, at_surf)

    call test_energy(pot1, pot2, at, log, at_surf, def_no = 2)

    call finalise(at_surf)

    call print("", file = log)

  end subroutine test_110_surface_e

  subroutine test_111_surface_e(pot1, pot2, at, log)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    type(Atoms) :: at_surf

    call print("# :::: TEST 111 SURFACE ENERGY", file = log)

    call create_111_surface(at, at_surf)

    call test_energy(pot1, pot2, at, log, at_surf, def_no = 2)

    call finalise(at_surf)

    call print("", file = log)

  end subroutine test_111_surface_e

  subroutine test_shear_x_along_z_e(pot1, pot2, at, log, disp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: disp
    type(Atoms) :: at_shear

    call print("# :::: TEST SHEAR X ALONG Z ENERGY", file = log)

    call shear_positions(at, at_shear, 1, 3, disp)

    call test_energy(pot1, pot2, at, log, at_shear)

    call finalise(at_shear)

    call print("", file = log)

  end subroutine test_shear_x_along_z_e

  subroutine test_shear_x_along_y_e(pot1, pot2, at, log, disp)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: disp
    type(Atoms) :: at_shear

    call print("# :::: TEST SHEAR X ALONG Y ENERGY", file = log)

    call shear_positions(at, at_shear, 1, 2, disp)

    call test_energy(pot1, pot2, at, log, at_shear)

    call finalise(at_shear)

    call print("", file = log)

  end subroutine test_shear_x_along_y_e

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
! ELASTIC CONSTANTS TEST SUBROUTINES
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine test_bulk_elastic_consts(pot1, pot2, at, log)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    type(MetaPotential) :: mpot1, mpot2
    real(dp) :: c(6,6), c0(6,6)

    call print("# :::: TEST BULK ELASTIC CONSTANTS", file = log)

    call print_bulk_structure(log, at)

    if (associated(pot1%ip)) call initialise(mpot1, 'Simple', pot1, mpi_obj = pot1%ip%mpi_glob)
    if (associated(pot1%tb)) call initialise(mpot1, 'Simple', pot1, mpi_obj = pot1%tb%mpi)
    if (associated(pot1%filepot)) call initialise(mpot1, 'Simple', pot1, mpi_obj = pot1%filepot%mpi)

    if (associated(pot2%ip)) call initialise(mpot2, 'Simple', pot2, mpi_obj = pot2%ip%mpi_glob)
    if (associated(pot2%tb)) call initialise(mpot2, 'Simple', pot2, mpi_obj = pot2%tb%mpi)
    if (associated(pot2%filepot)) call initialise(mpot2, 'Simple', pot2, mpi_obj = pot2%filepot%mpi)

    call calc_elastic_constants(mpot1, at, c = c, c0 = c0, relax_initial = .false.)

    call print("# :: Using Potential 1:", file = log)
    call print("# elastic constants (with relaxation):", file = log)
    call print(c, file = log)
    call print("# elastic constants (without relaxation)", file = log)
    call print(c0, file = log)

    call calc_elastic_constants(mpot2, at, c = c, c0 = c0, relax_initial = .false.)

    call print("# :: Using Potential 2:", file = log)
    call print("# elastic constants (with relaxation):", file = log)
    call print(c, file = log)
    call print("# elastic constants (without relaxation)", file = log)
    call print(c0, file = log)

    call finalise(mpot1)
    call finalise(mpot2)

    call print("", file = log)

  end subroutine test_bulk_elastic_consts

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
! CONVERGENCE TEST SUBROUTINES
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine converge_kpoints(pot, at, log, kpoints_n)

    type(Potential), intent(inout) :: pot
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    integer, intent(in) :: kpoints_n(:,:)
    type(Atoms) :: at_rand
    type(KPoints) :: kpoints_original
    type(Potential) :: pot_temporary
    integer :: i
    real(dp) :: energy_stats(3,size(kpoints_n, 2)), force_stats(4,at%N,size(kpoints_n, 2))

    call print("# :::: CONVERGE KPOINTS", file = log)

    if (current_verbosity() > NORMAL) then
       call print("# potential:", file = log)
       call print(pot, file = log)
    end if

    call randomise_positions(at, at_rand, 0.5_dp)

    call print_bulk_structure(log, at_rand)

    kpoints_original = pot%tb%tbsys%kpoints

    do i = 1, size(kpoints_n, 2)
       call print("# :: Number of Kpoints: " // kpoints_n(1,i) // " " // kpoints_n(2,i) // " " // kpoints_n(3,i), file = log)

       call finalise(pot%tb%tbsys%kpoints)
       call initialise(pot%tb%tbsys%kpoints, kpoints_n(:,i), monkhorst_pack = .true., mpi_obj = pot%tb%mpi)

       call finalise(pot_temporary)
       allocate(pot_temporary%tb)
       call finalise(pot_temporary%tb)
       pot_temporary%tb%init_args_str = pot%tb%init_args_str
       pot_temporary%tb%mpi = pot%tb%mpi
       call initialise(pot_temporary%tb%tbsys, pot%tb%tbsys, mpi_obj = pot%tb%mpi)

       call evaluate_energies(pot_temporary, at_rand, log, energy_stats(:,i))
       call evaluate_forces(pot_temporary, at_rand, log, force_stats(:,:,i))

       call finalise(pot_temporary)
    end do

    call print("# :: Overall Statistics:", file = log)
    call print_energy_statistics(log, energy_stats, args_name = 'kpoints n trial')
    call print_force_statistics(log, force_stats, atoms = (/1/), args_name = 'kpoints n trial')

    call finalise(pot%tb%tbsys%kpoints)
    pot%tb%tbsys%kpoints = kpoints_original
    call finalise(kpoints_original)

    call finalise(at_rand)

    call print("", file = log)

  end subroutine converge_kpoints

  subroutine converge_fermi_t(pot, at, log, fermi_t)

    type(Potential), intent(inout) :: pot
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: fermi_t(:)
    type(Atoms) :: at_rand
    logical :: has_default_fermi_t_original
    real(dp) :: default_fermi_t_original
    integer :: i
    real(dp) :: energy_stats(3,size(fermi_t)), force_stats(4,at%N,size(fermi_t))

    call print("# :::: CONVERGE FERMI TEMPERATURE", file = log)

    if (current_verbosity() > NORMAL) then
       call print("# potential:", file = log)
       call print(pot, file = log)
    end if

    call randomise_positions(at, at_rand, 0.5_dp)

    call print_bulk_structure(log, at_rand)

    has_default_fermi_t_original = pot%tb%tbsys%tbmodel%has_default_fermi_T
    default_fermi_t_original = pot%tb%tbsys%tbmodel%default_fermi_T

    do i = 1, size(fermi_t)
       call print("# :: Fermi Temperature: " // fermi_t(i), file = log)

       pot%tb%tbsys%tbmodel%has_default_fermi_T = .true.
       pot%tb%tbsys%tbmodel%default_fermi_T = fermi_t(i)

       call evaluate_energies(pot, at_rand, log, energy_stats(:,i))
       call evaluate_forces(pot, at_rand, log, force_stats(:,:,i))
    end do

    call print("# :: Overall Statistics:", file = log)
    call print_energy_statistics(log, energy_stats, args_name = 'fermi t', args_list = fermi_t)
    call print_force_statistics(log, force_stats, atoms = (/1/), args_name = 'fermi t', args_list = fermi_t)

    pot%tb%tbsys%tbmodel%has_default_fermi_T = has_default_fermi_t_original
    pot%tb%tbsys%tbmodel%default_fermi_T = default_fermi_t_original

    call finalise(at_rand)

    call print("", file = log)

  end subroutine converge_fermi_t

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
! OLD HOT MD SUBROUTINES
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine evaluate_forces(pot, at, log, forces)

    type(Potential), intent(inout) :: pot
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(out), optional :: forces(4,at%N)
    real(dp) :: force(3,at%N), force_norm(at%N)
    integer :: i

    if (associated(pot%tb)) then
       call setup_atoms(pot%tb, at)
    else
       call set_cutoff(at, cutoff(pot))
    end if

    call calc_connect(at)

    call calc(pot, at, f = force)

    do i = 1, at%N
       force_norm(i) = norm(force(:,i))
    end do

    call print("# forces:", file = log)
    call print("# N    f_x    f_y    f_z    |f|", file = log)
    do i = 1, at%N
       call print(i // "    " // force(1,i) // "    " // force(2,i) // "    " // force(3,i) // "    " // force_norm(i), file = log)
    end do

    if (present(forces)) then
       forces(1:3,:) = force(:,:)
       forces(4,:) = force_norm(:)
    end if

  end subroutine evaluate_forces

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
! ENERGY SUBROUTINES
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine test_energy(pot1, pot2, at_bulk, log, at_def, def_no)

    type(Potential), intent(inout) :: pot1, pot2
    type(Atoms), intent(inout) :: at_bulk
    type(Inoutput), intent(inout) :: log
    type(Atoms), intent(inout), optional :: at_def
    integer, intent(in), optional :: def_no

    if (present(at_def)) then
       call print("# relaxation tolerance:      " // default_e_rel_tol, file = log)
       call print("# relaxation max iterations: " // default_e_rel_iter, file = log)
       call print_bulk_structure(log, at_def)
    else
       call print_bulk_structure(log, at_bulk)
    end if

    call print("# :: Using Potential 1:", file = log)

    if (present(at_def)) then
       call calc_defect_energy(pot1, at_bulk, log, at_def, def_no = def_no)
    else
       call evaluate_energies(pot1, at_bulk, log)
    end if

    call print("# :: Using Potential 2:", file = log)

    if (present(at_def)) then
       call calc_defect_energy(pot2, at_bulk, log, at_def, def_no = def_no)
    else
       call evaluate_energies(pot2, at_bulk, log)
    end if

  end subroutine test_energy

  subroutine evaluate_energies(pot, at, log, energies)

    type(Potential), intent(inout) :: pot
    type(Atoms), intent(inout) :: at
    type(Inoutput), intent(inout) :: log
    real(dp), intent(out), optional :: energies(3)
    real(dp) :: bulk_e, atom_e, bind_e

    call calc_energy(pot, at, bulk_e, atom_e, bind_e)

    call print("bulk energy (per atom):    " // bulk_e, file = log)
    call print("atom energy (per atom):    " // atom_e, file = log)
    call print("binding energy (per atom): " // bind_e, file = log)

    if (present(energies)) then
       energies(1) = bulk_e
       energies(2) = atom_e
       energies(3) = bind_e
    end if

  end subroutine evaluate_energies

  subroutine calc_defect_energy(pot, at_bulk, log, at_def, def_no, def_e)

    type(Potential), intent(inout) :: pot
    type(Atoms), intent(inout) :: at_bulk, at_def
    type(Inoutput), intent(inout) :: log
    integer, intent(in), optional :: def_no
    real(dp), intent(out), optional :: def_e
    type(MetaPotential) :: mpot
    real(dp) :: def_bulk_e, bulk_e

    if (associated(pot%ip)) call initialise(mpot, 'Simple', pot, mpi_obj = pot%ip%mpi_glob)
    if (associated(pot%tb)) call initialise(mpot, 'Simple', pot, mpi_obj = pot%tb%mpi)
    if (associated(pot%filepot)) call initialise(mpot, 'Simple', pot, mpi_obj = pot%filepot%mpi)

    call do_relax(mpot, at_def, default_e_rel_tol, default_e_rel_iter, const_vol = .true.)

    call calc_energy(mpot%pot, at_def, def_bulk_e)
    call calc_energy(mpot%pot, at_bulk, bulk_e)

!    call calc_connect(at_def)
!    do i = 1, size(default_qw_order)
!       call calc_qw_at(at_def, default_qw_order(i))
!    end do

    call print_bulk_structure(log, at_def)

    call print("bulk energy (with defect, relaxed, per atom): " // def_bulk_e, file = log)
    call print("bulk energy (without defect, per atom):       " // bulk_e, file = log)

    if (present(def_e)) then
       if (present(def_no)) then
          def_e = ((def_bulk_e - bulk_e) * real(at_def%N, dp)) / real(def_no, dp)
       else
          def_e = (def_bulk_e - bulk_e) * real(at_def%N, dp)
       end if

       call print("defect energy:                                " // def_e, file = log)
    else
       if (present(def_no)) then
          call print("defect energy:                                " // (((def_bulk_e - bulk_e) * real(at_def%N, dp)) / real(def_no, dp)), file = log)
       else
          call print("defect energy:                                " // ((def_bulk_e - bulk_e) * real(at_def%N, dp)), file = log)
       end if
    end if

    call finalise(mpot)

  end subroutine calc_defect_energy

  subroutine do_relax(mpot, at, rel_tol, rel_iter, n_iter, const_vol)

    type(MetaPotential), intent(inout) :: mpot
    type(Atoms), intent(inout) :: at
    real(dp), intent(in) :: rel_tol
    integer, intent(in) :: rel_iter
    integer, intent(out), optional :: n_iter
    logical, intent(in), optional :: const_vol
    integer :: n_iter_temporary

    if (associated(mpot%pot%tb)) then
       call setup_atoms(mpot%pot%tb, at)
    else
       call set_cutoff(at, cutoff(mpot%pot))
    end if

    call calc_connect(at)

    n_iter_temporary = minim(mpot, at, 'cg', rel_tol, rel_iter, linminroutine = 'FAST_LINMIN', do_pos = .true., do_lat = (.not. const_vol))

    if (present(n_iter)) n_iter = n_iter_temporary

  end subroutine do_relax

  subroutine calc_energy(pot, at, bulk_e, atom_e, bind_e)

    type(Potential), intent(inout) :: pot
    type(Atoms), intent(inout) :: at
    real(dp), intent(out) :: bulk_e
    real(dp), intent(out), optional :: atom_e, bind_e
    type(Atoms) :: at_temporary
    real(dp) :: atom_e_temporary

    if (associated(pot%tb)) then
       call setup_atoms(pot%tb, at)
    else
       call set_cutoff(at, cutoff(pot))
    end if

    call calc_connect(at)

    call calc(pot, at, e = bulk_e)

    bulk_e = bulk_e / real(at%N, dp)

    if (present(atom_e) .or. present(bind_e)) then
       if (associated(pot%tb)) then
          ! 10.0_dp to fix NRL-TB, GSP crashes (should be 1.0_dp)
          call initialise(at_temporary, 1, reshape((/cutoff(pot) + 10.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, cutoff(pot) + 10.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, cutoff(pot) + 10.0_dp/), (/3, 3/)))
          at_temporary%pos(:,1) = (/0.0_dp, 0.0_dp, 0.0_dp/)
          at_temporary%Z(1) = at%Z(1)

          call setup_atoms(pot%tb, at_temporary)

          call calc_connect(at_temporary)

          call calc(pot, at_temporary, e = atom_e_temporary)

          call finalise(at_temporary)
       else
          atom_e_temporary = 0.0_dp
       end if

       if (present(atom_e)) atom_e = atom_e_temporary
       if (present(bind_e)) bind_e = bulk_e - atom_e_temporary
    end if

  end subroutine calc_energy

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
! PRINT SUBROUTINES
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine print_bulk_structure(log, at)

    type(Inoutput), intent(inout) :: log
    type(Atoms), intent(inout) :: at

    call print("# bulk structure:", file = log)

    call print_xyz(at, log, all_properties = .true.)

  end subroutine print_bulk_structure

  subroutine print_energy_statistics(log, energy_stats, args_name, args_list)

    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: energy_stats(:,:)
    character(len=*), intent(in), optional :: args_name
    real(dp), intent(in), optional :: args_list(:)
    integer :: i

    if (present(args_name)) then
       call print("# " // trim(args_name) // "    bulk_e    atom_e    bind_e", file = log)
    else
       call print("# n    bulk_e    atom_e    bind_e", file = log)
    end if
    do i = 1, size(energy_stats, 2)
       if (present(args_list)) then
          call print(args_list(i) // "    " // energy_stats(1,i) // "    " // energy_stats(2,i) // "    " // energy_stats(3,i), file = log)
       else
          call print(i // "    " // energy_stats(1,i) // "    " // energy_stats(2,i) // "    " // energy_stats(3,i), file = log)
       end if
    end do

  end subroutine print_energy_statistics

  subroutine print_force_statistics(log, force_stats, atoms, args_name, args_list)

    type(Inoutput), intent(inout) :: log
    real(dp), intent(in) :: force_stats(:,:,:)
    integer, intent(in), optional :: atoms(:)
    character(len=*), intent(in), optional :: args_name
    real(dp), intent(in), optional :: args_list(:)
    integer :: i, j

    if (present(atoms)) then
       do i = 1, size(atoms)
          call print("# atom: " // atoms(i), file = log)
          if (present(args_name)) then
             call print("# " // trim(args_name) // "    f_x    f_y    f_z    |f|", file = log)
          else
             call print("# n    f_x    f_y    f_z    |f|", file = log)
          end if
          do j = 1, size(force_stats, 3)
             if (present(args_list)) then
                call print(args_list(j) // "    " // force_stats(1,atoms(i),j) // "    " // force_stats(2,atoms(i),j) // "    " // force_stats(3,atoms(i),j) // "    " // force_stats(4,atoms(i),j), file = log)
             else
                call print(j // "    " // force_stats(1,atoms(i),j) // "    " // force_stats(2,atoms(i),j) // "    " // force_stats(3,atoms(i),j) // "    " // force_stats(4,atoms(i),j), file = log)
             end if
          end do
       end do
    else
       do i = 1, size(force_stats, 2)
          call print("# atom: " // i, file = log)
          if (present(args_name)) then
             call print("# " // trim(args_name) // "    f_x    f_y    f_z    |f|", file = log)
          else
             call print("# n    f_x    f_y    f_z    |f|", file = log)
          end if
          do j = 1, size(force_stats, 3)
             if (present(args_list)) then
                call print(args_list(j) // "    " // force_stats(1,i,j) // "    " // force_stats(2,i,j) // "    " // force_stats(3,i,j) // "    " // force_stats(4,i,j), file = log)
             else
                call print(j // "    " // force_stats(1,i,j) // "    " // force_stats(2,i,j) // "    " // force_stats(3,i,j) // "    " // force_stats(4,i,j), file = log)
             end if
          end do
       end do
    end if

  end subroutine print_force_statistics

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
! STRUCTURES SUBROUTINES
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !
  !% Create unit cell with a single vacancy
  !
  subroutine create_vacancy(at_bulk, at_vac, vac_position)

    type(Atoms), intent(in) :: at_bulk
    type(Atoms), intent(out) :: at_vac
    real(dp), intent(out), optional :: vac_position(3)

    at_vac = at_bulk

    if (present(vac_position)) vac_position = at_vac%pos(:,1)

    call remove_atoms(at_vac, 1)

  end subroutine create_vacancy

  !
  !% Create unit cell with a single interstitial
  !
  subroutine create_interstitial(at_bulk, at_int, pot, int_position)

    type(Atoms), intent(in) :: at_bulk
    type(Atoms), intent(out) :: at_int
    type(Potential), intent(in), optional :: pot
    real(dp), intent(out), optional :: int_position(3)

    at_int = at_bulk

    if (present(int_position)) int_position = (/0.0_dp, 0.0_dp, 0.0_dp/)

    call add_atoms(at_int, (/0.0_dp, 0.0_dp, 0.0_dp/), at_int%Z(1))

  end subroutine create_interstitial

  !
  !% Create rectengular unit cell with 100 surface
  !
  subroutine create_100_surface(at_bulk, at_surf, surf_indices)

    type(Atoms), intent(in) :: at_bulk
    type(Atoms), intent(out) :: at_surf
    integer, intent(out), optional :: surf_indices(3)

    at_surf = at_bulk

    at_surf%lattice(:,1) = 2.0_dp * at_surf%lattice(:,1)

    call matrix3x3_inverse(at_surf%lattice, at_surf%g)

    if (present(surf_indices)) surf_indices = (/2, 0, 0/)

  end subroutine create_100_surface

  !
  !% Create rectengular unit cell with 110 surface
  !
  subroutine create_110_surface(at_bulk, at_surf, surf_indices)

    type(Atoms), intent(in) :: at_bulk
    type(Atoms), intent(out) :: at_surf
    integer, intent(out), optional :: surf_indices(3)
    type(Atoms) :: at_temporary
    real(dp) :: lattice(3,3), rot_matrix(3,3), pos(3)
    integer :: N

    ! New lattice vectors:
    ! a = (  1  1 0 )
    ! b = ( -1  1 0 )
    ! c = (  0  0 1 )
    !
    ! To get rotation matrix solve:
    ! ( sqrt(2) )   (            )( 1 )  ( 0       )   (            )( -1 )
    ! ( 0       ) = ( rot_matrix )( 1 )  ( sqrt(2) ) = ( rot_matrix )(  1 )
    ! ( 0       )   (            )( 0 ), ( 0       )   (            )(  0 )
    !
    ! Rotation matrix:
    ! (  1/sqrt(2)   1/sqrt(2)    0 )
    ! ( -1/sqrt(2)   1/sqrt(2)    0 )
    ! (  0           0            1 )

    lattice(:,1) = sqrt(2.0_dp) * at_bulk%lattice(:,1)
    lattice(:,2) = sqrt(2.0_dp) * at_bulk%lattice(:,2)
    lattice(:,3) = at_bulk%lattice(:,3)

    rot_matrix = reshape((/1.0_dp / sqrt(2.0_dp), -1.0_dp / sqrt(2.0_dp), 0.0_dp, 1.0_dp / sqrt(2.0_dp), 1.0_dp / sqrt(2.0_dp), 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/), (/3, 3/))

    call initialise(at_surf, 0, lattice)

    call supercell(at_temporary, at_bulk, 2, 3, 1)

    do N = 1, at_temporary%N
       pos = rot_matrix .mult. at_temporary%pos(:,N)

       ! NUMERICAL_ZERO to fix inequality testing
       pos = pos + (/NUMERICAL_ZERO, NUMERICAL_ZERO, NUMERICAL_ZERO/)

       if ((pos(1) >= lattice(1,1)) .and. (pos(1) < (2.0_dp * lattice(1,1))) .and. (pos(2) >= 0.0_dp) .and. (pos(2) < lattice(2,2))) then
          ! NUMERICAL_ZERO to fix inequality testing
          pos = pos - (/NUMERICAL_ZERO, NUMERICAL_ZERO, NUMERICAL_ZERO/)

          pos(1) = pos(1) - lattice(1,1)

          call add_atoms(at_surf, pos, at_temporary%Z(N))
       end if
    end do

    call finalise(at_temporary)

    at_surf%lattice(:,1) = 2.0_dp * at_surf%lattice(:,1)

    call matrix3x3_inverse(at_surf%lattice, at_surf%g)

    if (present(surf_indices)) surf_indices = (/2, 0, 0/)

  end subroutine create_110_surface

  !
  !% Create rectengular unit cell with 111 surface
  !
  subroutine create_111_surface(at_bulk, at_surf, surf_indices)

    type(Atoms), intent(in) :: at_bulk
    type(Atoms), intent(out) :: at_surf
    integer, intent(out), optional :: surf_indices(3)
    type(Atoms) :: at_temporary
    real(dp) :: lattice(3,3), rot_matrix(3,3), pos(3)
    integer :: N

    ! New lattice vectors:
    ! a = (  1  1 1 )
    ! b = ( -1  1 0 )
    ! c = ( -1 -1 2 )
    !
    ! To get rotation matrix solve:
    ! ( sqrt(3) )   (            )( 1 )  ( 0       )   (            )( -1 )  ( 0       )   (            )( -1 )
    ! ( 0       ) = ( rot_matrix )( 1 )  ( sqrt(2) ) = ( rot_matrix )(  1 )  ( 0       ) = ( rot_matrix )( -1 )
    ! ( 0       )   (            )( 1 ), ( 0       )   (            )(  0 ), ( sqrt(6) )   (            )(  2 )
    !
    ! Rotation matrix:
    ! (  1/sqrt(3)    1/sqrt(3)   1/sqrt(3) )
    ! ( -1/sqrt(2)    1/sqrt(2)   0         )
    ! ( -1/sqrt(6)   -1/sqrt(6)   2/sqrt(6) )

    lattice(:,1) = sqrt(3.0_dp) * at_bulk%lattice(:,1)
    lattice(:,2) = sqrt(2.0_dp) * at_bulk%lattice(:,2)
    lattice(:,3) = sqrt(6.0_dp) * at_bulk%lattice(:,3)

    rot_matrix = reshape((/1.0_dp / sqrt(3.0_dp), -1.0_dp / sqrt(2.0_dp), -1.0_dp / sqrt(6.0_dp), 1.0_dp / sqrt(3.0_dp), 1.0_dp / sqrt(2.0_dp), -1.0_dp / sqrt(6.0_dp), 1.0_dp / sqrt(3.0_dp), 0.0_dp, 2.0_dp / sqrt(6.0_dp)/), (/3, 3/))

    call initialise(at_surf, 0, lattice)

    call supercell(at_temporary, at_bulk, 3, 4, 5)

    do N = 1, at_temporary%N
       pos = rot_matrix .mult. at_temporary%pos(:,N)

       ! NUMERICAL_ZERO to fix inequality testing
       pos = pos + (/NUMERICAL_ZERO, NUMERICAL_ZERO, NUMERICAL_ZERO/)

       if ((pos(1) >= (2.0_dp * lattice(1,1))) .and. (pos(1) < (3.0_dp * lattice(1,1))) .and. (pos(2) >= 0.0_dp) .and. (pos(2) < lattice(2,2)) .and. (pos(3) >= 0.0_dp) .and. (pos(3) < lattice(3,3))) then
          ! NUMERICAL_ZERO to fix inequality testing
          pos = pos - (/NUMERICAL_ZERO, NUMERICAL_ZERO, NUMERICAL_ZERO/)

          pos(1) = pos(1) - (2.0_dp * lattice(1,1))

          call add_atoms(at_surf, pos, at_temporary%Z(N))
       end if
    end do

    call finalise(at_temporary)

    at_surf%lattice(:,1) = 2.0_dp * at_surf%lattice(:,1)

    call matrix3x3_inverse(at_surf%lattice, at_surf%g)

    if (present(surf_indices)) surf_indices = (/2, 0, 0/)

  end subroutine create_111_surface

  !
  !% Shear atom positions in a given direction and along a given axis by quantities in the range displacement*(0,Lattice Constant)
  !
  subroutine shear_positions(at_orig, at_shear, direction, axis, displacement)

    type(Atoms), intent(in) :: at_orig
    type(Atoms), intent(out) :: at_shear
    integer, intent(in) :: direction, axis
    real(dp), intent(in) :: displacement
    type(Atoms) :: at_temporary
    real(dp) :: shear_matrix(3,3), pos(3)
    integer :: N

    shear_matrix = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/), (/3, 3/))

    shear_matrix(direction,axis) = displacement

    call initialise(at_shear, 0, at_orig%lattice)

    if (direction == 1) call supercell(at_temporary, at_orig, int(displacement) + 2, 1, 1)
    if (direction == 2) call supercell(at_temporary, at_orig, 1, int(displacement) + 2, 1)
    if (direction == 3) call supercell(at_temporary, at_orig, 1, 1, int(displacement) + 2)

    do N = 1, at_temporary%N
       pos = shear_matrix .mult. at_temporary%pos(:,N)

       ! NUMERICAL_ZERO to fix inequality testing
       pos(direction) = pos(direction) + NUMERICAL_ZERO

       if ((pos(direction) >= ((aint(displacement, dp) + 1.0_dp) * at_orig%lattice(direction,direction))) .and. (pos(direction) < ((aint(displacement, dp) + 2.0_dp) * at_orig%lattice(direction,direction)))) then
          ! NUMERICAL_ZERO to fix inequality testing
          pos(direction) = pos(direction) - NUMERICAL_ZERO

          pos(direction) = pos(direction) - ((aint(displacement, dp) + 1.0_dp) * at_orig%lattice(direction,direction))

          call add_atoms(at_shear, pos, at_temporary%Z(N))
       end if
    end do

  end subroutine shear_positions

  !
  !% Add to atom positions random quantities in the range displacement*(-Element Covalent Radius,Element Covalent Radius)
  !
  subroutine randomise_positions(at_orig, at_rand, displacement)

    type(Atoms), intent(in) :: at_orig
    type(Atoms), intent(out) :: at_rand
    real(dp), intent(in) :: displacement
    integer :: N

    at_rand = at_orig

    do N = 1, at_rand%N
       call randomise(at_rand%pos(:,N), 2.0_dp * displacement * ElementCovRad(at_rand%Z(N)))
    end do

  end subroutine randomise_positions

end module potential_test_module
