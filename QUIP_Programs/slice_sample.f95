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

#include "error.inc"

program slice_sample

use libAtoms_module
use Potential_module

implicit none

  type(Dictionary) :: cli_params
  type(Potential) :: pot
  type(Atoms) :: at
  real(dp) :: lattice_vec(9), atom_1(3), atom_2(3)
  integer :: Z_1, Z_2
  logical :: has_atom_2
  character(len=STRING_LENGTH) :: param_file, init_args
  integer :: n_configs, m_max, rng_seed, init_d
  real(dp) :: lattice_delta, atom_delta, e0, t
  integer :: n, d
  real(dp) :: lattice(3,3), theta
  real(dp), allocatable :: x(:), x_dash(:)

  call system_initialise()

  call enable_timing()

  call initialise(cli_params)
  call param_register(cli_params, "lattice", PARAM_MANDATORY, lattice_vec, help_string = "No help yet.")
  call param_register(cli_params, "atom_1", "0.0 0.0 0.0", atom_1, help_string = "No help yet.")
  call param_register(cli_params, "atom_2", "0.0 0.0 0.0", atom_2, help_string = "No help yet.", has_value_target=has_atom_2)
  call param_register(cli_params, "Z_1", PARAM_MANDATORY, Z_1, help_string = "No help yet.")
  call param_register(cli_params, "Z_2", "0", Z_2, help_string = "No help yet.", has_value_target=has_atom_2)
  call param_register(cli_params, 'param_file', 'quip_params.xml', param_file, help_string="No help yet.")
  call param_register(cli_params, 'init_args', '', init_args, help_string="No help yet.")
  call param_register(cli_params, "n_configs", "60", n_configs, help_string = "No help yet.")
  call param_register(cli_params, "e0", PARAM_MANDATORY, e0, help_string = "No help yet.")
  call param_register(cli_params, "temp", PARAM_MANDATORY, t, help_string = "No help yet.")
  call param_register(cli_params, "lattice_delta", PARAM_MANDATORY, lattice_delta, help_string = "No help yet.")
  call param_register(cli_params, "atom_delta", "0.0", atom_delta, help_string = "No help yet.")
  call param_register(cli_params, "m_max", PARAM_MANDATORY, m_max, help_string = "No help yet.")
  call param_register(cli_params, "rng_seed", "-1", rng_seed, help_string = "No help yet.")
  call param_register(cli_params, "init_d", "1", init_d, help_string = "No help yet.")

  if (.not. param_read_args(cli_params, task="slice_sample CLI arguments")) then
    call print("Usage: slice_sample lattice='{r r r r r r r r r}' [atom_1={0.0 0.0 0.0}] [atom_2='{r r r}']",PRINT_ALWAYS)
    call print("                    Z_1='i' [Z_2='i'] [param_file='file'] [init_args='str']", PRINT_ALWAYS)
    call print("                    [n_configs=60] e0='r' temp='r' lattice_delta='r' [atom_delta='r'] m_max=i", PRINT_ALWAYS)
    call print("                    [rng_seed=i] [init_d=1]", PRINT_ALWAYS)
    call system_abort("Confused by CLI arguments.")
  end if
  call finalise(cli_params)

  if (rng_seed >= 0) call system_set_random_seeds(rng_seed)

  if (trim(init_args) .ne. '') then
     if (trim(param_file) .ne. '') then
        call potential_filename_initialise(pot, args_str=init_args, param_filename=param_file)
     else
        call initialise(pot, args_str=init_args)
     end if
  else
     call system_abort("Potential not initialised.")
  end if

  lattice = reshape(lattice_vec, (/ 3, 3 /))

  call print("Input lattice:")
  call print(lattice)

  if (has_atom_2) atom_2 = atom_2 - atom_1
  atom_1 = (/ 0.0_dp, 0.0_dp, 0.0_dp /)

  if (lattice(2,1) .fne. 0.0_dp) then
     theta = - atan2(lattice(2,1), lattice(1,1))
     lattice = matmul(reshape((/ cos(theta), sin(theta), 0.0_dp, -sin(theta), cos(theta), 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp /), (/ 3, 3 /)), lattice)
     if (has_atom_2) atom_2 = matmul(reshape((/ cos(theta), sin(theta), 0.0_dp, -sin(theta), cos(theta), 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp /), (/ 3, 3 /)), atom_2)
  end if

  if (lattice(3,1) .fne. 0.0_dp) then
     theta = - atan2(lattice(3,1), lattice(1,1))
     lattice = matmul(reshape((/ cos(theta), 0.0_dp, -sin(theta), 0.0_dp, 1.0_dp, 0.0_dp, sin(theta), 0.0_dp, cos(theta) /), (/ 3, 3 /)), lattice)
     if (has_atom_2) atom_2 = matmul(reshape((/ cos(theta), 0.0_dp, -sin(theta), 0.0_dp, 1.0_dp, 0.0_dp, sin(theta), 0.0_dp, cos(theta) /), (/ 3, 3 /)), atom_2)
  end if

  if (lattice(3,2) .fne. 0.0_dp) then
     theta = - atan2(lattice(3,2), lattice(2,2))
     lattice = matmul(reshape((/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, cos(theta), sin(theta), 0.0_dp, -sin(theta), cos(theta) /), (/ 3, 3 /)), lattice)
     if (has_atom_2) atom_2 = matmul(reshape((/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, cos(theta), sin(theta), 0.0_dp, -sin(theta), cos(theta) /), (/ 3, 3 /)), atom_2)
  end if

  call print("Rotated lattice:")
  call print(lattice)

  if (has_atom_2) then
     allocate(x(0:9), x_dash(0:9))
  else
     allocate(x(0:6), x_dash(0:6))
  end if

  x(1:6) = (/ lattice(1,1), lattice(1,2), lattice(2,2), lattice(1,3), lattice(2,3), lattice(3,3) /)
  if (has_atom_2) x(7:9) = (/ atom_2(1), atom_2(2), atom_2(3) /)

  x(0) = density_function(x(1:))

  ! hack to output castep data
  if (associated(pot%simple%filepot)) then
     if (has_atom_2) then
        call system_command("tail -n 4 " // trim(pot%simple%filepot%filename) // "/" // trim(pot%simple%filepot%filename) // "_output >> slice_sample.xyz")
     else
        call system_command("tail -n 3 " // trim(pot%simple%filepot%filename) // "/" // trim(pot%simple%filepot%filename) // "_output >> slice_sample.xyz")
     end if
  else
     if (has_atom_2) then
        call initialise(at, 2, reshape((/ x(1), 0.0_dp, 0.0_dp, x(2), x(3), 0.0_dp, x(4), x(5), x(6) /), (/3, 3/)))
        at%pos(:,1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
        at%pos(:,2) = (/ x(7), x(8), x(9) /)
        call set_atoms(at, (/ Z_1, Z_2 /))
        call map_into_cell(at)
     else
        call initialise(at, 1, reshape((/ x(1), 0.0_dp, 0.0_dp, x(2), x(3), 0.0_dp, x(4), x(5), x(6) /), (/3, 3/)))
        at%pos(:,1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
        call set_atoms(at, Z_1)
     end if
     call write(at, 'stdout', prefix='AT')
     call finalise(at)
  endif

  call print("SLICE_SAMPLE " // x)

  n = 1
  d = init_d

  do
     x_dash = x
     x = x_new(x_dash(1:), d, lattice_delta, atom_delta, m_max, f_0 = x_dash(0))

     ! hack to output castep data
     if (associated(pot%simple%filepot)) then
        if (has_atom_2) then
           call system_command("tail -n 4 " // trim(pot%simple%filepot%filename) // "/" // trim(pot%simple%filepot%filename) // "_output >> slice_sample.xyz")
        else
           call system_command("tail -n 3 " // trim(pot%simple%filepot%filename) // "/" // trim(pot%simple%filepot%filename) // "_output >> slice_sample.xyz")
        end if
     else
        if (has_atom_2) then
           call initialise(at, 2, reshape((/ x(1), 0.0_dp, 0.0_dp, x(2), x(3), 0.0_dp, x(4), x(5), x(6) /), (/3, 3/)))
           at%pos(:,1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
           at%pos(:,2) = (/ x(7), x(8), x(9) /)
           call set_atoms(at, (/ Z_1, Z_2 /))
           call map_into_cell(at)
        else
           call initialise(at, 1, reshape((/ x(1), 0.0_dp, 0.0_dp, x(2), x(3), 0.0_dp, x(4), x(5), x(6) /), (/3, 3/)))
           at%pos(:,1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
           call set_atoms(at, Z_1)
        end if
        call write(at, 'stdout', prefix='AT')
        call finalise(at)
     endif

     call print("SLICE_SAMPLE " // x)

     if (n >= n_configs) then
        exit
     else
        n = n + 1
        d = d + 1
        if (has_atom_2) then
           if (d == 10) d = 1
        else
           if (d == 7) d = 1
        end if
     end if
  end do

  deallocate(x, x_dash)

  call finalise(pot)

  call system_finalise()

contains

  function density_function(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: density_function

    type(Atoms) :: at
    real(dp) :: E
    integer :: error

    if (size(x) == 9) then
       call initialise(at, 2, reshape((/ x(1), 0.0_dp, 0.0_dp, x(2), x(3), 0.0_dp, x(4), x(5), x(6) /), (/3, 3/)))
       at%pos(:,1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
       at%pos(:,2) = (/ x(7), x(8), x(9) /)
       call set_atoms(at, (/ Z_1, Z_2 /))
       call map_into_cell(at)
    else
       call initialise(at, 1, reshape((/ x(1), 0.0_dp, 0.0_dp, x(2), x(3), 0.0_dp, x(4), x(5), x(6) /), (/3, 3/)))
       at%pos(:,1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
       call set_atoms(at, Z_1)
    end if

    call set_cutoff(at, cutoff(pot))
    call calc_connect(at)
    call calc(pot, at, energy = E, error = error)

    if (error /= ERROR_NONE) E = e0 + 1000000000000.0_dp 

    call finalise(at)

    density_function = exp( - (E - e0) / (t * 0.00008617343_dp) )

  end function density_function

  function x_new(x_0, d, w_1, w_2, m, f_0)
    real(dp), intent(in) :: x_0(:)
    integer, intent(in) :: d
    real(dp), intent(in) :: w_1, w_2
    integer, intent(in) :: m
    real(dp), intent(in), optional :: f_0
    real(dp) :: x_new(0:size(x_0))

    real(dp) :: x(size(x_0))
    real(dp) :: f_slice
    real(dp) :: x_new_L, x_new_R, f_x_new_L, f_x_new_R, f_x_new
    logical :: increase_slice, f_x_new_L_up_to_date, f_x_new_R_up_to_date
    integer :: J, K

    x_new(0) = 0.0_dp
    x_new(1:size(x_0)) = x_0
    x = x_0

    if (present(f_0)) then
       f_slice = ran_uniform() * f_0
    else
       f_slice = ran_uniform() * density_function(x_0)
    end if

    if (d < 7) then
       x_new_L = x_0(d) - (w_1 * ran_uniform())
       x_new_R = x_new_L + w_1
    else
       x_new_L = x_0(d) - (w_2 * ran_uniform())
       x_new_R = x_new_L + w_2
    end if

    J = floor(real(m, dp) * ran_uniform())
    K = m - 1 - J

    f_x_new_L = 0.0_dp
    f_x_new_R = 0.0_dp
    f_x_new_L_up_to_date = .false.
    f_x_new_R_up_to_date = .false.

    do
       if (J > 0) then
          increase_slice = .false.

          if (.not. f_x_new_L_up_to_date) then
             x(d) = x_new_L
             f_x_new_L = density_function(x)
             f_x_new_L_up_to_date = .true.
          end if

          if (f_x_new_L > f_slice) increase_slice = .true.

          if (increase_slice) then
             if (d < 7) then
                x_new_L = x_new_L - w_1
             else
                x_new_L = x_new_L - w_2
             end if
             f_x_new_L_up_to_date = .false.

             J = J - 1
          else
             exit
          end if
       else
          exit
       end if
    end do

    do
       if (K > 0) then
          increase_slice = .false.

          if (.not. f_x_new_R_up_to_date) then
             x(d) = x_new_R
             f_x_new_R = density_function(x)
             f_x_new_R_up_to_date = .true.
          end if

          if (f_x_new_R > f_slice) increase_slice = .true.

          if (increase_slice) then
             if (d < 7) then
                x_new_R = x_new_R + w_1
             else
                x_new_R = x_new_R + w_2
             end if
             f_x_new_R_up_to_date = .false.

             K = K - 1
          else
             exit
          end if
       else
          exit
       end if
    end do

    do
       x_new(d) = x_new_L + (ran_uniform() * (x_new_R - x_new_L))

       f_x_new = density_function(x_new(1:))

       if (f_x_new > f_slice) then
          x_new(0) = f_x_new
          exit
       end if

       if (x_new(d) < x_0(d)) then
          x_new_L = x_new(d)
       else
          x_new_R = x_new(d)
       end if
    end do

  end function x_new

end program slice_sample
