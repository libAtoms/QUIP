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

! find largest spherical hole in a configuation starting from some initial guess
#include "error.inc"
module find_space_minim_mod
use libatoms_module
implicit none

real(dp), save :: rad, width

contains

  ! x > r   fo = 1-cos((x-r)/w * pi/2 + 3 pi /2)
  ! f'(r) = sin(3 pi /2) pi/(2 w) = -pi/(2 w)
  ! x <= r  fi = o - x pi/(2 w)
  ! fi(r) = fo(r)
  ! 1 = o - r pi/(2 w)
  ! o = 1 + r pi / (2 w)

function func(x,data)
  real(dp) :: x(:)
  character, optional :: data(:)
  real(dp) :: func

  integer :: i
  type(Atoms) :: at
  real(dp) :: r_mag

  if (.not. present(data)) call system_abort("find_space_minim_mod func needs data")
  if (size(x) /= 3) call system_abort("find_space_minim_mod func needs x(1:3)")

  at = transfer(data, at)
! call print("BOB 00", PRINT_ALWAYS)
  call atoms_repoint(at)
! call print("JOE 00", PRINT_ALWAYS)

  func = 0.0_dp
  do i=1, at%N
    r_mag = distance_min_image(at, i, x(1:3))
    if (r_mag < rad) then
      func = func + (1.0_dp/r_mag) * exp(-0.1_dp/(rad-r_mag))
    endif
  end do


end function func

function dfunc(x,data)
  real(dp) :: x(:)
  character, optional :: data(:)
  real(dp) :: dfunc(size(x))

  integer :: i
  type(Atoms) :: at
  real(dp) :: r(3), r_mag

  if (.not. present(data)) call system_abort("find_space_minim_mod func needs data")
  if (size(x) /= 3) call system_abort("find_space_minifind_space_minimeds x(1:3)")

  at = transfer(data, at)
! call print("BOB 10", PRINT_ALWAYS)
  call atoms_repoint(at)
! call print("JOE 10", PRINT_ALWAYS)

  dfunc = 0.0_dp
  do i=1, at%N
    r = diff_min_image(at, i, x(1:3))
    r_mag = norm(r)
    if (r_mag < rad) then
      dfunc = dfunc + ( (1.0_dp/r_mag) * exp(-0.1_dp/(rad-r_mag)) * (-0.1_dp/(rad-r_mag)**2) - (1.0_dp/r_mag**2) * exp(-0.1_dp/(rad-r_mag)) ) * r/r_mag
    endif
  end do


end function dfunc

subroutine apply_precond_dummy(x,g,P_g,data,error)
  real(dp) :: x(:), g(:), P_g(:)
  character, optional :: data(:)
  integer, optional :: error

  INIT_ERROR(error)
  P_g = g
end subroutine apply_precond_dummy

subroutine bothfunc(x,E,f,data,error)
  real(dp) :: x(:), E, f(:)
  character, optional :: data(:)
  integer,optional :: error

  integer :: i
  type(Atoms) :: at
  real(dp) :: r(3), r_mag

  INIT_ERROR(error)

  if (.not. present(data)) call system_abort("find_space_minim_mod func needs data")
  if (size(x) /= 3) call system_abort("find_space_minim_mod bothfunc needs x(1:3)")

  at = transfer(data, at)
! call print("BOB 00", PRINT_ALWAYS)
  call atoms_repoint(at)
! call print("JOE 00", PRINT_ALWAYS)

  f = 0.0_dp
  E = 0.0_dp
  do i=1, at%N
    r = diff_min_image(at, i, x(1:3))
    r_mag = norm(r)
    if (r_mag < rad) then
      E = E + (1.0_dp/r_mag) * exp(-0.1_dp/(rad-r_mag))
      f = f + ( (1.0_dp/r_mag) * exp(-0.1_dp/(rad-r_mag)) * (-0.1_dp/(rad-r_mag)**2) - (1.0_dp/r_mag**2) * exp(-0.1_dp/(rad-r_mag)) ) * r/r_mag
    endif
  end do

end subroutine bothfunc

subroutine hook(x,dx,E,done,do_print,data)
  use system_module
  real(dp), intent(in) ::x(:)
  real(dp), intent(in) ::dx(:)
  real(dp), intent(in) ::E
  logical, intent(out) :: done
  logical, optional, intent(in) :: do_print
  character(len=1),optional, intent(in) ::data(:)

  done = norm(dx) .feq. 0.0_dp

end subroutine hook

end module find_space_minim_mod

program find_space_minim
use libatoms_module
use find_space_minim_mod
implicit none

  real(dp) :: r(3), prev_r(3), orig_r(3)
  integer :: closest_list(100)
  type(Atoms) :: at, closest_at
  character(len=128) :: arg
  integer :: data_size
  character, allocatable :: data(:)
  integer :: i_r, n_iter
  real(dp) :: prev_val, final_val, prev_rad
  integer :: i_at
real(dp) :: vi, vf
  real(dp) :: expected_red

  call system_initialise()
  call read(at, "stdin")
  call get_cmd_arg(1, arg)
  read (unit=arg, fmt=*) r(1)
  call get_cmd_arg(2, arg)
  read (unit=arg, fmt=*) r(2)
  call get_cmd_arg(3, arg)
  read (unit=arg, fmt=*) r(3)

  do i_at=1, at%N
    at%pos(1:3,i_at) = at%pos(1:3,i_at) - r(1:3)
  end do
  call map_into_cell(at)
  orig_r = r
  r = 0.0_dp

  call find_closest(at, r, closest_list)
  call select(closest_at, at, list=closest_list)

  data_size = size(transfer(closest_at, data))
  allocate(data(data_size))
  data = transfer(closest_at,data)
  width = 0.2_dp

  i_r = 1
  prev_val = 0.0_dp
  final_val = 0.0_dp
  prev_r = 0.0_dp
  do while (final_val <= 1.0e-3 .and. i_r < 1000)
    prev_r = r
    prev_val = final_val
    rad = i_r*0.1_dp
    call verbosity_push(PRINT_SILENT)
    ! n_iter = minim(r, func, dfunc, 'cg', 1.0e-6_dp, 1000, 'NR_LINMIN', eps_guess=1e-3_dp, hook=hook, data=data)
    expected_red = 1.0e-2_dp
    n_iter = n_minim(r, bothfunc, .false., apply_precond_dummy, vi, vf, expected_red, 1000, 1e-8_dp, hook=hook, data=data)
    final_val = func(r, data)
    call verbosity_pop()
    call find_closest(at, r, closest_list(1:1))
    call print("rad " // rad// " n_iter " // n_iter // " relaxed val " // func(r, data) // " relaxed r" // r &
      // " radius " // distance_min_image(at,closest_list(1),r))
    i_r = i_r + 1
  end do

  r = r + orig_r
  prev_r = prev_r + orig_r
  do i_at=1, at%N
    at%pos(1:3,i_at) = at%pos(1:3,i_at) + orig_r(1:3)
  end do

  if (i_r == 2) then
    call print("WARNING: first try found no space at rad " // rad)
  else
    call find_closest(at, prev_r, closest_list(1:1))
    prev_rad = distance_min_image(at,closest_list(1),prev_r)

    call print("best r " // prev_r // " best_rad " // prev_rad // " best_val " // prev_val // " next_val " // final_val)
  endif
    
  call system_finalise()
end program
