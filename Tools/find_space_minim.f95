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
use iso_c_binding, only : c_ptr, c_loc, c_f_pointer
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
  type(c_ptr) :: at_ptr
  type(Atoms), pointer :: at
  real(dp) :: r_mag

  if (.not. present(data)) call system_abort("find_space_minim_mod func needs data")
  if (size(x) /= 3) call system_abort("find_space_minim_mod func needs x(1:3)")

  at_ptr = transfer(data, at_ptr)
  call c_f_pointer(at_ptr, at)
  call atoms_repoint(at)

  func = 0.0_dp
  do i=1, at%N
    r_mag = distance_min_image(at, i, x(1:3))
    if (r_mag < rad) then
      func = func + rep_E(r_mag)
    endif
  end do

end function func

function dfunc(x,data)
  real(dp) :: x(:)
  character, optional :: data(:)
  real(dp) :: dfunc(size(x))

  integer :: i
  type(c_ptr) :: at_ptr
  type(Atoms), pointer :: at
  real(dp) :: r(3), r_mag

  if (.not. present(data)) call system_abort("find_space_minim_mod func needs data")
  if (size(x) /= 3) call system_abort("find_space_minifind_space_minimeds x(1:3)")

  at_ptr = transfer(data, at_ptr)
  call c_f_pointer(at_ptr, at)
  call atoms_repoint(at)

  dfunc = 0.0_dp
  do i=1, at%N
    r = diff_min_image(at, i, x(1:3))
    r_mag = norm(r)
    if (r_mag < rad) then
      dfunc = dfunc + rep_dE(r, r_mag)
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
  type(c_ptr) :: at_ptr
  type(Atoms), pointer :: at
  real(dp) :: r(3), r_mag

  INIT_ERROR(error)

  if (.not. present(data)) call system_abort("find_space_minim_mod func needs data")
  if (size(x) /= 3) call system_abort("find_space_minim_mod bothfunc needs x(1:3)")

  at_ptr = transfer(data, at_ptr)
  call c_f_pointer(at_ptr, at)
  call atoms_repoint(at)

  f = 0.0_dp
  E = 0.0_dp
  do i=1, at%N
    r = diff_min_image(at, i, x(1:3))
    r_mag = norm(r)
    if (r_mag < rad) then
      E = E + rep_E(r_mag)
      f = f + rep_dE(r, r_mag)
    endif
  end do

end subroutine bothfunc

function rep_E(r_mag)
   real(dp) :: r_mag
   real(dp) :: rep_E

   real(dp) :: f_envelope

   if (r_mag > rad) then
      f_envelope = 0.0_dp
   else if (r_mag < rad-width) then
      f_envelope = 1.0_dp
   else
      f_envelope = 0.5_dp+0.5_dp*cos((r_mag-(rad-width))/width * 3.14159_dp)
   endif

   rep_E = (2.0_dp - r_mag/rad) * f_envelope

end function rep_E

function rep_dE(r, r_mag)
   real(dp) :: r(:), r_mag
   real(dp) :: rep_dE(size(r))

   real(dp) :: f_envelope, df_envelope

   if (r_mag > rad) then
      f_envelope = 0.0_dp
      df_envelope = 0.0_dp
   else if (r_mag < rad-width) then
      f_envelope = 1.0_dp
      df_envelope = 0.0_dp
   else
      f_envelope = 0.5_dp+0.5_dp*cos((r_mag-(rad-width))/width * 3.14159_dp)
      df_envelope = -0.5_dp*sin((r_mag-(rad-width))/width * 3.14159_dp) * (1.0_dp/width) * 3.14159_dp
   endif

   rep_dE = ((2.0_dp - r_mag/rad) * df_envelope  - 1.0_dp/rad * f_envelope) * r/r_mag
end function rep_dE

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
use find_space_minim_mod, only : func, bothfunc, apply_precond_dummy, hook, rad, width
use iso_c_binding, only : c_ptr, c_loc, c_f_pointer
implicit none

  real(dp) :: r(3), prev_r(3), orig_r(3)
  integer :: closest_list(100)
  type(Atoms) :: at
  type(Atoms), target :: closest_at
  type(c_ptr) :: closest_at_ptr
  character(len=128) :: arg
  integer :: data_size
  character, allocatable :: data(:)
  integer :: i_r, n_iter
  real(dp) :: prev_val, final_val, prev_rad
  integer :: i_at
real(dp) :: vi, vf
  real(dp) :: expected_red
  real(dp) :: rad_min, rad_increment

  call system_initialise()
  call read(at, "stdin")
  call get_cmd_arg(1, arg)
  read (unit=arg, fmt=*) r(1)
  call get_cmd_arg(2, arg)
  read (unit=arg, fmt=*) r(2)
  call get_cmd_arg(3, arg)
  read (unit=arg, fmt=*) r(3)
  call get_cmd_arg(4, arg)
  read (unit=arg, fmt=*) rad_min
  call get_cmd_arg(5, arg)
  read (unit=arg, fmt=*) rad_increment

  do i_at=1, at%N
    at%pos(1:3,i_at) = at%pos(1:3,i_at) - r(1:3)
  end do
  call map_into_cell(at)
  orig_r = r
  r = 0.0_dp

  if (at%N < 100) then
     call find_closest(at, r, closest_list(1:at%N))
  else
     call find_closest(at, r, closest_list)
  endif
  call select(closest_at, at, list=closest_list(1:at%N))

  closest_at_ptr = c_loc(closest_at)
  data_size = size(transfer(closest_at_ptr, data))
  allocate(data(data_size))
  data = transfer(closest_at_ptr,data)
  width = 0.2_dp

  i_r = 1
  prev_val = 0.0_dp
  final_val = 0.0_dp
  prev_r = 0.0_dp
  do while (final_val <= 1.0e-3 .and. i_r < 1000)
    prev_r = r
    prev_val = final_val
    rad = rad_min+(i_r-1)*rad_increment
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
