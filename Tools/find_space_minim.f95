! find largest spherical hole in a configuation starting from some initial guess
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
  real(dp) :: r_mag, offset

  if (.not. present(data)) call system_abort("find_space_minim_mod func needs data")
  if (size(x) /= 3) call system_abort("find_space_minim_mod func needs x(1:3)")

  at = transfer(data, at)
  call atoms_repoint(at)

  offset = 1 + rad*PI/(2.0_dp*width)

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
  real(dp) :: offset
  real(dp) :: r(3), r_mag

  if (.not. present(data)) call system_abort("find_space_minim_mod func needs data")
  if (size(x) /= 3) call system_abort("find_space_minifind_space_minimeds x(1:3)")

  at = transfer(data, at)
  call atoms_repoint(at)

  offset = 1 + rad*PI/(2.0_dp*width)

  dfunc = 0.0_dp
  do i=1, at%N
    r = diff_min_image(at, i, x(1:3))
    r_mag = norm(r)
    if (r_mag < rad) then
      dfunc = dfunc + ( (1.0_dp/r_mag) * exp(-0.1_dp/(rad-r_mag)) * (-0.1_dp/(rad-r_mag)**2) - (1.0_dp/r_mag**2) * exp(-0.1_dp/(rad-r_mag)) ) * r/r_mag
    endif
  end do


end function dfunc

subroutine hook(x,dx,E,done,do_print,data)
  use system_module
  real(dp)::x(:)
  real(dp)::dx(:)
  real(dp)::E
  logical :: done
  logical, optional:: do_print
  character,optional::data(:)

  done = norm(dx) .feq. 0.0_dp

end subroutine hook


subroutine find_closest(at, r, closest_list)
  type(Atoms), intent(in) :: at
  real(dp), intent(in) :: r(3)
  integer, intent(out) :: closest_list(:)

  integer :: n_p, i_p, j
  real(dp) :: prev_dist
  real(dp), allocatable :: r_a(:)
  real(dp), allocatable :: closest_r(:)

  n_p = size(closest_list)

  allocate(r_a(at%N))
  allocate(closest_r(n_p))

  if (at%N < n_p) call system_abort("not enought points ("//at%N//") in atoms object (need "//n_p//")")

  do i_p=1, at%N
    r_a(i_p) = distance_min_image(at, i_p, r)
  end do

  prev_dist = -1e38_dp
  do i_p=1, n_p
    closest_r(i_p) = 1.0e38_dp
    closest_list(i_p) = -1
    do j=1, at%N
      if (r_a(j) < closest_r(i_p) .and. r_a(j) >= prev_dist .and. .not. any(j == closest_list(1:i_p-1))) then
        closest_list(i_p) = j
        closest_r(i_p) = r_a(j)
      endif
    end do
    prev_dist = closest_r(i_p)
  end do

  deallocate(r_a, closest_r)

end subroutine find_closest

end module find_space_minim_mod

program find_space_minim
use libatoms_module
use find_space_minim_mod
implicit none

  real(dp) :: r(3), prev_r(3)
  integer :: closest_list(30)
  type(Atoms) :: at, closest_at
  character(len=128) :: arg
  integer :: data_size
  character, allocatable :: data(:)
  integer :: i_r, n_iter
  real(dp) :: prev_val, final_val, prev_rad

  call system_initialise()
  call read_xyz(at, "stdin")
  call get_cmd_arg(1, arg)
  read (unit=arg, fmt=*) r(1)
  call get_cmd_arg(2, arg)
  read (unit=arg, fmt=*) r(2)
  call get_cmd_arg(3, arg)
  read (unit=arg, fmt=*) r(3)

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
    call verbosity_push(SILENT)
    n_iter = minim(r, func, dfunc, 'cg', 1.0e-6_dp, 1000, 'NR_LINMIN', eps_guess=1e-3_dp, hook=hook, data=data)
    final_val = func(r, data)
    call verbosity_pop()
    call find_closest(at, r, closest_list(1:1))
    call print("rad " // rad// " n_iter " // n_iter // " relaxed val " // func(r, data) // " relaxed r" // r &
      // " radius " // distance_min_image(at,closest_list(1),r))
    i_r = i_r + 1
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
