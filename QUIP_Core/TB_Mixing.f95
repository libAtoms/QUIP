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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X TB_Mixing module
!X
!% do charge mixing for self-consistent TB models
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module tb_mixing_module

use System_module
use linearalgebra_module
implicit none

private

integer:: n_hist = 20

! Broyden
real(dp), allocatable:: u(:,:), dF(:,:), last_F(:), last_n(:)
real(dp), allocatable:: dn_mm1(:)
real(dp), allocatable:: a(:,:), beta(:,:), gamma(:,:), c(:,:), t(:,:), w(:)

! Ridders
real(dp):: xmin, fmin, xmax, fmax
real(dp):: bracket_step
logical:: need_bracket = .true.
real(dp)::  x3, f3, x4, f4
logical:: got_f3, got_f4

public :: do_mix_broyden
interface do_mix_broyden
  module procedure do_mix_broyden_array, do_mix_broyden_scalar
end interface do_mix_broyden

public :: do_mix_simple
interface do_mix_simple
  module procedure do_mix_simple_array, do_mix_simple_scalar
end interface do_mix_simple

public :: do_ridders_residual

contains

! From Johnson PRB vol. 38, p. 12807
subroutine do_mix_broyden_array(iter, n, F, np1, alpha, w0)
  integer, intent(in) :: iter
  real(dp), intent(in) :: n(:), F(:)
  real(dp), intent(out) :: np1(:)
  real(dp), intent(in) :: alpha, w0

  integer mi, i, j, ki, li, ni
  integer :: iter_offset
  integer :: effective_iter
  real(dp) :: dF_norm

  iter_offset = 0

  if (realloc_hist_dep_stuff()) then
    if (iter /= 1) then
      call print("WARNING: had to realloc hist_dep_stuff but iter /= 1")
      iter_offset = iter-1
    endif
  endif
  if (realloc_size_dep_stuff(size(n))) then
    if (iter /= 1) then
      call print("WARNING:had to realloc size_dep_stuff but iter /= 1")
      iter_offset = iter-1
    endif
  endif

  effective_iter = iter - iter_offset

  mi = min(effective_iter, n_hist)

  np1 = n + alpha * F

  if (effective_iter == 1) then
    last_F = F
    last_n = n
  else
      if (effective_iter .gt. n_hist) then
	  do i=1, mi-2
	      dF(:,i) = dF(:,i+1)
	      w(i) = w(i+1)
	      u(:,i) = u(:,i+1)
	  end do
	  do i=1, mi-2
	  do j=1, mi-2
	      a(i,j) = a(i+1,j+1)
	  end do
	  end do
	  do i=1, mi-1
	  do j=1, mi-1
	      c(i,j) = c(i,j+1)
	      gamma(i,j) = gamma(i+1,j)
	  end do
	  end do
      end if ! iter > n_hist

      ! calculate dF_(m-1)
      dF(:,mi-1) = F(:) - last_F(:)
      !!!dF_norm = sqrt(sum(dF**2))
      dF_norm = sqrt(sum(dF(:,mi-1)**2))
      dF(:,mi-1) = dF(:,mi-1) / dF_norm
      last_F = F

!! dF(:,mi-1)
!! w(mi-1)
!! a(1:mi-1,mi-1) a(mi-1,1:mi-1)
!! u(1:N_A,mi-1)
!! c(1:mi-1,mi)
!! gamma(mi,1:mi-1)

      ! calculate dn_(m-1)
      dn_mm1(:) = (n(:) - last_n(:)) / dF_norm
      last_n = n

      ! as per DD Johnson's suggestion in paper
      w(mi-1) = min(1.0_dp/sqrt(sum(F**2)), 1.0_dp)

      ! update a_ij
      do i=1, mi-1
	  a(i,mi-1) = w(i)*w(mi-1)*sum(dF(:,i)*dF(:,mi-1))
	  if (i .ne. mi-1) then
	      a(mi-1,i) = a(i,mi-1)
	  end if
      end do

      ! compute beta = (w0^2 + a ) ^ -1
      t = a
      do i=1, mi-1
	  t(i,i) = t(i,i) + w0**2
      end do
      call inverse(t(1:mi-1,1:mi-1), beta(1:mi-1,1:mi-1), .false.)

      ! compute u_(m-1)
      u(:,mi-1) = alpha*dF(:,mi-1) + dn_mm1(:)

      ! compute c_k^m
      do ki=1, mi-1
	  c(ki,mi) = w(ki)*sum(dF(:,ki)*F(:))
      end do
      do li=1, mi-1
	  gamma(mi,li) = 0.0_dp
	  do ki=1, mi-1
	      gamma(mi,li) = gamma(mi,li) + c(ki,mi)*beta(ki,li)
	  end do
      end do
      do ni=1, mi-1
	  np1(:) = np1(:) - w(ni)*gamma(mi,ni)*u(:,ni)
      end do
  endif

end subroutine do_mix_broyden_array

subroutine do_mix_broyden_scalar(iter, n, F, np1, alpha, w0)
  integer, intent(in) :: iter
  real(dp), intent(in) :: n, F
  real(dp), intent(out) :: np1
  real(dp), intent(in) :: alpha, w0

  real(dp) :: na(1), Fa(1), np1a(1)

  na = n
  Fa = F
  np1a = np1

  call do_mix_broyden(iter, na, Fa, np1a, alpha, w0)

  np1 = np1a(1)

end subroutine do_mix_broyden_scalar

subroutine do_mix_simple_array(iter, n, F, np1, alpha)
  integer, intent(in) :: iter
  real(dp), intent(in) :: n(:), F(:)
  real(dp), intent(out) :: np1(:)
  real(dp), intent(in) :: alpha

  np1 = n + alpha * F

end subroutine do_mix_simple_array

subroutine do_mix_simple_scalar(iter, n, F, np1, alpha)
  integer, intent(in) :: iter
  real(dp), intent(in) :: n, F
  real(dp), intent(out) :: np1
  real(dp), intent(in) :: alpha

  real(dp) :: na(1), Fa(1), np1a(1)

  na = n
  Fa = F
  np1a = np1

  call do_mix_simple(iter, na, Fa, np1a, alpha)

  np1 = np1a(1)

end subroutine do_mix_simple_scalar

function realloc_hist_dep_stuff()
  logical :: realloc_hist_dep_stuff

  realloc_hist_dep_stuff = .false.

  if (allocated(a)) then
    if (size(a,1) /= n_hist .or. size(a,2) /= n_hist) then
      deallocate(a)
    endif
  endif
  if (.not. allocated(a)) then
    allocate(a(n_hist,n_hist))
    realloc_hist_dep_stuff = .true.
    a = 0.0_dp
  endif

  if (allocated(beta)) then
    if (size(beta,1) /= n_hist .or. size(beta,2) /= n_hist) then
      deallocate(beta)
    endif
  endif
  if (.not. allocated(beta)) then
    allocate(beta(n_hist,n_hist))
    realloc_hist_dep_stuff = .true.
    beta = 0.0_dp
  endif

  if (allocated(gamma)) then
    if (size(gamma,1) /= n_hist .or. size(gamma,2) /= n_hist) then
      deallocate(gamma)
    endif
  endif
  if (.not. allocated(gamma)) then
    allocate(gamma(n_hist,n_hist))
    realloc_hist_dep_stuff = .true.
  endif
  gamma = 0.0_dp

  if (allocated(c)) then
    if (size(c,1) /= n_hist .or. size(c,2) /= n_hist) then
      deallocate(c)
    endif
  endif
  if (.not. allocated(c)) then
    allocate(c(n_hist,n_hist))
    realloc_hist_dep_stuff = .true.
    c = 0.0_dp
  endif

  if (allocated(t)) then
    if (size(t,1) /= n_hist .or. size(t,2) /= n_hist) then
      deallocate(t)
    endif
  endif
  if (.not. allocated(t)) then
    allocate(t(n_hist,n_hist))
    realloc_hist_dep_stuff = .true.
    t = 0.0_dp
  endif

  if (allocated(w)) then
    if (size(w) /= n_hist) then
      deallocate(w)
    endif
  endif
  if (.not. allocated(w)) then
    allocate(w(n_hist))
    realloc_hist_dep_stuff = .true.
    w = 0.0_dp
  endif

end function realloc_hist_dep_stuff

function realloc_size_dep_stuff(new_n)
  integer, intent(in) :: new_n
  logical :: realloc_size_dep_stuff

  realloc_size_dep_stuff = .false.

  if (allocated(u)) then
    if (size(u,1) /= new_n .or. size(u,2) /= n_hist) then
      deallocate(u)
    endif
  endif
  if (.not. allocated(u)) then
    allocate(u(new_n, n_hist))
    realloc_size_dep_stuff = .true.
    u = 0.0_dp
  endif

  if (allocated(dF)) then
    if (size(dF,1) /= new_n .or. size(dF,2) /= n_hist) then
      deallocate(dF)
    endif
  endif
  if (.not. allocated(dF)) then
    allocate(dF(new_n, n_hist))
    realloc_size_dep_stuff = .true.
    dF = 0.0_dp
  endif

  if (allocated(last_F)) then
    if (size(last_F) /= new_n) then
      deallocate(last_F)
    endif
  endif
  if (.not. allocated(last_F)) then
    allocate(last_F(new_n))
    realloc_size_dep_stuff = .true.
    last_F = 0.0_dp
  endif

  if (allocated(last_n)) then
    if (size(last_n) /= new_n) then
      deallocate(last_n)
    endif
  endif
  if (.not. allocated(last_n)) then
    allocate(last_n(new_n))
    realloc_size_dep_stuff = .true.
    last_n = 0.0_dp
  endif

  if (allocated(dn_mm1)) then
    if (size(dn_mm1) /= new_n) then
      deallocate(dn_mm1)
    endif
  endif
  if (.not. allocated(dn_mm1)) then
    allocate(dn_mm1(new_n))
    realloc_size_dep_stuff = .true.
    dn_mm1 = 0.0_dp
  endif

end function realloc_size_dep_stuff

!subroutine inverse (ns, n, a, b, err)
!  integer, intent(in) :: ns, n
!  real(dp), intent(inout) :: a(ns,ns)
!  real(dp), intent(out) :: b(ns,ns)
!  integer, intent(out) ::  err
!
!  real(dp) :: td(ns), ad(ns), bd(ns)
!
!  integer i, j, k
!
!  do i=1, n
!      if (dabs(a(i,i)) .lt. 1.0D-12) then
!	  print *, "Matrix has zero diagonal element ", i
!	  err = 1
!	  return
!      end if
!  end do
!
!  if (n .eq. 1) then
!      b(1,1) = 1.0_dp/a(1,1)
!      return
!  end if
!
!  b = 0.0_dp
!  do i=1, n
!      b(i,i) = 1.0_dp
!  end do
!
!  do i=1, n
!      do j=1, n
!	  td(j) = a(j,i) / a(i,i)
!      end do
!
!      td(i) = 0.0_dp
!      do k=1, n
!	  bd(k) = b(i,k)
!	  ad(k) = a(i,k)
!      end do
!
!      do k=1, n
!      do j=1, n
!	  b(j,k) = b(j,k) - (td(j)*bd(k))
!	  a(j,k) = a(j,k) - (td(j)*ad(k))
!      end do
!      end do
!
!  end do
!
!  do i=1, n
!  do j=1, n
!      b(j,i) = b(j,i) / a(j,j)
!  end do
!  end do
!
!  err = 0
!
!end subroutine

! From Numerical Recipes in C++, 2nd ed.
subroutine do_ridders_residual(iter, cur_x, cur_f, next_x)
  integer, intent(in) :: iter
  real(dp), intent(in) :: cur_x, cur_f
  real(dp), intent(out) :: next_x

  real(dp) t, num, denom

  if (iter == 1) need_bracket = .true.

  if (need_bracket) then
    if (iter == 1) then
      xmin = cur_x;
      fmin = cur_f;
      next_x = 0.99_dp*cur_x + 0.01_dp*(cur_f+cur_x);
      bracket_step = next_x - cur_x;
      return
    endif
    if (sign(1.0_dp,cur_f) == sign(1.0_dp,fmin)) then
      bracket_step = bracket_step*1.5
      if (abs(cur_f) < abs(fmin)) then
        xmin = cur_x
	fmin = cur_f
        next_x = cur_x + bracket_step
      else
        next_x = cur_x - bracket_step
      endif
      return
    else
      xmax = cur_x
      fmax = cur_f
      if (xmin > xmax) then
        t = xmin; xmin = xmax; xmax = t
        t = fmin; fmin = fmax; fmax = t
      endif
      x3 = (xmin+xmax)/2.0_dp
      next_x = x3
      need_bracket = .false.
      got_f4 = .false.
      got_f3 = .true.
      return
    endif
  endif

  if (got_f3) then
    f3 = cur_f
    num = sign(1.0_dp,fmin-fmax)*f3
    denom = sqrt(f3*f3 - fmin*fmax)
    if (denom == 0.0) then
      next_x = x3
      return
    endif
    x4 = x3 + (x3-xmin)*num/denom
    next_x = x4;
    got_f3 = .false.
    got_f4 = .true.
    return
  endif

  if (got_f4) then
    f4 = cur_f
    if (x3 < x4) then
      if (sign(1.0_dp,fmin) /= sign (1.0_dp,f3)) then
        xmax = x3
        fmax = f3
      else if (sign(1.0_dp,f3) /= sign(1.0_dp,f4)) then
        xmin = x3; fmin = f3;
        xmax = x4; fmax = f4;
      else if (sign(1.0_dp,f4) /= sign(1.0_dp,fmax)) then
        xmin = x4; fmin = f4;
      endif
    else
      if (sign(1.0_dp,fmin) /= sign (1.0_dp,f4)) then
        xmax = x4;
        fmax = f4;
      else if (sign(1.0_dp,f4) /= sign(1.0_dp,f3)) then
        xmin = x4; fmin = f4;
        xmax = x3; fmax = f3;
      else if (sign(1.0_dp,f3) /= sign(1.0_dp,fmax)) then
        xmin = x3; fmin = f3;
      endif
    endif

    x3 = (xmin+xmax)/2.0_dp;
    next_x = x3;
    got_f4 = .false.
    got_f3 = .true.
    return
  endif
end subroutine do_ridders_residual

end module tb_mixing_module
