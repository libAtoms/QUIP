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

!% A simple Gaussian process module for 1-D functions learned from samples of
!% their values

module gp_basic_module
use system_module
use linearalgebra_module
implicit none

public :: gp_basic, initialise, finalise, predict, predict_var, predict_grad

type gp_basic
   real(dp) :: l2, f_var
   integer  :: n_f, n_teach
   type(LA_Matrix), allocatable :: Cmat(:)
   real(dp), allocatable :: f_r(:), Cmat_inv_v(:), k(:), Cmat_inv_k(:)
   logical :: initialised = .false.
end type gp_basic

!% initialise (and teach) a gp_basic
interface initialise
   module procedure gp_basic_initialise
end interface initialise

!% finalise and deallocate a gp_basic
interface finalise
   module procedure gp_basic_finalise
end interface finalise

!% predict a function value from a gp
interface predict
   module procedure predict_r, predict_rr
end interface predict

!% predict a gradient of a function value from a gp
interface predict_grad
   module procedure predict_grad_r
end interface predict_grad

!% predict a function variance from a gp
interface predict_var
   module procedure predict_var_r, predict_var_rr
end interface predict_var

interface ff_kernel
   module procedure ff_kernel_r_r, ff_kernel_r_rr
end interface ff_kernel

interface df_kernel
   module procedure df_kernel_r_r, df_kernel_r_rr
end interface df_kernel

contains

subroutine gp_basic_initialise(self, f_r, f_v, f_n, len_scale, f_var, error)
   type(gp_basic), intent(inout) :: self !% object to store GP
   real(dp), intent(in) :: f_r(:), f_v(:), f_n(:) !% arrays of function positions, values, noise 
   real(dp), intent(in) :: len_scale, f_var !% length scale and variance prior for GP
   integer, optional, intent(out) :: error !% error status

   real(dp), allocatable :: Cmat(:,:)
   integer :: i

print *, "init 00"
   INIT_ERROR(error)

   call finalise(self)

   if (size(f_r) /= size(f_v) .or. size(f_r) /= size(f_n)) then
      RAISE_ERROR("size(f_r)="//size(f_r)//" size(f_v)="//size(f_v)//" size(f_n)="//size(f_n)//" not all equal", error)
   endif
   if (len_scale <= 0.0_dp) then
      RAISE_ERROR("invalid len_scale="//len_scale, error)
   endif
   if (f_var <= 0.0_dp) then
      RAISE_ERROR("invalid f_var="//f_var, error)
   endif

   self%l2 = len_scale**2
   self%f_var = f_var

   self%n_f = size(f_r)
   allocate(self%f_r(self%n_f))
   self%f_r(:) = f_r(:)

   self%n_teach = self%n_f

   allocate(Cmat(self%n_teach, self%n_teach))
   do i=1, self%n_f
      Cmat(i,:) = ff_kernel(self%f_r(i), self%f_r(:), self%f_var, self%l2)
      Cmat(i,i) = Cmat(i,i) + f_n(i)
   end do
   allocate(self%Cmat(1))
   call initialise(self%Cmat(1), Cmat)
   deallocate(Cmat)

   allocate(self%Cmat_inv_v(self%n_teach))
   call Matrix_QR_solve(self%Cmat(1), f_v, self%Cmat_inv_v)

   allocate(self%k(self%n_teach))
   allocate(self%Cmat_inv_k(self%n_teach))

   self%initialised = .true.
print *, "init 100"
end subroutine gp_basic_initialise

subroutine gp_basic_finalise(self)
   type(gp_basic), intent(inout) :: self !% object for GP
print *, "fin 00"

   if (self%initialised) then
      if (allocated(self%f_r)) deallocate(self%f_r)
      if (allocated(self%Cmat_inv_v)) deallocate(self%Cmat_inv_v)
      if (allocated(self%k)) deallocate(self%k)
      if (allocated(self%Cmat_inv_k)) deallocate(self%Cmat_inv_k)
      if (allocated(self%Cmat)) then
	 call finalise(self%Cmat(1))
	 deallocate(self%Cmat)
      endif
   endif
   self%l2 = 0.0_dp
   self%f_var = 0.0_dp
   self%n_f = 0
   self%n_teach = 0
   self%initialised = .false.
print *, "fin 100"
end subroutine gp_basic_finalise

function predict_r(self, r)
   type(gp_basic), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r !% position at which to predict value
   real(dp) :: predict_r

   if (.not. self%initialised) then
      predict_r = 0.0_dp
      return
   endif

   self%k(:) = ff_kernel(r, self%f_r(:), self%f_var, self%l2)
   predict_r = dot_product(self%k, self%Cmat_inv_v)
end function predict_r

function predict_rr(self, r)
   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   real(dp) :: predict_rr(size(r))

   integer :: i

   if (.not. self%initialised) then
      predict_rr = 0.0_dp
      return
   endif

   do i=1, size(r)
      self%k(:) = ff_kernel(r(i), self%f_r(:), self%f_var, self%l2)
      predict_rr(i) = dot_product(self%k, self%Cmat_inv_v)
   end do
end function predict_rr

function predict_var_r(self, r)
   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r
   real(dp) :: predict_var_r

   if (.not. self%initialised) then
      predict_var_r = 0.0_dp
      return
   endif

   self%k(:) = ff_kernel(r, self%f_r(:), self%f_var, self%l2)
   call Matrix_QR_Solve(self%Cmat(1), self%k, self%Cmat_inv_k)
   predict_var_r = self%f_var - dot_product(self%k, self%Cmat_inv_k)
end function predict_var_r

function predict_var_rr(self, r)
   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   real(dp) :: predict_var_rr(size(r))

   integer :: i

   if (.not. self%initialised) then
      predict_var_rr = 0.0_dp
      return
   endif

   do i=1, size(r)
      self%k(:) = ff_kernel(r(i), self%f_r(:), self%f_var, self%l2)
      call Matrix_QR_Solve(self%Cmat(1), self%k, self%Cmat_inv_k)
      predict_var_rr(i) = self%f_var - dot_product(self%Cmat_inv_k, self%Cmat_inv_k)
   end do
end function predict_var_rr

function predict_grad_r(self, r)
   type(gp_basic), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r !% position at which to predict value
   real(dp) :: predict_grad_r

   if (.not. self%initialised) then
      predict_grad_r = 0.0_dp
      return
   endif

   self%k(:) = df_kernel(r, self%f_r(:), self%f_var, self%l2)
   predict_grad_r = dot_product(self%k, self%Cmat_inv_v)
end function predict_grad_r


function ff_kernel_r_r(x1, x2, f_var, l2)
   real(dp), intent(in) :: x1, x2, f_var, l2
   real(dp) :: ff_kernel_r_r

   ff_kernel_r_r = f_var*exp(-0.5_dp*((x2-x1)**2)/l2)
end function ff_kernel_r_r

function ff_kernel_r_rr(x1, x2, f_var, l2)
   real(dp), intent(in) :: x1, x2(:), f_var, l2
   real(dp) :: ff_kernel_r_rr(size(x2))

   ff_kernel_r_rr(:) = f_var*exp(-0.5_dp*((x2(:)-x1)**2)/l2)
end function ff_kernel_r_rr

function df_kernel_r_r(x1, x2, f_var, l2)
   real(dp), intent(in) :: x1, x2, f_var, l2
   real(dp) :: df_kernel_r_r

   df_kernel_r_r = (x2-x1)/l2 * ff_kernel(x1, x2, f_var, l2)
end function df_kernel_r_r

function df_kernel_r_rr(x1, x2, f_var, l2)
   real(dp), intent(in) :: x1, x2(:), f_var, l2
   real(dp) :: df_kernel_r_rr(size(x2))

   df_kernel_r_rr(:) = (x2(:)-x1)/l2 * ff_kernel(x1, x2, f_var, l2)
end function df_kernel_r_rr

end module gp_basic_module
