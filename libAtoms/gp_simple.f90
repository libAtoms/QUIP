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

module gp_simple_module
use system_module
use linearalgebra_module
implicit none

public :: gp_simple, initialise, finalise, predict, predict_var

type gp_simple
   real(dp) :: l2, f_var
   integer  :: n_f, n_teach
   real(dp), allocatable :: Cmat(:,:), f_r(:), Cmat_inv_v(:), k(:), Cmat_inv_k(:)
end type gp_simple

!% initialise (and teach) a gp_simple
interface initialise
   module procedure gp_simple_initialise
end interface initialise

!% finalise and deallocate a gp_simple
interface finalise
   module procedure gp_simple_finalise
end interface finalise

!% predict a function value from a gp
interface predict
   module procedure predict_r, predict_rr
end interface predict

!% predict a function variance from a gp
interface predict_var
   module procedure predict_var_r, predict_var_rr
end interface predict_var

interface ff_kernel
   module procedure ff_kernel_r_r, ff_kernel_r_rr
end interface ff_kernel

contains

subroutine gp_simple_initialise(self, f_r, f_v, f_n, len_scale, f_var, error)
   type(gp_simple), intent(inout) :: self !% object to store GP
   real(dp), intent(in) :: f_r(:), f_v(:), f_n(:) !% arrays of function positions, values, noise 
   real(dp), intent(in) :: len_scale, f_var !% length scale and variance prior for GP
   integer, optional, intent(out) :: error !% error status

   integer :: i

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

   allocate(self%Cmat(self%n_teach, self%n_teach))

   do i=1, self%n_f
      self%Cmat(i,:) = ff_kernel(self%f_r(i), self%f_r(:), self%f_var, self%l2)
      self%Cmat(i,i) = self%Cmat(i,i) + f_n(i)
   end do

   allocate(self%Cmat_inv_v(self%n_teach))
   call symmetric_linear_solve(self%Cmat, f_v, self%Cmat_inv_v)

   allocate(self%k(self%n_teach))
   allocate(self%Cmat_inv_k(self%n_teach))
end subroutine gp_simple_initialise

subroutine gp_simple_finalise(self)
   type(gp_simple), intent(inout) :: self !% object for GP

   self%l2 = 0.0_dp
   self%f_var = 0.0_dp
   self%n_f = 0
   self%n_teach = 0
   if (allocated(self%Cmat)) deallocate(self%Cmat)
   if (allocated(self%f_r)) deallocate(self%f_r)
   if (allocated(self%Cmat_inv_v)) deallocate(self%Cmat_inv_v)
   if (allocated(self%k)) deallocate(self%k)
   if (allocated(self%Cmat_inv_k)) deallocate(self%Cmat_inv_k)
end subroutine gp_simple_finalise

function predict_r(self, r)
   type(gp_simple), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r !% position at which to predict value
   real(dp) :: predict_r

   self%k(:) = ff_kernel(r, self%f_r(:), self%f_var, self%l2)
   predict_r = dot_product(self%k, self%Cmat_inv_v)
end function predict_r

function predict_rr(self, r)
   type(gp_simple), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   real(dp) :: predict_rr(size(r))

   integer :: i

   do i=1, size(r)
      self%k(:) = ff_kernel(r(i), self%f_r(:), self%f_var, self%l2)
      predict_rr(i) = dot_product(self%k, self%Cmat_inv_v)
   end do
end function predict_rr

function predict_var_r(self, r)
   type(gp_simple), intent(inout) :: self
   real(dp), intent(in) :: r
   real(dp) :: predict_var_r

   self%k(:) = ff_kernel(r, self%f_r(:), self%f_var, self%l2)
   call symmetric_linear_solve(self%Cmat, self%k, self%Cmat_inv_k)
   predict_var_r = self%f_var - dot_product(self%k, self%Cmat_inv_k)
end function predict_var_r

function predict_var_rr(self, r)
   type(gp_simple), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   real(dp) :: predict_var_rr(size(r))

   integer :: i

   do i=1, size(r)
      self%k(:) = ff_kernel(r(i), self%f_r(:), self%f_var, self%l2)
      call symmetric_linear_solve(self%Cmat, self%k, self%Cmat_inv_k)
      predict_var_rr(i) = self%f_var - dot_product(self%Cmat_inv_k, self%Cmat_inv_k)
   end do
end function predict_var_rr



function ff_kernel_r_r(x1, x2, f_var, l2)
   real(dp), intent(in) :: x1, x2, f_var, l2
   real(dp) :: ff_kernel_r_r

   ff_kernel_r_r = f_var*exp(-0.5_dp*(x2-x1)**2/l2)
end function ff_kernel_r_r

function ff_kernel_r_rr(x1, x2, f_var, l2)
   real(dp), intent(in) :: x1, x2(:), f_var, l2
   real(dp) :: ff_kernel_r_rr(size(x2))

   ff_kernel_r_rr(:) = f_var*exp(-0.5_dp*(x2(:)-x1)**2/l2)
end function ff_kernel_r_rr

end module gp_simple_module
