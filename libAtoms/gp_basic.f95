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

public :: gp_basic, initialise, finalise, f_predict, f_predict_var, f_predict_grad
public :: SE_kernel_r_rr

type gp_basic
   real(dp) :: l2, f_var
   integer  :: n_f, n_g, n_teach
   type(LA_Matrix), allocatable :: Cmat(:)
   real(dp), allocatable :: f_r(:), g_r(:), Cmat_inv_v(:), k(:), Cmat_inv_k(:)
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

!% f_predict a function value from a gp
interface f_predict
   module procedure f_predict_r
end interface f_predict

!% f_predict a gradient of a function value from a gp
interface f_predict_grad
   module procedure f_predict_grad_r
end interface f_predict_grad

!% f_predict a function variance from a gp
interface f_predict_var
   module procedure f_predict_var_r
end interface f_predict_var

contains

subroutine gp_basic_initialise(self, len_scale, f_var, kernel, f_r, f_v, f_n, g_r, g_v, g_n, sparse_set, error)
   type(gp_basic), intent(inout) :: self !% object to store GP
   real(dp), intent(in) :: len_scale, f_var !% length scale and variance prior for GP
   interface
      subroutine kernel(x1, x2, f_var, l2, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1, x2(:), f_var, l2
	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
      end subroutine kernel
   end interface
   real(dp), optional, intent(in) :: f_r(:), f_v(:), f_n(:) !% arrays of function positions, values, noise 
   real(dp), optional, intent(in) :: g_r(:), g_v(:), g_n(:) !% arrays of function gradient positions, values, noise 
   integer, optional, intent(in) :: sparse_set(:) !% set of points to use for sparsifcation
   integer, optional, intent(out) :: error !% error status

   real(dp), allocatable :: Cmat(:,:), all_v(:)
   integer :: i, i_glob

   INIT_ERROR(error)

   call finalise(self)

   if (len_scale <= 0.0_dp) then
      RAISE_ERROR("invalid len_scale="//len_scale, error)
   endif
   if (f_var <= 0.0_dp) then
      RAISE_ERROR("invalid f_var="//f_var, error)
   endif

   i = count( (/ present(f_r), present(f_v), present(f_n) /) )
   if (i /= 0 .and. i /= 3) then
      RAISE_ERROR("got something other than all or none of f_r, f_v, f_n", error)
   endif
   i = count( (/ present(g_r), present(g_v), present(g_n) /) )
   if (i /= 0 .and. i /= 3) then
      RAISE_ERROR("got something other than all or none of g_r, g_v, g_n", error)
   endif

   if (present(f_r)) then
      if (size(f_r) /= size(f_v) .or. size(f_r) /= size(f_n)) then
	 RAISE_ERROR("size(f_r)="//size(f_r)//" size(f_v)="//size(f_v)//" size(f_n)="//size(f_n)//" not all equal", error)
      endif
   endif
   if (present(g_r)) then
      if (size(g_r) /= size(g_v) .or. size(g_r) /= size(g_n)) then
	 RAISE_ERROR("size(g_r)="//size(g_r)//" size(g_v)="//size(g_v)//" size(g_n)="//size(g_n)//" not all equal", error)
      endif
   endif

   self%l2 = len_scale**2
   self%f_var = f_var

   if (present(sparse_set)) then
      RAISE_ERROR("no sparsification implemented yet", error)
   endif

   if (present(f_r)) then
      self%n_f = size(f_r)
      allocate(self%f_r(self%n_f))
      self%f_r(:) = f_r(:)
   else
      self%n_f = 0
   endif
   if (present(g_r)) then
      self%n_g = size(g_r)
      allocate(self%g_r(self%n_g))
      self%g_r(:) = g_r(:)
   else
      self%n_g = 0
   endif

   if (self%n_f + self%n_g == 0) then
      RAISE_ERROR("no teaching data provided", error)
   endif

   self%n_teach = self%n_f + self%n_g

   allocate(Cmat(self%n_teach, self%n_teach))
   do i=1, self%n_f
      ! f on f
      call kernel(self%f_r(i), self%f_r(:), self%f_var, self%l2, f_f=Cmat(i,1:self%n_f))
      Cmat(i,i) = Cmat(i,i) + f_n(i)
      if (self%n_g > 0) then
	 ! f on g
	 call kernel(self%f_r(i), self%g_r(:), self%f_var, self%l2, f_g=Cmat(i,self%n_f+1:self%n_f+self%n_g))
      endif
   end do
   do i=1, self%n_g
      i_glob = self%n_f + i
      ! g on g
      call kernel(self%g_r(i), self%g_r(:), self%f_var, self%l2, g_g=Cmat(i_glob,self%n_f+1:self%n_f+self%n_g))
      Cmat(i_glob,i_glob) = Cmat(i_glob,i_glob) + g_n(i)
      if (self%n_f > 0) then
	 ! g on f
	 call kernel(self%g_r(i), self%f_r(:), self%f_var, self%l2, g_f=Cmat(i_glob,1:self%n_f))
      endif
   end do

   allocate(self%Cmat(1))
   call initialise(self%Cmat(1), Cmat)
   deallocate(Cmat)

   allocate(self%Cmat_inv_v(self%n_teach))
   allocate(all_v(self%n_teach))
   if (self%n_f > 0) all_v(1:self%n_f) = f_v(:)
   if (self%n_g > 0) all_v(self%n_f+1:self%n_f+self%n_g) = g_v(:)
   call Matrix_QR_solve(self%Cmat(1), all_v, self%Cmat_inv_v)
   deallocate(all_v)

   allocate(self%k(self%n_teach))
   allocate(self%Cmat_inv_k(self%n_teach))

   self%initialised = .true.
end subroutine gp_basic_initialise

subroutine gp_basic_finalise(self)
   type(gp_basic), intent(inout) :: self !% object for GP

   if (self%initialised) then
      if (allocated(self%f_r)) deallocate(self%f_r)
      if (allocated(self%g_r)) deallocate(self%g_r)
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
   self%n_g = 0
   self%n_teach = 0
   self%initialised = .false.
end subroutine gp_basic_finalise

function f_predict_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r !% position at which to f_predict value
   interface
      subroutine kernel(x1, x2, f_var, l2, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1, x2(:), f_var, l2
	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_r

   if (.not. self%initialised) then
      f_predict_r = 0.0_dp
      return
   endif

   if (self%n_f > 0) call kernel(r, self%f_r(1:self%n_f), self%f_var, self%l2, f_f=self%k(1:self%n_f))
   if (self%n_g > 0) call kernel(r, self%g_r(1:self%n_g), self%f_var, self%l2, f_g=self%k(self%n_f+1:self%n_f+self%n_g))
   f_predict_r = dot_product(self%k, self%Cmat_inv_v)
end function f_predict_r

function f_predict_var_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r
   interface
      subroutine kernel(x1, x2, f_var, l2, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1, x2(:), f_var, l2
	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_var_r

   if (.not. self%initialised) then
      f_predict_var_r = 0.0_dp
      return
   endif

   if (self%n_f > 0) call kernel(r, self%f_r(1:self%n_f), self%f_var, self%l2, f_f=self%k(1:self%n_f))
   if (self%n_g > 0) call kernel(r, self%g_r(1:self%n_g), self%f_var, self%l2, f_g=self%k(self%n_f+1:self%n_f+self%n_g))
   call Matrix_QR_Solve(self%Cmat(1), self%k, self%Cmat_inv_k)
   f_predict_var_r = self%f_var - dot_product(self%k, self%Cmat_inv_k)
end function f_predict_var_r

function f_predict_grad_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r !% position at which to f_predict value
   interface
      subroutine kernel(x1, x2, f_var, l2, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1, x2(:), f_var, l2
	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_grad_r

   if (.not. self%initialised) then
      f_predict_grad_r = 0.0_dp
      return
   endif

   if (self%n_f > 0) call kernel(r, self%f_r(1:self%n_f), self%f_var, self%l2, g_f=self%k(1:self%n_f))
   if (self%n_g > 0) call kernel(r, self%g_r(1:self%n_g), self%f_var, self%l2, g_g=self%k(self%n_f+1:self%n_f+self%n_g))
   f_predict_grad_r = dot_product(self%k, self%Cmat_inv_v)
end function f_predict_grad_r

subroutine SE_kernel_r_rr(x1, x2, f_var, l2, f_f, g_f, f_g, g_g)
   real(dp), intent(in) :: x1, x2(:), f_var, l2
   real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)

   if (present(f_f)) f_f(:) = f_var*exp(-0.5_dp*((x2(:)-x1)**2)/l2)
   if (present(g_f)) g_f(:) = (x2(:)-x1)/l2 * f_var*exp(-0.5_dp*((x2(:)-x1)**2)/l2)
   if (present(f_g)) f_g(:) = -(x2(:)-x1)/l2 * f_var*exp(-0.5_dp*((x2(:)-x1)**2)/l2)
   if (present(g_g)) g_g(:) = (1.0_dp-(x1-x2(:))**2/l2)/l2 * f_var*exp(-0.5_dp*((x2(:)-x1)**2)/l2)
end subroutine SE_kernel_r_rr

end module gp_basic_module
