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
   real(dp) :: len_scale_sq, f_var
   integer  :: m_f, m_g, m_teach
   type(LA_Matrix) :: Cmat, Kmm, noise_Kmm
   real(dp), allocatable :: f_r(:), g_r(:), Cmat_inv_v(:), k(:),  mat_inv_k(:), kp(:)
   logical :: sparsified = .false.
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
interface f_predict
   module procedure f_predict_r
end interface f_predict

!% predict a gradient of a function value from a gp
interface f_predict_grad
   module procedure f_predict_grad_r
end interface f_predict_grad

!% predict a function variance from a gp
interface f_predict_var
   module procedure f_predict_var_r
end interface f_predict_var

!% predict the variance of a gradient of a function from a gp
interface f_predict_grad_var
   module procedure f_predict_grad_var_r
end interface f_predict_grad_var

!% predict the gradient of a function variance from a gp
interface f_predict_var_grad
   module procedure f_predict_var_grad_r
end interface f_predict_var_grad

contains

subroutine gp_basic_initialise(self, len_scale, f_var, kernel, f_r, f_v, f_n, g_r, g_v, g_n, f_sparse_set, g_sparse_set, jitter, error)
   type(gp_basic), intent(inout) :: self !% object to store GP
   real(dp), intent(in) :: len_scale, f_var !% length scale and variance prior for GP
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1, x2(:), f_var, len_scale_sq
	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
      end subroutine kernel
   end interface
   real(dp), optional, intent(in) :: f_r(:), f_v(:), f_n(:) !% arrays of function positions, values, noise 
   real(dp), optional, intent(in) :: g_r(:), g_v(:), g_n(:) !% arrays of function gradient positions, values, noise 
   integer, optional, target, intent(in) :: f_sparse_set(:), g_sparse_set(:) !% sets of points to use for sparsifcation for values and gradients
   real(dp), intent(in), optional :: jitter
   integer, optional, intent(out) :: error !% error status

   integer, pointer :: u_f_sparse_set(:), u_g_sparse_set(:)
   real(dp), allocatable :: Cmat(:,:), Kmn(:,:), y(:), siginvsq_y(:), Kmn_siginvsq_y(:), siginvsq_Knm(:,:), Kmm(:,:)
   integer :: i, n_f, n_g, n_tot

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

   ! check and if needed point u_[fg]_sparse_set
   if (present(f_sparse_set)) then
      if (size(f_sparse_set) > 0 .and. .not. present(f_r)) then
	 RAISE_ERROR("got f_sparse_set but no f_r", error)
      endif
      if (any(f_sparse_set < 1) .or. any(f_sparse_set > size(f_r))) then
	 RAISE_ERROR("some of f_sparse_set out of range "//minval(f_sparse_set)//" "//maxval(f_sparse_set)//" "//size(f_r), error)
      endif
      u_f_sparse_set => f_sparse_set
      self%sparsified = .true.
   else
      if (present(f_r)) then
	 allocate(u_f_sparse_set(size(f_r)))
	 do i=1, size(f_r)
	    u_f_sparse_set(i) = i
	 end do
      end if
   endif
   if (present(g_sparse_set)) then
      if (size(g_sparse_set) > 0 .and. .not. present(g_r)) then
	 RAISE_ERROR("got g_sparse_set but no g_r", error)
      endif
      if (any(g_sparse_set < 1) .or. any(g_sparse_set > size(g_r))) then
	 RAISE_ERROR("some of g_sparse_set out of range "//minval(g_sparse_set)//" "//maxval(g_sparse_set)//" "//size(g_r), error)
      endif
      u_g_sparse_set => g_sparse_set
      self%sparsified = .true.
   else
      if (present(g_r)) then
	 allocate(u_g_sparse_set(size(g_r)))
	 do i=1, size(g_r)
	    u_g_sparse_set(i) = i
	 end do
      end if
   endif

   self%len_scale_sq = len_scale**2
   self%f_var = f_var

   if (present(f_r)) then
      self%m_f = size(u_f_sparse_set)
      allocate(self%f_r(self%m_f))
      self%f_r(:) = f_r(u_f_sparse_set)
      n_f = size(f_r)
   else
      self%m_f = 0
      n_f = 0
   endif
   if (present(g_r)) then
      self%m_g = size(u_g_sparse_set)
      allocate(self%g_r(self%m_g))
      self%g_r(:) = g_r(u_g_sparse_set)
      n_g = size(g_r)
   else
      self%m_g = 0
      n_g = 0
   endif

   if (self%m_f + self%m_g == 0) then
      RAISE_ERROR("no teaching data provided", error)
   endif

   self%m_teach = self%m_f + self%m_g
   n_tot = n_f + n_g

   allocate(Kmm(self%m_teach, self%m_teach))

   ! everyone needs Kmm
   call kernel_mat(self%m_f, self%f_r, self%m_g, self%g_r, &
                   self%m_f, self%f_r, self%m_g, self%g_r, &
                   self%f_var, self%len_scale_sq, kernel, Kmm)

   if (self%sparsified) then
      allocate(Kmn(self%m_teach, n_tot))
      allocate(Kmn_siginvsq_y(self%m_teach))
      allocate(siginvsq_Knm(n_tot, self%m_teach))
      allocate(siginvsq_y(n_tot))
      allocate(Cmat(self%m_teach, self%m_teach))

      call kernel_mat(self%m_f, self%f_r, self%m_g, self%g_r, & 
		      n_f, f_r, n_g, g_r, &
		      self%f_var, self%len_scale_sq, kernel, Kmn)

      ! we'll need Kmm^{-1} for variance
      call initialise(self%Kmm, Kmm)

      if (n_f > 0) siginvsq_y(1:n_f) = f_v(:)/f_n(:)
      if (n_g > 0) siginvsq_y(n_f+1:n_f+n_g) = g_v(:)/g_n(:)
      call matrix_product_sub(Kmn_siginvsq_y, Kmn, siginvsq_y)

      siginvsq_Knm = transpose(Kmn)
      do i=1, n_f
	 siginvsq_Knm(i,:) = siginvsq_Knm(i,:) / f_n(i)
      end do
      do i=1, n_g
	 siginvsq_Knm(n_f+i,:) = siginvsq_Knm(n_f+i,:) / g_n(i)
      end do

      Cmat = Kmm
      call matrix_product_sub(Cmat, Kmn, siginvsq_Knm, lhs_factor=1.0_dp, rhs_factor=1.0_dp)
      ! Cmat is now (K_{mm} + K_{mn} \Sigma^{-2} K_{nm})

      ! "jitter" (ABP e-mail 14 Aug)
      if (present(jitter)) then
	 do i=1, size(Cmat,1)
	    Cmat(i,i) = Cmat(i,i) + jitter
	 end do
      endif

      call initialise(self%Cmat, Cmat)
      ! self%Cmat is now (K_{mm} + K_{mn} \Sigma^{-2} K_{nm})


      allocate(self%Cmat_inv_v(self%m_teach))
      call Matrix_QR_Solve(self%Cmat, Kmn_siginvsq_y, self%Cmat_inv_v)

      deallocate(Cmat)
      deallocate(Kmn)
      deallocate(Kmn_siginvsq_y)
      deallocate(siginvsq_Knm)
      deallocate(siginvsq_y)
   else
      allocate(y(n_tot))

      ! we'll need self%noise_Kmm, which is Kmm shifted by noise, for variance
      do i=1, self%m_f
	 Kmm(i,i) = Kmm(i,i) + f_n(u_f_sparse_set(i))
      end do
      do i=1, self%m_g
	 Kmm(self%m_f+i,self%m_f+i) = Kmm(self%m_f+i,self%m_f+i) + g_n(u_g_sparse_set(i))
      end do
      call initialise(self%noise_Kmm, Kmm)

      if (n_f > 0) y(1:n_f) = f_v(:)
      if (n_g > 0) y(n_f+1:n_f+n_g) = g_v(:)

      allocate(self%Cmat_inv_v(self%m_teach))
      call Matrix_QR_Solve(self%noise_Kmm, y, self%Cmat_inv_v)

      deallocate(y)
   endif

   deallocate(Kmm)
   allocate(self%k(self%m_teach))
   allocate(self%kp(self%m_teach))
   allocate(self%mat_inv_k(self%m_teach))

   if (.not. present(f_sparse_set) .and. present(f_r)) deallocate(u_f_sparse_set)
   if (.not. present(g_sparse_set) .and. present(g_r)) deallocate(u_g_sparse_set)

   self%initialised = .true.
end subroutine gp_basic_initialise

subroutine kernel_mat(l_n_f, l_f_r, l_n_g, l_g_r, r_n_f, r_f_r, r_n_g, r_g_r, f_var, len_scale_sq, kernel, mat)
   integer :: l_n_f, l_n_g
   real(dp), optional :: l_f_r(:), l_g_r(:)
   integer :: r_n_f, r_n_g
   real(dp), optional :: r_f_r(:), r_g_r(:)
   real(dp) :: f_var, len_scale_sq
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1, x2(:), f_var, len_scale_sq
	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
      end subroutine kernel
   end interface
   real(dp) :: mat(:,:)

   integer :: i, i_glob

   do i=1, l_n_f
      ! f on f
      if (r_n_f > 0) call kernel(l_f_r(i), r_f_r(:), f_var, len_scale_sq, f_f=mat(i,1:r_n_f))
      ! f on g
      if (r_n_g > 0) call kernel(l_f_r(i), r_g_r(:), f_var, len_scale_sq, f_g=mat(i,r_n_f+1:r_n_f+r_n_g))
   end do
   do i=1, l_n_g
      i_glob = l_n_f + i
      ! g on f
      if (r_n_f > 0) call kernel(l_g_r(i), r_f_r(:), f_var, len_scale_sq, g_f=mat(i_glob,1:r_n_f))
      ! g on g
      if (r_n_g > 0) call kernel(l_g_r(i), r_g_r(:), f_var, len_scale_sq, g_g=mat(i_glob,r_n_f+1:r_n_f+r_n_g))
   end do
end subroutine

subroutine gp_basic_finalise(self)
   type(gp_basic), intent(inout) :: self !% object for GP

   if (self%initialised) then
      if (allocated(self%f_r)) deallocate(self%f_r)
      if (allocated(self%g_r)) deallocate(self%g_r)
      if (allocated(self%Cmat_inv_v)) deallocate(self%Cmat_inv_v)
      if (allocated(self%k)) deallocate(self%k)
      if (allocated(self%kp)) deallocate(self%kp)
      if (allocated(self%mat_inv_k)) deallocate(self%mat_inv_k)
      call finalise(self%Cmat)
      call finalise(self%noise_Kmm)
      call finalise(self%Kmm)
   endif
   self%len_scale_sq = 0.0_dp
   self%f_var = 0.0_dp
   self%m_f = 0
   self%m_g = 0
   self%m_teach = 0
   self%sparsified = .false.
   self%initialised = .false.
end subroutine gp_basic_finalise

function f_predict_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r !% position at which to f_predict value
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1, x2(:), f_var, len_scale_sq
	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_r

   if (.not. self%initialised) then
      f_predict_r = 0.0_dp
      return
   endif

   if (self%m_f > 0) call kernel(r, self%f_r(1:self%m_f), self%f_var, self%len_scale_sq, f_f=self%k(1:self%m_f))
   if (self%m_g > 0) call kernel(r, self%g_r(1:self%m_g), self%f_var, self%len_scale_sq, f_g=self%k(self%m_f+1:self%m_f+self%m_g))
   f_predict_r = dot_product(self%k, self%Cmat_inv_v)
end function f_predict_r

function f_predict_grad_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r !% position at which to f_predict value
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1, x2(:), f_var, len_scale_sq
	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_grad_r

   if (.not. self%initialised) then
      f_predict_grad_r = 0.0_dp
      return
   endif

   if (self%m_f > 0) call kernel(r, self%f_r(1:self%m_f), self%f_var, self%len_scale_sq, g_f=self%k(1:self%m_f))
   if (self%m_g > 0) call kernel(r, self%g_r(1:self%m_g), self%f_var, self%len_scale_sq, g_g=self%k(self%m_f+1:self%m_f+self%m_g))
   f_predict_grad_r = dot_product(self%k, self%Cmat_inv_v)
end function f_predict_grad_r

function f_predict_var_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1, x2(:), f_var, len_scale_sq
	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_var_r

   if (.not. self%initialised) then
      f_predict_var_r = 0.0_dp
      return
   endif

   if (self%m_f > 0) call kernel(r, self%f_r(1:self%m_f), self%f_var, self%len_scale_sq, f_f=self%k(1:self%m_f))
   if (self%m_g > 0) call kernel(r, self%g_r(1:self%m_g), self%f_var, self%len_scale_sq, f_g=self%k(self%m_f+1:self%m_f+self%m_g))

   f_predict_var_r = self%f_var
   if (self%sparsified) then
      call Matrix_QR_Solve(self%Kmm, self%k, self%mat_inv_k)
      f_predict_var_r = f_predict_var_r - dot_product(self%k, self%mat_inv_k)
      call Matrix_QR_Solve(self%Cmat, self%k, self%mat_inv_k)
      f_predict_var_r = f_predict_var_r + dot_product(self%k, self%mat_inv_k)
   else
      call Matrix_QR_Solve(self%noise_Kmm, self%k, self%mat_inv_k)
      f_predict_var_r = f_predict_var_r - dot_product(self%k, self%mat_inv_k)
   endif
end function f_predict_var_r

function f_predict_grad_var_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1, x2(:), f_var, len_scale_sq
	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_grad_var_r

   if (.not. self%initialised) then
      f_predict_grad_var_r = 0.0_dp
      return
   endif

   ! From Lockwood and Anitescu preprint ANL/MCS-P1808-1110 "Gradient-Enhanced Universal Kriging for Uncertainty Propagation"
   ! Eq. 2.22, not including last term which is relevant only for underlying polynomial basis (which we don't have)
   if (self%m_f > 0) call kernel(r, self%f_r(1:self%m_f), self%f_var, self%len_scale_sq, g_f=self%k(1:self%m_f))
   if (self%m_g > 0) call kernel(r, self%g_r(1:self%m_g), self%f_var, self%len_scale_sq, g_g=self%k(self%m_f+1:self%m_f+self%m_g))

   f_predict_grad_var_r = self%f_var/self%len_scale_sq 
   if (self%sparsified) then
      ! -k' K_{mm}^{-1} k'
      call Matrix_QR_Solve(self%Kmm, self%k, self%mat_inv_k)
      ! first term should be second derivative of covariance, but it's hard wired for now
      f_predict_grad_var_r = f_predict_grad_var_r - dot_product(self%k, self%mat_inv_k)
      ! -k' C^{-1} k'
      call Matrix_QR_Solve(self%Cmat, self%k, self%mat_inv_k)
      ! first term should be second derivative of covariance, but it's hard wired for now
      f_predict_grad_var_r = f_predict_grad_var_r + dot_product(self%k, self%mat_inv_k)
   else
      ! -k' C^{-1} k'
      call Matrix_QR_Solve(self%noise_Kmm, self%k, self%mat_inv_k)
      ! first term should be second derivative of covariance, but it's hard wired for now
      f_predict_grad_var_r = f_predict_grad_var_r - dot_product(self%k, self%mat_inv_k)
   endif
end function f_predict_grad_var_r

function f_predict_var_grad_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1, x2(:), f_var, len_scale_sq
	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_var_grad_r

   if (.not. self%initialised) then
      f_predict_var_grad_r = 0.0_dp
      return
   endif

   if (self%m_f > 0) call kernel(r, self%f_r(1:self%m_f), self%f_var, self%len_scale_sq, f_f=self%k(1:self%m_f))
   if (self%m_g > 0) call kernel(r, self%g_r(1:self%m_g), self%f_var, self%len_scale_sq, f_g=self%k(self%m_f+1:self%m_f+self%m_g))
   if (self%m_f > 0) call kernel(r, self%f_r(1:self%m_f), self%f_var, self%len_scale_sq, g_f=self%kp(1:self%m_f))
   if (self%m_g > 0) call kernel(r, self%g_r(1:self%m_g), self%f_var, self%len_scale_sq, g_g=self%kp(self%m_f+1:self%m_f+self%m_g))

   f_predict_var_grad_r = 0.0_dp
   if (self%sparsified) then
      call Matrix_QR_Solve(self%Kmm, self%k, self%mat_inv_k)
      f_predict_var_grad_r = f_predict_var_grad_r - dot_product(self%kp, self%mat_inv_k)
      call Matrix_QR_Solve(self%Kmm, self%kp, self%mat_inv_k)
      f_predict_var_grad_r = f_predict_var_grad_r - dot_product(self%k, self%mat_inv_k)
      call Matrix_QR_Solve(self%Cmat, self%k, self%mat_inv_k)
      f_predict_var_grad_r = f_predict_var_grad_r + dot_product(self%kp, self%mat_inv_k)
      call Matrix_QR_Solve(self%Cmat, self%kp, self%mat_inv_k)
      f_predict_var_grad_r = f_predict_var_grad_r + dot_product(self%k, self%mat_inv_k)
   else
      call Matrix_QR_Solve(self%noise_Kmm, self%k, self%mat_inv_k)
      f_predict_var_grad_r = f_predict_var_grad_r - dot_product(self%kp, self%mat_inv_k)
      call Matrix_QR_Solve(self%noise_Kmm, self%kp, self%mat_inv_k)
      f_predict_var_grad_r = f_predict_var_grad_r - dot_product(self%k, self%mat_inv_k)
   endif
end function f_predict_var_grad_r

subroutine SE_kernel_r_rr(x1, x2, f_var, len_scale_sq, f_f, g_f, f_g, g_g)
   real(dp), intent(in) :: x1, x2(:), f_var, len_scale_sq
   real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)

   if (present(f_f)) f_f(:) = f_var*exp(-0.5_dp*((x2(:)-x1)**2)/len_scale_sq)
   if (present(g_f)) g_f(:) = (x2(:)-x1)/len_scale_sq * f_var*exp(-0.5_dp*((x2(:)-x1)**2)/len_scale_sq)
   if (present(f_g)) f_g(:) = -(x2(:)-x1)/len_scale_sq * f_var*exp(-0.5_dp*((x2(:)-x1)**2)/len_scale_sq)
   if (present(g_g)) g_g(:) = (1.0_dp-(x1-x2(:))**2/len_scale_sq)/len_scale_sq * f_var*exp(-0.5_dp*((x2(:)-x1)**2)/len_scale_sq)
end subroutine SE_kernel_r_rr

end module gp_basic_module
