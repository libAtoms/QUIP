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

! public :: gp_basic, initialise, finalise, f_predict, f_predict_var, f_predict_grad
public :: SE_kernel_r_rr

type gp_basic
   real(dp), allocatable :: len_scale_sq(:)
   logical, allocatable :: periodic(:)
   real(dp) :: f_var
   integer  :: n_dof, m_f, m_g, m_teach
   type(LA_Matrix) :: Cmat, Kmm, noise_Kmm
   real(dp), allocatable :: f_r(:,:), g_r(:,:), Cmat_inv_v(:), k(:),  mat_inv_k(:), kp(:)
   logical :: sparsified = .false.
   logical :: initialised = .false.
end type gp_basic

!% initialise (and teach) a gp_basic
interface initialise
   module procedure gp_basic_initialise_nd
end interface initialise

!% finalise and deallocate a gp_basic
interface finalise
   module procedure gp_basic_finalise
end interface finalise

!% predict a function value from a gp
interface f_predict
   module procedure f_predict_r
end interface f_predict

!!% predict a gradient of a function value from a gp
!interface f_predict_grad
!   module procedure f_predict_grad_r
!end interface f_predict_grad

!% predict a function variance from a gp
interface f_predict_var
   module procedure f_predict_var_r
end interface f_predict_var

!!% predict the variance of a gradient of a function from a gp
!interface f_predict_grad_var
!   module procedure f_predict_grad_var_r
!end interface f_predict_grad_var
!
!!% predict the gradient of a function variance from a gp
!interface f_predict_var_grad
!   module procedure f_predict_var_grad_r
!end interface f_predict_var_grad

contains

subroutine gp_basic_initialise_nd(self, len_scale, periodic, f_var, kernel, f_r, f_v, f_n, g_r, g_v, g_n, f_sparse_set, g_sparse_set, jitter, error)
   type(gp_basic), intent(inout) :: self !% object to store GP
   real(dp), intent(in) :: len_scale(:)
   logical, intent(in) :: periodic(:)
   real(dp), intent(in) :: f_var !% length scale and variance prior for GP
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodic, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 logical, intent(in) :: periodic(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp), optional, intent(in) :: f_r(:,:), f_v(:), f_n(:) !% arrays of function positions, values, noise 
   real(dp), optional, intent(in) :: g_r(:,:), g_v(:,:), g_n(:,:) !% arrays of function gradient positions, values, noise 
   integer, optional, target, intent(in) :: f_sparse_set(:), g_sparse_set(:) !% sets of points to use for sparsifcation for values and gradients
   real(dp), intent(in), optional :: jitter
   integer, optional, intent(out) :: error !% error status

   integer, pointer :: u_f_sparse_set(:), u_g_sparse_set(:)
   real(dp), allocatable :: Cmat(:,:), Kmn(:,:), y(:), siginvsq_y(:), Kmn_siginvsq_y(:), siginvsq_Knm(:,:), Kmm(:,:)
   integer :: i, n_f, n_g, n_tot, i_glob, ii

   INIT_ERROR(error)

   call finalise(self)

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

   ! compare f_r to f_[gn]
   if (present(f_r)) then
      if (size(f_r,2) /= size(f_v) .or. size(f_r,2) /= size(f_n)) then
	 RAISE_ERROR("size(f_r,2)="//size(f_r,2)//" size(f_v)="//size(f_v)//" size(f_n)="//size(f_n)//" not all equal", error)
      endif
   endif
   ! compare g_r to g_[gn]
   if (present(g_r)) then
      if (any(shape(g_r) /= shape(g_v)) .or. any(shape(g_r) /= shape(g_n))) then
	 RAISE_ERROR("shape(g_r)="//shape(g_r)//" shape(g_v)="//shape(g_v)//" shape(g_n)="//shape(g_n)//" not all equal", error)
      endif
   endif
   ! compare f_r to g_r
   if (present(f_r) .and. present(g_r)) then
      if (size(f_r,1) /= size(g_r,1)) then
	 RAISE_ERROR("size(f_r,1)="//size(f_r,1)//" size(g_r,1)="//size(g_r,1)//" not equal", error)
      endif
   endif

   self%sparsified = .false.
   if (present(f_r)) then
      self%n_dof = size(f_r,1)
      n_f = size(f_r,2)
      ! check and if needed point u_[fg]_sparse_set
      if (present(f_sparse_set)) then
	 if (any(f_sparse_set < 1) .or. any(f_sparse_set > n_f)) then
	    RAISE_ERROR("some of f_sparse_set out of range "//minval(f_sparse_set)//" "//maxval(f_sparse_set)//" "//n_f, error)
	 endif
	 u_f_sparse_set => f_sparse_set
	 self%sparsified = .true.
      else
	 allocate(u_f_sparse_set(n_f))
	 do i=1, n_f
	    u_f_sparse_set(i) = i
	 end do
      endif
      self%m_f = size(u_f_sparse_set)
      allocate(self%f_r(self%n_dof,self%m_f))
      self%f_r(:,:) = f_r(:,u_f_sparse_set)
   else
      allocate(self%f_r(1,1))
      self%m_f = 0
      n_f = 0
   endif
   if (present(g_r)) then
      self%n_dof = size(g_r,1)
      n_g = size(g_r,2)
      ! check and if needed point u_[fg]_sparse_set
      if (present(g_sparse_set)) then
	 if (any(g_sparse_set < 1) .or. any(g_sparse_set > n_g)) then
	    RAISE_ERROR("some of g_sparse_set out of range "//minval(g_sparse_set)//" "//maxval(g_sparse_set)//" "//n_g, error)
	 endif
	 u_g_sparse_set => g_sparse_set
	 self%sparsified = .true.
      else
	 allocate(u_g_sparse_set(n_g))
	 do i=1, n_g
	    u_g_sparse_set(i) = i
	 end do
      endif
      self%m_g = size(u_g_sparse_set)
      allocate(self%g_r(self%n_dof, self%m_g))
      self%g_r(:,:) = g_r(:,u_g_sparse_set)
   else
      allocate(self%g_r(1,1))
      self%m_g = 0
      n_g = 0
   endif

   if (size(len_scale) /= self%n_dof) then
      RAISE_ERROR("size(len_scale)="//size(len_scale)//" /= n_dof="//self%n_dof, error)
   endif
   if (size(periodic) /= self%n_dof) then
      RAISE_ERROR("size(periodic)="//size(periodic)//" /= n_dof="//self%n_dof, error)
   endif

   allocate(self%len_scale_sq(self%n_dof))
   allocate(self%periodic(self%n_dof))
   self%len_scale_sq = len_scale**2
   self%periodic = periodic
   self%f_var = f_var

   if (n_f + n_g == 0) then
      RAISE_ERROR("no teaching data provided", error)
   endif
   if (self%m_f + self%m_g == 0) then
      RAISE_ERROR("no sparsified teaching data provided", error)
   endif

   self%m_teach = self%m_f + self%m_g*self%n_dof
   n_tot = n_f + n_g*self%n_dof

   allocate(Kmm(self%m_teach, self%m_teach))

   ! everyone needs Kmm
   call kernel_mat(self%m_f, self%f_r, self%m_g, self%g_r, &
                   self%m_f, self%f_r, self%m_g, self%g_r, &
                   self%f_var, self%len_scale_sq, self%periodic, kernel, Kmm)

   if (self%sparsified) then
      allocate(Kmn(self%m_teach, n_tot))
      allocate(Kmn_siginvsq_y(self%m_teach))
      allocate(siginvsq_Knm(n_tot, self%m_teach))
      allocate(siginvsq_y(n_tot))
      allocate(Cmat(self%m_teach, self%m_teach))

      call kernel_mat(self%m_f, self%f_r, self%m_g, self%g_r, & 
		      n_f, f_r, n_g, g_r, &
		      self%f_var, self%len_scale_sq, self%periodic, kernel, Kmn)

      ! we'll need Kmm^{-1} for variance
      call initialise(self%Kmm, Kmm)

      if (n_f > 0) siginvsq_y(1:n_f) = f_v(:)/f_n(:)
      if (n_g > 0) then
	 do i=1, n_g
	    siginvsq_y(n_f+(i-1)*self%n_dof+1:n_f+i*self%n_dof) = g_v(:,i)/g_n(:,i)
	 end do
      endif
      call matrix_product_sub(Kmn_siginvsq_y, Kmn, siginvsq_y)

      siginvsq_Knm = transpose(Kmn)
      do i=1, n_f
	 siginvsq_Knm(i,:) = siginvsq_Knm(i,:) / f_n(i)
      end do
      do i=1, n_g
	 do ii=1, self%n_dof
	    siginvsq_Knm(n_f+(i-1)*self%n_dof+ii,:) = siginvsq_Knm(n_f+(i-1)*self%n_dof+ii,:) / g_n(ii,i)
	 end do
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
      ! not sparsified
      allocate(y(n_tot))

      ! we'll need self%noise_Kmm, which is Kmm shifted by noise, for variance
      do i=1, self%m_f
	 Kmm(i,i) = Kmm(i,i) + f_n(i)
      end do
      do i=1, self%m_g
	 do ii=1, self%n_dof
	    i_glob = self%m_f + (i-1)*self%n_dof+ii
	    Kmm(i_glob, i_glob)= Kmm(i_glob, i_glob) + g_n(ii,i)
	 end do
      end do
      call initialise(self%noise_Kmm, Kmm)

      if (n_f > 0) y(1:n_f) = f_v(1:n_f)
      do i=1, n_g
	 i_glob = n_f + (i-1)*self%n_dof
	 y(i_glob+1:i_glob+self%n_dof) = g_v(1:self%n_dof,i)
      end do

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
end subroutine gp_basic_initialise_nd

subroutine gp_basic_finalise(self)
   type(gp_basic), intent(inout) :: self !% object for GP

   if (self%initialised) then
      if (allocated(self%f_r)) deallocate(self%f_r)
      if (allocated(self%g_r)) deallocate(self%g_r)
      if (allocated(self%Cmat_inv_v)) deallocate(self%Cmat_inv_v)
      if (allocated(self%k)) deallocate(self%k)
      if (allocated(self%kp)) deallocate(self%kp)
      if (allocated(self%mat_inv_k)) deallocate(self%mat_inv_k)
      if (allocated(self%len_scale_sq)) deallocate(self%len_scale_sq)
      if (allocated(self%periodic)) deallocate(self%periodic)
      call finalise(self%Cmat)
      call finalise(self%noise_Kmm)
      call finalise(self%Kmm)
   endif
   self%f_var = 0.0_dp
   self%m_f = 0
   self%m_g = 0
   self%m_teach = 0
   self%n_dof = 0
   self%sparsified = .false.
   self%initialised = .false.
end subroutine gp_basic_finalise


function f_predict_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r(:) !% position at which to f_predict value
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodic, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 logical, intent(in) :: periodic(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_r

   if (.not. self%initialised) then
      f_predict_r = 0.0_dp
      return
   endif

   call f_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, self%periodic, kernel, self%k)
   f_predict_r = dot_product(self%k, self%Cmat_inv_v)
end function f_predict_r

!function f_predict_grad_r(self, r, kernel)
!   type(gp_basic), intent(inout) :: self !% object for GP
!   real(dp), intent(in) :: r !% position at which to f_predict value
!   interface
!      subroutine kernel(x1, x2, f_var, len_scale_sq, f_f, g_f, f_g, g_g)
!	 use system_module, only : dp
!	 real(dp), intent(in) :: x1, x2(:), f_var, len_scale_sq
!	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
!      end subroutine kernel
!   end interface
!   real(dp) :: f_predict_grad_r
!
!   if (.not. self%initialised) then
!      f_predict_grad_r = 0.0_dp
!      return
!   endif
!
!   call g_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, kernel, self%k)
!   f_predict_grad_r = dot_product(self%k, self%Cmat_inv_v)
!end function f_predict_grad_r

function f_predict_var_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodic, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 logical, intent(in) :: periodic(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_var_r

   if (.not. self%initialised) then
      f_predict_var_r = 0.0_dp
      return
   endif

   call f_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, self%periodic, kernel, self%k)

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

!function f_predict_grad_var_r(self, r, kernel)
!   type(gp_basic), intent(inout) :: self
!   real(dp), intent(in) :: r
!   interface
!      subroutine kernel(x1, x2, f_var, len_scale_sq, f_f, g_f, f_g, g_g)
!	 use system_module, only : dp
!	 real(dp), intent(in) :: x1, x2(:), f_var, len_scale_sq
!	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
!      end subroutine kernel
!   end interface
!   real(dp) :: f_predict_grad_var_r
!
!   if (.not. self%initialised) then
!      f_predict_grad_var_r = 0.0_dp
!      return
!   endif
!
!   ! From Lockwood and Anitescu preprint ANL/MCS-P1808-1110 "Gradient-Enhanced Universal Kriging for Uncertainty Propagation"
!   ! Eq. 2.22, not including last term which is relevant only for underlying polynomial basis (which we don't have)
!   call g_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, kernel, self%k)
!
!   f_predict_grad_var_r = self%f_var/self%len_scale_sq 
!   if (self%sparsified) then
!      ! -k' K_{mm}^{-1} k'
!      call Matrix_QR_Solve(self%Kmm, self%k, self%mat_inv_k)
!      ! first term should be second derivative of covariance, but it's hard wired for now
!      f_predict_grad_var_r = f_predict_grad_var_r - dot_product(self%k, self%mat_inv_k)
!      ! -k' C^{-1} k'
!      call Matrix_QR_Solve(self%Cmat, self%k, self%mat_inv_k)
!      ! first term should be second derivative of covariance, but it's hard wired for now
!      f_predict_grad_var_r = f_predict_grad_var_r + dot_product(self%k, self%mat_inv_k)
!   else
!      ! -k' C^{-1} k'
!      call Matrix_QR_Solve(self%noise_Kmm, self%k, self%mat_inv_k)
!      ! first term should be second derivative of covariance, but it's hard wired for now
!      f_predict_grad_var_r = f_predict_grad_var_r - dot_product(self%k, self%mat_inv_k)
!   endif
!end function f_predict_grad_var_r
!
!function f_predict_var_grad_r(self, r, kernel)
!   type(gp_basic), intent(inout) :: self
!   real(dp), intent(in) :: r
!   interface
!      subroutine kernel(x1, x2, f_var, len_scale_sq, f_f, g_f, f_g, g_g)
!	 use system_module, only : dp
!	 real(dp), intent(in) :: x1, x2(:), f_var, len_scale_sq
!	 real(dp), optional :: f_f(:), g_f(:), f_g(:), g_g(:)
!      end subroutine kernel
!   end interface
!   real(dp) :: f_predict_var_grad_r
!
!   if (.not. self%initialised) then
!      f_predict_var_grad_r = 0.0_dp
!      return
!   endif
!
!   call f_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, kernel, self%k)
!   call g_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, kernel, self%kp)
!
!   f_predict_var_grad_r = 0.0_dp
!   if (self%sparsified) then
!      call Matrix_QR_Solve(self%Kmm, self%k, self%mat_inv_k)
!      f_predict_var_grad_r = f_predict_var_grad_r - dot_product(self%kp, self%mat_inv_k)
!      call Matrix_QR_Solve(self%Kmm, self%kp, self%mat_inv_k)
!      f_predict_var_grad_r = f_predict_var_grad_r - dot_product(self%k, self%mat_inv_k)
!      call Matrix_QR_Solve(self%Cmat, self%k, self%mat_inv_k)
!      f_predict_var_grad_r = f_predict_var_grad_r + dot_product(self%kp, self%mat_inv_k)
!      call Matrix_QR_Solve(self%Cmat, self%kp, self%mat_inv_k)
!      f_predict_var_grad_r = f_predict_var_grad_r + dot_product(self%k, self%mat_inv_k)
!   else
!      call Matrix_QR_Solve(self%noise_Kmm, self%k, self%mat_inv_k)
!      f_predict_var_grad_r = f_predict_var_grad_r - dot_product(self%kp, self%mat_inv_k)
!      call Matrix_QR_Solve(self%noise_Kmm, self%kp, self%mat_inv_k)
!      f_predict_var_grad_r = f_predict_var_grad_r - dot_product(self%k, self%mat_inv_k)
!   endif
!end function f_predict_var_grad_r

subroutine kernel_mat(l_n_f, l_f_r, l_n_g, l_g_r, r_n_f, r_f_r, r_n_g, r_g_r, f_var, len_scale_sq, periodic, kernel, mat)
   integer, intent(in) :: l_n_f, l_n_g
   real(dp), optional, intent(in) :: l_f_r(:,:), l_g_r(:,:)
   integer, intent(in) :: r_n_f, r_n_g
   real(dp), optional, intent(in) :: r_f_r(:,:), r_g_r(:,:)
   real(dp), intent(in) :: f_var, len_scale_sq(:)
   logical, intent(in) :: periodic(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodic, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 logical, intent(in) :: periodic(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp), intent(out) :: mat(:,:)

   integer :: i, i_glob, n_dof

   if (present(l_f_r)) n_dof=size(l_f_r,1)
   if (present(l_g_r)) n_dof=size(l_g_r,1)

   do i=1, l_n_f
      ! f on f
      if (r_n_f > 0) call kernel(l_f_r(1:n_dof,i), r_f_r(1:n_dof,:), f_var, len_scale_sq(1:n_dof), periodic(1:n_dof), f_f=mat(i,1:r_n_f))
      ! f on g
      if (r_n_g > 0) call kernel(l_f_r(1:n_dof,i), r_g_r(1:n_dof,:), f_var, len_scale_sq(1:n_dof), periodic(1:n_dof), f_g=mat(i,r_n_f+1:r_n_f+r_n_g*n_dof))
   end do
   do i=1, l_n_g
      i_glob = l_n_f + (i-1)*n_dof + 1
      ! g on f
      if (r_n_f > 0) call kernel(l_g_r(1:n_dof,i), r_f_r(1:n_dof,:), f_var, len_scale_sq(1:n_dof), periodic(1:n_dof), g_f=mat(i_glob:i_glob+n_dof-1,1:r_n_f))
      ! g on g
      if (r_n_g > 0) call kernel(l_g_r(1:n_dof,i), r_g_r(1:n_dof,:), f_var, len_scale_sq(1:n_dof), periodic(1:n_dof), g_g=mat(i_glob:i_glob+n_dof-1,r_n_f+1:r_n_f+r_n_g*n_dof))
   end do
end subroutine

subroutine f_kernel_vec(r, n_f, f_r, n_g, g_r, f_var, len_scale_sq, periodic, kernel, vec)
   real(dp), intent(in) :: r(:)
   integer, intent(in) :: n_f
   real(dp), intent(in) :: f_r(:,:)
   integer, intent(in) :: n_g
   real(dp), intent(in) :: g_r(:,:)
   real(dp), intent(in) :: f_var, len_scale_sq(:)
   logical, intent(in) :: periodic(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodic, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 logical, intent(in) :: periodic(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp), intent(out) :: vec(:)

   integer :: n_dof

   n_dof=size(r)

   if (n_f > 0) call kernel(r(1:n_dof), f_r(1:n_dof,1:n_f), f_var, len_scale_sq(1:n_dof), periodic(1:n_dof), f_f=vec(1:n_f))
   if (n_g > 0) call kernel(r(1:n_dof), g_r(1:n_dof,1:n_g), f_var, len_scale_sq(1:n_dof), periodic(1:n_dof), f_g=vec(n_f+1:n_f+n_g*n_dof))
end subroutine f_kernel_vec

subroutine g_kernel_vec(r, n_f, f_r, n_g, g_r, f_var, len_scale_sq, periodic, kernel, vec)
   real(dp), intent(in) :: r(:)
   integer, intent(in) :: n_f
   real(dp), intent(in) :: f_r(:,:)
   integer, intent(in) :: n_g
   real(dp), intent(in) :: g_r(:,:)
   real(dp), intent(in) :: f_var, len_scale_sq(:)
   logical, intent(in) :: periodic(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodic, f_f, g_f, f_g, g_g)
	 use system_module, only : dp
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 logical, intent(in) :: periodic(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp), intent(out) :: vec(:,:)

   integer n_dof

   n_dof = size(r)

   if (n_f > 0) call kernel(r(1:n_dof), f_r(1:n_dof,1:n_f), f_var, len_scale_sq(1:n_dof), periodic(1:n_dof), g_f=vec(1:n_dof,1:n_f))
   if (n_g > 0) call kernel(r(1:n_dof), g_r(1:n_dof,1:n_g), f_var, len_scale_sq(1:n_dof), periodic(1:n_dof), g_g=vec(1:n_dof,n_f+1:n_f+n_g*n_dof))
end subroutine g_kernel_vec

subroutine SE_kernel_r_rr(x1, x2, f_var, len_scale_sq, periodic, f_f, g_f, f_g, g_g)
   real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
   logical, intent(in) :: periodic(:)
   real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)

   real(dp), allocatable :: exp_arg(:), dexp_arg_i(:,:), ddexp_arg_i(:,:), dexp_arg_j(:,:)
   integer :: i, j, nv, n_dof


   n_dof = size(x1)
   nv = size(x2,2)

if (present(f_f)) then
if (size(f_f) /= nv) call system_abort("size(f_f) bad "//shape(f_f)//" "//nv//" "//n_dof)
endif
if (present(g_f)) then
if (size(g_f,1) /= n_dof .or. size(g_f,2) /= nv) call system_abort("size(g_f) bad "//shape(g_f)//" "//nv//" "//n_dof)
endif
if (present(f_g)) then
if (size(f_g,1) /= n_dof*nv) call system_abort("size(f_g) bad "//shape(f_g)//" "//nv//" "//n_dof)
endif
if (present(g_g)) then
if (size(g_g,1) /= n_dof .or. size(g_g,2) /= nv*n_dof) call system_abort("size(g_g) bad "//shape(f_g)//" "//nv//" "//n_dof)
endif

   allocate(exp_arg(nv))
   allocate(dexp_arg_i(n_dof,nv))
   allocate(ddexp_arg_i(n_dof,nv))
   allocate(dexp_arg_j(n_dof,nv))
   exp_arg = 0.0_dp
   do i=1, n_dof
      if (periodic(i)) then
	 exp_arg(:) = exp_arg(:) - 2.0_dp*sin((x2(i,:)-x1(i))/2.0_dp)**2/len_scale_sq(i)
      else
	 exp_arg(:) = exp_arg(:) - 0.5_dp*(x2(i,:)-x1(i))**2/len_scale_sq(i)
      endif
   end do

   if (present(f_f)) f_f(:) = f_var * exp(exp_arg(:))
   if (present(g_f) .or. present(g_g)) then
      do i=1, n_dof
	 if (periodic(i)) then
	    dexp_arg_i(i,:) = 2.0_dp*sin((x2(i,:)-x1(i))/2.0_dp)*cos((x2(i,:)-x1(i))/2.0_dp) / len_scale_sq(i)
	    if (present(g_g)) ddexp_arg_i(i,:) = (-sin((x2(i,:)-x1(i))/2.0_dp)**2+cos((x2(i,:)-x1(i))/2.0_dp)**2)/len_scale_sq(i)
	 else
	    dexp_arg_i(i,:) = (x2(i,:)-x1(i))/len_scale_sq(i)
	    if (present(g_g)) ddexp_arg_i(i,:) = 1.0_dp/len_scale_sq(i)
	 endif
      end do
   endif
   if (present(f_g) .or. present(g_g)) then
      do j=1, n_dof
	 if (periodic(j)) then
	    dexp_arg_j(j,:) = -2.0_dp*sin((x2(j,:)-x1(j))/2.0_dp)*cos((x2(j,:)-x1(j))/2.0_dp) / len_scale_sq(j)
	 else
	    dexp_arg_j(j,:) = -(x2(j,:)-x1(j))/len_scale_sq(j)
	 endif
      end do
   endif

   if (present(g_f)) then
      do i=1, n_dof
	 g_f(i,1:nv) = f_var * exp(exp_arg(:))*(dexp_arg_i(i,:))
      end do
   endif
   if (present(f_g)) then
      do j=1, n_dof
	 f_g(j:nv*n_dof:n_dof) = f_var * exp(exp_arg(:))*(dexp_arg_j(j,:))
      end do
   endif

   if (present(g_g)) then
      do i=1, n_dof
      do j=1, n_dof
	 if (i /= j) then
	    g_g(i,j:nv*n_dof:n_dof) = f_var * exp(exp_arg(:))*(dexp_arg_i(i,:))*(dexp_arg_j(j,:))
	 else
	    g_g(i,j:nv*n_dof:n_dof) = f_var * exp(exp_arg(:))*(ddexp_arg_i(i,:)-dexp_arg_i(i,:)**2)
	 endif
      end do
      end do
   end if

   deallocate(exp_arg)
   deallocate(dexp_arg_i)
   deallocate(ddexp_arg_i)
   deallocate(dexp_arg_j)

end subroutine SE_kernel_r_rr

end module gp_basic_module
