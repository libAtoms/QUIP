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

!% A simple Gaussian process module for n-D functions learned from samples of
!% their values or gradients, with sparsification

module gp_basic_module
implicit none
private

public :: gp_basic, initialise, finalise, f_predict, f_predict_var, f_predict_grad, f_predict_var_grad, f_predict_grad_var
public :: SE_kernel_r_rr

integer, parameter :: dp = selected_real_kind(15,307)

real(dp), parameter :: TWO_PI = 2.0_dp*3.14159265358979_dp

integer, parameter :: NOT_FACTORISED = 0
integer, parameter :: QR             = 2

type GP_Matrix
   real(dp), dimension(:,:), allocatable :: matrix, factor
   real(dp), dimension(:), allocatable :: s, tau
   integer :: n, m
   logical :: initialised = .false.
   logical :: equilibrated = .false.
   integer :: factorised = NOT_FACTORISED
end type GP_Matrix

type gp_basic
   real(dp), allocatable :: len_scale_sq(:)
   real(dp), allocatable :: periodicity(:)
   real(dp) :: f_var
   integer  :: n_dof, m_f, m_g, m_teach
   type(GP_Matrix) :: Cmat, Kmm, noise_Kmm
   real(dp), allocatable :: f_r(:,:), g_r(:,:), Cmat_inv_v(:), k(:), mat_inv_k(:)
   real(dp), allocatable :: k_grad(:,:), mat_inv_k_grad(:,:)
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

!% predict a gradient of a function value from a gp
interface f_predict_grad
   module procedure f_predict_grad_r
end interface f_predict_grad

!% predict a function variance from a gp
interface f_predict_var
   module procedure f_predict_var_r, f_predict_var_rr
end interface f_predict_var

!% predict the variance of a gradient of a function from a gp
interface f_predict_grad_var
   module procedure f_predict_grad_var_r
end interface f_predict_grad_var

!% predict the gradient of a function variance from a gp
interface f_predict_var_grad
   module procedure f_predict_var_grad_r
end interface f_predict_var_grad

interface Matrix_QR_Solve
   module procedure GP_Matrix_QR_Solve_Vector, GP_Matrix_QR_Solve_Matrix
endinterface Matrix_QR_Solve

contains

subroutine gp_basic_initialise_nd(self, len_scale, periodicity, f_var, kernel, f_r, f_v, f_n, g_r, g_v, g_n, f_sparse_set, &
                                  g_sparse_set, jitter)
   type(gp_basic), intent(inout) :: self !% object to store GP
   real(dp), intent(in) :: len_scale(:)
   real(dp), intent(in) :: periodicity(:)
   real(dp), intent(in) :: f_var !% length scale and variance prior for GP
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
	 integer, parameter :: dp = 8
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 real(dp), intent(in) :: periodicity(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp), optional, intent(in) :: f_r(:,:), f_v(:), f_n(:) !% arrays of function positions, values, noise 
   real(dp), optional, intent(in) :: g_r(:,:), g_v(:,:), g_n(:,:) !% arrays of function gradient positions, values, noise 
   integer, optional, target, intent(in) :: f_sparse_set(:), g_sparse_set(:) !% sets of points to use for sparsifcation for values and gradients
   real(dp), intent(in), optional :: jitter

   integer, pointer :: u_f_sparse_set(:), u_g_sparse_set(:)
   real(dp), allocatable :: Cmat(:,:), Kmn(:,:), y(:), siginvsq_y(:), Kmn_siginvsq_y(:), siginvsq_Knm(:,:), Kmm(:,:)
   integer :: i, n_f, n_g, n_tot, i_glob, ii

   call finalise(self)

   if (f_var <= 0.0_dp) then
      print *,"invalid f_var=",f_var
      stop
   endif

   i = count( (/ present(f_r), present(f_v), present(f_n) /) )
   if (i /= 0 .and. i /= 3) then
      print *,"got something other than all or none of f_r, f_v, f_n"
      stop
   endif
   i = count( (/ present(g_r), present(g_v), present(g_n) /) )
   if (i /= 0 .and. i /= 3) then
      print *,"got something other than all or none of g_r, g_v, g_n"
      stop
   endif

   ! compare f_r to f_[gn]
   if (present(f_r)) then
      if (size(f_r,2) /= size(f_v) .or. size(f_r,2) /= size(f_n)) then
	 print *,"size(f_r,2)=",size(f_r,2)," size(f_v)=",size(f_v)," size(f_n)=",size(f_n)," not all equal"
	 stop
      endif
   endif
   ! compare g_r to g_[gn]
   if (present(g_r)) then
      if (any(shape(g_r) /= shape(g_v)) .or. any(shape(g_r) /= shape(g_n))) then
	 print *,"shape(g_r)=",shape(g_r)," shape(g_v)=",shape(g_v)," shape(g_n)=",shape(g_n)," not all equal"
	 stop
      endif
   endif
   ! compare f_r to g_r
   if (present(f_r) .and. present(g_r)) then
      if (size(f_r,1) /= size(g_r,1)) then
	 print *,"size(f_r,1)=",size(f_r,1)," size(g_r,1)=",size(g_r,1)," not equal"
	 stop
      endif
   endif

   self%sparsified = .false.
   if (present(f_r)) then
      self%n_dof = size(f_r,1)
      n_f = size(f_r,2)
      ! check and if needed point u_[fg]_sparse_set
      if (present(f_sparse_set)) then
	 if (any(f_sparse_set < 1) .or. any(f_sparse_set > n_f)) then
	    print *,"some of f_sparse_set out of range ",minval(f_sparse_set)," ",maxval(f_sparse_set)," ",n_f
	    stop
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
	    print *,"some of g_sparse_set out of range ",minval(g_sparse_set)," ",maxval(g_sparse_set)," ",n_g
	    stop
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
      print *,"size(len_scale)=",size(len_scale)," , n_dof=",self%n_dof
      stop
   endif
   if (size(periodicity) /= self%n_dof) then
      print *,"size(periodicity)=",size(periodicity)," /= n_dof=",self%n_dof
      stop
   endif

   allocate(self%len_scale_sq(self%n_dof))
   allocate(self%periodicity(self%n_dof))
   self%len_scale_sq = len_scale**2
   self%periodicity = periodicity
   self%f_var = f_var

   if (n_f + n_g == 0) then
      print *,"no teaching data provided"
      stop
   endif
   if (self%m_f + self%m_g == 0) then
      print *,"no sparsified teaching data provided"
      stop
   endif

   self%m_teach = self%m_f + self%m_g*self%n_dof
   n_tot = n_f + n_g*self%n_dof

   allocate(Kmm(self%m_teach, self%m_teach))

   ! everyone needs Kmm
   call kernel_mat(self%m_f, self%f_r, self%m_g, self%g_r, &
                   self%m_f, self%f_r, self%m_g, self%g_r, &
                   self%f_var, self%len_scale_sq, self%periodicity, kernel, Kmm)

   if (self%sparsified) then
      allocate(Kmn(self%m_teach, n_tot))
      allocate(Kmn_siginvsq_y(self%m_teach))
      allocate(siginvsq_Knm(n_tot, self%m_teach))
      allocate(siginvsq_y(n_tot))
      allocate(Cmat(self%m_teach, self%m_teach))

      call kernel_mat(self%m_f, self%f_r, self%m_g, self%g_r, & 
		      n_f, f_r, n_g, g_r, &
		      self%f_var, self%len_scale_sq, self%periodicity, kernel, Kmn)

      ! "jitter" (ABP e-mail 14 Aug)
      if (present(jitter)) then
	 do i=1, size(Kmm,1)
	    Kmm(i,i) = Kmm(i,i) + jitter
	 end do
      endif
      ! we'll need Kmm^{-1} for variance
      call gp_matrix_initialise(self%Kmm, Kmm)

      if (n_f > 0) siginvsq_y(1:n_f) = f_v(:)/f_n(:)
      if (n_g > 0) then
	 do i=1, n_g
	    siginvsq_y(n_f+(i-1)*self%n_dof+1:n_f+i*self%n_dof) = g_v(:,i)/g_n(:,i)
	 end do
      endif
      call dgemv('N', size(Kmn,1), size(Kmn,2), 1.0_dp, Kmn(1,1), size(Kmn,1), &
	 siginvsq_y(1), 1, 0.0_dp, Kmn_siginvsq_y(1), 1)

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
      call dgemm('N', 'N', size(Cmat,1), size(Cmat,2), size(Kmn,2), &
	 1.0_dp, Kmn(1,1), size(Kmn,1), siginvsq_Knm(1,1), size(siginvsq_Knm,1), &
	 1.0_dp, Cmat(1,1), size(Cmat,1))
      ! Cmat is now (K_{mm} + K_{mn} \Sigma^{-2} K_{nm})

      call gp_matrix_initialise(self%Cmat, Cmat)
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
      call gp_matrix_initialise(self%noise_Kmm, Kmm)

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
   allocate(self%k_grad(self%n_dof, self%m_teach))
   allocate(self%mat_inv_k(self%m_teach))
   allocate(self%mat_inv_k_grad(self%m_teach, self%n_dof))

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
      if (allocated(self%k_grad)) deallocate(self%k_grad)
      if (allocated(self%mat_inv_k)) deallocate(self%mat_inv_k)
      if (allocated(self%mat_inv_k_grad)) deallocate(self%mat_inv_k_grad)
      if (allocated(self%len_scale_sq)) deallocate(self%len_scale_sq)
      if (allocated(self%periodicity)) deallocate(self%periodicity)
      call gp_matrix_finalise(self%Cmat)
      call gp_matrix_finalise(self%noise_Kmm)
      call gp_matrix_finalise(self%Kmm)
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
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
	 integer, parameter :: dp = 8
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 real(dp), intent(in) :: periodicity(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_r

   real(dp), external :: ddot

   if (.not. self%initialised) then
      f_predict_r = 0.0_dp
      return
   endif

   call f_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, self%periodicity, kernel, self%k)
   f_predict_r = ddot(size(self%k), self%k(1), 1, self%Cmat_inv_v(1), 1)
end function f_predict_r

function f_predict_grad_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r(:) !% position at which to f_predict value
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
	 integer, parameter :: dp = 8
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 real(dp), intent(in) :: periodicity(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_grad_r(size(r))

   if (.not. self%initialised) then
      f_predict_grad_r = 0.0_dp
      return
   endif

   call g_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, self%periodicity, &
      kernel, self%k_grad)
   call dgemv('N', size(self%k_grad,1), size(self%k_grad,2), 1.0_dp, self%k_grad(1,1), size(self%k_grad,1), &
      self%Cmat_inv_v(1), 1, 0.0_dp, f_predict_grad_r(1), 1)
end function f_predict_grad_r

function f_predict_var_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
	 integer, parameter :: dp = 8
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 real(dp), intent(in) :: periodicity(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_var_r

   real(dp), external :: ddot

   if (.not. self%initialised) then
      f_predict_var_r = 0.0_dp
      return
   endif

   call f_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, self%periodicity, kernel, self%k)
   ! print *, "kernel var size ", size(self%k), " num above 1% ",count(self%k > 0.01_dp*maxval(self%k))

   f_predict_var_r = self%f_var
   if (self%sparsified) then
      ! from Wlliams and Rasumussen, Gaussian Processes for Machine Learning, MIT Press, 2006 Eq. 8.27
      call Matrix_QR_Solve(self%Kmm, self%k, self%mat_inv_k)
      f_predict_var_r = f_predict_var_r - ddot(size(self%k), self%k(1), 1, self%mat_inv_k(1), 1)
      call Matrix_QR_Solve(self%Cmat, self%k, self%mat_inv_k)
      f_predict_var_r = f_predict_var_r + ddot(size(self%k), self%k(1), 1, self%mat_inv_k(1), 1)
   else
      call Matrix_QR_Solve(self%noise_Kmm, self%k, self%mat_inv_k)
      f_predict_var_r = f_predict_var_r - ddot(size(self%k), self%k(1), 1, self%mat_inv_k(1), 1)
   endif
end function f_predict_var_r

function f_predict_var_rr(self, rr, kernel)
   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: rr(:,:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
	 integer, parameter :: dp = 8
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 real(dp), intent(in) :: periodicity(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_var_rr(size(rr,2))

   real(dp), external :: ddot
   integer :: i
   real(dp), allocatable :: kk(:,:), mat_inv_kk(:,:)

   if (.not. self%initialised) then
      f_predict_var_rr = 0.0_dp
      return
   endif

   allocate(kk(size(self%k, 1), size(rr, 2)), mat_inv_kk(size(self%k,1), size(rr,2)))
   do i=1, size(rr, 2)
      call f_kernel_vec(rr(:,i), self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, self%periodicity, kernel, kk(:,i))
      ! print *, "kernel var size ", size(self%k), " num above 1% ",count(self%k > 0.01_dp*maxval(self%k))
   end do

   f_predict_var_rr = self%f_var
   if (self%sparsified) then
      ! from Wlliams and Rasumussen, Gaussian Processes for Machine Learning, MIT Press, 2006 Eq. 8.27
      call Matrix_QR_Solve(self%Kmm, kk, mat_inv_kk)
      do i=1, size(rr, 2)
	 f_predict_var_rr(i) = f_predict_var_rr(i) - ddot(size(kk,1), kk(1,i), 1, mat_inv_kk(1,i), 1)
      end do
      call Matrix_QR_Solve(self%Cmat, kk, mat_inv_kk)
      do i=1, size(rr, 2)
	 f_predict_var_rr(i) = f_predict_var_rr(i) + ddot(size(kk,1), kk(1,i), 1, mat_inv_kk(1,i), 1)
      end do
   else
      call Matrix_QR_Solve(self%noise_Kmm, kk, mat_inv_kk)
      do i=1, size(rr, 2)
	 f_predict_var_rr(i) = f_predict_var_rr(i) - ddot(size(kk,1), kk(1,i), 1, mat_inv_kk(1,i), 1)
      end do
   endif
   deallocate(kk, mat_inv_kk)
end function f_predict_var_rr

function f_predict_grad_var_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
	 integer, parameter :: dp = 8
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 real(dp), intent(in) :: periodicity(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_grad_var_r(size(r))

   integer :: ii

   real(dp), external :: ddot

   if (.not. self%initialised) then
      f_predict_grad_var_r = 0.0_dp
      return
   endif

   ! From Lockwood and Anitescu preprint ANL/MCS-P1808-1110 "Gradient-Enhanced Universal Kriging for Uncertainty Propagation"
   ! Eq. 2.22, not including last term which is relevant only for underlying polynomial basis (which we don't have)
   call g_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, self%periodicity, &
      kernel, self%k_grad)

   ! first term should be second derivative of covariance, but it's hard wired for now
   f_predict_grad_var_r = self%f_var/self%len_scale_sq 
   if (self%sparsified) then
      ! -k' K_{mm}^{-1} k'
      call Matrix_QR_Solve(self%Kmm, transpose(self%k_grad), self%mat_inv_k_grad)
      ! first term should be second derivative of covariance, but it's hard wired for now
      do ii=1, self%n_dof
	 f_predict_grad_var_r(ii) = f_predict_grad_var_r(ii) - ddot(size(self%k_grad,2), self%k_grad(ii,1), size(self%k_grad,1), &
											 self%mat_inv_k_grad(1,ii), 1)
      end do
      ! -k' C^{-1} k'
      call Matrix_QR_Solve(self%Cmat, transpose(self%k_grad), self%mat_inv_k_grad)
      do ii=1, self%n_dof
	 f_predict_grad_var_r(ii) = f_predict_grad_var_r(ii) + ddot(size(self%k_grad,2), self%k_grad(ii,1), size(self%k_grad,1), &
											 self%mat_inv_k_grad(1,ii), 1)
      end do
   else
      ! -k' C^{-1} k'
      call Matrix_QR_Solve(self%noise_Kmm, transpose(self%k_grad), self%mat_inv_k_grad)
      do ii=1, self%n_dof
	 f_predict_grad_var_r(ii) = f_predict_grad_var_r(ii) - ddot(size(self%k_grad,2), self%k_grad(ii,1), size(self%k_grad,1), &
										         self%mat_inv_k_grad(1,ii), 1)
      end do
   endif
end function f_predict_grad_var_r

function f_predict_var_grad_r(self, r, kernel)
   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
	 integer, parameter :: dp = 8
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 real(dp), intent(in) :: periodicity(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_var_grad_r(size(r))

   if (.not. self%initialised) then
      f_predict_var_grad_r = 0.0_dp
      return
   endif

   call f_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, self%periodicity, &
      kernel, self%k)
   call g_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, self%periodicity, &
      kernel, self%k_grad)

   f_predict_var_grad_r = 0.0_dp
   if (self%sparsified) then
      call Matrix_QR_Solve(self%Kmm, self%k, self%mat_inv_k)
      call dgemv('N', size(self%k_grad,1), size(self%k_grad,2), -2.0_dp, self%k_grad(1,1), size(self%k_grad,1), &
	 self%mat_inv_k(1), 1, 1.0_dp, f_predict_var_grad_r(1), 1)
      call Matrix_QR_Solve(self%Cmat, self%k, self%mat_inv_k)
      call dgemv('N', size(self%k_grad,1), size(self%k_grad,2), 2.0_dp, self%k_grad(1,1), size(self%k_grad,1), &
	 self%mat_inv_k(1), 1, 1.0_dp, f_predict_var_grad_r(1), 1)
   else
      call Matrix_QR_Solve(self%noise_Kmm, self%k, self%mat_inv_k)
      call dgemv('N', size(self%k_grad,1), size(self%k_grad,2), -2.0_dp, self%k_grad(1,1), size(self%k_grad,1), &
	 self%mat_inv_k(1), 1, 1.0_dp, f_predict_var_grad_r(1), 1)
   endif
end function f_predict_var_grad_r

subroutine kernel_mat(l_n_f, l_f_r, l_n_g, l_g_r, r_n_f, r_f_r, r_n_g, r_g_r, f_var, len_scale_sq, periodicity, kernel, mat)
   integer, intent(in) :: l_n_f, l_n_g
   real(dp), optional, intent(in) :: l_f_r(:,:), l_g_r(:,:)
   integer, intent(in) :: r_n_f, r_n_g
   real(dp), optional, intent(in) :: r_f_r(:,:), r_g_r(:,:)
   real(dp), intent(in) :: f_var, len_scale_sq(:)
   real(dp), intent(in) :: periodicity(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
	 integer, parameter :: dp = 8
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 real(dp), intent(in) :: periodicity(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp), intent(out) :: mat(:,:)

   integer :: i, i_glob, n_dof

   if (l_n_f > 0) n_dof=size(l_f_r,1)
   if (l_n_g > 0) n_dof=size(l_g_r,1)

   do i=1, l_n_f
      ! f on f
      if (r_n_f > 0) call kernel(l_f_r(1:n_dof,i), r_f_r(1:n_dof,:), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
				 f_f=mat(i,1:r_n_f))
      ! f on g
      if (r_n_g > 0) call kernel(l_f_r(1:n_dof,i), r_g_r(1:n_dof,:), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
				 f_g=mat(i,r_n_f+1:r_n_f+r_n_g*n_dof))
   end do
   do i=1, l_n_g
      i_glob = l_n_f + (i-1)*n_dof + 1
      ! g on f
      if (r_n_f > 0) call kernel(l_g_r(1:n_dof,i), r_f_r(1:n_dof,:), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
				 g_f=mat(i_glob:i_glob+n_dof-1,1:r_n_f))
      ! g on g
      if (r_n_g > 0) call kernel(l_g_r(1:n_dof,i), r_g_r(1:n_dof,:), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
				 g_g=mat(i_glob:i_glob+n_dof-1,r_n_f+1:r_n_f+r_n_g*n_dof))
   end do
end subroutine kernel_mat

subroutine f_kernel_vec(r, n_f, f_r, n_g, g_r, f_var, len_scale_sq, periodicity, kernel, vec)
   real(dp), intent(in) :: r(:)
   integer, intent(in) :: n_f
   real(dp), intent(in) :: f_r(:,:)
   integer, intent(in) :: n_g
   real(dp), intent(in) :: g_r(:,:)
   real(dp), intent(in) :: f_var, len_scale_sq(:)
   real(dp), intent(in) :: periodicity(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
	 integer, parameter :: dp = 8
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 real(dp), intent(in) :: periodicity(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp), intent(out) :: vec(:)

   integer :: n_dof

   n_dof=size(r)

   if (n_f > 0) call kernel(r(1:n_dof), f_r(1:n_dof,1:n_f), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
			    f_f=vec(1:n_f))
   if (n_g > 0) call kernel(r(1:n_dof), g_r(1:n_dof,1:n_g), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
			    f_g=vec(n_f+1:n_f+n_g*n_dof))
end subroutine f_kernel_vec

subroutine g_kernel_vec(r, n_f, f_r, n_g, g_r, f_var, len_scale_sq, periodicity, kernel, vec)
   real(dp), intent(in) :: r(:)
   integer, intent(in) :: n_f
   real(dp), intent(in) :: f_r(:,:)
   integer, intent(in) :: n_g
   real(dp), intent(in) :: g_r(:,:)
   real(dp), intent(in) :: f_var, len_scale_sq(:)
   real(dp), intent(in) :: periodicity(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
	 integer, parameter :: dp = 8
	 real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
	 real(dp), intent(in) :: periodicity(:)
	 real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp), intent(out) :: vec(:,:)

   integer n_dof

   n_dof = size(r)

   if (n_f > 0) call kernel(r(1:n_dof), f_r(1:n_dof,1:n_f), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
			    g_f=vec(1:n_dof,1:n_f))
   if (n_g > 0) call kernel(r(1:n_dof), g_r(1:n_dof,1:n_g), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
			    g_g=vec(1:n_dof,n_f+1:n_f+n_g*n_dof))
end subroutine g_kernel_vec

subroutine SE_kernel_r_rr(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
   real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
   real(dp), intent(in) :: periodicity(:)
   real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)

   real(dp), allocatable :: exp_arg(:), dexp_arg_i(:,:), ddexp_arg_i(:,:), dexp_arg_j(:,:)
   integer :: i, j, nv, n_dof


   n_dof = size(x1)
   nv = size(x2,2)

   allocate(exp_arg(nv))
   allocate(dexp_arg_i(n_dof,nv))
   allocate(ddexp_arg_i(n_dof,nv))
   allocate(dexp_arg_j(n_dof,nv))
   exp_arg = 0.0_dp
   do i=1, n_dof
      if (periodicity(i) /= 0.0_dp) then
	 exp_arg(:) = exp_arg(:) - 2.0_dp*sin(TWO_PI/periodicity(i)*(x2(i,:)-x1(i))/2.0_dp)**2/len_scale_sq(i)
      else
	 exp_arg(:) = exp_arg(:) - 0.5_dp*(x2(i,:)-x1(i))**2/len_scale_sq(i)
      endif
   end do

   if (present(f_f)) f_f(:) = f_var * exp(exp_arg(:))
   if (present(g_f) .or. present(g_g)) then
      do i=1, n_dof
	 if (periodicity(i) /= 0.0_dp) then
	    dexp_arg_i(i,:) = 2.0_dp*sin(TWO_PI/periodicity(i)*(x2(i,:)-x1(i))/2.0_dp)*cos(TWO_PI/periodicity(i)* &
	       (x2(i,:)-x1(i))/2.0_dp)*TWO_PI/periodicity(i) / len_scale_sq(i)
	    if (present(g_g)) ddexp_arg_i(i,:) = (-sin(TWO_PI/periodicity(i)*(x2(i,:)-x1(i))/2.0_dp)**2+ &
	       cos(TWO_PI/periodicity(i)*(x2(i,:)-x1(i))/2.0_dp)**2)*(TWO_PI/periodicity(i))**2/len_scale_sq(i)
	 else
	    dexp_arg_i(i,:) = (x2(i,:)-x1(i))/len_scale_sq(i)
	    if (present(g_g)) ddexp_arg_i(i,:) = 1.0_dp/len_scale_sq(i)
	 endif
      end do
   endif
   if (present(f_g) .or. present(g_g)) then
      do j=1, n_dof
	 if (periodicity(j) /= 0.0_dp) then
	    dexp_arg_j(j,:) = -2.0_dp*sin(TWO_PI/periodicity(j)*(x2(j,:)-x1(j))/2.0_dp)*cos(TWO_PI/periodicity(j)* &
	       (x2(j,:)-x1(j))/2.0_dp)*TWO_PI/periodicity(j) / len_scale_sq(j)
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


subroutine GP_Matrix_Initialise(this,matrix)

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(in) :: matrix

  if(this%initialised) call gp_matrix_finalise(this)

  this%n = size(matrix,1)
  this%m = size(matrix,2)
  allocate(this%matrix(this%n,this%m), this%factor(this%n,this%m), this%s(this%n),this%tau(this%m) )

  this%matrix = matrix
  this%initialised = .true.

end subroutine GP_Matrix_Initialise

subroutine GP_Matrix_Update(this,matrix)

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(in) :: matrix

  integer :: factorised

  factorised = this%factorised

  if(this%initialised) then
     if( all(shape(matrix) == (/this%n,this%m/)) ) then
	this%matrix = matrix
     else
	call gp_matrix_initialise(this,matrix)
     endif
  else
     call gp_matrix_initialise(this,matrix)
  endif

  select case(factorised)
  case(QR)
     call GP_Matrix_QR_Factorise(this)
  endselect

end subroutine GP_Matrix_Update

subroutine GP_Matrix_Finalise(this)

  type(GP_Matrix), intent(inout) :: this

  if(.not. this%initialised) return

  this%n = 0
  this%m = 0
  if(allocated(this%matrix) ) deallocate(this%matrix)
  if(allocated(this%factor) ) deallocate(this%factor)
  if(allocated(this%s) ) deallocate(this%s)
  if(allocated(this%tau) ) deallocate(this%tau)
  this%initialised = .false.
  this%equilibrated = .false.
  this%factorised = NOT_FACTORISED

end subroutine GP_Matrix_Finalise


subroutine GP_Matrix_QR_Factorise(this,q,r)
  type(GP_Matrix), intent(inout) :: this         
  real(dp), dimension(:,:), intent(out), optional :: q, r

  integer :: lwork
  real(dp), allocatable :: work(:)
  integer :: info

  this%factor = this%matrix

  allocate(work(1))
  lwork = -1
  call dgeqrf(this%n, this%m, this%factor, this%n, this%tau, work, lwork, info)
  lwork = nint(work(1))
  deallocate(work)

  allocate(work(lwork))
  call dgeqrf(this%n, this%m, this%factor, this%n, this%tau, work, lwork, info)
  deallocate(work)

  if( info /= 0 ) then
     print *,'GP_Matrix_QR_Factorise: ',(-info),'-th parameter had an illegal value.'
     stop
  endif

  this%factorised = QR

  if( present(q) .or. present(r) ) call GP_Matrix_GetQR(this,q,r)

end subroutine GP_Matrix_QR_Factorise

subroutine GP_Matrix_GetQR(this,q,r)
  type(GP_Matrix), intent(inout) :: this         
  real(dp), dimension(:,:), intent(out), optional :: q, r

  integer :: lwork
  real(dp), allocatable :: work(:)
  integer :: j, info

  if( this%factorised /= QR ) then
     print *,'GP_Matrix_GetQR: not QR-factorised, call GP_Matrix_QR_Factorise first.'
     stop
  endif

  if(present(q)) then
     if (size(q,1) /= this%n .or. size(q,2) /= this%m) then
        print *, "GT_Matrix_GetQRshape(q) ",shape(q),"does not match",this%n,this%m
     endif
     q = this%factor

     allocate(work(1))
     lwork = -1
     call dorgqr(this%n, this%m, this%m, q, this%n, this%tau, work, lwork, info)
     lwork = nint(work(1))
     deallocate(work)

     allocate(work(lwork))
     call dorgqr(this%n, this%m, this%m, q, this%n, this%tau, work, lwork, info)
     deallocate(work)
  endif

  if(present(r)) then
     if (size(r,1) /= this%n .or. size(r,2) /= this%m) then
        print *, "GP_Matrix_GetQR shape(r) ",shape(r),"does not match",this%n,this%m
     endif
     r = this%factor(1:this%m,1:this%m)
     do j = 1, this%m - 1
	r(j+1:this%m,j) = 0.0_dp
     enddo
  endif

  if( info /= 0 ) then
     print *,'GP_Matrix_GetQR: ',(info),'-th parameter had an illegal value.'
     stop
  endif

end subroutine GP_Matrix_GetQR

subroutine GP_Matrix_QR_Solve_Matrix(this,matrix,result)
  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(in) :: matrix
  real(dp), dimension(:,:), intent(out) :: result

  real(dp), dimension(:,:), allocatable :: my_result
  integer :: info, i, j, n, o

  integer :: lwork
  real(dp), allocatable :: work(:)

  if(this%factorised == NOT_FACTORISED) then
     call GP_Matrix_QR_Factorise(this)
  elseif(this%factorised /= QR) then
     print *,'GP_Matrix_QR_Solve_Matrix: matrix not QR-factorised'
     stop
  endif

  n = size(matrix,1)
  o = size(matrix,2)
  if (size(result,1) /= this%m .or. size(result,2) /= o) then
     print *, "GP_Matrix_GetQR shape(result) ",shape(result),"does not match",this%m,o
  endif

  if( n /= this%n ) then
     print *,'GP_Matrix_QR_Solve_Matrix: dimensions of Q and matrix do not match.'
     stop
  endif

  allocate(my_result(n,o))
  my_result = matrix
  lwork = -1
  allocate(work(1))
  call dormqr('L', 'T', this%n, o, this%m, this%factor, this%n, this%tau, my_result, this%n, work, lwork, info)
  lwork = nint(work(1))
  deallocate(work)

  allocate(work(lwork))
  call dormqr('L', 'T', this%n, o, this%m, this%factor, this%n, this%tau, my_result, this%n, work, lwork, info)
  deallocate(work)

  if( info /= 0 ) then
     print *,'GP_Matrix_QR_QR_Solve_Matrix: ',(-info),'-th parameter had an illegal value.'
     stop
  endif

  do i = 1, o
     do j = this%m, 2, -1
	my_result(j,i) = my_result(j,i)/this%factor(j,j)
	my_result(1:j-1,i) = my_result(1:j-1,i) - my_result(j,i)*this%factor(1:j-1,j)
     enddo
     my_result(1,i) = my_result(1,i) / this%factor(1,1)
  enddo

  result = my_result(1:this%m,:)
  deallocate(my_result)

end subroutine GP_Matrix_QR_Solve_Matrix

subroutine GP_Matrix_QR_Solve_Vector(this,vector,result)
  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:), intent(in) :: vector
  real(dp), dimension(:), intent(out) :: result

  real(dp), dimension(:,:), allocatable :: my_result
  integer :: n, m

  n = size(vector)
  m = size(result)

  allocate(my_result(m,1))

  call GP_Matrix_QR_Solve_Matrix(this,reshape(vector,(/n,1/)),my_result)
  result = my_result(:,1)

  deallocate(my_result)

end subroutine GP_Matrix_QR_Solve_Vector

end module gp_basic_module
