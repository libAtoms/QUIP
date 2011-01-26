! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   libAtoms+QUIP: atomistic simulation library
! HND X
! HND X   Portions of this code were written by
! HND X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! HND X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! HND X
! HND X   Copyright 2006-2010.
! HND X
! HND X   Not for distribution
! HND X
! HND X   Portions of this code were written by Noam Bernstein as part of
! HND X   his employment for the U.S. Government, and are not subject
! HND X   to copyright in the USA.
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   http://www.libatoms.org
! HND X
! HND X  Additional contributions by
! HND X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X Gaussian Process module
!X
!% Module for general GP function interpolations.
!% A gp object contains the training set (teaching points and function values),
!% important temporary matrices, vectors and parameters.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module gp_teach_module

   use libatoms_module
   use fox_wxml
   use FoX_sax, only: xml_t, dictionary_t, haskey, getvalue, parse, &
   open_xml_string, close_xml_t
   use gp_predict_module

   implicit none

   integer, parameter :: GP_TEACH_MEMORY_0 = 0
   integer, parameter :: GP_TEACH_MEMORY_1 = 1
   integer, parameter :: GP_TEACH_MEMORY_2 = 2

   private

   real(qp), parameter :: gp_jitter  = 1.0e-6_qp

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !
   !% Gaussian Process data type.
   !% Contains all the important data necessary to estimate function values,
   !% derivatives, optimising hyperparameters.
   !
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   type gp_sparse

      !% Data arrays. These are (re)allocated quite rarely, always bigger than they should be.
      !% Therefore adding new data points is easy and fast.
      real(qp), dimension(:,:), allocatable :: x, x_prime, x_sparse, k_mn, big_raw_k_mn, big_k_mn, k_mm, inverse_k_mm_k_mn
      real(qp), dimension(:), allocatable   :: y, lambda
      integer, dimension(:), allocatable    :: l, lf, ldf, xf, xdf, xz, xz_sparse, sp, target_type

      !% Hyperparameters.
      !% Error parameter, range parameter, sensitivity parameters.
      real(qp), dimension(:,:), allocatable       :: theta
      real(qp), dimension(:), allocatable :: delta, f0, sigma
      !% Vector sizes.
      !% d: dimensionality of input space
      !% n: number of teaching points
      !% m: number of teaching function values
      !% p: number of hyperparameters
      integer :: d, n, m, sr, nx, nxx, nxd, mf, mdf, nsp, n_target_type

      integer :: gp_teach_memory = GP_TEACH_MEMORY_0

      logical :: initialised = .false.

      character(len=1024) :: label, comment

   end type gp_sparse

   type gp_minimise
      type(gp_sparse), pointer :: minim_gp => null()
      real(qp), dimension(:,:), pointer :: theta_0 => null()
      logical :: do_sigma = .false., do_delta = .false., do_theta = .false., do_sparx = .false., do_f0 = .false., do_theta_fac = .false.
      integer :: li_sigma = 1, ui_sigma = 1, li_delta = 1, ui_delta = 1, li_theta = 1, ui_theta = 1, li_sparx = 1, ui_sparx = 1, li_f0 = 1, ui_f0 = 1,  li_theta_fac = 1, ui_theta_fac = 1
   end type gp_minimise

   type gp_simple_minimise
      type(gp), pointer :: minim_gp => null()
      real(qp), dimension(:,:), pointer :: theta_0 => null()
      logical :: do_sigma = .false., do_delta = .false., do_theta = .false., do_f0 = .false., do_theta_fac = .false.
      integer :: li_sigma = 1, ui_sigma = 1, li_delta = 1, ui_delta = 1, li_theta = 1, ui_theta = 1, li_f0 = 1, ui_f0 = 1,  li_theta_fac = 1, ui_theta_fac = 1
   end type gp_simple_minimise

   interface Initialise
      module procedure GP_Initialise, gp_simple_initialise
   end interface Initialise

   interface Finalise
      module procedure GP_Finalise_sparse, GP_Finalise_sparse_more
   end interface Finalise

   interface apply_l
      module procedure apply_l_matrix, apply_l_vector
   end interface apply_l

   interface covSEard
      module procedure covSEard_dp
#ifdef HAVE_QP
      module procedure covSEard_qp
#endif      
   end interface covSEard

   public :: gp_sparse, initialise, finalise, gp_update, &
   likelihood, test_gp_gradient, test_gp_simple_gradient, minimise_gp_gradient, minimise_gp_simple_gradient, &
   minimise_gp_ns, gp_sparsify, &
   gp_print_binary, gp_read_binary, minimise_gp_ns_new, fill_random_integer

   contains

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Initialise a gp object from scratch.
      !% Inputs needed:
      !% theta_in: hyperparameters
      !% y_in: function or derivative values (have to stick to one of them here)
      !% x_in: teaching points
      !% x_prime_in: optional, derivatives of teaching points
      !% l_in: if y values are sum of function values, this vector contains information about which points contribute.
      !% l_in is an integer array, same size as y. Each element represents the index of the \emph{last} contributing
      !% teaching vector in x_in. E.g.: if y(3) is the sum of the function to be learnt at the teaching points
      !% x_in(:,4), x_in(:,5), x_in(:,6) then l(2) = 3, l(3) = 6.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine GP_initialise(this, sparse)

         !% n: number of teaching points \\
         !% d: dimension of teaching vectors \\
         !% p: number of parameters \\

         type(gp), intent(inout)     :: this        !% gp to initialise
         type(gp_sparse), intent(in) :: sparse

         type(LA_matrix) :: LA_k_mm, LA_q_mm
#ifdef HAVE_QR
         integer :: i, j, n, m, info
         real(qp), dimension(:,:), allocatable :: a, k_mn_sq_inverse_lambda, factor_k_mm
         real(qp), dimension(:), allocatable :: y, alpha
#else
         real(qp), dimension(:,:), allocatable :: k_mn_inverse_lambda, k_mn_l_k_nm, inverse_q_mm, inverse_k_mm
         real(qp), dimension(:), allocatable :: alpha
         integer :: i, j, n, m, info
#endif         

         this%d = sparse%d
         this%n = sparse%sr
         this%nsp = sparse%nsp

         n = sparse%m
         m = sparse%sr

         allocate( this%x(this%d,this%n), this%alpha(this%n), this%theta(this%d,this%nsp), &
         & this%delta(this%nsp), this%f0(this%nsp), this%c(this%n,this%n) )
         allocate( this%xz(this%n), this%sp(this%nsp) )

         this%sigma = real(sparse%sigma(1),kind=dp) ! temporary
         this%delta = real(sparse%delta,kind=dp)
         this%f0    = real(sparse%f0   ,kind=dp)
         this%theta = real(sparse%theta,kind=dp)
         this%x = real(sparse%x_sparse,kind=dp)
         this%xz = sparse%xz_sparse
         this%sp = sparse%sp

	 call setup_x_div_theta(this)

#ifdef HAVE_QR
         allocate( k_mn_sq_inverse_lambda(m,n), factor_k_mm(m,m), a(n+m,m), y(n+m), alpha(m) )
         call matrix_product_vect_asdiagonal_sub(k_mn_sq_inverse_lambda,sparse%k_mn,sqrt(1.0_qp/sparse%lambda)) ! O(NM)

         call initialise(LA_k_mm,sparse%k_mm)
         call LA_Matrix_Factorise(LA_k_mm,factor_k_mm,info=info)
         if( info /= 0 ) call system_abort('GP_initialise: LA_k_mm')
         call finalise(LA_k_mm)

         do i = 1, m-1
            do j = i+1, m
               factor_k_mm(j,i) = 0.0_qp
            end do
         end do
         
         a(1:n,:) = transpose(k_mn_sq_inverse_lambda)
         a(n+1:,:) = factor_k_mm

         y = 0.0_dp
         y(1:n) = sparse%y*sqrt(1.0_qp/sparse%lambda)

         call initialise(LA_q_mm,a)
         call LA_Matrix_QR_Solve_Vector(LA_q_mm,y,alpha)
         this%alpha = real(alpha,dp)
         call finalise(LA_q_mm)
         deallocate( k_mn_sq_inverse_lambda, factor_k_mm, a, y, alpha)
         this%c = 0.0_dp
#else
         allocate( k_mn_inverse_lambda(m,n), k_mn_l_k_nm(m,m), inverse_q_mm(m,m), inverse_k_mm(m,m), alpha(m))
         call matrix_product_vect_asdiagonal_sub(k_mn_inverse_lambda,sparse%k_mn,1.0_qp/sparse%lambda) ! O(NM)
         k_mn_l_k_nm = matmul(k_mn_inverse_lambda,transpose(sparse%k_mn))

         call initialise(LA_q_mm,sparse%k_mm + k_mn_l_k_nm)
         call LA_Matrix_Inverse(LA_q_mm,inverse_q_mm,info=info)
         if( info /= 0 ) call system_abort('GP_initialise: LA_q_mm')

         !this%alpha = real( matmul(inverse_q_mm,matmul(k_mn_inverse_lambda, sparse%y)),kind=dp)          ! O(NM)
         call Matrix_Solve(LA_q_mm,matmul(k_mn_inverse_lambda, sparse%y),alpha)
         call finalise(LA_q_mm)
         this%alpha = real(alpha,kind=dp)

         call initialise(LA_k_mm,sparse%k_mm)
         call LA_Matrix_Inverse(LA_k_mm,inverse_k_mm,info=info)
         if( info /= 0 ) call system_abort('GP_initialise: LA_k_mm')

         call finalise(LA_k_mm)
         this%c = real(inverse_k_mm - inverse_q_mm,kind=dp)                                                            ! O(M^2)
         deallocate( k_mn_inverse_lambda, k_mn_l_k_nm, inverse_q_mm, inverse_k_mm, alpha )
#endif
         
         this%initialised = .true.

      end subroutine GP_initialise

      subroutine gp_simple_initialise(this,sigma_in, delta_in, theta_in, f0_in, y_in, x_in)

         type(gp), intent(inout)                        :: this
         real(dp), intent(in)                           :: sigma_in, delta_in, f0_in
         real(dp), dimension(:,:), intent(in)           :: theta_in    !% hyperparameters, dimension(p)
         real(dp), dimension(:), intent(in)             :: y_in        !% function (or derivative) value at teaching points, dimension(m)
         real(dp), dimension(:,:), intent(in)           :: x_in        !% teaching points, dimension(d,n)

         real(dp) :: mem_required, mem_total, mem_free

         if(this%initialised) call finalise(this)

         this%d = size(x_in,1)
         this%n = size(x_in,2) 
         this%nsp = 1

         if( this%n == 0 ) call system_abort('gp_simple_initialise: n = 0, nothing to teach.')

         call check_size('theta_in',theta_in,(/this%d,1/),'gp_simple_initialise')
         call check_size('y_in',y_in,(/this%n/),'gp_simple_initialise')

         mem_required = 3.0_dp*real(this%n,dp)**2 * real(dp,dp) / (1024.0_dp**3)
         call mem_info(mem_total,mem_free)
         mem_total = mem_total / (1024.0_dp**3)
         mem_free = mem_free / (1024.0_dp**3)

         call print('Memory required (approx.): '//mem_required//' GB')
         if( mem_required > mem_total ) call system_abort('Required memory ('//mem_required//' GB) exceeds available memory ('//mem_total//' GB).')

         allocate(this%x(this%d,this%n), this%alpha(this%n), this%y(this%n), this%x_div_theta(this%d,this%n), &
         this%c(this%n,this%n), this%theta(this%d,1), this%delta(1), this%f0(1), this%xz(this%n), this%sp(1))

         this%x = x_in
         this%sigma = sigma_in
         this%delta = delta_in
         this%theta = theta_in
         this%f0 = f0_in
         this%y = y_in

         this%xz = 0
         this%sp = 0

         call covariance_matrix_simple(this)

         this%initialised = .true.

      endsubroutine gp_simple_initialise

      subroutine GP_sparsify(this,r,sigma_in, delta_in, theta_in, y_in, yd_in, x_in, x_prime_in, xf_in, xdf_in, lf_in, ldf_in, xz_in,sp_in, f0_in, target_type_in, gp_teach_memory_in)
         type(gp_sparse), intent(inout)                 :: this        !% gp to initialise
         !integer, intent(in)                            :: sr
         integer, dimension(:), intent(in)              :: r
         real(dp), dimension(:), intent(in)             :: sigma_in
         real(dp), dimension(:), intent(in)             :: delta_in    !% hyperparameters, dimension(p)
         real(dp), dimension(:,:), intent(in)           :: theta_in    !% hyperparameters, dimension(p)
         real(dp), dimension(:), intent(in)             :: y_in        !% function (or derivative) value at teaching points, dimension(m)
         real(dp), dimension(:), intent(in)             :: yd_in        !% function (or derivative) value at teaching points, dimension(m)
         real(dp), dimension(:,:), intent(in)           :: x_in        !% teaching points, dimension(d,n)
         real(dp), dimension(:,:), intent(in)           :: x_prime_in  !% teaching points, dimension(d,n)
         integer, dimension(:), intent(in)              :: xf_in        !% how the sum is formed from teaching point(m)
         integer, dimension(:), intent(in)              :: xdf_in        !% how the sum is formed from teaching point(m)
         integer, dimension(:), intent(in)              :: lf_in        !% how the sum is formed from teaching point(m)
         integer, dimension(:), intent(in)              :: ldf_in        !% how the sum is formed from teaching point(m)
         integer, dimension(:), intent(in), optional    :: xz_in        !% how the sum is formed from teaching point(m)
         integer, dimension(:), intent(in), optional    :: sp_in        !% how the sum is formed from teaching point(m)
         real(dp), dimension(:), intent(in), optional   :: f0_in    !% hyperparameters, dimension(p)
         integer, dimension(:), intent(in), optional    :: target_type_in        !% what type of target is the function value
         integer, intent(in), optional                  :: gp_teach_memory_in

         integer :: i, lmf, nsp
         integer, dimension(:), allocatable :: my_target_type
         !integer, dimension(:), allocatable             :: r
         !integer, intent(in)                            :: sr

         if( this%initialised ) return
         call print_title('Sparse GP initialisation')

         this%d   = size(x_in,1) 
         this%nxx = size(x_in,2)
         this%nx  = size(xf_in)
         this%nxd = size(x_prime_in,2)
         this%m   = size(y_in)+size(yd_in)
         this%mf  = size(lf_in)
         this%mdf = size(ldf_in)
         this%sr  = size(r)
         this%n_target_type = size(sigma_in)
         
         if(present(sp_in)) then
            this%nsp = size(sp_in)
         else
            this%nsp = 1
         end if

         this%n = this%nx + this%nxd

         !allocate(r(sr))
         !call fill_random_integer(r,this%n)

         call print('Number of function values to sparsify: '//this%m)
         call print('Number of environments: '//this%nxx)
         call print('Sparse teaching points: '//this%sr)
         call print('Dimensions: '//this%d)

         if( (size(lf_in,1) + size(ldf_in,1)) /= this%m ) call system_abort('GP_Initialise: y_in and l_in do not conform')

         allocate(my_target_type(this%m))
         if( present( target_type_in ) ) then
            call check_size('target_type_in',target_type_in,(/this%m/),'GP_initialise')
            if( maxval(target_type_in) > this%n_target_type ) &
            & call system_abort('GP_initialise: target_type_in contains more types ('//maxval(target_type_in)//') than the dimensionality of sigma_in ('//this%n_target_type//')')
            if( minval(target_type_in) < 1 ) call system_abort('GP_initialise: target_type_in contains zeros')
            my_target_type = target_type_in
         else
            if( this%n_target_type < 2 ) &
            & call system_abort('GP_initialise: the dimensionality of sigma_in is '//this%n_target_type//', less than 2, cannot associate sigmas to different target function values')
            my_target_type(:this%mf) = 1
            my_target_type(this%mf+1:) = 2
         endif

         if(present(gp_teach_memory_in)) then
            this%gp_teach_memory = gp_teach_memory_in
         else
            this%gp_teach_memory = GP_TEACH_MEMORY_0
         endif

         call gp_allocate(this)

         this%sigma = real(sigma_in,kind=qp)
         this%delta = real(delta_in,kind=qp)
         this%theta = real(theta_in,kind=qp)
         if( this%mf>0 ) this%y(1:this%mf) = real(y_in,kind=qp)
         !this%y(this%mf+1:this%mdf) = yd_in
         if( this%mdf>0 ) this%y(this%mf+1:this%m) = real(yd_in,kind=qp)
         this%x = real(x_in,kind=qp)
         this%x_prime = real(x_prime_in,kind=qp)
         this%lf = lf_in
         this%ldf = ldf_in
         if(present(xz_in)) then
            this%xz = xz_in
         else
            this%xz = 0
         end if

         if( this%mf>0 ) then
            this%l(:this%mf) = lf_in
            lmf = lf_in(this%mf)
         else
            lmf = 0
         end if
         if( this%mdf>0 ) this%l(this%mf+1:) = ldf_in + lmf

         this%xf = xf_in
         this%xdf = xdf_in

         this%x_sparse = real(x_in(:,r),kind=qp)
         this%xz_sparse = this%xz(r)
         if(present(sp_in)) then
            this%sp = sp_in
         else
            this%sp = 1
         end if

         if(present(f0_in)) then
            this%f0 = f0_in
         else
            this%f0 = 0.0_qp
         end if
         this%target_type = my_target_type
         
         call covariance_matrix_sparse(this)

         this%initialised = .true.
         deallocate(my_target_type)

      end subroutine GP_sparsify

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Allocates data arrays, storing the $\mathbf{x}$ vectors or similarly sized other arrays.
      !% These are the `n' dimensional arrays.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_allocate(this)
         type(gp_sparse), intent(inout)       :: this

         real(dp) :: mem_required, mem_total, mem_free

         call gp_deallocate(this)

         if(this%gp_teach_memory == GP_TEACH_MEMORY_0) &
            mem_required = 0.0_dp

         if(this%gp_teach_memory == GP_TEACH_MEMORY_1) &
            mem_required = real(this%sr,dp) * real(this%n,dp) * real(dp,dp) / (1024.0_dp**3)

         if(this%gp_teach_memory == GP_TEACH_MEMORY_2) &
            mem_required = 2.0_dp * real(this%sr,dp) * real(this%n,dp) * real(dp,dp) / (1024.0_dp**3)

         call mem_info(mem_total,mem_free)
         mem_total = mem_total / (1024.0_dp**3)
         mem_free = mem_free / (1024.0_dp**3)

         call print('Memory required (approx.): '//mem_required//' GB')
         if( mem_required > mem_total ) call system_abort('Required memory ('//mem_required//' GB) exceeds available memory ('//mem_total//' GB).')


         allocate( this%sigma(this%n_target_type), this%theta(this%d,this%nsp), this%delta(this%nsp), this%f0(this%nsp), this%sp(this%nsp), this%target_type(this%m) )

         allocate( this%x(this%d,this%nxx), this%x_prime(this%d,this%nxd), this%xf(this%nx), this%xdf(this%nxd) )

         allocate( this%y(this%m), this%lf(this%mf), this%ldf(this%mdf), this%l(this%m), this%lambda(this%m), &
         this%x_sparse(this%d,this%sr), this%k_mn(this%sr,this%m), this%k_mm(this%sr,this%sr), this%inverse_k_mm_k_mn(this%sr,this%m))

         if(this%gp_teach_memory >= GP_TEACH_MEMORY_1) then
            allocate(this%big_k_mn(this%sr,this%n))
            this%big_k_mn = 0.0_qp
         endif

         if(this%gp_teach_memory >= GP_TEACH_MEMORY_2) then
            allocate(this%big_raw_k_mn(this%sr,this%n))
            this%big_raw_k_mn = 0.0_qp
         endif

         allocate(this%xz(this%nxx), this%xz_sparse(this%sr) )
         this%x     = 0.0_qp
         this%x_prime = 0.0_qp
         this%y     = 0.0_qp
         this%xf    = 0
         this%xdf   = 0
         this%lf    = 0
         this%ldf   = 0
         this%l     = 0
         this%x_sparse = 0.0_qp
         this%k_mn = 0.0_qp
         this%k_mm = 0.0_qp

         this%xz = 0
         this%xz_sparse = 0
         this%target_type = 0

      end subroutine gp_allocate

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Deallocates data arrays, storing the $\mathbf{x}$ vectors or similarly sized other arrays.
      !% These are the `n' dimensional arrays.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_deallocate(this)
         type(gp_sparse), intent(inout)       :: this

         if(allocated(this%x))            deallocate( this%x)
         if(allocated(this%x_prime))            deallocate( this%x_prime)

         if(allocated(this%y))        deallocate( this%y)
         if(allocated(this%lambda))   deallocate( this%lambda)
         if(allocated(this%lf))       deallocate( this%lf)
         if(allocated(this%ldf))      deallocate( this%ldf)
         if(allocated(this%l))        deallocate( this%l)
         if(allocated(this%xf))       deallocate( this%xf)
         if(allocated(this%xdf))      deallocate( this%xdf)
         if(allocated(this%x_sparse)) deallocate( this%x_sparse )
         if(allocated(this%k_mn))     deallocate( this%k_mn )
         if(allocated(this%k_mm))     deallocate( this%k_mm )
         if(allocated(this%big_k_mn)) deallocate( this%big_k_mn )
         if(allocated(this%big_raw_k_mn)) deallocate( this%big_raw_k_mn )
         if(allocated(this%inverse_k_mm_k_mn)) deallocate( this%inverse_k_mm_k_mn )

         if(allocated(this%xz)) deallocate( this%xz )
         if(allocated(this%xz_sparse)) deallocate( this%xz_sparse )
         if(allocated(this%sp)) deallocate( this%sp )
         if(allocated(this%sigma)) deallocate( this%sigma )
         if(allocated(this%theta)) deallocate( this%theta )
         if(allocated(this%delta)) deallocate( this%delta )
         if(allocated(this%f0)) deallocate( this%f0 )
         if(allocated(this%target_type)) deallocate( this%target_type )

      end subroutine gp_deallocate

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Changes one or more component in a gp object, updating temporary matrices and vectors.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_update(this, x_in, x_sparse_in, y_in, sigma_in, delta_in, theta_in, f0_in, do_covariance)

         type(gp_sparse), intent(inout)                 :: this       !% gp object to update
         real(qp), dimension(:,:), intent(in), optional :: x_in       !% update teaching points, dimension(d,n)
         real(qp), dimension(:,:), intent(in), optional :: x_sparse_in !% update derivatives of teaching points, dimension(d,n)
         real(qp), dimension(:), intent(in), optional   :: y_in       !% update function values, dimension(m)
         real(qp), dimension(:), intent(in), optional   :: sigma_in   !% update hyperparameters, dimension(p)
         real(qp), dimension(:), intent(in), optional   :: delta_in   !% update hyperparameters, dimension(p)
         real(qp), dimension(:,:), intent(in), optional :: theta_in   !% update hyperparameters, dimension(p)
         real(qp), dimension(:), intent(in), optional :: f0_in   !% update hyperparameters, dimension(p)
         logical, intent(in), optional :: do_covariance

         logical :: my_do_covariance

         if( .not. this%initialised ) return

         my_do_covariance = optional_default(.true.,do_covariance)

         if( present(x_in) ) then
             if( all( shape(this%x) /= shape(x_in) ) ) call system_abort('gp_update: array sizes do not conform')
             this%x = x_in
         end if

         if( present(x_sparse_in) ) then
             if( all( shape(this%x_sparse) /= shape(x_sparse_in) ) ) call system_abort('gp_update: array sizes do not conform')
             this%x_sparse = x_sparse_in
         end if

         if( present(y_in) ) then
             if( all( shape(this%y) /= shape(y_in) ) ) call system_abort('gp_update: array sizes do not conform')
             this%y = y_in
         end if

         if( present(sigma_in) ) then
            call check_size('sigma_in',sigma_in,shape(this%sigma),'gp_update')
            this%sigma = sigma_in
         endif

         if( present(delta_in) ) then
            call check_size('delta_in',delta_in,shape(this%delta),'gp_update')
            this%delta = delta_in
         endif

         if( present(f0_in) ) then
            call check_size('f0_in',f0_in,shape(this%f0),'gp_update')
            this%f0 = f0_in
         endif

         if( present(theta_in) ) then
             if( all( shape(this%theta) /= shape(theta_in) ) ) call system_abort('gp_update: array sizes do not conform')
             this%theta = theta_in
         end if

         if( my_do_covariance ) then
            if( present(x_in) .or. present(y_in) .or. present(x_sparse_in) .or. &
            & present(sigma_in) .or. present(delta_in) .or. present(theta_in) .or. present(f0_in) ) &
            & call covariance_matrix_sparse(this)
         end if

      end subroutine gp_update

      subroutine gp_update_simple(this, x_in, y_in, sigma_in, delta_in, theta_in, f0_in, do_covariance)
         type(gp), intent(inout)                        :: this       !% gp object to update
         real(qp), dimension(:,:), intent(in), optional :: x_in       !% update teaching points, dimension(d,n)
         real(qp), dimension(:), intent(in), optional   :: y_in       !% update function values, dimension(m)
         real(qp), intent(in), optional                 :: sigma_in   !% update hyperparameters, dimension(p)
         real(qp), dimension(:), intent(in), optional   :: delta_in   !% update hyperparameters, dimension(p)
         real(qp), dimension(:,:), intent(in), optional :: theta_in   !% update hyperparameters, dimension(p)
         real(qp), dimension(:), intent(in), optional   :: f0_in   !% update hyperparameters, dimension(p)
         logical, intent(in), optional :: do_covariance

         logical :: my_do_covariance

         if( .not. this%initialised ) return

         my_do_covariance = optional_default(.true.,do_covariance)

         if( present(x_in) ) then
             if( all( shape(this%x) /= shape(x_in) ) ) call system_abort('gp_update: array sizes do not conform')
             this%x = x_in
         end if

         if( present(y_in) ) then
             if( all( shape(this%y) /= shape(y_in) ) ) call system_abort('gp_update: array sizes do not conform')
             this%y = y_in
         end if

         if( present(sigma_in) ) then
            this%sigma = sigma_in
         endif

         if( present(delta_in) ) then
            call check_size('delta_in',delta_in,shape(this%delta),'gp_update')
            this%delta = delta_in
         endif

         if( present(f0_in) ) then
            call check_size('f0_in',f0_in,shape(this%f0),'gp_update')
            this%f0 = f0_in
         endif

         if( present(theta_in) ) then
             if( all( shape(this%theta) /= shape(theta_in) ) ) call system_abort('gp_update: array sizes do not conform')
             this%theta = theta_in
         end if

         if( my_do_covariance ) then
            if( present(x_in) .or. present(y_in) .or. &
            present(sigma_in) .or. present(delta_in) .or. present(theta_in) .or. present(f0_in) ) &
            call covariance_matrix_simple(this)
         end if

      endsubroutine gp_update_simple

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Finalise a gp object.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine GP_finalise_sparse(this)

         type(gp_sparse), intent(inout)              :: this !% gp object

         if( .not. this%initialised ) return

         call gp_deallocate(this)

         this%d = 0
         this%nxx =  0
         this%nx  = 0
         this%nxd = 0
         this%n = 0
         this%m = 0
         this%mf  = 0
         this%mdf = 0
         this%n_target_type = 0
         this%initialised = .false.

      end subroutine GP_finalise_sparse

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Finalise more gp objects.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine GP_finalise_sparse_more(this, this2, this3, this4, this5, this6, this7, this8, this9)

         type(gp_sparse), intent(inout)           :: this, this2
         type(gp_sparse), intent(inout), optional :: this3, this4, this5, this6, this7, this8, this9

         call GP_finalise_sparse(this)
         call GP_finalise_sparse(this2)

         if( present(this3) ) call GP_finalise_sparse(this3)
         if( present(this4) ) call GP_finalise_sparse(this4)
         if( present(this5) ) call GP_finalise_sparse(this5)
         if( present(this6) ) call GP_finalise_sparse(this6)
         if( present(this7) ) call GP_finalise_sparse(this7)
         if( present(this8) ) call GP_finalise_sparse(this8)
         if( present(this9) ) call GP_finalise_sparse(this9)

      end subroutine GP_finalise_sparse_more

      pure function covSEard_dp(delta,theta,xi,xj,dxj)

         ! Covariance function

         real(dp) :: covSEard_dp
         real(dp), intent(in) :: delta
         real(dp), dimension(:), intent(in) :: xi,xj,theta
         real(dp), dimension(:), intent(in), optional :: dxj
         real(dp), dimension(:), allocatable :: xixjtheta

         allocate(xixjtheta(size(xi)))
         xixjtheta = (xi-xj)/theta

         covSEard_dp = delta**2 * exp( - 0.5_dp * normsq(xixjtheta) ) 
         if(present(dxj)) covSEard_dp = covSEard_dp * dot_product(xixjtheta,dxj/theta)

         deallocate(xixjtheta)
      end function covSEard_dp

      pure function covSEard_qp(delta,theta,xi,xj,dxj)

         ! Covariance function

         real(qp) :: covSEard_qp
         real(qp), intent(in) :: delta
         real(qp), dimension(:), intent(in) :: xi,xj,theta
         real(qp), dimension(:), intent(in), optional :: dxj
         real(qp), dimension(:), allocatable :: xixjtheta

         allocate(xixjtheta(size(xi)))
         xixjtheta = (xi-xj)/theta

         covSEard_qp = delta**2 * exp( - 0.5_qp * normsq(xixjtheta) ) 
         if(present(dxj)) covSEard_qp = covSEard_qp * dot_product(xixjtheta,dxj/theta)
         deallocate(xixjtheta)
      end function covSEard_qp

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Updates covariance matrix and other important temporary arrays in a GP object.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine covariance_matrix_sparse(this)

         type(gp_sparse), intent(inout) :: this

         type(LA_Matrix) :: LA_k_mm

         integer :: i, j, lj, uj, j1, j2, id, jd, xj, xj1, xj2, k, Z_type, info
         real(qp) :: ep, big_raw_k_mn_ij, big_k_mn_ij
         real(qp), dimension(:,:), allocatable :: inverse_k_mm
         real(qp), dimension(:), allocatable   :: lknnl, diff_xijt
         real(qp), dimension(:,:), allocatable   :: theta2

         allocate(lknnl(this%m), inverse_k_mm(this%sr,this%sr), theta2(this%d,this%nsp)) ! , diff_xijt(this%d) )
         
         theta2 = 1.0_qp / this%theta**2
         this%k_mn = 0.0_qp

!$omp parallel do private(xj,k,Z_type,lj,i,big_raw_k_mn_ij,big_k_mn_ij)
         do j = 1, this%nx
            xj = this%xf(j)
            lj = count(j > this%l) + 1
            
            do k = 1, this%nsp; if( this%sp(k) == this%xz(xj) ) Z_type = k; end do

            do i = 1, this%sr
               if( this%xz(xj) == this%xz_sparse(i) ) then
                  big_raw_k_mn_ij = covSEard( this%delta(Z_type), this%theta(:,Z_type), this%x(:,xj), this%x_sparse(:,i) )
                  big_k_mn_ij = big_raw_k_mn_ij + this%f0(Z_type)**2 

                  if(this%gp_teach_memory >= GP_TEACH_MEMORY_2) this%big_raw_k_mn(i,j) = big_raw_k_mn_ij
                  if(this%gp_teach_memory >= GP_TEACH_MEMORY_1) this%big_k_mn(i,j) = big_k_mn_ij
!$omp critical                  
                  if(this%gp_teach_memory == GP_TEACH_MEMORY_0) this%k_mn(i,lj) = this%k_mn(i,lj) + big_k_mn_ij
!$omp end critical                  
               else
                  if(this%gp_teach_memory >= GP_TEACH_MEMORY_2) this%big_raw_k_mn(i,j) = 0.0_qp
                  if(this%gp_teach_memory >= GP_TEACH_MEMORY_1) this%big_k_mn(i,j) = 0.0_qp
               end if
              !& dot_product(( this%x_sparse(:,i) - this%x(:,j) )*theta2,this%x_prime(:,j))


            end do
         end do

!$omp parallel do private(jd,xj,k,Z_type,lj,i,big_raw_k_mn_ij,big_k_mn_ij)
         do j = 1, this%nxd
            jd = j + this%nx
            xj = this%xdf(j)
            lj = count(jd > this%l) + 1
            
            do k = 1, this%nsp; if( this%sp(k) == this%xz(xj) ) Z_type = k; end do
            
            do i = 1, this%sr
               if( this%xz(xj) == this%xz_sparse(i) ) then
                  big_raw_k_mn_ij = covSEard( this%delta(Z_type), this%theta(:,Z_type), this%x(:,xj), this%x_sparse(:,i) )
                  big_k_mn_ij = big_raw_k_mn_ij * &
                  dot_product(( this%x_sparse(:,i) - this%x(:,xj) )*theta2(:,Z_type),this%x_prime(:,j))

                  if(this%gp_teach_memory >= GP_TEACH_MEMORY_2) this%big_raw_k_mn(i,jd) = big_raw_k_mn_ij
                  if(this%gp_teach_memory >= GP_TEACH_MEMORY_1) this%big_k_mn(i,jd) = big_k_mn_ij
!$omp critical                  
                  if(this%gp_teach_memory == GP_TEACH_MEMORY_0) this%k_mn(i,lj) = this%k_mn(i,lj) + big_k_mn_ij
!$omp end critical                  

              else

                 if(this%gp_teach_memory >= GP_TEACH_MEMORY_2) this%big_raw_k_mn(i,jd) = 0.0_qp
                 if(this%gp_teach_memory >= GP_TEACH_MEMORY_1) this%big_k_mn(i,jd) = 0.0_qp

              end if
            end do
         end do

!$omp parallel do private(k,Z_type,i)
         do j = 1, this%sr
         
            do k = 1, this%nsp; if( this%sp(k) == this%xz_sparse(j) ) Z_type = k; end do

            do i = 1, this%sr
               if( this%xz_sparse(j) == this%xz_sparse(i) ) then
                  this%k_mm(i,j) = covSEard( this%delta(Z_type), this%theta(:,Z_type), this%x_sparse(:,j), this%x_sparse(:,i) ) + this%f0(Z_type)**2
               else
                  this%k_mm(i,j) = 0.0_qp
               end if
            end do
         end do

         ep = gp_jitter !*trace(this%k_mm)/this%sr
         do i = 1, this%sr
            this%k_mm(i,i) = this%k_mm(i,i) + ep
         end do

         lknnl = 0.0_qp

!$omp do private(i,lj,uj,j1,j2,xj1,xj2,k,Z_type)
         do i = 1, this%mf
            if(i==1) then 
               lj = 1
            else
               lj = this%lf(i-1)+1
            end if
            uj = this%lf(i)
            do j1 = lj, uj
               xj1 = this%xf(j1)

               do k = 1, this%nsp; if( this%sp(k) == this%xz(xj1) ) Z_type = k; end do

               do j2 = lj, uj
                  xj2 = this%xf(j2)
                  if( this%xz(xj1) == this%xz(xj2) ) &
                  & lknnl(i) = lknnl(i) + covSEard( this%delta(Z_type), this%theta(:,Z_type), this%x(:,xj1), this%x(:,xj2) )+ this%f0(Z_type)**2
               end do
            end do
         end do
!$omp end do 

!$omp parallel private(diff_xijt)
allocate(diff_xijt(this%d))
!$omp do private(i,id,lj,uj,j1,j2,xj1,xj2,k,Z_type)
         do i = 1, this%mdf
            id = i + this%mf
            if(i==1) then 
               lj = 1
            else
               lj = this%ldf(i-1)+1
            end if
            uj = this%ldf(i)
            do j1 = lj, uj
               xj1 = this%xdf(j1)
              
               do k = 1, this%nsp; if( this%sp(k) == this%xz(xj1) ) Z_type = k; end do
              
               do j2 = lj, uj
                  xj2 = this%xdf(j2)
                  if( this%xz(xj1) == this%xz(xj2) ) then
                     diff_xijt = (this%x(:,xj1) - this%x(:,xj2)) * theta2(:,Z_type)
                     lknnl(id) = lknnl(id) + covSEard( this%delta(Z_type), this%theta(:,Z_type), this%x(:,xj1), this%x(:,xj2) ) * &
                     & ( dot_product( this%x_prime(:,j1) * theta2(:,Z_type), this%x_prime(:,j2) ) - &
                     & dot_product( diff_xijt,this%x_prime(:,j1) ) * dot_product( diff_xijt,this%x_prime(:,j2) ) )
                  end if
               end do
            end do
         end do
!$omp end do 
deallocate(diff_xijt)
!$omp end parallel

         if(this%gp_teach_memory >= GP_TEACH_MEMORY_1) call apply_l(this%big_k_mn,this%l,this%k_mn)
         !call inverse(this%k_mm,inverse_k_mm)
         !inverse_k_mm = this%k_mm
         call initialise(LA_k_mm,this%k_mm)
         call Matrix_Solve(LA_k_mm,this%k_mn,this%inverse_k_mm_k_mn,info=info)
         if( info /= 0 ) call system_abort('covariance_matrix_sparse: LA_k_mm')
         call finalise(LA_k_mm)
         
!         call factorise_inverse_symmetric_matrix(this%k_mm,A_inverse=inverse_k_mm)            ! O(M^3)
!         call matrix_product_sub(this%inverse_k_mm_k_mn, inverse_k_mm, this%k_mn)
         do i = 1, this%m
            this%lambda(i) = lknnl(i) - dot_product( this%inverse_k_mm_k_mn(:,i), this%k_mn(:,i) ) + this%sigma(this%target_type(i))**2
         end do

         deallocate( inverse_k_mm, lknnl, theta2) !, diff_xijt )

      end subroutine covariance_matrix_sparse

      subroutine covariance_matrix_simple(this)
         type(gp), intent(inout) :: this

         integer :: i, j
         real(qp), dimension(:), allocatable :: xixjtheta
         real(qp), dimension(:,:), allocatable :: c
         real(qp), dimension(:), allocatable :: local_alpha

         allocate(local_alpha(this%N))
	 call setup_x_div_theta(this)
         allocate(xixjtheta(this%d),c(this%n,this%n))

         do i = 1, this%n
            do j = i + 1, this%n
               xixjtheta = this%x_div_theta(:,i) - this%x_div_theta(:,j)
               c(j,i) = exp(-0.5_dp * dot_product(xixjtheta,xixjtheta))
            enddo
            c(i,i) = 1.0_dp + this%sigma**2 / this%delta(1)**2
         enddo

         c = c * this%delta(1)**2

         do i = 1, this%n
            do j = 1, i - 1
               c(j,i) = c(i,j)
            enddo
         enddo

         call LA_Matrix_Update(this%LA_C_nn,c)

#ifdef HAVE_QR
         if(this%LA_C_nn%factorised /= QR) call LA_Matrix_QR_Factorise(this%LA_C_nn)
         call LA_Matrix_QR_Solve_Vector(this%LA_C_nn,real(this%y,qp),local_alpha)
#else
         if(this%LA_C_nn%factorised /= CHOLESKY) call LA_Matrix_Factorise(this%LA_C_nn)
         call LA_Matrix_Solve_Vector(this%LA_C_nn,real(this%y, qp),local_alpha)
#endif
         this%alpha = real(local_alpha, dp)
         deallocate(c,xixjtheta, local_alpha)
         
      endsubroutine covariance_matrix_simple
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% In case of sum of values as teaching points, apply the operation
      !% $ \mathbf{LCL}^\dagger $ or $ \mathbf{Lc} $, where $ \mathbf{L} $
      !% represents the linear combination matrix, here substituted by l integer array as we care about sums only.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine apply_l_matrix(x,l,lx)
         real(qp), dimension(:,:), intent(in)  :: x
         integer, dimension(:), intent(in)     :: l
         real(qp), dimension(:,:), intent(out) :: lx 

         integer :: j, uj, lj, m, n, sr

         m = size(l)
         if( size(x,1) /= size(lx,1) )  call system_abort('apply_l_matrix: x and lx do not conform')
         if( m /= size(lx,2) ) call system_abort('apply_llt_matrix: l and lx do not conform')

         sr = size(x,1)
         n = size(x,2)

         if( m /= n ) then
            do j = 1, m
               if( j == 1 ) then
                   lj = 1
               else
                   lj = l(j-1)+1
               end if
               uj = l(j)

               lx(:,j) = sum( x(:,lj:uj), dim=2 )
            end do
         else
            lx = x
         end if

      end subroutine apply_l_matrix

      subroutine apply_l_vector(x,l,lx)
         real(qp), dimension(:), intent(in)  :: x
         integer, dimension(:), intent(in)   :: l
         real(qp), dimension(:), intent(out) :: lx 

         integer :: i, m, n

         m = size(l)
         n = size(x)
         if( m /= size(lx) ) call system_abort('apply_llt_matrix: l and lx do not conform')

         if( m /= n ) then
             lx(1) = sum( x(1:l(1)) )
             do i = 2, m
                lx(i) = sum( x(l(i-1)+1:l(i)) )
             end do
         else
            lx = x
         end if
      end subroutine apply_l_vector

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Determines the log likelihood of a GP
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine likelihood(this,l,dl_dsigma, dl_ddelta,dl_dtheta,dl_dx,dl_df0,&
      & do_l, do_sigma, do_delta, do_theta, do_x, do_f0)

         type(gp_sparse), intent(inout)                :: this
         real(qp), intent(out), optional               :: l
         real(qp), dimension(:), intent(out), optional :: dl_dsigma
         real(qp), dimension(:), intent(out), optional :: dl_ddelta, dl_df0
         real(qp), dimension(:), intent(out), optional :: dl_dtheta, dl_dx
         logical, intent(in), optional :: do_l, do_sigma, do_delta, do_theta, do_x, do_f0

         type(LA_Matrix) :: LA_q_mm, LA_k_mm

         real(qp) :: det1, det2, det3, tr_dlambda_inverse_k, lknnl, &
         & diff_xijt_dot_x_prime_j1, diff_xijt_dot_x_prime_j2, big_raw_k_mn_ij

         real(qp), dimension(:,:), allocatable :: k_mn_inverse_lambda, k_mn_sq_inverse_lambda, k_mn_l_k_nm, inverse_mm, &
         & k_mn_ll_k_nm, k_mn_ll_k_nm_inverse_mm, d_big_k_mn, d_k_mn, d_k_mm, &
         & inverse_k_mm_k_mn_l_k_nm_inverse_k_mm, inverse_mm_k_mn_inverse_lambda, dk_mm_dx, inverse_k_mm_k_mn_l_k_nm, &
         & inverse_k_mm_k_mn_inverse_k_k_nm_inverse_k_mm, inverse_k_mm_k_mn_inverse_lambda, &
         & inverse_k_mm_k_mn_inverse_k, inverse_mm_k_mn_l_k_nm_inverse_k_mm, theta2, factor_k_mm, a, Q_q_mm, R_q_mm, Q_q_mm_sq_inverse_lambda, &
         & inverse_k_mm_k_mn_sq_inverse_lambda
         real(qp), dimension(:), allocatable   :: inverse_lambda, y_inverse_lambda, y_l_k_nm, &
         & y_l_k_nm_inverse_mm, y_inverse_k, lambda_dtheta, &
         & d_big_k_mn_dx, dk_mn_dx, y_inverse_c, y_inverse_c_k_nm_inverse_k_mm, y_l_k_nm_inverse_k_mm, &
         & dk_mm_inverse_k_mm_k_j, diff_xijt, diag_il_k_nm_inverse_mm_k_mn_il, diag_k_n_inverse_k_mm_inverse_k_mm_k_n, y, &
         & normsq_y_inverse_c, trace_inverse_c
         integer :: i, j, d, j1, j2, lj, uj, info, xj, xj1, xj2, jd, id, num_threads, k, Z_type

         logical :: my_do_l, my_do_sigma, my_do_delta, my_do_theta, my_do_x, my_do_f0

         if(this%gp_teach_memory == GP_TEACH_MEMORY_0) call system_abort('gp_teach_memory = '//this%gp_teach_memory// &
         '. You need to use at least '//GP_TEACH_MEMORY_1//' for likelihood optimization.')

         use_intrinsic_blas = .false.

         my_do_l = present(l) .and. optional_default(present(l), do_l)
         my_do_sigma = present(dl_dsigma) .and. optional_default(present(dl_dsigma), do_sigma)
         my_do_delta = present(dl_ddelta) .and. optional_default(present(dl_ddelta), do_delta)
         my_do_theta = present(dl_dtheta) .and. optional_default(present(dl_dtheta), do_theta)
         my_do_x = present(dl_dx) .and. optional_default(present(dl_dx), do_x)
         my_do_f0 = present(dl_df0) .and. optional_default(present(dl_df0), do_f0)

         if(my_do_l) l = 0.0_qp
         if(my_do_sigma) then
            call check_size('dl_dsigma',dl_dsigma,shape(this%sigma),'likelihood')
            dl_dsigma = 0.0_qp
         endif
         if(my_do_delta) dl_ddelta = 0.0_qp
         if(my_do_theta) dl_dtheta = 0.0_qp
         if(my_do_x) dl_dx = 0.0_qp
         if(my_do_f0) dl_df0 = 0.0_qp

         allocate( inverse_lambda(this%m), y_inverse_lambda(this%m), &
         & k_mn_inverse_lambda(this%sr,this%m), k_mn_l_k_nm(this%sr,this%sr), &
         & inverse_mm(this%sr,this%sr), y_l_k_nm(this%sr), &
         & y_l_k_nm_inverse_mm(this%sr),  k_mn_sq_inverse_lambda(this%sr,this%m), &
         & factor_k_mm(this%sr,this%sr), a(this%sr+this%m,this%sr), y(this%sr+this%m) )

         inverse_lambda = 1.0_qp / this%lambda                    ! O(N)
         y_inverse_lambda = this%y/this%lambda                 ! O(N)

         allocate( theta2(this%d,this%nsp) )
         theta2 = 1.0_qp / this%theta**2

         call matrix_product_vect_asdiagonal_sub(k_mn_inverse_lambda,this%k_mn,inverse_lambda) ! O(NM)
         call matrix_product_vect_asdiagonal_sub(k_mn_sq_inverse_lambda,this%k_mn,sqrt(1.0_qp/this%lambda)) ! O(NM)
!         call matrix_product_sub(k_mn_l_k_nm, k_mn_inverse_lambda, this%k_mn, m2_transpose = .true. ) ! O(NM^2)
         k_mn_l_k_nm = matmul( k_mn_inverse_lambda,transpose(this%k_mn) )
         call initialise(LA_k_mm,this%k_mm)
         call LA_Matrix_Factorise(LA_k_mm,factor_k_mm,info=info)
         if( info /= 0 ) call system_abort('likelihood: LA_k_mm')

         do i = 1, this%sr-1
            do j = i+1, this%sr
               factor_k_mm(j,i) = 0.0_qp
            end do
         end do

         a(1:this%m,:) = transpose(k_mn_sq_inverse_lambda)
         a(this%m+1:,:) = factor_k_mm

         y = 0.0_qp
         y(1:this%m) = this%y*sqrt(1.0_qp/this%lambda)

         call initialise(LA_q_mm,a)
         call LA_Matrix_QR_Solve_Vector(LA_q_mm,y,y_l_k_nm_inverse_mm)

!         call initialise(LA_mm,this%k_mm + k_mn_l_k_nm)
!         call LA_Matrix_Inverse(LA_mm,inverse_mm)
!         det1 = LA_Matrix_LogDet(LA_mm,info=info) / 2.0_qp
         det1 = LA_Matrix_LogDet(LA_q_mm,info=info) !/ 2.0_qp
         if( info /= 0 ) call system_abort('likelihood: LA_q_mm det1')

         y_l_k_nm = matmul( k_mn_inverse_lambda, this%y )         ! O(NM)
         !call Matrix_Solve(LA_mm,y_l_k_nm,y_l_k_nm_inverse_mm,info=info)
         !if( info /= 0 ) call system_abort('likelihood: LA_mm')

         det2 = LA_Matrix_LogDet(LA_k_mm,info=info) / 2.0_qp
         if( info /= 0 ) call system_abort('likelihood: LA_k_mm')

         det3 = 0.5_qp*sum(log(this%lambda))                      ! O(N)

         if( my_do_l ) &
           & l = -(det1-det2+det3) &
           & - 0.5_qp*dot_product(y_inverse_lambda,this%y) &      ! O(N)
           & + 0.5_qp*dot_product(y_l_k_nm,y_l_k_nm_inverse_mm) & ! O(M)
           & - 0.5_qp * this%m * log(2.0_qp * PI)

         if( my_do_sigma .or. my_do_delta .or. my_do_theta .or. my_do_x .or. my_do_f0 ) then
            allocate( y_inverse_k(this%m),y_inverse_c(this%m), inverse_mm_k_mn_inverse_lambda(this%sr,this%m) )
            allocate( Q_q_mm(this%m+this%sr,this%sr), R_q_mm(this%sr,this%sr), Q_q_mm_sq_inverse_lambda(this%sr,this%m) )
            y_inverse_k = matmul(y_l_k_nm_inverse_mm,k_mn_inverse_lambda) ! O(NM)
            y_inverse_c = y_inverse_lambda-y_inverse_k  ! O(N)
            !inverse_mm_k_mn_inverse_lambda = matmul( inverse_mm, k_mn_inverse_lambda ) ! O(NM^2)
            !call Matrix_Solve(LA_mm,k_mn_inverse_lambda,inverse_mm_k_mn_inverse_lambda,info=info)
            !if( info /= 0 ) call system_abort('likelihood: LA_mm')

            call LA_Matrix_GetQR(LA_q_mm,q=Q_q_mm,r=R_q_mm)
            y = 0.0_qp
            y(1:this%m) = sqrt(1.0_qp/this%lambda)
            call matrix_product_vect_asdiagonal_sub(Q_q_mm_sq_inverse_lambda,transpose(Q_q_mm(1:this%m,:)),sqrt(1.0_qp/this%lambda)) ! O(NM)

            inverse_mm_k_mn_inverse_lambda = Q_q_mm_sq_inverse_lambda
            do i = 1, this%m
               do j = this%sr, 2, -1
                  inverse_mm_k_mn_inverse_lambda(j,i) = inverse_mm_k_mn_inverse_lambda(j,i)/R_q_mm(j,j)
                  inverse_mm_k_mn_inverse_lambda(1:j-1,i) = inverse_mm_k_mn_inverse_lambda(1:j-1,i) - inverse_mm_k_mn_inverse_lambda(j,i)*R_q_mm(1:j-1,j)
               end do
               inverse_mm_k_mn_inverse_lambda(1,i) = inverse_mm_k_mn_inverse_lambda(1,i) / R_q_mm(1,1)
            end do
         end if

             !lambda_dtheta = (this%lambda - this%theta(1)**2) / this%theta(2) ! O(N)
         if( my_do_sigma .or. my_do_delta ) then
            allocate(diag_il_k_nm_inverse_mm_k_mn_il(this%m))
            allocate(normsq_y_inverse_c(this%n_target_type),trace_inverse_c(this%n_target_type))
            !allocate(k_mn_ll_k_nm(this%sr,this%sr),k_mn_ll_k_nm_inverse_mm(this%sr,this%sr))
            !call matrix_product_sub(k_mn_ll_k_nm, k_mn_inverse_lambda, k_mn_inverse_lambda, m2_transpose = .true. ) ! O(NM^2)
            !call Matrix_Solve(LA_mm,k_mn_ll_k_nm,k_mn_ll_k_nm_inverse_mm)
            !trace_inverse_c = sum(inverse_lambda) - trace( k_mn_ll_k_nm_inverse_mm ) ! O(N) + O(M^2)
            diag_il_k_nm_inverse_mm_k_mn_il = sum( inverse_mm_k_mn_inverse_lambda*k_mn_inverse_lambda, dim = 1 )
            do i = 1, this%n_target_type
               normsq_y_inverse_c(i) = sum( y_inverse_c**2, mask = (this%target_type==i) )
               trace_inverse_c(i) = sum(inverse_lambda, mask = (this%target_type==i) ) - sum(diag_il_k_nm_inverse_mm_k_mn_il, mask = (this%target_type==i) )
            enddo
               
            !normsq_y_inverse_c1 = normsq(y_inverse_c(:this%mf))      ! O(N)
            !normsq_y_inverse_c2 = normsq(y_inverse_c(this%mf+1:))      ! O(N)

            !trace_inverse_c1 = sum(inverse_lambda(:this%mf)) &
            !& - sum(diag_il_k_nm_inverse_mm_k_mn_il(:this%mf))
            !trace_inverse_c2 = sum(inverse_lambda(this%mf+1:)) &
            !& - sum(diag_il_k_nm_inverse_mm_k_mn_il(this%mf+1:))

            !trace_inverse_c1 = sum(inverse_lambda(:this%mf)) &
            !& - sum(inverse_mm_k_mn_inverse_lambda(:,:this%mf)*k_mn_inverse_lambda(:,:this%mf)) ! O(N) + O(NM)
            !trace_inverse_c2 = sum(inverse_lambda(this%mf+1:)) &
            !& - sum(inverse_mm_k_mn_inverse_lambda(:,this%mf+1:)*k_mn_inverse_lambda(:,this%mf+1:)) ! O(N) + O(NM)
         end if

         if( my_do_delta .or. my_do_theta .or. my_do_x .or. my_do_f0 ) then
            allocate( inverse_k_mm_k_mn_l_k_nm(this%sr,this%sr), &
            & inverse_k_mm_k_mn_inverse_k_k_nm_inverse_k_mm(this%sr,this%sr), &
            & y_inverse_c_k_nm_inverse_k_mm(this%sr), y_l_k_nm_inverse_k_mm(this%sr), &
            & inverse_mm_k_mn_l_k_nm_inverse_k_mm(this%sr,this%sr), inverse_k_mm_k_mn_sq_inverse_lambda(this%sr,this%sr+this%m) )
             
            call Matrix_Solve(LA_k_mm,k_mn_l_k_nm,inverse_k_mm_k_mn_l_k_nm,info=info)
            if( info /= 0 ) call system_abort('likelihood: LA_k_mm')
             
            call Matrix_Solve(LA_k_mm,y_l_k_nm,y_l_k_nm_inverse_k_mm,info=info)
            if( info /= 0 ) call system_abort('likelihood: LA_k_mm')

            call Matrix_Solve(LA_k_mm,k_mn_sq_inverse_lambda,inverse_k_mm_k_mn_sq_inverse_lambda(:,1:this%m),info=info)
            if( info /= 0 ) call system_abort('likelihood: LA_k_mm')
            inverse_k_mm_k_mn_sq_inverse_lambda(:,this%m+1:) = 0.0_qp

            y_inverse_c_k_nm_inverse_k_mm = y_l_k_nm_inverse_k_mm - &
            & matmul( inverse_k_mm_k_mn_l_k_nm, y_l_k_nm_inverse_mm )

            !call Matrix_Solve(LA_mm,transpose(inverse_k_mm_k_mn_l_k_nm), &
            !& inverse_mm_k_mn_l_k_nm_inverse_k_mm, info=info )
            !if( info /= 0 ) call system_abort('likelihood: LA_mm')

            call LA_Matrix_QR_Solve_Matrix(LA_q_mm,transpose(inverse_k_mm_k_mn_sq_inverse_lambda),inverse_mm_k_mn_l_k_nm_inverse_k_mm)

            inverse_k_mm_k_mn_inverse_k_k_nm_inverse_k_mm = matmul( inverse_k_mm_k_mn_l_k_nm, &
            & inverse_mm_k_mn_l_k_nm_inverse_k_mm )
         end if

         if( my_do_sigma ) then
            dl_dsigma = this%sigma * ( normsq_y_inverse_c - trace_inverse_c )
         end if

         if( my_do_delta .or. my_do_f0 .or. my_do_theta .or. my_do_x ) then
            allocate( inverse_k_mm_k_mn_l_k_nm_inverse_k_mm(this%sr,this%sr), &
            & inverse_k_mm_k_mn_inverse_lambda(this%sr,this%m), inverse_k_mm_k_mn_inverse_k(this%sr,this%m) )
            
            call Matrix_Solve(LA_k_mm,k_mn_inverse_lambda,inverse_k_mm_k_mn_inverse_lambda,info=info)
            if( info /= 0 ) call system_abort('likelihood: LA_k_mm')
            
            call Matrix_Solve(LA_k_mm,transpose(inverse_k_mm_k_mn_l_k_nm), &
            & inverse_k_mm_k_mn_l_k_nm_inverse_k_mm, info=info )
            if( info /= 0 ) call system_abort('likelihood: LA_k_mm')

            inverse_k_mm_k_mn_inverse_k = matmul( inverse_k_mm_k_mn_l_k_nm, &   ! O(M^3)
            & inverse_mm_k_mn_inverse_lambda ) ! O(NM^2)
         end if

         if( my_do_delta ) then
            dl_ddelta = 0.0_qp

!$omp parallel private(d_big_k_mn,d_k_mm,d_k_mn,diff_xijt,lambda_dtheta,lknnl,lj,uj,j1,j2,diff_xijt_dot_x_prime_j1,diff_xijt_dot_x_prime_j2,tr_dlambda_inverse_k,id,xj1,xj2)
            allocate( d_big_k_mn(this%sr,this%n), d_k_mm(this%sr,this%sr), d_k_mn(this%sr,this%m), &
            & lambda_dtheta(this%m), diff_xijt(this%d) )

!$omp do private(j,xj,jd,i,k)
             do k = 1, this%nsp
                do j = 1, this%nx
                   xj = this%xf(j)
                   do i = 1, this%sr
                      if( (this%sp(k) == this%xz(xj)) .and. (this%xz(xj) == this%xz_sparse(i)) ) then
                         if(this%gp_teach_memory == GP_TEACH_MEMORY_2) then
                            d_big_k_mn(i,j) = this%big_raw_k_mn(i,j) / this%delta(k)
                         elseif(this%gp_teach_memory == GP_TEACH_MEMORY_1) then
                            d_big_k_mn(i,j) = (this%big_k_mn(i,j) - this%f0(k)) / this%delta(k)
                         endif
                      else
                         d_big_k_mn(i,j) = 0.0_qp
                      end if
                   end do
                end do
                do j = 1, this%nxd
                   xj = this%xdf(j)
                   jd = j + this%nx
                   do i = 1, this%sr
                      if( (this%sp(k) == this%xz(xj)) .and. (this%xz(xj) == this%xz_sparse(i)) ) then
                         d_big_k_mn(i,jd) = this%big_k_mn(i,jd) / this%delta(k)
                      else
                         d_big_k_mn(i,jd) = 0.0_qp
                      end if
                   end do
                end do
                do j = 1, this%sr
                   do i = 1, this%sr
                      if( (this%sp(k) == this%xz_sparse(j) ) .and. (this%xz_sparse(j) == this%xz_sparse(i)) ) then
                         if(i==j) then
                            d_k_mm(i,j) = (this%k_mm(i,j) -this%f0(k)**2 -gp_jitter)/ this%delta(k)
                         else
                            d_k_mm(i,j) = (this%k_mm(i,j) -this%f0(k)**2 )/ this%delta(k)
                         end if

                      else
                         d_k_mm(i,j) = 0.0_qp
                      end if
                   end do
                end do
                call apply_l(d_big_k_mn,this%l,d_k_mn)

                do i = 1, this%mf
                   lknnl = 0.0_qp
                   if(i==1) then
                      lj = 1
                   else
                      lj = this%lf(i-1)+1
                   end if
                   uj = this%lf(i)
                   do j1 = lj, uj
                      xj1 = this%xf(j1)
                      do j2 = lj, uj
                         xj2 = this%xf(j2)
                         if( (this%sp(k) == this%xz(xj1)) .and. (this%xz(xj1) == this%xz(xj2)) ) &
                         & lknnl = lknnl + covSEard( this%delta(k), this%theta(:,k), this%x(:,xj1), this%x(:,xj2) ) / this%delta(k)
                      end do
                   end do
                   lambda_dtheta(i) = lknnl - 2.0_qp * dot_product( d_k_mn(:,i), this%inverse_k_mm_k_mn(:,i) ) + &
                   & dot_product( this%inverse_k_mm_k_mn(:,i), matmul(d_k_mm,this%inverse_k_mm_k_mn(:,i)) )
                end do

                do i = 1, this%mdf
                   id = i + this%mf
                   lknnl = 0.0_qp
                   if(i==1) then
                      lj = 1
                   else
                      lj = this%ldf(i-1)+1
                   end if
                   uj = this%ldf(i)
                   do j1 = lj, uj
                      xj1 = this%xdf(j1)
                      do j2 = lj, uj
                         xj2 = this%xdf(j2)
                         if( (this%sp(k) == this%xz(xj1)) .and. (this%xz(xj1) == this%xz(xj2)) ) then
                            diff_xijt = (this%x(:,xj1) - this%x(:,xj2)) * theta2(:,k)
                            lknnl = lknnl + covSEard( this%delta(k), this%theta(:,k), this%x(:,xj1), this%x(:,xj2) ) * &
                            & ( dot_product( this%x_prime(:,j1) * theta2(:,k), this%x_prime(:,j2) ) - &
                            & dot_product( diff_xijt,this%x_prime(:,j1) ) * dot_product( diff_xijt,this%x_prime(:,j2) ) )/ this%delta(k)
                         end if
                      end do
                   end do
                   lambda_dtheta(id) = lknnl - 2.0_qp * dot_product( d_k_mn(:,id), this%inverse_k_mm_k_mn(:,id) ) + &
                   & dot_product( this%inverse_k_mm_k_mn(:,id), matmul(d_k_mm,this%inverse_k_mm_k_mn(:,id)) )
                end do
                tr_dlambda_inverse_k = 0.0_qp
                do i = 1, this%m
                   tr_dlambda_inverse_k = tr_dlambda_inverse_k + &
                   & lambda_dtheta(i) * dot_product( k_mn_inverse_lambda(:,i),inverse_mm_k_mn_inverse_lambda(:,i) )
                end do
                dl_ddelta(k) = 2.0_qp * dot_product(y_inverse_c_k_nm_inverse_k_mm, matmul(d_k_mn,y_inverse_c) ) - &
                & dot_product(y_inverse_c_k_nm_inverse_k_mm, matmul(d_k_mm,y_inverse_c_k_nm_inverse_k_mm) ) + &
                & dot_product( y_inverse_c*lambda_dtheta, y_inverse_c) - &
                & 2.0_qp * sum( d_k_mn*inverse_k_mm_k_mn_inverse_lambda) + &
                & sum( d_k_mm * inverse_k_mm_k_mn_l_k_nm_inverse_k_mm ) - sum(lambda_dtheta * inverse_lambda) + &
                & 2.0_qp * sum( d_k_mn * inverse_k_mm_k_mn_inverse_k ) - &
                & sum( d_k_mm * inverse_k_mm_k_mn_inverse_k_k_nm_inverse_k_mm ) + tr_dlambda_inverse_k

            end do
!$omp end do
            deallocate( d_big_k_mn, d_k_mm, d_k_mn, lambda_dtheta, diff_xijt )
!$omp end parallel            
         end if

         if( my_do_f0 ) then
            dl_df0 = 0.0_qp

!$omp parallel private(d_big_k_mn,d_k_mm,d_k_mn,diff_xijt,lambda_dtheta,lknnl,lj,uj,j1,j2,diff_xijt_dot_x_prime_j1,diff_xijt_dot_x_prime_j2,tr_dlambda_inverse_k,id,xj1,xj2)
            allocate( d_big_k_mn(this%sr,this%n), d_k_mm(this%sr,this%sr), d_k_mn(this%sr,this%m), &
            & lambda_dtheta(this%m), diff_xijt(this%d) )

!$omp do private(j,xj,jd,i,k)
                do k = 1, this%nsp
                do j = 1, this%nx
                   xj = this%xf(j)
                   do i = 1, this%sr
                      if( (this%sp(k) == this%xz(xj)) .and. (this%xz(xj) == this%xz_sparse(i)) ) then
                         d_big_k_mn(i,j) = this%f0(k)
                      else
                         d_big_k_mn(i,j) = 0.0_qp
                      end if
                   end do
                end do
                do j = 1, this%nxd
                   xj = this%xdf(j)
                   jd = j + this%nx
                   do i = 1, this%sr
                      d_big_k_mn(i,jd) = 0.0_qp
                   end do
                end do
                do j = 1, this%sr
                   do i = 1, this%sr
                      if( (this%sp(k) == this%xz_sparse(j) ) .and. (this%xz_sparse(j) == this%xz_sparse(i)) ) then
                         d_k_mm(i,j) = this%f0(k)
                      else
                         d_k_mm(i,j) = 0.0_qp
                      end if
                   end do
                end do
                call apply_l(d_big_k_mn,this%l,d_k_mn)

                do i = 1, this%mf
                   lknnl = 0.0_qp
                   if(i==1) then
                      lj = 1
                   else
                      lj = this%lf(i-1)+1
                   end if
                   uj = this%lf(i)
                   do j1 = lj, uj
                      xj1 = this%xf(j1)
                      do j2 = lj, uj
                         xj2 = this%xf(j2)
                         if( (this%sp(k) == this%xz(xj1)) .and. (this%xz(xj1) == this%xz(xj2)) ) &
                         & lknnl = lknnl + this%f0(k)
                      end do
                   end do
                   lambda_dtheta(i) = lknnl - 2.0_qp * dot_product( d_k_mn(:,i), this%inverse_k_mm_k_mn(:,i) ) + &
                   & dot_product( this%inverse_k_mm_k_mn(:,i), matmul(d_k_mm,this%inverse_k_mm_k_mn(:,i)) )
                end do

                do i = 1, this%mdf
                   id = i + this%mf
                   lknnl = 0.0_qp
                   if(i==1) then
                      lj = 1
                   else
                      lj = this%ldf(i-1)+1
                   end if
                   uj = this%ldf(i)
                   do j1 = lj, uj
                      xj1 = this%xdf(j1)
                      do j2 = lj, uj
                         xj2 = this%xdf(j2)
                         if( (this%sp(k) == this%xz(xj1)) .and. (this%xz(xj1) == this%xz(xj2)) ) then
                            lknnl = 0.0_qp
                         end if
                      end do
                   end do
                   lambda_dtheta(id) = lknnl - 2.0_qp * dot_product( d_k_mn(:,id), this%inverse_k_mm_k_mn(:,id) ) + &
                   & dot_product( this%inverse_k_mm_k_mn(:,id), matmul(d_k_mm,this%inverse_k_mm_k_mn(:,id)) )
                end do
                tr_dlambda_inverse_k = 0.0_qp
                do i = 1, this%m
                   tr_dlambda_inverse_k = tr_dlambda_inverse_k + &
                   & lambda_dtheta(i) * dot_product( k_mn_inverse_lambda(:,i),inverse_mm_k_mn_inverse_lambda(:,i) )
                end do
                dl_df0(k) = 2.0_qp * dot_product(y_inverse_c_k_nm_inverse_k_mm, matmul(d_k_mn,y_inverse_c) ) - &
                & dot_product(y_inverse_c_k_nm_inverse_k_mm, matmul(d_k_mm,y_inverse_c_k_nm_inverse_k_mm) ) + &
                & dot_product( y_inverse_c*lambda_dtheta, y_inverse_c) - &
                & 2.0_qp * sum( d_k_mn*inverse_k_mm_k_mn_inverse_lambda) + &
                & sum( d_k_mm * inverse_k_mm_k_mn_l_k_nm_inverse_k_mm ) - sum(lambda_dtheta * inverse_lambda) + &
                & 2.0_qp * sum( d_k_mn * inverse_k_mm_k_mn_inverse_k ) - &
                & sum( d_k_mm * inverse_k_mm_k_mn_inverse_k_k_nm_inverse_k_mm ) + tr_dlambda_inverse_k

            end do
!$omp end do
            deallocate( d_big_k_mn, d_k_mm, d_k_mn, lambda_dtheta, diff_xijt )
!$omp end parallel            
         end if
         if( my_do_theta ) then
            num_threads = 1
!$omp parallel
!$omp master
!$ num_threads = omp_get_num_threads()
!$omp end master
!$omp end parallel

!$omp parallel private(d_big_k_mn,d_k_mm,d_k_mn,diff_xijt,lambda_dtheta,lknnl,lj,uj,j1,j2,diff_xijt_dot_x_prime_j1,diff_xijt_dot_x_prime_j2,tr_dlambda_inverse_k,id,xj1,xj2)
            allocate( d_big_k_mn(this%sr,this%n), d_k_mm(this%sr,this%sr), d_k_mn(this%sr,this%m), &
            & lambda_dtheta(this%m), diff_xijt(this%d) )

!$omp do private(d,j,xj,jd,i,k)
             do d = 1, this%d
                do k = 1, this%nsp
                do j = 1, this%nx
                   xj = this%xf(j)
                   do i = 1, this%sr
                      if( (this%sp(k) == this%xz(xj)) .and. (this%xz(xj) == this%xz_sparse(i)) ) then
                         if(this%gp_teach_memory == GP_TEACH_MEMORY_2) then
                            d_big_k_mn(i,j) = this%big_raw_k_mn(i,j) * (this%x(d,xj) - this%x_sparse(d,i))**2
                         elseif(this%gp_teach_memory == GP_TEACH_MEMORY_1) then
                            d_big_k_mn(i,j) = (this%big_k_mn(i,j) - this%f0(k)) * (this%x(d,xj) - this%x_sparse(d,i))**2
                         endif
                      else
                         d_big_k_mn(i,j) = 0.0_qp
                      end if
                   end do
                end do
                do j = 1, this%nxd
                   xj = this%xdf(j)
                   jd = j + this%nx
                   do i = 1, this%sr
                      if( (this%sp(k) == this%xz(xj)) .and. (this%xz(xj) == this%xz_sparse(i)) ) then
                         if(this%gp_teach_memory == GP_TEACH_MEMORY_2) then
                            d_big_k_mn(i,jd) = this%big_k_mn(i,jd) * (this%x(d,xj) - this%x_sparse(d,i))**2 - &
                            2.0_qp * this%big_raw_k_mn(i,jd) * ( this%x_sparse(d,i) - this%x(d,xj) ) * this%x_prime(d,j) 
                         elseif(this%gp_teach_memory == GP_TEACH_MEMORY_1) then
                            big_raw_k_mn_ij = this%big_k_mn(i,jd) / &
                            dot_product(( this%x_sparse(:,i) - this%x(:,xj) )*theta2(:,k),this%x_prime(:,j))

                            d_big_k_mn(i,jd) = this%big_k_mn(i,jd) * (this%x(d,xj) - this%x_sparse(d,i))**2 - &
                            2.0_qp * big_raw_k_mn_ij * ( this%x_sparse(d,i) - this%x(d,xj) ) * this%x_prime(d,j) 
                         endif
                      else
                         d_big_k_mn(i,jd) = 0.0_qp
                      end if
                   end do
                end do
                do j = 1, this%sr
                   do i = 1, this%sr
                      if( (this%sp(k) == this%xz_sparse(j) ) .and. (this%xz_sparse(j) == this%xz_sparse(i)) ) then
                         if(i==j) then
                            d_k_mm(i,j) = (this%k_mm(i,j)-this%f0(k)**2-gp_jitter) * ( this%x_sparse(d,j) - this%x_sparse(d,i))**2
                         else
                            d_k_mm(i,j) = (this%k_mm(i,j)-this%f0(k)**2) * ( this%x_sparse(d,j) - this%x_sparse(d,i))**2
                         endif
                      else
                         d_k_mm(i,j) = 0.0_qp
                      endif
                   enddo
                enddo
                call apply_l(d_big_k_mn,this%l,d_k_mn)

                do i = 1, this%mf
                   lknnl = 0.0_qp
                   if(i==1) then
                      lj = 1
                   else
                      lj = this%lf(i-1)+1
                   endif
                   uj = this%lf(i)
                   do j1 = lj, uj
                      xj1 = this%xf(j1)
                      do j2 = lj, uj
                         xj2 = this%xf(j2)
                         if( (this%sp(k) == this%xz(xj1)) .and. (this%xz(xj1) == this%xz(xj2)) ) &
                         & lknnl = lknnl + covSEard( this%delta(k), this%theta(:,k), this%x(:,xj1), this%x(:,xj2) ) * &
                         & ( this%x(d,xj1) - this%x(d,xj2))**2
                      enddo
                   enddo
                   lambda_dtheta(i) = lknnl - 2.0_qp * dot_product( d_k_mn(:,i), this%inverse_k_mm_k_mn(:,i) ) + &
                   & dot_product( this%inverse_k_mm_k_mn(:,i), matmul(d_k_mm,this%inverse_k_mm_k_mn(:,i)) )
                enddo

                do i = 1, this%mdf
                   id = i + this%mf
                   lknnl = 0.0_qp
                   if(i==1) then
                      lj = 1
                   else
                      lj = this%ldf(i-1)+1
                   endif
                   uj = this%ldf(i)
                   do j1 = lj, uj
                      xj1 = this%xdf(j1)
                      do j2 = lj, uj
                         xj2 = this%xdf(j2)
                         if( (this%sp(k) == this%xz(xj1)) .and. (this%xz(xj1) == this%xz(xj2)) ) then
                            diff_xijt = (this%x(:,xj1) - this%x(:,xj2)) * theta2(:,k)
                            diff_xijt_dot_x_prime_j1 = dot_product( diff_xijt,this%x_prime(:,j1) )
                            diff_xijt_dot_x_prime_j2 = dot_product( diff_xijt,this%x_prime(:,j2) )
                            lknnl = lknnl + covSEard( this%delta(k), this%theta(:,k), this%x(:,xj1), this%x(:,xj2) ) * &
                            & ( ( dot_product( this%x_prime(:,j1) * theta2(:,k), this%x_prime(:,j2) ) - &
                            & diff_xijt_dot_x_prime_j1*diff_xijt_dot_x_prime_j2 ) * ( this%x(d,xj1) - this%x(d,xj2))**2 + 2*(- &
                            & this%x_prime(d,j1)*this%x_prime(d,j2) + &
                            & diff_xijt_dot_x_prime_j1 * (this%x(d,xj1) - this%x(d,xj2)) * this%x_prime(d,j2) + &
                            & diff_xijt_dot_x_prime_j2 * (this%x(d,xj1) - this%x(d,xj2)) * this%x_prime(d,j1)) )
                         endif
                      enddo
                   enddo
                   lambda_dtheta(id) = lknnl - 2.0_qp * dot_product( d_k_mn(:,id), this%inverse_k_mm_k_mn(:,id) ) + &
                   & dot_product( this%inverse_k_mm_k_mn(:,id), matmul(d_k_mm,this%inverse_k_mm_k_mn(:,id)) )
                enddo
                tr_dlambda_inverse_k = 0.0_qp
                do i = 1, this%m
                   tr_dlambda_inverse_k = tr_dlambda_inverse_k + &
                   & lambda_dtheta(i) * dot_product( k_mn_inverse_lambda(:,i),inverse_mm_k_mn_inverse_lambda(:,i) )
                enddo
                dl_dtheta(d+(k-1)*this%d) = 2.0_qp * dot_product(y_inverse_c_k_nm_inverse_k_mm, matmul(d_k_mn,y_inverse_c) ) - &
                & dot_product(y_inverse_c_k_nm_inverse_k_mm, matmul(d_k_mm,y_inverse_c_k_nm_inverse_k_mm) ) + &
                & dot_product( y_inverse_c*lambda_dtheta, y_inverse_c) - &
                & 2.0_qp * sum( d_k_mn*inverse_k_mm_k_mn_inverse_lambda) + &
                & sum( d_k_mm * inverse_k_mm_k_mn_l_k_nm_inverse_k_mm ) - sum(lambda_dtheta * inverse_lambda) + &
                & 2.0_qp * sum( d_k_mn * inverse_k_mm_k_mn_inverse_k ) - &
                & sum( d_k_mm * inverse_k_mm_k_mn_inverse_k_k_nm_inverse_k_mm ) + tr_dlambda_inverse_k
                dl_dtheta(d+(k-1)*this%d) = 0.5_qp * dl_dtheta(d+(k-1)*this%d) / this%theta(d,k)**3
             enddo
            enddo
!$omp end do
            deallocate( d_big_k_mn, d_k_mm, d_k_mn, lambda_dtheta, diff_xijt )
!$omp end parallel            
         endif

         if( my_do_x ) then

!$omp parallel private(d_big_k_mn_dx,dk_mn_dx,dk_mm_dx,dk_mm_inverse_k_mm_k_j,lambda_dtheta,tr_dlambda_inverse_k)
            allocate( d_big_k_mn_dx(this%n), dk_mn_dx(this%m), &
            & dk_mm_dx(this%sr,this%sr), dk_mm_inverse_k_mm_k_j(this%sr), lambda_dtheta(this%m))
            
!$omp do private(i,d,j,xj,jd,Z_type,k) 
             do i = 1, this%sr
                
                do k = 1, this%nsp; if( this%sp(k) == this%xz_sparse(i) ) Z_type = k; enddo
                
                do d = 1, this%d
                   do j = 1, this%nx
                      xj = this%xf(j)
                      if( this%xz(xj) == this%xz_sparse(i) ) then
                         if(this%gp_teach_memory == GP_TEACH_MEMORY_2) then
                            d_big_k_mn_dx(j) = this%big_raw_k_mn(i,j) * &
                            ( this%x(d,xj) - this%x_sparse(d,i))
                         elseif(this%gp_teach_memory == GP_TEACH_MEMORY_1) then 
                            d_big_k_mn_dx(j) = (this%big_k_mn(i,j) - this%f0(Z_type)**2) * &
                            ( this%x(d,xj) - this%x_sparse(d,i))
                         endif
                      else
                         d_big_k_mn_dx(j) = 0.0_qp
                      endif
                   enddo

                   do j = 1, this%nxd
                      jd = j + this%nx
                      xj = this%xdf(j)
                      if( this%xz(xj) == this%xz_sparse(i) ) then
                         if(this%gp_teach_memory == GP_TEACH_MEMORY_2) then
                            d_big_k_mn_dx(jd) = this%big_k_mn(i,jd) * &
                            ( this%x(d,xj) - this%x_sparse(d,i) ) + &
                            this%big_raw_k_mn(i,jd) * this%x_prime(d,j)
                         elseif(this%gp_teach_memory == GP_TEACH_MEMORY_1) then
                            big_raw_k_mn_ij = this%big_k_mn(i,jd) / &
                            dot_product(( this%x_sparse(:,i) - this%x(:,xj) )*theta2(:,Z_type),this%x_prime(:,j))

                            d_big_k_mn_dx(jd) = this%big_k_mn(i,jd) * &
                            ( this%x(d,xj) - this%x_sparse(d,i) ) + &
                            big_raw_k_mn_ij * this%x_prime(d,j)
                         endif
                      else
                         d_big_k_mn_dx(jd) = 0.0_qp
                      endif
                   enddo
                   call apply_l(d_big_k_mn_dx,this%l,dk_mn_dx)
                   dk_mm_dx = 0.0_qp
                   do j = 1, this%sr
                      if( this%xz_sparse(j) == this%xz_sparse(i) ) then
                         if(i==j) then
                            dk_mm_dx(j,i) = (this%k_mm(j,i)-this%f0(Z_type)**2-gp_jitter) * ( this%x_sparse(d,j) - this%x_sparse(d,i) )
                         else
                            dk_mm_dx(j,i) = (this%k_mm(j,i)-this%f0(Z_type)**2) * ( this%x_sparse(d,j) - this%x_sparse(d,i) )
                         endif
                         dk_mm_dx(i,j) = dk_mm_dx(j,i)
                      else
                         dk_mm_dx(j,i) = 0.0_qp
                         dk_mm_dx(i,j) = 0.0_qp
                      endif
                   enddo
                   do j = 1, this%m
                      dk_mm_inverse_k_mm_k_j = dk_mm_dx(:,i) * this%inverse_k_mm_k_mn(i,j)
                      dk_mm_inverse_k_mm_k_j(i) = dot_product( dk_mm_dx(:,i), this%inverse_k_mm_k_mn(:,j) )
                      lambda_dtheta(j) = -2.0_qp * dk_mn_dx(j) * this%inverse_k_mm_k_mn(i,j) + &
                      & dot_product( this%inverse_k_mm_k_mn(:,j), dk_mm_inverse_k_mm_k_j )
                   enddo
                   tr_dlambda_inverse_k = 0.0_qp
                   do j = 1, this%m
                      tr_dlambda_inverse_k = tr_dlambda_inverse_k + &
                      & lambda_dtheta(j) * dot_product( k_mn_inverse_lambda(:,j),inverse_mm_k_mn_inverse_lambda(:,j) )
                   enddo

                   dl_dx((i-1)*this%d+d) = &
                   & 2.0_qp * dot_product(y_inverse_c,dk_mn_dx) * y_inverse_c_k_nm_inverse_k_mm(i) &
                   & - dot_product(y_inverse_c_k_nm_inverse_k_mm, matmul(dk_mm_dx,y_inverse_c_k_nm_inverse_k_mm) ) &
                   & + dot_product( y_inverse_c*lambda_dtheta, y_inverse_c) &
                   & - 2.0_qp * dot_product(dk_mn_dx,(inverse_k_mm_k_mn_inverse_lambda(i,:)-inverse_k_mm_k_mn_inverse_k(i,:))) & !#1
                   & + sum( dk_mm_dx * (inverse_k_mm_k_mn_l_k_nm_inverse_k_mm-inverse_k_mm_k_mn_inverse_k_k_nm_inverse_k_mm) ) &
                   & - sum(lambda_dtheta * inverse_lambda) &  !#2
                   & + tr_dlambda_inverse_k
                   dl_dx((i-1)*this%d+d) = 0.5_qp * dl_dx((i-1)*this%d+d) / this%theta(d,Z_type)**2 
                enddo
             enddo
!$omp end do
            deallocate( d_big_k_mn_dx, dk_mn_dx, dk_mm_dx, dk_mm_inverse_k_mm_k_j, lambda_dtheta)
!$omp end parallel
         endif
                   !& 2.0_qp * dot_product(inverse_k_mm(:,i),matmul(k_mn_inverse_lambda,dk_mn_dx)) + & !#1
                   !& 2.0_qp * sum( (inverse_k_mm(:,i) .outer. dk_mn_dx)*k_mn_inverse_lambda) + & !#1
                   !& 2.0_qp * sum( dk_mn_dx * inverse_k_mm_k_mn_inverse_k(i,:) ) - & !#2

         call finalise(LA_q_mm)
         call finalise(LA_k_mm)

         if(allocated(inverse_lambda)) deallocate(inverse_lambda)
         if(allocated(y_inverse_lambda)) deallocate(y_inverse_lambda)
         if(allocated(k_mn_inverse_lambda)) deallocate(k_mn_inverse_lambda)
         if(allocated(k_mn_l_k_nm)) deallocate(k_mn_l_k_nm)
         if(allocated(inverse_mm)) deallocate(inverse_mm)
         if(allocated(y_l_k_nm)) deallocate(y_l_k_nm)
         if(allocated(y_l_k_nm_inverse_mm)) deallocate(y_l_k_nm_inverse_mm)

         if(allocated(y_inverse_k)) deallocate(y_inverse_k)
         if(allocated(y_inverse_c)) deallocate(y_inverse_c)

         if(allocated(k_mn_ll_k_nm)) deallocate(k_mn_ll_k_nm)
         if(allocated(k_mn_ll_k_nm_inverse_mm)) deallocate(k_mn_ll_k_nm_inverse_mm)

         if(allocated(lambda_dtheta)) deallocate(lambda_dtheta)
         if(allocated(inverse_mm_k_mn_inverse_lambda)) deallocate(inverse_mm_k_mn_inverse_lambda)
         if(allocated(y_inverse_c_k_nm_inverse_k_mm)) deallocate(y_inverse_c_k_nm_inverse_k_mm)
         if(allocated(inverse_k_mm_k_mn_l_k_nm)) deallocate(inverse_k_mm_k_mn_l_k_nm)
         if(allocated(inverse_k_mm_k_mn_l_k_nm_inverse_k_mm)) deallocate(inverse_k_mm_k_mn_l_k_nm_inverse_k_mm)
         if(allocated(inverse_k_mm_k_mn_inverse_k)) deallocate(inverse_k_mm_k_mn_inverse_k)
         if(allocated(inverse_k_mm_k_mn_inverse_k_k_nm_inverse_k_mm)) &
         & deallocate(inverse_k_mm_k_mn_inverse_k_k_nm_inverse_k_mm)

         if(allocated(inverse_mm_k_mn_l_k_nm_inverse_k_mm)) deallocate(inverse_mm_k_mn_l_k_nm_inverse_k_mm)
         if(allocated(y_l_k_nm_inverse_k_mm)) deallocate(y_l_k_nm_inverse_k_mm)

         if(allocated(d_big_k_mn)) deallocate(d_big_k_mn)
         if(allocated(d_k_mm)) deallocate(d_k_mm)
         if(allocated(d_k_mn)) deallocate(d_k_mn)
         if(allocated(inverse_k_mm_k_mn_inverse_lambda)) deallocate(inverse_k_mm_k_mn_inverse_lambda)
         if(allocated(diff_xijt)) deallocate(diff_xijt)
         if(allocated(theta2)) deallocate(theta2)

         if(allocated(d_big_k_mn_dx)) deallocate(d_big_k_mn_dx)
         if(allocated(dk_mn_dx)) deallocate(dk_mn_dx)
         if(allocated(dk_mm_dx)) deallocate(dk_mm_dx)
         if(allocated(dk_mm_inverse_k_mm_k_j)) deallocate(dk_mm_inverse_k_mm_k_j )
         if(allocated(diag_il_k_nm_inverse_mm_k_mn_il)) deallocate(diag_il_k_nm_inverse_mm_k_mn_il )
         if(allocated(diag_k_n_inverse_k_mm_inverse_k_mm_k_n)) deallocate(diag_k_n_inverse_k_mm_inverse_k_mm_k_n )
         if(allocated(k_mn_sq_inverse_lambda)) deallocate(k_mn_sq_inverse_lambda )
         if(allocated(factor_k_mm)) deallocate(factor_k_mm )
         if(allocated(a)) deallocate(a)
         if(allocated(y)) deallocate(y)
         if(allocated(Q_q_mm)) deallocate(Q_q_mm)
         if(allocated(R_q_mm)) deallocate(R_q_mm)
         if(allocated(Q_q_mm_sq_inverse_lambda)) deallocate(Q_q_mm_sq_inverse_lambda)
         if(allocated(inverse_k_mm_k_mn_sq_inverse_lambda)) deallocate(inverse_k_mm_k_mn_sq_inverse_lambda)
         if(allocated(normsq_y_inverse_c)) deallocate(normsq_y_inverse_c)
         if(allocated(trace_inverse_c)) deallocate(trace_inverse_c)

      end subroutine likelihood
      
      subroutine likelihood_simple(this,l,dl_dsigma, dl_ddelta,dl_dtheta,dl_df0, &
      & do_l, do_sigma, do_delta, do_theta, do_f0)

         type(gp), intent(inout)                       :: this
         real(qp), intent(out), optional               :: l
         real(qp), intent(out), optional               :: dl_dsigma, dl_ddelta, dl_df0
         real(qp), dimension(:), intent(out), optional :: dl_dtheta
         logical, intent(in), optional :: do_l, do_sigma, do_delta, do_theta, do_f0

         logical :: my_do_l, my_do_sigma, my_do_delta, my_do_theta, my_do_f0
         real(qp), dimension(:,:), allocatable :: c_inverse, outer_alpha_minus_c_inverse, k, dk
#ifdef HAVE_QR
         real(qp), dimension(:,:), allocatable :: one
#endif
         integer :: i, j, d

         my_do_l = present(l) .and. optional_default(present(l), do_l)
         my_do_sigma = present(dl_dsigma) .and. optional_default(present(dl_dsigma), do_sigma)
         my_do_delta = present(dl_ddelta) .and. optional_default(present(dl_ddelta), do_delta)
         my_do_theta = present(dl_dtheta) .and. optional_default(present(dl_dtheta), do_theta)
         my_do_f0 = present(dl_df0) .and. optional_default(present(dl_df0), do_f0)

         if(my_do_l) l = 0.0_qp

         if( my_do_l ) &
         l = - 0.5_qp * LA_Matrix_LogDet(this%LA_C_nn)  &
             - 0.5_qp * dot_product( this%y, this%alpha )  &
             - 0.5_qp * this%n * log(2.0_qp * PI)

         if( my_do_sigma .or. my_do_delta .or. my_do_theta .or. my_do_f0 ) then
#ifdef HAVE_QR
            allocate(c_inverse(this%n,this%n),one(this%n,this%n))
            one = 0.0_dp
            do i = 1, this%n        
               one(i,i) = 1.0_dp        
            enddo
            call LA_Matrix_QR_Solve_Matrix(this%LA_C_nn, one, c_inverse)     
            deallocate(one)
#else
            allocate(c_inverse(this%n,this%n))
            call LA_Matrix_Inverse(this%LA_C_nn, c_inverse)
#endif
         endif

         if(my_do_delta .or. my_do_theta) then
            allocate(outer_alpha_minus_c_inverse(this%n,this%n),k(this%n,this%n))
            k = this%LA_C_nn%matrix
            do j = 1, this%n
               do i = 1, this%n
                  outer_alpha_minus_c_inverse(i,j) = this%alpha(i) * this%alpha(j) - c_inverse(i,j)
               enddo
               k(j,j) = k(j,j) - this%sigma**2
            enddo
         endif

         if( my_do_sigma ) &
         dl_dsigma = this%sigma * ( dot_product(this%alpha,this%alpha) - trace(c_inverse) )

         if( my_do_delta ) &
         dl_ddelta = sum( outer_alpha_minus_c_inverse * k ) / this%delta(1)

         if( my_do_theta ) then
            dl_dtheta = 0.0_qp
            allocate(dk(this%n,this%n))
            do d = 1, this%d
               do j = 1, this%n
                  do i = 1, this%n
                     dk(i,j) = (this%x(d,i)-this%x(d,j))**2
                  enddo
               enddo
               dl_dtheta(d) = 0.5_dp * sum( outer_alpha_minus_c_inverse * dk * k) / this%theta(d,1)**3
            enddo
         endif


         if( my_do_f0 ) dl_df0 = 0.0_qp

         if(allocated(c_inverse)) deallocate(c_inverse)
         if(allocated(outer_alpha_minus_c_inverse)) deallocate(outer_alpha_minus_c_inverse)
         if(allocated(k)) deallocate(k)
         if(allocated(dk)) deallocate(dk)

      endsubroutine likelihood_simple

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Maximises the log likelihood of a GP in the space of the hyperparameters
      !% using nested sampling.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      function minimise_gp_ns(this,theta_min,theta_max,N_live,max_iter,tighten)

         type(gp_sparse), intent(inout) :: this                               !% GP
         real(dp), dimension(:), intent(inout) :: theta_min,theta_max  !% range of hyperparameters, dimension(p)
         integer, optional, intent(in)         :: N_live               !% number of living points
         integer, optional, intent(in)         :: max_iter             !% maximum number of iterations
         logical, optional, intent(in)         :: tighten              !% tighten the range of hyperparameters

         logical :: minimise_gp_ns

         real(dp), dimension(:), allocatable :: x

         if( .not. this%initialised ) call system_abort('gp_minimise_ns: initialise first')
         !if( ( size(theta_min) /= this%p ) .or. ( size(theta_max) /= this%p ) ) &
         !& call system_abort('gp_minimise_ns: size mismatch')

         !allocate( x(3+this%d+this%sr*this%d) )
         !x(1:2) = real(this%sigma,kind=dp)
         !x(3) = real(this%delta,kind=dp)
         !x(4:3+this%d) = real(this%theta,kind=dp)
         !x(4+this%d:) = real(reshape(this%x_sparse,(/this%sr*this%d/)),kind=dp)

         allocate( x(this%d*this%nsp) )
         x = real(reshape(this%theta,(/this%d*this%nsp/)),kind=dp)

         !minimise_gp_ns = ns(x,l,theta_min,theta_max,N_live,max_iter,tighten)

         deallocate( x )

         contains

            function l(x_in)

               real(dp), dimension(:), intent(in) :: x_in
               real(dp)                           :: l
               real(qp), dimension(:), allocatable :: x
               real(qp) :: l_qp

               allocate(x(size(x_in)))
               x = real(x_in,kind=qp)

               call gp_update(this, theta_in=reshape(x,(/this%d,this%nsp/)) )
               call likelihood(this,l=l_qp)
               l = real(l_qp,kind=dp)
               deallocate(x)

            end function l

      end function minimise_gp_ns

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Tests the derivative of the log likelihood with test_gradient.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      function test_gp_gradient(this,sigma,delta,theta,sparx,f0,theta_fac)

         type(gp_sparse), intent(inout), target :: this                 !% GP
         logical, intent(in), optional :: sigma, delta, theta, sparx, f0, theta_fac
         logical                 :: test_gp_gradient
         real(dp), dimension(:), allocatable :: xx_dp
         real(qp), dimension(:), allocatable :: xx
         integer :: n, li, ui

         type(gp_minimise) :: am
         character, dimension(:), allocatable :: am_data
         character, dimension(1) :: am_mold
         integer :: am_data_size
         real(qp), dimension(:,:), allocatable, target :: theta_0

         am_data_size = size(transfer(am,am_mold))
         allocate(am_data(am_data_size))

         am%minim_gp => this

         am%do_sigma = optional_default(.true.,sigma)
         am%do_delta = optional_default(.true.,delta)
         am%do_theta = optional_default(.true.,theta)
         am%do_sparx = optional_default(.true.,sparx)
         am%do_f0 = optional_default(.true.,f0)
         am%do_theta_fac = optional_default(.false.,theta_fac)

         if(am%do_theta .and. am%do_theta_fac) then
            call print_warning('test_gp_gradient called with both theta and theta_fac active. Switching theta_fac to false.')
            am%do_theta_fac = .false.
         endif

         if(am%do_theta_fac) then
            allocate(theta_0(this%d,this%nsp))
            theta_0 = this%theta
            am%theta_0 => theta_0
         endif

         test_gp_gradient = .false.
         if( (.not.am%do_sigma) .and. (.not.am%do_delta) .and. (.not.am%do_theta) .and. (.not.am%do_sparx) .and. (.not.am%do_f0) .and. (.not.am%do_theta_fac)  ) return

         n = 0
         li = 0
         ui = 0
         if( am%do_sigma ) n = n + this%n_target_type
         if( am%do_delta ) n = n + this%nsp
         if( am%do_theta ) n = n + this%d*this%nsp
         if( am%do_sparx ) n = n + this%sr*this%d
         if( am%do_f0 ) n = n + this%nsp
         if( am%do_theta_fac ) n = n + this%nsp

         allocate( xx(n), xx_dp(n) )

         if(am%do_sigma) then
            li = li + 1
            ui = li + this%n_target_type - 1
            am%li_sigma = li
            am%ui_sigma = ui
            xx(li:ui) = this%sigma
         endif
         if(am%do_delta) then
            li = ui + 1
            ui = li + this%nsp - 1
            am%li_delta = li
            am%ui_delta = ui
            xx(li:ui) = this%delta
         endif
         if(am%do_theta) then
            li = ui + 1
            ui = li + this%d*this%nsp - 1
            am%li_theta = li
            am%ui_theta = ui
            xx(li:ui) = reshape(this%theta,(/this%d*this%nsp/))
         endif
         if(am%do_sparx) then
            li = ui + 1
            ui = li + this%sr*this%d - 1
            am%li_sparx = li
            am%ui_sparx = ui
            xx(li:ui) = reshape(this%x_sparse,(/this%sr*this%d/))
         endif
         if(am%do_f0) then
            li = ui + 1
            ui = li + this%nsp - 1
            am%li_f0 = li
            am%ui_f0 = ui
            xx(li:ui) = this%f0
         endif
         if(am%do_theta_fac) then
            li = ui + 1
            ui = li + this%nsp - 1
            am%li_theta_fac = li
            am%ui_theta_fac = ui
            xx(li:ui) = 1.0_qp
         endif

         xx_dp = real(xx,kind=dp)
         am_data = transfer(am, am_data)

         test_gp_gradient = test_gradient(xx_dp, likelihood_function, likelihood_gradient, data=am_data)

         deallocate(xx,xx_dp,am_data)
         am%minim_gp => null()
         am%theta_0 => null()
         if(allocated(theta_0)) deallocate(theta_0)

      end function test_gp_gradient

      function test_gp_simple_gradient(this,sigma,delta,theta,f0,theta_fac)

         type(gp), intent(inout), target :: this                 !% GP
         logical, intent(in), optional :: sigma, delta, theta, f0, theta_fac
         logical                 :: test_gp_simple_gradient
         real(dp), dimension(:), allocatable :: xx_dp
         real(qp), dimension(:), allocatable :: xx
         integer :: n, li, ui

         type(gp_simple_minimise) :: am
         character, dimension(:), allocatable :: am_data
         character, dimension(1) :: am_mold
         integer :: am_data_size
         real(qp), dimension(:,:), allocatable, target :: theta_0

         am_data_size = size(transfer(am,am_mold))
         allocate(am_data(am_data_size))

         am%minim_gp => this

         am%do_sigma = optional_default(.true.,sigma)
         am%do_delta = optional_default(.true.,delta)
         am%do_theta = optional_default(.true.,theta)
         am%do_f0 = optional_default(.true.,f0)
         am%do_theta_fac = optional_default(.false.,theta_fac)

         if(am%do_theta .and. am%do_theta_fac) then
            call print_warning('test_gp_gradient called with both theta and theta_fac active. Switching theta_fac to false.')
            am%do_theta_fac = .false.
         endif

         if(am%do_theta_fac) then
            allocate(theta_0(this%d,this%nsp))
            theta_0 = this%theta
            am%theta_0 => theta_0
         endif

         test_gp_simple_gradient = .false.
         if( (.not.am%do_sigma) .and. (.not.am%do_delta) .and. (.not.am%do_theta) .and. (.not.am%do_f0) .and. (.not.am%do_theta_fac)  ) return

         n = 0
         li = 0
         ui = 0
         if( am%do_sigma ) n = n + 1
         if( am%do_delta ) n = n + this%nsp
         if( am%do_theta ) n = n + this%d*this%nsp
         if( am%do_f0 ) n = n + this%nsp
         if( am%do_theta_fac ) n = n + this%nsp

         allocate( xx(n), xx_dp(n) )

         if(am%do_sigma) then
            li = li + 1
            ui = li
            am%li_sigma = li
            am%ui_sigma = ui
            xx(li:ui) = this%sigma
         endif
         if(am%do_delta) then
            li = ui + 1
            ui = li + this%nsp - 1
            am%li_delta = li
            am%ui_delta = ui
            xx(li:ui) = this%delta
         endif
         if(am%do_theta) then
            li = ui + 1
            ui = li + this%d*this%nsp - 1
            am%li_theta = li
            am%ui_theta = ui
            xx(li:ui) = reshape(this%theta,(/this%d*this%nsp/))
         endif
         if(am%do_f0) then
            li = ui + 1
            ui = li + this%nsp - 1
            am%li_f0 = li
            am%ui_f0 = ui
            xx(li:ui) = this%f0
         endif
         if(am%do_theta_fac) then
            li = ui + 1
            ui = li + this%nsp - 1
            am%li_theta_fac = li
            am%ui_theta_fac = ui
            xx(li:ui) = 1.0_qp
         endif

         xx_dp = real(xx,kind=dp)
         am_data = transfer(am, am_data)

         test_gp_simple_gradient = test_gradient(xx_dp, likelihood_simple_function, likelihood_simple_gradient, data=am_data)

         deallocate(xx,xx_dp,am_data)
         am%minim_gp => null()
         am%theta_0 => null()
         if(allocated(theta_0)) deallocate(theta_0)

      end function test_gp_simple_gradient

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Maximises the log likelihood in the space of hyperparameters
      !% using minim (conjugate gradients).
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      function minimise_gp_gradient(this,convergence_tol,max_steps,always_do_test_gradient,sigma,delta,theta,sparx,f0,theta_fac)

         type(gp_sparse), intent(inout), target :: this                       !% GP
         integer, intent(in), optional :: max_steps            !% maximum number of steps, default: 100
         real(dp), intent(in), optional :: convergence_tol     !% convergence tolerance, default: 0.1
         logical, intent(in), optional :: always_do_test_gradient
         logical, intent(in), optional :: sigma, delta, theta, sparx, f0, theta_fac
         integer                 :: minimise_gp_gradient
         real(dp), dimension(:), allocatable :: xx_dp
         real(qp), dimension(:), allocatable :: xx
         integer :: n, li, ui

         integer :: do_max_steps
         real(dp) :: do_convergence_tol
         type(gp_minimise) :: am
         character, dimension(:), allocatable :: am_data
         character, dimension(1) :: am_mold
         integer :: am_data_size
         real(qp), dimension(:,:), allocatable, target :: theta_0

         do_max_steps = optional_default(100,max_steps)
         do_convergence_tol = optional_default(0.0001_dp,convergence_tol)

         am_data_size = size(transfer(am,am_mold))
         allocate(am_data(am_data_size))

         am%minim_gp => this

         am%do_sigma = optional_default(.true.,sigma)
         am%do_delta = optional_default(.true.,delta)
         am%do_theta = optional_default(.true.,theta)
         am%do_sparx = optional_default(.true.,sparx)
         am%do_f0 = optional_default(.true.,f0)
         am%do_theta_fac = optional_default(.false.,theta_fac)

         if(am%do_theta .and. am%do_theta_fac) then
            call print_warning('minimise_gp_gradient called with both theta and theta_fac active. Switching theta_fac to false.')
            am%do_theta_fac = .false.
         endif

         if(am%do_theta_fac) then
            allocate(theta_0(this%d,this%nsp))
            theta_0 = this%theta
            am%theta_0 => theta_0
         endif

         minimise_gp_gradient = 0
         if( (.not.am%do_sigma) .and. (.not.am%do_delta) .and. (.not.am%do_theta) .and. (.not.am%do_sparx) .and. (.not.am%do_f0) .and. (.not.am%do_theta_fac) ) return

         n = 0
         li = 0
         ui = 0
         if( am%do_sigma ) n = n + this%n_target_type
         if( am%do_delta ) n = n + this%nsp
         if( am%do_theta ) n = n + this%d*this%nsp
         if( am%do_sparx ) n = n + this%sr*this%d
         if( am%do_f0 ) n = n + this%nsp
         if( am%do_theta_fac ) n = n + this%nsp

         allocate( xx(n), xx_dp(n) )

         if(am%do_sigma) then
            li = li + 1
            ui = li + this%n_target_type - 1
            am%li_sigma = li
            am%ui_sigma = ui
            xx(li:ui) = this%sigma
         endif
         if(am%do_delta) then
            li = ui + 1
            ui = li + this%nsp - 1
            am%li_delta = li
            am%ui_delta = ui
            xx(li:ui) = this%delta
         endif
         if(am%do_theta) then
            li = ui + 1
            ui = li + this%d*this%nsp - 1
            am%li_theta = li
            am%ui_theta = ui
            xx(li:ui) = reshape(this%theta,(/this%d*this%nsp/))
         endif
         if(am%do_sparx) then
            li = ui + 1
            ui = li + this%sr*this%d - 1
            am%li_sparx = li
            am%ui_sparx = ui
            xx(li:ui) = reshape(this%x_sparse,(/this%sr*this%d/))
         endif
         if(am%do_f0) then
            li = ui + 1
            ui = li + this%nsp - 1
            am%li_f0 = li
            am%ui_f0 = ui
            xx(li:ui) = this%f0
         endif
         if(am%do_theta_fac) then
            li = ui + 1
            ui = li + this%nsp - 1
            am%li_theta_fac = li
            am%ui_theta_fac = ui
            xx(li:ui) = 1.0_qp
         endif

         xx_dp = real(xx,kind=dp)
         am_data = transfer(am, am_data)

         minimise_gp_gradient = minim(xx_dp,likelihood_function,likelihood_gradient,&
         & method='cg',convergence_tol=do_convergence_tol,max_steps=do_max_steps,&
         & hook = save_likelihood_parameters, hook_print_interval=1,data=am_data, always_do_test_gradient=always_do_test_gradient)

         deallocate(xx,xx_dp,am_data)
         am%minim_gp => null()
         am%theta_0 => null()
         if(allocated(theta_0)) deallocate(theta_0)

      end function minimise_gp_gradient

      function minimise_gp_simple_gradient(this,convergence_tol,max_steps,always_do_test_gradient,sigma,delta,theta,f0,theta_fac)

         type(gp), intent(inout), target :: this                       !% GP
         integer, intent(in), optional :: max_steps            !% maximum number of steps, default: 100
         real(dp), intent(in), optional :: convergence_tol     !% convergence tolerance, default: 0.1
         logical, intent(in), optional :: always_do_test_gradient
         logical, intent(in), optional :: sigma, delta, theta, f0, theta_fac
         integer                 :: minimise_gp_simple_gradient
         real(dp), dimension(:), allocatable :: xx_dp
         real(qp), dimension(:), allocatable :: xx
         integer :: n, li, ui

         integer :: do_max_steps
         real(dp) :: do_convergence_tol
         type(gp_simple_minimise) :: am
         character, dimension(:), allocatable :: am_data
         character, dimension(1) :: am_mold
         integer :: am_data_size
         real(qp), dimension(:,:), allocatable, target :: theta_0

         do_max_steps = optional_default(100,max_steps)
         do_convergence_tol = optional_default(0.0001_dp,convergence_tol)

         am_data_size = size(transfer(am,am_mold))
         allocate(am_data(am_data_size))

         am%minim_gp => this

         am%do_sigma = optional_default(.true.,sigma)
         am%do_delta = optional_default(.true.,delta)
         am%do_theta = optional_default(.true.,theta)
         am%do_f0 = optional_default(.true.,f0)
         am%do_theta_fac = optional_default(.false.,theta_fac)

         if(am%do_theta .and. am%do_theta_fac) then
            call print_warning('minimise_gp_gradient called with both theta and theta_fac active. Switching theta_fac to false.')
            am%do_theta_fac = .false.
         endif

         if(am%do_theta_fac) then
            allocate(theta_0(this%d,this%nsp))
            theta_0 = this%theta
            am%theta_0 => theta_0
         endif

         minimise_gp_simple_gradient = 0
         if( (.not.am%do_sigma) .and. (.not.am%do_delta) .and. (.not.am%do_theta) .and. (.not.am%do_f0) .and. (.not.am%do_theta_fac) ) return

         n = 0
         li = 0
         ui = 0
         if( am%do_sigma ) n = n + 1
         if( am%do_delta ) n = n + this%nsp
         if( am%do_theta ) n = n + this%d*this%nsp
         if( am%do_f0 ) n = n + this%nsp
         if( am%do_theta_fac ) n = n + this%nsp

         allocate( xx(n), xx_dp(n) )

         if(am%do_sigma) then
            li = li + 1
            ui = li
            am%li_sigma = li
            am%ui_sigma = ui
            xx(li:ui) = this%sigma
         endif
         if(am%do_delta) then
            li = ui + 1
            ui = li + this%nsp - 1
            am%li_delta = li
            am%ui_delta = ui
            xx(li:ui) = this%delta
         endif
         if(am%do_theta) then
            li = ui + 1
            ui = li + this%d*this%nsp - 1
            am%li_theta = li
            am%ui_theta = ui
            xx(li:ui) = reshape(this%theta,(/this%d*this%nsp/))
         endif
         if(am%do_f0) then
            li = ui + 1
            ui = li + this%nsp - 1
            am%li_f0 = li
            am%ui_f0 = ui
            xx(li:ui) = this%f0
         endif
         if(am%do_theta_fac) then
            li = ui + 1
            ui = li + this%nsp - 1
            am%li_theta_fac = li
            am%ui_theta_fac = ui
            xx(li:ui) = 1.0_qp
         endif

         xx_dp = real(xx,kind=dp)
         am_data = transfer(am, am_data)

         minimise_gp_simple_gradient = minim(xx_dp,likelihood_simple_function,likelihood_simple_gradient,&
         & method='cg',convergence_tol=do_convergence_tol,max_steps=do_max_steps,&
         & hook = save_likelihood_parameters, hook_print_interval=1,data=am_data, always_do_test_gradient=always_do_test_gradient)

         deallocate(xx,xx_dp,am_data)
         am%minim_gp => null()
         am%theta_0 => null()
         if(allocated(theta_0)) deallocate(theta_0)

      end function minimise_gp_simple_gradient

      subroutine minimise_gp_ns_new(this,N_live,max_steps)
         type(gp_sparse), intent(inout) :: this                       !% GP
         integer, optional, intent(in)         :: N_live              !% number of living points
         integer, optional, intent(in)         :: max_steps            !% number of living points

         integer :: i, my_N_live, my_max_steps, stat
         integer, dimension(1) :: ml
         integer, dimension(:,:), allocatable :: r
         real(qp) :: lklhd_old, lklhd_new
         real(qp), dimension(:), allocatable :: lklhd

         my_N_live = optional_default(20,N_live)
         my_max_steps = optional_default(100,max_steps)

         allocate(r(this%sr,my_N_live), lklhd(my_N_live) )

         do i = 1, my_N_live
            call fill_random_integer(r(:,i),this%nxx)
            this%x_sparse = this%x(:,r(:,i))
!            call randomise(this%x_sparse,0.2_dp)
            call covariance_matrix_sparse(this)
!            stat = minimise_gp_gradient(this,max_steps=5,sigma=.false.,delta=.false.,sparx=.false.)
            call likelihood(this, lklhd(i))
         enddo

         do i = 1, my_max_steps
            ml = minloc(lklhd)
            lklhd_old = lklhd(ml(1))
            do
               call fill_random_integer(r(:,ml(1)),this%nxx)
               this%x_sparse = this%x(:,r(:,ml(1)))
!               call randomise(this%x_sparse,0.2_dp)
               call covariance_matrix_sparse(this)
!               stat = minimise_gp_gradient(this,max_steps=5,sparx=.false.)
               call likelihood(this, lklhd_new)
               if(lklhd_new > lklhd_old) exit
            enddo
            lklhd(ml(1)) = lklhd_new
            !call print('Iteration: '//i//', function value dropped: '//lklhd_old,verbosity=PRINT_NORMAL)
         enddo

         ml = maxloc(lklhd)
         this%x_sparse = this%x(:,r(:,ml(1)))

         call covariance_matrix_sparse(this)
         stat = minimise_gp_gradient(this,max_steps=15,sigma=.false.,delta=.false.,theta=.true.,sparx=.false.)
         stat = minimise_gp_gradient(this,max_steps=10,sigma=.true.,delta=.true.,theta=.false.,sparx=.false.)
         stat = minimise_gp_gradient(this,max_steps=15,sigma=.false.,delta=.false.,theta=.false.,sparx=.true.)
         stat = minimise_gp_gradient(this,max_steps=10,sigma=.true.,delta=.true.,theta=.false.,sparx=.false.)

         deallocate(r, lklhd)
      end subroutine minimise_gp_ns_new

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Nested sampling routine.
      !% Maximises a function.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      function ns(x_inout,func,x_min,x_max,N_live,max_iter,tighten)

        real(dp), dimension(:), intent(inout) :: x_inout         !% x value to start at
        real(dp), dimension(:), intent(inout) :: x_min, x_max    !% boundary for x
        integer, optional, intent(in)         :: N_live          !% number of living points
        integer, optional, intent(in)         :: max_iter        !% number of maximum iterations
        logical, optional, intent(in)         :: tighten         !% whether to adapt the boundaries
        logical :: ns

        interface
           function func(x)
             use system_module
             real(dp) :: func
             real(dp), dimension(:), intent(in) :: x
           end function func
        end interface

        integer               :: iter, i, do_N_live, do_max_iter
        integer, dimension(1) :: ml
        logical :: do_tighten

        real(dp) :: e_old, e_new, e_old_prev
        real(dp), dimension(:), allocatable   :: e, x_temp
        real(dp), dimension(:,:), allocatable :: x_hold

        call print_title('Nested sampling',verbosity=PRINT_NORMAL)
        call print('Number of variables is '//size(x_inout),verbosity=PRINT_NORMAL)

        ns = .false.
        e_old_prev = huge(1.0_dp)

        do_N_live = optional_default(100,N_live)
        do_max_iter = optional_default(1000,max_iter)
        do_tighten = optional_default(.false.,tighten)

        allocate( e(do_N_live), x_hold(size(x_inout),do_N_live), x_temp(size(x_inout)) )

        x_hold(:,1) = x_inout
        e(1) = func( x_hold(:,1) )

        do i = 2, do_N_live
           call get_random_vector( x_hold(:,i), x_min, x_max )
           e(i) = func( x_hold(:,i) )
        enddo
        
        do iter = 1, do_max_iter
           ml = minloc(e)
           e_old = e(ml(1))

           do
              call get_random_vector(x_temp,x_min,x_max)
              e_new = func(x_temp)
              if (e_new > e_old) exit
           enddo

           e(ml(1)) = e_new
           x_hold(:,ml(1)) = x_temp(:)

           if( do_tighten ) then
              x_min = minval( x_hold, dim=2)
              x_max = maxval( x_hold, dim=2)
           endif

           call print('Iteration: '//iter//', function value dropped: '//e_old,verbosity=PRINT_NORMAL)
           call print('Point dropped: '//x_temp(:),verbosity=PRINT_NORMAL)

           e_old_prev = e_old
        enddo

        ml = maxloc(e)
        x_inout(:) = x_hold(:,ml(1))

        deallocate( e, x_hold, x_temp )

        ns = .true.

        contains

          subroutine get_random_vector(x,x_min,x_max)

             real(dp), dimension(:), intent(out) :: x
             real(dp), dimension(:), intent(in) :: x_min, x_max

             call random_number(x)

             x = x * (x_max-x_min) + x_min

          end subroutine get_random_vector

      end function ns

      function likelihood_function(x_in, am_data) result(l)

         real(dp), dimension(:) :: x_in
         real(dp)               :: l
         character,optional     :: am_data(:)
         real(qp), dimension(:), allocatable :: x
         real(qp) :: l_qp

         type(gp_minimise) :: am

         l = 0.0_dp
         allocate(x(size(x_in)))
         x = real(x_in,kind=qp)

         am = transfer(am_data,am)

         if( am%do_sigma ) call gp_update( am%minim_gp, sigma_in=x(am%li_sigma:am%ui_sigma), do_covariance = .false. )
         if( am%do_delta ) call gp_update( am%minim_gp, delta_in=x(am%li_delta:am%ui_delta), do_covariance = .false. )
         if( am%do_theta ) call gp_update( am%minim_gp, theta_in= &
         & reshape(x(am%li_theta:am%ui_theta),(/am%minim_gp%d,am%minim_gp%nsp/)), &
         & do_covariance = .false. )
         if( am%do_sparx ) call gp_update( am%minim_gp, &
         & x_sparse_in=reshape(x(am%li_sparx:am%ui_sparx),(/am%minim_gp%d,am%minim_gp%sr/)), do_covariance = .false. )
         if( am%do_f0 ) call gp_update( am%minim_gp, f0_in=x(am%li_f0:am%ui_f0), do_covariance = .false. )
         if( am%do_theta_fac ) call gp_update( am%minim_gp, theta_in= &
         & spread(x(am%li_theta_fac:am%ui_theta_fac),dim=1,ncopies=am%minim_gp%d)*am%theta_0, do_covariance = .false. )

         call covariance_matrix_sparse(am%minim_gp)

         am_data = transfer(am, am_data)

         call likelihood(am%minim_gp,l=l_qp)
         l = -real(l_qp,kind=dp)
         deallocate(x)

      end function likelihood_function

      function likelihood_simple_function(x_in, am_data) result(l)

         real(dp), dimension(:) :: x_in
         real(dp)               :: l
         character,optional     :: am_data(:)
         real(qp), dimension(:), allocatable :: x
         real(qp) :: l_qp

         type(gp_simple_minimise) :: am

         l = 0.0_dp
         allocate(x(size(x_in)))
         x = real(x_in,kind=qp)

         am = transfer(am_data,am)

         if( am%do_sigma ) call gp_update_simple( am%minim_gp, sigma_in=x(am%li_sigma), do_covariance = .false. )
         if( am%do_delta ) call gp_update_simple( am%minim_gp, delta_in=x(am%li_delta:am%ui_delta), do_covariance = .false. )
         if( am%do_theta ) call gp_update_simple( am%minim_gp, theta_in= &
         reshape(x(am%li_theta:am%ui_theta),(/am%minim_gp%d,am%minim_gp%nsp/)), &
         do_covariance = .false. )
         if( am%do_f0 ) call gp_update_simple( am%minim_gp, f0_in=x(am%li_f0:am%ui_f0), do_covariance = .false. )
         if( am%do_theta_fac ) call gp_update_simple( am%minim_gp, theta_in= &
         spread(x(am%li_theta_fac:am%ui_theta_fac),dim=1,ncopies=am%minim_gp%d)*am%theta_0, do_covariance = .false. )

         call covariance_matrix_simple(am%minim_gp)

         am_data = transfer(am, am_data)

         call likelihood_simple(am%minim_gp,l=l_qp)
         l = -real(l_qp,kind=dp)
         deallocate(x)

      end function likelihood_simple_function

      function likelihood_gradient(x_in,am_data) result(dl_out)

         real(dp), dimension(:)          :: x_in
         real(dp), dimension(size(x_in)) :: dl_out
         character,optional              :: am_data(:)

         real(qp), dimension(:), allocatable :: x, dl, dl_dtheta
         type(gp_minimise) :: am

         allocate(x(size(x_in)),dl(size(x_in)))

         x = real(x_in,kind=qp)
         dl = 0.0_qp

         am = transfer(am_data,am)

         if( am%do_sigma ) call gp_update( am%minim_gp, sigma_in=x(am%li_sigma:am%ui_sigma), do_covariance = .false. )
         if( am%do_delta ) call gp_update( am%minim_gp, delta_in=x(am%li_delta:am%ui_delta), do_covariance = .false. )
         if( am%do_theta ) call gp_update( am%minim_gp, theta_in= &
         & reshape(x(am%li_theta:am%ui_theta),(/am%minim_gp%d,am%minim_gp%nsp/)), &
         & do_covariance = .false. )
         if( am%do_sparx ) call gp_update( am%minim_gp, &
         & x_sparse_in=reshape(x(am%li_sparx:am%ui_sparx),(/am%minim_gp%d,am%minim_gp%sr/)), do_covariance = .false. )
         if( am%do_f0 ) call gp_update( am%minim_gp, f0_in=x(am%li_f0:am%ui_f0), do_covariance = .false. )
         if( am%do_theta_fac ) call gp_update( am%minim_gp, theta_in= &
         & spread(x(am%li_theta_fac:am%ui_theta_fac),dim=1,ncopies=am%minim_gp%d)*am%theta_0, do_covariance = .false. )

         if(am%do_theta .or. am%do_theta_fac) allocate(dl_dtheta(am%minim_gp%d*am%minim_gp%nsp))

         call covariance_matrix_sparse(am%minim_gp)

         am_data = transfer(am, am_data)

         call likelihood(am%minim_gp,dl_dsigma=dl(am%li_sigma:am%ui_sigma),do_sigma = am%do_sigma, &
         & dl_ddelta=dl(am%li_delta:am%ui_delta), do_delta = am%do_delta, &
         & dl_dtheta=dl_dtheta, do_theta = (am%do_theta .or.am%do_theta_fac) , &
         & dl_dx = dl(am%li_sparx:am%ui_sparx), do_x = am%do_sparx, &
         & dl_df0=dl(am%li_f0:am%ui_f0), do_f0 = am%do_f0 )

         if( am%do_theta ) dl(am%li_theta:am%ui_theta) = dl_dtheta
         if( am%do_theta_fac ) dl(am%li_theta_fac:am%ui_theta_fac) = sum(reshape(dl_dtheta,(/am%minim_gp%d,am%minim_gp%nsp/))*am%theta_0,dim=1)

         dl_out = -real(dl,kind=dp)
         deallocate(x,dl)
         if(allocated(dl_dtheta)) deallocate(dl_dtheta)

      end function likelihood_gradient

      function likelihood_simple_gradient(x_in,am_data) result(dl_out)

         real(dp), dimension(:)          :: x_in
         real(dp), dimension(size(x_in)) :: dl_out
         character,optional              :: am_data(:)

         real(qp), dimension(:), allocatable :: x, dl, dl_dtheta
         type(gp_simple_minimise) :: am

         allocate(x(size(x_in)),dl(size(x_in)))

         x = real(x_in,kind=qp)
         dl = 0.0_qp

         am = transfer(am_data,am)

         if( am%do_sigma ) call gp_update_simple( am%minim_gp, sigma_in=x(am%li_sigma), do_covariance = .false. )
         if( am%do_delta ) call gp_update_simple( am%minim_gp, delta_in=x(am%li_delta:am%ui_delta), do_covariance = .false. )
         if( am%do_theta ) call gp_update_simple( am%minim_gp, theta_in= &
         reshape(x(am%li_theta:am%ui_theta),(/am%minim_gp%d,am%minim_gp%nsp/)), &
         do_covariance = .false. )
         if( am%do_f0 ) call gp_update_simple( am%minim_gp, f0_in=x(am%li_f0:am%ui_f0), do_covariance = .false. )
         if( am%do_theta_fac ) call gp_update_simple( am%minim_gp, theta_in= &
         spread(x(am%li_theta_fac:am%ui_theta_fac),dim=1,ncopies=am%minim_gp%d)*am%theta_0, do_covariance = .false. )

         if(am%do_theta .or. am%do_theta_fac) allocate(dl_dtheta(am%minim_gp%d*am%minim_gp%nsp))

         call covariance_matrix_simple(am%minim_gp)

         am_data = transfer(am, am_data)

         call likelihood_simple(am%minim_gp,dl_dsigma=dl(am%li_sigma),do_sigma = am%do_sigma, &
         dl_ddelta=dl(am%li_delta), do_delta = am%do_delta, &
         dl_dtheta=dl_dtheta, do_theta = (am%do_theta .or.am%do_theta_fac) , &
         dl_df0=dl(am%li_f0), do_f0 = am%do_f0 )

         if( am%do_theta ) dl(am%li_theta:am%ui_theta) = dl_dtheta
         if( am%do_theta_fac ) dl(am%li_theta_fac:am%ui_theta_fac) = sum(reshape(dl_dtheta,(/am%minim_gp%d,am%minim_gp%nsp/))*am%theta_0,dim=1)

         dl_out = -real(dl,kind=dp)
         deallocate(x,dl)
         if(allocated(dl_dtheta)) deallocate(dl_dtheta)

      end function likelihood_simple_gradient

      subroutine save_likelihood_parameters(x,dx,e,done,do_print,data)
         real(dp), dimension(:), intent(in) :: x
         real(dp), dimension(:), intent(in) :: dx
         real(dp), intent(in) :: e
         logical, intent(out) :: done
         logical, optional, intent(in) :: do_print
         character,optional, intent(in) :: data(:)

         if( present(do_print) ) then
             if( do_print ) call print('SAVE HYPERS:'//x)
         endif
         done = .false.

      end subroutine save_likelihood_parameters

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Prints out GP in binary format to a file.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_print_binary(this,filename)

         type(gp), intent(in) :: this
         character(len=*), intent(in) :: filename

         integer :: i, j

         if( .not. this%initialised ) call system_abort('gp_print_binary: gp not initialised')

         open(unit=666, file=trim(filename), form='unformatted')

         write(666) this%d, this%nsp, this%n, this%sigma

         do i = 1, this%nsp
            write(666) this%delta(i)
         enddo

         do j = 1, this%nsp
            do i = 1, this%d
               write(666) this%theta(i,j)
            enddo
         enddo

         do i = 1, this%n
            do j = 1, this%d
               write(666) this%x(j,i)
            enddo
         enddo

         do i = 1, this%n
            do j = 1, this%n
               write(666) this%c(j,i)
            enddo
         enddo

         do i = 1, this%n
            write(666) this%alpha(i)
         enddo

         do i = 1, this%n
            write(666) this%xz(i)
         enddo

         do i = 1, this%nsp
            write(666) this%sp(i)
         enddo

         do i = 1, this%nsp
            write(666) this%f0(i)
         enddo

         write(666) this%comment

         close(unit=666)
         
      end subroutine gp_print_binary

      subroutine gp_read_binary(this,filename)

         type(gp), intent(inout) :: this
         character(len=*), intent(in) :: filename

         integer :: i, j, stat
         logical :: file_exist

         inquire(file=filename,exist=file_exist)

         if( .not. file_exist ) call system_abort('gp_read: '//trim(filename)//' does not exist')

         if( this%initialised ) call finalise(this)

         open(unit=666, file=trim(filename), form='unformatted', action='read')

         read(666) this%d, this%nsp, this%n, this%sigma

         allocate(this%theta(this%d,this%nsp), this%delta(this%nsp), this%x(this%d,this%n), &
         & this%c(this%n,this%n), this%alpha(this%n), this%xz(this%n), this%sp(this%nsp), this%f0(this%nsp))

         do i = 1, this%nsp
            read(666) this%delta(i)
         end do

         do j = 1, this%nsp
            do i = 1, this%d
               read(666) this%theta(i,j)
            enddo
         enddo

         do i = 1, this%n
            do j = 1, this%d
               read(666) this%x(j,i)
            enddo
         enddo

	 call setup_x_div_theta(this)

         do i = 1, this%n
            do j = 1, this%n
               read(666) this%c(j,i)
            enddo
         enddo

         do i = 1, this%n
            read(666) this%alpha(i)
         enddo

         do i = 1, this%n
            read(666,iostat=stat) this%xz(i)
            if( stat/=0 ) then
               this%xz = 6
               backspace(unit=666)
               backspace(unit=666)
               exit
            endif
         enddo

         do i = 1, this%nsp
            read(666,iostat=stat) this%sp(i)
            if( stat/=0 ) then
               this%sp = 6
               backspace(unit=666)
               backspace(unit=666)
               exit
            endif
         enddo
         do i = 1, this%nsp
            read(666) this%f0(i)
         enddo

         read(666,iostat=stat) this%comment
         if( stat/=0 ) this%comment = ''
         close(unit=666)
         
         this%initialised = .true.
      end subroutine gp_read_binary

      subroutine fill_random_integer(r,n,b)
         integer, dimension(:), intent(out) :: r
         integer, intent(in) :: n
         integer, dimension(:), intent(in), optional :: b

         integer :: i, i0, irnd, sr
         real(qp) :: rnd

         sr = size(r)
         if( sr > n ) call system_abort('fill_random_integer: cannot fill array, too short')
         r = 0
         i0 = 1
         if( present(b) ) then
            if(size(b) > sr) call system_abort('fill_random_integer: cannot fill array, optinal b too long')
            r(1:size(b)) = b
            i0 = size(b) + 1
         endif

         do i = i0, sr
            do 
               call random_number(rnd)
               irnd = ceiling(rnd*n)
               if( all(r /= irnd) ) exit
            enddo
            r(i) = irnd
         enddo

      end subroutine fill_random_integer

end module gp_teach_module
