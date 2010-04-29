! HEADER_CHECK

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X Gaussian Process module
!X
!% Module for general GP function interpolations.
!% A gp object contains the training set (teaching points and function values),
!% important temporary matrices, vectors and parameters.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module gp_module

   use libatoms_module

   implicit none

   integer, parameter :: gp_n_max_length = 1000
   integer, parameter :: gp_m_max_length = 1000
   integer, parameter :: gp_n_increment  = 1000
   integer, parameter :: gp_m_increment  = 1000

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !
   !% Gaussian Process data type.
   !% Contains all the important data necessary to estimate function values,
   !% derivatives, optimising hyperparameters.
   !
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   type gp

      !% Data arrays. These are (re)allocated quite rarely, always bigger than they should be.
      !% Therefore adding new data points is easy and fast.
      real(dp), dimension(:,:), allocatable :: x_data, x_prime_data, c_data, c_factorised_data, lcl_data
      real(dp), dimension(:), allocatable   :: y_data, alpha_data
      integer, dimension(:), allocatable    :: l_data
      logical, dimension(:), allocatable    :: is_prime_data

      !% Pointer arrays. These point to the data arrays, they are reassociated every time a new
      !% teaching point is included.
      real(dp), dimension(:,:), pointer :: x, x_prime, c, c_factorised, lcl
      real(dp), dimension(:), pointer   :: y, alpha
      integer, dimension(:), pointer    :: l
      logical, dimension(:), pointer    :: is_prime

      !% Hyperparameters.
      !% Error parameter, range parameter, sensitivity parameters.
      real(dp), dimension(:), allocatable       :: theta

      !% Vector sizes.
      !% d: dimensionality of input space
      !% n: number of teaching points
      !% m: number of teaching function values
      !% p: number of hyperparameters
      integer :: d, n, m, p
      integer :: n_max_length = gp_n_max_length
      integer :: m_max_length = gp_m_max_length

      logical :: initialised = .false.

      !% Which covariance function to use.
      character(len=1024) :: cov_func_string, label

   endtype gp

   interface Initialise
     module procedure GP_Initialise
   endinterface Initialise

   interface Finalise
     module procedure GP_Finalise, GP_Finalise_more
   endinterface Finalise

   interface GP_Add_Point
     module procedure gp_add_single_point, gp_add_multiple_points
   endinterface GP_Add_Point

   interface covSEiso
     module procedure covSEiso_diag, covSEiso_gen
   endinterface covSEiso

   interface covSEard
     module procedure covSEard_diag, covSEard_gen
   endinterface covSEard

   interface apply_llt
     module procedure apply_llt_matrix, apply_llt_vector, apply_llt_gp
   endinterface apply_llt

   private
   public :: gp, initialise, finalise, gp_update, gp_mean, gp_variance, gp_predict, &
           & gp_add_point, likelihood, likelihood_gradient, test_gp_gradient, minimise_gp_gradient, &
           & minimise_gp_ns, gp_print_binary, gp_read_binary

   contains

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Initialise a gp object from scratch.
      !% Inputs needed:
      !% cov_func_in: 'iso' or 'ard', depending on which covariance function to use.
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

      subroutine GP_initialise(this, cov_func_in, theta_in, y_in, x_in, x_prime_in, l_in)

         !% n: number of teaching points \\
         !% d: dimension of teaching vectors \\
         !% p: number of parameters \\

         type(gp), intent(inout)                        :: this        !% gp to initialise
         character(len=*), intent(in)                   :: cov_func_in !% covariance function to use (iso or ard)
         real(dp), dimension(:), intent(in)             :: theta_in    !% hyperparameters, dimension(p)
         real(dp), dimension(:), intent(in)             :: y_in        !% function (or derivative) value at teaching points, dimension(m)
         real(dp), dimension(:,:), intent(in)           :: x_in        !% teaching points, dimension(d,n)
         real(dp), dimension(:,:), intent(in), optional :: x_prime_in  !% derivative of teaching points, dimension(d,n)
         integer, dimension(:), intent(in), optional    :: l_in        !% how the sum is formed from teaching point(m)

         integer :: i

         if( this%initialised ) return

         this%d = size(x_in,1) 
         this%n = size(x_in,2)
         this%m = size(y_in)
         this%p = size(theta_in)

         if(present(x_prime_in)) then
            if( (size(x_prime_in,1) /= this%d) .or. (size(x_prime_in,2) /= this%n) ) &
            & call system_abort('GP_Initialise: x_in and x_prime_in do not conform')
         endif

         if( present(l_in) ) then
             if( size(l_in,1) /= this%m ) call system_abort('GP_Initialise: y_in and l_in do not conform')
         endif

         if( (this%m /= this%n) .and. .not. present(l_in) ) &
         & call system_abort('GP_Initialise: if values y_in are linear combinations of values at teaching points, &
         & there has to be an l_in to show how the linear combination is done.')

         if( this%p /= num_parameters(this%d,cov_func_in) ) &
         & call system_abort('GP_initialise: number of parameters incorrect')

         call gp_allocate(this)
         call gp_pointers(this)

         this%cov_func_string = trim(cov_func_in)
         this%theta = theta_in
         this%y = y_in
         this%x = x_in
         if(present(x_prime_in)) then
            this%x_prime = x_prime_in
            this%is_prime = .true.
         endif
         if(present(l_in)) then
            this%l = l_in
         else
            forall(i=1:this%m) this%l(i) = i
         endif

         call covariance_matrix(this)

         this%initialised = .true.

      endsubroutine GP_initialise

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Adds a teaching point to an initialised gp object.
      !% This subroutine expects a singe value (or derivative) and a single teaching point and
      !% derivative, if necessary.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_add_single_point(this, y_in, x_in, x_prime_in)

         type(gp), intent(inout)                      :: this       !% gp object
         real(dp), intent(in)                         :: y_in       !% function (or derivative) value at teaching point
         real(dp), dimension(:), intent(in)           :: x_in       !% teaching vector, dimension(d)
         real(dp), dimension(:), intent(in), optional :: x_prime_in !% derivative of teaching vector, dimension(d)

         if( .not. this%initialised ) &
         & call system_abort('gp_add_single_point: GP not initialised')

         if(present(x_prime_in) ) then
            call gp_add_multiple_points(this, y_in, reshape(x_in,(/this%d,1/)),reshape(x_prime_in,(/this%d,1/)))
         else
            call gp_add_multiple_points(this, y_in, reshape(x_in,(/this%d,1/)))
         endif

      endsubroutine gp_add_single_point

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Adds multiple teaching points to an initialised gp object.
      !% Here a sum of a function at different teaching points are expected.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_add_multiple_points(this, y_in, x_in, x_prime_in)

         type(gp), intent(inout)                        :: this       !% gp object
         real(dp), intent(in)                           :: y_in       !% sum of function values (or derivatives) at teaching points
         real(dp), dimension(:,:), intent(in)           :: x_in       !% teaching vectors, dimension(d,n)
         real(dp), dimension(:,:), intent(in), optional :: x_prime_in !% derivatives of teaching vectors, dimension(d,n)

         integer :: n, m, n_old

         if( .not. this%initialised ) &
         & call system_abort('gp_add_multiple_points: GP not initialised')

         if( size(x_in,1) /= this%d ) &
         & call system_abort('gp_add_multiple_points: gp requires '//this%d//' dimension teaching vector, &
         & given: '//size(x_in,1))

         if( present(x_prime_in) ) then
             if( size(x_prime_in,1) /= this%d ) &
             & call system_abort('gp_add_multiple_point: gp requires '//this%d//' dimension gradient of teaching vector, &
             & given: '//size(x_prime_in,1))
         endif

         n_old = this%n
         n = n_old + size(x_in,2) 
         m = this%m + 1

         call gp_extend(this,n,m)

         this%x(:,n_old+1:) = x_in(:,:)
         if(present(x_prime_in)) then
            this%x_prime(:,n_old+1:) = x_prime_in(:,:)
            this%is_prime(n_old+1:) = .true.
         endif
         this%y(m) = y_in
         this%l(m) = n

         call covariance_matrix_extend(this)

      endsubroutine gp_add_multiple_points

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Extends gp. If necessary, reallocates data arrays, otherwise reassociate pointers only.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_extend(this,n,m)

         type(gp), intent(inout)       :: this    !% GP object
         integer, intent(in)           :: n, m    !% new dimensions

         real(dp), dimension(:,:), allocatable :: x_data, x_prime_data, c_data, c_factorised_data, lcl_data
         real(dp), dimension(:), allocatable   :: y_data, alpha_data
         integer, dimension(:), allocatable    :: l_data
         logical, dimension(:), allocatable    :: is_prime_data

         if( this%n_max_length < n ) then
             allocate( x_data(this%d,this%n), x_prime_data(this%d,this%n), c_data(this%n,this%n), &
             & is_prime_data(this%n) )
             x_data        = this%x
             x_prime_data  = this%x_prime
             c_data        = this%c
             is_prime_data = this%is_prime

             call gp_deallocate_n(this)
             this%n_max_length = this%n_max_length + gp_n_increment
             call gp_allocate_n(this)

             this%x_data(:,:this%n)       = x_data
             this%x_prime_data(:,:this%n) = x_prime_data
             this%c_data(:this%n,:this%n) = c_data
             this%is_prime_data(:this%n)  = is_prime_data

             deallocate( x_data, x_prime_data, c_data, is_prime_data )
         endif

         if( this%m_max_length < m ) then
             allocate( y_data(this%m), l_data(this%m), c_factorised_data(this%m,this%m), alpha_data(this%m) )
             if( allocated(this%lcl_data) ) allocate(lcl_data(this%m,this%m) )

             y_data  = this%y
             l_data  = this%l
             c_factorised_data = this%c_factorised
             alpha_data = this%alpha
             if( allocated(this%lcl_data) ) lcl_data = this%lcl

             call gp_deallocate_m(this)
             this%m_max_length = this%m_max_length + gp_m_increment
             call gp_allocate_m(this)

             this%y_data(:this%m) = y_data
             this%l_data(:this%m) = l_data
             this%c_factorised_data(:this%m,:this%m) = c_factorised_data
             this%alpha_data(:this%m) = alpha_data
             if( allocated(lcl_data) ) this%lcl_data(:this%m,:this%m) = lcl_data

             deallocate( y_data, l_data, c_factorised_data, alpha_data )
             if( allocated(lcl_data) ) deallocate( lcl_data )
         endif

         if( (this%m==this%n).and.(m/=n) ) then
             allocate( this%lcl_data(this%m_max_length,this%m_max_length) )
             this%lcl_data = 0.0_dp
             this%lcl_data(:this%m,:this%m) = this%c_data(:this%m,:this%m)
         endif

         this%n = n
         this%m = m

         call gp_pointers(this)

      endsubroutine gp_extend

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Reassociates pointer arrays.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_pointers(this)
         type(gp), target, intent(inout)       :: this    !% GP object

         this%x        => this%x_data(:,:this%n)
         this%x_prime  => this%x_prime_data(:,:this%n)
         this%c        => this%c_data(:this%n,:this%n)
         this%is_prime => this%is_prime_data(:this%n)

         this%y => this%y_data(:this%m)
         this%l => this%l_data(:this%m)
         this%c_factorised => this%c_factorised_data(:this%m,:this%m)
         this%alpha => this%alpha_data(:this%m)

         if( this%m == this%n ) then
             this%lcl => this%c
         else
             this%lcl => this%lcl_data(:this%m,:this%m)
         endif

      endsubroutine gp_pointers

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Nullifies arrays.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_nullify(this)
         type(gp), intent(inout)       :: this !% GP object

         this%x        => null()
         this%x_prime  => null()
         this%c        => null()
         this%is_prime => null()

         this%y => null()
         this%l => null()
         this%c_factorised => null()
         this%alpha => null()
         this%lcl => null()

      endsubroutine gp_nullify

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Allocates data arrays, storing the $\mathbf{x}$ vectors or similarly sized other arrays.
      !% These are the `n' dimensional arrays.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_allocate_n(this)
         type(gp), intent(inout)       :: this

         call gp_deallocate_n(this)

         if( this%n > this%n_max_length ) &
         & this%n_max_length = ceiling(real(this%n,dp) / real(this%n_max_length,dp)) * this%n_max_length

         allocate( this%x_data(this%d,this%n_max_length), this%x_prime_data(this%d,this%n_max_length), &
         & this%c_data(this%n_max_length,this%n_max_length), this%is_prime_data(this%n_max_length) )

         this%x_data        = 0.0_dp
         this%x_prime_data  = 0.0_dp
         this%c_data        = 0.0_dp
         this%is_prime_data = .false.

      endsubroutine gp_allocate_n

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Allocates data arrays, storing the function values or similarly sized other arrays.
      !% These are the `m' dimensional arrays.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_allocate_m(this)
         type(gp), intent(inout)       :: this

         call gp_deallocate_m(this)

         if( this%m > this%m_max_length ) &
         & this%m_max_length = ceiling(real(this%m,dp) / real(this%m_max_length,dp)) * this%m_max_length

         allocate( this%y_data(this%m_max_length), this%l_data(this%m_max_length), &
         & this%c_factorised_data(this%m_max_length,this%m_max_length), this%alpha_data(this%m_max_length) )

         this%y_data        = 0.0_dp
         this%l_data        = 0
         this%c_factorised_data = 0.0_dp
         this%alpha_data = 0.0_dp

         if( this%n /= this%m ) allocate( this%lcl_data(this%m_max_length,this%m_max_length) )

      endsubroutine gp_allocate_m

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Deallocates data arrays, storing the $\mathbf{x}$ vectors or similarly sized other arrays.
      !% These are the `n' dimensional arrays.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_deallocate_n(this)
         type(gp), intent(inout)       :: this

         if(allocated(this%x_data))            deallocate( this%x_data )
         if(allocated(this%x_prime_data))      deallocate( this%x_prime_data )
         if(allocated(this%c_data))            deallocate( this%c_data )
         if(allocated(this%is_prime_data))     deallocate( this%is_prime_data )

      endsubroutine gp_deallocate_n

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Deallocates data arrays, storing the function values or similarly sized other arrays.
      !% These are the `m' dimensional arrays.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_deallocate_m(this)
         type(gp), intent(inout)       :: this

         if(allocated(this%y_data))            deallocate( this%y_data )
         if(allocated(this%l_data))            deallocate( this%l_data )
         if(allocated(this%c_factorised_data)) deallocate( this%c_factorised_data )
         if(allocated(this%alpha_data))        deallocate( this%alpha_data )
         if(allocated(this%lcl_data))          deallocate( this%lcl_data )

      endsubroutine gp_deallocate_m

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Wrapper to allocate every data array.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_allocate(this)
         type(gp), intent(inout)       :: this

         call gp_deallocate(this)

         allocate( this%theta(this%p) )

         call gp_allocate_n(this)
         call gp_allocate_m(this)

      endsubroutine gp_allocate

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Wrapper to deallocate every data array.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_deallocate(this)

         type(gp), intent(inout) :: this

         if(allocated(this%theta)) deallocate( this%theta )
         call gp_deallocate_n(this)
         call gp_deallocate_m(this)

      endsubroutine gp_deallocate

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Returns the number of hyperparameters given the dimensions of teaching vectors
      !% and covariance function.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      function num_parameters(d,cov_func_string)

         integer, intent(in) :: d
         character(len=*), intent(in) :: cov_func_string
         integer :: num_parameters

         selectcase(trim(cov_func_string))
           case('ard')
             num_parameters = d + 2
           case('iso')
             num_parameters = 3
           case default
             call system_abort('num_parameters: '//trim(cov_func_string)//' not known')
         endselect

      endfunction num_parameters 

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Changes one or more component in a gp object, updating temporary matrices and vectors.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_update(this, x_in, x_prime_in, y_in, theta_in)

         type(gp), intent(inout)                        :: this       !% gp object to update
         real(dp), dimension(:,:), intent(in), optional :: x_in       !% update teaching points, dimension(d,n)
         real(dp), dimension(:,:), intent(in), optional :: x_prime_in !% update derivatives of teaching points, dimension(d,n)
         real(dp), dimension(:), intent(in), optional   :: y_in       !% update function values, dimension(m)
         real(dp), dimension(:), intent(in), optional   :: theta_in   !% update hyperparameters, dimension(p)

         if( .not. this%initialised ) return

         if( present(x_in) ) then
             if( all( shape(this%x) /= shape(x_in) ) ) call system_abort('gp_update: array sizes do not conform')
             this%x = x_in
         endif

         if( present(x_prime_in) ) then
             if( all( shape(this%x_prime) /= shape(x_prime_in) ) ) call system_abort('gp_update: array sizes do not conform')
             this%x_prime = x_prime_in
         endif

         if( present(y_in) ) then
             if( all( shape(this%y) /= shape(y_in) ) ) call system_abort('gp_update: array sizes do not conform')
             this%y = y_in
         endif

         if( present(theta_in) ) then
             if( all( shape(this%theta) /= shape(theta_in) ) ) call system_abort('gp_update: array sizes do not conform')
             this%theta = theta_in
         endif

         if( present(x_in) .or. present(x_prime_in) .or. present(y_in) .or. present(theta_in) ) then
            call covariance_matrix(this)
         endif
         
      endsubroutine gp_update

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Finalise a gp object.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine GP_finalise(this)

         type(gp), intent(inout)              :: this !% gp object

         if( .not. this%initialised ) return

         call gp_deallocate(this)
         call gp_nullify(this)

         this%d = 0
         this%n = 0
         this%m = 0
         this%p = 0 
         this%cov_func_string = ''
         this%initialised = .false.

      endsubroutine GP_finalise

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Finalise more gp objects.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine GP_finalise_more(this, this2, this3, this4, this5, this6, this7, this8, this9)

         type(gp), intent(inout)           :: this, this2
         type(gp), intent(inout), optional :: this3, this4, this5, this6, this7, this8, this9

         call GP_finalise(this)
         call GP_finalise(this2)

         if( present(this3) ) call GP_finalise(this3)
         if( present(this4) ) call GP_finalise(this4)
         if( present(this5) ) call GP_finalise(this5)
         if( present(this6) ) call GP_finalise(this6)
         if( present(this7) ) call GP_finalise(this7)
         if( present(this8) ) call GP_finalise(this8)
         if( present(this9) ) call GP_finalise(this9)

      endsubroutine GP_finalise_more

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Predict a function value or derivative using GP.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      function gp_mean(gp_in,x_star,x_prime_star)

         real(dp) :: gp_mean                                          !% output, predicted value at test point
         type(gp), intent(in)                         :: gp_in        !% GP
         real(dp), dimension(:), intent(in)           :: x_star       !% test point, dimension(d)
         real(dp), dimension(:), intent(in), optional :: x_prime_star !% derivative of test point, dimension(d)

         call gp_predict(gp_data = gp_in, mean = gp_mean, x_star = x_star, x_prime_star = x_prime_star )

      endfunction gp_mean

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Predict the variance of a function value using GP.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      function gp_variance(gp_in,x_star,x_prime_star)

         real(dp) :: gp_variance                         !% output, predicted variance at test point
         type(gp), intent(in)                         :: gp_in        !% GP
         real(dp), dimension(:), intent(in)           :: x_star       !% test point, dimension(d)
         real(dp), dimension(:), intent(in), optional :: x_prime_star !% derivative of test point, dimension(d)

         call gp_predict(gp_data = gp_in, variance = gp_variance, x_star = x_star, x_prime_star = x_prime_star )

      endfunction gp_variance

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Predict the function value and its variance using GP.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_predict(gp_data, mean, variance, x_star, x_prime_star)

         type(gp), intent(in)               :: gp_data                !% GP
         real(dp), intent(out), optional    :: mean, variance         !% output, predicted value and variance at test point
         real(dp), dimension(:), intent(in)           :: x_star       !% test point, dimension(d)
         real(dp), dimension(:), intent(in), optional :: x_prime_star !% derivative of test point, dimension(d)

         real(dp), dimension(:), allocatable :: my_x_prime_star, k, lk, v

         integer :: i

         if( .not. gp_data%initialised ) &
         & call system_abort('gp_predict: not initialised, call gp_initialise first')

         if( size(x_star) /= gp_data%d ) call system_abort('gp_predict: array sizes do not conform')

         if( present(x_prime_star) ) then
             if( size(x_prime_star) /= gp_data%d ) call system_abort('gp_predict: array sizes do not conform')
         endif

         allocate( k(gp_data%n), lk(gp_data%m), v(gp_data%m) )
         allocate( my_x_prime_star(gp_data%d) )
         my_x_prime_star = 0.0_dp
         if( present(x_prime_star) ) my_x_prime_star = x_prime_star

         selectcase(gp_data%cov_func_string)
            case('ard')   
!$omp parallel do
               do i = 1, gp_data%n
                  k(i) = covSEard( x_star, my_x_prime_star, present(x_prime_star), &
                  & gp_data%x(:,i), gp_data%x_prime(:,i), gp_data%is_prime(i), gp_data%theta )
               enddo
            case('iso')   
!$omp parallel do
               do i = 1, gp_data%n
                  k(i) = covSEiso(  x_star, my_x_prime_star, present(x_prime_star), &
                  & gp_data%x(:,i), gp_data%x_prime(:,i), gp_data%is_prime(i), gp_data%theta )
               enddo
            case default
               call system_abort('gp_predict: '//trim(gp_data%cov_func_string)//' not known')
         endselect

         call apply_llt(k,gp_data%l,lk)

         if( present(mean) ) mean = dot_product( lk, gp_data%alpha )

         if( present(variance) ) then

             call Matrix_BackSubstitute(gp_data%c_factorised,v,lk)

             selectcase(gp_data%cov_func_string)
                case('ard')   
                   variance = sqrt( covSEard_diag( my_x_prime_star, present(x_prime_star), gp_data%theta) &
                            & + gp_data%theta(1)**2 - dot_product(lk,v))
                case('iso')   
                   variance =  sqrt( covSEiso_diag( my_x_prime_star, present(x_prime_star), gp_data%theta) &
                            & + gp_data%theta(1)**2 - dot_product(lk,v))
                case default
                   call system_abort('gp_predict: '//trim(gp_data%cov_func_string)//' not known')
             endselect

         endif

         deallocate( k, lk, v )
         deallocate( my_x_prime_star )

      endsubroutine gp_predict

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      ! Squared Exponential covariance function with one characteristic lengthscale.
      ! Covariance between two (xi and xj) vectors, with derivatives dxi and dxj,
      ! whether they are derivatives or not (lxi, lxj).
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      pure function covSEiso_gen(xi,dxi,lxi,xj,dxj,lxj,theta)

         ! Covariance function

         real(dp) :: covSEiso_gen
         real(dp), dimension(:), intent(in) :: xi,xj,dxi,dxj,theta
         logical, intent(in) :: lxi, lxj

         real(dp), dimension(:), allocatable :: diff_xij
         integer :: n

         n = size(xi)

         allocate( diff_xij(n) )

         diff_xij  = xi-xj

         if( .not.(lxi.or.lxj) ) then    
              covSEiso_gen = theta(2)**2 * exp( -0.5_dp * dot_product( diff_xij,diff_xij ) / theta(3)**2 )
         elseif( lxi .and. (.not.lxj) ) then
              covSEiso_gen = - theta(2)**2 * exp( -0.5_dp * dot_product( diff_xij,diff_xij ) / theta(3)**2 ) * &
                       & dot_product( diff_xij,dxi ) / theta(3)**2
         elseif( (.not.lxi) .and. lxj ) then
              covSEiso_gen = theta(2)**2 * exp( -0.5_dp * dot_product( diff_xij,diff_xij ) / theta(3)**2 ) * &
                       & dot_product( diff_xij,dxj ) / theta(3)**2
         else 
              covSEiso_gen = theta(2)**2 * exp( -0.5_dp * dot_product( diff_xij,diff_xij ) / theta(3)**2 ) * &
                       & ( dot_product(dxi,dxj) / theta(3)**2 - &
                       & dot_product(diff_xij,dxi) * dot_product(diff_xij,dxj) / theta(3)**4 )
         endif

         deallocate( diff_xij )

      endfunction covSEiso_gen 

      pure function covSEiso_diag(dxi,lxi,theta)

         ! Covariance function

         real(dp) :: covSEiso_diag                                ! function value
         real(dp), dimension(:), intent(in) :: dxi,theta   ! input variables
         logical, intent(in) :: lxi

         if( .not.lxi ) then    
              covSEiso_diag = theta(2)**2 
         else 
              covSEiso_diag = theta(2)**2 * norm2(dxi) / theta(3)**2 
         endif

      endfunction covSEiso_diag  

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      ! Squared Exponential covariance function with one characteristic lengthscale
      ! for each dimension.
      ! Covariance between two (xi and xj) vectors, with derivatives dxi and dxj,
      ! whether they are derivatives or not (lxi, lxj).
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      pure function covSEard_gen(xi,dxi,lxi,xj,dxj,lxj,theta)

         ! Covariance function

         real(dp) :: covSEard_gen                            
         real(dp), dimension(:), intent(in) :: xi,xj,dxi,dxj,theta
         logical, intent(in) :: lxi, lxj

         real(dp), dimension(:), allocatable :: diff_xij, diff_xijt, theta2
         integer :: n

         n = size(xi)

         allocate( diff_xij(n), diff_xijt(n), theta2(n) )

         diff_xij  = xi-xj
         theta2    = theta(3:)**2
         diff_xijt = diff_xij / theta2

         if( .not.(lxi.or.lxj) ) then    
              covSEard_gen = theta(2)**2 * exp( -0.5_dp * dot_product( diff_xijt,diff_xij ) )
         elseif( lxi .and. (.not.lxj) ) then
              covSEard_gen = - theta(2)**2 * exp( -0.5_dp * dot_product( diff_xijt,diff_xij ) ) * &
                       & dot_product( diff_xijt,dxi )
         elseif( (.not.lxi) .and. lxj ) then
              covSEard_gen = theta(2)**2 * exp( -0.5_dp * dot_product( diff_xijt,diff_xij ) ) * &
                       & dot_product( diff_xijt,dxj )
         else 
              covSEard_gen = theta(2)**2 * exp( -0.5_dp * dot_product( diff_xijt,diff_xij ) ) * &
                       & ( dot_product(dxi/theta2,dxj) - dot_product(diff_xijt,dxi) * dot_product(diff_xijt,dxj) )
         endif

         deallocate( diff_xij, diff_xijt, theta2 )

      endfunction covSEard_gen  

      pure function covSEard_diag(dxi,lxi,theta)

         ! Covariance function

         real(dp) :: covSEard_diag                                ! function value
         real(dp), dimension(:), intent(in) :: dxi,theta   ! input variables
         logical, intent(in) :: lxi

         if( .not. lxi ) then    
              covSEard_diag = theta(2)**2 
         else 
              covSEard_diag = theta(2)**2 * norm2(dxi/theta(3:)**2) 
         endif

      endfunction covSEard_diag  

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Updates covariance matrix and other important temporary arrays in a GP object.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine covariance_matrix(this)

         type(gp), intent(inout) :: this

         integer :: i, j
         real(dp) :: sig2

         sig2 = this%theta(1)**2
         
         selectcase(this%cov_func_string)
            case('ard')   
!$omp parallel do
              do i = 1, this%n
                 do j = 1, this%n
                    this%c(j,i) = covSEard( this%x(:,j), this%x_prime(:,j), this%is_prime(j), &
                    & this%x(:,i), this%x_prime(:,i), this%is_prime(i), this%theta ) 
                 enddo
              enddo
            case('iso')   
!$omp parallel do
              do i = 1, this%n
                 do j = 1, this%n
                    this%c(j,i) = covSEiso( this%x(:,j), this%x_prime(:,j), this%is_prime(j), &
                    & this%x(:,i), this%x_prime(:,i), this%is_prime(i), this%theta )
                 enddo
              enddo
            case default
               call system_abort('gp_predict: '//trim(this%cov_func_string)//' not known')
         endselect

         call apply_llt(this)
         call add_xidentity(this%lcl,sig2)
         !call factorise_symmetric_matrix(this%lcl, this%c_factorised)
         call Matrix_CholFactorise(this%lcl, this%c_factorised)

         call Matrix_BackSubstitute(this%c_factorised,this%alpha,this%y)

      endsubroutine covariance_matrix

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Updates covariance matrix and other important temporary arrays in a GP object
      !% if a teaching point was added. Faster that covariance_matrix, as it uses
      !% the already calculated parts of objects.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine covariance_matrix_extend(this)

         type(gp), intent(inout) :: this

         integer :: i, j
         real(dp) :: sig2
         real(dp), dimension(:), allocatable :: k

         sig2 = this%theta(1)**2
         allocate( k(this%m-1) )
         
         selectcase(this%cov_func_string)
            case('ard')   
!$omp parallel do
              do i = this%l(this%m-1) + 1, this%n
                 do j = 1, i
                    this%c(j,i) = covSEard( this%x(:,j), this%x_prime(:,j), this%is_prime(j), &
                    & this%x(:,i), this%x_prime(:,i), this%is_prime(i), this%theta ) 
                 enddo
              enddo
            case('iso')   
!$omp parallel do
              do i = this%l(this%m-1) + 1, this%n
                 do j = 1, i
                    this%c(j,i) = covSEiso( this%x(:,j), this%x_prime(:,j), this%is_prime(j), &
                    & this%x(:,i), this%x_prime(:,i), this%is_prime(i), this%theta )
                 enddo
              enddo
            case default
               call system_abort('gp_predict: '//trim(this%cov_func_string)//' not known')
         endselect

         do i = 1, this%n
            do j = max( this%l(this%m-1) + 1, i+1) , this%n
               this%c(j,i) = this%c(i,j)
            enddo
         enddo

         call apply_llt(this,is_last=.true.)

         this%lcl(this%m,this%m) = this%lcl(this%m,this%m) + this%theta(1)**2

         call Matrix_Solve_Upper_Triangular(this%c_factorised(:this%m-1,:this%m-1),k,this%lcl(:this%m-1,this%m))

         this%c_factorised(:this%m-1,this%m) = k
         this%c_factorised(this%m,this%m) = sqrt( this%lcl(this%m,this%m) - dot_product(k,k) )

         call Matrix_BackSubstitute(this%c_factorised,this%alpha,this%y)

         deallocate( k )

      endsubroutine covariance_matrix_extend

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% In case of sum of values as teaching points, apply the operation
      !% $ \mathbf{LCL}^\dagger $ or $ \mathbf{Lc} $, where $ \mathbf{L} $
      !% represents the linear combination matrix, here substituted by l integer array as we care about sums only.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine apply_llt_matrix(x,l,lxl,is_last)
         real(dp), dimension(:,:), intent(in)  :: x
         integer, dimension(:), intent(in)     :: l
         real(dp), dimension(:,:), intent(out) :: lxl 
         logical, intent(in), optional :: is_last

         integer :: i, j, ui, uj, li, lj, m, n, first_i
         logical :: my_is_last

         my_is_last = optional_default(.false.,is_last)

         m = size(l)
         if( size(x,1) /= size(x,2) )  call system_abort('apply_llt_matrix: matrix has to be square')
         if( size(lxl,1) /= size(lxl,2) )  call system_abort('apply_llt_matrix: matrix has to be square')
         if( m /= size(lxl,1) ) call system_abort('apply_llt_matrix: l and lxl do not conform')

         n = size(x,1)

         first_i = 1
         if( my_is_last ) first_i = m

         if( m /= n ) then
             do i = first_i, m
                if( i == 1 ) then
                    li = 1
                else
                    li = l(i-1)+1
                endif
                ui = l(i)

                do j = 1, i
                   if( j == 1 ) then
                       lj = 1
                   else
                       lj = l(j-1)+1
                   endif
                   uj = l(j)

                   lxl(j,i) = sum( x(lj:uj,li:ui) )
                enddo
             enddo

             do i = 1, m
                do j = max(i+1,first_i), m
                   lxl(j,i) = lxl(i,j)
                enddo
             enddo
         else
            lxl = x
         endif

      endsubroutine apply_llt_matrix

      subroutine apply_llt_vector(x,l,lx)
         real(dp), dimension(:), intent(in)  :: x
         integer, dimension(:), intent(in)   :: l
         real(dp), dimension(:), intent(out) :: lx 

         integer :: i, m, n

         m = size(l)
         n = size(x)
         if( m /= size(lx) ) call system_abort('apply_llt_matrix: l and lx do not conform')

         if( m /= n ) then
             lx(1) = sum( x(1:l(1)) )
             do i = 2, m
                lx(i) = sum( x(l(i-1)+1:l(i)) )
             enddo
         else
            lx = x
         endif
      endsubroutine apply_llt_vector

      subroutine apply_llt_gp(this,is_last)
         type(gp), intent(inout)       :: this
         logical, intent(in), optional :: is_last

         call apply_llt_matrix(this%c,this%l,this%lcl,is_last)

      endsubroutine apply_llt_gp

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Determines the log likelihood of a GP
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      function likelihood(this)

         type(gp), intent(inout) :: this
         real(dp) :: likelihood
         real(dp) :: determinant

         integer :: i

         determinant = 0.0_dp

         do i = 1, this%m
            determinant = determinant + log( this%c_factorised(i,i) )
         enddo

         likelihood = - determinant &
           & - 0.5_dp * dot_product( this%y, this%alpha )  &
           & - 0.5_dp * this%m * log(2.0_dp * PI)

      endfunction likelihood

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Determines the derivative of the log likelihood of a GP with respect to the hyperparameters
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      function likelihood_gradient(this)

         type(gp), intent(inout) :: this
         real(dp), dimension(this%p) :: likelihood_gradient

         real(dp) :: sig2
         real(dp), dimension(:,:), allocatable :: c_inverse, outer_alpha_minus_c_inverse, c, xc, dc, lcl
         integer                :: i, j, d

         sig2 = this%theta(1)**2

         allocate( c_inverse( this%m, this%m ), outer_alpha_minus_c_inverse( this%m, this%m ), &
         & c( this%n, this%n ), xc( this%n, this%n ), dc( this%n, this%n ) , lcl( this%m, this%m ) )

         c = 0.0_dp
         xc = 0.0_dp
         dc = 0.0_dp

         call Matrix_Factorised_Inverse(this%c_factorised,c_inverse)

!$omp parallel do
         do j = 1, this%m
            do i = 1, this%m
               outer_alpha_minus_c_inverse(i,j) = this%alpha(i) * this%alpha(j) - c_inverse(i,j)
            enddo
         enddo

         lcl = this%lcl
         forall(i = 1:this%m) lcl(i,i) = lcl(i,i) - sig2

         likelihood_gradient = 0.0_dp

         likelihood_gradient(1) = &
         & -this%theta(1) * trace(c_inverse) + this%theta(1) * norm2( this%alpha )

         likelihood_gradient(2) = sum( outer_alpha_minus_c_inverse * lcl ) / this%theta(2)  

         selectcase(this%cov_func_string)
            case('ard')
!$omp parallel do
              do i = 1, this%n
                 do j = 1, this%n
                    c(j,i) = covSEard(this%x(:,j), this%x_prime(:,j), .false., &
                    & this%x(:,i), this%x_prime(:,i), .false., this%theta )
                    if( this%is_prime(i).and.this%is_prime(j) ) xc(j,i) = &
                    & dot_product( (this%x(:,i)-this%x(:,j))/this%theta(3:)**2, this%x_prime(:,j) )
                 enddo
              enddo
            case('iso')
!$omp parallel do
              do i = 1, this%n
                 do j = 1, this%n
                    c(j,i) = covSEiso(this%x(:,j), this%x_prime(:,j), .false., &
                    & this%x(:,i), this%x_prime(:,i), .false., this%theta )
                    if( this%is_prime(i).and.this%is_prime(j) ) xc(j,i) = &
                    & dot_product( (this%x(:,i)-this%x(:,j)), this%x_prime(:,j) )/this%theta(3)
                 enddo
              enddo
         endselect

         selectcase(this%cov_func_string)
            case('ard')
              do d = 1, this%d
!$omp parallel do
                 do i = 1, this%n
                    do j = 1, this%n
                       dc(j,i) = (this%x(d,i)-this%x(d,j))**2 * this%c(j,i)
                       if(this%is_prime(i).and.this%is_prime(j)) then
                          dc(j,i) = dc(j,i) + 2.0_dp * c(j,i) * ( -this%x_prime(d,i) * this%x_prime(d,j) + &
                          & (this%x(d,i) - this%x(d,j)) * ( this%x_prime(d,i)*xc(j,i) - this%x_prime(d,j)*xc(i,j)) )
                       elseif(this%is_prime(i).and..not.this%is_prime(j)) then
                          dc(j,i) = dc(j,i) + 2.0_dp * ( this%x(d,i) - this%x(d,j) ) * this%x_prime(d,i) * c(j,i)
                       elseif(this%is_prime(j).and..not.this%is_prime(i)) then
                          dc(j,i) = dc(j,i) - 2.0_dp * ( this%x(d,j) - this%x(d,i) ) * this%x_prime(d,j) * c(j,i)
                       endif
                    enddo
                 enddo
                 call apply_llt(dc,this%l,lcl)
                 likelihood_gradient(d+2) = 0.5_dp * sum( outer_alpha_minus_c_inverse * &
                 & lcl ) / this%theta(d+2)**3
              enddo
            case('iso')
!$omp parallel do
              do i = 1, this%n
                 do j = 1, this%n
                    dc(j,i) = norm2(this%x(:,i)-this%x(:,j)) * this%c(j,i)
                    if(this%is_prime(i).and.this%is_prime(j)) then
                       dc(j,i) = dc(j,i) + 2.0_dp * c(j,i) * ( - dot_product( this%x_prime(:,i) ,this%x_prime(:,j) ) + &
                       & 4.0_dp * xc(j,i) * xc(i,j) )
                    elseif(this%is_prime(i).xor.this%is_prime(j)) then
                       dc(j,i) = dc(j,i) - 2.0_dp * this%c(j,i)
                    endif
                 enddo
              enddo
              call apply_llt(dc,this%l,lcl)
              likelihood_gradient(d+2) = 0.5_dp * sum( outer_alpha_minus_c_inverse * &
              & lcl ) / this%theta(3)**3
            case default
              call system_abort('likelihood_gradient: '//trim(this%cov_func_string)//' not known')
         endselect

         deallocate( c_inverse, outer_alpha_minus_c_inverse, c, xc, dc , lcl )

      endfunction likelihood_gradient

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

         write(666) this%d, this%n, this%m, this%p, this%n_max_length, this%m_max_length, this%cov_func_string

         do i = 1, this%p
            write(666) this%theta(i)
         enddo

         do i = 1, this%n
            do j = 1, this%d
               write(666) this%x(j,i)
            enddo
         enddo
         do i = 1, this%n
            do j = 1, this%d
               write(666) this%x_prime(j,i)
            enddo
         enddo

         do i = 1, this%n
            do j = 1, this%n
               write(666) this%c(j,i)
            enddo
         enddo
         do i = 1, this%m
            do j = 1, this%m
               write(666) this%c_factorised(j,i)
            enddo
         enddo

         do i = 1, this%m
            do j = 1, this%m
               write(666) this%lcl(j,i)
            enddo
         enddo

         do i = 1, this%m
            write(666) this%y(i)
         enddo
         do i = 1, this%m
            write(666) this%alpha(i)
         enddo
         do i = 1, this%m
            write(666) this%l(i)
         enddo
         do i = 1, this%n
            write(666) this%is_prime(i)
         enddo

         close(unit=666)
         
      endsubroutine gp_print_binary

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Reads GP in binary format from a file.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_read_binary(this,filename)

         type(gp), intent(inout) :: this
         character(len=*), intent(in) :: filename

         integer :: i, j
         logical :: file_exist

         inquire(file=filename,exist=file_exist)

         if( .not. file_exist ) call system_abort('gp_read: '//trim(filename)//' does not exist')

         if( this%initialised ) call finalise(this)

         open(unit=666, file=filename, form='unformatted', action='read')

         read(666) this%d, this%n, this%m, this%p, this%n_max_length, this%m_max_length, this%cov_func_string

         call gp_allocate(this)
         call gp_pointers(this)

         do i = 1, this%p
            read(666) this%theta(i)
         enddo

         do i = 1, this%n
            do j = 1, this%d
               read(666) this%x(j,i)
            enddo
         enddo
         do i = 1, this%n
            do j = 1, this%d
               read(666) this%x_prime(j,i)
            enddo
         enddo

         do i = 1, this%n
            do j = 1, this%n
               read(666) this%c(j,i)
            enddo
         enddo
         do i = 1, this%m
            do j = 1, this%m
               read(666) this%c_factorised(j,i)
            enddo
         enddo

         do i = 1, this%m
            do j = 1, this%m
               read(666) this%lcl(j,i)
            enddo
         enddo

         do i = 1, this%m
            read(666) this%y(i)
         enddo
         do i = 1, this%m
            read(666) this%alpha(i)
         enddo
         do i = 1, this%m
            read(666) this%l(i)
         enddo
         do i = 1, this%n
            read(666) this%is_prime(i)
         enddo

         close(unit=666)
         
         this%initialised = .true.

      endsubroutine gp_read_binary

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Maximises the log likelihood of a GP in the space of the hyperparameters
      !% using nested sampling.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      function minimise_gp_ns(this,theta_min,theta_max,N_live,max_iter,tighten)

         type(gp), intent(inout) :: this                               !% GP
         real(dp), dimension(:), intent(inout) :: theta_min,theta_max  !% range of hyperparameters, dimension(p)
         integer, optional, intent(in)         :: N_live               !% number of living points
         integer, optional, intent(in)         :: max_iter             !% maximum number of iterations
         logical, optional, intent(in)         :: tighten              !% tighten the range of hyperparameters

         logical :: minimise_gp_ns

         real(dp), dimension(:), allocatable :: theta

         if( .not. this%initialised ) call system_abort('gp_minimise_ns: initialise first')
         if( ( size(theta_min) /= this%p ) .or. ( size(theta_max) /= this%p ) ) &
         & call system_abort('gp_minimise_ns: size mismatch')

         allocate( theta(this%p) )
         theta = abs(this%theta)

         minimise_gp_ns = ns(theta,l,theta_min,theta_max,N_live,max_iter,tighten)

         deallocate( theta )

         contains

            function l(theta)

               real(dp), dimension(:), intent(in) :: theta
               real(dp)                           :: l

               call gp_update(this,theta_in=theta)
               l = likelihood(this)

            endfunction l

      endfunction minimise_gp_ns

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Tests the derivative of the log likelihood with test_gradient.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      function test_gp_gradient(this)

         type(gp), intent(inout) :: this                 !% GP
         logical                 :: test_gp_gradient

         test_gp_gradient = test_gradient(this%theta, l, dl)

         contains

            function l(theta, data)

               real(dp), dimension(:) :: theta
               real(dp)               :: l
               character,optional     :: data(:)

               call gp_update(this,theta_in=theta)
               l = likelihood(this)

            endfunction l

            function dl(theta,data)

               real(dp), dimension(:)           :: theta
               real(dp), dimension(size(theta)) :: dl
               character,optional               :: data(:)

               call gp_update(this,theta_in=theta)
               dl = likelihood_gradient(this)

            endfunction dl

      endfunction test_gp_gradient

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Maximises the log likelihood in the space of hyperparameters
      !% using minim (conjugate gradients).
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      function minimise_gp_gradient(this,convergence_tol,max_steps)

         type(gp), intent(inout) :: this                       !% GP
         integer, intent(in), optional :: max_steps            !% maximum number of steps, default: 100
         real(dp), intent(in), optional :: convergence_tol     !% convergence tolerance, default: 0.1
         integer                 :: minimise_gp_gradient      

         integer :: do_max_steps
         real(dp) :: do_convergence_tol

         do_max_steps = optional_default(100,max_steps)
         do_convergence_tol = optional_default(0.1_dp,convergence_tol)

         minimise_gp_gradient = minim(this%theta,l,dl,method='cg',convergence_tol=do_convergence_tol, &
         & max_steps=do_max_steps, hook_print_interval=1)

         contains

            function l(theta, data)

               real(dp), dimension(:) :: theta
               real(dp)               :: l
               character,optional     :: data(:)

               call gp_update(this,theta_in=theta)
               l = - likelihood(this)

            endfunction l

            function dl(theta,data)

               real(dp), dimension(:)           :: theta
               real(dp), dimension(size(theta)) :: dl
               character,optional               :: data(:)

               call gp_update(this,theta_in=theta)
               dl = - likelihood_gradient(this)

            endfunction dl

      endfunction minimise_gp_gradient

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
           endfunction func
        endinterface

        integer               :: iter, i, do_N_live, do_max_iter
        integer, dimension(1) :: ml
        logical :: do_tighten

        real(dp) :: e_old, e_new, e_old_prev
        real(dp), dimension(:), allocatable   :: e, x_temp
        real(dp), dimension(:,:), allocatable :: x_hold

        call print_title('Nested sampling',verbosity=NORMAL)
        call print('Number of variables is '//size(x_inout),verbosity=NORMAL)
        call random_seed()

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

           call print('Iteration: '//iter//', function value dropped: '//e_old,verbosity=NORMAL)
           call print('Point dropped: '//x_temp(:),verbosity=NERD)

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

          endsubroutine get_random_vector

      endfunction ns

endmodule gp_module
