module gp_predict_module

   use libatoms_module
   use fox_wxml
   use FoX_sax, only: xml_t, dictionary_t, haskey, getvalue, parse, &
   open_xml_string, close_xml_t

   implicit none
   private

   type gp

      !% Data arrays. These are (re)allocated quite rarely, always bigger than they should be.
      !% Therefore adding new data points is easy and fast.
      type(LA_Matrix) :: LA_C_nn
      real(dp), dimension(:,:), allocatable :: x, c
      real(dp), dimension(:), allocatable   :: alpha, y
      real(dp), dimension(:,:), allocatable :: x_div_theta

      !% Hyperparameters.
      !% Error parameter, range parameter, sensitivity parameters.
      real(dp), dimension(:,:), allocatable       :: theta
      real(dp) :: sigma
      real(dp), dimension(:), allocatable :: delta, f0

      integer, dimension(:), allocatable :: xz, sp

      integer, dimension(:,:), allocatable :: Z_index

      !% Vector sizes.
      !% d: dimensionality of input space
      !% n: number of teaching points
      !% m: number of teaching function values
      !% p: number of hyperparameters
      integer :: d, n, nsp

      logical :: initialised = .false.

      character(len=1024) :: label = ""
      character(len=10000) :: comment

   end type gp

   integer :: parse_cur_type, parse_cur_sparse_point
   type(extendable_str), private, save :: parse_cur_data
   logical, private :: parse_in_gp, parse_matched_label
   type (gp), pointer :: parse_gp

   interface finalise
      module procedure gp_finalise
   endinterface finalise

   interface gp_read_xml
      module procedure gp_read_xml_t_xml, gp_read_params_xml
   endinterface gp_read_xml

   public :: gp, finalise, gp_predict, gp_precompute_covariance, &
      gp_print_xml, gp_read_xml, setup_x_div_theta, GP_sort_data_by_Z

   contains

      subroutine GP_finalise(this)

         type(gp), intent(inout)              :: this !% gp object

         if( .not. this%initialised ) return

         if(allocated(this%theta)) deallocate(this%theta)
         if(allocated(this%delta)) deallocate(this%delta)
         if(allocated(this%x)) deallocate(this%x)
         if(allocated(this%x_div_theta)) deallocate(this%x_div_theta)
         if(allocated(this%c)) deallocate(this%c)
         if(allocated(this%alpha)) deallocate(this%alpha)
         if(allocated(this%y)) deallocate(this%y)
         if(allocated(this%xz)) deallocate(this%xz)
         if(allocated(this%sp)) deallocate(this%sp)
         if(allocated(this%f0)) deallocate(this%f0)

         if(allocated(this%Z_index)) deallocate(this%Z_index)

         call finalise(this%LA_C_nn)
            
         this%d = 0
         this%n = 0
         this%sigma = 0.0_dp
         this%initialised = .false.

      end subroutine GP_finalise


      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !
      !% Predict the function value and its variance using GP.
      !
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine gp_precompute_covariance(gp_data,x_star,c,Z,mpi)
         type(gp), intent(in)                        :: gp_data
         real(dp), dimension(:,:), intent(in)        :: x_star                 !% test point, dimension(d)
         real(dp), dimension(:,:), intent(out)       :: c
         integer, dimension(:), intent(in), optional :: Z
         type(mpi_context), intent(in), optional     :: mpi

         integer :: i, j, Z_type, sp, m
         real(dp), dimension(gp_data%d,gp_data%nsp) :: inv_theta
         real(dp), dimension(gp_data%d) :: x_star_div_theta, xixjtheta
         integer, dimension(:), allocatable :: my_Z

         allocate(my_Z(size(x_star,2)))
         if( size(x_star,1) /= gp_data%d ) call system_abort('gp_precompute_covariance: first dimension of x_star is '//size(x_star,1)//', not '//gp_data%d)
         m = size(x_star,2)

         if(present(Z)) then
            call check_size('Z',Z,(/m/), 'gp_precompute_covariance')
            my_Z = Z
         else
            my_Z = 0
         endif

         call check_size('c',c,(/gp_data%n,m/), 'gp_precompute_covariance')

         inv_theta = 1.0_dp / gp_data%theta
         c = 0.0_dp

!$omp do private(sp, Z_type,x_star_div_theta,j,xixjtheta)
         do i = 1, m
            if(present(mpi)) then
               if (mpi%active) then
                  if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
               endif
            endif
            ! Determine what type of atom we have
            do sp = 1, gp_data%nsp; if( gp_data%sp(sp) == my_Z(i) ) Z_type = sp; enddo

            x_star_div_theta = x_star(:,i) * inv_theta(:,Z_type)
            do j = 1, gp_data%n
               if( my_Z(i) == gp_data%xz(j) ) then
                  xixjtheta = gp_data%x_div_theta(:,j)-x_star_div_theta
                  c(j,i) = gp_data%delta(Z_type)**2 * exp( - 0.5_dp * dot_product(xixjtheta,xixjtheta) )
               endif
            enddo
         enddo
!$omp end do 

         if(present(mpi)) call sum_in_place(mpi,c)

         deallocate(my_Z)
      endsubroutine gp_precompute_covariance

      subroutine gp_predict(gp_data, mean, variance, x_star, x_prime_star, Z, c_in, xixjtheta_in, new_x_star, use_dgemv)

         type(gp), intent(in)               :: gp_data                !% GP
         real(dp), intent(out), optional    :: mean, variance         !% output, predicted value and variance at test point
         real(dp), dimension(:), intent(in) :: x_star                 !% test point, dimension(d)
         real(dp), dimension(:), intent(in), optional :: x_prime_star !% test point, dimension(d)
         integer, intent(in), optional :: Z
         real(dp), dimension(:), intent(in), optional :: c_in 
#ifdef GP_PREDICT_QP
         real(qp), dimension(:,:), intent(inout), target, optional :: xixjtheta_in
#else
         real(dp), dimension(:,:), intent(inout), target, optional :: xixjtheta_in
#endif
         logical, intent(in), optional :: new_x_star, use_dgemv

#ifdef GP_PREDICT_QP
         real(qp), dimension(:), allocatable :: k, c
         real(qp), dimension(:,:), pointer :: xixjtheta !, tmp
	 real(qp) :: x_prime_star_div_theta(size(x_star)), x_star_div_theta(size(x_star))
#else
         real(dp), dimension(:), allocatable :: k, c
         real(dp), dimension(:,:), pointer :: xixjtheta !, tmp
	 real(dp) :: x_prime_star_div_theta(size(x_star)), x_star_div_theta(size(x_star))
#endif
         real(dp) :: kappa

	 logical :: my_new_x_star, my_use_dgemv
         integer :: i, do_Z, Z_type

	 allocate(k(gp_data%n))
	 allocate(c(gp_data%n))
	 if (present(xixjtheta_in)) then
	    if (size(xixjtheta_in,1) /= gp_data%d) call system_abort("BAD")
	    if (size(xixjtheta_in,2) /= gp_data%n) call system_abort("BAD")
	    xixjtheta => xixjtheta_in
	 else
	    allocate(xixjtheta(gp_data%d, gp_data%n))
	 endif

	 my_new_x_star = optional_default(.true., new_x_star)
	 my_use_dgemv = optional_default(.false., use_dgemv)

         !if( .not. gp_data%initialised ) &
         !& call system_abort('gp_predict: not initialised, call gp_initialise first')

         !if( size(x_star) /= gp_data%d ) call system_abort('gp_predict: array sizes do not conform')

         do_Z = optional_default(gp_data%xz(1),Z)

         do i = 1, gp_data%nsp; if( gp_data%sp(i) == do_Z ) Z_type = i; end do
	 x_star_div_theta = x_star / gp_data%theta(:,Z_type)
         
	 ! not needed - all values to be used will be set before they're used
         ! xixjtheta = 0.0_dp 
         c = 0.0_dp
         do i = 1, gp_data%n
            if( do_Z == gp_data%xz(i) ) then
               if (my_new_x_star) xixjtheta(:,i) = gp_data%x_div_theta(:,i)-x_star_div_theta
               if(.not. present(c_in)) c(i) = gp_data%delta(Z_type)**2 * exp( - 0.5_dp * dot_product(xixjtheta(:,i),xixjtheta(:,i)) )
            endif
         enddo
         if(present(c_in)) c = real(c_in, qp)

         if(present(x_prime_star)) then
	    x_prime_star_div_theta = x_prime_star / gp_data%theta(:,Z_type)
	    if (my_use_dgemv) then
	       k = 0.0_dp
	       call dgemv('T',gp_data%d,gp_data%Z_index(2,do_Z)-gp_data%Z_index(1,do_Z)+1,1.0_dp, &
		  xixjtheta(1,gp_data%Z_index(1,do_Z)), gp_data%d, x_prime_star_div_theta, 1, 0.0_dp, k(gp_data%Z_index(1,do_Z)), 1)
	       k(gp_data%Z_index(1,do_Z):gp_data%Z_index(2,do_Z)) = k(gp_data%Z_index(1,do_Z):gp_data%Z_index(2,do_Z)) * &
								    c(gp_data%Z_index(1,do_Z):gp_data%Z_index(2,do_Z))
	    else
	       do i = 1, gp_data%n
		  !k(i) = covSEard( gp_data%delta, gp_data%theta, gp_data%x(:,i), x_star, x_prime_star )
		  if( do_Z == gp_data%xz(i) ) then
#ifdef GP_PREDICT_QP
		     !xixjtheta = real((gp_data%x(:,i)-x_star)/gp_data%theta(:,Z_type),qp)
		     !k(i) = gp_data%delta(Z_type)**2 * exp( - 0.5_qp * dot_product(xixjtheta,xixjtheta) ) * &
		     !dot_product(xixjtheta,real(x_prime_star/gp_data%theta(:,Z_type),qp))
#else
		     !xixjtheta = (gp_data%x_div_theta(:,i)-x_star_div_theta)
		     !k(i) = gp_data%delta(Z_type)**2 * exp( - 0.5_dp * dot_product(xixjtheta,xixjtheta) ) * &
		     !& dot_product(xixjtheta,x_prime_star_div_theta)
#endif
		     k(i) = c(i) * dot_product(xixjtheta(:,i),x_prime_star_div_theta)
		  else
		     k(i) = 0.0_dp
		  end if
		  
		  !k(i) = exp( - 0.5_qp * dot_product(xixjtheta,xixjtheta) ) * &
		  !& dot_product(xixjtheta,real(x_prime_star,qp)/real(gp_data%theta,qp))
		  !xixjtheta = real((real(gp_data%x(:,i),qp)-real(x_star,qp))/real(gp_data%theta,qp),dp)

		  !k(i) = gp_data%delta**2 * exp( - 0.5_dp * normsq(xixjtheta) ) * dot_product(xixjtheta,x_prime_star*tmp)
		  !xixjtheta = (real(gp_data%x(:,i),qp)-real(x_star,qp))/real(gp_data%theta,qp)
	       end do
	    endif ! use_dgemv
         else
            do i = 1, gp_data%n
!               xixjtheta = (gp_data%x(:,i)-x_star)*tmp
               if( do_Z == gp_data%xz(i) ) then
#ifdef GP_PREDICT_QP
                  !xixjtheta = real((gp_data%x(:,i)-x_star)/gp_data%theta(:,Z_type),qp)
                  !k(i) = gp_data%delta(Z_type)**2 * exp( - 0.5_qp * dot_product(xixjtheta,xixjtheta) ) + gp_data%f0(Z_type)**2
#else
                  !xixjtheta = (gp_data%x_div_theta(:,i)-x_star_div_theta)
                  !k(i) = gp_data%delta(Z_type)**2 * exp( - 0.5_dp * dot_product(xixjtheta,xixjtheta) ) + gp_data%f0(Z_type)**2
#endif
                  k(i) = c(i) + gp_data%f0(Z_type)**2
               else
                  k(i) = 0.0_dp
               end if
               !xixjtheta = (real(gp_data%x(:,i),qp)-real(x_star,qp))/real(gp_data%theta,qp)
               !k(i) = exp( - 0.5_qp * dot_product(xixjtheta,xixjtheta) )
            end do
         end if
#ifdef GP_PREDICT_QP
         if( present(mean) ) mean = dot_product( k, gp_data%alpha )
#else
         if( present(mean) ) mean = real(dot_product( k, real(gp_data%alpha,qp) ), dp)
#endif
 !        if( present(variance) ) then
 !           if( present(x_prime_star) ) then
 !              kappa = gp_data%sigma(2)**2 + gp_data%delta**2 * normsq(x_prime_star*gp_data%theta)
 !           else
 !              kappa = gp_data%sigma(1)**2 + gp_data%delta**2
 !           end if
 !           variance = sqrt( kappa - dot_product(k,matmul(gp_data%c,k)) )
 !        end if
         !& variance = sqrt( gp_data%delta + gp_data%sigma**2 - dot_product(k,matmul(gp_data%c,k)) )

!         deallocate( k )

         deallocate(k, c)
         if (.not. present(xixjtheta_in)) deallocate(xixjtheta)
      end subroutine gp_predict

      subroutine GP_sort_data_by_Z(this)
	 type(gp), intent(inout) :: this

	 integer :: i
	 integer, allocatable :: perm(:)
	 real(dp), allocatable :: xz(:)

	 allocate(perm(this%n))
	 allocate(xz(this%n))

	 do i=1, this%n
	    perm(i) = i
	 end do

	 ! do sort
	 xz = this%xz
	 call sort_array(xz, perm)

	 ! apply permutation to get sorted arrays
	 this%xz(:) = this%xz(perm(:))
	 this%alpha(:) = this%alpha(perm(:))
	 do i=1, this%d
	    this%x(i,:) = this%x(i,perm(:))
	 end do

	 ! find beginning and end indices of each xz value
	 allocate(this%Z_index(2,maxval(this%xz)))
	 this%Z_index = 0
	 this%Z_index(1, this%xz(1)) = 1
	 do i=2, this%n
	    if (this%xz(i) < this%xz(i-1)) then ! improperly sorted - error
	       call system_abort("GP_sort_data_by_Z got unsorted xz array for indexing")
	    endif
	    if (this%xz(i) > this%xz(i-1)) then ! sorted, new Z
	       this%Z_index(2,this%xz(i-1)) = i-1 ! last Z ended on last i
	       this%Z_index(1,this%xz(i)) = i ! this Z starts on this i
	    endif
	 end do
	 this%Z_index(2,this%xz(this%n)) = this%n

	 deallocate(perm)
	 deallocate(xz)

      end subroutine GP_sort_data_by_Z

      subroutine gp_print_xml(this,xf,label)

         type(gp), intent(in) :: this
         type(xmlf_t), intent(inout) :: xf
         character(len=*), intent(in), optional :: label

         integer :: i

         if( .not. this%initialised ) call system_abort('gp_print_xml: gp not initialised')

         call xml_NewElement(xf,"GP_data")
         call xml_AddAttribute(xf,"dimensions", ""//this%d)
         call xml_AddAttribute(xf,"n_species", ""//this%nsp)
         call xml_AddAttribute(xf,"n_sparse_x", ""//this%n)
         if(present(label)) call xml_AddAttribute(xf,"label", trim(label))

         do i = 1, this%nsp
            call xml_NewElement(xf,"per_species_data")
            call xml_AddAttribute(xf,"type", ""//i)
            call xml_AddAttribute(xf,"signal_variance", ""//this%delta(i))
            call xml_AddAttribute(xf,"atomic_num",""//this%sp(i))
            call xml_AddAttribute(xf,"function_offset",""//this%f0(i))
            call xml_NewElement(xf,"theta")
            call xml_AddCharacters(xf, ""//this%theta(:,i)//" ")
            call xml_EndElement(xf,"theta")
            call xml_EndElement(xf,"per_species_data")
         enddo

         do i = 1, this%n
            call xml_NewElement(xf,"sparse_x")
            call xml_AddAttribute(xf,"i", ""//i)
            call xml_AddAttribute(xf,"coefficient", ""//this%alpha(i))
            call xml_AddAttribute(xf,"sparse_x_atomic_number", ""//this%xz(i))
            call xml_AddCharacters(xf, ""//this%x(:,i)//" ")
            call xml_EndElement(xf,"sparse_x")
         enddo

         call xml_EndElement(xf,"GP_data")
         
      end subroutine gp_print_xml

      subroutine gp_read_xml_t_xml(this,xp,label)
         type(gp), intent(inout), target :: this
         type(xml_t), intent(inout) :: xp
         character(len=*), intent(in), optional :: label

         if( this%initialised ) call finalise(this)

         parse_in_gp = .false.
         parse_matched_label = .false.
         parse_gp => this
         call initialise(parse_cur_data)

         parse_gp%label = optional_default("",label)

         call parse(xp, &
         characters_handler = GP_characters_handler, &
         startElement_handler = GP_startElement_handler, &
         endElement_handler = GP_endElement_handler)

         call finalise(parse_cur_data)

	 call setup_x_div_theta(this)
         this%initialised = .true.

      endsubroutine gp_read_xml_t_xml

      subroutine gp_read_params_xml(this,params_str,label)

         type(gp), intent(inout), target :: this
         character(len=*), intent(in) :: params_str
         character(len=*), intent(in), optional :: label

         type(xml_t) :: xp

         call open_xml_string(xp, params_str)
         call gp_read_xml(this,xp,label)
         call close_xml_t(xp)

      end subroutine gp_read_params_xml

      subroutine GP_startElement_handler(URI, localname, name, attributes)
         character(len=*), intent(in)   :: URI
         character(len=*), intent(in)   :: localname
         character(len=*), intent(in)   :: name
         type(dictionary_t), intent(in) :: attributes

         integer             :: status, ti, xi
         character(len=1024) :: value

         if(name == 'GP_data') then ! new GP_data
            if(parse_in_gp) &
            call system_abort("GP_startElement_handler entered GP_data with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

            if(parse_matched_label) return ! we already found an exact match for this label

            call GP_FoX_get_value(attributes, 'label', value, status)
            if (status /= 0) value = ''

            if(len(trim(parse_gp%label)) > 0) then ! we were passed in a label
               if(trim(value) == trim(parse_gp%label)) then
                  parse_matched_label = .true.
                  parse_in_gp = .true.
               else ! no match
                  parse_in_gp = .false.
               endif
            else ! no label passed in
               parse_in_gp = .true.
            endif

            if(parse_in_gp) then
               if(parse_gp%initialised) call finalise(parse_gp)

               call GP_FoX_get_value(attributes, 'dimensions', value, status)
               if (status == 0) then
                  read (value,*) parse_gp%d
               else
                  call system_abort("GP_startElement_handler did not find the dimensions attribute.")
               endif

               call GP_FoX_get_value(attributes, 'n_species', value, status)
               if (status == 0) then
                  read (value,*) parse_gp%nsp
               else
                  call system_abort("GP_startElement_handler did not find the n_species attribute.")
               endif

               call GP_FoX_get_value(attributes, 'n_sparse_x', value, status)
               if (status == 0) then
                  read (value,*) parse_gp%n
               else
                  call system_abort("GP_startElement_handler did not find the n_sparse_x attribute.")
               endif

               allocate(parse_gp%theta(parse_gp%d,parse_gp%nsp), parse_gp%delta(parse_gp%nsp), parse_gp%x(parse_gp%d,parse_gp%n), &
               & parse_gp%alpha(parse_gp%n), parse_gp%xz(parse_gp%n), parse_gp%sp(parse_gp%nsp), parse_gp%f0(parse_gp%nsp))
            endif

         elseif(parse_in_gp .and. name == 'per_species_data') then
            call GP_FoX_get_value(attributes, 'type', value, status)
            if (status == 0) then
               read (value,*) ti
               if(ti > parse_gp%nsp) call system_abort("GP_startElement_handler got invalid type "//ti)
            else
               call system_abort("GP_startElement_handler did not find the type attribute.")
            endif

            call GP_FoX_get_value(attributes, 'signal_variance', value, status)
            if (status == 0) then
               read (value,*) parse_gp%delta(ti)
            else
               call system_abort("GP_startElement_handler did not find the signal_variance attribute.")
            endif

            call GP_FoX_get_value(attributes, 'atomic_num', value, status)
            if (status == 0) then
               read (value,*) parse_gp%sp(ti)
            else
               call system_abort("GP_startElement_handler did not find the atomic_num attribute.")
            endif

            call GP_FoX_get_value(attributes, 'function_offset', value, status)
            if (status == 0) then
               read (value,*) parse_gp%f0(ti)
            else
               call system_abort("GP_startElement_handler did not find the function_offset attribute.")
            endif

            parse_cur_type = ti

         elseif(parse_in_gp .and. name == 'sparse_x') then

            call GP_FoX_get_value(attributes, 'i', value, status)
            if (status == 0) then
               read (value,*) xi
               if(xi > parse_gp%n) call system_abort("GP_startElement_handler got sparse point out of range ("//xi//")")
            else
               call system_abort("GP_startElement_handler did not find the i attribute.")
            endif

            call GP_FoX_get_value(attributes, 'coefficient', value, status)
            if (status == 0) then
               read (value,*) parse_gp%alpha(xi)
            else
               call system_abort("GP_startElement_handler did not find the coefficient attribute.")
            endif

            call GP_FoX_get_value(attributes, 'sparse_x_atomic_number', value, status)
            if (status == 0) then
               read (value,*) parse_gp%xz(xi)
            else
               call system_abort("GP_startElement_handler did not find the sparse_x_atomic_number attribute.")
            endif

            parse_cur_sparse_point = xi

            call zero(parse_cur_data)

         elseif(parse_in_gp .and. name == 'theta') then
            call zero(parse_cur_data)
         endif

      endsubroutine GP_startElement_handler

      subroutine GP_endElement_handler(URI, localname, name)
         character(len=*), intent(in)   :: URI
         character(len=*), intent(in)   :: localname
         character(len=*), intent(in)   :: name
      
         character(len=100*parse_gp%d) :: val
      
         if(parse_in_gp) then
            if(name == 'GP_data') then
               parse_in_gp = .false.
            elseif(name == 'sparse_x') then
               val = string(parse_cur_data)
               read(val,*) parse_gp%x(:,parse_cur_sparse_point)
            elseif(name == 'theta') then
               val = string(parse_cur_data)
               read(val,*) parse_gp%theta(:,parse_cur_type)
            endif
         endif
            
      endsubroutine GP_endElement_handler

      subroutine GP_characters_handler(in)
         character(len=*), intent(in) :: in
      
         if(parse_in_gp) then
           call concat(parse_cur_data, in, keep_lf=.false.)
         endif
      end subroutine GP_characters_handler

      subroutine setup_x_div_theta(this)
	 type(gp), intent(inout) :: this

	 integer :: i, j, Z_type

	 if (allocated(this%x_div_theta)) deallocate(this%x_div_theta)
	 allocate(this%x_div_theta(this%d,this%n))
	 do i = 1, this%n

            do j = 1, this%nsp
               if( this%sp(j) == this%xz(i) ) Z_type = j
            enddo

            this%x_div_theta(:,i) = this%x(:,i) / this%theta(:,Z_type)
	 enddo
      end subroutine

      subroutine GP_FoX_get_value(attributes, key, val, status)
        type(dictionary_t), intent(in) :: attributes
        character(len=*), intent(in) :: key
        character(len=*), intent(inout) :: val
        integer, intent(out), optional :: status
      
        if (HasKey(attributes,key)) then
          val = GetValue(attributes, trim(key))
          if (present(status)) status = 0
        else
          val = ""
          if (present(status)) status = 1
        endif
      end subroutine GP_FoX_get_value

end module gp_predict_module
