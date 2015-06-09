module real_space_covariance_module

   use libatoms_module

   implicit none

   private

   !% 24-point quaternion grid from centers of the cells of a
   !% truncated-cubic tetracontaoctachoron (48-cell), coverage angle
   !% alpha=62.8 degrees. See http://charles.karney.info/orientation
   real(dp), dimension(4,24) :: quat_grid_c48u1 = reshape( &
      (/1.        ,  0.        ,  0.        ,  0.        , &
        0.        ,  1.        ,  0.        ,  0.        , &
        0.        ,  0.        ,  1.        ,  0.        , &
        0.        ,  0.        ,  0.        ,  1.        , &
        0.5       ,  0.5       ,  0.5       ,  0.5       , &
        0.5       ,  0.5       ,  0.5       , -0.5       , &
        0.5       ,  0.5       , -0.5       ,  0.5       , &
        0.5       ,  0.5       , -0.5       , -0.5       , &
        0.5       , -0.5       ,  0.5       ,  0.5       , &
        0.5       , -0.5       ,  0.5       , -0.5       , &
        0.5       , -0.5       , -0.5       ,  0.5       , &
        0.5       , -0.5       , -0.5       , -0.5       , &
        0.70710678,  0.70710678,  0.        ,  0.        , &
        0.70710678, -0.70710678,  0.        ,  0.        , &
        0.70710678,  0.        ,  0.70710678,  0.        , &
        0.70710678,  0.        , -0.70710678,  0.        , &
        0.70710678,  0.        ,  0.        ,  0.70710678, &
        0.70710678,  0.        ,  0.        , -0.70710678, &
        0.        ,  0.70710678,  0.70710678,  0.        , &
        0.        ,  0.70710678, -0.70710678,  0.        , &
        0.        ,  0.70710678,  0.        ,  0.70710678, &
        0.        ,  0.70710678,  0.        , -0.70710678, &
        0.        ,  0.        ,  0.70710678,  0.70710678, &
        0.        ,  0.        ,  0.70710678, -0.70710678 /), (/4, 24/))

   public RealSpaceCovariance
   type RealSpaceCovariance
      integer n, n_grid_minim
      type(Atoms), dimension(:), allocatable :: data
      real(dp) :: sigma_alpha, sigma_beta, sigma_covariance, sigma_error, cutoff, transition_width, covariance_tol
      real(dp), allocatable, dimension(:) :: k
      real(dp), allocatable, dimension(:,:) :: covariance_matrix, cov_inv, q_min, rot_force, &
           distance_matrix, force_covariance_matrix, force_difference_matrix
      real(dp), allocatable, dimension(:,:,:) :: alignment
   end type RealSpaceCovariance

   type minimise_overlap
      integer tag1, tag2
      type(atoms), pointer :: at1, at2
      type(RealSpaceCovariance), pointer :: rscov
   endtype minimise_overlap

   type optimise_likelihood
      logical use_qr
      type(RealSpaceCovariance), pointer :: rscov
   endtype optimise_likelihood

   public initialise
   interface initialise
      module procedure realspacecovariance_initialise
   end interface initialise

   public finalise
   interface finalise
      module procedure realspacecovariance_finalise
   end interface finalise

   public teach
   interface teach
      module procedure realspacecovariance_teach
   end interface teach

   public predict
   interface predict
      module procedure realspacecovariance_predict
   end interface predict

   public predict_2
   interface predict_2
      module procedure realspacecovariance_predict_2
   end interface predict_2

   public predict_3
   interface predict_3
      module procedure realspacecovariance_predict_3
   end interface predict_3

   public covariance
   interface covariance
      module procedure realspacecovariance_covariance
   end interface covariance

   public align
   interface align
      module procedure realspacecovariance_align_pair
      module procedure realspacecovariance_align_with_data_set
      module procedure realspacecovariance_align_data_sets
   end interface align

   public set_n
   interface set_n
      module procedure realspacecovariance_set_n
   end interface set_n

   public centre_on_atom_1


   public minimise_overlap
   public overlap_func, doverlap_func
   public random_quaternion, rotmat_quaternion, drotmat_dquaternion

   contains

   subroutine realspacecovariance_finalise(this)
     type(realspacecovariance), intent(inout) :: this
     integer i

     if (allocated(this%data)) then
        do i=1,size(this%data)
           call finalise(this%data(i))
        end do
        deallocate(this%data) 
     end if
     this%n = 0
     if (allocated(this%covariance_matrix)) deallocate(this%covariance_matrix)
     if (allocated(this%cov_inv)) deallocate(this%cov_inv)
     if (allocated(this%distance_matrix)) deallocate(this%distance_matrix)
     if (allocated(this%force_covariance_matrix)) deallocate(this%force_covariance_matrix)
     if (allocated(this%force_difference_matrix)) deallocate(this%force_difference_matrix)
     if (allocated(this%q_min)) deallocate(this%q_min)
     if (allocated(this%rot_force)) deallocate(this%rot_force)
     if (allocated(this%k)) deallocate(this%k)
     if (allocated(this%alignment)) deallocate(this%alignment)

   end subroutine realspacecovariance_finalise

   subroutine realspacecovariance_initialise(this, n, n_grid_minim, sigma_alpha, sigma_beta, &
        sigma_covariance, sigma_error, cutoff_in, transition_width, covariance_tol)
     type(realspacecovariance), intent(inout) :: this
     integer, intent(in) :: n, n_grid_minim
     real(dp), intent(in) :: sigma_alpha, sigma_beta, sigma_covariance, sigma_error, cutoff_in, transition_width, covariance_tol

     integer i, j
     real(dp) tmp_lattice(3,3)

     call finalise(this)

     this%n = n
     this%n_grid_minim = n_grid_minim
     this%sigma_alpha = sigma_alpha
     this%sigma_beta = sigma_beta
     this%sigma_covariance = sigma_covariance
     this%sigma_error = sigma_error
     this%cutoff = cutoff_in
     this%transition_width = transition_width
     this%covariance_tol = covariance_tol

     allocate(this%data(n))
     allocate(this%covariance_matrix(size(this%data),size(this%data)), this%cov_inv(size(this%data),size(this%data)))
     allocate(this%distance_matrix(size(this%data),size(this%data)), this%force_covariance_matrix(size(this%data),size(this%data)) )
     allocate(this%force_difference_matrix(size(this%data),size(this%data)) )
     allocate(this%q_min(4,size(this%data)))
     allocate(this%rot_force(3,size(this%data)))
     allocate(this%k(size(this%data)))
     allocate(this%alignment(4,size(this%data),size(this%data)))

     tmp_lattice(:,:) = 0.0_dp
     do i=1,size(this%data)
        call initialise(this%data(i), 0, tmp_lattice)
        this%q_min(:,i) = (/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /)
        do j=1,size(this%data)
           this%alignment(:,i,j) = (/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /)
        end do
     end do

   end subroutine realspacecovariance_initialise

   subroutine realspacecovariance_set_n(this, n)
     type(realspacecovariance), intent(inout) :: this
     integer, intent(in) :: n
     
     call finalise(this)
     call initialise(this, n, this%n_grid_minim, this%sigma_alpha, this%sigma_beta, this%sigma_covariance, &
          this%sigma_error, this%cutoff, this%transition_width, this%covariance_tol)

   end subroutine realspacecovariance_set_n


   ! Return uniform random unit quaternion. From
   ! K. Shoesmith, Uniform Random Rotations, in D. Kirk (ed.)
   ! Graphics Gems III, pages 124-132, Academic Press (1992).
   function random_quaternion()
     real(dp), dimension(4) :: random_quaternion

     real(dp) :: rand0, rand1, rand2, r1, r2, t1, t2

     rand0 = ran_uniform()
     rand1 = ran_uniform()
     rand2 = ran_uniform()

     r1 = sqrt(1.0_dp - rand0)
     r2 = sqrt(rand0)
     t1 = PI*2.0_dp * rand1
     t2 = PI*2.0_dp * rand2
     random_quaternion = (/cos(t2)*r2, sin(t1)*r1, &
                           cos(t1)*r1, sin(t2)*r2/)

   end function random_quaternion

   subroutine overlap_rho_rot_rho(this,at1,at2,q,overlap,doverlap)
      type(RealSpaceCovariance), intent(in) :: this
      type(Atoms), intent(in) :: at1, at2
      real(dp), dimension(4), intent(in), optional :: q
      real(dp), dimension(4), intent(out), optional :: doverlap
      real(dp), intent(out), optional :: overlap

      real(dp) :: diff_ij(3), rot(3,3), drot(3,3,4), overlap_ij, r_i, r_j, fcut_i, fcut_j, sigma_i, sigma_j
      real(dp), dimension(:,:), allocatable :: pos
      real(dp), dimension(:,:,:), allocatable :: dpos
      integer :: i, j, k
 
      if(.not. present(overlap) .and. .not.present(doverlap)) return

      allocate(pos(3,at2%N))

      if (present(q)) then
         rot = rotmat_quaternion(q)
         pos = matmul(rot,at2%pos)
      else
         pos = at2%pos
      end if

      if(present(overlap)) overlap = 0.0_dp

      if(present(doverlap)) then
         allocate(dpos(3,at2%N,4))
         if (present(q)) then
            drot = drotmat_dquaternion(q)
         else
            drot = 0.0_dp
         end if

         do i = 1, 4
            dpos(:,:,i) = matmul(drot(:,:,i),at2%pos)
         enddo
         doverlap = 0.0_dp
      endif

      do i = 2, at1%N
         r_i = norm(at1%pos(:,i))
         if (r_i > this%cutoff) cycle
         fcut_i = coordination_function(r_i, this%cutoff, this%transition_width)
         sigma_i = this%sigma_alpha*r_i**2 + this%sigma_beta

         do j = 2, at2%N
            r_j = norm(at2%pos(:,j))
            if (r_j > this%cutoff) cycle
            fcut_j = coordination_function(r_j, this%cutoff, this%transition_width)
            sigma_j = this%sigma_alpha*r_j**2 + this%sigma_beta

            diff_ij = at1%pos(:,i)-pos(:,j)
            overlap_ij = sigma_i*sigma_j/sqrt(sigma_i**2 + sigma_j**2) * &
                 exp(-sum( diff_ij**2 ) /(2.0_dp*(sigma_i**2 + sigma_j**2)))*fcut_i*fcut_j

            !write (*,*) i, j, r_i, r_j, sigma_i, sigma_j, overlap_ij

            if(present(overlap)) overlap = overlap + overlap_ij

            if(present(doverlap)) then
               do k = 1, 4
                  doverlap(k) = doverlap(k) + dot_product(diff_ij,dpos(:,j,k))*overlap_ij / (sigma_i**2 + sigma_j**2)
               enddo
            endif

         enddo
      enddo
 
      if(allocated(dpos)) deallocate(dpos)
      if(allocated(pos)) deallocate(pos)
 
   endsubroutine overlap_rho_rot_rho

   function overlap_func(x,am_data)
      real(dp), dimension(:) :: x
      character(len=1), dimension(:), optional :: am_data

      real(dp) :: overlap_func
      type(minimise_overlap) :: am

      am = transfer(am_data,am)
      call overlap_rho_rot_rho(am%rscov,am%at1,am%at2,x,overlap=overlap_func)

      overlap_func = -overlap_func
   endfunction overlap_func

   function doverlap_func(x,am_data)
      real(dp), dimension(:) :: x
      character(len=1), dimension(:), optional :: am_data

      real(dp), dimension(size(x)) :: doverlap_func
      type(minimise_overlap) :: am

      am = transfer(am_data,am)
      call overlap_rho_rot_rho(am%rscov,am%at1,am%at2,x,doverlap=doverlap_func)

      doverlap_func = -doverlap_func
   endfunction doverlap_func

   function sum_overlap_func(x,am_data)
      real(dp), dimension(:) :: x
      character(len=1), dimension(:), optional :: am_data

      integer i, j
      integer tag_i, tag_j
      real(dp) :: time_i, time_j
      real(dp) :: overlap, sum_overlap_func, weight
      type(minimise_overlap) :: am

      am = transfer(am_data,am)

      if (am%tag2 == 0) then
         if (.not. get_value(am%at1%params, 'time', time_j)) then
            call system_abort('cannot find time in at1')
         end if
      end if

      sum_overlap_func = 0.0_dp
      do i=1,size(am%rscov%data)
         if (.not. get_value(am%rscov%data(i)%params, 'tag', tag_i)) then
            call system_abort('cannot find tag in rscov%data('//i//')')
         end if
         if (.not. get_value(am%rscov%data(i)%params, 'time', time_i)) then
            call system_abort('cannot find time in rscov%data('//i//')')
         end if
         if (tag_i /= am%tag1) cycle

         if (am%tag2 /= 0) then
            do j=1,size(am%rscov%data)
               if(.not. get_value(am%rscov%data(j)%params, 'tag', tag_j)) then
                  call system_abort('cannot find tag in rscov%data('//j//')')
               end if
               if (tag_j /= am%tag2) cycle

               call overlap_rho_rot_rho(am%rscov,am%rscov%data(i),am%rscov%data(j),x,overlap=overlap)               
               weight = 1.0_dp !exp(-abs(time_i-time_j)/am%rscov%tau)
               sum_overlap_func = sum_overlap_func + weight*overlap
            end do
         else
            call overlap_rho_rot_rho(am%rscov,am%rscov%data(i),am%at1,x,overlap=overlap)
            weight = 1.0_dp !exp(-abs(time_i-time_j)/am%rscov%tau)
            sum_overlap_func = sum_overlap_func + weight*overlap
         end if
      end do
      sum_overlap_func = -sum_overlap_func

   end function sum_overlap_func

   function sum_doverlap_func(x,am_data)
      real(dp), dimension(:) :: x
      character(len=1), dimension(:), optional :: am_data

      integer i, j
      integer :: tag_i, tag_j
      real(dp) :: time_i, time_j, weight
      real(dp) :: doverlap(size(x)), sum_doverlap_func(size(x))
      type(minimise_overlap) :: am

      am = transfer(am_data,am)

      if (am%tag2 == 0) then
         if (.not. get_value(am%at1%params, 'time', time_j)) then
            call system_abort('cannot find time in at1')
         end if
      end if

      sum_doverlap_func = 0.0_dp
      do i=1,size(am%rscov%data)
         if(.not. get_value(am%rscov%data(i)%params, 'tag', tag_i)) then
            call system_abort('cannot find tag in rscov%data('//i//')')
         end if
         if (.not. get_value(am%rscov%data(i)%params, 'time', time_i)) then
            call system_abort('cannot find time in rscov%data('//i//')')
         end if
         if (tag_i /= am%tag1) cycle

         if (am%tag2 /= 0) then
            do j=1,size(am%rscov%data)
               if(.not. get_value(am%rscov%data(j)%params, 'tag', tag_j)) then
                  call system_abort('cannot find tag in rscov%data('//j//')')
               end if
               if (.not. get_value(am%rscov%data(j)%params, 'time', time_j)) then
                  call system_abort('cannot find time in rscov%data('//i//')')
               end if
               if (tag_j /= am%tag2) cycle

               call overlap_rho_rot_rho(am%rscov,am%rscov%data(i),am%rscov%data(j),x,doverlap=doverlap)               
               weight = 1.0_dp !exp(-abs(time_i-time_j)/am%rscov%tau)
               sum_doverlap_func = sum_doverlap_func + weight*doverlap
            end do
         else
            call overlap_rho_rot_rho(am%rscov,am%rscov%data(i),am%at1,x,doverlap=doverlap)
            weight = 1.0_dp !exp(-abs(time_i-time_j)/am%rscov%tau)
            sum_doverlap_func = sum_doverlap_func + weight*doverlap
         end if
      end do
      sum_doverlap_func = -sum_doverlap_func

   end function sum_doverlap_func

   function rotmat(a)
      real(dp), dimension(3), intent(in) :: a
      real(dp), dimension(3,3) :: rotmat

      real(dp) :: c1, c2, c3, s1, s2, s3

      c1 = cos(a(1))
      c2 = cos(a(2))
      c3 = cos(a(3))
      s1 = sin(a(1))
      s2 = sin(a(2))
      s3 = sin(a(3))

      rotmat(1,1) = c1*c2*c3 - s1*s3;  rotmat(1,2) = -c1*s2; rotmat(1,3) = c2*c1*s3 + c3*s1
      rotmat(2,1) = c3*s2;             rotmat(2,2) = c2;     rotmat(2,3) = s3*s2
      rotmat(3,1) = -c1*s3 - c3*c2*s1; rotmat(3,2) = s2*s1;  rotmat(3,3) = c1*c3 - c2*s1*s3

   endfunction rotmat


   !% Find the quaterion corresponding to the rotation which maximise the overlap
   !% between the atomic densities of at1 and at2, starting minimisations at each
   !% quaterion in the (4,n) array qs
   function realspacecovariance_align_pair(this, at1, at2, qs, overlap) result(q)
     type(RealSpaceCovariance), intent(in), target :: this
     type(Atoms), intent(in), target :: at1, at2
     real(dp), intent(in) :: qs(:,:)
     real(dp), optional, intent(out) :: overlap
     real(dp) q(4)

     real(dp) :: this_overlap, max_overlap, this_q(4)
     type(minimise_overlap) :: am
     integer :: am_data_size, n_minim, i
     character(len=1), dimension(1) :: am_mold
     character(len=1), allocatable :: am_data(:)

     am_data_size = size(transfer(am,am_mold))
     allocate(am_data(am_data_size))
     am%rscov => this
     am%at1 => at1
     am%at2 => at2   
     am_data = transfer(am,am_mold)
     max_overlap = -huge(1.0_dp)

     do i=1,size(qs,2)
        this_q = qs(:,i)
        n_minim = minim(this_q, overlap_func,doverlap_func, method='cg', &
             convergence_tol=this%covariance_tol, max_steps = 100, data=am_data)

        this_overlap = -overlap_func(this_q,am_data=am_data)
        if (this_overlap > max_overlap) then
           max_overlap = this_overlap
           q = this_q
        end if
     end do

     q = q/norm(q)
     if (present(overlap)) overlap = max_overlap
     deallocate(am_data)
     
   end function realspacecovariance_align_pair
   
   !% Align the configuration 'at' with the set of configurations in 'this%data'
   !% which have tag equal to 'tag'. Objective function is sum of overlaps.
   !% A single minimisation is performed, starting from quaterions in 'qs'.
   function realspacecovariance_align_with_data_set(this, at, tag, qs, overlap) result(q)
     type(RealSpaceCovariance), intent(in), target :: this
     type(Atoms), intent(in), target :: at
     integer, intent(in) :: tag
     real(dp), intent(in) :: qs(:,:)
     real(dp), optional, intent(out) :: overlap
     real(dp) :: q(4)

     real(dp) :: this_overlap, max_overlap, this_q(4)
     type(minimise_overlap) :: am
     integer :: am_data_size, n_minim, i
     character(len=1), dimension(1) :: am_mold
     character(len=1), allocatable :: am_data(:)

     am_data_size = size(transfer(am,am_mold))
     allocate(am_data(am_data_size))
     am%rscov => this
     am%at1 => at
     am%tag1 = tag
     am%tag2 = 0
     am_data = transfer(am, am_mold)
     max_overlap = -huge(1.0_dp)

     call n_test_gradient(q,sum_overlap_func,sum_doverlap_func,data=am_data)
     print*, test_gradient(q,overlap_func,doverlap_func,data=am_data)

     do i=1,size(qs,2)
        this_q = qs(:,i)
        n_minim = minim(this_q, sum_overlap_func, sum_doverlap_func, method='cg', &
             convergence_tol=this%covariance_tol, max_steps = 100, data=am_data)
        this_overlap = -sum_overlap_func(this_q,am_data=am_data)
        if (this_overlap > max_overlap) then
           max_overlap = this_overlap
           q = this_q
        end if
     end do

     q = q/norm(q)
     if (present(overlap)) overlap = max_overlap
     deallocate(am_data)   
     
   end function realspacecovariance_align_with_data_set


   !% Align two sets of configurations, identified by the tags 'tag1' and 'tag2'
   function realspacecovariance_align_data_sets(this, tag1, tag2, qs, overlap) result(q)
     type(RealSpaceCovariance), intent(in), target :: this
     integer, intent(in) :: tag1, tag2
     real(dp), intent(in) :: qs(:,:)
     real(dp), optional, intent(out) :: overlap
     real(dp) q(4)

     real(dp) :: this_overlap, max_overlap, this_q(4)
     type(minimise_overlap) :: am
     integer :: am_data_size, n_minim, i
     character(len=1), dimension(1) :: am_mold
     character(len=1), allocatable :: am_data(:)

     am_data_size = size(transfer(am,am_mold))
     allocate(am_data(am_data_size))
     am%rscov => this
     am%tag1 = tag1
     am%tag2 = tag2
     am_data = transfer(am, am_mold)
     max_overlap = -huge(1.0_dp)

     

     do i=1,size(qs,2)
        this_q = qs(:,i)
        n_minim = minim( this_q, sum_overlap_func, sum_doverlap_func, method='cg', &
             convergence_tol=this%covariance_tol, max_steps = 100, data=am_data)
        this_overlap = -sum_overlap_func(this_q,am_data=am_data)
        if (this_overlap > max_overlap) then
           max_overlap = this_overlap
           q = this_q
        end if
     end do

     q = q/norm(q)
     if (present(overlap)) overlap = -sum_overlap_func(q,am_data=am_data)     
     deallocate(am_data)     

   end function realspacecovariance_align_data_sets

   function realspacecovariance_covariance(this, at1, at2, q, distance, force_covariance, force_difference, do_align) result(covariance)
     type(realspacecovariance), intent(inout) :: this
     type(Atoms), intent(in) :: at1, at2
     real(dp), intent(inout) :: q(4)
     real(dp), intent(out), optional :: distance, force_covariance, force_difference
     logical, intent(in), optional :: do_align

     real(dp) :: covariance
     real(dp) :: d2, rho1_2, rho2_2, rho_12, rot(3,3)
     real(dp), dimension(:,:), pointer :: force1_ptr, force2_ptr
     logical my_do_align
     
     my_do_align = optional_default(.true., do_align)
     if (my_do_align) then
        q = align(this, at1, at2, reshape(q, (/1,4/)))
     end if
           
     if(present(force_covariance) .or. present(force_difference)) then
        rot = rotmat_quaternion(q)
        call assign_property_pointer(at1, 'force', ptr=force1_ptr)
        call assign_property_pointer(at2, 'force', ptr=force2_ptr)
        if (present(force_difference)) &
             force_difference = sqrt(sum((force1_ptr(:,1) - matmul(rot, force2_ptr(:,1)))**2))
        if (present(force_covariance)) &
             force_covariance = dot_product( force1_ptr(:,1), matmul(rot,force2_ptr(:,1)) )
     endif

     call overlap_rho_rot_rho(this,at1,at1,overlap=rho1_2)
     call overlap_rho_rot_rho(this,at2,at2,overlap=rho2_2)
     call overlap_rho_rot_rho(this,at1,at2,q=q, overlap=rho_12)
     d2 = rho1_2 + rho2_2 - 2.0_dp*rho_12
     if(present(distance)) distance = d2
     covariance = exp(-0.5_dp*d2/this%sigma_covariance**2)
     
   end function realspacecovariance_covariance

   subroutine realspacecovariance_teach(this, do_optimise, use_qr)
     type(realspacecovariance), intent(inout), target :: this
     logical, intent(in) :: do_optimise, use_qr

     integer :: am_data_size
     character(len=1), dimension(1) :: am_mold
     character(len=1), allocatable :: am_data(:)
     real(dp) :: x_sigma(2)

     type(optimise_likelihood) :: am_likelihood
     integer i,j,n

     do i=1,size(this%data)
        call centre_on_atom_1(this%data(i))
     end do

     call system_timer('covariance')
     call verbosity_push(PRINT_SILENT)
     !$omp parallel default(none) shared(this) private(j)
     !$omp do
     do i=1,size(this%data)
        do j=i,size(this%data)
           this%covariance_matrix(i,j) = covariance(this, this%data(i), this%data(j), &
                q=this%alignment(:,i,j), &
                distance=this%distance_matrix(i,j), force_covariance=this%force_covariance_matrix(i,j), &
                force_difference=this%force_difference_matrix(i,j))
         
           this%covariance_matrix(j,i) = this%covariance_matrix(i,j)
           this%distance_matrix(j,i) = this%distance_matrix(i,j)
           this%force_covariance_matrix(j,i) = this%force_covariance_matrix(i,j)
           this%force_difference_matrix(j,i) = this%force_difference_matrix(i,j)
           this%alignment(:,j,i) = (/ this%alignment(1,i,j), -this%alignment(2:4,j,i) /) !quaterion conjugate
        end do
        this%covariance_matrix(i,i) = this%covariance_matrix(i,i) + this%sigma_error**2
     end do
     !$omp end parallel
     call system_timer('covariance')
     call verbosity_pop()

     mainlog%prefix = 'DISTANCE '
     call print(this%distance_matrix)

     mainlog%prefix = 'FORCEDIF '
     call print(this%force_difference_matrix)

     mainlog%prefix = 'FORCECOV '
     call print(this%force_covariance_matrix)

     mainlog%prefix = 'COV1 '
     call print(this%covariance_matrix)
     
     call inverse(this%covariance_matrix, this%cov_inv)
     mainlog%prefix = 'INV1 '
     call print(this%cov_inv)
     mainlog%prefix = ''

     if (do_optimise) then
        am_likelihood%use_qr = use_qr
        am_likelihood%rscov => this
        am_data_size = size(transfer(am_likelihood,am_mold))
        allocate(am_data(am_data_size))
        am_data = transfer(am_likelihood,am_mold)
     
        ! Testing the implementation of likelihood and dlikelihood.
        if(.false.) then
           call print('starting test gradient')
           call verbosity_push(PRINT_VERBOSE)
           write (*,*) test_gradient((/this%sigma_error,this%sigma_covariance/),likelihood,dlikelihood,data=am_data)
           call verbosity_pop()
        endif
        call print_title('likelihood optimisation')

        if (use_qr) then
           call print('doing likelihood optimisation with QR decomposition')
        else
           call print('doing likelihood optimisation with Cholesky decomposition')
        end if
        x_sigma = (/this%sigma_error,this%sigma_covariance/)
        n = minim( x_sigma, likelihood,dlikelihood, method='cg', convergence_tol=this%covariance_tol, max_steps = 100, data=am_data)

        deallocate(am_data)

        this%sigma_error = x_sigma(1)
        this%sigma_covariance = x_sigma(2)
        call print('Optimisied hypers. ERROR = '//this%sigma_error//' COVARIANCE = '//this%sigma_covariance)

        ! update covariance matrix
        this%covariance_matrix = exp(-0.5_dp * this%distance_matrix/this%sigma_covariance**2)
        do i = 1, size(this%covariance_matrix,1)
           this%covariance_matrix(i,i) = this%covariance_matrix(i,i) + this%sigma_error**2
        enddo
     end if

     mainlog%prefix = 'COV2 '
     call print(this%covariance_matrix)
     
     call inverse(this%covariance_matrix, this%cov_inv)
     mainlog%prefix = 'INV2 '
     call print(this%cov_inv)
     mainlog%prefix = ''
     
   end subroutine realspacecovariance_teach

   subroutine realspacecovariance_predict(this, test, force, error)
     type(realspacecovariance), intent(inout), target :: this
     type(Atoms), intent(inout), target :: test
     real(dp), intent(out) :: force(3)
     real(dp), intent(out) :: error

     type(Atoms), save :: rot_at
     real(dp) :: rot(3,3), d2, rho1_2, rho2_2, rho_12, kappa, unit_q(4)
     real(dp), pointer, dimension(:,:) :: rot_force_ptr
     integer i
     
     !$omp threadprivate(rot_at)

     call centre_on_atom_1(test)
     call overlap_rho_rot_rho(this,test,test,overlap=rho2_2)

     force(:) = 0.0_dp
     this%rot_force(:,:) = 0.0_dp

     call verbosity_push(PRINT_SILENT)
     !$omp parallel default(none) shared(test, this, rho2_2) private(rot, rot_force_ptr, d2, rho1_2, rho_12)
     !$omp do
     do i=1,size(this%data)
        this%q_min(:,i) = align(this, test, this%data(i), reshape(this%q_min(:,i), (/1,4/)))
        rot = rotmat_quaternion(this%q_min(:,i))

        call atoms_copy_without_connect(rot_at, this%data(i))
        rot_at%pos = matmul(rot,rot_at%pos)
        call assign_property_pointer(rot_at, 'force', ptr=rot_force_ptr)
        rot_force_ptr(:,:) = matmul(rot, rot_force_ptr)
        this%rot_force(:,i) = rot_force_ptr(:,1)

        call overlap_rho_rot_rho(this,rot_at,rot_at,overlap=rho1_2)
        call overlap_rho_rot_rho(this,rot_at,test,overlap=rho_12)
        d2 = rho1_2 + rho2_2 - 2.0_dp*rho_12

        this%k(i) = exp(-0.5_dp*d2/this%sigma_covariance**2)
     end do
     !$omp end parallel
     call verbosity_pop()

     force = matmul(this%rot_force, matmul(this%cov_inv, this%k))

     unit_q = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
     kappa = covariance(this, test, test, q=unit_q)
     error = sqrt(kappa - min(kappa, (this%k .dot. matmul(this%cov_inv, this%k))))
        
   end subroutine realspacecovariance_predict

   subroutine realspacecovariance_predict_2(this, test, force, error)
     type(realspacecovariance), intent(inout), target :: this
     type(Atoms), intent(inout), target :: test
     real(dp), intent(out) :: force(3)
     real(dp), intent(out) :: error
     
     type(Atoms), save :: rot_at
     real(dp) :: rot(3,3), d2, rho1_2, rho2_2, rho_12, kappa, unit_q(4)
     real(dp), pointer, dimension(:,:) :: rot_force_ptr
     integer i,j
     
     !$omp threadprivate(rot_at)

     call centre_on_atom_1(test)
     call overlap_rho_rot_rho(this,test,test,overlap=rho2_2)
     
     force(:) = 0.0_dp
     this%rot_force(:,:) = 0.0_dp

     call verbosity_push(PRINT_SILENT)

     !$omp parallel default(none) shared(test, this)
     !$omp do
     do i=1,size(this%data)
        this%q_min(:,i) = align(this, test, this%data(i), reshape(this%q_min(:,i), (/1,4/)))
     end do
     !$omp end parallel

     call verbosity_pop()

     !$omp parallel default(none) shared(this, test, rho2_2) private(j, rot, rot_force_ptr, d2, rho1_2, rho_12)
     !$omp do
     do i=1,size(this%data)
        ! rotate data(i) into frame of test config
        rot = rotmat_quaternion(this%q_min(:,i))
        call atoms_copy_without_connect(rot_at, this%data(i))
        rot_at%pos = matmul(rot,rot_at%pos)

        do j=i+1,size(this%data)
           
           ! compute covariance when data(j) is also rotated into frame of test config
           this%covariance_matrix(i,j) = covariance(this, rot_at, this%data(j), &
                distance=this%distance_matrix(i,j), q=this%q_min(:,j))

           write (*,*) i, j, this%distance_matrix(i,j)

           this%distance_matrix(j,i) = this%distance_matrix(i,j)
           this%covariance_matrix(j,i) = this%covariance_matrix(i,j)

        end do
        this%distance_matrix(i,i) = 0.0_dp
        this%covariance_matrix(i,i) = 1.0_dp + this%sigma_error**2

        call overlap_rho_rot_rho(this,rot_at,rot_at,overlap=rho1_2)
        call overlap_rho_rot_rho(this,rot_at,test,overlap=rho_12)
        d2 = rho1_2 + rho2_2 - 2.0_dp*rho_12
        this%k(i) =  exp(-0.5_dp*d2/this%sigma_covariance**2)

        ! store rotated force
        call assign_property_pointer(rot_at, 'force', ptr=rot_force_ptr)
        rot_force_ptr(:,:) = matmul(rot, rot_force_ptr)
        this%rot_force(:,i) = rot_force_ptr(:,1)

     end do
     !$omp end parallel

     call inverse(this%covariance_matrix, this%cov_inv)
     call print('cov=')
     call print(this%covariance_matrix)
     call print('inv=')
     call print(this%cov_inv)

     force = matmul(this%rot_force, matmul(this%cov_inv, this%k))
     unit_q = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
     kappa = covariance(this, test, test, q=unit_q)
     error = sqrt(kappa - min(kappa, (this%k .dot. matmul(this%cov_inv, this%k))))
        
   end subroutine realspacecovariance_predict_2

   subroutine realspacecovariance_predict_3(this, test, crystal, k, force, error)
     type(realspacecovariance), intent(inout), target :: this
     type(Atoms), intent(inout), target :: test, crystal
     integer, intent(in) :: k
     real(dp), intent(out) :: force(3)
     real(dp), intent(out) :: error
     
     type(Atoms), save :: rot_at
     type(Quaternion) :: QQ_Ci, QQ_Cj, QQ_Ti, QQ_Tj, QQ_Tk, QQ_ij
     real(dp) :: rot(3,3), d2, rho1_2, rho2_2, rho_12, kappa, unit_q(4), q_Tk(4), q_Ti(4), q_Tj(4), q_ij(4)
     real(dp), pointer, dimension(:,:) :: rot_force_ptr
     integer i,j
     
     !$omp threadprivate(rot_at)

     call centre_on_atom_1(test)
     call overlap_rho_rot_rho(this,test,test,overlap=rho2_2)
     
     force(:) = 0.0_dp
     this%rot_force(:,:) = 0.0_dp

     call verbosity_push(PRINT_SILENT)

     ! (1) align all learning set to cold crystal

     !$omp parallel default(none) shared(test, crystal, this)
     !$omp do
     do i=1,size(this%data)
        this%q_min(:,i) = align(this, crystal, this%data(i), reshape(this%q_min(:,i), (/1,4/)))
     end do
     !$omp end parallel

     call verbosity_pop()

     ! (2) align test config to one representative element k of the learning set
     q_Tk = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
     q_Tk = align(this, test, this%data(k), reshape(q_Tk, (/1,4/)))
     call print('q_Tk = '//q_Tk)
     QQ_Tk = q_Tk

     !$omp parallel default(none) shared(this, test, rho2_2, QQ_Tk) private(j, rot, rot_force_ptr, d2, rho1_2, rho_12, QQ_Ci, QQ_Ti, q_Ti, QQ_Cj, QQ_Tj, q_Tj, QQ_ij, q_ij)
     !$omp do
     do i=1,size(this%data)
        
        ! q_T,i = conj(q_C,k) * q_C,i * q_T,k

        QQ_Ci = this%q_min(:,i)
        QQ_Ti = (.conj. QQ_Tk) * QQ_Ci * QQ_Tk
        q_Ti = QQ_Ti

        ! rotate data(i) into frame of test config
        rot = rotmat_quaternion(q_Ti)
        call atoms_copy_without_connect(rot_at, this%data(i))
        rot_at%pos = matmul(rot,rot_at%pos)

        do j=i+1,size(this%data)

           ! q_T,j = conj(q_C,k) * q_C,j * q_T,k
           ! q_i,j = conj(q_T,i) * q_T,j

           QQ_Cj = this%q_min(:,j)
           QQ_Tj = (.conj. QQ_Tk) * QQ_Cj * QQ_Tk
           q_Tj = QQ_Tj

           QQ_ij = (.conj. QQ_Ti) * QQ_Tj
           q_ij = QQ_ij

           write (*,*) i, j, q_ij

           this%covariance_matrix(i,j) = covariance(this, rot_at, this%data(j), &
                distance=this%distance_matrix(i,j), q=q_ij, do_align=.false.)

           this%distance_matrix(j,i) = this%distance_matrix(i,j)
           this%covariance_matrix(j,i) = this%covariance_matrix(i,j)

        end do
        this%distance_matrix(i,i) = 0.0_dp
        this%covariance_matrix(i,i) = 1.0_dp + this%sigma_error**2

        call overlap_rho_rot_rho(this,rot_at,rot_at,overlap=rho1_2)
        call overlap_rho_rot_rho(this,rot_at,test,overlap=rho_12)
        d2 = rho1_2 + rho2_2 - 2.0_dp*rho_12
        this%k(i) =  exp(-0.5_dp*d2/this%sigma_covariance**2)

        ! store rotated force
        call assign_property_pointer(rot_at, 'force', ptr=rot_force_ptr)
        rot_force_ptr(:,:) = matmul(rot, rot_force_ptr)
        this%rot_force(:,i) = rot_force_ptr(:,1)

     end do
     !$omp end parallel

     call inverse(this%covariance_matrix, this%cov_inv)
     call print('cov=')
     call print(this%covariance_matrix)
     call print('inv=')
     call print(this%cov_inv)

     force = matmul(this%rot_force, matmul(this%cov_inv, this%k))
     unit_q = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
     kappa = covariance(this, test, test, q=unit_q, do_align=.false.)
     error = sqrt(kappa - min(kappa, (this%k .dot. matmul(this%cov_inv, this%k))))
        
   end subroutine realspacecovariance_predict_3


   function rotmat_quaternion(q)
      real(dp), dimension(4), intent(in) :: q
      real(dp), dimension(3,3) :: rotmat_quaternion

      real(dp), dimension(4) :: q_in          

      q_in = q / sqrt(dot_product(q,q))

      rotmat_quaternion(1,1) = 1.0_dp - 2.0_dp*q_in(3)**2 - 2.0_dp*q_in(4)**2
      rotmat_quaternion(1,2) = 2.0_dp*q_in(2)*q_in(3) - 2.0_dp*q_in(1)*q_in(4)
      rotmat_quaternion(1,3) = 2.0_dp*q_in(1)*q_in(3) + 2.0_dp*q_in(2)*q_in(4)
      rotmat_quaternion(2,1) = 2.0_dp*q_in(2)*q_in(3) + 2.0_dp*q_in(1)*q_in(4)
      rotmat_quaternion(2,2) = 1.0_dp - 2.0_dp*q_in(2)**2 - 2.0_dp*q_in(4)**2
      rotmat_quaternion(2,3) = -2.0_dp*q_in(1)*q_in(2) + 2.0_dp*q_in(3)*q_in(4)
      rotmat_quaternion(3,1) = -2.0_dp*q_in(1)*q_in(3) + 2.0_dp*q_in(2)*q_in(4)
      rotmat_quaternion(3,2) = 2.0_dp*q_in(1)*q_in(2) + 2.0_dp*q_in(3)*q_in(4)
      rotmat_quaternion(3,3) = 1.0_dp - 2.0_dp*q_in(2)**2 - 2.0_dp*q_in(3)**2

   endfunction rotmat_quaternion

   function drotmat_dquaternion(q)
      real(dp), dimension(4), intent(in) :: q
      real(dp), dimension(3,3,4) :: drotmat_dquaternion

      drotmat_dquaternion(1,1,1) = 4*q(1)*(q(3)**2 + q(4)**2)
      drotmat_dquaternion(1,2,1) = -4*q(1)*q(2)*q(3) + 2*q(1)**2*q(4) - 2*q(4)*(q(2)**2 + q(3)**2 + q(4)**2)
      drotmat_dquaternion(1,3,1) = 2*(-(q(1)**2*q(3)) - 2*q(1)*q(2)*q(4) + q(3)*(q(2)**2 + q(3)**2 + q(4)**2))
      drotmat_dquaternion(2,1,1) = 2*(-2*q(1)*q(2)*q(3) - q(1)**2*q(4) + q(4)*(q(2)**2 + q(3)**2 + q(4)**2))
      drotmat_dquaternion(2,2,1) = 4*q(1)*(q(2)**2 + q(4)**2)
      drotmat_dquaternion(2,3,1) = 2*q(1)**2*q(2) - 4*q(1)*q(3)*q(4) - 2*q(2)*(q(2)**2 + q(3)**2 + q(4)**2)
      drotmat_dquaternion(3,1,1) = 2*q(1)**2*q(3) - 4*q(1)*q(2)*q(4) - 2*q(3)*(q(2)**2 + q(3)**2 + q(4)**2)
      drotmat_dquaternion(3,2,1) = 2*(-(q(1)**2*q(2)) - 2*q(1)*q(3)*q(4) + q(2)*(q(2)**2 + q(3)**2 + q(4)**2))
      drotmat_dquaternion(3,3,1) = 4*q(1)*(q(2)**2 + q(3)**2)

      drotmat_dquaternion(1,1,2) = 4*q(2)*(q(3)**2 + q(4)**2)
      drotmat_dquaternion(1,2,2) = 2*(q(1)**2*q(3) + 2*q(1)*q(2)*q(4) + q(3)*(-q(2)**2 + q(3)**2 + q(4)**2))
      drotmat_dquaternion(1,3,2) = 2*(-2*q(1)*q(2)*q(3) + (q(1)**2 - q(2)**2 + q(3)**2)*q(4) + q(4)**3)
      drotmat_dquaternion(2,1,2) = 2*(q(1)**2*q(3) - 2*q(1)*q(2)*q(4) + q(3)*(-q(2)**2 + q(3)**2 + q(4)**2))
      drotmat_dquaternion(2,2,2) = -4*q(2)*(q(1)**2 + q(3)**2)
      drotmat_dquaternion(2,3,2) = -2*(q(1)**3 + 2*q(2)*q(3)*q(4) + q(1)*(-q(2)**2 + q(3)**2 + q(4)**2))
      drotmat_dquaternion(3,1,2) = 2*(2*q(1)*q(2)*q(3) + (q(1)**2 - q(2)**2 + q(3)**2)*q(4) + q(4)**3)
      drotmat_dquaternion(3,2,2) = 2*(q(1)**3 - 2*q(2)*q(3)*q(4) + q(1)*(-q(2)**2 + q(3)**2 + q(4)**2))
      drotmat_dquaternion(3,3,2) = -4*q(2)*(q(1)**2 + q(4)**2)

      drotmat_dquaternion(1,1,3) = -4*(q(1)**2 + q(2)**2)*q(3)
      drotmat_dquaternion(1,2,3) = 2*(q(1)**2*q(2) + 2*q(1)*q(3)*q(4) + q(2)*(q(2)**2 - q(3)**2 + q(4)**2))
      drotmat_dquaternion(1,3,3) = 2*(q(1)**3 - 2*q(2)*q(3)*q(4) + q(1)*(q(2)**2 - q(3)**2 + q(4)**2))
      drotmat_dquaternion(2,1,3) = 2*(q(1)**2*q(2) - 2*q(1)*q(3)*q(4) + q(2)*(q(2)**2 - q(3)**2 + q(4)**2))
      drotmat_dquaternion(2,2,3) = 4*q(3)*(q(2)**2 + q(4)**2)
      drotmat_dquaternion(2,3,3) = 2*(2*q(1)*q(2)*q(3) + (q(1)**2 + q(2)**2 - q(3)**2)*q(4) + q(4)**3)
      drotmat_dquaternion(3,1,3) = -2*(q(1)**3 + 2*q(2)*q(3)*q(4) + q(1)*(q(2)**2 - q(3)**2 + q(4)**2))
      drotmat_dquaternion(3,2,3) = 2*(-2*q(1)*q(2)*q(3) + (q(1)**2 + q(2)**2 - q(3)**2)*q(4) + q(4)**3)
      drotmat_dquaternion(3,3,3) = -4*q(3)*(q(1)**2 + q(4)**2)

      drotmat_dquaternion(1,1,4) = -4*(q(1)**2 + q(2)**2)*q(4)
      drotmat_dquaternion(1,2,4) = -2*(q(1)**3 + 2*q(2)*q(3)*q(4) + q(1)*(q(2)**2 + q(3)**2 - q(4)**2))
      drotmat_dquaternion(1,3,4) = 2*q(2)*(q(1)**2 + q(2)**2 + q(3)**2) - 4*q(1)*q(3)*q(4) - 2*q(2)*q(4)**2
      drotmat_dquaternion(2,1,4) = 2*(q(1)**3 - 2*q(2)*q(3)*q(4) + q(1)*(q(2)**2 + q(3)**2 - q(4)**2))
      drotmat_dquaternion(2,2,4) = -4*(q(1)**2 + q(3)**2)*q(4)
      drotmat_dquaternion(2,3,4) = 2*q(3)*(q(1)**2 + q(2)**2 + q(3)**2) + 4*q(1)*q(2)*q(4) - 2*q(3)*q(4)**2
      drotmat_dquaternion(3,1,4) = 2*q(2)*(q(1)**2 + q(2)**2 + q(3)**2) + 4*q(1)*q(3)*q(4) - 2*q(2)*q(4)**2
      drotmat_dquaternion(3,2,4) = 2*q(3)*(q(1)**2 + q(2)**2 + q(3)**2) - 4*q(1)*q(2)*q(4) - 2*q(3)*q(4)**2
      drotmat_dquaternion(3,3,4) = 4*(q(2)**2 + q(3)**2)*q(4)

      drotmat_dquaternion = drotmat_dquaternion / dot_product(q,q)**2

   endfunction drotmat_dquaternion

   function likelihood(x,am_data)
      real(dp), dimension(:) :: x
      character(len=1), dimension(:), optional :: am_data
      real(dp) :: likelihood

      type(optimise_likelihood) :: am
      real(dp) :: sigma_error, sigma_covariance
      real(dp), dimension(:,:), allocatable :: ICF
      integer :: i
      type(LA_Matrix) :: LA_covariance

      am = transfer(am_data,am)
      sigma_error = x(1)
      sigma_covariance = x(2) 

      am%rscov%covariance_matrix = exp(-0.5_dp * am%rscov%distance_matrix/sigma_covariance**2)
      do i = 1, size(am%rscov%covariance_matrix,1)
         am%rscov%covariance_matrix(i,i) = am%rscov%covariance_matrix(i,i) + sigma_error**2
      enddo

      LA_covariance = am%rscov%covariance_matrix
      allocate(ICF(size(am%rscov%covariance_matrix,1),size(am%rscov%covariance_matrix,1)))

      if (am%use_qr) then
         call Matrix_QR_Solve(LA_covariance,am%rscov%force_covariance_matrix,ICF)
      else
         call Matrix_Solve(LA_covariance,am%rscov%force_covariance_matrix,ICF)
      end if

      likelihood = 3 * 0.5_dp * LA_Matrix_LogDet(LA_covariance) + 0.5_dp * trace(ICF)

      call finalise(LA_covariance)
      if(allocated(ICF)) deallocate(ICF)

   endfunction likelihood

   function dlikelihood(x,am_data)
      real(dp), dimension(:) :: x
      character(len=1), dimension(:), optional :: am_data
      real(dp), dimension(size(x)) :: dlikelihood

      type(optimise_likelihood) :: am
      real(dp) :: sigma_error, sigma_covariance
      real(dp), dimension(:,:), allocatable :: ICF, IC, ICDCDS, DCDS
      integer :: i, n
      type(LA_Matrix) :: LA_covariance

      am = transfer(am_data,am)
      sigma_error = x(1)
      sigma_covariance = x(2) 

      n = size(am%rscov%covariance_matrix,1)
      allocate(ICF(n,n),DCDS(n,n),IC(n,n),ICDCDS(n,n))

      am%rscov%covariance_matrix = exp(-0.5_dp * am%rscov%distance_matrix/sigma_covariance**2)
      DCDS = am%rscov%covariance_matrix * am%rscov%distance_matrix / sigma_covariance**3

      do i = 1, n
         am%rscov%covariance_matrix(i,i) = am%rscov%covariance_matrix(i,i) + sigma_error**2
      enddo
      
      LA_covariance = am%rscov%covariance_matrix

      if (am%use_qr) then
         call Matrix_QR_Solve(LA_covariance,am%rscov%force_covariance_matrix,ICF)
         call Matrix_QR_Solve(LA_covariance,DCDS,ICDCDS)
         call LA_Matrix_QR_Inverse(LA_covariance,IC)
      else
         call Matrix_Solve(LA_covariance,am%rscov%force_covariance_matrix,ICF)
         call Matrix_Solve(LA_covariance,DCDS,ICDCDS)
         call LA_Matrix_Inverse(LA_covariance,IC)
      end if

      dlikelihood(1) = sigma_error * (3*trace(IC) - sum(IC*ICF) )
      dlikelihood(2) = + 0.5_dp * 3*trace(ICDCDS) - 0.5_dp * sum(ICDCDS*transpose(ICF))

      call finalise(LA_covariance)
      if(allocated(ICF)) deallocate(ICF)
      if(allocated(DCDS)) deallocate(DCDS)
      if(allocated(IC)) deallocate(IC)
      if(allocated(ICDCDS)) deallocate(ICDCDS)

   endfunction dlikelihood

   subroutine centre_on_atom_1(at)
     type(Atoms), intent(inout) :: at

     integer i

     ! make all positions relative to first atom in file
     do i = 2, at%N
        at%pos(:,i) = diff_min_image(at,1,i)
     enddo
     at%pos(:,1) = 0.0_dp

   end subroutine centre_on_atom_1

endmodule real_space_covariance_module

