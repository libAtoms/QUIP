program force_gaussian_prediction

   use libatoms_module
   implicit none

   type(Atoms)                                  :: at_in
   type(CInOutput)                              :: in
   type(Dictionary)                             :: params
   real(dp)                                     :: r_cut, r_min, m_min, m_max, z_len, theta, thresh
   real(dp), parameter                          :: TOL_REAL = 1e-7_dp, SCALE_IVS = 100.0_dp
   integer                                      :: i,j, k, n, t, at_n, n_loop, preci, r_mesh, m_mesh
   real(dp)                                     :: force(3), z_inv(3,3)
   real(dp), dimension(:), allocatable          :: r_grid, m_grid, sigma, covariance_predict, force_proj_ivs_predict
   real(dp), dimension(:,:,:), allocatable      :: z_matrix_normalised, z_matrix
   real(dp), dimension(:,:), allocatable        :: force_proj_ivs, covariance, inv_covariance, distance_ivs
   real(dp), dimension(:,:), allocatable        :: z_matrix_normalised_predict, z_matrix_predict
   real(dp), dimension(:,:), pointer            :: force_ptr
   logical                                      :: spherical_cluster_teach, spherical_cluster_predict 
   character(STRING_LENGTH)                     :: teaching_file, grid_file, test_file

   call system_initialise(enable_timing=.true.)

   call initialise(params)
   call param_register(params, 'theta',  '2.5', theta, "the correlation length of the data")
   call param_register(params, 'r_min',  '0.5', r_min, "the minimum radius for the spherical atomic environment")
   call param_register(params, 'r_cut',  '8.0', r_cut, "the cutoff radius for the spherical atomic environment")
   call param_register(params, 'm_min',  '1.0', m_min, "the minimum m for calculating the atomic environment")
   call param_register(params, 'm_max',  '5.0', m_max, "the maxmium m for calculating the atomic environment")
   call param_register(params, 'thresh', '1.0', thresh, "the threshold for doing the Sigular Value Decompostion of the Covraince Matrix")
   call param_register(params, 'preci',  '6',   preci,  "the screening accuracy on the edge atoms")
   call param_register(params, 'r_mesh', '6',   r_mesh, "grid finess of r0")
   call param_register(params, 'm_mesh', '6',   m_mesh, "grid finess of m")
   call param_register(params, 'spherical_cluster_teach', 'T', spherical_cluster_teach, "only the central atom in the cluster are considered when doing teaching")
   call param_register(params, 'spherical_cluster_predict', 'T', spherical_cluster_predict, "only the central atom in the cluster are considered when doing predicting")
   call param_register(params, 'teaching_file', 'data.xyz', teaching_file, "file to read teaching configurations from")
   call param_register(params, 'grid_file', 'grid.xyz', grid_file, "file to generate the proper pairs of r0 and m")
   call param_register(params, 'test_file', 'test.xyz', test_file, "file to read the testing configurations from")

   if (.not. param_read_args(params, task="fgp command line params")) then
      call print("Usage: fgp [options]")
      call system_abort("Confused by command line arguments")
   end if
   call print_title('params')
   call param_print(params)
   call finalise(params)

   call print_title('Testing the grid')
   ! generate the r0_m grids
   call initialise(in, grid_file, INPUT)
   call read(in, at_in)
   ! including the tranlational symmetry images
   call set_cutoff(at_in, r_cut)
   call calc_connect(at_in)
 
   call  grid_m_r0(at_in, r_mesh, m_mesh, r_min, r_cut, m_min, m_max, preci, k, r_grid, m_grid)
   call  finalise(in)
   call print('The number of valid (r0, m):  '//k)

  ! to know the total number of teaching confs: n
  call initialise(in, teaching_file, INPUT)
  n = 0  
  do i=1,  in%n_frame
        call read(in, at_in)
        if (spherical_cluster_teach) then
            n = n + 1
        else
            n = n +  at_in%N
        endif   
  enddo
  call finalise(in)
  call print("The total number of teaching configurations:    "//n)
 
  ! the main work starts
  call initialise(in, teaching_file, INPUT)
  call print('Got '//in%n_frame//' input configurations.')

  allocate(z_matrix(k,3,n))   
  allocate(z_matrix_normalised(k,3, n))
  allocate(force_proj_ivs(k, n))

 t = 0 ! counter for loop the atomic configurations
 do i=1, in%n_frame
        call read(in, at_in)
        call assign_property_pointer(at_in, 'force', ptr=force_ptr)
        call print('the frame: '//i)

        ! Including the symmetry images
        call set_cutoff(at_in, r_cut)
        call calc_connect(at_in)

     if (spherical_cluster_teach) then        
        n_loop = 1
     else
        n_loop = at_in%N   
     endif

     do at_n =1,  n_loop
          t = t + 1
          force = force_ptr(:,at_n)
          call print("The atom:"//at_n)
           
        do j=1, k
           z_matrix(j, :, t) = internal_vector(at_in, r_grid(j), m_grid(j), at_n)*SCALE_IVS
           z_len = norm(z_matrix(j,:, t))
           if (z_len < TOL_REAL)    z_len=1.0_dp
           z_matrix_normalised(j,:,t) = z_matrix(j,:, t)/z_len
           write(*,*) "internal vectors (", r_grid(j), m_grid(j), ")",  z_matrix(j,1,t), z_matrix(j,2,t), z_matrix(j,3,t) 
        enddo

        ! store the force info. for every configuration
        do j=1, k
            force_proj_ivs(j,t) = dot_product(z_matrix_normalised(j,:,t), force) 
            call print("Force projection on IV"//j//" is:  "//force_proj_ivs(j,t))
        enddo   
 
      enddo  ! loop over the atoms
 end do     ! loop over the frames

 call finalise(in)

 call print_title('starting the teaching process')

 allocate(distance_ivs(n,n))
 allocate(sigma(k))
 do t = 1, k 
   do i = 1, n
        do j=1, n
              distance_ivs(i,j) = distance_bent_space(z_matrix(t,:,i), z_matrix(t,:,j), z_matrix_normalised(:,:,i), z_matrix_normalised(:,:,j)) 
        enddo
   enddo
   sigma(t) = maxval(distance_ivs(:,:))/theta 
   if  (abs(sigma(t))<TOL_REAL)  sigma(t)=TOL_REAL
 enddo
 call print('sigma is:  '//sigma)

! establish the Covariance Matrix
  allocate(covariance(n,n))
  do i = 1, n
      do j=1, n
             covariance(i,j) = cov(z_matrix(:,:,i), z_matrix(:,:,j), z_matrix_normalised(:,:,i), z_matrix_normalised(:,:,j), sigma, k) 
      enddo
  enddo
  
  do i=1, n
     write(*,*) "Covariance:", covariance(i,2), force_proj_ivs(3,i)-force_proj_ivs(3,2)
  enddo

 allocate(inv_covariance(n,n))

! before inverse operations on the covariance matrix, do SVD
! Sigular Value Decomposition (SVD): A = U*SIGMA*VT
 inv_covariance = inverse_svd_threshold(covariance, n, thresh)
 write(*,*) "MAX and MIN components of inverse_covariance:", maxval(inv_covariance(2,:)), minval(inv_covariance(2,:))

 deallocate(covariance)
 deallocate(distance_ivs)

 call print_title('starting the predicting process')

 allocate(covariance_predict(n))
 allocate(z_matrix_predict(k,3))
 allocate(z_matrix_normalised_predict(k,3))
 allocate(force_proj_ivs_predict(k))

 call initialise(in, test_file, INPUT)
 do i=1, in%n_frame 
   call read(in, at_in)
   call assign_property_pointer(at_in, 'force', ptr=force_ptr)
   call set_cutoff(at_in, r_cut)
   call calc_connect(at_in)
   if (spherical_cluster_predict) then
        n_loop = 1
   else
        n_loop = at_in%N
   endif

    do at_n=1, n_loop
       do j= 1, k
            z_matrix_predict(j,:) = internal_vector(at_in, r_grid(j), m_grid(j), at_n)*SCALE_IVS
            write(*,*) "internal vectors (",r_grid(j), m_grid(j),"):   ",z_matrix_predict(j,1), z_matrix_predict(j,2), z_matrix_predict(j,3)
            z_len = norm(z_matrix_predict(j,:))
            if (z_len < TOL_REAL)    z_len=1.0_dp
            z_matrix_normalised_predict(j,:) = z_matrix_predict(j,:)/z_len
       enddo

       do t= 1,n
             covariance_predict(t) = cov(z_matrix_predict, z_matrix(:,:,t), z_matrix_normalised_predict(:,:), z_matrix_normalised(:,:,t), sigma, k)
       enddo
   
       force_proj_ivs_predict(:) = matmul(covariance_predict, matmul(inv_covariance, transpose(force_proj_ivs(:,:)) )) 
       write(*,*) "The predicted force in IVs space:", force_proj_ivs_predict
       ! using least-squares to restore the force in the External Cartesian Space  
       z_inv = inverse_svd_threshold(matmul(transpose(z_matrix_normalised_predict), z_matrix_normalised_predict), 3, thresh)  
       force = z_inv .mult. transpose(z_matrix_normalised_predict) .mult. force_proj_ivs_predict
       call print("the force in external space:"//force//"  the original force:"//force_ptr(:, at_n))  
       call print("max error :    "//maxval(abs(force_ptr(:,at_n)-force)))
       call print("predicted error :"//(1.0_dp - covariance_predict .dot. matmul(inv_covariance, covariance_predict)))
    enddo  ! loop over postions
  enddo    ! loop over frames

 call finalise(in)

 deallocate(m_grid)
 deallocate(r_grid)
 deallocate(force_proj_ivs_predict)
 deallocate(z_matrix_normalised)
 deallocate(z_matrix_normalised_predict)
 deallocate(z_matrix)
 deallocate(z_matrix_predict)
 deallocate(sigma)
 deallocate(force_proj_ivs)
 deallocate(inv_covariance) 
 deallocate(covariance_predict)
 call system_finalise()

 
 contains 

 function internal_vector(at, r0, m, at_n)

      real(dp)                            :: internal_vector(3), delta_r(3), delta_r_len
      real(dp), intent(in)                :: r0, m
      type(Atoms), intent(in)             :: at
      integer                             :: i, j, at_n

      internal_vector = 0.0_dp
      do i=1, atoms_n_neighbours(at, at_n)
          j = atoms_neighbour(at, at_n, i, distance=delta_r_len, diff=delta_r) 
          internal_vector = internal_vector + (delta_r *exp(-((delta_r_len/r0)**m)) /delta_r_len)
      enddo
 endfunction  internal_vector


 function internal_vector_dm(at, r0, m, at_n)

      real(dp)                            :: internal_vector_dm, internal_vector_prim(3), delta_r(3), delta_r_len
      real(dp), intent(in)                :: r0, m
      type(Atoms), intent(in)             :: at
      integer                             :: i, j, at_n

      internal_vector_prim = 0.0_dp
      do i=1, atoms_n_neighbours(at, at_n)
          j = atoms_neighbour(at, at_n, i, distance=delta_r_len, diff=delta_r)
          internal_vector_prim = internal_vector_prim + (delta_r *exp(-((delta_r_len/r0)**m) *(-m) *log(delta_r_len/r0) ) /delta_r_len)
      enddo
      internal_vector_dm = norm(internal_vector_prim)
 endfunction  internal_vector_dm


 function cutoff_m(r_cut, r0, preci)
 
     real(dp)          :: r0, r_cut, cutoff_m
     integer           :: preci

     cutoff_m = log(real(preci)*log(10.0_dp)) / log(r_cut/r0)
 endfunction cutoff_m


 subroutine grid_m_r0(at, r_mesh, m_mesh, r_min, r_cut, m_min, m_max, preci, k, r_grid, m_grid)

  type(Atoms), intent(in), optional                        :: at
  real(dp), dimension(:), intent(out), allocatable         :: r_grid, m_grid
  real(dp), dimension(:), allocatable                      :: r_point, m_point 
  real(dp), intent(in)                                     :: r_cut, r_min, m_min, m_max
  real(dp)                                                 :: ivs(3)
  integer, intent(in)                                      :: r_mesh, m_mesh
  integer                                                  :: i, j
  integer, intent(in)                                      :: preci
  integer, intent(out)                                     :: k                             


  allocate(r_point(r_mesh))
  allocate(m_point(m_mesh))
 
  do i=1, r_mesh                                                                             ! r_mesh, the mesh size for r
     r_point(i) = r_min + real(i-1)*(r_cut - r_min)/real(r_mesh)
  enddo 

  do i=1, m_mesh
    m_point(i) = m_min + real(i-1)*(m_max - m_min)/real(m_mesh)
  enddo

  k = 0
  do i = 1, m_mesh     ! to know the value of k
     do j = 1, r_mesh
         ivs  = internal_vector(at, r_point(j),  m_point(i), 1)*SCALE_IVS
         if ((m_point(i) >= cutoff_m(r_cut, r_point(j), preci)) .and. (dot_product(ivs, ivs) > 0.0004_dp) )   k = k + 1
     enddo
  enddo

  allocate(r_grid(k))
  allocate(m_grid(k))

  k = 0
  do i=1, m_mesh 
     do j = 1, r_mesh 
         ivs  = internal_vector(at, r_point(j), m_point(i), 1)*SCALE_IVS
         write(*,*) "IVs using ",  r_point(j), m_point(i),"------", ivs                             
         write(*,*) "IVs_dr:   ",  internal_vector_dm(at, r_point(j), m_point(i), 1)*SCALE_IVS           
         write(*,*) "cutoff m:",   cutoff_m(r_cut, r_point(j), preci)
         if ((m_point(i) >= cutoff_m(r_cut, r_point(j), preci)) .and. (dot_product(ivs, ivs) > 0.0004_dp) )   then
             k = k + 1
             r_grid(k) = r_point(j) 
             m_grid(k) = m_point(i)  
         endif
     enddo
  enddo

  deallocate(r_point)
  deallocate(m_point)

 end subroutine grid_m_r0


 function cov(z_matrix1, z_matrix2, bent_space1, bent_space2, sigma, k)

    real(dp), intent(in)                  ::        z_matrix1(:,:), z_matrix2(:,:), bent_space1(:,:), bent_space2(:,:), sigma(:)
    real(dp)                              ::        cov, d_sq
    integer                               ::        i
    integer, intent(in)                   ::        k
 
    d_sq = 0.0_dp    
    do i=1, k
          d_sq = d_sq + (distance_bent_space(z_matrix1(i,:), z_matrix2(i,:), bent_space1, bent_space2) **2) /(sigma(i)**2) /2.0_dp
    enddo
    cov = exp(-1.0_dp * d_sq)

 endfunction cov


 function inverse_svd_threshold(in_matrix, n, thresh)
   real(dp), intent(in)                  ::         in_matrix(:,:), thresh
   integer,  intent(in)                  ::         n
   real(dp)                              ::         w(n), sigmainv(n,n), u(n,n), vt(n,n), inverse_svd_threshold(n,n)
   real(dp), allocatable                 ::         work(:)
   real(dp), parameter                   ::         TOL_SVD = 1e-13_dp
   integer                               ::         info, i, lwork, j
 
   call print('entering inverse_svd_threshold')

   if (n <= 3) then
      call print('in_matrix')
      call print(in_matrix)
   end if

   lwork = -1
   allocate(work(1))
   call DGESVD('A','A',n,n,in_matrix,n,w,u,n,vt,n, work, lwork, info)

   lwork = work(1) 
   deallocate(work)
   allocate(work(lwork))

   call DGESVD('A','A',n,n,in_matrix,n,w,u,n,vt,n, work, lwork, info)
   write (*,*) "DGESVD finished with exit code", info
   deallocate(work)
   if (info /= 0) then
          if (info < 0) then
                 write(line,'(a,i0)')'SVD: Problem with argument ',-info
                 call system_abort(line)
          else
                 call system_abort('SVD: DBDSQR (called from DGESVD) did not converge')
          end if
   end if

   sigmainv = 0.0_dp
   j = 0

   do i=1, n
       if (w(i) < thresh*TOL_SVD) then
            sigmainv(i,i) = 0.0_dp
            j = j + 1
       else
           sigmainv(i,i) = 1.0_dp/w(i)
       end if
   end do
   call print("the number of zero singular values:"//j)
   inverse_svd_threshold = transpose(vt) .mult. sigmainv .mult. transpose(u)

 endfunction inverse_svd_threshold

 
 function  distance_bent_space(vect1, vect2, bent_space1, bent_space2)

    real(dp), intent(in)                    :: vect1(3), vect2(3), bent_space1(:,:), bent_space2(:,:)
    real(dp)                                :: distance_bent_space

    distance_bent_space = norm( (bent_space1 .mult. vect1) - (bent_space2 .mult. vect2))  

 endfunction distance_bent_space

end program force_gaussian_prediction
