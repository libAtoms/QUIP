program force_gaussian_prediction

   use libatoms_module
   implicit none

   type(Atoms)                 :: at_in
   type(CInOutput)             :: in
   type(Dictionary)            :: params
   real(dp)                    :: r_cut, r_min, m_min, m_max, feature_len, theta, thresh, sigma_error, cutoff_len_ivs, sigma_factor, dist_primitive, distance_ivs_stati
   real(dp), dimension(:), allocatable          :: max_value, mean_value, deviation_value
   real(dp), dimension(:), allocatable          :: r_grid, m_grid, sigma, theta_array, covariance_pred, force_proj_ivs_pred, force_proj_target, distance_confs
   real(dp), parameter                          :: TOL_REAL=1e-7_dp, SCALE_IVS=100.0_dp
   real(dp)                                     :: force(3), feature_inv(3,3)
   real(dp), dimension(:,:,:), allocatable      :: feature_matrix_normalised, feature_matrix
   real(dp), dimension(:,:), allocatable        :: force_proj_ivs, force_proj_ivs_select, covariance, inv_covariance, distance_ivs
   real(dp), dimension(:,:), allocatable        :: feature_matrix_normalised_pred, feature_matrix_pred
   real(dp), dimension(:,:), pointer            :: force_ptr, force_ptr_mm
   integer                                      :: i,j, k, n, t, at_n, n_loop, preci, r_mesh, m_mesh, add_vector, n_relavant_confs
   integer, dimension(:), allocatable           :: distance_index
   logical                                      :: spherical_cluster_teach, spherical_cluster_pred, do_gp, fix_sigma,print_dist_stati
   character(STRING_LENGTH)                     :: teaching_file, grid_file, test_file
 
 

   call system_initialise(enable_timing=.true.)

   call initialise(params)
   call param_register(params, 'theta',  '2.5', theta, "the correlation length of the data")
   call param_register(params, 'r_min',  '0.5', r_min, "the minimum radius for the spherical atomic environment")
   call param_register(params, 'r_cut',  '8.0', r_cut, "the cutoff radius for the spherical atomic environment")
   call param_register(params, 'm_min',  '1.0', m_min, "the minimum m for calculating the atomic environment")
   call param_register(params, 'm_max',  '5.0', m_max, "the maxmium m for calculating the atomic environment")
   call param_register(params, 'cutoff_len_ivs', '0.04', cutoff_len_ivs, "the cutoff lenght for IVs to be considered valid when generating the grid")
   call param_register(params, 'thresh', '1.0', thresh, "the threshold for doing the Sigular Value Decompostion of the Covariance Matrix")
   call param_register(params, 'preci',  '6',   preci,  "the screening accuracy on the edge atoms")
   call param_register(params, 'add_vector', '0', add_vector, "the number of vectors you like to add into the internal_vector representation")
   call param_register(params, 'sigma_error', '0.01', sigma_error, "the noise assumed on the teaching data")
   call param_register(params, 'r_mesh', '6',   r_mesh, "grid finess of r0")
   call param_register(params, 'm_mesh', '6',   m_mesh, "grid finess of m")
   call param_register(params, 'do_gp',  'F', do_gp, "true for doing a gaussian processes, instead of SVD")
   call param_register(params, 'n_relavant_confs', '200', n_relavant_confs, "the number of relavant confs you would like to do machine learning with")
   call param_register(params, 'print_dist_stati', 'F', print_dist_stati, "set true to print out the distance on every single dimension of the IVs space")
   call param_register(params, 'spherical_cluster_teach', 'T', spherical_cluster_teach, "only the first atom in the cluster are considered when doing teaching")
   call param_register(params, 'spherical_cluster_pred', 'T', spherical_cluster_pred, "only the first atom in the cluster are considered when doing predicting")
   call param_register(params, 'fix_sigma',  'F', fix_sigma, "true, if you want manually input sigma")
   call param_register(params, 'sigma_factor', '1.0', sigma_factor, "the factor multiplied by sigma when using fix_sigma")
   call param_register(params, 'teaching_file', 'data.xyz', teaching_file, "file to read teaching configurations from")
   call param_register(params, 'grid_file', 'grid.xyz', grid_file, "file to generate the proper pairs of (r0, m)")
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

  if (add_vector > 0)   k = k + add_vector   
  call print('The number of valid  Internal Vectors :  '//k)
   
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
  call print('Got '//in%n_frame//' input frames.')

  allocate(feature_matrix(k,3,n))   
  allocate(feature_matrix_normalised(k,3, n))
  allocate(force_proj_ivs(k, n))

 t = 0 ! counter for loop the atomic configurations
 do i=1, in%n_frame
        call read(in, at_in)
        call assign_property_pointer(at_in, 'force', ptr=force_ptr)
        call assign_property_pointer(at_in, 'mm_force', ptr=force_ptr_mm)
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
          call print("atomic force in real space : "//norm(force))
          call print("The atom:"//at_n)

        do j=1, k-add_vector
           feature_matrix(j,:,t) = internal_vector(at_in, r_grid(j), m_grid(j), at_n)*SCALE_IVS
           feature_len = norm(feature_matrix(j,:, t))

           if (feature_len < TOL_REAL)  then
                   feature_len=1.0_dp
                   call print("WARNNING: TEACHING, encountered the decimal limit in getting the unit direction of IVs") 
           endif

           feature_matrix_normalised(j,:,t) = feature_matrix(j,:, t)/feature_len
           call print("internal vectors ( "//r_grid(j)//" "//m_grid(j)//" ) "//feature_matrix(j,1,t)//" "//feature_matrix(j,2,t)//" "//feature_matrix(j,3,t)) 
        enddo

        if (add_vector >0 ) then
           do j=k-add_vector+1, k
                feature_matrix(j,:,t) = force_ptr_mm(:,at_n)
                feature_len = norm(feature_matrix(j,:,t))

                if (feature_len < TOL_REAL)  then
                    feature_len=1.0_dp
                    call print("WARNNING: TEACHING, encountered the decimal limit in getting the unit direction of IVs")
                endif

                feature_matrix_normalised(j,:,t) = feature_matrix(j,:,t)/feature_len
                call print("added internal vectors :"//feature_matrix(j,1,t)//"  "//feature_matrix(j,2,t)//" "//feature_matrix(j,3,t))
           enddo
        endif

        ! store the force info for every configuration
        do j=1, k
             force_proj_ivs(j,t) = dot_product(feature_matrix_normalised(j,:,t), force) 
             call print("Force projection on IV"//j//" is:  "//force_proj_ivs(j,t))
        enddo   
     enddo   ! loop over the atoms
 enddo       ! loop over the frames

 call finalise(in)

 call print_title('starting the teaching process')

 allocate(distance_ivs(n,n))
 allocate(sigma(k))
 allocate(theta_array(k))

! the parameters for doing statistics on the distance matrix
call system_timer('Getting Hyper-parameters from the DATABASE')

if (fix_sigma) then
   if (.not. print_dist_stati) then
     open(unit=1, file='sigma.dat')
     read(1, *) sigma
     close(unit=1) 
     sigma = sigma * sigma_factor
   endif

   if (print_dist_stati) then
       allocate(mean_value(k))
       allocate(max_value(k))
       allocate(deviation_value(k))
       
       dist_primitive=0.0_dp
       do t = 1, k
          do i = 1, n
            do j=1, n
               distance_ivs(i,j) = distance_in_bent_space(feature_matrix(t,:,i), feature_matrix(t,:,j), feature_matrix_normalised(:,:,i), feature_matrix_normalised(:,:,j))
               if (i>j)  call print("Dimensional Distance : "//distance_ivs(i, j)//"  dimension : "//t)
            enddo   ! i
          enddo     ! j

          call matrix_statistics(distance_ivs, max_value(t), mean_value(t), deviation_value(t), n)
          call print("Statistics max value: "//max_value(t)//" mean value: "//mean_value(t)//" deviation_value: "//deviation_value(t))
          dist_primitive = dist_primitive + (mean_value(t)/deviation_value(t)/2.0_dp)**2 
       enddo

       dist_primitive = sqrt(dist_primitive / real(k,dp)) 
       call print("primitive_distance : "//dist_primitive)

       ! normalised distance with hyperparameters derived from statistics analysis      
       do i=1,n
          do j=1,n
          distance_ivs_stati = 0.0_dp   
            do t=1, k                 
              distance_ivs_stati = distance_ivs_stati + &
                  (distance_in_bent_space(feature_matrix(t,:,i), feature_matrix(t,:,j), feature_matrix_normalised(:,:,i), feature_matrix_normalised(:,:,j))/2.0_dp/deviation_value(t))**2  
            enddo
            distance_ivs_stati = sqrt(distance_ivs_stati/real(k,dp))          ! to take the square root and get a normalised distance 
            if (i>j)   call print("Normalisd Dimensional Distance : "//distance_ivs_stati//"  Normalised Value : "//distance_ivs_stati/dist_primitive)
          enddo  
       enddo  

       ! to write the derived sigma vector into file "sigma.dat" for later use
       sigma = deviation_value * sigma_factor   ! sigma derived from the statistical deviation on each single dimension of the IVs space
       open(unit=1, file='sigma.dat')
       write(1, *) sigma
       close(unit=1)
       deallocate(mean_value)
       deallocate(max_value)
       deallocate(deviation_value)
   endif                                                            ! do statistics on the distance matrix
else 
  theta_array=theta
  sigma(t) = maxval(distance_ivs(:,:))/theta_array(t)
  if  (abs(sigma(t))<TOL_REAL) then
              sigma(t)=TOL_REAL
              call print("WARNNING: SIGMA, encountered the decimal limit in getting Sigma from Theta")
  endif

  if (print_dist_stati) then
       allocate(mean_value(k))
       allocate(max_value(k))
       allocate(deviation_value(k))

   dist_primitive=0.0_dp
   do t = 1, k 
     do i = 1, n
         do j=1, n
         distance_ivs(i,j) = distance_in_bent_space(feature_matrix(t,:,i), feature_matrix(t,:,j), feature_matrix_normalised(:,:,i), feature_matrix_normalised(:,:,j)) 
         if (print_dist_stati) then  
            if (i>j)  call print("Dimensional Distance : "//distance_ivs(i, j)//"  dimension : "//t)
         endif
        enddo
    enddo

    ! statistically tackling the distance matrix
    call matrix_statistics(distance_ivs, max_value(t), mean_value(t), deviation_value(t), n)
    call print("Statistics max value: "//max_value(t)//" mean value: "//mean_value(t)//" deviation_value: "//deviation_value(t))
    dist_primitive = dist_primitive + ( mean_value(t)/deviation_value(t)/2.0_dp )**2
  enddo ! loop over the k dimension of ivs space

  dist_primitive = sqrt(dist_primitive/real(k,dp))
  call print("primitive distance : "//dist_primitive)

  ! normalised distance with statistical parameters       
  do i=1,n
       do j=1,n
         distance_ivs_stati = 0.0_dp
         do t=1, k
            distance_ivs_stati = distance_ivs_stati + &
               (distance_in_bent_space(feature_matrix(t,:,i), feature_matrix(t,:,j), feature_matrix_normalised(:,:,i), feature_matrix_normalised(:,:,j))/2.0_dp/deviation_value(t))**2
         enddo
         distance_ivs_stati = sqrt(distance_ivs_stati/real(k,dp))
         if (i>j)   call print("Normalisd Dimensional Distance : "//distance_ivs_stati//"  Normalised Value : "//distance_ivs_stati/dist_primitive)
       enddo
  enddo
  deallocate(mean_value)
  deallocate(max_value)
  deallocate(deviation_value)
 endif         ! doing print_dist_stati   
endif          ! doing fix_sigma or using theta

call print('sigma is:    '//sigma)
deallocate(distance_ivs)
deallocate(theta_array)

call system_timer('Getting hyper-parameters from the DATABASE')


call print_title('starting the predicting process')

 allocate(covariance_pred(n_relavant_confs))
 allocate(feature_matrix_pred(k,3))
 allocate(feature_matrix_normalised_pred(k,3))
 allocate(force_proj_ivs_pred(k))
 allocate(force_proj_target(k))
 allocate(covariance(n_relavant_confs,n_relavant_confs))
 allocate(inv_covariance(n_relavant_confs,n_relavant_confs))
 allocate(distance_confs(n_relavant_confs))
 allocate(distance_index(n_relavant_confs))
 allocate(force_proj_ivs_select(k, n_relavant_confs))
 

 call initialise(in, test_file, INPUT)
 do i=1, in%n_frame 
   call read(in, at_in)
   call assign_property_pointer(at_in, 'force', ptr=force_ptr)
   call assign_property_pointer(at_in, 'mm_force', ptr=force_ptr_mm)
   call set_cutoff(at_in, r_cut)
   call calc_connect(at_in)
   if (spherical_cluster_pred) then
        n_loop = 1
   else
        n_loop = at_in%N
   endif

   do at_n=1, n_loop

       do j= 1, k-add_vector
            feature_matrix_pred(j,:) = internal_vector(at_in, r_grid(j), m_grid(j), at_n)*SCALE_IVS
            call print("internal vectors ( "//r_grid(j)//" "//m_grid(j)//" ):   "//feature_matrix_pred(j,1)//"  "//feature_matrix_pred(j,2)//"  "//feature_matrix_pred(j,3))
            feature_len = norm(feature_matrix_pred(j,:))

            if (feature_len < TOL_REAL)  then
                   feature_len=1.0_dp
                   call print("WARNNING: PREDICTION, encountered the decimal limit in getting the unit direction of IVs")
            endif  

            feature_matrix_normalised_pred(j,:) = feature_matrix_pred(j,:)/feature_len
      enddo

      if (add_vector > 0) then
           do j = k-add_vector+1, k
              feature_matrix_pred(j,:)=force_ptr_mm(:, at_n)
              feature_len = norm(feature_matrix_pred(j,:))

              if (feature_len < TOL_REAL) then
                   feature_len=1.0_dp
                   call print("WARNNING: PREDICTION, encountered the decimal limit in getting the unit direction of IVs")
              endif

              feature_matrix_normalised_pred(j,:) = feature_matrix_pred(j,:)/feature_len
              call print("added internal vectors : "//feature_matrix_pred(j,1)//"  "//feature_matrix_pred(j,2)//"  "//feature_matrix_pred(j,3))
           enddo
      endif

      call system_timer('Sorting the DATABASE')

      ! do the sorting and selection
      call sorting_configuration(feature_matrix_pred, feature_matrix, feature_matrix_normalised_pred, feature_matrix_normalised, distance_confs, distance_index)
      call print("Min and Max DISTANCE with Index after Sorting: "//distance_confs(1)//" and "// &
                  distance_confs(n_relavant_confs)//"  the INDEX: "//distance_index(1)//" and "//distance_index(n_relavant_confs))

      call system_timer('Sorting the DATABASE')

       
      ! to establish the Covariance Matrix for DATABASE
      do t = 1, n_relavant_confs
         do j=1, n_relavant_confs
              covariance(t,j) = cov(feature_matrix(:,:,distance_index(t)), feature_matrix(:,:,distance_index(j)), &
                          feature_matrix_normalised(:,:,distance_index(t)), feature_matrix_normalised(:,:,distance_index(j)), sigma=sigma)
         enddo
      enddo

     !  doing Gaussian Processes, adding an noise "sigma_error" for the teaching data
     if (do_gp) then
     do t=1, n_relavant_confs
          covariance(t,t) = covariance(t,t) + sigma_error**2
     enddo
     endif


     if (do_gp) then
        call inverse(covariance, inv_covariance)
     else
        ! To Do Sigular Value Decomposition (SVD): A = U*SIGMA*VT
        inv_covariance = inverse_svd_threshold(covariance, n_relavant_confs, thresh)
     endif
     call print("MAX and MIN components of inverse_covariance : "//maxval(inv_covariance(2,:))//" "//minval(inv_covariance(2,:)))

      do t=1, n_relavant_confs
           force_proj_ivs_select(:,t)=force_proj_ivs(:,distance_index(t))
      enddo
      
      ! the prediction covariance

      do t= 1, n_relavant_confs
            covariance_pred(t) = cov(feature_matrix_pred, feature_matrix(:,:,distance_index(t)), feature_matrix_normalised_pred, &
                                                                                     feature_matrix_normalised(:,:,distance_index(t)), sigma=sigma)
      enddo


     force_proj_ivs_pred(:) = matmul(covariance_pred, matmul(inv_covariance, transpose(force_proj_ivs_select(:,:)) )) 

     do j=1, k
            force_proj_target(j) = dot_product(feature_matrix_normalised_pred(j,:), force_ptr(:,at_n))
     enddo   

   
     do j=1, k 
            call print("Force in IV space"//j//": "//force_proj_ivs_pred(j)//": "//force_proj_target(j)//": "//abs(force_proj_ivs_pred(j)-force_proj_target(j)))
            ! predicted force: real force: absolute difference 
     enddo

     ! using least-squares to restore the target force in the External Cartesian Space  
     call inverse(matmul(transpose(feature_matrix_normalised_pred), feature_matrix_normalised_pred), feature_inv)  
     force = feature_inv .mult. transpose(feature_matrix_normalised_pred) .mult. force_proj_ivs_pred
     call print("force in external space:"//force)
     call print("the original force:"//force_ptr(:, at_n))
     call print("max error :    "//maxval(abs(force_ptr(:,at_n)-force))//" norm  error :  "//norm(force_ptr(:,at_n)-force))

     if (do_gp) then
           call print("predicted error :"//( 1.0_dp + sigma_error - covariance_pred .dot. matmul(inv_covariance, covariance_pred)))
     endif  
 
    enddo  ! loop over postions
  enddo    ! loop over frames

 call finalise(in)

 deallocate(m_grid)
 deallocate(r_grid)
 deallocate(force_proj_ivs_pred)
 deallocate(feature_matrix_normalised)
 deallocate(feature_matrix_normalised_pred)
 deallocate(feature_matrix)
 deallocate(feature_matrix_pred)
 deallocate(sigma)
 deallocate(force_proj_ivs)
 deallocate(force_proj_ivs_select)
 deallocate(force_proj_target)
 deallocate(inv_covariance) 
 deallocate(covariance)
 deallocate(covariance_pred)
 deallocate(distance_confs)
 deallocate(distance_index)

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

     cutoff_m = log(real(preci,dp)*log(10.0_dp)) / log(r_cut/r0)
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
     r_point(i) = r_min + real(i-1,dp)*(r_cut - r_min)/real(r_mesh,dp)
  enddo 

  do i=1, m_mesh
    m_point(i) = m_min + real(i-1,dp)*(m_max - m_min)/real(m_mesh,dp)
  enddo

  k = 0
  do i = 1, m_mesh     ! to know the value of k
     do j = 1, r_mesh
         ivs  = internal_vector(at, r_point(j),  m_point(i), 1)*SCALE_IVS
         if ((m_point(i) >= cutoff_m(r_cut, r_point(j), preci)) .and. (norm(ivs) > cutoff_len_ivs) )  k = k + 1
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
         if ((m_point(i) >= cutoff_m(r_cut, r_point(j), preci)) .and. (norm(ivs) > cutoff_len_ivs))   then
             k = k + 1
             r_grid(k) = r_point(j) 
             m_grid(k) = m_point(i)  
         endif
     enddo
  enddo

  deallocate(r_point)
  deallocate(m_point)

 end subroutine grid_m_r0


 function cov(feature_matrix1, feature_matrix2, bent_space1, bent_space2, sigma, distance)

    real(dp), intent(in)                  ::        feature_matrix1(:,:), feature_matrix2(:,:), bent_space1(:,:), bent_space2(:,:) 
    real(dp), intent(in), optional        ::        sigma(:)
    real(dp), intent(out), optional       ::        distance
    real(dp)                              ::        cov, d_sq
    integer                               ::        i, k

    k=size(feature_matrix1(1,:))    
    if (present(sigma)) then
       d_sq = 0.0_dp    
       do i=1, k
               d_sq = d_sq + (distance_in_bent_space(feature_matrix1(i,:), feature_matrix2(i,:), bent_space1, bent_space2) **2) /(2.0_dp*sigma(i))**2
       enddo
    else                       ! if no sigma is provided, sigma=1.0_dp is assumed for every dimension
       d_sq = 0.0_dp
       do i=1, k
             d_sq = d_sq + (distance_in_bent_space(feature_matrix1(i,:), feature_matrix2(i,:), bent_space1, bent_space2) **2) /(2.0_dp*1.0_dp)**2
       enddo
    endif
    d_sq = d_sq/real(k, dp)      ! normalised by the dimensionality of the Internal Space

    if (present(distance))  distance=sqrt(d_sq)     
    cov = exp(-1.0_dp*d_sq)

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
   call print("DGESVD finished with exit code "//info)
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

 function  distance_in_bent_space(vect1, vect2, bent_space1, bent_space2)

    real(dp), intent(in)                    :: vect1(3), vect2(3), bent_space1(:,:), bent_space2(:,:)
    real(dp)                                :: distance_in_bent_space

    distance_in_bent_space = norm( (bent_space1 .mult. vect1) - (bent_space2 .mult. vect2))  
 
 endfunction distance_in_bent_space


 subroutine matrix_statistics(in_matrix, max_value, mean_value, deviation_value, n)

   real(dp), intent(in)                    ::  in_matrix(:,:)
   integer, intent(in)                     ::  n
   integer                                 ::  i, j, t
   real(dp), intent(out)                   ::  max_value, mean_value, deviation_value
   real(dp), allocatable                   ::  data_array(:)

   allocate(data_array(n*(n-1)/2))

    t=0
    do i=1, n
       do j=1,n
          if (i>j) then
             t=t+1
             data_array(t) = in_matrix(i,j)  
          endif
       enddo
    enddo
    mean_value = sum(data_array)/real(size(data_array),dp)
    max_value  = maxval(data_array)
    
    deviation_value = 0.0_dp
    do i=1, size(data_array)
       deviation_value = deviation_value + ( (data_array(i) - mean_value)**2 )
    enddo
    deviation_value = sqrt(deviation_value/real(size(data_array),dp))

    deallocate(data_array)

 end subroutine  matrix_statistics

 subroutine normal_squared_matrix(in_matrix, out_matrix) 

  real(dp), intent(in)                      ::  in_matrix(:,:)
  real(dp), intent(inout)                   ::  out_matrix(:,:)
  integer                                   ::  i 
   
  do i = 1, size(in_matrix(1,:))
       out_matrix(:,i) = in_matrix(:,i)/norm(in_matrix(:,i))
  enddo

 end subroutine normal_squared_matrix


  subroutine sorting_configuration(matrix_predict, matrix_data, matrix_predict_norm, matrix_data_norm, distance_confs, distance_index)
 
    real(dp), intent(in)                      ::  matrix_predict(:,:), matrix_data(:,:,:), matrix_predict_norm(:,:), matrix_data_norm(:,:,:)
    real(dp), intent(inout)                   ::  distance_confs(:)
    integer, intent(inout)                    ::  distance_index(:)
    real(dp)                                  ::  cov_tmp  
    integer                                   ::  i
    
    do i=1, size(matrix_data(1,1,:)) 
        cov_tmp=cov(matrix_predict, matrix_data(:,:,i), matrix_predict_norm(:,:), matrix_data_norm(:,:,i), distance=distance_confs(i))   
    enddo

    call insertion_sort(distance_confs, idx=distance_index)

  end subroutine sorting_configuration

end program force_gaussian_prediction
