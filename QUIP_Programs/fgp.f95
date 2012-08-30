module likelihood_optimisation_module

  use libAtoms_module
  implicit none

  type optimise_likelihood
      logical use_qr
      real(dp), dimension(:,:), pointer :: distance_matrix, covariance_matrix, force_covariance_matrix
  endtype optimise_likelihood

contains

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

      am%covariance_matrix = exp(-0.5_dp * am%distance_matrix/sigma_covariance**2)
      do i = 1, size(am%covariance_matrix,1)
         am%covariance_matrix(i,i) = am%covariance_matrix(i,i) + sigma_error**2
      enddo

      LA_covariance = am%covariance_matrix
      allocate(ICF(size(am%covariance_matrix,1),size(am%covariance_matrix,1)))

      if (am%use_qr) then
         call Matrix_QR_Solve(LA_covariance,am%force_covariance_matrix,ICF)
      else
         call Matrix_Solve(LA_covariance,am%force_covariance_matrix,ICF)
      end if

      likelihood = 1.5_dp * LA_Matrix_LogDet(LA_covariance) + 0.5_dp * trace(ICF)
!      write(*,*) "DEBUG Likelihood value: ", likelihood

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

      n = size(am%covariance_matrix,1)
      allocate(ICF(n,n),DCDS(n,n),IC(n,n),ICDCDS(n,n))

      am%covariance_matrix = exp(-0.5_dp * am%distance_matrix/sigma_covariance**2)
      DCDS = am%covariance_matrix * am%distance_matrix / sigma_covariance**3

      do i = 1, n
         am%covariance_matrix(i,i) = am%covariance_matrix(i,i) + sigma_error**2
      enddo
      
      LA_covariance = am%covariance_matrix

      if (am%use_qr) then
         call Matrix_QR_Solve(LA_covariance,am%force_covariance_matrix,ICF)
         call Matrix_QR_Solve(LA_covariance,DCDS,ICDCDS)
         call LA_Matrix_QR_Inverse(LA_covariance,IC)
      else
         call Matrix_Solve(LA_covariance,am%force_covariance_matrix,ICF)
         call Matrix_Solve(LA_covariance,DCDS,ICDCDS)
         call LA_Matrix_Inverse(LA_covariance,IC)
      end if

      dlikelihood(1) = sigma_error * (3*trace(IC) - sum(IC*ICF) )
      dlikelihood(2) = 1.5_dp * trace(ICDCDS) - 0.5_dp * sum(ICDCDS*transpose(ICF))
 !    write(*,*) "DEBUG dlikelihood", dlikelihood(1)**2, dlikelihood(2)**2
      
      call finalise(LA_covariance)
      if(allocated(ICF)) deallocate(ICF)
      if(allocated(DCDS)) deallocate(DCDS)
      if(allocated(IC)) deallocate(IC)
      if(allocated(ICDCDS)) deallocate(ICDCDS)

   endfunction dlikelihood

end module likelihood_optimisation_module


program force_gaussian_prediction

   use libatoms_module
   use likelihood_optimisation_module
   implicit none

   type(Atoms)           :: at_in
   type(CInOutput)       :: in
   type(Dictionary)      :: params
   real(dp)              :: r_cut, r_min, m_min, m_max, feature_len, theta, thresh, sigma_error, cutoff_len_ivs, sigma_covariance, dist_primitive, distance_ivs_stati
   real(dp), dimension(:), allocatable           :: max_value, mean_value, deviation_value
   real(dp), dimension(:), allocatable           :: r_grid, m_grid, sigma, theta_array, covariance_pred, force_proj_ivs_pred, force_proj_target, distance_confs
   real(dp), parameter                           :: TOL_REAL=1e-7_dp, SCALE_IVS=100.0_dp
   real(dp)                                      :: force(3), feature_inner_matrix(3,3), feature_inv(3,3), kappa
   real(dp), dimension(:,:,:), allocatable       :: feature_matrix_normalised, feature_matrix
   real(dp), dimension(:,:), allocatable         :: force_proj_ivs, force_proj_ivs_select, inv_covariance, distance_ivs, force_pred
   real(dp), dimension(:,:), allocatable, target :: covariance, distance_matrix, force_covariance_matrix, covariance_tiny
   real(dp), dimension(:,:), allocatable         :: feature_matrix_normalised_pred, feature_matrix_pred
   real(dp), dimension(:,:), pointer             :: force_ptr, force_ptr_mm, force_ptr_tff
   integer                                       :: i,j, k, n, t, n_center_atom, n_loop, preci, r_mesh, m_mesh 
   integer                                       :: add_vector, n_relevant_confs, func_type, selector, temp_integer, ii
   integer	                 		 :: local_ml_optim_size
   integer, dimension(:), allocatable            :: distance_index
   logical                                       :: spherical_cluster_teach, spherical_cluster_pred, do_gp, fix_sigma,print_dist_stati, least_sq, fixed_iv, print_ft_matrix_test,  print_at, linear_vectors, do_svd
   character(STRING_LENGTH)                      :: teaching_file, grid_file, test_file, out_file, iv_params_file
   
 
   call system_initialise(enable_timing=.true.)

   call initialise(params)
   call param_register(params, 'theta',  '2.5', theta, "the correlation length of the data")
   call param_register(params, 'r_min',  '0.5', r_min, "the minimum radius for the spherical atomic environment")
   call param_register(params, 'r_cut',  '8.0', r_cut, "the cutoff radius for the spherical atomic environment")
   call param_register(params, 'm_min',  '1.0', m_min, "the minimum m for calculating the atomic environment")
   call param_register(params, 'm_max',  '5.0', m_max, "the maximum m for calculating the atomic environment")
   call param_register(params, 'cutoff_len_ivs', '0.04', cutoff_len_ivs, "the cutoff lenght for IVs to be considered valid when generating the grid")
   call param_register(params, 'thresh', '1.0', thresh, "the threshold for doing the Singular Value Decompostion of the Covariance Matrix")
   call param_register(params, 'preci',  '6',   preci,  "the screening accuracy on the edge atoms")
   call param_register(params, 'add_vector', '0', add_vector, "the number of vectors you like to add into the internal_vector representation")
   call param_register(params, 'sigma_error', '0.01', sigma_error, "the noise assumed on the teaching data")
   call param_register(params, 'r_mesh', '6',   r_mesh, "grid finess of r0")
   call param_register(params, 'm_mesh', '6',   m_mesh, "grid finess of m")
   call param_register(params, 'func_type', '0', func_type, "which kernel function is used to build the covariance matrix")
   call param_register(params, 'do_gp',  'F', do_gp, "true for doing a gaussian processes")
   call param_register(params, 'n_relevant_confs', '200', n_relevant_confs, "the number of relevant confs you would like to do machine learning with")
   call param_register(params, 'print_dist_stati', 'F', print_dist_stati, "set true to print out the distance on every single dimension of the IVs space")
   call param_register(params, 'spherical_cluster_teach', 'T', spherical_cluster_teach, "only the first atom in the cluster are considered when doing teaching")
   call param_register(params, 'spherical_cluster_pred', 'T', spherical_cluster_pred, "only the first atom in the cluster are considered when doing predicting")
   call param_register(params, 'fix_sigma',  'F', fix_sigma, "true, if you want manually input sigma")
   call param_register(params, 'least_sq',    'T', least_sq, "if true, the internal force components will be tranformed to real force using least squares")
   call param_register(params, 'do_svd', 'F', do_svd, "if true, doing inverting by SVD")
   call param_register(params, 'sigma_covariance', '1.0', sigma_covariance, "the factor multiplied by sigma when using fix_sigma")
   call param_register(params, 'teaching_file', 'data.xyz', teaching_file, "file to read teaching configurations from")
   call param_register(params, 'grid_file', 'grid.xyz', grid_file, "file to generate the proper pairs of (r0, m)")
   call param_register(params, 'test_file', 'test.xyz', test_file, "file to read the testing configurations from")
   call param_register(params, 'out_file', 'out.xyz', out_file, "output configuration file containing the predicted forces")
   call param_register(params, 'print_at', 'F', print_at, "true for print out the testing configuration with predicted forces")
   call param_register(params, 'print_ft_matrix_test', 'F', print_ft_matrix_test, "print the feature matrix of the testing conf if true")
   call param_register(params, 'fixed_iv', 'F', fixed_iv, "does a file containing the internal vector parameters exist?")
   call param_register(params, 'iv_params_file', 'iv_params.csv', iv_params_file, "location of internal vectors parameters file (if necessary)")
   call param_register(params, 'local_ml_optim_size', '0', local_ml_optim_size, "Optimise sigma_error and sigma_covariance using local_ml_optim_size LS confs. If 0, no optimisation is performed")
   
   if (.not. param_read_args(params, task="fgp command line params")) then
      call print("Usage: fgp [options]")
      call system_abort("Confused by command line arguments")
   end if

   call print_title('params')
   call param_print(params)
   call finalise(params)

   if (fixed_iv) then
      call load_iv_params(iv_params_file, r_grid, m_grid, k) 
   else   
      call print_title('Testing the grid')
      ! generate the r0_m grids
      call initialise(in, grid_file, INPUT)
      call read(in, at_in)
      ! including the tranlational symmetry images
      call set_cutoff(at_in, r_cut)
      call calc_connect(at_in)
      call  grid_m_r0(at_in, r_mesh, m_mesh, r_min, r_cut, m_min, m_max, preci, k, r_grid, m_grid, cutoff_len_ivs)
      call  finalise(in)
   end if
   
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

  call system_timer('Preparing Information from the DATABASE')
 
  ! the main work starts
  call initialise(in, teaching_file, INPUT)
  call print('Got '//in%n_frame//' input frames.')

  allocate(feature_matrix(k, 3, n))   
  allocate(feature_matrix_normalised(k, 3, n))
  allocate(force_proj_ivs(k, n))

 t = 0 ! counter for loop the atomic configurations
 do i=1, in%n_frame
        call read(in, at_in)
        call assign_property_pointer(at_in, 'force', ptr=force_ptr)
        if (add_vector > 0) call assign_property_pointer(at_in, 'mm_force', ptr=force_ptr_mm) 
        ! the first vector is labelled mm_force in the input file, and is usually the SW force
        if (add_vector > 1) call assign_property_pointer(at_in, 'tersoff_f', ptr=force_ptr_tff) 
        ! the second vector is labelled tersoff_f
        ! call print('the frame: '//i)

        ! Including the symmetry images
        call set_cutoff(at_in, r_cut)
        call calc_connect(at_in)

     if (spherical_cluster_teach) then        
        n_loop = 1
     else
        n_loop = at_in%N   
     endif

    do n_center_atom =1,  n_loop
          t = t + 1
          force = force_ptr(:,n_center_atom)
          !call print("atomic force in real space : "//norm(force))
          !call print("The atom:"//n_center_atom)

        do j=1, k-add_vector
           feature_matrix(j,:,t) = internal_vector(at_in, r_grid(j), m_grid(j), n_center_atom)*SCALE_IVS
           feature_len = norm(feature_matrix(j,:, t))

           if (feature_len < TOL_REAL)  then
                   feature_len=1.0_dp
                   call print("WARNING: TEACHING, encountered the decimal limit in getting the unit direction of IVs") 
           endif

           feature_matrix_normalised(j,:,t) = feature_matrix(j,:, t)/feature_len
           !call print("internal vectors ( "//r_grid(j)//" "//m_grid(j)//" ) "//feature_matrix(j,1,t)//" "//feature_matrix(j,2,t)//" "//feature_matrix(j,3,t)) 
        enddo

        if (add_vector > 0 ) then
           do j=k-add_vector+1, k
              temp_integer=k-j
              select case (temp_integer)
              case (0)
                 feature_matrix(j,:,t) = force_ptr_mm(:,n_center_atom)
              case (1)
                 feature_matrix(j,:,t) = force_ptr_tff(:,n_center_atom)
              end select
              
              feature_len = norm(feature_matrix(j,:,t))
              
              if (feature_len < TOL_REAL)  then
                 feature_len=1.0_dp
                 call print("WARNING: TEACHING, encountered the decimal limit in getting the unit direction of IVs")
              endif
              
              feature_matrix_normalised(j,:,t) = feature_matrix(j,:,t)/feature_len
              !call print("added internal vectors :"//feature_matrix(j,1,t)//"  "//feature_matrix(j,2,t)//" "//feature_matrix(j,3,t))
           enddo
        endif

        ! store the force info for every configuration
        do j=1, k
             force_proj_ivs(j,t) = dot_product(feature_matrix_normalised(j,:,t), force) 
             !call print("Force projection on IV"//j//" is:  "//force_proj_ivs(j,t))
        enddo   
     enddo   ! loop over the atoms
 enddo       ! loop over the frames

 call finalise(in)

 call system_timer('Preparing Information from the DATABASE')

 call print_title('starting the teaching process')

 allocate(distance_ivs(n,n))
 allocate(sigma(k))
 allocate(theta_array(k))
 allocate(mean_value(k))
 allocate(max_value(k))
 allocate(deviation_value(k))

! the parameters for doing statistics on the distance matrix
call system_timer('Getting Hyper-parameters from the DATABASE')

if (fix_sigma) then
   if (.not. print_dist_stati) then
      open(unit=1, file='sigma.dat')
      read(1, *) sigma
      close(unit=1) 
   endif
   
   if (print_dist_stati) then
   
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
            if (i>j)   call print("Normalised Dimensional Distance : "//distance_ivs_stati//"  Normalised Value : "//distance_ivs_stati/dist_primitive)
         enddo
      enddo

      ! to write the derived sigma vector into file "sigma.dat" for later use
      sigma = deviation_value           ! sigma derived from the statistical deviation on each single dimension of the IVs space
      open(unit=1, file='sigma.dat')
      write(1,*) sigma
      close(unit=1)
   endif                                !  statistics on the distance matrix

else 
   theta_array=theta
   sigma(t) = maxval(distance_ivs(:,:))/theta_array(t)
   if  (abs(sigma(t))<TOL_REAL) then
      sigma(t)=TOL_REAL
      call print("WARNING: SIGMA, encountered the decimal limit in getting Sigma from Theta")
   endif
      
endif          ! doing fix_sigma or using sigma derived from theta

 call print('sigma is:    '//sigma)

deallocate(distance_ivs)
deallocate(theta_array)
deallocate(mean_value)
deallocate(max_value)
deallocate(deviation_value)

call system_timer('Getting hyper-parameters from the DATABASE')


call print_title('starting the predicting process')

if (n_relevant_confs > n)  call system_abort("Trying to get more confs than the database has")

allocate(covariance_pred(n_relevant_confs))
allocate(feature_matrix_pred(k,3))
allocate(feature_matrix_normalised_pred(k,3))
allocate(force_proj_ivs_pred(k))
allocate(force_proj_target(k))
allocate(covariance(n_relevant_confs,n_relevant_confs)) 
if (local_ml_optim_size > 0) then 
   allocate(distance_matrix(local_ml_optim_size,local_ml_optim_size), force_covariance_matrix(local_ml_optim_size,local_ml_optim_size), covariance_tiny(local_ml_optim_size,local_ml_optim_size))
endif
allocate(inv_covariance(n_relevant_confs,n_relevant_confs))
allocate(distance_confs(n))
allocate(distance_index(n))
allocate(force_proj_ivs_select(k, n_relevant_confs))

call initialise(in, test_file, INPUT)
do i=1, in%n_frame 
   call read(in, at_in)
   call assign_property_pointer(at_in, 'force', ptr=force_ptr)
   if (add_vector > 0) call assign_property_pointer(at_in, 'mm_force', ptr=force_ptr_mm)
   if (add_vector > 1) call assign_property_pointer(at_in, 'tersoff_f', ptr=force_ptr_tff)

   call set_cutoff(at_in, r_cut)
   call calc_connect(at_in)
   if (spherical_cluster_pred) then
      n_loop = 1
   else
      n_loop = at_in%N
   endif

   do n_center_atom=1, n_loop
     do j= 1, k-add_vector
        feature_matrix_pred(j,:) = internal_vector(at_in, r_grid(j), m_grid(j), n_center_atom)*SCALE_IVS
        call print("internal vectors ( "//r_grid(j)//" "//m_grid(j)//" ):   "//feature_matrix_pred(j,1)//"  "//feature_matrix_pred(j,2)//"  "//feature_matrix_pred(j,3))
        feature_len = norm(feature_matrix_pred(j,:))

        if (feature_len < TOL_REAL)  then
            feature_len=1.0_dp
            call print("WARNING: PREDICTION, encountered the decimal limit in getting the unit direction of IVs")
        endif

        feature_matrix_normalised_pred(j,:) = feature_matrix_pred(j,:)/feature_len
     enddo

      if (add_vector > 0) then
         do j = k-add_vector+1, k
            temp_integer=k-j
            select case (temp_integer)
            case (0)
               feature_matrix_pred(j,:)=force_ptr_mm(:, n_center_atom)
            case (1)
               feature_matrix_pred(j,:)=force_ptr_tff(:, n_center_atom)
            end select
            feature_len = norm(feature_matrix_pred(j,:))

            if (feature_len < TOL_REAL) then
               feature_len=1.0_dp
               call print("WARNING: PREDICTION, encountered the decimal limit in getting the unit direction of IVs")
            endif
            
            feature_matrix_normalised_pred(j,:) = feature_matrix_pred(j,:)/feature_len
            call print("added internal vectors : "//feature_matrix_pred(j,1)//"  "//feature_matrix_pred(j,2)//"  "//feature_matrix_pred(j,3))
         enddo
      endif


     if (print_ft_matrix_test) then
       do j=1, k 
          do ii=1, k   
            call print("ATOM : "//n_center_atom//" feature_matrix: ("//j//" "//ii//") "//dot_product(feature_matrix_pred(j,:), feature_matrix_normalised_pred(ii,:)) )
          enddo   
       enddo
     endif

      
      call system_timer('Sorting the DATABASE')
      
      ! do the sorting and selection
      call sorting_configuration(feature_matrix_pred, feature_matrix, feature_matrix_normalised_pred, feature_matrix_normalised, sigma, distance_confs, distance_index)
      ! call print("Min and Max DISTANCE with Index after Sorting: "//distance_confs(1)//" and "// &
      !     distance_confs(n_relevant_confs)//"  the INDEX: "//distance_index(1)//" and "//distance_index(n_relevant_confs))
   
      if (distance_confs(1)>0.2_dp) then
                sigma_covariance = sigma_covariance*distance_confs(1)/0.2_dp     ! distance scaling
                write(*,*) "distance are scaled with a factor of 0.2"
      endif

      if (.false.) then        ! massive output, only for testing use
         call print("Detail on the teaching procedure")
         do ii = 1, n_relevant_confs
            call print("Index: "//distance_index(ii)//" Distance: "//distance_confs(ii))
            call print("Force in IVs space: "//force_proj_ivs(:,distance_index(ii)))
         enddo
      endif
      call system_timer('Sorting the DATABASE')

      if (local_ml_optim_size > 0) call print("Optimise: "//local_ml_optim_size)
      !!! ML optimisation here: first calculates the necessary input matrices, then optimises
      if (local_ml_optim_size > 0) then
         call print("Local Maximum Likelihood selected")
         do t = 1, local_ml_optim_size  ! loop to generate the force_cov_mat and the distance_matrix of the restricted subset of the LS
            do j = t, local_ml_optim_size
               force_covariance_matrix(t,j) = dot_product(force_proj_ivs(:,distance_index(t)), force_proj_ivs(:,distance_index(j))) 
               ! check if possible, may be a wreck because IVs make absolute representation of the force impossible 
               force_covariance_matrix(j,t) = force_covariance_matrix(t,j)
               covariance_tiny(t,j) = cov(feature_matrix(:,:,distance_index(t)), feature_matrix(:,:,distance_index(j)), &
                    feature_matrix_normalised(:,:,distance_index(t)), feature_matrix_normalised(:,:,distance_index(j)), sigma, sigma_covariance, func_type=func_type, &
                    distance=distance_matrix(t,j))
               covariance_tiny(j,t) = covariance_tiny(t,j)
               distance_matrix(j,t) = distance_matrix(t,j)
            enddo
         enddo
         call do_optimise_likelihood(sigma_error, sigma_covariance) 
      endif
      
      
      ! Generate covariance matrix
      do t = 1, n_relevant_confs  ! loop to generate the covariance matrix of the learning set
         covariance(t,t) = cov(feature_matrix(:,:,distance_index(t)), feature_matrix(:,:,distance_index(t)), &
                feature_matrix_normalised(:,:,distance_index(t)), feature_matrix_normalised(:,:,distance_index(t)), sigma, sigma_covariance, func_type=func_type)
         if (do_gp) covariance(t,t) = covariance(t,t) + sigma_error**2
         do j = t+1, n_relevant_confs ! above diagonal
            covariance(t,j) = cov(feature_matrix(:,:,distance_index(t)), feature_matrix(:,:,distance_index(j)), &
                 feature_matrix_normalised(:,:,distance_index(t)), feature_matrix_normalised(:,:,distance_index(j)), sigma, sigma_covariance, func_type=func_type)
            covariance(j,t) = covariance(t,j)  
        enddo
      enddo
       
      call system_timer('Inverting the Covariance Matrix')
      
      if (do_gp) then
         call inverse(covariance, inv_covariance)
      else
         ! To Do Sigular Value Decomposition (SVD): A = U*SIGMA*VT
         inv_covariance = inverse_svd_threshold(covariance, n_relevant_confs, thresh)
      endif
      
      call system_timer('Inverting the Covariance Matrix')
      

      do t=1, n_relevant_confs
         force_proj_ivs_select(:,t)=force_proj_ivs(:,distance_index(t))
      enddo
      
      ! the test configuration covariance vector
      
      do t= 1, n_relevant_confs
        covariance_pred(t) = cov(feature_matrix_pred, feature_matrix(:,:,distance_index(t)), feature_matrix_normalised_pred, &
                                         feature_matrix_normalised(:,:,distance_index(t)), sigma, sigma_covariance, func_type=func_type)
      enddo
      
      ! the predicted force projected on each of the internal directions
      force_proj_ivs_pred(:) = matmul(covariance_pred, matmul(inv_covariance, transpose(force_proj_ivs_select(:,:)) )) 

      do j=1, k
         force_proj_target(j) = dot_product(feature_matrix_normalised_pred(j,:), force_ptr(:,n_center_atom))
      enddo
      
      
      do j=1, k 
         call print("Force in IV space"//j//": "//force_proj_ivs_pred(j)//": "//force_proj_target(j)//": "//abs(force_proj_ivs_pred(j)-force_proj_target(j)))  !//&
             ! " : "//(abs(force_proj_ivs_pred(j)-force_proj_target(j))/abs(force_proj_target(j))))   
         ! predicted force: real force: absolute difference 
      enddo
      
      ! using least-squares to restore the target force in the External Cartesian Space
      feature_inner_matrix=matmul(transpose(feature_matrix_normalised_pred), feature_matrix_normalised_pred)

      if (do_svd) then
            write(*,*)  "feature inner matrix :", feature_inner_matrix
            feature_inv = inverse_svd_threshold(feature_inner_matrix, 3, thresh)
      else
            write(*,*)  "feature_inner_matrix :", feature_inner_matrix
            call inverse(feature_inner_matrix, feature_inv) 
      endif

      if (least_sq) then 
         force = feature_inv .mult. transpose(feature_matrix_normalised_pred) .mult. force_proj_ivs_pred
      else
         call real_expection_sampling_components(feature_matrix_normalised_pred, force_proj_ivs_pred, force)        
      endif 


      ! a trick to cope with symmetry problem 
      call internal_vector_linearity(feature_matrix_normalised_pred, linear_vectors, 0.1_dp)
      if (linear_vectors) then              
            force=(transpose(feature_matrix_normalised_pred) .mult. force_proj_ivs_pred )/real(k, dp)
      endif

      call print("force in external space:"//force)
      call print("the original force:"//force_ptr(:, n_center_atom))
      call print("max error :    "//maxval(abs(force_ptr(:,n_center_atom)-force))//" norm  error :  "//norm(force_ptr(:,n_center_atom)-force))

      if (print_at) force_ptr(:,n_center_atom)=force     

      kappa = cov(feature_matrix_pred, feature_matrix_pred, feature_matrix_normalised_pred, feature_matrix_normalised_pred, sigma, sigma_covariance, func_type=func_type) + sigma_error**2
    
      call print("predicted error : "//sqrt(abs(kappa - covariance_pred .dot. matmul(inv_covariance, covariance_pred))))
     
   enddo  ! loop over positions

   if (print_at) then 
      call write(at_in, out_file, append=.true.)      ! write the output configurations with predicted force
   endif

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
if(allocated(force_covariance_matrix)) deallocate(force_covariance_matrix)
if(allocated(covariance_tiny)) deallocate(covariance_tiny)
if(allocated(distance_matrix)) deallocate(distance_matrix)

call system_finalise()
 
contains 
   
  function internal_vector(at, r0, m, n_center_atom)

      real(dp)                            :: internal_vector(3), delta_r(3), delta_r_len
      real(dp), intent(in)                :: r0, m
      type(Atoms), intent(in)             :: at
      integer                             :: i, j, n_center_atom

      internal_vector = 0.0_dp
      do i=1, atoms_n_neighbours(at, n_center_atom)
          j = atoms_neighbour(at, n_center_atom, i, distance=delta_r_len, diff=delta_r) 
          internal_vector = internal_vector + (delta_r *exp(-((delta_r_len/r0)**m)) /delta_r_len)
      enddo
 endfunction  internal_vector


 function internal_vector_dm(at, r0, m, n_center_atom)

      real(dp)                            :: internal_vector_dm, internal_vector_prim(3), delta_r(3), delta_r_len, alpha
      real(dp), intent(in)                :: r0, m
      type(Atoms), intent(in)             :: at
      integer                             :: i, j, n_center_atom

      internal_vector_prim = 0.0_dp
      do i=1, atoms_n_neighbours(at, n_center_atom)
          j = atoms_neighbour(at, n_center_atom, i, distance=delta_r_len, diff=delta_r)
          alpha = delta_r_len / r0
          internal_vector_prim = internal_vector_prim - (delta_r * exp(-(alpha**m)) *(alpha**m) *log(alpha) /delta_r_len)
      enddo
      internal_vector_dm = norm(internal_vector_prim)
 endfunction  internal_vector_dm


 function cutoff_m(r_cut, r0, preci)
 
     real(dp)          :: r0, r_cut, cutoff_m
     integer           :: preci

     cutoff_m = log(real(preci,dp)*log(10.0_dp)) / log(r_cut/r0)
 endfunction cutoff_m


 subroutine load_iv_params(iv_params_file, r_grid, m_grid, k)
   ! reads r and m from iv_params file. The first line must contain the number
   ! of internal vector, and the following lines should be formatted in two columns with commas as spacers (a common csv file)
   character(STRING_LENGTH), intent(in)                     :: iv_params_file                           
   real(dp), dimension(:), intent(out), allocatable         :: r_grid, m_grid
   integer                                                  :: i
   integer, intent(out)                                     :: k 
   
   open (unit=22, file=iv_params_file, status='old', action='read')
   read(22,*), k
   write(*,*) "Number of iv: ", k
   allocate(r_grid(k))
   allocate(m_grid(k))
   
   write(*,*) "Reading internal vectors parameters from file ", iv_params_file
   
   do i=1, k
      read(22,*) r_grid(i), m_grid(i)
      write(*,*) "Vector", i, ":", r_grid(i), m_grid(i)
   end do

   close(22)
   
 end subroutine load_iv_params

 subroutine grid_m_r0(at, r_mesh, m_mesh, r_min, r_cut, m_min, m_max, preci, k, r_grid, m_grid, cutoff_len_ivs)

  type(Atoms), intent(in), optional                        :: at
  real(dp), dimension(:), intent(out), allocatable         :: r_grid, m_grid
  real(dp), dimension(:), allocatable                      :: r_point, m_point 
  real(dp), intent(in)                                     :: r_cut, r_min, m_min, m_max, cutoff_len_ivs
  real(dp)                                                 :: ivs(3)
  real(dp), parameter                                      :: SCALE_IVS = 100.0_dp
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


 function cov(feature_matrix1, feature_matrix2, bent_space1, bent_space2, iv_weights, sigma_cov, distance, func_type)

    real(dp), intent(in)                  ::        feature_matrix1(:,:), feature_matrix2(:,:), bent_space1(:,:), bent_space2(:,:), iv_weights(:), sigma_cov
    real(dp), intent(out), optional       ::        distance
    real(dp)                              ::        cov, d_sq
    integer, intent(in), optional         ::        func_type
    integer                               ::        i, k
    
    k = size(feature_matrix1(:,1))    
    d_sq = 0.0_dp    
    do i=1, k
       d_sq = d_sq + (distance_in_bent_space(feature_matrix1(i,:), feature_matrix2(i,:), bent_space1, bent_space2))**2 / (2.0_dp * (iv_weights(i)**2))
    enddo
 
    ! normalised with the dimensionality k of the Internal Space
    d_sq = d_sq/real(k,dp)                            

    if (present(distance))  distance = d_sq  ! N.B. distance squared

    if (present(func_type)) then
         selector=func_type
    else
         selector=0
    endif

    SELECT CASE (selector)
    CASE (0)
       cov = exp(-0.5_dp*d_sq/sigma_cov**2)
    CASE (1)
       cov = sqrt(d_sq)
    CASE (2)
       cov = d_sq
    CASE (3)
       cov = (sqrt(d_sq))**3
    CASE (4)
       cov = d_sq**2
    END SELECT

 endfunction cov

 function inverse_svd_threshold(in_matrix, n, thresh)

   real(dp), intent(in)                  ::         in_matrix(:,:), thresh
   integer,  intent(in)                  ::         n
   real(dp), allocatable                 ::         w(:), sigmainv(:,:), u(:,:), vt(:,:), inverse_svd_threshold(:,:)
   real(dp), allocatable                 ::         work(:)
   real(dp), parameter                   ::         TOL_SVD = 1e-13_dp
   integer                               ::         info, i, lwork, j
 
   call print('entering inverse_svd_threshold')

   allocate(w(n), sigmainv(n,n), u(n,n), vt(n,n))
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
   call print("the number of zero singular values : "//j)
   inverse_svd_threshold = transpose(vt) .mult. sigmainv .mult. transpose(u)
   deallocate(w)
   deallocate(sigmainv)
   deallocate(u)
   deallocate(vt)

 endfunction inverse_svd_threshold

 function  distance_in_bent_space(vect1, vect2, bent_space1, bent_space2)

    real(dp), intent(in)                    :: vect1(3), vect2(3), bent_space1(:,:), bent_space2(:,:)
    real(dp)                                :: distance_in_bent_space

    distance_in_bent_space = norm( (bent_space1 .mult. vect1) - (bent_space2 .mult. vect2))  
 
 endfunction distance_in_bent_space

! function eval_lhood(covariance_matrix,FCM)
!   real(dp), dimension(:,:), intent(in) :: covariance_matrix, FCM
!   real(dp) :: eval_lhood   
!   real(dp), dimension(:,:), allocatable :: ICF
!   type(LA_Matrix) :: LA_covariance
      
!   LA_covariance = covariance_matrix
!   allocate(ICF(size(FCM,1),size(FCM,1)))
   
!   call Matrix_QR_Solve(LA_covariance, FCM, ICF)
   !write(*,*) "DEBUG ICF: ", ICF
      
!   eval_lhood = 1.5_dp * LA_Matrix_LogDet(LA_covariance) + 0.5_dp * trace(ICF)
   
!   call finalise(LA_covariance)
!   if(allocated(ICF)) deallocate(ICF)
   
! endfunction eval_lhood


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


 subroutine sorting_configuration(matrix_predict, matrix_data, matrix_predict_norm, matrix_data_norm, sigma, distance_confs, distance_index)
 
   ! returns a list of indexes corresponding to the configurations ordered according to their distance to target 
   real(dp), intent(in)                        ::  matrix_predict(:,:), matrix_data(:,:,:), matrix_predict_norm(:,:), matrix_data_norm(:,:,:), sigma(:)
   real(dp), intent(inout)                     ::  distance_confs(:)
   integer,  intent(inout)                     ::  distance_index(:)
   real(dp)                                    ::  cov_tmp  
   integer                                     ::  i
 
   do i=1, size(matrix_data(1,1,:)) 
      cov_tmp = cov(matrix_predict, matrix_data(:,:,i), matrix_predict_norm, matrix_data_norm(:,:,i), sigma, 1.0_dp, distance=distance_confs(i))   
   enddo

    call insertion_sort(distance_confs, idx=distance_index)

 end subroutine sorting_configuration


 subroutine  real_expection_sampling_components(ivs_direction_matrix, internal_component_matrix, target_force)

    real(dp), intent(in)                        ::  ivs_direction_matrix(:,:), internal_component_matrix(:)
    real(dp), intent(out)                       ::  target_force(3)
    real(dp)                                    ::  converting_matrix(3,3), converting_matrix_inv(3,3), force_value_matrix(3), const
    integer                                     ::  i, j, t, n_counter, s

    n_counter=0

    target_force=0.0_dp

    s=size(internal_component_matrix)
    write(*,*) s

    do i=1, s
       do j=1, s
          do t=1, s
             if ((i>j) .and. (j>t)) then
               converting_matrix(:,1)=ivs_direction_matrix(i,:)
               converting_matrix(:,2)=ivs_direction_matrix(j,:)
               converting_matrix(:,3)=ivs_direction_matrix(t,:)
               force_value_matrix(1)=internal_component_matrix(i)
               force_value_matrix(2)=internal_component_matrix(j)
               force_value_matrix(3)=internal_component_matrix(t)

               ! a criteria needed here to decide whether or not to do inversing 
               const=0.01_dp

               if ((norm(converting_matrix(:,1) .cross. converting_matrix(:,2)) > const) .and. (norm(converting_matrix(:,2) &
                           .cross. converting_matrix(:,3)) > const) .and. (norm(converting_matrix(:,3) .cross. converting_matrix(:,1)) > const)) then
                  write(*,*) "the converting matrix", converting_matrix
                  converting_matrix_inv=inverse_svd_threshold(converting_matrix, 3, 10.0_dp)
                  target_force=transpose(converting_matrix_inv) .mult. force_value_matrix
                  write(*,*) "the inverse matrix : ", converting_matrix_inv
                  write(*,*) "the number of sequence of ivs:", i, j, t
                  if (.true.) call print("Target Force Distribution : "//target_force(1)//" "//target_force(2)//" "//target_force(3))
                  n_counter = n_counter + 1
                endif     ! orientational condition
             endif        ! condition of i<j<t
           enddo
       enddo
     enddo

  target_force=target_force / real(n_counter,dp)

 end subroutine real_expection_sampling_components


 subroutine internal_vector_linearity(in_matrix, linear_vectors, TOL)

   real(dp), intent(in)                        ::  in_matrix(:,:), TOL
   real(dp)                                    ::  x_sin
   logical, intent(out)                        ::  linear_vectors
   integer                                     ::  k_size, i, j
   k_size = size(in_matrix(:,1))

   x_sin =0.0_dp
   linear_vectors=.true.

   do i=1, k_size-1
      do j=i+1, k_size
         x_sin = norm(in_matrix(i,:) .cross. in_matrix(j,:))
         if ( 1.0_dp-TOL > x_sin > TOL) linear_vectors=.false.
      enddo
   enddo

   write(*,*) "the threshold for being linear : ", TOL, linear_vectors
  
  end subroutine internal_vector_linearity
 

 subroutine do_optimise_likelihood(sigma_error, sigma_covariance)

   type(optimise_likelihood)           :: am_likelihood 
   logical                             :: use_qr
   integer                             :: am_data_size
   character(len=1), dimension(1)      :: am_mold
   character(len=1), allocatable       :: am_data(:)
   real(dp)                            :: x_sigma(2), covariance_tol
   real(dp), intent(inout)             :: sigma_error, sigma_covariance
   
   covariance_tol = 1.0e-2
   use_qr=.true.

   am_likelihood%use_qr = use_qr ! filling the am_likelihood container
   am_likelihood%distance_matrix => distance_matrix
   am_likelihood%covariance_matrix => covariance_tiny
   am_likelihood%force_covariance_matrix => force_covariance_matrix

   am_data_size = size(transfer(am_likelihood,am_mold))
   allocate(am_data(am_data_size))
   am_data = transfer(am_likelihood,am_mold) ! pouring the content of am_likelihood into am_data
     
   ! Testing the implementation of likelihood and dlikelihood.
   if(.false.) then
      call print('starting test gradient')
      call verbosity_push(PRINT_VERBOSE)
      write (*,*) test_gradient((/sigma_error,sigma_covariance/),likelihood,dlikelihood,data=am_data)
      call verbosity_pop()
   endif
   call print_title('likelihood optimisation')

   if (use_qr) then
      call print('doing likelihood optimisation with QR decomposition')
   else
      call print('doing likelihood optimisation with Cholesky decomposition')
   end if
   x_sigma = (/sigma_error, sigma_covariance/)    ! starting points for the values
   ! n = minim( x_sigma, likelihood, dlikelihood, method='cg', convergence_tol=covariance_tol, max_steps = 100, data=am_data)
   ! x_sigma is the input vector, lh is the target function to minimise (using dlh), using conjugate gradient, with that tolerance, 
   ! with that max number of steps, and using the data contained in am_data
   write(*,*) "Fixed LIKELIHOOD", likelihood(x_sigma, am_data), dlikelihood(x_sigma, am_data)
   deallocate(am_data)

   sigma_error = x_sigma(1)  
   sigma_covariance = x_sigma(2)   ! return the values of sigma_error and sigma_covariance after the optimisation
   call print('Optimised hypers. ERROR = '//sigma_error//' COVARIANCE = '//sigma_covariance)
   
 end subroutine do_optimise_likelihood

end program force_gaussian_prediction
