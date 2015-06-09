! This Fortran95 Code/Module performs Machine Learning Force prediction based on XYZ structural Files;
! Sorting/Selecting Algorithms are used to construct sub training database during heavy machine learning task;
! Database can be updated according to an error_threshold, e.g. pred_err; 
! Internal Vector representation is used for describing the atomic structures; 
! Any issues or bugs regarding this code please contact Zhenwei Li (zhenwei.li@unibas.ch)
! Reference article : Molecular Dynamics with On-the-Fly Machine Learning of Quantum-Mechanical Forces, Zhenwei Li, James Kermode, Alessandro De Vita, PRL, 114, 096405 (2015) 
program force_gaussian_prediction

   use potential_module
   use libatoms_module
   use force_machine_learning_module

   implicit none

   type(Atoms)         :: at_in
   type(Potential)     :: pot
   type(CInOutput)     :: in
   type(Dictionary)    :: params
   real(dp)            :: r_cut, feature_len, thresh, sigma_error, error_ls, error_gp, error_gp_sq, time, error_ls_max
   real(dp)            ::  sigma_covariance, sigma_cov, dist_shift_factor, dim_tol, error_frame, force_magntd_pred, dist
   real(dp), dimension(:), allocatable           :: r_grid, m_grid, sigma, covariance_pred,force_magntd_data,force_magntd_data_select 
   real(dp), dimension(:), allocatable           :: force_proj_test, force_proj_ivs_pred, distance_confs_unsorted
   real(dp), parameter                           :: TOL_REAL=1e-7_dp, SCALE_IVS=100.0_dp
   real(dp)                                      :: force(3), variance_xyz(3), feature_inner_matrix(3,3), feature_inv(3,3), kappa, tmp
   real(dp), dimension(:,:,:), allocatable       :: ivs_set_norm, ivs_set, ivs_set_pred, nearest_data
   real(dp), dimension(:,:), allocatable         :: force_proj_ivs, force_proj_ivs_select, inv_covariance, out_u, out_vt, geometry_matrix
   real(dp), dimension(:,:), allocatable, target :: covariance, distance_matrix, distance_matrix_dij, force_covariance_matrix, covariance_tiny
   real(dp), dimension(:,:), allocatable         :: ivs_set_norm_pred, ivs_set_norm_pred_t
   real(dp), dimension(:,:), pointer             :: force_ptr, force_ptr_mm, force_ptr_tff, pred_err
   integer                                       :: i,j, k, n, n_0, t,  n_center_atom, n_loop, error 
   integer                                       :: add_vector, n_teach, func_type, weight_index, temp_integer, ii, do_atz_index
   integer	                 		 :: local_ml_optim_size, dimension_mark(3)
   integer, dimension(:), allocatable            :: distance_index, mark_zero_ivs, species
   logical                         ::  top_teach, top_pred, do_gp, do_sigma, least_sq, print_at, do_svd, print_data,do_teach, do_read, print_verbosity, update_data, do_insert, do_dynamic_x, write_data_dij, read_data_dij
   character(STRING_LENGTH)                      :: teach_file, test_file, out_file, data_dir
   
   call system_initialise(enable_timing=.true.)

   call initialise(params)
   call param_register(params, 'r_cut',  '8.0', r_cut, "the cutoff radius for the spherical atomic environment")
   call param_register(params, 'thresh', '10.0', thresh, "the threshold for doing the Singular Value Decompostion of the Covariance Matrix")
   call param_register(params, 'weight_index', '0', weight_index, "Type of weight function in deriving IVs, zero stands for the form of $exp(r/rcut)^m$")
   call param_register(params, 'sigma_error', '0.05', sigma_error, "the noise assumed on the teaching data")
   call param_register(params, 'do_atz_index', '0', do_atz_index, "0: charge-weight disabled; 1: using at%Z, 2: using ionic charge instead")
   call param_register(params, 'func_type', '0', func_type, "which kernel function is used to build the covariance matrix")
   call param_register(params, 'do_gp',  'T', do_gp, "true for doing a gaussian processes")
   call param_register(params, 'n_teach', '1000', n_teach, "the number of relevant confs you would like to do machine learning with")
   call param_register(params, 'do_sigma', 'F', do_sigma, "set true to print out the distance on every single dimension of the IVs space")
   call param_register(params, 'do_dynamic_x', 'F', do_dynamic_x, "true for updating sigma.dat every test configuration based on nearest data points")
   call param_register(params, 'top_teach', 'T', top_teach, "only the first configuration is considered when doing teaching")
   call param_register(params, 'top_pred', 'T', top_pred, "only the first configuration is considered when doing predicting")
   call param_register(params, 'least_sq',    'T', least_sq, "if true, the internal force components will be tranformed to real force using least squares")
   call param_register(params, 'do_svd', 'F', do_svd, "if true, doing inverting the feature_inner_matrix by SVD")
   call param_register(params, 'sigma_cov', '1.0', sigma_cov, "correlation length for pairs of data points") 
   call param_register(params, 'dist_shift_factor', '0.2', dist_shift_factor, "distance shifted with a given factor")
   call param_register(params, 'dim_tol', '0.005', dim_tol, "threshold for determing the dimensionality of the group of IVs")
   call param_register(params, 'teach_file', 'data.xyz', teach_file, "file to read teaching configurations from")
   call param_register(params, 'test_file', 'test.xyz', test_file, "file to read the testing configurations from")
   call param_register(params, 'print_data', 'F', print_data, "if true, print out the information for the database" )
   call param_register(params, 'data_dir', './info/', data_dir, "the directory where data files locate")
   call param_register(params, 'do_teach', 'T', do_teach, "if false, the teaching part will be skipped")
   call param_register(params, 'do_read', 'F', do_read, "if true, reading data file from the database")
   call param_register(params, 'do_insert', 'F', do_insert, "INsert any addiitional vector as force_plus")
   call param_register(params, 'print_verbosity', 'F', print_verbosity, "if true, print out verbosity stuff")
   call param_register(params, 'out_file', 'out.xyz', out_file, "output configuration file containing the predicted forces")
   call param_register(params, 'print_at', 'F', print_at, "true for print out the testing configuration with predicted forces")
   call param_register(params, 'update_data', 'F', update_data, "whether or not updating the database")
   call param_register(params, 'read_data_dij', 'F', read_data_dij, 'whether or not writing the distance file')
   call param_register(params, 'write_data_dij', 'F', write_data_dij, 'whether or not reading the pair distnace matrix') 
   call param_register(params, 'local_ml_optim_size', '0', local_ml_optim_size, "Optimise sigma_error and sigma_cov using local_ml_optim_size LS confs. If 0, no optimisation is performed")
   
   if (.not. param_read_args(params, task="fgp command line params")) then
      call print("Usage: fgp [options]")
      call system_abort("Confused by command line arguments")
   end if

   call print_title('params')
   call param_print(params)
   call finalise(params)

if (do_teach) then  
  !
  call load_iv_params(data_dir, r_grid, m_grid, species, k, add_vector, n_data=n_0) 
  ! the current k is number of defined internal vectors 
  if (add_vector > 0)   k=k+add_vector
  call print('The number of valid Feature Vectors (inc. f_mm): '//k)
  allocate(sigma(k)) ! knowing k, sigma canbe assigned
   
  ! to know the total number of teaching confs: n
  call initialise(in, teach_file, INPUT)
  n = 0  
  do i=1,  in%n_frame
        call read(in, at_in)
        if (top_teach) then
            n = n + 1
        else
            n = n +  at_in%N
        endif   
  enddo
  call finalise(in)
  call print("The total number of teaching configurations:    "//n)


  call print_title('starting the teaching process')
  call system_timer('Preparing Information from the DATABASE')
 
  ! the main work starts
  call initialise(in, teach_file, INPUT)
  call print('Got '//in%n_frame//' input frames.')

  allocate(ivs_set(k,3,n))   
  allocate(ivs_set_norm(k,3,n))
  allocate(force_proj_ivs(k, n))
  allocate(force_magntd_data(n))

 t=0 ! counter for loop the atomic configurations
 do i=1, in%n_frame
        call read(in, at_in)
        call assign_property_pointer(at_in, 'force', ptr=force_ptr)
        if (add_vector > 0) then 
           if (do_insert) then
              call assign_property_pointer(at_in, 'force_plus', ptr=force_ptr_mm) !treat the plus Force as mm
           else
             call Potential_Filename_Initialise(pot, args_str='IP SW', param_filename='SW.xml')
             call calc(pot, at_in, args_str='force=force_mm',  error=error)
             call assign_property_pointer(at_in, 'force_mm', ptr=force_ptr_mm) 
             call finalise(pot)
           endif
        endif
        ! the first vector is labelled force_mm in the input file, and is usually the SW force
        if (add_vector > 1) call assign_property_pointer(at_in, 'force_i', ptr=force_ptr_tff) 
        ! the second vector is labelled force_i
        ! call print('the frame: '//i)

        ! Including the symmetry images
        call set_cutoff(at_in, r_cut)
        call calc_connect(at_in)

     if (top_teach) then        
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
           ivs_set(j,:,t) = internal_vector(at_in, r_grid(j), m_grid(j), species(j), n_center_atom, weight_index, do_atz_index) *SCALE_IVS
           feature_len = norm(ivs_set(j,:, t))
           if (feature_len < TOL_REAL)  then
                   feature_len=1.0_dp
                   call print("WARNING: TEACHING, encountered the numerical limit in getting the unit direction of IVs") 
           endif

           ivs_set_norm(j,:,t) = ivs_set(j,:, t)/feature_len
           if (print_verbosity) call print("internal vectors ( "//r_grid(j)//" "//m_grid(j)//" "//species(j)//" )"//ivs_set(j,1,t)//" "//ivs_set(j,2,t)//" "//ivs_set(j,3,t)) 
        enddo

        if (add_vector > 0 ) then
           do j=k-add_vector+1, k
              temp_integer=k-j
              select case (temp_integer)
              case (0)
                 ivs_set(j,:,t) = force_ptr_mm(:,n_center_atom)
              case (1)
                 ivs_set(j,:,t) = force_ptr_tff(:,n_center_atom)
              end select
              
              feature_len = norm(ivs_set(j,:,t))
              if (feature_len < TOL_REAL)  then
                 feature_len=1.0_dp
                 call print("WARNING: TEACHING, encountered the numerical limit in getting the unit direction of IVs")
              endif
              
              ivs_set_norm(j,:,t) = ivs_set(j,:,t)/feature_len
              !call print("added internal vectors :"//ivs_set(j,1,t)//"  "//ivs_set(j,2,t)//" "//ivs_set(j,3,t))
           enddo
        endif

        ! store the force info for every configuration
        force_magntd_data(t) = norm(force) 
        do j=1, k
             force_proj_ivs(j,t) = dot_product(ivs_set_norm(j,:,t), force) 
        enddo 
         
     enddo   ! loop over the atoms
enddo       ! loop over the frames
  
call finalise(in)
call system_timer('Preparing Information from the DATABASE')

if (do_sigma) then
   call sigma_factor_from_single_dimension(ivs_set(:,:,:), sigma, print_verbosity, TOL_REAL, data_dir)
endif

!---------------------------------------------------------------------------------------------------------------------------------
!  	WRITTING OUT THE DATAFILES  
!---------------------------------------------------------------------------------------------------------------------------------

if (print_data) then 
  call system_timer('Writing Data File: Grid, Force and IV')
  ! for updating database, initial value taken into account
  call  write_iv_params(data_dir, r_grid, m_grid, species, k-add_vector, add_vector, n)
  if (.not. update_data ) call write_distance(0.0_dp, data_dir)  
  !
  OPEN(2,status='replace',file=trim(data_dir)//'Force.dat', form='UNFORMATTED')
  OPEN(3,status='replace',file=trim(data_dir)//'IV.dat', form='UNFORMATTED')

	do t=1, n
	    write(2)  force_proj_ivs(:, t), force_magntd_data(t)
	    do i=1,3
	      write(3)  ivs_set(:,i,t)
	    enddo
	enddo ! for n configuration to augument the existing DATABASE
        close(2)
        close(3)
        call system_timer('Writing Data File: Grid, Force and IV')

endif  ! print_data or not
end if ! do_teach, otherwise skipped


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! if do_read, teaching information to be collected from data files 
if (do_read) then

    call print_title('READING TEACHING INFORMATION')
    call system_timer('Collecting Information From Files')
    call load_iv_params(data_dir, r_grid, m_grid, species, k, add_vector, n_data=n)
    k = k + add_vector  ! adding the Classical Force Vectors

    allocate(force_proj_ivs(k,n))
    allocate(ivs_set(k,3,n))
    allocate(ivs_set_norm(k,3,n))
    allocate(force_magntd_data(n))
    allocate(sigma(k))

    OPEN(2,status='old',file=trim(data_dir)//'Force.dat',form='UNFORMATTED')
    OPEN(3,status='old',file=trim(data_dir)//'IV.dat',form='UNFORMATTED')

    do t=1, n
       read(2) force_proj_ivs(:, t), force_magntd_data(t)
       !write(*,*) "%", force_proj_ivs(:, t), force_magntd_data(t)
       do i=1, 3
          read(3) ivs_set(:,i,t)
       enddo
       !write(*,*) "$", ivs_set(:,:,t)
    enddo
    close(2)
    close(3)

  ! calculating the normalised vectors
  do t=1, n
    do j=1, k
     feature_len = norm(ivs_set(j,:, t))
     if (feature_len < TOL_REAL)  then
             feature_len=1.0_dp
             call print("WARNING: TEACHING, encountered the numerical limit in getting the unit direction of IVs")
     endif
     ivs_set_norm(j,:,t)=ivs_set(j,:, t)/feature_len
    enddo
 enddo
  
    ! reading sigma file
    open(unit=1, status='old', file=trim(data_dir)//'sigma.dat', form='formatted')
    read(1, *) sigma
    close(unit=1)
    call system_timer('Collecting Information From Files')
    call print('sigma is:    '//sigma)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ((do_teach) .and. (.not. do_read)) then 
    open(unit=1, status='old', file=trim(data_dir)//'sigma.dat', form='formatted')
    read(1, *) sigma
    close(unit=1)
    call system_timer('Collecting Information From Files')
    call print('sigma is:    '//sigma)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! the following part can be changed to a subroutine
call print_title('starting the predicting process')
!
if (n_teach > n) then
  n_teach = n
  call print("the Actual Number of Configs for Prediction: "//n_teach)
  call print("WARNNING: Trying to select more configurations than the database has")
endif
!
allocate(covariance_pred(n_teach))
allocate(force_proj_ivs_pred(k))
allocate(mark_zero_ivs(k))
allocate(geometry_matrix(3,k))
allocate(covariance(n_teach,n_teach)) 
allocate(distance_confs_unsorted(n), distance_index(n))
if (local_ml_optim_size > 0) then 
allocate(distance_matrix(local_ml_optim_size,local_ml_optim_size), force_covariance_matrix(local_ml_optim_size,local_ml_optim_size), covariance_tiny(local_ml_optim_size,local_ml_optim_size))
endif
allocate(inv_covariance(n_teach,n_teach))
allocate(force_proj_ivs_select(k, n_teach))
allocate(force_magntd_data_select(n_teach))
allocate(ivs_set_norm_pred(k,3))
if (read_data_dij) allocate(distance_matrix_dij(n, n))

if (read_data_dij) then
          call system_timer('Read Pair Distance and Bulid Covariance')
          call read_distance_file(distance_matrix_dij(:,:), data_dir)
          !write(*,*)  'Distance Read', distance_matrix_dij
          call system_timer('Read Pair Distance and Bulid Covariance')
endif

call initialise(in, test_file, INPUT)
do i=1, in%n_frame 
   call read(in, at_in)
   if (.not. get_value(at_in%params, 'time', time))  then
     call print("TIME of FRAME : UNKOWN")
   else
     call print("TIME of FRAME :"//time)
   endif
   call add_property(at_in, 'force', 0.0_dp, n_cols=3, overwrite=.false.) ! how do you know dont want overwrite
   call assign_property_pointer(at_in, 'force', ptr=force_ptr)

   ! for writing out the predicted error
   call add_property(at_in, 'pred_err', 0.0_dp, n_cols=3, overwrite=.true.)
   call assign_property_pointer(at_in, 'pred_err', ptr=pred_err)

   if (add_vector > 0) then 
           if (do_insert) then
              call assign_property_pointer(at_in, 'force_plus', ptr=force_ptr_mm)
           else
             call Potential_Filename_Initialise(pot, args_str='IP SW', param_filename='SW.xml')
             call calc(pot, at_in, args_str='force=force_mm',  error=error)
             call assign_property_pointer(at_in, 'force_mm', ptr=force_ptr_mm)
             call finalise(pot)
           endif
   endif
   if (add_vector>1) call assign_property_pointer(at_in, 'force_i', ptr=force_ptr_tff)

   call set_cutoff(at_in, r_cut)
   call calc_connect(at_in)
   if (top_pred) then
      n_loop = 1
   else
      n_loop = at_in%N
   endif

   allocate(ivs_set_pred(k,3,n_loop))

   error_frame=0.0_dp           ! for calculating the average error of the system.

  do n_center_atom=1, n_loop    ! loop over atoms in the Frame, make prediction (incl. error)

     do j=1, k                  ! reset the mark for zero internal vectors
         mark_zero_ivs(j)=1
     enddo

     do j= 1, k-add_vector
        ivs_set_pred(j,:,n_center_atom) = internal_vector(at_in, r_grid(j), m_grid(j), species(j), n_center_atom, weight_index, do_atz_index) *SCALE_IVS
        if (print_verbosity) call print("internal vectors ( "//r_grid(j)//" "//m_grid(j)//" "//species(j)//" ):  "//ivs_set_pred(j,1,n_center_atom)//"  "//ivs_set_pred(j,2,n_center_atom)//"  "//ivs_set_pred(j,3,n_center_atom))
       feature_len = norm(ivs_set_pred(j,:,n_center_atom))

        if (feature_len < TOL_REAL)  then
            feature_len=1.0_dp
            mark_zero_ivs(j)=0
            call print("WARNING: PREDICTION, encountered the numerical limit in getting the unit direction of IVs")
        endif

        ivs_set_norm_pred(j,:) = ivs_set_pred(j,:,n_center_atom)/feature_len
     enddo

     if (add_vector>0) then
         do j = k-add_vector+1, k
            temp_integer=k-j
            select case (temp_integer)
            case (0)
               ivs_set_pred(j,:,n_center_atom)=force_ptr_mm(:, n_center_atom)
            case (1)
               ivs_set_pred(j,:,n_center_atom)=force_ptr_tff(:, n_center_atom)
            end select
            feature_len = norm(ivs_set_pred(j,:,n_center_atom))

            if (feature_len < TOL_REAL) then
               feature_len=1.0_dp
               mark_zero_ivs(j)=0
               call print("WARNING: PREDICTION, encountered the numerical limit in getting the unit direction of IVs")
            endif
            
            ivs_set_norm_pred(j,:) = ivs_set_pred(j,:,n_center_atom)/feature_len
            if (print_verbosity) call print("added internal vectors : "//ivs_set_pred(j,1,n_center_atom)//"  "//ivs_set_pred(j,2,n_center_atom)//"  "//ivs_set_pred(j,3,n_center_atom))
         enddo
      endif


     if (print_verbosity) then
       do j=1, k 
          do ii=1, k   
           if(print_verbosity) then
              call print("ATOM : "//n_center_atom//" ivs_set: ("//j//" "//ii//") "//dot_product(ivs_set_pred(j,:,n_center_atom), ivs_set_norm_pred(ii,:)) )
           endif
          enddo   
       enddo
     endif

      
      ! initialise distance_index(:)
      do t=1, n
         distance_index(t) =t 
      enddo

        ! do the sorting and selection
      call sorting_configuration(ivs_set_pred(:,:,n_center_atom), ivs_set, ivs_set_norm_pred, ivs_set_norm, sigma, distance_confs_unsorted, distance_index)
      call print("The First 20 Closest data Configurations are :"//distance_index(:20))
      call print("Min and Max DISTANCE with Index after Sorting: "//distance_confs_unsorted(distance_index(1))//" and "// &
	                                    distance_confs_unsorted(distance_index(n_teach))//"  the INDEX: "//distance_index(1)//" and "//distance_index(n_teach))
       !!!!!!!!!!!!!
    if (do_dynamic_x) then 
        allocate(nearest_data(k,3,n_teach))
        do ii=1, n_teach
           nearest_data(:,:,ii)=ivs_set(:,:,distance_index(ii))
        enddo             
        call sigma_factor_from_single_dimension(nearest_data, sigma, print_verbosity, TOL_REAL, data_dir)
        ! update 'sigma.dat'
        deallocate(nearest_data) 
     endif     
   
     if ( distance_confs_unsorted( distance_index(1) ) > dist_shift_factor) then  
                ! distance_confs_unsorted: the pair distance betwen test and database.
                ! sigma_cov: the original input for covariance length 
                sigma_covariance = sigma_cov*distance_confs_unsorted(distance_index(1))/dist_shift_factor    
                ! this intrinsically changes the value of SIGMA_COVARIANCE and enters next loop
                call print("sigma_covariance are shifted by a factor of :"//distance_confs_unsorted(distance_index(1))/dist_shift_factor)
                call print("Modified sigma_covariance :"//sigma_covariance)
     else
               sigma_covariance = sigma_cov
               ! the sigma_covariance used in GP
     endif

     if (.false.) then
         ! massive output, only for testing use
         call print("Detail on the teaching information")
         do ii = 1, n_teach
            call print("Index: "//distance_index(ii)//" Distance: "//distance_confs_unsorted(ii))
            call print("Force in IVs space: "//force_proj_ivs(:,distance_index(ii)))
         enddo
     endif
   
     ! max marginal likelihood 
      if (local_ml_optim_size > 0) call print("Optimise: "//local_ml_optim_size)
      !!! ML optimisation here: first calculates the necessary input matrices, then optimises
      if (local_ml_optim_size > 0) then
         call print("Local Maximum Likelihood selected")
         do t = 1, local_ml_optim_size  ! loop to generate the force_cov_mat and the distance_matrix of the restricted subset of the LS
            do j = t, local_ml_optim_size
               force_covariance_matrix(t,j) = dot_product(force_proj_ivs(:,distance_index(t)), force_proj_ivs(:,distance_index(j))) 
               ! check if possible, may be a wreck because IVs make absolute representation of the force impossible 
               force_covariance_matrix(j,t) = force_covariance_matrix(t,j)
               covariance_tiny(t,j) = cov(ivs_set(:,:,distance_index(t)), ivs_set(:,:,distance_index(j)), &
                    ivs_set_norm(:,:,distance_index(t)), ivs_set_norm(:,:,distance_index(j)), sigma, sigma_covariance, func_type=func_type, &
                    distance=distance_matrix(t,j))
               covariance_tiny(j,t) = covariance_tiny(t,j)
               distance_matrix(j,t) = distance_matrix(t,j)
            enddo
         enddo
         call do_optimise_likelihood(sigma_error, sigma_covariance) 
      endif

      if (read_data_dij) then
          do t = 1, n_teach  ! loop to generate the covariance matrix of the learning set
                 dist=distance_matrix_dij(distance_index(t), distance_index(t))
                 covariance(t,t) = cov_dij(sigma_covariance, dist, func_type)
           !      write(*,*) 'diagoanl value', covariance(t,t)
                 if (do_gp) covariance(t,t) = covariance(t,t) + sigma_error**2
                 do j = t+1, n_teach ! above diagonal
                    dist=distance_matrix_dij(distance_index(t),distance_index(j))
                    covariance(t,j) = cov_dij(sigma_covariance, dist, func_type)
                    covariance(j,t) = covariance(t,j)
                 enddo
          enddo
      else
	      call system_timer('Evaluating the Cov')
	      !Generate covariance matrix  
	      !$omp parallel default(none) shared(covariance, n_teach, ivs_set, ivs_set_norm, distance_index, sigma, sigma_error, sigma_covariance, do_gp, func_type), private(j)
	      !$omp do 
	      do t = 1, n_teach  ! loop to generate the covariance matrix of the learning set
		 covariance(t,t) = cov(ivs_set(:,:,distance_index(t)), ivs_set(:,:,distance_index(t)), &
			ivs_set_norm(:,:,distance_index(t)), ivs_set_norm(:,:,distance_index(t)), sigma, sigma_covariance, func_type=func_type)
		 if (do_gp) covariance(t,t) = covariance(t,t) + sigma_error**2
		 do j = t+1, n_teach ! above diagonal
		    covariance(t,j) = cov(ivs_set(:,:,distance_index(t)), ivs_set(:,:,distance_index(j)), &
			 ivs_set_norm(:,:,distance_index(t)), ivs_set_norm(:,:,distance_index(j)), sigma, sigma_covariance, func_type=func_type)
		    covariance(j,t) = covariance(t,j)  
		enddo
	       enddo
	      !$omp end parallel
	       call system_timer('Evaluating the Cov')
      endif ! not read_data_dij

     call system_timer('Inverting the Covariance Matrix')
      if (do_gp) then
         call inverse(covariance, inv_covariance)
      else
         ! To Do Sigular Value Decomposition (SVD): A = U*SIGMA*VT
         call inverse_svd_threshold(covariance, thresh, result_inv=inv_covariance)
      endif
      call system_timer('Inverting the Covariance Matrix')
      

      do t=1, n_teach
         force_proj_ivs_select(:,t)=force_proj_ivs(:,distance_index(t))
         force_magntd_data_select(t) = force_magntd_data(distance_index(t))  
      enddo
      
      ! the test configuration covariance vector
      do t= 1, n_teach
        covariance_pred(t) = cov(ivs_set_pred(:,:,n_center_atom), ivs_set(:,:,distance_index(t)), ivs_set_norm_pred, &
                                         ivs_set_norm(:,:,distance_index(t)), sigma, sigma_covariance, func_type=func_type)
      enddo
      
      ! the predicted force components on each of the internal directions
      force_proj_ivs_pred(:) = matmul(covariance_pred, matmul(inv_covariance, transpose(force_proj_ivs_select(:,:)) )) 
        
      ! Predict the Force Magnitude
      force_magntd_pred =  dot_product(covariance_pred, matmul(inv_covariance, force_magntd_data_select))

      ! we need to mannually set the force components to be zero, projected on the zero internal-vector directions
      do j=1, k
         force_proj_ivs_pred(j)=force_proj_ivs_pred(j)*real(mark_zero_ivs(j), dp) 
      enddo

     if (print_verbosity) then 
         do j=1, k
               call print("ivs_set_norm_pred"//ivs_set_norm_pred(j,:))
         enddo
     endif

     allocate(out_u(3,3), out_vt(3,3))      ! to obtain the transformation matrix with the principal axis
     call inverse_svd_threshold(transpose(ivs_set_norm_pred) .mult. ivs_set_norm_pred, thresh, u_out=out_u, vt_out=out_vt)

     allocate(ivs_set_norm_pred_t(k,3))

     call internal_dimension_mark(ivs_set_pred(:,:,n_center_atom) .mult. out_u, dim_tol, dimension_mark)

     ivs_set_norm_pred_t = ivs_set_norm_pred .mult. out_u
 
  ! ivs_set_pred_t contains the Internal-directions after PCA transformation.
  if (print_verbosity) then
       do j=1, k
           write(*,*) "ivs_set_norm_pred_t : ", ivs_set_norm_pred_t(j, :)
       enddo
  endif

  select case (sum(dimension_mark(:)))
    case(0)
         force= 0.0_dp
    case(1)
         force=(transpose(ivs_set_norm_pred) .mult. force_proj_ivs_pred )/real(sum(mark_zero_ivs),dp)
    case(2) 
         feature_inner_matrix=transpose(ivs_set_norm_pred_t) .mult. ivs_set_norm_pred_t
         call rank_lowered_inverse(feature_inner_matrix, feature_inv) 
         force = out_u .mult. feature_inv .mult. transpose(ivs_set_norm_pred_t) .mult. force_proj_ivs_pred 
    case(3)

         feature_inner_matrix=transpose(ivs_set_norm_pred) .mult. ivs_set_norm_pred
         if (do_svd) then   
             if (print_verbosity) write(*,*)  "feature inner matrix :", feature_inner_matrix
             call inverse_svd_threshold(feature_inner_matrix, thresh, result_inv=feature_inv)
             if (print_verbosity) write(*,*)  "inverse feature inner matrix :", feature_inv
         else
             if (print_verbosity) write(*,*)  "feature_inner_matrix :", feature_inner_matrix
             call inverse(feature_inner_matrix, feature_inv)
             if (print_verbosity) write(*,*)  "inverse feature inner matrix :", feature_inv
         endif

         if (least_sq) then
            force = feature_inv .mult. transpose(ivs_set_norm_pred) .mult. force_proj_ivs_pred
         else
            call real_expection_sampling_components(ivs_set_norm_pred, force_proj_ivs_pred, force)
         endif
   end select

  deallocate(out_u)
  deallocate(out_vt)
  deallocate(ivs_set_norm_pred_t)

   call print("Predicted Force Magnitude: "//force_magntd_pred)
   call print("force in external space : "//force//" Magnitude: "//norm(force))
   call print("the original force:"//force_ptr(:, n_center_atom)//" Magnitude: "//norm(force_ptr(:, n_center_atom)))
   call print("max error :    "//maxval(abs(force_ptr(:,n_center_atom)-force))//" norm  error :  "//norm(force_ptr(:,n_center_atom)-force))
   if (norm(force_ptr(:,n_center_atom)) > TOL_REAL) then 
           call print("Relative Error : "//norm(force_ptr(:,n_center_atom)-force)/norm(force_ptr(:,n_center_atom)) )
   else
           call print("Relative Error : "//" ZERO Magnititude ")
   endif
   ! measure the consistence between the least-square inference and the components obtained from GP processes
   if (least_sq) then
       error_ls=0.0_dp
       do j=1, k
            error_ls = error_ls + (dot_product(ivs_set_norm_pred(j,:), force(:)) - force_proj_ivs_pred(j))**2
       enddo
       error_ls = sqrt(error_ls / real(k,dp))
       call print("deviation of least-square inference : "//error_ls)
   error_ls_max=maxval(abs( matmul(ivs_set_norm_pred, force) - force_proj_ivs_pred(:) ))
   call print("Max deviation Component of LS : "//error_ls_max)
   endif 

   kappa = cov(ivs_set_pred(:,:,n_center_atom), ivs_set_pred(:,:,n_center_atom), ivs_set_norm_pred, ivs_set_norm_pred, sigma, sigma_covariance, func_type=func_type) + sigma_error**2
   ! the passing of  uncertainty from  Gaussian Processes to the Force Vector
   error_gp_sq = kappa - covariance_pred .dot. matmul(inv_covariance, covariance_pred)
   error_gp = sqrt(abs(error_gp_sq))
   geometry_matrix = feature_inv .mult. transpose(ivs_set_norm_pred) 

   do j=1, 3  ! Deviation in External Cartesian Coordinate 
      variance_xyz(j)= sqrt(dot_product(geometry_matrix(j,:), geometry_matrix(j,:))) * (error_gp) 
      !call print("Predicted Error Factor :"//sum(abs(geometry_matrix(j,:)))//" Dimension :"//j)
   enddo

   call print("GP Uncertainty in Cartesian Coordinate: "//variance_xyz//" GP variance: "//error_gp//" variance Squre :"//error_gp_sq)   
   pred_err(1, n_center_atom)=norm(variance_xyz) 
   pred_err(2, n_center_atom)=abs(force_magntd_pred - norm(force))
   pred_err(3, n_center_atom)=error_ls
  if (update_data) then 
          call system_timer('Updating the DATABSE : dist.dat force.dat ivs.dat')
          allocate(force_proj_test(k)) ! the projection of real force on test conf
          do j = 1, k 
               force_proj_test(j) = dot_product(force_ptr(:,n_center_atom), ivs_set_norm_pred(j,:) )
          enddo
          call  write_iv_params(data_dir, r_grid, m_grid, species, k-add_vector, add_vector, n+1)
          if (write_data_dij) call write_distance(distance_confs_unsorted, data_dir) 
          ! append distance file
          if (print_verbosity) then  
                  write(*,*) 'Distance', distance_confs_unsorted(:)
          endif
          OPEN(2,status='old',file=trim(data_dir)//'Force.dat', form='UNFORMATTED', position='APPEND')
          OPEN(3,status='old',file=trim(data_dir)//'IV.dat', form='UNFORMATTED',  position='APPEND')
          write(2)  force_proj_test(:), norm(force_ptr(:,n_center_atom))
          do j=1, 3
                write(3)  ivs_set_pred(:,j,n_center_atom)
          enddo
	  close(2)
	  close(3)
          deallocate(force_proj_test)
          call system_timer('Updating the DATABSE : dist.dat force.dat ivs.dat')
   endif
   if (print_at) force_ptr(:,n_center_atom)=force
 enddo  ! loop over positions

  deallocate(ivs_set_pred)

   error_frame = maxval( pred_err)
   call print("Error Bar of the Frame :"//error_frame)

   if (print_at) then 
       ! write the output configurations with predicted force, we dont want appending practically
       call write(at_in, out_file, append=.False.) 
   endif

enddo                         ! loop over frames

call finalise(in)

deallocate(ivs_set_norm_pred)   ! because it is defined within the loop

deallocate(geometry_matrix)
deallocate(r_grid, m_grid, species, sigma)
deallocate(force_proj_ivs_pred)
deallocate(ivs_set_norm)
deallocate(force_magntd_data, force_magntd_data_select)
deallocate(ivs_set)
deallocate(mark_zero_ivs)
deallocate(force_proj_ivs)
deallocate(force_proj_ivs_select)
deallocate(covariance, inv_covariance)
deallocate(covariance_pred)
deallocate(distance_index, distance_confs_unsorted)
if(allocated(force_covariance_matrix)) deallocate(force_covariance_matrix)
if(allocated(covariance_tiny)) deallocate(covariance_tiny)
if(allocated(distance_matrix)) deallocate(distance_matrix)
if (allocated(distance_matrix_dij)) deallocate(distance_matrix_dij)
call system_finalise()
call print('-------------END OF FGP----------------')

end program force_gaussian_prediction
