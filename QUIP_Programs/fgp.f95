program force_gaussian_prediction
   use libatoms_module
   implicit none

   type(Atoms)                                  :: at_in, at
   type(CInOutput)                              :: in
   type(Dictionary)                             :: params
   real(dp)                                     :: r_cut, r_min, m_min, m_max, feature_len, theta, thresh, sigma_error, cutoff_len_ivs
   real(dp), parameter                          :: TOL_REAL=1e-7_dp, SCALE_IVS=100.0_dp
   integer                                      :: i,j, k, n, t, at_n, n_loop, preci, r_mesh, m_mesh, add_vector, n_relavant_confs
   integer, dimension(:), allocatable           :: distance_index
   real(dp)                                     :: force(3), feature_inv(3,3)
   real(dp), dimension(:), allocatable          :: r_grid, m_grid, sigma, theta_array, covariance_pred, force_proj_ivs_pred, force_proj_tmp, distance_confs
   real(dp), dimension(:,:,:), allocatable      :: feature_matr_normalised, feature_matr
   real(dp), dimension(:,:), allocatable        :: force_proj_ivs, covariance, inv_covariance, distance_ivs
   real(dp), dimension(:,:), allocatable        :: feature_matr_normalised_pred, feature_matr_pred
   real(dp), dimension(:,:), pointer            :: force_ptr, force_ptr_mm
   logical                                      :: spherical_cluster_teach, spherical_cluster_pred, do_gp, fix_sigma, do_conf_filter
   character(STRING_LENGTH)                     :: teaching_file, grid_file, test_file

   call system_initialise(enable_timing=.true.)

   call initialise(params)
   call param_register(params, 'theta',  '2.5', theta, "the correlation length of the data")
   call param_register(params, 'r_min',  '0.5', r_min, "the minimum radius for the spherical atomic environment")
   call param_register(params, 'r_cut',  '8.0', r_cut, "the cutoff radius for the spherical atomic environment")
   call param_register(params, 'm_min',  '1.0', m_min, "the minimum m for calculating the atomic environment")
   call param_register(params, 'm_max',  '5.0', m_max, "the maxmium m for calculating the atomic environment")
   call param_register(params, 'n_relavant_confs', '200', n_relavant_confs, "the number of relavant confs you'd like to teach")  
   call param_register(params, 'cutoff_len_ivs', '0.2', cutoff_len_ivs, "the cutoff lenght for IVs to be considered valid when generating the grid")
   call param_register(params, 'thresh', '1.0', thresh, "the threshold for doing the Sigular Value Decompostion of the Covariance Matrix")
   call param_register(params, 'preci',  '6',   preci,  "the screening accuracy on the edge atoms")
   call param_register(params, 'add_vector', '0', add_vector, "the number of vectors you like to add into the internal_vector representation")
   call param_register(params, 'sigma_error', '0.01', sigma_error, "the noise assumed on the teaching data")
   call param_register(params, 'r_mesh', '6',   r_mesh, "grid finess of r0")
   call param_register(params, 'm_mesh', '6',   m_mesh, "grid finess of m")
   call param_register(params, 'do_gp',  'F', do_gp, "true for doing a gaussian processes, instead of SVD")
   call param_register(params, 'do_conf_filter', 'F', do_conf_filter, "true for doing configuratio filtering")
   call param_register(params, 'spherical_cluster_teach', 'T', spherical_cluster_teach, "only the first atom in the cluster are considered when doing teaching")
   call param_register(params, 'spherical_cluster_pred', 'T', spherical_cluster_pred, "only the first atom in the cluster are considered when doing predicting")
   call param_register(params, 'fix_sigma',  'F', fix_sigma, "true, if you want manually input sigma")
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

  allocate(feature_matr(k,3,n))   
  allocate(feature_matr_normalised(k,3, n))
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
          write(*,*) "atomic force in real space:", norm(force)
          call print("The atom:"//at_n)

        do j=1, k-add_vector
           feature_matr(j, :, t) = internal_vector(at_in, r_grid(j), m_grid(j), at_n)*SCALE_IVS
           feature_len = norm(feature_matr(j,:, t))

           if (feature_len < TOL_REAL)   then
                   feature_len=1.0_dp
                   call print("WARNNING: TEACHING, encountered the decimal limit in getting the unit direction of IVs") 
           endif

           feature_matr_normalised(j,:,t) = feature_matr(j,:, t)/feature_len
           write(*,*) "internal vectors (", r_grid(j), m_grid(j), ")", feature_matr(j,1,t), feature_matr(j,2,t), feature_matr(j,3,t) 
        enddo

        if (add_vector >0 ) then
           do j=k-add_vector+1, k
                feature_matr(j,:,t) = force_ptr_mm(:,at_n)
                feature_len = norm(feature_matr(j,:,t))

                if (feature_len < TOL_REAL)  then
                    feature_len=1.0_dp
                    call print("WARNNING: TEACHING, encountered the decimal limit in getting the unit direction of IVs")
                endif

                feature_matr_normalised(j,:,t) = feature_matr(j,:,t) / feature_len
                write(*,*) "internal vectors (", "                )", feature_matr(j,1,t), feature_matr(j,2,t), feature_matr(j,3,t)
           enddo
        endif

        ! store the force info. for every configuration
        do j=1, k
            force_proj_ivs(j,t) = dot_product(feature_matr_normalised(j,:,t), force) 
            call print("Force projection on IV"//j//" is:  "//force_proj_ivs(j,t))
        enddo   
     enddo   ! loop over the atoms
 enddo       ! loop over the frames

 call finalise(in)

 call print_title('starting the teaching process')

 allocate(distance_ivs(n,n))
 allocate(sigma(k))
 allocate(theta_array(k))

if (fix_sigma) then
     open(unit=1, file='sigma.dat')
     read(1, *) sigma
     close(unit=1)
else 
  theta_array=theta
 do t = 1, k 
    do i = 1, n
        do j=1, n
              distance_ivs(i,j) = distance_bent_space(feature_matr(t,:,i), feature_matr(t,:,j), feature_matr_normalised(:,:,i), feature_matr_normalised(:,:,j)) 
        enddo
    enddo
    sigma(t) = maxval(distance_ivs(:,:))/theta_array(t) 
    if  (abs(sigma(t))<TOL_REAL) then
              sigma(t)=TOL_REAL
              call print("WARNNING: SIGMA, encountered the decimal limit in the getting Sigma from Theta")
    endif
 enddo
endif
 call print('sigma is:    '//sigma)

  ! to establish the Covariance Matrix
  allocate(covariance(n,n))
  do i = 1, n
      do j=1, n
            covariance(i,j) = cov(feature_matr(:,:,i), feature_matr(:,:,j), feature_matr_normalised(:,:,i), feature_matr_normalised(:,:,j), sigma, k) 
      enddo
  enddo

  ! if doing Gaussian Processes, adding an noise "sigma_error" for the teaching data
  if (do_gp) then 
  do i=1, n
      covariance(i,i) = covariance(i,i) + sigma_error 
  enddo
  endif
  
  do i=1, n
     write(*,*) "Covariance:", covariance(i,2), force_proj_ivs(3,i)-force_proj_ivs(3,2)
  enddo

 allocate(inv_covariance(n,n))

 if (do_gp) then
     call inverse(covariance, inv_covariance)
 else
     ! To Do Sigular Value Decomposition (SVD): A = U*SIGMA*VT
     inv_covariance = inverse_svd_threshold(covariance, n, thresh)
 endif
 
 write(*,*) "MAX and MIN components of inverse_covariance:", maxval(inv_covariance(2,:)), minval(inv_covariance(2,:))

 deallocate(covariance)
 deallocate(distance_ivs)

 call print_title('starting the predicting process')

 allocate(covariance_pred(n))
 allocate(feature_matr_pred(k,3))
 allocate(feature_matr_normalised_pred(k,3))
 allocate(force_proj_ivs_pred(k))
 allocate(force_proj_tmp(k))

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
            feature_matr_pred(j,:) = internal_vector(at_in, r_grid(j), m_grid(j), at_n)*SCALE_IVS
            write(*,*) "internal vectors (",r_grid(j), m_grid(j),"):   ",feature_matr_pred(j,1), feature_matr_pred(j,2), feature_matr_pred(j,3)
            feature_len = norm(feature_matr_pred(j,:))

            if (feature_len < TOL_REAL)  then
                   feature_len=1.0_dp
                   call print("WARNNING: PREDICTION, encountered the decimal limit in getting the unit direction of IVs")
            endif  

            feature_matr_normalised_pred(j,:) = feature_matr_pred(j,:)/feature_len
      enddo

      if (add_vector > 0) then
           do j = k-add_vector+1, k
              feature_matr_pred(j,:) = force_ptr_mm(:, at_n)
              feature_len = norm(feature_matr_pred(j,:))

              if (feature_len < TOL_REAL) then
                   feature_len=1.0_dp
                   call print("WARNNING: PREDICTION, encountered the decimal limit in getting the unit direction of IVs")
              endif

              feature_matr_normalised_pred(j,:) = feature_matr_pred(j,:)/feature_len
              write(*,*) "internal vectors (","                         ):",feature_matr_pred(j,1), feature_matr_pred(j,2), feature_matr_pred(j,3)
           enddo
      endif

      allocate(distance_confs(n))
      allocate(distance_index(n))

      do t= 1,n
         covariance_pred(t) = cov(feature_matr_pred, feature_matr(:,:,t), feature_matr_normalised_pred(:,:), feature_matr_normalised(:,:,t), sigma, k, distance = distance_confs(t))
      enddo
     
  
      call insertion_sort(distance_confs, idx=distance_index)   
      write(*,*) "DISTANCE:", distance_confs(1), distance_confs(n), distance_index(1), distance_index(n)
      deallocate(distance_confs)

      if (do_conf_filter) then
          open(unit=100, file='data_embed.xyz', status='replace')
          do t=1, n_relavant_confs
              call read(at, teaching_file, frame=distance_index(t)-1)
              call write(at, 'data_embed.xyz', append=.true.)
          enddo
          call finalise(at)
          close(100)
       endif 
 
      force_proj_ivs_pred(:) = matmul(covariance_pred, matmul(inv_covariance, transpose(force_proj_ivs(:,:)) )) 

      do j=1, k
            force_proj_tmp(j) = dot_product(feature_matr_normalised_pred(j,:), force_ptr(:,at_n))
      enddo   

   
      do j=1, k 
            call print("Force in IV space"//j//": "//force_proj_ivs_pred(j)//": "//force_proj_tmp(j)//": "//abs(force_proj_ivs_pred(j)-force_proj_tmp(j)))
            !  predicted force: real force: absolute difference 
      enddo

       ! using least-squares to restore the force in the External Cartesian Space  
       call inverse(matmul(transpose(feature_matr_normalised_pred), feature_matr_normalised_pred), feature_inv)  
       force = feature_inv .mult. transpose(feature_matr_normalised_pred) .mult. force_proj_ivs_pred
       call print("force in external space:"//force)
       call print("the original force:"//force_ptr(:, at_n))
       call print("max error :    "//maxval(abs(force_ptr(:,at_n)-force)))

       if (do_gp) then
           call print("predicted error :"//( 1.0_dp + sigma_error - covariance_pred .dot. matmul(inv_covariance, covariance_pred)))
       endif  
 
    enddo  ! loop over postions
  enddo    ! loop over frames

 call finalise(in)

 deallocate(m_grid)
 deallocate(r_grid)
 deallocate(force_proj_ivs_pred)
 deallocate(feature_matr_normalised)
 deallocate(feature_matr_normalised_pred)
 deallocate(feature_matr)
 deallocate(feature_matr_pred)
 deallocate(sigma)
 deallocate(theta_array)
 deallocate(distance_index)
 deallocate(force_proj_ivs)
 deallocate(force_proj_tmp)
 deallocate(inv_covariance) 
 deallocate(covariance_pred)
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


 function cov(feature_matr1, feature_matr2, bent_space1, bent_space2, sigma, k, distance)

    real(dp), intent(in)                  ::        feature_matr1(:,:), feature_matr2(:,:), bent_space1(:,:), bent_space2(:,:), sigma(:)
    real(dp), intent(out), optional        ::        distance
    real(dp)                              ::        cov, d_sq
    integer                               ::        i
    integer, intent(in)                   ::        k
 
    d_sq = 0.0_dp    
    do i=1, k
          d_sq = d_sq + (distance_bent_space(feature_matr1(i,:), feature_matr2(i,:), bent_space1, bent_space2) **2) /(sigma(i)**2) /2.0_dp
    enddo

    if (present(distance))  distance=sqrt(d_sq)
    cov = exp(-1.0_dp * d_sq)

 endfunction cov


 function inverse_svd_threshold(in_matr, n, thresh)
   real(dp), intent(in)                  ::         in_matr(:,:), thresh
   integer,  intent(in)                  ::         n
   real(dp)                              ::         w(n), sigmainv(n,n), u(n,n), vt(n,n), inverse_svd_threshold(n,n)
   real(dp), allocatable                 ::         work(:)
   real(dp), parameter                   ::         TOL_SVD = 1e-13_dp
   integer                               ::         info, i, lwork, j
 
   call print('entering inverse_svd_threshold')

   if (n <= 3) then
      call print('in_matr')
      call print(in_matr)
   end if

   lwork = -1
   allocate(work(1))
   call DGESVD('A','A',n,n,in_matr,n,w,u,n,vt,n, work, lwork, info)

   lwork = work(1) 
   deallocate(work)
   allocate(work(lwork))

   call DGESVD('A','A',n,n,in_matr,n,w,u,n,vt,n, work, lwork, info)
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
