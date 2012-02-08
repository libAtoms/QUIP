program test_cova

   use libatoms_module
   use real_space_covariance_module

   implicit none


   type(RealSpaceCovariance), target :: rscov
   type(CInOutput) :: in
   type(Atoms) :: at
   real(dp), pointer :: force_ptr(:,:)
   real(dp) :: force(3), error, r0
   type(Dictionary) :: params
   character(STRING_LENGTH) :: teaching_file, test_file
   logical :: do_optimise, use_qr, test_data, do_teach, print_covariance, test_data_2
   integer n_grid_minim
   real(dp) :: sigma_alpha, sigma_beta, sigma_covariance, sigma_error, cutoff_in, transition_width, covariance_tol

   type(minimise_overlap) :: am
   integer :: i, j, am_data_size
   character(len=1), dimension(1) :: am_mold
   character(len=1), allocatable :: am_data(:)
   real(dp) :: q(4), cov, d2, force_covariance, force_difference

   call system_initialise(enable_timing=.false.)

   call initialise(params)
   call param_register(params, 'n_grid_minim', '5', n_grid_minim, "start minimisations from the n_grid_minim biggest overlap grid points")
   call param_register(params, 'sigma_alpha', '0.1', r0, "Set sigma_rho via equation sigma_rho = sigma_alpha*r**2 + sigma_beta")
   call param_register(params, 'sigma_beta', '0.0', r0, "Set sigma_rho via equation sigma_rho = sigma_alpha*r**2 + sigma_beta")
   call param_register(params, 'sigma_covariance', '0.5', sigma_covariance, "sigma used for construction of covariance")
   call param_register(params, 'sigma_error', '0.1', sigma_error, "expected error on forces")
   call param_register(params, 'cutoff', '5.0', cutoff_in, "cutoff radius used for rho")
   call param_register(params, 'transition_width', '0.5', transition_width, "transition width used for rho")
   call param_register(params, 'covariance_tol', '0.001', covariance_tol, "Covariance tolerance used for minimizations")
   call param_register(params, 'teaching_file', 'data.xyz', teaching_file, "file to read teaching configurations from")
   call param_register(params, 'test_file', 'test.xyz', test_file, "file to read test configuration from")
   call param_register(params, 'do_teach', 'T', do_teach, "do teaching? (default T)")
   call param_register(params, 'do_optimise', 'F', do_optimise, "do likelihood optimisation in the space of hyperparameters")
   call param_register(params, 'print_covariance', 'F', print_covariance, "print covariance of each config with frame 0")
   call param_register(params, 'use_qr', 'F', use_qr, "Use QR decomposition instead of Cholesky for likelihood optimisation")
   call param_register(params, 'test_data', 'T', test_data, "Include teaching data in test set (default T)")
   call param_register(params, 'test_data_2', 'F', test_data_2, "Test data using method #2 (default F)")

   if (.not. param_read_args(params, task="test_cova command line params")) then
      call print("Usage: test_cova [options]")
      call system_abort("Confused by command line arguments")
   end if

   call print_title('params')
   call param_print(params)
   call finalise(params)

   call print_title('teaching')
   call initialise(in, teaching_file, INPUT)
   call print('Got '//in%n_frame//' teaching configurations.')

   call initialise(rscov, in%n_frame, n_grid_minim, sigma_alpha, sigma_beta, sigma_covariance, sigma_error, cutoff_in, transition_width, covariance_tol)
   do i=1,in%n_frame
      call read(in, rscov%data(i))
      call centre_on_atom_1(rscov%data(i))
   end do
   call finalise(in)

   if (print_covariance) then
      mainlog%prefix = 'COVARIANCE'
      do i=1,rscov%n
         call verbosity_push(PRINT_SILENT)
         cov = covariance(rscov, rscov%data(1), rscov%data(i), &
              q=rscov%q_min(:,i), &
              distance=d2, force_covariance=force_covariance, &
              force_difference=force_difference)
         call verbosity_pop()
         call print(i//' '//cov//'  '//d2//' '//force_covariance//' '//force_difference)
      end do
      mainlog%prefix = ''
   end if

   if (do_teach) then
      call teach(rscov, do_optimise, use_qr)
   end if

   ! Testing the implementation of overlap_func and doverlap_func. 
   if(.false.) then
      am%at1 => rscov%data(1)   
      am%at2 => rscov%data(1)   
      am%rscov => rscov
      am_data_size = size(transfer(am,am_mold))
      allocate(am_data(am_data_size))
      am_data = transfer(am,am_mold)
      q = random_quaternion()
      !q = (/.99619_dp, .08716_dp, + 0.00000_dp, 0.00000_dp/)
      call print('starting test gradient')
      call verbosity_push(PRINT_NERD)
      call n_test_gradient(q,overlap_func,doverlap_func,data=am_data)
      print*, test_gradient(q,overlap_func,doverlap_func,data=am_data)
      call verbosity_pop()
   endif

   if (test_data) then
      mainlog%prefix = 'DATA'
      do i=1,size(rscov%data)
         call print_title('predict data('//i//')')
         
         call system_timer('predict_force')
         call predict(rscov, rscov%data(i), force, error)
         call system_timer('predict_force')
         call assign_property_pointer(rscov%data(i), 'force', force_ptr)
         call print('RESULT '//force_ptr(:,1)//' '//force//' '//' '//maxval(abs(force_ptr(:,1)-force))//' '//error)
      end do
      mainlog%prefix = ''
   end if

   if (trim(test_file) /= '') then
      call print_title('Test config')

      mainlog%prefix = 'TEST'
      call initialise(in, test_file, INPUT)
      do i=1,in%n_frame
         call read(in, at)
         call centre_on_atom_1(at)
         call system_timer('predict_force')
         call predict(rscov, at, force, error)
         call system_timer('predict_force')
         call assign_property_pointer(at, 'force', force_ptr)
         do j=1,size(rscov%data)
            call print('RESULT '//j//' q='//rscov%q_min(:,j)//' k='//rscov%k(j))
         end do
         call print('RESULT '//force_ptr(:,1)//' '//force//' '//maxval(abs(force_ptr(:,1)-force))//' '//error)
      end do
      mainlog%prefix = ''
      call finalise(in)
   end if

   if (test_data_2) then
      mainlog%prefix = 'DATA'
      do i=1,size(rscov%data)
         call print_title('predict_2 data('//i//')')
         
         call system_timer('predict_force')
         call predict_2(rscov, rscov%data(i), force, error)
         call system_timer('predict_force')
         call assign_property_pointer(rscov%data(i), 'force', force_ptr)
         call print('RESULT_2 '//force_ptr(:,1)//' '//force//' '//' '//maxval(abs(force_ptr(:,1)-force))//' '//error)
      end do
      mainlog%prefix = ''
      
      call print_title('Test config - method 2')

      mainlog%prefix = 'TEST'
      call initialise(in, test_file, INPUT)
      do i=1,in%n_frame
         call read(in, at)
         call centre_on_atom_1(at)
         call system_timer('predict_force')
         call predict_2(rscov, at, force, error)
         call system_timer('predict_force')
         call assign_property_pointer(at, 'force', force_ptr)
         do j=1,size(rscov%data)
            call print('RESULT_2 '//j//' q='//rscov%q_min(:,j)//' k='//rscov%k(j))
         end do
         call print('RESULT_2 '//force_ptr(:,1)//' '//force//' '//maxval(abs(force_ptr(:,1)-force))//' '//error)
      end do
      mainlog%prefix = ''
      call finalise(in)      

   end if
   call system_finalise()
   
   call finalise(rscov)

end program test_cova
