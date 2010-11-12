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

program teach_sparse_program

  use libatoms_module
  use descriptors_module
  use gp_predict_module
  use gp_teach_module
  use clustering_module
  use teach_sparse_mod

  implicit none

  integer, parameter :: SPARSE_LENGTH = 10000
  integer, parameter :: SPARSE_N_FIELDS = 2000
  integer, parameter :: THETA_LENGTH = 10000
  real(dp), parameter :: THETA_MIN = 0.000000001

  type(teach_sparse) :: main_teach_sparse
  type(inoutput) :: bispectrum_inout, theta_inout, sparse_inout
  type(Dictionary) :: params
  type(gp_sparse) :: gp_sp

  character(len=FIELD_LENGTH) :: verbosity
  character(len=FIELD_LENGTH) :: qw_cutoff_string, qw_cutoff_f_string, qw_cutoff_r1_string, &
  theta_file, sparse_file, bispectrum_file, m_sparse_in_type_string, tmp
  real(dp) :: mem_required, mem_total, mem_free
  logical :: has_e0, has_f0, has_theta_file, has_sparse_file, has_bispectrum_file, test_gp_gradient_result

  character(len=FIELD_LENGTH), dimension(99) :: qw_cutoff_fields, qw_cutoff_f_fields, qw_cutoff_r1_fields, m_sparse_in_type_fields
  character(len=SPARSE_LENGTH) :: sparse_string
  character(len=FIELD_LENGTH), dimension(:), allocatable :: sparse_string_array
  character(len=THETA_LENGTH) :: theta_string
  character(len=FIELD_LENGTH), dimension(:), allocatable :: theta_string_array
  integer :: i, j, k, l, o, dd, dt, m_sparse_in_type_num_fields, i_default, m_total
  character(len=FIELD_LENGTH) :: gp_file

  call system_initialise(verbosity=PRINT_NORMAL)
  call initialise(params)
  call param_register(params, 'at_file', PARAM_MANDATORY, main_teach_sparse%at_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'm', '50', main_teach_sparse%m, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'r_cut', '2.75', main_teach_sparse%r_cut, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'j_max', '4', main_teach_sparse%j_max, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'z0', '0.0', main_teach_sparse%z0, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'coordinates', 'bispectrum', main_teach_sparse%coordinates, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'l_max', '6', main_teach_sparse%qw_l_max, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'cutoff', '', qw_cutoff_string, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'cutoff_f', '', qw_cutoff_f_string, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'cutoff_r1', '', qw_cutoff_r1_string, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'no_q', 'F', main_teach_sparse%qw_no_q, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'no_w', 'F', main_teach_sparse%qw_no_w, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'e0', '0.0', main_teach_sparse%e0, has_value_target = has_e0, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'f0', '0.0', main_teach_sparse%f0, has_value_target = has_f0, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'sgm', '0.1 0.1 0.1', main_teach_sparse%sgm, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'dlt', '1.0', main_teach_sparse%dlt, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'theta_file', '', theta_file, has_value_target = has_theta_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'sparse_file', '', sparse_file, has_value_target = has_sparse_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'theta_fac', '1.5', main_teach_sparse%theta_fac, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'do_sigma', 'F', main_teach_sparse%do_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'do_delta', 'F', main_teach_sparse%do_delta, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'do_theta', 'F', main_teach_sparse%do_theta, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'do_sparx', 'F', main_teach_sparse%do_sparx, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'do_f0', 'F', main_teach_sparse%do_f0, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'do_theta_fac', 'F', main_teach_sparse%do_theta_fac, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'do_cluster', 'F', main_teach_sparse%do_cluster, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'do_pivot', 'F', main_teach_sparse%do_pivot, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'min_steps', '10', main_teach_sparse%min_steps, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'min_save', '0', main_teach_sparse%min_save, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'do_test_gp_gradient', 'F', main_teach_sparse%do_test_gp_gradient, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'bispectrum_file', '', bispectrum_file, has_value_target = has_bispectrum_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'ip_args', '', main_teach_sparse%ip_args, has_value_target = main_teach_sparse%do_core, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'energy_property_name', 'energy', main_teach_sparse%energy_property_name, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'force_property_name', 'force', main_teach_sparse%force_property_name, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'virial_property_name', 'virial', main_teach_sparse%virial_property_name, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'm_sparse_in_type','',m_sparse_in_type_string,has_value_target = main_teach_sparse%has_m_sparse_in_type, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'do_sparse', 'T', main_teach_sparse%do_sparse, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'do_pca', 'F', main_teach_sparse%do_pca, help_string='PCA analysis is performed on input data')
  call param_register(params, 'verbosity', 'NORMAL', verbosity, help_string="No help yet.  This source file was $LastChangedBy$")


  if (.not. param_read_args(params, do_check = .true., command_line=main_teach_sparse%command_line)) then
     call print("Usage: teach_sparse [at_file=file] [m=50] &
     & [r_cut=2.75] [j_max=4] [z0=0.0] [coordinates={bispectrum,qw,water_monomer,water_dimer}] [l_max=6] [cutoff={:}] [cutoff_f={:}] [cutoff_r1={:}] [no_q] [no_w] &
     & [e0=0.0] [f0=avg] [sgm={0.1 0.1 0.1}] [dlt=1.0] [theta_file=file] [sparse_file=file] [theta_fac=3.0] &
     & [do_sigma=F] [do_delta=F] [do_theta=F] [do_sparx=F] [do_f0=F] [do_theta_fac=F] &
     & [do_cluster=F] [do_pivot=F] [min_steps=10] [min_save=0] &
     & [do_test_gp_gradient=F] [bispectrum_file=file] [ip_args={}] &
     & [energy_property_name=energy] [force_property_name=force] [virial_property_name=virial] [m_sparse_in_type={default:200:crystal:300}] [do_sparse=T] [verbosity=NORMAL]")
     call system_abort('Exit: Mandatory argument(s) missing...')
  endif
  call finalise(params)

  select case(verbosity)
    case ("NORMAL")
      call verbosity_push(PRINT_NORMAL)
    case ("VERBOSE")
      call verbosity_push(PRINT_VERBOSE)
    case ("NERD")
      call verbosity_push(PRINT_NERD)
    case ("ANAL")
      call verbosity_push(PRINT_ANAL)
    case default
      call system_abort("confused by verbosity " // trim(verbosity))
  end select

  if( count( (/has_e0,has_f0/) ) > 1 ) &
  & call print('Warning - you have specified both e0 and f0 - careful!')

  if( count( (/has_sparse_file,main_teach_sparse%do_cluster,main_teach_sparse%do_pivot/) ) > 1 ) &
  & call system_abort('There has been more than one method specified for sparsification.')
     
  if( main_teach_sparse%has_m_sparse_in_type ) then
     call parse_string(m_sparse_in_type_string,':',m_sparse_in_type_fields,m_sparse_in_type_num_fields)

     ! find "default" if present
     i_default = 0
     do i = 1, m_sparse_in_type_num_fields, 2
        if( lower_case(trim(m_sparse_in_type_fields(i))) == "default" ) i_default = i
     enddo


     if( i_default > 1 ) then
        ! default is somewhere else, swap it to the first place
        tmp = m_sparse_in_type_fields(1)
        m_sparse_in_type_fields(1) = m_sparse_in_type_fields(i_default)
        m_sparse_in_type_fields(i_default) = tmp

        tmp = m_sparse_in_type_fields(2)
        m_sparse_in_type_fields(2) = m_sparse_in_type_fields(i_default+1)
        m_sparse_in_type_fields(i_default+1) = tmp

        i_default = 1
     endif

     if( i_default == 0 ) then
        ! no default present in the string
        allocate(main_teach_sparse%m_sparse_in_type(0 : m_sparse_in_type_num_fields/2))
     else
        ! default is present, it's the first
        allocate(main_teach_sparse%m_sparse_in_type(0 : m_sparse_in_type_num_fields/2-1))
     endif

     m_total = 0
     do i = i_default + 1, m_sparse_in_type_num_fields/2 ! don't include the default for now
        main_teach_sparse%m_sparse_in_type(i-i_default)%type = m_sparse_in_type_fields(2*(i-1)+1)
        main_teach_sparse%m_sparse_in_type(i-i_default)%m = string_to_int(m_sparse_in_type_fields(2*i))
        m_total = m_total + main_teach_sparse%m_sparse_in_type(i-i_default)%m
     enddo

     main_teach_sparse%m_sparse_in_type(0)%type = 'default'

     if( i_default == 0 ) then
        if( m_total > main_teach_sparse%m ) then
           call print_warning('Total number of sparse points for different types '//m_total//' &
           is greater than number of sparse points '//main_teach_sparse%m//', latter is overwritten.')
           main_teach_sparse%m_sparse_in_type(0)%m = 0
           main_teach_sparse%m = m_total
        else
           main_teach_sparse%m_sparse_in_type(0)%m = main_teach_sparse%m - m_total
        endif
     else
        main_teach_sparse%m_sparse_in_type(0)%m = string_to_int(m_sparse_in_type_fields(2))
        m_total = m_total + main_teach_sparse%m_sparse_in_type(0)%m
        if( m_total > main_teach_sparse%m ) then
           call print_warning('Total number of sparse points for different types '//m_total//' &
           is greater than number of sparse points '//main_teach_sparse%m//', latter is overwritten.')
           main_teach_sparse%m = m_total
        else
           main_teach_sparse%m_sparse_in_type(0)%m = main_teach_sparse%m_sparse_in_type(0)%m + main_teach_sparse%m - m_total
        endif
     endif
     call print('Sparse points per pre-defined types of configurations')
     do i = lbound(main_teach_sparse%m_sparse_in_type,dim=1), ubound(main_teach_sparse%m_sparse_in_type,dim=1)
        call print(""//trim(main_teach_sparse%m_sparse_in_type(i)%type)//"   "//main_teach_sparse%m_sparse_in_type(i)%m)
     enddo
  endif



  main_teach_sparse%coordinates = lower_case(trim(main_teach_sparse%coordinates))
  select case(trim(main_teach_sparse%coordinates))
  case('water_monomer')
     main_teach_sparse%d = 3
  case('water_dimer')
     main_teach_sparse%d = WATER_DIMER_D
  case('qw')
     main_teach_sparse%qw_cutoff = 0.0_dp
     main_teach_sparse%qw_cutoff_f = 0
     main_teach_sparse%qw_cutoff_r1 = 0.0_dp
     call parse_string(qw_cutoff_string, ':', qw_cutoff_fields, main_teach_sparse%qw_f_n)
     call parse_string(qw_cutoff_f_string, ':', qw_cutoff_f_fields, main_teach_sparse%qw_f_n)
     call parse_string(qw_cutoff_r1_string, ':', qw_cutoff_r1_fields, main_teach_sparse%qw_f_n)
     do i = 1, main_teach_sparse%qw_f_n
        main_teach_sparse%qw_cutoff(i) = string_to_real(qw_cutoff_fields(i))
        main_teach_sparse%qw_cutoff_f(i) = string_to_int(qw_cutoff_f_fields(i))
        main_teach_sparse%qw_cutoff_r1(i) = string_to_real(qw_cutoff_r1_fields(i))
     enddo
     main_teach_sparse%r_cut = maxval(main_teach_sparse%qw_cutoff(1:main_teach_sparse%qw_f_n))
  case('bispectrum')
     main_teach_sparse%z0 = max(1.0_dp,main_teach_sparse%z0) * main_teach_sparse%r_cut/(PI-0.02_dp)
     main_teach_sparse%d = j_max2d(main_teach_sparse%j_max)
  case default
     call system_abort('Unknown coordinates '//trim(main_teach_sparse%coordinates))
  endselect

  if(main_teach_sparse%do_core) &
     call read(main_teach_sparse%quip_string, "quip_params.xml",keep_lf=.true.)

  call teach_n_from_xyz(main_teach_sparse)

  if (.not. has_f0 .and. .not. has_e0) then
     call e0_avg_from_xyz(main_teach_sparse)
  end if

  call w_Z_from_xyz(main_teach_sparse)

  call teach_data_from_xyz(main_teach_sparse)

  if( has_bispectrum_file ) then
     call initialise(bispectrum_inout,bispectrum_file,action=OUTPUT)
     do i = 1, size(main_teach_sparse%x,2)
        write(bispectrum_inout%unit,"("//main_teach_sparse%d//"f16.8)") main_teach_sparse%x(:,i)
     enddo
     call finalise(bispectrum_inout)
  endif

  allocate(main_teach_sparse%dlta(main_teach_sparse%n_species))

  main_teach_sparse%dlta = main_teach_sparse%dlt 

  if(main_teach_sparse%do_pca) then
     allocate(main_teach_sparse%pca_matrix(main_teach_sparse%d,main_teach_sparse%d), &
     main_teach_sparse%pca_mean(main_teach_sparse%d))
     call pca(main_teach_sparse%x,main_teach_sparse%pca_mean,main_teach_sparse%pca_matrix)

     do i = 1, main_teach_sparse%nn
        main_teach_sparse%x(:,i) = matmul(main_teach_sparse%x(:,i)-main_teach_sparse%pca_mean, &
        main_teach_sparse%pca_matrix)
     enddo

     do i = 1, main_teach_sparse%n
        main_teach_sparse%xd(:,i) = matmul(main_teach_sparse%xd(:,i),main_teach_sparse%pca_matrix)
     enddo
  endif

  if(main_teach_sparse%do_sparse) then
     if( has_sparse_file ) then
        allocate(sparse_string_array(SPARSE_N_FIELDS))
        call initialise(sparse_inout,sparse_file)
        read(sparse_inout%unit,'(a)') sparse_string
        call parse_string(sparse_string,' ',sparse_string_array,main_teach_sparse%m)
        allocate(main_teach_sparse%r(main_teach_sparse%m))
        do i = 1, main_teach_sparse%m
           main_teach_sparse%r(i) = string_to_int(sparse_string_array(i))
        enddo
        deallocate(sparse_string_array)
        call finalise(sparse_inout)
     elseif(main_teach_sparse%do_cluster) then
        allocate(main_teach_sparse%r(main_teach_sparse%m))
        call bisect_kmedoids(main_teach_sparse%x,main_teach_sparse%m,med=main_teach_sparse%r,theta_fac=main_teach_sparse%theta_fac)
     elseif(main_teach_sparse%do_pivot) then
        allocate(main_teach_sparse%r(main_teach_sparse%m))
        call pivot(main_teach_sparse%x, main_teach_sparse%r,theta_fac=main_teach_sparse%theta_fac)
     else
        allocate(main_teach_sparse%r(main_teach_sparse%m))
        call fill_random_integer(main_teach_sparse%r,size(main_teach_sparse%x,2))
     endif
     call sort_array(main_teach_sparse%r)

     call print('')
     call print('Atomic environments used in sparsification')
     call print(main_teach_sparse%r)
     call print('')
  else
     allocate(main_teach_sparse%r(main_teach_sparse%nn))
     main_teach_sparse%r = (/ (i,i=1,main_teach_sparse%nn) /)
  endif

  allocate(main_teach_sparse%theta(main_teach_sparse%d,main_teach_sparse%n_species))

  if( has_theta_file ) then
     allocate(theta_string_array(main_teach_sparse%d))
     call initialise(theta_inout,theta_file)
     read(theta_inout%unit,'(a)') theta_string
     call parse_string(theta_string,' ',theta_string_array,dt)
     if(main_teach_sparse%d /= dt) call system_abort('File '//trim(theta_file)//' does not contain the right number of hyperparameters')
     do k = 1, main_teach_sparse%n_species  
        do dd = 1, main_teach_sparse%d
           main_teach_sparse%theta(dd,k) = string_to_real(theta_string_array(dd+(k-1)*main_teach_sparse%d))
        enddo
     enddo
     deallocate(theta_string_array)
     call finalise(theta_inout)
  else
     do k = 1, main_teach_sparse%n_species
        do dd = 1, main_teach_sparse%d
!           theta(dd,k) = ( maxval(x(dd,:),mask=(xz(:)==species_Z(k))) - minval(x(dd,:),mask=(xz(:)==species_Z(k))) )
           main_teach_sparse%theta(dd,k) = ( maxval(main_teach_sparse%x(dd,main_teach_sparse%r),&
           mask=(main_teach_sparse%xz(main_teach_sparse%r)==main_teach_sparse%species_Z(k))) &
           - minval(main_teach_sparse%x(dd,main_teach_sparse%r),&
           mask=(main_teach_sparse%xz(main_teach_sparse%r)==main_teach_sparse%species_Z(k))) )
!           theta(dd) = sqrt( & !take square root
!                          & sum( x(dd,:)**2 ) / size(x(dd,:)) - &
!                          & (sum( x(dd,:) ) / size(x(dd,:)))**2 )

           if( main_teach_sparse%theta(dd,k) >= THETA_MIN ) then
              main_teach_sparse%theta(dd,k) = main_teach_sparse%theta_fac*main_teach_sparse%theta(dd,k)
           else
              main_teach_sparse%theta(dd,k) = 1.0_dp
           endif
        enddo
     enddo
  endif
  ! Stop execution if required memory is greater than the available memory. 
  ! The biggest arrays allocated are 2*sr*(nx+nxd), where sr is the
  ! number of sparse points, nx and nxd are the number of bispectra and partial
  ! derivatives.
  if(main_teach_sparse%do_sparse) then
     mem_required = real(size(main_teach_sparse%r),dp) * (real(size(main_teach_sparse%xf),dp) &
     + real(size(main_teach_sparse%xdf),dp)) * real(dp,dp) / (1024.0_dp**3)
#ifdef SPEEDOPT
     mem_required = 2.0_dp * mem_required
#endif                  
  else
     mem_required = 3.0_dp*real(size(main_teach_sparse%x,2),dp)**2 * real(dp,dp) / (1024.0_dp**3)
  endif

  call mem_info(mem_total,mem_free)
  mem_total = mem_total / (1024.0_dp**3)
  mem_free = mem_free / (1024.0_dp**3)

  call print('Memory required (approx.): '//mem_required//' GB')
  if( mem_required > mem_total ) call system_abort('Required memory ('//mem_required//' GB) exceeds available memory ('//mem_total//' GB).')

  if(main_teach_sparse%do_sparse) then
     call gp_sparsify(gp_sp,main_teach_sparse%r,&
     main_teach_sparse%sgm,main_teach_sparse%dlta,main_teach_sparse%theta,&
     main_teach_sparse%yf,main_teach_sparse%ydf,main_teach_sparse%x,main_teach_sparse%xd,&
     main_teach_sparse%xf,main_teach_sparse%xdf,main_teach_sparse%lf,main_teach_sparse%ldf,&
     main_teach_sparse%xz,main_teach_sparse%species_Z,(/(main_teach_sparse%f0,i=1,main_teach_sparse%n_species)/),&
     main_teach_sparse%target_type)

     call print('')
     call print('theta')
     do l = 1, size(gp_sp%theta, 2)
        do o = 1, size(gp_sp%theta, 1)
           call print(real(gp_sp%theta(o,l),kind=dp))
        enddo
     enddo
     call print('delta')
     call print(real(gp_sp%delta,kind=dp))
     call print('sigma')
     call print(real(gp_sp%sigma,kind=dp))
     call print('')

     call enable_timing()
     
     if( main_teach_sparse%do_test_gp_gradient ) then
        call verbosity_push(PRINT_NERD)
        test_gp_gradient_result = test_gp_gradient(gp_sp,&
        sigma=main_teach_sparse%do_sigma,delta=main_teach_sparse%do_delta,&
        theta=main_teach_sparse%do_theta,sparx=main_teach_sparse%do_sparx,&
        f0=main_teach_sparse%do_f0,theta_fac=main_teach_sparse%do_theta_fac)
        call verbosity_pop()
     endif

     ! Conjugate gradient minimiser's counter starts at 1, and stops when it reaches min_steps,
     ! so if min_steps is equal to 1, no iterations are made!
    
     if (main_teach_sparse%min_save == 0) main_teach_sparse%min_save = main_teach_sparse%min_steps

     k = 0
     do i = 1, ((main_teach_sparse%min_steps / main_teach_sparse%min_save) + 1)
        if (k == main_teach_sparse%min_steps) exit

        if ((main_teach_sparse%min_steps - k) >= main_teach_sparse%min_save) then
           j = minimise_gp_gradient(gp_sp,max_steps=(main_teach_sparse%min_save + 1),&
           sigma=main_teach_sparse%do_sigma,delta=main_teach_sparse%do_delta,&
           theta=main_teach_sparse%do_theta,sparx=main_teach_sparse%do_sparx,&
           f0=main_teach_sparse%do_f0,theta_fac=main_teach_sparse%do_theta_fac)
           k = k + main_teach_sparse%min_save
        elseif ((main_teach_sparse%min_steps - k) < main_teach_sparse%min_save) then
           j = minimise_gp_gradient(gp_sp,max_steps=(main_teach_sparse%min_steps - k + 1),&
           sigma=main_teach_sparse%do_sigma,delta=main_teach_sparse%do_delta,&
           theta=main_teach_sparse%do_theta,sparx=main_teach_sparse%do_sparx,&
           f0=main_teach_sparse%do_f0,theta_fac=main_teach_sparse%do_theta_fac)

           k = main_teach_sparse%min_steps
        endif

        call print('')
        call print(k // ' iterations completed:')
        call print('theta')
        do l = 1, size(gp_sp%theta, 2)
           do o = 1, size(gp_sp%theta, 1)
              call print(real(gp_sp%theta(o,l),kind=dp))
           enddo
        enddo
        call print('delta')
        call print(real(gp_sp%delta,kind=dp))
        call print('sigma')
        call print(real(gp_sp%sigma,kind=dp))
        call print('f0')
        call print(real(gp_sp%f0,kind=dp))
        call print('')

        call initialise(main_teach_sparse%my_gp,gp_sp)

        gp_file = 'gp_'//main_teach_sparse%m//'_'//k//'.xml'

        call teach_sparse_print_xml(main_teach_sparse,gp_file)

        call system_command('ln -fs '//trim(gp_file)//' gp.xml')

        call finalise(main_teach_sparse%my_gp)
     enddo
  else
     call initialise(main_teach_sparse%my_gp, main_teach_sparse%sgm(1), &
     main_teach_sparse%dlta(1), main_teach_sparse%theta, main_teach_sparse%f0, &
     main_teach_sparse%yf, main_teach_sparse%x)

     if( main_teach_sparse%do_test_gp_gradient ) then
        call verbosity_push(PRINT_NERD)
        test_gp_gradient_result = test_gp_simple_gradient(main_teach_sparse%my_gp,&
        sigma=main_teach_sparse%do_sigma,delta=main_teach_sparse%do_delta,&
        theta=main_teach_sparse%do_theta, &
        f0=main_teach_sparse%do_f0,theta_fac=main_teach_sparse%do_theta_fac)
        call verbosity_pop()
     endif
     j = minimise_gp_simple_gradient(main_teach_sparse%my_gp,max_steps=main_teach_sparse%min_steps ,&
         sigma=main_teach_sparse%do_sigma,delta=main_teach_sparse%do_delta,&
         theta=main_teach_sparse%do_theta,&
         f0=main_teach_sparse%do_f0,theta_fac=main_teach_sparse%do_theta_fac)

     gp_file = 'gp_'//main_teach_sparse%my_gp%n//'.xml'

     call teach_sparse_print_xml(main_teach_sparse,gp_file)
     call system_command('ln -fs '//trim(gp_file)//' gp.xml')
     call finalise(main_teach_sparse%my_gp)
  endif


  call print("model parameters:")
  call print("r_cut     = "//main_teach_sparse%r_cut)
  select case(trim(main_teach_sparse%coordinates))
  case('water_monomer','water_dimer')
  case('qw')
     call print("l_max     = "//main_teach_sparse%qw_l_max)
     call print("cutoff    = "//qw_cutoff_string)
     call print("cutoff_f  = "//qw_cutoff_f_string)
     call print("cutoff_r1 = "//qw_cutoff_r1_string)
     call print("q         = "//(.not. main_teach_sparse%qw_no_q))
     call print("w         = "//(.not. main_teach_sparse%qw_no_w))
  case('bispectrum')
     call print("j_max     = "//main_teach_sparse%j_max)
     call print("z0        = "//main_teach_sparse%z0)
  case default
     call system_abort('Unknown coordinates '//trim(main_teach_sparse%coordinates))
  endselect

  call print("n_species = "//main_teach_sparse%n_species)
  call print("species_Z = "//main_teach_sparse%species_Z)
  call print("w         = "//main_teach_sparse%w_Z(main_teach_sparse%species_Z))
  call print("e0        = "//main_teach_sparse%e0)

  call finalise(gp_sp)

  call system_finalise()

end program teach_sparse_program
