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

  character(len=STRING_LENGTH) :: verbosity
  character(len=STRING_LENGTH) :: qw_cutoff_string, qw_cutoff_f_string, qw_cutoff_r1_string, &
  theta_file, sparse_file, bispectrum_file, config_type_hypers_string, tmp
  real(dp) :: mem_required, mem_total, mem_free, z0_ratio
  logical :: has_e0, has_f0, has_theta_file, has_sparse_file, has_bispectrum_file, test_gp_gradient_result

  character(len=STRING_LENGTH), dimension(99) :: qw_cutoff_fields, qw_cutoff_f_fields, qw_cutoff_r1_fields, config_type_hypers_fields
  character(len=SPARSE_LENGTH) :: sparse_string
  character(len=STRING_LENGTH), dimension(:), allocatable :: sparse_string_array
  character(len=THETA_LENGTH) :: theta_string
  character(len=STRING_LENGTH), dimension(:), allocatable :: theta_string_array
  integer :: i, j, k, l, o, dd, dt, config_type_hypers_num_fields, i_default, m_total, li, ui, gp_teach_memory
  integer, dimension(:), allocatable :: type_indices, sparse_points
  integer, dimension(:,:), allocatable :: permutation
  character(len=STRING_LENGTH) :: gp_file
  real(dp), dimension(:), allocatable :: theta

  call system_initialise(verbosity=PRINT_NORMAL)
  call initialise(params)
  call param_register(params, 'at_file', PARAM_MANDATORY, main_teach_sparse%at_file, help_string="XYZ file with teaching configurations")
  call param_register(params, 'm', '50', main_teach_sparse%m, help_string="Number of sparse environments")
  call param_register(params, 'r_cut', '2.75', main_teach_sparse%r_cut, help_string="Cutoff for GAP environment")
  call param_register(params, 'j_max', '4', main_teach_sparse%j_max, help_string="Max of expansion of bispectrum, i.e. resulution")
  call param_register(params, 'z0_ratio', '0.0', z0_ratio, help_string="Developer option. Ratio of radius of 4D projection sphere times PI and the cutoff.")
  call param_register(params, 'coordinates', 'bispectrum', main_teach_sparse%coordinates, help_string="{bispectrum,qw,water_monomer,water_dimer} representation to use")
  call param_register(params, 'l_max', '6', main_teach_sparse%qw_l_max, help_string="Maximum sphearical harmonic expansion")
  call param_register(params, 'cutoff', '', qw_cutoff_string, help_string="Cutoffs for the radial basis functions")
  call param_register(params, 'cutoff_f', '', qw_cutoff_f_string, help_string="Types of the radial basis functions: 1 - Fermi, 2 - Gaussian, 3 - Skewed Gaussian")
  call param_register(params, 'cutoff_r1', '', qw_cutoff_r1_string, help_string="Characteristic width of the cutoff function")
  call param_register(params, 'no_q', 'F', main_teach_sparse%qw_no_q, help_string="Turn off Qs")
  call param_register(params, 'no_w', 'F', main_teach_sparse%qw_no_w, help_string="Turn off Ws")

  call param_register(params, 'cosnx_l_max', '6', main_teach_sparse%cosnx_l_max, help_string="Maximum cos(l theta) expansion")
  call param_register(params, 'cosnx_n_max', '3', main_teach_sparse%cosnx_n_max, help_string="Maximum number of radial basis functions")
  call param_register(params, 'use_rdf', 'F', main_teach_sparse%use_rdf, help_string="Whether to optimise the radial basis set with respect of the distribution of distances")

  call param_register(params, 'e0', '0.0', main_teach_sparse%e0, has_value_target = has_e0, help_string="Value to be subtracted from energies before fittint (and added back on after prediction)")
  call param_register(params, 'f0', '0.0', main_teach_sparse%f0, has_value_target = has_f0, help_string="Mean of the gp predictor, defaults to data mean. another useful value is 0.0")
  call param_register(params, 'sgm', '0.1 0.1 0.1', main_teach_sparse%sgm, help_string="error in [energies forces virials]")
  call param_register(params, 'dlt', '1.0', main_teach_sparse%dlt, help_string="Range of GP")
  call param_register(params, 'theta_file', '', theta_file, has_value_target = has_theta_file, help_string="Width of Gaussians from a file. There should be as many real numbers as the number of descriptors, in a single line")
  call param_register(params, 'sparse_file', '', sparse_file, has_value_target = has_sparse_file, help_string="Sparse environments from a file. Integers, in single line")
  call param_register(params, 'theta_fac', '1.5', main_teach_sparse%theta_fac, help_string="Width of Gaussians, determined from multiplying the range of each descriptor by theta_fac")
  call param_register(params, 'do_sigma', 'F', main_teach_sparse%do_sigma, help_string="Likelihood optimization with respect of sigma")
  call param_register(params, 'do_delta', 'F', main_teach_sparse%do_delta, help_string="Likelihood optimization with respect of delta")
  call param_register(params, 'do_theta', 'F', main_teach_sparse%do_theta, help_string="Likelihood optimization with respect of theta")
  call param_register(params, 'do_sparx', 'F', main_teach_sparse%do_sparx, help_string="Likelihood optimization with respect of sparse points")
  call param_register(params, 'do_f0', 'F', main_teach_sparse%do_f0, help_string="Likelihood optimization with respect of f0 (offset)")
  call param_register(params, 'do_theta_fac', 'F', main_teach_sparse%do_theta_fac, help_string="Likelihood optimization with respect of theta_fac")
  call param_register(params, 'do_cluster', 'F', main_teach_sparse%do_cluster, help_string="Do k-means clustering of data points to initialise the sparsifier")
  call param_register(params, 'do_pivot', 'F', main_teach_sparse%do_pivot, help_string="Initialise the sparsifier by looking at the covariance matrix")
  call param_register(params, 'min_steps', '10', main_teach_sparse%min_steps, help_string="How many iteration steps in the likelihood optimization")
  call param_register(params, 'min_save', '0', main_teach_sparse%min_save, help_string="How often print out a GAP xml file during the optimization")
  call param_register(params, 'do_test_gp_gradient', 'F', main_teach_sparse%do_test_gp_gradient, help_string="Developer option. Test likelihood gradient")
  call param_register(params, 'bispectrum_file', '', bispectrum_file, has_value_target = has_bispectrum_file, help_string="Developer option. Print out descriptors of each atom in a text file.")
  call param_register(params, 'ip_args', '', main_teach_sparse%ip_args, has_value_target = main_teach_sparse%do_core, help_string=" QUIP init string for a potential to subtract from data (and added back after prediction)")
  call param_register(params, 'energy_parameter_name', 'energy', main_teach_sparse%energy_parameter_name, help_string="Name of energy property in the at_file that describe the data")
  call param_register(params, 'force_parameter_name', 'force', main_teach_sparse%force_parameter_name, help_string="Name of force property in the at_file that describe the data")
  call param_register(params, 'virial_parameter_name', 'virial', main_teach_sparse%virial_parameter_name, help_string="Name of virial property in the at_file that describe the data")
  call param_register(params, 'config_type_parameter_name', 'config_type', main_teach_sparse%config_type_parameter_name, help_string="Identifier of property determining the type of input data in the at_file")
  call param_register(params, 'mask_name', '', main_teach_sparse%mask_name, help_string="Name of logical property to use as a mask for atomic enviornments")
  call param_register(params, 'config_type_hypers','',config_type_hypers_string,has_value_target = main_teach_sparse%has_config_type_hypers, help_string="How many sparse points to choose for each type of data and what sigma parameters to associate with the type. E.g. {liquid:200:0.01:0.2:0.02:crystal:300:0.001:0.1:0.001} will select 300 sparse points from crystal type data and 200 from liquid")
  call param_register(params, 'do_sparse', 'T', main_teach_sparse%do_sparse, help_string="Do sparsification or regular GP. Latter: no derivative information is used")
  call param_register(params, 'do_pca', 'F', main_teach_sparse%do_pca, help_string='PCA analysis is performed on input data')
  call param_register(params, 'mark_sparse_atoms', '', main_teach_sparse%mark_sparse_atoms, has_value_target = main_teach_sparse%do_mark_sparse_atoms, help_string="Reprints the original xyz file after sparsification process. sparse propery added, true for atoms associated with a sparse point.")
  call param_register(params, 'verbosity', 'NORMAL', verbosity, help_string="Verbosity control")


  if (.not. param_read_args(params, command_line=main_teach_sparse%command_line)) then
     call print("Usage: teach_sparse [at_file=file] [m=50] &
     [r_cut=2.75] [j_max=4] [z0_ratio=0.0] [coordinates={bispectrum,qw,water_monomer,water_dimer,cosnx}] [l_max=6] [cutoff={:}] [cutoff_f={:}] [cutoff_r1={:}] [no_q] [no_w] &
     [cosnx_l_max=6] [cosnx_n_max=3] [use_rdf=T] &
     [e0=0.0] [f0=avg] [sgm={0.1 0.1 0.1}] [dlt=1.0] [theta_file=file] [sparse_file=file] [theta_fac=3.0] &
     [do_sigma=F] [do_delta=F] [do_theta=F] [do_sparx=F] [do_f0=F] [do_theta_fac=F] &
     [do_cluster=F] [do_pivot=F] [min_steps=10] [min_save=0] &
     [do_test_gp_gradient=F] [bispectrum_file=file] [ip_args={}] &
     [energy_parameter_name=energy] [force_parameter_name=force] [virial_parameter_name=virial] &
     [config_type_hypers={liquid:200:0.01:0.2:0.02:crystal:300:0.001:0.1:0.001}] [do_sparse=T] [verbosity=NORMAL]")
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
     
  if( main_teach_sparse%has_config_type_hypers ) then
     call parse_string(config_type_hypers_string,':',config_type_hypers_fields,config_type_hypers_num_fields)

     ! find "default" if present
     i_default = 0
     do i = 1, config_type_hypers_num_fields, 5
        if( lower_case(trim(config_type_hypers_fields(i))) == "default" ) i_default = i
     enddo


     if( i_default > 1 ) then
        ! default is somewhere else, swap it to the first place

        do i = 1, 5
           tmp = config_type_hypers_fields(i)
           config_type_hypers_fields(i) = config_type_hypers_fields(i_default+i-1)
           config_type_hypers_fields(i_default+i-1) = tmp
        enddo

        i_default = 1
     endif

     if( i_default == 0 ) then
        ! no default present in the string
        allocate(main_teach_sparse%config_type_hypers(0 : config_type_hypers_num_fields/5))
     else
        ! default is present, it's the first
        allocate(main_teach_sparse%config_type_hypers(0 : config_type_hypers_num_fields/5-1))
     endif

     m_total = 0
     do i = i_default + 1, config_type_hypers_num_fields/5 ! don't include the default for now
        main_teach_sparse%config_type_hypers(i-i_default)%type = config_type_hypers_fields(5*(i-1)+1)
        main_teach_sparse%config_type_hypers(i-i_default)%m = string_to_int(config_type_hypers_fields(5*(i-1)+2))
        main_teach_sparse%config_type_hypers(i-i_default)%sgm(1) = string_to_real(config_type_hypers_fields(5*(i-1)+3))
        main_teach_sparse%config_type_hypers(i-i_default)%sgm(2) = string_to_real(config_type_hypers_fields(5*(i-1)+4))
        main_teach_sparse%config_type_hypers(i-i_default)%sgm(3) = string_to_real(config_type_hypers_fields(5*(i-1)+5))
        m_total = m_total + main_teach_sparse%config_type_hypers(i-i_default)%m
     enddo

     main_teach_sparse%config_type_hypers(0)%type = 'default'
     main_teach_sparse%config_type_hypers(0)%sgm = main_teach_sparse%sgm

     if( i_default == 0 ) then
        if( m_total > main_teach_sparse%m ) then
           call print_warning('Total number of sparse points for different types '//m_total//' &
           is greater than number of sparse points '//main_teach_sparse%m//', latter is overwritten.')
           main_teach_sparse%config_type_hypers(0)%m = 0
           main_teach_sparse%m = m_total
        else
           main_teach_sparse%config_type_hypers(0)%m = main_teach_sparse%m - m_total
        endif
     else
        main_teach_sparse%config_type_hypers(0)%m = string_to_int(config_type_hypers_fields(2))
        m_total = m_total + main_teach_sparse%config_type_hypers(0)%m
        if( m_total > main_teach_sparse%m ) then
           call print_warning('Total number of sparse points for different types '//m_total//' &
           is greater than number of sparse points '//main_teach_sparse%m//', latter is overwritten.')
           main_teach_sparse%m = m_total
        else
           main_teach_sparse%config_type_hypers(0)%m = main_teach_sparse%config_type_hypers(0)%m + main_teach_sparse%m - m_total
        endif
     endif
     call print('Sparse points and target errors per pre-defined types of configurations')
     do i = lbound(main_teach_sparse%config_type_hypers,dim=1), ubound(main_teach_sparse%config_type_hypers,dim=1)
        call print(""//trim(main_teach_sparse%config_type_hypers(i)%type)//"   "//main_teach_sparse%config_type_hypers(i)%m//"  "//main_teach_sparse%config_type_hypers(i)%sgm)
     enddo
  else
     allocate(main_teach_sparse%config_type_hypers(0:0))
     main_teach_sparse%config_type_hypers(0)%type = "default"
     main_teach_sparse%config_type_hypers(0)%m = main_teach_sparse%m
     main_teach_sparse%config_type_hypers(0)%sgm = main_teach_sparse%sgm
  endif

  allocate(main_teach_sparse%sigma(3*size(main_teach_sparse%config_type_hypers)))
  k = 0
  do i = lbound(main_teach_sparse%config_type_hypers,dim=1), ubound(main_teach_sparse%config_type_hypers,dim=1)
     do j = 1, 3
        k = k+1
        main_teach_sparse%sigma(k) = main_teach_sparse%config_type_hypers(i)%sgm(j)
     enddo
  enddo

  main_teach_sparse%coordinates = lower_case(trim(main_teach_sparse%coordinates))
  select case(trim(main_teach_sparse%coordinates))
  case('hf_dimer')
     main_teach_sparse%d = 6
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
     main_teach_sparse%z0 = max(1.0_dp,z0_ratio) * main_teach_sparse%r_cut/(PI-0.02_dp)
     main_teach_sparse%d = j_max2d(main_teach_sparse%j_max)
  case('cosnx')
     main_teach_sparse%d = (main_teach_sparse%cosnx_l_max+1)*main_teach_sparse%cosnx_n_max
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

  if( has_bispectrum_file ) then
     call initialise(bispectrum_inout,bispectrum_file,action=OUTPUT)
     do i = 1, size(main_teach_sparse%x,2)
        write(bispectrum_inout%unit,"("//main_teach_sparse%d//"f16.8)") main_teach_sparse%x(:,i)
     enddo
     call finalise(bispectrum_inout)
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
     else

        allocate(main_teach_sparse%r(main_teach_sparse%m))
        li = 0
        ui = 0
        do i = lbound(main_teach_sparse%config_type_hypers,dim=1), ubound(main_teach_sparse%config_type_hypers,dim=1)
           if( main_teach_sparse%config_type_hypers(i)%m == 0 ) cycle
           allocate( type_indices(count(main_teach_sparse%config_type == i)), sparse_points(main_teach_sparse%config_type_hypers(i)%m) )
           type_indices = find(main_teach_sparse%config_type == i)
           call print("Choosing "//main_teach_sparse%config_type_hypers(i)%m//" sparse points for "//&
           trim(main_teach_sparse%config_type_hypers(i)%type)//" type configurations of "//size(type_indices)//" atomic environments")

           if(main_teach_sparse%do_cluster) then
              call bisect_kmedoids(main_teach_sparse%x(:,type_indices), main_teach_sparse%config_type_hypers(i)%m, med=sparse_points, theta_fac=main_teach_sparse%theta_fac)
           elseif(main_teach_sparse%do_pivot) then
              call pivot(main_teach_sparse%x(:,type_indices), sparse_points,theta_fac=main_teach_sparse%theta_fac)
           else
              call fill_random_integer(sparse_points,size(type_indices))
           endif

           li = ui + 1
           ui = ui + main_teach_sparse%config_type_hypers(i)%m
           main_teach_sparse%r(li:ui) = type_indices(sparse_points)

           deallocate(type_indices, sparse_points)
        enddo

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

  call print_sparse(main_teach_sparse)

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

  if(main_teach_sparse%do_sparse) then

     if( any( (/main_teach_sparse%do_sigma, main_teach_sparse%do_delta, main_teach_sparse%do_theta, &
     main_teach_sparse%do_sparx, main_teach_sparse%do_f0, main_teach_sparse%do_theta_fac /) ) ) then
        gp_teach_memory = 2
     else
        gp_teach_memory = 0
     endif

     call enable_timing()
     call system_timer('GP sparsify')


     select case(trim(main_teach_sparse%coordinates))
     case('hf_dimer', 'water_dimer')
        if(trim(main_teach_sparse%coordinates) .eq. 'hf_dimer') then
           allocate(permutation(6,2))
           permutation(:,1) = (/1, 2, 3, 4, 5, 6/) ! original order
           permutation(:,2) = (/4, 3, 2, 1, 5, 6/) ! swap the two HF molecules
        end if
        if(trim(main_teach_sparse%coordinates) .eq. 'water_dimer') then
           allocate(permutation(15,8))
           permutation(:,1) = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15/) ! original order
           permutation(:,2) = (/1, 3, 2, 4, 5, 7, 6, 8, 9, 10, 13, 14, 11, 12, 15/) ! swap Hs on monomer A
           permutation(:,3) = (/1, 2, 3, 5, 4, 6, 7, 9, 8, 10, 12, 11, 14, 13, 15/) ! swap Hs on monomer B
           permutation(:,4) = (/1, 3, 2, 5, 4, 7, 6, 9, 8, 10, 14, 13, 12, 11, 15/) ! swap Hs on both monomers
           permutation(:,5) = (/1, 8, 9, 6, 7, 4, 5, 2, 3, 15, 11, 13, 12, 14, 10/) ! swap monomers A and B
           permutation(:,6) = (/1, 9, 8, 6, 7, 5, 4, 2, 3, 15, 12, 14, 11, 13, 10/) ! swap monomers and Hs on monomer A
           permutation(:,7) = (/1, 8, 9, 7, 6, 4, 5, 3, 2, 15, 13, 11, 14, 12, 10/) ! swap monomers and Hs on monomer B
           permutation(:,8) = (/1, 9, 8, 7, 6, 5, 4, 3, 2, 15, 14, 12, 13, 11, 10/) ! swap monomers and Hs on both monomers
        end if

        allocate(theta(size(main_teach_sparse%theta,1)))
        theta = 0.0_dp
        do i = 1, size(permutation,2)
           theta = theta + main_teach_sparse%theta(permutation(:,i),1)
        enddo
        main_teach_sparse%theta(:,1) = theta / size(permutation,2)
        deallocate(theta)

        call gp_sparsify(gp_sp,main_teach_sparse%r,&
        main_teach_sparse%sigma,main_teach_sparse%dlta,main_teach_sparse%theta,&
        main_teach_sparse%yf,main_teach_sparse%ydf,main_teach_sparse%x,main_teach_sparse%xd,&
        main_teach_sparse%xf,main_teach_sparse%xdf,main_teach_sparse%lf,main_teach_sparse%ldf,&
        main_teach_sparse%xz,main_teach_sparse%species_Z,(/(main_teach_sparse%f0,i=1,main_teach_sparse%n_species)/),&
        main_teach_sparse%target_type, gp_teach_memory_in=gp_teach_memory,permutation_in = permutation)
     case default
        call gp_sparsify(gp_sp,main_teach_sparse%r,&
        main_teach_sparse%sigma,main_teach_sparse%dlta,main_teach_sparse%theta,&
        main_teach_sparse%yf,main_teach_sparse%ydf,main_teach_sparse%x,main_teach_sparse%xd,&
        main_teach_sparse%xf,main_teach_sparse%xdf,main_teach_sparse%lf,main_teach_sparse%ldf,&
        main_teach_sparse%xz,main_teach_sparse%species_Z,(/(main_teach_sparse%f0,i=1,main_teach_sparse%n_species)/),&
        main_teach_sparse%target_type, gp_teach_memory_in=gp_teach_memory)
     endselect

     call system_timer('GP sparsify')

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

        call system_timer('GP initialise')
        call initialise(main_teach_sparse%my_gp,gp_sp)
        call system_timer('GP initialise')

        gp_file = 'gp_'//main_teach_sparse%m//'_'//k//'.xml'

        call teach_sparse_print_xml(main_teach_sparse,gp_file)

        call system_command('ln -fs '//trim(gp_file)//' gp.xml')

        call finalise(main_teach_sparse%my_gp)
     enddo
  else
     select case(trim(main_teach_sparse%coordinates))
     case('hf_dimer', 'water_dimer')
        if(trim(main_teach_sparse%coordinates) .eq. 'hf_dimer') then
           allocate(permutation(6,2))
           permutation(:,1) = (/1, 2, 3, 4, 5, 6/) ! original order
           permutation(:,2) = (/4, 3, 2, 1, 5, 6/) ! swap the two HF molecules
        end if
        if(trim(main_teach_sparse%coordinates) .eq. 'water_dimer') then
           allocate(permutation(15,8))
           permutation(:,1) = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15/) ! original order
           permutation(:,2) = (/1, 3, 2, 4, 5, 7, 6, 8, 9, 10, 13, 14, 11, 12, 15/) ! swap Hs on monomer A
           permutation(:,3) = (/1, 2, 3, 5, 4, 6, 7, 9, 8, 10, 12, 11, 14, 13, 15/) ! swap Hs on monomer B
           permutation(:,4) = (/1, 3, 2, 5, 4, 7, 6, 9, 8, 10, 14, 13, 12, 11, 15/) ! swap Hs on both monomers
           permutation(:,5) = (/1, 8, 9, 6, 7, 4, 5, 2, 3, 15, 11, 13, 12, 14, 10/) ! swap monomers A and B
           permutation(:,6) = (/1, 9, 8, 6, 7, 5, 4, 2, 3, 15, 12, 14, 11, 13, 10/) ! swap monomers and Hs on monomer A
           permutation(:,7) = (/1, 8, 9, 7, 6, 4, 5, 3, 2, 15, 13, 11, 14, 12, 10/) ! swap monomers and Hs on monomer B
           permutation(:,8) = (/1, 9, 8, 7, 6, 5, 4, 3, 2, 15, 14, 12, 13, 11, 10/) ! swap monomers and Hs on both monomers
        end if


        allocate(theta(size(main_teach_sparse%theta,1)))
        theta = 0.0_dp
        do i = 1, size(permutation,2)
           theta = theta + main_teach_sparse%theta(permutation(:,i),1)
        enddo
        main_teach_sparse%theta(:,1) = theta / size(permutation,2)
        deallocate(theta)

        call initialise(main_teach_sparse%my_gp, main_teach_sparse%sgm(1), &
        main_teach_sparse%dlta(1), main_teach_sparse%theta, main_teach_sparse%f0, &
        main_teach_sparse%yf, main_teach_sparse%x, permutation_in = permutation )

        !print*, main_teach_sparse%my_gp%c
        !print*, main_teach_sparse%yf
     case default
        call initialise(main_teach_sparse%my_gp, main_teach_sparse%sgm(1), &
        main_teach_sparse%dlta(1), main_teach_sparse%theta, main_teach_sparse%f0, &
        main_teach_sparse%yf, main_teach_sparse%x )
     endselect


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
     call finalise(main_teach_sparse%my_gp)
     call system_command('ln -fs '//trim(gp_file)//' gp.xml')
  endif


  call print("model parameters:")
  call print("r_cut     = "//main_teach_sparse%r_cut)
  select case(trim(main_teach_sparse%coordinates))
  case('water_monomer','water_dimer','hf_dimer')
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
  case('cosnx')
     call print("l_max     = "//main_teach_sparse%cosnx_l_max)
     call print("n_max     = "//main_teach_sparse%cosnx_n_max)
     call print("use_rdf   = "//main_teach_sparse%use_rdf)
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
