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

  use teach_sparse_module
  use libatoms_module
  use bispectrum_module
  use gp_sparse_module
  use clustering_module

  implicit none

  integer, parameter :: SPARSE_LENGTH = 10000
  integer, parameter :: SPARSE_N_FIELDS = 2000
  integer, parameter :: THETA_LENGTH = 10000
  real(dp), parameter :: THETA_MIN = 0.000000001

  type(teach_sparse) :: main_teach_sparse
  type(inoutput) :: bispectrum_inout, theta_inout, sparse_inout
  type(Dictionary) :: params
  type(gp_sparse) :: gp_sp

  character(len=FIELD_LENGTH) :: qw_cutoff_string, qw_cutoff_f_string, qw_cutoff_r1_string, &
  theta_file, sparse_file, z_eff_string, bispectrum_file
  real(dp) :: mem_required, mem_total, mem_free
  logical :: has_e0, has_f0, has_theta_file, has_sparse_file, has_bispectrum_file, test_gp_gradient_result

  character(len=FIELD_LENGTH), dimension(232) :: z_eff_fields
  integer :: num_z_eff_fields
  character(len=FIELD_LENGTH), dimension(99) :: qw_cutoff_fields, qw_cutoff_f_fields, qw_cutoff_r1_fields
  character(len=SPARSE_LENGTH) :: sparse_string
  character(len=FIELD_LENGTH), dimension(:), allocatable :: sparse_string_array
  character(len=THETA_LENGTH) :: theta_string
  character(len=FIELD_LENGTH), dimension(:), allocatable :: theta_string_array
  integer :: i, j, k, l, o, dd, dt
  character(len=FIELD_LENGTH) :: gp_file

  call system_initialise(verbosity=PRINT_NORMAL)
  call initialise(params)
  call param_register(params, 'at_file', PARAM_MANDATORY, main_teach_sparse%at_file)
  call param_register(params, 'm', '50', main_teach_sparse%m)
  call param_register(params, 'r_cut', '2.75', main_teach_sparse%r_cut)
  call param_register(params, 'j_max', '4', main_teach_sparse%j_max)
  call param_register(params, 'z0', '0.0', main_teach_sparse%z0)
  call param_register(params, 'qw_so3', 'F', main_teach_sparse%do_qw_so3)
  call param_register(params, 'l_max', '6', main_teach_sparse%qw_l_max)
  call param_register(params, 'cutoff', '', qw_cutoff_string)
  call param_register(params, 'cutoff_f', '', qw_cutoff_f_string)
  call param_register(params, 'cutoff_r1', '', qw_cutoff_r1_string)
  call param_register(params, 'no_q', 'F', main_teach_sparse%qw_no_q)
  call param_register(params, 'no_w', 'F', main_teach_sparse%qw_no_w)
  call param_register(params, 'e0', '0.0', main_teach_sparse%e0, has_e0)
  call param_register(params, 'f0', '0.0', main_teach_sparse%f0, has_f0)
  call param_register(params, 'sgm', '0.1 0.1 0.1', main_teach_sparse%sgm)
  call param_register(params, 'dlt', '1.0', main_teach_sparse%dlt)
  call param_register(params, 'theta_file', '', theta_file, has_theta_file)
  call param_register(params, 'sparse_file', '', sparse_file, has_sparse_file)
  call param_register(params, 'theta_fac', '1.5', main_teach_sparse%theta_fac)
  call param_register(params, 'do_sigma', 'F', main_teach_sparse%do_sigma)
  call param_register(params, 'do_delta', 'F', main_teach_sparse%do_delta)
  call param_register(params, 'do_theta', 'F', main_teach_sparse%do_theta)
  call param_register(params, 'do_sparx', 'F', main_teach_sparse%do_sparx)
  call param_register(params, 'do_f0', 'F', main_teach_sparse%do_f0)
  call param_register(params, 'do_theta_fac', 'F', main_teach_sparse%do_theta_fac)
  call param_register(params, 'do_cluster', 'F', main_teach_sparse%do_cluster)
  call param_register(params, 'do_pivot', 'F', main_teach_sparse%do_pivot)
  call param_register(params, 'min_steps', '10', main_teach_sparse%min_steps)
  call param_register(params, 'min_save', '0', main_teach_sparse%min_save)
  call param_register(params, 'z_eff', '', z_eff_string,main_teach_sparse%do_ewald)
  call param_register(params, 'do_test_gp_gradient', 'F', main_teach_sparse%do_test_gp_gradient)
  call param_register(params, 'bispectrum_file', '', bispectrum_file, has_bispectrum_file)
  call param_register(params, 'ip_args', '', main_teach_sparse%ip_args, main_teach_sparse%do_core)
  call param_register(params, 'do_ewald_corr', 'T', main_teach_sparse%do_ewald_corr)
  call param_register(params, 'energy_property_name', 'energy', main_teach_sparse%energy_property_name)
  call param_register(params, 'force_property_name', 'force', main_teach_sparse%force_property_name)
  call param_register(params, 'virial_property_name', 'virial', main_teach_sparse%virial_property_name)

  if (.not. param_read_args(params, do_check = .true., command_line=main_teach_sparse%command_line)) then
     call print("Usage: teach_sparse [at_file=file] [m=50] &
     & [r_cut=2.75] [j_max=4] [z0=0.0] [qw_so3] [l_max=6] [cutoff={:}] [cutoff_f={:}] [cutoff_r1={:}] [no_q] [no_w] &
     & [e0=0.0] [f0=avg] [sgm={0.1 0.1 0.1}] [dlt=1.0] [theta_file=file] [sparse_file=file] [theta_fac=3.0] &
     & [do_sigma=F] [do_delta=F] [do_theta=F] [do_sparx=F] [do_f0=F] [do_theta_fac=F] &
     & [do_cluster=F] [do_pivot=F] [min_steps=10] [min_save=0] [z_eff={Ga:1.0:N:-1.0}] &
     & [do_test_gp_gradient=F] [bispectrum_file=file] [ip_args={}] [do_ewald_corr=F] &
     & [energy_property_name=energy] [force_property_name=force] [virial_property_name=virial]")
     call system_abort('Exit: Mandatory argument(s) missing...')
  endif
  call finalise(params)

  if( count( (/has_e0,has_f0/) ) > 1 ) &
  & call print('Warning - you have specified both e0 and f0 - careful!')

  if( count( (/has_sparse_file,main_teach_sparse%do_cluster,main_teach_sparse%do_pivot/) ) > 1 ) &
  & call system_abort('There has been more than one method specified for sparsification.')
     
  main_teach_sparse%z_eff = 0.0_dp
  if(main_teach_sparse%do_ewald) then
     call parse_string(z_eff_string,':',z_eff_fields,num_z_eff_fields)
     do i = 1, num_z_eff_fields, 2
        j = atomic_number_from_symbol(z_eff_fields(i))
        if(j < 1 .or. j > 116) call system_abort("Invalid atomic number "//j//" parsed from "//z_eff_fields(i))
        main_teach_sparse%z_eff(j) = string_to_real(z_eff_fields(i+1))
     enddo
  endif
  main_teach_sparse%do_ewald_corr = main_teach_sparse%do_ewald .and. main_teach_sparse%do_ewald_corr

  if (main_teach_sparse%do_qw_so3) then
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
  else
     main_teach_sparse%z0 = max(1.0_dp,main_teach_sparse%z0) * main_teach_sparse%r_cut/(PI-0.02_dp)
     main_teach_sparse%d = j_max2d(main_teach_sparse%j_max)
  endif

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
     do i = 1, main_teach_sparse%nn/3
        write(bispectrum_inout%unit,"("//main_teach_sparse%d//"f16.8)") main_teach_sparse%x(:,i)
     enddo
     call finalise(bispectrum_inout)
  endif

  allocate(main_teach_sparse%dlta(main_teach_sparse%n_species))

  main_teach_sparse%dlta = main_teach_sparse%dlt 

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
  mem_required = 2.0_dp * real(size(main_teach_sparse%r),dp) * (real(size(main_teach_sparse%xf),dp) &
  + real(size(main_teach_sparse%xdf),dp)) * real(dp,dp) / (1024.0_dp**3)
  call mem_info(mem_total,mem_free)
  mem_total = mem_total / (1024.0_dp**3)
  mem_free = mem_free / (1024.0_dp**3)

  call print('Memory required (approx.): '//mem_required//' GB')
  if( mem_required > mem_free ) call system_abort('Required memory ('//mem_required//' GB) exceeds available memory ('//mem_free//' GB).')

  call gp_sparsify(gp_sp,main_teach_sparse%r,&
  main_teach_sparse%sgm,main_teach_sparse%dlta,main_teach_sparse%theta,&
  main_teach_sparse%yf,main_teach_sparse%ydf,main_teach_sparse%x,main_teach_sparse%xd,&
  main_teach_sparse%xf,main_teach_sparse%xdf,main_teach_sparse%lf,main_teach_sparse%ldf,&
  main_teach_sparse%xz,main_teach_sparse%species_Z,(/(main_teach_sparse%f0,i=1,main_teach_sparse%n_species)/),&
  main_teach_sparse%target_type)

  !deallocate(x,xd,xf,xdf,yf,ydf,lf,ldf)

  call print('')
  call print('theta')
  do l = 1, size(gp_sp%theta, 2)
     do o = 1, size(gp_sp%theta, 1)
        call print(real(gp_sp%theta(o,l),kind=dp))
     enddo
  enddo
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
     call print('delta')
     call print(real(gp_sp%delta,kind=dp))
     call print('theta')
     do l = 1, size(gp_sp%theta, 2)
        do o = 1, size(gp_sp%theta, 1)
           call print(real(gp_sp%theta(o,l),kind=dp))
        enddo
     enddo
     call print('f0')
     call print(real(gp_sp%f0,kind=dp))
     call print('')

     call initialise(main_teach_sparse%my_gp,gp_sp)

     gp_file = 'gp_'//main_teach_sparse%m//'_'//k//'.xml'

     call teach_sparse_print_xml(main_teach_sparse,gp_file)

     call system_command('ln -fs '//trim(gp_file)//' gp.xml')

     call finalise(main_teach_sparse%my_gp)
  enddo

  call print("model parameters:")
  call print("r_cut     = "//main_teach_sparse%r_cut)
  if (main_teach_sparse%do_qw_so3) then
     call print("l_max     = "//main_teach_sparse%qw_l_max)
     call print("cutoff    = "//qw_cutoff_string)
     call print("cutoff_f  = "//qw_cutoff_f_string)
     call print("cutoff_r1 = "//qw_cutoff_r1_string)
     call print("q         = "//(.not. main_teach_sparse%qw_no_q))
     call print("w         = "//(.not. main_teach_sparse%qw_no_w))
  else
     call print("j_max     = "//main_teach_sparse%j_max)
     call print("z0        = "//main_teach_sparse%z0)
  endif
  call print("n_species = "//main_teach_sparse%n_species)
  call print("species_Z = "//main_teach_sparse%species_Z)
  call print("w         = "//main_teach_sparse%w_Z(main_teach_sparse%species_Z))
  call print("z_eff     = "//main_teach_sparse%z_eff(main_teach_sparse%species_Z))
  call print("do_ewald  = "//main_teach_sparse%do_ewald)
  call print("do_ewald_corr  = "//main_teach_sparse%do_ewald_corr)
  call print("e0        = "//main_teach_sparse%e0)

  call finalise(gp_sp)

  call system_finalise()

end program teach_sparse_program
