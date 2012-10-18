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
  use gp_teach_module
  use clustering_module
  use teach_sparse_mod

  implicit none

  type(teach_sparse) :: main_teach_sparse
  type(Dictionary) :: params

  character(len=STRING_LENGTH) :: verbosity
  character(len=STRING_LENGTH) :: descriptor_str, sparse_method_str, covariance_type_str, core_param_file

  logical :: has_e0

  integer :: i_coordinate
  character(len=STRING_LENGTH) :: gp_file
  logical :: has_theta_uniform

  call system_initialise(verbosity=PRINT_NORMAL)
  call enable_timing()

  call initialise(params)
  
  call param_register(params, 'at_file', PARAM_MANDATORY, main_teach_sparse%at_file, help_string="XYZ file with teaching configurations")

  call param_register(params, 'descriptor_str', PARAM_MANDATORY, descriptor_str, help_string="Initialisation string for descriptors")

  call param_register(params, 'e0', '0.0', main_teach_sparse%e0, has_value_target = has_e0, &
       help_string="Value to be subtracted from energies before fitting (and added back on after prediction)")
  
  call param_register(params, 'default_sigma', PARAM_MANDATORY, main_teach_sparse%default_sigma, help_string="error in [energies forces virials]")

  call param_register(params, 'sparse_jitter', PARAM_MANDATORY, main_teach_sparse%sparse_jitter, help_string="intrisic error of atomic/bond energy")
  
  call param_register(params, 'core_param_file', 'quip_params.xml', core_param_file, &
       help_string=" QUIP XML file for a potential to subtract from data (and added back after prediction)")
  
  call param_register(params, 'ip_args', '', main_teach_sparse%ip_args, has_value_target = main_teach_sparse%do_core, &
       help_string=" QUIP init string for a potential to subtract from data (and added back after prediction)")
  
  call param_register(params, 'energy_parameter_name', 'energy', main_teach_sparse%energy_parameter_name, &
       help_string="Name of energy property in the at_file that describe the data")
  
  call param_register(params, 'force_parameter_name', 'force', main_teach_sparse%force_parameter_name, &
       help_string="Name of force property in the at_file that describe the data")
  
  call param_register(params, 'virial_parameter_name', 'virial', main_teach_sparse%virial_parameter_name, &
       help_string="Name of virial property in the at_file that describe the data")
  
  call param_register(params, 'config_type_parameter_name', 'config_type', main_teach_sparse%config_type_parameter_name, &
       help_string="Identifier of property determining the type of input data in the at_file")
  
  call param_register(params, 'config_type_sigma','',main_teach_sparse%config_type_sigma_string, &
       has_value_target = main_teach_sparse%has_config_type_sigma, &
       help_string="What sigma values to choose for each type of data. Format: {type1:energy:force:virial:type2:energy:force:virial}")
  
  call param_register(params, 'do_sparse', 'T', main_teach_sparse%do_sparse, &
       help_string="Do sparsification or regular GP.")

  call param_register(params, 'sparseX_separate_file', 'T', main_teach_sparse%sparseX_separate_file, &
       help_string="Save sparse coordinates data in separate file")
  
  call param_register(params, 'gp_file', 'gp_new.xml', gp_file, help_string="output XML file")

  call param_register(params, 'verbosity', 'NORMAL', verbosity, help_string="Verbosity control")

  if (.not. param_read_args(params, command_line=main_teach_sparse%command_line)) then
     call print("Usage: teach_sparse at_file=file &
     descriptor_str=<string> default_sigma={<float> <float> <float>} sparse_jitter=<float> [e0]  [ip_args={}] &
     [energy_parameter_name=energy] [force_parameter_name=force] [virial_parameter_name=virial] &
     [config_type_parameter_name=config_type] &
     [config_type_sigma={type1:energy:force:virial:type2:energy:force:virial}] [do_sparse=T] [verbosity=NORMAL]")
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

  call print_title('Gaussian Approximation Potentials - Database fitting')
  call print('')
  call print('Initial parsing of command line arguments finished.')

  call split_string(descriptor_str,':','{}',main_teach_sparse%descriptor_str(:),main_teach_sparse%n_coordinate,matching=.true.)

  call print('Found '//main_teach_sparse%n_coordinate//' descriptors.')

  allocate(main_teach_sparse%delta(main_teach_sparse%n_coordinate))
  allocate(main_teach_sparse%f0(main_teach_sparse%n_coordinate))
  allocate(main_teach_sparse%n_sparseX(main_teach_sparse%n_coordinate))
  allocate(main_teach_sparse%config_type_n_sparseX_string(main_teach_sparse%n_coordinate))
  allocate(main_teach_sparse%theta_fac_string(main_teach_sparse%n_coordinate))
  allocate(main_teach_sparse%theta_file(main_teach_sparse%n_coordinate))
  allocate(main_teach_sparse%theta_uniform(main_teach_sparse%n_coordinate))
  allocate(main_teach_sparse%sparse_file(main_teach_sparse%n_coordinate))
  allocate(main_teach_sparse%mark_sparse_atoms(main_teach_sparse%n_coordinate))
  allocate(main_teach_sparse%sparse_method(main_teach_sparse%n_coordinate))
  allocate(main_teach_sparse%add_species(main_teach_sparse%n_coordinate))
  allocate(main_teach_sparse%covariance_type(main_teach_sparse%n_coordinate))
  allocate(main_teach_sparse%theta(main_teach_sparse%n_coordinate))
  allocate(main_teach_sparse%zeta(main_teach_sparse%n_coordinate))

  do i_coordinate = 1, main_teach_sparse%n_coordinate
     call initialise(params)

     call param_register(params, 'delta', '1.0', main_teach_sparse%delta(i_coordinate), help_string="Signal variance")

     call param_register(params, 'f0', '0.0', main_teach_sparse%f0(i_coordinate), help_string="Signal mean")

     call param_register(params, 'n_sparseX', '0', main_teach_sparse%n_sparseX(i_coordinate), &
     help_string="Number of sparse teaching points")

     call param_register(params, 'config_type_n_sparseX', '', main_teach_sparse%config_type_n_sparseX_string(i_coordinate), &
     help_string="Number of sparse teaching points in each config type. Format: {type1:50:type2:100}")

     call param_register(params, 'sparse_method', 'RANDOM', sparse_method_str, &
     help_string="Sparsification method. RANDOM(default), PIVOT, CLUSTER, KMEANS")

     call param_register(params, 'theta_fac', '1.0', main_teach_sparse%theta_fac_string(i_coordinate), &
     help_string="Width of Gaussians, determined from multiplying the range of each descriptor by theta_fac. &
     Can be uniform or different for each direction. For multiple theta_fac separate each value by whitespaces.")

     call param_register(params, 'theta_uniform', '0.0', main_teach_sparse%theta_uniform(i_coordinate), &
     help_string="Width of Gaussians, same in each dimension.", has_value_target = has_theta_uniform)

     call param_register(params, 'theta_file', '', main_teach_sparse%theta_file(i_coordinate), &
     help_string="Width of Gaussians from a file. There should be as many real numbers as the number of descriptors, in a single line")

     call param_register(params, 'sparse_file', '', main_teach_sparse%sparse_file(i_coordinate), &
     help_string="Sparse environments from a file. Integers, in single line")

     call param_register(params, 'mark_sparse_atoms', 'F', main_teach_sparse%mark_sparse_atoms(i_coordinate), &
     help_string="Reprints the original xyz file after sparsification process. sparse propery added, true for atoms associated with a sparse point.")

     call param_register(params, 'add_species', 'F', main_teach_sparse%add_species(i_coordinate), &
     help_string="Create species-specific descriptor, using the descriptor string as a template.")

     call param_register(params, 'covariance_type', 'ARD', covariance_type_str, &
        help_string="Type of covariance function to use. Available: ARD, DOT_PRODUCT, BOND_REAL_SPACE")

     call param_register(params, 'theta', '1.0', main_teach_sparse%theta(i_coordinate), &
     help_string="Width of Gaussians for use with real space covariance.")

     call param_register(params, 'zeta', '1.0', main_teach_sparse%zeta(i_coordinate), &
     help_string="Exponent of soap type dot product covariance kernel")

     if (.not. param_read_line(params, main_teach_sparse%descriptor_str(i_coordinate), ignore_unknown=.true., task='main program descriptor_str('//i_coordinate//')')) then
        call system_abort("main program failed to parse descriptor_str("//i_coordinate//")='"//trim(main_teach_sparse%descriptor_str(i_coordinate))//"'")
     endif
     call finalise(params)

     select case(lower_case(trim(sparse_method_str)))
     case('random')
        main_teach_sparse%sparse_method(i_coordinate) = GP_SPARSE_RANDOM
     case('pivot')
        main_teach_sparse%sparse_method(i_coordinate) = GP_SPARSE_PIVOT
     case('cluster')
        main_teach_sparse%sparse_method(i_coordinate) = GP_SPARSE_CLUSTER
     case('uniform')
        main_teach_sparse%sparse_method(i_coordinate) = GP_SPARSE_UNIFORM
     case('kmeans')
        main_teach_sparse%sparse_method(i_coordinate) = GP_SPARSE_KMEANS
     case('covariance')
        main_teach_sparse%sparse_method(i_coordinate) = GP_SPARSE_COVARIANCE
     case default
        call system_abort("unknown sparse method "//trim(sparse_method_str))
     endselect

     select case(lower_case(trim(covariance_type_str)))
     case('none')
        call system_abort("covariance type cannot be"//trim(covariance_type_str))
        main_teach_sparse%covariance_type(i_coordinate) = COVARIANCE_NONE
     case('ard')
        main_teach_sparse%covariance_type(i_coordinate) = COVARIANCE_ARD
     case('dot_product')
        main_teach_sparse%covariance_type(i_coordinate) = COVARIANCE_DOT_PRODUCT
     case('BOND_REAL_SPACE')
        main_teach_sparse%covariance_type(i_coordinate) = COVARIANCE_BOND_REAL_SPACE
     case default
        call system_abort("unknown covariance type"//trim(covariance_type_str)//". Available: ARD, DOT_PRODUCT, BOND_REAL_SPACE")
     endselect


     if(.not. has_theta_uniform) main_teach_sparse%theta_uniform(i_coordinate) = 0.0_dp
  enddo

  call print('Descriptors have been parsed')

  if(main_teach_sparse%do_core) call read(main_teach_sparse%quip_string, trim(core_param_file), keep_lf=.true.)

  call read_descriptors(main_teach_sparse) ! initialises descriptors from the descriptor_str and sets max_cutoff according to that.
  call read_teach_xyz(main_teach_sparse)   ! reads in xyz into an array of atoms objects. sets cutoff and does calc_connect on each frame
  call print('XYZ file read')

  call get_species_xyz(main_teach_sparse)  ! counts the number of species present in the xyz file.
  call add_multispecies_descriptors(main_teach_sparse)

  call parse_config_type_sigma(main_teach_sparse)
  call parse_config_type_n_sparseX(main_teach_sparse)

  if(any(main_teach_sparse%add_species)) then ! descriptor_str might have changed. reinitialises descriptors from the descriptor_str and sets max_cutoff according to that.
     call read_descriptors(main_teach_sparse)
  endif
  call print('Multispecies support added where requested')

  call teach_n_from_xyz(main_teach_sparse) ! counts number of energies, forces, virials. computes number of descriptors and gradients.

  if ( .not. has_e0) then
     call e0_avg_from_xyz(main_teach_sparse) ! calculates the average atomic energy so it can be subtracted later.
  end if

  call teach_data_from_xyz(main_teach_sparse) ! converts atomic neighbourhoods (bond neighbourhoods etc.) do descriptors, and feeds those to the GP
  call print('Cartesian coordinates transformed to descriptors')

  if(main_teach_sparse%do_sparse) then

     call system_timer('GP sparsify')

     call gp_covariance_sparse(main_teach_sparse%my_gp)
     call initialise(main_teach_sparse%gp_sp,main_teach_sparse%my_gp)
     call teach_sparse_print_xml(main_teach_sparse,gp_file,main_teach_sparse%sparseX_separate_file)

     call system_timer('GP sparsify')

  else

  endif

  call system_finalise()

end program teach_sparse_program
