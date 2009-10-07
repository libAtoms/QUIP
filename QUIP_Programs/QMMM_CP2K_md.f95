!last modified: 2008-08-30
!reads atoms & QM list & runs CP2K & writes movie xyz with QM flags
!filepot_program can be e.g. /Users/csilla/QUIP/build.darwin_x86_64_g95/cp2k_filepot

program qmmm_md

  use libatoms_module
  use cp2k_driver_module
  use quip_module

!  use atoms_module,            only: atoms, finalise, &
!                                     read_xyz, print_xyz, &
!                                     add_property, &
!                                     map_into_cell, calc_dists, cell_volume, &
!                                     set_cutoff, &
!!                                     assignment(=), &
!                                     calc_connect, diff_min_image, &
!                                     distance_min_image
!  use constraints_module,      only: print
!  use cp2k_driver_module,      only: get_qm_list, &
!                                     extend_qmlist, &
!                                     calc_topology, delete_metal_connects, &
!                                     write_psf_file, &
!                                     create_centred_qmcore, &
!                                     read_qmlist, &
!                                     quip_combine_forces, &
!                                     spline_pot, spline_force
!  use dictionary_module,       only: dictionary, initialise, finalise, &
!                                     get_value, set_value
!  use dynamicalsystem_module,  only: dynamicalsystem, &
!                                     initialise, finalise, &
!                                     add_thermostat, &
!                                     constrain_bond_diff, &
!                                     advance_verlet1, &
!                                     advance_verlet2, &
!                                     ds_print_status, &
!                                     rescale_velo
!  use linearalgebra_module,    only: norm, operator(.dot.), operator(.feq.)
!  use metapotential_module,    only: metapotential, finalise
!  use paramreader_module,      only: param_register, param_read_args, &
!                                     FIELD_LENGTH, STRING_LENGTH
!  use potential_module,        only: potential, initialise, finalise, &
!                                     calc
!  use system_module,           only: dp, inoutput, initialise, finalise, &
!                                     INPUT, OUTPUT, INOUT, &
!                                     system_timer, &
!                                     system_abort, &
!                                     system_initialise, &
!                                     system_finalise, &
!                                     system_reseed_rng, &
!                                     print, print_title, &
!                                     operator(//), round, &
!                                     NORMAL, ANAL, NERD
!  use table_module,            only: table, finalise, &
!                                     append, allocate
!  use thermostat_module,       only: NONE, LANGEVIN, NOSE_HOOVER, &
!                                     NOSE_HOOVER_LANGEVIN, &
!!                                     assignment(=), &
!                                     nose_hoover_mass, print
!  use units_module,            only: BOLTZMANN_K, PI
!

  implicit none

  integer, parameter :: TOPOLOGY_NO_PSF = 0
  integer, parameter :: TOPOLOGY_CP2K_ONLY_STEP0 = -1
  integer, parameter :: TOPOLOGY_DRIVER_ONLY_STEP0 = -2
  integer, parameter :: TOPOLOGY_USE_EXISTING_PSF = -3

  type(DynamicalSystem)               :: ds
  type(Potential)                     :: CP2K_potential
  type(MetaPotential)                 :: my_metapotential
  real(dp)                            :: energy,check,TI_force, TI_corr
  real(dp), dimension(:,:), allocatable :: f,f0,f1, add_force
  type(Table)                         :: embedlist
  type(Table)                         :: constraints
  integer                             :: i, n, N_constraints
  type(Inoutput)                      :: xyz
  type(Atoms)                         :: my_atoms
  character(len=STRING_LENGTH)                  :: args_str
  logical                             :: list_changed
  logical                             :: list_changed1
  logical                             :: empty_QM_core
  integer                             :: l
  character(len=FIELD_LENGTH)         :: PSF_Print
  real(dp)                            :: Q_QM, Q_MM, ndof_QM, ndof_MM, temp
  type(spline_pot)                    :: my_spline
  type(Table)                         :: intrares_impropers

  type(Dictionary)            :: params_in
  character(len=FIELD_LENGTH),dimension(5) :: Run_Type_array               !_MM_, QS, QMMM_EXTENDED or QMMM_CORE
  character(len=FIELD_LENGTH),dimension(4) :: Print_PSF_array               !_MM_, QS, QMMM_EXTENDED or QMMM_CORE
  character(len=FIELD_LENGTH) :: Run_Type1               !_MM_, QS, QMMM_EXTENDED or QMMM_CORE
  character(len=FIELD_LENGTH) :: Run_Type2               !_NONE_, MM, or QMMM_CORE
  integer                     :: IO_Rate                 !print coordinates at every n-th step
  integer                     :: Thermostat_Type         !_0_ none, 1 Langevin
  character(len=FIELD_LENGTH) :: Print_PSF               !_NO_PSF_, DRIVER_AT0, CP2K_AT0, EVERY_#,USE_EXIST
  integer                     :: Topology_Print          !_0_ never, -1 let CP2K print one at the 0th step and use that
                                                         !           -2 print one at the 0th time step and use that
                                                         !n>0 every n-th step
!  integer                     :: Topology_Calculate      !_0_ never
  real(dp)                    :: Time_Step
  real(dp)                    :: Equilib_Time
  real(dp)                    :: Run_Time
  real(dp)                    :: Inner_QM_Radius         !for hysteretic quantum
  real(dp)                    :: Outer_QM_Radius         !          region selection
  real(dp)                    :: Inner_Core_Radius
  real(dp)                    :: Outer_Core_Radius
  real(dp)                    :: Connect_Cutoff
  real(dp)                    :: Simulation_Temperature
  character(len=FIELD_LENGTH) :: coord_file
  character(len=FIELD_LENGTH) :: new_coord_file          !output XYZ file
  character(len=FIELD_LENGTH) :: backup_coord_file          !output XYZ file
  integer :: backup_i
  character(len=FIELD_LENGTH) :: qm_list_filename        !QM list file with a strange format
  character(len=FIELD_LENGTH) :: Residue_Library
  logical                     :: Use_Constraints         !_0_ no, 1 yes
  character(len=FIELD_LENGTH) :: constraint_file
  integer                     :: Charge
  real(dp)                    :: Tau
  logical                     :: Buffer_general, do_general
  logical                     :: Continue_it
  real(dp)                    :: avg_time
  integer                     :: Seed
  logical                     :: origin_centre
  integer                     :: atom_centre
  logical                     :: qm_zone_cartesian_centre
  real(dp)                    :: atom_centre_hysteresis, qm_zone_centre(3)
  character(len=FIELD_LENGTH) :: print_prop
  logical                     :: Delete_Metal_Connections
  real(dp)                    :: nneightol
  real(dp)                    :: spline_from
  real(dp)                    :: spline_to
  real(dp)                    :: spline_dpot
  logical                     :: use_spline
  character(len=FIELD_LENGTH) :: cp2k_calc_args               ! other args to calc(cp2k,...)
  character(len=FIELD_LENGTH) :: filepot_program
  integer :: max_n_steps
real(dp) :: pot
logical :: momentum_conservation
!pointers
integer, pointer :: qm_flag_p(:), thermostat_region_p(:)
real(dp), pointer :: pot_p(:)


    call system_initialise(verbosity=normal,enable_timing=.true.)
    call system_timer('program')

      call initialise(params_in)
      call param_register(params_in, 'Run_Type1', 'MM', Run_Type1)
      call param_register(params_in, 'Run_Type2', 'NONE', Run_Type2)
      call param_register(params_in, 'IO_Rate', '1', IO_Rate)
      call param_register(params_in, 'Thermostat_Type', '0', Thermostat_Type)
      call param_register(params_in, 'Print_PSF', 'NO_PSF', Print_PSF)
!      call param_register(params_in, 'Topology_Calculate', '0', Topology_Calculate)
      call param_register(params_in, 'Time_Step', '0.5', Time_Step)
      call param_register(params_in, 'Equilib_Time', '0.0', Equilib_Time)
      call param_register(params_in, 'Run_Time', '0.5', Run_Time)
      call param_register(params_in, 'Inner_QM_Radius', '0.0', Inner_QM_Radius)
      call param_register(params_in, 'Outer_QM_Radius', '0.0', Outer_QM_Radius)
      call param_register(params_in, 'Inner_Core_Radius', '0.0', Inner_Core_Radius)
      call param_register(params_in, 'Outer_Core_Radius', '0.0', Outer_Core_Radius)
      call param_register(params_in, 'Connect_Cutoff', '1.2', Connect_cutoff)
      call param_register(params_in, 'Simulation_Temperature', '300.0', Simulation_Temperature)
      call param_register(params_in, 'coord_file', 'coord.xyz',coord_file) 
      call param_register(params_in, 'new_coord_file', 'movie.xyz',new_coord_file) 
      call param_register(params_in, 'qm_list_filename', 'qmlist.dat', qm_list_filename)
      call param_register(params_in, 'Residue_Library', 'all_res.CHARMM.lib',Residue_Library) 
      call param_register(params_in, 'Use_Constraints', 'F', Use_Constraints)
      call param_register(params_in, 'constraint_file', 'constraints.dat',constraint_file) 
      call param_register(params_in, 'Charge', '0', Charge)
      call param_register(params_in, 'Tau', '500.0', Tau)
      call param_register(params_in, 'Buffer_general', 'F', Buffer_general)
      call param_register(params_in, 'Continue', 'F', Continue_it)
      call param_register(params_in, 'avg_time', '100.0', avg_time)
      call param_register(params_in, 'Seed', '-1', Seed)
      call param_register(params_in, 'origin_centre', 'F', origin_centre)
      call param_register(params_in, 'atom_centre', '0', atom_centre)
      call param_register(params_in, 'atom_centre_hysteresis', '0.75', atom_centre_hysteresis)
      call param_register(params_in, 'print_prop', 'all', print_prop)
      call param_register(params_in, 'nneightol', '1.2', nneightol)
      call param_register(params_in, 'Delete_Metal_Connections', 'T', Delete_Metal_Connections)
      call param_register(params_in, 'spline_from', '0.0', spline_from)
      call param_register(params_in, 'spline_to', '0.0', spline_to)
      call param_register(params_in, 'spline_dpot', '0.0', spline_dpot)
      call param_register(params_in, 'use_spline', 'F', use_spline)
      call param_register(params_in, 'max_n_steps', '-1', max_n_steps)
      cp2k_calc_args=''
      call param_register(params_in, 'cp2k_calc_args', '', cp2k_calc_args)
      call param_register(params_in, 'filepot_program', param_mandatory, filepot_program)
      call param_register(params_in, 'momentum_conservation', 'T', momentum_conservation)

      if (.not. param_read_args(params_in, do_check = .true.)) then
        call system_abort('could not parse argument line')
      end if

      if (origin_centre .and. atom_centre /= 0) &
        call system_abort("parsed both origin_centre="//origin_centre//" and non-zero atom_centre="//atom_centre)
      qm_zone_cartesian_centre = origin_centre .or. (atom_centre /= 0)

      if (Seed.gt.0) call system_reseed_rng(Seed)
!      call hello_world(seed, common_seed)

!check different run types
      Run_Type_array(1) ='QS'
      Run_Type_array(2) ='QMMM_EXTENDED'
      Run_Type_array(3) ='QMMM_CORE'
      Run_Type_array(4) ='MM'
      Run_Type_array(5) ='NONE'
      if (.not.any(Run_Type1.eq.Run_Type_array(1:4))) &
         call system_abort('Run_Type1 must be one of "QS", "MM", "QMMM_CORE", "QMMM_EXTENDED"')
      if (.not.any(Run_Type2.eq.Run_Type_array(3:5))) &
         call system_abort('Run_Type1 must be one of "NONE", "MM", "QMMM_CORE"')
      if ( (trim(Run_Type1).eq.trim(Run_Type2) .or. &
            any(trim(Run_Type1).eq.(/'MM','QS'/))) .and. &
          trim(Run_Type2).ne.'NONE' ) then
         Run_Type2 = 'NONE'
         call print('RunType2 set to NONE')
      endif
      if ((trim(Run_Type1)).eq.'QMMM_EXTENDED' .and..not.any(trim(Run_Type2).eq.Run_Type_array(3:5))) call system_abort('Run_Type1 must be higher level of accuracy than Run_Type2')
      if ((trim(Run_Type1)).eq.'QMMM_CORE' .and..not.any(trim(Run_Type2).eq.Run_Type_array(4:5))) call system_abort('Run_Type1 must be higher level of accuracy than Run_Type2')

!check PSF printing
      Print_PSF_array(1) = 'NO_PSF'
      Print_PSF_array(2) = 'DRIVER_AT_0'
      Print_PSF_array(3) = 'CP2K_AT_0'
      Print_PSF_array(4) = 'USE_EXIST'
      if (.not.any(trim(Print_PSF).eq.Print_PSF_array(1:4))) then
         if (.not.Print_PSF(1:6).eq.'EVERY_') call system_abort('Print_PSF must be one of "NO_PSF","DRIVER_AT_0","CP2K_AT_0","USE_EXIST","EVERY_#')
      endif
      if (Print_PSF.eq.'NO_PSF') Topology_Print=TOPOLOGY_NO_PSF !0
      if (Print_PSF.eq.'DRIVER_AT_0') Topology_Print=TOPOLOGY_DRIVER_ONLY_STEP0 !-2
      if (Print_PSF.eq.'CP2K_AT_0') Topology_Print=TOPOLOGY_CP2K_ONLY_STEP0 !-1
      if (Print_PSF.eq.'USE_EXIST') Topology_Print=TOPOLOGY_USE_EXISTING_PSF !-3
      if (Print_PSF(1:6).eq.'EVERY_') read(Print_PSF(7:len(Print_PSF)),*) Topology_Print !#>0

      call finalise(params_in)

      call print('Run parameters:')
      call print('  filepot_program '//trim(filepot_program))
      call print('  Run_Type1 '//Run_Type1)
      call print('  Run_Type2 '//Run_Type2)
      call print('  IO_Rate '//IO_Rate)
      if (Thermostat_Type.eq.1) then
         call print('  Thermostat_Type '//'Langevin')
         call print('  Tau '//Tau)
      else
         if (Thermostat_Type.eq.2) then
            call print('  Thermostat_Type '//'Steve')
         endif
      endif
      call print('  Print_PSF '//Print_PSF)
      call print('  nneightol '//nneightol)
      call print('  Time_Step '//round(Time_Step,3))
      call print('  Equilib_Time '//round(Equilib_Time,3))
      call print('  Run_Time '//round(Run_Time,3))
      call print('  Inner_QM_Radius '//round(Inner_QM_Radius,3))
      call print('  Outer_QM_Radius '//round(Outer_QM_Radius,3))
      call print('! - not used any more -  Connect_Cutoff '//round(Connect_cutoff,3))
      call print('  Simulation_Temperature '//round(Simulation_Temperature,3))
      call print('  coord_file '//coord_file) 
      call print('  new_coord_file '//new_coord_file) 
      if (qm_zone_cartesian_centre) then
         if (origin_centre) then
           call print('  QM core is centred around origin')
         else
           call print('  QM core is centred around atom ' // atom_centre)
           call print('    hysteresis ' // atom_centre_hysteresis)
         endif
         call print('  Inner_Core_Radius '//round(Inner_Core_Radius,3))
         call print('  Outer_Core_Radius '//round(Outer_Core_Radius,3))
         call print('  use_spline '//use_spline)
         if (use_spline) then
            call print('  spline_from '//spline_from)
            call print('  spline_to '//spline_to)
            call print('  spline_dpot '//spline_dpot)
           ! initialise spline
            my_spline%from = spline_from
            my_spline%to = spline_to
            my_spline%dpot = spline_dpot
         endif
      else
         call print('  qm_list_filename '//qm_list_filename)
      endif
      call print('  Residue_Library '//Residue_Library) 
      call print('  Use_Constraints '//Use_Constraints)
      call print('  constraint_file '//constraint_file) 
      call print('  Charge '//Charge)
      call print('  Buffer_general '//Buffer_general)
      call print('  Continue '//Continue_it)
      call print('  avg_time '//avg_time)
      call print('  Seed '//Seed)
      call print('  Properties to print '//trim(print_prop))
      call print('---------------------------------------')
      call print('')

! starts here

    if (is_file_readable(trim(new_coord_file))) then
      call print("WARNING: new_coord_file " // trim(new_coord_file) // " exists, backing it up")
      backup_i=1
      backup_coord_file=trim(new_coord_file)//".backup_"//backup_i
      do while (is_file_readable(trim(backup_coord_file)))
	backup_i = backup_i + 1
	backup_coord_file=trim(new_coord_file)//".backup_"//backup_i
      end do
      call print("WARNING:      to backup_coord_file " // trim(backup_coord_file))
      call system("cp "//trim(new_coord_file)//" "//trim(backup_coord_file))
    endif

    call initialise(xyz,new_coord_file,action=OUTPUT)
  
    !read in the coordinates, create atoms object with connectivities
    call print('Reading in the coordinates from file '//trim(coord_file)//'...')
    call read_xyz(my_atoms,coord_file)
    !use if there are CP2K velocities in the coord file, converts CP2K units to libAtoms units
    !call velocity_conversion(my_atoms)

    if (origin_centre) then
      qm_zone_centre = 0.0_dp
      call print("origin_centred qm_zone, initial qm_zone_centre " // qm_zone_centre)
    else if (atom_centre /= 0) then
      if (atom_centre < 0) call system_abort("atom_centre="//atom_centre//" < 0")
      if (atom_centre > my_atoms%N) call system_abort("atom_centre="//atom_centre//" > atoms%N="//my_atoms%N)
      if (.not. get_value(my_atoms%params,'qm_zone_centre', qm_zone_centre)) then
        qm_zone_centre = my_atoms%pos(:,atom_centre)
        call print("atom_centred qm zone, initial qm_zone_centre from atom pos " // qm_zone_centre)
      else
        call print("atom_centred qm zone, initial qm_zone_centre from input file value " // qm_zone_centre)
      endif
    endif

    N_constraints = 0
    if (Use_Constraints) then
       call print('Reading in the constraints from file '//trim(constraint_file)//'...')
       call read_constraints_bond_diff(my_atoms,constraints,trim(constraint_file))
       N_constraints = constraints%N
    endif

    call initialise(ds,my_atoms,constraints=N_constraints)
    if (Continue_it) then
      if (get_value(my_atoms%params,'Time',ds%t)) then
	  call print('Found Time in atoms%params'//ds%t)
      endif
    endif

    ds%avg_time = avg_time
    if (Thermostat_Type.eq.1) then
       call add_thermostat(ds,type=LANGEVIN,T=Simulation_Temperature,tau=Tau)
       call print('Langevin Thermostat added')
    else
       if (Thermostat_Type.eq.2) then
         !ndof is the estimated number of atoms in the QM zone (= core + buffer)
          call print('Cell_volume: '//cell_volume(ds%atoms))
          call print('Number of atoms: '//ds%atoms%N)
          call print('Estimated volume of QM zone: '//(4._dp * ( (Inner_Core_Radius + Outer_Core_Radius + Inner_QM_Radius + Outer_QM_Radius) * 0.5_dp )**3._dp * PI / 3._dp ))
          ndof_QM = nint(4._dp * ( (Inner_Core_Radius + Outer_Core_Radius + Inner_QM_Radius + Outer_QM_Radius) * 0.5_dp )**3._dp * PI / 3._dp * real(ds%atoms%N,dp) / cell_volume(ds%atoms)) * 3
          call print('density: '//(real(ds%atoms%N,dp)/cell_volume(ds%atoms)))
          ndof_MM = 3._dp * ds%atoms%N - ndof_QM
          call print('Estimated ndof QM: '//ndof_QM)
          call print('Estimated ndof MM: '//ndof_MM)
          Q_QM = nose_hoover_mass(Ndof=ndof_QM,T=Simulation_Temperature,tau=74._dp)
          Q_MM = nose_hoover_mass(Ndof=ndof_MM,T=Simulation_Temperature,tau=74._dp)
          call add_thermostat(ds,type=NOSE_HOOVER,T=Simulation_Temperature,tau=Tau,Q=Q_QM)
          call add_thermostat(ds,type=NOSE_HOOVER_LANGEVIN,T=Simulation_Temperature,tau=Tau,Q=Q_MM)
       endif
    endif
    call finalise(my_atoms)
    call add_property(ds%atoms,'pot',0._dp) ! always do this, it's just 0 if spline isn't active - no need to change print_props
    call write_cp2k_info(ds) ! to compare temperature with CP2K

!SET CONSTRAINTS

    if (Use_Constraints) then
       do i=1,constraints%N
          call Constrain_Bond_Diff(ds,constraints%int(1,i),constraints%int(2,i),constraints%int(3,i))
          call print('Set constraint: '//constraints%int(1,i)//' -- '//constraints%int(2,i)//' -- '//constraints%int(3,i))
!          call Constrain_Bond(ds,constraints%int(1,i),constraints%int(2,i))
!          call print('Set constraint: '//constraints%int(1,i)//' -- '//constraints%int(2,i))
          check = distance_min_image(ds%atoms,constraints%int(1,i),constraints%int(2,i)) - &
                  distance_min_image(ds%atoms,constraints%int(2,i),constraints%int(3,i))
          call print('constrained bond length diff: '//round(check,10))
       enddo
    endif

!SET VELOCITIES

    if (.not.Continue_it) then
       call rescale_velo(ds,Simulation_Temperature)
    endif
    !use next lines if there's no velocity in the coord file, creates also an XYZ with CP2K unit velocities
    ! call velocity_conversion_rev(ds%atoms)
    ! call print_xyz(ds%atoms,'cp2k_coord_vel.xyz',all_properties=.true.,real_format='f17.10')
    ! call velocity_conversion(ds%atoms)

! CALC. CONNECTIONS

!    call set_cutoff(ds%atoms,Connect_Cutoff)
    call set_cutoff(ds%atoms,0._dp)
    call calc_connect(ds%atoms)
    if (Delete_Metal_Connections) call delete_metal_connects(ds%atoms)

! QM LIST + THERMOSTATTING

   !QM CORE
    if ((trim(Run_Type1).eq.'QMMM_CORE') .or. &
        (trim(Run_Type1).eq.'QMMM_EXTENDED')) then
       if (.not.Continue_it) then
          if (qm_zone_cartesian_centre) then
             call add_property(ds%atoms,'QM_flag',0)
             call map_into_cell(ds%atoms)
             call calc_dists(ds%atoms)
             if (atom_centre > 0) call update_qm_zone_centre(ds%atoms, atom_centre, atom_centre_hysteresis, qm_zone_centre)
             call create_centred_qmcore(ds%atoms,Inner_Core_Radius,Outer_Core_Radius,origin=qm_zone_centre,list_changed=list_changed1)

              if (list_changed1) then
                 call print('Core has changed')
                 ! do nothing: both core and buffer belong to the QM of QM/MM
              endif
          else
             call read_qmlist(ds%atoms,qm_list_filename)
          endif
          call print('QM_flag property added')
       endif

   !QM BUFFER + THERMOSTATTING
      ! set general/heavy atom selection before QM region selection
       call set_value(ds%atoms%params,'Buffer_general',Buffer_general)
       call print('set Buffer_general into ds%atoms%params')
       if (get_value(ds%atoms%params,'Buffer_general',do_general)) then
           call print('Found Buffer_general in atoms%params'//do_general)
           buffer_general=do_general
       else
           call print('Not found Buffer_general in atoms%params')
           buffer_general=.false.
       endif

       if (trim(Run_Type1).eq.'QMMM_EXTENDED') then
          list_changed = extend_qmlist(ds%atoms,Inner_QM_Radius,Outer_QM_Radius)
          call set_value(ds%atoms%params,'QM_list_changed',list_changed)
          if (Thermostat_Type.eq.2) then !match thermostat_region to QM_flag property
             if (.not.(assign_pointer(ds%atoms, "QM_flag", qm_flag_p))) &
                call system_abort("couldn't find QM_flag property")
             if (.not.(assign_pointer(ds%atoms, "thermostat_region", thermostat_region_p))) &
                call system_abort("couldn't find thermostat_region property")
             do i=1,ds%atoms%N
                select case(qm_flag_p(i))
                   case(0)
                     thermostat_region_p(i) = 2
                   case(1,2)
                     thermostat_region_p(i) = 1
                   case default
                     call system_abort('Unknown QM_flag '//qm_flag_p(i))
                end select
             enddo
          else
            ! all atoms are in thermostat 1 by default
          endif
       endif
       if (.not.Continue_it) then
         if (qm_zone_cartesian_centre) then
               call center_atoms(ds%atoms,p=qm_zone_centre)
               call finalise(embedlist)
            endif
         else
         ! ev
         !center the whole system around the first QM atom, otherwise CP2K has difficulties with finding
         !the centre of QM box => QM atoms will be on the edges, nothing in the middle!
           call get_qm_list(ds%atoms,1,embedlist)
           call center_atoms(ds%atoms,center_i=embedlist%int(1,1))
           call finalise(embedlist)
       endif
    endif
    call map_into_cell(ds%atoms)
    call calc_dists(ds%atoms)

!TOPOLOGY

   ! topology calculation
    if (trim(Run_Type1).ne.'QS') then
!       if (.not.Continue_it) then
          call set_value(ds%atoms%params,'Library',trim(Residue_Library))
          temp = ds%atoms%nneightol
          ds%atoms%nneightol = nneightol
          call calc_topology(ds%atoms,do_CHARMM=.true.,intrares_impropers=intrares_impropers)
          call check_topology(ds%atoms)
          !call write_psf_file(ds%atoms,psf_file='psf.CHARMM.psf',run_type=Run_Type_1,intrares_impropers=intrares_impropers)
          call write_psf_file(ds%atoms,psf_file='psf.CHARMM.psf',run_type_string=trim(Run_Type1),intrares_impropers=intrares_impropers)
          ds%atoms%nneightol = temp
!       endif
    endif
    if (trim(Run_Type1).eq.'QS') then
       call set_value(ds%atoms%params,'Charge',Charge)
    endif

!PRINTING

     !----------------------------------------------------
    call set_value(ds%atoms%params,'Time',ds%t)
    if (trim(print_prop).eq.'all') then
        call print_xyz(ds%atoms,xyz,all_properties=.true.,real_format='f17.10')
    else
        call print_xyz(ds%atoms,xyz,properties=trim(print_prop),real_format='f17.10')
    endif
     !----------------------------------------------------


    !call initialise(CP2K_potential,'wrapper=.true.')
    if ((trim(Run_Type1).eq.'QS').or.(trim(Run_Type1).eq.'MM')) then
       call initialise(CP2K_potential,'FilePot command='//trim(filepot_program)//' property_list=pos min_cutoff=0.0')
    else
       call initialise(CP2K_potential,'FilePot command='//trim(filepot_program)//' property_list=pos:QM_flag min_cutoff=0.0')
    endif
    call initialise(my_metapotential,args_str='Simple=T',pot=CP2K_potential)

    !allocate force lists
    allocate(f0(3,ds%N),f1(3,ds%N),f(3,ds%N))

!FIRST STEP - first step of velocity verlet

    call system_timer('step')

     n = 0

    !force

     PSF_Print = 'DRIVER_PRINT_AND_SAVE' !in case every #th step
     if (Topology_Print.eq.TOPOLOGY_DRIVER_ONLY_STEP0) PSF_Print = 'USE_EXISTING_PSF' !DRIVER_PRINT_AND_SAVE ! we have just printed one
     if (Topology_Print.eq.TOPOLOGY_CP2K_ONLY_STEP0) PSF_Print = 'CP2K_PRINT_AND_SAVE'
     if (Topology_Print.eq.TOPOLOGY_NO_PSF) PSF_Print='NO_PSF'
     if (Topology_Print.eq.TOPOLOGY_USE_EXISTING_PSF) PSF_Print='USE_EXISTING_PSF'
     if (Topology_Print.gt.0 .or. Topology_Print.eq.(-2)) PSF_Print = 'DRIVER_PRINT_AND_SAVE'    !generate PSF (at every n-th step) and then use it
!added now
     if (Topology_Print.eq.TOPOLOGY_DRIVER_ONLY_STEP0) PSF_Print = 'USE_EXISTING_PSF' !DRIVER_PRINT_AND_SAVE ! we have just printed one
!added now
     empty_QM_core = .false.
     if (qm_zone_cartesian_centre) then
        if (.not.(assign_pointer(ds%atoms, "QM_flag", qm_flag_p))) &
           call system_abort("couldn't find QM_flag property")
        if (.not.any(qm_flag_p(1:ds%atoms%N).eq.1)) empty_QM_core = .true.
     endif
     if ((qm_zone_cartesian_centre).and.empty_QM_core) then
        call print('Empty QM core. MM run will be performed instead of QM/MM.')
        args_str=trim(cp2k_calc_args) // ' Run_Type=MM PSF_Print='// &
        trim(PSF_Print)
     else
        args_str=trim(cp2k_calc_args) // ' Run_Type='//trim(Run_Type1)// &
          ' PSF_Print='//trim(PSF_Print)
     endif
     if (.not.((trim(Run_Type1).eq.'QS').or.(trim(Run_Type1).eq.'MM'))) then
        call print_qm_region(ds%atoms)
     endif
     call calc(my_metapotential,ds%atoms,e=energy,f=f1,args_str=trim(args_str))

     PSF_Print = 'USE_EXISTING_PSF' !every case but
     if (Topology_Print.eq.TOPOLOGY_NO_PSF) PSF_Print='NO_PSF'

    !spline force calculation, if needed
     if ((qm_zone_cartesian_centre).and.use_spline) then
        allocate(add_force(1:3,1:ds%atoms%N))
	call verbosity_push_decrement()
	  call print('Force due to added spline potential (eV/A):')
	  call print('atom     F(x)     F(y)     F(z)')
          if (.not.(assign_pointer(ds%atoms, "pot", pot_p))) &
             call system_abort("couldn't find pot property")
	  do i = 1, ds%atoms%N
	     add_force(1:3,i) = spline_force(ds%atoms,i,my_spline, pot=pot)
	     pot_p(i) = pot
	     call print('  '//i//'    '//round(add_force(1,i),5)//'  '//round(add_force(2,i),5)//'  '//round(add_force(3,i),5))
	  enddo
	call verbosity_pop()
        call print('Sum of the forces: '//sum(add_force(1,1:ds%atoms%N))//' '//sum(add_force(2,1:ds%atoms%N))//' '//sum(add_force(3,1:ds%atoms%N)))
     endif

    !second force calculation, only for QMMM
     if (trim(Run_Type2).ne.'NONE') then
        if (Topology_Print.ne.0) PSF_Print = 'USE_EXISTING_PSF'  !use existing PSF
        if ((qm_zone_cartesian_centre).and.empty_QM_core) then
	   args_str = trim(cp2k_calc_args) // ' Run_Type=MM PSF_Print='//trim(PSF_Print)
        else
	   args_str = trim(cp2k_calc_args) // ' Run_Type='//trim(Run_Type2)//' PSF_Print='//trim(PSF_Print)
        endif
        if (.not.((trim(Run_Type1).eq.'QS').or.(trim(Run_Type1).eq.'MM'))) then
           call print_qm_region(ds%atoms)
        endif
        call calc(my_metapotential,ds%atoms,e=energy,f=f0,args_str=trim(args_str))
        if ((qm_zone_cartesian_centre).and.use_spline) then
           f1 = f1 + add_force
           f0 = f0 + add_force
           deallocate(add_force)
        endif
        if ((qm_zone_cartesian_centre).and.empty_QM_core) then !no 2nd run, f = f1
           f = sum0(f1)
        else
           if (momentum_conservation) then
              call QUIP_combine_forces(f1,f0,f,ds%atoms)
           else
              call abrupt_force_mixing(f1,f0,f,ds%atoms)
              f1=f
              f = sum0(f1)
           endif
        endif
     else !no 2nd run, f = f1
        if ((qm_zone_cartesian_centre).and.use_spline) then
           f = sum0(f1+add_force)
           deallocate(add_force)
        else
           f = sum0(f1)
        endif
     endif

    !PRINT DS,CONSTRAINT
     call ds_print_status(ds, 'E',energy)
     call print(ds%thermostat)
     if (Use_Constraints) then
        call print(ds%constraint)
        do i=1,constraints%N
           TI_force = force_on_collective_variable(ds%atoms,(/f(1:3,constraints%int(1,i)),f(1:3,constraints%int(2,i)),f(1:3,constraints%int(3,i))/),constraints%int(1:3,i), TI_corr, check)
           call print('constrained bond length diff: '//round(check,10))
           call print('force on colvar '//i//' :'//round(TI_force,10)//' '//round(TI_corr,10))
        enddo
     endif

    !advance verlet1
     call advance_verlet1(ds, Time_Step, f)
!call print('after advance verlet:')
!do l=1,ds%atoms%N
!call print('hello!ds: atom '//l//': '//ds%atoms%pos(1,l)//ds%atoms%velo(2,l)//ds%atoms%acc(3,l))
!enddo
!call print('ds%1: '//ds%atoms%pos(1,1)//' '//ds%atoms%pos(2,1)//' '//ds%atoms%pos(3,1))
!call print('ds%5: '//ds%atoms%pos(1,5)//' '//ds%atoms%pos(2,5)//' '//ds%atoms%pos(3,5))
!call print('ds%6: '//ds%atoms%pos(1,6)//' '//ds%atoms%pos(2,6)//' '//ds%atoms%pos(3,6))
    !PRINT XYZ
     call set_value(ds%atoms%params,'Time',ds%t)
     if (trim(print_prop).eq.'all') then
         call print_xyz(ds%atoms,xyz,all_properties=.true.,real_format='f17.10')
     else
         call print_xyz(ds%atoms,xyz,properties=trim(print_prop),real_format='f17.10')
     endif

    call system_timer('step')

!LOOP - force calc, then VV-2, VV-1

  do while (ds%t < (Equilib_Time + Run_Time) .and. ((max_n_steps < 0) .or. (n < (max_n_steps-1))))

    call system_timer('step')
     n = n + 1

   !QM CORE + BUFFER UPDATE + THERMOSTAT REASSIGNMENT
     if (trim(Run_Type1).eq.'QMMM_EXTENDED') then
        if (qm_zone_cartesian_centre) then
           if (atom_centre > 0) call update_qm_zone_centre(ds%atoms, atom_centre, atom_centre_hysteresis, qm_zone_centre)
           call create_centred_qmcore(ds%atoms,Inner_Core_Radius,Outer_Core_Radius,origin=qm_zone_centre,list_changed=list_changed1)
           if (list_changed1) then
              call print('Core has changed')
              ! do nothing: both core and buffer belong to the QM of QM/MM
           endif
        endif
        list_changed = extend_qmlist(ds%atoms,R_inner=Inner_QM_Radius,R_outer=Outer_QM_Radius)
        call set_value(ds%atoms%params,'QM_list_changed',list_changed)
!        if (list_changed) then
!           call print('QM list changed at time '//round(ds%t,2)//'. Recalculate topology.')
!           temp = ds%atoms%nneightol
!           ds%atoms%nneightol = nneightol
!           call calc_connect(ds%atoms)
!           if (Delete_Metal_Connections) call delete_metal_connects(ds%atoms)
!           call calc_topology(ds%atoms,do_CHARMM=.true.,intrares_impropers=intrares_impropers)
!           call check_topology(ds%atoms)
!!           call write_psf_file(ds%atoms,psf_file='psf.CHARMM.psf',run_type=Run_Type_1,intrares_impropers=intrares_impropers)
!           call write_psf_file(ds%atoms,psf_file='psf.CHARMM.psf',run_type_string=Run_Type1,intrares_impropers=intrares_impropers)
!           ds%atoms%nneightol = temp
!        endif
        if (Thermostat_Type.eq.2) then !match thermostat_region to QM_flag property
           if (.not.(assign_pointer(ds%atoms, "QM_flag", qm_flag_p))) &
              call system_abort("couldn't find QM_flag property")
           if (.not.(assign_pointer(ds%atoms, "thermostat_region", thermostat_region_p))) &
              call system_abort("couldn't find thermostat_region property")
           do i=1,ds%atoms%N
              select case(qm_flag_p(i))
                 case(0)
                   thermostat_region_p(i) = 2
                 case(1,2)
                   thermostat_region_p(i) = 1
                 case default
                   call system_abort('Unknown QM_flag '//qm_flag_p(i))
              end select
           enddo
        else
          ! all is in thermostat 1 by default
        endif
     endif

   !FORCE
     if (Topology_Print.gt.0) then !every #th step
        if (mod(n,Topology_Print).eq.0) then    !recalc connectivity & generate PSF (at every n-th step) and then use it
           PSF_Print = 'DRIVER_PRINT_AND_SAVE'
!added now
     if (Topology_Print.eq.TOPOLOGY_DRIVER_ONLY_STEP0) PSF_Print = 'USE_EXISTING_PSF' !DRIVER_PRINT_AND_SAVE ! we have just printed one
!added now
           if (trim(Run_Type2).ne.'QS') then
              call print_title('Recalculate Connectivity & Topology')
              call map_into_cell(ds%atoms)
              call calc_dists(ds%atoms)
              call calc_topology(ds%atoms,do_CHARMM=.true.)
              call print('CHARMM parameters added')
           endif
        endif
     endif

     empty_QM_core = .false.
     if (qm_zone_cartesian_centre) then
        if (.not.(assign_pointer(ds%atoms, "QM_flag", qm_flag_p))) &
           call system_abort("couldn't find QM_flag property")
        if (.not.any(qm_flag_p(1:ds%atoms%N).eq.1)) empty_QM_core = .true.
     endif
     if ((qm_zone_cartesian_centre).and.empty_QM_core) then
        call print('Empty QM core. MM run will be performed instead of QM/MM.')
	args_str = trim(cp2k_calc_args) // ' Run_Type=MM PSF_Print='//trim(PSF_Print)
     else
	args_str = trim(cp2k_calc_args) // ' Run_Type='//trim(Run_Type1)//' PSF_Print='//trim(PSF_Print)
     endif
     if (.not.((trim(Run_Type1).eq.'QS').or.(trim(Run_Type1).eq.'MM'))) then
        call print_qm_region(ds%atoms)
     endif
     call calc(my_metapotential,ds%atoms,e=energy,f=f1,args_str=trim(args_str))

    !spline force calculation, if needed
     if ((qm_zone_cartesian_centre).and.use_spline) then
	call verbosity_push_decrement()
	  call print('Force due to added spline potential (eV/A):')
	  call print('atom     F(x)     F(y)     F(z)')
	  allocate(add_force(1:3,1:ds%atoms%N))
          if (.not.(assign_pointer(ds%atoms, "pot", pot_p))) &
             call system_abort("couldn't find pot property")
	  do i = 1, ds%atoms%N
	     add_force(1:3,i) = spline_force(ds%atoms,i,my_spline, pot=pot)
	     pot_p(i) = pot
	     call print('  '//i//'    '//round(add_force(1,i),5)//'  '//round(add_force(2,i),5)//'  '//round(add_force(3,i),5))
	  enddo
	call verbosity_pop()
        call print('Sum of the forces: '//sum(add_force(1,1:ds%atoms%N))//' '//sum(add_force(2,1:ds%atoms%N))//' '//sum(add_force(3,1:ds%atoms%N)))
     endif

   !FORCE 2 + FORCE MIXING
     if (trim(Run_Type2).ne.'NONE') then
        if (Topology_Print.ne.TOPOLOGY_NO_PSF) PSF_Print = 'USE_EXISTING_PSF'  !use existing PSF
        if ((qm_zone_cartesian_centre).and.empty_QM_core) then
           args_str = trim(cp2k_calc_args) // ' Run_Type=MM PSF_Print='//trim(PSF_Print)
        else
           args_str = trim(cp2k_calc_args) // ' Run_Type='//trim(Run_Type2)//' PSF_Print='//trim(PSF_Print)
        endif
        if (.not.((trim(Run_Type1).eq.'QS').or.(trim(Run_Type1).eq.'MM'))) then
	   call print_qm_region(ds%atoms)
        endif
        call calc(my_metapotential,ds%atoms,e=energy,f=f0,args_str=trim(args_str))
        if ((qm_zone_cartesian_centre).and.use_spline) then
           f1 = f1 + add_force
           f0 = f0 + add_force
           deallocate(add_force)
        endif
        if ((qm_zone_cartesian_centre).and.empty_QM_core) then !no 2nd run, f = f1
           f = sum0(f1)
        else
           if (momentum_conservation) then
              call QUIP_combine_forces(f1,f0,f,ds%atoms)
           else
              call abrupt_force_mixing(f1,f0,f,ds%atoms)
              f1=f
              f = sum0(f1)
           endif
        endif
     else !no 2nd run, f = f1
        if ((qm_zone_cartesian_centre).and.use_spline) then
           f = sum0(f1+add_force)
           deallocate(add_force)
        else
           f = sum0(f1)
        endif
     endif

    !advance verlet2
     call advance_verlet2(ds, Time_Step, f)

   !PRINT DS,THERMOSTAT,CONSTRAINT,XYZ
     if (ds%t < Equilib_Time) then
        call ds_print_status(ds, 'E',energy)
     else
        call ds_print_status(ds, 'I',energy)
     end if
     call print(ds%thermostat)

    !CONSTRAINT
     if (Use_Constraints) then
        call print(ds%constraint)
        do i=1,constraints%N
           TI_force = force_on_collective_variable(ds%atoms,(/f(1:3,constraints%int(1,i)),f(1:3,constraints%int(2,i)),f(1:3,constraints%int(3,i))/),constraints%int(1:3,i), TI_corr, check)
           call print('constrained bond length diff: '//round(check,10))
           call print('force on colvar '//i//' :'//round(TI_force,10)//' '//round(TI_corr,10))
        enddo
     endif

    !XYZ
     if (mod(n,IO_Rate)==0) then
        call set_value(ds%atoms%params,'Time',ds%t)
        if (trim(print_prop).eq.'all') then
            call print_xyz(ds%atoms,xyz,all_properties=.true.,real_format='f17.10')
        else
            call print_xyz(ds%atoms,xyz,properties=trim(print_prop),real_format='f17.10')
        endif
     end if
     
    !advance verlet1
     call advance_verlet1(ds, Time_Step, f)

     call write_cp2k_info(ds)
    call system_timer('step')

  enddo

  deallocate(f,f0,f1)

  call finalise(ds)
  call finalise(xyz)
  call finalise(my_metapotential)
  call finalise(CP2K_potential)

  call print_title('THE')
  call print('Finished. CP2K is now having a rest, since deserved it. Bye-Bye!')
  call print_title('END')

    call system_timer('program')
  call system_finalise

contains

subroutine write_cp2k_info(ds)

  type(DynamicalSystem), intent(in) :: ds

  print *,'CP2K E_COM = ',( (sum(ds%atoms%mass(1:ds%atoms%N)*ds%atoms%velo(1,1:ds%atoms%N)))**2 + &
          (sum(ds%atoms%mass(1:ds%atoms%N)*ds%atoms%velo(2,1:ds%atoms%N)))**2 + &
          (sum(ds%atoms%mass(1:ds%atoms%N)*ds%atoms%velo(3,1:ds%atoms%N)))**2 ) / &
        (sum(ds%atoms%mass(1:ds%atoms%N)))
  print *,'CP2K v_COM = ',( sum(ds%atoms%mass(1:ds%atoms%N)*ds%atoms%velo(1,1:ds%atoms%N)) ) / (sum(ds%atoms%mass(1:ds%atoms%N))), &
                   ( sum(ds%atoms%mass(1:ds%atoms%N)*ds%atoms%velo(2,1:ds%atoms%N)) ) / (sum(ds%atoms%mass(1:ds%atoms%N))), &
                   ( sum(ds%atoms%mass(1:ds%atoms%N)*ds%atoms%velo(3,1:ds%atoms%N)) ) / (sum(ds%atoms%mass(1:ds%atoms%N)))
  print *,'CP2K Temp.:', 1/(285*BOLTZMANN_K) * ((sum(ds%atoms%mass(1:ds%atoms%N)*ds%atoms%velo(1,1:ds%atoms%N)**2) + &
                                                   sum(ds%atoms%mass(1:ds%atoms%N)*ds%atoms%velo(2,1:ds%atoms%N)**2) + &
                                                   sum(ds%atoms%mass(1:ds%atoms%N)*ds%atoms%velo(3,1:ds%atoms%N)**2)) - &
        ( (sum(ds%atoms%mass(1:ds%atoms%N)*ds%atoms%velo(1,1:ds%atoms%N)))**2 + &
          (sum(ds%atoms%mass(1:ds%atoms%N)*ds%atoms%velo(2,1:ds%atoms%N)))**2 + &
          (sum(ds%atoms%mass(1:ds%atoms%N)*ds%atoms%velo(3,1:ds%atoms%N)))**2 ) / &
          (sum(ds%atoms%mass(1:ds%atoms%N))) )

end subroutine

  subroutine center_atoms(at,center_i,p)

    type(atoms), intent(inout) :: at
    integer, intent(in), optional        :: center_i
    real(dp), intent(in), optional        :: p(3)
    integer                    :: i
    real(dp)                   :: shift(3)


    if (present(center_i)) then
      if (center_i.gt.at%N.or.center_i.lt.1) &
        call system_abort('center_atoms: no atom '//center_i//'in atoms: 1 - '//at%N)
      shift = 0.0_dp - at%pos(1:3,center_i)
    else
      if (.not. present(p)) call system_abort("center_atoms called with neither center_i nor p(3)")
      shift = -p
    endif

    do i=1,at%N
       at%pos(1:3,i) = at%pos(1:3,i) + shift
    enddo

  end subroutine center_atoms

  subroutine read_constraints_bond_diff(my_atoms,constraints,constraint_file)

    type(Atoms), intent(in) :: my_atoms
    type(Table), intent(out) :: constraints
    character(len=*), intent(in) :: constraint_file
    type(InOutput) :: cons
    integer :: n, num_constraints, cons_a, cons_b,cons_c

       call initialise(cons,trim(constraint_file),action=INPUT)
       read (cons%unit,*) num_constraints
       call allocate(constraints,3,0,0,0,num_constraints)
       do n=1,num_constraints
          cons_a = 0
          cons_b = 0
          cons_c = 0
          read (cons%unit,*) cons_a,cons_b, cons_c
          call append(constraints,(/cons_a,cons_b,cons_c/))
       enddo
       if (constraints%N.ne.num_constraints) call system_abort('read_constraint: Something wrong with the constraints file')
    if (any(constraints%int(1:2,1:constraints%N).gt.my_atoms%N).or.any(constraints%int(1:2,1:constraints%N).lt.1)) &
       call system_abort("read_constraints: Constraint atom(s) is <1 or >"//my_atoms%N)
       call finalise(cons)

  end subroutine read_constraints_bond_diff

  function force_on_collective_variable(my_atoms,frc,at123,F_corr,lambda) result(F_lambda)

    type(Atoms),            intent(in)  :: my_atoms
    real(dp), dimension(9), intent(in)  :: frc
    integer,  dimension(3), intent(in)  :: at123
    real(dp), optional,     intent(out) :: F_corr         ! the metric tensor correction
    real(dp), optional,     intent(out) :: lambda         ! the metric tensor correction
    real(dp)                            :: F_lambda       ! the force on the bondlength difference

    real(dp) :: d              ! the bondlength difference
    real(dp) :: DD             ! the bondlength sum
    real(dp) :: d_12(3), d12   ! the 1>2 bond vector and its norm
    real(dp) :: d_23(3), d23   ! the 2>3 bond vector and its norm
    integer  :: a1, a2, a3     ! the 3 atoms a1-a2-a3
    real(dp) :: F_1(3), F_2(3), F_3(3)   ! the force on  the 3 atoms

     a1 = at123(1)
     a2 = at123(2)
     a3 = at123(3)

     F_1(1:3) = frc(1:3)
     F_2(1:3) = frc(4:6)
     F_3(1:3) = frc(7:9)

     d_12 = diff_min_image(my_atoms,a1,a2)
     d_23 = diff_min_image(my_atoms,a2,a3)
     d12 = distance_min_image(my_atoms,a1,a2)
     d23 = distance_min_image(my_atoms,a2,a3)

     if (abs(sqrt(dot_product(d_12,d_12))-d12).gt.0.000001_dp .or. &
         abs(sqrt(dot_product(d_23,d_23))-d23).gt.0.000001_dp) then
        call print(sqrt(dot_product(d_12,d_12))//' '//d12)
        call print(sqrt(dot_product(d_23,d_23))//' '//d23)
        call system_abort('wrong realpos')
     endif

!    ! calc. F_lambda from bondlength forces - not good
!     F_d12 = dot_product((F_1(1:3)*m_2-F_2(1:3)*m_1),(-d_12(1:3))) / ((m_1+m_2)*d12)
!     F_d23 = dot_product((F_2(1:3)*m_3-F_3(1:3)*m_2),(-d_23(1:3))) / ((m_2+m_3)*d23)
!     F_lambda = F_d12 - F_d23

     F_lambda = dot_product(F_1(1:3),-d_12(1:3)) / (2._dp*d12) - &
                dot_product(F_3(1:3),d_23(1:3)) / (2._dp*d23)

     !calc. metric tensor correction
     d = d12 - d23
     DD = d12 + d23
     if (present(F_corr)) &
        F_corr = 4 * BOLTZMANN_K * 300 * DD / (DD*DD-d*d)

     if (present(lambda)) lambda = d

  end function force_on_collective_variable

! for QM/MM and MM runs, to check water topology
  subroutine check_topology(my_atoms)

    type(Atoms), intent(in)      :: my_atoms
  
    integer                                 :: i, N
    logical                                 :: do_mm
    integer,                        pointer :: qm_flag_p(:)
    character(TABLE_STRING_LENGTH), pointer :: atom_res_name_p(:)

    do_mm = .false.

    if (.not.(assign_pointer(my_atoms, "QM_flag", qm_flag_p))) then ! MM RUN
       do_mm = .true.
    end if

    if (.not.(assign_pointer(my_atoms, "atom_res_name", atom_res_name_p))) &
       call system_abort("couldn't find atom_res_name property")

    if (do_mm) then
       do i=1, my_atoms%N
          if ( ('H3O'.eq.trim(atom_res_name_p(i))) .or. &
               ('HYD'.eq.trim(atom_res_name_p(i))) .or. &
               ('HWP'.eq.trim(atom_res_name_p(i))) ) then
             call system_abort('wrong topology calculated')
          endif
       enddo
    else
       N = 0
       do i=1,my_atoms%N
          if ( .not.(qm_flag_p(i).eq.1) .and. &
!          if ( .not.any((qm_flag_p(i).eq.(/1,2/))) .and. &
               any((/'H3O','HYD','HWP'/).eq.trim(atom_res_name_p(i)))) then
            N = N + 1
            call print('ERROR: classical or buffer atom '//i//'has atom_res_name '//trim(atom_res_name_p(i)))
          endif
       enddo
       if (N.gt.0) call system_abort('wrong topology calculated')
    endif

  end subroutine check_topology

! momentum conservation, do not care about masses
  function sum0(force) result(force0)

    real(dp), dimension(:,:), intent(in) :: force
    real(dp), allocatable, dimension(:,:) :: force0
    integer :: i
    real(dp) :: sumF(3)

    allocate(force0(size(force,1),size(force,2)))

    do i = 1, size(force,2)
       sumF(1) = sum(force(1,1:size(force,2)))
       sumF(2) = sum(force(2,1:size(force,2)))
       sumF(3) = sum(force(3,1:size(force,2)))
    enddo
    if ((sumF(1).feq.0.0_dp).and.(sumF(2).feq.0.0_dp).and.(sumF(3).feq.0.0_dp)) then
       call print('Sum of the forces are zero.')
       force0(1:3,1:size(force,2)) = force(1:3,1:size(force,2))
    endif

    call print('Sum of the forces was '//sumF(1:3))
    sumF = sumF / size(force,2)

    do i = 1, size(force,2)
       force0(1:3,i) = force(1:3,i) - sumF(1:3)
    enddo

    do i = 1, size(force0,2)
       sumF(1) = sum(force0(1,1:size(force0,2)))
       sumF(2) = sum(force0(2,1:size(force0,2)))
       sumF(3) = sum(force0(3,1:size(force0,2)))
    enddo
    call print('Sum of the forces after mom.cons.: '//sumF(1:3))

  end function sum0

  subroutine update_qm_zone_centre(at, at_i, hyster, p)
    type(Atoms), intent(inout) :: at
    integer, intent(in) :: at_i
    real(dp), intent(in) :: hyster
    real(dp), intent(inout) :: p(3)

    real(dp) :: dist

    dist = distance_min_image(at, at_i, p)
    if (dist > hyster) then
      call print("updating qm_zone_centre, prev pos was " // p // " new pos of atom " // at_i //" is " // at%pos(:,at_i))
      p = at%pos(:,at_i)
      call set_value(at%params, 'qm_zone_centre', p)
    endif
  end subroutine update_qm_zone_centre


end program qmmm_md
