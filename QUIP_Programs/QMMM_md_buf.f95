!adaptive hybrid QMMM MD program.
!reads atoms & QM list & runs CP2K & writes movie xyz with QM flags
!uses the metapotential with momentum conservation, adding the corrective force only to the QM core and buffer atoms
!filepot_program can be e.g. /Users/csilla/QUIP/build.darwin_x86_64_g95/cp2k_filepot
!TODO:
!use the buffer carving of the metapotential force mixing routines.

program qmmm_md

  use libatoms_module
  use cp2k_driver_template_module
  use cp2k_driver_module
  use quip_module

  implicit none

  integer, parameter :: TOPOLOGY_NO_PSF = 0
  integer, parameter :: TOPOLOGY_CP2K_ONLY_STEP0 = -1
  integer, parameter :: TOPOLOGY_DRIVER_ONLY_STEP0 = -2
  integer, parameter :: TOPOLOGY_USE_EXISTING_PSF = -3

  type(DynamicalSystem)               :: ds
  type(Atoms)                         :: my_atoms
  type(Table)                         :: constraints
  integer                             :: i, n, N_constraints
  character(len=FIELD_LENGTH)         :: Run_Type_array(5)               !_MM_, QS, QMMM_EXTENDED or QMMM_CORE
  integer                             :: l

  !Force calc.
  type(Potential)                     :: CP2K_MM_potential
  type(Potential)                     :: CP2K_QM_potential
  type(MetaPotential)                 :: my_metapotential
  character(len=STRING_LENGTH)        :: args_str
  character(len=STRING_LENGTH)        :: qm_args_str, mm_args_str
  real(dp)                            :: energy,check,TI_force, TI_corr
  real(dp), dimension(:,:), allocatable :: f,f0,f1, add_force

  !QM list generation
  type(Table)                         :: embedlist
  logical                             :: list_changed
  logical                             :: list_changed1
  logical                             :: empty_QM_core
  integer, pointer                    :: qm_flag_p(:)
  integer, pointer                    :: hybrid(:)
  integer, pointer                    :: hybrid_mark_p(:)
  integer, pointer                    :: cluster_mark_p(:)

  !Thermostat
  real(dp)                            :: Q_QM, Q_MM, ndof_QM, ndof_MM, temp
  real(dp)                            :: Q_QM_heavy, Q_QM_H
  integer, pointer                    :: thermostat_region_p(:)

  !Spline
  type(spline_pot)                    :: my_spline
  real(dp)                            :: pot
  real(dp), pointer                   :: pot_p(:)

  !Output XYZ
  type(Inoutput)                      :: xyz
!  type(CInoutput)                      :: xyz
  character(len=FIELD_LENGTH)         :: backup_coord_file          !output XYZ file
  integer                             :: backup_i

  !Topology
  character(len=FIELD_LENGTH)         :: Print_PSF_array(4)               !_MM_, QS, QMMM_EXTENDED or QMMM_CORE
  character(len=FIELD_LENGTH)         :: PSF_Print
  integer                             :: Topology_Print          !_0_ never, -1 let CP2K print one at the 0th step and use that
                                                         !           -2 print one at the 0th time step and use that
                                                         !n>0 every n-th step
  type(Table)                         :: intrares_impropers

  !Input parameters
  type(Dictionary)            :: params_in
  character(len=FIELD_LENGTH) :: Run_Type1               !_MM_, QS, QMMM_EXTENDED or QMMM_CORE
  character(len=FIELD_LENGTH) :: Run_Type2               !_NONE_, MM, or QMMM_CORE
  integer                     :: IO_Rate                 !print coordinates at every n-th step
  integer                     :: Thermostat_Type         !_0_ none, 1 Langevin
  character(len=FIELD_LENGTH) :: Print_PSF               !_NO_PSF_, DRIVER_AT0, CP2K_AT0, EVERY_#,USE_EXIST
  real(dp)                    :: Time_Step
  real(dp)                    :: Equilib_Time
  real(dp)                    :: Run_Time
  real(dp)                    :: Inner_Buffer_Radius         !for hysteretic quantum
  real(dp)                    :: Outer_Buffer_Radius         !          region selection
  real(dp)                    :: Inner_QM_Region_Radius
  real(dp)                    :: Outer_QM_Region_Radius
  real(dp)                    :: Connect_Cutoff
  real(dp)                    :: Simulation_Temperature
  character(len=FIELD_LENGTH) :: coord_file
  character(len=FIELD_LENGTH) :: new_coord_file          !output XYZ file
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
  character(len=FIELD_LENGTH) :: print_prop
  real(dp)                    :: nneightol
  logical                     :: Delete_Metal_Connections
  real(dp)                    :: spline_from
  real(dp)                    :: spline_to
  real(dp)                    :: spline_dpot
  logical                     :: use_spline
  integer                     :: max_n_steps
  character(len=FIELD_LENGTH) :: cp2k_calc_args               ! other args to calc(cp2k,...)
  character(len=FIELD_LENGTH) :: filepot_program
  logical                     :: do_carve_cluster
type(inoutput) :: csilla_out
real(dp) :: my_origin(3)
logical :: have_silica_potential
real(dp) :: use_cutoff

!    call system_initialise(verbosity=ANAL,enable_timing=.true.)
    call system_initialise(verbosity=NERD,enable_timing=.true.)
!    call system_initialise(verbosity=NORMAL,enable_timing=.true.)
    call system_timer('program')

    !INPUT
      call initialise(params_in)
      call param_register(params_in, 'Run_Type1', 'MM', Run_Type1)
      call param_register(params_in, 'Run_Type2', 'NONE', Run_Type2)
      call param_register(params_in, 'IO_Rate', '1', IO_Rate)
      call param_register(params_in, 'Thermostat_Type', '0', Thermostat_Type)
      call param_register(params_in, 'Print_PSF', 'NO_PSF', Print_PSF)
      call param_register(params_in, 'Time_Step', '0.5', Time_Step)
      call param_register(params_in, 'Equilib_Time', '0.0', Equilib_Time)
      call param_register(params_in, 'Run_Time', '0.5', Run_Time)
      call param_register(params_in, 'Inner_Buffer_Radius', '0.0', Inner_Buffer_Radius)
      call param_register(params_in, 'Outer_Buffer_Radius', '0.0', Outer_Buffer_Radius)
      call param_register(params_in, 'Inner_QM_Region_Radius', '0.0', Inner_QM_Region_Radius)
      call param_register(params_in, 'Outer_QM_Region_Radius', '0.0', Outer_QM_Region_Radius)
      call param_register(params_in, 'Connect_Cutoff', '0.0', Connect_cutoff)
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
      call param_register(params_in, 'carve_cluster', 'F', do_carve_cluster)
      call param_register(params_in, 'origin', '(/0.0 0.0 0.0/)', my_origin)
      call param_register(params_in, 'have_silica_potential', 'F', have_silica_potential)

      if (.not. param_read_args(params_in, do_check = .true.)) then
        call system_abort('could not parse argument line')
      end if

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

    !PRINT INPUT PARAMETERS
      call print('Run parameters:')
      call print('  filepot_program '//trim(filepot_program))
      call print('  Run_Type1 '//Run_Type1)
      call print('  Run_Type2 '//Run_Type2)
      call print('  IO_Rate '//IO_Rate)
      if (Thermostat_Type.eq.1) then
         call print('  Thermostat_Type '//'Langevin')
         call print('  Tau '//Tau)
      elseif (Thermostat_Type.eq.2) then
         call print('  Thermostat_Type: QM core & buffer atoms in the 1st thermostat')
         call print('                   classical O & H in the 3rd thermostat')
         call system_abort('no more used. Use Thermostat_Type=5 instead.')
      elseif (Thermostat_Type.eq.5) then
         call print('  Thermostat_Type: QM core & buffer heavy atoms in the 1st thermostat')
         call print('                   QM core & buffer H in the 2nd thermostat')
         call print('                   classical O & H in the 3rd thermostat')
      endif
      call print('  Print_PSF '//Print_PSF)
      call print('  nneightol '//nneightol)
      call print('  Time_Step '//round(Time_Step,3))
      call print('  Equilib_Time '//round(Equilib_Time,3))
      call print('  Run_Time '//round(Run_Time,3))
      call print('  Inner_Buffer_Radius '//round(Inner_Buffer_Radius,3))
      call print('  Outer_Buffer_Radius '//round(Outer_Buffer_Radius,3))
      call print('! - not used any more -  Connect_Cutoff '//round(Connect_cutoff,3))
      call print('  Simulation_Temperature '//round(Simulation_Temperature,3))
      call print('  coord_file '//coord_file) 
      call print('  new_coord_file '//new_coord_file) 
      if (origin_centre) then
         call print('  QM core is centred around origin')
         call print('  Inner_QM_Region_Radius '//round(Inner_QM_Region_Radius,3))
         call print('  Outer_QM_Region_Radius '//round(Outer_QM_Region_Radius,3))
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
      call print('  carve_cluster '//do_carve_cluster)
      call print('  my_origin '//my_origin)
      call print('  have_silica_potential '//have_silica_potential)
      call print('---------------------------------------')
      call print('')

! STARTS HERE

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
  
  !READ COORDINATES AND CONSTRAINTS

    call print('Reading in the coordinates from file '//trim(coord_file)//'...')
    call read_xyz(my_atoms,coord_file)
!    call read(my_atoms,coord_file)

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

  !THERMOSTAT

    ds%avg_time = avg_time
    if (Thermostat_Type.eq.1) then
       call add_thermostat(ds,type=LANGEVIN,T=Simulation_Temperature,tau=Tau)
       call print('Langevin Thermostat added')
    elseif (Thermostat_Type.eq.2) then
      !ndof is the estimated number of atoms in the QM zone (= core + buffer)
       call print('Cell_volume: '//cell_volume(ds%atoms))
       call print('Number of atoms: '//ds%atoms%N)
       call print('Estimated volume of QM zone: '//(4._dp * ( (Inner_QM_Region_Radius + Outer_QM_Region_Radius + Inner_Buffer_Radius + Outer_Buffer_Radius) * 0.5_dp )**3._dp * PI / 3._dp ))
       ndof_QM = nint(4._dp * ( (Inner_QM_Region_Radius + Outer_QM_Region_Radius + Inner_Buffer_Radius + Outer_Buffer_Radius) * 0.5_dp )**3._dp * PI / 3._dp * real(ds%atoms%N,dp) / cell_volume(ds%atoms)) * 3
       call print('density: '//(real(ds%atoms%N,dp)/cell_volume(ds%atoms)))
       ndof_MM = 3._dp * ds%atoms%N - ndof_QM
       call print('Estimated ndof QM: '//ndof_QM)
       call print('Estimated ndof MM: '//ndof_MM)
       Q_QM = nose_hoover_mass(Ndof=ndof_QM,T=Simulation_Temperature,tau=74._dp)
       Q_MM = nose_hoover_mass(Ndof=ndof_MM,T=Simulation_Temperature,tau=74._dp)
       call add_thermostat(ds,type=NOSE_HOOVER,T=Simulation_Temperature,Q=Q_QM,gamma=0._dp)
       call add_thermostat(ds,type=NOSE_HOOVER_LANGEVIN,T=Simulation_Temperature,tau=Tau,Q=Q_MM)
    elseif (Thermostat_Type.eq.5) then
      !ndof is the estimated number of atoms in the QM zone (= core + buffer)
       call print('Cell_volume: '//cell_volume(ds%atoms))
       call print('Number of atoms: '//ds%atoms%N)
       call print('Estimated volume of QM zone: '//(4._dp * ( (Inner_QM_Region_Radius + Outer_QM_Region_Radius + Inner_Buffer_Radius + Outer_Buffer_Radius) * 0.5_dp )**3._dp * PI / 3._dp ))
       ndof_QM = nint(4._dp * ( (Inner_QM_Region_Radius + Outer_QM_Region_Radius + Inner_Buffer_Radius + Outer_Buffer_Radius) * 0.5_dp )**3._dp * PI / 3._dp * real(ds%atoms%N,dp) / cell_volume(ds%atoms)) * 3
       call print('density: '//(real(ds%atoms%N,dp)/cell_volume(ds%atoms)))
       ndof_MM = 3._dp * ds%atoms%N - ndof_QM
       call print('Estimated ndof QM: '//ndof_QM)
       call print('     out of which ndof QM heavy atom: '//(ndof_QM/3.0_dp))
       call print('     out of which ndof QM H: '//(ndof_QM*2.0_dp/3.0_dp))
       call print('Estimated ndof MM: '//ndof_MM)
       Q_QM_heavy = nose_hoover_mass(Ndof=(ndof_QM/3.0_dp),T=Simulation_Temperature,tau=74._dp)
       Q_QM_H = nose_hoover_mass(Ndof=(ndof_QM*2.0_dp/3.0_dp),T=Simulation_Temperature,tau=74._dp)
       Q_MM = nose_hoover_mass(Ndof=ndof_MM,T=Simulation_Temperature,tau=74._dp)
       call add_thermostat(ds,type=NOSE_HOOVER,T=Simulation_Temperature,Q=Q_QM_heavy,gamma=0._dp)
       call add_thermostat(ds,type=NOSE_HOOVER,T=Simulation_Temperature,Q=Q_QM_H,gamma=0._dp)
       call add_thermostat(ds,type=NOSE_HOOVER_LANGEVIN,T=Simulation_Temperature,tau=Tau,Q=Q_MM)
    endif
    call finalise(my_atoms)
    call add_property(ds%atoms,'pot',0._dp) ! always do this, it's just 0 if spline isn't active - no need to change print_props

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

  !CALC. CONNECTIONS

!    call set_cutoff(ds%atoms,Connect_Cutoff)
!    call set_cutoff(ds%atoms,0._dp) !use the covalent radii to determine bonds
    use_cutoff = max(nneightol, Outer_Buffer_Radius)
    use_cutoff = max(use_cutoff, Outer_QM_Region_Radius)
    if (have_silica_potential) then
	use_cutoff = max(SILICA_2BODY_CUTOFF, Outer_Buffer_Radius)
    endif
    call set_cutoff(ds%atoms,use_cutoff)
    call calc_connect(ds%atoms)
    if (Delete_Metal_Connections) call delete_metal_connects(ds%atoms)

  !READ / CREATE QM LIST + THERMOSTATTING

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

   !QM CORE
    if ((trim(Run_Type1).eq.'QMMM_CORE') .or. &
        (trim(Run_Type1).eq.'QMMM_EXTENDED')) then
       if (.not.Continue_it) then
          if (origin_centre) then
!             call add_property(ds%atoms,'QM_flag',0)
call add_property(ds%atoms,'hybrid',HYBRID_NO_MARK)
call add_property(ds%atoms,'hybrid_mark',HYBRID_NO_MARK)
call add_property(ds%atoms,'old_hybrid_mark',HYBRID_NO_MARK)
call add_property(ds%atoms,'cluster_mark',HYBRID_NO_MARK)
call add_property(ds%atoms,'old_cluster_mark',HYBRID_NO_MARK)
call add_property(ds%atoms, 'cut_bonds', 0, n_cols=4) !MAX_CUT_BONDS)
             call map_into_cell(ds%atoms)
             call calc_dists(ds%atoms)
call print('hi')
             call create_centred_qmcore(ds%atoms,Inner_QM_Region_Radius,Outer_QM_Region_Radius,origin=my_origin,list_changed=list_changed1)
!!!!!!!!!
!call initialise(csilla_out,filename='csillaQM.xyz',ACTION=OUTPUT,append=.true.)
!call print_xyz(ds%atoms,xyzfile=csilla_out,properties="species:pos:hybrid:hybrid_mark")
!call finalise(csilla_out)
!!!!!!!!!
!stop

              if (list_changed1) then
                 call print('Core has changed')
                 ! do nothing: both core and buffer belong to the QM of QM/MM
!                call set_value(ds%atoms%params,'QM_core_changed',list_changed1)
              endif
          else
call add_property(ds%atoms,'hybrid',HYBRID_NO_MARK)
call add_property(ds%atoms,'hybrid_mark',HYBRID_NO_MARK)
call add_property(ds%atoms,'old_hybrid_mark',HYBRID_NO_MARK)
call add_property(ds%atoms,'cluster_mark',HYBRID_NO_MARK)
call add_property(ds%atoms,'old_cluster_mark',HYBRID_NO_MARK)
call add_property(ds%atoms, 'cut_bonds', 0, n_cols=4) !MAX_CUT_BONDS)
             call read_qmlist(ds%atoms,qm_list_filename)
                if (.not.(assign_pointer(ds%atoms, "hybrid_mark", hybrid_mark_p))) call system_abort('??')
                if (.not.(assign_pointer(ds%atoms, "cluster_mark", cluster_mark_p))) call system_abort('??')
                call print('hybrid_mark'//count(hybrid_mark_p.eq.1))
                cluster_mark_p = hybrid_mark_p
                call print('cluster_mark'//count(cluster_mark_p.eq.1))
          endif
          call print('hybrid, hybrid_mark and old_hybrid_mark properties added')
       endif


       if (origin_centre) then
         ! no centering needed: QM is centred around the origin
       else
         ! ev
         !center the whole system around the first QM atom, otherwise CP2K has difficulties with finding
         !the centre of QM box => QM atoms will be on the edges, nothing in the middle!
          if (.not.Continue_it) then
             call get_qm_list(ds%atoms,1,embedlist,'hybrid_mark')
             if (embedlist%N.eq.0) call system_abort('empty embedlist')
             call center_atoms(ds%atoms,embedlist%int(1,1))
             call finalise(embedlist)
          endif
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
	  call map_into_cell(ds%atoms)
	  call calc_dists(ds%atoms)
          call create_CHARMM(ds%atoms,do_CHARMM=.true.,intrares_impropers=intrares_impropers)
!          call calc_topology(ds%atoms,do_CHARMM=.true.,intrares_impropers=intrares_impropers)
          call check_topology(ds%atoms)
          call write_psf_file(ds%atoms,psf_file='quip_cp2k.psf',run_type_string=trim(Run_Type1),intrares_impropers=intrares_impropers,add_silica_23body=have_silica_potential)
          ds%atoms%nneightol = temp
!       endif
    endif

  !CHARGE

    if (trim(Run_Type1).eq.'QS') then
       call set_value(ds%atoms%params,'Charge',Charge)
    endif

  !PRINTING
     !----------------------------------------------------
    call set_value(ds%atoms%params,'Time',ds%t)
    if (trim(print_prop).eq.'all') then
        call print_xyz(ds%atoms,xyz,all_properties=.true.,real_format='f17.10')
        !call write(this=ds%atoms,cio=xyz,real_format='%17.10f')
    else
        call print_xyz(ds%atoms,xyz,properties=trim(print_prop),real_format='f17.10')
        !call write(this=ds%atoms,cio=xyz,real_format='%17.10f')
        !call write(this=ds%atoms,cio=xyz,properties=trim(print_prop),real_format='%17.10f')
    endif
     !----------------------------------------------------

  !INIT. METAPOTENTIAL

    !only QMMM_EXTENDED for the moment **************
    !if (trim(Run_Type1).ne.'QMMM_EXTENDED') call system_abort('ONLY QMMM_EXTENDED')

    !call initialise(CP2K_potential,'wrapper=.true.')
    if ((trim(Run_Type1).eq.'QS').or.(trim(Run_Type1).eq.'MM')) then
       call initialise(CP2K_QM_potential,'FilePot command='//trim(filepot_program)//' property_list=pos min_cutoff=0.0 Run_Type='//trim(Run_Type1)//' '//trim(cp2k_calc_args))
       call initialise(my_metapotential,args_str='Simple=T',pot=CP2K_QM_potential)
    elseif ((trim(Run_Type2).eq.'NONE')) then
       if  ((trim(Run_Type1).eq.'QMMM_EXTENDED')) call system_abort('not yet implemented')
!       call initialise(CP2K_QM_potential,'FilePot command='//trim(filepot_program)//' property_list=pos:cluster_mark:old_cluster_mark:cut_bonds min_cutoff=0.0 clean_up_files=F Run_Type='//trim(Run_Type1))
       call initialise(CP2K_QM_potential,'FilePot command='//trim(filepot_program)//' property_list=pos:cluster_mark:old_cluster_mark min_cutoff=0.0 clean_up_files=F Run_Type='//trim(Run_Type1))
       call initialise(my_metapotential,args_str='Simple=T',pot=CP2K_QM_potential)
    else
!       call initialise(CP2K_QM_potential,'FilePot command='//trim(filepot_program)//' property_list=pos:hybrid_mark:cut_bonds min_cutoff=0.0 clean_up_files=F Run_Type='//trim(Run_Type1))
       call initialise(CP2K_QM_potential,'FilePot command='//trim(filepot_program)//' property_list=pos:cluster_mark:old_cluster_mark:cut_bonds min_cutoff=0.0 clean_up_files=F Run_Type='//trim(Run_Type1))
!       call initialise(CP2K_QM_potential,'FilePot command='//trim(filepot_program)//' property_list=pos:hybrid_mark min_cutoff=0.0 clean_up_files=F Run_Type='//trim(Run_Type1))
!       call initialise(CP2K_MM_potential,'FilePot command='//trim(filepot_program)//' property_list=pos:hybrid_mark min_cutoff=0.0 clean_up_files=F Run_Type='//trim(Run_Type2)//' PSF_Print=USE_EXISTING_PSF') !every case but
       call initialise(CP2K_MM_potential,'FilePot command='//trim(filepot_program)//' property_list=pos:cluster_mark min_cutoff=0.0 clean_up_files=F Run_Type='//trim(Run_Type2)//' PSF_Print=USE_EXISTING_PSF') !every case but
       call initialise(my_metapotential,args_str='ForceMixing=T use_buffer_for_fitting=T add_cut_H_in_fitlist=T'// &
            ' method=conserve_momentum conserve_momentum_weight_method=mass calc_weights=T'// &
            ' min_images_only=F nneighb_only=F lotf_nneighb_only=F fit_hops=1 hysteretic_buffer=T'// &
            ' buffer_inner_radius='//Inner_Buffer_Radius// &
            ' buffer_outer_radius='//Outer_Buffer_Radius// &
            ' single_cluster=T little_clusters=F carve_cluster='//do_carve_cluster &
!next line is for playing with silica carving
!          //' even_electrons=T terminate=T cluster_same_lattice=T termination_clash_check=T' &
            //' construct_buffer_use_only_heavy_atoms='//(.not.(buffer_general)) &
            , pot=CP2K_MM_potential, pot2=CP2K_QM_potential)
    endif

    !allocate force lists
    allocate(f0(3,ds%N),f1(3,ds%N),f(3,ds%N))

!FIRST STEP - first step of velocity verlet

    call system_timer('step')

     n = 0

  !FORCE

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
     if (origin_centre) then
        if (.not.(assign_pointer(ds%atoms, "hybrid_mark", qm_flag_p))) &
           call system_abort("couldn't find hybrid_mark property")
        if (.not.any(qm_flag_p(1:ds%atoms%N).eq.1)) empty_QM_core = .true.
     endif
     if (.not.((trim(Run_Type1).eq.'QS').or.(trim(Run_Type1).eq.'MM'))) then
        call print_qm_region(ds%atoms)
     endif

     if ((trim(Run_Type1).eq.'QS').or.(trim(Run_Type1).eq.'MM')) then
        call print(trim(Run_Type1)//' run will be performed with simple metapotential.')
        args_str=trim(cp2k_calc_args) // &
          ' Run_Type='//trim(Run_Type1)// &
          ' PSF_Print='//trim(PSF_Print)
        call calc(CP2K_QM_potential,ds%atoms,f=f1,args_str=trim(args_str))
     elseif ((trim(Run_Type1).eq.'QMMM_CORE').or.(trim(Run_Type1).eq.'NONE')) then
        call print(trim(Run_Type1)//' run will be performed with simple metapotential.')
        args_str=trim(cp2k_calc_args) // &
          ' Run_Type='//trim(Run_Type1)// &
          ' PSF_Print='//trim(PSF_Print)
        call calc(CP2K_QM_potential,ds%atoms,f=f1,args_str=trim(args_str))
     elseif (origin_centre.and.empty_QM_core) then !only MM, no force mixing
        call print('Empty QM core. MM run will be performed instead of QM/MM.')
        args_str=trim(cp2k_calc_args) // &
          ' Run_Type=MM'// &
          ' PSF_Print='//trim(PSF_Print)
        call print('ARGS_STR | '//trim(args_str))
        call calc(CP2K_MM_potential,ds%atoms,f=f1,args_str=trim(args_str))
     else ! force mixing
        qm_args_str=trim(cp2k_calc_args) // &
          ' Run_Type='//trim(Run_Type1)// &
          ' PSF_Print='//trim(PSF_Print)// &
          ' single_cluster=T carve_cluster='//do_carve_cluster//' cluster_nneighb_only=F termination_clash_check=T terminate=T even_electrons=F'// &
          ' clean_up_files=F' !// &
!next 2lines are for playing with silica carving
!            ' single_cluster=T little_clusters=F carve_cluster='//do_carve_cluster// &
!            ' even_electrons=T terminate=T cluster_same_lattice=T termination_clash_check=T'
        mm_args_str=trim(cp2k_calc_args) // &
          ' Run_Type='//trim(Run_Type2)// &
          ' PSF_Print='//trim(PSF_Print)// &
          ' clean_up_files=F'
        args_str='qm_args_str={'//trim(qm_args_str)// &
!        args_str='qm_args_str={'//trim(mm_args_str)// &
          '} mm_args_str={'//trim(mm_args_str)//'}'
        call print('ARGS_STR | '//trim(args_str))
!        call create_hybrid_mark_and_weight_region1_from_QM_flag(ds%atoms)
!!!!!!!!!!
!call initialise(csilla_out,filename='csilla.xyz',ACTION=OUTPUT,append=.true.)
!call print_xyz(ds%atoms,xyzfile=csilla_out,properties="species:pos:hybrid:hybrid_mark:weight_region1")
!call finalise(csilla_out)
!!!!!!!!!!
!        call set_cutoff(ds%atoms,Inner_Buffer_Radius)
!        call calc_connect(ds%atoms)
        call calc(my_metapotential,ds%atoms,f=f1,args_str=trim(args_str))
!        call set_cutoff(ds%atoms,0._dp)
!        call calc_connect(ds%atoms)
     endif
     energy=0._dp !no energy

    !spline force calculation, if needed
     if (origin_centre.and.use_spline) then
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
        f = sumBUFFER(f1+add_force,ds%atoms,my_spline)
        deallocate(add_force)
     else
        f = sum0(f1,ds%atoms)
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

do i=1,ds%atoms%N
     call print('FFF '//f(1,i)//' '//f(2,i)//' '//f(3,i))
enddo

  !THERMOSTATTING now - hybrid_mark was updated only in calc
       if (trim(Run_Type1).eq.'QMMM_EXTENDED') then
          if (Thermostat_Type.eq.2) then !match thermostat_region to cluster_mark property
!             if (.not.(assign_pointer(ds%atoms, "hybrid_mark", qm_flag_p))) &
             if (.not.(assign_pointer(ds%atoms, "cluster_mark", qm_flag_p))) &
!                call system_abort("couldn't find hybrid_mark property")
                call system_abort("couldn't find cluster_mark property")
             if (.not.(assign_pointer(ds%atoms, "thermostat_region", thermostat_region_p))) &
                call system_abort("couldn't find thermostat_region property")
             do i=1,ds%atoms%N
                select case(qm_flag_p(i))
                   case(HYBRID_NO_MARK, HYBRID_TERM_MARK)
                     thermostat_region_p(i) = 2
                   case(HYBRID_ACTIVE_MARK, HYBRID_BUFFER_MARK)
                     thermostat_region_p(i) = 1
                   case default
!                     call system_abort('Unknown hybrid_mark '//qm_flag_p(i))
                     call system_abort('Unknown cluster_mark '//qm_flag_p(i))
                end select
             enddo
          elseif (Thermostat_Type.eq.5) then !match thermostat_region to cluster_mark property
             if (.not.(assign_pointer(ds%atoms, "cluster_mark", qm_flag_p))) &
                call system_abort("couldn't find cluster_mark property")
             if (.not.(assign_pointer(ds%atoms, "thermostat_region", thermostat_region_p))) &
                call system_abort("couldn't find thermostat_region property")
             do i=1,ds%atoms%N
                select case(qm_flag_p(i))
                   case(HYBRID_NO_MARK, HYBRID_TERM_MARK)
                     thermostat_region_p(i) = 3
                   case(HYBRID_ACTIVE_MARK, HYBRID_BUFFER_MARK)
                     if (ds%atoms%Z(i).eq.1) then !QM or buffer H
                        thermostat_region_p(i) = 2
                     else !QM or buffer heavy atom
                        thermostat_region_p(i) = 1
                     endif
                   case default
                     call system_abort('Unknown cluster_mark '//qm_flag_p(i))
                end select
             enddo
          else
            ! all atoms are in thermostat 1 by default
          endif
       endif

  !ADVANCE VERLET 1

     call advance_verlet1(ds, Time_Step, f)

  !PRINT XYZ

     call set_value(ds%atoms%params,'Time',ds%t)
     if (trim(print_prop).eq.'all') then
         call print_xyz(ds%atoms,xyz,all_properties=.true.,real_format='f17.10')
!         call write(ds%atoms,xyz,real_format='%17.10f')
     else
         call print_xyz(ds%atoms,xyz,properties=trim(print_prop),real_format='f17.10')
!         call write(ds%atoms,xyz,real_format='%17.10f')
         !call write(ds%atoms,xyz,properties=trim(print_prop),real_format='%17.10f')
     endif

    call system_timer('step')

!LOOP - force calc, then VV-2, VV-1

  do while (ds%t < (Equilib_Time + Run_Time) .and. ((max_n_steps < 0) .or. (n < (max_n_steps-1))))

    call system_timer('step')
     n = n + 1

  !QM CORE + BUFFER UPDATE + THERMOSTAT REASSIGNMENT

     if (trim(Run_Type1).eq.'QMMM_EXTENDED') then
        if (origin_centre) then
           call create_centred_qmcore(ds%atoms,Inner_QM_Region_Radius,Outer_QM_Region_Radius,origin=my_origin,list_changed=list_changed1)
           if (list_changed1) then
              call print('Core has changed')
!             call set_value(ds%atoms%params,'QM_core_changed',list_changed)
              ! do nothing: both core and buffer belong to the QM of QM/MM
           endif
        endif
     endif

  !FORCE

     if (Topology_Print.gt.0) then !every #th step
        if (mod(n,Topology_Print).eq.0) then    !recalc connectivity & generate PSF (at every n-th step) and then use it
           PSF_Print = 'DRIVER_PRINT_AND_SAVE'
           if (Topology_Print.eq.TOPOLOGY_DRIVER_ONLY_STEP0) PSF_Print = 'USE_EXISTING_PSF' !DRIVER_PRINT_AND_SAVE ! we have just printed one
           if (trim(Run_Type1).ne.'QS') then
              call print_title('Recalculate Connectivity & Topology')
              call map_into_cell(ds%atoms)
	      call calc_dists(ds%atoms)
              ! call calc_dists(ds%atoms)
              call create_CHARMM(ds%atoms,do_CHARMM=.true.,intrares_impropers=intrares_impropers)
!              call calc_topology(ds%atoms,do_CHARMM=.true.)
              call print('CHARMM parameters added')
           endif
        endif
     endif

     empty_QM_core = .false.
     if (origin_centre) then
        if (.not.(assign_pointer(ds%atoms, "hybrid_mark", qm_flag_p))) &
           call system_abort("couldn't find hybrid_mark property")
        if (.not.any(qm_flag_p(1:ds%atoms%N).eq.1)) empty_QM_core = .true.
     endif
     if (.not.((trim(Run_Type1).eq.'QS').or.(trim(Run_Type1).eq.'MM'))) then
        call print_qm_region(ds%atoms)
     endif

     if ((trim(Run_Type1).eq.'QS').or.(trim(Run_Type1).eq.'MM')) then
        call print(trim(Run_Type1)//' run will be performed with simple metapotential.')
        args_str=trim(cp2k_calc_args) // &
          ' Run_Type='//trim(Run_Type1)// &
          ' PSF_Print='//trim(PSF_Print)
        call calc(CP2K_QM_potential,ds%atoms,f=f1,args_str=trim(args_str))
     elseif ((trim(Run_Type1).eq.'QMMM_CORE').or.(trim(Run_Type1).eq.'NONE')) then
        call print(trim(Run_Type1)//' run will be performed with simple metapotential.')
        args_str=trim(cp2k_calc_args) // &
          ' Run_Type='//trim(Run_Type1)// &
          ' PSF_Print='//trim(PSF_Print)
        call calc(CP2K_QM_potential,ds%atoms,f=f1,args_str=trim(args_str))
     elseif (origin_centre.and.empty_QM_core) then !only MM, no force mixing
        call print('Empty QM core. MM run will be performed instead of QM/MM.')
        args_str=trim(cp2k_calc_args) // &
          ' Run_Type=MM'// &
          ' PSF_Print='//trim(PSF_Print)
        call print('ARGS_STR | '//trim(args_str))
        call calc(CP2K_MM_potential,ds%atoms,f=f1,args_str=trim(args_str))
     else ! force mixing
        qm_args_str=trim(cp2k_calc_args) // &
          ' Run_Type='//trim(Run_Type1)// &
          ' PSF_Print='//trim(PSF_Print)// &
          ' single_cluster=T carve_cluster='//do_carve_cluster//' cluster_nneighb_only=F termination_clash_check=T terminate=T even_electrons=F'// &
          ' clean_up_files=F' !// &
!next 2lines are for playing with silica carving
!            ' single_cluster=T little_clusters=F carve_cluster='//do_carve_cluster &
!          //' even_electrons=T terminate=T cluster_same_lattice=T termination_clash_check=T'
        mm_args_str=trim(cp2k_calc_args) // &
          ' Run_Type='//trim(Run_Type2)// &
          ' PSF_Print='//trim(PSF_Print)// &
          ' clean_up_files=F'
        args_str='qm_args_str={'//trim(qm_args_str)// &
!        args_str='qm_args_str={'//trim(mm_args_str)// &
          '} mm_args_str={'//trim(mm_args_str)//'}'
        call print('ARGS_STR | '//trim(args_str))
!        call set_cutoff(ds%atoms,Inner_Buffer_Radius)
!        call calc_connect(ds%atoms)
        call calc(my_metapotential,ds%atoms,f=f1,args_str=trim(args_str))
!        call set_cutoff(ds%atoms,0._dp)
!        call calc_connect(ds%atoms)
     endif
     energy=0._dp !no energy

    !SPLINE force calculation, if needed
     if (origin_centre.and.use_spline) then
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
        f = sumBUFFER(f1+add_force,ds%atoms,my_spline)
        deallocate(add_force)
     else
        f = sum0(f1,ds%atoms)
     endif

do i=1,ds%atoms%N
     call print('FFF '//f(1,i)//' '//f(2,i)//' '//f(3,i))
enddo

  !THERMOSTATTING now - hybrid_mark was updated only in calc
       if (trim(Run_Type1).eq.'QMMM_EXTENDED') then
          if (Thermostat_Type.eq.2) then !match thermostat_region to QM_flag property
!             if (.not.(assign_pointer(ds%atoms, "hybrid_mark", qm_flag_p))) &
             if (.not.(assign_pointer(ds%atoms, "cluster_mark", qm_flag_p))) &
!                call system_abort("couldn't find hybrid_mark property")
                call system_abort("couldn't find cluster_mark property")
             if (.not.(assign_pointer(ds%atoms, "thermostat_region", thermostat_region_p))) &
                call system_abort("couldn't find thermostat_region property")
             do i=1,ds%atoms%N
                select case(qm_flag_p(i))
                   case(HYBRID_NO_MARK, HYBRID_TERM_MARK)
                     thermostat_region_p(i) = 2
                   case(HYBRID_ACTIVE_MARK, HYBRID_BUFFER_MARK)
                     thermostat_region_p(i) = 1
                   case default
!                     call system_abort('Unknown hybrid_mark '//qm_flag_p(i))
                     call system_abort('Unknown cluster_mark '//qm_flag_p(i))
                end select
             enddo
          elseif (Thermostat_Type.eq.5) then !match thermostat_region to cluster_mark property
             if (.not.(assign_pointer(ds%atoms, "cluster_mark", qm_flag_p))) &
                call system_abort("couldn't find cluster_mark property")
             if (.not.(assign_pointer(ds%atoms, "thermostat_region", thermostat_region_p))) &
                call system_abort("couldn't find thermostat_region property")
             do i=1,ds%atoms%N
                select case(qm_flag_p(i))
                   case(HYBRID_NO_MARK, HYBRID_TERM_MARK)
                     thermostat_region_p(i) = 3
                   case(HYBRID_ACTIVE_MARK, HYBRID_BUFFER_MARK)
                     if (ds%atoms%Z(i).eq.1) then !QM or buffer H
                        thermostat_region_p(i) = 2
                     else !QM or buffer heavy atom
                        thermostat_region_p(i) = 1
                     endif
                   case default
                     call system_abort('Unknown cluster_mark '//qm_flag_p(i))
                end select
             enddo
          else
            ! all atoms are in thermostat 1 by default
          endif
       endif

  !ADVANCE VERLET 2

     call advance_verlet2(ds, Time_Step, f)

  !PRINT DS,THERMOSTAT,CONSTRAINT,XYZ

     if (ds%t < Equilib_Time) then
        call ds_print_status(ds, 'E',energy)
     else
        call ds_print_status(ds, 'I',energy)
     end if

     !Thermostat
     call print(ds%thermostat)

     !Constraint
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
!            call write(ds%atoms,xyz,real_format='%17.10f')
        else
            call print_xyz(ds%atoms,xyz,properties=trim(print_prop),real_format='f17.10')
!            call write(ds%atoms,xyz,real_format='%17.10f')
            !call write(ds%atoms,xyz,properties=trim(print_prop),real_format='%17.10f')
        endif
     end if
     
  !ADVANCE VERLET 1

     call advance_verlet1(ds, Time_Step, f)

     call system_timer('step')

  enddo

  deallocate(f,f0,f1)

  call finalise(ds)
  call finalise(xyz)
  call finalise(my_metapotential)
  call finalise(CP2K_QM_potential)
  if (.not.((trim(Run_Type1).eq.'QS').or.(trim(Run_Type1).eq.'MM'))) call finalise(CP2K_MM_potential)

  call print_title('THE')
  call print('Finished. CP2K is now having a rest, since deserved it. Bye-Bye!')
  call print_title('END')

  call system_timer('program')
  call system_finalise

contains

  subroutine center_atoms(at,this)

    type(atoms), intent(inout) :: at
    integer, intent(in)        :: this
    integer                    :: i
    real(dp)                   :: shift(3)

    if (this.gt.at%N.or.this.lt.1) call system_abort('center_atoms: no atom '//this//'in atoms: 1 - '//at%N)

    shift = 0.0_dp - at%pos(1:3,this)

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

    if (.not.(assign_pointer(my_atoms, "hybrid_mark", qm_flag_p))) then ! MM RUN
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

! momentum conservation over all atoms, mass weighted
  function sum0(force,at) result(force0)

    real(dp), dimension(:,:), intent(in) :: force
    type(Atoms),              intent(in) :: at
    real(dp) :: force0(size(force,1),size(force,2))
    integer :: i
    real(dp) :: sumF(3), sum_weight

    sumF = sum(force,2)
    call print('Sum of the forces is '//sumF(1:3))

    if ((sumF(1) .feq. 0.0_dp) .and.  (sumF(2) .feq. 0.0_dp) .and.  (sumF(3) .feq. 0.0_dp)) then
       call print('Sum of the forces is zero.')
       force0 = force
    else
      !F_corr weighted by element mass
      sum_weight = sum(ElementMass(at%Z(1:at%N)))

      force0(1,:) = force(1,:) - sumF(1) * ElementMass(at%Z(:)) / sum_weight
      force0(2,:) = force(2,:) - sumF(2) * ElementMass(at%Z(:)) / sum_weight
      force0(3,:) = force(3,:) - sumF(3) * ElementMass(at%Z(:)) / sum_weight
    endif

    sumF = sum(force0,2)
    call print('Sum of the forces after mom.cons.: '//sumF(1:3))

  end function sum0

! momentum conservation over atoms with QM flag 1 or 2, mass weighted
  function sumBUFFER(force,at,my_spline) result(force0)

    real(dp), dimension(:,:), intent(in) :: force
    type(Atoms),              intent(in) :: at
    type(spline_pot),         intent(in) :: my_spline
    real(dp) :: force0(size(force,1),size(force,2))
    real(dp) :: F_corr(size(force,1),size(force,2))
    integer :: i
    real(dp) :: sumF(3), sum_weight
    integer, pointer :: qm_flag_p(:)

!    if (.not. assign_pointer(at, 'hybrid_mark', qm_flag_p)) &
    if (.not. assign_pointer(at, 'cluster_mark', qm_flag_p)) &
!         call system_abort('MetaPotential_FM_Calc: hybrid_mark property missing')
         call system_abort('MetaPotential_FM_Calc: cluster_mark property missing')
    if ( .not. any(qm_flag_p(:).eq.(/1,2/)) ) then !no buffer or core: maybe empty_qm_core?
        force0 = sum0(force,at)
        return
    endif

    !sum F
    sumF = sum(force,2)
    call print('Sum of the forces is '//sumF(1:3))

    if ((sumF(1) .feq. 0.0_dp) .and.  (sumF(2) .feq. 0.0_dp) .and.  (sumF(3) .feq. 0.0_dp)) then
       call print('Sum of the forces is zero.')
       force0 = force
    else
      !F_corr weighted by element mass
      sum_weight = 0._dp
      F_corr = 0._dp
      do i=1, at%N
         if ( any(qm_flag_p(i).eq.(/1,2/)) ) then
            F_corr(1:3,i) = ElementMass(at%Z(i))
            sum_weight = sum_weight + ElementMass(at%Z(i))
         endif
      enddo
      if (sum_weight .feq. 0._dp) call system_abort('sum_buffer: 0 element masses? the presence of core or buffer atoms has already been checked.')

      force0(1,:) = force(1,:) - sumF(1) * F_corr(1,:) / sum_weight
      force0(2,:) = force(2,:) - sumF(2) * F_corr(2,:) / sum_weight
      force0(3,:) = force(3,:) - sumF(3) * F_corr(3,:) / sum_weight
    endif

    sumF = sum(force0,2)
    call print('Sum of the forces after mom.cons.: '//sumF(1:3))

  end function sumBUFFER

end program qmmm_md
