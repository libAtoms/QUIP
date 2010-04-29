! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   Portions of this code were written by
! H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! H0 X
! H0 X   Copyright 2006-2010.
! H0 X
! H0 X   These portions of the source code are released under the GNU General
! H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
! H0 X
! H0 X   If you would like to license the source code under different terms,
! H0 X   please contact Gabor Csanyi, gabor@csanyi.net
! H0 X
! H0 X   Portions of this code were written by Noam Bernstein as part of
! H0 X   his employment for the U.S. Government, and are not subject
! H0 X   to copyright in the USA.
! H0 X
! H0 X
! H0 X   When using this software, please cite the following reference:
! H0 X
! H0 X   http://www.libatoms.org
! H0 X
! H0 X  Additional contributions by
! H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! H0 X
! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!annealing program. reprints psf file at each step. if you do not want this, change DRIVER_PRINT_AND_SAVE to USE_EXISTING_PSF in the args_strings.
!from qmmm_md. Uses CP2K filepot.
!annealing_time_step: when to step with the temperature
!max_temp: when reached, keep the temperature, do not increase it
!with_restraint: add harmonic restraint in space on the heavy atoms
!harmonic: the strength of the harmonic potential E = harmonic * (_x_ - _x_ref_)
!Temp_step: temperature step

program qmmm_md

  use libatoms_module
  use cp2k_driver_module
  use quip_module

  implicit none

  integer, parameter :: TOPOLOGY_NO_PSF = 0
  integer, parameter :: TOPOLOGY_CP2K_ONLY_STEP0 = -1
  integer, parameter :: TOPOLOGY_DRIVER_ONLY_STEP0 = -2
  integer, parameter :: TOPOLOGY_USE_EXISTING_PSF = -3

  type(DynamicalSystem)               :: ds
  type(Potential)                     :: CP2K_potential
  type(MetaPotential)                 :: my_metapotential
  real(dp)                            :: energy
  real(dp), dimension(:,:), allocatable :: f,f0,f1
  integer                             :: N_constraints = 0
  type(Inoutput)                      :: xyz
  type(Atoms)                         :: my_atoms, reference
  character(len=STRING_LENGTH)        :: args_str
  integer                             :: n, k
  real(dp)                            :: temp
  type(Table)                         :: intrares_impropers

  type(Dictionary)            :: params_in
  integer                     :: IO_Rate                 !print coordinates at every n-th step
  integer                     :: Thermostat_Type         !_0_ none, 1 Langevin
  integer                     :: Topology_Print          !_0_ never, -1 let CP2K print one at the 0th step and use that
                                                         !           -2 print one at the 0th time step and use that
                                                         !n>0 every n-th step
!  integer                     :: Topology_Calculate      !_0_ never
  real(dp)                    :: Time_Step
  real(dp)                    :: Equilib_Time
  real(dp)                    :: Run_Time
  real(dp)                    :: Simulation_Temperature
  character(len=FIELD_LENGTH) :: coord_file
  character(len=FIELD_LENGTH) :: new_coord_file          !output XYZ file
  character(len=FIELD_LENGTH) :: backup_coord_file          !output XYZ file
  integer :: backup_i
  character(len=FIELD_LENGTH) :: Residue_Library
  integer                     :: Charge
  real(dp)                    :: Tau
  real(dp)                    :: avg_time
  integer                     :: Seed
  character(len=FIELD_LENGTH) :: print_prop
  logical                     :: Delete_Metal_Connections
  real(dp)                    :: nneightol
  character(len=FIELD_LENGTH) :: cp2k_calc_args               ! other args to calc(cp2k,...)
  character(len=FIELD_LENGTH) :: filepot_program
  integer :: max_n_steps

real(dp) :: annealing_time_step
real(dp) :: Temp_step, max_temp
logical :: with_restraints
real(dp) :: harmonic !depth of the harmonic restraint potential
real(dp) :: QQ

    call system_initialise(verbosity=normal,enable_timing=.true.)
    call system_timer('program')

      call initialise(params_in)
      call param_register(params_in, 'IO_Rate', '1', IO_Rate)
      call param_register(params_in, 'Thermostat_Type', '0', Thermostat_Type)
!      call param_register(params_in, 'Topology_Calculate', '0', Topology_Calculate)
      call param_register(params_in, 'Time_Step', '0.5', Time_Step)
      call param_register(params_in, 'Temp_step', '0.5', Temp_step)
      call param_register(params_in, 'max_temp', '0.5', max_temp)
      call param_register(params_in, 'annealing_time_step', '0.5', annealing_time_step)
      call param_register(params_in, 'harmonic', '0.5', harmonic)
      call param_register(params_in, 'Equilib_Time', '0.0', Equilib_Time)
      call param_register(params_in, 'Run_Time', '0.5', Run_Time)
      call param_register(params_in, 'coord_file', 'coord.xyz',coord_file) 
      call param_register(params_in, 'new_coord_file', 'movie.xyz',new_coord_file) 
      call param_register(params_in, 'Residue_Library', 'all_res.CHARMM.lib',Residue_Library) 
      call param_register(params_in, 'Charge', '0', Charge)
      call param_register(params_in, 'Tau', '500.0', Tau)
      call param_register(params_in, 'avg_time', '100.0', avg_time)
      call param_register(params_in, 'Seed', '-1', Seed)
      call param_register(params_in, 'print_prop', 'all', print_prop)
      call param_register(params_in, 'nneightol', '1.2', nneightol)
      call param_register(params_in, 'Delete_Metal_Connections', 'T', Delete_Metal_Connections)
      call param_register(params_in, 'with_restraints', 'T', with_restraints)
      call param_register(params_in, 'max_n_steps', '-1', max_n_steps)
      cp2k_calc_args=''
      call param_register(params_in, 'cp2k_calc_args', '', cp2k_calc_args)
      call param_register(params_in, 'filepot_program', param_mandatory, filepot_program)

      if (.not. param_read_args(params_in, do_check = .true.)) then
        call system_abort('could not parse argument line')
      end if

      if (Seed.gt.0) call system_reseed_rng(Seed)
!      call hello_world(seed, common_seed)

      call finalise(params_in)

      call print('Run parameters:')
      call print('  filepot_program '//trim(filepot_program))
      call print('  Run_Type '//'MM')
      call print('  IO_Rate '//IO_Rate)
      if (Thermostat_Type.eq.1) then
         call print('  Thermostat_Type '//'Langevin')
         call print('  Tau '//Tau)
      else
         if (Thermostat_Type.eq.2) then
            call print('  Thermostat_Type '//'Steve')
         endif
      endif
      call print('  nneightol '//nneightol)
      call print('  Time_Step '//round(Time_Step,3))
      call print('  Temp_step '//round(Temp_step,3))
      call print('  annealing_time_step '//round(annealing_time_step,3))
      call print('  with_restraints '//with_restraints)
      call print('  harmonic '//round(harmonic,3))
      call print('  Equilib_Time '//round(Equilib_Time,3))
      call print('  Run_Time '//round(Run_Time,3))
      call print('  coord_file '//coord_file) 
      call print('  new_coord_file '//new_coord_file) 
      call print('  Residue_Library '//Residue_Library) 
      call print('  Charge '//Charge)
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
    reference = my_atoms

    call initialise(ds,my_atoms,constraints=N_constraints)

    Simulation_Temperature = Temp_step
    ds%avg_time = avg_time
    if (Thermostat_Type.eq.1) then
       QQ = nose_hoover_mass(Ndof=3._dp*ds%atoms%N-6,T=Simulation_Temperature,tau=74._dp)
       call add_thermostat(ds,type=NOSE_HOOVER_LANGEVIN,T=Simulation_Temperature,tau=Tau,Q=QQ)
!       call add_thermostat(ds,type=LANGEVIN,T=Simulation_Temperature,tau=Tau)
       call print('Langevin Thermostat added')
    endif
    call finalise(my_atoms)

!SET CONSTRAINTS

!SET VELOCITIES

       call rescale_velo(ds,Simulation_Temperature)

! CALC. CONNECTIONS

    call set_cutoff(ds%atoms,0._dp)
    call calc_connect(ds%atoms)
    if (Delete_Metal_Connections) call delete_metal_connects(ds%atoms)

! QM LIST + THERMOSTATTING

    call map_into_cell(ds%atoms)
    call calc_dists(ds%atoms)

!TOPOLOGY

   ! topology calculation
          call set_value(ds%atoms%params,'Library',trim(Residue_Library))
          temp = ds%atoms%nneightol
          ds%atoms%nneightol = nneightol
          call calc_topology(ds%atoms,do_CHARMM=.true.,intrares_impropers=intrares_impropers)
          call check_topology(ds%atoms)
          call write_psf_file(ds%atoms,psf_file='psf.CHARMM.psf',run_type_string=trim('MM'),intrares_impropers=intrares_impropers)
          ds%atoms%nneightol = temp

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
       call initialise(CP2K_potential,'FilePot command='//trim(filepot_program)//' property_list=pos min_cutoff=0.0')
    call initialise(my_metapotential,args_str='Simple=T',pot=CP2K_potential)

    !allocate force lists
    allocate(f0(3,ds%N),f1(3,ds%N),f(3,ds%N))

!FIRST STEP - first step of velocity verlet

    call system_timer('step')

     n = 0
     k = 0

    !force

        args_str=trim(cp2k_calc_args) // ' Run_Type='//trim('MM')// &
          ' PSF_Print='//trim('DRIVER_PRINT_AND_SAVE')
     call calc(my_metapotential,ds%atoms,e=energy,f=f1,args_str=trim(args_str))

     if (with_restraints) then
        call restrain_to_ref(ds%atoms,reference,harmonic,f0)
        f = sum0(f1+f0)
     else
        f = sum0(f1)
     endif

    !PRINT DS,CONSTRAINT
     call ds_print_status(ds, 'E',energy)
     call print(ds%thermostat)

    !advance verlet1
     call advance_verlet1(ds, Time_Step, f)
    !PRINT XYZ
     call set_value(ds%atoms%params,'Time',ds%t)
     if (trim(print_prop).eq.'all') then
         call print_xyz(ds%atoms,xyz,all_properties=.true.,real_format='f17.10')
     else
         call print_xyz(ds%atoms,xyz,properties=trim(print_prop),real_format='f17.10')
     endif

    call system_timer('step')

!LOOP - force calc, then VV-2, VV-1

  do while (ds%t < Run_Time) !(Equilib_Time + Run_Time) .and. ((max_n_steps < 0) .or. (n < (max_n_steps-1))))

    call system_timer('step')
     n = n + 1
     k = k + 1

    if (abs(real(k,dp)*Time_Step-annealing_time_step).lt.epsilon(0._dp)) then
       if (abs(Simulation_Temperature-max_temp).gt.epsilon(0._dp) .and. Simulation_Temperature.lt.max_temp) then
          call print('Rescaling velocities, stepping in temperature at time '//ds%t//' from '//Simulation_Temperature// &
                     ' K  to '//(Simulation_Temperature+Temp_step)//' K.')
          k = 0
          Simulation_Temperature = Simulation_Temperature + Temp_step
          call rescale_velo(ds,Simulation_Temperature)
          ds%thermostat(1)%T = Simulation_Temperature
       endif
    endif

   !FORCE
     if (Topology_Print.gt.0) then !every #th step
        if (mod(n,Topology_Print).eq.0) then    !recalc connectivity & generate PSF (at every n-th step) and then use it
              call print_title('Recalculate Connectivity & Topology')
              call map_into_cell(ds%atoms)
              call calc_dists(ds%atoms)
              call calc_topology(ds%atoms,do_CHARMM=.true.)
              call print('CHARMM parameters added')
        endif
     endif

        args_str=trim(cp2k_calc_args) // ' Run_Type='//trim('MM')// &
          ' PSF_Print='//trim('DRIVER_PRINT_AND_SAVE')
     call calc(my_metapotential,ds%atoms,e=energy,f=f1,args_str=trim(args_str))

     if (with_restraints) then
        call restrain_to_ref(ds%atoms,reference,harmonic,f0)
call print("BOB 20")
        f = sum0(f1+f0)
     else
call print("BOB 30")
        f = sum0(f1)
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

! calc restraint potential only on heavy atoms
  subroutine restrain_to_ref(at,ref,harm,ff)

    type(Atoms), intent(in) :: at,ref
    real(dp), intent(in) :: harm
    real(dp), dimension(:,:),intent(out) :: ff
    integer :: i

!    allocate(force0(size(force,1),size(force,2)))

    ff = 0.0_dp
    do i = 1, at%N
       if (at%Z(i).eq.1) cycle
       !extra potential is E(x,y,z) = harmonic * [(x-x_ref)**2.0 + (y-y_ref)**2.0 + (z-z_ref)**2.0 ]
       ff(1,i) = - harm * 2.0_dp * (distance_min_image(at,(/at%pos(1,i),0._dp,0._dp/),(/ref%pos(1,i),0._dp,0._dp/)))
       ff(2,i) = - harm * 2.0_dp * (distance_min_image(at,(/0._dp,at%pos(2,i),0._dp/),(/0._dp,ref%pos(2,i),0._dp/)))
       ff(3,i) = - harm * 2.0_dp * (distance_min_image(at,(/0._dp,0._dp,at%pos(3,i)/),(/0._dp,0._dp,ref%pos(3,i)/)))
    enddo

  end subroutine restrain_to_ref

! momentum conservation, do not care about masses
  function sum0(force) result(force0)

    real(dp), dimension(:,:), intent(in) :: force
    real(dp) :: force0(size(force,1),size(force,2))
    integer :: i
    real(dp) :: sumF(3)

    sumF = sum(force,2)
    call print('Sum of the forces is '//sumF(1:3))

    if ((sumF(1) .feq. 0.0_dp) .and.  (sumF(2) .feq. 0.0_dp) .and.  (sumF(3) .feq. 0.0_dp)) then
       call print('Sum of the forces is zero.')
       force0 = force
    else
      sumF = sumF / size(force,2)

      force0(1,:) = force(1,:) - sumF(1)
      force0(2,:) = force(2,:) - sumF(2)
      force0(3,:) = force(3,:) - sumF(3)
    endif

    sumF = sum(force0,2)
    call print('Sum of the forces after mom.cons.: '//sumF(1:3))

  end function sum0

end program qmmm_md
