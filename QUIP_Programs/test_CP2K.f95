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

!last modified: 2008-08-30
!reads atoms & QM list & runs CP2K & writes movie xyz with QM flags

program test_CP2K

  use libatoms_module
  use cp2k_driver_module
!  use quantumselection_module
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
  type(Table)                         :: embedlist
  integer                             :: n
  type(Inoutput)                      :: xyz
  type(Atoms)                         :: my_atoms
  character(len=STRING_LENGTH)                  :: args_str
  logical                             :: list_changed
  integer                             :: PSF_Print

  type(Dictionary)            :: params_in
  character(len=FIELD_LENGTH),dimension(5) :: Run_Type_array               !_MM_, QS, QMMM_EXT or QMMM_CORE
  character(len=FIELD_LENGTH),dimension(4) :: Print_PSF_array               !_MM_, QS, QMMM_EXT or QMMM_CORE
  character(len=FIELD_LENGTH) :: Run_Type1               !_MM_, QS, QMMM_EXT or QMMM_CORE
  character(len=FIELD_LENGTH) :: Run_Type2               !_NONE_, MM, or QMMM_CORE
  integer                     :: Run_Type_1              !type of force calculation: MM_RUN, QS_RUN, QMMM_RUN_CORE or QMMM_RUN_EXTENDED
  integer                     :: Run_Type_2              !type of second force calculation: NONE_RUN, MM_RUN or QMMM_RUN_CORE
  integer                     :: nproc
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
  real(dp)                    :: Connect_Cutoff
  real(dp)                    :: Simulation_Temperature
  character(len=FIELD_LENGTH) :: coord_file
  character(len=FIELD_LENGTH) :: new_coord_file          !output XYZ file
  character(len=FIELD_LENGTH) :: qm_list_filename        !QM list file with a strange format
  character(len=FIELD_LENGTH) :: Residue_Library

    call system_initialise(verbosity=normal)

      call initialise(params_in)
      call param_register(params_in, 'Run_Type1', 'MM', Run_Type1)
      call param_register(params_in, 'Run_Type2', 'NONE', Run_Type2)
      call param_register(params_in, 'nproc', '1', nproc)
      call param_register(params_in, 'IO_Rate', '1', IO_Rate)
      call param_register(params_in, 'Thermostat_Type', '0', Thermostat_Type)
      call param_register(params_in, 'Print_PSF', 'NO_PSF', Print_PSF)
!      call param_register(params_in, 'Topology_Calculate', '0', Topology_Calculate)
      call param_register(params_in, 'Time_Step', '0.5', Time_Step)
      call param_register(params_in, 'Equilib_Time', '0.0', Equilib_Time)
      call param_register(params_in, 'Run_Time', '0.5', Run_Time)
      call param_register(params_in, 'Inner_QM_Radius', '0.0', Inner_QM_Radius)
      call param_register(params_in, 'Outer_QM_Radius', '0.0', Outer_QM_Radius)
      call param_register(params_in, 'Connect_Cutoff', '1.2', Connect_cutoff)
      call param_register(params_in, 'Simulation_Temperature', '300.0', Simulation_Temperature)
      call param_register(params_in, 'coord_file', 'coord.xyz',coord_file) 
      call param_register(params_in, 'new_coord_file', 'movie.xyz',new_coord_file) 
      call param_register(params_in, 'qm_list_filename', 'qmlist.dat', qm_list_filename)
      call param_register(params_in, 'Residue_Library', 'all_res.CHARMM.lib',Residue_Library) 

      if (.not. param_read_args(params_in, do_check = .true.)) then
        call system_abort('could not parse argument line')
      end if

!check different run types
      Run_Type_array(1) ='QS'
      Run_Type_array(2) ='QMMM_EXT'
      Run_Type_array(3) ='QMMM_CORE'
      Run_Type_array(4) ='MM'
      Run_Type_array(5) ='NONE'
      if (.not.any(Run_Type1.eq.Run_Type_array(1:4))) &
         call system_abort('Run_Type1 must be one of "QS", "MM", "QMMM_CORE", "QMMM_EXT"')
      if (.not.any(Run_Type2.eq.Run_Type_array(3:5))) &
         call system_abort('Run_Type1 must be one of "NONE", "MM", "QMMM_CORE"')
      if ( (trim(Run_Type1).eq.trim(Run_Type2) .or. &
            any(trim(Run_Type1).eq.(/'MM','QS'/))) .and. &
          trim(Run_Type2).ne.'NONE' ) then
         Run_Type2 = 'NONE'
         call print('RunType2 set to NONE')
      endif
      if ((trim(Run_Type1)).eq.'QMMM_EXT' .and..not.any(trim(Run_Type2).eq.Run_Type_array(3:5))) call system_abort('Run_Type1 must be higher level of accuracy than Run_Type2')
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
      call print('  Run_Type1 '//Run_Type1)
      call print('  Run_Type2 '//Run_Type2)
      call print('#not useable yet#  nproc '//nproc)
      call print('  IO_Rate '//IO_Rate)
      call print('  Thermostat_Type '//Thermostat_Type)
      call print('  Print_PSF '//Print_PSF)
      call print('  Time_Step '//round(Time_Step,3))
      call print('  Equilib_Time '//round(Equilib_Time,3))
      call print('  Run_Time '//round(Run_Time,3))
      call print('  Inner_QM_Radius '//round(Inner_QM_Radius,3))
      call print('  Outer_QM_Radius '//round(Outer_QM_Radius,3))
      call print('  Connect_Cutoff '//round(Connect_cutoff,3))
      call print('  Simulation_Temperature '//round(Simulation_Temperature,3))
      call print('  coord_file '//coord_file) 
      call print('  new_coord_file '//new_coord_file) 
      call print('  qm_list_filename '//qm_list_filename)
      call print('  Residue_Library '//Residue_Library) 
      call print('---------------------------------------')
      call print('')

    if (trim(Run_Type1).eq.'QS') Run_Type_1 = QS_RUN
    if (trim(Run_Type1).eq.'MM') Run_Type_1 = MM_RUN
    if (trim(Run_Type1).eq.'QMMM_CORE') Run_Type_1 = QMMM_RUN_CORE
    if (trim(Run_Type1).eq.'QMMM_EXT') Run_Type_1 = QMMM_RUN_EXTENDED
    Run_Type_2 = NONE_RUN
    if (trim(Run_Type2).eq.'MM') Run_Type_2 = MM_RUN
    if (trim(Run_Type2).eq.'QMMM_CORE') Run_Type_2 = QMMM_RUN_CORE

    call initialise(xyz,new_coord_file,action=OUTPUT)
  
    !read in the coordinates, create atoms object with connectivities
    call print('Reading in the coordinates from file '//trim(coord_file)//'...')
    call read_xyz(my_atoms,coord_file)
!!!!!!!use if there are CP2K velocities in the coord file, converts CP2K units to libAtoms units
!    call velocity_conversion(my_atoms)
!!!!!!!

    call initialise(ds,my_atoms,constraints=0)
    ds%avg_time = 10._dp
!    ds%avg_time = 50._dp
    if (Thermostat_Type.eq.1) call add_thermostat(ds,type=LANGEVIN,T=Simulation_Temperature,tau=500.0_dp)
    call finalise(my_atoms)
    call write_cp2k_info(ds) ! to compare temperature with CP2K

!!!!!!!use next lines if there's no velocity in the coord file, creates also an XYZ with CP2K unit velocities
    call rescale_velo(ds,Simulation_Temperature)
    call velocity_conversion_rev(ds%atoms)
    call print_xyz(ds%atoms,'cp2k_coord_vel.xyz',all_properties=.true.,real_format='f17.10')
    call velocity_conversion(ds%atoms)
!!!!!!!
   
    call set_cutoff(ds%atoms,Connect_Cutoff)
    call calc_connect(ds%atoms)
    call print('Number of bonds found: '//num_of_bonds(ds%atoms))
  !!!only in case of water!! otherwise ignore warnings
    call check_neighbour_numbers(ds%atoms)
    if (num_of_bonds(ds%atoms)*3.ne.ds%atoms%N*2) call print('calc_connect calculated wrong number of bonds!!')
  !!!

    !allocate force lists
    allocate(f0(3,ds%N),f1(3,ds%N),f(3,ds%N))

   !read in initial QM list for QMMM run
    if (any(Run_Type_1.eq.(/QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
       call allocate(embedlist,4,0,0,0,20)
       call read_qmlist(ds%atoms,embedlist,qm_list_filename)
       call print('QM_flag property added')

       if (Run_Type_1.eq.QMMM_RUN_EXTENDED) list_changed = extend_qmlist(ds%atoms,Inner_QM_Radius,Outer_QM_Radius)

      !center the whole system around the first QM atom, otherwise CP2K has difficulties with finding
      !the centre of QM box => QM atoms will be on the edges, nothing in the middle!
       call center_atoms(ds%atoms,embedlist%int(1,1))    !1=atomic number,2-4=shift  of  1st QM atom
    endif
    call map_into_cell(ds%atoms)

    if (Run_Type_1.ne.QS_RUN) then
       call set_value(ds%atoms%params,'Library',trim(Residue_Library))
       call calc_topology(ds%atoms,do_CHARMM=.true.)
    endif

     !----------------------------------------------------
        call print_xyz(ds%atoms,xyz,all_properties=.true.,real_format='f17.10')
!        call print_xyz(ds%atoms,xyz,properties='species:pos:velo:QM_flag',real_format='f17.10')
     !----------------------------------------------------

    call initialise(CP2K_potential,'wrapper=.true.')
!    call initialise(my_metapotential,'Simple',CP2K_potential)

!FIRST STEP - first step of velocity verlet

     n = 0

    !force

     PSF_Print = DRIVER_PRINT_AND_SAVE !in case every #th step
     if (Topology_Print.eq.TOPOLOGY_DRIVER_ONLY_STEP0) PSF_Print = DRIVER_PRINT_AND_SAVE
     if (Topology_Print.eq.TOPOLOGY_CP2K_ONLY_STEP0) PSF_Print = CP2K_PRINT_AND_SAVE
     if (Topology_Print.eq.TOPOLOGY_NO_PSF) PSF_Print=NO_PSF
     if (Topology_Print.eq.TOPOLOGY_USE_EXISTING_PSF) PSF_Print=USE_EXISTING_PSF
     if (Topology_Print.gt.0 .or. Topology_Print.eq.(-2)) PSF_Print = 1    !generate PSF (at every n-th step) and then use it

     write (args_str,'(a,i0,a,i0,a,i0,a)') 'Run_Type=',Run_Type_1,' nproc=',nproc,' PSF_Print=',PSF_Print,' cp2k_program=cp2k_serial'
call print(args_str)
     call calc(CP2K_potential,ds%atoms,e=energy,f=f1,args_str=args_str)

     PSF_Print = USE_EXISTING_PSF !every case but
     if (Topology_Print.eq.TOPOLOGY_NO_PSF) PSF_Print=NO_PSF

    !second force calculation, only for QMMM
     if (Run_Type_2.ne.NONE_RUN) then
        if (Topology_Print.ne.0) PSF_Print = 2  !use existing PSF
        write (args_str,'(a,i0,a,i0,a,i0,a)') 'Run_Type=',Run_Type_2,' nproc=',nproc,' PSF_Print=',PSF_Print,' cp2k_program=cp2k_serial'
        call calc(CP2K_potential,ds%atoms,e=energy,f=f0,args_str=args_str)
        call QUIP_combine_forces(f1,f0,f,embedlist,ds%atoms)
     else
        f=f1
     endif

     call ds_print_status(ds, 'E',energy)

    !advance verlet1
     call advance_verlet1(ds, Time_Step, f)

     call print_xyz(ds%atoms,xyz,all_properties=.true.,real_format='f17.10')
!     call print_xyz(ds%atoms,xyz,properties='species:pos:velo:QM_flag',real_format='f17.10')

!LOOP - force calc, then VV-2, VV-1

  do while (ds%t < (Equilib_Time + Run_Time))

     n = n + 1

     if (Run_Type_1.eq.QMMM_RUN_EXTENDED) list_changed = extend_qmlist(ds%atoms,Inner_QM_Radius,Outer_QM_Radius)
     if (list_changed) call print('QM list changed at time '//round(ds%t,2))
    !force
     if (Topology_Print.gt.0) then !every #th step
        if (mod(n,Topology_Print).eq.0) then    !recalc connectivity & generate PSF (at every n-th step) and then use it
           PSF_Print = DRIVER_PRINT_AND_SAVE
           if (Run_Type_2.ne.QS_RUN) then
              call print_title('Recalculate Connectivity & Topology')
              call calc_topology(ds%atoms,do_CHARMM=.true.)
              call print('CHARMM parameters added')
           endif
        endif
     endif

     write (args_str,'(a,i0,a,i0,a,i0,a)') 'Run_Type=',Run_Type_1,' nproc=',nproc,' PSF_Print=',PSF_Print,' cp2k_program=cp2k_serial'
     call calc(CP2K_potential,ds%atoms,e=energy,f=f1,args_str=args_str)

    !second force calculation
     if (Run_Type_2.ne.NONE_RUN) then
        if (Topology_Print.ne.TOPOLOGY_NO_PSF) PSF_Print = USE_EXISTING_PSF  !use existing PSF
        write (args_str,'(a,i0,a,i0,a,i0,a)') 'Run_Type=',Run_Type_2,' nproc=',nproc,' PSF_Print=',PSF_Print,' cp2k_program=cp2k_serial'
        call calc(CP2K_potential,ds%atoms,e=energy,f=f0,args_str=args_str)
        call QUIP_combine_forces(f1,f0,f,embedlist,ds%atoms)
     else
        f=f1
     endif

    !advance verlet2
     call advance_verlet2(ds, Time_Step, f)
     if (ds%t < Equilib_Time) then
        call ds_print_status(ds, 'E',energy)
     else
        call ds_print_status(ds, 'I',energy)
     end if

    !print
     if (mod(n,IO_Rate)==0) then
        call print_xyz(ds%atoms,xyz,all_properties=.true.,real_format='f17.10')
!       call print_xyz(ds%atoms,xyz,properties='species:pos:velo:QM_flag')
     end if
     
    !advance verlet1
     call advance_verlet1(ds, Time_Step, f)

     call write_cp2k_info(ds)

  enddo

  deallocate(f,f0,f1)

  call finalise(ds)
  call finalise(embedlist)
  call finalise(xyz)
  call finalise(my_metapotential)
  call finalise(CP2K_potential)

  call print_title('THE')
  call print('Finished. CP2K is now having a rest, since deserved it. Bye-Bye!')
  call print_title('END')

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

end program test_CP2K
