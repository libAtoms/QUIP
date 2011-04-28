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

!adaptive hybrid QMMM MD program.
!reads atoms & QM list & runs CP2K & writes movie xyz with QM flags
!uses the potential with momentum conservation, adding the corrective force only to the QM core and buffer atoms
!filepot_program can be e.g. /Users/csilla/QUIP/build.darwin_x86_64_g95/cp2k_filepot
!TODO:
!use the buffer carving of the potential force mixing routines.
#include "error.inc"

program qmmm_md

  use libatoms_module
  use cp2k_driver_template_module
  use potential_module
  use restraints_constraints_xml_module

  implicit none

  type(DynamicalSystem)               :: ds
  type(Atoms)                         :: my_atoms
  integer                             :: i, n
  character(len=FIELD_LENGTH)         :: Run_Type_array(5)               !_MM_, QS, QMMM_EXTENDED or QMMM_CORE

  !Force calc.
  type(Potential)                     :: cp2k_fast_pot, cp2k_slow_pot 
  type(Potential)                     :: evbsub_cp2k_fast_pot, evbsub_cp2k_slow_pot 
  type(Potential)                     :: pot, empty_qm_pot
  type(Potential)                     :: evbsub_pot, evbsub_empty_qm_pot
  real(dp)                            :: energy,check,TI_force, TI_corr
  real(dp), dimension(:,:), allocatable :: f,f0,f1, add_force

  !QM list generation
  !logical                             :: list_changed1
  !integer, pointer                    :: hybrid_mark_p(:)
  !integer, pointer                    :: cluster_mark_p(:)

  !Thermostat
  real(dp)                            :: temp

  !Spline

  type spline_pot
    real(dp) :: from
    real(dp) :: to
    real(dp) :: dpot
  end type spline_pot

  type(spline_pot)                    :: my_spline
  real(dp)                            :: spline_pot_val
  real(dp), pointer                   :: spline_pot_val_p(:)

  !Output XYZ
  type(CInoutput)                      :: traj_xyz, latest_xyz
!  type(CInoutput)                      :: xyz
  character(len=FIELD_LENGTH)         :: backup_coord_file          !output XYZ file
  integer                             :: backup_i

  !Topology
  character(len=FIELD_LENGTH)         :: driver_PSF_Print
  integer                             :: Topology_Print_rate     !_-1_ never, 0 print one at the 0th time step and use that
                                                                 ! n>0 print at 0th and then every n-th step
  type(Table)                         :: intrares_impropers

  !Input parameters
  type(Dictionary)            :: params_in
  character(len=FIELD_LENGTH) :: Run_Type1               !_MM_, QS, QMMM_EXTENDED or QMMM_CORE
  character(len=FIELD_LENGTH) :: Run_Type2               !_NONE_, MM, or QMMM_CORE
  integer                     :: IO_Rate                 !print coordinates at every n-th step
  integer                     :: Thermostat_Type         !_0_ none, 1 Langevin
  real(dp)                    :: Thermostat_7_rs(2)
  character(len=FIELD_LENGTH) :: PSF_Print               !_NO_PSF_, DRIVER_AT_0, DRIVER_EVERY_#, USE_EXISTING_PSF
  real(dp)                    :: Time_Step
  real(dp)                    :: Equilib_Time
  real(dp)                    :: Run_Time
  real(dp)                    :: Inner_Buffer_Radius         !for hysteretic quantum
  real(dp)                    :: Outer_Buffer_Radius         !          region selection
  real(dp)                    :: Inner_QM_Region_Radius, Inner_QM_Core_Region_Radius
  real(dp)                    :: Outer_QM_Region_Radius, Outer_QM_Core_Region_Radius
  real(dp)                    :: Connect_Cutoff
  real(dp)                    :: Simulation_Temperature
  character(len=FIELD_LENGTH) :: coord_file
  character(len=FIELD_LENGTH) :: latest_coord_file          !output XYZ file
  character(len=FIELD_LENGTH) :: traj_file          !output XYZ file
  character(len=FIELD_LENGTH) :: qm_list_filename, qm_core_list_filename        !QM list file with a strange format
  type(Table)                 :: qm_seed_extended, qm_seed_core
  character(len=FIELD_LENGTH) :: Residue_Library
  character(len=FIELD_LENGTH) :: restraint_constraint_xml_file
  type(Extendable_Str)        :: restraint_constraint_xml_es
  integer                     :: Charge
  real(dp)                    :: Langevin_Tau, Nose_Hoover_Tau
  logical                     :: Buffer_general, do_general
  logical                     :: Continue_it
  logical                     :: reinitialise_qm_region
  real(dp)                    :: avg_time
  integer                     :: Seed
  logical                     :: qm_region_pt_ctr, qm_core_region_pt_ctr
  integer                     :: qm_region_atom_ctr, qm_core_region_atom_ctr
  character(len=FIELD_LENGTH) :: print_prop
  logical                     :: print_forces
  logical                     :: print_forces_at0
  real(dp)                    :: nneightol
  logical                     :: Delete_Metal_Connections
  real(dp)                    :: spline_from
  real(dp)                    :: spline_to
  real(dp)                    :: spline_dpot
  logical                     :: use_spline
  integer                     :: max_n_steps
  integer                     :: update_QM_region_interval
  character(len=FIELD_LENGTH) :: cp2k_calc_args               ! other args to calc(cp2k,...)
  character(len=FIELD_LENGTH) :: evb_args_str                 ! args to calc(EVB(cp2k))
  character(len=FIELD_LENGTH) :: extra_calc_args                 ! extra args to mm and qm calcs (like FilePot_log)
  character(len=FIELD_LENGTH) :: filepot_program

  logical                     :: do_carve_cluster
  logical                     :: use_create_cluster_info_for_core
  character(len=FIELD_LENGTH) :: core_create_cluster_info_args
  real(dp) :: qm_region_ctr(3), qm_core_region_ctr(3)
  real(dp) :: use_cutoff
  real(dp) :: H_extra_heat_ctr(3), H_extra_heat_r(2), H_extra_heat_velo_factor
  real(dp) :: nH_extra_heat_ctr(3), nH_extra_heat_r(2), nH_extra_heat_velo_factor
  logical :: EVB

  logical :: distance_ramp
  real(dp) :: distance_ramp_inner_radius, distance_ramp_outer_radius
  character(len=128) :: weight_interpolation

  real(dp) :: max_move_since_calc_connect
  real(dp) :: calc_connect_buffer
logical :: have_silica_potential
  integer :: stat

  logical have_QMMM_core, have_QMMM_extended

  integer :: error=ERROR_NONE

  real(dp), allocatable :: restraint_stuff(:,:)
  real(dp) :: r
  integer :: tmp_run_dir_i, copy_latest_every, copy

!    call system_initialise(verbosity=PRINT_ANAL,enable_timing=.true.)
!    call system_initialise(verbosity=PRINT_NERD,enable_timing=.true.)
    call system_initialise(verbosity=PRINT_VERBOSE,enable_timing=.true.)
!    call system_initialise(verbosity=PRINT_NORMAL,enable_timing=.true.)
    call system_timer('program')

    !INPUT
      call initialise(params_in)
      call param_register(params_in, 'Run_Type1', 'MM', Run_Type1, help_string="Type of run1 (or only run if Run_Type2=NONE) QS, MM, QMMM_CORE, QMMM_EXTENDED")
      call param_register(params_in, 'Run_Type2', 'NONE', Run_Type2, help_string="Type of run1 (or only run if Run_Type2=NONE) QS, MM, QMMM_CORE, QMMM_EXTENDED")
      call param_register(params_in, 'IO_Rate', '1', IO_Rate, help_string="Frequency of output printing (every IO_Rate steps).")
      call param_register(params_in, 'Thermostat_Type', '0', Thermostat_Type, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'Thermostat_7_rs', '0.0 0.0', Thermostat_7_rs, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'PSF_Print', 'NO_PSF', PSF_Print, help_string="when to print PSF file: NO_PSF, EVERY_#, DRIVER_PRINT_AT_0, USE_EXISTING_PSF")
      call param_register(params_in, 'Time_Step', '0.5', Time_Step, help_string="Time step in fs.")
      call param_register(params_in, 'Equilib_Time', '0.0', Equilib_Time, help_string="Equilibration time in fs.")
      call param_register(params_in, 'Run_Time', '0.5', Run_Time, help_string="Run time in fs.")
      call param_register(params_in, 'Inner_Buffer_Radius', '0.0', Inner_Buffer_Radius, help_string="The inner radius of the hysteretic buffer of the QM region.")
      call param_register(params_in, 'Outer_Buffer_Radius', '0.0', Outer_Buffer_Radius, help_string="The outer radius of the hysteretic buffer of the QM region.")
      call param_register(params_in, 'Inner_QM_Region_Radius', '0.0', Inner_QM_Region_Radius, help_string="The inner radius of the QM region (around an atom list or a point in space).")
      call param_register(params_in, 'Outer_QM_Region_Radius', '0.0', Outer_QM_Region_Radius, help_string="The outer radius of the QM region (around an atom list or a point in space).")
      call param_register(params_in, 'Inner_QM_Core_Region_Radius', '0.0', Inner_QM_Core_Region_Radius, help_string="The inner radius of the core QM region (around an atom list or a point in space).")
      call param_register(params_in, 'Outer_QM_Core_Region_Radius', '0.0', Outer_QM_Core_Region_Radius, help_string="The outer radius of the core QM region (around an atom list or a point in space).")
      call param_register(params_in, 'Connect_Cutoff', '0.0', Connect_cutoff, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'Simulation_Temperature', '300.0', Simulation_Temperature, help_string="Simulation temperature in K.")
      call param_register(params_in, 'coord_file', 'coord.xyz',coord_file, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'latest_coord_file', 'latest.xyz',latest_coord_file, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'traj_file', 'movie.xyz',traj_file, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'Residue_Library', 'all_res.CHARMM.lib',Residue_Library, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'restraint_constraint_xml_file', '', restraint_constraint_xml_file, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'Charge', '0', Charge, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'Langevin_Tau', '500.0', Langevin_Tau, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'Nose_Hoover_Tau', '74.0', Nose_Hoover_Tau, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'Buffer_general', 'F', Buffer_general, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'Continue', 'F', Continue_it, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'reinitialise_qm_region', 'F', reinitialise_qm_region, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'avg_time', '10.0', avg_time, help_string="Averaging time for the calculation of averaged coordinates to smooth out bond vibrations for the topology calculation")
      call param_register(params_in, 'Seed', '-1', Seed, help_string="No help yet.  This source file was $LastChangedBy$")

      call param_register(params_in, 'qm_list_filename', '', qm_list_filename,  &
	help_string="Filename containing list of atoms for (extended) QM region (in funny format).  Mutually exclusive with qm_region_pt_ctr and qm_region_atom_ctr")
      call param_register(params_in, 'qm_region_pt_ctr', 'F', qm_region_pt_ctr, &
	help_string="Use Cartesian point to center (extended) QM region.  Mutually exclusive with qm_list_filename and qm_region_atom_ctr")
      call param_register(params_in, 'qm_region_ctr', '(/0.0 0.0 0.0/)', qm_region_ctr, help_string="Centre of the (extended) QM region in space.")
      call param_register(params_in, 'qm_region_atom_ctr', '0', qm_region_atom_ctr, help_string="Atom to center (extended) QM region about.  Mutually exclusive with qm_region_pt_ctr and qm_list_filename")

      call param_register(params_in, 'qm_core_list_filename', '', qm_core_list_filename, &
	help_string="Filename containing list of atoms for core QM region (in funny format).  Mutually exclusive with qm_core_region_pt_ctr and qm_core_region_atom_ctr")
      call param_register(params_in, 'qm_core_region_pt_ctr', 'F', qm_core_region_pt_ctr, &
	help_string="Use Cartesian point to center core QM region.  Mutually exclusive with qm_core_list_filename and qm_core_region_atom_ctr")
      call param_register(params_in, 'qm_core_region_ctr', '(/0.0 0.0 0.0/)', qm_core_region_ctr, help_string="Centre of the QM core region in space.")
      call param_register(params_in, 'qm_core_region_atom_ctr', '0', qm_core_region_atom_ctr, &
	help_string="Atom to center core QM region about.  Mutually exclusive with qm_core_region_pt_ctr and qm_core_list_filename")

      call param_register(params_in, 'print_prop', 'all', print_prop, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'print_forces', 'T', print_forces, help_string="whether to print forces at time>0 on the mainlog")
      call param_register(params_in, 'print_forces_at0', 'T', print_forces_at0, help_string="whether to print forces at time=0 on the mainlog")
      call param_register(params_in, 'nneightol', '1.2', nneightol, help_string="Nearest neighbour tolerance factor.  Two atoms are bonded if d(a1,a2)<[CovRad(a1)+CovRad(a2)]*Neighbour_Tolerance.")
      call param_register(params_in, 'Delete_Metal_Connections', 'T', Delete_Metal_Connections, help_string="Whether not to take into account bonds of non(C,H,O,N,P) elements.")
      call param_register(params_in, 'spline_from', '0.0', spline_from, help_string="inner radius of the external spline potential")
      call param_register(params_in, 'spline_to', '0.0', spline_to, help_string="outer radius of the external spline potential")
      call param_register(params_in, 'spline_dpot', '0.0', spline_dpot, help_string="depth of the external spline potential")
      call param_register(params_in, 'use_spline', 'F', use_spline, help_string="whether to use an external spline potential within a spherical shell around (/0. 0. 0./).")
      call param_register(params_in, 'max_n_steps', '-1', max_n_steps, help_string="Maximum number of time steps to run")
      call param_register(params_in, 'update_QM_region_interval', '1', update_QM_region_interval, help_string="interval between updates of QM region")
      call param_register(params_in, 'cp2k_calc_args', '', cp2k_calc_args, help_string="calc args for cp2k_filepot")
      call param_register(params_in, 'evb_args_str', '', evb_args_str, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'extra_calc_args', '', extra_calc_args, help_string="Extra args for mm and qm calc(), for things like FilePot_log, and arguments for buffer create_cluster_info")
      call param_register(params_in, 'filepot_program', param_mandatory, filepot_program, help_string="Filepot program.")
      call param_register(params_in, 'carve_cluster', 'F', do_carve_cluster, help_string="If true, carve cluster instead of doing QM/MM.  Untested")
      call param_register(params_in, 'use_create_cluster_info_for_core', 'F', use_create_cluster_info_for_core, help_string="If true, use cluster_create_inf() to create core, which does things like residue completing heuristics")
      call param_register(params_in, 'core_create_cluster_info_args', '', core_create_cluster_info_args, help_string="Arguments to pass to create_cluster_info for QM core, if used by update_QM_region")
      call param_register(params_in, 'calc_connect_buffer', '0.2', calc_connect_buffer, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'have_silica_potential', 'F', have_silica_potential, help_string="Whether there is a silica unit in the system to be treated with the Danny potential (implemented in CP2K).")
      call param_register(params_in, 'EVB', 'F', EVB, help_string="Whether to use the EVB MM potential instead of a simple MM.")

      call param_register(params_in, 'distance_ramp', 'F', distance_ramp, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'distance_ramp_inner_radius', '3.0', distance_ramp_inner_radius, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'distance_ramp_outer_radius', '4.0', distance_ramp_outer_radius, help_string="No help yet.  This source file was $LastChangedBy$")

      call param_register(params_in, 'H_extra_heat_ctr', '0.0 0.0 0.0', H_extra_heat_ctr, help_string="spherical centre of the spherical shell where extra heating is introduced on hydrogen atoms")
      call param_register(params_in, 'H_extra_heat_r', '0.0 -1.0', H_extra_heat_r, help_string="inner and outer radii of the spherical shell where extra heating is introduced on hydrogen atoms")
      call param_register(params_in, 'H_extra_heat_velo_factor', '1.0', H_extra_heat_velo_factor, help_string="velocities of hydrogen atoms will be multiplied by this every time step")
      call param_register(params_in, 'nH_extra_heat_ctr', '0.0 0.0 0.0', nH_extra_heat_ctr, help_string="spherical centre of the spherical shell where extra heating is introduced on nonhydrogen atoms")
      call param_register(params_in, 'nH_extra_heat_r', '0.0 -1.0', nH_extra_heat_r, help_string="inner and outer radii of the spherical shell where extra heating is introduced on nonhydrogen atoms")
      call param_register(params_in, 'nH_extra_heat_velo_factor', '1.0', nH_extra_heat_velo_factor, help_string="velocities of nonhydrogen atoms will be multiplied by this every time step")
      call param_register(params_in, 'tmp_run_dir_i', '-1', tmp_run_dir_i, help_string="if >0, the cp2k run directory will be /tmp/cp2k_run_$tmp_run_dir_i$, and all input files are also copied here when first called")
      call param_register(params_in, 'copy_latest_every', '-1', copy_latest_every, help_string="if >0, the copy the latest.xyz and the wfn files from the run directory /tmp/cp2k_run_$tmp_run_dir_i$ back to the home at this frequency -- only active when $tmp_run_dir_i$>0")

      if (.not. param_read_args(params_in)) then
        call system_abort('could not parse argument line')
      end if

      call finalise(params_in)

    !CHECK INPUT PARAMETERS

      if (Run_Type1(1:13) == 'QMMM_EXTENDED' .or. Run_Type2(1:13) == 'QMMM_EXTENDED') then
	if (count((/qm_region_pt_ctr, qm_region_atom_ctr /= 0, len_trim(qm_list_filename) /= 0 /)) /= 1) then
	    call system_abort("Doing Run_Type1="//trim(Run_Type1)//" Run_Type2="//trim(Run_Type2)//&
			      ", need exactly one of qm_region_pt_ctr, qm_region_atom_ctr="//qm_region_atom_ctr// &
			      "/=0, len_trim(qm_list_filename='"//trim(qm_list_filename)//"') /= 0")
	endif
      endif
      if (Run_Type1(1:9) == 'QMMM_CORE' .or. Run_Type2(1:9) == 'QMMM_CORE') then
	if (count((/qm_core_region_pt_ctr, qm_core_region_atom_ctr /= 0, len_trim(qm_core_list_filename) /= 0 /)) /= 1) then
	    call system_abort("Doing Run_Type1="//trim(Run_Type1)//" Run_Type2="//trim(Run_Type2)//&
			      ", need exactly one of qm_core_region_pt_ctr, qm_core_region_atom_ctr="//qm_core_region_atom_ctr// &
			      "/=0, len_trim(qm_core_list_filename='"//trim(qm_core_list_filename)//"') /= 0")
	endif
      endif

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
      if ((trim(Run_Type1)).eq.'QMMM_EXTENDED' .and..not.any(trim(Run_Type2).eq.Run_Type_array(3:5))) &
	call system_abort('Run_Type1 must be higher level of accuracy than Run_Type2')
      if ((trim(Run_Type1)).eq.'QMMM_CORE' .and..not.any(trim(Run_Type2).eq.Run_Type_array(4:5))) &
	call system_abort('Run_Type1 must be higher level of accuracy than Run_Type2')

!check PSF printing
      if (trim(PSF_Print) == 'NO_PSF' .or. trim(PSF_Print) == 'USE_EXISTING_PSF') then
	 topology_print_rate=-1
	 driver_PSF_Print=PSF_Print
      else if (trim(PSF_Print) == 'DRIVER_AT_0') then
	 topology_print_rate=0
      else if (len(PSF_Print) > 13) then
	 if (PSF_Print(1:13) == 'DRIVER_EVERY_') then
	    read(unit=PSF_Print(14:len_trim(PSF_Print)),fmt=*,iostat=stat) topology_print_rate
	    if (stat /= 0) &
	       call system_abort("PSF_Print='"//trim(PSF_Print)//"' unable to parse N from DRIVER_EVERY_N '"// &
	 	 PSF_Print(14:len_trim(PSF_Print))//"'")
	 else
	   call system_abort("Unknown PSF_Print '"//trim(PSF_Print)//"'")
	 endif
      else
	 call system_abort("Unknown PSF_Print '"//trim(PSF_Print)//"'")
      endif

    !PRINT INPUT PARAMETERS
      call print('Run parameters:')
      call print('  filepot_program '//trim(filepot_program))
      call print('  Run_Type1 '//Run_Type1)
      call print('  Run_Type2 '//Run_Type2)
      call print('  EVB (MMpot) '//EVB)
      call print('  IO_Rate '//IO_Rate)
      if (Thermostat_Type.eq.0) then
         call print('  Thermostat_Type 0: None, NVE')
      elseif (Thermostat_Type.eq.1) then
         call print('  Thermostat_Type 1: Langevin everywhere')
         call print('  Langevin_Tau '//Langevin_Tau)
      elseif (Thermostat_Type.eq.2) then
         call print('  Thermostat_Type 1: Nose-Hoover everywhere')
         call print('  Nose_Hoover_Tau '//Nose_Hoover_Tau)
      elseif (Thermostat_Type.eq.3) then
         call print('  Thermostat_Type 3: separate Nose-Hoover for each atom')
         call print('  Nose_Hoover_Tau '//Nose_Hoover_Tau)
      elseif (Thermostat_Type.eq.4) then
         call print('  Thermostat_Type 4: separate Nose-Hoover-Langevin for each atom')
         call print('  Langevin_Tau ' // Langevin_Tau // ' Nose_Hoover_Tau '//Nose_Hoover_Tau)
      elseif (Thermostat_Type.eq.5) then
         call print('  Thermostat_Type 5: QM core & buffer heavy atoms in the 1st thermostat')
         call print('                     QM core & buffer H in the 2nd thermostat')
         call print('                     classical O & H in the 3rd thermostat')
         call print('  Langevin_Tau '//Langevin_Tau // ', Nose_Hoover_Tau ' // Nose_Hoover_Tau)
      elseif (Thermostat_Type.eq.6) then
         call print('  Thermostat_Type 6: 3 regions (QM, buffer, MM) x  2 types (H, heavy)')
         call print('                   each with its own Nose-Hoover thermostat')
         call print('  Nose_Hoover_Tau ' // Nose_Hoover_Tau)
      elseif (Thermostat_Type.eq.7) then
         call print('  Thermostat_Type 7: 3 regions, by radius,')
         call print('                   each with its own Nose-Hoover thermostat')
         call print('  Nose_Hoover_Tau ' // Nose_Hoover_Tau)
      else
	 call print('  Thermostat_Type '//thermostat_type//' unknown')
      endif
      call print('  PSF_Print '//PSF_Print)
      call print('  nneightol '//nneightol)
      call print('  Time_Step '//round(Time_Step,3))
      call print('  Equilib_Time '//round(Equilib_Time,3))
      call print('  Run_Time '//round(Run_Time,3))
      call print('  max_n_steps '//max_n_steps)
      call print('  Inner_Buffer_Radius '//round(Inner_Buffer_Radius,3))
      call print('  Outer_Buffer_Radius '//round(Outer_Buffer_Radius,3))
      call print('! - not used any more -  Connect_Cutoff '//round(Connect_cutoff,3))
      call print('  Simulation_Temperature '//round(Simulation_Temperature,3))
      call print('  coord_file '//coord_file) 
      call print('  latest_coord_file '//latest_coord_file) 
      call print('  traj_file '//traj_file) 
      if (len_trim(qm_list_filename) /= 0) then
         call print('  qm_list_filename '//trim(qm_list_filename))
      else if (qm_region_pt_ctr .or. qm_region_atom_ctr /= 0) then
	 if (qm_region_pt_ctr) then
	   call print('  QM extended region is centred around qm_region_ctr= '//qm_region_ctr)
	 else
	   call print('  QM extended region is centred around atom '//qm_region_atom_ctr)
	   qm_region_pt_ctr = .true.
	 endif
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
      endif
      if (len_trim(qm_core_list_filename) /= 0) then
         call print('  qm_core_list_filename '//trim(qm_core_list_filename))
      else if (qm_core_region_pt_ctr .or. qm_core_region_atom_ctr /= 0) then
	 if (qm_core_region_pt_ctr) then
	   call print('  QM core is centred around qm_core_region_ctr= '//qm_core_region_ctr)
	 else
	   call print('  QM core is centred around atom '//qm_core_region_atom_ctr)
	   qm_core_region_pt_ctr = .true.
	 endif
         call print('  Inner_QM_Core_Region_Radius '//round(Inner_QM_Core_Region_Radius,3))
         call print('  Outer_QM_Core_Region_Radius '//round(Outer_QM_Core_Region_Radius,3))
      endif
      call print('  Residue_Library '//Residue_Library) 
      call print('  restraint_constraint_xml_file '//restraint_constraint_xml_file) 
      call print('  Charge '//Charge)
      call print('  Buffer_general '//Buffer_general)
      call print('  Continue '//Continue_it)
      call print('  reinitialise_qm_region '//reinitialise_qm_region)
      call print('  avg_time '//avg_time)
      call print('  Seed '//Seed)
      call print('  Properties to print '//trim(print_prop))
      call print('  Print forces? '//print_forces)
      call print('  Print forces at t=0? '//print_forces_at0)
      call print('  carve_cluster '//do_carve_cluster)
      call print('  have_silica_potential '//have_silica_potential)
      call print('  evb_args_str '//trim(evb_args_str))
      call print('  cp2k_calc_args '//trim(cp2k_calc_args))
      call print('  filepot_program '//trim(filepot_program))
      call print('  tmp_run_dir_i '//tmp_run_dir_i)
      call print('  copy_latest_every '//copy_latest_every)
      call print('---------------------------------------')
      call print('')

! STARTS HERE

    if (is_file_readable(trim(traj_file))) then
      call print("WARNING: traj_file " // trim(traj_file) // " exists, backing it up")
      backup_i=1
      backup_coord_file=trim(traj_file)//".backup_"//backup_i
      do while (is_file_readable(trim(backup_coord_file)))
	backup_i = backup_i + 1
	backup_coord_file=trim(traj_file)//".backup_"//backup_i
      end do
      call print("WARNING:      to backup_coord_file " // trim(backup_coord_file))
      call system("cp "//trim(traj_file)//" "//trim(backup_coord_file))
    endif

    call initialise(traj_xyz,traj_file,action=OUTPUT)
  
  !READ COORDINATES

    call print('Reading in the coordinates from file '//trim(coord_file)//'...')
    call read(my_atoms,coord_file)

    call initialise(ds,my_atoms)

    if (Continue_it) then
      if (get_value(my_atoms%params,'Time',ds%t)) then
	  call print('Found Time in atoms%params'//ds%t)
      endif
    endif

    if (len_trim(restraint_constraint_xml_file) > 0) then
       call initialise(restraint_constraint_xml_es)
       call read(restraint_constraint_xml_es, restraint_constraint_xml_file, convert_to_string=.true.)
       call init_restraints_constraints(ds, string(restraint_constraint_xml_es))
       call finalise(restraint_constraint_xml_es)
       if (ds%Nrestraints > 0) allocate(restraint_stuff(5,ds%Nrestraints))
    endif

    ds%avg_time = avg_time

  !THERMOSTAT
    call add_thermostats(ds, thermostat_type, Simulation_Temperature)
    call set_thermostat_masses(ds%atoms, Thermostat_Type, Simulation_Temperature, &
      Inner_QM_Region_Radius, Outer_QM_Region_Radius, Inner_Buffer_Radius, Outer_Buffer_Radius, Thermostat_7_rs)

    call finalise(my_atoms)
    call add_property(ds%atoms,'pot',0._dp) ! always do this, it's just 0 if spline isn't active - no need to change print_props

    if (.not. has_property(ds%atoms, 'force')) call add_property(ds%atoms, 'force', 0.0_dp, n_cols=3)

  !PRINT CONSTRAINTS AND RESTRAINTS

   if (ds%Nconstraints > 0) then
      call print("Constraints:")
      call print(ds%constraint)
   end if
   if (ds%Nrestraints > 0) then
      call print("Restraints:")
      call print(ds%restraint)
   end if

  !SET VELOCITIES

    if (.not.Continue_it) then
       call rescale_velo(ds,Simulation_Temperature)
    endif

  !CALC. CONNECTIONS

!    call set_cutoff(ds%atoms,Connect_Cutoff)
!    call set_cutoff(ds%atoms,0._dp) !use the covalent radii to determine bonds
    use_cutoff = max(nneightol, Outer_Buffer_Radius)
    use_cutoff = max(use_cutoff, Outer_QM_Region_Radius)
    if (distance_ramp) use_cutoff = max(use_cutoff, distance_ramp_outer_radius)
    if (have_silica_potential) then
	use_cutoff = max(SILICA_2BODY_CUTOFF, Outer_Buffer_Radius)
    endif
    call set_cutoff(ds%atoms,use_cutoff+calc_connect_buffer)
    call calc_connect(ds%atoms)
    if (Delete_Metal_Connections) call delete_metal_connects(ds%atoms)

  !READ / CREATE QM LIST + THERMOSTATTING

   !QM BUFFER + THERMOSTATTING
      ! set general/heavy atom selection before QM region selection
!       call set_value(ds%atoms%params,'Buffer_general',Buffer_general)
!       call print('set Buffer_general into ds%atoms%params')
!       if (get_value(ds%atoms%params,'Buffer_general',do_general)) then
!           call print('Found Buffer_general in atoms%params'//do_general)
!           buffer_general=do_general
!       else
!           call print('Not found Buffer_general in atoms%params')
!           buffer_general=.false.
!       endif
    do_general=Buffer_general

    ds%atoms%nneightol = nneightol
    ! topology calculation
    if (trim(Run_Type1).ne.'QS') then
       if (.not.Continue_it) then
	  ! residues are needed here for heuristics based on them (e.g. don't break residue, don't break molecule)
	  ! in update_QM_region() (and also buffer cluster creating in potential_simple_calc())
          call set_value(ds%atoms%params,'Library',trim(Residue_Library))
	  call map_into_cell(ds%atoms)
	  call calc_dists(ds%atoms)
          call create_residue_labels_arb_pos(ds%atoms,do_CHARMM=.true.,intrares_impropers=intrares_impropers,have_silica_potential=have_silica_potential)
          call check_topology(ds%atoms)
       endif
    endif

   !QM CORE
    have_QMMM_core=.false.
    if (trim(Run_Type1) == 'QMMM_CORE' .or. trim(Run_Type2) == 'QMMM_CORE') have_QMMM_core =.true.
    have_QMMM_extended=.false.
    if (trim(Run_Type1) == 'QMMM_EXTENDED' .or. trim(Run_Type2) == 'QMMM_EXTENDED') have_QMMM_extended =.true.

    if (have_QMMM_core .or. have_QMMM_extended) then

       !even if continue, we need the QM seed. -- no marks overwritten!
       if (Continue_it .and. len_trim(qm_list_filename)/=0) call read_qmlist(ds%atoms,qm_list_filename,qm_seed_extended)
       if (Continue_it .and. len_trim(qm_core_list_filename)/=0) call read_qmlist(ds%atoms,qm_core_list_filename,qm_seed_core)

       !(re)initialise QM region -- marks overwritten!
       if (reinitialise_qm_region .or. .not. Continue_it) then
	  if (have_QMMM_extended) call update_QM_region(ds%atoms, '_extended', qm_region_pt_ctr, qm_region_ctr, qm_region_atom_ctr, &
	    qm_list_filename, Inner_QM_Region_Radius, Outer_QM_Region_Radius, buffer_general, use_create_cluster_info_for_core, &
	    core_create_cluster_info_args, .true., qm_seed_extended)
	  if (have_QMMM_core) call update_QM_region(ds%atoms, '_core', qm_core_region_pt_ctr, qm_core_region_ctr, qm_core_region_atom_ctr, &
	    qm_core_list_filename, Inner_QM_Core_Region_Radius, Outer_QM_Core_Region_Radius, buffer_general, use_create_cluster_info_for_core, &
	    core_create_cluster_info_args, .true., qm_seed_core)

          call print('hybrid, hybrid_mark and old_hybrid_mark properties added')
       endif ! .not. Continue_it

    endif ! have_QMMM_core || have_QMMM_extended
    call map_into_cell(ds%atoms)
    call calc_dists(ds%atoms)

  !TOPOLOGY

  !CHARGE

    if (trim(Run_Type1).eq.'QS') then
       call set_value(ds%atoms%params,'Charge',Charge)
    endif

  !INIT. POTENTIAL

    !only QMMM_EXTENDED for the moment **************
    !if (trim(Run_Type1).ne.'QMMM_EXTENDED') call system_abort('ONLY QMMM_EXTENDED')

    ! set up pot
    if (trim(Run_Type2) == 'NONE') then ! no force mixing
       if (EVB.and.trim(Run_Type1)=='MM') then
	  call print("SETUP_POT EVBSUB_POT")
          call setup_pot(evbsub_pot, Run_Type1, filepot_program, tmp_run_dir_i)
          call initialise(pot,args_str='EVB=T',  &
             !' topology_suffix1=_EVB1 topology_suffix2=_EVB2'// &
             !' form_bond={1 2} break_bond={2 6} ', &
             pot1=evbsub_pot)
	  call print("INITIALISE_POT EVB")
       else
	  call print("SETUP_POT POT")
          call setup_pot(pot, Run_Type1, filepot_program, tmp_run_dir_i)
       endif
       ! set up mm only pot, in case we need it for empty QM core
       if (EVB) then
	  call print("SETUP_POT EVBSUB_EMPTY_QM_POT")
          call setup_pot(evbsub_empty_qm_pot, 'MM', filepot_program, tmp_run_dir_i)
          call initialise(empty_qm_pot,args_str='EVB=T',  &
             !' topology_suffix1=_EVB1 topology_suffix2=_EVB2'// &
             !' form_bond={1 2} break_bond={2 6} ', &
             pot1=evbsub_empty_qm_pot)
	  call print("INITIALISE_POT EVB")
       else
	  call print("SETUP_POT EMPTY_QM_POT")
          call setup_pot(empty_qm_pot, 'MM', filepot_program, tmp_run_dir_i)
       endif
    else ! doing force mixing
       if (EVB.and.trim(Run_Type1)=='MM') then
          call print("WARNING: Force Mixing with MM as the slow potential!")
	  call print("SETUP_POT EVBSUB_CP2K_SLOW_POT")
          call setup_pot(evbsub_cp2k_slow_pot, Run_Type1, filepot_program, tmp_run_dir_i)
          call initialise(cp2k_slow_pot,args_str='EVB=T',  &
             !' topology_suffix1=_EVB1 topology_suffix2=_EVB2'// &
             !' form_bond={1 2} break_bond={2 6} ', &
             pot1=evbsub_cp2k_slow_pot)
	  call print("INITIALISE_POT EVB")
       else
	  call print("SETUP_POT CP2K_SLOW_POT")
          call setup_pot(cp2k_slow_pot, Run_Type1, filepot_program, tmp_run_dir_i)
       endif
       if (EVB.and.trim(Run_Type2)=='MM') then
	  call print("SETUP_POT EVBSUB_CP2K_FAST_POT")
          call setup_pot(evbsub_cp2k_fast_pot, Run_Type2, filepot_program, tmp_run_dir_i)
          call initialise(cp2k_fast_pot,args_str='EVB=T',  &
             !' topology_suffix1=_EVB1 topology_suffix2=_EVB2'// &
             !' form_bond={1 2} break_bond={2 6} ', &
             pot1=evbsub_cp2k_fast_pot)
	  call print("INITIALISE_POT EVB")
       else
	  call print("SETUP_POT CP2K_FAST_POT")
          call setup_pot(cp2k_fast_pot, Run_Type2, filepot_program, tmp_run_dir_i)
       endif
       if (distance_ramp) then
	 if (.not. qm_region_pt_ctr) call system_abort("Distance ramp needs qm_region_pt_ctr (or qm_region_atom_ctr)")
	 weight_interpolation='distance_ramp'
       else
	 weight_interpolation='hop_ramp'
       endif
       call initialise(pot,args_str='ForceMixing=T use_buffer_for_fitting=T add_cut_H_in_fitlist=T'// &
	  ' method=conserve_momentum conserve_momentum_weight_method=mass calc_weights=T'// &
	  ' min_images_only=F nneighb_only=F lotf_nneighb_only=F fit_hops=1 hysteretic_buffer=T'// &
	  ' hysteretic_buffer_inner_radius='//Inner_Buffer_Radius// &
	  ' hysteretic_buffer_outer_radius='//Outer_Buffer_Radius// &
	  ' weight_interpolation='//trim(weight_interpolation)// &
	  ' distance_ramp_inner_radius='//distance_ramp_inner_radius//' distance_ramp_outer_radius='//distance_ramp_outer_radius// &
	  ' single_cluster=T little_clusters=F carve_cluster='//do_carve_cluster &
!next line is for playing with silica carving
!          //' even_electrons=T terminate=T cluster_same_lattice=T termination_clash_check=T' &
	  //' construct_buffer_use_only_heavy_atoms='//(.not.(buffer_general)), &
	  pot1=cp2k_fast_pot, pot2=cp2k_slow_pot)
       call print('INITIALISE_POT ForceMixing=T use_buffer_for_fitting=T add_cut_H_in_fitlist=T'// &
	  ' method=conserve_momentum conserve_momentum_weight_method=mass calc_weights=T'// &
	  ' min_images_only=F nneighb_only=F lotf_nneighb_only=F fit_hops=1 hysteretic_buffer=T'// &
	  ' hysteretic_buffer_inner_radius='//Inner_Buffer_Radius// &
	  ' hysteretic_buffer_outer_radius='//Outer_Buffer_Radius// &
	  ' weight_interpolation='//trim(weight_interpolation)// &
	  ' distance_ramp_inner_radius='//distance_ramp_inner_radius//' distance_ramp_outer_radius='//distance_ramp_outer_radius// &
	  ' single_cluster=T little_clusters=F carve_cluster='//do_carve_cluster // &
	  ' construct_buffer_use_only_heavy_atoms='//(.not.(buffer_general)))

       ! if Run_Type2 = QMMM_CORE, we'll crash if QM core is ever empty
       if (trim(Run_Type2) == 'MM') then
          if (EVB) then
	     call print("SETUP_POT EVBSUB_EMPTY_QM_POT")
             call setup_pot(evbsub_empty_qm_pot, Run_Type2, filepot_program, tmp_run_dir_i)
             call initialise(empty_qm_pot,args_str='EVB=T',  &
                !' topology_suffix1=_EVB1 topology_suffix2=_EVB2'// &
                !' form_bond={1 2} break_bond={2 6} ', &
                pot1=evbsub_empty_qm_pot)
	     call print("INITIALISE_POT EVB")
          else
	     call print("SETUP_POT EMPTY_QM_POT")
             call setup_pot(empty_qm_pot, Run_Type2, filepot_program, tmp_run_dir_i)
          endif
       endif
    endif

    !allocate force lists
    allocate(f0(3,ds%N),f1(3,ds%N),f(3,ds%N))

!FIRST STEP - first step of velocity verlet

    call system_timer('step')

     n = 0

  !FORCE

     if (topology_print_rate >= 0) driver_PSF_Print='DRIVER_PRINT_AND_SAVE'
call print("MAIN CALL CALC EVB")
call print(pot)
     call do_calc_call(pot, empty_qm_pot, ds%atoms, Run_Type1, Run_Type2, EVB, qm_region_pt_ctr, qm_core_region_pt_ctr, &
       distance_ramp, qm_region_ctr, qm_core_region_ctr, cp2k_calc_args, evb_args_str, extra_calc_args, do_carve_cluster, driver_PSF_Print, f1, energy)
call print("MAIN CALLED CALC EVB")
     if (topology_print_rate >= 0) driver_PSF_Print='USE_EXISTING_PSF'

     !spline force calculation, if needed
     if (qm_region_pt_ctr.and.use_spline) then
        allocate(add_force(1:3,1:ds%atoms%N))
	call verbosity_push_decrement()
	  call print('Force due to added spline potential (eV/A):')
	  call print('atom     F(x)     F(y)     F(z)')
          if (.not.(assign_pointer(ds%atoms, "pot", spline_pot_val_p))) &
             call system_abort("couldn't find pot property")
	  do i = 1, ds%atoms%N
	     add_force(1:3,i) = spline_force(ds%atoms,i,my_spline, pot=spline_pot_val)
	     spline_pot_val_p(i) = spline_pot_val
	     call print('  '//i//'    '//round(add_force(1,i),5)//'  '//round(add_force(2,i),5)//'  '//round(add_force(3,i),5))
	  enddo
	call verbosity_pop()
        call print('Sum of the forces: '//sum(add_force(1,1:ds%atoms%N))//' '//sum(add_force(2,1:ds%atoms%N))//' '//sum(add_force(3,1:ds%atoms%N)))
        f = sumBUFFER(f1+add_force,ds%atoms)
        deallocate(add_force)
     else
        f = sum0(f1,ds%atoms)
     endif

     if (ds%Nrestraints > 0) call get_restraint_stuff(ds, restraint_stuff)

  !THERMOSTATTING now - hybrid_mark was updated only in calc
     call set_thermostat_regions(ds%atoms, Thermostat_Type, Thermostat_7_rs, qm_region_ctr)
     call set_thermostat_masses(ds%atoms, Thermostat_Type, Simulation_Temperature, &
	 Inner_QM_Region_Radius, Outer_QM_Region_Radius, Inner_Buffer_Radius, Outer_Buffer_Radius, Thermostat_7_rs)

  !PRINT DS,CONSTRAINT
     call ds_print_status(ds, 'E',energy)
     if (ds%Nrestraints > 0) call print_restraint_stuff(ds, restraint_stuff, 'E')
     call print(ds%thermostat)

    call set_value(ds%atoms%params,'Time',ds%t)

     if (ds%Nconstraints > 0) then
        call print(ds%constraint)
        do i=1,ds%Nconstraints
           call print("TI"//i//" Time[fs]  Constraint_value  Lambda_R  Lambda_V  Z_xi")
           call print("TI"//i//" "//ds%t//" "//ds%constraint(i)%C//" "// &
                      ds%constraint(i)%lambdaR//" "// &
                      ds%constraint(i)%lambdaV//" "// &
                      ds%constraint(i)%Z_coll)
	   if (ds%constraint(i)%N /= 3) then
	      call print("WARNING: constraint " // i // " does not involve 3 atoms, ignoring", PRINT_ALWAYS)
	   else
	      TI_force = force_on_collective_variable(ds%atoms, (/ f(1:3,ds%constraint(i)%atom(1)), &
								   f(1:3,ds%constraint(i)%atom(2)), &
								   f(1:3,ds%constraint(i)%atom(3)) /), ds%constraint(i)%atom(1:3), TI_corr, check)
	      call print('constrained bond length diff: '//round(check,10))
	      call print('force on colvar '//i//' :'//round(TI_force,10)//' '//round(TI_corr,10))
	   endif
        enddo
        call print_cv(ds,Time_Step,.true.)
     endif

  !PRINTING
     !----------------------------------------------------
    if (trim(print_prop).eq.'all') then
        call write(ds%atoms,traj_xyz,real_format='%17.10f', error=error)
	HANDLE_ERROR(error)
    else
        call write(ds%atoms,traj_xyz,properties=trim(print_prop),real_format='%17.10f',error=error)
	HANDLE_ERROR(error)
    endif
    call initialise(latest_xyz,trim(latest_coord_file)//".new",action=OUTPUT)
    call write(ds%atoms,latest_xyz,real_format='%17.10f',error=error)
    HANDLE_ERROR(error)
    call finalise(latest_xyz)
    call system("mv "//trim(latest_coord_file)//".new "//trim(latest_coord_file))
     !----------------------------------------------------

    if (print_forces_at0) then
       do i=1,ds%atoms%N
	  call print('FFF '//f(1,i)//' '//f(2,i)//' '//f(3,i))
       enddo
    endif

  !ADVANCE VERLET 1

     call advance_verlet1(ds, Time_Step)

  !PRINT XYZ

     call set_value(ds%atoms%params,'Time',ds%t)
     if (trim(print_prop).eq.'all') then
         call write(ds%atoms,traj_xyz,real_format='%17.10f',error=error)
	 HANDLE_ERROR(error)
!         call write(ds%atoms,xyz,real_format='%17.10f')
     else
         call write(ds%atoms,traj_xyz,properties=trim(print_prop),real_format='%17.10f', error=error)
	 HANDLE_ERROR(error)
!         call write(ds%atoms,xyz,real_format='%17.10f')
         !call write(ds%atoms,xyz,properties=trim(print_prop),real_format='%17.10f')
     endif
     call initialise(latest_xyz,trim(latest_coord_file)//".new",action=OUTPUT)
     call write(ds%atoms,latest_xyz,real_format='%17.10f', error=error)
     HANDLE_ERROR(error)
     call finalise(latest_xyz)
     call system("mv "//trim(latest_coord_file)//".new "//trim(latest_coord_file))

     !save latest.xyz and wfn files at every $copy_latest_every$ step, starting now
     if (tmp_run_dir_i>0 .and. copy_latest_every>0) then
        copy=1
     !   call system("cp "//trim(latest_coord_file)//" .")
     !   call system("for i in /tmp/cp2k_run_"//tmp_run_dir_i//"/wfn.restart.wfn* ; do cp ${i} . ; done")
     endif

    call system_timer('step')

!LOOP - force calc, then VV-2, VV-1

  ! keep track of how far atoms could have moved
  max_move_since_calc_connect = 0.0_dp
  do while (ds%t < (Equilib_Time + Run_Time) .and. ((max_n_steps < 0) .or. (n < (max_n_steps-1))))

    max_move_since_calc_connect = max_move_since_calc_connect + Time_Step* maxval(abs(ds%atoms%velo))
    ! max distance atom could have moved is about max_move_since_calc_connect
    ! 2.0 time that is max that interatomic distance could have changed
    ! another factor of 1.5 to account for inaccuracy in max_move_since_calc_connect (because atoms really move a distance dependent on both velo and higher derivatives like acc)
    if (max_move_since_calc_connect*2.0_dp*1.5_dp >= calc_connect_buffer) then
      call system_timer("calc_connect")
      call calc_connect(ds%atoms)
      if (Delete_Metal_Connections) call delete_metal_connects(ds%atoms)
      max_move_since_calc_connect = 0.0_dp
      call system_timer("calc_connect")
    else
      call calc_dists(ds%atoms)
    endif

    call system_timer('step')
     n = n + 1

  !QM CORE + BUFFER UPDATE + THERMOSTAT REASSIGNMENT

     if (mod(n,update_QM_region_interval) == 0) then
     !NB should we really recalculate residue labels, so heuristics that keep residues together function here?  probably
       if (have_QMMM_extended) call update_QM_region(ds%atoms, '_extended', qm_region_pt_ctr, qm_region_ctr, &
	 qm_region_atom_ctr, qm_list_filename, Inner_QM_Region_Radius, Outer_QM_Region_Radius, buffer_general, &
	 use_create_cluster_info_for_core, core_create_cluster_info_args, .false., qm_seed_extended)
       if (have_QMMM_core) call update_QM_region(ds%atoms, '_core', qm_core_region_pt_ctr, qm_core_region_ctr, &
	 qm_core_region_atom_ctr, qm_core_list_filename, Inner_QM_Core_Region_Radius, Outer_QM_Core_Region_Radius, buffer_general, &
	 use_create_cluster_info_for_core, core_create_cluster_info_args, .false., qm_seed_core)
     endif

  !FORCE

     if (Topology_Print_rate > 0) then !every #th step
        if (mod(n,Topology_Print_rate) == 0) then    !recalc connectivity & generate PSF (at every n-th step) and then use it
           driver_PSF_Print = 'DRIVER_PRINT_AND_SAVE'
	else
           driver_PSF_Print = 'USE_EXISTING_PSF'
        endif
     endif
     call do_calc_call(pot, empty_qm_pot, ds%atoms, Run_Type1, Run_Type2, EVB, qm_region_pt_ctr, qm_core_region_pt_ctr, &
       distance_ramp, qm_region_ctr, qm_core_region_ctr, cp2k_calc_args, evb_args_str, extra_calc_args, do_carve_cluster, driver_PSF_Print, f1, energy)

    !SPLINE force calculation, if needed
     if (qm_region_pt_ctr.and.use_spline) then
	call verbosity_push_decrement()
	  call print('Force due to added spline potential (eV/A):')
	  call print('atom     F(x)     F(y)     F(z)')
	  allocate(add_force(1:3,1:ds%atoms%N))
          if (.not.(assign_pointer(ds%atoms, "pot", spline_pot_val_p))) &
             call system_abort("couldn't find pot property")
	  do i = 1, ds%atoms%N
	     add_force(1:3,i) = spline_force(ds%atoms,i,my_spline, pot=spline_pot_val)
	     spline_pot_val_p(i) = spline_pot_val
	     call print('  '//i//'    '//round(add_force(1,i),5)//'  '//round(add_force(2,i),5)//'  '//round(add_force(3,i),5))
	  enddo
	call verbosity_pop()
        call print('Sum of the forces: '//sum(add_force(1,1:ds%atoms%N))//' '//sum(add_force(2,1:ds%atoms%N))//' '//sum(add_force(3,1:ds%atoms%N)))
        f = sumBUFFER(f1+add_force,ds%atoms)
        deallocate(add_force)
     else
        f = sum0(f1,ds%atoms)
     endif

     if (print_forces) then
       do i=1,ds%atoms%N
	   call print('FFF '//f(1,i)//' '//f(2,i)//' '//f(3,i))
       enddo
     endif

  !THERMOSTATTING now - hybrid_mark was updated only in calc
       call set_thermostat_regions(ds%atoms, Thermostat_Type, Thermostat_7_rs, qm_region_ctr)
       call set_thermostat_masses(ds%atoms, Thermostat_Type, Simulation_Temperature, &
	 Inner_QM_Region_Radius, Outer_QM_Region_Radius, Inner_Buffer_Radius, Outer_Buffer_Radius, Thermostat_7_rs)

  !ADVANCE VERLET 2

     call advance_verlet2(ds, Time_Step, f)

  !RESTRAINTS
     if (ds%Nrestraints > 0) call get_restraint_stuff(ds, restraint_stuff)

  !PRINT DS,THERMOSTAT,CONSTRAINT,XYZ

     if (ds%t < Equilib_Time) then
        call ds_print_status(ds, 'E',energy)
        if (ds%Nrestraints > 0) call print_restraint_stuff(ds, restraint_stuff, 'E')
     else
        call ds_print_status(ds, 'I',energy)
	if (ds%Nrestraints > 0) call print_restraint_stuff(ds, restraint_stuff, 'I')
     end if

     !Thermostat (only if not per-atom)
     if (thermostat_type /= 3 .and. thermostat_type /= 4) then
       call print(ds%thermostat)
     endif

     !Constraint
     if (ds%Nconstraints > 0) then
        call print(ds%constraint)
        do i=1,ds%Nconstraints
           call print("TI"//i//" "//ds%t//" "//ds%constraint(i)%C//" "// &
                      ds%constraint(i)%lambdaR//" "// &
                      ds%constraint(i)%lambdaV//" "// &
                      ds%constraint(i)%Z_coll)
	   if (ds%constraint(i)%N /= 3) then
	      call print("WARNING: constraint " // i // " does not involve 3 atoms, ignoring", PRINT_ALWAYS)
	   else
	      TI_force = force_on_collective_variable(ds%atoms, (/ f(1:3,ds%constraint(i)%atom(1)), &
								   f(1:3,ds%constraint(i)%atom(2)), &
								   f(1:3,ds%constraint(i)%atom(3)) /), ds%constraint(i)%atom(1:3), TI_corr, check)
	      call print('constrained bond length diff: '//round(check,10))
	      call print('force on colvar '//i//' :'//round(TI_force,10)//' '//round(TI_corr,10))
	   endif
        enddo
        call print_cv(ds,Time_Step,.false.)
     endif

     !XYZ
     if (mod(n,IO_Rate)==0) then
        call set_value(ds%atoms%params,'Time',ds%t)
        if (trim(print_prop).eq.'all') then
            call write(ds%atoms,traj_xyz,real_format='%17.10f',error=error)
	    HANDLE_ERROR(error)
!            call write(ds%atoms,xyz,real_format='%17.10f')
        else
            call write(ds%atoms,traj_xyz,properties=trim(print_prop),real_format='%17.10f',error=error)
	    HANDLE_ERROR(error)
!            call write(ds%atoms,xyz,real_format='%17.10f')
            !call write(ds%atoms,xyz,properties=trim(print_prop),real_format='%17.10f')
        endif
        call initialise(latest_xyz,trim(latest_coord_file)//".new",action=OUTPUT)
        call write(ds%atoms,latest_xyz,real_format='%17.10f',error=error)
	HANDLE_ERROR(error)
        call finalise(latest_xyz)
        call system("mv "//trim(latest_coord_file)//".new "//trim(latest_coord_file))
     end if

     !save latest.xyz and wfn files from /tmp at every $copy_latest_every$ step
     if (tmp_run_dir_i>0 .and. copy_latest_every>0) then
       if (copy==copy_latest_every) then
         copy=0
         call system("cp "//trim(latest_coord_file)//" .")
         call system("for i in /tmp/cp2k_run_"//tmp_run_dir_i//"/wfn.restart.wfn* ; do cp ${i} . ; done")
       endif
       copy=copy+1
     endif

     if (H_extra_heat_r(2) >= H_extra_heat_r(1)) then
       do i=1, ds%atoms%N
	if (ds%atoms%Z(i) == 1) then
	  r = distance_min_image(ds%atoms, i, H_extra_heat_ctr)
	  if (r >= H_extra_heat_r(1) .and. r <= H_extra_heat_r(2)) ds%atoms%velo(:,i) = H_extra_heat_velo_factor*ds%atoms%velo(:,i)
	endif
       end do
     endif
     if (nH_extra_heat_r(2) >= nH_extra_heat_r(1)) then
       do i=1, ds%atoms%N
	if (ds%atoms%Z(i) /= 1) then
	  r = distance_min_image(ds%atoms, i, nH_extra_heat_ctr)
	  if (r >= nH_extra_heat_r(1) .and. r <= nH_extra_heat_r(2)) ds%atoms%velo(:,i) = nH_extra_heat_velo_factor*ds%atoms%velo(:,i)
	endif
       end do
     endif
     
  !ADVANCE VERLET 1

     call advance_verlet1(ds, Time_Step)

     call system_timer('step')

  enddo

  deallocate(f,f0,f1)

  call finalise(qm_seed_extended)
  call finalise(qm_seed_core)
  call finalise(ds)
  call finalise(traj_xyz)

  call finalise(pot)
  call finalise(empty_qm_pot)
  call finalise(CP2K_slow_pot)
  call finalise(CP2K_fast_pot)

  call print_title('THE')
  call print('Finished. CP2K is now having a rest, since deserved it. Bye-Bye!')
  call print_title('END')

  call system_timer('program')
  call system_finalise

contains

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
        F_corr = - 4 * BOLTZMANN_K * 300 * d / (DD*DD-d*d)

     if (present(lambda)) lambda = d

  end function force_on_collective_variable

! for QM/MM and MM runs, to check water topology
  subroutine check_topology(my_atoms)

    type(Atoms), intent(in)      :: my_atoms
  
    integer                                 :: i, N
    logical                                 :: is_bad_H
    integer,                        pointer :: qm_flag_p(:), qm_core_flag_p(:)
    character(len=1), pointer :: atom_res_name_p(:,:)


    if (.not.(assign_pointer(my_atoms, "hybrid_mark_extended", qm_flag_p))) nullify(qm_flag_p)
    if (.not.(assign_pointer(my_atoms, "hybrid_mark_core", qm_core_flag_p))) nullify(qm_core_flag_p)

    if (.not.(assign_pointer(my_atoms, "atom_res_name", atom_res_name_p))) &
       call system_abort("couldn't find atom_res_name property")

     N = 0
     do i=1,my_atoms%N
	is_bad_H = .false.
	if (any((/'H3O','HYD','HWP'/).eq.trim(a2s(atom_res_name_p(:,i))))) then
	  if (associated(qm_flag_p)) then
	    if (qm_flag_p(i) /= 1) is_bad_H = .true.
	  endif
	  if (associated(qm_core_flag_p)) then
	    if (qm_core_flag_p(i) /= 1) is_bad_H = .true.
	  endif
	  if (.not. associated(qm_flag_p) .and. .not. associated(qm_core_flag_p)) is_bad_H = .true.
	endif
	if (is_bad_H) then
	  N = N + 1
	  call print('ERROR: classical or buffer atom '//i//'has atom_res_name '//trim(a2s(atom_res_name_p(:,i))))
	endif
     enddo
     if (N > 0) call system_abort('Wrong topology calculated: found H3O, HYD, or HWP in MM or buffer atom.')

  end subroutine check_topology

! momentum conservation over all atoms, mass weighted
  function sum0(force,at) result(force0)

    real(dp), dimension(:,:), intent(in) :: force
    type(Atoms),              intent(in) :: at
    real(dp) :: force0(size(force,1),size(force,2))
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
  function sumBUFFER(force,at) result(force0)

    real(dp), dimension(:,:), intent(in) :: force
    type(Atoms),              intent(in) :: at
    real(dp) :: force0(size(force,1),size(force,2))
    real(dp) :: F_corr(size(force,1),size(force,2))
    integer :: i
    real(dp) :: sumF(3), sum_weight
    integer, pointer :: qm_flag_p(:)

    if (.not. assign_pointer(at, 'cluster_mark_extended', qm_flag_p)) &
         call system_abort('Potential_FM_Calc: cluster_mark property missing')
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

  !% Calculates the force on the $i$th atom due to an external potential that has
  !% the form of a spline, $my_spline$.
  !
  function spline_force(at, i, my_spline, pot) result(force)

    type(Atoms), intent(in)  :: at
    integer, intent(in) :: i
    type(spline_pot), intent(in) :: my_spline
    real(dp), optional, intent(out) :: pot
    real(dp), dimension(3) :: force

    real(dp) :: dist, factor

! dist = distance from the origin
! spline: f(dist) = spline%dpot/(spline%from-spline%to)**3._dp * (dist-spline%to)**2._dp * (3._dp*spline%from - spline%to - 2._dp*dist)
! f(spline%from) = spline%dpot; f(spline%to) = 0.
! f`(spline%from) = 0.; f`(spline%to) = 0.
! force = - grad f(dist)

    dist = distance_min_image(at,i,(/0._dp,0._dp,0._dp/))
    if (dist.ge.my_spline%to .or. dist.le.my_spline%from) then
       force = 0._dp
       if (present(pot)) then
          if (dist.ge.my_spline%to) pot = 0._dp
          if (dist.le.my_spline%from) pot = my_spline%dpot
       endif
    else
       factor = 1._dp
      ! force on H should be 1/16 times the force on O, to get the same acc.
       if (at%Z(i).eq.1) factor = ElementMass(1)/ElementMass(8)
       force = - factor * 2._dp *my_spline%dpot / ((my_spline%from - my_spline%to)**3._dp * dist) * &
               ( (3._dp * my_spline%from - my_spline%to - 2._dp * dist) * (dist - my_spline%to) &
               - (dist - my_spline%to)**2._dp ) * at%pos(1:3,i)
       if (present(pot)) pot = my_spline%dpot/(my_spline%from-my_spline%to)**3._dp * (dist-my_spline%to)**2._dp * (3._dp*my_spline%from - my_spline%to - 2._dp*dist)
    endif

  end function spline_force

  !% Reads the QM list from a file
  !% Will not overwrite properties!
  !
  subroutine read_qmlist(my_atoms,qmlistfilename,qm_list,verbose)

    type(Atoms),       intent(inout) :: my_atoms
    character(*),      intent(in)    :: qmlistfilename
    type(Table),       intent(out)   :: qm_list
    logical, optional, intent(in)    :: verbose

    type(Inoutput)                   :: qmlistfile
    integer                          :: n,num_qm_atoms,qmatom,status
    character(80)                    :: title,testline
    logical                          :: my_verbose
    character(20), dimension(10)     :: fields
    integer                          :: num_fields


    my_verbose = optional_default(.false., verbose)

    if (my_verbose) call print('In Read_QM_list:')
    call print('Reading the QM list from file '//trim(qmlistfilename)//'...')

    call initialise(qmlistfile,filename=trim(qmlistfilename),action=INPUT)
    title = read_line(qmlistfile,status)
    if (status > 0) then
       call system_abort('read_qmlist: Error reading from '//qmlistfile%filename)
    else if (status < 0) then
       call system_abort('read_qmlist: End of file when reading from '//qmlistfile%filename)
    end if

    call parse_line(qmlistfile,' ',fields,num_fields,status)
    if (status > 0) then
       call system_abort('read_qmlist: Error reading from '//qmlistfile%filename)
    else if (status < 0) then
       call system_abort('read_qmlist: End of file when reading from '//qmlistfile%filename)
    end if

    num_qm_atoms = string_to_int(fields(1))
    if (num_qm_atoms.gt.my_atoms%N) call print('WARNING! read_qmlist: more QM atoms then atoms in the atoms object, possible redundant QM list file',verbosity=PRINT_ALWAYS)
    call print('Number of QM atoms: '//num_qm_atoms)
    call allocate(qm_list,4,0,0,0,num_qm_atoms)      !1 int, 0 reals, 0 str, 0 log, num_qm_atoms entries

    do while (status==0)
       testline = read_line(qmlistfile,status)
       !print *,testline
       if (testline(1:4)=='list') exit
    enddo
   ! Reading and storing QM list...
    do n=1,num_qm_atoms
       call parse_line(qmlistfile,' ',fields,num_fields,status)
       qmatom = string_to_int(fields(1))
       if (my_verbose) call print(n//'th quantum atom is: '//qmatom)
       call append(qm_list,(/qmatom,0,0,0/))
    enddo

    call finalise(qmlistfile)

    if (qm_list%N/=num_qm_atoms) call system_abort('read_qmlist: Something wrong with the QM list file')
    if (any(int_part(qm_list,1).gt.my_atoms%N).or.any(int_part(qm_list,1).lt.1)) &
       call system_abort('read_qmlist: at least 1 QM atom is out of range')
    if ((size(int_part(qm_list,1)).gt.my_atoms%N).or.(size(int_part(qm_list,1)).lt.1)) &
       call system_abort("read_qmlist: QM atoms' number is <1 or >"//my_atoms%N)

    if (my_verbose) call print('Finished. '//qm_list%N//' QM atoms have been read successfully.')

  end subroutine read_qmlist

  !%Prints the quantum region mark with the cluster_mark property (/=0).
  !%If file is given, into file, otherwise to the standard io.
  !
  subroutine print_qm_region(at, file)
    type(Atoms), intent(inout) :: at
    type(CInoutput), optional :: file

    integer, pointer :: cluster_mark_p(:)
    type(Atoms) :: cluster

    if (.not.(assign_pointer(at, "cluster_mark_extended", cluster_mark_p))) &
      call system_abort("print_qm_region couldn't find cluster_mark property")

    call select(cluster, at, mask=(cluster_mark_p /= 0))
    if (present(file)) then
      call write(cluster, file, prefix="QM_REGION")
    else
      call write(cluster, "stdout", prefix="QM_REGION")
    end if
  end subroutine print_qm_region

  subroutine do_calc_call(pot, empty_qm_pot, at, Run_Type1, Run_Type2, EVB_MM, qm_region_pt_ctr, qm_core_region_pt_ctr, &
			  distance_ramp, qm_region_ctr, qm_core_region_ctr, cp2k_calc_args, evb_args_str, extra_calc_args, &
			  do_carve_cluster, driver_PSF_Print, f1, energy)
     type(Potential), intent(inout) :: pot, empty_qm_pot
     type(Atoms), intent(inout) :: at
     logical, intent(in) :: EVB_MM, qm_region_pt_ctr, qm_core_region_pt_ctr, distance_ramp, do_carve_cluster
     real(dp), intent(in) :: qm_region_ctr(3), qm_core_region_ctr(3)
     character(len=*), intent(in) :: Run_Type1, Run_Type2, cp2k_calc_args, evb_args_str, extra_calc_args, driver_PSF_Print
     real(dp), intent(inout) :: f1(:,:)
     real(dp), intent(out) :: energy

     integer, pointer :: qm_flag_p(:)
     character(len=STRING_LENGTH)        :: slow_args_str, fast_args_str, args_str
     logical :: empty_QM_core, empty_QM_extended
     real(dp), pointer :: force_p(:,:)

     empty_QM_core = .false.
     if (qm_core_region_pt_ctr) then
        if (.not.(assign_pointer(at, "hybrid_mark_core", qm_flag_p))) &
           call system_abort("couldn't find hybrid_mark_core property")
        if (.not.any(qm_flag_p(1:at%N).eq.HYBRID_ACTIVE_MARK)) empty_QM_core = .true.
	if (empty_QM_core .and. trim(Run_Type2) == 'QMMM_CORE') &
	  call system_abort("Can't handle Run_Type2=QMMM_CORE but QM core appears empty")
     endif
     empty_QM_extended = .false.
     if (qm_region_pt_ctr) then
        if (.not.(assign_pointer(at, "hybrid_mark_extended", qm_flag_p))) &
           call system_abort("couldn't find hybrid_mark_extended property")
        if (.not.any(qm_flag_p(1:at%N).eq.HYBRID_ACTIVE_MARK)) empty_QM_extended = .true.
     endif
     if (.not.((trim(Run_Type1).eq.'QS').or.(trim(Run_Type1).eq.'MM'))) then
        call print_qm_region(at)
     endif

     if (trim(Run_Type2) == 'NONE' .or. (qm_region_pt_ctr .and. empty_QM_extended)) then ! no force mixing
        call print(trim(Run_Type1)//' run will be performed with simple potential.')
        if (EVB_MM .and. &
           !only instead of MM
           ((trim(Run_Type1) == "MM" .and. trim(Run_Type2) == 'NONE') .or. &
            (Run_Type1(1:4) == 'QMMM' .and. (qm_region_pt_ctr .and. empty_QM_extended) .and. trim(Run_Type2)=="MM"))) then
           !Add args_str into mm_args_str of EVB pot
           args_str=trim(evb_args_str)//' mm_args_str="'// &
             trim(cp2k_calc_args) // &
             ' Run_Type=MM'// & !only MM
             ' PSF_Print='//trim(driver_PSF_Print)//'"'
        else
           args_str=trim(cp2k_calc_args) // &
	     ' '//trim(cp2k_driver_run_type_args(Run_Type1))// &
             ' PSF_Print='//trim(driver_PSF_Print) !// &
!	     ' clean_up_files=F'
        endif
        call print('ARGS_STR | '//trim(args_str))
	if (Run_Type1(1:4) == 'QMMM') then
	  if ( qm_region_pt_ctr .and. empty_QM_extended) then
	    call print('WARNING: Empty QM extended region. MM run will be performed instead of QM/MM.', PRINT_ALWAYS)
	    call print("CALC_POT '"//trim(args_str)//"'")
	    call calc(empty_qm_pot,at,energy=energy,force=f1,args_str=trim(args_str))
	    if (assign_pointer(at, 'force', force_p)) force_p = f1
	    return
	  else
	    args_str = trim(args_str) // &
	      ' single_cluster=T carve_cluster='//do_carve_cluster//' cluster_nneighb_only=T ' // &
	      ' termination_clash_check=T terminate=T even_electrons=F auto_centre'
	  endif
          !cluster_mark_postfix to save cluster_mark & old_cluster_mark under different name for QMMM_extended & QMMM_core
          if (trim(Run_Type1) == 'QMMM_EXTENDED') then
            args_str = trim(args_str) // ' hybrid_mark_postfix=_extended'
          elseif (trim(Run_Type1) == 'QMMM_CORE') then
            args_str = trim(args_str) // ' hybrid_mark_postfix=_core'
          endif
	endif
	if (Run_Type1(1:4) == 'QMMM' .and. trim(Run_Type2)=="NONE") then
          !potential simple will not calculate energy with single_cluster=T
	  call print("CALC_POT '"//trim(args_str)//"'")
	  call calc(pot,at,force=f1,args_str=trim(args_str))
          energy=0.0_dp
	else
	  call print("CALC_POT '"//trim(args_str)//"'")
	  call calc(pot,at,energy=energy,force=f1,args_str=trim(args_str))
        endif
     else ! do force mixing

       if (EVB_MM .and. trim(Run_Type1)=="MM") then
          call print("WARNING: Force Mixing with MM as the slow potential!")
          !Add slow_args_str into mm_args_str of EVB pot
          slow_args_str=trim(evb_args_str)//' mm_args_str="'// &
             trim(cp2k_calc_args) // ' ' // trim(cp2k_driver_run_type_args(Run_Type1)) // ' PSF_Print='//trim(driver_PSF_print)//'"'
       else
          slow_args_str=trim(cp2k_calc_args) // ' '//trim(cp2k_driver_run_type_args(Run_Type1)) //' PSF_Print='//trim(driver_PSF_print) !//' clean_up_files=F'
       endif
       if (Run_Type1(1:4) == 'QMMM' .and. .not. (qm_region_pt_ctr .and. empty_QM_core)) then
	 slow_args_str = trim(slow_args_str) // &
           ' single_cluster=T carve_cluster='//do_carve_cluster//' cluster_nneighb_only=T ' // &
	   ' termination_clash_check=T terminate=T even_electrons=F auto_centre'
         !cluster_mark_postfix to save cluster_mark & old_cluster_mark under different name for QMMM_extended & QMMM_core
         if (trim(Run_Type1) == 'QMMM_EXTENDED') then
           slow_args_str = trim(slow_args_str) // ' cluster_mark_postfix=_extended hybrid_mark_postfix=_extended'
         elseif (trim(Run_Type1) == 'QMMM_CORE') then
           slow_args_str = trim(slow_args_str) // ' cluster_mark_postfix=_core hybrid_mark_postfix=_core'
         endif
       endif

       if (EVB_MM .and. trim(Run_Type2)=="MM") then
          !Add fast_args_str into mm_args_str of EVB pot
          fast_args_str=trim(evb_args_str)//' mm_args_str="'// &
             trim(cp2k_calc_args) // &
             ' Run_Type=MM'// & !only MM
             ' PSF_Print='//trim(driver_PSF_Print)//'"'
       else
          fast_args_str=trim(cp2k_calc_args) // ' '// trim(cp2k_driver_run_type_args(Run_Type2))//' PSF_Print='//trim(driver_PSF_print) !//' clean_up_files=F'
       endif
       if (Run_Type2(1:4) == 'QMMM' .and. .not. (qm_region_pt_ctr .and. empty_QM_core)) then
	 fast_args_str = trim(fast_args_str) // &
           ' single_cluster=T carve_cluster='//do_carve_cluster//' cluster_nneighb_only=T ' // &
	   ' termination_clash_check=T terminate=T even_electrons=F auto_centre'
         !cluster_mark_postfix to save cluster_mark & old_cluster_mark under different name for QMMM_extended & QMMM_core
         if (trim(Run_Type2) == 'QMMM_EXTENDED') then
           fast_args_str = trim(fast_args_str) // ' hybrid_mark_postfix=_extended'
         elseif (trim(Run_Type2) == 'QMMM_CORE') then
           fast_args_str = trim(fast_args_str) // ' hybrid_mark_postfix=_core'
         endif
       endif

       args_str='qm_args_str={'//trim(slow_args_str)//' '//trim(extra_calc_args)//'} mm_args_str={'//trim(fast_args_str)//' '//trim(extra_calc_args)//'} hybrid_mark_postfix=_extended'
       if (distance_ramp) then
	 args_str = trim(args_str) // ' distance_ramp_centre='//qm_region_ctr
       endif
       call print('ARGS_STR | '//trim(args_str))
       if (qm_region_pt_ctr .and. empty_QM_extended) then
	 call system_abort("Got to end of force mixing part of do_calc_call, but (qm_region_ptr_ctr .and. empty_QM_extended) is true.  Should never happen!")
!	 if (trim(Run_Type2) /= 'MM') &
!	   call system_abort("Doing force mixing, but Run_Type2='"//trim(Run_Type2)//"' /= MM")
!	 call calc(empty_qm_pot,at,force=f1,args_str=trim(fast_args_str))
       else
	 call print("CALC_POT '"//trim(args_str)//"'")
	 call calc(pot,at,force=f1,args_str=trim(args_str))
       endif
       energy=0._dp !no energy
     endif ! do force mixing

     ! save forces if force property is there
     if (assign_pointer(at, 'force', force_p)) force_p = f1

  end subroutine do_calc_call

  function cp2k_driver_run_type_args(run_type)
    character(len=*), intent(in) :: run_type
    character(len=300) :: cp2k_driver_run_type_args

    if (run_type(1:2) == 'MM') then
      cp2k_driver_run_type_args = "Run_Type=MM"
    else if (run_type(1:2) == 'QS') then
      cp2k_driver_run_type_args = "Run_Type=QS"
    else ! QMMM
      cp2k_driver_run_type_args = "Run_Type=QMMM"
      if (run_type(1:9) == 'QMMM_CORE') then
	cp2k_driver_run_type_args = trim(cp2k_driver_run_type_args)//' use_buffer=F qm_name_postfix=_core'
      else
	cp2k_driver_run_type_args = trim(cp2k_driver_run_type_args)//' use_buffer=T qm_name_postfix=_extended'
      endif
    endif
  end function cp2k_driver_run_type_args

  subroutine setup_pot(pot, Run_Type, filepot_program, tmp_run_dir_i)
    type(Potential), intent(inout) :: pot
    character(len=*), intent(in) :: Run_Type, filepot_program
    integer, intent(in) :: tmp_run_dir_i
    character(len=FIELD_LENGTH) :: filename

    if (tmp_run_dir_i > 0) then
       filename="/tmp/cp2k_run_"//tmp_run_dir_i//"/filepot"
    else
       filename="filepot"
    endif
    if (trim(Run_Type) == 'QS') then
       call print('INITIALISE_POT QS FilePot filename='//trim(filename)//' command='//trim(filepot_program)//' property_list=species:pos min_cutoff=0.0')
       call initialise(pot,'FilePot filename='//trim(filename)//' command='//trim(filepot_program)//' property_list=species:pos min_cutoff=0.0')
    else if (trim(Run_Type) == 'MM') then
       call print('INITIALISE_POT MM FilePot filename='//trim(filename)//' command='//trim(filepot_program)//' property_list=species:pos:avgpos:mol_id:atom_res_number min_cutoff=0.0')
       call initialise(pot,'FilePot filename='//trim(filename)//' command='//trim(filepot_program)//' property_list=species:pos:avgpos:mol_id:atom_res_number min_cutoff=0.0')
    else if (trim(Run_Type) == 'QMMM_CORE') then
       call print('INITIALISE_POT QMMM_CORE FilePot filename='//trim(filename)//' command='//trim(filepot_program)//' property_list=species:pos:avgpos:atom_charge:mol_id:atom_res_number:cluster_mark_core:old_cluster_mark_core:cut_bonds_core min_cutoff=0.0')
       call initialise(pot,'FilePot filename='//trim(filename)//' command='//trim(filepot_program)//' property_list=species:pos:avgpos:atom_charge:mol_id:atom_res_number:cluster_mark_core:old_cluster_mark_core:cut_bonds_core min_cutoff=0.0')
    else if (trim(Run_Type) == 'QMMM_EXTENDED') then
       call print('INITIALISE_POT QMMM_EXTENDED FilePot filename='//trim(filename)//' command='//trim(filepot_program)//' property_list=species:pos:avgpos:atom_charge:mol_id:atom_res_number:cluster_mark_extended:old_cluster_mark_extended:cut_bonds_extended min_cutoff=0.0')
       call initialise(pot,'FilePot filename='//trim(filename)//' command='//trim(filepot_program)//' property_list=species:pos:avgpos:atom_charge:mol_id:atom_res_number:cluster_mark_extended:old_cluster_mark_extended:cut_bonds_extended min_cutoff=0.0')
    else
       call system_abort("Run_Type='"//trim(Run_Type)//"' not supported")
    endif
  end subroutine setup_pot

  subroutine add_thermostats(ds, thermostat_type, T)
    type(DynamicalSystem), intent(inout) :: ds
    integer, intent(in) :: thermostat_type
    real(dp), intent(in) :: T

    select case(Thermostat_Type)
      case(0)
	call print("No thermostat, NVE!!", PRINT_ALWAYS)
      case(1)
	call add_thermostat(ds,type=LANGEVIN,T=T,tau=Langevin_Tau)
	call print('Added single Langevin Thermostat')
      case(2)
	call add_thermostat(ds,type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp)
	call print('Added single Nose-Hoover Thermostat')
      case(3)
	do i=1, ds%atoms%N
	  call add_thermostat(ds,type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp)
	end do
	ds%print_thermostat_temps = .false.
	call print("Added 1 Nose-Hoover thermostat for each atom")
      case(4)
	call add_thermostat(ds,type=NOSE_HOOVER_LANGEVIN,T=T,Q=1.0_dp, tau=Langevin_Tau, n=ds%atoms%N)
	! do i=1, ds%atoms%N
	  ! call add_thermostat(ds,type=NOSE_HOOVER_LANGEVIN,T=T,Q=1.0_dp, tau=Langevin_Tau)
	! end do
	ds%print_thermostat_temps = .false.
	call print("Added 1 Nose-Hoover-Langevin thermostat for each atom")
      case(5)
	call add_thermostat(ds, type=NOSE_HOOVER,T=T, Q=1.0_dp, gamma=0.0_dp) ! heavy QM+buffer
	call add_thermostat(ds, type=NOSE_HOOVER,T=T, Q=1.0_dp, gamma=0.0_dp) ! H QM+buffer
	call add_thermostat(ds, type=NOSE_HOOVER_LANGEVIN,T=T,tau=Langevin_Tau,Q=1.0_dp) ! MM
	call print("Added 1 Nose-Hoover for QM+buffer heavy, 1 Nose-Hoover for QM+buffer H, and 1 Nose-Hoover-Langevin for MM")
      case(6)
	call add_thermostat(ds, type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp) ! heavy QM
	call add_thermostat(ds, type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp) ! H QM
	call add_thermostat(ds, type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp) ! heavy buffer
	call add_thermostat(ds, type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp) ! H buffer
	call add_thermostat(ds, type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp) ! heavy MM
	call add_thermostat(ds, type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp) ! H MM
	call print("Added 6 Nose-Hoover thermostats, 3 regions (QM, Buffer, MM) x 2 kinds (H, heavy)")
      case(7)
	call add_thermostat(ds, type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp) ! heavy QM
	call add_thermostat(ds, type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp) ! H QM
	call add_thermostat(ds, type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp) ! heavy buffer
	call add_thermostat(ds, type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp) ! H buffer
	call add_thermostat(ds, type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp) ! heavy MM
	call add_thermostat(ds, type=NOSE_HOOVER,T=T,Q=1.0_dp, gamma=0.0_dp) ! H MM
	call print("Added 6 Nose-Hoover thermostats, 3 regions (QM, Buffer, MM) x 2 kinds (H, heavy)")
      case default
	call system_abort("add_thermostats: Unknown Thermostat_Type="//Thermostat_type)
    end select
  end subroutine add_thermostats

  subroutine set_thermostat_masses(at, Thermostat_type, T, Inner_QM_Region_Radius, Outer_QM_Region_Radius, Inner_Buffer_Radius, Outer_Buffer_Radius, Thermostat_7_rs)
    type(Atoms), intent(in) :: at
    integer, intent(in) :: thermostat_type
    real(dp) :: T
    real(dp), intent(in) :: Inner_QM_Region_Radius, Outer_QM_Region_Radius, Inner_Buffer_Radius, Outer_Buffer_Radius
    real(dp), intent(in) :: Thermostat_7_rs(2)

    integer, pointer :: cluster_mark_p(:)

    real(dp) :: QM_vol, Buffer_vol, cell_vol
    integer :: Ndof_QM_H, Ndof_QM_heavy
    integer :: Ndof_Buffer_H, Ndof_Buffer_heavy
    integer :: Ndof_MM_H, Ndof_MM_heavy
    real(dp) :: n_H, n_heavy

    if (Thermostat_type == 7) then
      QM_vol = 4.0_dp/3.0_dp*PI*Thermostat_7_rs(1)**3
      Buffer_vol = 4.0_dp/3.0_dp*PI*Thermostat_7_rs(2)**3 - QM_vol
      cell_vol = cell_volume(at)
      n_H = count(at%Z == 1)
      n_heavy = at%N - n_H
      Ndof_QM_H = 3*int(QM_vol * n_H/cell_vol)
      Ndof_QM_heavy = 3*int(QM_vol * n_heavy/cell_vol)
      Ndof_Buffer_H = 3*int(Buffer_vol * n_H/cell_vol)
      Ndof_Buffer_heavy = 3*int(Buffer_vol * n_heavy/cell_vol)
      Ndof_MM_H = 3*n_H - ndof_QM_H - ndof_Buffer_H
      Ndof_MM_heavy = 3*n_heavy - ndof_QM_heavy - ndof_Buffer_heavy
      call print("Thermostat_type == 7, estimating NDOFs from volumes QM,Buffer,cell " // QM_vol // " " // &
	Buffer_vol//" "//cell_vol // " density H,heavy" // (n_H/cell_vol) // " " // (n_heavy/cell_vol))
    else
      if (.not. assign_pointer(at, "cluster_mark_extended", cluster_mark_p)) then
	QM_vol = 4.0_dp/3.0_dp*PI*((Inner_QM_Region_Radius + Outer_QM_Region_Radius)/2.0_dp)**3
	Buffer_vol = 4.0_dp/3.0_dp*PI*((Inner_QM_Region_Radius+Inner_Buffer_Radius + Outer_QM_Region_Radius+Outer_Buffer_Radius)/2.0_dp)**3 - QM_vol
	cell_vol = cell_volume(at)
	n_H = count(at%Z == 1)
	n_heavy = at%N - n_H
	Ndof_QM_H = 3*int(QM_vol * n_H/cell_vol)
	Ndof_QM_heavy = 3*int(QM_vol * n_heavy/cell_vol)
	Ndof_Buffer_H = 3*int(Buffer_vol * n_H/cell_vol)
	Ndof_Buffer_heavy = 3*int(Buffer_vol * n_heavy/cell_vol)
	Ndof_MM_H = 3*n_H - ndof_QM_H - ndof_Buffer_H
	Ndof_MM_heavy = 3*n_heavy - ndof_QM_heavy - ndof_Buffer_heavy
	call print("no cluster_mark_extended, estimating NDOFs from volumes QM,Buffer,cell " // QM_vol // " " // &
	  Buffer_vol//" "//cell_vol // " density H,heavy" // (n_H/cell_vol) // " " // (n_heavy/cell_vol))
      else
	Ndof_QM_H = 3*count(cluster_mark_p == HYBRID_ACTIVE_MARK .and. at%Z == 1)
	Ndof_QM_heavy = 3*count(cluster_mark_p == HYBRID_ACTIVE_MARK) - Ndof_QM_H
	Ndof_Buffer_H = 3*count(cluster_mark_p == HYBRID_BUFFER_MARK .and. at%Z == 1)
	Ndof_Buffer_heavy = 3*count(cluster_mark_p == HYBRID_BUFFER_MARK) - Ndof_Buffer_H
	n_H = count(at%Z == 1)
	n_heavy = at%N - n_H
	Ndof_MM_H = 3*n_H - (Ndof_QM_H+Ndof_Buffer_H)
	Ndof_MM_heavy = 3*n_heavy - (Ndof_QM_heavy+Ndof_Buffer_heavy)
      endif
    endif

    call print("set_thermostat_masses got NDOFs (H, heavy) QM " // Ndof_QM_H // " " // Ndof_QM_heavy // &
      " Buffer " // Ndof_Buffer_H // " " // Ndof_Buffer_heavy // " MM " // Ndof_MM_H // " " // Ndof_MM_heavy)

    select case(Thermostat_Type)
      case (0, 1)
	continue
      case (2)
	ds%thermostat(1)%Q = nose_hoover_mass(Ndof=3*at%N, T=T, tau=Nose_Hoover_tau)
      case(3, 4)
	do i=1, ds%atoms%N
	  ds%thermostat(i)%Q = nose_hoover_mass(Ndof=3, T=T, tau=Nose_Hoover_tau)
	end do
      case(5)
	ds%thermostat(1)%Q = nose_hoover_mass(Ndof=Ndof_QM_heavy+Ndof_Buffer_heavy, T=T, tau=Nose_Hoover_tau)
	ds%thermostat(2)%Q = nose_hoover_mass(Ndof=Ndof_QM_H+Ndof_Buffer_H, T=T, tau=Nose_Hoover_tau)
	ds%thermostat(3)%Q = nose_hoover_mass(Ndof=Ndof_MM_H+Ndof_MM_heavy, T=T, tau=Nose_Hoover_tau)
      case(6)
	ds%thermostat(1)%Q = nose_hoover_mass(Ndof=Ndof_QM_heavy, T=T, tau=Nose_Hoover_tau)
	ds%thermostat(2)%Q = nose_hoover_mass(Ndof=Ndof_QM_H, T=T, tau=Nose_Hoover_tau)
	ds%thermostat(3)%Q = nose_hoover_mass(Ndof=Ndof_Buffer_heavy, T=T, tau=Nose_Hoover_tau)
	ds%thermostat(4)%Q = nose_hoover_mass(Ndof=Ndof_Buffer_H, T=T, tau=Nose_Hoover_tau)
	ds%thermostat(5)%Q = nose_hoover_mass(Ndof=Ndof_MM_heavy, T=T, tau=Nose_Hoover_tau)
	ds%thermostat(6)%Q = nose_hoover_mass(Ndof=Ndof_MM_H, T=T, tau=Nose_Hoover_tau)
      case(7)
	ds%thermostat(1)%Q = nose_hoover_mass(Ndof=Ndof_QM_heavy, T=T, tau=Nose_Hoover_tau)
	ds%thermostat(2)%Q = nose_hoover_mass(Ndof=Ndof_QM_H, T=T, tau=Nose_Hoover_tau)
	ds%thermostat(3)%Q = nose_hoover_mass(Ndof=Ndof_Buffer_heavy, T=T, tau=Nose_Hoover_tau)
	ds%thermostat(4)%Q = nose_hoover_mass(Ndof=Ndof_Buffer_H, T=T, tau=Nose_Hoover_tau)
	ds%thermostat(5)%Q = nose_hoover_mass(Ndof=Ndof_MM_heavy, T=T, tau=Nose_Hoover_tau)
	ds%thermostat(6)%Q = nose_hoover_mass(Ndof=Ndof_MM_H, T=T, tau=Nose_Hoover_tau)
      case default
	call system_abort("set_thermostat_masses: Unknown thermostat_type="//thermostat_type//" in set_thermostat_masses")
    end select
  end subroutine set_thermostat_masses

  subroutine set_thermostat_regions(at, thermostat_type, Thermostat_7_rs, Thermostat_7_centre)
    type(Atoms), intent(inout) :: at
    integer, intent(in) :: thermostat_type
    real(dp), intent(in) :: thermostat_7_rs(2), Thermostat_7_centre(3)

    integer, pointer :: cluster_mark_p(:)
    real(dp) :: r

    select case(Thermostat_Type)
      case(0) 
	  continue
      case(1, 2)
	  at%thermostat_region = 1
      case(3, 4)
	  do i=1, ds%atoms%N
	    at%thermostat_region(i) = i
	  end do
      case (5)
	if (.not. assign_pointer(at, 'cluster_mark_extended', cluster_mark_p)) &
	  call system_abort("set_thermostat_region failed to find cluster_mark_extended for thermostat_type="//thermostat_type)
	at%thermostat_region = 3 ! MM
	where ((cluster_mark_p /= HYBRID_NO_MARK .and. cluster_mark_p /= HYBRID_TERM_MARK) .and. at%Z /= 1) ! QM+Buffer heave
	  at%thermostat_region = 1
	end where
	where ((cluster_mark_p /= HYBRID_NO_MARK .and. cluster_mark_p /= HYBRID_TERM_MARK) .and. at%Z == 1) ! QM+Buffer H
	  at%thermostat_region = 2
	end where
      case(6)
	if (.not. assign_pointer(at, 'cluster_mark_extended', cluster_mark_p)) &
	  call system_abort("set_thermostat_region failed to find cluster_mark for thermostat_type="//thermostat_type)
	where ((cluster_mark_p == HYBRID_ACTIVE_MARK) .and. at%Z /= 1) ! QM heavy
	  at%thermostat_region = 1
	end where
	where ((cluster_mark_p == HYBRID_ACTIVE_MARK) .and. at%Z == 1) ! QM H
	  at%thermostat_region = 2
	end where
	where ((cluster_mark_p == HYBRID_BUFFER_MARK) .and. at%Z /= 1) ! Buffer heavy
	  at%thermostat_region = 3
	end where
	where ((cluster_mark_p == HYBRID_BUFFER_MARK) .and. at%Z == 1) ! Buffer H
	  at%thermostat_region = 4
	end where
	where ((cluster_mark_p /= HYBRID_ACTIVE_MARK .and. cluster_mark_p /= HYBRID_BUFFER_MARK) .and. at%Z /= 1) ! MM heavy
	  at%thermostat_region = 5
	end where
	where ((cluster_mark_p /= HYBRID_ACTIVE_MARK .and. cluster_mark_p /= HYBRID_BUFFER_MARK) .and. at%Z == 1) ! MM H
	  at%thermostat_region = 6
	end where
      case (7)
	call print("Thermostat=7, center at " // thermostat_7_centre, PRINT_ALWAYS)
	do i=1, at%N
	  r = distance_min_image(at, i, thermostat_7_centre)
	  if (r < Thermostat_7_rs(1)) then
	    if (at%Z(i) /= 1) then
	      at%thermostat_region(i) = 1
	    else
	      at%thermostat_region(i) = 2
	    endif
	  else if (r < Thermostat_7_rs(2)) then
	    if (at%Z(i) /= 1) then
	      at%thermostat_region(i) = 3
	    else
	      at%thermostat_region(i) = 4
	    endif
	  else
	    if (at%Z(i) /= 1) then
	      at%thermostat_region(i) = 5
	    else
	      at%thermostat_region(i) = 6
	    endif
	  endif
	end do
      case default
	call system_abort("set_thermostat_regions: Unknown thermostat_type="//thermostat_type//" in set_thermostat_masses")
    end select
  end subroutine set_thermostat_regions

   subroutine print_restraint_stuff(ds, restraint_stuff, suffix)
     type(DynamicalSystem), intent(in) :: ds
     real(dp), intent(in) :: restraint_stuff(:,:)
     character(len=*), intent(in) :: suffix
     logical, save :: firstcall = .true.

     if (firstcall) then
       call print("#R t  target_v C E dE/dcoll dE/dk ....")
     endif

     call print('R'//trim(suffix) // " " // ds%t // " " // reshape( restraint_stuff, (/ size(restraint_stuff) /) ))
   end subroutine print_restraint_stuff

   subroutine get_restraint_stuff(ds, restraint_stuff)
     type(DynamicalSystem), intent(in) :: ds
     real(dp), intent(out) :: restraint_stuff(:,:)
 
     integer i_r
 
     do i_r = 1, ds%Nrestraints
        restraint_stuff(1,i_r) = ds%restraint(i_r)%target_v
        restraint_stuff(2,i_r) = ds%restraint(i_r)%C
        restraint_stuff(3,i_r) = ds%restraint(i_r)%E
        restraint_stuff(4,i_r) = -ds%restraint(i_r)%dE_dcoll
        restraint_stuff(5,i_r) = -ds%restraint(i_r)%dE_dk
     end do
   end subroutine get_restraint_stuff

   !% print coord, velo and force of constraint atoms in a cv file
   subroutine print_cv(ds,time_step,first)
     type(dynamicalsystem), intent(in) :: ds
     real(dp), intent(in) :: time_step
     logical, intent(in) :: first
     real(dp), pointer :: v(:,:),f(:,:)
     character(len=FIELD_LENGTH) :: output_line
     integer :: i,j
     nullify(v)
     if (.not. assign_pointer(ds%atoms, 'velo', v)) call system_abort("Could not find velo property.")
     nullify(f)
     if (.not. assign_pointer(ds%atoms, 'frc', f)) then
        if (.not. assign_pointer(ds%atoms, 'force', f)) call system_abort("Could not find frc or force property.")
     endif
     do j=1,ds%Nconstraints
        if (first) call print("CV"//j//"| "//"# Step Nr.   Time[fs]  Atom Nr.  Pos_x[A]     Pos_y[A]     Pos_z[A]    Vel_x[A/fs]  Vel_y[A/fs]  Vel_z[A/fs]  Frc_x[eV/A]  Frc_y[eV/A]  Frc_z[eV/A]")
        do i=1,ds%constraint(j)%N
           write (output_line,FMT="(I8,F13.3,I8,9F13.8)") &
               nint(ds%t/time_step), ds%t, ds%constraint(j)%atom(i), &
               ds%atoms%pos(1:3,ds%constraint(j)%atom(i)), &
               (v(1:3,ds%constraint(j)%atom(i))), &
               (f(1:3,ds%constraint(j)%atom(i)))
           call print("CV"//j//"| "//trim(output_line))
        enddo
     enddo
   end subroutine print_cv

  subroutine update_QM_region(ds_atoms, mark_postfix, qm_region_pt_ctr, qm_region_ctr, qm_region_atom_ctr, qm_list_filename, inner_radius, outer_radius, buffer_general, use_create_cluster_info_for_core, create_cluster_info_args, first_time, qm_seed)
    type(Atoms), intent(inout) :: ds_atoms
    character(len=*), intent(in) :: mark_postfix
    logical, intent(in) :: qm_region_pt_ctr
    real(dp), intent(inout) :: qm_region_ctr(3)
    integer, intent(in) :: qm_region_atom_ctr
    character(len=*), intent(in) :: qm_list_filename
    real(dp), intent(in) :: inner_radius, outer_radius
    logical, intent(in) :: buffer_general, use_create_cluster_info_for_core
    character(len=*), intent(in) :: create_cluster_info_args
    logical, intent(in) :: first_time
    type(table), intent(inout) :: qm_seed
    integer, pointer :: hybrid_p(:), hybrid_mark_p(:), cluster_mark_p(:)

    logical :: list_changed1

    if (first_time) then
      call add_property(ds_atoms,'hybrid'//trim(mark_postfix),HYBRID_NO_MARK)
      call add_property(ds_atoms,'hybrid_mark'//trim(mark_postfix),HYBRID_NO_MARK)
      call add_property(ds_atoms,'old_hybrid_mark'//trim(mark_postfix),HYBRID_NO_MARK)
      call add_property(ds_atoms,'cluster_mark'//trim(mark_postfix),HYBRID_NO_MARK)
      call add_property(ds_atoms,'old_cluster_mark'//trim(mark_postfix),HYBRID_NO_MARK)
      call add_property(ds_atoms,'cut_bonds'//trim(mark_postfix), 0, n_cols=4) !MAX_CUT_BONDS)
    endif

    if (qm_region_pt_ctr) then
       if (first_time) then
	 call map_into_cell(ds_atoms)
	 call calc_dists(ds_atoms)
       endif
       if (qm_region_atom_ctr /= 0) qm_region_ctr = ds_atoms%pos(:,qm_region_atom_ctr)
       !update hybrid_*, hybrid_mark_*, old_hybrid_mark_*
       call create_pos_or_list_centred_hybrid_region(ds_atoms,inner_radius,outer_radius, &
	 origin=qm_region_ctr,add_only_heavy_atoms=(.not. buffer_general), &
	 use_create_cluster_info=use_create_cluster_info_for_core,create_cluster_info_args=create_cluster_info_args,&
	 list_changed=list_changed1, mark_postfix=trim(mark_postfix))
       if (.not.(assign_pointer(ds_atoms, "hybrid_mark"//trim(mark_postfix), hybrid_mark_p))) call system_abort('??')
       if (list_changed1) then
	  call print('Region with postfix "'//trim(mark_postfix)//'" has changed')
       endif
    else ! not qm_region_pt_ctr
       !QM region centred around an atom list

       if (first_time) then !reinitialise QM region

         !read QM seed
	 call read_qmlist(ds_atoms,qm_list_filename,qm_seed)

         if (.not. (assign_pointer(ds_atoms, 'hybrid'//trim(mark_postfix), hybrid_p))) call system_abort("??")
	 if (.not.(assign_pointer(ds_atoms, "hybrid_mark"//trim(mark_postfix), hybrid_mark_p))) call system_abort('??')
	 if (.not.(assign_pointer(ds_atoms, "cluster_mark"//trim(mark_postfix), cluster_mark_p))) call system_abort('??')

          !call print("MARKSX0 "//count(hybrid_p(1:ds_atoms%N)==1)//" + "//count(hybrid_p(1:ds_atoms%N)==2)//" hybrid atoms",PRINT_ANAL)
          !call print("MARKSX0 "//count(hybrid_mark_p(1:ds_atoms%N)==1)//" + "//count(hybrid_mark_p(1:ds_atoms%N)==2)//" hybrid_mark atoms",PRINT_ANAL)
          !call print("MARKSX0 "//count(cluster_mark_p(1:ds_atoms%N)==1)//" + "//count(cluster_mark_p(1:ds_atoms%N)==2)//" cluster_mark atoms",PRINT_ANAL)

         !(re)initialise the hysteretic qm region with (re)initialising hybrid_*, hybrid_mark_* and cluster_mark_*
         hybrid_p(1:ds_atoms%N) = HYBRID_NO_MARK
         hybrid_p(int_part(qm_seed,1)) = HYBRID_ACTIVE_MARK
         hybrid_mark_p(1:ds_atoms%n) = hybrid_p(1:ds_atoms%n)
	 call print('hybrid_mark '//count(hybrid_mark_p.eq.1))
	 cluster_mark_p(1:ds_atoms%N) = hybrid_mark_p(1:ds_atoms%N)
	 call print('cluster_mark '//count(cluster_mark_p.eq.1))

          !call print("MARKSX1 "//count(hybrid_p(1:ds_atoms%N)==1)//" + "//count(hybrid_p(1:ds_atoms%N)==2)//" hybrid atoms",PRINT_ANAL)
          !call print("MARKSX1 "//count(hybrid_mark_p(1:ds_atoms%N)==1)//" + "//count(hybrid_mark_p(1:ds_atoms%N)==2)//" hybrid_mark atoms",PRINT_ANAL)
          !call print("MARKSX1 "//count(cluster_mark_p(1:ds_atoms%N)==1)//" + "//count(cluster_mark_p(1:ds_atoms%N)==2)//" cluster_mark atoms",PRINT_ANAL)

       endif

        !if (.not. (assign_pointer(ds_atoms, 'hybrid'//trim(mark_postfix), hybrid_p))) call system_abort("??")
        !if (.not.(assign_pointer(ds_atoms, "hybrid_mark"//trim(mark_postfix), hybrid_mark_p))) call system_abort('??')
        !if (.not.(assign_pointer(ds_atoms, "cluster_mark"//trim(mark_postfix), cluster_mark_p))) call system_abort('??')
        !call print("MARKS0 "//count(hybrid_p(1:ds_atoms%N)==1)//" + "//count(hybrid_p(1:ds_atoms%N)==2)//" hybrid atoms",PRINT_ANAL)
        !call print("MARKS0 "//count(hybrid_mark_p(1:ds_atoms%N)==1)//" + "//count(hybrid_mark_p(1:ds_atoms%N)==2)//" hybrid_mark atoms",PRINT_ANAL)
        !call print("MARKS0 "//count(cluster_mark_p(1:ds_atoms%N)==1)//" + "//count(cluster_mark_p(1:ds_atoms%N)==2)//" cluster_mark atoms",PRINT_ANAL)

       !extend QM region around seed atoms
       !update hybrid_*, hybrid_mark_*, old_hybrid_mark_*
       call create_pos_or_list_centred_hybrid_region(ds_atoms,inner_radius,outer_radius,atomlist=qm_seed, &
	 add_only_heavy_atoms=(.not. buffer_general),nneighb_only=.false.,min_images_only=.true., &
	 use_create_cluster_info=use_create_cluster_info_for_core,create_cluster_info_args=create_cluster_info_args, &
	 list_changed=list_changed1, mark_postfix=trim(mark_postfix))

        !if (.not. (assign_pointer(ds_atoms, 'hybrid'//trim(mark_postfix), hybrid_p))) call system_abort("??")
        !if (.not.(assign_pointer(ds_atoms, "hybrid_mark"//trim(mark_postfix), hybrid_mark_p))) call system_abort('??')
        !if (.not.(assign_pointer(ds_atoms, "cluster_mark"//trim(mark_postfix), cluster_mark_p))) call system_abort('??')
        !call print("MARKS1 "//count(hybrid_p(1:ds_atoms%N)==1)//" + "//count(hybrid_p(1:ds_atoms%N)==2)//" hybrid atoms",PRINT_ANAL)
        !call print("MARKS1 "//count(hybrid_mark_p(1:ds_atoms%N)==1)//" + "//count(hybrid_mark_p(1:ds_atoms%N)==2)//" hybrid_mark atoms",PRINT_ANAL)
        !call print("MARKS1 "//count(cluster_mark_p(1:ds_atoms%N)==1)//" + "//count(cluster_mark_p(1:ds_atoms%N)==2)//" cluster_mark atoms",PRINT_ANAL)

    endif ! qm_region_pt_ctr
  end subroutine update_QM_region

end program qmmm_md
