!CP2K driver

!to use external basis set or potential list add to atoms%params BASIS_SET_list=... or POTENTIAL_list=...,
!   or pass to calc in args_str (see below)
!the format is for H and O (for example):
!2                  !number of lines
!1 DZVP-GTH-PADE    !atomic number and basis set
!8 DZVP-GTH-PADE
!
!add HAVE_CP2K=1 to the Makefile.programs
!
!arguments to be passed to go_cp2k in the args_str:
!    Run_Type, default: 'MM', also understands 'QS', 'QMMM_CORE', and 'QMMM_EXTENDED'
!    PSF_Print, default: 'NO_PSF', also understands 'CP2K_PRINT_AND_SAVE', 'DRIVER_PRINT_AND_SAVE', 'USE_EXISTING_PSF'
!    cp2k_program, default:'cp2k.sopt'
!    nproc, default: 1      ! runs mpirun -np $nproc$ if the cp2k program is stored in CP2K_POPT variable, otherwise uses the one in args_str
!    basis_set_file: default none, name of file (instead of atoms%params value BASIS_SET_list)
!    potential_file: default none, name of file (instead of atoms%params value POTENTIAL_list)
!
!contains the following subroutines and functions:
!
!S  create_centred_qmcore(my_atoms,R_inner,R_outer,origin)
!F  extend_qmlist(my_atoms,R_inner,R_outer) result(list_changed)
!S  construct_buffer_RADIUS(my_atoms,core,radius,buffer,use_avgpos,verbosity)
!F  check_qmlist_change(old_qmlist,new_qmlist) result(list_changed)
!F  num_of_bonds(at) result(bonds)
!S  param_initialise(this,at,run_type,nproc,use_cubic_cell)
!S  get_qm_list_int(my_atoms,qmflag,qmlist)
!S  get_qm_list_int_rec(my_atoms,qmflag,qmlist, do_recursive)
!S  get_qm_list_array(my_atoms,qmflag,qmlist)
!F  get_property(my_atoms,prop) result(prop_index)
!S  go_cp2k(my_atoms,forces,energy,args_str)
!S  read_qmlist(my_atoms,qm_list,qmlistfilename,verbose)
!S  check_neighbour_numbers(at)
!S  calc_topology(at,do_CHARMM)
!S  create_CHARMM(at,do_CHARMM)
!S  next_motif(library,res_name,pdb_name,motif,atom_names,atom_charges,do_CHARMM)
!S  write_cp2k_input_files(my_atoms,qm_list,param,run_type,PSF_print)
!S  write_cp2k_pdb_file(at,pdb_file,run_type)
!S  write_cp2k_input_file(my_atoms,qm_list,param,run_type,PSF_print)
!S  read_cp2k_forces(forces,energy,param,my_atoms,run_type,verbose)
!S  read_convert_back_pos(forces,param,my_atoms,run_type,verbose)
!F  real_feq2(x,y) result(feq)
!F  matrix_feq2(matrix1,matrix2) result (feq)
!S  QUIP_combine_forces(qmmm_forces,mm_forces,combined_forces,qmlist,my_atoms)
!F  spline_force(at, i) result(force)
!S  combine_forces(qmmm_forces,mm_forces,combined_forces,qmlist,my_atoms)
!S  energy_conversion(energy)
!S  force_conversion(force)
!S  velocity_conversion(at)
!S  velocity_conversion_rev(at)

module cp2k_driver_module

!  use libatoms_module

  use atoms_module,            only: atoms, initialise, finalise, &
                                     read_xyz, add_property, &
                                     map_into_cell, calc_dists, &
                                     set_cutoff, set_cutoff_minimum, &
                                     calc_connect, DEFAULT_NNEIGHTOL, &
                                     distance_min_image, &
                                     read_line, parse_line, &
                                     assignment(=), atoms_n_neighbours, remove_bond, &
				     assign_pointer, print_xyz
  use clusters_module,         only: bfs_step,&
                                     construct_buffer, &
                                     select_hysteretic_quantum_region, &
                                     add_cut_hydrogens
  use dictionary_module,       only: dictionary, initialise, finalise, &
                                     get_value, set_value, &
                                     value_len, read_string
  use linearalgebra_module,    only: find_in_array, find, &
                                     print, operator(.feq.), &
                                     check_size
  use paramreader_module,      only: param_register, param_read_line, &
                                     FIELD_LENGTH
  use periodictable_module,    only: ElementName, ElementCovRad, ElementMass
  use structures_module,       only: find_motif
  use system_module,           only: dp, inoutput, initialise, finalise, &
                                     INPUT, OUTPUT, INOUT, &
                                     system_timer, system_command, system_abort, &
                                     optional_default, &
                                     print, print_title, &
                                     string_to_int, string_to_real, round, &
                                     parse_string, read_line, &
                                     operator(//), &
                                     NORMAL, ANAL, NERD, ERROR, &
				     verbosity_push_decrement, verbosity_pop, current_verbosity, &
				     mainlog
  use table_module,            only: table, initialise, finalise, &
                                     append, allocate, delete, &
                                     int_part, TABLE_STRING_LENGTH
  use topology_module,         only: write_psf_file, create_CHARMM, &
                                     NONE_RUN, QS_RUN, MM_RUN, &
                                     QMMM_RUN_CORE, QMMM_RUN_EXTENDED
  use units_module,            only: HARTREE, BOHR, HBAR

!  use quantumselection_module
  implicit none

  private :: &
            param_initialise, &
            finalise, &
            construct_buffer_RADIUS, &
            get_qm_list_array, &
            get_qm_list_int, &
            get_qm_list_int_rec, &
            read_cp2k_forces, &
            read_convert_back_pos, &
            real_feq2, &
            matrix_feq2, &
            write_cp2k_input_files, &
            write_cp2k_input_file
!            combine_forces, &

  public :: create_centred_qmcore, &
            extend_qmlist, &
            check_qmlist_change, &
            num_of_bonds, &
            get_qm_list, &
            get_property, &
            go_cp2k, &
            read_qmlist, &
            check_neighbour_numbers, &
            calc_topology, &
            delete_metal_connects, &
            write_cp2k_pdb_file, &
            QUIP_combine_forces, &
            spline_force, &
            energy_conversion, &
            force_conversion, &
            velocity_conversion, &
            velocity_conversion_rev

!parameters for Run_Type
! see in ../libAtoms/Topology.f95
!parameters for PSF_Print
  integer, parameter :: CP2K_PRINT_AND_SAVE = -1
  integer, parameter :: NO_PSF = 0
  integer, parameter :: DRIVER_PRINT_AND_SAVE = 1
  integer, parameter :: USE_EXISTING_PSF = 2

  integer,  parameter, private :: MAX_CHAR_LENGTH = 200 ! for the system_command
  real(dp), parameter :: DEFAULT_QM_BUFFER_WIDTH = 4.5_dp

  real(dp), parameter :: EPSILON_ZERO = 1.e-8_dp  !to check if xyz atoms equal

  type spline_pot
    real(dp) :: from
    real(dp) :: to
    real(dp) :: dpot
  end type spline_pot

!!!!! ============================================================================================ !!!!!
!!!!!                                                                                              !!!!!
!!!!! \/ \/ \/ \/ \/ \/ ----------    parameters type  to run CP2K    ---------- \/ \/ \/ \/ \/ \/ !!!!!
!!!!!                                                                                              !!!!!
!!!!! ============================================================================================ !!!!!

  type param_dft_nml
    character(len=value_len) :: file
    real(dp)      :: mgrid_cutoff != 300.0_dp
    character(20) :: scf_guess != 'atomic'
    character(20) :: xc_functional !!= 'PADE'
  end type param_dft_nml

  type param_mm_forcefield_nml
    character(20) :: parmtype != 'CHM'
    character(20) :: parm_file_name != 'charmm.pot'
    real(dp)      :: charge_OT != -0.834_dp
    real(dp)      :: charge_HT != 0.417_dp
  end type param_mm_forcefield_nml

  type param_mm_ewald_nml
    character(20) :: ewald_type != 'ewald'
    real(dp)      :: ewald_alpha != 0.44_dp
    integer       :: ewald_gmax != 25
  end type param_mm_ewald_nml

  type param_qmmm_nml
    real(dp),dimension(3) :: qmmm_cell != (/8.0_dp,8.0_dp,8.0_dp/)
    character(20)         :: qmmm_cell_unit != 'ANGSTROM'
    logical               :: reuse_wfn
    character(20)         :: ecoupl != 'GAUSS'
    real(dp)              :: radius_H != 0.44_dp
    real(dp)              :: radius_O != 0.78_dp
    real(dp)              :: radius_C != 0.80_dp
    real(dp)              :: radius_N != 0.80_dp
    real(dp)              :: radius_S != 0.80_dp
    real(dp)              :: radius_Cl != 0.80_dp
    real(dp)              :: radius_Na != 0.80_dp
    real(dp)              :: radius_Si !=0.80_dp
  end type param_qmmm_nml

  type param_global_nml
    character(20) :: project != 'cp2k'
    character(20) :: runtype != 'MD'
    character(len=value_len) :: cell_file
  end type param_global_nml

  type param_md_nml
    character(20) :: ensemble != 'NVE'
    integer       :: steps != 0
    real(dp)      :: timestep != 0.5_dp
    real(dp)      :: temperature != 300.0_dp
  end type param_md_nml

  type working_env
    character(80) :: working_directory    !directory for temporary files used by cp2k
    character(80) :: pdb_file
    character(80) :: psf_file
    character(80) :: wfn_file
    character(80) :: cp2k_input_filename !='cp2k_input.inp'
    character(80) :: force_file != 'cp2k-frc-1.xyz'
    character(80) :: pos_file != 'cp2k-pos-1.xyz'
    character(80) :: wrk_filename != 'wrk.dat'
    character(800) :: cp2k_program != 'cp2k_serial' or 'cp2k_popt'. long, for possible long path
  end type working_env

  type param_cp2k

    type(Dictionary)              :: basis_set
    type(Dictionary)              :: potential
    type(param_dft_nml)           :: dft
    type(param_mm_forcefield_nml) :: mm_forcefield
    type(param_mm_ewald_nml)      :: mm_ewald
    type(param_qmmm_nml)          :: qmmm
    type(param_global_nml)        :: global
    type(param_md_nml)            :: md
    type(working_env)             :: wenv

    integer                       :: N
    logical                       :: initialised

  end type param_cp2k

  interface finalise
     module procedure param_finalise
  end interface finalise

  interface get_qm_list
     module procedure get_qm_list_int, get_qm_list_int_rec, get_qm_list_array
  end interface get_qm_list

contains

  !!!!! ============================================================================================ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! \/ \/ \/ \/ \/ \/ ---------- initialise parameters  to run CP2K ---------- \/ \/ \/ \/ \/ \/ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! ============================================================================================ !!!!!

  subroutine param_initialise(this,at,run_type,nproc,cp2k_program,basis_set_file,potential_file,&
    dft_file,cell_file,use_cubic_cell)

  type(param_cp2k),  intent(out) :: this
  type(atoms),       intent(inout) :: at
  integer,           intent(in)  :: run_type
  integer,           intent(in)  :: nproc
  character(len=*),  intent(in)  :: cp2k_program
  character(len=*),  intent(in)  :: basis_set_file, potential_file, dft_file, cell_file
  logical, optional, intent(in)  :: use_cubic_cell

  real(dp), dimension(3)         :: QM_maxdist
  type(Table)                    :: fitlist
  character(len=20)              :: dir
  integer                        :: i,j
  logical                        :: ex
  logical                        :: cubic
  character(len=800)             :: env_program_name
  integer                        :: name_len, status
  real(dp),dimension(3)          :: old_QM_cell
  logical                        :: QM_list_changed
  character(len=value_len)       :: bs_pot_string
  character(len=value_len)       :: basis_set_file_val, potential_file_val

  call system_timer('param_init')

  if (.not.any(run_type.eq.(/QS_RUN,MM_RUN,QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) &
     call system_abort('param_initialise Run type is not recognized: '//run_type)
  call print('Initializing parameters to run CP2K...')
  cubic = optional_default(.false.,use_cubic_cell)

  if (any(run_type.eq.(/QS_RUN,QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
    this%dft%mgrid_cutoff = 280.0_dp
    if (run_type.eq.QS_RUN .and. any(at%Z(1:at%N).eq.17)) this%dft%mgrid_cutoff = 300._dp
    this%dft%scf_guess = 'atomic'
    this%dft%xc_functional = 'PADE' !can be B3LYP,PBE0,BLYP,BP,PADE,PBE,TPSS,HTCH120,OLYP,NO_SHORTCUT,NONE
  endif
  if (any(run_type.eq.(/MM_RUN,QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
    this%mm_forcefield%parmtype = 'CHM'
    this%mm_forcefield%parm_file_name = 'charmm.pot'
    this%mm_forcefield%charge_OT = -0.834_dp
    this%mm_forcefield%charge_HT = 0.417_dp
  endif
    this%mm_ewald%ewald_type = 'ewald'
    this%mm_ewald%ewald_alpha = 0.35_dp
    this%mm_ewald%ewald_gmax = int(max(at%lattice(1,1),at%lattice(2,2),at%lattice(3,3))/2._dp)*2+1

  if (any(run_type.eq.(/QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
   ! let's have 3 Angstroms on any side of the cell
    call print('INFO: The size of the QM cell is either the MM cell itself, or it will have at least 3-3 Angstrom around the QM atoms.')
    call print('WARNING! Please check if your cell is centered around the QM region!')
    call print('WARNING! CP2K centering algorithm fails if QM atoms are not all in the')
    call print('WARNING! 0,0,0 cell. If you have checked it, please ignore this message.')
    this%qmmm%qmmm_cell = 0._dp
    call get_qm_list_int_rec(at,run_type,fitlist, do_recursive=.true.)  ! with QM_flag 1 or 2
    if (any(at%Z(int_part(fitlist,1)).gt.8)) this%dft%mgrid_cutoff = 300._dp

! whole numbers => to use wavefunction extrapolation later
    QM_maxdist = 0
    do i=1,fitlist%N
      do j=i+1,fitlist%N
         QM_maxdist(1) = max(QM_maxdist(1),distance_min_image(at,(/at%pos(1,i),0._dp,0._dp/),(/at%pos(1,j),0._dp,0._dp/)))
         QM_maxdist(2) = max(QM_maxdist(2),distance_min_image(at,(/0._dp,at%pos(2,i),0._dp/),(/0._dp,at%pos(2,j),0._dp/)))
         QM_maxdist(3) = max(QM_maxdist(3),distance_min_image(at,(/0._dp,0._dp,at%pos(3,i)/),(/0._dp,0._dp,at%pos(3,j)/)))
!         call print('QM_maxdist: '//round(QM_maxdist(1),3)//' '//round(QM_maxdist(2),3)//' '//round(QM_maxdist(3),3))
      enddo
    enddo

    if (cubic) then
       QM_maxdist(1) = maxval(QM_maxdist(1:3))
       QM_maxdist(2) = QM_maxdist(1)
       QM_maxdist(3) = QM_maxdist(1)
       call print('cubic cell')
    endif

    this%qmmm%qmmm_cell(1) = min(real(ceiling(QM_maxdist(1)))+6._dp,at%lattice(1,1))
    this%qmmm%qmmm_cell(2) = min(real(ceiling(QM_maxdist(2)))+6._dp,at%lattice(2,2))
    this%qmmm%qmmm_cell(3) = min(real(ceiling(QM_maxdist(3)))+6._dp,at%lattice(3,3))

    call print('Lattice:')
    call print('      A '//round(at%lattice(1,1),6)//' '//round(at%lattice(2,1),6)//' '//round(at%lattice(3,1),6))
    call print('      B '//round(at%lattice(1,2),6)//' '//round(at%lattice(2,2),6)//' '//round(at%lattice(3,2),6))
    call print('      C '//round(at%lattice(1,3),6)//' '//round(at%lattice(2,3),6)//' '//round(at%lattice(3,3),6))
    call print('QM cell:')
    call print('      A '//round(this%qmmm%qmmm_cell(1),6)//' '//round(0._dp,6)//' '//round(0._dp,6))
    call print('      B '//round(0._dp,6)//' '//round(this%qmmm%qmmm_cell(2),6)//' '//round(0._dp,6))
    call print('      C '//round(0._dp,6)//' '//round(0._dp,6)//' '//round(this%qmmm%qmmm_cell(3),6))

   ! save the QM cell size, if the same as last time, might reuse wavefunction
    ex = .false.
    this%qmmm%reuse_wfn = .false.
    ex = get_value(at%params,'QM_cell',old_QM_cell)
    if (ex) then
       if (old_QM_cell(1:3) .feq. this%qmmm%qmmm_cell(1:3)) then
          this%qmmm%reuse_wfn = .true.
          if (Run_Type.eq.QMMM_RUN_CORE) call print('QM cell size have not changed. Reuse previous wavefunction.')
!          call print('QM_cell size has not changed in at%params, will reuse wfn'//old_QM_cell(1:3)//this%qmmm%qmmm_cell(1:3)//' if QM list has not changed')
       else
          this%qmmm%reuse_wfn = .false.
!          call print('QM_cell size changed in at%params, will not reuse wfn'//old_QM_cell(1:3)//this%qmmm%qmmm_cell(1:3))
          call set_value(at%params,'QM_cell',this%qmmm%qmmm_cell(1:3))
       endif
    else !for the first step there is no QM_cell, no wavefunction
       this%qmmm%reuse_wfn = .false.
!       call print('did not find QM_cell property in at%params, will not reuse wfn')
       call set_value(at%params,'QM_cell',this%qmmm%qmmm_cell(1:3))
    endif
   ! only in case of extended run also QM buffer zone change counts
    if (Run_Type.eq.QMMM_RUN_EXTENDED) then
       ex = .false.
       ex = get_value(at%params,'QM_list_changed',QM_list_changed)
       if (ex) then
       this%qmmm%reuse_wfn = this%qmmm%reuse_wfn .and. (.not.QM_list_changed)
!       if (QM_list_changed) then
!          call print('QM list changed, will not use wfn')
!       else
!          call print('QM list has not changed, will use wfn if QM cell is the same')
!       endif
       else
          this%qmmm%reuse_wfn = .false.
       endif
       if (this%qmmm%reuse_wfn) then
          this%dft%scf_guess = 'RESTART'
          call print('QM list and QM cell size have not changed. Reuse previous wavefunction.')
       endif
    endif

!    call print('QM coordinates:')
!    do i=1,fitlist%N
!    call print('QM ATOM '//round(at%pos(1,fitlist%int(1,i)),6)//' '//round(at%pos(2,fitlist%int(1,i)),6)//' '//round(at%pos(3,fitlist%int(1,i)),6))
!    enddo

    this%qmmm%qmmm_cell_unit = 'ANGSTROM'
    this%qmmm%ecoupl = 'GAUSS'
! ONLY 0.44, 0.78 & 0.80 values can be used at the moment!
    this%qmmm%radius_H  = 0.44_dp
    this%qmmm%radius_O  = 0.78_dp
    this%qmmm%radius_C  = 0.80_dp !0.77_dp
    this%qmmm%radius_N  = 0.80_dp !0.75_dp
    this%qmmm%radius_S  = 0.80_dp !0.75_dp
    this%qmmm%radius_Cl = 0.80_dp !0.75_dp
    this%qmmm%radius_Na = 0.80_dp !0.75_dp
  endif

  if (any(run_type.eq.(/QS_RUN,QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then

    ! read basis_set/potentials from file
    call initialise(this%basis_set)
    basis_set_file_val = basis_set_file
    if (len(trim(basis_set_file_val)) == 0) then ! if no file was passed in calc_str, try to get it from config header
      if (get_value(at%params,'BASIS_SET_list',basis_set_file_val)) then ! if file was defined in config header, set filename
         call print("Using basis_set_file from BASIS_SET_list property of config")
      endif
    endif
    if (len(trim(basis_set_file_val)) /= 0) then ! if a filename has been defined anywhere
      call print("Reading basis set from '"//trim(basis_set_file_val)//"'")
      call read_list_file(this,basis_set_file_val,basis_set=.true.)
    endif

    this%dft%file = trim(dft_file)

    call initialise(this%potential)
    potential_file_val = potential_file
    if (len(trim(potential_file_val)) == 0) then ! if no file was passed in calc_str, try to get it from config header
      if (get_value(at%params,'POTENTIAL_list',potential_file_val)) then ! if file was defined in config header, set filename
         call print("Using potential_file from POTENTIAL_list property of config")
      endif
    endif
    if (len(trim(potential_file_val)) /= 0) then ! if a filename has been defined anywhere
      call print("Reading potential from '"//trim(potential_file_val)//"'")
      call read_list_file(this,potential_file_val,potential=.true.)
    endif

! add default basis sets and potentials:
    if (.not.get_value(this%basis_set,'H',bs_pot_string).and.find_in_array(at%Z(1:at%N),1).gt.0) then
       call print('Added default basis set for H: DZVP-GTH-PADE')
       call set_value(this%basis_set,'H','DZVP-GTH-PADE')
    endif
    if (.not.get_value(this%potential,'H',bs_pot_string).and.find_in_array(at%Z(1:at%N),1).gt.0) then
       call print('Added default potential for H: GTH-PADE-q1')
       call set_value(this%potential,'H','GTH-PADE-q1')
    endif
    if (.not.get_value(this%basis_set,'O',bs_pot_string).and.find_in_array(at%Z(1:at%N),8).gt.0) then
       call print('Added default basis set for O: DZVP-GTH-PADE')
       call set_value(this%basis_set,'O','DZVP-GTH-PADE')
    endif
    if (.not.get_value(this%potential,'O',bs_pot_string).and.find_in_array(at%Z(1:at%N),8).gt.0) then
       call print('Added default potential for O: GTH-PADE-q6')
       call set_value(this%potential,'O','GTH-PADE-q6')
    endif
    if (.not.get_value(this%basis_set,'C',bs_pot_string).and.find_in_array(at%Z(1:at%N),6).gt.0) then
       call print('Added default basis set for C: DZVP-GTH-PADE')
       call set_value(this%basis_set,'C','DZVP-GTH-PADE')
    endif
    if (.not.get_value(this%potential,'C',bs_pot_string).and.find_in_array(at%Z(1:at%N),6).gt.0) then
       call print('Added default potential for C: GTH-PADE-q4')
       call set_value(this%potential,'C','GTH-PADE-q4')
    endif
    if (.not.get_value(this%basis_set,'N',bs_pot_string).and.find_in_array(at%Z(1:at%N),7).gt.0) then
       call print('Added default basis set for N: DZVP-GTH-PADE')
       call set_value(this%basis_set,'N','DZVP-GTH-PADE')
    endif
    if (.not.get_value(this%potential,'N',bs_pot_string).and.find_in_array(at%Z(1:at%N),7).gt.0) then
       call print('Added default potential for N: GTH-PADE-q5')
       call set_value(this%potential,'N','GTH-PADE-q5')
    endif
    if (.not.get_value(this%basis_set,'S',bs_pot_string).and.find_in_array(at%Z(1:at%N),16).gt.0) then
       call print('Added default basis set for S: DZVP-GTH-PADE')
       call set_value(this%basis_set,'S','DZVP-GTH-PADE')
    endif
    if (.not.get_value(this%potential,'S',bs_pot_string).and.find_in_array(at%Z(1:at%N),16).gt.0) then
       call print('Added default potential for S: GTH-PADE-q6')
       call set_value(this%potential,'S','GTH-PADE-q6')
    endif
    if (.not.get_value(this%basis_set,'Cl',bs_pot_string).and.find_in_array(at%Z(1:at%N),17).gt.0) then
       call print('Added default basis set for Cl: DZVP-GTH-PADE')
       call set_value(this%basis_set,'Cl','DZVP-GTH-PADE')
    endif
    if (.not.get_value(this%potential,'Cl',bs_pot_string).and.find_in_array(at%Z(1:at%N),17).gt.0) then
       call print('Added default potential for Cl: GTH-PADE-q7')
       call set_value(this%potential,'Cl','GTH-PADE-q7')
    endif
    if (.not.get_value(this%basis_set,'F',bs_pot_string).and.find_in_array(at%Z(1:at%N),9).gt.0) then
       call print('Added default basis set for F: DZVP-GTH-PADE')
       call set_value(this%basis_set,'F','DZVP-GTH-PADE')
    endif
    if (.not.get_value(this%potential,'F',bs_pot_string).and.find_in_array(at%Z(1:at%N),9).gt.0) then
       call print('Added default potential for F: GTH-PADE-q7')
       call set_value(this%potential,'F','GTH-PADE-q7')
    endif
    if (.not.get_value(this%basis_set,'Na',bs_pot_string).and.find_in_array(at%Z(1:at%N),11).gt.0) then
       call print('Added default basis set for Na: DZVP-GTH-PADE')
       call set_value(this%basis_set,'Na','DZVP-GTH-PADE')
    endif
    if (.not.get_value(this%potential,'Na',bs_pot_string).and.find_in_array(at%Z(1:at%N),11).gt.0) then
       call print('Added default potential for Na: GTH-PADE-q1')
       call set_value(this%potential,'Na','GTH-PADE-q1')
    endif
    if (.not.get_value(this%basis_set,'Fe',bs_pot_string).and.find_in_array(at%Z(1:at%N),26).gt.0) then
       call print('Added default basis set for Fe: DZVP-GTH-PADE')
       call set_value(this%basis_set,'Fe','DZVP-GTH-PADE')
    endif
    if (.not.get_value(this%potential,'Fe',bs_pot_string).and.find_in_array(at%Z(1:at%N),26).gt.0) then
       call print('Added default potential for Fe: GTH-PADE-q16')
       call set_value(this%potential,'Fe','GTH-PADE-q16')
    endif
    if (.not.get_value(this%basis_set,'Si',bs_pot_string).and.find_in_array(at%Z(1:at%N),14).gt.0) then
       call print('Added default basis set for Si: DZVP-GTH-PADE')
       call set_value(this%basis_set,'Si','DZVP-GTH-PADE')
    endif
    if (.not.get_value(this%potential,'Si',bs_pot_string).and.find_in_array(at%Z(1:at%N),14).gt.0) then
       call print('Added default potential for Si: GTH-PADE-q4')
       call set_value(this%potential,'Si','GTH-PADE-q4')
    endif
  endif

    this%global%project = 'cp2k'
    this%global%runtype = 'MD'
    this%global%cell_file = trim(cell_file)

    this%md%ensemble = 'NVE'
    this%md%steps = 0
    this%md%timestep = 0.5_dp
    this%md%temperature = 300.0_dp

! use numbers in ascending order
!    this%wenv%working_directory = 'cp2k_run_'//ran_string(8)    !directory for temporary files used by cp2k
    ex = .true.
    i=1
    do while (ex)
       write (dir,'(a,i0)') 'cp2k_run_',i
       inquire(file=(trim(dir)//'/cp2k_input.inp'),exist=ex)
       i = i + 1
    enddo

!added status optionally to system_command!
    call system_command('mkdir '//trim(dir),status=status)
    call print('system_command status: '//status)
    if (status /= 0) then
      call system_abort('Failed to mkdir '//trim(dir))
    endif

    call print('the working directory: '//trim(dir))
    this%wenv%working_directory = trim(dir)
    this%wenv%pdb_file = 'pdb.CHARMM.pdb'
    this%wenv%psf_file = 'psf.CHARMM.psf'
    this%wenv%wfn_file = 'wfn.restart.wfn'
    this%wenv%cp2k_input_filename = 'cp2k_input.inp'
    this%wenv%force_file = trim(this%global%project)//'-frc-1.xyz'
    this%wenv%pos_file = trim(this%global%project)//'-pos-1.xyz'
    this%wenv%wrk_filename    = 'wrk.dat'

    this%N=at%N

   ! get program name from environment variable
    call print('use cp2k_program passed in args_str: '//trim(cp2k_program))
    this%wenv%cp2k_program = trim(cp2k_program)

  call system_timer('param_init')

  end subroutine param_initialise

  subroutine param_finalise(this)

  type(param_cp2k),  intent(out) :: this

    call finalise(this%basis_set)
    call finalise(this%potential)
    this%dft%file = ''
    this%dft%mgrid_cutoff = 0._dp
    this%dft%scf_guess = ''
    this%dft%xc_functional = ''
    this%mm_forcefield%parmtype = ''
    this%mm_forcefield%parm_file_name = ''
    this%mm_forcefield%charge_OT = 0._dp
    this%mm_forcefield%charge_HT = 0._dp
    this%mm_ewald%ewald_type = ''
    this%mm_ewald%ewald_alpha = 0._dp
    this%mm_ewald%ewald_gmax = 0
    this%qmmm%qmmm_cell = 0._dp
    this%qmmm%qmmm_cell_unit = ''
    this%qmmm%reuse_wfn = .false.
    this%qmmm%ecoupl = ''
    this%qmmm%radius_H = 0._dp
    this%qmmm%radius_O = 0._dp
    this%qmmm%radius_C = 0._dp
    this%qmmm%radius_N = 0._dp
    this%qmmm%radius_S = 0._dp
    this%qmmm%radius_Cl = 0._dp
    this%qmmm%radius_Na = 0._dp
    this%qmmm%radius_Si = 0._dp
    this%global%cell_file = ''
    this%global%project = ''
    this%global%runtype = ''
    this%md%ensemble = ''
    this%md%steps = 0
    this%md%timestep = 0._dp
    this%md%temperature = 0._dp
    this%wenv%working_directory = ''
    this%wenv%pdb_file = ''
    this%wenv%psf_file = ''
    this%wenv%wfn_file = ''
    this%wenv%cp2k_input_filename = ''
    this%wenv%force_file = ''
    this%wenv%pos_file = ''
    this%wenv%wrk_filename = ''
    this%wenv%cp2k_program = ''
    this%N = 0

    this%initialised = .false.

  end subroutine param_finalise

 ! updates the core QM flags (=1), where an atom leaves becomes part of buffer (QM_flag=2)
  subroutine create_centred_qmcore(my_atoms,R_inner,R_outer,origin,list_changed)

    type(Atoms),        intent(inout) :: my_atoms
    real(dp),           intent(in)    :: R_inner
    real(dp),           intent(in)    :: R_outer
    real(dp), optional, intent(in)    :: origin(3)

    type(Atoms) :: atoms_for_add_cut_hydrogens
    logical     :: list_changed
    type(Table) :: inner_list, outer_list, core, core1, ext_qmlist
    integer     :: qm_flag_index, i
    real(dp)    :: my_origin(3)

    my_origin = optional_default((/0._dp,0._dp,0._dp/),origin)

    call map_into_cell(my_atoms)

    call allocate(core,4,0,0,0,0)
    call get_qm_list_int(my_atoms,1,core1)
    call get_qm_list_int(my_atoms,1,core)
    call get_qm_list_int_rec(my_atoms,2,ext_qmlist, do_recursive=.true.)

    call construct_buffer_origin(my_atoms,R_inner,inner_list,my_origin)
    call construct_buffer_origin(my_atoms,R_outer,outer_list,my_origin)

    call select_hysteretic_quantum_region(my_atoms,inner_list,outer_list,core)
    call finalise(inner_list)
    call finalise(outer_list)

!TO BE OPTIMIZED : add avgpos to add_cut_hydrogen
   ! add cut hydrogens, according to avgpos
    atoms_for_add_cut_hydrogens = my_atoms
    atoms_for_add_cut_hydrogens%oldpos = my_atoms%avgpos
    atoms_for_add_cut_hydrogens%avgpos = my_atoms%avgpos
    atoms_for_add_cut_hydrogens%pos = my_atoms%avgpos
    call set_cutoff(atoms_for_add_cut_hydrogens,DEFAULT_NNEIGHTOL)
    call calc_connect(atoms_for_add_cut_hydrogens)
    call add_cut_hydrogens(atoms_for_add_cut_hydrogens,core)
    !call print('Atoms in hysteretic quantum region after adding the cut hydrogens:')
    !do i=1,core%N
    !   call print(core%int(1,i))
    !enddo
    call finalise(atoms_for_add_cut_hydrogens)

   ! check changes in QM list and set the new QM list
    list_changed = check_qmlist_change(old_qmlist=core1,new_qmlist=core)
    if (list_changed)  call print('QM list around the origin  has changed')

   ! update QM_flag of my_atoms
    qm_flag_index = get_property(my_atoms,'QM_flag')
    my_atoms%data%int(qm_flag_index,1:my_atoms%N) = 0
    my_atoms%data%int(qm_flag_index,int_part(ext_qmlist,1)) = 2
    my_atoms%data%int(qm_flag_index,int_part(core,1)) = 1

    call finalise(core)
    call finalise(core1)
    call finalise(ext_qmlist)

  end subroutine create_centred_qmcore

  subroutine print_qm_region(at, file)
    type(Atoms), intent(inout) :: at
    type(Inoutput), optional :: file

    integer, pointer :: qm_flag(:)

    if (.not.(assign_pointer(at, "QM_flag", qm_flag))) &
      call system_abort("print_qm_region couldn't find QM_flag property")

    if (present(file)) then
      file%prefix="QM_REGION"
      call print_xyz(at, file, mask=(qm_flag /= 0))
      file%prefix=""
    else
      mainlog%prefix="QM_REGION"
      call print_xyz(at, mainlog, mask=(qm_flag /= 0))
      mainlog%prefix=""
    end if
  end subroutine print_qm_region

  subroutine construct_buffer_origin(my_atoms,Radius,list,origin)

    type(Atoms),        intent(inout) :: my_atoms
    real(dp),           intent(in)    :: Radius
    type(Table),        intent(out)   :: list
    real(dp), optional, intent(in)    :: origin(3)

    integer  :: i
    real(dp) :: my_origin(3)

    my_origin = optional_default((/0._dp,0._dp,0._dp/),origin)

    call initialise(list,4,0,0,0,0)

    do i = 1,my_atoms%N
       if (my_atoms%Z(i).eq.1) cycle
       if (distance_min_image(my_atoms,i,my_origin(1:3)).lt.Radius) call append(list,(/i,0,0,0/))
    enddo

  end subroutine construct_buffer_origin

 ! to read BASIS_SET and POTENTIAL from external file 
  subroutine read_list_file(this,filename,basis_set,potential)

    type(param_cp2k),  intent(inout) :: this
    character(len=*),  intent(in)    :: filename
    logical, optional, intent(in)    :: basis_set,potential

    type(InOutput)                  :: listfile
    integer                         :: num_list, n, Z
    character(80), dimension(2)     :: fields
    character(len=80)               :: bs_pot
    integer                         :: num_fields
    logical                         :: do_basis_set
    integer                         :: status

    if (present(basis_set).and.present(potential)) then
       if (basis_set.eqv.potential) call system_abort('read_list_file: exactly one of basis_set and potential should be true')
    endif
    if (present(basis_set)) do_basis_set = basis_set
    if (present(potential)) do_basis_set = .not.potential

    call initialise(listfile,filename=trim(filename),action=INPUT)

    call parse_line(listfile,' ',fields,num_fields,status)
    if (status > 0) then
       call system_abort('read_list_file: Error reading from '//listfile%filename)
    else if (status < 0) then
       call system_abort('read_list_file: End of file when reading from '//listfile%filename)
    end if

    num_list = string_to_int(fields(1))
    if (num_list.le.0) call system_abort('read_list_file: number of basis/potential should be >0:'//num_list)
!    if (do_basis_set) then
!        call initialise(this%basis_set)
!       call allocate(this%basis_set,1,0,1,0,num_list)      !1 int, 0 reals, 0 str, 0 log, num_qm_atoms entries
!    else
!        call initialise(this%potential)
!       call allocate(this%basis_set,1,0,1,0,num_list)      !1 int, 0 reals, 0 str, 0 log, num_qm_atoms entries
!    endif

   ! Reading and storing list...
    do n=1,num_list
       call parse_line(listfile,' ',fields,num_fields,status)
       if (status.ne.0) call system_abort('wrong file '//listfile%filename)
       if (num_fields.lt.2) call system_abort('wrong file format '//listfile%filename)
       Z = string_to_int(fields(1))
       bs_pot = fields(2)
!       if (len(trim(bs_pot)).gt.TABLE_STRING_LENGTH) call system_abort('too long (>10) name: '//trim(bs_pot))
       if (do_basis_set) then
          call print('Use for '//ElementName(Z)//' basis set: '//trim(bs_pot))
          call set_value(this%basis_set,ElementName(Z),bs_pot)
!          call append(this%basis_set,(/Z/))
!          call append(this%basis_set,bs_pot)
       else
          call print('Use for '//ElementName(Z)//' potential: '//trim(bs_pot))
          call set_value(this%potential,ElementName(Z),bs_pot)
!          call append(this%potential,(/Z/))
!          call append(this%potential,bs_pot)
       endif
    enddo

    call finalise(listfile)

  end subroutine read_list_file


  !!!!! ============================================================================================ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! \/ \/ \/ \/ \/ \/ ------------  main subroutine  to run CP2K  ------------ \/ \/ \/ \/ \/ \/ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! ============================================================================================ !!!!!

  subroutine go_cp2k(my_atoms,forces,energy,args_str)

  type(Atoms),                        intent(inout)  :: my_atoms
  real(dp), dimension(:,:), optional, intent(out) :: forces
  real(dp),                 optional, intent(out) :: energy
  character(len=*),         optional, intent(in)  :: args_str

  type(param_cp2k)                      :: param
  type(Table)                           :: qmlist
  real(dp), allocatable, dimension(:,:) :: CP2K_forces
  real(dp)                              :: CP2K_energy

  type(Dictionary)                      :: params
  integer                               :: run_type
  integer                               :: nproc, PSF_print
  character(len=FIELD_LENGTH)           :: cp2k_program
  character(len=FIELD_LENGTH)           :: run_type_str, psf_print_str
  character(len=FIELD_LENGTH)           :: fileroot_str
  character(len=FIELD_LENGTH)           :: basis_set_file, potential_file, dft_file, cell_file

  character(max_char_length)            :: run_command='', &
                                           fin_command=''
  integer                               :: status
  logical                               :: ex


    call system_timer('go_cp2k_start')
    call print_title('Setting up the CP2K run')

    if (.not.present(energy).and..not.present(forces)) call system_abort('go_cp2k: nothing to be calculated. Neither energy, nor forces.')

   ! Read and print run parameters
    fileroot_str=''
    basis_set_file=''
    potential_file=''
    dft_file=''
    cell_file=''
    call initialise(params)
    call param_register(params, 'Run_Type', 'MM', run_type_str)
    call param_register(params, 'nproc', '1', nproc)
    call param_register(params, 'PSF_Print', 'NO_PSF', psf_print_str)
    call param_register(params, 'cp2k_program','cp2k.sopt', cp2k_program)
    call param_register(params, 'root', '', fileroot_str)
    call param_register(params, 'basis_set_file', '', basis_set_file)
    call param_register(params, 'potential_file', '', potential_file)
    call param_register(params, 'dft_file', '', dft_file)
    call param_register(params, 'cell_file', '', cell_file)

    if (present(args_str)) then
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.)) &
      call system_abort("Potential_Initialise_str failed to parse args_str='"//trim(args_str)//"'")
    else
      call system_abort("Potential_Initialise_str failed to parse args_str, no args_str")
    endif
    call finalise(params)

    if (len_trim(fileroot_str) > 0) then
      if (len_trim(basis_set_file) == 0) basis_set_file=trim(fileroot_str)//'.basis_set'
      if (len_trim(potential_file) == 0) potential_file=trim(fileroot_str)//'.potential'
      if (len_trim(dft_file) == 0) dft_file=trim(fileroot_str)//'.dft'
      if (len_trim(cell_file) == 0) cell_file=trim(fileroot_str)//'.cell'
    end if

    select case(trim(psf_print_str))
      case('CP2K_PRINT_AND_SAVE') 
        PSF_print = CP2K_PRINT_AND_SAVE
      case('NO_PSF') 
        PSF_print = NO_PSF
      case('DRIVER_PRINT_AND_SAVE') 
        PSF_print = DRIVER_PRINT_AND_SAVE
      case('USE_EXISTING_PSF') 
        PSF_print = USE_EXISTING_PSF
      case default
        call system_abort("go_cp2k got psf_print_str='"//trim(psf_print_str)//"'" // &
          ", only know about NO_PSF, CP2K_PRINT_AND_SAVE, DRIVER_PRINT_AND_SAVE, USE_EXISTING_PSF")
    end select

    select case(trim(run_type_str))
      case('MM') 
        run_type = MM_RUN
      case('QS')
        run_type = QS_RUN
      case('QMMM_CORE')
        run_type = QMMM_RUN_CORE
      case('QMMM_EXTENDED')
        run_type = QMMM_RUN_EXTENDED
      case default
        call system_abort("go_cp2k got run_type_str='"//trim(run_type_str)//"'" // &
          ", only know about MM, QS, QMMM_CORE, QMMM_EXTENDED")
    end select

   if (.not.any(run_type.eq.(/QS_RUN,MM_RUN,QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) &
      call system_abort('go_cp2k: Run type is not recognized: '//run_type)
   if (.not.any(PSF_print.eq.(/NO_PSF,CP2K_PRINT_AND_SAVE,DRIVER_PRINT_AND_SAVE,USE_EXISTING_PSF/))) &
      call system_abort('go_cp2k: PSF print is not recognized: '//PSF_print)

   call print ('Run parameters: ')
   if (run_type.eq.QS_RUN) call print ('Run type: QS')
   if (run_type.eq.MM_RUN) call print ('Run type: MM')
   if (run_type.eq.QMMM_RUN_CORE) call print ('Run type: QM/MM run, QM: quantum core')
   if (run_type.eq.QMMM_RUN_EXTENDED) call print ('Run type: QM/MM run, QM: extended quantum region')
   call print ('Number of processors used: '//nproc)
   if (PSF_Print.eq.NO_PSF) call print ('PSF printing: No PSF')
   if (PSF_Print.eq.CP2K_PRINT_AND_SAVE) call print ('PSF printing: CP2K prints PSF')
   if (PSF_Print.eq.DRIVER_PRINT_AND_SAVE) call print ('PSF printing: driver prints PSF')
   if (PSF_Print.eq.USE_EXISTING_PSF) call print ('PSF printing: Use existing PSF')

  ! If QM/MM, get QM list
    if (any(run_type.eq.(/QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
       call get_qm_list_int_rec(my_atoms,run_type,qmlist, do_recursive=.true.)
      !check if all the atoms are QM, QS used instead
       if (qmlist%N.eq.my_atoms%N) then
          call print('WARNING: go_cp2k: QM/MM calculation was requested, but all the atoms are within the quantum zone: fully QS will be run!')
          run_type = QS_RUN
          call print ('Run type: QS')
       endif
    endif

  ! set params, CHARMM formats and create working directory
    call param_initialise(param,my_atoms,run_type,nproc,cp2k_program,basis_set_file,potential_file,dft_file,cell_file)

  ! check CHARMM topology
  ! Write the coordinate file and the QM/MM or MM input file
    call write_cp2k_input_files(my_atoms,qmlist,param=param,run_type=run_type,PSF_print=PSF_print)

  ! Run cp2k MM serial or QM/MM parallel
    call print('Running CP2K...')

    run_command = 'cd '//trim(param%wenv%working_directory)//';'//trim(param%wenv%cp2k_program)//' '//trim(param%wenv%cp2k_input_filename)//' > cp2k_output.out'

    call print(run_command)
!added status optionally to system_command.
    call system_timer('go_cp2k_start')
    call system_timer('system_command(run_command)')
    call system_command(run_command,status=status)
    call system_timer('system_command(run_command)')
    call system_timer('go_cp2k_end')
    call print_title('...CP2K...')
    call print('grep -i warning '//trim(param%wenv%working_directory)//'/cp2k_output.out')
    call system_command('grep -i warning '//trim(param%wenv%working_directory)//'/cp2k_output.out')
    call system_command('grep -i error '//trim(param%wenv%working_directory)//'/cp2k_output.out')
    if (status.ne.0) call system_abort('CP2K aborted. See output file '//trim(param%wenv%working_directory)//'/cp2k_output.out')

  ! Read energy and forces - also re-reordering!
!    call print('file='//trim(param%wenv%working_directory)//'/'//trim(param%wenv%force_file)//',exist=ex')
!    inquire(file=trim(param%wenv%working_directory)//'/'//trim(param%wenv%force_file),exist=ex)
!    if (ex) call system_abort('CP2K aborted before printing forces. Check the output file '//trim(param%wenv%working_directory)//'/cp2k_output.out')
    call read_cp2k_forces(CP2K_forces,CP2K_energy,param,my_atoms,run_type=run_type)

    if (present(energy)) energy = CP2K_energy
    if (present(forces)) forces = CP2K_forces

   ! Save PSF if needed
    if (PSF_print.eq.CP2K_PRINT_AND_SAVE) then
!       inquire(file=trim(param%wenv%working_directory)//'/'//trim(param%global%project)//'-dump-1.psf',exist=ex)
!       if (ex) &
        call system_command('cp '//trim(param%wenv%working_directory)//'/'//trim(param%global%project)//'-dump-1.psf '//trim(param%wenv%psf_file))
    endif

   ! Save Wfn if needed
    if (any(run_type.eq.(/QS_RUN,QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
!       inquire(file=trim(param%wenv%working_directory)//'/'//trim(param%global%project)//'-RESTART.wfn',exist=ex)
!       if (ex) &
        call system_command('cp '//trim(param%wenv%working_directory)//'/'//trim(param%global%project)//'-RESTART.wfn '//trim(param%wenv%wfn_file))
    endif

  ! remove all unnecessary files
    fin_command = 'rm -rf '//trim(param%wenv%working_directory)
    call system_command(fin_command,status=status)

    call finalise(param)    

    call print('CP2K finished. Go back to main program.')
    call system_timer('go_cp2k_end')

  end subroutine go_cp2k

  !!!!! ============================================================================================ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! \/ \/ \/ \/ \/ \/ -----------------  QM list operations  ----------------- \/ \/ \/ \/ \/ \/ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! ============================================================================================ !!!!!

  function extend_qmlist(my_atoms,R_inner,R_outer) result(list_changed)

  type(Atoms),           intent(inout) :: my_atoms
  real(dp),    optional, intent(in)    :: R_inner, R_outer  !! to select hysteretic quantum region

  real(dp)    :: inner, outer
  type(Table) :: inner_buffer, outer_buffer
  type(Table) :: core, ext_qmlist, old_ext_qmlist
  type(Atoms) :: fake, atoms_for_add_cut_hydrogens
  integer     :: qm_flag_index, i
  logical     :: list_changed,buffer_general,do_general

  ! checking the optional inputs
    if (present(R_inner)) then
       inner = R_inner
       if (present(R_outer)) then
          outer = R_outer
          call print('Updating QM_flag properties according to R_inner = '//round(R_inner,3)//', R_outer = '//round(R_outer,3))
       else
          call print('Warning: extend_qmlist: R_inner present, but not R_outer! Use R_outer=R_inner!')
          outer = R_inner
       endif
    else
       if (present(R_outer)) then
          call print('Warning: extend_qmlist: R_outer present, but not R_inner! Use R_inner=R_outer!')
          inner = R_outer
          outer = R_outer
       else
          call print('Warning: extend_qmlist: R_inner and R_outer not present! qm_list will not be extended!')
          return
       endif
    endif

    call get_qm_list_int(my_atoms,1,core)
    call get_qm_list_int_rec(my_atoms,2,ext_qmlist, do_recursive=.true.)
    call get_qm_list_int_rec(my_atoms,2,old_ext_qmlist, do_recursive=.true.)

!TO BE OPTIMIZED:
  ! extend the qm region > ext_qmlist
  !this should go without fake atoms and with use_avgpos! but not the cluster making, only the add_cut_hydrogen!
    fake=my_atoms
!    cutoff_old = my_atoms%cutoff
    call set_cutoff_minimum(fake,outer*(maxval(ElementCovRad(fake%Z(1:fake%N)))/ElementCovRad(8))) ! to find everything like O
    call calc_connect(fake)

    if (get_value(my_atoms%params,'Buffer_general',do_general)) then
        call print('Found Buffer_general in atoms%params'//do_general)
        buffer_general=do_general
    else
        buffer_general=.false.
    endif

    if (buffer_general) then
       call print('Not element specific buffer selection')
       call construct_buffer(fake,core,inner,inner_buffer,use_avgpos=.false.,verbosity=NORMAL)
       !call print('Atoms in inner buffer:')
       !do i=1,inner_buffer%N
       !   call print(inner_buffer%int(1,i))
       !enddo
       call construct_buffer(fake,core,outer,outer_buffer,use_avgpos=.false.,verbosity=NORMAL)
       !call print('Atoms in outer buffer:')
       !do i=1,outer_buffer%N
       !   call print(outer_buffer%int(1,i))
       !enddo
       !call print('Atoms in the quantum region before updating:')
       !do i=1,ext_qmlist%N
       !   call print(ext_qmlist%int(1,i))
       !enddo
    else
       call print('Element specific buffer selection, only heavy atoms count')
       call construct_buffer_RADIUS(fake,core,inner,inner_buffer,use_avgpos=.false.,verbosity=NORMAL)
       !call print('Atoms in inner buffer:')
       !do i=1,inner_buffer%N
       !   call print(inner_buffer%int(1,i))
       !enddo
       call construct_buffer_RADIUS(fake,core,outer,outer_buffer,use_avgpos=.false.,verbosity=NORMAL)
       !call print('Atoms in outer buffer:')
       !do i=1,outer_buffer%N
       !   call print(outer_buffer%int(1,i))
       !enddo
       !call print('Atoms in the quantum region before updating:')
       !do i=1,ext_qmlist%N
       !   call print(ext_qmlist%int(1,i))
       !enddo
    endif
    call select_hysteretic_quantum_region(fake,inner_buffer,outer_buffer,ext_qmlist,verbosity=NORMAL)
    !call print('Atoms in hysteretic quantum region:')
    !do i=1,ext_qmlist%N
    !   call print(ext_qmlist%int(1,i))
    !enddo
    call finalise(fake)
    call finalise(inner_buffer)
    call finalise(outer_buffer)

   ! add cut hydrogens, according to avgpos
    atoms_for_add_cut_hydrogens = my_atoms
    atoms_for_add_cut_hydrogens%oldpos = my_atoms%avgpos
    atoms_for_add_cut_hydrogens%avgpos = my_atoms%avgpos
    atoms_for_add_cut_hydrogens%pos = my_atoms%avgpos
    call set_cutoff(atoms_for_add_cut_hydrogens,0._dp)
    call calc_connect(atoms_for_add_cut_hydrogens)
    call add_cut_hydrogens(atoms_for_add_cut_hydrogens,ext_qmlist)
    !call print('Atoms in hysteretic quantum region after adding the cut hydrogens:')
    !do i=1,ext_qmlist%N
    !   call print(ext_qmlist%int(1,i))
    !enddo
    call finalise(atoms_for_add_cut_hydrogens)

   ! check changes in QM list and set the new QM list
    list_changed = check_qmlist_change(old_qmlist=old_ext_qmlist,new_qmlist=ext_qmlist)

   ! update QM_flag of my_atoms
    qm_flag_index = get_property(my_atoms,'QM_flag')
    my_atoms%data%int(qm_flag_index,1:my_atoms%N) = 0
    my_atoms%data%int(qm_flag_index,int_part(ext_qmlist,1)) = 2
    my_atoms%data%int(qm_flag_index,int_part(core,1)) = 1

    call finalise(core)
    call finalise(ext_qmlist)
    call finalise(old_ext_qmlist)

  end function extend_qmlist

  subroutine construct_buffer_RADIUS(my_atoms,core,radius,buffer,use_avgpos,verbosity)

    type(Atoms),       intent(in)  :: my_atoms
    real(dp),          intent(in)  :: radius
    type(Table),       intent(in)  :: core
    type(Table),       intent(out) :: buffer
    logical, optional, intent(in)  :: use_avgpos
    integer, optional, intent(in)  :: verbosity

    integer :: i,j
    logical :: delete_it

    call construct_buffer(my_atoms,core,(radius*maxval(ElementCovRad(my_atoms%Z(1:my_atoms%N)))/ElementCovRad(8)),buffer,use_avgpos=use_avgpos,verbosity=verbosity)
    !call print('Atoms in buffer:')
    !do i=1,buffer%N
    !   call print(buffer%int(1,i))
    !enddo

    i = 1
    do while (i.le.buffer%N)
       delete_it = .true.
       if (any(buffer%int(1,i).eq.core%int(1,1:core%N))) then ! keep core
          delete_it = .false.
       else
          if (my_atoms%Z(buffer%int(1,i)).ne.1) then ! delete Hs from buffer
             do j=1,core%N
                if (my_atoms%Z(core%int(1,j)).eq.1) cycle ! only consider nonH -- nonH distances
                if (distance_min_image(my_atoms,buffer%int(1,i),core%int(1,j)).le. &
                   radius / (ElementCovrad(8) + ElementCovRad(8)) * (ElementCovRad(my_atoms%Z(buffer%int(1,i))) + ElementCovRad(my_atoms%Z(core%int(1,j))))) then
                   delete_it = .false.
                   !call print('Atom '//ElementName(my_atoms%Z(buffer%int(1,i)))//buffer%int(1,i)//' is within element specific buffer distance '//round(distance_min_image(my_atoms,buffer%int(1,i),core%int(1,j)),3)//' < '//round(radius / (ElementCovrad(8) + ElementCovRad(8)) * (ElementCovRad(my_atoms%Z(buffer%int(1,i))) + ElementCovRad(my_atoms%Z(core%int(1,j)))),3))
                endif
             enddo
          endif
       endif
       if (delete_it) then
          !call print('delete atom '//ElementName(my_atoms%Z(buffer%int(1,i)))//buffer%int(1,i)//'from buffer')
          call delete(buffer,i)
          i = i - 1
       endif
       i = i + 1
    enddo

  end subroutine construct_buffer_RADIUS

  subroutine get_qm_list_int(my_atoms,qmflag,qmlist)

    type(Atoms), intent(in)  :: my_atoms
    integer,     intent(in)  :: qmflag
    type(Table), intent(out) :: qmlist
  
    integer              :: i
    integer              :: qm_flag_index

    qm_flag_index = get_property(my_atoms,'QM_flag')

    call initialise(qmlist,4,0,0,0,0)      !1 int, 0 reals, 0 str, 0 log, num_qm_atoms entries
    do i=1,my_atoms%N
       if (my_atoms%data%int(qm_flag_index,i).eq.qmflag) &
          call append(qmlist,(/i,0,0,0/))
    enddo

    if (qmlist%N.eq.0) call print('Empty QM list with QM_flag '//qmflag)

  end subroutine get_qm_list_int

  subroutine get_qm_list_int_rec(my_atoms,qmflag,qmlist,do_recursive)

    type(Atoms), intent(in)  :: my_atoms
    integer,     intent(in)  :: qmflag
    type(Table), intent(out) :: qmlist
    logical, intent(in) :: do_recursive
  
    integer              :: i
    integer              :: qm_flag_index
    logical              :: my_do_recursive

!    my_do_recursive = optional_default(.true.,do_recursive)
    qm_flag_index = get_property(my_atoms,'QM_flag')

    call initialise(qmlist,4,0,0,0,0)      !1 int, 0 reals, 0 str, 0 log, num_qm_atoms entries
!    if (my_do_recursive) then
    if (do_recursive) then
       do i=1,my_atoms%N
          if (my_atoms%data%int(qm_flag_index,i).gt.0.and. &
             my_atoms%data%int(qm_flag_index,i).le.qmflag) &
               call append(qmlist,(/i,0,0,0/))
       enddo
    else
       do i=1,my_atoms%N
          if (my_atoms%data%int(qm_flag_index,i).eq.qmflag) &
             call append(qmlist,(/i,0,0,0/))
       enddo
    endif

    if (qmlist%N.eq.0) call print('Empty QM list with QM_flag '//qmflag)

  end subroutine get_qm_list_int_rec

  subroutine get_qm_list_array(my_atoms,qmflags,qmlist)

    type(Atoms), intent(in)  :: my_atoms
    integer,     intent(in)  :: qmflags(:)
    type(Table), intent(out) :: qmlist
  
    integer              :: i
    integer              :: qm_flag_index

    qm_flag_index = get_property(my_atoms,'QM_flag')

    call initialise(qmlist,4,0,0,0,0)      !1 int, 0 reals, 0 str, 0 log, num_qm_atoms entries
    do i=1,my_atoms%N
       if (any(my_atoms%data%int(qm_flag_index,i).eq.qmflags(1:size(qmflags)))) &
          call append(qmlist,(/i,0,0,0/))
    enddo

    if (qmlist%N.eq.0) call print('Empty QM list with QM_flag '//qmflags(1:size(qmflags)))

  end subroutine get_qm_list_array

  function get_property(my_atoms,prop) result(prop_index)

    type(Atoms), intent(in)      :: my_atoms
    character(len=*), intent(in) :: prop
  
    integer,dimension(3) :: pos_indices
    integer              :: prop_index

    if (get_value(my_atoms%properties,trim(prop),pos_indices)) then
       prop_index = pos_indices(2)
    else
       call system_abort('get_property: No '//trim(prop)//' property assigned to the Atoms object!')
    end if

  end function get_property


  !!!!! ============================================================================================ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! \/ \/ \/ \/ \/ \/ ----------  reading in the QM list from file  ---------- \/ \/ \/ \/ \/ \/ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! ============================================================================================ !!!!!

  subroutine read_qmlist(my_atoms,qmlistfilename,verbose)

    type(Atoms),       intent(inout) :: my_atoms
    character(*),      intent(in)    :: qmlistfilename
    logical, optional, intent(in)    :: verbose

    type(table)                      :: qm_list
    type(Inoutput)                   :: qmlistfile
    integer                          :: n,num_qm_atoms,qmatom,status
    character(80)                    :: title,testline
    logical                          :: my_verbose
    integer                          :: qm_flag_index
    character(20), dimension(10)     :: fields
    integer                          :: num_fields


    my_verbose = .false.
    if (present(verbose)) my_verbose = verbose

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
    if (num_qm_atoms.gt.my_atoms%N) call print('WARNING! read_qmlist: more QM atoms then atoms in the atoms object, possible redundant QM list file')
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

    call add_property(my_atoms,'QM_flag',0)
    qm_flag_index = get_property(my_atoms,'QM_flag')
    my_atoms%data%int(qm_flag_index,1:my_atoms%N) = 0
    my_atoms%data%int(qm_flag_index,int_part(qm_list,1)) = 1

    if (my_verbose) call print('Finished. '//qm_list%N//' QM atoms have been read successfully.')
    call finalise(qm_list)

  end subroutine read_qmlist


  !!!!! ============================================================================================ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! \/ \/ \/ \/ \/ \/ --------  create the CHARMM params to run CP2K  -------- \/ \/ \/ \/ \/ \/ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! ============================================================================================ !!!!!

  !calculates number of bonds generated during connectivity generation
  !bonds2 checks if symmetric
  function num_of_bonds(at) result(bonds)

    type(atoms),intent(in) :: at
    integer                :: bonds
    integer                :: i,bonds2

    bonds  = 0
    bonds2 = 0
    
    do i=1,at%N
       if (associated(at%connect%neighbour1(i)%t)) bonds  = bonds  + at%connect%neighbour1(i)%t%N
       if (associated(at%connect%neighbour2(i)%t)) bonds2 = bonds2 + at%connect%neighbour2(i)%t%N
    enddo

    if (bonds.ne.bonds2) call system_abort('num_of_bonds: connectivities of atoms are not symmetric.')

  end function num_of_bonds

!prints the number of neighbours, if not: 1 for H, 2 for O.
!the two H of a flexible water molecule can get too close => correct it!
  subroutine check_neighbour_numbers(at)

    type(atoms),intent(in) :: at
    integer                :: i, n_neigh

    do i=1,at%N

       n_neigh = atoms_n_neighbours(at, i)
       if ( (at%Z(i).eq.1) .AND. (n_neigh.ne.1) ) then
!          do n=1,atoms_n_neighbours(at,i)
!             j=atoms_neighbour(at,i,n)
!          enddo
          call print('Warning: '//ElementName(at%Z(i))//' atom '//i//' position '//round(at%pos(1,i),5)//' '//round(at%pos(2,i),5)//' '//round(at%pos(3,i),5)//' has unusual number of neighbours: '//(n_neigh))
       endif
       if ( (at%Z(i).eq.8) .AND. (n_neigh.ne.2) ) &
          call print('Warning: '//ElementName(at%Z(i))//' atom '//i//' position '//round(at%pos(1,i),5)//' '//round(at%pos(2,i),5)//' '//round(at%pos(3,i),5)//' has unusual number of neighbours: '//(n_neigh))

       if ( .not.any(at%Z(i).eq.(/1,8/)) ) &
          call print('Info: '//ElementName(at%Z(i))//' atom '//i//' position '//round(at%pos(1,i),5)//' '//round(at%pos(2,i),5)//' '//round(at%pos(3,i),5)//' has number of neighbours: '//(n_neigh))

    enddo

  end subroutine check_neighbour_numbers

  function check_qmlist_change(old_qmlist,new_qmlist) result(list_changed)

    type(Table), intent(in) :: old_qmlist, new_qmlist
    integer :: i
    logical :: list_changed

    list_changed = .false.
    if (old_qmlist%N.ne.new_qmlist%N) then
       call print ('fitlist has changed: new number of qm atoms: '//new_qmlist%N//' instead of '//old_qmlist%N)
       list_changed = .true.
    else
       if (any(int_part(old_qmlist,2).ne.int_part(new_qmlist,2))) then
           do i=1,old_qmlist%N
              if (.not.find_in_array(int_part(old_qmlist,1),(new_qmlist%int(1,i))).gt.0) then
                 call print('fitlist has changed: atom '//new_qmlist%int(1,i)//' has entered QM zone')
                 list_changed = .true.
              endif
              if (.not.find_in_array(int_part(new_qmlist,1),(old_qmlist%int(1,i))).gt.0) then
                 call print('fitlist has changed: atom '//old_qmlist%int(1,i)//' has left QM zone')
                 list_changed = .true.
              endif
           enddo
       endif
    endif

  end function check_qmlist_change

  !!!!! ============================================================================================ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! \/ \/ \/ \/ \/ \/ ------------  CHARMM  atom type generation  ------------ \/ \/ \/ \/ \/ \/ !!!!!
  !!!!! \/ \/ \/ \/ \/ \/ ------------    with average coordinates    ------------ \/ \/ \/ \/ \/ \/ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! ============================================================================================ !!!!!
! topology calculation using average coordinates
  subroutine calc_topology(at,do_CHARMM,intrares_impropers)

    type(Atoms),       intent(inout) :: at
    logical, optional, intent(in)    :: do_CHARMM
    type(Table), optional, intent(out) :: intrares_impropers
    type(Atoms)                      :: atoms_for_find_motif
    integer                          :: pos_indices(3), &
                                        atom_type_index, &
                                        atom_res_name_index, &
                                        atom_mol_name_index, &
                                        atom_res_number_index, &
                                        atom_charge_index, &
                                        atom_type_index2, &
                                        atom_res_name_index2, &
                                        atom_mol_name_index2, &
                                        atom_res_number_index2, &
                                        atom_charge_index2
    logical :: my_do_CHARMM

  call system_timer('calc_topology')
   my_do_CHARMM = optional_default(.true.,do_CHARMM)
   !this should be only the 3 calls with use_avgpos!
    call print('Creating CHARMM format...')

    if (.not.get_value(at%properties,'oldpos',pos_indices)) &
       call system_abort('calc_topology: atoms object does not have oldpos property')
    if (.not.get_value(at%properties,'avgpos',pos_indices)) &
       call system_abort('calc_topology: atoms object does not have avgpos property')

    atoms_for_find_motif = at
    atoms_for_find_motif%oldpos = at%avgpos
    atoms_for_find_motif%avgpos = at%avgpos
    atoms_for_find_motif%pos = at%avgpos
!     call set_cutoff(atoms_for_find_motif,DEFAULT_NNEIGHTOL)
!    call set_cutoff(atoms_for_find_motif,at%cutoff)
call set_cutoff(atoms_for_find_motif, 0._dp)
    call calc_connect(atoms_for_find_motif)
    call delete_metal_connects(atoms_for_find_motif)
    call map_into_cell(atoms_for_find_motif)
    call calc_dists(atoms_for_find_motif)
    if (my_do_CHARMM.and.present(intrares_impropers)) then
       call create_CHARMM(atoms_for_find_motif,do_CHARMM=my_do_CHARMM,intrares_impropers=intrares_impropers)
    else
       call create_CHARMM(atoms_for_find_motif,do_CHARMM=my_do_CHARMM)
    endif

    atoms_for_find_motif%pos = atoms_for_find_motif%oldpos

   ! copy data to at
    call add_property(at,'atom_type',repeat(' ',TABLE_STRING_LENGTH))
    call add_property(at,'atom_res_name',repeat(' ',TABLE_STRING_LENGTH))
    call add_property(at,'atom_mol_name',repeat(' ',TABLE_STRING_LENGTH))
    call add_property(at,'atom_res_number',0)
    call add_property(at,'atom_charge',0._dp)

    atom_type_index = get_property(at,'atom_type')
    atom_type_index2 = get_property(atoms_for_find_motif,'atom_type')
    atom_res_name_index = get_property(at,'atom_res_name')
    atom_res_name_index2 = get_property(atoms_for_find_motif,'atom_res_name')
    atom_mol_name_index = get_property(at,'atom_mol_name')
    atom_mol_name_index2 = get_property(atoms_for_find_motif,'atom_mol_name')
    atom_res_number_index = get_property(at,'atom_res_number')
    atom_res_number_index2 = get_property(atoms_for_find_motif,'atom_res_number')
    atom_charge_index = get_property(at,'atom_charge')
    atom_charge_index2 = get_property(atoms_for_find_motif,'atom_charge')

    at%data%str(atom_type_index,1:at%N) = atoms_for_find_motif%data%str(atom_type_index2,1:at%N)
    at%data%str(atom_res_name_index,1:at%N) = atoms_for_find_motif%data%str(atom_res_name_index2,1:at%N)
    at%data%str(atom_mol_name_index,1:at%N) = atoms_for_find_motif%data%str(atom_mol_name_index2,1:at%N)
    at%data%int(atom_res_number_index,1:at%N) = atoms_for_find_motif%data%int(atom_res_number_index2,1:at%N)
    at%data%real(atom_charge_index,1:at%N) = atoms_for_find_motif%data%real(atom_charge_index2,1:at%N)

    if (any(at%data%int(atom_res_number_index,1:at%N).le.0)) &
       call system_abort('calc_topology: atom_res_number is not >0 for every atom')
    if (any(at%data%str(atom_type_index,1:at%N).eq.'X')) &
       call system_abort('calc_topology: atom_type is not saved for at least one atom')
    if (any(at%data%str(atom_res_name_index,1:at%N).eq.'X')) &
       call system_abort('calc_topology: atom_res_name is not saved for at least one atom')
    call print('CHARMM parameters added and stored in atoms object')
    call finalise(atoms_for_find_motif)
  call system_timer('calc_topology')

  end subroutine calc_topology

   ! remove bonds for metal ions - everything but H, C, N, O, Si, P, S, Cl
  subroutine delete_metal_connects(my_atoms)

    type(Atoms), intent(inout) :: my_atoms
    integer :: i, j

    do i = 1, my_atoms%N
       if (any(my_atoms%Z(i).eq.(/1,6,7,8,14,15,16,17/))) cycle
!       call print('found metal atom '//i)
       do j = 1, my_atoms%N
          if (i.ne.j) call delete_bond(my_atoms, i, j)
       enddo
    enddo

  end subroutine

   ! removes a bond between i and j, if the bond is present
  subroutine delete_bond(my_atoms, i, j)
    type(Atoms), intent(inout) :: my_atoms
    integer, intent(in) :: i, j

    integer :: ii, jj
!!    integer :: kk, k, ll, change
!!    integer, allocatable :: bond_table(:,:)

    if (i.eq.j) return

    if (i.lt.j) then
       ii = i
       jj = j
    else
       jj = i
       ii = j
    endif

    call remove_bond(my_atoms%connect, ii, jj)

!!     allocate (bond_table(4,my_atoms%connect%neighbour1(ii)%N))
!!     bond_table = int_part(my_atoms%connect%neighbour1(ii))
!!     kk = find_in_array(bond_table(1,:),jj)
!!     if (kk.gt.0) then
!! !       call print('found bond to delete for atoms '//ii//' '//jj)
!!        call delete(my_atoms%connect%neighbour1(ii),kk)
!!     endif
!!     deallocate(bond_table)
!! 
!!    ! if I delete a bond from neighbour1, I should update in neighbour2 that
!!    ! it's the $(n-1)$-th (or the last is now takes the place of the deleted) neighbour from now on
!!     if (kk.gt.0) then
!!        do k = kk, my_atoms%connect%neighbour1(ii)%N     ! k-th neighbour of ii
!!           ll = my_atoms%connect%neighbour1(ii)%int(1,k)     ! is ll
!!           allocate (bond_table(4,my_atoms%connect%neighbour2(ll)%N))
!!           bond_table = int_part(my_atoms%connect%neighbour2(ll))
!!           change = find_in_array(bond_table(1,:),ii)    ! ll has an account for ii
!!           if (change.eq.0) call system_abort('Found nonsymmetrical connectivity.')
!!           my_atoms%connect%neighbour2(ll)%int(2,change) = k ! this account is updated now
!!           deallocate (bond_table)
!!        enddo
!!     endif
!! 
!!     allocate (bond_table(4,my_atoms%connect%neighbour2(jj)%N))
!!     bond_table = int_part(my_atoms%connect%neighbour2(jj))
!!     kk = find_in_array(bond_table(1,:),ii)
!!     if (kk.gt.0) then
!! !       call print('found bond to delete for atoms '//jj//' '//ii)
!!        call delete(my_atoms%connect%neighbour2(jj),kk)
!!     endif
!!     deallocate(bond_table)

  end subroutine


  !!!!! ============================================================================================ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! \/ \/ \/ \/ \/ \/ ------------  input file printing for CP2K  ------------ \/ \/ \/ \/ \/ \/ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! ============================================================================================ !!!!!

  subroutine write_cp2k_input_files(my_atoms,qm_list,param,run_type,PSF_print)

    type(Atoms),       intent(in) :: my_atoms
    type(Table),       intent(in) :: qm_list
    type(param_cp2k),  intent(in) :: param
    integer,           intent(in) :: run_type, &
                                     PSF_print

    call system_timer('write_inputs')

    if (.not.any(run_type.eq.(/QS_RUN,MM_RUN,QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) &
       call system_abort('write_cp2k_input_files: Run type is not recognized: '//run_type)
    if (.not.any(PSF_print.eq.(/NO_PSF,CP2K_PRINT_AND_SAVE,DRIVER_PRINT_AND_SAVE,USE_EXISTING_PSF/))) &
       call system_abort('write_cp2k_input_files: PSF print is not recognized: '//PSF_print)

    call print('Writing the CP2K input file(s)')
    if (run_type.ne.QS_RUN) then
       call write_cp2k_pdb_file(my_atoms,trim(param%wenv%working_directory)//'/'//trim(param%wenv%pdb_file),run_type)
       if (PSF_print.eq.DRIVER_PRINT_AND_SAVE) then
          call write_psf_file(my_atoms,param%wenv%psf_file,run_type)
       endif
    endif
    call write_cp2k_input_file(my_atoms,qm_list,param,run_type,PSF_print)
    call print ('Finished. Input files are ready for the run.')
    if (PSF_print.eq.USE_EXISTING_PSF) call print('Use already existing PSF file '//param%wenv%psf_file)

    call system_timer('write_inputs')

  end subroutine write_cp2k_input_files

  !writes modified PDB format for accurate coordinates and charges, for CP2K
  !use CHARGE_EXTENDED keyword i.e. reads charges from the last column of PDB file 
  !ATOM      1  CT3 ALA A   1       0.76700000   0.80100000  13.31100000  0.00  0.00     ALA   C  -0.27
  !ATOM      2   HA ALA A   1       0.07400000  -0.06000000  13.18800000  0.00  0.00     ALA   H   0.09
  !ATOM      3   HA ALA A   1       0.17600000   1.74100000  13.29800000  0.00  0.00     ALA   H   0.09
  subroutine write_cp2k_pdb_file(at,pdb_file,run_type)

    character(len=*),  intent(in)  :: pdb_file
    type(atoms),       intent(in)  :: at
    integer,           intent(in)  :: run_type

!    character(*), parameter  :: pdb_format = '(a6,i5,1x,a4,1x,a4,1x,i4,1x,3x,3f13.8,2f6.2,10x,a2,2x,f6.4)'
    character(*), parameter  :: pdb_format = '(a6,i5,1x,a4,1x,a4,i5,1x,3x,3f13.8,2f6.2,5x,a4,1x,a2,2x,f7.4)'

    type(Inoutput)           :: pdb
    character(103)           :: sor
    integer                  :: mm
    character(4)             :: QM_prefix_atom_mol_name
    integer                  :: qm_flag_index, &
                                atom_type_index, &
                                atom_res_name_index, &
                                atom_mol_name_index, &
                                atom_res_number_index, &
                                atom_charge_index

  !Brookhaven PDB format
  !    character(86)                        :: sor
  !    character(*), parameter               :: pdb_format = '(a6,i5,1x,a4,1x,a4,1x,i4,1x,3x,3f8.3,2f6.2,10x,a2,2x,f7.3)'

  !CP2K modified PDB format 
  !       sor(1:6)   = 'ATOM  '
  !       sor(7:11)  = mm
  !       sor(13:16) = this%atom_type(mm)
  !       sor(18:21) = this%res_name(mm)
  !!      sor(22:22) = ' A'                   !these two are now
  !       sor(23:26) = residue_number(mm)     !   merged to handle >9999 residues
  !       sor(31:43) = at%pos(1,mm)
  !       sor(44:56) = at%pos(2,mm)
  !       sor(57:69) = at%pos(3,mm)
  !       sor(70:75) = '  0.00'
  !       sor(76:81) = '  0.00'
  !!      sor(87:90) = 'MOL1'
  !       sor(92:93) = ElementName(at%Z(mm))
  !       sor(96:)   = this%atom_charge(mm)

    call initialise(pdb,trim(pdb_file),action=OUTPUT)
    call print('   PDB file: '//trim(pdb%filename))
!    call print('REMARK'//at%N,file=pdb)

    if (any(run_type.eq.(/QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
       qm_flag_index = get_property(at,'QM_flag')
    endif
    atom_type_index = get_property(at,'atom_type')
    atom_res_name_index = get_property(at,'atom_res_name')
    atom_mol_name_index = get_property(at,'atom_mol_name')
    atom_res_number_index = get_property(at,'atom_res_number')
    atom_charge_index = get_property(at,'atom_charge')

    do mm=1,at%N
      ! CP2K needs different name for QM molecules, if use isolated atoms
       sor = ''
       QM_prefix_atom_mol_name = ''
       QM_prefix_atom_mol_name = trim(at%data%str(atom_mol_name_index,mm))
       if (any(run_type.eq.(/QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
          if (at%data%int(qm_flag_index,mm).ge.QMMM_RUN_CORE .and. at%data%int(qm_flag_index,mm).le.run_type) then
             QM_prefix_atom_mol_name = 'QM'//trim(at%data%str(atom_mol_name_index,mm))
!             call print('QM molecule '//QM_prefix_atom_mol_name)
          endif
       endif
!             call print('molecule '//QM_prefix_atom_mol_name)
!       call print('writing PDB file: atom type '//at%data%str(atom_type_index,mm))
       write(sor,pdb_format) 'ATOM  ',mm,at%data%str(atom_type_index,mm),at%data%str(atom_res_name_index,mm),at%data%int(atom_res_number_index,mm), &
                             at%pos(1:3,mm),0._dp,0._dp,QM_prefix_atom_mol_name,ElementName(at%Z(mm)),at%data%real(atom_charge_index,mm)
!                             at%pos(1:3,mm),0._dp,0._dp,this%atom_res_name(mm),ElementName(at%Z(mm)),this%atom_charge(mm)
       call print(sor,file=pdb)
    enddo
    call print('END',file=pdb)
    call finalise(pdb)

  end subroutine write_cp2k_pdb_file


  subroutine write_cp2k_input_file(my_atoms,qm_list,param,run_type,PSF_print)

    type(Atoms),       intent(in) :: my_atoms
    type(Table),       intent(in) :: qm_list
    type(param_cp2k),  intent(in) :: param
    integer,           intent(in) :: run_type, &
                                     PSF_print

    type(Inoutput)                   :: input_file 
    integer,allocatable,dimension(:) :: list   !in order to separate QM atom indices for H and O
    integer                          :: i,j,counter
    type(Table)                      :: qm_list_check
    integer                          :: qm_flag_index, &
                                        atom_type_index, &
                                        atom_res_name_index, &
                                        atom_charge_index
    logical                          :: ex, lsd
    integer                          :: charge
    character(len=value_len)         :: basis_set, potential
    integer                          :: stat
    type(inoutput)                   :: dft_in_io, cell_in_io
    character(len=1024)              :: line

    if (.not.any(run_type.eq.(/QS_RUN,MM_RUN,QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) &
       call system_abort('write_cp2k_input_file: Run type is not recognized: '//run_type)

    if (any(run_type.eq.(/MM_RUN,QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
       if (any(run_type.eq.(/QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
          qm_flag_index = get_property(my_atoms,'QM_flag')
       endif
       atom_type_index = get_property(my_atoms,'atom_type')
       atom_res_name_index = get_property(my_atoms,'atom_res_name')
       atom_charge_index = get_property(my_atoms,'atom_charge')
    endif

    if (any(run_type.eq.(/QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
       call get_qm_list_int_rec(my_atoms,run_type,qm_list_check, do_recursive=.true.)
       if (qm_list_check%N.ne.qm_list%N) call system_abort ('qm_list is different from Atoms QM_flag!!')
       do i = 1, qm_list_check%N
          if (.not.any(qm_list_check%int(1,i).eq.qm_list%int(1,1:qm_list%N))) call system_abort ('qm_list is different from Atoms QM_flag!!')
       enddo
       do i = 1, qm_list%N
          if (.not.any(qm_list%int(1,i).eq.qm_list_check%int(1,1:qm_list_check%N))) call system_abort ('qm_list is different from Atoms QM_flag!!')
       enddo
    endif

!    call print_title('Writing CP2K input file')

    call initialise(input_file,trim(param%wenv%working_directory)//'/'//trim(param%wenv%cp2k_input_filename),action=OUTPUT)
    call print('   inp file: '//trim(input_file%filename))

  !  input_file%default_real_precision=16

    call print('&FORCE_EVAL',file=input_file)
    if (run_type.eq.MM_RUN) then
       call print('  METHOD Fist',file=input_file)
    else
       if (run_type.eq.QS_RUN) then
        call print('  METHOD Quickstep',file=input_file)
       else !QMMM_RUN_CORE or QMMM_RUN_EXTENDED
        call print('  METHOD QMMM',file=input_file)
       endif
    endif

    if (any(run_type.eq.(/QS_RUN,QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
       call print('  &DFT',file=input_file)
       if (len(trim(param%dft%file)) > 0) then
          call initialise(dft_in_io,trim(param%dft%file), INPUT)
          stat = 0
          do while (stat == 0)
            line = read_line(dft_in_io, stat)
            if (stat == 0) call print(trim(line), file=input_file)
          end do
          call finalise(dft_in_io)
          if (run_type.eq.QS_RUN) then
             call print('    WFN_RESTART_FILE_NAME ../'//trim(param%wenv%wfn_file),file=input_file)
          else
             if (param%qmmm%reuse_wfn) then
                call print('    WFN_RESTART_FILE_NAME ../'//trim(param%wenv%wfn_file),file=input_file)
             else
                call print('#    WFN_RESTART_FILE_NAME ../'//trim(param%wenv%wfn_file),file=input_file)
             endif
          endif
          if (run_type == QS_RUN) then
            ! check if the system is charged, given Charge=... in the comment line
             if (get_value(my_atoms%params,'Charge',charge)) then
               call print('    CHARGE '//charge,file=input_file)
             endif
             lsd = .false.
             if (get_value(my_atoms%params,'LSD',lsd)) then
               if (lsd) call print('    LSD ',file=input_file)
             endif
          else
             if (mod(nint(sum(my_atoms%data%real(atom_charge_index,int_part(qm_list,1)))),2).ne.0) call print('    LSD',file=input_file)
             if (nint(sum(my_atoms%data%real(atom_charge_index,int_part(qm_list,1)))).ne.0) &
                call print('    CHARGE '//nint(sum(my_atoms%data%real(atom_charge_index,int_part(qm_list,1)))),file=input_file)
          endif
       else
         ! call print('    &PRINT',file=input_file)
         ! call print('      &LOCALIZATION',file=input_file)
         ! call print('        &LOCALIZE',file=input_file)
         ! call print('        &END LOCALIZE',file=input_file)
         ! call print('      &END LOCALIZATION',file=input_file)
         ! call print('    &END PRINT',file=input_file)

          call print('    BASIS_SET_FILE_NAME ../BASIS_SET',file=input_file)
          call print('    POTENTIAL_FILE_NAME ../POTENTIAL',file=input_file)
          if (run_type.eq.QS_RUN) then
             call print('    WFN_RESTART_FILE_NAME ../'//trim(param%wenv%wfn_file),file=input_file)
          else
             if (param%qmmm%reuse_wfn) then
                call print('    WFN_RESTART_FILE_NAME ../'//trim(param%wenv%wfn_file),file=input_file)
             else
                call print('#    WFN_RESTART_FILE_NAME ../'//trim(param%wenv%wfn_file),file=input_file)
             endif
          endif
          call print('    &MGRID',file=input_file)
          call print('      COMMENSURATE',file=input_file)
          call print('      CUTOFF '//round(param%dft%mgrid_cutoff,1),file=input_file)   !some cutoff !Ry
          call print('    NGRIDS 5',file=input_file)
          call print('    &END MGRID',file=input_file)
          call print('    &QS',file=input_file)
          call print('      EXTRAPOLATION PS',file=input_file)
          call print('      EXTRAPOLATION_ORDER 1',file=input_file)
          call print('    &END QS',file=input_file)
          call print('    &SCF',file=input_file)
          if (run_type.eq.QS_RUN) then
             call print('      SCF_GUESS RESTART',file=input_file)
          else
  !           if (param%qmmm%wfn_reuse) then
  !              call print('      SCF_GUESS RESTART',file=input_file)
  !           else
                call print('      SCF_GUESS '//trim(param%dft%scf_guess),file=input_file)
  !           endif
          endif
          call print('      EPS_SCF 1.0E-6',file=input_file)
          call print('      MAX_SCF 200',file=input_file)
          call print('      &OUTER_SCF',file=input_file)
          call print('         EPS_SCF 1.0E-6',file=input_file)
          call print('         MAX_SCF 200',file=input_file)
          call print('      &END',file=input_file)
          call print('      &OT',file=input_file)
          call print('         MINIMIZER CG',file=input_file)
          call print('         PRECONDITIONER FULL_ALL',file=input_file)
  !        call print('#         PRECONDITIONER SPARSE_DIAG',file=input_file)
  !        call print('#         ENERGY_GAP 0.001',file=input_file)
          call print('      &END OT',file=input_file)
          call print('    &END SCF',file=input_file)
          call print('    &XC',file=input_file)
          call print('      &XC_FUNCTIONAL '//trim(param%dft%xc_functional),file=input_file)
          call print('      &END XC_FUNCTIONAL',file=input_file)
          call print('    &END XC',file=input_file)
          ! in case of charged QM system use unrestricted KS and specify charge
          if (run_type.eq.QS_RUN) then
             call print('    &POISSON',file=input_file)
             call print('      &EWALD',file=input_file)
             call print('        EWALD_TYPE '//trim(param%mm_ewald%ewald_type),file=input_file)
             call print('        EWALD_ACCURACY 1.E-6',file=input_file)
             call print('        ALPHA '//round(param%mm_ewald%ewald_alpha,2),file=input_file)
             call print('        GMAX '//param%mm_ewald%ewald_gmax,file=input_file)
             call print('      &END EWALD',file=input_file)
             call print('      POISSON_SOLVER PERIODIC',file=input_file)
             call print('      PERIODIC XYZ',file=input_file)
             call print('      &MULTIPOLE',file=input_file)
             call print('        RCUT 10',file=input_file)
             call print('        EWALD_PRECISION 1.E-6',file=input_file)
             call print('        NGRIDS 50 50 50',file=input_file)
             call print('        &INTERPOLATOR',file=input_file)
  !           call print('#          AINT_PRECOND SPL3_NOPBC_AINT1',file=input_file)
  !           call print('#          PRECOND SPL3_NOPBC_PRECOND3',file=input_file)
  !           call print('#          EPS_X 1.E-15',file=input_file)
  !           call print('#          EPS_R 1.E-15',file=input_file)
  !           call print('#          MAX_ITER 200',file=input_file)
             call print('        &END INTERPOLATOR',file=input_file)
             call print('      &END MULTIPOLE',file=input_file)
             call print('    &END POISSON',file=input_file)
            ! check if the system is charged, given Charge=... in the comment line
             ex = .false.
             ex = get_value(my_atoms%params,'Charge',charge)
             if (ex) call print('    CHARGE '//charge,file=input_file)
             ex = .false.
             lsd = .false.
             ex = get_value(my_atoms%params,'LSD',lsd)
             if (ex .and. lsd) call print('    LSD ',file=input_file)
          else
             if (mod(nint(sum(my_atoms%data%real(atom_charge_index,int_part(qm_list,1)))),2).ne.0) call print('    LSD',file=input_file)
             if (nint(sum(my_atoms%data%real(atom_charge_index,int_part(qm_list,1)))).ne.0) &
                call print('    CHARGE '//nint(sum(my_atoms%data%real(atom_charge_index,int_part(qm_list,1)))),file=input_file)
          endif
        end if ! dft file
      call print('  &END DFT',file=input_file)
   endif

  if (any(run_type.eq.(/MM_RUN,QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
    call print('  &MM',file=input_file)
    call print('    &FORCEFIELD',file=input_file)
    call print('      PARMTYPE '//trim(param%mm_forcefield%parmtype),file=input_file)
    call print('      PARM_FILE_NAME ../'//trim(param%mm_forcefield%parm_file_name),file=input_file) ! charmm.pot
#ifdef HAVE_DANNY
    call print('      DANNY T',file=input_file)
    call print('      DANNY_CUTOFF 5.5',file=input_file)
    call print('      &SPLINE',file=input_file)
    call print('        EMAX_SPLINE 0.5',file=input_file)
    call print('        EMAX_ACCURACY 0.1',file=input_file)
    call print('      &END SPLINE',file=input_file)
#endif

!!for old CP2K, QM/MM
!    if (any(my_atoms%data%str(atom_type_index,1:my_atoms%N).eq.'CCL') .and. any(my_atoms%data%str(atom_type_index,1:my_atoms%N).eq.'CLA')) then
!       call print('      &SPLINE',file=input_file)
!       call print('        EMAX_SPLINE 20.',file=input_file)
!       call print('      &END SPLINE',file=input_file)
!    endif
!!for old CP2K, QM/MM

    if (any(my_atoms%data%str(atom_type_index,1:my_atoms%N).eq.'HFL') .or. any(my_atoms%data%str(atom_type_index,1:my_atoms%N).eq.'OFL')) then
       call print('      &BEND',file=input_file)
       call print('        ATOMS HFL OFL HFL',file=input_file)
       call print('        K 0.074',file=input_file)
       call print('        THETA0 1.87',file=input_file)
       call print('      &END BEND',file=input_file)
       call print('      &BOND',file=input_file)
       call print('        ATOMS OFL HFL',file=input_file)
       call print('        K 0.27',file=input_file)
       call print('        R0 1.84',file=input_file)
       call print('      &END BOND',file=input_file)
!       call print('#      &CHARGE',file=input_file)
!       call print('#        ATOM OFL',file=input_file)
!       call print('#        CHARGE -0.8476',file=input_file)
!       call print('#      &END CHARGE',file=input_file)
!       call print('#      &CHARGE',file=input_file)
!       call print('#        ATOM HFL',file=input_file)
!       call print('#        CHARGE 0.4238',file=input_file)
!       call print('#      &END CHARGE',file=input_file)
       call print('      &NONBONDED',file=input_file)
       call print('        &LENNARD-JONES',file=input_file)
       call print('          atoms OFL OFL',file=input_file)
       call print('          EPSILON 78.198',file=input_file)
       call print('          SIGMA 3.166',file=input_file)
       call print('          RCUT 11.4',file=input_file)
       call print('        &END LENNARD-JONES',file=input_file)
       call print('        &LENNARD-JONES',file=input_file)
       call print('          atoms OFL HFL',file=input_file)
       call print('          EPSILON 0.0',file=input_file)
       call print('          SIGMA 3.6705',file=input_file)
       call print('          RCUT 11.4',file=input_file)
       call print('        &END LENNARD-JONES',file=input_file)
       call print('        &LENNARD-JONES',file=input_file)
       call print('          atoms HFL HFL',file=input_file)
       call print('          EPSILON 0.0',file=input_file)
       call print('          SIGMA 3.30523',file=input_file)
       call print('          RCUT 11.4',file=input_file)
       call print('        &END LENNARD-JONES',file=input_file)
       call print('      &END NONBONDED',file=input_file)
    endif
    call print('    &END FORCEFIELD',file=input_file)
    call print('    &POISSON',file=input_file)
    call print('      &EWALD',file=input_file)
    call print('        EWALD_TYPE '//trim(param%mm_ewald%ewald_type),file=input_file)
    call print('        EWALD_ACCURACY 1.E-6',file=input_file)
    call print('        ALPHA '//round(param%mm_ewald%ewald_alpha,2),file=input_file)
    call print('        GMAX '//param%mm_ewald%ewald_gmax,file=input_file)
    call print('      &END EWALD',file=input_file)
    call print('    &END POISSON',file=input_file)
    call print('  &END MM',file=input_file)
  endif

    if (any(run_type.eq.(/QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
       call print('  &QMMM',file=input_file)
       call print('    MM_POTENTIAL_FILE_NAME ../MM_POTENTIAL',file=input_file)
       call print('    &CELL',file=input_file)
       call print('      ABC '//round(param%qmmm%qmmm_cell(1),6)//' '//round(param%qmmm%qmmm_cell(2),6)//' '//round(param%qmmm%qmmm_cell(3),6),file=input_file) !some cell size for the QM region
!       call print('      UNIT '//trim(param%qmmm%qmmm_cell_unit),file=input_file)   !in Angstroms
       call print('      PERIODIC XYZ',file=input_file)   !in Angstroms
       call print('    &END CELL',file=input_file)
       call print('    ECOUPL '//trim(param%qmmm%ecoupl),file=input_file) !QMMM electrostatic coupling, can be NONE,COULOMB,GAUSS,S-WAVE
     
!       call print('Writing decoupling section...')
     
       call print('    &PERIODIC',file=input_file)
!       call print('#      GMAX 1.',file=input_file)
!       call print('#      REPLICA -1',file=input_file)
!       call print('#      NGRIDS 50 50 50',file=input_file)
       call print('      &MULTIPOLE',file=input_file)
       call print('        RCUT 10',file=input_file)
       call print('        EWALD_PRECISION 1.E-6',file=input_file)
       call print('        NGRIDS 50 50 50',file=input_file)
       call print('        &INTERPOLATOR',file=input_file)
!       call print('#          AINT_PRECOND SPL3_NOPBC_AINT1',file=input_file)
!       call print('#          PRECOND SPL3_NOPBC_PRECOND3',file=input_file)
!       call print('#          EPS_X 1.E-15',file=input_file)
!       call print('#          EPS_R 1.E-15',file=input_file)
!       call print('#          MAX_ITER 200',file=input_file)
       call print('        &END INTERPOLATOR',file=input_file)
       call print('      &END MULTIPOLE',file=input_file)
       call print('    &END PERIODIC',file=input_file)

       !print non default COVALENT RADII: only 0.44, 0.78, 0.80 (=default) available, use 0.44 for H, 0.78 for O
       do i=1,my_atoms%N-1
          if (any(my_atoms%data%str(atom_type_index,i).eq.my_atoms%data%str(atom_type_index,i+1:my_atoms%N))) cycle
          if (any(my_atoms%Z(i).eq.(/1,8/))) then
             call print('    &MM_KIND '//trim(my_atoms%data%str(atom_type_index,i)),file=input_file)
             if (my_atoms%Z(i).eq.1) then
                call print('      RADIUS '//round(param%qmmm%radius_H,3),file=input_file) !radius of the atomic kind
             else
                call print('      RADIUS '//round(param%qmmm%radius_O,3),file=input_file) !radius of the atomic kind
             endif
             call print('    &END MM_KIND',file=input_file)
!          else
!             call print(trim(my_atoms%data%str(atom_type_index,i))//'is not H or O')
          endif
       enddo
       if (any(my_atoms%Z(my_atoms%N).eq.(/1,8/))) then
          call print('    &MM_KIND '//trim(my_atoms%data%str(atom_type_index,my_atoms%N)),file=input_file)
          if (my_atoms%Z(my_atoms%N).eq.1) then
             call print('      RADIUS '//round(param%qmmm%radius_H,3),file=input_file) !radius of the atomic kind
          else
             call print('      RADIUS '//round(param%qmmm%radius_O,3),file=input_file) !radius of the atomic kind
          endif
          call print('    &END MM_KIND',file=input_file)
       endif
       call print('#    otherwise use default: 0.800',file=input_file)
     
!       call print('Writing QM list...')
       allocate(list(qm_list%N))
       list=int_part(qm_list,1)
       counter = 0

       !print QM list
       do j=1,116
         if (any(my_atoms%Z(list(1:size(list))).eq.j)) then
            call print('    &QM_KIND '//ElementName(j),file=input_file)
            do i=1,qm_list%N
               if (my_atoms%Z(list(i)).eq.j) then
                  call print('      MM_INDEX '//list(i),file=input_file)
                  counter = counter + 1
               endif
            enddo
            call print('    &END QM_KIND',file=input_file)
         endif
       enddo

       if (size(list).ne.counter) call system_abort('cp2k_write_input: Number of QM list atoms '//size(list)// &
                                                    ' differs from the written atom numbers '//counter)
     
       deallocate(list)
     
       call print('  &END QMMM',file=input_file)
    endif

!    call print(' ...lattice parameters...')
    call print('  &SUBSYS',file=input_file)
  
    call print('    &CELL',file=input_file)
    call print('      A '//round(my_atoms%lattice(1,1),6)//' '//round(my_atoms%lattice(2,1),6)//' '//round(my_atoms%lattice(3,1),6),file=input_file)
    call print('      B '//round(my_atoms%lattice(1,2),6)//' '//round(my_atoms%lattice(2,2),6)//' '//round(my_atoms%lattice(3,2),6),file=input_file)
    call print('      C '//round(my_atoms%lattice(1,3),6)//' '//round(my_atoms%lattice(2,3),6)//' '//round(my_atoms%lattice(3,3),6),file=input_file)
!    call print('      UNIT ANGSTROM',file=input_file)         !in Angstroms
    if (len(trim(param%global%cell_file)) > 0) then
      call initialise(cell_in_io,trim(param%global%cell_file), INPUT)
      stat = 0
      do while (stat == 0)
        line = read_line(cell_in_io, stat)
        if (stat == 0) call print(trim(line), file=input_file)
      end do
      call finalise(cell_in_io)
    else
      call print('      PERIODIC XYZ',file=input_file)   !in Angstroms
    endif
    call print('    &END CELL',file=input_file)

!    call print('...coordinates...')
  if (run_type.eq.QS_RUN) then
    call print('    &COORD',file=input_file)
    do i = 1, my_atoms%N
    call print('  '//ElementName(my_atoms%Z(i))//'  '//round(my_atoms%pos(1,i),10)//'  '//round(my_atoms%pos(2,i),10)//'  '//round(my_atoms%pos(3,i),10),file=input_file)
    enddo
    call print('    &END COORD',file=input_file)
    call print('#    &VELOCITY',file=input_file)
    call print('#    &END VELOCITY',file=input_file)
    call print('    &TOPOLOGY',file=input_file)
    call print('      CENTER_COORDINATES',file=input_file)
    call print('    &END TOPOLOGY',file=input_file)
  else
    call print('#    &COORD',file=input_file)
    call print('#    &END COORD',file=input_file)
    call print('#    &VELOCITY',file=input_file)
    call print('#    &END VELOCITY',file=input_file)
    call print('    &TOPOLOGY',file=input_file)
    call print('      &DUMP_PSF',file=input_file)
    call print('      &END DUMP_PSF',file=input_file)
    call print('      &DUMP_PDB',file=input_file)
    call print('      &END DUMP_PDB',file=input_file)
!    call print('      CENTER_COORDINATES',file=input_file)
    call print('      &GENERATE',file=input_file)
    call print('        BONDLENGTH_MAX     5.0',file=input_file)
    call print('        REORDER F',file=input_file)
    call print('        CREATE_MOLECULES F',file=input_file)

   ! print isolated atom list!
    call print('        &ISOLATED_ATOMS',file=input_file)
    do i = 1, my_atoms%N
!       if ((my_atoms%connect%neighbour1(i)%N + my_atoms%connect%neighbour2(i)%N).lt.1) then
!          call print('          LIST '//i,file=input_file)
!          call print('Written isolated atom '//my_atoms%Z(i)//' with atom number '//i//': has no neighbour')
!          cycle
!       endif
! needs different name for the QM molecules, say, only the atom name could be both the res_name and the mol_name!
       if (any(run_type.eq.(/QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
          if (any(i.eq.int_part(qm_list,1))) then
             call print('          LIST '//i,file=input_file)
!             call print('Written isolated atom '//my_atoms%Z(i)//' with atom number '//i//': QM atom')
             cycle
          endif
       endif
       if (any(trim(my_atoms%data%str(atom_res_name_index,i)).eq.(/'SOD','CLA'/))) then !'MCL'?
          call print('          LIST '//i,file=input_file)
!          call print('Written isolated atom '//my_atoms%Z(i)//' with atom number '//i//': methyl-chloride or 1 atom molecule: '//trim(my_atoms%data%str(atom_res_name_index,i)))
       cycle
       endif
    enddo
    call print('        &END ISOLATED_ATOMS',file=input_file)

    call print('      &END GENERATE',file=input_file)
    call print('      CHARGE_EXTENDED',file=input_file)
    call print('      COORD_FILE_NAME '//param%wenv%pdb_file,file=input_file)
    call print('      COORDINATE PDB',file=input_file)
    if (any(PSF_print.eq.(/DRIVER_PRINT_AND_SAVE,USE_EXISTING_PSF/))) then
       call print('      CONN_FILE_NAME ../'//trim(param%wenv%psf_file),file=input_file)
       call print('      CONN_FILE_FORMAT PSF',file=input_file)
    endif
    call print('    &END TOPOLOGY',file=input_file)
  endif
!    call print(' ...basis sets and potentials...')

    if (any(run_type.eq.(/QS_RUN,QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
       do j = 1,100
          if (any(my_atoms%Z(1:my_atoms%N).eq.j)) then
             call print('    &KIND '//ElementName(j),file=input_file)
             if (get_value(param%basis_set,ElementName(j),basis_set)) then
                call print('      BASIS_SET '//trim(basis_set),file=input_file)
             else
                call system_abort('There is no defined basis set for element '//ElementName(j))
             endif
             if (get_value(param%potential,ElementName(j),potential)) then
                call print('      POTENTIAL '//trim(potential),file=input_file)
             else
                call system_abort('There is no defined potential for element '//ElementName(j))
             endif
             call print('    &END KIND',file=input_file)
          endif
       enddo
    endif
    call print('  &END SUBSYS',file=input_file)
    call print('&END FORCE_EVAL',file=input_file)
     
!    call print(' ...any other internal variables...')

    call print('&GLOBAL',file=input_file)
    call print('  PROJECT '//trim(param%global%project),file=input_file)
    call print('  RUN_TYPE '//trim(param%global%runtype),file=input_file)
    call print('  PRINT_LEVEL LOW',file=input_file)
    call print('&END GLOBAL',file=input_file)
    call print('&MOTION',file=input_file)
    call print('  &PRINT',file=input_file)
    call print('    &FORCES',file=input_file)
    call print('      FORMAT XMOL',file=input_file)
    call print('    &END FORCES',file=input_file)
    call print('  &END PRINT',file=input_file)
    call print('  &MD',file=input_file)
    call print('    ENSEMBLE '//trim(param%md%ensemble),file=input_file)  !not needed,only for input check
    call print('    STEPS '//param%md%steps,file=input_file)      !only calculates forces and writes them
    call print('    TIMESTEP '//round(param%md%timestep,2),file=input_file)  !not needed, only for input check
    call print('    TEMPERATURE '//round(param%md%temperature,2),file=input_file)  !TEMPERATURE
    call print('    TEMP_KIND',file=input_file)
    call print('    &PRINT',file=input_file)
    call print('      &TEMP_KIND',file=input_file)
    call print('      &END TEMP_KIND',file=input_file)
    call print('    &END PRINT',file=input_file)
    call print('  &END MD',file=input_file)
    call print('&END MOTION',file=input_file)
  
    call finalise(input_file)

  end subroutine write_cp2k_input_file


  !!!!! ============================================================================================ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! \/ \/ \/ \/ \/ \/ ----- read forces and energy from CP2K output file ----- \/ \/ \/ \/ \/ \/ !!!!!
  !!!!!                                                                                              !!!!!
  !!!!! ============================================================================================ !!!!!

  subroutine read_cp2k_forces(forces,energy,param,my_atoms,run_type,verbose)

    real(dp),allocatable,dimension(:,:),intent(out) :: forces
    real(dp),                           intent(out) :: energy
    type(param_cp2k),                   intent(in)  :: param
    type(Atoms),                        intent(in)  :: my_atoms
    integer,                            intent(in)  :: run_type
    logical, optional,                  intent(in)  :: verbose
 
    integer :: i,j,num_atoms,m,status
    type(Inoutput) :: force_file
    logical :: my_verbose
    character(len=1024) :: commentline
    character(len=2048), dimension(10) :: fields
    integer :: num_fields
    type(Dictionary) :: params
    logical :: energy_found

    my_verbose = .false.
    if (present(verbose)) my_verbose = verbose

    call initialise(force_file,trim(param%wenv%working_directory)//'/'//trim(param%wenv%force_file),action=INPUT,isformatted=.true.)
    call print('Reading forces from file '//trim(force_file%filename)//'...')

    call parse_line(force_file,' ',fields,num_fields,status)
    if (status > 0) then
       call system_abort('read_cp2k_forces: Error reading from '//force_file%filename)
    else if (status < 0) then
       call system_abort('read_cp2k_forces: End of file when reading from '//force_file%filename)
    end if

    num_atoms = string_to_int(fields(1))
    if (num_atoms.ne.my_atoms%N) call system_abort('read_cp2k_forces: different number of atoms in the force file and in the my_atoms object')
    if (my_verbose) &
       call print('Reading force coordinates of '//num_atoms//' atoms')
    allocate(forces(3,num_atoms))

    !read the energy from the command line
    commentline = read_line(force_file)
    call read_string(params, commentline)
    energy_found = get_value(params, 'E',energy)
    if (.not.energy_found) call system_abort('read_cp2k_forces: no energy found in '//force_file%filename//', update CP2K you are using!!')

    !read forces
    do i=1,num_atoms
       call parse_line(force_file,' ',fields,num_fields)
       if (num_fields < 4) exit
       do j = 2, 4
          forces((j-1),i) = string_to_real(fields(j))
       end do
!       call print('read force for atom '//i//': '//round(forces(1,i),3)//' '//round(forces(2,i),3)//' '//round(forces(3,i),3))
    enddo
    call energy_conversion(energy)
    call force_conversion(forces)
    call finalise(force_file)

   ! re-reorder atoms/forces in case CP2K reordered them
    call read_convert_back_pos(forces,param,my_atoms,run_type=run_type,verbose=my_verbose)

    call print('')
    call print('The energy of the system: '//energy)
    call verbosity_push_decrement()
      call print('The forces acting on each atom (eV/A):')
      call print('atom     F(x)     F(y)     F(z)')
      do m=1,size(forces,2)
        call print('  '//m//'    '//forces(1,m)//'  '//forces(2,m)//'  '//forces(3,m))
      enddo
    call verbosity_pop()
    call print('Sum of the forces: '//sum(forces(1,1:num_atoms))//' '//sum(forces(2,1:num_atoms))//' '//sum(forces(3,1:num_atoms)))

!    if (any(abs(forces(1:3,1:my_atoms%N)).gt.(4.21648784E-02*HARTREE/BOHR))) &
!       call print('at least one too large force')

  end subroutine read_cp2k_forces

!recommended:
!correct forces according to new *-pos-* file, no need to hack cp2k!!
  subroutine read_convert_back_pos(forces,param,my_atoms,run_type,verbose)

    real(dp), allocatable, dimension(:,:), intent(inout) :: forces
    type(param_cp2k),                      intent(in)  :: param
    type(Atoms),                           intent(in)  :: my_atoms
    integer,                               intent(in)  :: run_type
    logical, optional,                     intent(in)  :: verbose

    type(InOutput)                        :: pos
    type(Atoms)                           :: reordered
    type(Table)                           :: qm_list
    real(dp), dimension(3)                :: shift
    integer, allocatable, dimension(:)    :: oldpos
    real(dp), allocatable, dimension(:,:) :: tforces
    integer                               :: i, j, status
    logical                               :: my_verbose

    if (.not.any(run_type.eq.(/QS_RUN,MM_RUN,QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) &
       call system_abort('read_convert_back_matrix: Run type is not recognized: '//run_type)

    my_verbose = optional_default(.false.,verbose)
    if (my_verbose) &
       call print_title('Reading matrix of new atomic coordinates from CP2K output')

   ! Getting the oldpos matrix from *-pos-* file...
    call initialise(pos,filename=trim(param%wenv%working_directory)//'/'//trim(param%wenv%pos_file))

    call print('Reading possibly reordered atoms from file '//pos%filename)
    call read_xyz(reordered,pos,status=status)
    if (status.ne.0) call system_abort('no atoms could be read from file '//pos%filename)
    if (reordered%N.ne.my_atoms%N) call system_abort('different number of atoms in the original and reordered atoms object')

    allocate(oldpos(my_atoms%N))
    oldpos = 0

    !match the atoms in the 2 atoms object to obtain the oldpos list

   ! shifted cell in case of QMMM (cp2k/src/topology_coordinate_util.F)
    if (any(run_type.eq.(/QMMM_RUN_CORE,QMMM_RUN_EXTENDED/))) then
       call get_qm_list_int_rec(my_atoms,run_type,qm_list,do_recursive=.true.)
       do i=1,3
         shift(i) = 0.5_dp * param%qmmm%qmmm_cell(i) &
                    - ( MINVAL(my_atoms%pos(i,int_part(qm_list,1))) + &
                        MAXVAL(my_atoms%pos(i,int_part(qm_list,1))) ) * 0.5_dp
       enddo
    else
       shift = 0._dp
    endif

    ! try to get reordering using default shift
    do i=1,my_atoms%N
       do j=1,reordered%N
!          if ((my_atoms%pos(1:3,i)+shift).feq.reordered%pos(1:3,j)) then
          if (matrix_feq2(my_atoms%pos(1:3,i)+shift,reordered%pos(1:3,j))) then
!             call print('atom '//i//' became atom '//j)
             oldpos(i) = j
             exit
          endif
       enddo
    enddo

    ! if failed, try reordering again with shift of a/2 b/2 c/2 (in case user set TOPOLOGY%CENTER_COORDINATES)
    if (any(oldpos.eq.0) .and. run_type .eq. QS_RUN) then
      shift = sum(my_atoms%lattice(:,:),2)/2.0_dp - &
              (minval(my_atoms%pos(:,:),2)+maxval(my_atoms%pos(:,:),2))/2.0_dp
      do i=1,my_atoms%N
         do j=1,reordered%N
  !          if ((my_atoms%pos(1:3,i)+shift).feq.reordered%pos(1:3,j)) then
            if (matrix_feq2(my_atoms%pos(1:3,i)+shift,reordered%pos(1:3,j))) then
  !             call print('atom '//i//' became atom '//j)
               oldpos(i) = j
               exit
            endif
        enddo
      enddo
    endif

    if (any(oldpos.eq.0)) then
      call system_abort('could not match the 2 atoms objects')
    endif

    call finalise(pos)
   ! End of getting the oldpos matrix from *-pos-* file...

    allocate(tforces(3,my_atoms%N))
    tforces = 0.
    do i=1,my_atoms%N
    tforces(1:3,oldpos(i)) = forces(1:3,i)
    enddo

    forces = tforces

    deallocate(tforces)
    deallocate(oldpos)
  end subroutine read_convert_back_pos

  subroutine QUIP_combine_forces(qmmm_forces,mm_forces,combined_forces,my_atoms)
  ! use momentum conservation DOES NOT ASSUME that SUM(F(MM)) = 0._dp:
  ! for  fitlist:   F_i = F(MM)    - m_i * SUMcore(F_i) / SUMfitlist(m_i)
  ! for QM  core:   F_i = F(QM/MM) - m_i * SUMcore(F_i) / SUMfitlist(m_i)
  ! for MM atoms:   F_i = F(MM)
  ! inner qm_list should be a reordered qm_list !!
  ! wrk1 should be written into file !!

    real(dp), allocatable, dimension(:,:), intent(in)    :: qmmm_forces
    real(dp), allocatable, dimension(:,:), intent(in)    :: mm_forces
    real(dp), allocatable, dimension(:,:), intent(out)   :: combined_forces
    type(Atoms),                           intent(in)    :: my_atoms

    integer               :: i_atom,m
    real(dp),dimension(3) :: force_corr
    real(dp)              :: mass_corr
    type(Table)           :: core, fitlist

    if (size(qmmm_forces,2).ne.my_atoms%N.or.size(mm_forces,2).ne.my_atoms%N) call system_abort('combine_forces: the size of QMMM and/or MM forces is different: different number of atoms')

    call get_qm_list_int(my_atoms,1,core)
    call get_qm_list_int_rec(my_atoms,2,fitlist, do_recursive=.true.)

    force_corr = 0._dp
    mass_corr  = 0._dp
    allocate(combined_forces(3,my_atoms%N))
    combined_forces = 0._dp

    force_corr = sum(mm_forces,2)    !this can be deleted if SUM(F(MM)) = 0._dp
    force_corr = force_corr + sum((qmmm_forces(1:3,int_part(core,1))-mm_forces(1:3,int_part(core,1))),2)
    mass_corr = sum(ElementMass(my_atoms%Z(int_part(fitlist,1))))

    force_corr = force_corr / mass_corr
    call print('Force mixing: the force correction is '//force_corr(1)//' '//force_corr(2)//' '//force_corr(3))

    combined_forces(1:3,1:my_atoms%N) = mm_forces(1:3,1:my_atoms%N)
    combined_forces(1:3,int_part(core,1)) = qmmm_forces(1:3,int_part(core,1))
    do i_atom=1,my_atoms%N
       if (any(int_part(fitlist,1).eq.i_atom)) combined_forces(1:3,i_atom) = combined_forces(1:3,i_atom) - ElementMass(my_atoms%Z(i_atom)) * force_corr 
    enddo

    call print('')
    call print('No energy after force mixing.')
    call verbosity_push_decrement()
      call print('The forces acting on each atom (eV/A):')
      call print('atom     F(x)     F(y)     F(z)')
      do m=1,size(combined_forces,2)
      call print('  '//m//'    '//combined_forces(1,m)//'  '//combined_forces(2,m)//'  '//combined_forces(3,m))
      enddo
    call verbosity_pop()
    call print('Sum of the combined forces: '//sum(combined_forces(1,1:my_atoms%N))//' '//sum(combined_forces(2,1:my_atoms%N))//' '//sum(combined_forces(3,1:my_atoms%N)))

    
  end subroutine QUIP_combine_forces

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

!  subroutine combine_forces(qmmm_forces,mm_forces,combined_forces,qmlist,my_atoms)
!  ! use momentum conservation:
!  ! for QM atoms:   F_corr,i = m_i * SUMall(F_i) / SUMqm(m_i)
!  ! for MM atoms:   F_corr,i = 0
!  ! inner qm_list should be a reordered qm_list !!
!  ! wrk1 should be written into file !!
!
!    real(dp), allocatable, dimension(:,:), intent(in)    :: qmmm_forces
!    real(dp), allocatable, dimension(:,:), intent(in)    :: mm_forces
!    real(dp), allocatable, dimension(:,:), intent(out)   :: combined_forces
!    type(Table),                           intent(in)    :: qmlist
!    type(Atoms),                           intent(in)    :: my_atoms
!
!    integer               :: i_atom
!    real(dp),dimension(3) :: force_corr
!    real(dp)              :: mass_corr
!
!    if (size(qmmm_forces,2).ne.my_atoms%N.or.size(mm_forces,2).ne.my_atoms%N) call system_abort('combine_forces: the size of QMMM and/or MM forces is different: different number of atoms')
!    if (any(int_part(qmlist,1).gt.my_atoms%N).or.any(int_part(qmlist,1).lt.1)) call system_abort('combine_forces: QM atom out of range 1-'//my_atoms%N)
!
!    force_corr = 0._dp
!    mass_corr  = 0._dp
!    allocate(combined_forces(3,my_atoms%N))
!    combined_forces = 0._dp
!
!    combined_forces(1:3,1:my_atoms%N) = mm_forces(1:3,1:my_atoms%N)
!    combined_forces(1:3,int_part(qmlist,1)) = qmmm_forces(1:3,int_part(qmlist,1))
!    force_corr = sum(combined_forces,2)
!    mass_corr = sum(ElementMass(my_atoms%Z(int_part(qmlist,1))))
!
!    force_corr = force_corr / mass_corr
!
!    do i_atom=1,my_atoms%N
!       if (any(int_part(qmlist,1).eq.i_atom)) combined_forces(1:3,i_atom) = qmmm_forces(1:3,i_atom) - ElementMass(my_atoms%Z(i_atom)) * force_corr 
!    enddo
!    
!  end subroutine combine_forces

  subroutine energy_conversion(energy)

    real(dp), intent(inout) :: energy

    energy = energy * HARTREE

  end subroutine energy_conversion

  subroutine force_conversion(force)

    real(dp), dimension(:,:), intent(inout) :: force
    integer                    :: i

    if (size(force,2).le.0) call system_abort('force_conversion: no forces!')
    if (size(force,1).ne.3) call system_abort('force_conversion: '//size(force,1)//' components of force instead of 3')

    do i=1,size(force,2)
       force(1:3,i) = force(1:3,i) * HARTREE / BOHR
    enddo

  end subroutine force_conversion

  subroutine velocity_conversion(at)

    type(Atoms), intent(inout) :: at
    integer                    :: i

    if (.not.associated(at%velo)) call system_abort('velocity_conversion: no velocities!')

    do i=1,at%N
       at%velo(1:3,i) = at%velo(1:3,i) * (HARTREE * BOHR / HBAR)
    enddo

  end subroutine velocity_conversion

  subroutine velocity_conversion_rev(at)

    type(Atoms), intent(inout) :: at
    integer                    :: i

    if (.not.associated(at%velo)) call system_abort('velocity_conversion: no velocities!')

    do i=1,at%N
       at%velo(1:3,i) = at%velo(1:3,i) / (HARTREE * BOHR / HBAR)
    enddo

  end subroutine velocity_conversion_rev

   ! a .feq. b
   !
   ! floating point logical comparison
   ! used to compare atoms objects
   function real_feq2(x,y) result(feq)

     real(dp), intent(in) :: x, y
     logical              :: feq

     if (abs(x-y) > EPSILON_ZERO) then
        feq = .false.
     else
        feq = .true.
     end if
     
   end function real_feq2

   ! m1 .feq. m2
   !
   ! matrix floating point comparison 
   ! used to compare atoms objects
   function matrix_feq2(matrix1,matrix2) result (feq)
     real(dp),intent(in), dimension(:) :: matrix1
     real(dp),intent(in), dimension(:) :: matrix2

     integer::j
     logical::feq
     
     call check_size('Matrix2',matrix2,shape(matrix1),'Matrix_FEQ')
     
     feq =.true.
     do j=1,size(matrix1)
           if (.not.real_feq2(matrix1(j),matrix2(j))) then
              feq=.false.
              return
           end if
     end do
     
   end function matrix_feq2

end module cp2k_driver_module
