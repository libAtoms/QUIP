module cp2k_driver_template_module
use libatoms_module
implicit none

private

public :: do_cp2k_calc


contains

  subroutine do_cp2k_calc(at, f, e, args_str)
    type(Atoms), intent(inout) :: at
    real(dp), intent(out) :: f(:,:), e
    character(len=*), intent(in) :: args_str

    type(Dictionary) :: cli
    character(len=FIELD_LENGTH) :: run_type, cp2k_template_file, psf_print, cp2k_program
    logical :: clean_up_files, save_output_files
    integer :: max_n_tries
    real(dp) :: max_force_warning
    real(dp) :: qm_vacuum
    logical :: try_reuse_wfn

    character(len=128) :: method

    type(Inoutput) :: template_io
    integer :: template_n_lines
    character(len=1024), allocatable :: cp2k_template_a(:)

    character(len=1024) :: run_dir

    integer :: run_type_i

    type(Table) :: qm_list
    integer, allocatable :: qm_list_a(:)
    integer :: counter

    integer :: charge
    logical :: do_lsd

    logical :: can_reuse_wfn, qm_list_changed

    logical :: use_QM, use_MM, use_QMMM

    integer, pointer :: isolated_atom(:)
    integer i, atno, insert_pos
    real(dp) :: cur_qmmm_qm_abc(3), old_qmmm_qm_abc(3)

    type(Atoms) :: at_cp2k

    character(len=TABLE_STRING_LENGTH), pointer :: from_str(:), to_str(:)
    character(len=TABLE_STRING_LENGTH) :: dummy_s
    real(dp), pointer :: from_dp(:), to_dp(:)
    integer, pointer :: from_i(:), to_i(:)

    integer, pointer :: old_cluster_mark_p(:), cluster_mark_p(:)
    logical :: dummy, have_silica_potential
    type(Table) :: intrares_impropers

    call system_timer('do_cp2k_calc')

    call initialise(cli)
      run_type = ""
      call param_register(cli, "Run_Type", PARAM_MANDATORY, run_type)
      cp2k_template_file = ""
      call param_register(cli, "cp2k_template_file", "cp2k_input.template", cp2k_template_file)
      psf_print = ""
      call param_register(cli, "PSF_print", "NO_PSF", psf_print)
      cp2k_program = ""
      call param_register(cli, "cp2k_program", PARAM_MANDATORY, cp2k_program)
      call param_register(cli, "clean_up_files", "T", clean_up_files)
      call param_register(cli, "save_output_files", "T", save_output_files)
      call param_register(cli, "max_n_tries", "2", max_n_tries)
      call param_register(cli, "max_force_warning", "2.0", max_force_warning)
      call param_register(cli, "qm_vacuum", "6.0", qm_vacuum)
      call param_register(cli, "try_reuse_wfn", "T", try_reuse_wfn)
      call param_register(cli, 'have_silica_potential', 'F', have_silica_potential) !if yes, use 2.8A SILICA_CUTOFF for the connectivities
      if (.not.param_read_line(cli, args_str, do_check=.true.,ignore_unknown=.false.,task='cp2k_filepot_template args_str')) &
	call system_abort('could not parse argument line')
    call finalise(cli)

    call print("do_cp2k_calc command line arguments")
    call print("  Run_Type " // Run_Type)
    call print("  cp2k_template_file " // cp2k_template_file)
    call print("  PSF_print " // PSF_print)
    call print("  clean_up_files " // clean_up_files)
    call print("  save_output_files " // save_output_files)
    call print("  max_n_tries " // max_n_tries)
    call print("  max_force_warning " // max_force_warning)
    call print("  qm_vacuum " // qm_vacuum)
    call print("  try_reuse_wfn " // try_reuse_wfn)
    call print('  have_silica_potential '//have_silica_potential)

    ! read template file
    call initialise(template_io, trim(cp2k_template_file), INPUT)
    call read_file(template_io, cp2k_template_a, template_n_lines)
    call finalise(template_io)

    call prefix_cp2k_input_sections(cp2k_template_a(1:template_n_lines))

    if ( (trim(psf_print) /= 'NO_PSF') .and. &
         (trim(psf_print) /= 'CP2K_PRINT_AND_SAVE') .and. &
         (trim(psf_print) /= 'DRIVER_PRINT_AND_SAVE') .and. &
         (trim(psf_print) /= 'USE_EXISTING_PSF')) &
      call system_abort("Unknown value for psf_print '"//trim(psf_print)//"'")

    ! parse run_type
    use_QM = .false.
    use_MM = .false.
    use_QMMM = .false.
    select case(trim(run_type))
      case("QS")
	use_QM=.true.
	method="QS"
	run_type_i = QS_RUN
      case("MM")
	use_MM=.true.
	method="Fist"
	run_type_i = MM_RUN
      case("QMMM_CORE")
	use_QM=.true.
	use_MM=.true.
	use_QMMM=.true.
	method="QMMM"
	run_type_i = QMMM_RUN_CORE
      case("QMMM_EXTENDED")
	use_QM=.true.
	use_MM=.true.
	use_QMMM=.true.
	method="QMMM"
	run_type_i = QMMM_RUN_EXTENDED
      case default
	call system_abort("Unknown run_type "//trim(run_type))
    end select

    ! prepare CHARMM params if necessary
    if (use_MM) then
      if (have_silica_potential) then
        call set_cutoff(at,SILICA_2body_CUTOFF)
      else
        call set_cutoff(at,0._dp)
      endif
      call calc_connect(at)

      call map_into_cell(at)
      call calc_dists(at)

      call create_CHARMM(at,do_CHARMM=.true.,intrares_impropers=intrares_impropers)
    endif

    ! put in method
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "", "&FORCE_EVAL")
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL METHOD "//trim(method), after_line = insert_pos, n_l = template_n_lines)

    ! get qm_list
    if (use_QMMM) then
      call get_qm_list(at, run_type_i, qm_list, do_union=.true.)
      allocate(qm_list_a(qm_list%N))
      qm_list_a = int_part(qm_list,1)
    endif

    if (qm_list%N == at%N) then
      call print("WARNING: requested '"//trim(run_type)//"' but all atoms are in QM region, doing full QM run instead", ERROR)
      run_type='QS'
      use_QM = .true.
      use_MM = .false.
      use_QMMM = .false.
      method = 'QS'
    endif

    can_reuse_wfn = .true.

    ! put in things needed for QMMM
    if (use_QMMM) then

      if (trim(run_type) == "QMMM_CORE") then
	run_type_i = QMMM_RUN_CORE
      else if (trim(run_type) == "QMMM_EXTENDED") then
	run_type_i = QMMM_RUN_EXTENDED
      else
	call system_abort("Unknown run_type '"//trim(run_type)//"' with use_QMMM true")
      endif

      insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&QMMM", "&CELL")
      call print('INFO: The size of the QM cell is either the MM cell itself, or it will have at least '//(qm_vacuum/2.0_dp)// &
			' Angstrom around the QM atoms.')
      call print('WARNING! Please check if your cell is centered around the QM region!',ERROR)
      call print('WARNING! CP2K centering algorithm fails if QM atoms are not all in the',ERROR)
      call print('WARNING! 0,0,0 cell. If you have checked it, please ignore this message.',ERROR)
      cur_qmmm_qm_abc = qmmm_qm_abc(at, run_type_i, qm_vacuum)
      call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&QMMM&CELL ABC " // cur_qmmm_qm_abc, after_line=insert_pos, n_l=template_n_lines); insert_pos = insert_pos + 1
      call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&QMMM&CELL PERIODIC XYZ", after_line=insert_pos, n_l=template_n_lines); insert_pos = insert_pos + 1

      if (get_value(at%params, "QM_cell", old_qmmm_qm_abc)) then
	if (cur_qmmm_qm_abc .fne. old_qmmm_qm_abc) can_reuse_wfn = .false.
      endif

      !check if QM list changed: compare cluster_mark and old_cluster_mark
!      if (get_value(at%params, "QM_list_changed", qm_list_changed)) then
!       if (qm_list_changed) can_reuse_wfn = .false.
!      endif
       if (.not.has_property(at, 'cluster_mark')) call system_abort('no cluster_mark found in atoms object')
       if (.not.has_property(at, 'old_cluster_mark')) call system_abort('no old_cluster_mark found in atoms object')
       dummy = assign_pointer(at, 'old_cluster_mark', old_cluster_mark_p)
       dummy = assign_pointer(at, 'cluster_mark', cluster_mark_p)

       qm_list_changed = .false.
       do i=1,at%N
          !only hybrid_no_mark matters
          if (old_cluster_mark_p(i).ne.cluster_mark_p(i) .and. &
              any((/old_cluster_mark_p(i),cluster_mark_p(i)/).eq.HYBRID_NO_MARK)) then
              qm_list_changed = .true.
          endif
       enddo
       call set_value(at%params,'QM_list_changed',qm_list_changed)
       call print('set_value QM_list_changed '//qm_list_changed)

       if (qm_list_changed) can_reuse_wfn = .false.

      counter = 0
      do atno=minval(at%Z), maxval(at%Z)
	if (any(at%Z(qm_list_a) == atno)) then
	  insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&QMMM", "&QM_KIND "//ElementName(atno))
	  do i=1, qm_list%N
	    if (at%Z(qm_list_a(i)) == atno) then
	      call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&QMMM&QM_KIND-"//trim(ElementName(atno))// &
							" MM_INDEX "//qm_list_a(i), after_line = insert_pos, n_l = template_n_lines)
	      insert_pos = insert_pos + 1
	      counter = counter + 1
	    endif
	  end do
	end if
      end do
      if (size(qm_list_a) /= counter) &
	call system_abort("Number of QM list atoms " // size(qm_list_a) // " doesn't match number of QM_KIND atoms " // counter)
    endif

    ! put in things needed for QM
    if (use_QM) then
      if (try_reuse_wfn .and. can_reuse_wfn) then 
	insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL", "&DFT")
	call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&DFT WFN_RESTART_FILE_NAME ../wfn.restart.wfn", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
	insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&DFT", "&SCF")
	call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&DFT&SCF SCF_GUESS RESTART", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
      endif
      call calc_charge_lsd(at, qm_list, charge, do_lsd)
      insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL", "&DFT")
      call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&DFT CHARGE "//charge, after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
      if (do_lsd) call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&DFT LSD ", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    endif

    ! put in unit cell
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL", "&SUBSYS")

    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS", "&CELL")
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&CELL A " // at%lattice(:,1), after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&CELL B " // at%lattice(:,2), after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&CELL C " // at%lattice(:,3), after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1

    ! put in topology
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS", "&TOPOLOGY")
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS&TOPOLOGY", "&DUMP_PSF")
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS&TOPOLOGY", "&DUMP_PDB")
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS&TOPOLOGY", "&GENERATE")
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE REORDER F", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE CREATE_MOLECULES F", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE", "&ISOLATED_ATOMS")
    if (use_QMMM) then
      do i=1, size(qm_list_a)
	call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE&ISOLATED_ATOMS LIST " // qm_list_a(i), after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
      end do
    endif
    if (assign_pointer(at, "isolated_atom", isolated_atom)) then
      do i=1, at%N
	if (isolated_atom(i) /= 0) then
	  call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE&ISOLATED_ATOMS LIST " // i, after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
	endif
      end do
    endif

    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS", "&TOPOLOGY")
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY COORD_FILE_NAME quip_cp2k.xyz", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY COORDINATE EXYZ", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    if (trim(psf_print) == "DRIVER_PRINT_AND_SAVE" .or. trim(psf_print) == "USE_EXISTING_PSF") then
      call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY CONN_FILE_NAME ../quip_cp2k.psf", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
      call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY CONN_FILE_FORMAT PSF", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    endif

    ! put in global stuff to run a single force evalution, print out appropriate things
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "", "&GLOBAL")
    call insert_cp2k_input_line(cp2k_template_a, "&GLOBAL   PROJECT quip", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    call insert_cp2k_input_line(cp2k_template_a, "&GLOBAL   RUN_TYPE MD", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1

    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "", "&MOTION")
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&MOTION", "&PRINT")
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&MOTION&PRINT", "&FORCES")
    call insert_cp2k_input_line(cp2k_template_a, "&MOTION&PRINT&FORCES       FORMAT XMOL", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&MOTION", "&MD")
    call insert_cp2k_input_line(cp2k_template_a, "&MOTION&MD     ENSEMBLE NVE", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    call insert_cp2k_input_line(cp2k_template_a, "&MOTION&MD     STEPS 0", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1

    ! write psf file if necessary

    if (run_type /= "QS") then
      if (trim(psf_print) == "DRIVER_PRINT_AND_SAVE") then
	call write_psf_file(at, "quip_cp2k.psf", run_type_string=trim(run_type),intrares_impropers=intrares_impropers,add_silica_23body=have_silica_potential)
      endif
    endif

    run_dir = make_run_directory()

    call write_cp2k_input_file(cp2k_template_a(1:template_n_lines), trim(run_dir)//'/cp2k_input.inp')

    ! prepare xyz file for input to cp2k

!    call atoms_copy_without_connect(at_cp2k, at, properties="species:pos")
!    dummy_s = ""
!    ! hopefully no need for this once topology.f95 and cp2k patch agree on property names
!!    if (run_type /= "QS") then
!      call add_property(at_cp2k, "atmname", dummy_s, n_cols=1)
!      call add_property(at_cp2k, "molname", dummy_s, n_cols=1)
!      call add_property(at_cp2k, "resname", dummy_s, n_cols=1)
!      call add_property(at_cp2k, "resid", 0, n_cols=1)
!      call add_property(at_cp2k, "atm_charge", 0.0_dp, n_cols=1)
!      if (.not. assign_pointer(at_cp2k, "atmname", to_str)) &
!        call system_abort("impossible failure to set pointer to at_cp2k%atmname")
!      if (.not. assign_pointer(at, "atom_type", from_str)) &
!        call system_abort("failed to set pointer to at%atom_type")
!      to_str = from_str
!      if (.not. assign_pointer(at_cp2k, "molname", to_str)) &
!        call system_abort("impossible failure to set pointer to at_cp2k%molname")
!      if (.not. assign_pointer(at, "atom_mol_name", from_str)) &
!        call system_abort("failed to set pointer to at%atom_mol_name")
!      to_str = from_str
!      if (.not. assign_pointer(at_cp2k, "resname", to_str)) &
!        call system_abort("impossible failure to set pointer to at_cp2k%resname")
!      if (.not. assign_pointer(at, "atom_res_name", from_str)) &
!        call system_abort("failed to set pointer to at%atom_res_name")
!      to_str = from_str
!      if (.not. assign_pointer(at_cp2k, "resid", to_i)) &
!        call system_abort("impossible failure to set pointer to at_cp2k%resid")
!      if (.not. assign_pointer(at, "atom_res_number", from_i)) &
!        call system_abort("failed to set pointer to at%atom_res_number")
!      to_i = from_i
!      if (.not. assign_pointer(at_cp2k, "atm_charge", to_dp)) &
!        call system_abort("impossible failure to set pointer to at_cp2k%atm_charge")
!      if (.not. assign_pointer(at, "atom_charge", from_dp)) &
!        call system_abort("failed to set pointer to at%atom_charge")
!      to_dp = from_dp
!!    endif
!    call print_xyz(at_cp2k, trim(run_dir)//'/quip_cp2k.xyz', all_properties=.true.)
    call print_xyz(at, trim(run_dir)//'/quip_cp2k.xyz', all_properties=.true.)

    ! actually run cp2k

    call run_cp2k_program(trim(cp2k_program), trim(run_dir), max_n_tries)

    ! parse output
    call read_energy_forces(at, qm_list, cur_qmmm_qm_abc, trim(run_dir), "quip", e, f)

    if (maxval(abs(f)) > max_force_warning) &
      call print('WARNING cp2k forces max component ' // maxval(abs(f)) // ' at ' // maxloc(abs(f)) // &
		 ' exceeds warning threshold ' // max_force_warning, ERROR)

    ! save output

    if (trim(psf_print) == "CP2K_PRINT_AND_SAVE") &
      call system_command('cp '//trim(run_dir)//'/quip-dump-1.psf quip_cp2k.psf')

    if (use_QM) &
      call system_command('cp '//trim(run_dir)//'/quip-RESTART.wfn wfn.restart.wfn')

    if (save_output_files) then
      call system_command('cat '//trim(run_dir)//'/cp2k_input.inp'// &
        ' >> cp2k_input_log; echo "##############" >> cp2k_input_log;' // &
        ' cat '//trim(run_dir)//'/quip-frc-1.xyz'// &
        ' >> cp2k_force_file_log; echo "##############" >> cp2k_force_file_log;' // &
        ' cat '//trim(run_dir)//'/cp2k_output.out >> cp2k_output_log; echo "##############" >> cp2k_output_log')
    endif

    ! clean up

    if (clean_up_files) call system_command('rm -rf '//trim(run_dir))

    call system_timer('do_cp2k_calc')

  end subroutine do_cp2k_calc

  function find_make_cp2k_input_section(l_a, n_l, base_sec, new_sec) result(line_n)
    character(len=*), allocatable, intent(inout) :: l_a(:)
    integer, intent(inout) :: n_l
    character(len=*), intent(in) :: base_sec, new_sec
    integer :: line_n

    integer :: i, pamp, pspc
    character(len=1024) :: sec, word, arg, base_sec_root, base_sec_tail, new_sec_end

    line_n = 0

    do i=1, n_l
      call split_cp2k_input_line(trim(l_a(i)), sec, word, arg)

      if (trim(sec) == trim(base_sec) .and. trim(word) == trim(new_sec)) then
	line_n = i
	return
      endif
    end do

    if (len_trim(base_sec) == 0) then
      i = n_l
      call insert_cp2k_input_line(l_a, " "//trim(new_sec), after_line=i, n_l=n_l); i = i + 1
      pspc = index(trim(new_sec)," ")
      if (pspc == 0) then
	new_sec_end = "&END " // new_sec(2:len_trim(new_sec))
      else
	new_sec_end = "&END " // new_sec(2:pspc-1)
      endif
      call insert_cp2k_input_line(l_a, " " // trim(new_sec_end), after_line=i, n_l=n_l); i = i + 1
      line_n = n_l-1
      return
    endif

    pamp = index(base_sec,"&", back=.true.)
    if (pamp <= 1) then
      base_sec_root = ""
      base_sec_tail = trim(base_sec)
    else
      base_sec_root = base_sec(1:pamp-1)
      base_sec_tail = base_sec(pamp:len_trim(base_sec))
    endif


    do i=1, n_l
      call split_cp2k_input_line(trim(l_a(i)), sec, word, arg)
      if (trim(sec) == trim(base_sec_root) .and. trim(word) == trim(base_sec_tail)) then
	call insert_cp2k_input_line(l_a, trim(base_sec)//" "//trim(new_sec), after_line=i, n_l=n_l)
	pspc = index(trim(new_sec)," ")
	if (pspc == 0) then
	  new_sec_end = "&END " // new_sec(2:len_trim(new_sec))
	else
	  new_sec_end = "&END " // new_sec(2:pspc-1)
	endif
	call insert_cp2k_input_line(l_a, trim(base_sec)//" " // trim(new_sec_end), after_line=i+1, n_l=n_l)
	line_n = i+1
	return
      endif
    end do

    if (line_n == 0) &
      call system_abort("Could not find or make section '"//trim(new_sec)//" in base section '"//trim(base_sec))

  end function find_make_cp2k_input_section

  subroutine read_energy_forces(at, qm_list, cur_qmmm_qm_abc, run_dir, proj, e, f)
    type(Atoms), intent(in) :: at
    type(Table), intent(in) :: qm_list
    real(dp), intent(in) :: cur_qmmm_qm_abc(3)
    character(len=*), intent(in) :: run_dir, proj
    real(dp), intent(out) :: e, f(:,:)

    type(Atoms) :: f_xyz, p_xyz
    integer :: m

    call read_xyz(f_xyz, trim(run_dir)//'/'//trim(proj)//'-frc-1.xyz')
    call read_xyz(p_xyz, trim(run_dir)//'/'//trim(proj)//'-pos-1.xyz')

    if (.not. get_value(f_xyz%params, "E", e)) &
      call system_abort('read_energy_forces failed to find E value in '//trim(run_dir)//'/quip-frc-1.xyz file')
    f = f_xyz%pos

    e = e * HARTREE
    f  = f * HARTREE/BOHR 
    call reorder_if_necessary(at, qm_list, cur_qmmm_qm_abc, p_xyz%pos, f)

    call print('')
    call print('The energy of the system: '//e)
    call verbosity_push_decrement()
      call print('The forces acting on each atom (eV/A):')
      call print('atom     F(x)     F(y)     F(z)')
      do m=1,size(f,2)
        call print('  '//m//'    '//f(1,m)//'  '//f(2,m)//'  '//f(3,m))
      enddo
    call verbosity_pop()
    call print('Sum of the forces: '//sum(f,2))

  end subroutine read_energy_forces

  subroutine reorder_if_necessary(at, qm_list, qmmm_qm_abc, new_p, new_f)
    type(Atoms), intent(in) :: at
    type(Table), intent(in) :: qm_list
    real(dp), intent(in) :: qmmm_qm_abc(3)
    real(dp), intent(in) :: new_p(:,:)
    real(dp), intent(inout) :: new_f(:,:)

    real(dp) :: shift(3)
    integer, allocatable :: oldpos(:)
    integer :: i, j

    ! shifted cell in case of QMMM (cp2k/src/toplogy_coordinate_util.F)
    if (qm_list%N > 0) then
      do i=1,3
	shift(i) = 0.5_dp * qmmm_qm_abc(i) - (minval(at%pos(i,int_part(qm_list,1))) + maxval(at%pos(i,int_part(qm_list,1))))*0.5_dp
      end do
    else
      shift = 0.0_dp
    endif

    allocate(oldpos(at%N))

    do i=1, at%N
      do j=1, size(new_p,2)
	if (norm(at%pos(:,i)+shift-new_p(:,j)) <= 1.0e-4_dp) then
	  oldpos(i) = j
	  cycle
	endif
      end do
    end do

    if (any(oldpos == 0)) then
      ! try again with shift of a/2 b/2 c/2 in case TOPOLOGY%CENTER_COORDINATES is set
      shift = sum(at%lattice(:,:),2)/2.0_dp - &
	      (minval(at%pos(:,:),2)+maxval(at%pos(:,:),2))/2.0_dp
      do i=1, at%N
	do j=1, size(new_p,2)
	  if (norm(at%pos(:,i) + shift - new_p(:,j)) <= 1.0e-4_dp) oldpos(i) = j
	end do
      end do
    endif

    if (any(oldpos == 0)) &
      call system_abort("Could not match orig and read in atom objects")

    new_f(1,oldpos(:)) = new_f(1,:)
    new_f(2,oldpos(:)) = new_f(2,:)
    new_f(3,oldpos(:)) = new_f(3,:)

  end subroutine reorder_if_necessary

  subroutine run_cp2k_program(cp2k_program, run_dir, max_n_tries)
    character(len=*), intent(in) :: cp2k_program, run_dir
    integer, intent(in) :: max_n_tries

    integer :: n_tries
    logical :: converged
    character(len=1024) :: cp2k_run_command
    integer :: stat, error_stat

    n_tries = 0
    converged = .false.

    do while (.not. converged .and. (n_tries < max_n_tries))
      n_tries = n_tries + 1

      cp2k_run_command = 'cd ' // trim(run_dir)//'; '//trim(cp2k_program)//' cp2k_input.inp >> cp2k_output.out'
      call print("Doing '"//trim(cp2k_run_command)//"'")
      call system_timer('cp2k_run_command')
      call system_command(trim(cp2k_run_command), status=stat)
      call system_timer('cp2k_run_command')
      call print('grep -i warning '//trim(run_dir)//'/cp2k_output.out', ERROR)
      call system_command("fgrep -i 'warning' "//trim(run_dir)//"/cp2k_output.out")
      call system_command("fgrep -i 'error' "//trim(run_dir)//"/cp2k_output.out", status=error_stat)
      if (stat /= 0) &
	call system_abort('cp2k_run_command has non zero return status ' // stat //'. check output file '//trim(run_dir)//'/cp2k_output.out')
      if (error_stat == 0) &
	call system_abort('cp2k_run_command generated ERROR message in output file '//trim(run_dir)//'/cp2k_output.out')

      call system_command('egrep "FORCE_EVAL.* QS " '//trim(run_dir)//'/cp2k_output.out',status=stat)
      if (stat == 0) then ! QS or QMMM run
	call system_command('grep "FAILED to converge" '//trim(run_dir)//'/cp2k_output.out',status=stat)
	if (stat == 0) then
	  call print("WARNING: cp2k_driver failed to converge, trying again",ERROR)
	  converged = .false.
	else
	  call system_command('grep "SCF run converged" '//trim(run_dir)//'/cp2k_output.out',status=stat)
	  if (stat == 0) then
	    converged = .true.
	  else
	    call print("WARNING: cp2k_driver couldn't find definitive sign of convergence or failure to converge in output file, trying again",ERROR)
	    converged = .false.
	  endif
	end if
      else ! MM run
	converged = .true.
      endif
    end do

    if (.not. converged) &
      call system_abort('cp2k failed to converge after n_tries='//n_tries//'. see output file '//trim(run_dir)//'/cp2k_output.out')

  end subroutine run_cp2k_program

  function make_run_directory() result(dir)
    integer i
    character(len=1024) :: dir

    logical :: exists
    integer stat

    exists = .true.
    i = 0
    do while (exists)
      i = i + 1
      dir = "cp2k_run_"//i
      inquire(file=trim(dir)//'/cp2k_input.inp', exist=exists)
    end do
    call system_command("mkdir "//trim(dir), status=stat)

    if (stat /= 0) &
      call system_abort("Failed to mkdir "//trim(dir)//" status " // stat)

  end function make_run_directory

  subroutine write_cp2k_input_file(l_a, filename)
    character(len=*), intent(in) :: l_a(:)
    character(len=*), intent(in) :: filename

    integer :: i, pspc
    type(Inoutput) :: io

    call initialise(io, trim(filename), OUTPUT)
    do i=1, size(l_a)
      pspc = index(l_a(i), " ")
      call print(l_a(i)(pspc+1:len_trim(l_a(i))), file=io)
    end do
    call finalise(io)

  end subroutine write_cp2k_input_file

  subroutine split_cp2k_input_line(l, sec, word, arg)
    character(len=*) :: l, sec, word, arg

    character(len=len(l)) :: t_l
    integer :: pspc

    t_l = l

    pspc = index(trim(t_l), " ")
    if (pspc == 0 .or. pspc == 1) then ! top level
      sec = ""
    else
      sec = t_l(1:pspc-1)
    endif

    t_l = adjustl(t_l(pspc+1:len_trim(t_l)))

    pspc = index(trim(t_l), " ")
    if (pspc == 0) then ! no arg
      word = trim(t_l)
      arg = ""
    else
      word = t_l(1:pspc-1)
      arg = t_l(pspc+1:len_trim(t_l))
    endif

    sec = adjustl(sec)
    word = adjustl(word)
    arg = adjustl(arg)

  end subroutine split_cp2k_input_line

  subroutine insert_cp2k_input_line(l_a, l, after_line, n_l)
    character(len=*), allocatable, intent(inout) :: l_a(:)
    character(len=*), intent(in) :: l
    integer, intent(in) :: after_line
    integer, intent(inout) :: n_l

    integer :: i

    if (n_l+1 > size(l_a)) call extend_char_array(l_a)

    do i=n_l, after_line+1, -1
      l_a(i+1) = l_a(i)
    end do
    l_a(after_line+1) = trim(l)
    n_l = n_l + 1
  end subroutine insert_cp2k_input_line

  subroutine prefix_cp2k_input_sections(l_a)
    character(len=*), intent(inout) :: l_a(:)

    integer :: pamp
    character(len=1024) :: section_str, new_section_str
    integer :: i, j

    section_str = ""
    do i=1, size(l_a)
      if (index(l_a(i),"&END") /= 0) then
	pamp = index(section_str, "&", .true.)
	new_section_str = section_str(1:pamp-1)
	section_str = new_section_str
      else if (index(l_a(i),"&") /= 0) then
	pamp = index(l_a(i), "&")
	new_section_str = trim(section_str)
	do j=pamp, len_trim(l_a(i))
	  if (l_a(i)(j:j) == " ") then
	    new_section_str = trim(new_section_str) // "-"
	  else
	    new_section_str = trim(new_section_str) // l_a(i)(j:j)
	  endif
	end do
      endif
      l_a(i) = (trim(section_str) //" "//trim(l_a(i)))
      section_str = new_section_str
    end do
  end subroutine

  subroutine get_qm_list(at,qmflag,qm_list,do_union,int_property)
    type(Atoms), intent(in)  :: at
    integer,     intent(in)  :: qmflag
    type(Table), intent(out) :: qm_list
    logical, intent(in), optional :: do_union
    character(len=*), optional, intent(in) :: int_property
  
    integer :: i
    integer, pointer :: QM_flag(:)
    logical              :: my_do_union
    integer :: qmflag_min
    character(STRING_LENGTH) :: my_int_property

    my_do_union = optional_default(.false.,do_union)

    my_int_property = ''
    if (present(int_property)) then
       my_int_property = trim(int_property)
    else
       my_int_property = "cluster_mark"
    endif
    if (.not.(assign_pointer(at, trim(my_int_property), QM_flag))) &
      call system_abort("get_qm_list couldn't find "//trim(my_int_property)//" field")

    if (my_do_union) then
      qmflag_min = 1
    else
      qmflag_min = qmflag
    endif

    call initialise(qm_list,4,0,0,0,0)      !1 int, 0 reals, 0 str, 0 log, num_qm_atoms entries
    do i=1,at%N
      if (QM_flag(i) >= qmflag_min .and.  QM_flag(i) <= qmflag) call append(qm_list,(/i,0,0,0/))
    enddo

    if (qm_list%N.eq.0) call print('Empty QM list with QM_flag '//qmflag//' and my_do_union ' // my_do_union ,ERROR)
  end subroutine get_qm_list

  function qmmm_qm_abc(at, run_type_i, qm_vacuum)
    type(Atoms), intent(in) :: at
    integer, intent(in) :: run_type_i
    real(dp), intent(in) :: qm_vacuum
    real(dp) :: qmmm_qm_abc(3)

    real(dp) :: qm_maxdist(3)
    integer i, j
    type(Table) :: qm_list

    call get_qm_list(at,run_type_i,qm_list, do_union=.true.)

    qm_maxdist = 0.0_dp
    do i=1, qm_list%N
    do j=1, qm_list%N
      qm_maxdist(1) = max(qm_maxdist(1), at%pos(1,i)-at%pos(1,j))
      qm_maxdist(2) = max(qm_maxdist(2), at%pos(2,i)-at%pos(2,j))
      qm_maxdist(3) = max(qm_maxdist(3), at%pos(3,i)-at%pos(3,j))
    end do
    end do

    qmmm_qm_abc(1) = min(real(ceiling(qm_maxdist(1)))+qm_vacuum,at%lattice(1,1))
    qmmm_qm_abc(2) = min(real(ceiling(qm_maxdist(2)))+qm_vacuum,at%lattice(2,2))
    qmmm_qm_abc(3) = min(real(ceiling(qm_maxdist(3)))+qm_vacuum,at%lattice(3,3))

  end function

  subroutine calc_charge_lsd(at, qm_list, charge, do_lsd)
    type(Atoms), intent(in) :: at
    type(Table), intent(in) :: qm_list
    integer, intent(out) :: charge
    logical, intent(out) :: do_lsd

    real(dp), pointer :: atom_charge(:)

    if (qm_list%N > 0) then
      if (.not. assign_pointer(at, "atom_charge", atom_charge)) &
	call system_abort("calc_charge_lsd could not find atom_charge")
      charge = nint(sum(atom_charge(int_part(qm_list,1))))
      do_lsd = (mod(charge,2) /= 0)
    else
      if (get_value(at%params, 'Charge', charge)) call print("Using Charge " // charge)
      do_lsd = .false.
      if (get_value(at%params, 'LSD', do_lsd)) call print("Using do_lsd " // do_lsd)
    endif

  end subroutine


end module cp2k_driver_template_module
