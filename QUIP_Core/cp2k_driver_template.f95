module cp2k_driver_template_module
use libatoms_module
implicit none

private

public :: go_cp2k_template


contains

  subroutine go_cp2k_template(at, f, e, args_str)
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
    integer :: template_n_lines, additions_n_lines, combined_n_lines
    character(len=1024), allocatable :: cp2k_template_a(:)
    character(len=1024), allocatable :: cp2k_template_additions(:)
    character(len=1024), allocatable :: cp2k_template_combined(:)

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
    integer i, atno
    real(dp) :: cur_qmmm_qm_abc(3), old_qmmm_qm_abc(3)

    call system_timer('go_cp2k_template')

    call initialise(cli)
      run_type = ""
      call param_register(cli, "Run_Type", PARAM_MANDATORY, run_type)
      cp2k_template_file = ""
      call param_register(cli, "cp2k_template_file", "cp2k_template", cp2k_template_file)
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
      if (.not.param_read_line(cli, args_str, do_check=.true.,ignore_unknown=.false.,task='cp2k_filepot_template args_str')) &
	call system_abort('could not parse argument line')
    call finalise(cli)

    ! read template file
    call initialise(template_io, trim(cp2k_template_file), INPUT)
    call read_file(template_io, cp2k_template_a, template_n_lines)
    call finalise(template_io)

    call prefix_sections(cp2k_template_a(1:template_n_lines))

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
      call set_cutoff(at,0._dp)
      call calc_connect(at)

      call map_into_cell(at)
      call calc_dists(at)

      call create_CHARMM(at,do_CHARMM=.true.)
    endif

    additions_n_lines = 0

    ! put in method
    call append_line(cp2k_template_additions, "&FORCE_EVAL METHOD "//trim(method), additions_n_lines)

    ! get qm_list
    if (use_QMMM) then
      call get_qm_list(at, run_type_i, qm_list, do_union=.true.)
      allocate(qm_list_a(qm_list%N))
      qm_list_a = int_part(qm_list,1)
    endif

    if (qm_list%N == at%N) then
      call print("WARNING: requested '"//trim(run_type)//"' but all atoms are in QM region, doing full QM run instead")
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

      call append_line(cp2k_template_additions, "&FORCE_EVAL&QMMM &CELL", additions_n_lines)
      call print('INFO: The size of the QM cell is either the MM cell itself, or it will have at least '//(qm_vacuum/2.0_dp)// &
			' Angstrom around the QM atoms.')
      call print('WARNING! Please check if your cell is centered around the QM region!',ERROR)
      call print('WARNING! CP2K centering algorithm fails if QM atoms are not all in the',ERROR)
      call print('WARNING! 0,0,0 cell. If you have checked it, please ignore this message.',ERROR)
      cur_qmmm_qm_abc = qmmm_qm_abc(at, run_type_i, qm_vacuum)
      call append_line(cp2k_template_additions, "&FORCE_EVAL&QMMM&CELL ABC " // cur_qmmm_qm_abc, additions_n_lines)
      call append_line(cp2k_template_additions, "&FORCE_EVAL&QMMM&CELL PERIODIC XYZ", additions_n_lines)
      call append_line(cp2k_template_additions, "&FORCE_EVAL&QMMM &END CELL", additions_n_lines)

      if (get_value(at%params, "QM_cell", old_qmmm_qm_abc)) then
	if (cur_qmmm_qm_abc .fne. old_qmmm_qm_abc) can_reuse_wfn = .false.
      endif

      if (get_value(at%params, "QM_list_changed", qm_list_changed)) then
	if (qm_list_changed) can_reuse_wfn = .false.
      endif

call print("qm_list_a " // qm_list_a, ERROR)
call print("at%Z " // at%Z, ERROR)
call print("at%Z min max " // minval(at%Z) // " " //  maxval(at%Z), ERROR)

      counter = 0
      do atno=minval(at%Z), maxval(at%Z)
call print("looking for QM_KIND Z="//atno, ERROR)
	if (any(at%Z(qm_list_a) == atno)) then
call print("    found at least one atom", ERROR)
	  call append_line(cp2k_template_additions, "&FORCE_EVAL&QMMM &QM_KIND "//ElementName(atno), additions_n_lines)
	  do i=1, qm_list%N
	    if (at%Z(qm_list_a(i)) == atno) then
call print("    adding atom " // qm_list_a(i), ERROR)
	      call append_line(cp2k_template_additions, "&FORCE_EVAL&QMMM&QM_KIND-"//trim(ElementName(atno))// &
							" MM_INDEX "//qm_list_a(i), additions_n_lines)
	      counter = counter + 1
	    endif
	  end do
	  call append_line(cp2k_template_additions, "&FORCE_EVAL&QMMM &END QM_KIND", additions_n_lines)
	end if
      end do
      if (size(qm_list_a) /= counter) &
	call system_abort("Number of QM list atoms " // size(qm_list_a) // " doesn't match number of QM_KIND atoms " // counter)
    endif

    ! put in things needed for QM
    if (use_QM) then
      if (try_reuse_wfn .and. can_reuse_wfn) then 
	call append_line(cp2k_template_additions, "&FORCE_EVAL&DFT WFN_RESTART_FILE_NAME ../wfn.restart.wfn", additions_n_lines)
	call append_line(cp2k_template_additions, "&FORCE_EVAL&DFT &SCF", additions_n_lines)
	call append_line(cp2k_template_additions, "&FORCE_EVAL&DFT&SCF SCF_GUESS RESTART", additions_n_lines)
	call append_line(cp2k_template_additions, "&FORCE_EVAL&DFT &END SCF", additions_n_lines)
      endif
      call calc_charge_lsd(at, qm_list, charge, do_lsd)
      call append_line(cp2k_template_additions, "&FORCE_EVAL&DFT CHARGE "//charge, additions_n_lines)
      if (do_lsd) call append_line(cp2k_template_additions, "&FORCE_EVAL&DFT LSD", additions_n_lines)
    endif

    ! put in unit cell
    call append_line(cp2k_template_additions, "&FORCE_EVAL &SUBSYS", additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS &CELL", additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&CELL A " // at%lattice(:,1), additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&CELL B " // at%lattice(:,2), additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&CELL C " // at%lattice(:,3), additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS &END CELL", additions_n_lines)

    ! put in topology
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS &TOPOLOGY", additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY &DUMP_PSF", additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY &END DUMP_PSF", additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY &DUMP_PDB", additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY &END DUMP_PDB", additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY &GENERATE", additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE REORDER F", additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE CREATE_MOLECULES F", additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE &ISOLATED_ATOMS", additions_n_lines)
    if (use_QMMM) then
      do i=1, size(qm_list_a)
	call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE&ISOLATED_ATOMS " // qm_list_a(i), additions_n_lines)
      end do
    endif
    if (assign_pointer(at, "isolated_atom", isolated_atom)) then
      do i=1, at%N
	if (isolated_atom(i) /= 0) then
	  call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE&ISOLATED_ATOMS " // i, additions_n_lines)
	endif
      end do
    endif
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE &END ISOLATED_ATOMS", additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY &END GENERATE", additions_n_lines)

    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY COORD_FILE_NAME quip_cp2k.exyz", additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY COORDINATE EXYZ", additions_n_lines)
    if (trim(psf_print) == "DRIVER_PRINT_AND_SAVE" .or. trim(psf_print) == "USE_EXISTING_PSF") then
      call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY CONN_FILE_NAME quip_cp2k.psf", additions_n_lines)
      call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS&TOPOLOGY CONN_FILE_FORMAT PSF", additions_n_lines)
    endif
    call append_line(cp2k_template_additions, "&FORCE_EVAL&SUBSYS &END TOPOLOGY", additions_n_lines)
    call append_line(cp2k_template_additions, "&FORCE_EVAL &END SUBSYS", additions_n_lines)

    ! put in global stuff to run a single force evalution, print out appropriate things
    call append_line(cp2k_template_additions, " &GLOBAL", additions_n_lines)
    call append_line(cp2k_template_additions, "&GLOBAL   PROJECT quip", additions_n_lines)
    call append_line(cp2k_template_additions, "&GLOBAL   RUN_TYPE MD", additions_n_lines)
    call append_line(cp2k_template_additions, " &END GLOBAL", additions_n_lines)
    call append_line(cp2k_template_additions, " &MOTION", additions_n_lines)
    call append_line(cp2k_template_additions, "&MOTION   &PRINT", additions_n_lines)
    call append_line(cp2k_template_additions, "&MOTION   &END PRINT", additions_n_lines)
    call append_line(cp2k_template_additions, "&MOTION&PRINT     &FORCES", additions_n_lines)
    call append_line(cp2k_template_additions, "&MOTION&PRINT     &END FORCES", additions_n_lines)
    call append_line(cp2k_template_additions, "&MOTION&PRINT&FORCES       FORMAT XMOL", additions_n_lines)
    call append_line(cp2k_template_additions, "&MOTION   &MD", additions_n_lines)
    call append_line(cp2k_template_additions, "&MOTION   &END MD", additions_n_lines)
    call append_line(cp2k_template_additions, "&MOTION&MD     ENSEMBLE NVE", additions_n_lines)
    call append_line(cp2k_template_additions, "&MOTION&MD     STEPS 0", additions_n_lines)
    call append_line(cp2k_template_additions, " &END MOTION", additions_n_lines)

call print("cp2k_template_additions")
do i=1, additions_n_lines
  call print(trim(cp2k_template_additions(i)), ERROR)
end do

    ! merge template and QUIP stuff
    call merge_templates(cp2k_template_a(1:template_n_lines), cp2k_template_additions(1:additions_n_lines), &
      cp2k_template_combined, combined_n_lines)

call print("cp2k_template_combined")
do i=1, combined_n_lines
  call print(trim(cp2k_template_combined(i)), ERROR)
end do

    if (run_type /= "QS") then
      if (trim(psf_print) == "DRIVER_PRINT_AND_SAVE") then
	call write_psf_file(at, "quip_cp2k.psf", run_type_string=trim(run_type))
      endif
    endif

    run_dir = make_run_directory()

    call write_cp2k_input_file(cp2k_template_combined(1:combined_n_lines), trim(run_dir)//'/cp2k_input.inp')

    call run_cp2k(trim(cp2k_program), trim(run_dir), max_n_tries)

    !! parse output
    call read_energy_forces(at, qm_list, cur_qmmm_qm_abc, trim(run_dir), "quip", e, f)

    if (maxval(abs(f)) > max_force_warning) &
      call print('WARNING cp2k forces max component ' // maxval(abs(f)) // ' at ' // maxloc(abs(f)) // &
		 ' exceeds warning threshold ' // max_force_warning, ERROR)

    !! save output

    if (trim(psf_print) == "CP2K_PRINT_AND_SAVE") &
      call system_command('cp '//trim(run_dir)//'/quip-dump-1.psf quip_cp2k.psf')

    if (use_QM) &
      call system_command('cp '//trim(run_dir)//'quip-RESTART.wfn wfn.restart.wfn')

    if (save_output_files) then
      call system_command('cat '//trim(run_dir)//'/cp2k_input.inp'// &
        ' >> cp2k_input_log; echo "##############" >> cp2k_input_log;' // &
        ' cat '//trim(run_dir)//'/quip-frc-1.xyz'// &
        ' >> cp2k_force_file_log; echo "##############" >> cp2k_force_file_log;' // &
        ' cat '//trim(run_dir)//'/cp2k_output.out >> cp2k_output_log; echo "##############" >> cp2k_output_log')
    endif

    if (clean_up_files) call system_command('rm -rf '//trim(run_dir))

    call system_timer('go_cp2k_template')

  end subroutine go_cp2k_template

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
	if (norm(at%pos(:,i)+shift-new_p(:,j)) .feq. 0.0_dp) oldpos(i) = j
      end do
    end do

    ! try again with shift of a/2 b/2 c/2 in case TOPOLOGY%CENTER_COORDINATES is set
    if (any(oldpos == 0)) then
      shift = sum(at%lattice(:,:),2)/2.0_dp - &
	      (minval(at%pos(:,:),2)+maxval(at%pos(:,:),2))/2.0_dp
      do i=1, at%N
	do j=1, size(new_p,2)
	  if (norm(at%pos(:,i) + shift- new_p(:,j)) .feq. 0.0_dp) oldpos(i) = j
	end do
      end do
    endif

    if (any(oldpos == 0)) &
      call system_abort("Could not match 2 atom objects")

    new_f(1,oldpos(:)) = new_f(1,:)
    new_f(2,oldpos(:)) = new_f(2,:)
    new_f(3,oldpos(:)) = new_f(3,:)

  end subroutine reorder_if_necessary

  subroutine run_cp2k(cp2k_program, run_dir, max_n_tries)
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

  end subroutine run_cp2k

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

  subroutine merge_templates(base, additions, combined, combined_n_lines)
    character(len=*), intent(in) :: base(:), additions(:)
    character(len=*), allocatable, intent(inout) :: combined(:)
    integer, intent(out) :: combined_n_lines

    integer :: i, j, sec_end
    logical :: already_present
    character(len=1024) :: sec_c, word_c, arg_c, sec_a, word_a, arg_a

    ! copy base template into combined
    combined_n_lines = 0
    do i=1, size(base)
      call append_line(combined, trim(base(i)), combined_n_lines)
    end do

    ! for each line in additions
    do i=1, size(additions)
call print("merging line '"//trim(additions(i)), ERROR)
      call split_cp2k_input_line(trim(additions(i)), sec_a, word_a, arg_a)
      already_present = .false.
      do j=1, combined_n_lines
	call split_cp2k_input_line(trim(combined(j)), sec_c, word_c, arg_c)
	if (trim(sec_a) == trim(sec_c) .and. trim(word_a) == trim(word_c)) then
	  if (trim(arg_a) == trim(arg_c)) then 
	    already_present = .true.
!	  else
!	    if (word_a(1:1) /= "&") &
!	      call system_abort("ERROR: QUIP inserting line '"//trim(additions(i))//"' already present with different argument in template '"//trim(combined(j)))
	  endif
	endif
      end do
      if (already_present) cycle

      if (len_trim(sec_a) == 0) then ! top level
	call append_line(combined, trim(additions(i)), combined_n_lines)
      else ! within an section
	sec_end = find_section_end(combined(1:combined_n_lines), sec_a)
	if (sec_end <= 0) &
	  call system_abort("merge_templates failed to find section '"//trim(sec_a)//"'")

	call insert_line(combined, trim(additions(i)), after_line=sec_end, n_l=combined_n_lines)
      endif
    end do

  end subroutine merge_templates

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

! call print("splitting '"//trim(l) // " into '" // trim(sec) // "' '" // trim(word)//"' '"//trim(arg)//"'", ERROR)
  end subroutine split_cp2k_input_line

  function find_section_end(l_a, sec) result(line)
    character(len=*), intent(in) :: l_a(:)
    character(len=*) :: sec
    integer :: line

    character(len=1024) :: sec_base, sec_tail, l_sec, l_word
    integer :: i, j, pamp, pspc

call print("find_section_end sec '"//trim(sec)//"'", ERROR)
    pamp = index(trim(sec), "&", back=.true.)
    if (pamp == 1) then
      sec_base=""
    else
      sec_base = sec(1:pamp-1)
    endif
    sec_tail=sec(pamp:len_trim(sec))

call print(" sec_base '"//trim(sec_base)//"'", ERROR)
call print(" sec_tail '"//trim(sec_tail)//"'", ERROR)

    line = 0

    do i=size(l_a), 1, -1
call print("check line '"//trim(l_a(i))//"'")
      if (l_a(i)(1:len_trim(sec)) == trim(sec)) then
call print("found last line that matches this section " // i // " '"//trim(l_a(i))//"'", ERROR)
	line = i
	return
      endif
      pspc = index(l_a(i), " ")
      l_sec = l_a(i)(1:pspc-1)
      l_word = trim(adjustl(l_a(i)(pspc+1:len_trim(l_a(i)))))
      do j=1, len_trim(l_word)
	if (l_word(j:j) == " ") l_word(j:j) = "-"
      end do
call print ("this line's section is '"//trim(l_sec)//"' and word '"//trim(l_word), ERROR)
      if (trim(l_sec) == trim(sec_base) .and. trim(l_word) == trim(sec_tail)) then
call print("found section beginning", ERROR)
	line = i
	return
      endif
    end do

  end function find_section_end

  subroutine insert_line(l_a, l, after_line, n_l)
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
  end subroutine insert_line

  subroutine prefix_sections(l_a)
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

  subroutine append_line(l_a, l, n_l)
    character(len=*), allocatable, intent(inout) :: l_a(:)
    character(len=*), intent(in) :: l
    integer, intent(inout) :: n_l

    if (.not. allocated(l_a)) call extend_char_array(l_a)
    if (n_l+1 > size(l_a)) call extend_char_array(l_a)

    l_a(n_l+1) = trim(l)
    n_l = n_l + 1

  end subroutine append_line

  subroutine get_qm_list(at,qmflag,qm_list,do_union)
    type(Atoms), intent(in)  :: at
    integer,     intent(in)  :: qmflag
    type(Table), intent(out) :: qm_list
    logical, intent(in), optional :: do_union
  
    integer :: i
    integer, pointer :: QM_flag(:)
    logical              :: my_do_union
    integer :: qmflag_min

    my_do_union = optional_default(.false.,do_union)

    if (.not.(assign_pointer(at, "QM_flag", QM_flag))) &
      call system_abort("get_qm_list couldn't find QM_flag field")

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
      qm_maxdist(1) = max(qm_maxdist(1), distance_min_image(at, (/ at%pos(1,i), 0.0_dp, 0.0_dp /), (/at%pos(1,j), 0.0_dp, 0.0_dp /) ))
      qm_maxdist(2) = max(qm_maxdist(2), distance_min_image(at, (/ 0.0_dp, at%pos(2,i), 0.0_dp /), (/0.0_dp, at%pos(2,j), 0.0_dp /)))
      qm_maxdist(3) = max(qm_maxdist(3), distance_min_image(at, (/ 0.0_dp, 0.0_dp, at%pos(3,i) /), (/ 0.0_dp, 0.0_dp, at%pos(3,j) /)))
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
