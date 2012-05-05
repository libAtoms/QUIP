#include "error.inc"
program vasp_driver
use libatoms_module
implicit none

   type(Atoms) :: at

   character(len=STRING_LENGTH) :: infile, outfile
   character(len=STRING_LENGTH) :: args_str
   character(len=STRING_LENGTH)  :: arg
   integer :: i, index_insert
   type(CInOutput) :: outfile_io

   integer :: error

   call system_initialise(verbosity=PRINT_SILENT, enable_timing=.true.)
   call verbosity_push(PRINT_NORMAL)
   call system_timer('vasp_filepot')

   if (cmd_arg_count() < 2) then
      call system_abort("Usage: vasp_filepot infile outfile [ other_cli_arg1=v1 other_cli_arg2=v2 ...]")
   endif

   call get_cmd_arg(1, infile)
   call get_cmd_arg(2, outfile)

    args_str = ""
    if (cmd_arg_count() > 2) then
      do i=3, cmd_arg_count()
	call get_cmd_arg(i, arg)
        !add {} if there is space in the arg
        if (index(trim(arg)," ").ne.0) then
            index_insert = index(trim(arg),"=")
            arg(index_insert+1:len_trim(arg)+2) = "{"//arg(index_insert+1:len_trim(arg))//"}"
            call print('arg: '//trim(arg),PRINT_SILENT)
        endif
	args_str = trim(args_str) // " " // trim(arg)
      end do
    endif

    call read(at,infile)

    call print('vasp_filepot args_str '//trim(args_str), PRINT_ALWAYS)

    call do_vasp_calc(at, trim(args_str), error)
    HANDLE_ERROR(error)

    call initialise(outfile_io,outfile,action=OUTPUT)
    call write(at, outfile_io, real_format='%20.13f')
    call finalise(outfile_io)

    call finalise(at)

    call system_timer('vasp_filepot')
    mainlog%prefix=""
    call verbosity_push(PRINT_SILENT)
    call system_finalise()

contains

subroutine do_vasp_calc(at, args_str, error)
   type(Atoms), intent(inout) :: at
   character(len=*), intent(in) :: args_str
   integer, intent(out), optional :: error

   character(len=STRING_LENGTH) :: incar_template_file, kpoints_file, potcar_files, vasp_path, run_suffix, verbosity_str
   character(len=STRING_LENGTH) :: calc_energy, calc_force, calc_virial, calc_local_energy
   logical :: do_calc_energy, do_calc_force, do_calc_virial, do_calc_local_energy, persistent, persistent_n_atoms_exist
   logical :: clean_up_files, ignore_convergence, no_use_WAVECAR, force_constant_basis
   integer :: persistent_restart_interval, persistent_n_atoms
   real(dp) :: persistent_max_wait_time
   type(Dictionary) :: incar_dict, cli

   integer :: persistent_run_i
   logical :: persistent_already_started = .false.
   type(inoutput) :: persistent_io

   integer :: run_dir_i, force_run_dir_i

   integer :: i_charge
   real(dp) :: nelect, nelect_of_elem(size(ElementName)), r_charge

   character(len=STRING_LENGTH) :: run_dir
   type(Inoutput) :: io

   integer :: at_i, stat
   integer, pointer :: sort_index_p(:)
   logical :: converged, have_wavecar

   INIT_ERROR(error)

   call initialise(cli)
   call param_register(cli, "verbosity", "NORMAL", verbosity_str, help_string="verbosity level")
   call param_register(cli, "INCAR_template", "INCAR.template", incar_template_file, help_string="name of file with template for main input (INCAR) file")
   call param_register(cli, "kpoints_file", "KPOINTS", kpoints_file, help_string="name of k-points file (KPOINTS)")
   call param_register(cli, "potcar_files", PARAM_MANDATORY, potcar_files, help_string="mapping from atomic numbers to POTCAR files, Z_1 POTCAR_FILE_1 [ Z_2 POTCAR_FILE_2 ... ]")
   call param_register(cli, "vasp", PARAM_MANDATORY, vasp_path, help_string="vasp executable")
   call param_register(cli, "energy", "", calc_energy, help_string="if set, calculate energy and put it in parameter named this")
   call param_register(cli, "force", "", calc_force, help_string="if set, calculate forces and put them in property named this")
   call param_register(cli, "virial", "", calc_virial, help_string="if set, calculate virial and put it in parameter named this")
   call param_register(cli, "local_energy", "", calc_local_energy, help_string="not supported")
   call param_register(cli, "clean_up_files", "T", clean_up_files, help_string="if true, delete run directory when done")
   call param_register(cli, "ignore_convergence", "F", ignore_convergence, help_string="If true, accept result even if failed to converge SCF")
   call param_register(cli, "no_use_WAVECAR", "F", no_use_WAVECAR, help_string="If true, don't reuse previous wavefunction file WAVECAR")
   call param_register(cli, "force_constant_basis", "F", force_constant_basis, help_string="If true, do a restart with the same basis")
   call param_register(cli, "run_suffix", "run", run_suffix, help_string="suffix to label persistent state from run to run")
   call param_register(cli, "persistent", "F", persistent, help_string="If true, use persistent connection to a long-running vasp process")
   call param_register(cli, "persistent_restart_interval", "100", persistent_restart_interval, help_string="If persistent, restart vasp this often")
   call param_register(cli, "persistent_max_wait_time", "600", persistent_max_wait_time, help_string="If persistent, never wait more than this long (s) for calculation to finish")
   if (.not. param_read_line(cli, args_str, ignore_unknown=.true.,task='do_vasp_calc args_str')) then
      call print("Args:  verbosity INCAR_template kpoints_file potcar_files vasp energy force", PRINT_ALWAYS)
      call print("       virial local_energy clean_up_files ignore_convergence no_use_WAVECAR", PRINT_ALWAYS)
      call print("       force_constant_basis run_suffix", PRINT_ALWAYS)
      RAISE_ERROR('do_vasp_calc failed to parse args_str="'//trim(args_str)//'"', error)
   endif
   call finalise(cli)

   call verbosity_push(verbosity_of_str(trim(verbosity_str)))

   do_calc_energy = len_trim(calc_energy) > 0
   do_calc_force = len_trim(calc_force) > 0
   do_calc_virial = len_trim(calc_virial) > 0
   do_calc_local_energy = len_trim(calc_local_energy) > 0

   call print("do_vasp_calc using args:", PRINT_VERBOSE)
   call print("  INCAR_template " // trim(INCAR_template_file), PRINT_VERBOSE)
   call print("  kpoints_file " // trim(kpoints_file), PRINT_VERBOSE)
   call print("  potcar_files " // trim(potcar_files), PRINT_VERBOSE)
   call print("  vasp " // trim(vasp_path), PRINT_VERBOSE)
   call print("  energy " // do_calc_energy, PRINT_VERBOSE)
   call print("  force " // do_calc_force, PRINT_VERBOSE)
   call print("  virial " // do_calc_virial, PRINT_VERBOSE)
   call print("  local_energy " // do_calc_local_energy, PRINT_VERBOSE)
   call print("  clean_up_files " // clean_up_files, PRINT_VERBOSE)
   call print("  ignore_convergence " // ignore_convergence, PRINT_VERBOSE)
   call print("  no_use_WAVECAR " // no_use_WAVECAR, PRINT_VERBOSE)
   call print("  force_constant_basis " // force_constant_basis, PRINT_VERBOSE)
   call print("  run_suffix " // run_suffix, PRINT_VERBOSE)
   call print("  persistent " // persistent, PRINT_VERBOSE)
   if (persistent) then
      call print("  persistent_restart_interval " // persistent_restart_interval, PRINT_VERBOSE)
      call print("  persistent_max_wait_time " // persistent_max_wait_time, PRINT_VERBOSE)
   endif

   call system_timer('vasp_filepot/prep')
   ! check for local energy
   if (do_calc_local_energy) then
      RAISE_ERROR("do_vasp_calc can't do local_energy", error)
   endif

   call print("reading template", PRINT_VERBOSE)
   ! read the template incar
   call read_vasp_incar_dict(incar_dict, trim(incar_template_file), error)
   PASS_ERROR(error)

   call print("setting INCAR dictionary parameters", PRINT_VERBOSE)
   ! check for stress, set ISIF
   if (.not. do_calc_energy .and. .not. do_calc_force .and. .not. do_calc_virial) then
      return
   else if (.not. do_calc_virial) then
      call set_value(incar_dict, "ISIF", 0)
   else
      if (persistent) then
	 call set_value(incar_dict, "ISIF", 3)
      else
	 call set_value(incar_dict, "ISIF", 2)
      endif
   endif

   force_run_dir_i = -1
   if (persistent) then ! make run_dir manually
      inquire(file="persistent_run_i", exist=persistent_already_started)
      run_dir = "vasp_run_0"
      force_run_dir_i = 0
      if (.not. persistent_already_started) then
	 call system_command("mkdir "//trim(run_dir))
	 call initialise(persistent_io, "persistent_run_i", action=OUTPUT)
	 call print("-1", file=persistent_io)
	 call finalise(persistent_io)
      endif

      inquire(file="persistent_n_atoms", exist=persistent_n_atoms_exist)
      if (persistent_n_atoms_exist) then
         call initialise(persistent_io, "persistent_n_atoms", action=INPUT)
         line = read_line(persistent_io)
         read (unit=line, fmt=*, iostat=stat) persistent_n_atoms
         call finalise(persistent_io)
      else
         persistent_n_atoms = 0
      end if
      call print('got old persistent_n_atoms = '//persistent_n_atoms, PRINT_VERBOSE)
      call print('writing new persistent_n_atoms = '//at%n, PRINT_VERBOSE)
      call initialise(persistent_io, "persistent_n_atoms", action=OUTPUT)
      call print(at%n, file=persistent_io)
      call finalise(persistent_io)
   endif
   run_dir = make_run_directory("vasp_run", force_run_dir_i, run_dir_i)

   if (persistent) then
   ! set dynamics to REFTRAJ, number of ionic steps to 1000000
      call set_value(incar_dict, "NSW", 1000000)
      call set_value(incar_dict, "IBRION", 12)
   else
   ! set dynamics to MD, number of ionic steps to 0
      call set_value(incar_dict, "NSW", 0)
      call set_value(incar_dict, "IBRION", 0)
   endif
   if (no_use_WAVECAR) then
      call set_value(incar_dict, "ISTART", 0)
   else
      inquire(file="WAVECAR."//trim(run_suffix), exist=have_wavecar)
      if (have_wavecar) then
	 call print("copying old WAVECAR.run_suffix to run_dir/WAVECAR", PRINT_VERBOSE)
	 call system_command("cp WAVECAR."//trim(run_suffix)//" "//trim(run_dir)//"/WAVECAR", status=stat)
	 if (stat /= 0) then
	    RAISE_ERROR("do_vasp_calc copying file WAVECAR."//trim(run_suffix)//" into "//trim(run_dir), error)
	 endif
      endif
      if (force_constant_basis) then
	 ! 3 is faster, but not tested
	 call set_value(incar_dict, "ISTART", 2)
      else
	 call set_value(incar_dict, "ISTART", 1)
      endif
   endif

   call print("doing sorting", PRINT_VERBOSE)
   call add_property(at, "sort_index", 0, n_cols=1, ptr=sort_index_p, error=error)
   PASS_ERROR(error)
   do at_i=1, at%N
      sort_index_p(at_i) = at_i
   end do
   call atoms_sort(at, 'Z', error=error)
   PASS_ERROR(error)

   if (persistent) then
      call print("reading persistent_run_i", PRINT_VERBOSE)
      call initialise(persistent_io, "persistent_run_i", action=INPUT)
      line = read_line(persistent_io)
      read (unit=line, fmt=*, iostat=stat) persistent_run_i
      call finalise(persistent_io)
      call print("got persistent_run_i="//persistent_run_i, PRINT_VERBOSE)
      persistent_run_i = persistent_run_i + 1
      call initialise(persistent_io, "persistent_run_i", action=OUTPUT)
      call print(persistent_run_i, file=persistent_io)
      call finalise(persistent_io)
   endif

   if (.not. persistent_already_started) then
      call print("writing POTCAR", PRINT_VERBOSE)
      call write_vasp_potcar(at, trim(run_dir), trim(potcar_files), nelect_of_elem, error=error)
      PASS_ERROR(error)
   endif

   ! total number of electrons for each element from POTCAR files
   nelect = sum(nelect_of_elem(at%Z))
   ! get Charge from xyz file
   r_charge = 0.0_dp
   if (get_value(at%params, 'Charge', r_charge)) then
      nelect = nelect - r_charge
   else if (get_value(at%params, 'Charge', i_charge)) then
      r_charge = i_charge
      nelect = nelect - i_charge
   endif
   if (r_charge .fne. 0.0_dp) then ! set total number of electrons
      call set_value(incar_dict, "NELECT", nelect)
   endif

   if (.not. persistent_already_started) then
      ! write INCAR
      call print("writing INCAR", PRINT_VERBOSE)
      call initialise(io, trim(run_dir)//"/INCAR", action=OUTPUT, error=error)
      PASS_ERROR(error)
      call print(write_string(incar_dict, entry_sep=quip_new_line), file=io)
      call finalise(io)

      ! write KPOINTS if needed
      if (trim(kpoints_file) /= "_NONE_") then
	 call print("writing KPOINTS", PRINT_VERBOSE)
	 call system_command("cp "//trim(kpoints_file)//" "//trim(run_dir)//"/KPOINTS", status=stat)
	 if (stat /= 0) then
	    RAISE_ERROR("do_vasp_calc copying kpoints_file "//trim(kpoints_file)//" into "//trim(run_dir), error)
	 endif
      endif
   endif

   if (.not. persistent_already_started .or. mod(persistent_run_i, persistent_restart_interval) == 0 .or. at%n /= persistent_n_atoms) then
      call print("writing POSCAR", PRINT_VERBOSE)
      call write_vasp_poscar(at, trim(run_dir), error=error)
      PASS_ERROR(error)
   endif
   call system_timer('vasp_filepot/prep')

   if (persistent) then ! restart VASP if needed, write to REFTRAJCAR if needed, initiate calc
      call system_timer('vasp_filepot/start_restart')
      if (mod(persistent_run_i, persistent_restart_interval) == 0 .or. at%n /= persistent_n_atoms) then ! start or restart VASP
	 call print("starting or restarting VASP", PRINT_VERBOSE)
	 if (persistent_already_started) then ! tell currently running vasp to quit
	    call initialise(persistent_io, trim(run_dir)//"/REFTRAJCAR", action=OUTPUT)
	    call print("0", file=persistent_io)
	    call finalise(persistent_io)
	    call initialise(persistent_io, trim(run_dir)//"/REFTRAJ_READY", action=OUTPUT)
	    call print("-", file=persistent_io)
	    call finalise(persistent_io)
	    ! wait for signal that it's finished
	    call system_command("echo 'before wait'; ps auxw | egrep 'vasp|mpi'")
	    call wait_for_file_to_exist(trim(run_dir)//"/vasp.done", max_wait_time=600.0_dp, error=error)
	    call system_command("echo 'after wait'; ls -l "//trim(run_dir))
	    call system_command("echo 'after wait'; ps auxw | egrep 'vasp|mpi'")
	    PASS_ERROR(error)
	    ! clean up
	    call unlink(trim(run_dir)//"/vasp.done")
	    call system_command("cd "//trim(run_dir)//"; rm -f REFTRAJCAR REFTRAJ_READY OUTCAR CONTCAR OSZICAR")
	    call system_command("cat vasp.stdout."//run_dir_i//" >> vasp.stdout")
	 endif
	 ! start vasp in background
	 call          print("cd "//trim(run_dir)//"; bash -c "//'"'//"("//'\"'//trim(vasp_path)//'\"'//" > ../vasp.stdout."//run_dir_i//"; echo done > vasp.done )"//'"'//" 2>&1 &")
	 call system_command("cd "//trim(run_dir)//"; bash -c "//'"'//"("//'\"'//trim(vasp_path)//'\"'//" > ../vasp.stdout."//run_dir_i//"; echo done > vasp.done )"//'"'//" 2>&1 &", status=stat)
	 if (stat /= 0) then
	    RAISE_ERROR("error running "//trim(vasp_path)//", status="//stat, error)
	 endif
	 persistent_already_started = .true.
      endif
      call system_timer('vasp_filepot/start_restart')
      call system_timer('vasp_filepot/write_reftrajcar')
      if (mod(persistent_run_i,persistent_restart_interval) /= 0 .and. at%n == persistent_n_atoms) then ! need to write current coords to REFTRAJCAR
	 call print("writing current coords to REFTRAJCAR", PRINT_VERBOSE)
	 call initialise(persistent_io, trim(run_dir)//"/REFTRAJCAR", action=OUTPUT)
	 call print(at%N, file=persistent_io)
	 call print(at%lattice(:,1), file=persistent_io)
	 call print(at%lattice(:,2), file=persistent_io)
	 call print(at%lattice(:,3), file=persistent_io)
	 do i=1, at%N
	    call print(matmul(at%g, at%pos(:,i)), file=persistent_io)
	 end do
	 call finalise(persistent_io)
	 call print("touching REFTRAJ_READY", PRINT_VERBOSE)
	 call initialise(persistent_io, trim(run_dir)//"/REFTRAJ_READY", action=OUTPUT)
	 call print("-", file=persistent_io)
	 call finalise(persistent_io)
      endif
      call system_timer('vasp_filepot/write_reftrajcar')
      call system_timer('vasp_filepot/wait_reftraj_step_done')
      call print("waiting for output REFTRAJ_STEP_DONE", PRINT_VERBOSE)
      ! wait for output to be ready
      call wait_for_file_to_exist(trim(run_dir)//"/REFTRAJ_STEP_DONE", persistent_max_wait_time, cycle_time=0.1_dp, error=error)
      PASS_ERROR(error)
      call system_timer('vasp_filepot/wait_reftraj_step_done')
   else ! not persistent, just run vasp
      call system_timer('vasp_filepot/run_vasp')
      call print("running vasp", PRINT_VERBOSE)
      call print("cd "//trim(run_dir)//"; bash -c "//'"'//trim(vasp_path)//'"'//" > ../vasp.stdout."//run_dir_i//" 2>&1")
      call system_command("cd "//trim(run_dir)//"; bash -c "//'"'//trim(vasp_path)//'"'//" > ../vasp.stdout."//run_dir_i//" 2>&1", status=stat)
      if (stat /= 0) then
	 RAISE_ERROR("error running "//trim(vasp_path)//", status="//stat, error)
      endif
      call system_command("cat vasp.stdout."//run_dir_i//" >> vasp.stdout")
      call system_timer('vasp_filepot/run_vasp')
   endif

   ! look for errors in stdout and runtime in OUTCAR
   call system_command("fgrep -i 'error' vasp.stdout."//run_dir_i)
   call system_command("fgrep -i 'total cpu time used' "//trim(run_dir)//"/OUTCAR | tail -1; fgrep -c 'LOOP:' "//trim(run_dir)//"/OUTCAR | tail -1")

   call system_timer('vasp_filepot/read_output')
   if (persistent) then ! read REFTRAJ output
      call print("reading output REFTRAJ_OUTPUT", PRINT_VERBOSE)
      ! read output
      call read_vasp_output_persistent(at, run_dir, do_calc_energy, do_calc_force, do_calc_virial, error)
      PASS_ERROR(error)
      converged = .true. ! should really find a way of telling if converged for persistent connections
      ! clean up
      call unlink(trim(run_dir)//"/REFTRAJ_OUTPUT")
      call unlink(trim(run_dir)//"/REFTRAJ_STEP_DONE")
   else ! read regular VASP output (OUTCAR)
      call print("reading vasp output", PRINT_VERBOSE)
      call read_vasp_output(trim(run_dir), do_calc_energy, do_calc_force, do_calc_virial, converged, error=error)
      PASS_ERROR(error)
   endif
   call system_timer('vasp_filepot/read_output')

   call system_timer('vasp_filepot/end')
   if (.not. converged) then
      call print("WARNING: failed to find sign of convergence in OUTCAR", verbosity=PRINT_ALWAYS)
      if (.not. ignore_convergence) then
	 RAISE_ERROR("failed to find sign of convergence in OUTCAR", error)
      endif
   endif

   inquire(file=trim(run_dir)//"/WAVECAR", exist=have_wavecar)
   if (have_wavecar) then
      call print("copying run_dir/WAVECAR to WAVECAR.run_suffix", PRINT_VERBOSE)
      call system_command("cp "//trim(run_dir)//"/WAVECAR WAVECAR."//trim(run_suffix), status=stat)
      if (stat /= 0) then
	 RAISE_ERROR("do_vasp_calc copying file "//trim(run_dir)//"WAVECAR into ./"//trim(run_dir)//"WAVECAR."//trim(run_suffix), error)
      endif
   endif

   call print("doing unsorting", PRINT_VERBOSE)
   call atoms_sort(at, 'sort_index', error=error)
   PASS_ERROR(error)

   if (clean_up_files .and. .not. persistent) then
      call print("cleaning up", PRINT_VERBOSE)
      call system_command("rm -rf "//trim(run_dir))
   endif

   call system_timer('vasp_filepot/end')

   call verbosity_pop()

end subroutine do_vasp_calc

subroutine write_vasp_potcar(at, run_dir, potcar_files, nelect_of_elem, error)
   type(Atoms), intent(in) :: at
   character(len=*), intent(in) :: run_dir, potcar_files
   real(dp), intent(out) :: nelect_of_elem(:)
   integer, intent(out), optional :: error

   integer :: i, Z_i, potcar_files_n_fields
   character(len=STRING_LENGTH) :: potcar_files_fields(128), potcar_files_a(128)
   integer, allocatable :: uniq_Z(:)

   integer :: stat, t_Z
   real(dp) :: t_nelect
   character(len=STRING_LENGTH) :: line
   type(inoutput) :: io

   INIT_ERROR(error)

   call split_string_simple(trim(potcar_files), potcar_files_fields, potcar_files_n_fields, " ", error)
   PASS_ERROR(error)

   if (potcar_files_n_fields <= 0 .or. mod(potcar_files_n_fields,2) /= 0) then
      RAISE_ERROR("write_vasp_potcar didn't find even number of fields (Z1 potcar1 [Z2 potcar2 ...]) in potcar_files='"//trim(potcar_files)//"'", error)
   endif

   call uniq(at%Z, uniq_Z)

   potcar_files_a=""
   do i=1, potcar_files_n_fields, 2
      read(unit=potcar_files_fields(i),fmt=*) Z_i
      potcar_files_a(Z_i) = trim(potcar_files_fields(i+1))
   end do

   call system_command("rm -f nelect_summary; touch nelect_summary")
   call system_command("rm -f "//trim(run_dir)//"/POTCAR")
   do Z_i=1, size(uniq_Z)
      if (uniq_Z(Z_i) == 0) exit
      if (len_trim(potcar_files_a(uniq_Z(Z_i))) == 0) then
	 RAISE_ERROR("write_vasp_potcar could not find a potcar for Z="//uniq_Z(Z_i), error)
      endif
      call system_command("cat "//trim(potcar_files_a(uniq_Z(Z_i)))//" >> "//trim(run_dir)//"/POTCAR")
      call system_command("(echo -n '"//uniq_Z(Z_i)//" '; cat "//trim(potcar_files_a(uniq_Z(Z_i)))//" | fgrep ZVAL | sed -e 's/.*ZVAL[ ]*=[ ]*//' -e 's/ .*//') >> nelect_summary")
   end do

   ! read back nelect_summary
   call initialise(io, "nelect_summary")
   do Z_i=1, size(uniq_Z)
      line = read_line(io, stat)
      read (unit=line, fmt=*) t_Z, t_nelect
      if (t_Z <= 0 .or. t_Z > size(nelect_of_elem)) then
	 RAISE_ERROR("write_vasp_potcar got value out of range reading nelect_summary value t_Z="//t_Z, error)
      endif
      if (t_Z /= uniq_Z(Z_i)) then
	 RAISE_ERROR("write_vasp_potcar got mismatch reading nelect_summary value t_Z="//t_Z//" uniq_Z("//Z_i//")="//uniq_Z(Z_i), error)
      endif
      nelect_of_elem(t_Z) = t_nelect
   end do
   call finalise(io)

end subroutine write_vasp_potcar

subroutine write_vasp_poscar(at, run_dir, error)
   type(Atoms), intent(inout) :: at
   character(len=*), intent(in) :: run_dir
   integer, intent(out), optional :: error

   integer, allocatable :: uniq_Z(:), n_Z(:)
   integer :: i
   logical :: swapped_a1_a2
   real(dp) :: vol, t_lat(3)

   type(inoutput) :: io

   INIT_ERROR(error)

   vol = scalar_triple_product(at%lattice(:,1), at%lattice(:,2), at%lattice(:,3))
   if (vol < 0.0_dp) then
      t_lat(:) = at%lattice(:,1)
      at%lattice(:,1) = at%lattice(:,2)
      at%lattice(:,2) = t_lat(:)
      swapped_a1_a2 = .true.
   else
      swapped_a1_a2 = .false.
   endif

   call initialise(io, trim(run_dir)//"/POSCAR",action=OUTPUT)
   if (swapped_a1_a2) then
      call print("QUIP VASP run "//trim(run_dir)//" swapped a1 and a2", file=io)
   else
      call print("QUIP VASP run "//trim(run_dir), file=io)
   endif
   call print("1.00 ! scale factor", file=io)

   call print(at%lattice(:,1), file=io)
   call print(at%lattice(:,2), file=io)
   call print(at%lattice(:,3), file=io)

   call uniq(at%Z, uniq_Z)
   allocate(n_Z(size(uniq_Z)))
   do i=1, size(uniq_Z)
      n_Z(i) = count(at%Z == uniq_Z(i))
   end do
   call print(""//n_Z//" ! Zs = "//uniq_Z, file=io)
   call print("Selective Dynamics", file=io)
   call print("Cartesian", file=io)
   do i=1, at%N
      call print(at%pos(:,i)//" T T T", file=io)
   end do
   call finalise(io)

end subroutine write_vasp_poscar

subroutine read_vasp_output(run_dir, do_calc_energy, do_calc_force, do_calc_virial, converged, error)
   character(len=*), intent(in) :: run_dir
   logical, intent(in) :: do_calc_energy, do_calc_force, do_calc_virial
   logical, intent(out) :: converged
   integer, intent(out), optional :: ERROR

   type(inoutput) :: outcar_io
   integer :: stat
   character(len=STRING_LENGTH) :: line
   character(len=STRING_LENGTH)  :: fields(100), t_s
   integer :: n_fields, line_i

   real(dp) :: energy, virial(3,3), t_pos(3)
   real(dp), pointer :: force_p(:,:)
   integer :: force_start_line_i, virial_start_line_i

   INIT_ERROR(error)

   if (do_calc_force) then
      if (.not. assign_pointer(at, "force", force_p)) then
	 call add_property(at, "force", 0.0_dp, n_cols=3, ptr2=force_p, error=error)
	 PASS_ERROR(error)
      endif
   endif

   converged = .false.
   force_start_line_i = -1
   call initialise(outcar_io, trim(run_dir)//"/OUTCAR")
   stat = 0
   line_i = 1
   do while (stat == 0)
      line = read_line(outcar_io, stat)
      if (do_calc_energy) then
	 if (index(trim(line),"free  energy") > 0) then ! found energy
	    call split_string_simple(trim(line), fields, n_fields, " ", error=error)
	    PASS_ERROR(error)
	    if (n_fields /= 6) then
	       RAISE_ERROR("read_vasp_output confused by energy line# "//line_i//" in OUTCAR line '"//trim(line)//"'", error)
	    endif
	    if (fields(6) /= "eV") then
	       RAISE_ERROR("read_vasp_output confused by units energy line# "//line_i//" in OUTCAR line '"//trim(line)//"'", error)
	    endif
	    read(unit=fields(5),fmt=*) energy
	    call set_param_value(at, 'energy', energy)
	 endif
      endif
      if (do_calc_force) then
	 if (force_start_line_i > 0 .and. line_i >= force_start_line_i .and. line_i <= force_start_line_i+at%N-1) then
	    read(unit=line,fmt=*) t_pos, force_p(:,line_i-force_start_line_i+1)
	 else
	    if (index(trim(line),"TOTAL-FORCE") > 0) then ! found force
	       call split_string_simple(trim(line), fields, n_fields, " ", error=error)
	       if (fields(3) /= "(eV/Angst)") then
		  RAISE_ERROR("read_vasp_output confused by units force header line# "//line_i//" in OUTCAR line '"//trim(line)//"'", error)
	       endif
	       force_start_line_i = line_i+2
	    endif
	 endif
      endif
      if (do_calc_virial) then
	 if (virial_start_line_i > 0 .and. line_i == virial_start_line_i) then
	    read(unit=line,fmt=*) t_s, virial(1,1), virial(2,2), virial(3,3), virial(1,2), virial(2,3), virial(1,3)
	    virial(2,1) = virial(1,2)
	    virial(3,2) = virial(2,3)
	    virial(3,1) = virial(1,3)
	    call set_param_value(at, 'virial', virial)
	 else
	    if (index(trim(line),"FORCE on cell") > 0) then ! found virial
	       call split_string_simple(trim(line), fields, n_fields, " ", error=error)
	       if (fields(9) /= "(eV):") then
		  RAISE_ERROR("read_vasp_output confused by units virial header line# "//line_i//" in OUTCAR line '"//trim(line)//"'", error)
	       endif
	       virial_start_line_i = line_i+13
	    endif
	 endif
      endif
      if (.not. converged) then
	 if (index(trim(line),"aborting loop because EDIFF is reached") > 0) then ! found convergence of SCF
	    converged = .true.
	 endif
      endif
      ! if (stat == 0) call print("GOT OUTPUT "//trim(line))
      line_i = line_i + 1
   end do
   call finalise(outcar_io)

end subroutine read_vasp_output

subroutine read_vasp_output_persistent(at, run_dir, do_calc_energy, do_calc_force, do_calc_virial, error)
   type(atoms), intent(inout) :: at
   character(len=*), intent(in) :: run_dir
   logical, intent(in) :: do_calc_energy, do_calc_force, do_calc_virial
   integer, optional, intent(out) :: error

   real(dp) :: energy, virial(3,3), force_t(3)
   real(dp), pointer :: force_p(:,:)
   type(inoutput) :: reftraj_output_io
   character(len=STRING_LENGTH) :: line
   integer :: ni

   INIT_ERROR(error)

   if (do_calc_force) then
      if (.not. assign_pointer(at, "force", force_p)) then
	 call add_property(at, "force", 0.0_dp, n_cols=3, ptr2=force_p, error=error)
	 PASS_ERROR(error)
      endif
   endif

   call initialise(reftraj_output_io, trim(run_dir)//"/REFTRAJ_OUTPUT", action=INPUT)

   line = read_line(reftraj_output_io); read(unit=line, fmt=*) ni
   if (ni /= at%N) then
      RAISE_ERROR("read_vasp_output_persistent got NI == "//ni//" != at%N == "//at%N, error)
   endif

   ! read energy
   line = read_line(reftraj_output_io); read(unit=line, fmt=*) energy
   if (do_calc_energy) call set_param_value(at, 'energy', energy)
   ! read force
   do i=1, at%N
      line = read_line(reftraj_output_io); read(unit=line, fmt=*) force_t
      if (do_calc_force) force_p(:,i) = force_t(:)
   end do
   ! read virial
   if (do_calc_virial) then
      line = read_line(reftraj_output_io); read(unit=line, fmt=*) virial(1,1), virial(2,2), virial(3,3), virial(1,2), virial(2,3), virial(3,1)
      virial(2,1) = virial(1,2)
      virial(3,2) = virial(2,3)
      virial(1,3) = virial(3,1)
      call set_param_value(at, "virial", virial)
   endif

   call finalise(reftraj_output_io)
end subroutine read_vasp_output_persistent

subroutine read_vasp_incar_dict(incar_dict, incar_template_file, error)
   type(Dictionary), intent(inout) :: incar_dict
   character(len=*), intent(in) :: incar_template_file
   integer, intent(out), optional :: error

   type(Inoutput) :: incar_io
   character(len=STRING_LENGTH), allocatable :: incar_a(:), incar_line_fields(:), incar_field_fields(:)
   integer :: incar_n_lines, incar_line_n_fields, incar_field_n_fields
   integer :: i, j
   integer :: comment_pos

   INIT_ERROR(error)

   call initialise(incar_io, trim(incar_template_file), INPUT)
   call read_file(incar_io, incar_a, incar_n_lines)
   call finalise(incar_io)

   call initialise(incar_dict)
   allocate(incar_line_fields(100))
   allocate(incar_field_fields(100))
   ! loop over lines
   do i=1, incar_n_lines
      comment_pos = index(trim(incar_a(i)), '#')
      if (comment_pos > 1) then
	 incar_a(i) = incar_a(i)(1:comment_pos-1)
      elseif (comment_pos == 1) then
	 incar_a(i) = ""
      endif
      call split_string(trim(incar_a(i)), ";", "''"//'""', incar_line_fields, incar_line_n_fields, .true.)
      ! loop over ';' separated fields
      do j=1, incar_line_n_fields
	 call split_string(incar_line_fields(j), "=", "''"//'""', incar_field_fields, incar_field_n_fields, .true.)
	 if (incar_field_n_fields > 2) then
	    RAISE_ERROR("read_vasp_incar_dict got more than one '=' in field "//trim(incar_line_fields(j)), error)
	 endif
	 call set_value(incar_dict, trim(adjustl(incar_field_fields(1))), trim(adjustl(incar_field_fields(2))))
      end do
   end do

   deallocate(incar_line_fields)
   deallocate(incar_field_fields)

end subroutine read_vasp_incar_dict

end program vasp_driver
