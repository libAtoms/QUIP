! Reads a sequence of configs, optionally with decimation, min_time, max_time, and sorting
!  by Time as well as eliminating duplicates (basically a driver for atoms_ll_read_xyz_filename())

program clean_traj
use libatoms_module
implicit none

  type(Dictionary) :: cli_params
  character(len=FIELD_LENGTH) :: infilename
  logical :: infile_is_list
  character(len=FIELD_LENGTH) :: outfilename
  type(CInOutput) :: outfile
  integer :: decimation
  real(dp) :: min_time, max_time
  logical :: sort_Time, no_Time_dups, quiet, no_compute_index

  type(Atoms_ll) :: structure_ll
  type(Atoms_ll_entry), pointer :: entry

  call system_initialise(NORMAL)

  call initialise(cli_params)
  call param_register(cli_params, 'infile', param_mandatory, infilename)
  call param_register(cli_params, 'infile_is_list', 'F', infile_is_list)
  call param_register(cli_params, 'outfile', 'stdout', outfilename)
  call param_register(cli_params, 'decimation', '1', decimation)
  call param_register(cli_params, 'min_time', '-1.0', min_time)
  call param_register(cli_params, 'max_time', '-1.0', max_time)
  call param_register(cli_params, 'sort_Time', 'F', sort_Time)
  call param_register(cli_params, 'no_Time_dups', 'F', no_Time_dups)
  call param_register(cli_params, 'quiet', 'F', quiet)
  if (.not. param_read_args(cli_params, do_check = .true.)) then
      call print_usage()
    call system_abort('could not parse argument line')
  end if
  call finalise(cli_params)

  call print("infile " // trim(infilename) // " infile_is_list " // infile_is_list)
  call print("outfilename " // trim(outfilename))
  call print("decimation " // decimation // " min_time " // min_time // " max_time " // max_time)
  call print("sort_Time " // sort_Time // " no_Time_dups " // no_Time_dups // " quiet " // quiet)

  no_compute_index=.false.
  if ((decimation == 1) .and. (.not. sort_Time) .and. (.not. no_Time_dups)) no_compute_index=.true.

  call read_xyz(structure_ll, infilename, infile_is_list, decimation, min_time, max_time, sort_Time, no_Time_dups, quiet, no_compute_index)

  call initialise(outfile, outfilename, action=OUTPUT)
  entry => structure_ll%first
  do while (associated(entry))
    call write(outfile, entry%at)
    entry => entry%next
  end do
  call finalise(outfile)

  call system_finalise()

contains

  subroutine print_usage()

    character(len=1024) my_exec_name

    if (EXEC_NAME == "<UNKNOWN>") then
      my_exec_name=""
    else
      my_exec_name=EXEC_NAME
    endif

    call print("Usage: " // trim(my_exec_name)//" infile=filename [infile_is_list=logical(F)]", ERROR)
    call print("       outfile=filename [decimation=n(1)]", ERROR)
    call print("       [min_time=t(-1.0)] [max_time=t(-1.0)]", ERROR)
    call print("       [sort_Time(F)] [no_Time_dups(F)] [quiet(F)]", ERROR)
  end subroutine print_usage

end program clean_traj
