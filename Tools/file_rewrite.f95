! Reads a sequence of configs and writes them back, mostly to convert from xyz to netcdf and back

program file_rewrite
use libatoms_module
implicit none

  type(Dictionary) :: cli_params
  character(len=FIELD_LENGTH) :: infilename
  character(len=FIELD_LENGTH) :: outfilename
  logical :: netcdf4
  type(CInOutput) :: infile, outfile
  type(Atoms) :: at
  integer i, status

  call system_initialise(NORMAL)

  call initialise(cli_params)
  call param_register(cli_params, 'infile', 'stdin', infilename)
  call param_register(cli_params, 'outfile', 'stdout', outfilename)
  call param_register(cli_params, 'netcdf4', 'F', netcdf4)
  if (.not. param_read_args(cli_params, do_check = .true.)) then
      call print_usage()
    call system_abort('could not parse argument line')
  end if
  call finalise(cli_params)

  call print("infile " // trim(infilename) // " outfile " // trim(outfilename) // " netcdf4 " // netcdf4)

  call initialise(infile, infilename, action=INPUT)
  call initialise(outfile, outfilename, action=OUTPUT, netcdf4=netcdf4)

  i = 1
  call read(infile, at, status=status)
  do while (status == 0) 
    call write(outfile, at)
    if (mod(i,100) == 0) write (*,'(I1,$)') mod(i/100,10)
    call read(infile, at, status=status)
    i = i + 1
  end do

  call finalise(outfile)
  call finalise(infile)
  call system_finalise()

contains

  subroutine print_usage()

    character(len=1024) my_exec_name

    if (EXEC_NAME == "<UNKNOWN>") then
      my_exec_name=""
    else
      my_exec_name=EXEC_NAME
    endif

    call print("Usage: " // trim(my_exec_name)//" infile=filename(stdin) outfile=filename(stdout) [netcdf4 (for output)]", ERROR)
  end subroutine print_usage

end program file_rewrite
