! print 1 out of every n frames
program decimate
use libatoms_module
implicit none
  type(Atoms) :: at
  integer :: i, n, stat
  type(dictionary) :: cli_params
  character(len=FIELD_LENGTH) :: filename
  type(inoutput) :: xyzfile

  call system_initialise(verbosity=SILENT)
  call verbosity_push(NORMAL)

  call initialise(cli_params)
  call param_register(cli_params,"n",PARAM_MANDATORY, n)
  call param_register(cli_params,"file","stdin", filename)
  if (.not. param_read_args(cli_params, do_check = .true.)) then
    call system_abort("Usage: decimate n=(1) [file=(stdin)]")
  endif
  call finalise(cli_params)

  call initialise(xyzfile, filename, INPUT)

  call read_xyz(at, xyzfile, status=stat)
  call print_xyz(at, mainlog, all_properties=.true.)
  do while (stat == 0)
    do i=1, n-1
      call read_xyz(xyzfile, status=stat)
      if (stat /= 0) exit
    end do
    if (stat == 0) then
      call read_xyz(at, xyzfile, status=stat)
      if (stat == 0) call print_xyz(at, mainlog, all_properties=.true.)
    endif
  end do

  call verbosity_pop()
  call system_finalise()
end program decimate
