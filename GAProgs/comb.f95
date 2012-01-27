program comb

  use libatoms_module

  implicit none

  type(atoms) :: at
  type(cinoutput) :: infile, outfile
  type(Dictionary) :: params

  integer  :: i, n
  
  character(len=STRING_LENGTH) :: in_file, out_file

  call system_initialise(verbosity=PRINT_SILENT)

  call initialise(params)
  call param_register(params, 'in_file', PARAM_MANDATORY, in_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'out_file', 'stdout', out_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'n', '20', n, help_string="No help yet.  This source file was $LastChangedBy$")

  if (.not. param_read_args(params)) then
     call verbosity_push(PRINT_NORMAL)
     call print("Usage: comb [in_file=file] [out_file=stdout] [n=20]")
     call system_abort('Exit: Mandatory argument(s) missing...')
  endif
  call finalise(params)

  call initialise(infile,trim(in_file))
  call initialise(outfile,trim(out_file),action=OUTPUT)

  do i = n, infile%n_frame, n
     call read(infile,at,frame=i-1)
     call write(at,outfile)
     call finalise(at)
  enddo

  call finalise(infile)
  call finalise(outfile)

end program comb
