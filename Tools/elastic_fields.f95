program elastic

  use libAtoms_module
 
  implicit none

  type(Atoms) :: at
  type(Inoutput) :: infile, outfile
  type(Dictionary) :: params
  integer :: i, n, iostat, nargs
  character(len=2048) :: comment, arg1, arg2, missing_params
  real(dp) :: a, C11, C12, C44, cutoff, nneightol

  call system_initialise(SILENT)
  call verbosity_push(NORMAL)

  call param_register(params, 'a', PARAM_MANDATORY, a)
  call param_register(params, 'C11', PARAM_MANDATORY, C11)
  call param_register(params, 'C12', PARAM_MANDATORY, C12)
  call param_register(params, 'C44', PARAM_MANDATORY, C44)
  call param_register(params, 'cutoff', '5.0', cutoff)
  call param_register(params, 'nneightol', '1.2', nneightol)

  nargs = cmd_arg_count()

  if (nargs < 2) then
     call print('Usage: elastic_fields <infile> <outfile> [params]')
     call print('')
     call print('Parameters and default values are:')
     call param_print(params)
     call verbosity_push(SILENT)
     call system_finalise()
     stop
  end if

  call get_cmd_arg(1, arg1)
  call get_cmd_arg(2, arg2)

  call initialise(infile, arg1)
  call initialise(outfile, arg2, action=OUTPUT)


  if (nargs > 2) then
     if (.not. param_read_args(params, (/ (i, i=3,nargs ) /))) &
          call system_abort('Error parsing command line')
  end if

  if (.not. param_check(params, missing_params)) &
       call system_abort('Missing mandatory parameters: '//missing_params)

  call print('Parameters:')
  call param_print(params)
  call print('')


  ! Loop over all frames in input file
  n = 0
  do
     call Read_xyz(at, infile, comment, iostat)
     if (iostat /= 0) exit !EgOF

     call print('Frame '//n)

     call atoms_set_cutoff(at, cutoff)
     at%nneightol = nneightol

     call elastic_fields(at, a, C11, C12, C44)

     call print_xyz(at, outfile, &
          properties='pos:S_xx_sub1:S_yy_sub1:S_zz_sub1:S_yz:S_xz:S_xy:'//&
                     'Sig_xx:Sig_yy:Sig_zz:Sig_yz:Sig_xz:Sig_xy:'//&
                     'SigEval1:SigEval2:SigEval3:SigEvec1:SigEvec2:SigEvec3:'//&
                     'von_mises_stress:von_mises_strain')
     n = n + 1

  end do
  call print('Done '//n//' frames!')

  call verbosity_push(SILENT)
  call system_finalise

end program elastic
