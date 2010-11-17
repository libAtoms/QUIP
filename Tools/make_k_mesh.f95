program make_k_mesh
use libatoms_module
use tb_kpoints_module
implicit none
  type(KPoints) :: kp
  integer :: i
  type(Dictionary) :: cli_params
  character(len=FIELD_LENGTH) :: outfile
  type(inoutput) :: out_io
  integer :: mesh(3)
  logical :: monkhorst_pack

  call system_initialise(verbosity=PRINT_SILENT)

  call initialise(cli_params)
  call param_register(cli_params, "outfile", "stdout", outfile, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "mesh", "1 1 1", mesh, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "monkhorst_pack", "F", monkhorst_pack, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli_params, ignore_unknown=.false.)) then
    call print("Usage: make_k_mesh outfile=filename(stdout) mesh='n1 n2 n3'(1 1 1 ) monkhorst_pack=L(F)", PRINT_ALWAYS)
    call system_abort("Failed to parse command line arguments")
  endif

  if (any(mesh <= 0)) &
    call system_abort("Invalid mesh " // mesh)
  call initialise(kp, mesh, monkhorst_pack)

  call initialise(out_io, trim(outfile))
  call print('<KPoints N="'//kp%N//'" >', file=out_io)
  do i=1, kp%N
    call print('<point weight="'//kp%weights(i)//'" > '//kp%k_pts(:,i)//' </point>', file=out_io)
  end do
  call print('</KPoints>', file=out_io)
  call finalise(out_io)
  call finalise(kp)

  call system_finalise()

end program make_k_mesh
