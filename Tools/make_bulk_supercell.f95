program make_bulk_supercell
use libatoms_module
implicit none

  character(len=FIELD_LENGTH) :: struct, outfile, Z_values_str
  real(dp) :: vol_per_atom, vol_per_unit_cell
  integer :: repeat(3)
  type(Dictionary) :: cli_params

  type(Atoms) :: dup_cell

  call system_initialise(verbosity=PRINT_SILENT)

  call initialise(cli_params)
  call param_register(cli_params,"struct", PARAM_MANDATORY, struct, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"outfile", "stdout", outfile, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"vol_per_atom", "-1.0", vol_per_atom, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"vol_per_unit_cell", "-1.0", vol_per_unit_cell, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"repeat", "1 1 1", repeat, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"Z_values", "", Z_values_str, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli_params)) then
    call print("Usage: make_bulk_supercell struct=[struct_name] outfile=[filename](stdout)", PRINT_ALWAYS)
    call print("       [ vol_per_atom=volume | vol_per_unit_cell=volume ] [ repeat='n1 n2 n3'(1 1 1) ]", PRINT_ALWAYS)
    call print("       [ Z_values='Z1 Z2 ...' ]", PRINT_ALWAYS)
    call print("In addition, struct names that do not begin with . or / will be searched for", PRINT_ALWAYS)
    call print("  in $QUIP_DIR/structures/ or $HOME/share/quip_structures/, in that order", PRINT_ALWAYS)
    call system_abort("Failed to parse command line arguments")
  endif
  call finalise(cli_params)

  dup_cell = structure_from_file(struct, vol_per_atom, vol_per_unit_cell, repeat, Z_values_str)
  call verbosity_push(PRINT_NORMAL)
  call write(dup_cell, outfile)
  call verbosity_pop()

  call system_finalise()

end program
