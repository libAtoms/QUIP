! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   Portions of this code were written by
! H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! H0 X
! H0 X   Copyright 2006-2010.
! H0 X
! H0 X   These portions of the source code are released under the GNU General
! H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
! H0 X
! H0 X   If you would like to license the source code under different terms,
! H0 X   please contact Gabor Csanyi, gabor@csanyi.net
! H0 X
! H0 X   Portions of this code were written by Noam Bernstein as part of
! H0 X   his employment for the U.S. Government, and are not subject
! H0 X   to copyright in the USA.
! H0 X
! H0 X
! H0 X   When using this software, please cite the following reference:
! H0 X
! H0 X   http://www.libatoms.org
! H0 X
! H0 X  Additional contributions by
! H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! H0 X
! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

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
  real(dp) :: scale_lattice, scale_pos
  character(len=FIELD_LENGTH) :: properties

  type(Atoms_ll) :: structure_ll
  type(Atoms_ll_entry), pointer :: entry

  call system_initialise(PRINT_NORMAL)

  call initialise(cli_params)
  call param_register(cli_params, 'infile', param_mandatory, infilename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'infile_is_list', 'F', infile_is_list, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'outfile', 'stdout', outfilename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'decimation', '1', decimation, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'min_time', '-1.0', min_time, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'max_time', '-1.0', max_time, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'sort_Time', 'F', sort_Time, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'no_Time_dups', 'F', no_Time_dups, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'quiet', 'F', quiet, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'scale_lattice', '1.0', scale_lattice, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'scale_pos', '1.0', scale_pos, help_string="No help yet.  This source file was $LastChangedBy$")
  properties = ""
  call param_register(cli_params, 'properties', 'species:pos:Z', properties, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli_params, do_check = .true.)) then
      call print_usage()
    call system_abort('could not parse argument line')
  end if
  call finalise(cli_params)

  call print("infile " // trim(infilename) // " infile_is_list " // infile_is_list)
  call print("outfilename " // trim(outfilename))
  call print("decimation " // decimation // " min_time " // min_time // " max_time " // max_time)
  call print("sort_Time " // sort_Time // " no_Time_dups " // no_Time_dups // " quiet " // quiet)
  call print("scale_lattice " // scale_lattice // " scale_pos " // scale_pos)

  no_compute_index=.false.
  if ((.not. sort_Time) .and. (.not. no_Time_dups)) no_compute_index=.true.

  call read_xyz(structure_ll, infilename, infile_is_list, decimation, min_time, max_time, sort_Time, no_Time_dups, quiet, no_compute_index, properties=properties)

  call initialise(outfile, outfilename, action=OUTPUT)
  entry => structure_ll%first
  do while (associated(entry))
    if (scale_lattice /= 1.0_dp) entry%at%lattice = scale_lattice * entry%at%lattice
    if (scale_pos /= 1.0_dp) entry%at%pos = scale_pos * entry%at%pos
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

    call print("Usage: " // trim(my_exec_name)//" infile=filename [infile_is_list=logical(F)]", PRINT_ALWAYS)
    call print("       outfile=filename [decimation=n(1)]", PRINT_ALWAYS)
    call print("       [min_time=t(-1.0)] [max_time=t(-1.0)]", PRINT_ALWAYS)
    call print("       [sort_Time(F)] [no_Time_dups(F)] [quiet(F)]", PRINT_ALWAYS)
    call print("       [scale_lattice=scale(1.0)] [scale_pos=scale(1.0)]", PRINT_ALWAYS)
    call print("       [properties=prop1:prop2...(species:pos:Z,use blank for all props.)]", PRINT_ALWAYS)
  end subroutine print_usage

end program clean_traj
