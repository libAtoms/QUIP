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

!X
!X FilePot Module
!X
!% FilePot is a potential that computes things by writing atomic config to a
!% file, running a command, and reading its output
!%
!% it takes an argument string on initialization with one mandatory parameter
!%>   command=path_to_command
!% and three optional parameters
!%>   property_list=prop1:T1:N1:prop2:T2:N2...
!% which defaults to 'pos', 
!%>    filename=<string>
!% which defaults to 'filepot', and
!%>   min_cutoff=cutoff
!% which default to zero. If min_cutoff is non zero and the cell is narrower
!% than $2*min_cutoff$ in any direction then it will be replicated before
!% being written to the file. The forces are taken from the primitive cell
!% and the energy is reduced by a factor of the number of repeated copies.
!% 
!% The command takes 2 arguments, the names of the input and the output files.
!% command output is in extended xyz form.
!%
!% energy and virial (both optional) are passed via the comment, labeled as
!%     'energy=E' and 'virial="vxx vxy vxz vyx vyy vyz vzx vzy vzz"'.
!%
!%  per atoms data is at least atomic type and optionally 
!%     a local energy (labeled 'local_e:R:1') 
!%  and 
!%  forces (labeled 'force:R:3')
!%
!% right now (14/2/2008) the atoms_xyz reader requires the 1st 3 columns after
!%   the atomic type to be the position.
!% 
!% If you ask for some quantity from FilePot_Calc and it's not in the output file, it
!% returns an error status or crashes (if err isn't present).
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! some day implement calculation queues
! how to do this? one possibility:
!  add args_str to FilePot_Calc, and indicate
!   queued_force=i to indicate this atoms structure should be queued for the calculation
!    of forces on atom i (or a list, or a range?)
! or
!   process_queue, which would process the queue and fill in all the forces
#include "error.inc"
module FilePot_module

use error_module
use system_module, only : dp, print, PRINT_NORMAL, PRINT_ALWAYS, PRINT_VERBOSE, inoutput, current_verbosity, system_command, optional_default, parse_string, operator(//)
use mpi_context_module
use dictionary_module
use paramreader_module
use linearalgebra_module
use connection_module
use atoms_types_module
use atoms_module
use structures_module
use CInOutput_module

implicit none
private


public :: FilePot_type
type FilePot_type
  character(len=STRING_LENGTH) :: command
  character(len=STRING_LENGTH) :: property_list
  character(len=STRING_LENGTH) :: read_extra_property_list
  character(len=STRING_LENGTH) :: read_extra_param_list
  character(len=STRING_LENGTH) :: property_list_prefixes
  character(len=STRING_LENGTH) :: filename
  real(dp)            :: min_cutoff

  character(len=STRING_LENGTH) :: init_args_str
  type(MPI_context) :: mpi

end type FilePot_type

public :: Initialise
interface Initialise
  module procedure FilePot_Initialise
end interface Initialise

public :: Finalise
interface Finalise
  module procedure FilePot_Finalise
end interface Finalise

public :: cutoff
interface cutoff
   module procedure FilePot_cutoff
end interface

public :: Wipe
interface Wipe
  module procedure FilePot_Wipe
end interface Wipe

public :: Print
interface Print
  module procedure FilePot_Print
end interface Print

public :: calc
interface calc
  module procedure FilePot_Calc
end interface

contains


subroutine FilePot_Initialise(this, args_str, mpi, error)
  type(FilePot_type), intent(inout) :: this
  character(len=*), intent(in) :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  type(Dictionary) ::  params
  character(len=STRING_LENGTH) :: command, property_list, read_extra_property_list, &
       read_extra_param_list, property_list_prefixes, filename
  real(dp) :: min_cutoff
  real(dp) :: r_scale, E_scale
  logical :: do_rescale_r, do_rescale_E

  INIT_ERROR(error)

  this%init_args_str = args_str

  call initialise(params)
  call param_register(params, 'command', PARAM_MANDATORY, command, help_string="system command to execute that should read the structure file, run the model and deposit the output file")
  call param_register(params, 'property_list', 'species:pos', property_list, help_string="list of properties to print with the structure file")
  call param_register(params, 'read_extra_property_list', '', read_extra_property_list, help_string="names of extra properties to read from filepot.out files")
  call param_register(params, 'property_list_prefixes', '', property_list_prefixes, help_string="list of prefixes to which run_suffix will be applied during calc()")
  call param_register(params, 'read_extra_param_list', 'QM_cell', read_extra_param_list, help_string="list of extra params (comment line in XYZ) to read from filepot.out files. Default is 'QM_cell'")
  call param_register(params, 'filename', 'filepot', filename, help_string="seed name for directory and structure files to be used")
  call param_register(params, 'min_cutoff', '0.0', min_cutoff, help_string="if the unit cell does not fit into this cutoff, it is periodically replicated so that it does")
  call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
  call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")

  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='filepot_initialise args_str')) then
    RAISE_ERROR("FilePot_initialise failed to parse args_str='"//trim(args_str)//"'", error)
  endif
  call finalise(params)
  if (do_rescale_r .or. do_rescale_E) then
     RAISE_ERROR("FilePot_Initialise: rescaling of potential with r_scale and E_scale not yet implemented!", error)
  end if

  this%command = command
  this%property_list = property_list
  this%read_extra_property_list = read_extra_property_list
  this%read_extra_param_list = read_extra_param_list
  this%property_list_prefixes = property_list_prefixes
  this%min_cutoff = min_cutoff
  this%filename = filename
  if (present(mpi)) this%mpi = mpi

end subroutine FilePot_Initialise

subroutine FilePot_Finalise(this)
  type(FilePot_type), intent(inout) :: this

  call wipe(this)

end subroutine FilePot_Finalise

subroutine FilePot_Wipe(this)
  type(FilePot_type), intent(inout) :: this

  this%command=""
  this%property_list=""
  this%read_extra_property_list=""
  this%read_extra_param_list=""
  this%min_cutoff = 0.0_dp
  this%filename = ""

end subroutine FilePot_Wipe

function FilePot_cutoff(this)
  type(FilePot_type), intent(in) :: this
  real(dp) :: FilePot_cutoff
  FilePot_cutoff = 0.0_dp ! return zero, because FilePot does its own connection calculation
end function FilePot_cutoff

subroutine FilePot_Print(this, file)
  type(FilePot_type),    intent(in)           :: this
  type(Inoutput), intent(inout),optional,target:: file

  if (current_verbosity() < PRINT_NORMAL) return

  call print("FilePot: command='"//trim(this%command)// &
       "' filename='"//trim(this%filename)//&
       "' property_list='"//trim(this%property_list)//&
       "' read_extra_property_list='"//trim(this%read_extra_property_list)//&
       "' read_extra_param_list='"//trim(this%read_extra_param_list)//&
       "' property_list_prefixes='"//trim(this%property_list_prefixes)//&
       "' min_cutoff="//this%min_cutoff,file=file)

end subroutine FilePot_Print

subroutine FilePot_Calc(this, at, energy, local_e, forces, virial, local_virial, args_str, error)
  type(FilePot_type), intent(inout) :: this
  type(Atoms), intent(inout) :: at
  real(dp), intent(out), optional :: energy
  real(dp), intent(out), target, optional :: local_e(:)
  real(dp), intent(out), optional :: forces(:,:), local_virial(:,:)
  real(dp), intent(out), optional :: virial(3,3)
  character(len=*), intent(in), optional :: args_str
  integer, intent(out), optional :: error

  character(len=STRING_LENGTH)  :: xyzfile, outfile, filename, run_suffix
  character(len=STRING_LENGTH) :: my_args_str
  integer :: nx, ny, nz, i
  type(Atoms) :: sup
  integer :: status, n_properties, my_err
  character(len=STRING_LENGTH) :: read_extra_property_list, read_extra_param_list, property_list, tmp_properties_array(100)
  type(Dictionary) :: cli
  logical :: FilePot_log, filename_override
  
  INIT_ERROR(error)

  if (present(energy)) energy = 0.0_dp
  if (present(local_e)) local_e = 0.0_dp
  if (present(forces)) forces = 0.0_dp
  if (present(virial)) virial = 0.0_dp
  if (present(local_virial)) local_virial = 0.0_dp
  my_args_str = ''
  if (present(args_str)) my_args_str = args_str

  call initialise(cli)
  call param_register(cli, "FilePot_log", "T", FilePot_log, help_string="if True, save logfile of all the filepot.xyz and filepot.out")
  call param_register(cli, "read_extra_property_list", trim(this%read_extra_property_list), read_extra_property_list, help_string="extra properties to read from filepot.out. Overrides init_args version.")
  call param_register(cli, "read_extra_param_list", trim(this%read_extra_param_list), read_extra_param_list, help_string="extra params to read from filepot.out. Overrides init_args version.")
  call param_register(cli, "run_suffix", '', run_suffix, help_string="suffix to apply to property names in this%property_list_prefixes")
  call param_register(cli, 'filename', 'filepot', filename, has_value_target=filename_override, help_string="seed name for directory and structure files to be used")

  if (.not. param_read_line(cli, my_args_str, ignore_unknown=.true.,task='filepot_calc args_str')) then
    RAISE_ERROR("FilePot_calc failed to parse args_str='"//trim(args_str)//"'",error)
  endif
  call finalise(cli)

  if(filename_override) then
     this%filename = filename
  end if

  ! Run external command either if MPI object is not active, or if it is active and we're the
  ! master process. Function does not return on any node until external command is finished.

  if (.not. this%mpi%active .or.  (this%mpi%active .and. this%mpi%my_proc == 0)) then
     
     if(this%mpi%active) then
        xyzfile=(trim(this%filename)//"."//this%mpi%my_proc//".xyz")
        outfile=(trim(this%filename)//"."//this%mpi%my_proc//".out")
     else
        xyzfile=(trim(this%filename)//".xyz")
        outfile=(trim(this%filename)//".out")
     end if

     call print("FilePot: filename seed=`"//trim(this%filename)//"'", PRINT_VERBOSE)
     call print("FilePot: outfile=`"//trim(outfile)//"'", PRINT_VERBOSE)
     call print("FilePot: xyzfile=`"//trim(xyzfile)//"'", PRINT_VERBOSE)
     call system_command("rm -f "//trim(outfile), status=status)
     if (status /= 0) call print("WARNING: FilePot_calc failed to delete outfile="//trim(outfile)//" before running filepot command", PRINT_ALWAYS)
     call system_command("rm -f "//trim(xyzfile)//".idx", status=status)
     if (status /= 0) call print("WARNING: FilePot_calc failed to delete index file="//trim(xyzfile)//".idx before running filepot command", PRINT_ALWAYS)
     call system_command("rm -f "//trim(outfile)//".idx", status=status)
     if (status /= 0) call print("WARNING: FilePot_calc failed to delete index file="//trim(outfile)//".idx before running filepot command", PRINT_ALWAYS)

     ! Do we need to replicate cell to exceed min_cutoff ?
     if (this%min_cutoff .fne. 0.0_dp) then
        call fit_box_in_cell(this%min_cutoff, this%min_cutoff, this%min_cutoff, at%lattice, nx, ny, nz)
     else
        nx = 1; ny = 1; nz = 1
     end if

     property_list = this%property_list
     if (len_trim(this%property_list_prefixes) /= 0) then
        call parse_string(this%property_list_prefixes, ':', tmp_properties_array, n_properties, error=error)
        do i=1, n_properties
           property_list = trim(property_list)//':'//trim(tmp_properties_array(i))//run_suffix
        end do
     end if

     if (nx /= 1 .or. ny /= 1 .or. nz /= 1) then
        call Print('FilePot: replicating cell '//nx//'x'//ny//'x'//nz//' times.')
        call supercell(sup, at, nx, ny, nz)
        call write(sup, xyzfile, properties=property_list)
     else
        call write(at, xyzfile, properties=property_list)
     end if

     if (FilePot_log) then
       if (nx /= 1 .or. ny /= 1 .or. nz /= 1) then
          call write(sup, "FilePot_pos_log.xyz", properties=property_list,append=.true.)
        else
          call write(at, "FilePot_pos_log.xyz", properties=property_list,append=.true.)
        endif
     endif

!     call print("FilePot: invoking external command "//trim(this%command)//" "//' '//trim(xyzfile)//" "// &
!          trim(outfile)//" on "//at%N//" atoms...")
     call print("FilePot: invoking external command "//trim(this%command)//' '//trim(xyzfile)//" "// &
          trim(outfile)//" "//trim(my_args_str)//" on "//at%N//" atoms...")

     ! call the external command here
!     call system_command(trim(this%command)//" "//trim(xyzfile)//" "//trim(outfile),status=status)
     call system_command(trim(this%command)//' '//trim(xyzfile)//" "//trim(outfile)//" "//trim(my_args_str),status=status)
     call print("FilePot: got status " // status // " from external command")

     ! read back output from external command
     call filepot_read_output(outfile, at, nx, ny, nz, energy, local_e, forces, virial, local_virial, &
          read_extra_property_list, read_extra_param_list, run_suffix, filepot_log=FilePot_log, error=error)
     PASS_ERROR_WITH_INFO("Filepot_Calc reading output", error)
  end if

  if (this%mpi%active) then
     ! Share results with other nodes
     if (present(energy))  call bcast(this%mpi, energy)
     if (present(local_e)) call bcast(this%mpi, local_e)

     if (present(forces))  call bcast(this%mpi, forces)
     if (present(virial))  call bcast(this%mpi, virial)
     if (present(local_virial))  call bcast(this%mpi, local_virial)

     call bcast(this%mpi, my_err)
  end if

end subroutine FilePot_calc

subroutine filepot_read_output(outfile, at, nx, ny, nz, energy, local_e, forces, virial, local_virial, &
     read_extra_property_list, read_extra_param_list, run_suffix, filepot_log, error)
  character(len=*), intent(in) :: outfile
  type(Atoms), intent(inout) :: at
  integer, intent(in) :: nx, ny, nz
  real(dp), intent(out), optional :: energy
  real(dp), intent(out), target, optional :: local_e(:)
  real(dp), intent(out), optional :: forces(:,:), local_virial(:,:)
  real(dp), intent(out), optional :: virial(3,3)
  character(len=*), intent(in) :: read_extra_property_list, read_extra_param_list, run_suffix
  logical, intent(in), optional :: filepot_log
  integer, intent(out), optional :: error

  character(STRING_LENGTH) :: tmp_params_array(100), copy_keys(100)
  integer :: i, n_params, n_copy
  type(atoms) :: at_out, primitive
  integer, pointer :: Z_p(:)
  real(dp) :: virial_1d(9)
  real(dp), pointer :: local_e_p(:), forces_p(:,:), local_virial_p(:,:)
  logical :: my_filepot_log

  INIT_ERROR(error)

  my_filepot_log = optional_default(.false., filepot_log)

  call print('Filepot: reading back results from file '//trim(outfile))
  call read(at_out, outfile, error=error)
  PASS_ERROR(error)

  if (nx /= 1 .or. ny /= 1 .or. nz /= 1) then
     ! Discard atoms outside the primitive cell
     call select(primitive, at_out, list=(/ (i, i=1,at_out%N/(nx*ny*nz) ) /))
     at_out = primitive
     call finalise(primitive)
  end if

  call print("FilePot: the following keys were found in the atoms structure:")
  call print_keys(at_out%params)

  if (at_out%N /= at%N) then
     RAISE_ERROR("filepot_read_output in '"//trim(outfile)//"' got N="//at_out%N//" /= at%N="//at%N, error)
  endif

  if (.not. assign_pointer(at_out,'Z',Z_p)) then
     RAISE_ERROR("filepot_read_output in '"//trim(outfile)//"' couldn't assign pointer for field Z", error)
  endif
  do i=1, at%N
    if (at%Z(i) /= Z_p(i)) then
      RAISE_ERROR("filepot_read_output in '"//trim(outfile)//"' got Z("//i//")="//at_out%Z(i)//" /= at%Z("//i//")="//at%Z(i), error)
    endif
  end do

  if (present(energy)) then
    if (.not. get_value(at_out%params,'energy',energy)) then
      RAISE_ERROR("filepot_read_output needed energy, but couldn't find energy in '"//trim(outfile)//"'", error)
    endif
    ! If cell was repeated, reduce energy by appropriate factor
    ! to give energy of primitive cell.
    if (nx /= 1 .or. ny /= 1 .or. nz /= 1) &
         energy = energy/(nx*ny*nz)
  endif

  if (present(virial)) then
     if (nx /= 1 .or. ny /= 1 .or. nz /= 1) then
          RAISE_ERROR("filepot_read_output: don't know how to rescale virial for repicated system", error)
     endif

    if ( get_value(at_out%params,'virial',virial_1d)) then
      virial(:,1) = virial_1d(1:3)
      virial(:,2) = virial_1d(4:6)
      virial(:,3) = virial_1d(7:9)
    elseif( .not. get_value(at_out%params,'virial',virial) ) then
      RAISE_ERROR("filepot_read_output needed virial, but couldn't find virial in '"//trim(outfile)//"'", error)
    endif
  endif

  if (present(local_e)) then
    if (.not. assign_pointer(at_out, 'local_e', local_e_p)) then
	RAISE_ERROR("filepot_read_output needed local_e, but couldn't find local_e in '"//trim(outfile)//"'", error)
    endif
    local_e = local_e_p
  endif

  if (present(forces)) then
    if (.not. assign_pointer(at_out, 'force', forces_p)) then
	RAISE_ERROR("filepot_read_output needed forces, but couldn't find force in '"//trim(outfile)//"'", error)
    endif
    forces = forces_p
  endif

  if (present(local_virial)) then
    if (.not. assign_pointer(at_out, 'local_virial', local_virial_p)) then
	RAISE_ERROR("filepot_read_output needed local_virial, but couldn't find local_virial in '"//trim(outfile)//"'", error)
    endif
    local_virial = local_virial_p
  endif

  if (len_trim(read_extra_property_list) > 0) then
     call copy_properties(at, at_out, trim(read_extra_property_list))
  endif
  
  if (len_trim(read_extra_param_list) > 0) then
     call parse_string(read_extra_param_list, ':', tmp_params_array, n_params, error=error)
     PASS_ERROR(error)

     n_copy = 0
     do i=1,n_params
        if (has_key(at_out%params, trim(tmp_params_array(i)))) then
           n_copy = n_copy + 1
           copy_keys(n_copy) = tmp_params_array(i)
           call print("FilePot copying param key "//trim(copy_keys(n_copy)), PRINT_VERBOSE)
        else if  (has_key(at_out%params, trim(tmp_params_array(i))//trim(run_suffix))) then
           n_copy = n_copy + 1
           copy_keys(n_copy) =  trim(tmp_params_array(i))//trim(run_suffix)
           call print("FilePot copying param key "//trim(copy_keys(n_copy)), PRINT_VERBOSE)
        end if
     end do
     
     call subset(at_out%params, copy_keys(1:n_copy), at%params, out_no_initialise=.true.)
  end if

  if (my_filepot_log) then
     call write(at_out, "FilePot_out_log.xyz", append=.true.)
  endif

  call finalise(at_out)

end subroutine filepot_read_output

end module FilePot_module
