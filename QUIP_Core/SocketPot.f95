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
!X SocketPot Module
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"
module SocketPot_module

use error_module
use system_module
use mpi_context_module
use dictionary_module
use paramreader_module
use atoms_types_module
use atoms_module
use SocketTools_module

implicit none
private

public :: SocketPot_type
type SocketPot_type
   character(len=STRING_LENGTH) :: init_args_str
   character(len=STRING_LENGTH) :: ip
   character(len=STRING_LENGTH) :: property_list
   character(len=STRING_LENGTH) :: property_list_prefixes
   character(len=STRING_LENGTH) :: read_extra_property_list
   character(len=STRING_LENGTH) :: read_extra_param_list
   integer :: port, client_id, label, last_label, buffsize
   type(MPI_context) :: mpi
end type SocketPot_type

public :: Initialise
interface Initialise
  module procedure SocketPot_Initialise
end interface Initialise

public :: Finalise
interface Finalise
  module procedure SocketPot_Finalise
end interface Finalise

public :: cutoff
interface cutoff
   module procedure SocketPot_cutoff
end interface

public :: Wipe
interface Wipe
  module procedure SocketPot_Wipe
end interface Wipe

public :: Print
interface Print
  module procedure SocketPot_Print
end interface Print

public :: calc
interface calc
  module procedure SocketPot_Calc
end interface

contains


subroutine SocketPot_Initialise(this, args_str, mpi, error)
  type(SocketPot_type), intent(inout) :: this
  character(len=*), intent(in) :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  type(Dictionary) :: params

  INIT_ERROR(error)

  call finalise(this)

  call initialise(params)
  call param_register(params, 'server_ip', '', this%ip, help_string="IP address to send configs to")
  call param_register(params, 'server_port', '0', this%port, help_string="TCP port number to send configs to. Default 0.")
  call param_register(params, 'client_id', '0', this%client_id, help_string="Identity of this client. Default 0.")
  call param_register(params, 'buffsize', '1000000', this%buffsize, help_string='Size of recv buffer in bytes. Default 100000')
  call param_register(params, 'property_list', 'species:pos', this%property_list, help_string="list of properties to send with the structure")
  call param_register(params, 'read_extra_property_list', '', this%read_extra_property_list, help_string="names of extra properties to read back")
  call param_register(params, 'read_extra_param_list', 'QM_cell', this%read_extra_param_list, help_string="list of extra params (comment line in XYZ) to read back. Default is 'QM_cell'")
  call param_register(params, 'property_list_prefixes', '', this%property_list_prefixes, help_string="list of prefixes to which run_suffix will be applied during calc()")

  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='SocketPot_initialise')) then
       RAISE_ERROR('SocketPot_Initialise failed to parse args_str="'//trim(args_str)//"'", error)
  endif
  call finalise(params)

  this%label = 0
  this%init_args_str = args_str
  if (present(mpi)) this%mpi = mpi

end subroutine SocketPot_Initialise

subroutine SocketPot_Finalise(this)
  type(SocketPot_type), intent(inout) :: this

  call wipe(this)

end subroutine SocketPot_Finalise

subroutine SocketPot_Wipe(this)
  type(SocketPot_type), intent(inout) :: this

  this%ip = ''
  this%port = 0
  this%client_id = 0
  this%buffsize = 1000000
  this%property_list=""
  this%read_extra_property_list=""
  this%read_extra_param_list=""
  this%property_list_prefixes=""

end subroutine SocketPot_Wipe

function SocketPot_cutoff(this)
  type(SocketPot_type), intent(in) :: this
  real(dp) :: SocketPot_cutoff
  SocketPot_cutoff = 0.0_dp ! return zero, because SocketPot does its own connection calculation
end function SocketPot_cutoff

subroutine SocketPot_Print(this, file)
  type(SocketPot_type),    intent(in)           :: this
  type(Inoutput), intent(inout),optional,target:: file

  if (current_verbosity() < PRINT_NORMAL) return

  call print('SocketPot: connecting to QUIP server on host '//trim(this%ip)//':'//this%port//' with client_id '//this%client_id, file=file)
  call print(" buffsize="//this%buffsize//&
             " property_list='"//trim(this%property_list)//&
             "' read_extra_property_list='"//trim(this%read_extra_property_list)//&
             "' read_extra_param_list='"//trim(this%read_extra_param_list)//&
             "' property_list_prefixes='"//trim(this%property_list_prefixes)//"'", file=file)

end subroutine SocketPot_Print

subroutine SocketPot_Calc(this, at, energy, local_e, forces, virial, local_virial, args_str, error)
  type(SocketPot_type), intent(inout) :: this
  type(Atoms), intent(inout) :: at
  real(dp), intent(out), optional :: energy
  real(dp), intent(out), target, optional :: local_e(:)
  real(dp), intent(out), optional :: forces(:,:), local_virial(:,:)
  real(dp), intent(out), optional :: virial(3,3)
  character(len=*), intent(in), optional :: args_str
  integer, intent(out), optional :: error

  type(Atoms) :: at_copy
  type(Dictionary) :: cli
  integer n_params, n_copy, i, n_properties, label
  real(dp), pointer :: local_e_ptr(:), force_ptr(:,:), local_virial_ptr(:,:)
  logical :: calc_energy, calc_local_e, calc_force, calc_virial, calc_local_virial
  character(STRING_LENGTH) :: my_args_str, tmp_params_array(100), copy_keys(100), read_extra_property_list, &
       read_extra_param_list, run_suffix, property_list, tmp_properties_array(100)

  INIT_ERROR(error)

  my_args_str = ''
  if (present(args_str)) my_args_str = args_str

  call param_register(cli, "read_extra_property_list", trim(this%read_extra_property_list), read_extra_property_list, help_string="extra properties to read from filepot.out. Overrides init_args version.")
  call param_register(cli, "read_extra_param_list", trim(this%read_extra_param_list), read_extra_param_list, help_string="extra params to read from filepot.out. Overrides init_args version.")
  call param_register(cli, "run_suffix", '', run_suffix, help_string="suffix to apply to property names in this%property_list_prefixes")

  if (.not. param_read_line(cli, my_args_str, ignore_unknown=.true.,task='filepot_calc args_str')) then
     RAISE_ERROR("FilePot_calc failed to parse args_str='"//trim(args_str)//"'",error)
  endif
  call finalise(cli)

  if (.not. this%mpi%active .or. (this%mpi%active .and. this%mpi%my_proc == 0)) then
     calc_energy = .false.
     calc_local_e = .false.
     calc_force = .false.
     calc_virial = .false.
     calc_local_virial = .false.
     if (present(energy)) then
        energy = 0.0_dp
        calc_energy = .true.
     end if
     if (present(local_e)) then
        local_e = 0.0_dp
        calc_local_e = .true.
     end if
     if (present(forces)) then
        forces = 0.0_dp
        calc_force = .true.
     end if
     if (present(virial)) then
        virial = 0.0_dp
        calc_virial = .true.
     end if
     if (present(local_virial)) then
        local_virial = 0.0_dp
        calc_local_virial = .true.
     end if

     at_copy = at ! work with a copy of original Atoms
     call set_value(at_copy%params, 'calc_energy', calc_energy)
     call set_value(at_copy%params, 'calc_local_e', calc_local_e)
     call set_value(at_copy%params, 'calc_virial', calc_virial)
     call set_value(at_copy%params, 'calc_force', calc_force)
     call set_value(at_copy%params, 'calc_local_virial', calc_local_virial)
     call set_value(at_copy%params, 'label', this%label)
     if (present(args_str)) then
        call set_value(at_copy%params, 'calc_args_str', args_str)
     end if

     property_list = this%property_list
     if (len_trim(this%property_list_prefixes) /= 0) then
        call parse_string(this%property_list_prefixes, ':', tmp_properties_array, n_properties, error=error)
        do i=1, n_properties
           property_list = trim(property_list)//':'//trim(tmp_properties_array(i))//run_suffix
        end do
     end if

     call system_timer('socket_send')
     call socket_send_xyz(this%ip, this%port, this%client_id, at_copy, properties=property_list)
     call system_timer('socket_send')

     this%label = this%label + 1

     call system_timer('socket_recv')
     call socket_recv_xyz(this%ip, this%port, this%client_id, this%buffsize, at_copy)
     call system_timer('socket_recv')

     if (.not. get_value(at_copy%params, 'label', label)) then
        RAISE_ERROR('SocketPot_Calc: missing label param in received atoms', error)
     end if
     if (label /= this%label) then
        RAISE_ERROR('SocketPot_Calc: mismatch between labels expected ('//this%label//') and received ('//label//').', error)
     end if

     if (present(energy)) then
        if (.not. get_value(at_copy%params, 'energy', energy)) then
           call print('WARNING SocketPot_calc: "energy" requested but not returned by callback')
        endif
     end if

     if (present(local_e)) then
        if (.not. assign_pointer(at_copy, 'local_e', local_e_ptr)) then
           call print('WARNING SocketPot_calc: "local_e" requested but not returned by callback')
        else
           local_e(:) = local_e_ptr
        end if
     end if

     if (present(forces)) then
        if (.not. assign_pointer(at_copy, 'force', force_ptr)) then
           call print('WARNING SocketPot_calc: "forces" requested but not returned by callback')
        else
           forces(:,:) = force_ptr
        end if
     end if

     if (present(virial)) then
        if (.not. get_value(at_copy%params, 'virial', virial)) then
           call print('WARNING SocketPot_calc: "virial" requested but not returned by callback')
        end if
     end if

     if (present(local_virial)) then
        if (.not. assign_pointer(at_copy, 'local_virial', local_virial_ptr)) then
           call print('WARNING SocketPot_calc: "local_virial" requested but not returned by callback')
        else
           local_virial(:,:) = local_virial_ptr
        end if
     end if

     if (len_trim(read_extra_property_list) > 0) then
        call copy_properties(at, at_copy, trim(read_extra_property_list))
     endif

     if (len_trim(read_extra_param_list) > 0) then
        call parse_string(read_extra_param_list, ':', tmp_params_array, n_params, error=error)
        PASS_ERROR(error)

        n_copy = 0
        do i=1,n_params
           if (has_key(at_copy%params, trim(tmp_params_array(i)))) then
              n_copy = n_copy + 1
              copy_keys(n_copy) = tmp_params_array(i)
              call print("SocketPot copying param key "//trim(copy_keys(n_copy)), PRINT_VERBOSE)
           else if  (has_key(at_copy%params, trim(tmp_params_array(i))//trim(run_suffix))) then
              n_copy = n_copy + 1
              copy_keys(n_copy) =  trim(tmp_params_array(i))//trim(run_suffix)
              call print("SocketPot copying param key "//trim(copy_keys(n_copy)), PRINT_VERBOSE)
           end if
        end do

        call subset(at_copy%params, copy_keys(1:n_copy), at%params, out_no_initialise=.true.)
     end if
  end if

  if (this%mpi%active) then
     ! Share results with other nodes
     if (present(energy))  call bcast(this%mpi, energy)
     if (present(local_e)) call bcast(this%mpi, local_e)

     if (present(forces))  call bcast(this%mpi, forces)
     if (present(virial))  call bcast(this%mpi, virial)
     if (present(local_virial))  call bcast(this%mpi, local_virial)
  end if

  call finalise(at_copy)

end subroutine SocketPot_calc

end module SocketPot_module
