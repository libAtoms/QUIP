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
!X Callbackpot Module
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"
module Callbackpot_module

use libatoms_module
use mpi_context_module

implicit none
private

integer, parameter :: MAX_CALLBACKS = 200
integer :: n_callbacks = 0

public :: Callbackpot_type
type CallbackPot_type
   character(len=STRING_LENGTH) :: init_args_str
   character(len=STRING_LENGTH) :: label
   integer :: callback_id
   type(MPI_context) :: mpi
end type CallbackPot_type

public :: Initialise
interface Initialise
  module procedure Callbackpot_Initialise
end interface Initialise

public :: Finalise
interface Finalise
  module procedure Callbackpot_Finalise
end interface Finalise

public :: cutoff
interface cutoff
   module procedure Callbackpot_cutoff
end interface

public :: Wipe
interface Wipe
  module procedure Callbackpot_Wipe
end interface Wipe

public :: Print
interface Print
  module procedure Callbackpot_Print
end interface Print

public :: calc
interface calc
  module procedure Callbackpot_Calc
end interface

public :: set_callback
interface set_callback
   module procedure callbackpot_set_callback
end interface


interface 
   subroutine register_callbackpot_sub(sub)
     interface 
        subroutine sub(at)
          integer, intent(in) :: at(12)
        end subroutine sub
     end interface
   end subroutine register_callbackpot_sub
end interface
  
interface 
   subroutine call_callbackpot_sub(i,at)
     integer, intent(in) :: at(12)
     integer, intent(in) :: i
   end subroutine call_callbackpot_sub
end interface


contains


subroutine Callbackpot_Initialise(this, args_str, mpi, error)
  type(Callbackpot_type), intent(inout) :: this
  character(len=*), intent(in) :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  type(Dictionary) :: params

  INIT_ERROR(error)

  call finalise(this)

  call initialise(params)
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='CallbackPot_initialise')) then
       RAISE_ERROR('CallbackPot_Initialise failed to parse args_str="'//trim(args_str)//"'", error)
  endif
  call finalise(params)

  this%init_args_str = args_str
  if (present(mpi)) this%mpi = mpi

end subroutine Callbackpot_Initialise

subroutine callbackpot_set_callback(this, callback, error)
  type(Callbackpot_type), intent(inout) :: this
  interface
     subroutine callback(at)
       integer, intent(in) :: at(12)
     end subroutine callback
  end interface
  integer, intent(out), optional :: error

  INIT_ERROR(error)

  if (n_callbacks >= MAX_CALLBACKS) then
       RAISE_ERROR('CallbackPot_Initialise: Too many registered callback routines', error)
  endif
  this%callback_id = n_callbacks
  n_callbacks = n_callbacks + 1
  call register_callbackpot_sub(callback)

end subroutine Callbackpot_Set_callback

subroutine Callbackpot_Finalise(this)
  type(Callbackpot_type), intent(inout) :: this

  call wipe(this)

end subroutine Callbackpot_Finalise

subroutine Callbackpot_Wipe(this)
  type(Callbackpot_type), intent(inout) :: this

  this%callback_id = -1

end subroutine Callbackpot_Wipe

function Callbackpot_cutoff(this)
  type(Callbackpot_type), intent(in) :: this
  real(dp) :: Callbackpot_cutoff
  Callbackpot_cutoff = 0.0_dp ! return zero, because Callbackpot does its own connection calculation
end function Callbackpot_cutoff

subroutine Callbackpot_Print(this, file)
  type(Callbackpot_type),    intent(in)           :: this
  type(Inoutput), intent(inout),optional,target:: file

  if (current_verbosity() < PRINT_NORMAL) return

  call print("Callbackpot: callback_id="//this%callback_id)
  call print("Callbackpot: label="//this%label)

end subroutine Callbackpot_Print

subroutine Callbackpot_Calc(this, at, energy, local_e, forces, virial, local_virial, args_str, error)
  type atoms_ptr_type
     type(atoms), pointer :: p
  end type atoms_ptr_type
  type(Callbackpot_type), intent(inout) :: this
  type(Atoms), intent(inout) :: at
  real(dp), intent(out), optional :: energy
  real(dp), intent(out), target, optional :: local_e(:)
  real(dp), intent(out), optional :: forces(:,:), local_virial(:,:)
  real(dp), intent(out), optional :: virial(3,3)
  character(len=*), intent(in), optional :: args_str
  integer, intent(out), optional :: error

  type(Atoms), target :: at_copy
  type(atoms_ptr_type) :: at_ptr
  integer :: at_ptr_i(12)
  real(dp), pointer :: local_e_ptr(:), force_ptr(:,:), local_virial_ptr(:,:)
  logical :: calc_energy, calc_local_e, calc_force, calc_virial, calc_local_virial

  INIT_ERROR(error)

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

  if (this%callback_id < 0) then
     RAISE_ERROR('callbackpot_calc: callback_id < 0', error)
  endif

  at_copy = at
  call set_value(at_copy%params, 'calc_energy', calc_energy)
  call set_value(at_copy%params, 'calc_local_e', calc_local_e)
  call set_value(at_copy%params, 'calc_virial', calc_virial)
  call set_value(at_copy%params, 'calc_force', calc_force)
  call set_value(at_copy%params, 'calc_local_virial', calc_local_virial)
  call set_value(at_copy%params, 'callback_id', this%callback_id)
  call set_value(at_copy%params, 'label', this%label)
  if (present(args_str)) then
     call set_value(at_copy%params, 'calc_args_str', args_str)
  end if

  at_ptr%p => at_copy
  at_ptr_i = transfer(at_ptr, at_ptr_i)
  call call_callbackpot_sub(this%callback_id, at_ptr_i)

  if (present(energy)) then
     if (.not. get_value(at_copy%params, 'energy', energy)) then
          call print('WARNING Callbackpot_calc: "energy" requested but not returned by callback')
     endif
  end if

  if (present(local_e)) then
     if (.not. assign_pointer(at_copy, 'local_e', local_e_ptr)) then
        call print('WARNING Callbackpot_calc: "local_e" requested but not returned by callback')
     else
        local_e(:) = local_e_ptr
     end if
  end if

  if (present(forces)) then
     if (.not. assign_pointer(at_copy, 'force', force_ptr)) then
        call print('WARNING Callbackpot_calc: "forces" requested but not returned by callback')
     else
        forces(:,:) = force_ptr
     end if
  end if

  if (present(virial)) then
     if (.not. get_value(at_copy%params, 'virial', virial)) then
          call print('WARNING Callbackpot_calc: "virial" requested but not returned by callback')
     end if
  end if

  if (present(local_virial)) then
     if (.not. assign_pointer(at_copy, 'local_virial', local_virial_ptr)) then
        call print('WARNING Callbackpot_calc: "local_virial" requested but not returned by callback')
     else
        local_virial(:,:) = local_virial_ptr
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

end subroutine Callbackpot_calc

end module Callbackpot_module
