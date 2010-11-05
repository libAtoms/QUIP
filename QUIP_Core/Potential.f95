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
!% Potential module to provide interface to all energy/force/virial 
!% calculating schemes, including actual calculations, as well as
!% abstract hybrid schemes such as LOTF, Force Mixing, ONIOM,
!% and the local energy scheme.
!%
!% It can also be used to geometry optimise an Atoms structure,
!% using the 'Minim' interface.
!%
!% Potential_LOTF also has to live in this module so that it can call
!% Potential_minim() without any circular depenendancies.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#ifdef HAVE_LOCAL_E_MIX
#define HAVE_HYBRID
#endif
#ifdef HAVE_ONIOM
#define HAVE_HYBRID
#endif

#include "error.inc"
module Potential_module

  use libAtoms_module
  use QUIP_Common_module
  use MPI_context_module
  use Potential_simple_module
#ifdef HAVE_TB
  use TB_module
#endif
  use adjustablepotential_module, only: adjustable_potential_init, adjustable_potential_optimise, &
       adjustable_potential_force, adjustable_potential_finalise

  implicit none

#ifndef POTENTIAL_NO_DEFAULT_PRIVATE
  private
#endif

  !*************************************************************************
  !*
  !*  Potential header
  !*
  !*************************************************************************

  public :: Potential
  type Potential
     type(MPI_context) :: mpi
     character(len=FIELD_LENGTH) :: init_args_pot1, init_args_pot2, xml_label, xml_init_args

     logical :: is_simple = .false.
     type(Potential_simple) :: simple

     type(Potential), pointer :: l_mpot1 => null() ! local copies of potentials, if they are created by potential_initialise from an args string
     type(Potential), pointer :: l_mpot2 => null() ! local copies of potentials, if they are created by potential_initialise from an args string

     logical :: is_sum = .false.
     type(Potential_Sum), pointer :: sum => null()

     logical :: is_forcemixing = .false.
     type(Potential_FM), pointer :: forcemixing => null()

     logical :: is_evb = .false.
     type(Potential_EVB), pointer :: evb => null()

#ifdef HAVE_LOCAL_E_MIX
     logical :: is_local_e_mix = .false.
     type(Potential_Local_E_Mix), pointer :: local_e_mix => null()
#endif

#ifdef HAVE_ONIOM
     logical :: is_oniom = .false.
     type(Potential_ONIOM), pointer :: oniom => null()
#endif
  end type Potential

  public :: Initialise, Potential_Filename_Initialise
  interface Initialise
     module procedure Potential_Initialise, Potential_Initialise_inoutput
  end interface

  public :: Finalise
  interface Finalise
     module procedure Potential_Finalise
  end interface

  public :: Setup_Parallel
  interface Setup_Parallel
     module procedure Potential_Setup_Parallel
  end interface

  public :: Print
  interface Print
     module procedure Potential_Print
  end interface

  public :: Cutoff
  interface Cutoff
     module procedure Potential_Cutoff
  end interface

  public :: Calc
  interface Calc
     module procedure Potential_Calc
  end interface

  public :: Minim
  interface Minim
     module procedure Potential_Minim
  end interface

  public :: test_gradient
  interface test_gradient
     module procedure pot_test_gradient
  end interface

  public :: n_test_gradient
  interface n_test_gradient
     module procedure pot_n_test_gradient
  end interface

  public :: set_callback
  interface set_callback
     module procedure potential_set_callback
  end interface

  public :: potential_minimise
  type Potential_minimise
     real(dp) :: minim_pos_lat_preconditioner = 1.0_dp

     real(dp) :: minim_save_lat(3,3)
     character(len=1024) :: minim_args_str
     logical :: minim_do_pos, minim_do_lat
     integer :: minim_n_eval_e, minim_n_eval_f, minim_n_eval_ef
     real(dp) :: pos_lat_preconditioner_factor
     type(Atoms), pointer :: minim_at => null()
     real(dp), pointer :: last_connect_x(:) => null()
     type(Potential), pointer :: minim_pot => null()
     type(CInoutput), pointer :: minim_cinoutput_movie => null()
     logical, dimension(3,3) :: lattice_fix = .false.
     real(dp), dimension(3,3) :: external_pressure = 0.0_dp
  end type Potential_minimise

  type(Potential), pointer :: parse_pot
  logical, save :: parse_in_pot, parse_in_pot_done, parse_matched_label

  public :: run
  interface run
     module procedure dynamicalsystem_run
  end interface run

#include "Potential_Sum_header.f95"
#include "Potential_ForceMixing_header.f95"
#include "Potential_EVB_header.f95"

#ifdef HAVE_LOCAL_E_MIX
#include "Potential_Local_E_Mix_header.f95"
#endif
#ifdef HAVE_ONIOM
#include "Potential_ONIOM_header.f95"
#endif

  contains

  !*************************************************************************
  !*
  !*  Potential routines
  !*
  !*************************************************************************

recursive subroutine potential_Filename_Initialise(this, args_str, param_filename, bulk_scale, mpi_obj, error)
  type(Potential), intent(inout) :: this
  character(len=*), intent(in) :: args_str !% Valid arguments are 'Sum', 'ForceMixing', 'EVB', 'Local_E_Mix' and 'ONIOM', and any type of simple_potential
  character(len=*), intent(in) :: param_filename !% name of xml parameter file for potential initializers
  type(Atoms), optional, intent(inout) :: bulk_scale !% optional bulk structure for calculating space and E rescaling
  type(MPI_Context), intent(in), optional :: mpi_obj
  integer, intent(out), optional :: error

  type(inoutput) :: io

  INIT_ERROR(error)

  call initialise(io, param_filename, INPUT, master_only=.true.)
  call initialise(this, args_str, io, bulk_scale, mpi_obj, error=error)
  PASS_ERROR(error)
  call finalise(io)

end subroutine potential_Filename_Initialise

subroutine potential_initialise_inoutput(this, args_str, io_obj, bulk_scale, mpi_obj, error)
  type(Potential), intent(inout) :: this
  character(len=*), intent(in) :: args_str !% Valid arguments are 'Sum', 'ForceMixing', 'EVB', 'Local_E_Mix' and 'ONIOM', and any type of simple_potential
  type(InOutput), intent(in) :: io_obj !% name of xml parameter inoutput for potential initializers
  type(Atoms), optional, intent(inout) :: bulk_scale !% optional bulk structure for calculating space and E rescaling
  type(MPI_Context), intent(in), optional :: mpi_obj
  integer, intent(out), optional :: error

  type(extendable_str) :: es

  INIT_ERROR(error)

  call initialise(es)
  if (present(mpi_obj)) then
    call read(es, io_obj%unit, convert_to_string=.true., mpi_comm=mpi_obj%communicator)
  else
    call read(es, io_obj%unit, convert_to_string=.true.)
  endif

  call initialise(this, args_str, param_str=string(es), bulk_scale=bulk_scale, mpi_obj=mpi_obj, error=error)
  PASS_ERROR(error)
  call finalise(es)

end subroutine potential_initialise_inoutput

recursive subroutine potential_initialise(this, args_str, pot1, pot2, param_str, bulk_scale, mpi_obj, error)
  type(Potential), intent(inout) :: this
  character(len=*), intent(in) :: args_str !% Valid arguments are 'Sum', 'ForceMixing', 'EVB', 'Local_E_Mix' and 'ONIOM', and any type of simple_potential
  type(Potential), optional, intent(in), target :: pot1 !% Optional first Potential upon which this Potential is based
  type(Potential), optional, intent(in), target :: pot2 !% Optional second potential
  character(len=*), optional, intent(in) :: param_str !% contents of xml parameter file for potential initializers, if needed
  type(Atoms), optional, intent(inout) :: bulk_scale !% optional bulk structure for calculating space and E rescaling
  type(MPI_Context), optional, intent(in) :: mpi_obj
  integer, intent(out), optional :: error

  type(Potential), pointer :: u_pot1, u_pot2
  type(Dictionary) :: params

  logical :: has_xml_label
  character(len=FIELD_LENGTH) :: my_args_str

  INIT_ERROR(error)

  call finalise(this)

  call initialise(params)
  call param_register(params, 'xml_label', '', this%xml_label, has_value_target=has_xml_label)
  if(.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Initialise args_str')) then
    RAISE_ERROR("Potential_initialise failed to parse args_str='"//trim(args_str)//"'", error)
  endif
  call finalise(params)

  if(has_xml_label) then
     call Potential_read_params_xml(this, param_str)
     my_args_str = trim(this%xml_init_args)
  else
     my_args_str = trim(args_str)
  endif

  call initialise(params)
  call param_register(params, 'init_args_pot1', '', this%init_args_pot1)
  call param_register(params, 'init_args_pot2', '', this%init_args_pot2)
  call param_register(params, 'Sum', 'false', this%is_sum)
  call param_register(params, 'ForceMixing', 'false', this%is_forcemixing)
  call param_register(params, 'EVB', 'false', this%is_evb)
#ifdef HAVE_LOCAL_E_MIX
  call param_register(params, 'Local_E_Mix', 'false', this%is_local_e_mix)
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
  call param_register(params, 'ONIOM', 'false', this%is_oniom)
#endif /* HAVE_ONIOM */

  if(.not. param_read_line(params, my_args_str, ignore_unknown=.true.,task='Potential_Initialise args_str')) then
    RAISE_ERROR("Potential_initialise failed to parse args_str='"//trim(my_args_str)//"'", error)
  endif
  call finalise(params)

  if (present(pot1) .and. len_trim(this%init_args_pot1) > 0) then
    RAISE_ERROR("Potential_initialise got both pot and args_str with init_args_pot1 passed in, conflict", error)
   endif
  if (present(pot2) .and. len_trim(this%init_args_pot2) > 0) then
    RAISE_ERROR("Potential_initialise got both pot2 and args_str with init_args_pot2 passed in, conflict", error)
  endif

  if (len_trim(this%init_args_pot1) > 0) then
    allocate(this%l_mpot1)
    call initialise(this%l_mpot1, args_str=this%init_args_pot1, param_str=param_str, bulk_scale=bulk_scale, mpi_obj=mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing pot1", error)
    u_pot1 => this%l_mpot1
  else
    u_pot1 => pot1
  endif
  if (len_trim(this%init_args_pot2) > 0) then
    allocate(this%l_mpot2)
    call initialise(this%l_mpot2, args_str=this%init_args_pot2, param_str=param_str, bulk_scale=bulk_scale, mpi_obj=mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing pot2", error)
    u_pot2 => this%l_mpot2
  else
    if (present(pot2)) then
      u_pot2 => pot2
    else
      nullify(u_pot2)
    endif
  endif

  if (present(mpi_obj)) this%mpi = mpi_obj

  this%is_simple = .not. any( (/ this%is_sum, this%is_forcemixing, this%is_evb &
#ifdef HAVE_LOCAL_E_MIX
  , this%is_local_e_mix &
#endif
#ifdef HAVE_ONIOM
  , this%is_oniom  &
#endif
  /) )

  if (count( (/this%is_simple, this%is_sum, this%is_forcemixing, this%is_evb &
#ifdef HAVE_LOCAL_E_MIX
          , this%is_local_e_mix &
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
          , this%is_oniom &
#endif /* HAVE_ONIOM */
          /) ) /= 1) then
     RAISE_ERROR("Potential_initialise found too few or two many Potential types args_str='"//trim(my_args_str)//"'", error)
  end if

  if (this%is_simple) then
    if (present(bulk_scale)) call print("Potential_initialise Simple ignoring bulk_scale passed in", PRINT_ALWAYS)

    call initialise(this%simple, my_args_str, param_str, mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing pot", error)
  else if (this%is_sum) then
    if (present(bulk_scale)) call print("Potential_initialise Sum ignoring bulk_scale passed in", PRINT_ALWAYS)

    if (.not. associated(u_pot2)) then
      RAISE_ERROR('Potential_initialise: two potentials needs for sum potential', error)
    endif

    allocate(this%sum)
    call initialise(this%sum, my_args_str, u_pot1, u_pot2, mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing sum", error)

  else if (this%is_forcemixing) then

    allocate(this%forcemixing)
    if(associated(u_pot2)) then
       call initialise(this%forcemixing, my_args_str, u_pot1, u_pot2, bulk_scale, mpi_obj, error=error)
    else
       ! if only one pot is given, assign it to QM. this is useful for time-embedding LOTF
       call initialise(this%forcemixing, my_args_str, qmpot=u_pot1, reference_bulk=bulk_scale, mpi=mpi_obj, error=error)
    endif
    PASS_ERROR_WITH_INFO("Initializing forcemixing", error)

  else if (this%is_evb) then
    if (present(bulk_scale)) call print("Potential_initialise EVB ignoring bulk_scale passed in", PRINT_ALWAYS)

    !1 is enough, everything is the same apart from the passed topology
    if (associated(u_pot2)) then
      call print('WARNING! Potential_initialise EVB ignoring u_pot2 (it will use u_pot1 twice)', PRINT_ALWAYS)
    endif

    allocate(this%evb)
    call initialise(this%evb, my_args_str, u_pot1, mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing evb", error)

#ifdef HAVE_LOCAL_E_MIX
  else if (this%is_local_e_mix) then
    if (.not. associated(u_pot2)) then
      RAISE_ERROR('Potential_initialise: two potentials needs for local_e_mix potential', error)
   endif

    allocate(this%local_e_mix)
    call initialise(this%local_e_mix, my_args_str, u_pot1, u_pot2, bulk_scale, mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing local_e_mix", error)
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
  else if (this%is_oniom) then
    if (.not. associated(u_pot2)) then
      RAISE_ERROR('Potential_initialise: two potentials needs for oniom potential', error)
    endif

    allocate(this%oniom)
    call initialise(this%oniom, my_args_str, u_pot1, u_pot2, bulk_scale, mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing oniom", error)
#endif /* HAVE_ONIOM */
  end if

  end subroutine potential_initialise

  recursive subroutine potential_finalise(this, error)
    type(Potential), intent(inout) :: this
    integer, intent(out), optional :: error

  INIT_ERROR(error)

    if (this%is_simple) then
       call finalise(this%simple)
    else if (this%is_sum) then
       call finalise(this%sum)
       deallocate(this%sum)
    else if (this%is_forcemixing) then
       call finalise(this%forcemixing)
       deallocate(this%forcemixing)
    else if (this%is_evb) then
       call finalise(this%evb)
       deallocate(this%evb)
#ifdef HAVE_LOCAL_E_MIX
    else if (this%is_local_e_mix) then
       call finalise(this%local_e_mix)
       deallocate(this%local_e_mix)
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
    else if (this%is_oniom) then
       call finalise(this%oniom)
       deallocate(this%oniom)
#endif /* HAVE_ONIOM */
    end if

    if (associated(this%l_mpot1)) call finalise(this%l_mpot1)
    if (associated(this%l_mpot2)) call finalise(this%l_mpot2)
    nullify(this%l_mpot1)
    nullify(this%l_mpot2)

    this%is_simple = .false.
    this%is_sum = .false.
    this%is_forcemixing = .false.
    this%is_evb = .false.
#ifdef HAVE_LOCAL_E_MIX
    this%is_local_e_mix = .false.
#endif
#ifdef HAVE_ONIOM
    this%is_oniom = .false.
#endif

  end subroutine potential_finalise

  recursive subroutine potential_calc(this, at, energy, force, virial, local_energy, local_virial, args_str, error)
    type(Potential), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    real(dp), intent(out), optional :: energy
    real(dp), intent(out), optional :: force(:,:)
    real(dp), intent(out), optional :: virial(3,3)
    real(dp), intent(out), optional :: local_virial(:,:)
    real(dp), intent(out), optional :: local_energy(:)
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    real(dp), pointer :: at_force_ptr(:,:), at_local_energy_ptr(:), at_local_virial_ptr(:,:)
    type(Dictionary) :: params
    character(len=STRING_LENGTH) :: calc_energy, calc_force, calc_virial, calc_local_energy, calc_local_virial, extra_args_str
    character(len=STRING_LENGTH) :: use_calc_energy, use_calc_force, use_calc_virial, use_calc_local_energy, use_calc_local_virial

    integer i
    real(dp) :: effective_cutoff

    INIT_ERROR(error)

    if (cutoff(this) > 0.0_dp) then
       ! For Potentials which need connectivity information, ensure Atoms cutoff is >= Potential cutoff

       if (.not. at%connect%initialised) then
          call print('Potential_calc: setting Atoms cutoff to Potential cutoff ('//cutoff(this)//')', PRINT_VERBOSE)
          call set_cutoff(at, cutoff(this))
          call calc_connect(at)
       end if

       if (at%use_uniform_cutoff) then
          if (at%cutoff < cutoff(this)) then
             RAISE_ERROR('IP_calc: cutoff of Atoms object ('//at%cutoff//') < Potential cutoff ('//cutoff(this)//')', error)
          end if
       else
          ! Find effective cutoff 
          effective_cutoff = 0.0_dp
          do i = 1, at%N
             if (at%Z(i) < 1 .or. at%Z(i) > 116) then
                RAISE_ERROR('IP_calc: bad atomic number i='//i//' z='//at%z(i), error)
             end if
             if (ElementCovRad(at%Z(i)) > effective_cutoff) effective_cutoff = ElementCovRad(at%Z(i))
          end do
          effective_cutoff = (2.0_dp * effective_cutoff) * at%cutoff
          if (effective_cutoff < cutoff(this)) then
             RAISE_ERROR('IP_calc: effective cutoff of Atoms object ('//effective_cutoff//') < Potential cutoff ('//cutoff(this)//')', error)
          end if
       end if
    end if

    call initialise(params)
    call param_register(params, "energy", "", calc_energy)
    call param_register(params, "virial", "", calc_virial)
    call param_register(params, "force", "", calc_force)
    call param_register(params, "local_energy", "", calc_local_energy)
    call param_register(params, "local_virial", "", calc_local_virial)
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Calc args_str')) then
       RAISE_ERROR('Potential_Calc failed to parse args_str="'//trim(args_str)//'"', error)
    endif
    call finalise(params)

    extra_args_str = ""

    ! create property/param names and possibly storage
    use_calc_energy = trim(calc_energy)
    if (present(energy) .and. len_trim(calc_energy) == 0) then ! have optional and no args_str - make new name for param
       use_calc_energy = "energy"
       do while (lookup_entry_i(at%params, trim(calc_energy)) > 0)
	 use_calc_energy = "T"//trim(use_calc_energy)
       end do
       extra_args_str = trim(extra_args_str) // " energy="//trim(use_calc_energy)
    endif

    use_calc_force = trim(calc_force)

    if (present(force)) then
       if (len_trim(calc_force) == 0) then ! have force optional but not args_str - add property from pointer in new property name
	  use_calc_force = "force"
	  do while (has_property(at, trim(use_calc_force)))
	    use_calc_force = "T"//trim(use_calc_force)
	  end do
	  extra_args_str = trim(extra_args_str) // " force="//trim(use_calc_force)
	  call add_property_from_pointer(at, trim(use_calc_force), force)
       else ! has force optional _and_ args_str - add property with its own storage
	  call add_property(at, trim(calc_force), 0.0_dp, n_cols=3, ptr2=at_force_ptr)
       endif
    else if (len_trim(calc_force) > 0) then ! no optional, have args_str - add property with its own storage
       call add_property(at, trim(calc_force), 0.0_dp, n_cols=3, ptr2=at_force_ptr)
    endif

    use_calc_virial = trim(calc_virial) ! have optional and no args_str - make new name for param
    if (present(virial) .and. len_trim(calc_virial) == 0) then
       use_calc_virial = "virial"
       do while (lookup_entry_i(at%params, trim(calc_virial)) > 0)
	 use_calc_virial = "T"//trim(use_calc_virial)
       end do
       extra_args_str = trim(extra_args_str) // " virial="//trim(use_calc_virial)
    endif

    use_calc_local_energy = trim(calc_local_energy)
    if (present(local_energy)) then
       if (len_trim(calc_local_energy) == 0) then ! have local_energy optional but not args_str - add property from pointer in new property name
	  use_calc_local_energy = "local_energy"
	  do while (has_property(at, trim(use_calc_local_energy)))
	    use_calc_local_energy = "T"//trim(use_calc_local_energy)
	  end do
	  extra_args_str = trim(extra_args_str) // " local_energy="//trim(use_calc_local_energy)
	  call add_property_from_pointer(at, trim(use_calc_local_energy), local_energy)
       else ! has local_energy optional _and_ args_str - add property with its own storage
	  call add_property(at, trim(calc_local_energy), 0.0_dp, n_cols=1, ptr=at_local_energy_ptr)
       endif
    else if (len_trim(calc_local_energy) > 0) then ! no optional, have args_str - add property with its own storage
       call add_property(at, trim(calc_local_energy), 0.0_dp, n_cols=1, ptr=at_local_energy_ptr)
    endif

    use_calc_local_virial = trim(calc_local_virial)
    if (present(local_virial)) then
       if (len_trim(calc_local_virial) == 0) then ! have local_virial optional but not args_str - add property from pointer in new property name
	  use_calc_local_virial = "local_virial"
	  do while (has_property(at, trim(use_calc_local_virial)))
	    use_calc_local_virial = "T"//trim(use_calc_local_virial)
	  end do
	  extra_args_str = trim(extra_args_str) // " local_virial="//trim(use_calc_local_virial)
	  call add_property_from_pointer(at, trim(use_calc_local_virial), local_virial)
       else ! has local_virial optional _and_ args_str - add property with its own storage
	  call add_property(at, trim(calc_local_virial), 0.0_dp, n_cols=9, ptr2=at_local_virial_ptr)
       endif
    else if (len_trim(calc_local_virial) > 0) then ! no optional, have args_str - add property with its own storage
       call add_property(at, trim(calc_local_virial), 0.0_dp, n_cols=9, ptr2=at_local_virial_ptr)
    endif

    ! do actual calculation using args_str
    if (this%is_simple) then
       call Calc(this%simple, at, trim(args_str)//" "//trim(extra_args_str), error=error)
       PASS_ERROR(error)
    else if (this%is_sum) then
       call Calc(this%sum, at, trim(args_str)//" "//trim(extra_args_str), error=error)
       PASS_ERROR(error)
    else if (this%is_forcemixing) then
       call Calc(this%forcemixing, at, trim(args_str)//" "//trim(extra_args_str), error=error)
       PASS_ERROR(error)
    else if (this%is_evb) then
       call Calc(this%evb, at, trim(args_str)//" "//trim(extra_args_str), error=error)
       PASS_ERROR(error)
#ifdef HAVE_LOCAL_E_MIX
    else if (this%is_local_e_mix) then
      call Calc(this%local_e_mix, at, trim(args_str)//" "//trim(extra_args_str), error=error)
      PASS_ERROR(error)
#endif /* HAVE_local_e_mix */
#ifdef HAVE_ONIOM
    else if (this%is_oniom) then
      call Calc(this%oniom, at, trim(args_str)//" "//trim(extra_args_str), error=error)
      PASS_ERROR(error)
#endif /* HAVE_ONIOM */
    else
      RAISE_ERROR('Potential_Calc: no potential type is set',error)
    endif

    ! get values from property/param if necessary, and remove temporary if only optional was passed in
    if (present(energy)) then
      call get_param_value(at, trim(use_calc_energy), energy)
      if (len_trim(calc_energy) == 0) call remove_value(at%params, trim(use_calc_energy))
    endif
    if (present(force)) then
      if (len_trim(calc_force) == 0) then ! property should be a pointer to force optional
	call remove_property(at, trim(use_calc_force))
      else
        force = at_force_ptr
      endif
    endif
    if (present(virial)) then
      call get_param_value(at, trim(use_calc_virial), virial)
      if (len_trim(calc_virial) == 0) call remove_value(at%params, trim(use_calc_virial))
    endif
    if (present(local_energy)) then
      if (len_trim(calc_local_energy) == 0) then ! property should be pointer to local_energy optional
	call remove_property(at, trim(use_calc_local_energy))
      else
	local_energy = at_local_energy_ptr
      endif
    endif
    if (present(local_virial)) then
      if (len_trim(calc_local_virial) == 0) then ! property should be pointer to local_virial optional
	call remove_property(at, trim(use_calc_local_virial))
      else
	local_virial = at_local_virial_ptr
      endif
    endif

  end subroutine potential_calc

  subroutine Potential_setup_parallel(this, at, args_str, error)
    type(Potential), intent(inout) :: this
    type(Atoms), intent(inout) :: at     !% The atoms structure to compute energy and forces
    character(len=*), intent(in) :: args_str
    integer, intent(out), optional :: error

    INIT_ERROR(error)

    if(this%is_simple) then
       call setup_parallel(this%simple, at, args_str)
    endif
  end subroutine Potential_setup_parallel

  recursive subroutine potential_print(this, file, error)
    type(Potential), intent(inout) :: this
    type(Inoutput), intent(inout), optional :: file
    integer, intent(out), optional :: error

    INIT_ERROR(error)

    if (this%is_simple) then
       call Print('Potential containing potential')
       call Print(this%simple, file=file)
    else if (this%is_sum) then
       call Print(this%sum, file=file)
    else if (this%is_forcemixing) then
       call Print(this%forcemixing, file=file)
    else if (this%is_evb) then
       call Print(this%evb, file=file)
#ifdef HAVE_LOCAL_E_MIX
    else if (this%is_local_e_mix) then
       call Print(this%local_e_mix, file=file)
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
    else if (this%is_oniom) then
       call Print(this%oniom, file=file)
#endif /* HAVE_oniom */
    else
       RAISE_ERROR('Potential_Print: no potential type is set', error)
    end if
  end subroutine potential_print


  recursive function potential_cutoff(this, error)
    type(Potential), intent(in) :: this
    integer, intent(out), optional :: error
    real(dp) :: potential_cutoff

    INIT_ERROR(error)

    if (this%is_simple) then
       potential_cutoff = cutoff(this%simple)
    else if (this%is_sum) then
       potential_cutoff = cutoff(this%sum)
    else if (this%is_forcemixing) then
       potential_cutoff = cutoff(this%forcemixing)
    else if (this%is_evb) then
       potential_cutoff = cutoff(this%evb)
#ifdef HAVE_LOCAL_E_MIX
    else if (this%is_local_e_mix) then
       potential_cutoff = cutoff(this%local_e_mix)
#endif /* HAVE_local_e_mix */
#ifdef HAVE_ONIOM
    else if (this%is_oniom) then
       potential_cutoff = cutoff(this%oniom)
#endif /* HAVE_ONIOM */
    else
       RAISE_ERROR('Potential_Cutoff: no potential type if set', error)
    end if
  end function potential_cutoff



  !*************************************************************************
  !*
  !*  Minimisation routines
  !*
  !*************************************************************************


  function potential_minim(this, at, method, convergence_tol, max_steps, linminroutine, do_print, print_inoutput, print_cinoutput, &
       do_pos, do_lat, args_str, eps_guess, use_n_minim, use_fire, lattice_fix, external_pressure, use_precond, hook_print_interval, error)
    type(Atoms), intent(inout), target :: at !% starting configuration
    type(Potential), intent(inout), target :: this !% potential to evaluate energy/forces with
    character(*), intent(in)    :: method !% passed to minim()
    real(dp),     intent(in)    :: convergence_tol !% Minimisation is treated as converged once $|\mathbf{\nabla}f|^2 <$
                                                    !% 'convergence_tol'. 
    integer,      intent(in)    :: max_steps  !% Maximum number of steps
    character(*), intent(in),optional    :: linminroutine !% Name of the line minisation routine to use, passed to base minim()
    logical, optional :: do_print !% if true, print configurations using minim's hook()
    type(inoutput), intent(inout), optional, target :: print_inoutput !% inoutput object to print configs to, needed if do_print is true
    type(cinoutput), intent(inout), optional, target :: print_cinoutput !% cinoutput object to print configs to, needed if do_print is true
    logical, optional :: do_pos, do_lat !% do relaxation w.r.t. positions and/or lattice (is neither is included, do both)
    character(len=*), intent(in), optional :: args_str !% arguments to pass to calc()
    real(dp), intent(in), optional :: eps_guess !% eps_guess argument to pass to minim
    logical, intent(in), optional :: use_n_minim !% if true, use n_minim instead of minim
    logical, intent(in), optional :: use_fire   !% if true, use fire_minim instead of minim
    logical, dimension(3,3), optional :: lattice_fix !% Mask to fix some components of lattice. Defaults to all false.
    real(dp), dimension(3,3), optional :: external_pressure
    logical, intent(in), optional :: use_precond
    integer, intent(in), optional :: hook_print_interval !% how often to print xyz from hook function
    integer, intent(out), optional :: error !% set to 1 if an error occurred during minimisation
    integer::potential_minim

    logical :: my_use_precond
    integer n_iter, n_iter_tot
    real(dp), allocatable :: x(:)
    real(dp) :: deform_grad(3,3)
    logical my_do_print
    logical done
    real(dp) :: my_eps_guess
    logical my_use_n_minim, my_use_fire

    real(dp) :: initial_E, final_E, mass
    type(potential_minimise) am
    integer :: am_data_size
    character, allocatable :: am_data(:)
    character, dimension(1) :: am_mold
    integer :: status

    INIT_ERROR(error)

    my_use_precond = optional_default(.false., use_precond)

    am_data_size = size(transfer(am, am_mold))
    allocate(am_data(am_data_size))

    my_use_n_minim = optional_default(.false., use_n_minim)
    my_use_fire = optional_default(.false., use_fire)

    my_eps_guess = optional_default(1.0e-2_dp/at%N, eps_guess)
    if (my_eps_guess .feq. 0.0_dp) my_eps_guess = 1.0e-2_dp/at%N

    call calc_connect(at)
    am%minim_at => at
    am%pos_lat_preconditioner_factor = am%minim_pos_lat_preconditioner*am%minim_at%N

    if (present(args_str)) then
      am%minim_args_str = args_str
    else
      am%minim_args_str = ""
    endif
    am%minim_pot => this

    if (.not.present(do_pos) .and. .not. present(do_lat)) then
      am%minim_do_pos = .true.
      am%minim_do_lat = .true.
    else
      am%minim_do_pos = optional_default(.false., do_pos)
      am%minim_do_lat = optional_default(.false., do_lat)
    endif

    am%lattice_fix = .false.
    if (present(lattice_fix)) am%lattice_fix = lattice_fix

    am%external_pressure = 0.0_dp
    if (present(external_pressure)) then
       am%external_pressure = external_pressure
       if( (am%external_pressure(1,1) .fne. am%external_pressure(2,2)) .or. &
       & (am%external_pressure(1,1) .fne. am%external_pressure(3,3)) .or. &
       & (am%external_pressure(1,2) .fne. 0.0_dp) .or. &
       & (am%external_pressure(1,3) .fne. 0.0_dp) .or. &
       & (am%external_pressure(2,1) .fne. 0.0_dp) .or. &
       & (am%external_pressure(2,3) .fne. 0.0_dp) .or. &
       & (am%external_pressure(3,1) .fne. 0.0_dp) .or. &
       & (am%external_pressure(3,2) .fne. 0.0_dp) ) then
          if( .not. my_use_fire ) then
             call print_warning('Anisotrpic pressure is being used. Switching to fire_minim.')
             my_use_fire = .true.
          endif
       endif
    endif

    my_do_print = optional_default(.false., do_print)

    if (my_do_print .and. .not. present(print_inoutput) .and. .not. present(print_cinoutput)) &
         call system_abort("potential_minim: do_print is true, but no print_inoutput or print_cinoutput present")

    if (my_do_print) then
       if (present(print_cinoutput)) then
          am%minim_cinoutput_movie => print_cinoutput
          if (this%is_forcemixing) &
               this%forcemixing%minim_cinoutput_movie => print_cinoutput
#ifdef HAVE_LOCAL_E_MIX
          if (this%is_local_e_mix) &
               this%local_e_mix%minim_cinoutput_movie => print_cinoutput
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
          if (this%is_oniom) &
               this%oniom%minim_cinoutput_movie => print_cinoutput
#endif /* HAVE_ONIOM */
       end if
    else
      nullify(am%minim_cinoutput_movie)
    endif

    am%minim_n_eval_e = 0
    am%minim_n_eval_f = 0

    allocate(x(9+am%minim_at%N*3))
    allocate(am%last_connect_x(size(x)))
    am%last_connect_x=1.0e38_dp
    deform_grad = 0.0_dp; call add_identity(deform_grad)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor)

    am_data = transfer(am, am_data)
    if (my_use_n_minim) then
       n_iter = n_minim(x, both_func, my_use_precond, apply_precond_func, initial_E, final_E, my_eps_guess, max_steps, convergence_tol, print_hook, &
            hook_print_interval=hook_print_interval, data=am_data, error=error)
       PASS_ERROR(error)
    else if (my_use_fire) then
       if (has_property(at, 'mass')) then
          mass = at%mass(1)
       else
          mass = ElementMass(at%Z(1))
       end if
       n_iter = fire_minim(x, mass, dummy_energy_func, gradient_func, 1.0_dp, convergence_tol, max_steps, &
            print_hook, hook_print_interval=hook_print_interval, data=am_data, status=status)
    else
       n_iter = minim(x, energy_func, gradient_func, method, convergence_tol, max_steps, linminroutine, &
            print_hook, hook_print_interval=hook_print_interval, eps_guess=my_eps_guess, data=am_data, status=status)
    endif
    call print("minim relax w.r.t. both n_iter " // n_iter, PRINT_VERBOSE)
    call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)
    call calc_connect(am%minim_at)
    n_iter_tot = n_iter
    done = .true.
    deallocate(am%last_connect_x)
    deallocate(x)
    call print("MINIM_N_EVAL E " // am%minim_n_eval_e // " F " // am%minim_n_eval_f // &
      " EF " // am%minim_n_eval_ef, PRINT_VERBOSE)

    potential_minim = n_iter_tot

    deallocate(am_data)

  end function potential_minim

  subroutine undo_travel(at)
    type(Atoms), intent(inout) :: at

!    if (any(at%travel /= 0)) then
!      at%pos = at%pos + (at%lattice .mult. at%travel)
!      at%travel = 0
!    endif
end subroutine undo_travel

  ! test the gradient given a potential, string to pass to calc, and atoms structure
  function pot_test_gradient(pot, at, do_pos, do_lat, args_str, dir_field)
    type(Atoms), intent(inout), target :: at
    type(Potential), target, intent(inout) :: pot
    logical, intent(in), optional :: do_pos, do_lat
    character(len=*), intent(in), optional :: args_str
    character(len=*), intent(in), optional :: dir_field
    logical pot_test_gradient

    real(dp) :: dg(3,3)
    real(dp), allocatable :: x(:), dir_x(:)
    real(dp), pointer :: dir_p(:,:)
    integer :: am_data_size
    type(potential_minimise) :: am
    character, allocatable :: am_data(:)
    character, dimension(1) :: am_mold

    am_data_size = size(transfer(am, am_mold))
    allocate(am_data(am_data_size))

    pot_test_gradient = .true.

    if (present(args_str)) then
      am%minim_args_str = args_str
    else
      am%minim_args_str = ""
    endif
    am%minim_pot => pot

    ! neither do_pos not do_lat is specified, so do both by default
    if (.not. present(do_pos) .and. .not. present(do_lat)) then
      am%minim_do_pos= .true.
      am%minim_do_lat= .true.
    else ! something is specified, so do exactly that
      am%minim_do_pos = optional_default(.false., do_pos)
      am%minim_do_lat = optional_default(.false., do_lat)
    endif

    am%pos_lat_preconditioner_factor = am%minim_pos_lat_preconditioner*at%N

    if (present(dir_field)) then
      if (.not. assign_pointer(at, trim(dir_field), dir_p)) &
	call system_abort("pot_test_gradient got dir_field '"//trim(dir_field)//"' but can't assign pointer to it")
      call print("pot_test_gradient using '"//trim(dir_field)//"' field for gradient test direction")
    else
      nullify(dir_p)
    endif

    am%minim_at => at
    allocate(x(at%N*3+9))
    allocate(am%last_connect_x(at%N*3+9))
    am%last_connect_x=1.0e38_dp
    dg = 0.0_dp; call add_identity(dg)
    call pack_pos_dg(at%pos, dg, x, am%pos_lat_preconditioner_factor)
    if (associated(dir_p)) then
      allocate(dir_x(at%N*3+9))
      dg = 0.0_dp
      call pack_pos_dg(dir_p, dg, dir_x, am%pos_lat_preconditioner_factor)
    endif
    am_data = transfer(am, am_data)
    if (associated(dir_p)) then
      pot_test_gradient = test_gradient(x, energy_func, gradient_func, dir=dir_x, data=am_data)
    else
      pot_test_gradient = test_gradient(x, energy_func, gradient_func, data=am_data)
    endif
    deallocate(x)
    if (allocated(dir_x)) deallocate(dir_x)
    deallocate(am%last_connect_x)

    deallocate(am_data)

  end function pot_test_gradient

  ! test the gradient given a potential, string to pass to calc, and atoms structure
  subroutine pot_n_test_gradient(pot, at, do_pos, do_lat, args_str, dir_field)
    type(Potential), target, intent(inout) :: pot
    type(Atoms), intent(inout), target :: at
    logical, intent(in), optional :: do_pos, do_lat
    character(len=*), intent(in), optional :: args_str
    character(len=*), intent(in), optional :: dir_field

    real(dp), allocatable :: x(:), dir_x(:)
    real(dp), pointer :: dir_p(:,:)
    real(dp) :: dg(3,3)
    integer :: am_data_size
    type(potential_minimise) :: am
    character, allocatable :: am_data(:)
    character, dimension(1) :: am_mold

    am_data_size = size(transfer(am, am_mold))
    allocate(am_data(am_data_size))

    if (present(args_str)) then
      am%minim_args_str = args_str
    else
      am%minim_args_str = ""
    endif
    am%minim_pot => pot

    ! neither do_pos not do_lat is specified, so do both by default
    if (.not. present(do_pos) .and. .not. present(do_lat)) then
      am%minim_do_pos= .true.
      am%minim_do_lat= .true.
    else ! something is specified, so do exactly that
      am%minim_do_pos = optional_default(.false., do_pos)
      am%minim_do_lat = optional_default(.false., do_lat)
    endif

    am%pos_lat_preconditioner_factor = am%minim_pos_lat_preconditioner*at%N

    if (present(dir_field)) then
      if (.not. assign_pointer(at, trim(dir_field), dir_p)) &
	call system_abort("pot_test_gradient got dir_field '"//trim(dir_field)//"' but can't assign pointer to it")
      call print("pot_test_gradient using '"//trim(dir_field)//"' field for gradient test direction")
    else
      nullify(dir_p)
    endif

    am%minim_at => at
    allocate(x(at%N*3+9))
    allocate(am%last_connect_x(at%N*3+9))
    am%last_connect_x=1.0e38_dp
    dg = 0.0_dp; call add_identity(dg)
    call pack_pos_dg(at%pos, dg, x, am%pos_lat_preconditioner_factor)
    if (associated(dir_p)) then
      allocate(dir_x(at%N*3+9))
      dg = 0.0_dp
      call pack_pos_dg(dir_p, dg, dir_x, am%pos_lat_preconditioner_factor)
    endif
    am_data = transfer(am, am_data)
    if (associated(dir_p)) then
      call n_test_gradient(x, energy_func, gradient_func, data=am_data, dir=dir_x)
    else
      call n_test_gradient(x, energy_func, gradient_func, data=am_data)
    endif
    deallocate(x)
    if (allocated(dir_x)) deallocate(dir_x)
    deallocate(am%last_connect_x)

    deallocate(am_data)

  end subroutine pot_n_test_gradient

  ! hook for minim to print configuration during position relaxation
  subroutine print_hook(x,dx,E,done,do_print,am_data)
    real(dp), intent(in) :: x(:), dx(:)
    real(dp), intent(in) :: E
    logical, intent(out) :: done
    logical, optional, intent(in) :: do_print
    character, intent(in), optional :: am_data(:)

    real(dp), allocatable :: f(:,:)
    real(dp), pointer :: f_p(:,:)
    integer, pointer, dimension(:) :: hybrid_mark
    real(dp) :: virial(3,3), deform_grad(3,3), deform_grad_inv(3,3)
    type(potential_minimise) :: am
    logical :: my_do_print

    if (.not. present(am_data)) call system_abort("potential_minimise print_hook must have am_data")
    my_do_print = optional_default(.true., do_print)

    am = transfer(am_data, am)
    ! beware of transfer and pointers !!!
    call atoms_repoint(am%minim_at)

    if (associated(am%minim_cinoutput_movie)) then
      if (size(x) /= am%minim_at%N*3+9) call system_abort("Called gradient_func() with size mismatch " // &
        size(x) // " /= " // am%minim_at%N // "*3+9")

      call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
      call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

      allocate(f(3,am%minim_at%N))
      f = reshape(dx(10:), (/3,am%minim_at%N/))
      f = transpose(deform_grad) .mult. f

      virial = reshape(dx(1:9), (/3,3/) )
      call inverse(deform_grad, deform_grad_inv)
      virial = virial .mult. transpose(deform_grad_inv)

      ! Set atom parameters to print energy, maxforce and maxvirial
      call set_value(am%minim_at%params, 'E', E)
      if (am%minim_do_pos) then
         call set_value(am%minim_at%params, 'MaxForce', maxval(norm(f,1)))
         call set_value(am%minim_at%params, 'df2', norm2(reshape(f,(/3*am%minim_at%N/))))


         if (am%minim_pot%is_forcemixing) then
            if (assign_pointer(am%minim_at, 'hybrid_mark', hybrid_mark)) then 
               if (.not. am%minim_pot%forcemixing%minimise_mm) then
                  call set_value(am%minim_at%params, 'QM_MaxForce', &
                       maxval(abs(f(:,find(hybrid_mark == HYBRID_ACTIVE_MARK)))))
            
                  call set_value(am%minim_at%params, 'QM_df2', &
                       norm2(reshape(f(:,find(hybrid_mark == HYBRID_ACTIVE_MARK)),&
                       (/3*count(hybrid_mark == HYBRID_ACTIVE_MARK)/))))

               end if
            end if
         end if

      end if

      if (am%minim_do_lat) &
	call set_value(am%minim_at%params, 'MaxVirial', maxval(abs(virial)))

      if (assign_pointer(am%minim_at, "forces", f_p)) f_p = -f
      if (assign_pointer(am%minim_at, "force", f_p)) f_p = -f ! compatibility with crack program

      if (my_do_print) then
         if (associated(am%minim_cinoutput_movie)) then
            if (.not. am%minim_pot%mpi%active .or. (am%minim_pot%mpi%active .and. am%minim_pot%mpi%my_proc == 0)) &
                 call write(am%minim_cinoutput_movie, am%minim_at)
         end if
      end if

      deallocate(f)
      call fix_atoms_deform_grad(deform_grad, am%minim_at,am)
    endif
    done = .false.

  end subroutine print_hook

  function dummy_energy_func(xx, am_data)
    real(dp) :: xx(:)
    character(len=1), optional :: am_data(:)
    real(dp) :: dummy_energy_func

    dummy_energy_func = 0.0_dp

  end function dummy_energy_func


  ! compute energy
  function energy_func(x, am_data)
    real(dp) :: x(:)
    character(len=1), optional :: am_data(:)
    real(dp) :: energy_func

    real(dp) :: max_atom_rij_change
    real(dp) :: deform_grad(3,3)
    type(potential_minimise)  :: am

    call system_timer("energy_func")

    if (.not. present(am_data)) call system_abort("potential_minimise energy_func must have am_data")
    am = transfer(am_data, am)
    call atoms_repoint(am%minim_at)

    am%minim_n_eval_e = am%minim_n_eval_e + 1

    if (size(x) /= am%minim_at%N*3+9) call system_abort("Called energy_func() with size mismatch " // &
      size(x) // " /= " // am%minim_at%N // "*3+9")

    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    max_atom_rij_change = max_rij_change(am%last_connect_x, x, cutoff(am%minim_pot), &
      1.0_dp/am%pos_lat_preconditioner_factor)

    call print("energy_func got x " // x, PRINT_NERD)

    call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

    ! Safety factor of 1.1, just in case
    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    if (1.1*max_atom_rij_change >= am%minim_at%cutoff - cutoff(am%minim_pot)) then
      call print("gradient_func: Do calc_connect, atoms moved " // max_atom_rij_change // &
        "*1.1 >= buffer " // (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_connect(am%minim_at)
      am%last_connect_x = x
    else
      call print("gradient_func: Do calc_dists, atoms moved " // max_atom_rij_change // &
        " *1.1 < buffer " // (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_dists(am%minim_at)
    end if

    call print("energy_func using am%minim_at", PRINT_NERD)
    if (current_verbosity() >= PRINT_NERD) call write(am%minim_at, 'stdout')

    call calc(am%minim_pot, am%minim_at, energy = energy_func, args_str = am%minim_args_str)
    call print ("energy_func got energy " // energy_func, PRINT_NERD)

    energy_func = energy_func + cell_volume(am%minim_at)*trace(am%external_pressure) / 3.0_dp
    call print ("energy_func got enthalpy " // energy_func, PRINT_NERD)

    call fix_atoms_deform_grad(deform_grad, am%minim_at, am)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor)

    am_data = transfer(am, am_data)

    call system_timer("energy_func")

  end function energy_func

  ! compute gradient (forces and virial)
  ! x is vectorized version of atomic positions
  ! result is vectorized version of forces
  function gradient_func(x, am_data)
    real(dp) :: x(:)
    character(len=1), optional :: am_data(:)
    real(dp) :: gradient_func(size(x))

    real(dp) :: deform_grad(3,3), virial(3,3), deform_grad_inv(3,3)
    real(dp), allocatable :: f(:,:)
    real(dp) :: max_atom_rij_change
    integer :: t_i(2), i
    integer, pointer, dimension(:) :: move_mask, fixed_pot

    type(potential_minimise) :: am

    call system_timer("gradient_func")

    if (.not. present(am_data)) call system_abort("potential_minimise gradient_func must have am_data")
    am = transfer(am_data, am)
    ! beware of transfer and pointers !!!
    call atoms_repoint(am%minim_at)

    am%minim_n_eval_e = am%minim_n_eval_e + 1

    am%minim_n_eval_f = am%minim_n_eval_f + 1

    if (size(x) /= am%minim_at%N*3+9) call system_abort("Called gradient_func() with size mismatch " // &
      size(x) // " /= " // am%minim_at%N // "*3+9")

    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    max_atom_rij_change = max_rij_change(am%last_connect_x, x, cutoff(am%minim_pot), &
      1.0_dp/am%pos_lat_preconditioner_factor)

    if (current_verbosity() >= PRINT_NERD) call print("gradient_func got x " // x, PRINT_NERD)

    call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

    ! Safety factor of 1.1, just in case
    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    if (1.1*max_atom_rij_change >= am%minim_at%cutoff - cutoff(am%minim_pot)) then
      call print("gradient_func: Do calc_connect, atoms moved " // max_atom_rij_change // &
        "*1.1 >= buffer " // (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_connect(am%minim_at)
      am%last_connect_x = x
    else
      call print("gradient_func: Do calc_dists, atoms moved " // max_atom_rij_change // &
        " *1.1 < buffer " // (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_dists(am%minim_at)
    end if

    if (current_verbosity() >= PRINT_NERD) call add_property(am%minim_at, "force", 0.0_dp, 3)

    allocate(f(3,am%minim_at%N))
    f = 0.0_dp
    virial = 0.0_dp
    if (am%minim_do_pos .and. am%minim_do_lat) then
      call calc(am%minim_pot, am%minim_at, force = f, virial = virial, args_str = am%minim_args_str)
    else if (am%minim_do_pos) then
      call calc(am%minim_pot, am%minim_at, force = f, args_str = am%minim_args_str)
    else
      call calc(am%minim_pot, am%minim_at, virial = virial, args_str = am%minim_args_str)
    endif

    call print("gradient_func used am%minim_at, got forces", PRINT_NERD)
    if (current_verbosity() >= PRINT_NERD) call write(am%minim_at,'stdout')

    ! zero forces if fixed by potential
    if (am%minim_do_pos .and. assign_pointer(am%minim_at, "fixed_pot", fixed_pot)) then
      do i=1, am%minim_at%N
	if (fixed_pot(i) /= 0) f(:,i) = 0.0_dp
      end do
    endif

    ! Zero force on any fixed atoms
    if (assign_pointer(am%minim_at, 'move_mask', move_mask)) then
       do i=1,am%minim_at%N
          if (move_mask(i) == 0) f(:,i) = 0.0_dp
       end do
    end if

    call print ("gradient_func got f", PRINT_NERD)
    call print(f, PRINT_NERD)
    call print ("gradient_func got virial", PRINT_NERD)
    call print(virial, PRINT_NERD)

    virial = virial - am%external_pressure*cell_volume(am%minim_at)
    call print ("gradient_func got virial, external pressure subtracted", PRINT_NERD)
    call print(virial, PRINT_NERD)

    f = transpose(deform_grad) .mult. f

    call inverse(deform_grad, deform_grad_inv)
    virial = virial .mult. transpose(deform_grad_inv)

    ! Zero out components corresponding to fixed lattice elements
    virial = merge(0.0_dp, virial, am%lattice_fix)

    call pack_pos_dg(-f, -virial, gradient_func, 1.0_dp/am%pos_lat_preconditioner_factor)
    if (current_verbosity() >= PRINT_NERD) then
      call print ("gradient_func packed as", PRINT_NERD)
      call print (gradient_func, PRINT_NERD)
    endif

    t_i = maxloc(abs(f))
    call print ("gradient_func got max force component " // f(t_i(1),t_i(2)) // " on atom " // t_i(2) // &
      " max virial component " // maxval(abs(virial)))

    deallocate(f)

    call fix_atoms_deform_grad(deform_grad, am%minim_at,am)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor)
    am_data = transfer(am, am_data)

    call system_timer("gradient_func")

  end function gradient_func

! Cristoph Ortner's simple Hessian preconditioner, dense matrix inverse right now
  subroutine apply_precond_func(x,g,P_g,am_data,error)
    real(dp) :: x(:), g(:), P_g(:)
    character(len=1), optional :: am_data(:)
    integer, optional :: error

    type(potential_minimise) :: am
    real(dp) :: deform_grad(3,3)
    real(dp), allocatable :: P(:,:)
    real(dp) :: c0, c1
    integer :: i, jj, j
    real(dp) :: max_atom_rij_change
    integer :: n_nearest

    INIT_ERROR(error)

    am = transfer(am_data, am)
    call atoms_repoint(am%minim_at)

    max_atom_rij_change = max_rij_change(am%last_connect_x, x, cutoff(am%minim_pot), &
      1.0_dp/am%pos_lat_preconditioner_factor)
max_atom_rij_change = 1.038_dp

    call print("apply_precond_func got x " // x, PRINT_NERD)

    call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

    ! Safety factor of 1.1, just in case
    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    if (1.1*max_atom_rij_change >= am%minim_at%cutoff - cutoff(am%minim_pot)) then
      call print("apply_precond_func: Do calc_connect, atoms moved " // max_atom_rij_change // "*1.1 >= buffer " // &
        (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_connect(am%minim_at)
      am%last_connect_x = x
    else
      call print("apply_precond_func: Do calc_dists, atoms moved " // max_atom_rij_change // " *1.1 < buffer " // &
        (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_dists(am%minim_at)
    end if

    ! identity
    ! c0 = 1.0_dp
    ! c1 = 0.0_dp
    ! simplest expression of connectivity
    c0 = 0.01_dp
    c1 = 1.0_dp

    allocate(P(3*am%minim_at%N,3*am%minim_at%N))
    P = 0.0_dp
    do i=1, am%minim_at%N
      n_nearest = 0
      do jj=1, atoms_n_neighbours(am%minim_at, i)
	if (is_nearest_neighbour(am%minim_at, i, jj)) then
	  n_nearest = n_nearest + 1
	  j = atoms_neighbour(am%minim_at, i, jj)
	  P(3*(i-1)+1,3*(j-1)+1) = -c1
	  P(3*(i-1)+2,3*(j-1)+2) = -c1
	  P(3*(i-1)+3,3*(j-1)+3) = -c1
	endif
      end do
      P(3*(i-1)+1,3*(i-1)+1) = c0+c1*n_nearest
      P(3*(i-1)+2,3*(i-1)+2) = c0+c1*n_nearest
      P(3*(i-1)+3,3*(i-1)+3) = c0+c1*n_nearest
    end do

    P_g(1:9) = g(1:9)
    call symmetric_linear_solve(P, g(10:10+am%minim_at%N*3-1), P_g(10:10+am%minim_at%N*3-1))

  end subroutine apply_precond_func

  ! compute energy and gradient (forces and virial)
  ! x is vectorized version of atomic positions
  ! result is vectorized version of forces
  subroutine both_func(x,val,grad,am_data,error)
    real(dp) :: x(:)
    real(dp) :: val
    real(dp) :: grad(:)
    character(len=1), optional :: am_data(:)
    integer, optional :: error

    real(dp) :: deform_grad(3,3), virial(3,3), deform_grad_inv(3,3)
    real(dp), allocatable :: f(:,:)
    real(dp) :: max_atom_rij_change
    integer :: t_i(2), i
    integer, pointer, dimension(:) :: move_mask, fixed_pot

    type(potential_minimise) :: am

    INIT_ERROR(error)

    call system_timer("both_func")

    if (.not. present(am_data)) call system_abort("potential_minimise both_func must have am_data")
    am = transfer(am_data, am)
    ! beware of transfer and pointers !!!
    call atoms_repoint(am%minim_at)

    am%minim_n_eval_ef = am%minim_n_eval_ef + 1

    if (size(x) /= am%minim_at%N*3+9) call system_abort("Called both_func() with size mismatch " // &
      size(x) // " /= " // am%minim_at%N // "*3+9")

    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    max_atom_rij_change = max_rij_change(am%last_connect_x, x, cutoff(am%minim_pot), &
      1.0_dp/am%pos_lat_preconditioner_factor)
max_atom_rij_change = 1.038_dp

    call print("both_func got x " // x, PRINT_NERD)

    call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

    ! Safety factor of 1.1, just in case
    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    if (1.1*max_atom_rij_change >= am%minim_at%cutoff - cutoff(am%minim_pot)) then
      call print("both_func: Do calc_connect, atoms moved " // max_atom_rij_change // "*1.1 >= buffer " // &
        (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_connect(am%minim_at)
      am%last_connect_x = x
    else
      call print("both_func: Do calc_dists, atoms moved " // max_atom_rij_change // " *1.1 < buffer " // &
        (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_dists(am%minim_at)
    end if

    call print("both_func using am%minim_at", PRINT_NERD)
    if (current_verbosity() >= PRINT_NERD) call write(am%minim_at,'stdout')

    allocate(f(3,am%minim_at%N))
    f = 0.0_dp
    virial = 0.0_dp
    if (am%minim_do_pos .and. am%minim_do_lat) then
      call calc(am%minim_pot, am%minim_at, energy = val, force = f, virial = virial, args_str = am%minim_args_str)
    else if (am%minim_do_pos) then
      call calc(am%minim_pot, am%minim_at, energy = val, force = f, args_str = am%minim_args_str)
    else
      call calc(am%minim_pot, am%minim_at, energy = val, virial = virial, args_str = am%minim_args_str)
    endif

    call print ("both_func got f", PRINT_NERD)
    call print(f, PRINT_NERD)
    call print ("both_func got virial", PRINT_NERD)
    call print(virial, PRINT_NERD)

    virial = virial - am%external_pressure*cell_volume(am%minim_at)
    call print ("both_func got virial, external pressure subtracted", PRINT_NERD)
    call print(virial, PRINT_NERD)

    val = val + cell_volume(am%minim_at)*trace(am%external_pressure) / 3.0_dp
    call print ("both_func got enthalpy " // val, PRINT_NERD)
    ! zero forces if fixed by potential
    if (am%minim_do_pos .and. assign_pointer(am%minim_at, "fixed_pot", fixed_pot)) then
      do i=1, am%minim_at%N
	if (fixed_pot(i) /= 0) f(:,i) = 0.0_dp
      end do
    endif

    ! Zero force on any fixed atoms
    if (assign_pointer(am%minim_at, 'move_mask', move_mask)) then
       do i=1,am%minim_at%N
          if (move_mask(i) == 0) f(:,i) = 0.0_dp
       end do
    end if

    f = transpose(deform_grad) .mult. f

    call inverse(deform_grad, deform_grad_inv)
    virial = virial .mult. transpose(deform_grad_inv)

    call pack_pos_dg(-f, -virial, grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call print ("both_func gradient packed as", PRINT_NERD)
    call print (grad, PRINT_NERD)

    t_i = maxloc(abs(f))
    call print ("both_func got max force component " // f(t_i(1),t_i(2)) // " on atom " // t_i(2) // &
      " max virial component " // maxval(abs(virial)))

    deallocate(f)

    call fix_atoms_deform_grad(deform_grad, am%minim_at,am)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor)

    am_data = transfer(am, am_data)

    call system_timer("both_func")

  end subroutine both_func

  function max_rij_change(last_connect_x, x, r_cut, lat_factor)
    real(dp), intent(inout) :: last_connect_x(:)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: r_cut, lat_factor
    real(dp) :: max_rij_change

    real(dp) :: F0(3,3), F0_evals(3), dF(3,3), dF_evals(3)

    call print("max_rij_change, size(last_connect_x) " // size(last_connect_x) // &
      " size(x) " // size(x), PRINT_NERD)
    if (size(last_connect_x) == size(x)) then
      F0 = lat_factor*reshape(last_connect_x(1:9), (/ 3,3 /))
      dF = lat_factor*reshape(x(1:9), (/ 3,3 /)) - F0
      call diagonalise(matmul(F0,transpose(F0)), F0_evals)
      call diagonalise(matmul(dF,transpose(dF)), dF_evals)
      F0_evals = sqrt(F0_evals)
      dF_evals = sqrt(dF_evals)
      call print("max_rij_change = " // maxval(abs(F0_evals)) //"*2.0_dp*"// &
        maxval(abs(last_connect_x-x))*sqrt(3.0_dp)  // &
        " + " // maxval(abs(dF_evals))//"*"//r_cut, PRINT_VERBOSE)
      max_rij_change = maxval(abs(F0_evals))*2.0_dp*maxval(abs(last_connect_x-x))*sqrt(3.0_dp) + &
                 maxval(abs(dF_evals))*r_cut
    else
      call system_abort("max_rij_change: size(last_connect_x)="//size(last_connect_x)// &
        " /= size(x)="//size(x))
    endif

  end function max_rij_change

  ! prepare atoms structure given a deformation gradient matrix
  subroutine prep_atoms_deform_grad(deform_grad, at, am)
    real(dp), intent(in) :: deform_grad(3,3)
    type(Atoms), intent(inout) :: at
    type(potential_minimise), intent(inout) :: am

    am%minim_save_lat = at%lattice
    call set_lattice(at, deform_grad .mult. am%minim_save_lat, scale_positions=.true.)
  end subroutine prep_atoms_deform_grad

  ! remove effect of deformation gradient applied by prep_atoms_deform_grad
  subroutine fix_atoms_deform_grad(deform_grad, at,am)
    real(dp), intent(in) :: deform_grad(3,3)
    type(Atoms), intent(inout) :: at
    type(potential_minimise), intent(inout) :: am

    real(dp) :: deform_grad_inv(3,3)

    call inverse(deform_grad, deform_grad_inv)
    at%pos = deform_grad_inv .mult. at%pos
    call set_lattice(at, am%minim_save_lat, scale_positions=.false.)

  end subroutine fix_atoms_deform_grad

  ! unpack a deformation grad and a pos array from 1-D array into a atoms%pos 2-D array
  subroutine unpack_pos_dg(xx, at_N, at_pos, dg, lat_factor)
    real(dp), intent(in) :: xx(:)
    integer :: at_N
    real(dp), intent(inout) :: at_pos(:,:)
    real(dp) :: dg(3,3)
    real(dp), intent(in) :: lat_factor

    if (3*at_N+9 /= size(xx)) call system_abort("Called unpack_pos with mismatching sizes x " // size(xx) // " at " // at_N)

    dg = lat_factor*reshape(xx(1:9), (/ 3,3 /) )
    at_pos = reshape(xx(10:), (/ 3, at_N /) )

end subroutine unpack_pos_dg

  ! pack a 3xN 2-D array into a 1-D array
  subroutine pack_pos_dg(x2d, dg2d, x, lat_factor)
    real(dp), intent(in) :: x2d(:,:)
    real(dp), intent(in)  :: dg2d(3,3)
    real(dp), intent(out) :: x(:)
    real(dp), intent(in)  :: lat_factor

    if (size(x2d,1) /= 3) call system_abort("Called pack with mismatching size(x2d,1) " // &
      size(x2d,1) // " != 3")
    
    if (size(dg2d,1) /= 3 .or. size(dg2d,2) /= 3) &
      call system_abort("Called pack with mismatching size(dg2d,1) " // shape(dg2d) // " != 3x3")
    if (size(x) /= (size(x2d)+size(dg2d))) &
      call system_abort("Called pack with mismatching size x " // size(x) // " x2d " // size(x2d) // &
      " dg2d " // size(dg2d))

    x(1:9) = lat_factor*reshape(dg2d, (/ 9 /) )

      x(10:) = reshape(x2d, (/ size(x2d) /) )

end subroutine pack_pos_dg

  subroutine Potential_read_params_xml(this, param_str)
     type(Potential), intent(inout), target :: this
     character(len=*), intent(in) :: param_str

     type(xml_t) :: fxml

     if (len(trim(param_str)) <= 0) &
     call system_abort('Potential_read_params_xml: invalid param_str length '//len(trim(param_str)) )

     parse_in_pot = .false.
     parse_in_pot_done = .false.
     parse_matched_label = .false.
     parse_pot => this

     call open_xml_string(fxml, param_str)

     call parse(fxml,  &
       startElement_handler = Potential_startElement_handler, &
       endElement_handler = Potential_endElement_handler)
     call close_xml_t(fxml)

     if(.not. parse_in_pot_done) &
     call system_abort('Potential_read_params_xml: could not initialise potential from xml_label.')

  endsubroutine Potential_read_params_xml

  subroutine Potential_startElement_handler(URI, localname, name, attributes)
     character(len=*), intent(in)   :: URI
     character(len=*), intent(in)   :: localname
     character(len=*), intent(in)   :: name
     type(dictionary_t), intent(in) :: attributes
   
     integer :: status
     character(len=FIELD_LENGTH) :: value

     if(name == 'Potential') then ! new Potential stanza

        if(parse_in_pot) &
           call system_abort("Potential_startElement_handler entered GAP_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

        if(parse_matched_label) return ! we already found an exact match for this label

        call QUIP_FoX_get_value(attributes, 'label', value, status)
        if(status /= 0) value = ''

        if(len(trim(parse_pot%xml_label)) > 0) then ! we were passed in a label
           if(trim(value) == trim(parse_pot%xml_label)) then ! exact match
              parse_matched_label = .true.
              parse_in_pot = .true.
           else ! no match
              parse_in_pot = .false.
           endif
        else ! no label passed in
           call system_abort("Potential_startElement_handler: no label passed in")
        endif

        call QUIP_FoX_get_value(attributes, 'init_args', value, status)
        if(status == 0) then
           read (value, '(a)') parse_pot%xml_init_args
        else
           call system_abort("Potential_startElement_handler: no init_args attribute found")
        endif
     endif

  endsubroutine Potential_startElement_handler

  subroutine Potential_endElement_handler(URI, localname, name)
     character(len=*), intent(in)   :: URI
     character(len=*), intent(in)   :: localname
     character(len=*), intent(in)   :: name

     if(parse_in_pot) then
        if(name == 'Potential') then
           parse_in_pot = .false.
           parse_in_pot_done = .true.
        endif
     endif

  endsubroutine Potential_endElement_handler

  subroutine potential_set_callback(this, callback)
    type(Potential), intent(inout) :: this
    interface
       subroutine callback(at)
         integer, intent(in) :: at(12)
       end subroutine callback
    end interface
    
    if (this%is_simple) then
       call set_callback(this%simple, callback)
    else
       call system_abort('potential_set_callback() only implemented for simple Potentials.')
    end if

  end subroutine potential_set_callback


#include "Potential_Sum_routines.f95"
#include "Potential_ForceMixing_routines.f95"
#include "Potential_EVB_routines.f95"

#ifdef HAVE_LOCAL_E_MIX
#include "Potential_Local_E_Mix_routines.f95"
#endif
#ifdef HAVE_ONIOM
#include "Potential_ONIOM_routines.f95"
#endif

#include "Potential_Hybrid_utils.f95"

  subroutine DynamicalSystem_run(this, pot, dt, n_steps, hook, hook_interval, write_interval, connect_interval, trajectory, args_str, error)
    type atoms_ptr_type
       type(atoms), pointer :: p
    end type atoms_ptr_type
    type(DynamicalSystem), intent(inout), target :: this
    type(Potential), intent(inout) :: pot
    real(dp), intent(in) :: dt
    integer, intent(in) :: n_steps
    integer, intent(in), optional :: hook_interval, write_interval, connect_interval
    type(CInOutput), intent(inout), optional :: trajectory
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: error
    interface
       subroutine hook()
       end subroutine hook
    end interface    
    
    integer :: n, my_hook_interval, my_write_interval, my_connect_interval
    real(dp) :: e
    real(dp), pointer, dimension(:,:) :: f
    character(len=1024) :: my_args_str
    type(Dictionary) :: params

    INIT_ERROR(error)

    my_hook_interval = optional_default(1, hook_interval)
    my_write_interval = optional_default(1, write_interval)
    my_connect_interval = optional_default(1, connect_interval)
    my_args_str = optional_default("", args_str)
    call initialise(params)
    call read_string(params, my_args_str)
    if (.not. has_key(params, 'energy')) call set_value(params, 'energy')
    if (.not. has_key(params, 'force'))  call set_value(params, 'force')
    my_args_str = write_string(params)
    call finalise(params)

    call calc_connect(this%atoms)
    call calc(pot, this%atoms, args_str=my_args_str, error=error)
    PASS_ERROR(error)
    call set_value(this%atoms%params, 'time', this%t)
    if (.not. get_value(this%atoms%params, 'energy', e)) &
         call system_abort("dynamicalsystem_run failed to get energy")
    if (.not. assign_pointer(this%atoms, 'force', f)) &
         call system_abort("dynamicalsystem_run failed to get forces")
    call ds_print_status(this, epot=e)
    call hook()
    if (present(trajectory)) call write(trajectory, this%atoms)

    do n=1,n_steps
       call advance_verlet1(this, dt)
       call calc(pot, this%atoms, args_str=my_args_str, error=error)
       PASS_ERROR(error)
       call advance_verlet2(this, dt, f)
       if (.not. get_value(this%atoms%params, 'energy', e)) &
            call system_abort("dynamicalsystem_run failed to get energy")
       if (.not. assign_pointer(this%atoms, 'force', f)) &
            call system_abort("dynamicalsystem_run failed to get forces")
       call ds_print_status(this, epot=e)
       call set_value(this%atoms%params, 'time', this%t)

       if (mod(n,my_hook_interval) == 0) call hook()
       if (present(trajectory) .and. mod(n,my_write_interval) == 0) call write(trajectory, this%atoms)
       if (mod(n,my_connect_interval) == 0) call calc_connect(this%atoms)
    end do

  end subroutine DynamicalSystem_run

end module Potential_module
