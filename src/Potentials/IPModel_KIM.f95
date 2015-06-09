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
!X IPModel_KIM module  
!X
!% Module for KIM (Knowledgebase of Interatomic potential Models)
!%
!% The IPModel_KIM object contains all the parameters read from a
!% 'KIM_params' XML stanza. (?)
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"
#include "KIM_API_status.h"



module IPModel_KIM_module

use error_module
use system_module, only : dp, inoutput, quip_new_line, split_string_simple
use periodictable_module
use extendable_str_module
use linearalgebra_module
use dictionary_module
use paramreader_module
use atoms_module

use mpi_context_module
use QUIP_Common_module
use KIM_API

implicit none

private 

include 'IPModel_interface.h'

public :: IPModel_KIM
type IPModel_KIM
  integer(kind=kim_intptr) :: pkim = 0

  real(dp) :: cutoff = 0.0_dp    !% Cutoff for computing connection.

  character(len=STRING_LENGTH) KIM_Model_name, KIM_Test_name

  logical :: kim_is_initialised = .false.

  integer, pointer :: atomTypeOfZ(:) => null()
  character(len=512), pointer :: model_optional_outputs(:) => null()

  logical :: kim_model_has_energy, kim_model_has_forces, kim_model_has_particleenergy, kim_model_has_virial 

end type IPModel_KIM

!NB in an ideal world these should not be global variables, but that requires
!NB passing in a neighObject explicitly into the iterator, rather than just a pkim
!NB maybe something can be done with a KIM buffer?
integer :: kim_iterator_current_i = -1
type(Atoms), pointer :: kim_at

interface Initialise
  module procedure IPModel_KIM_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_KIM_Finalise
end interface Finalise

interface Print
  module procedure IPModel_KIM_Print
end interface Print

interface Calc
  module procedure IPModel_KIM_Calc
end interface Calc

contains

subroutine IPModel_KIM_Initialise_str(this, args_str, param_str)
  type(IPModel_KIM), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  integer :: n_atom_types, kim_at_type_code, Z_i, at_type_i, j, k
  character(len=KIM_KEY_STRING_LENGTH) atom_types_list(1); pointer(p_atom_types_list, atom_types_list)

  type(Dictionary) :: params
  integer :: kim_error
  type(Extendable_Str) :: test_kim_es


  ! Some day it would be nice if this could be cleaner with F2003 pointers
  real(dp) :: cutoffstub; pointer(pcutoff,cutoffstub);

  call Finalise(this)

  ! get params from QUIP, specifically KIM model and test names
  call initialise(params)
  call param_register(params, 'KIM_Model_name', PARAM_MANDATORY, this%KIM_Model_name, help_string="Name of KIM Model to initialise")
  call param_register(params, 'KIM_Test_name', '-', this%KIM_Test_name, help_string="Name of KIM Test to initialise. If '-', autogenerate test.kim file")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_KIM_Initialise_str args_str')) then
    call system_abort("IPModel_KIM_Initialise_str failed to parse args_str="//trim(args_str))
  endif
  call finalise(params)

#ifdef KIM_NO_AUTOGENERATE_TEST_KIM
   if (trim(this%KIM_Test_name) == '-') then
     call system_abort("IPModel_KIM_Initialise_str: no support for KIM_Test_name=- autogenerate test.kim file")
   endif
#endif

#ifndef KIM_NO_AUTOGENERATE_TEST_KIM
  if (trim(this%KIM_Test_name) == '-') then
      ! create KIM test file string from model kim string
      call write_test_kim_file_from_model(this, test_kim_es)

      ! KIM routines to initialise pkim based on this%KIM_Test_name string
      kim_error = KIM_API_string_init_f(this%pkim, string(test_kim_es)//char(0), this%KIM_Model_name)
      if (KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_init_str_testname_f failed", kim_error) /= KIM_STATUS_OK) then
	call print(test_kim_es)
	call system_abort("IPModel_KIM_Initialise_str: kim_api_init_str_testname_f with Test='"//trim(this%KIM_Test_name)// &
	  "' and Model='"//trim(this%KIM_Model_name)//"' failed")
      endif

   else
#endif
    ! KIM routines to initialise pkim based on this%KIM_Test_name
    kim_error = KIM_API_init_f(this%pkim, this%KIM_Test_name, this%KIM_Model_name)
    if (KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_init_f failed", kim_error) /= KIM_STATUS_OK) &
      call system_abort("IPModel_KIM_Initialise_str: kim_api_init_f with Test='"//trim(this%KIM_Test_name)// &
	"' and Model='"//trim(this%KIM_Model_name)//"' failed")
#ifndef KIM_NO_AUTOGENERATE_TEST_KIM
  endif
#endif

  this%kim_is_initialised = .true.
  kim_iterator_current_i = -1
  nullify(kim_at)

  ! create array for conversion of atomic numbers to KIM atom type codes
  p_atom_types_list = KIM_API_get_model_partcl_typs_f(this%pkim, n_atom_types, kim_error)

  allocate(this%atomTypeOfZ(size(ElementName)))
  this%atomTypeOfZ(:) = -1
  do at_type_i=1, n_atom_types
    ! replace null termination with space padding (C -> Fortran)
    do j=1, KIM_KEY_STRING_LENGTH
	if (atom_types_list(at_type_i)(j:j) == char(0)) then
	    do k=j, KIM_KEY_STRING_LENGTH
		atom_types_list(at_type_i)(k:k) = ' '
	    end do
	    exit
	endif
    end do
    kim_at_type_code = KIM_API_get_partcl_type_code_f(this%pkim, trim(atom_types_list(at_type_i)), kim_error)
    if (KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_get_partcl_type_code_f failed", kim_error) /= KIM_STATUS_OK) &
       call system_abort("failed to find name of atom type "//at_type_i//" '"//trim(atom_types_list(at_type_i))//"' in pkim")
    Z_i = find_in_array(ElementName, trim(atom_types_list(at_type_i)))
    if (Z_i < 1) then
       call system_abort("failed to find atomic number for name of atom type "//at_type_i//" '"//trim(atom_types_list(at_type_i))//"'")
    endif
    Z_i = Z_i - 1
    this%atomTypeOfZ(Z_i) = kim_at_type_code
  end do

  ! register the neighbor list iterator
  kim_error = kim_api_set_data_f(this%pkim, "get_neigh", int(1,8), loc(quip_neighbour_iterator))
  if (KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set_data_f failed", kim_error) /= KIM_STATUS_OK) &
    call system_abort("IPModel_KIM_Initialise_str failed to register quip_neighbour_iterator in kim")

  ! register the cutoff
  kim_error = kim_api_set_data_f(this%pkim, "cutoff", int(1,8), loc(this%cutoff))
  if (KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set_data_f failed", kim_error) /= KIM_STATUS_OK) &
    call system_abort("IPModel_KIM_Initialise_str failed to register cutoff in kim")

  ! initialize the model
  kim_error = KIM_API_model_init_f(this%pkim)
  if (KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_model_init_f failed", kim_error) /= KIM_STATUS_OK) &
    call system_abort("IPModel_KIM_Initialise_str: kim_api_model_init_f failed")

end subroutine IPModel_KIM_Initialise_str

function quip_neighbour_iterator(pkim, iterator_mode, request, atom, nneigh, p_neigh_list, p_neigh_rij) result(outval)
  integer(kind=kim_intptr) :: pkim
  integer :: iterator_mode, request, atom, nneigh

  !NB 512 should be replaced with a KIM constant
  integer, save :: neigh_list(512)
  integer :: neigh_list_stub(1); pointer (p_neigh_list, neigh_list_stub)

  !NB 512 should be replaced with a KIM constant
  real(dp), save :: neigh_rij(3,512)
  real(dp) :: neigh_rij_stub(3,1); pointer (p_neigh_rij, neigh_rij_stub)

  integer :: outval

  integer :: ji

  p_neigh_rij = loc(neigh_rij(1,1))
  p_neigh_list = loc(neigh_list(1))

  if (iterator_mode == 1) then ! locator mode
    if (request < 1 .or. request > kim_at%N) then
      outval = KIM_STATUS_NEIGH_INVALID_REQUEST
      return
    endif
    atom = request
  else if (iterator_mode == 0) then ! iterator mode
    if (request == 0) then
       kim_iterator_current_i = 1
       outval = KIM_STATUS_NEIGH_ITER_INIT_OK
       return
    else if (request == 1) then
      if (kim_iterator_current_i > kim_at%N) then
        outval = KIM_STATUS_NEIGH_ITER_PAST_END
        return
      endif
      atom = kim_iterator_current_i
      kim_iterator_current_i = kim_iterator_current_i + 1
    else ! other iterator requests
      outval = KIM_STATUS_NEIGH_INVALID_REQUEST
      return
    endif
  else ! other mode
    outval = KIM_STATUS_NEIGH_INVALID_MODE
    return
  endif

  nneigh = n_neighbours(kim_at, atom)
  do ji=1, nneigh
      neigh_list(ji) = neighbour(kim_at, atom, ji, diff = neigh_rij(1:3,ji))
  end do

  outval = KIM_STATUS_OK

end function quip_neighbour_iterator

subroutine IPModel_KIM_Finalise(this)
  type(IPModel_KIM), intent(inout) :: this

  integer :: kim_error

  !! BAD - bug in gfortran 4.5 makes (apparently) this%kim_is_initialised is not set correctly for new structures, so we can't trust it
  if (this%kim_is_initialised) then
    kim_error = kim_api_model_destroy_f(this%pkim)
    call kim_api_free_f(this%pkim, kim_error)
  endif

  this%KIM_Test_name = ''
  this%KIM_Model_name = ''
  this%cutoff = 0.0_dp
  this%kim_is_initialised = .false.
  kim_iterator_current_i = -1
  if (this%pkim /= 0) then
     if (associated(this%AtomTypeOfZ)) deallocate(this%AtomTypeOfZ)
     if (associated(this%model_optional_outputs)) deallocate(this%model_optional_outputs)
  endif
  this%pkim = 0
  nullify(kim_at)
end subroutine IPModel_KIM_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator: this routine computes energy, forces and the virial.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_KIM_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_KIM), intent(inout) :: this
  type(Atoms), intent(inout), target :: at
  real(dp), intent(out), target, optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), target, optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial
  character(len=*), intent(in), optional      :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  character(len=STRING_LENGTH) :: extra_calcs_list(10)
  integer :: n_extra_calcs
  integer*8 :: ppos, pf, pe, pN, pNlocal
  integer*8 :: Nstub

  real(dp) :: virial_sym(6)
  integer :: i
  integer :: numberAtomTypes
  integer, allocatable :: atomTypes(:)
  integer kim_error
  integer :: bad_i(1)

  real(dp), target :: t_e
  real(dp), allocatable, target :: t_local_e(:), t_f(:,:)
  real(dp), pointer :: p_e, p_local_e(:), p_f(:,:)

  INIT_ERROR(error)

  ! check for unsupported output quantities
  if (present(local_virial)) then
    RAISE_ERROR("IPModel_KIM_Calc can't handle local_virial yet", error)
  endif
  if (present(args_str)) then
    if (len_trim(args_str) > 0) then
      n_extra_calcs = parse_extra_calcs(args_str, extra_calcs_list)
      if (n_extra_calcs > 0) then
	RAISE_ERROR("No extra calcs supported, found "//n_extra_calcs//" in args_str " // trim(args_str), error)
      endif
    endif
  endif

  ! allocate atomTypes array for conversion from atomic number to atom type
  allocate(atomTypes(at%N))
  atomTypes(1:at%N) = this%atomTypeOfZ(at%Z(1:at%N))
  if (any(atomTypes < 0)) then
    bad_i = minloc(atomTypes)
    RAISE_ERROR("Some Z mapped to an invalid atom type, e.g. atom "//bad_i(1)//" Z "//at%Z(bad_i(1)), error)
  endif
  numberAtomTypes = maxval(this%atomTypeOfZ)

  ! register config input information
  call kim_api_setm_data_f(this%pkim, kim_error, &
    "numberOfParticles", int(1,8), loc(at%N), KIM_L(.true.), &
    "coordinates", int(3*at%N,8), loc(at%pos(1,1)), KIM_L(.true.), &
    "numberParticleTypes", int(numberAtomTypes,8), loc(numberAtomTypes), KIM_L(.true.), &
    "particleTypes", int(at%N,8), loc(atomTypes(1)), KIM_L(.true.))
  ! deal with ghost atoms later (numberOfContributingAtoms)
  if (KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_setm_data_f for inputs failed", kim_error) /= KIM_STATUS_OK) then
    RAISE_ERROR("Failed to setm_data for numberOfParticles, coordinates, particleTypes fields in pkim", error)
  endif

  ! it would be nice to know if model outputs are optional or not, so we don't have to allocate them
  ! if we don't need to
  if (present(e)) then
     if (.not. this%kim_model_has_energy) &
       call system_abort("IPModel_KIM_calc got energy, but kim model does not return energy")
     p_e => e
  else
     p_e => t_e
  endif
  if (present(f)) then
     if (.not. this%kim_model_has_forces) &
       call system_abort("IPModel_KIM_calc got forces, but kim model does not return forces")
     p_f => f
  else
     if (this%kim_model_has_forces) then
	allocate(t_f(3, at%N))
     else
	allocate(t_f(1, 1))
     endif
     p_f => t_f
  endif
  if (present(local_e)) then
     if (.not. this%kim_model_has_particleenergy) &
       call system_abort("IPModel_KIM_calc got local_e, but kim model does not return particleEnergy")
     p_f => f
  else
     if (this%kim_model_has_particleenergy) then
	allocate(t_local_e(at%N))
     else
	allocate(t_local_e(1))
     endif
     p_local_e => t_local_e
  endif
  ! virial is already dealt with because virial_sym needs to be defined separately from virial (6 vs. 3x3)

  ! register config output information
  call kim_api_setm_data_f(this%pkim, kim_error, &
    "energy", int(1,8), loc(p_e), KIM_L(.true.), &
    "particleEnergy", int(at%N,8), loc(p_local_e(1)), KIM_L(.true.), &
    "forces", int(3*at%N,8), loc(p_f(1,1)), KIM_L(.true.), &
    "virial", int(6,8), loc(virial_sym(1)), KIM_L(.true.))
  if (KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_setm_data_f for output failed", kim_error) /= KIM_STATUS_OK) then
    RAISE_ERROR("Failed to setm_data for energy, particleEnergy, forces, virial  fields in pkim", error)
  endif

  ! set compute for config output information
  call kim_api_setm_compute_f(this%pkim, kim_error, &
    "energy", KIM_COMPUTE_L(present(e)), KIM_L(.true.), &
    "particleEnergy", KIM_COMPUTE_L(present(local_e)), KIM_L(.true.), &
    "forces", KIM_COMPUTE_L(present(f)), KIM_L(.true.), &
    "virial", KIM_COMPUTE_L(present(virial)), KIM_L(.true.))
  if (KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_setm_compute_f for outputs failed", kim_error) /= KIM_STATUS_OK) then
    RAISE_ERROR("Failed to setm_compute for energy, particleEnergy, forces, virial field in pkim", error)
  endif

  ! set data for neighbor stuff

  kim_at => at

  kim_error = kim_api_set_data_f(this%pkim, "neighObject", int(1,8), loc(at%N))
  if (KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set_data_f for neighObject failed", kim_error) /= KIM_STATUS_OK) then
    RAISE_ERROR("Failed to create neighObject field in pkim", error)
  endif

  kim_error = kim_api_model_compute_f(this%pkim)
  if (KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_model_compute_f failed", kim_error) /= KIM_STATUS_OK) then
    RAISE_ERROR("Failed to compute with KIM Model", error)
  endif

  if (present(virial)) then
    virial(1,1) = -virial_sym(1)
    virial(2,2) = -virial_sym(2)
    virial(3,3) = -virial_sym(3)
    virial(2,3) = -virial_sym(4)
    virial(3,2) = -virial_sym(4)
    virial(1,3) = -virial_sym(5)
    virial(3,1) = -virial_sym(5)
    virial(1,2) = -virial_sym(6)
    virial(2,1) = -virial_sym(6)
  endif

  deallocate(atomTypes)
  if (allocated(t_f)) deallocate(t_f)
  if (allocated(t_local_e)) deallocate(t_local_e)

end subroutine IPModel_KIM_Calc

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of KIM parameters: number of different types, cutoff radius, atomic numbers, etc.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_KIM_Print (this, file)
  type(IPModel_KIM), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_KIM : KIM", file=file)
  call Print("IPModel_KIM : Test name = " // trim(this%KIM_Test_name) // " Model name = " // &
    trim(this%KIM_Model_Name) // " cutoff = " // this%cutoff, file=file)

end subroutine IPModel_KIM_Print

function parse_extra_calcs(args_str, extra_calcs_list) result(n_extra_calcs)
  character(len=*), intent(in) :: args_str
  character(len=*), intent(out) :: extra_calcs_list(:)
  integer :: n_extra_calcs

  character(len=STRING_LENGTH) :: extra_calcs_str
  type(Dictionary) :: params

  n_extra_calcs = 0
  call initialise(params)
  call param_register(params, "extra_calcs", "", extra_calcs_str, help_string="No help yet.  This source file was $LastChangedBy$")
  if (param_read_line(params, args_str, ignore_unknown=.true.,task='parse_extra_calcs')) then
    if (len_trim(extra_calcs_str) > 0) then
      call split_string_simple(extra_calcs_str, extra_calcs_list, n_extra_calcs, ":")
    end if
  end if
  call finalise(params)

end function parse_extra_calcs

subroutine write_test_kim_file_from_model(this, test_kim_es)
    type(IPModel_KIM), intent(inout) :: this
    type(Extendable_Str) :: test_kim_es

    integer kim_error
    integer :: i, str_len
    character(len=1) :: model_str_stub(1); pointer(p_model_str, model_str_stub)

    integer :: sec_start, sec_end
    character(len=80) :: cur_sec_name
    integer :: t_index
    integer(kind=kim_intptr) :: t_pkim = 0
    character(len=1000) :: test_model_output


    p_model_str = kim_api_get_model_kim_str_f(trim(this%KIM_model_name), str_len, kim_error)
    if (KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_get_model_kim_str failed", kim_error) /= KIM_STATUS_OK) &
      call system_abort("Failed to get model .kim string")

    call initialise(test_kim_es)
    do i=1, str_len
      call concat(test_kim_es, model_str_stub(i), no_trim=.true.)
    end do

    ! keep SUPPORTED_ATOM/PARTICLE_TYPES - will fail at compute type if it doesn't match

    ! get rid of MODEL_PARAMETERS section
    call find_section(test_kim_es, "MODEL_PARAMETERS", sec_start, sec_end)
    if (sec_start > 0) then
       call substr_replace(test_kim_es, sec_start, sec_end, "")
    endif

    ! set CONVENTIONS section to QUIP's restrictions
    call find_section(test_kim_es, "CONVENTIONS", sec_start, sec_end)
    if (sec_start > 0) then
       call substr_replace(test_kim_es, sec_start, sec_end, &
	   "CONVENTIONS:"//quip_new_line// &
	   "  OneBasedLists flag"//quip_new_line// &
	   "  Neigh_BothAccess flag"//quip_new_line// &
	   "  NEIGH_RVEC_F flag"//quip_new_line//quip_new_line// &
	   "################################################################################")
    else
    	call system_abort("write_test_kim_file_from_model failed to find CONVENTIONS section")
    endif

    ! set MODEL_INPUT section to what QUIP supports
    call find_section(test_kim_es, "MODEL_INPUT", sec_start, sec_end)
    if (sec_start > 0) then
       call substr_replace(test_kim_es, sec_start, sec_end, &
	    "MODEL_INPUT:"//quip_new_line// &
	    "# Name                      Type         Unit                Shape              Requirements"//quip_new_line// &
	    "numberOfParticles           integer      none                []"//quip_new_line// &
	    "numberParticleTypes         integer      none                []"//quip_new_line// &
	    "particleTypes               integer      none                [numberOfParticles]"//quip_new_line// &
	    "coordinates                 real*8       length              [numberOfParticles,3]"//quip_new_line// &
	    "get_neigh                   method       none                []"//quip_new_line// &
	    "neighObject                 pointer      none                []"//quip_new_line// &
	    "#######################################################################################################")
    else
    	call system_abort("write_test_kim_file_from_model failed to find MODEL_INPUT section")
    endif

    ! promise to use every model output
    kim_error = KIM_API_model_info_f(t_pkim, this%KIM_Model_name)
    if (KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_model_info_f failed", kim_error) /= KIM_STATUS_OK) &
      call system_abort("Failed to get info for kim model "//trim(this%KIM_model_name))

    t_index = KIM_API_get_index_f(t_pkim, "energy", kim_error)
    this%kim_model_has_energy = (kim_error == KIM_STATUS_OK)
    t_index = KIM_API_get_index_f(t_pkim, "forces", kim_error)
    this%kim_model_has_forces = (kim_error == KIM_STATUS_OK)
    t_index = KIM_API_get_index_f(t_pkim, "particleEnergy", kim_error)
    this%kim_model_has_particleenergy = (kim_error == KIM_STATUS_OK)
    t_index = KIM_API_get_index_f(t_pkim, "virial", kim_error)
    this%kim_model_has_virial = (kim_error == KIM_STATUS_OK)

    test_model_output = &
	 "destroy                     method       none                []"//quip_new_line// &
	 "compute                     method       none                []"//quip_new_line// &
	 "cutoff                      real*8       length              []"//quip_new_line
    if (this%kim_model_has_energy) &
      test_model_output = trim(test_model_output)// &
	 "energy                      real*8       energy              []"//quip_new_line
    if (this%kim_model_has_forces) &
      test_model_output = trim(test_model_output)// &
	 "forces                      real*8       force      [numberOfParticles,3]"//quip_new_line
    if (this%kim_model_has_particleenergy) &
      test_model_output = trim(test_model_output)// &
	 "particleEnergy              real*8       energy     [numberOfParticles]"//quip_new_line
    if (this%kim_model_has_virial) &
      test_model_output = trim(test_model_output)// &
	 "virial                      real*8       energy               [6]"//quip_new_line

    call find_section(test_kim_es, "MODEL_OUTPUT", sec_start, sec_end)
    if (sec_start > 0) then
       call substr_replace(test_kim_es, sec_start, sec_end, &
	    "MODEL_OUTPUT:"//quip_new_line// &
	    "# Name                      Type         Unit                Shape        "//quip_new_line// &
	    trim(test_model_output)// &
	    "#######################################################################################################")
    else
    	call system_abort("write_test_kim_file_from_model failed to find MODEL_OUTPUT section")
    endif

end subroutine write_test_kim_file_from_model

subroutine find_section(test_kim_es, sec_name, sec_start, sec_end)
    type(Extendable_Str) :: test_kim_es
    character(len=*), intent(in) :: sec_name
    integer, intent(out) :: sec_start, sec_end

    integer :: cur_sec_loc, next_sec_loc
    character(len=80) :: cur_sec_name

    sec_start = -1
    sec_end = -1

    cur_sec_loc = find_section_label(test_kim_es, 1, cur_sec_name)
    do while (cur_sec_loc > 0)
      if (trim(cur_sec_name) == trim(sec_name)) then
	 sec_start = cur_sec_loc
	 next_sec_loc = find_section_label(test_kim_es, cur_sec_loc+len_trim(cur_sec_name)+1)
	 if (next_sec_loc <= 0) then
	    sec_end = test_kim_es%len
	 else
	    sec_end = next_sec_loc - 2
	 endif
	 return
      endif
      cur_sec_loc = find_section_label(test_kim_es, cur_sec_loc+len_trim(cur_sec_name)+1, cur_sec_name)
    end do 

end subroutine find_section

function find_section_label(test_kim_es, start, sec_name) result(sec_loc)
   type(Extendable_Str) :: test_kim_es
   integer, intent(in) :: start
   character(len=*), optional :: sec_name
   integer :: sec_loc

   integer :: i
   logical :: after_newline, in_section, is_newline, is_section_char, is_section_end
   character(len=80) :: t_sec_name

   after_newline = (start == 1)

   sec_loc = -1
   in_section = .false.
   do i=start, test_kim_es%len
      is_newline = (substr(test_kim_es,i,i) == quip_new_line)
      is_section_char = (scan(substr(test_kim_es,i,i),"ABCDEFGHIJKLMNOPQRSTUVWXYZ_") > 0)
      is_section_end = (substr(test_kim_es,i,i) == ':')
      if (after_newline) then
	 if (is_section_char) then
	    in_section = .true.
	    t_sec_name = substr(test_kim_es, i, i)
	 endif
      else ! not after newline
	 if (in_section) then
	    if (is_section_end) then
	       sec_loc = i-len_trim(t_sec_name)
	       if (present(sec_name)) sec_name = trim(t_sec_name)
	       return
	    else if (is_section_char) then
	       t_sec_name = trim(t_sec_name) // substr(test_kim_es, i, i)
	    else
	       sec_loc = -1
	       in_section = .false.
	    endif
	 endif
      endif
      after_newline = is_newline
   end do

   sec_loc = -1
   in_section = .false.
 
end function find_section_label

function KIM_L(l)
   logical, intent(in) :: l
   integer :: KIM_L

   if (l) then
      KIM_L = 1
   else
      KIM_L = 0
   endif
end function KIM_L

function KIM_COMPUTE_L(l)
   logical, intent(in) :: l
   integer :: KIM_COMPUTE_L

   if (l) then
      KIM_COMPUTE_L = KIM_COMPUTE_TRUE
   else
      KIM_COMPUTE_L = KIM_COMPUTE_FALSE
   endif
end function KIM_COMPUTE_L

end module IPModel_KIM_module
