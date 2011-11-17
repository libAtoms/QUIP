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
#include "KIMstatus.h"

module IPModel_KIM_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module
use KIMservice

implicit none

private 

include 'IPModel_interface.h'

public :: IPModel_KIM
type IPModel_KIM
  integer(kind=kim_intptr) :: pkim = 0

  real(dp) :: cutoff = 0.0_dp    !% Cutoff for computing connection.

  character(len=FIELD_LENGTH) KIM_Model_name, KIM_Test_name

  logical :: kim_is_initialised = .false.

  integer, allocatable :: atomTypeOfZ(:)

end type IPModel_KIM

!NB in an ideal world these should not be global variables, but that requires
!NB passing in a neighObject explicitly into the iterator, rather than just a pkim
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
  character(len=KEY_CHAR_LENGTH) atom_types_list(1); pointer(p_atom_types_list, atom_types_list)

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
      call write_test_kim_file(test_kim_es, this%KIM_Test_name, this%KIM_Model_name)

      ! KIM routines to initialise pkim based on this%KIM_Test_name string
      kim_error = KIM_API_init_str_testname_f(this%pkim, string(test_kim_es)//char(0), this%KIM_Model_name)
      if (kim_error /= KIM_STATUS_OK) then
	call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_init_str_testname_f failed", kim_error)
	call system_abort("IPModel_KIM_Initialise_str: kim_api_init_str_testname_f with Test='"//trim(this%KIM_Test_name)// &
	  "' and Model='"//trim(this%KIM_Model_name)//"' failed")
      endif
   else
#endif
    ! KIM routines to initialise pkim based on this%KIM_Test_name
    kim_error = KIM_API_init_f(this%pkim, this%KIM_Test_name, this%KIM_Model_name)
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_init_f failed", kim_error)
      call system_abort("IPModel_KIM_Initialise_str: kim_api_init_f with Test='"//trim(this%KIM_Test_name)// &
	"' and Model='"//trim(this%KIM_Model_name)//"' failed")
    endif
#ifndef KIM_NO_AUTOGENERATE_TEST_KIM
  endif
#endif

  this%kim_is_initialised = .true.
  kim_iterator_current_i = -1
  nullify(kim_at)

  ! create array for conversion of atomic numbers to KIM atom type codes
  p_atom_types_list = KIM_API_get_ListAtomTypes_f(this%pkim, n_atom_types, kim_error)

  allocate(this%atomTypeOfZ(size(ElementName)))
  this%atomTypeOfZ(:) = -1
  do at_type_i=1, n_atom_types
    do j=1, KEY_CHAR_LENGTH
	if (atom_types_list(at_type_i)(j:j) == char(0)) then
	    do k=j, KEY_CHAR_LENGTH
		atom_types_list(at_type_i)(k:k) = ' '
	    end do
	    exit
	endif
    end do
    kim_at_type_code = kim_api_get_atypecode_f(this%pkim, trim(atom_types_list(at_type_i)), kim_error)
    if (kim_error /= KIM_STATUS_OK) then
       call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_get_atypecode_f failed", kim_error)
       call system_abort("failed to find name of atom type "//at_type_i//" '"//trim(atom_types_list(at_type_i))//"' in pkim")
    endif
    Z_i = find_in_array(ElementName, trim(atom_types_list(at_type_i)))
    if (Z_i < 1) then
       call system_abort("failed to match name of atom type "//at_type_i//" '"//trim(atom_types_list(at_type_i))//"'")
    endif
    Z_i = Z_i - 1
    this%atomTypeOfZ(Z_i) = kim_at_type_code
  end do

  call print("atomTypeOfZ")
  call print(this%atomTypeOfZ)

  ! register the neighbor list iterator
  kim_error = kim_api_set_data_f(this%pkim, "get_full_neigh", int(1,8), loc(quip_neighbour_iterator))
  if(kim_error /= KIM_STATUS_OK) then
    call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set_data_f failed", kim_error)
    call system_abort("IPModel_KIM_Initialise_str failed to register quip_neighbour_iterator in kim")
  endif

  kim_error = kim_api_set_data_f(this%pkim, "cutoff", int(1,8), loc(this%cutoff))
  if(kim_error /= KIM_STATUS_OK) then
    call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set_data_f failed", kim_error)
    call system_abort("IPModel_KIM_Initialise_str failed to register cutoff in kim")
  endif

  ! initialize the model
  kim_error = KIM_API_model_init_f(this%pkim)
  if (kim_error /= KIM_STATUS_OK) then
    call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_model_init_f failed", kim_error)
    call system_abort("IPModel_KIM_Initialise_str: kim_api_model_init_f failed")
  endif

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

!!NB call print("quip_neighbour_iterator")
!!NB call print("got mode " // iterator_mode // " request " // request)
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

!!NB call print("atom is "//atom)
  nneigh = atoms_n_neighbours(kim_at, atom)
!!NB call print("nneigh is  "//nneigh)
  do ji=1, nneigh
!!NB call print("trying to do atom ji "//ji)
!!NB call print("test memory " // neigh_rij(1:3,ji))
      neigh_list(ji) = atoms_neighbour(kim_at, atom, ji, diff = neigh_rij(1:3,ji))
!!NB call print("   got "//neigh_list(ji) // " rij " // neigh_rij(1:3,ji))
  end do
  !!NB neigh_rij(1:3,1:nneigh) = -neigh_rij(1:3,1:nneigh)

  outval = KIM_STATUS_OK

end function quip_neighbour_iterator

subroutine IPModel_KIM_Finalise(this)
  type(IPModel_KIM), intent(inout) :: this

  integer :: kim_error

  !! BAD - bug in gfortran 4.5 makes (apparently) this%kim_is_initialised is not set correctly for new structures, so we can't trust it
  !! if (this%kim_is_initialised) then
  !!   call kim_api_model_destroy_f(this%pkim, kim_error)
  !!   call kim_api_free_f(this%pkim, kim_error)
  !! endif

  this%KIM_Test_name = ''
  this%KIM_Model_name = ''
  this%cutoff = 0.0_dp
  this%kim_is_initialised = .false.
  kim_iterator_current_i = -1
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
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial
  character(len=*), intent(in), optional      :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  character(len=FIELD_LENGTH) :: extra_calcs_list(10)
  integer :: n_extra_calcs
  integer*8 :: ppos, pf, pe, pN, pNlocal
  integer*8 :: Nstub

  real(dp) :: virial_sym(6)
  integer :: i
  integer, allocatable :: atomTypes(:)
  integer kim_error
  integer :: bad_i(1)

  INIT_ERROR(error)

  ! register pos, atom types
  kim_error = kim_api_set_data_f(this%pkim, "coordinates", int(3*at%N,8), loc(at%pos(1,1)))
  if (kim_error /= KIM_STATUS_OK) then
    call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set_data_f for pos failed", kim_error)
    RAISE_ERROR("Failed to register coordinates field in pkim", error)
  endif

  allocate(atomTypes(at%N))
  atomTypes(1:at%N) = this%atomTypeOfZ(at%Z(1:at%N))
  if (any(atomTypes < 0)) then
    bad_i = minloc(atomTypes)
    RAISE_ERROR("Some Z mapped to an invalid atom type, e.g. atom "//bad_i(1)//" Z "//at%Z(bad_i(i)), error)
  endif

  kim_error = kim_api_set_data_f(this%pkim, "atomTypes", int(at%N,8), loc(atomTypes(1)))
  if (kim_error /= KIM_STATUS_OK) then
    call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set_data_f for atom types failed", kim_error)
    RAISE_ERROR("Failed to register atomTypes field in pkim", error)
  endif

  if (present(e)) then
    kim_error = kim_api_set_data_f(this%pkim, "energy", int(1,8), loc(e))
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set_data_f for energy failed", kim_error)
      RAISE_ERROR("Failed to create energy field in pkim", error)
    endif
    ! call set_compute energy...
    call kim_api_set2_compute_f(this%pkim,"energy", kim_error)
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set2_compute_f for energy failed", kim_error)
      RAISE_ERROR("Failed to create energy field in pkim", error)
    endif
  else
    call kim_api_set2_donotcompute_f(this%pkim,"energy", kim_error)
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set2_compute_f for energy failed", kim_error)
      RAISE_ERROR("Failed to create energy field in pkim", error)
    endif
  endif

  if (present(local_e)) then
    kim_error = kim_api_set_data_f(this%pkim, "energyPerAtom", int(at%N,8), loc(local_e(1)))
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set_data_f for energyPerAtom failed", kim_error)
      RAISE_ERROR("Failed to create energyPerAtom field in pkim", error)
    endif
    ! call set_compute energyPerAtom...
    call kim_api_set2_compute_f(this%pkim,"energyPerAtom", kim_error)
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set2_compute_f for energyPerAtom failed", kim_error)
      RAISE_ERROR("Failed to create energyPerAtom field in pkim", error)
    endif
  else
    call kim_api_set2_donotcompute_f(this%pkim,"energyPerAtom", kim_error)
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set2_compute_f for energyPerAtom failed", kim_error)
      RAISE_ERROR("Failed to create energyPerAtom field in pkim", error)
    endif
  endif

  if (present(f)) then
    kim_error = kim_api_set_data_f(this%pkim, "forces", int(3*at%N,8), loc(f(1,1)))
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set_data_f for forces failed", kim_error)
      RAISE_ERROR("Failed to create forces field in pkim", error)
    endif
    ! call set_compute forces...
    call kim_api_set2_compute_f(this%pkim,"forces", kim_error)
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set2_compute_f for forces failed", kim_error)
      RAISE_ERROR("Failed to create forces field in pkim", error)
    endif
  else
    call kim_api_set2_donotcompute_f(this%pkim,"forces", kim_error)
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set2_compute_f for forces failed", kim_error)
      RAISE_ERROR("Failed to create forces field in pkim", error)
    endif
  endif

  if (present(virial)) then
    kim_error = kim_api_set_data_f(this%pkim, "virialGlobal", int(6,8), loc(virial_sym(1)))
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set_data_f for virial failed", kim_error)
      RAISE_ERROR("Failed to create virial field in pkim", error)
    endif
    ! call set_compute virial...
    call kim_api_set2_compute_f(this%pkim,"virial", kim_error)
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set2_compute_f for virial failed", kim_error)
      RAISE_ERROR("Failed to create virial field in pkim", error)
    endif
  else
    call kim_api_set2_donotcompute_f(this%pkim,"virial", kim_error)
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set2_compute_f for virial failed", kim_error)
      RAISE_ERROR("Failed to create virial field in pkim", error)
    endif
  endif

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

  kim_error = kim_api_set_data_f(this%pkim, "numberOfAtoms", int(1,8), loc(at%N))
  if (kim_error /= KIM_STATUS_OK) then
    call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set_data_f for numberOfAtoms failed", kim_error)
    RAISE_ERROR("Failed to create numberOfAtoms field in pkim", error)
  endif
  ! deal with ghost atoms later (numberOfContributingAtoms)

  kim_at => at

  kim_error = kim_api_set_data_f(this%pkim, "neighObject", int(1,8), loc(at%N))
  if (kim_error /= KIM_STATUS_OK) then
    call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_set_data_f for neighObject failed", kim_error)
    RAISE_ERROR("Failed to create neighObject field in pkim", error)
  endif

  call kim_api_model_compute_f(this%pkim, kim_error)
  if (kim_error /= KIM_STATUS_OK) then
    call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_model_compute_f failed", kim_error)
    RAISE_ERROR("Failed to compute with KIM Model", error)
  endif

  !!!!!NB  Need to convert from virial_sym to virial
  if (present(virial)) then
    virial(1,1) = virial_sym(1)
    virial(2,2) = virial_sym(2)
    virial(3,3) = virial_sym(3)
    virial(2,3) = virial_sym(4)
    virial(3,2) = virial_sym(4)
    virial(1,3) = virial_sym(5)
    virial(3,1) = virial_sym(5)
    virial(1,2) = virial_sym(6)
    virial(2,1) = virial_sym(6)
  endif

  deallocate(atomTypes)

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

  character(len=FIELD_LENGTH) :: extra_calcs_str
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

subroutine write_test_kim_file(test_kim_es, test_name, model_name)
    type(Extendable_Str) :: test_kim_es
    character(len=*), intent(in) :: test_name, model_name

    integer kim_error
    integer :: str_len
    character(len=1) :: model_str_stub(1); pointer(p_model_str, model_str_stub)

    integer i

    p_model_str = kim_api_get_model_kim_str_f(trim(model_name), str_len, kim_error)
    if (kim_error /= KIM_STATUS_OK) then
      call KIM_API_report_error_f(__LINE__, __FILE__, "KIM_API_model_kim_str failed", kim_error)
      call system_abort("Failed to get model .kim string")
    endif

    call initialise(test_kim_es)
    do i=1, str_len
      !! call print("TEST KIM str "//model_str_stub(i))
      call concat(test_kim_es, model_str_stub(i), no_trim=.true.)
    end do

    !! call print(test_kim_es)

    !NB need to rewrite kim_es MODEL_NAME into TEST_NAME

end subroutine write_test_kim_file

end module IPModel_KIM_module
