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
  integer*8 :: pkim

  real(dp) :: cutoff = 0.0_dp    !% Cutoff for computing connection.

  character(len=FIELD_LENGTH) KIM_Model_name, KIM_TEST_name
  integer :: iter_current_i

  logical :: kim_is_initialised = .false.

end type IPModel_KIM

type KIM_iterator_data
   integer :: current_i
   type(Atoms), pointer :: at
end type KIM_iterator_data

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

  type(Dictionary) :: params
  real(dp) :: cutoffstub; pointer(pcutoff,cutoffstub);
  integer*8 :: one = 1

  call Finalise(this)

  call initialise(params)
  call param_register(params, 'KIM_Model_name', PARAM_MANDATORY, this%KIM_Model_name, help_string="Name of KIM Model to initialise")
  call param_register(params, 'KIM_TEST_name', PARAM_MANDATORY, this%KIM_TEST_name, help_string="Name of KIM TEST to initialise")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_KIM_Initialise_str args_str')) then
    call system_abort("IPModel_KIM_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  ! KIM routines to initialise pkim based on this%KIM_Model_name string
  ! .
  if (kim_api_init_f(this%pkim, this%KIM_TEST_name, this%KIM_Model_name).ne.1) then
    call system_abort("IPModel_KIM_Initialise_str: kim_api_init_f with Test='"//trim(this%KIM_TEST_name)// &
      "' and Model='"//trim(this%KIM_Model_name)//"' failed")
  endif
  this%kim_is_initialised = .true.

  ! this will not be hardwired some day - kim_api_init(this%pkim,"QUIP_Model_SW_v0_f") for example
  ! call sample_01_lj_quip_cpp_init(this%pkim)
  call QUIP_Model_SW_v0_f_init(this%pkim)

  if(kim_api_set_data_f(this%pkim,"neighIterator",int(1,8),loc(quip_neighbour_iterator)).ne.1) then
    call system_abort("IPModel_KIM_Initialise_str failed to register quip_neighbour_iterator in kim")
  endif

  ! would be nice to get rid of stub and combine these 2 functions - requires F2003 pointers?
  pcutoff = kim_api_get_data_f(this%pkim,"cutoff")
  this%cutoff = cutoffstub! value of cutoff from KIM Model

end subroutine IPModel_KIM_Initialise_str

subroutine quip_neighbour_iterator(pneigh_obj, pneighbours, n_neighbours, restart, pr, prij)
  use iso_c_binding
  ! intent(in) neigh_obj, restart, r
  ! intent(out) neighbours, n_neighbours, rij
  integer*8, pointer :: pneigh_obj
  integer*8 :: pneighbours
  integer :: n_neighbours, restart
  integer*8 :: pr, prij
  type(c_ptr) :: ppneigh_obj

  type(KIM_iterator_data), pointer :: ki
  integer :: neighbours(1)
  real(dp) :: r(3,512), rij(3,512)
  pointer(pneighbours, neighbours)
  pointer(pr, r)
  pointer(prij, rij)

  integer :: cur_n_neighbours, i, ji

print *, "quip_neighbour_iterator got pneigh_obj ", pneigh_obj, c_loc(pneigh_obj)

  ppneigh_obj = c_loc(pneigh_obj)
  call c_f_pointer(ppneigh_obj, ki)

  ! atoms_n_neighbours(at, i)
  ! atoms_neighbours(at, i, ji, rij)

  if (restart == 1) then
    ki%current_i = 1
    return 1
  else
    ki%current_i = ki%current_i + 1
  endif
  i = ki%current_i

  cur_n_neighbours = atoms_n_neighbours(ki%at, i)
  if (i == ki%at%N) then ! some day this will change to a simple flag
    n_neighbours = -1
  else
    n_neighbours = cur_n_neighbours
  endif

  neighbours(1) = cur_n_neighbours+2
  neighbours(2) = i-1
  do ji=1, cur_n_neighbours
    neighbours(ji+2) = atoms_neighbour(ki%at, i, ji, diff=rij(1:3,ji)) -1
  end do
  return 
end subroutine quip_neighbour_iterator

subroutine IPModel_KIM_Finalise(this)
  type(IPModel_KIM), intent(inout) :: this

  !! BAD - bug in gfortran 4.5 makes (apparently) this%kim_is_initialised is not set correctly for new structures, so we can't trust it
  !! if (this%kim_is_initialised) call kim_api_free_f(this%pkim)
  this%KIM_TEST_name = ''
  this%KIM_Model_name = ''
  this%cutoff = 0.0_dp
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

  character(len=STRING_LENGTH) :: extra_calcs_list(10)
  integer :: n_extra_calcs
  type(KIM_iterator_data) :: ki
  integer*8 :: ppos, pf, pe, pN, pNlocal
  integer*8 :: Nstub

  integer :: i
  integer, allocatable :: ZofType(:)

  INIT_ERROR(error)

  ! would be nice to compress these 2 into 1 line, hide int(3*at%N,8) stuff
  ppos = loc(at%pos(1,1))
  if (kim_api_set_data_f(this%pkim, "coordinates", int(3*at%N,8), ppos) /= 1) then
    RAISE_ERROR("Failed to create coordinates field in pkim", error)
  endif

  allocate(ZofType(maxval(at%Z)))
  do i=1, maxval(at%Z)
        ZofType(i) = i
  end do
  if (kim_api_set_data_f(this%pkim, "ZofType", int(maxval(at%Z),8), loc(ZofType(1))) /= 1) then
    RAISE_ERROR("Failed to create coordinates field in pkim", error)
  endif
  if (kim_api_set_data_f(this%pkim, "atomTypes", int(at%N,8), loc(at%Z(1))) /= 1) then
    RAISE_ERROR("Failed to create coordinates field in pkim", error)
  endif

  if (present(e)) then
    pe = loc(e)
    if (kim_api_set_data_f(this%pkim, "energy", int(1,8), pe) /= 1) then
      RAISE_ERROR("Failed to create energy field in pkim", error)
    endif
    ! call set_compute energy...
    call kim_api_set2_compute_f(this%pkim,"energy")
  else
    call kim_api_set2_uncompute_f(this%pkim,"energy")
  endif

  if (present(f)) then
    pf = loc(f(1,1))
    if (kim_api_set_data_f(this%pkim, "forces", int(3*at%N,8), pf) /= 1) then
      RAISE_ERROR("Failed to create forces field in pkim", error)
    endif
    ! call set_compute energy...
    call kim_api_set2_compute_f(this%pkim,"forces")
  else
    call kim_api_set2_uncompute_f(this%pkim,"forces")
  endif

  Nstub = at%N; pN = loc(Nstub)
  if (kim_api_set_data_f(this%pkim, "numberOfAtoms", int(1,8), pN) /= 1) then
    RAISE_ERROR("Failed to create forces field in pkim", error)
  endif
  pNlocal = loc(at%N) ! we'll deal with ghost atoms eventually
  if (kim_api_set_data_f(this%pkim, "numberOfLocals", int(1,8), pNlocal) /= 1) then
    RAISE_ERROR("Failed to create forces field in pkim", error)
  endif

  if (present(local_e)) then
    RAISE_ERROR("IPModel_KIM_Calc can't handle local_e yet", error)
  endif
  if (present(virial)) then
    RAISE_ERROR("IPModel_KIM_Calc can't handle virial yet", error)
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

  ki%at => at
print *, "IPModel_KIM_Calc registering ki in neighObject ", loc(ki)
  if (kim_api_set_data_f(this%pkim, "neighObject", int(1,8), loc(ki)) /= 1) then
    RAISE_ERROR("Failed to create forces field in pkim", error)
  endif

print *, "IPModel_KIM_Calc calling kim_api_model_compute"
  call kim_api_model_compute(this%pkim)

  deallocate(ZofType)

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
  call Print("IPModel_KIM : Test name = " // trim(this%KIM_Test_Name) // " Model name = " // &
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

end module IPModel_KIM_module
