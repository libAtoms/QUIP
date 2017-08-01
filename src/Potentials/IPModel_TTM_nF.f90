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
!X IPModel_TTM_nF
!X
!% TTM_nF module
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_TTM_nF_module

use error_module
use system_module, only : dp, inoutput, print, verbosity_push_decrement, verbosity_pop, operator(//), lower_case
use dictionary_module
use units_module, only : KCAL_MOL
use linearalgebra_module, only: is_diagonal
use topology_module, only: find_water_monomer
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module

use mpi_context_module
use QUIP_Common_module


use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE

implicit none

#ifdef HAVE_TTM_NF

integer(C_INT), protected, bind(C,name="POT_QTIP4PF") :: POT_QTIP4PF
integer(C_INT), protected, bind(C,name="POT_TTM2F") :: POT_TTM2F
integer(C_INT), protected, bind(C,name="POT_TTM3F") :: POT_TTM3F
integer(C_INT), protected, bind(C,name="POT_TTM4F") :: POT_TTM4F

interface
   subroutine ttm_nf_energy(pot_type,n,pos,energy) bind(C)
      use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE
      integer(C_INT), value :: pot_type
      integer(C_SIZE_T), value :: n
      real(C_DOUBLE), dimension(*) :: pos
      real(C_DOUBLE) :: energy
   endsubroutine ttm_nf_energy
endinterface

interface
   subroutine ttm_nf_energy_gradient(pot_type,n,pos,energy,gradient) bind(C)
      use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE
      integer(C_INT), value :: pot_type
      integer(C_SIZE_T), value :: n
      real(C_DOUBLE), dimension(*) :: pos
      real(C_DOUBLE) :: energy
      real(C_DOUBLE), dimension(*) :: gradient
   endsubroutine ttm_nf_energy_gradient
endinterface

#endif

private

include 'IPModel_interface.h'

public :: IPModel_TTM_nF
type IPModel_TTM_nF
  !integer :: n_types = 0
  !integer, allocatable :: atomic_num(:), type_of_atomic_num(:)

  real(dp) :: cutoff = 0.0_dp
  integer :: potential_type

  character(len=STRING_LENGTH) :: label

end type IPModel_TTM_nF

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_TTM_nF), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_TTM_nF_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_TTM_nF_Finalise
end interface Finalise

interface Print
  module procedure IPModel_TTM_nF_Print
end interface Print

interface Calc
  module procedure IPModel_TTM_nF_Calc
end interface Calc

contains

subroutine IPModel_TTM_nF_Initialise_str(this, args_str, param_str)
  type(IPModel_TTM_nF), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  character(len=STRING_LENGTH) :: potential_type

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'potential_type', 'TTM4F', potential_type, help_string="Type of potential. Allowed values: QTIP4PF, TTM2F, TTM3F, TTM4F$")

  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_TTM_nF_Initialise_str args_str')) then
    call system_abort("IPModel_TTM_nF_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  this%cutoff = 2.0_dp
  select case(lower_case(trim(potential_type)))
  case("qtip4pf")
     this%potential_type = POT_QTIP4PF
  case("ttm2f")
     this%potential_type = POT_TTM2F
  case("ttm3f")
     this%potential_type = POT_TTM3F
  case("ttm4f")
     this%potential_type = POT_TTM4F
  case default
     call system_abort("IPModel_TTM_nF_Initialise_str: unknown potential_type "//potential_type)
  endselect


  !  Add initialisation code here

end subroutine IPModel_TTM_nF_Initialise_str

subroutine IPModel_TTM_nF_Finalise(this)
  type(IPModel_TTM_nF), intent(inout) :: this

  ! Add finalisation code here
  this%cutoff = 0.0_dp

  this%label = ''
end subroutine IPModel_TTM_nF_Finalise

subroutine IPModel_TTM_nF_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_TTM_nF), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   integer(C_SIZE_T) :: nWater
   integer :: i, a
   integer, dimension(:,:), allocatable :: water_monomer_index
   real(dp) :: energy
   real(dp), dimension(3) :: lattice
   real(dp), dimension(:), allocatable :: pos, gradient

   INIT_ERROR(error)
#ifndef HAVE_TTM_NF
   RAISE_ERROR('IPModel_TTM_nF_Calc - not linked to the TTM_nF library',error)
#endif

   if (present(e)) e = 0.0_dp
   
   if (present(local_e)) then
      RAISE_ERROR('IPModel_TTM_nF_Calc - local energies not implemented',error)
      call check_size('Local_E',local_e,(/at%N/),'IPModel_TTM_nF_Calc', error)
      local_e = 0.0_dp
   endif
   
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_TTM_nF_Calc', error)
      f = 0.0_dp
   end if
   
   if (present(virial)) then
      RAISE_ERROR('IPModel_TTM_nF_Calc - virials not implemented',error)
      virial = 0.0_dp
   endif

   if (present(local_virial)) then
      RAISE_ERROR('IPModel_TTM_nF_Calc - local virials not implemented',error)
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_TTM_nF_Calc', error)
      local_virial = 0.0_dp
   endif

   nWater = count(at%Z==8)
   allocate(water_monomer_index(3,nWater))
   call find_water_monomer(at,water_monomer_index,error=error)

   allocate(pos(3*at%N),gradient(3*at%N))

   do i = 1, nWater
      do a = 1, 3
         pos(a + (i-1)*9)  = at%pos(a,water_monomer_index(1,i))
         pos(a + (i-1)*9 + 3)  = at%pos(a,water_monomer_index(2,i))
         pos(a + (i-1)*9 + 6)  = at%pos(a,water_monomer_index(3,i))
      enddo
   enddo
 
#ifdef HAVE_TTM_NF
   if(present(f)) then
      call ttm_nf_energy_gradient(this%potential_type,nWater,pos,energy,gradient)
   endif

   if(present(e)) then
      call ttm_nf_energy(this%potential_type,nWater,pos,energy)
   endif
#endif

   if (present(e)) e = energy * KCAL_MOL

   if(present(f)) then
      do i = 1, nWater
         do a = 1, 3
            f(a,water_monomer_index(1,i)) = - gradient(a + (i-1)*9)
            f(a,water_monomer_index(2,i)) = - gradient(a + (i-1)*9 + 3)
            f(a,water_monomer_index(3,i)) = - gradient(a + (i-1)*9 + 6)
         enddo
      enddo
       f = f * KCAL_MOL
   endif

   if(allocated(water_monomer_index)) deallocate(water_monomer_index)
   if(allocated(pos)) deallocate(pos)
   if(allocated(gradient)) deallocate(gradient)

end subroutine IPModel_TTM_nF_Calc


subroutine IPModel_TTM_nF_Print(this, file)
  type(IPModel_TTM_nF), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti

  call Print("IPModel_TTM_nF : TTM_nF Potential", file=file)

  if (this%potential_type == POT_QTIP4PF) then
     call print("IPModel_TTM_nF : type is QTIP4PF")
  elseif( this%potential_type == POT_TTM2F ) then
     call print("IPModel_TTM_nF : type is TTM2F")
  elseif( this%potential_type == POT_TTM3F ) then
     call print("IPModel_TTM_nF : type is TTM3F")
  elseif( this%potential_type == POT_TTM4F ) then
     call print("IPModel_TTM_nF : type is TTM4F")
  else
     call system_abort("IPModel_TTM_nF_Print: unknown potential_type "//this%potential_type)
  endif

end subroutine IPModel_TTM_nF_Print

subroutine IPModel_TTM_nF_read_params_xml(this, param_str)
  type(IPModel_TTM_nF), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type(xml_t) :: fxml

  if (len(trim(param_str)) <= 0) return

  parse_in_ip = .false.
  parse_matched_label = .false.
  parse_ip => this

  call open_xml_string(fxml, param_str)
  call parse(fxml,  &
    startElement_handler = IPModel_startElement_handler, &
    endElement_handler = IPModel_endElement_handler)
  call close_xml_t(fxml)

end subroutine IPModel_TTM_nF_read_params_xml

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=STRING_LENGTH) :: value

  if (name == 'TTM_nF_params') then ! new TTM_nF stanza

    if (parse_matched_label) return ! we already found an exact match for this label

    call QUIP_FoX_get_value(attributes, 'label', value, status)
    if (status /= 0) value = ''

    if (len(trim(parse_ip%label)) > 0) then ! we were passed in a label
      if (value == parse_ip%label) then ! exact match
        parse_matched_label = .true.
        parse_in_ip = .true.
      else ! no match
        parse_in_ip = .false.
      endif
    else ! no label passed in
      parse_in_ip = .true.
    endif

    if(parse_in_ip) then
       call finalise(parse_ip)
    endif
  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'TTM_nF_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_TTM_nF_module
