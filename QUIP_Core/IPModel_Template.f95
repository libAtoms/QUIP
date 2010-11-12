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
!X IPModel_Template
!X
!% Template module for use when implementing a new Interatomic Potential
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_Template_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_Template
type IPModel_Template
  integer :: n_types = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)

  real(dp) :: cutoff = 0.0_dp

  character(len=FIELD_LENGTH) :: label

end type IPModel_Template

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_Template), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_Template_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Template_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Template_Print
end interface Print

interface Calc
  module procedure IPModel_Template_Calc
end interface Calc

contains

subroutine IPModel_Template_Initialise_str(this, args_str, param_str)
  type(IPModel_Template), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_Template_Initialise_str args_str')) then
    call system_abort("IPModel_Template_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_Template_read_params_xml(this, param_str)

  !  Add initialisation code here

end subroutine IPModel_Template_Initialise_str

subroutine IPModel_Template_Finalise(this)
  type(IPModel_Template), intent(inout) :: this

  ! Add finalisation code here

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_Template_Finalise


subroutine IPModel_Template_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_Template), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   ! Add calc() code here

   INIT_ERROR(error)

   if (present(e)) e = 0.0_dp
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_Template_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_Template_Calc', error)
      f = 0.0_dp
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_Template_Calc', error)
      local_virial = 0.0_dp
   endif

   RAISE_ERROR('IPModel_Calc - not implemented',error)

end subroutine IPModel_Template_Calc


subroutine IPModel_Template_Print(this, file)
  type(IPModel_Template), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti

  call Print("IPModel_Template : Template Potential", file=file)
  call Print("IPModel_Template : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_Template : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    call Print ("IPModel_Template : " // &
        "cutoff " // this%cutoff , &
         file=file)
    call verbosity_pop()
  end do

end subroutine IPModel_Template_Print

subroutine IPModel_Template_read_params_xml(this, param_str)
  type(IPModel_Template), intent(inout), target :: this
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

  if (this%n_types == 0) then
    call system_abort("IPModel_Template_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_Template_read_params_xml

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
  character(len=FIELD_LENGTH) :: value

  integer ti

  if (name == 'Template_params') then ! new Template stanza

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

    if (parse_in_ip) then
      if (parse_ip%n_types /= 0) then
        call finalise(parse_ip)
      endif

      call QUIP_FoX_get_value(attributes, 'n_types', value, status)
      if (status == 0) then
        read (value, *), parse_ip%n_types
      else
        call system_abort("Can't find n_types in Template_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      call QUIP_FoX_get_value(attributes, "cutoff", value, status)
      if (status /= 0) call system_abort ("IPModel_Template_read_params_xml cannot find cutoff")
      read (value, *) parse_ip%cutoff
    endif


  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_Template_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_Template_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
        parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do


  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'Template_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_Template_module
