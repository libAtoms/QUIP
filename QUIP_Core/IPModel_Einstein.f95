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
!X IPModel_Einstein
!X
!% Einstein module
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_Einstein_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_Einstein
type IPModel_Einstein
  integer :: n_types = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)
  real(dp), dimension(:), allocatable :: spring_constant

  type(atoms) :: ref
  character(len=FIELD_LENGTH) :: ref_file

  real(dp) :: cutoff = 0.0_dp

  character(len=FIELD_LENGTH) :: label

end type IPModel_Einstein

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_Einstein), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_Einstein_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Einstein_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Einstein_Print
end interface Print

interface Calc
  module procedure IPModel_Einstein_Calc
end interface Calc

contains

subroutine IPModel_Einstein_Initialise_str(this, args_str, param_str)
  type(IPModel_Einstein), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_Einstein_Initialise_str args_str')) then
    call system_abort("IPModel_Einstein_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_Einstein_read_params_xml(this, param_str)

  call read(this%ref,trim(this%ref_file))
  !  Add initialisation code here

end subroutine IPModel_Einstein_Initialise_str

subroutine IPModel_Einstein_Finalise(this)
  type(IPModel_Einstein), intent(inout) :: this

  ! Add finalisation code here

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%spring_constant)) deallocate(this%spring_constant)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_Einstein_Finalise


subroutine IPModel_Einstein_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_Einstein), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   real(dp), dimension(3) :: diff
   real(dp), dimension(3,3) :: virial_i
   real(dp) :: ki, ei
   integer :: i, ti
   ! Add calc() code here

   type(Dictionary) :: params
   logical, dimension(:), pointer :: atom_mask_pointer
   logical :: has_atom_mask_name
   character(FIELD_LENGTH) :: atom_mask_name
   real(dp) :: r_scale, E_scale
   logical :: do_rescale_r, do_rescale_E

   INIT_ERROR(error)

   if (present(e)) e = 0.0_dp
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_Einstein_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_Einstein_Calc', error)
      f = 0.0_dp
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_Einstein_Calc', error)
      local_virial = 0.0_dp
      RAISE_ERROR("IPModel_Einstein_Calc: local_virial calculation requested but not supported yet.", error)
   endif

   atom_mask_pointer => null()
   if(present(args_str)) then
      call initialise(params)
      call param_register(params, 'atom_mask_name', 'NONE',atom_mask_name,has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
      call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")

      if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='IPModel_Einstein_Calc args_str')) &
      call system_abort("IPModel_Einstein_Calc failed to parse args_str='"//trim(args_str)//"'")
      call finalise(params)
 
 
      if( has_atom_mask_name ) then
         if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) &
         call system_abort("IPModel_Einstein_Calc did not find "//trim(atom_mask_name)//" propery in the atoms object.")
      else
         atom_mask_pointer => null()
      endif
      if (do_rescale_r .or. do_rescale_E) then
         RAISE_ERROR("IPModel_Einstein_Calc: rescaling of potential with r_scale and E_scale not yet implemented!", error)
      end if
   endif


   if(at%N /= this%ref%N) call system_abort("IPModel_Einstein_Calc: number of atoms "//at%N//" does not match the number of atoms "//this%ref%N//" in the &
   reference structure")

   do i = 1, at%N

      if(associated(atom_mask_pointer)) then
         if(.not. atom_mask_pointer(i)) cycle
      endif

      ti = get_type(this%type_of_atomic_num, at%Z(i))
      ki = this%spring_constant(ti)

      diff = diff_min_image(at,this%ref%pos(:,i),at%pos(:,i))
      if(present(local_e) .or. present(e)) ei = 0.5_dp * ki * sum(diff**2)

      if(present(local_e)) local_e(i) = ei
      if(present(e)) e = e + ei

      if(present(f)) f(:,i) = f(:,i) - ki*diff
      if(present(virial) .or. present(local_virial)) virial_i = ki*(diff .outer. at%pos(:,i))

      if(present(local_virial)) local_virial(:,i) = local_virial(:,i) - reshape(virial_i,(/9/))
      if(present(virial)) virial = virial - virial_i
   enddo

end subroutine IPModel_Einstein_Calc


subroutine IPModel_Einstein_Print(this, file)
  type(IPModel_Einstein), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti

  call Print("IPModel_Einstein : Einstein Potential", file=file)
  call Print("IPModel_Einstein : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_Einstein : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call Print ("IPModel_Einstein : type " // ti // " spring_constant " // this%spring_constant(ti), file=file)
    call verbosity_push_decrement()
    call Print ("IPModel_Einstein : " // &
        "cutoff " // this%cutoff , &
         file=file)
    call verbosity_pop()
  end do

end subroutine IPModel_Einstein_Print

subroutine IPModel_Einstein_read_params_xml(this, param_str)
  type(IPModel_Einstein), intent(inout), target :: this
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
    call system_abort("IPModel_Einstein_read_params_xml parsed file, but n_types = 0")
  endif
end subroutine IPModel_Einstein_read_params_xml

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

  if (name == 'Einstein_params') then ! new Einstein stanza

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
        call system_abort("Can't find n_types in Einstein_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0
      allocate(parse_ip%spring_constant(parse_ip%n_types))
      parse_ip%spring_constant = 0.0_dp

      call QUIP_FoX_get_value(attributes, "cutoff", value, status)
      if (status /= 0) call system_abort ("IPModel_Einstein_read_params_xml cannot find cutoff")
      read (value, *) parse_ip%cutoff

      call QUIP_FoX_get_value(attributes, "ref_file", value, status)
      if (status /= 0) call system_abort ("IPModel_Einstein_read_params_xml cannot find ref_file")
      parse_ip%ref_file = trim(value)
    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_Einstein_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_Einstein_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    call QUIP_FoX_get_value(attributes, "spring_constant", value, status)
    if (status /= 0) call system_abort ("IPModel_Einstein_read_params_xml cannot find spring_constant")
    read (value, *) parse_ip%spring_constant(ti)

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
    if (name == 'Einstein_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_Einstein_module
