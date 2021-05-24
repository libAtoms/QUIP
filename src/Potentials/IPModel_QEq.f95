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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X IPModel_Qeq
!X
!% Charge Equlibration module. Calculates electrostatic interactions between charged
!% species. Supported methods:
!% Ewald: Ewald summation technique
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module IPModel_QEq_module

use error_module
use system_module, only : dp, inoutput, print, lower_case, verbosity_push_decrement, verbosity_pop, operator(//), split_string, string_to_int
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module
use units_module

use mpi_context_module
use QUIP_Common_module

use IPEwald_module


implicit none
private

include 'IPModel_interface.h'

public :: IPModel_Qeq
type IPModel_Qeq

  integer :: n_types = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)
  real(dp), dimension(:), allocatable :: en, hardness

  real(dp) :: cutoff = 0.0_dp

  real(dp) :: ewald_error
  real(dp) :: smooth_QEq_cutoff

  character(len=STRING_LENGTH) :: label


endtype IPModel_QEq

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_QEq), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_QEq_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_QEq_Finalise
end interface Finalise

interface Print
  module procedure IPModel_QEq_Print
end interface Print

interface Calc
  module procedure IPModel_QEq_Calc
end interface Calc

contains

subroutine IPModel_QEq_Initialise_str(this, args_str, param_str)
  type(IPModel_QEq), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params
  character(len=STRING_LENGTH) :: method_str
  logical :: has_method

  call Finalise(this)
  call initialise(params)
  this%label=''
  method_str=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_QEq_Initialise_str args_str')) then
    call system_abort("IPModel_QEq_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_QEq_read_params_xml(this, param_str)

  !  Add initialisation code here

end subroutine IPModel_QEq_Initialise_str

subroutine IPModel_QEq_Finalise(this)
  type(IPModel_QEq), intent(inout) :: this

  ! Add finalisation code here

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%en)) deallocate(this%en)
  if (allocated(this%hardness)) deallocate(this%hardness)

  this%n_types = 0
  this%label = ''
  this%cutoff = 0.0_dp

end subroutine IPModel_QEq_Finalise

recursive subroutine IPModel_QEq_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   use minimization_module

   type(IPModel_QEq), intent(inout):: this
   type(Atoms), intent(inout)  :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   type(Dictionary) :: params
   real(dp), dimension(:), allocatable, target :: my_en, my_hardness, x, b
   real(dp), dimension(:,:), allocatable :: J_matrix
   real(dp), dimension(:), pointer :: en, hardness, charge
   real(dp) :: r_scale, E_scale
   logical :: do_rescale_r, do_rescale_E, real_space, reciprocal_space

   integer :: i, n_it
   character(len=STRING_LENGTH) :: charge_property_name, en_property_name, hardness_property_name, atom_mask_name, source_mask_name
   type(LA_Matrix) :: LA_J_matrix

   INIT_ERROR(error)

   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_QEq_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_QEq_Calc', error)
      f = 0.0_dp
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_QEq_Calc', error)
      local_virial = 0.0_dp
   endif

   if (present(args_str)) then
      call initialise(params)
      call param_register(params, 'en_property_name', 'en', en_property_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'hardness_property_name', 'hardness', hardness_property_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'charge_property_name', 'charge', charge_property_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'atom_mask_name', '', atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'source_mask_name', '', source_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Rescaling factor for distances. Default 1.0.")
      call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Rescaling factor for energy. Default 1.0.")
      call param_register(params, 'real_space', 'T',real_space, help_string="Do real-space contribution only")
      call param_register(params, 'reciprocal_space', 'T',reciprocal_space, help_string="Do reciprocal-space contribution only")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_QEq_Calc args_str')) then
         RAISE_ERROR("IPModel_QEq_Calc failed to parse args_str="//trim(args_str), error)
      endif
      call finalise(params)
      if (do_rescale_r ) then
         RAISE_ERROR("IPModel_QEq_Calc: rescaling of potential with r_scale not yet implemented!", error)
      end if
   else
      en_property_name = 'en'
      hardness_property_name = 'hardness'
      charge_property_name = 'charge'
   endif

   if(has_property(at,en_property_name)) then
      if(.not. assign_pointer(at, en_property_name, en)) then
         RAISE_ERROR('IPModel_QEq_Calc failed to assign pointer to '//trim(en_property_name)//' property', error)
      endif
   else
      allocate(my_en(at%N))
      en => my_en
      en = 0.0_dp
      do i = 1, at%N
         en(i) = this%en(this%type_of_atomic_num(at%Z(i))) 
      enddo
   endif

   if(has_property(at,hardness_property_name)) then
      if(.not. assign_pointer(at, hardness_property_name, en)) then
         RAISE_ERROR('IPModel_QEq_Calc failed to assign pointer to '//trim(hardness_property_name)//' property', error)
      endif
   else
      allocate(my_hardness(at%N))
      hardness => my_hardness
      hardness = 0.0_dp
      do i = 1, at%N
         hardness(i) = this%hardness(this%type_of_atomic_num(at%Z(i))) 
      enddo
   endif

   if(has_property(at,charge_property_name)) then
      if(.not. assign_pointer(at, charge_property_name, charge)) then
         RAISE_ERROR('IPModel_QEq_Calc failed to assign pointer to '//trim(charge_property_name)//' property', error)
      endif
   else
      call add_property(at,charge_property_name,0.0_dp,n_cols=1,ptr=charge,error=error)
      charge = [-1.0_dp,0.25_dp,0.25_dp,0.25_dp,0.25_dp]
   endif

   allocate(x(at%N+1),b(at%N+1),J_matrix(at%N+1,at%N+1))
   !call Direct_Coulomb_calc(at, charge, e=e, d2e_dq2 = J_matrix(1:at%N,1:at%N) )
   call Ewald_calc(at, charge, ewald_error=this%ewald_error, use_ewald_cutoff=.false., smooth_coulomb_cutoff=this%smooth_qeq_cutoff, &
   d2e_dq2 = J_matrix(1:at%N,1:at%N), real_space = real_space, reciprocal_space=reciprocal_space, error=error )

   J_matrix(at%N+1,1:at%N) = -1.0_dp 
   J_matrix(1:at%N,at%N+1) = -1.0_dp 
   J_matrix(at%N+1,at%N+1) = 0.0_dp
   b(at%N+1) = 0.0_dp
   do i = 1, at%N
      b(i) = -en(i)
      J_matrix(i,i) = J_matrix(i,i) + hardness(i)
   enddo

   call initialise(LA_J_matrix,J_matrix)
   call matrix_qr_solve(LA_J_matrix,b,x)
   charge = x(1:at%N)

   call Ewald_calc(at,charge,e=e,f=f,virial=virial,ewald_error=this%ewald_error, use_ewald_cutoff=.false., smooth_coulomb_cutoff=this%smooth_qeq_cutoff,error=error)
   !if(present(e)) then
   !   e = 0.5_dp*dot_product(charge,matmul(J_matrix(1:at%N,1:at%N),charge)) + dot_product(charge,en)
   !endif
   !x(1:at%N) = charge
   !x(at%N+1) = 1.0_dp

   !n_it = minim(x,QEq_e,dQEq_e,"cg",1.0e-6_dp,100)
   !charge = x(1:at%N)

   !call print("charge "//charge)
   charge => null()
   en => null()
   hardness => null()
   if(allocated(my_en)) deallocate(my_en)
   if(allocated(my_hardness)) deallocate(my_hardness)

   if (do_rescale_E) then
      if (present(e)) e = e*E_scale
      if (present(local_e)) local_e = local_e*E_scale
      if (present(f)) f = f*E_scale
      if (present(virial)) virial=virial*E_scale
      if (present(local_virial)) local_virial=local_virial*E_scale
   end if

   contains

      function QEq_e(q,data)
         real(dp), dimension(:) :: q
         character(len=1),optional::data(:)
         real(dp) :: QEq_e

         !call Ewald_calc(at, q, e=e, ewald_error=this%ewald_error, use_ewald_cutoff=.false., smooth_coulomb_cutoff=this%smooth_coulomb_cutoff, de_dq=de_dq,error=error)
         !call Ewald_calc(at, q(1:at%N), e=e, ewald_error=this%ewald_error, use_ewald_cutoff=.false., smooth_coulomb_cutoff=this%smooth_qeq_cutoff)
         call Direct_Coulomb_calc(at, q(1:at%N), e=e )

         !QEq_e = e + dot_product(en,q(1:at%N)) + 0.5_dp*dot_product(hardness,q(1:at%N)**2) + q(at%N+1)**2 * sum(q(1:at%N))**2
         QEq_e = e + dot_product(en,q(1:at%N)) + 0.5_dp*dot_product(hardness,q(1:at%N)**2) + q(at%N+1) * sum(q(1:at%N))
         print*,"e",e

      endfunction QEq_e

      function dQEq_e(q,data)
         real(dp), dimension(:) :: q
         character(len=1),optional::data(:)
         real(dp), dimension(size(q)) :: dQEq_e
         real(dp), dimension(:), allocatable :: de_dq

         !call Ewald_calc(at, q, e=e, ewald_error=this%ewald_error, use_ewald_cutoff=.false., smooth_coulomb_cutoff=this%smooth_coulomb_cutoff, de_dq=de_dq,error=error)
         allocate(de_dq(at%N))
         !call Ewald_calc(at, q(1:at%N), ewald_error=this%ewald_error, use_ewald_cutoff=.false., de_dq=de_dq, smooth_coulomb_cutoff=this%smooth_qeq_cutoff)
         call Direct_Coulomb_calc(at, q(1:at%N), de_dq=de_dq )
         !dQEq_e(1:at%N) = de_dq + en + hardness*q(1:at%N) + 2*q(at%N+1)**2 * sum(q(1:at%N))
         !dQEq_e(at%N+1) = 2*q(at%N+1)*sum(q(1:at%N))**2
         dQEq_e(1:at%N) = de_dq + en + hardness*q(1:at%N) + q(at%N+1)
         dQEq_e(at%N+1) = sum(q(1:at%N))
         deallocate(de_dq)

      endfunction dQEq_e

end subroutine IPModel_QEq_Calc


subroutine IPModel_QEq_Print(this, file)
  type(IPModel_QEq), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti , tj

  call Print("IPModel_QEq : QEq Potential", file=file)
  call Print("IPModel_QEq : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_QEq : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    call Print ("IPModel_QEq : " // &
        "cutoff " // this%cutoff , &
         file=file)
    call verbosity_pop()
  end do

end subroutine IPModel_QEq_Print

subroutine IPModel_QEq_read_params_xml(this, param_str)
  type(IPModel_QEq), intent(inout), target :: this
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

  if (this%n_types == 0 ) then
    call system_abort("IPModel_QEq_read_params_xml parsed file, but n_types = 0 ")
  endif

end subroutine IPModel_QEq_read_params_xml

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

  character(len=STRING_LENGTH) :: value, signature_string, pos_type_str
  character(len=STRING_LENGTH), dimension(99) :: signature_fields

  logical :: energy_shift, linear_force_shift
  integer :: ti, tj, i ,n_atoms
  real(dp) :: alpha_au

  if (name == 'QEq_params') then ! new QEq stanza

    if (parse_matched_label) return ! we already found an exact match for this label

    call QUIP_FoX_get_value(attributes, 'label', value, status)
    if (status /= 0) value = ''

    if (len(trim(parse_ip%label)) > 0) then ! we were passed in a label
      if (trim(value) == trim(parse_ip%label)) then ! exact match
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
        read (value, *) parse_ip%n_types
      else
        call system_abort("Can't find n_types in QEq_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      allocate(parse_ip%en(parse_ip%n_types))
      parse_ip%en = 0.0_dp

      allocate(parse_ip%hardness(parse_ip%n_types))
      parse_ip%hardness = 0.0_dp

      call QUIP_FoX_get_value(attributes, "cutoff", value, status)
      if (status /= 0) call system_abort ("IPModel_QEq_read_params_xml cannot find cutoff")
      read (value, *) parse_ip%cutoff

      call QUIP_FoX_get_value(attributes, "ewald_error", value, status)
      if (status == 0) then
         read (value, *) parse_ip%ewald_error
      else
         parse_ip%ewald_error = 1.0e-6_dp
      endif

      call QUIP_FoX_get_value(attributes, "smooth_QEq_cutoff", value, status)
      if (status == 0) then
         read (value, *) parse_ip%smooth_QEq_cutoff
      else
         parse_ip%smooth_QEq_cutoff = 0.0_dp
      endif

    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_QEq_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_QEq_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    call QUIP_FoX_get_value(attributes, "en", value, status)
    if (status /= 0) call system_abort ("IPModel_QEq_read_params_xml cannot find electronegativity")
    read (value, *) parse_ip%en(ti)

    call QUIP_FoX_get_value(attributes, "hardness", value, status)
    if (status /= 0) call system_abort ("IPModel_QEq_read_params_xml cannot find hardness")
    read (value, *) parse_ip%hardness(ti)

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
    if (name == 'QEq_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_QEq_module
