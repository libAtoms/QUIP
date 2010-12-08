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
!X IPModel_Sutton_Chen
!X
!% Sutton_Chen module for use when implementing a new Interatomic Potential
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_Sutton_Chen_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_Sutton_Chen
type IPModel_Sutton_Chen
  integer :: n_types = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)

  real(dp), dimension(:), allocatable :: a, epsilon, c, v_cutoff, rho_cutoff, &
  dv_cutoff, drho_cutoff
  integer, dimension(:), allocatable  :: m, n
  real(dp) :: cutoff = 0.0_dp

  character(len=FIELD_LENGTH) :: label

end type IPModel_Sutton_Chen

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_Sutton_Chen), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_Sutton_Chen_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Sutton_Chen_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Sutton_Chen_Print
end interface Print

interface Calc
  module procedure IPModel_Sutton_Chen_Calc
end interface Calc

contains

subroutine IPModel_Sutton_Chen_Initialise_str(this, args_str, param_str)
  type(IPModel_Sutton_Chen), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  integer :: ti

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_Sutton_Chen_Initialise_str args_str')) then
    call system_abort("IPModel_Sutton_Chen_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_Sutton_Chen_read_params_xml(this, param_str)

  do ti = 1, this%n_types
     this%v_cutoff(ti) = ( this%a(ti) / this%cutoff )**this%n(ti)
     this%dv_cutoff(ti) = - this%n(ti) * ( this%a(ti) / this%cutoff )**this%n(ti) / this%cutoff
     this%rho_cutoff(ti) = ( this%a(ti) / this%cutoff )**this%m(ti)
     this%drho_cutoff(ti) = - this%m(ti) * ( this%a(ti) / this%cutoff )**this%m(ti) / this%cutoff
  enddo

end subroutine IPModel_Sutton_Chen_Initialise_str

subroutine IPModel_Sutton_Chen_Finalise(this)
  type(IPModel_Sutton_Chen), intent(inout) :: this

  ! Add finalisation code here

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%a)) deallocate(this%a)
  if (allocated(this%epsilon)) deallocate(this%epsilon)
  if (allocated(this%c)) deallocate(this%c)
  if (allocated(this%m)) deallocate(this%m)
  if (allocated(this%n)) deallocate(this%n)
  if (allocated(this%v_cutoff)) deallocate(this%v_cutoff)
  if (allocated(this%dv_cutoff)) deallocate(this%dv_cutoff)
  if (allocated(this%rho_cutoff)) deallocate(this%rho_cutoff)
  if (allocated(this%drho_cutoff)) deallocate(this%drho_cutoff)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_Sutton_Chen_Finalise


subroutine IPModel_Sutton_Chen_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_Sutton_Chen), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   integer :: i, j, n, ti, tj
   real(dp) :: r_ij, rho_i, v_i, a_over_r, drho_i_drij, dv_i_drij, df_i, e_i
   real(dp), dimension(3) :: drho_i_dri, dv_i_dri, u_ij
   real(dp), dimension(3,3) :: drho_i_drij_outer_rij, dv_i_drij_outer_rij, virial_i

   type(Dictionary) :: params
   logical, dimension(:), pointer :: atom_mask_pointer
   logical :: has_atom_mask_name
   character(FIELD_LENGTH) :: atom_mask_name
   real(dp) :: r_scale, E_scale
   logical :: do_rescale_r, do_rescale_E


   INIT_ERROR(error)

   if (present(e)) e = 0.0_dp
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_Sutton_Chen_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_Sutton_Chen_Calc', error)
      f = 0.0_dp
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_Sutton_Chen_Calc', error)
      local_virial = 0.0_dp
   endif

   atom_mask_pointer => null()
   if(present(args_str)) then
      call initialise(params)
      call param_register(params, 'atom_mask_name', 'NONE',atom_mask_name,has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
      call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")
      if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='IPModel_Sutton_Chen_Calc args_str')) then
         RAISE_ERROR("IPModel_Sutton_Chen_Calc failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)
 
 
      if( has_atom_mask_name ) then
         if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) then
            RAISE_ERROR("IPModel_Sutton_Chen_Calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
         endif
      else
         atom_mask_pointer => null()
      endif
      if (do_rescale_r .or. do_rescale_E) then
         RAISE_ERROR("IPModel_Sutton_Chen_Calc: rescaling of potential with r_scale and E_scale not yet implemented!", error)
      end if
   endif

   do i = 1, at%N
      if(associated(atom_mask_pointer)) then
         if(.not. atom_mask_pointer(i)) cycle
      endif

      rho_i = 0.0_dp ! Local density from neighbours
      v_i   = 0.0_dp

      drho_i_dri = 0.0_dp
      dv_i_dri   = 0.0_dp

      drho_i_drij_outer_rij = 0.0_dp
      dv_i_drij_outer_rij   = 0.0_dp

      ! Get the type of this species from its atomic number
      ti = get_type(this%type_of_atomic_num, at%Z(i))

      do n = 1, atoms_n_neighbours(at, i)
         j = atoms_neighbour(at, i, n, distance=r_ij, cosines=u_ij, max_dist=this%cutoff)
         if( j <= 0 ) cycle
         tj = get_type(this%type_of_atomic_num, at%Z(j))
         a_over_r = this%a(tj) / r_ij 

         rho_i = rho_i + ( a_over_r**this%m(tj) - this%rho_cutoff(tj) )
         v_i = v_i + ( a_over_r**this%n(tj) - this%v_cutoff(tj) )

         if( present(f) .or. present(virial) .or. present(local_virial) ) then
            drho_i_drij = - this%m(tj) * a_over_r**this%m(tj) / r_ij
            dv_i_drij = -this%n(tj) * a_over_r**this%n(tj) / r_ij

            if( present(f) ) then
               drho_i_dri = drho_i_dri + drho_i_drij * u_ij
               dv_i_dri = dv_i_dri + dv_i_drij * u_ij
            endif

            if( present(virial) .or. present(local_virial) ) then
               drho_i_drij_outer_rij = drho_i_drij_outer_rij + drho_i_drij*(u_ij .outer. u_ij)*r_ij
               dv_i_drij_outer_rij = dv_i_drij_outer_rij + dv_i_drij*(u_ij .outer. u_ij)*r_ij
            endif
         endif
      enddo

      df_i = - 0.5_dp * this%c(ti) / sqrt(rho_i)

      if(present(f)) then
         f(:,i) = f(:,i) + this%epsilon(ti) * ( drho_i_dri * df_i + dv_i_dri )
         do n = 1, atoms_n_neighbours(at, i)
            j = atoms_neighbour(at, i, n, distance=r_ij, cosines=u_ij, max_dist=this%cutoff)
            if( j <= 0 ) cycle

            tj = get_type(this%type_of_atomic_num, at%Z(j))

            a_over_r = this%a(tj) / r_ij
            drho_i_drij = - this%m(tj) * a_over_r**this%m(tj) / r_ij

            f(:,j) = f(:,j) - this%epsilon(tj) * df_i * drho_i_drij * u_ij
         enddo
      endif

      if( present(e) .or. present(local_e) ) then
         e_i = this%epsilon(ti) * ( 0.5_dp * v_i - this%c(ti) * sqrt(rho_i) )
         if( present(e) ) e = e + e_i
         if( present(local_e) ) local_e(i) = e_i
      endif

      if( present(virial) .or. present(local_virial) ) then
         virial_i = this%epsilon(ti) * ( - 0.5_dp * dv_i_drij_outer_rij - &
         df_i * drho_i_drij_outer_rij )

         if(present(virial)) virial = virial + virial_i
         if(present(local_virial)) local_virial(:,i) = local_virial(:,i) + reshape(virial_i,(/9/))
      endif
   enddo ! i

   atom_mask_pointer => null()

end subroutine IPModel_Sutton_Chen_Calc


subroutine IPModel_Sutton_Chen_Print(this, file)
  type(IPModel_Sutton_Chen), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti

  call Print("IPModel_Sutton_Chen : Sutton_Chen Potential", file=file)
  call Print("IPModel_Sutton_Chen : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_Sutton_Chen : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    call Print ("IPModel_Sutton_Chen : " // &
        "cutoff " // this%cutoff , &
         file=file)
    call verbosity_pop()
  end do

end subroutine IPModel_Sutton_Chen_Print

subroutine IPModel_Sutton_Chen_read_params_xml(this, param_str)
  type(IPModel_Sutton_Chen), intent(inout), target :: this
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
    call system_abort("IPModel_Sutton_Chen_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_Sutton_Chen_read_params_xml

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

  if (name == 'Sutton_Chen_params') then ! new Sutton_Chen stanza

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
        call system_abort("Can't find n_types in Sutton_Chen_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0
      allocate(parse_ip%a(parse_ip%n_types))
      parse_ip%a = 0.0_dp
      allocate(parse_ip%epsilon(parse_ip%n_types))
      parse_ip%epsilon = 0.0_dp
      allocate(parse_ip%c(parse_ip%n_types))
      parse_ip%c = 0.0_dp
      allocate(parse_ip%m(parse_ip%n_types))
      parse_ip%m = 0
      allocate(parse_ip%n(parse_ip%n_types))
      parse_ip%n = 0
      allocate(parse_ip%v_cutoff(parse_ip%n_types))
      parse_ip%v_cutoff = 0.0_dp
      allocate(parse_ip%dv_cutoff(parse_ip%n_types))
      parse_ip%dv_cutoff = 0.0_dp
      allocate(parse_ip%rho_cutoff(parse_ip%n_types))
      parse_ip%rho_cutoff = 0.0_dp
      allocate(parse_ip%drho_cutoff(parse_ip%n_types))
      parse_ip%drho_cutoff = 0.0_dp

      call QUIP_FoX_get_value(attributes, "cutoff", value, status)
      if (status /= 0) call system_abort ("IPModel_Sutton_Chen_read_params_xml cannot find cutoff")
      read (value, *) parse_ip%cutoff
    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_Sutton_Chen_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_Sutton_Chen_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    call QUIP_FoX_get_value(attributes, "a", value, status)
    if (status /= 0) call system_abort ("IPModel_Sutton_Chen_read_params_xml cannot find a")
    read (value, *) parse_ip%a(ti)

    call QUIP_FoX_get_value(attributes, "epsilon", value, status)
    if (status /= 0) call system_abort ("IPModel_Sutton_Chen_read_params_xml cannot find epsilon")
    read (value, *) parse_ip%epsilon(ti)

    call QUIP_FoX_get_value(attributes, "c", value, status)
    if (status /= 0) call system_abort ("IPModel_Sutton_Chen_read_params_xml cannot find c")
    read (value, *) parse_ip%c(ti)

    call QUIP_FoX_get_value(attributes, "m", value, status)
    if (status /= 0) call system_abort ("IPModel_Sutton_Chen_read_params_xml cannot find m")
    read (value, *) parse_ip%m(ti)

    call QUIP_FoX_get_value(attributes, "n", value, status)
    if (status /= 0) call system_abort ("IPModel_Sutton_Chen_read_params_xml cannot find n")
    read (value, *) parse_ip%n(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
        parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    enddo

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'Sutton_Chen_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_Sutton_Chen_module
