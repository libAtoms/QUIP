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
!X IPModel_FB module  
!%
!% Module for Flikkema \& Bromley interatomic potential
!%
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_FB_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none

private 

include 'IPModel_interface.h'

public :: IPModel_FB
type IPModel_FB
  integer :: n_types = 0         !% Number of atomic types. 
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)  !% Atomic number dimensioned as \texttt{n_types}. 

  real(dp) :: cutoff = 0.0_dp    !% Cutoff for computing connection.

  real(dp), dimension(:,:), allocatable :: A, B, C

  character(len=STRING_LENGTH) label

end type IPModel_FB

logical, private :: parse_in_ip, parse_matched_label
integer :: parse_cur_type_i, parse_cur_type_j
type(IPModel_FB), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_FB_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_FB_Finalise
end interface Finalise

interface Print
  module procedure IPModel_FB_Print
end interface Print

interface Calc
  module procedure IPModel_FB_Calc
end interface Calc

contains

subroutine IPModel_FB_Initialise_str(this, args_str, param_str)
  type(IPModel_FB), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label = ''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_FB_Initialise_str args_str')) then
    call system_abort("IPModel_FB_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_FB_read_params_xml(this, param_str)

end subroutine IPModel_FB_Initialise_str

subroutine IPModel_FB_Finalise(this)
  type(IPModel_FB), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if (allocated(this%A)) deallocate(this%A)
  if (allocated(this%B)) deallocate(this%B)
  if (allocated(this%C)) deallocate(this%C)

  this%n_types = 0
  this%cutoff = 0.0_dp
  this%label = ''
end subroutine IPModel_FB_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator: this routine computes energy, forces and the virial.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_FB_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_FB), intent(in) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(out), optional :: e !% \texttt{e} = System total energy
  real(dp), dimension(:), intent(out), optional :: local_e !% \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), dimension(3,3), intent(out), optional :: virial   !% Virial
  character(len=*), optional      :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  integer :: i, j, n, ti, tj
  real(dp) :: r_ij, de, de_dr
  real(dp), dimension(3) :: u_ij
  real(dp), dimension(3,3) :: virial_i

  type(Dictionary) :: params
  logical, dimension(:), pointer :: atom_mask_pointer
  logical :: has_atom_mask_name
  character(STRING_LENGTH) :: atom_mask_name
  real(dp) :: r_scale, E_scale
  logical :: do_rescale_r, do_rescale_E

  type(Atoms) :: at_ew

  INIT_ERROR(error)

  if (present(e)) e = 0.0_dp
  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_FB_Calc', error)
     local_e = 0.0_dp
  endif
  if (present(f)) then
     call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_FB_Calc', error)
     f = 0.0_dp
  end if
  if (present(virial)) virial = 0.0_dp
  if (present(local_virial)) then
     call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_FB_Calc', error)
     local_virial = 0.0_dp
  endif

  atom_mask_pointer => null()
  if(present(args_str)) then
     call initialise(params)
     call param_register(params, 'atom_mask_name', 'NONE',atom_mask_name,has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
     call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")

     if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='IPModel_FB_Calc args_str')) then
        RAISE_ERROR("IPModel_FB_Calc failed to parse args_str='"//trim(args_str)//"'",error)
     endif
     call finalise(params)


     if( has_atom_mask_name ) then
        if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) then
           RAISE_ERROR("IPModel_FB_Calc did not find "//trim(atom_mask_name)//" property in the atoms object.",error)
        endif
     else
        atom_mask_pointer => null()
     endif
     if (do_rescale_r .or. do_rescale_E) then
        RAISE_ERROR("IPModel_FB_Calc: rescaling of potential with r_scale and E_scale not yet implemented!", error)
     end if
  endif

  do i = 1, at%N
    if (present(mpi)) then
       if (mpi%active) then
	 if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
       endif
    endif

    if(associated(atom_mask_pointer)) then
       if(.not. atom_mask_pointer(i)) cycle
    endif

    ti = get_type(this%type_of_atomic_num, at%Z(i))
    do n = 1, atoms_n_neighbours(at, i)
      j = atoms_neighbour(at, i, n, distance=r_ij, cosines = u_ij)
      if (r_ij .feq. 0.0_dp) cycle
      tj = get_type(this%type_of_atomic_num, at%Z(j))
      if (present(e) .or. present(local_e)) then
	de = 0.5_dp * (this%A(ti,tj) * exp( - r_ij / this%B(ti,tj) ) - this%C(ti,tj) / r_ij**6)
	if (present(local_e)) local_e(i) = local_e(i) + de
	if (present(e)) e = e + de
      endif
      if (present(f) .or. present(virial)) then
	de_dr = 0.5_dp * (-this%A(ti,tj) * exp( - r_ij / this%B(ti,tj) ) / this%B(ti,tj) + &
              & 6.0_dp * this%C(ti,tj) / r_ij**7 )
	if (present(f)) then
	  f(:,i) = f(:,i) + de_dr*u_ij
	  f(:,j) = f(:,j) - de_dr*u_ij
	endif
	if (present(virial) .or. present(local_virial)) virial_i = de_dr*(u_ij .outer. u_ij)*r_ij
        if (present(virial))  virial = virial - virial_i
        if (present(local_virial)) local_virial(:,i) = local_virial(:,i) - reshape(virial_i,(/9/))
      endif
    end do
  end do

  if (present(mpi)) then
     if (present(e)) e = sum(mpi, e)
     if (present(local_e)) call sum_in_place(mpi, local_e)
     if (present(virial)) call sum_in_place(mpi, virial)
     if (present(f)) call sum_in_place(mpi, f)
  endif

end subroutine IPModel_FB_Calc

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions.
!% An example for XML stanza is given below.
!%
!%> 
!%> <FB_potential>
!%>    <!-- Flikkema & Bromley  -->
!%>    <Potential label="FB_Potential" init_args="Sum init_args_pot1={IP FB} init_args_pot2={IP Coulomb}"/>
!%> 
!%>    <Coulomb_params n_types="2" cutoff="6.0" method="ewald" label="default">
!%>       <per_type_data type="1" atomic_num="8" charge="-1.2"/>
!%>       <per_type_data type="2" atomic_num="14" charge="2.4"/>
!%>    </Coulomb_params>
!%> 
!%>    <FB_params n_types="2" label="default">
!%>       <per_type_data type="1" atomic_num="8"  />
!%>       <per_type_data type="2" atomic_num="14" />
!%>       <per_pair_data atomic_num_i="8" atomic_num_j="8" A="1428.406" B="0.358" C="41.374" r_cut="6.0" />
!%>       <per_pair_data atomic_num_i="8" atomic_num_j="14" A="10454.202" B="0.208" C="63.047" r_cut="6.0" />
!%>       <per_pair_data atomic_num_i="14" atomic_num_j="14" A="79502.113" B="0.201" C="446.780" r_cut="6.0" />
!%>    </FB_params>
!%> 
!%> </FB_potential>
!%> 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=1024) :: value

  integer :: ti, Zi, Zj

  if (name == 'FB_params') then ! new FB stanza

    if (parse_in_ip) &
      call system_abort("IPModel_startElement_handler entered FB_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

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
	call system_abort("Can't find n_types in FB_params")
      endif

      call QUIP_FoX_get_value(attributes, 'cutoff', value, status)
      if (status == 0) then
	read (value, *), parse_ip%cutoff
      else
	call system_abort("Can't find cutoff in FB_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0
      allocate(parse_ip%A(parse_ip%n_types,parse_ip%n_types))
      parse_ip%A = 0.0_dp
      allocate(parse_ip%B(parse_ip%n_types,parse_ip%n_types))
      parse_ip%B = 0.0_dp
      allocate(parse_ip%C(parse_ip%n_types,parse_ip%n_types))
      parse_ip%C = 0.0_dp
    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
	parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "atomic_num_i", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find atomic_num_i")
    read (value, *) Zi

    call QUIP_FoX_get_value(attributes, "atomic_num_j", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find atomic_num_i")
    read (value, *) Zj

    parse_cur_type_i = get_type(parse_ip%type_of_atomic_num,Zi)
    parse_cur_type_j = get_type(parse_ip%type_of_atomic_num,Zj)

    call QUIP_FoX_get_value(attributes, "A", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find A")
    read (value, *) parse_ip%A(parse_cur_type_i,parse_cur_type_j)

    call QUIP_FoX_get_value(attributes, "B", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find B")
    read (value, *) parse_ip%B(parse_cur_type_i,parse_cur_type_j)

    call QUIP_FoX_get_value(attributes, "C", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find C")
    read (value, *) parse_ip%C(parse_cur_type_i,parse_cur_type_j)

    if (parse_cur_type_i /= parse_cur_type_j) then
        parse_ip%A(parse_cur_type_j,parse_cur_type_i) = parse_ip%A(parse_cur_type_i,parse_cur_type_j)
        parse_ip%B(parse_cur_type_j,parse_cur_type_i) = parse_ip%B(parse_cur_type_i,parse_cur_type_j)
        parse_ip%C(parse_cur_type_j,parse_cur_type_i) = parse_ip%C(parse_cur_type_i,parse_cur_type_j)
    end if

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if (name == 'FB_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_FB_read_params_xml(this, param_str)
  type(IPModel_FB), intent(inout), target :: this
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
    call system_abort("IPModel_FB_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_FB_read_params_xml


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of FB parameters: number of different types, cutoff radius, atomic numbers, etc.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_FB_Print (this, file)
  type(IPModel_FB), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_FB : Flikkema Bromley", file=file)
  call Print("IPModel_FB : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call verbosity_push_decrement()
    do tj=1, this%n_types
      call Print ("IPModel_FB : interaction " // ti // " " // tj // " A " // this%A(ti,tj) // " B " // &
      & this%B(ti,tj) // " C " // this%C(ti,tj), file=file)
    end do
    call verbosity_pop()
  end do

end subroutine IPModel_FB_Print

end module IPModel_FB_module

