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
!X IPModel_RS module  
!X
!% Module for a Repulsive Step pair potential.
!% \begin{equation} 
!%   \nonumber
!%     V(r) = \epsilon \left[  \left( \frac{\sigma}{r} \right)^{14} + \frac{1}{2} \left( 1 - \tanh (k (r - \sigma_1)) \right) \right] 
!% \end{equation} 
!% From J. Chem. Phys. 129, 064512 (2008)
!%
!% The IPModel_RS object contains all the parameters read from a
!% 'RS_params' XML stanza.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_RS_module

use error_module
use system_module, only : dp, inoutput, print, verbosity_push_decrement, verbosity_pop, split_string_simple, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use units_module, only : PI
use atoms_types_module
use atoms_module

use mpi_context_module
use QUIP_Common_module

implicit none

private 

include 'IPModel_interface.h'

public :: IPModel_RS
type IPModel_RS
  integer :: n_types = 0         !% Number of atomic types. 
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)  !% Atomic number dimensioned as \texttt{n_types}. 

  real(dp) :: cutoff = 0.0_dp    !% Cutoff for computing connection.

  real(dp), allocatable :: eps(:,:), sigma(:,:), sigma1(:,:), k(:,:)

  character(len=STRING_LENGTH) label

end type IPModel_RS

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_RS), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_RS_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_RS_Finalise
end interface Finalise

interface Print
  module procedure IPModel_RS_Print
end interface Print

interface Calc
  module procedure IPModel_RS_Calc
end interface Calc

contains

subroutine IPModel_RS_Initialise_str(this, args_str, param_str)
  type(IPModel_RS), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label = ''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_RS_Initialise_str args_str')) then
    call system_abort("IPModel_RS_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_RS_read_params_xml(this, param_str)

end subroutine IPModel_RS_Initialise_str

subroutine IPModel_RS_Finalise(this)
  type(IPModel_RS), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if (allocated(this%sigma)) deallocate(this%sigma)
  if (allocated(this%sigma1)) deallocate(this%sigma1)
  if (allocated(this%eps)) deallocate(this%eps)
  if (allocated(this%k)) deallocate(this%k)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_RS_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator: this routine computes energy, forces and the virial.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_RS_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_RS), intent(inout) :: this
  type(Atoms), intent(inout) :: at
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial
  character(len=*), intent(in), optional      :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  real(dp), pointer :: w_e(:)
  integer i, ji, j, ti, tj
  real(dp) :: dr(3), dr_mag
  real(dp) :: de, de_dr

  type(Dictionary)                :: params
  logical :: has_atom_mask_name
  character(STRING_LENGTH) :: atom_mask_name
  real(dp) :: r_scale, E_scale
  logical :: do_rescale_r, do_rescale_E

  INIT_ERROR(error)

  if (present(e)) e = 0.0_dp
  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_RS_Calc', error)
     local_e = 0.0_dp
  endif
  if (present(f)) then 
     call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_RS_Calc', error)
     f = 0.0_dp
  end if
  if (present(virial)) virial = 0.0_dp
  if (present(local_virial)) then
     call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_RS_Calc', error)
     local_virial = 0.0_dp
     RAISE_ERROR("IPModel_RS_Calc: local_virial calculation requested but not supported yet.", error)
  endif

  if (present(args_str)) then
    call initialise(params)
    call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
    call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")

    if(.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_RS_Calc args_str')) then
       RAISE_ERROR("IPModel_RS_Calc failed to parse args_str='"//trim(args_str)//"'",error)
    endif
    call finalise(params)
    if(has_atom_mask_name) then
       RAISE_ERROR('IPModel_RS_Calc: atom_mask_name found, but not supported', error)
    endif
    if (do_rescale_r .or. do_rescale_E) then
       RAISE_ERROR("IPModel_RS_Calc: rescaling of potential with r_scale and E_scale not yet implemented!", error)
    end if

  endif ! present(args_str)

  if (.not. assign_pointer(at, "weight", w_e)) nullify(w_e)

  do i = 1, at%N

    if (present(mpi)) then
       if (mpi%active) then
	 if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
       endif
    endif

    ti = get_type(this%type_of_atomic_num, at%Z(i))

    do ji = 1, n_neighbours(at, i)
      j = neighbour(at, i, ji, dr_mag, cosines = dr)

      if (dr_mag .feq. 0.0_dp) cycle
      if ((i < j)) cycle

      tj = get_type(this%type_of_atomic_num, at%Z(j))

      if (present(e) .or. present(local_e)) then
	de = IPModel_RS_pairenergy(this, ti, tj, dr_mag)

	if (present(local_e)) then
	  local_e(i) = local_e(i) + 0.5_dp*de
          if(i/=j) local_e(j) = local_e(j) + 0.5_dp*de
	endif
	if (present(e)) then
	  if (associated(w_e)) then
	    de = de*0.5_dp*(w_e(i)+w_e(j))
	  endif
          if(i==j) then
             e = e + 0.5_dp*de
          else
             e = e + de
          endif
	endif
      endif
      if (present(f) .or. present(virial)) then
	de_dr = IPModel_RS_pairenergy_deriv(this, ti, tj, dr_mag)
	if (associated(w_e)) then
	  de_dr = de_dr*0.5_dp*(w_e(i)+w_e(j))
	endif
	if (present(f)) then
	  f(:,i) = f(:,i) + de_dr*dr
	  if(i/=j) f(:,j) = f(:,j) - de_dr*dr
	endif
	if (present(virial)) then
          if(i==j) then
             virial = virial - 0.5_dp*de_dr*(dr .outer. dr)*dr_mag
          else
             virial = virial - de_dr*(dr .outer. dr)*dr_mag
          endif
	endif
      endif
    end do
  end do

  if (present(mpi)) then
     if (present(e)) e = sum(mpi, e)
     if (present(local_e)) call sum_in_place(mpi, local_e)
     if (present(virial)) call sum_in_place(mpi, virial)
     if (present(f)) call sum_in_place(mpi, f)
  endif

end subroutine IPModel_RS_Calc

!% This routine computes the two-body term for a pair of atoms  separated by a distance r.
function IPModel_RS_pairenergy(this, ti, tj, r)
  type(IPModel_RS), intent(in) :: this
  integer, intent(in) :: ti, tj   !% Atomic types.
  real(dp), intent(in) :: r       !% Distance.
  real(dp) :: IPModel_RS_pairenergy

  if ((r .feq. 0.0_dp) .or. (r > this%cutoff)) then
    IPModel_RS_pairenergy = 0.0_dp
    return
  endif

  IPModel_RS_pairenergy = this%eps(ti,tj) * ( &
     ( this%sigma(ti,tj) / r )**14 + &
     0.5_dp * ( 1.0_dp - tanh( this%k(ti,tj) * ( r - this%sigma1(ti,tj) ) ) ) &
     )

end function IPModel_RS_pairenergy

!% Derivative of the two-body term.
function IPModel_RS_pairenergy_deriv(this, ti, tj, r)
  type(IPModel_RS), intent(in) :: this
  integer, intent(in) :: ti, tj   !% Atomic types.
  real(dp), intent(in) :: r       !% Distance.
  real(dp) :: IPModel_RS_pairenergy_deriv

  if ((r .feq. 0.0_dp) .or. (r > this%cutoff)) then
    IPModel_RS_pairenergy_deriv = 0.0_dp
    return
  endif

  IPModel_RS_pairenergy_deriv = - this%eps(ti,tj) * ( &
     14.0_dp * ( this%sigma(ti,tj) / r )**14 / r + &
     0.5_dp * this%k(ti,tj) / cosh( this%k(ti,tj) * ( r - this%sigma1(ti,tj) ) )**2 &
     )
  
end function IPModel_RS_pairenergy_deriv

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions.
!% An example for XML stanza is given below, please notice that
!% they are simply dummy parameters for testing purposes, with no physical meaning.
!%
!%> <RS_params n_types="2" cutoff="10.0" label="default">
!%> <per_type_data type="1" atomic_num="29" />
!%> <per_type_data type="2" atomic_num="79" />
!%> <per_pair_data type1="1" type2="1" sigma="1.0" eps="1.0" 
!%>       sigma1="1.35" k="10" />
!%> <per_pair_data type1="1" type2="2" sigma="1.0" eps="1.0" 
!%>       sigma1="1.35" k="10" />
!%> <per_pair_data type1="2" type2="2" sigma="1.0" eps="1.0" 
!%>       sigma1="1.35" k="10" />
!%> </RS_params>
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=1024) :: value

  integer :: ti, tj

  if (name == 'RS_params') then ! new RS stanza

    if (parse_in_ip) &
      call system_abort("IPModel_startElement_handler entered RS_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

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
	read (value, *) parse_ip%n_types
      else
	call system_abort("Can't find n_types in RS_params")
      endif

      call QUIP_FoX_get_value(attributes, 'cutoff', value, status)
      if (status == 0) then
	read (value, *) parse_ip%cutoff
      else
	call system_abort("Can't find cutoff in RS_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0
      allocate(parse_ip%sigma(parse_ip%n_types,parse_ip%n_types))
      parse_ip%sigma = 1.0_dp
      allocate(parse_ip%eps(parse_ip%n_types,parse_ip%n_types))
      parse_ip%eps = 0.0_dp
      allocate(parse_ip%sigma1(parse_ip%n_types,parse_ip%n_types))
      parse_ip%sigma1 = 0.0_dp
      allocate(parse_ip%k(parse_ip%n_types,parse_ip%n_types))
      parse_ip%k = 0.0_dp
    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find type")
    read (value, *) ti

    if (ti < 1) call system_abort("IPModel_RS_read_params_xml got per_type_data type="//ti//" < 1")
    if (ti > parse_ip%n_types) call system_abort("IPModel_RS_read_params_xml got per_type_data type="//ti//" > n_types="//parse_ip%n_types)

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
	parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "type1", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find type1")
    read (value, *) ti
    call QUIP_FoX_get_value(attributes, "type2", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find type2")
    read (value, *) tj

    if (ti < 1) call system_abort("IPModel_RS_read_params_xml got per_type_data type1="//ti//" < 1")
    if (ti > parse_ip%n_types) call system_abort("IPModel_RS_read_params_xml got per_pair_data type1="//ti//" > n_types="//parse_ip%n_types)
    if (tj < 1) call system_abort("IPModel_RS_read_params_xml got per_type_data type2="//tj//" < 1")
    if (tj > parse_ip%n_types) call system_abort("IPModel_RS_read_params_xml got per_pair_data type2="//tj//" > n_types="//parse_ip%n_types)

    call QUIP_FoX_get_value(attributes, "sigma", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find sigma")
    read (value, *) parse_ip%sigma(ti,tj)
    call QUIP_FoX_get_value(attributes, "eps", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find eps6")
    read (value, *) parse_ip%eps(ti,tj)
    call QUIP_FoX_get_value(attributes, "sigma1", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find sigma1")
    read (value, *) parse_ip%sigma1(ti,tj)
    call QUIP_FoX_get_value(attributes, "k", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find k")
    read (value, *) parse_ip%k(ti,tj)

    if (ti /= tj) then
      parse_ip%eps(tj,ti) = parse_ip%eps(ti,tj)
      parse_ip%sigma(tj,ti) = parse_ip%sigma(ti,tj)
      parse_ip%sigma1(tj,ti) = parse_ip%sigma1(ti,tj)
      parse_ip%k(tj,ti) = parse_ip%k(ti,tj)
    end if

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if (name == 'RS_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_RS_read_params_xml(this, param_str)
  type(IPModel_RS), intent(inout), target :: this
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
    call system_abort("IPModel_RS_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_RS_read_params_xml


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of RS parameters: number of different types, cutoff radius, atomic numbers, etc.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_RS_Print (this, file)
  type(IPModel_RS), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_RS : Repulsive Step", file=file)
  call Print("IPModel_RS : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_RS : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    do tj=1, this%n_types
      call Print ("IPModel_RS : interaction " // ti // " " // tj // " sigma " // this%sigma(ti,tj) // " eps " // &
	this%eps(ti,tj) // " sigma1 " // this%sigma1(ti,tj) // " k " // this%k(ti,tj), file=file)
    end do
    call verbosity_pop()
  end do

end subroutine IPModel_RS_Print

end module IPModel_RS_module
