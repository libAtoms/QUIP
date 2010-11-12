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
!X IPModel_Morse module  
!X
!% Module for Morse pair potential.
!% \begin{equation} 
!%   \nonumber
!%     V(r) = D \left( \exp(-2 \alpha (r-r_0)) - 2 \exp( -\alpha (r-r_0)) \right)
!% \end{equation} 
!%
!% The IPModel_Morse object contains all the parameters read from a
!% 'Morse_params' XML stanza.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_Morse_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none

private 

include 'IPModel_interface.h'

public :: IPModel_Morse
type IPModel_Morse
  integer :: n_types = 0         !% Number of atomic types. 
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)  !% Atomic number dimensioned as \texttt{n_types}. 

  real(dp) :: cutoff = 0.0_dp    !% Cutoff for computing connection.

  real(dp), allocatable :: D(:,:), alpha(:,:), r0(:,:), cutoff_a(:,:) !% IP parameters.

  character(len=FIELD_LENGTH) label

end type IPModel_Morse

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_Morse), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_Morse_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Morse_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Morse_Print
end interface Print

interface Calc
  module procedure IPModel_Morse_Calc
end interface Calc

contains

subroutine IPModel_Morse_Initialise_str(this, args_str, param_str)
  type(IPModel_Morse), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label = ''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_Morse_Initialise_str args_str')) then
    call system_abort("IPModel_Morse_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_Morse_read_params_xml(this, param_str)

  this%cutoff = maxval(this%cutoff_a)

end subroutine IPModel_Morse_Initialise_str

subroutine IPModel_Morse_Finalise(this)
  type(IPModel_Morse), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if (allocated(this%r0)) deallocate(this%r0)
  if (allocated(this%D)) deallocate(this%D)
  if (allocated(this%alpha)) deallocate(this%alpha)
  if (allocated(this%cutoff_a)) deallocate(this%cutoff_a)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_Morse_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator: this routine computes energy, forces and the virial.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_Morse_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_Morse), intent(inout) :: this
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
  logical :: i_is_min_image

  integer :: i_calc, n_extra_calcs
  character(len=20) :: extra_calcs_list(10)

  logical :: do_flux = .false.
  real(dp), pointer :: velo(:,:)
  real(dp) :: flux(3)

  type(Dictionary)                :: params
  logical :: has_atom_mask_name
  character(FIELD_LENGTH) :: atom_mask_name

   INIT_ERROR(error)

  if (present(e)) e = 0.0_dp

  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_Morse_Calc', error)
     local_e = 0.0_dp
  endif

  if (present(f)) then 
     call check_size('Force',f,(/3,at%N/),'IPModel_Morse_Calc', error)
     f = 0.0_dp
  end if

  if (present(virial)) virial = 0.0_dp

  if (present(local_virial)) then
     call check_size('Local_virial',local_virial,(/9,at%N/),'IPModel_Morse_Calc', error)
     local_virial = 0.0_dp
     RAISE_ERROR("IPModel_Morse_Calc: local_virial calculation requested but not supported yet.", error)
  endif

  if (present(args_str)) then
    if (len_trim(args_str) > 0) then
      n_extra_calcs = parse_extra_calcs(args_str, extra_calcs_list)
      if (n_extra_calcs > 0) then
	do i_calc=1, n_extra_calcs
	  select case(trim(extra_calcs_list(i_calc)))
	    case("flux")
	      if (.not. assign_pointer(at, "velo", velo)) &
		call system_abort("IPModel_Morse_Calc Flux calculation requires velo field")
	      do_flux = .true.
	      flux = 0.0_dp
	    case default
	      call system_abort("Unsupported extra_calc '"//trim(extra_calcs_list(i_calc))//"'")
	  end select
	end do
      endif ! n_extra_calcs
    endif ! len_trim(args_str)
    call initialise(params)
    call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")

    if(.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_Morse_Calc args_str')) then
       RAISE_ERROR("IPModel_Morse_Calc failed to parse args_str='"//trim(args_str)//"'",error)
    endif
    call finalise(params)
    if(has_atom_mask_name) then
       RAISE_ERROR('IPModel_Morse_Calc: atom_mask_name found, but not supported', error)
    endif
  endif ! present(args_str)

  if (.not. assign_pointer(at, "weight", w_e)) nullify(w_e)

  do i = 1, at%N
    i_is_min_image = at%connect%is_min_image(i)

    if (present(mpi)) then
       if (mpi%active) then
	 if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
       endif
    endif

    do ji = 1, atoms_n_neighbours(at, i)
      j = atoms_neighbour(at, i, ji, dr_mag, cosines = dr)
      
      if (dr_mag .feq. 0.0_dp) cycle
      if ((i < j) .and. i_is_min_image) cycle

      ti = get_type(this%type_of_atomic_num, at%Z(i))
      tj = get_type(this%type_of_atomic_num, at%Z(j))

      if (present(e) .or. present(local_e)) then
	de = IPModel_Morse_pairenergy(this, ti, tj, dr_mag)
	if (present(local_e)) then
	  local_e(i) = local_e(i) + 0.5_dp*de
          if(i_is_min_image) local_e(j) = local_e(j) + 0.5_dp*de
	endif
	if (present(e)) then
	  if (associated(w_e)) then
	    de = de*0.5_dp*(w_e(i)+w_e(j))
	  endif
          if(i_is_min_image) then
             e = e + de
          else
             e = e + 0.5_dp*de
          endif
	endif
      endif
      if (present(f) .or. present(virial) .or. do_flux) then
	de_dr = IPModel_Morse_pairenergy_deriv(this, ti, tj, dr_mag)
	if (associated(w_e)) then
	  de_dr = de_dr*0.5_dp*(w_e(i)+w_e(j))
	endif
	if (present(f)) then
	  f(:,i) = f(:,i) + de_dr*dr
	  if(i_is_min_image) f(:,j) = f(:,j) - de_dr*dr
	endif
	if (do_flux) then
	  ! -0.5 (v_i + v_j) . F_ij * dr_ij
	  flux = flux - 0.5_dp*sum((velo(:,i)+velo(:,j))*(de_dr*dr))*(dr*dr_mag)
	endif
	if (present(virial)) then
	  if(i_is_min_image) then
             virial = virial - de_dr*(dr .outer. dr)*dr_mag
          else
             virial = virial - 0.5_dp*de_dr*(dr .outer. dr)*dr_mag
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
  if (do_flux) then
    flux = flux / cell_volume(at)
    if (present(mpi)) call sum_in_place(mpi, flux)
    call set_value(at%params, "Flux", flux)
  endif

end subroutine IPModel_Morse_Calc

!% This routine computes the two-body term for a pair of atoms  separated by a distance r.
function IPModel_Morse_pairenergy(this, ti, tj, r)
  type(IPModel_Morse), intent(in) :: this
  integer, intent(in) :: ti, tj   !% Atomic types.
  real(dp), intent(in) :: r       !% Distance.
  real(dp) :: IPModel_Morse_pairenergy

  real(dp) :: texp

  if ((r .feq. 0.0_dp) .or. (r > this%cutoff_a(ti,tj))) then
    IPModel_Morse_pairenergy = 0.0
    return
  endif

  texp = exp(-this%alpha(ti,tj)*(r-this%r0(ti,tj)))

  IPModel_Morse_pairenergy = this%D(ti,tj) * ( texp*texp - 2.0_dp*texp )
end function IPModel_Morse_pairenergy

!% Derivative of the two-body term.
function IPModel_Morse_pairenergy_deriv(this, ti, tj, r)
  type(IPModel_Morse), intent(in) :: this
  integer, intent(in) :: ti, tj   !% Atomic types.
  real(dp), intent(in) :: r       !% Distance.
  real(dp) :: IPModel_Morse_pairenergy_deriv

  real(dp) :: texp, texp_d

  if ((r .feq. 0.0_dp) .or. (r > this%cutoff_a(ti,tj))) then
    IPModel_Morse_pairenergy_deriv = 0.0
    return
  endif

  texp = exp(-this%alpha(ti,tj)*(r-this%r0(ti,tj)))
  texp_d = -this%alpha(ti,tj) * texp 

  IPModel_Morse_pairenergy_deriv = this%D(ti,tj) * ( 2.0_dp*texp*texp_d - 2.0_dp * texp_d)
end function IPModel_Morse_pairenergy_deriv


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions.
!% An example for XML stanza is given below.  Please notice that
!% these are simply dummy parameters for testing purposes, with no physical meaning.
!%
!%> <Morse_params n_types="2" label="default">
!%>   <pair atnum_i="6" atnum_j="6"   r0="2.0" D="1.0" alpha="0.1" cutoff="4.0" />
!%>   <pair atnum_i="6" atnum_j="14"  r0="2.5" D="2.0" alpha="0.15" cutoff="5.0" />
!%>   <pair atnum_i="14" atnum_j="14" r0="3.0" D="3.0" alpha="0.2" cutoff="6.0" />
!%> </Morse_params>
!%
!% cutoff defaults to r0 - log(1e-8) / alpha
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=1024) :: value

  integer atnum_i, atnum_j, ti, tj, ti_a(1)

  if (name == 'Morse_params') then ! new Morse stanza
    if (parse_in_ip) &
      call system_abort("IPModel_startElement_handler entered Morse_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

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
	call system_abort("Can't find n_types in Morse_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0
      allocate(parse_ip%r0(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%D(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%alpha(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%cutoff_a(parse_ip%n_types,parse_ip%n_types))

      parse_ip%r0 = 0.0_dp
      parse_ip%D = 0.0_dp
      parse_ip%alpha = 0.0_dp

    endif ! parse_in_ip
  elseif (parse_in_ip .and. name == 'pair') then

    call QUIP_FoX_get_value(attributes, "atnum_i", value, status)
    if (status /= 0) call system_abort ("IPModel_Morse_read_params_xml cannot find atnum_i")
    read (value, *) atnum_i
    call QUIP_FoX_get_value(attributes, "atnum_j", value, status)
    if (status /= 0) call system_abort ("IPModel_Morse_read_params_xml cannot find atnum_j")
    read (value, *) atnum_j

    if (all(parse_ip%atomic_num /= atnum_i)) then
      ti_a = minloc(parse_ip%atomic_num)
      parse_ip%atomic_num(ti_a(1)) = atnum_i
    endif
    if (all(parse_ip%atomic_num /= atnum_j)) then
      ti_a = minloc(parse_ip%atomic_num)
      parse_ip%atomic_num(ti_a(1)) = atnum_j
    endif
    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
	parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

    ti = parse_ip%type_of_atomic_num(atnum_i)
    tj = parse_ip%type_of_atomic_num(atnum_j)

    call QUIP_FoX_get_value(attributes, "r0", value, status)
    if (status /= 0) call system_abort ("IPModel_Morse_read_params_xml cannot find r0")
    read (value, *) parse_ip%r0(ti,tj)
    call QUIP_FoX_get_value(attributes, "D", value, status)
    if (status /= 0) call system_abort ("IPModel_Morse_read_params_xml cannot find D")
    read (value, *) parse_ip%D(ti,tj)
    call QUIP_FoX_get_value(attributes, "alpha", value, status)
    if (status /= 0) call system_abort ("IPModel_Morse_read_params_xml cannot find alpha")
    read (value, *) parse_ip%alpha(ti,tj)

    parse_ip%cutoff_a(ti,tj) = -1.0_dp
    call QUIP_FoX_get_value(attributes, "cutoff", value, status)
    if (status == 0) read (value, *) parse_ip%cutoff_a(ti,tj)

    if (ti /= tj) then
      parse_ip%r0(tj,ti) = parse_ip%r0(ti,tj)
      parse_ip%D(tj,ti) = parse_ip%D(ti,tj)
      parse_ip%alpha(tj,ti) = parse_ip%alpha(ti,tj)
      parse_ip%cutoff_a(tj,ti) = parse_ip%cutoff_a(ti,tj)
    endif

  endif ! parse_in_ip .and. name = 'Morse'

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if (name == 'Morse_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_Morse_read_params_xml(this, param_str)
  type(IPModel_Morse), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type(xml_t) :: fxml
  integer :: ti, tj

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
    call system_abort("IPModel_Morse_read_params_xml parsed file, but n_types = 0")
  endif

  ! default cutoff:
  ! exp(-alpha (r-r0)) < 1e-8
  ! -alpha (r-r0) < log(1e-8)
  ! r - r0 < log(1e-8) / -alpha
  ! r < r0 + log(1e-8) / -alpha
  ! r < r0 + 13.8 / alpha

  do ti=1, this%n_types
  do tj=1, this%n_types
    if (this%cutoff_a(ti,tj) < 0.0_dp) this%cutoff_a(ti,tj) = this%r0(ti,tj) - log(1.0e-8_dp) / this%alpha(ti,tj)
  end do
  end do
  this%cutoff = maxval(this%cutoff_a)

end subroutine IPModel_Morse_read_params_xml



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of Morse parameters: number of different types, cutoff radius, atomic numbers, etc.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_Morse_Print (this, file)
  type(IPModel_Morse), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_Morse : Morse", file=file)
  call Print("IPModel_Morse : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_Morse : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    do tj=1, this%n_types
      call Print ("IPModel_Morse : interaction " // ti // " " // tj // " r0 " // this%r0(ti,tj) // " D " // &
	this%D(ti,tj) // " alpha " // this%alpha(ti,tj) // " cutoff_a " // this%cutoff_a(ti,tj), file=file)
    end do
    call verbosity_pop()
  end do

end subroutine IPModel_Morse_Print

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

end module IPModel_Morse_module
