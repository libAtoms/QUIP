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
! H0 X   Public License, version 2, https://www.gnu.org/copyleft/gpl.html
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
! H0 X   https://www.libatoms.org
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
use system_module, only : dp, inoutput, print, verbosity_push_decrement, verbosity_pop, split_string_simple, string_to_numerical, operator(//)
use extendable_str_module  
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
type shoulder_params_t
   integer :: n = 0
   real(dp), dimension(:), allocatable :: sigma1, lambda, k
   real(dp) :: lambda0 = 0.0_dp
endtype shoulder_params_t

type IPModel_RS
  integer :: n_types = 0         !% Number of atomic types.
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)  !% Atomic number dimensioned as \texttt{n_types}.

  real(dp) :: cutoff = 0.0_dp    !% Cutoff for computing connection.

  real(dp), allocatable :: eps(:,:), sigma(:,:)
  type(shoulder_params_t), dimension(:,:), allocatable :: shoulder_params

  character(len=STRING_LENGTH) label

end type IPModel_RS

logical, private :: parse_in_ip, parse_matched_label
integer, private :: parse_ti, parse_tj
type(extendable_str), private, save :: parse_element_data
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
  integer :: i, j

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if (allocated(this%sigma)) deallocate(this%sigma)
  if (allocated(this%eps)) deallocate(this%eps)
  if (allocated(this%shoulder_params)) then
     do i = 1, this%n_types
        do j = 1, this%n_types
           if(allocated(this%shoulder_params(j,i)%sigma1)) deallocate(this%shoulder_params(j,i)%sigma1)
           if(allocated(this%shoulder_params(j,i)%lambda)) deallocate(this%shoulder_params(j,i)%lambda)
           if(allocated(this%shoulder_params(j,i)%k)) deallocate(this%shoulder_params(j,i)%k)
           this%shoulder_params(j,i)%n = 0
           this%shoulder_params(j,i)%lambda0 = 0.0_dp
        enddo
     enddo
     deallocate(this%shoulder_params)
  endif

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
  real(dp) :: de, de_dr, virial_i(3,3)

  type(Dictionary)                :: params
  logical :: has_atom_mask_name
  character(STRING_LENGTH) :: atom_mask_name
  logical, dimension(:), pointer :: atom_mask_pointer
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
  endif

  atom_mask_pointer => null()
  if (present(args_str)) then
    call initialise(params)
    call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, help_string="Name of logical property used to mask atoms")
    call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
    call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")

    if(.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_RS_Calc args_str')) then
       RAISE_ERROR("IPModel_RS_Calc failed to parse args_str='"//trim(args_str)//"'",error)
    endif
    call finalise(params)
    if(has_atom_mask_name) then
        if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) &
        call system_abort("IPModel_RS_Calc did not find "//trim(atom_mask_name)//" property in the atoms object.")
    else
        atom_mask_pointer => null()
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

    if(associated(atom_mask_pointer)) then
       if(.not. atom_mask_pointer(i)) cycle
    endif

    ti = get_type(this%type_of_atomic_num, at%Z(i))

    do ji = 1, n_neighbours(at, i)
      j = neighbour(at, i, ji, dr_mag, cosines = dr)

      if (dr_mag .feq. 0.0_dp) cycle
      !if (j <= 0) cycle

      tj = get_type(this%type_of_atomic_num, at%Z(j))

      if (present(e) .or. present(local_e)) then
        de = 0.5_dp * IPModel_RS_pairenergy(this, ti, tj, dr_mag)

        if (present(local_e)) then
          local_e(i) = local_e(i) + de
        endif
        if (present(e)) then
          if (associated(w_e)) then
            de = de*0.5_dp*(w_e(i)+w_e(j))
          endif
          e = e + de
        endif
      endif
      if (present(f) .or. present(virial)) then
        de_dr = 0.5_dp*IPModel_RS_pairenergy_deriv(this, ti, tj, dr_mag)
        if (associated(w_e)) then
          de_dr = de_dr*0.5_dp*(w_e(i)+w_e(j))
        endif
        if (present(f)) then
          f(:,i) = f(:,i) + de_dr*dr
          f(:,j) = f(:,j) - de_dr*dr
        endif

        if (present(virial) .or. present(local_virial) ) virial_i = de_dr*(dr .outer. dr)*dr_mag
        if (present(virial)) virial = virial - virial_i
        if (present(local_virial)) local_virial(:,i) = local_virial(:,i) - reshape(virial_i, (/9/))
      endif
    end do
  end do
  if (present(mpi)) then
     if (present(e)) e = sum(mpi, e)
     if (present(local_e)) call sum_in_place(mpi, local_e)
     if (present(virial)) call sum_in_place(mpi, virial)
     if (present(local_virial)) call sum_in_place(mpi, local_virial)
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
     sum( this%shoulder_params(ti,tj)%lambda * tanh( this%shoulder_params(ti,tj)%k * ( r - this%shoulder_params(ti,tj)%sigma1 ) ) ) + &
     this%shoulder_params(ti,tj)%lambda0 )

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

  IPModel_RS_pairenergy_deriv = this%eps(ti,tj) * ( &
     - 14.0_dp * ( this%sigma(ti,tj) / r )**14 / r + &
     sum( this%shoulder_params(ti,tj)%lambda * &
        this%shoulder_params(ti,tj)%k / cosh( this%shoulder_params(ti,tj)%k * ( r - this%shoulder_params(ti,tj)%sigma1 ) )**2 ) )

end function IPModel_RS_pairenergy_deriv

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!% XML param reader functions.
!% An example for XML stanza is given below, please notice that
!% they are simply dummy parameters for testing purposes, with no physical meaning.
!%
!% <quip>
!% <RS_params n_types="1" cutoff="10.0" label="default">
!% <per_type_data type="1" atomic_num="29" />
!% <per_pair_data type1="1" type2="1" sigma="1.0" eps="1.0" n="1">
!%    <lambda>-0.5</lambda>
!%    <sigma1>1.35</sigma1>
!%    <k>10</k>
!% </per_pair_data>
!% </RS_params>
!% <RS_params n_types="1" cutoff="10.0" label="attractive">
!% <per_type_data type="1" atomic_num="29" />
!% <per_pair_data type1="1" type2="1" sigma="1.0" eps="1.0" n="2">
!%    <lambda>-0.8 0.3</lambda>
!%    <sigma1>1.35 1.80</sigma1>
!%    <k>10 10</k>
!% </per_pair_data>
!% </RS_params>
!% </quip>
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=1024) :: value

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
      allocate(parse_ip%shoulder_params(parse_ip%n_types,parse_ip%n_types))
    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find type")
    read (value, *) parse_ti

    if (parse_ti < 1) call system_abort("IPModel_RS_read_params_xml got per_type_data type="//parse_ti//" < 1")
    if (parse_ti > parse_ip%n_types) call system_abort("IPModel_RS_read_params_xml got per_type_data type="//parse_ti//" > n_types="//parse_ip%n_types)

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(parse_ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do parse_ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(parse_ti) > 0) &
        parse_ip%type_of_atomic_num(parse_ip%atomic_num(parse_ti)) = parse_ti
    end do

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "type1", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find type1")
    read (value, *) parse_ti
    call QUIP_FoX_get_value(attributes, "type2", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find type2")
    read (value, *) parse_tj

    if (parse_ti < 1) call system_abort("IPModel_RS_read_params_xml got per_type_data type1="//parse_ti//" < 1")
    if (parse_ti > parse_ip%n_types) call system_abort("IPModel_RS_read_params_xml got per_pair_data type1="//parse_ti//" > n_types="//parse_ip%n_types)
    if (parse_tj < 1) call system_abort("IPModel_RS_read_params_xml got per_type_data type2="//parse_tj//" < 1")
    if (parse_tj > parse_ip%n_types) call system_abort("IPModel_RS_read_params_xml got per_pair_data type2="//parse_tj//" > n_types="//parse_ip%n_types)

    call QUIP_FoX_get_value(attributes, "sigma", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find sigma")
    read (value, *) parse_ip%sigma(parse_ti,parse_tj)
    call QUIP_FoX_get_value(attributes, "eps", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find eps6")
    read (value, *) parse_ip%eps(parse_ti,parse_tj)



    call QUIP_FoX_get_value(attributes, "n", value, status)
    if (status /= 0) call system_abort ("IPModel_RS_read_params_xml cannot find n")
    read (value, *) parse_ip%shoulder_params(parse_ti,parse_tj)%n
    allocate( parse_ip%shoulder_params(parse_ti,parse_tj)%sigma1( parse_ip%shoulder_params(parse_ti,parse_tj)%n ), &
       parse_ip%shoulder_params(parse_ti,parse_tj)%lambda( parse_ip%shoulder_params(parse_ti,parse_tj)%n ), &
       parse_ip%shoulder_params(parse_ti,parse_tj)%k( parse_ip%shoulder_params(parse_ti,parse_tj)%n ) )

    if (parse_ti /= parse_tj) then
      parse_ip%eps(parse_tj,parse_ti) = parse_ip%eps(parse_ti,parse_tj)
      parse_ip%sigma(parse_tj,parse_ti) = parse_ip%sigma(parse_ti,parse_tj)
      parse_ip%shoulder_params(parse_tj,parse_ti)%n = parse_ip%shoulder_params(parse_ti,parse_tj)%n

      allocate( parse_ip%shoulder_params(parse_tj,parse_ti)%sigma1( parse_ip%shoulder_params(parse_tj,parse_ti)%n ), &
         parse_ip%shoulder_params(parse_tj,parse_ti)%lambda( parse_ip%shoulder_params(parse_tj,parse_ti)%n ), &
         parse_ip%shoulder_params(parse_tj,parse_ti)%k( parse_ip%shoulder_params(parse_tj,parse_ti)%n ) )
    end if

  elseif (parse_in_ip .and. name == 'lambda') then
     call zero(parse_element_data)
  elseif (parse_in_ip .and. name == 'sigma1') then
     call zero(parse_element_data)
  elseif (parse_in_ip .and. name == 'k') then
     call zero(parse_element_data)
  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'RS_params') then
       parse_in_ip = .false.
    elseif( name == "lambda" ) then
       call string_to_numerical( string(parse_element_data), parse_ip%shoulder_params(parse_ti,parse_tj)%lambda )
       parse_ip%shoulder_params(parse_ti,parse_tj)%lambda0 = -sum(parse_ip%shoulder_params(parse_ti,parse_tj)%lambda)
       if( parse_ti /= parse_tj ) then
          parse_ip%shoulder_params(parse_tj,parse_ti)%lambda = parse_ip%shoulder_params(parse_ti,parse_tj)%lambda
          parse_ip%shoulder_params(parse_tj,parse_ti)%lambda0  = parse_ip%shoulder_params(parse_ti,parse_tj)%lambda0
       endif
    elseif( name == "sigma1" ) then
       call string_to_numerical( string(parse_element_data), parse_ip%shoulder_params(parse_ti,parse_tj)%sigma1 )
       if( parse_ti /= parse_tj ) parse_ip%shoulder_params(parse_tj,parse_ti)%sigma1 = parse_ip%shoulder_params(parse_ti,parse_tj)%sigma1
    elseif( name == "k" ) then
       call string_to_numerical( string(parse_element_data), parse_ip%shoulder_params(parse_ti,parse_tj)%k )
       if( parse_ti /= parse_tj ) parse_ip%shoulder_params(parse_tj,parse_ti)%k = parse_ip%shoulder_params(parse_ti,parse_tj)%k
    end if

  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_characters_handler(in)
   character(len=*), intent(in) :: in

   if( parse_in_ip ) then
      call concat(parse_element_data, in, keep_lf=.false.,lf_to_whitespace=.true.)
   endif
endsubroutine IPModel_characters_handler

subroutine IPModel_RS_read_params_xml(this, param_str)
  type(IPModel_RS), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type(xml_t) :: fxml

  if (len(trim(param_str)) <= 0) return

  parse_in_ip = .false.
  parse_matched_label = .false.
  parse_ip => this
  call initialise(parse_element_data)

  call open_xml_string(fxml, param_str)
  call parse(fxml,  &
    characters_handler = IPModel_characters_handler, &
    startElement_handler = IPModel_startElement_handler, &
    endElement_handler = IPModel_endElement_handler)
  call close_xml_t(fxml)

  call finalise(parse_element_data)

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
  call Print("IPModel_RS : label = "//this%label, file=file)
  call Print("IPModel_RS : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_RS : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    do tj=1, this%n_types
      call Print ("IPModel_RS : interaction " // ti // " " // tj // " sigma " // this%sigma(ti,tj) // " eps " // &
        this%eps(ti,tj) // " sigma1 " // this%shoulder_params(ti,tj)%sigma1 // " k " // this%shoulder_params(ti,tj)%k // &
        " lambda " // this%shoulder_params(ti,tj)%lambda, file=file)
    end do
    call verbosity_pop()
  end do

end subroutine IPModel_RS_Print

end module IPModel_RS_module
