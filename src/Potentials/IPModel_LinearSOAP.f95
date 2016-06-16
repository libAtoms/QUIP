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
!X IPModel_LinearSOAP module  
!X
!% Module for a Linear SOAP potential.
!% \begin{equation} 
!%   \nonumber
!%     E = \sum_i \beta \cdot p_{ss'nn'l,i}
!% \end{equation} 
!% where $p_{ss'nn'l,i}$ is the SOAP sphericial power spectrum of atom i, and
!% $\beta$ is a vector of coefficients. 
!%
!% The IPModel_LinearSOAP object contains all the parameters read from a
!% 'LinearSOAP_params' XML stanza.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_LinearSOAP_module

use error_module
use system_module, only : dp, inoutput, print, verbosity_push_decrement, verbosity_pop, split_string_simple, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module
#ifdef HAVE_GAP
use descriptors_module
#endif 
use mpi_context_module
use QUIP_Common_module
use extendable_str_module

implicit none

private 

include 'IPModel_interface.h'

public :: IPModel_LinearSOAP
type IPModel_LinearSOAP
  integer :: n_types = 0         !% Number of atomic types. 
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)  !% Atomic number dimensioned as \texttt{n_types}. 

  real(dp) :: cutoff = 0.0_dp     !% Cutoff for computing connection.
  real(dp) :: atom_sigma = 0.0_dp !% atom_sigma to be passed to SOAP
  integer  :: l_max = 0
  integer  :: n_max = 0
  integer  :: pissnl_dimension !% actual dimension of fit. 
  real(dp), allocatable :: e0(:)     !% energy shift per atom
  real(dp), allocatable :: beta(:,:) !% pissnl coefficients. size is the maximum of the pissnl_dimensions over all species
  
  character(len=STRING_LENGTH) label

end type IPModel_LinearSOAP

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_LinearSOAP), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_LinearSOAP_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_LinearSOAP_Finalise
end interface Finalise

interface Print
  module procedure IPModel_LinearSOAP_Print
end interface Print

interface Calc
  module procedure IPModel_LinearSOAP_Calc
end interface Calc

contains

subroutine IPModel_LinearSOAP_Initialise_str(this, args_str, param_str)
  type(IPModel_LinearSOAP), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

! now initialise the potential
#ifndef HAVE_GAP
  call system_abort('IPModel_LinearSOAP_Initialise_str: must be compiled with HAVE_GAP')
#else

  call initialise(params)
  this%label = ''
  call param_register(params, 'label', '', this%label, help_string="XML label of the potential")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_LinearSOAP_Initialise_str args_str')) then
    call system_abort("IPModel_LinearSOAP_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_LinearSOAP_read_params_xml(this, param_str)
#endif
end subroutine IPModel_LinearSOAP_Initialise_str

subroutine IPModel_LinearSOAP_Finalise(this)
#ifdef HAVE_GAP
  type(IPModel_LinearSOAP), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if (allocated(this%beta)) deallocate(this%beta)
  if (allocated(this%e0)) deallocate(this%e0)

  this%n_types = 0
  this%label = ''

#endif
end subroutine IPModel_LinearSOAP_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator: this routine computes energy, forces and the virial.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_LinearSOAP_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_LinearSOAP), intent(inout) :: this
  type(Atoms), intent(inout) :: at
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial
  character(len=*), intent(in), optional      :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

#ifdef HAVE_GAP 
  integer i, atomtype, d

  type(Dictionary)                :: params
  logical :: has_atom_mask_name
  character(STRING_LENGTH) :: atom_mask_name
  real(dp) :: r_scale, E_scale, e_i
  logical :: do_rescale_r, do_rescale_E

  logical, dimension(:), allocatable :: mpi_local_mask
  real(dp), dimension(:), allocatable   :: local_e_in
  logical, dimension(:), pointer :: atom_mask_pointer

  type(descriptor), dimension(:), allocatable :: my_descriptor
  type(descriptor_data) :: my_descriptor_data
  type(extendable_str) :: Zstring
  type(extendable_str) :: my_args_str
  

  
  INIT_ERROR(error)

  if (present(e)) e = 0.0_dp
  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_LinearSOAP_Calc', error)
     local_e = 0.0_dp
  endif
  if (present(f)) then 
     call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_LinearSOAP_Calc', error)
     f = 0.0_dp
     RAISE_ERROR("IPModel_LinearSOAP_Calc: force calculation requested but not supported yet.", error)
  end if
  if (present(virial)) then
     virial = 0.0_dp
     RAISE_ERROR("IPModel_LinearSOAP_Calc: virial calculation requested but not supported yet.", error)
  end if
  if (present(local_virial)) then
     call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_LinearSOAP_Calc', error)
     local_virial = 0.0_dp
     RAISE_ERROR("IPModel_LinearSOAP_Calc: local_virial calculation requested but not supported yet.", error)
  endif


  allocate(local_e_in(at%N))
  local_e_in = 0.0_dp

  atom_mask_pointer => null()
  has_atom_mask_name = .false.
  atom_mask_name = ""

  
  if (present(args_str)) then
    call initialise(params)
    call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
    call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")

    if(.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_LinearSOAP_Calc args_str')) then
       RAISE_ERROR("IPModel_LinearSOAP_Calc failed to parse args_str='"//trim(args_str)//"'",error)
    endif
    call finalise(params)
    
    if(has_atom_mask_name) then
       if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) &
        call system_abort("IPModel_GAP_Calc did not find "//trim(atom_mask_name)//" property in the atoms object.")
    else
        atom_mask_pointer => null()
    endif

    if (do_rescale_r .or. do_rescale_E) then
       RAISE_ERROR("IPModel_LinearSOAP_Calc: rescaling of potential with r_scale and E_scale not yet implemented!", error)
    end if

    my_args_str = trim(args_str)
  else
     call initialise(my_args_str)
  endif ! present(args_str)



  allocate(my_descriptor(this%n_types))

  
  ! compute descriptors for all atoms in the cell

  call initialise(Zstring)
  call concat(Zstring," species_Z={")
  do atomtype=1,this%n_types
     call concat(Zstring, " "//this%atomic_num(atomtype))
  end do
  call concat(Zstring, "} ")

  if( present(mpi) ) then
     if(mpi%active) then
        if(has_atom_mask_name) then
           RAISE_ERROR("IPModel_GAP: atom_mask_name "//trim(atom_mask_name)//" present while running MPI version. &
              The use of atom_mask_name is intended for serial-compiled code called from an external parallel code, such as LAMMPS",error)
        endif

        if( has_property(at,"mpi_local_mask") ) then
           RAISE_ERROR("IPModel_GAP: mpi_local_mask property already present", error)
        endif

        allocate(mpi_local_mask(at%N))
        call add_property_from_pointer(at,'mpi_local_mask',mpi_local_mask,error=error)

        call concat(my_args_str," atom_mask_name=mpi_local_mask")
     endif
  endif

  
  do atomtype=1,this%n_types
     call initialise(my_descriptor(atomtype),"soap cutoff="//this%cutoff//" n_max="//this%n_max//" l_max="//this%l_max//" atom_sigma="//this%atom_sigma//" n_species="//this%n_types//" Z="//this%atomic_num(atomtype)//Zstring//" "//string(my_args_str))
   
     if(mpi%active) call descriptor_MPI_setup(my_descriptor(atomtype),at,mpi,mpi_local_mask,error)

     d = descriptor_dimensions(my_descriptor(atomtype))

     ! compute descriptors for this atom type center
     call calc(my_descriptor(atomtype),at,my_descriptor_data, &
          do_descriptor=.true.,do_grad_descriptor=.false., &
          args_str=trim(string(my_args_str)), error=error)
     PASS_ERROR(error)

     ! this loop is over the found descriptors, i.e. all atoms of the current type
     do i=1,size(my_descriptor_data%x)

        !beta . pissnl
        e_i = sum(this%beta(:,atomtype) * my_descriptor_data%x(i)%data(:))
        
        local_e_in(my_descriptor_data%x(i)%ci(1)) = e_i
     end do

     call finalise(my_descriptor_data)
  end do

  !
  ! collect data and cleanup
  
  if(present(e)) e = sum(local_e_in)
  if(present(local_e)) local_e = local_e_in

  if (present(mpi)) then
     if( mpi%active ) then
        if (present(e)) e = sum(mpi, e)
        if (present(local_e)) call sum_in_place(mpi, local_e)
        if (present(virial)) call sum_in_place(mpi, virial)
        if (present(f)) call sum_in_place(mpi, f)

        call remove_property(at,'mpi_local_mask', error=error) 
        deallocate(mpi_local_mask)
     end if
  endif

  ! cleanup
  if(allocated(local_e_in)) deallocate(local_e_in)
  
#endif 
  
end subroutine IPModel_LinearSOAP_Calc



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions.
!% An example for XML stanza is given below, please notice that
!% they are simply dummy parameters for testing purposes, with no physical meaning.
!%
!%> <LinearSOAP_params n_types="2" label="default">
!%> <params l_max="12" n_max="12" cutoff="7.0" atom_sigma="0.5" pissnl_dimension="100"/>
!%> <per_type_data atomtype="1" atomic_num="29" e0="-10.1" />
!%> <per_type_data atomtype="2" atomic_num="79" e0="-12.2" />
!%> <beta atomtype="1" i="1" value="0.1" />
!%> <beta atomtype="1" i="2" value="0.1" />
!%> <beta atomtype="1" i="3" value="0.1" />
!%> </LinearSOAP_params>
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=1024) :: value

  integer :: ti, atomtype, i
  real(dp) :: v

  if (name == 'LinearSOAP_params') then ! new LinearSOAP stanza

    if (parse_in_ip) &
      call system_abort("IPModel_startElement_handler entered LinearSOAP_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

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
	call system_abort("Can't find n_types in LinearSOAP_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0
      allocate(parse_ip%e0(parse_ip%n_types))
      parse_ip%e0 = 0.0_dp

      if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
      allocate(parse_ip%type_of_atomic_num(108))
      parse_ip%type_of_atomic_num = 0
   endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "atomtype", value, status)
    if (status /= 0) call system_abort ("IPModel_LinearSOAP_read_params_xml cannot find atomtype")
    read (value, *) ti

    if (ti < 1) call system_abort("IPModel_LinearSOAP_read_params_xml got per_type_data type="//ti//" < 1")
    if (ti > parse_ip%n_types) call system_abort("IPModel_LinearSOAP_read_params_xml got per_type_data type="//ti//" > n_types="//parse_ip%n_types)

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_LinearSOAP_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti

    call QUIP_FoX_get_value(attributes, "e0", value, status)
    if(status /= 0)  call system_abort ("IPModel_LinearSOAP_read_params_xml cannot find e0")
    read (value, *) parse_ip%e0(ti)


  elseif (parse_in_ip .and. name == 'params') then

    call QUIP_FoX_get_value(attributes, "l_max", value, status)
    if (status /= 0) call system_abort ("IPModel_LinearSOAP_read_params_xml cannot find l_max")
    read (value, *) parse_ip%l_max

    call QUIP_FoX_get_value(attributes, "n_max", value, status)
    if (status /= 0) call system_abort ("IPModel_LinearSOAP_read_params_xml cannot find n_max")
    read (value, *) parse_ip%n_max

    call QUIP_FoX_get_value(attributes, "cutoff", value, status)
    if (status /= 0) call system_abort ("IPModel_LinearSOAP_read_params_xml cannot find cutoff")
    read (value, *) parse_ip%cutoff

    call QUIP_FoX_get_value(attributes, "atom_sigma", value, status)
    if (status /= 0) call system_abort ("IPModel_LinearSOAP_read_params_xml cannot find atom_sigma")
    read (value, *) parse_ip%atom_sigma

    call QUIP_FoX_get_value(attributes, "pissnl_dimension", value, status)
    if (status /= 0) call system_abort ("IPModel_LinearSOAP_read_params_xml cannot find pissnl_dimension")
    read (value, *) parse_ip%pissnl_dimension

    if(allocated(parse_ip%beta)) deallocate(parse_ip%beta)
    allocate(parse_ip%beta(parse_ip%pissnl_dimension, parse_ip%n_types))

  elseif (parse_in_ip .and. name == 'beta') then

     if(.not. allocated(parse_ip%beta)) call system_abort ("IPModel_LinearSOAP_read_params_xml encountered <beta> stanza but beta array is not yet allcoated")

     call QUIP_FoX_get_value(attributes, "atomtype", value, status)
     if (status /= 0) call system_abort ("IPModel_LinearSOAP_read_params_xml cannot find atomtype in <beta> stanza")
     read (value, *) atomtype

     call QUIP_FoX_get_value(attributes, "i", value, status)
     if (status /= 0) call system_abort ("IPModel_LinearSOAP_read_params_xml cannot find i in <beta> stanza")
     read (value, *) i

     call QUIP_FoX_get_value(attributes, "value", value, status)
     if (status /= 0) call system_abort ("IPModel_LinearSOAP_read_params_xml cannot find value in <beta> stanza")
     read (value, *) v

     if (i > parse_ip%pissnl_dimension) call system_abort("IPModel_LinearSOAP_read_params_xml: i > pissnl_dimension in <beta> stanza")
     if (atomtype > parse_ip%n_types) call system_abort("IPModel_LinearSOAP_read_params_xml: atomtype > n_types in <beta> stanza")
     
     parse_ip%beta(i,atomtype) = v

     
  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if (name == 'LinearSOAP_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_LinearSOAP_read_params_xml(this, param_str)
  type(IPModel_LinearSOAP), intent(inout), target :: this
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
    call system_abort("IPModel_LinearSOAP_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_LinearSOAP_read_params_xml


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of LinearSOAP parameters: number of different types, cutoff radius, atomic numbers, etc.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_LinearSOAP_Print (this, file)
  type(IPModel_LinearSOAP), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti
  
  call Print("IPModel_LinearSOAP ", file=file)
  call Print("IPModel_LinearSOAP : n_types = " // this%n_types // " cutoff = " // this%cutoff // " l_max = " // this%l_max // " n_max = " // this%n_max // " atom_sigma = " // this%atom_sigma // " pissnl_dimension = " // this%pissnl_dimension, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_LinearSOAP : type " // ti // " atomic_num " // this%atomic_num(ti) // " e0 = " // this%e0(ti), file=file)
  end do

end subroutine IPModel_LinearSOAP_Print


end module IPModel_LinearSOAP_module
