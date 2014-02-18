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
!X IPModel_Coulomb
!X
!% Coulomb module. Calculates electrostatic interactions between charged
!% species. Supported methods:
!% Direct: the $1/r$ potential
!% Yukawa: Yukawa-screened electrostatic interactions
!% Ewald: Ewald summation technique
!% DSF: Damped Shifted Force Coulomb potential. The interaction is damped by the
!% error function and the potential is force-shifted so both the potential and its
!% derivative goes smoothly to zero at the cutoff. Reference: JCP, 124, 234104 (2006)
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module IPModel_Coulomb_module

use error_module
use system_module, only : dp, inoutput, print, lower_case, verbosity_push_decrement, verbosity_pop, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module

use mpi_context_module
use QUIP_Common_module

use Yukawa_module
use IPEwald_module
use Ewald_module

implicit none
private

include 'IPModel_interface.h'

integer, parameter :: IPCoulomb_Method_Direct = 1
integer, parameter :: IPCoulomb_Method_Yukawa = 2
integer, parameter :: IPCoulomb_Method_Ewald  = 3
integer, parameter :: IPCoulomb_Method_DSF    = 4
integer, parameter :: IPCoulomb_Method_Ewald_NB  = 5
integer, parameter :: IPCoulomb_Method_Multipole_Moments = 6

integer, parameter :: Multipole_Position_Atomic = 1
integer, parameter :: Multipole_Position_Centre_of_Mass = 2

integer, parameter :: Multipole_Moments_Method_Fixed = 1
integer, parameter :: Multipole_Moments_Method_GAP = 2
integer, parameter :: Multipole_Moments_Method_Partridge_Schwenke = 3

public :: IPModel_Coulomb
type IPModel_Coulomb

  integer :: n_types = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)
  real(dp), dimension(:), allocatable :: charge
  integer :: method = 0

  real(dp) :: cutoff = 0.0_dp

  real(dp) :: yukawa_alpha = 0.0_dp
  real(dp) :: yukawa_smooth_length = 0.0_dp
  real(dp) :: yukawa_grid_size = 0.0_dp
  logical :: yukawa_pseudise = .false.

  real(dp) :: ewald_error
  real(dp) :: smooth_coulomb_cutoff

  real(dp) :: dsf_alpha = 0.0_dp

  character(len=STRING_LENGTH) :: label

endtype IPModel_Coulomb

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_Coulomb), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_Coulomb_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Coulomb_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Coulomb_Print
end interface Print

interface Calc
  module procedure IPModel_Coulomb_Calc
end interface Calc

contains

subroutine IPModel_Coulomb_Initialise_str(this, args_str, param_str)
  type(IPModel_Coulomb), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params
  character(len=STRING_LENGTH) :: method_str

  call Finalise(this)

  call initialise(params)
  this%label=''
  method_str=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'method', '', method_str, help_string="If present, method for Coulomb calculation.  Will be overridden by xml parameters if present")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_Coulomb_Initialise_str args_str')) then
    call system_abort("IPModel_Coulomb_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  if (trim(method_str) /= "") then
      select case(lower_case(trim(method_str)))
	 case("direct")
	    this%method = IPCoulomb_Method_Direct
	 case("yukawa")
	    this%method = IPCoulomb_Method_Yukawa
	 case("ewald")
	    this%method = IPCoulomb_Method_Ewald
	 case("ewald_nb")
	    this%method = IPCoulomb_Method_Ewald_NB
	 case("dsf")
	    this%method = IPCoulomb_Method_DSF
	 case("multipole_moments")
	    this%method = IPCoulomb_Method_Multipole_Moments
	 case default
	    call system_abort ("IPModel_Coulomb_Initialise_str: method "//trim(method_str)//" unknown")
      end select
  end if

  call IPModel_Coulomb_read_params_xml(this, param_str)

  this%multipoles%initialised=.false.  
  if (this%multipoles%n_monomer_types /= 0) then
    ! check for ambiguities / errors in site specification of the monomers 
     call print("specification of multipoles based on chemical environment")
     this%multipoles%initialised=.true.
  end if
  !  Add initialisation code here

end subroutine IPModel_Coulomb_Initialise_str

subroutine IPModel_Coulomb_Finalise(this)
  type(IPModel_Coulomb), intent(inout) :: this

  ! Add finalisation code here

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%charge)) deallocate(this%charge)

  this%n_types = 0
  this%label = ''
  this%cutoff = 0.0_dp
  this%method = 0

end subroutine IPModel_Coulomb_Finalise

recursive subroutine IPModel_Coulomb_Calc(this, at_in, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_Coulomb), intent(inout):: this
   type(Atoms), intent(inout),target  :: at_in
   type(Atoms), pointer               :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   type(Dictionary) :: params
   real(dp), dimension(:), allocatable, target :: my_charge
   real(dp), dimension(:), pointer :: charge
   real(dp), dimension(:,:), allocatable :: dummy_force
   real(dp) :: r_scale, E_scale, e_pre_calc
   logical :: do_rescale_r, do_rescale_E, do_pairwise_by_Z,do_e, do_f, intermolecular_only

   real(dp), pointer :: local_e_by_Z(:,:), local_e_contrib(:)
   integer, allocatable :: Z_s(:), Z_u(:)
   integer :: n_uniq_Zs
   type(Atoms),target :: at_dummy

   real(dp), allocatable :: gamma_mat(:,:)

   integer :: i, i_Z, intramolecular_factor

   character(len=STRING_LENGTH) :: charge_property_name, atom_mask_name, source_mask_name

   INIT_ERROR(error)
   do_f=present(f)
   do_e=present(e)

   at => at_in
   if (.not. has_property(at,"dummy_charge")) then
     if (do_e) e = 0.0_dp
     if (present(local_e)) then
        call check_size('Local_E',local_e,(/at%N/),'IPModel_Coulomb_Calc', error)
        local_e = 0.0_dp
     endif
     if (do_f) then
        call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_Coulomb_Calc', error)
        f = 0.0_dp
     end if
     if (present(virial)) virial = 0.0_dp
     if (present(local_virial)) then
        call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_Coulomb_Calc', error)
        local_virial = 0.0_dp
     endif
   else
     if (do_f) then
        call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_Coulomb_Calc', error)
     end if
   end if

   if (present(args_str)) then
      call initialise(params)
      call param_register(params, 'charge_property_name', 'charge', charge_property_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'atom_mask_name', '', atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'source_mask_name', '', source_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Rescaling factor for distances. Default 1.0.")
      call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Rescaling factor for energy. Default 1.0.")
      call param_register(params, 'pairwise_by_Z', 'F',do_pairwise_by_Z, help_string="If true, calculate pairwise contributions to local_e broken down by Z")
      call param_register(params, 'intermolecular_only', 'F',intermolecular_only, help_string="If true, ignore interactions between multipoles on the same molecule. Default F")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_Coulomb_Calc args_str')) then
         RAISE_ERROR("IPModel_Coulomb_Calc failed to parse args_str="//trim(args_str), error)
      endif
      call finalise(params)
      if (do_rescale_r .or. do_rescale_E) then
         RAISE_ERROR("IPModel_Coulomb_Calc: rescaling of potential with r_scale and E_scale not yet implemented!", error)
      end if
   else
      charge_property_name = 'charge'
   endif

   intramolecular_factor = 0
   if (intermolecular_only) then ! figure out whether pre_calc needs to subtract(-1), or ignore(0) interactions between sites on same molecule
     intramolecular_factor = -1 
   end if

   if(has_property(at,charge_property_name)) then
      if(.not. assign_pointer(at, charge_property_name, charge)) then
         RAISE_ERROR('IPModel_Coulomb_Calc failed to assign pointer to '//trim(charge_property_name)//' property', error)
      endif
   else if (this%multipoles%initialised) then
      call Multipole_Moments_Pre_Calc(at,this%multipoles,intramolecular_factor,do_f=do_f,e=e,error=error) ! specify positions and values of multipole moments, their derivatives w/r/t atomic positions
      if (do_e) e_pre_calc = e
      call Multipole_Moments_Make_Dummy_Atoms(at_dummy,at,this%multipoles,error) ! if multipoles are all charges then we can use other summation method by making a dummy atoms object
      at => at_dummy
      if (do_f) then                                                             
        if (this%method == IPCoulomb_Method_Multipole_Moments) then
          call Multipole_Moments_Calc(at, this%multipoles,this%multipoles%intermolecular_only, e=e, local_e=local_e, f=f, virial=virial,local_virial=local_virial, error = error)               
        else
          allocate(dummy_force(3,at%N))                                            ! on recursive call "dummy_charge" is present so this whole block gets skipped 
          dummy_force=0.0_dp
          call IPModel_Coulomb_Calc(this, at, e=e, f=dummy_force, virial=virial, args_str=args_str//" charge_property_name = dummy_charge", mpi=mpi, error=error)                                      
          call Multipole_Moments_Forces_From_Dummy_Atoms(at,this%multipoles,dummy_force,f,error)
          deallocate(dummy_force)
        end if
      else
        call IPModel_Coulomb_Calc(this, at, e=e, virial=virial, args_str=args_str//" charge_property_name = dummy_charge", mpi=mpi, error=error)
      end if
      if (do_e) e = e + e_pre_calc 
      at => at_in 
      return
   else
      allocate(my_charge(at%N))
      charge => my_charge
      charge = 0.0_dp
      do i = 1, at%N
         charge(i) = this%charge(this%type_of_atomic_num(at%Z(i))) 
      enddo
   endif

   selectcase(this%method)
   case(IPCoulomb_Method_Direct)
      call Direct_Coulomb_calc(at, charge, e=e, f=f, virial=virial, error = error)
   case(IPCoulomb_Method_Yukawa)
      call yukawa_charges(at, charge, this%cutoff, this%yukawa_alpha, this%yukawa_smooth_length, &
      e, local_e, f, virial, &
      mpi=mpi, atom_mask_name=atom_mask_name, source_mask_name=source_mask_name, type_of_atomic_num=this%type_of_atomic_num, &
      pseudise=this%yukawa_pseudise, grid_size=this%yukawa_grid_size, error=error)
   case(IPCoulomb_Method_Ewald)
      call Ewald_calc(at, charge, e, f, virial, ewald_error=this%ewald_error, use_ewald_cutoff=.false., smooth_coulomb_cutoff=this%smooth_coulomb_cutoff, error=error)
   case(IPCoulomb_Method_Ewald_NB)
      if (present(f) .or. present(virial) .or. present(local_virial)) then
	 RAISE_ERROR("IPModel_Coulomb_Calc: method ewald_nb doesn't have F or V implemented yet", error)
      endif
      allocate(gamma_mat(at%n,at%n))
      gamma_mat = 0.0_dp
      call add_madelung_matrix(at%N, at%lattice(:,1), at%lattice(:,2), at%lattice(:,3), at%pos, gamma_mat, redo_lattice=.true.)
call print("gamma_mat")
call print(gamma_mat)
      if (present(e)) e = 0.5_dp*sum(charge*matmul(gamma_mat, charge))
      if (present(local_e)) local_e = 0.5_dp*charge*matmul(gamma_mat, charge)
      if (do_pairwise_by_Z) then
	 allocate(Z_s(at%n))
	 Z_s = at%Z
	 call sort_array(Z_s)
	 call uniq(Z_s, Z_u)
	 deallocate(Z_s)

call print("local_pot "//matmul(gamma_mat, charge))

	 n_uniq_Zs = size(Z_u)
	 if (has_property(at, "local_e_pairwise_by_Z")) call remove_property(at, "local_e_pairwise_by_Z")
	 call add_property(at, "local_e_pairwise_by_Z", 0.0_dp, n_cols=n_uniq_Zs, ptr2 = local_e_by_Z)
	 allocate(local_e_contrib(at%N))
	 do i=1, at%N
	    local_e_contrib = 0.5_dp * gamma_mat(i,:) * charge(:) * charge(i)
call print("local_e_contrib "//i //" "//local_e_contrib)
	    do i_Z = 1, n_uniq_Zs
	       local_e_by_Z(i_Z,i) = sum(local_e_contrib, mask=(at%Z == Z_u(i_Z)))
	    end do
	 end do
	 deallocate(local_e_contrib)

      endif
      deallocate(gamma_mat)
   case(IPCoulomb_Method_DSF)
      call DSF_Coulomb_calc(at, charge, this%DSF_alpha, e=e, local_e=local_e, f=f, virial=virial, cutoff=this%cutoff, error = error)

   case(IPCoulomb_Method_Multipole_Moments)
      RAISE_ERROR("IPModel_Coulomb_Calc : something went wrong, you shouldn't have got here ",error)
   case default
      RAISE_ERROR("IPModel_Coulomb_Calc: unknown method", error)
   endselect

   charge => null()
   if(allocated(my_charge)) deallocate(my_charge)

end subroutine IPModel_Coulomb_Calc


subroutine IPModel_Coulomb_Print(this, file)
  type(IPModel_Coulomb), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti , tj

  call Print("IPModel_Coulomb : Coulomb Potential", file=file)
  select case(this%method)
  case(IPCoulomb_Method_Direct)
     call Print("IPModel_Coulomb method: Direct")
  case(IPCoulomb_Method_Yukawa)
     call Print("IPModel_Coulomb method: Yukawa")
  case(IPCoulomb_Method_Ewald)
     call Print("IPModel_Coulomb method: Ewald")
  case(IPCoulomb_Method_Ewald_NB)
     call Print("IPModel_Coulomb method: Ewald_NB")
  case(IPCoulomb_Method_DSF)
     call Print("IPModel_Coulomb method: Damped Shifted Force Coulomb")
  case(IPCoulomb_Method_Multipole_Moments)
     call Print("IPModel_Coulomb method: Direct sum over multipole interactions")
  case default
     call system_abort ("IPModel_Coulomb: method identifier "//this%method//" unknown")
  endselect

  call Print("IPModel_Coulomb : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_Coulomb : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    call Print ("IPModel_Coulomb : " // &
        "cutoff " // this%cutoff , &
         file=file)
    call verbosity_pop()
  end do
  if (this%multipoles%initialised) then
    do ti=1,this%multipoles%n_monomer_types
      call print("IPModel_Coulomb : monomer " // ti // " signature " // this%multipoles%monomer_types(ti)%signature//" moments method : "//this%multipoles%monomer_types(ti)%moments_method, file=file)
      do tj=1,size(this%multipoles%monomer_types(ti)%site_types)
        call print("IPModel_Coulomb :   site " // tj // " pos_type " // this%multipoles%monomer_types(ti)%site_types(tj)%pos_type, file=file)
      end do
    end do
  end if

end subroutine IPModel_Coulomb_Print

subroutine IPModel_Coulomb_read_params_xml(this, param_str)
  type(IPModel_Coulomb), intent(inout), target :: this
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
    call system_abort("IPModel_Coulomb_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_Coulomb_read_params_xml

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

  integer ti

  if (name == 'Coulomb_params') then ! new Coulomb stanza

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
        call system_abort("Can't find n_types in Coulomb_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      allocate(parse_ip%charge(parse_ip%n_types))
      parse_ip%charge = 0.0_dp

      call QUIP_FoX_get_value(attributes, 'n_monomer_types', value, status)
      if (status /= 0) then
        value='0'
      endif
      read (value, *), parse_ip%multipoles%n_monomer_types
      allocate(parse_ip%multipoles%monomer_types(parse_ip%multipoles%n_monomer_types))

      

      call QUIP_FoX_get_value(attributes, "cutoff", value, status)
      if (status /= 0) call system_abort ("IPModel_Coulomb_read_params_xml cannot find cutoff")
      read (value, *) parse_ip%cutoff

      call QUIP_FoX_get_value(attributes, "method", value, status)
      if (status /= 0) call system_abort ("IPModel_Coulomb_read_params_xml cannot find method")
      select case(lower_case(trim(value)))
      case("direct")
         parse_ip%method = IPCoulomb_Method_Direct
      case("yukawa")
         parse_ip%method = IPCoulomb_Method_Yukawa
      case("ewald")
         parse_ip%method = IPCoulomb_Method_Ewald
      case("ewald_nb")
         parse_ip%method = IPCoulomb_Method_Ewald_NB
      case("dsf")
         parse_ip%method = IPCoulomb_Method_DSF
      case("multipole_moments")
         parse_ip%method = IPCoulomb_Method_Multipole_Moments
      case default
         call system_abort ("IPModel_Coulomb_read_params_xml: method "//trim(value)//" unknown")
      endselect

      call QUIP_FoX_get_value(attributes, "yukawa_alpha", value, status)
      if (status /= 0) then
         if( parse_ip%method == IPCoulomb_Method_Yukawa ) call system_abort("IPModel_Coulomb_read_params_xml: Yukawa method requested but no yukawa_alpha parameter found.")
      else
         read (value, *) parse_ip%yukawa_alpha
      endif

      call QUIP_FoX_get_value(attributes, "yukawa_smooth_length", value, status)
      if (status /= 0) then
         if( parse_ip%method == IPCoulomb_Method_Yukawa ) call system_abort("IPModel_Coulomb_read_params_xml: Yukawa method requested but no yukawa_smooth_length parameter found.")
      else
         read (value, *) parse_ip%yukawa_smooth_length
      endif

      call QUIP_FoX_get_value(attributes, "yukawa_pseudise", value, status)
      if (status == 0) then
         read (value, *) parse_ip%yukawa_pseudise
      else
         parse_ip%yukawa_pseudise = .false.
      endif

      call QUIP_FoX_get_value(attributes, "yukawa_grid_size", value, status)
      if (status == 0) then
         read (value, *) parse_ip%yukawa_grid_size
      else
         parse_ip%yukawa_grid_size = 0.0_dp
      endif

      call QUIP_FoX_get_value(attributes, "ewald_error", value, status)
      if (status == 0) then
         read (value, *) parse_ip%ewald_error
      else
         parse_ip%ewald_error = 1.0e-6_dp
      endif

      call QUIP_FoX_get_value(attributes, "smooth_coulomb_cutoff", value, status)
      if (status == 0) then
         read (value, *) parse_ip%smooth_coulomb_cutoff
      else
         parse_ip%smooth_coulomb_cutoff = 0.0_dp
      endif

      call QUIP_FoX_get_value(attributes, "dsf_alpha", value, status)
      if (status /= 0) then
         if( parse_ip%method == IPCoulomb_Method_DSF ) call system_abort("IPModel_Coulomb_read_params_xml: Damped Shifted Force method requested but no dsf_alpha parameter found.")
      else
         read (value, *) parse_ip%dsf_alpha
      endif

    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_Coulomb_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_Coulomb_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    call QUIP_FoX_get_value(attributes, "charge", value, status)
    if (status /= 0) call system_abort ("IPModel_Coulomb_read_params_xml cannot find charge")
    read (value, *) parse_ip%charge(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
        parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  elseif (parse_in_ip .and. name == 'monomer') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_Coulomb_read_params_xml cannot find type")
    read (value, *) ti
    if (ti < 1) call system_abort("IPModel_Coulomb_read_params_xml got monomer type="//ti//" < 1")
    if (ti > parse_ip%multipoles%n_monomer_types) call system_abort("IPModel_Coulomb_read_params_xml got monomer type="//ti//" > n_types="//parse_ip%multipoles%n_monomer_types)

    call QUIP_FoX_get_value(attributes, "signature", value, status)
    if (status /= 0) call system_abort ("IPModel_Coulomb_read_params_xml cannot find signature")

    read (value, *) signature_string
    call split_string(signature_string,'_','{}',signature_fields(:),n_atoms,matching=.true.) !read the signature
    allocate(parse_ip%multipoles%monomer_types(ti)%signature(n_atoms))
    allocate(parse_ip%multipoles%monomer_types(ti)%masses(n_atoms))

    do i=1,n_atoms
      parse_ip%multipoles%monomer_types(ti)%signature(i) = string_to_int(signature_fields(i)) ! pass it to the monomer type
    end do

    call QUIP_FoX_get_value(attributes, "monomer_cutoff", value, status)
    if (size(parse_ip%multipoles%monomer_types(ti)%signature) > 1) then
      if (status /= 0 ) then
        call system_abort ("IPModel_Coulomb_read_params_xml cannot find monomer_cutoff")
      else
        read (value, *) parse_ip%multipoles%monomer_types(ti)%monomer_cutoff
      end if
    else
        parse_ip%multipoles%monomer_types(ti)%monomer_cutoff = 0.01_dp
    end if

    call QUIP_FoX_get_value(attributes, "n_sites", value, status)
    if (status /= 0) call system_abort ("IPModel_Coulomb_read_params_xml cannot find number of sites")
    allocate(parse_ip%multipoles%monomer_types(ti)%site_types(string_to_int(value)))

    call QUIP_FoX_get_value(attributes, "moments_method", value, status)
    if (status /= 0) call system_abort ("IPModel_Coulomb_read_params_xml cannot find moments_method")
    read (value, *) moments_method_str
    if (trim(moments_method_str) /= "") then
        select case(lower_case(trim(moments_method_str)))
  	 case("fixed")
           parse_ip%multipoles%monomer_types(ti)%moments_method = Multipole_Moments_Method_Fixed
  	 case("gap")
           parse_ip%multipoles%monomer_types(ti)%moments_method = Multipole_Moments_Method_GAP
  	 case("partridge_schwenke")
           parse_ip%multipoles%monomer_types(ti)%moments_method = Multipole_Moments_Method_Partridge_Schwenke
      case default
         call system_abort ("IPModel_Coulomb_read_params_xml: moments_method "//trim(value)//" unknown")
        end select
    end if

    call QUIP_FoX_get_value(attributes, "step", value, status)
    if (status /= 0) value = '1D-08'  
    read (value, *) parse_ip%multipoles%monomer_types(ti)%step

  elseif (parse_in_ip .and. name == 'per_site_data') then

    call QUIP_FoX_get_value(attributes, "monomer", value, status)
    if (status /= 0) call system_abort ("IPModel_Coulomb_read_params_xml cannot find which monomer this site belongs to, specify with monomer='n' n=0,1,2...")
    read (value, *) ti
    if (ti < 1) call system_abort("IPModel_Coulomb_read_params_xml got monomer type="//ti//" < 1")
    if (.not. allocated(parse_ip%multipoles%monomer_types(ti)%signature)) call system_abort("IPModel_Coulomb_read_params_xml, site specified but monomer not initialised")

    call QUIP_FoX_get_value(attributes, "site", value, status)
    if (status /= 0) call system_abort ("IPModel_Coulomb_read_params_xml cannot find site number (e.g. site='1') ")

    read (value, *) tj

    if (tj < 1) call system_abort("IPModel_Coulomb_read_params_xml got monomer type="//tj//" < 1")
    if (tj > size(parse_ip%multipoles%monomer_types(ti)%site_types)) call system_abort("IPModel_Coulomb_read_params_xml got site type="//tj//" > n_sites="//size(parse_ip%multipoles%monomer_types(ti)%sites))
    call QUIP_FoX_get_value(attributes, "charge", value, status)
    if (status == 0) then
      read (value, *) parse_ip%multipoles%monomer_types(ti)%site_types(tj)%charge    
      parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d = 1
    else
      parse_ip%multipoles%monomer_types(ti)%site_types(tj)%charge = 0.0_dp
    end if

    call QUIP_FoX_get_value(attributes, "pos_type", value, status)
    if (status /= 0) call system_abort ("IPModel_Coulomb_read_params_xml cannot find site pos_type ")    
    read (value, *) pos_type_str

    if (trim(pos_type_str) /= "") then
        select case(lower_case(trim(pos_type_str)))
  	 case("atom")
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%pos_type = Multipole_Position_Atomic
  	 case("com")
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%pos_type = Multipole_Position_Centre_of_Mass
         case default
           call system_abort ("IPModel_Coulomb_read_params_xml: pos_type "//trim(value)//" unknown")         
        end select
    end if

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) value='0'    
    read (value, *) parse_ip%multipoles%monomer_types(ti)%site_types(tj)%atomic_number

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'Coulomb_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_Coulomb_module
