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
!X IPModel_Multipoles
!X
!% Multipoles interatomic potential: 
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_Multipoles_module

use error_module
!use system_module, only : dp, inoutput, print, lower_case, verbosity_push_decrement, verbosity_pop, operator(//), split_string, string_to_int,PRINT_ANAL,reallocate
use system_module
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module
use Multipoles_module
use units_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

integer, parameter :: IPMultipoles_Method_Direct = 1
integer, parameter :: IPMultipoles_Method_Yukawa = 2
integer, parameter :: IPMultipoles_Method_Ewald  = 3

integer, parameter :: Multipole_Position_Atomic = 1                
integer, parameter :: Multipole_Position_Centre_of_Mass = 2        
integer, parameter :: Multipole_Position_M_Site = 3                
                                                                   
integer, parameter :: Charge_Method_None = 0                       
integer, parameter :: Charge_Method_Fixed = 1                      
integer, parameter :: Charge_Method_GAP = 2                        
integer, parameter :: Charge_Method_Partridge_Schwenke = 3         
                                                                   
integer, parameter :: Dipole_Method_None = 0                       
integer, parameter :: Dipole_Method_Partridge_Schwenke = 1         
integer, parameter :: Dipole_Method_GAP = 2                        
                                                                   
integer, parameter :: Polarisation_Method_None = 0                 
integer, parameter :: Polarisation_Method_FPI = 1                  
integer, parameter :: Polarisation_Method_GMRES = 2                
integer, parameter :: Polarisation_Method_QR = 3                   
                                                                   
integer, parameter :: Damping_None = 0                             
integer, parameter :: Damping_Exp = 1                              
integer, parameter :: Damping_Erf = 2                              
integer, parameter :: Damping_Erf_Uniform = 3                      
                                                                   
integer, parameter :: Screening_None = 0                           
integer, parameter :: Screening_Yukawa = 1                         
integer, parameter :: Screening_Erfc_Uniform = 2

public :: IPModel_Multipoles
type IPModel_Multipoles

  type(Multipole_Moments) :: multipoles

  integer :: n_monomer_types = 0

  integer :: method = 0
  integer :: damping = 0
  integer :: screening = 0
  integer :: polarisation = 0

  real(dp) :: cutoff = 0.0_dp

  real(dp) :: ewald_error
  real(dp) :: smooth_coulomb_cutoff

  character(len=STRING_LENGTH) :: label

end type IPModel_Multipoles

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_Multipoles), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_Multipoles_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Multipoles_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Multipoles_Print
end interface Print

interface Calc
  module procedure IPModel_Multipoles_Calc
end interface Calc

contains

subroutine IPModel_Multipoles_Initialise_str(this, args_str, param_str, error)
  type(IPModel_Multipoles), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(Dictionary) :: params
  integer, optional, intent(out):: error


  INIT_ERROR(error)
  call Finalise(this)
  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if(.not. param_read_line(params, args_str, ignore_unknown=.true., task='IPModel_Multipoles_Initialise args_str')) then
     RAISE_ERROR("IPModel_Multipoles_Init failed to parse args_str='"//trim(args_str)//"'", error)
  end if

  call finalise(params)

  call IPModel_Multipoles_read_params_xml(this, param_str)

end subroutine IPModel_Multipoles_Initialise_str

subroutine IPModel_Multipoles_Finalise(this)
  type(IPModel_Multipoles), intent(inout) :: this

  ! Add finalisation code here
  call finalise(this%multipoles)

end subroutine IPModel_Multipoles_Finalise


subroutine IPModel_Multipoles_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_Multipoles), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   type(Atoms)                            :: dummy_atoms
   real(dp),dimension(:,:), allocatable   :: pol_matrix
   type(Ewald_arrays)                     :: ewald
   type(Dictionary)                       :: params
   real(dp)                               :: r_scale, E_scale, pol_energy
   logical                                :: do_rescale_r, do_rescale_E,do_e, do_f, strict

   logical :: my_use_ewald_cutoff

   INIT_ERROR(error)
   do_f=present(f)
   do_e=present(e)
   if(do_e)e=0.0_dp
   if (do_f) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_Coulomb_Calc', error)
      f = 0.0_dp
   end if
   if (present(args_str)) then
      call initialise(params)
      call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Rescaling factor for distances. Not supported in multipole calcs.")
      call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Rescaling factor for energy. Default 1.0.")
      call param_register(params, 'strict', 'T',strict, help_string="If true, make sure every atom in system belongs to one of the monomers in this potential, default T")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_Multipoles_Calc args_str')) then
         RAISE_ERROR("IPModel_Multipoles_Calc failed to parse args_str="//trim(args_str), error)
      endif
      call finalise(params)
      if (do_rescale_r ) then
         RAISE_ERROR("IPModel_Multipoles_Calc: rescaling of potential with r_scale not yet implemented!", error)
      end if
   endif

   call multipole_sites_setup(at,this%multipoles,dummy_atoms,do_f,strict) ! multipoles includes exclude_list

   if (this%method == IPMultipoles_Method_Ewald) then
     call ewald_setup(dummy_atoms,this%multipoles,ewald,this%ewald_error) ! changes this%multipoles%cutoff
   end if

   if (this%polarisation /= Polarisation_Method_None) then
     call electrostatics_calc(dummy_atoms,this%multipoles,ewald,do_field=.true.)
     call build_polarisation_matrix(dummy_atoms,this%multipoles,pol_matrix,ewald) ! A^-1-T
     call calc_induced_dipoles(pol_matrix,this%multipoles,this%polarisation,pol_energy) ! this updates the dipoles on any polarisable sites
     if(do_e) e = pol_energy
   end if

   call electrostatics_calc(dummy_atoms,this%multipoles,ewald,e=e,do_force=.true.)
   if (present(f)) then
     call atomic_forces_from_sites(this%multipoles%sites,f=f) ! no support for virial or local stuff yet
   end if

   ! clean up
   if(allocated(pol_matrix)) deallocate(pol_matrix)
   call finalise(dummy_atoms)

end subroutine IPModel_Multipoles_Calc


subroutine IPModel_Multipoles_Print(this, file)
  type(IPModel_Multipoles), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti , tj

  call Print("IPModel_Multipoles : Multipoles Potential", file=file)
  select case(this%method)
  case(IPMultipoles_Method_Direct)
     call Print("IPModel_Multipoles method: Direct")
  case(IPMultipoles_Method_Yukawa)
     call Print("IPModel_Multipoles method: Yukawa")
  case(IPMultipoles_Method_Ewald)
     call Print("IPModel_Multipoles method: Ewald")
  case default
     call system_abort ("IPModel_Multipoles: method identifier "//this%method//" unknown")
  endselect

  call Print("IPModel_Multipoles : n_monomer_types = " // this%n_monomer_types // " cutoff = " // this%cutoff, file=file)

  call print("Electrostatics Options for multipole interactions",file=file)
  select case(this%multipoles%calc_opts%damping)
  case(Damping_None)
     call Print("IPModel_Multipoles damping :  No damping, bare electrostatic interactions at short range",file=file)
  case(Damping_Erf)
     call Print("IPModel_Multipoles damping :  Gaussian damping of electrostatic interactions at short range",file=file)
  case(Damping_Exp)
     call Print("IPModel_Multipoles damping :  Exponential damping of electrostatic interactions at short range",file=file)
     call Print("IPModel_Multipoles damping :  Exponential damping rank : "//this%multipoles%calc_opts%damp_exp_order,file=file)
     call Print("IPModel_Multipoles damping :  Exponential damping scale : "//this%multipoles%calc_opts%damp_exp_scale,file=file)
  case default
     call Print("IPModel_Multipoles damping :  Unknown damping scheme",file=file)
  endselect

  select case(this%multipoles%calc_opts%screening)
  case(Screening_None)
     call Print("IPModel_Multipoles screening :  No screening, including all long range electrostatic interactions",file=file)
  case(Screening_Yukawa)
     call Print("IPModel_Multipoles screening :  Yukawa screening of long range electrostatic interactions",file=file)
     call Print("IPModel_Multipoles screening :  yukawa alpha : "//this%multipoles%calc_opts%yukawa_alpha,file=file)
     call Print("IPModel_Multipoles screening :  yukawa smooth_length : "//this%multipoles%calc_opts%yukawa_smooth_length,file=file)
  case default
     call Print("IPModel_Multipoles screening :  Unknown screening scheme",file=file)
  endselect

  do ti=1,this%n_monomer_types
    call print("IPModel_Multipoles : monomer " // ti // " signature " // this%multipoles%monomer_types(ti)%signature, file=file)
    do tj=1,size(this%multipoles%monomer_types(ti)%site_types)
      call print("IPModel_Multipoles :   site " // tj // " d " // this%multipoles%monomer_types(ti)%site_types(tj)%d, file=file)
      call print("IPModel_Multipoles :   site " // tj //" charge method : "//this%multipoles%monomer_types(ti)%site_types(tj)%charge_method, file=file)
      call print("IPModel_Multipoles :   site " // tj //" dipole method : "//this%multipoles%monomer_types(ti)%site_types(tj)%dipole_method, file=file)
      call print("IPModel_Multipoles :   site " // tj // " pos_type " // this%multipoles%monomer_types(ti)%site_types(tj)%pos_type, file=file)
    end do
  end do


end subroutine IPModel_Multipoles_Print

subroutine IPModel_Multipoles_read_params_xml(this, param_str)
  type(IPModel_Multipoles), intent(inout), target :: this
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

  if (this%n_monomer_types == 0) then
    call system_abort("IPModel_Multipoles_read_params_xml parsed file, but n_monomer_types = 0")
  endif

end subroutine IPModel_Multipoles_read_params_xml

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

  character(len=STRING_LENGTH) :: value, signature_string, pos_type_str, moments_method_str
  character(len=STRING_LENGTH), dimension(99) :: signature_fields

  logical :: energy_shift, linear_force_shift
  integer :: ti, tj, i , j, n_atoms,n_sites, n_pairs,cursor
  real(dp) :: alpha_au

  if (name == 'Multipoles_params') then ! new Multipoles stanza

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

      if (parse_ip%n_monomer_types /= 0) then
        call finalise(parse_ip)
      endif

      call QUIP_FoX_get_value(attributes, 'n_monomer_types', value, status)
      if (status /= 0) then
        value='0'
      endif
      read (value, *), parse_ip%n_monomer_types

      allocate(parse_ip%multipoles%monomer_types(parse_ip%n_monomer_types))

      call QUIP_FoX_get_value(attributes, "cutoff", value, status)
      if (status /= 0) call system_abort ("IPModel_Multipoles_read_params_xml cannot find cutoff")
      read (value, *) parse_ip%cutoff
      parse_ip%multipoles%cutoff = parse_ip%cutoff

      call QUIP_FoX_get_value(attributes, 'dipole_tolerance', value, status)
      if (status /= 0) then
        value='0.0D0'
      endif
      read (value, *), parse_ip%multipoles%dipole_tolerance

      call QUIP_FoX_get_value(attributes, 'intermolecular_only', value, status)
      if (status == 0) then 
        read (value, *), parse_ip%multipoles%intermolecular_only
      endif


      call QUIP_FoX_get_value(attributes, "method", value, status)
      if (status /= 0) call system_abort ("IPModel_Multipoles_read_params_xml cannot find method")
      select case(lower_case(trim(value)))
      case("direct")
         parse_ip%method = IPMultipoles_Method_Direct
      case("yukawa")
         parse_ip%method = IPMultipoles_Method_Yukawa
      case("ewald")
         parse_ip%method = IPMultipoles_Method_Ewald
      case default
         call system_abort ("IPModel_Multipoles_read_params_xml: method "//trim(value)//" unknown")
      endselect

      call QUIP_FoX_get_value(attributes, "polarisation", value, status)
      if (status /= 0)  value = "none"
      select case(lower_case(trim(value)))
      case("fpi")
         parse_ip%polarisation = Polarisation_Method_FPI
      case("gmres")
         parse_ip%polarisation = Polarisation_Method_GMRES
      case("qr")
         parse_ip%polarisation = Polarisation_Method_QR
      case("none")
         parse_ip%polarisation = Polarisation_Method_None
      case default
         parse_ip%polarisation = Polarisation_Method_None
      endselect

      call QUIP_FoX_get_value(attributes, "damping", value, status)
      if (status /= 0)  value = "none"
      select case(lower_case(trim(value)))
      case("exp")
         parse_ip%multipoles%calc_opts%damping = Damping_Exp
      case("erf")
         parse_ip%multipoles%calc_opts%damping = Damping_Erf
      case("none")
         parse_ip%multipoles%calc_opts%damping = Damping_None
      case default
         parse_ip%multipoles%calc_opts%damping = Damping_None
      endselect


      call QUIP_FoX_get_value(attributes, "damp_exp_scale", value, status) ! global param to adjust strength of damping, just a convenience param for damping_exp
      if (status == 0) then
         if (parse_ip%multipoles%calc_opts%damping == Damping_Exp) then
           read (value, *) parse_ip%multipoles%calc_opts%damp_exp_scale
         else
           call system_abort("IPModel_Multipoles_read_params_xml: damp_exp_scale specified but exponential damping not selected")
         end if
      else
         parse_ip%multipoles%calc_opts%damp_exp_scale = 1.0_dp
      endif

      call QUIP_FoX_get_value(attributes, "damp_exp_order", value, status) ! order of exponential damping, typically 3 or 4, default 3 like in FX water model
      if (status == 0) then
         if (parse_ip%multipoles%calc_opts%damping == Damping_Exp) then
           read (value, *)  parse_ip%multipoles%calc_opts%damp_exp_order
         else
           call system_abort("IPModel_Multipoles_read_params_xml: damp_exp_order specified but exponential damping not selected")
         end if
      else
         parse_ip%multipoles%calc_opts%damp_exp_order = 3
      endif

      call QUIP_FoX_get_value(attributes, "yukawa_alpha", value, status)
      if (status /= 0) then
         if( parse_ip%method == IPMultipoles_Method_Yukawa ) call system_abort("IPModel_Multipoles_read_params_xml: Yukawa method requested but no yukawa_alpha parameter found.")
      else
         read (value, *) parse_ip%multipoles%calc_opts%yukawa_alpha
         parse_ip%screening = Screening_Yukawa
         parse_ip%multipoles%calc_opts%screening = Screening_Yukawa
      endif

      call QUIP_FoX_get_value(attributes, "yukawa_smooth_length", value, status)
      if (status /= 0) then
         if( parse_ip%screening == Screening_Yukawa ) call system_abort("IPModel_Multipoles_read_params_xml: Yukawa method requested but no yukawa_smooth_length parameter found.")
      else
         read (value, *) parse_ip%multipoles%calc_opts%yukawa_smooth_length
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

    endif

  elseif (parse_in_ip .and. name == 'monomer') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_Multipoles_read_params_xml cannot find type")
    read (value, *) ti
    if (ti < 1) call system_abort("IPModel_Multipoles_read_params_xml got monomer type="//ti//" < 1")
    if (ti > parse_ip%n_monomer_types) call system_abort("IPModel_Multipoles_read_params_xml got monomer type="//ti//" > n_monomer_types="//parse_ip%n_monomer_types)

    call QUIP_FoX_get_value(attributes, "signature", value, status)
    if (status /= 0) call system_abort ("IPModel_Multipoles_read_params_xml cannot find signature")

    read (value, *) signature_string
    call split_string(signature_string,'_','{}',signature_fields(:),n_atoms,matching=.true.) !read the signature
    allocate(parse_ip%multipoles%monomer_types(ti)%signature(n_atoms))
    allocate(parse_ip%multipoles%monomer_types(ti)%masses(n_atoms))

    do i=1,n_atoms
      parse_ip%multipoles%monomer_types(ti)%signature(i) = string_to_int(signature_fields(i)) ! pass it to the monomer type
    end do

    !call QUIP_FoX_get_value(attributes, "exclude_pairs", value, status)

    call QUIP_FoX_get_value(attributes, "monomer_cutoff", value, status)
    if (size(parse_ip%multipoles%monomer_types(ti)%signature) > 1) then
      if (status /= 0 ) then
        call system_abort ("IPModel_Multipoles_read_params_xml cannot find monomer_cutoff")
      else
        read (value, *) parse_ip%multipoles%monomer_types(ti)%monomer_cutoff
      end if
    else
        parse_ip%multipoles%monomer_types(ti)%monomer_cutoff = 0.01_dp
    end if

    call QUIP_FoX_get_value(attributes, "n_sites", value, status)
    if (status /= 0) call system_abort ("IPModel_Multipoles_read_params_xml cannot find number of sites")
    n_sites=string_to_int(value)
    allocate(parse_ip%multipoles%monomer_types(ti)%site_types(n_sites))    
    if(parse_ip%multipoles%intermolecular_only) then
      n_pairs=(n_sites*(n_sites-1))/2
      call reallocate(parse_ip%multipoles%monomer_types(ti)%excluded_pairs,2,n_pairs,zero=.true.)
      cursor=0
      do i=1,n_sites-1
        do j=i+1,n_sites
          cursor=cursor+1
          parse_ip%multipoles%monomer_types(ti)%excluded_pairs(1,cursor)=i
          parse_ip%multipoles%monomer_types(ti)%excluded_pairs(2,cursor)=j
        end do
      end do          
    end if
      
    call QUIP_FoX_get_value(attributes, "step", value, status)
    if (status /= 0) value = '1D-08'  
    read (value, *) parse_ip%multipoles%monomer_types(ti)%step

    call QUIP_FoX_get_value(attributes, "gammaM", value, status)
    if (status == 0) read (value, *) parse_ip%multipoles%monomer_types(ti)%gammaM


  elseif (parse_in_ip .and. name == 'per_site_data') then

    call QUIP_FoX_get_value(attributes, "monomer", value, status)
    if (status /= 0) call system_abort ("IPModel_Multipoles_read_params_xml cannot find which monomer this site belongs to, specify with monomer='n' n=0,1,2...")
    read (value, *) ti
    if (ti < 1) call system_abort("IPModel_Multipoles_read_params_xml got monomer type="//ti//" < 1")
    if (.not. allocated(parse_ip%multipoles%monomer_types(ti)%signature)) call system_abort("IPModel_Multipoles_read_params_xml, site specified but monomer not initialised")

    call QUIP_FoX_get_value(attributes, "site", value, status)
    if (status /= 0) call system_abort ("IPModel_Multipoles_read_params_xml cannot find site number (e.g. site='1') ")

    read (value, *) tj

    if (tj < 1) call system_abort("IPModel_Multipoles_read_params_xml got site type="//tj//" < 1")
    if (tj > size(parse_ip%multipoles%monomer_types(ti)%site_types)) call system_abort("IPModel_Multipoles_read_params_xml got site type="//tj//" > n_sites="//size(parse_ip%multipoles%monomer_types(ti)%site_types))
    call QUIP_FoX_get_value(attributes, "charge", value, status)
    if (status == 0) then
      read (value, *) parse_ip%multipoles%monomer_types(ti)%site_types(tj)%charge    
      parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d = 1
    else
      parse_ip%multipoles%monomer_types(ti)%site_types(tj)%charge = 0.0_dp
    end if

    call QUIP_FoX_get_value(attributes, "pos_type", value, status)
    if (status /= 0) call system_abort ("IPModel_Multipoles_read_params_xml cannot find site pos_type ")    
    read (value, *) pos_type_str

    if (trim(pos_type_str) /= "") then
        select case(lower_case(trim(pos_type_str)))
  	 case("atom")
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%pos_type = Multipole_Position_Atomic
  	 case("com")
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%pos_type = Multipole_Position_Centre_of_Mass
  	 case("m-site")
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%pos_type = Multipole_Position_M_Site
         case default
           call system_abort ("IPModel_Multipoles_read_params_xml: pos_type "//trim(value)//" unknown")         
        end select
    end if

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) value='0'    
    read (value, *) parse_ip%multipoles%monomer_types(ti)%site_types(tj)%atomic_number

    parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d = 0

    parse_ip%multipoles%monomer_types(ti)%site_types(tj)%charge_method = Charge_Method_None
    call QUIP_FoX_get_value(attributes, "charge_method", value, status)
    if (status == 0) then
      read (value, *) moments_method_str
      select case(lower_case(trim(moments_method_str)))
         case("none")
           continue
  	 case("fixed")
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%charge_method = Charge_Method_Fixed
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d =  parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d + 1
  	 case("gap")
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%charge_method = Charge_Method_GAP
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d =  parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d + 1
  	 case("partridge_schwenke")
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%charge_method = Charge_Method_Partridge_Schwenke
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d =  parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d + 1
      case default
         call system_abort ("IPModel_Multipoles_read_params_xml: moments_method "//trim(value)//" unknown")
        end select
    end if

    call QUIP_FoX_get_value(attributes, "dipole_method", value, status)
    parse_ip%multipoles%monomer_types(ti)%site_types(tj)%dipole_method = Dipole_Method_None
    if (status == 0) then
      read (value, *) moments_method_str
      select case(lower_case(trim(moments_method_str)))
         case("none")
           continue
  	 case("gap")
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%dipole_method = Dipole_Method_GAP
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d =  parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d + 3
  	 case("partridge_schwenke")
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%dipole_method = Dipole_Method_Partridge_Schwenke
           parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d =  parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d + 3
      case default
         call system_abort ("IPModel_Multipoles_read_params_xml: moments_method "//trim(value)//" unknown")
        end select
    end if

    call QUIP_FoX_get_value(attributes, "pol_alpha", value, status)
    if (status /= 0) then 
      parse_ip%multipoles%monomer_types(ti)%site_types(tj)%polarisable=.false.
    else
      parse_ip%multipoles%monomer_types(ti)%site_types(tj)%polarisable=.true.
      if (parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d .lt. 3 ) parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d = parse_ip%multipoles%monomer_types(ti)%site_types(tj)%d + 3
    end if
    read (value, *) alpha_au
    parse_ip%multipoles%monomer_types(ti)%site_types(tj)%alpha=alpha_au*CUBIC_BOHR

    ! by default, polarisable sites are damped and unpolarisable ones are not, can override default behaviour by setting damp_rad
    ! if damp_rad=0.0, then no damping occurs. If want to damp a non-polarisable site, have to provide a radius manually
    parse_ip%multipoles%monomer_types(ti)%site_types(tj)%damped = parse_ip%multipoles%monomer_types(ti)%site_types(tj)%polarisable

    call QUIP_FoX_get_value(attributes, "damp_rad", value, status)
    if (status /= 0) then 
      if (parse_ip%multipoles%calc_opts%damping==Damping_Erf) then
        parse_ip%multipoles%monomer_types(ti)%site_types(tj)%damp_rad = (alpha_au*sqrt(2.0_dp/PI)/3.0_dp)**(1.0_dp/3.0_dp)*BOHR
      else if (parse_ip%multipoles%calc_opts%damping==Damping_Exp) then
        parse_ip%multipoles%monomer_types(ti)%site_types(tj)%damp_rad = (alpha_au)**(1.0_dp/3.0_dp)*BOHR
      end if
    else
      if (parse_ip%multipoles%calc_opts%damping==Damping_None) then
        call system_abort("damp_rad specified but no damping type selected, please specify damping=exp or damping=erf")
      else
        read (value, *) parse_ip%multipoles%monomer_types(ti)%site_types(tj)%damp_rad ! override default damp_rad value
        parse_ip%multipoles%monomer_types(ti)%site_types(tj)%damped = .true.
      end if
    end if

    if (.not. parse_ip%multipoles%monomer_types(ti)%site_types(tj)%damped .and. parse_ip%multipoles%calc_opts%damping/=Damping_None ) then
      call system_abort ("IPModel_Multipoles_read_params_xml: damping requested but no radii specified for non-polarisable sites ")
    end if

    parse_ip%multipoles%monomer_types(ti)%site_types(tj)%initialised = .true.

  endif


end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'Multipoles_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler


end module IPModel_Multipoles_module
