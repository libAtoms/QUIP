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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X TBModel_Bowler module
!X
!% Calculate energies using Bowler tight-binding model
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module TBModel_Bowler_module

use system_module, only : dp, print, inoutput, system_abort, verbosity_push_decrement, verbosity_pop, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use spline_module
use atoms_module

use TB_Common_module
use QUIP_Common_module

implicit none
private

integer, parameter :: max_n_orb_sets = 2

include 'TBModel_interface.h'

public :: TBModel_Bowler
type TBModel_Bowler
  integer :: n_types = 0
  character(len=STRING_LENGTH) label

  real(dp) :: cutoff = 0.0_dp
  logical :: is_orthogonal = .true.

  integer, allocatable :: type_of_atomic_num(:)
  integer, allocatable :: n_orbs(:), n_elecs(:), n_orb_sets(:), orb_set_type(:,:)
  integer, allocatable :: atomic_num(:)

  !! Bowler parameters
  real(dp), allocatable :: H_coeff(:,:,:), Vrep(:,:)
  real(dp), allocatable :: r0(:,:), rc(:,:), n(:,:), nc(:,:), dc(:,:), m(:,:), mc(:,:)
  real(dp) :: tailx0

  real(dp), allocatable :: E(:,:)

  type(spline), allocatable :: H_tail_spline(:,:)
  type(spline), allocatable :: Vrep_tail_spline(:,:)

  logical :: has_default_fermi_E = .false., has_default_fermi_T = .true., has_default_band_width = .false., has_default_k_density=.false.
  real(dp) :: default_fermi_E, default_fermi_T = 0.001_dp, default_band_width, default_k_density

end type

integer, private :: parse_cur_type
logical, private :: parse_in_tbm, parse_matched_label
type(TBModel_Bowler), private, pointer :: parse_tbm

interface Initialise
  module procedure TBModel_Bowler_Initialise_str
end interface Initialise

interface Finalise
  module procedure TBModel_Bowler_Finalise
end interface Finalise

interface Print
  module procedure TBModel_Bowler_Print
end interface Print

interface n_orbs_of_Z
  module procedure TBModel_Bowler_n_orbs_of_Z
end interface n_orbs_of_Z

interface n_orb_sets_of_Z
  module procedure TBModel_Bowler_n_orb_sets_of_Z
end interface n_orb_sets_of_Z

interface n_orbs_of_orb_set_of_Z
  module procedure TBModel_Bowler_n_orbs_of_orb_set_of_Z
end interface n_orbs_of_orb_set_of_Z

interface orb_type_of_orb_set_of_Z
  module procedure TBModel_Bowler_orb_type_of_orb_set_of_Z
end interface orb_type_of_orb_set_of_Z

interface n_elecs_of_Z
  module procedure TBModel_Bowler_n_elecs_of_Z
end interface n_elecs_of_Z

interface get_HS_blocks
  module procedure TBModel_Bowler_get_HS_blocks
end interface get_HS_blocks

interface get_dHS_masks
  module procedure TBModel_Bowler_get_dHS_masks
end interface get_dHS_masks

interface get_dHS_blocks
  module procedure TBModel_Bowler_get_dHS_blocks
end interface get_dHS_blocks

interface get_local_rep_E
  module procedure TBModel_Bowler_get_local_rep_E
end interface get_local_rep_E

interface get_local_rep_E_force
  module procedure TBModel_Bowler_get_local_rep_E_force
end interface get_local_rep_E_force

interface get_local_rep_E_virial
  module procedure TBModel_Bowler_get_local_rep_E_virial
end interface get_local_rep_E_virial

interface calc_H_coeff
  module procedure TBModel_Bowler_calc_H_coeff
end interface calc_H_coeff

interface calc_H_coeff_deriv
  module procedure TBModel_Bowler_calc_H_coeff_deriv
end interface calc_H_coeff_deriv

contains

subroutine TBModel_Bowler_Initialise_str(this, args_str, param_str)
  type(TBModel_Bowler), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='TBModel_Bowler_Initialise_str args_str')) then
    call system_abort("TBModel_Bowler_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call TBModel_Bowler_read_params_xml(this, param_str)
end subroutine TBModel_Bowler_Initialise_str

subroutine TBModel_Bowler_Finalise(this)
  type(TBModel_Bowler), intent(inout) :: this

  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%n_orbs)) deallocate(this%n_orbs)
  if (allocated(this%n_elecs)) deallocate(this%n_elecs)
  if (allocated(this%n_orb_sets)) deallocate(this%n_orb_sets)
  if (allocated(this%orb_set_type)) deallocate(this%orb_set_type)
  if (allocated(this%atomic_num)) deallocate(this%atomic_num)

  if (allocated(this%H_coeff)) deallocate(this%H_coeff)
  if (allocated(this%Vrep)) deallocate(this%Vrep)

  if (allocated(this%r0)) deallocate(this%r0)
  if (allocated(this%rc)) deallocate(this%rc)
  if (allocated(this%n)) deallocate(this%n)
  if (allocated(this%nc)) deallocate(this%nc)
  if (allocated(this%dc)) deallocate(this%dc)
  if (allocated(this%m)) deallocate(this%m)
  if (allocated(this%mc)) deallocate(this%mc)

  if (allocated(this%E)) deallocate(this%E)

  if (allocated(this%H_tail_spline)) deallocate(this%H_tail_spline)
  if (allocated(this%Vrep_tail_spline)) deallocate(this%Vrep_tail_spline)

  this%n_types = 0
  this%label = ''
end subroutine TBModel_Bowler_Finalise

subroutine TBM_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  character(len=1024) :: value
  integer status
  integer ti, tj, Zi, Zj

  if (name == 'Bowler_params') then ! new Bowler stanza

    if (parse_in_tbm) &
      call system_abort("IPModel_startElement_handler entered Bowler_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

    if (parse_matched_label) return ! we already found an exact match for this label

    call QUIP_FoX_get_value(attributes, 'label', value, status)
    if (status /= 0) value = ''

    if (len(trim(parse_tbm%label)) > 0) then ! we were passed in a label
      if (value == parse_tbm%label) then ! exact match
        parse_matched_label = .true.
        parse_in_tbm = .true.
      else ! no match
        parse_in_tbm = .false.
      endif
    else ! no label passed in
      parse_in_tbm = .true.
    endif

    if (parse_in_tbm) then
      if (parse_tbm%n_types /= 0) then
        call finalise(parse_tbm)
      endif

      parse_cur_type = 1

      call QUIP_FoX_get_value(attributes, "n_types", value, status);
      if (status == 0) read (value, *) parse_tbm%n_types
      call QUIP_FoX_get_value(attributes, "tailx0", value, status);
      if (status == 0) read (value, *) parse_tbm%tailx0
      call QUIP_FoX_get_value(attributes, "cutoff", value, status);
      if (status == 0) read (value, *) parse_tbm%cutoff

      allocate(parse_tbm%atomic_num(parse_tbm%n_types))
      parse_tbm%atomic_num = 0
      allocate(parse_tbm%n_orbs(parse_tbm%n_types))
      allocate(parse_tbm%n_elecs(parse_tbm%n_types))
      allocate(parse_tbm%n_orb_sets(parse_tbm%n_types))
      allocate(parse_tbm%orb_set_type(max_n_orb_sets,parse_tbm%n_types))
      allocate(parse_tbm%E(max_n_orb_sets,parse_tbm%n_types))

      allocate(parse_tbm%H_coeff(4, parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%Vrep(parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%rc(parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%r0(parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%n(parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%nc(parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%dc(parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%m(parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%mc(parse_tbm%n_types, parse_tbm%n_types))

      parse_tbm%H_coeff = 0.0_dp
      parse_tbm%Vrep    = 0.0_dp
      parse_tbm%rc  = 1.0_dp
      parse_tbm%r0  = 0.0_dp
      parse_tbm%n   = 0.0_dp
      parse_tbm%nc  = 0.0_dp
      parse_tbm%dc  = 1.0_dp
      parse_tbm%m   = 0.0_dp
      parse_tbm%mc  = 0.0_dp
    endif

  elseif (parse_in_tbm .and. name == 'per_type_data') then
    if (parse_cur_type > parse_tbm%n_types) &
      call system_abort('Too many types defined in Bowler_params')

    call QUIP_FoX_get_value(attributes, "Z", value, status);
    if (status == 0) read (value, *) parse_tbm%atomic_num(parse_cur_type)

    call QUIP_FoX_get_value(attributes, "n_orbs", value, status);
    if (status == 0) read (value, *) parse_tbm%n_orbs(parse_cur_type)
    if (parse_tbm%n_orbs(parse_cur_type) == 1) then
      parse_tbm%n_orb_sets(parse_cur_type) = 1
      parse_tbm%orb_set_type(1,parse_cur_type) = ORB_S
      call QUIP_FoX_get_value(attributes, "E_s", value, status);
      if (status == 0) read (value, *) parse_tbm%E(ORB_S,parse_cur_type)
    else if (parse_tbm%n_orbs(parse_cur_type) == 4) then
      parse_tbm%n_orb_sets(parse_cur_type) = 2
      parse_tbm%orb_set_type(1,parse_cur_type) = ORB_S
      parse_tbm%orb_set_type(2,parse_cur_type) = ORB_P
      call QUIP_FoX_get_value(attributes, "E_s", value, status);
      if (status == 0) read (value, *) parse_tbm%E(ORB_S,parse_cur_type)
      call QUIP_FoX_get_value(attributes, "E_p", value, status);
      if (status == 0) read (value, *) parse_tbm%E(ORB_P,parse_cur_type)
    else
      call system_abort("TBModel_Bowler_read_params_xml can only do Bowler with s or sp")
    endif

    call QUIP_FoX_get_value(attributes, "n_elecs", value, status);
    if (status == 0) read (value, *) parse_tbm%n_elecs(parse_cur_type)

    if (allocated(parse_tbm%type_of_atomic_num)) deallocate(parse_tbm%type_of_atomic_num)
    allocate(parse_tbm%type_of_atomic_num(maxval(parse_tbm%atomic_num(:))))
    parse_tbm%type_of_atomic_num(:) = 0
    do ti=1, parse_tbm%n_types
      if (parse_Tbm%atomic_num(ti) > 0) &
	parse_tbm%type_of_atomic_num(parse_tbm%atomic_num(ti)) = ti
    end do

    parse_cur_type = parse_cur_type + 1

  elseif (parse_in_tbm .and. name == 'defaults') then

    call QUIP_FoX_get_value(attributes, "fermi_e", value, status)
    if (status == 0) then
      parse_tbm%has_default_fermi_e = .true.
      read (value, *) parse_tbm%default_fermi_e
    endif
    call QUIP_FoX_get_value(attributes, "fermi_T", value, status)
    if (status == 0) then
      parse_tbm%has_default_fermi_T = .true.
      read (value, *) parse_tbm%default_fermi_T
    endif
    call QUIP_FoX_get_value(attributes, "band_width", value, status)
    if (status == 0) then
      parse_tbm%has_default_band_width = .true.
      read (value, *) parse_tbm%default_band_width
    endif
    call QUIP_FoX_get_value(attributes, "k_density", value, status)
    if (status == 0) then
      parse_tbm%has_default_k_density = .true.
      read (value, *) parse_tbm%default_k_density
    endif

  elseif (parse_in_tbm .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "Z1", value, status);
    if (status == 0) read (value, *) Zi
    call QUIP_FoX_get_value(attributes, "Z2", value, status);
    if (status == 0) read (value, *) Zj

    ti = get_type(parse_tbm%type_of_atomic_num,Zi)
    tj = get_type(parse_tbm%type_of_atomic_num,Zj)

    call QUIP_FoX_get_value(attributes, "rc", value, status);
    if (status == 0) read (value, *) parse_tbm%rc(ti,tj)
    call QUIP_FoX_get_value(attributes, "r0", value, status);
    if (status == 0) read (value, *) parse_tbm%r0(ti,tj)
    call QUIP_FoX_get_value(attributes, "n", value, status);
    if (status == 0) read (value, *) parse_tbm%n(ti,tj)
    call QUIP_FoX_get_value(attributes, "nc", value, status);
    if (status == 0) read (value, *) parse_tbm%nc(ti,tj)
    call QUIP_FoX_get_value(attributes, "dc", value, status);
    if (status == 0) read (value, *) parse_tbm%dc(ti,tj)
    call QUIP_FoX_get_value(attributes, "m", value, status);
    if (status == 0) read (value, *) parse_tbm%m(ti,tj)
    call QUIP_FoX_get_value(attributes, "mc", value, status);
    if (status == 0) read (value, *) parse_tbm%mc(ti,tj)

    call QUIP_FoX_get_value(attributes, "H_sss", value, status);
    if (status == 0) read (value, *) parse_tbm%H_coeff(SK_SSS,ti,tj)
    call QUIP_FoX_get_value(attributes, "H_sps", value, status);
    if (status == 0) read (value, *) parse_tbm%H_coeff(SK_SPS,ti,tj)
    call QUIP_FoX_get_value(attributes, "H_pps", value, status);
    if (status == 0) read (value, *) parse_tbm%H_coeff(SK_PPS,ti,tj)
    call QUIP_FoX_get_value(attributes, "H_ppp", value, status);
    if (status == 0) read (value, *) parse_tbm%H_coeff(SK_PPP,ti,tj)
    call QUIP_FoX_get_value(attributes, "Vrep", value, status);
    if (status == 0) read (value, *) parse_tbm%Vrep(ti,tj)

    if (ti /= tj) then
      parse_tbm%rc(tj,ti) = parse_tbm%rc(ti,tj)
      parse_tbm%r0(tj,ti) = parse_tbm%r0(ti,tj)
      parse_tbm%n(tj,ti) = parse_tbm%n(ti,tj)
      parse_tbm%nc(tj,ti) = parse_tbm%nc(ti,tj)
      parse_tbm%dc(tj,ti) = parse_tbm%dc(ti,tj)
      parse_tbm%m(tj,ti) = parse_tbm%m(ti,tj)
      parse_tbm%mc(tj,ti) = parse_tbm%mc(ti,tj)
      parse_tbm%H_coeff(SK_SSS:SK_PPP,tj,ti) = parse_tbm%H_coeff(SK_SSS:SK_PPP,ti,tj)
      parse_tbm%Vrep(tj,ti) = parse_tbm%Vrep(ti,tj)

      call QUIP_FoX_get_value(attributes, "H_pss", value, status);
      if (status == 0) read (value, *) parse_tbm%H_coeff(SK_SPS,tj,ti)
    end if
  endif

end subroutine TBM_startElement_handler

subroutine TBM_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (name == 'Bowler_params') then
    parse_in_tbm = .false.
  endif

end subroutine TBM_endElement_handler

subroutine TBModel_Bowler_read_params_xml(this, param_str)
  type(TBModel_Bowler), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type(xml_t) :: fxml

  if (len(trim(param_str)) <= 0) return

  parse_in_tbm = .false.
  parse_matched_label = .false.
  parse_tbm => this

  call open_xml_string(fxml, param_str)

  call parse(fxml, &
    startElement_handler = TBM_startElement_handler, &
    endElement_handler = TBM_endElement_handler)

  call close_xml_t(fxml)

  if (this%n_types <= 0) call system_abort("TBModel_Bowler_read_params_xml couldn't find any types defined")

  call TBModel_Bowler_fix_tails(this)

end subroutine TBModel_Bowler_read_params_xml
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  type(inoutput), pointer :: in
!
!  integer ti, tj, Zi, Zj, tt
!
!  integer status
!  type(dictionary_t) :: attributes
!  type(xml_t) :: fxml
!  character(len=1024) :: value
!  integer :: ndata
!  character(len=10240) :: pcdata
!
!  if (present(from_io)) then
!    in => from_io
!  else
!    allocate(in)
!    call Initialise(in, trim(bowler_default_file), INPUT)
!  endif
!
!  call open_xmlfile(in%unit, fxml, status)
!! call enable_debug(sax=.true.)
!  if (status /= 0) call system_abort('TBModel_Bowler_read_params_xml Cannot open xml file')
!
!  call rewind_xmlfile(fxml)
!
!  call get_node(fxml, path='//Bowler_params', attributes=attributes, status=status)
!  if (status /= 0) call system_abort('TBModel_Bowler_read_params_xml Cannot find /Bowler_params/header')
!  call QUIP_FoX_get_value(attributes, "name", value, status);
!  if (status == 0) read (value, *) this%name
!
!  call QUIP_FoX_get_value(attributes, "n_types", value, status);
!  if (status == 0) read (value, *) this%n_types
!  call QUIP_FoX_get_value(attributes, "tailx0", value, status);
!  if (status == 0) read (value, *) this%tailx0
!  call QUIP_FoX_get_value(attributes, "cutoff", value, status);
!  if (status == 0) read (value, *) this%cutoff
!
!  allocate(this%atomic_num(this%n_types))
!  allocate(this%n_orbs(this%n_types))
!  allocate(this%n_elecs(this%n_types))
!  allocate(this%n_orb_sets(this%n_types))
!  allocate(this%orb_set_type(max_n_orb_sets,this%n_types))
!  allocate(this%E(max_n_orb_sets,this%n_types))
!
!  call rewind_xmlfile(fxml)
!
!  !! read bowler params
!  do ti=1, this%n_types
!    call get_node(fxml, path='//Bowler_params/per_type_data', attributes=attributes, status=status)
!    if (status /= 0) call system_abort("Couldn't find per_type_data for " // str(ti))
!
!    call QUIP_FoX_get_value(attributes, "Z", value, status);
!    if (status == 0) read (value, *) this%atomic_num(ti)
!
!    call QUIP_FoX_get_value(attributes, "n_orbs", value, status);
!    if (status == 0) read (value, *) this%n_orbs(ti)
!    if (this%n_orbs(ti) == 1) then
!      this%n_orb_sets(ti) = 1
!      this%orb_set_type(1,ti) = ORB_S
!      call QUIP_FoX_get_value(attributes, "E_s", value, status);
!      if (status == 0) read (value, *) this%E(ORB_S,ti)
!    else if (this%n_orbs(ti) == 4) then
!      this%n_orb_sets(ti) = 2
!      this%orb_set_type(1,ti) = ORB_S
!      this%orb_set_type(2,ti) = ORB_P
!      call QUIP_FoX_get_value(attributes, "E_s", value, status);
!      if (status == 0) read (value, *) this%E(ORB_S,ti)
!      call QUIP_FoX_get_value(attributes, "E_p", value, status);
!      if (status == 0) read (value, *) this%E(ORB_P,ti)
!    else
!      call system_abort("TBModel_Bowler_read_params_xml can only do Bowler with s or sp")
!    endif
!
!    call QUIP_FoX_get_value(attributes, "n_elecs", value, status);
!    if (status == 0) read (value, *) this%n_elecs(ti)
!
!  end do
!
!  allocate(this%type_of_atomic_num(maxval(this%atomic_num(:))))
!  this%type_of_atomic_num(:) = 0
!  do ti=1, this%n_types
!    this%type_of_atomic_num(this%atomic_num(ti)) = ti
!  end do
!
!  allocate(this%H_coeff(4, this%n_types, this%n_types))
!  allocate(this%Vrep(this%n_types, this%n_types))
!  allocate(this%rc(this%n_types, this%n_types))
!  allocate(this%r0(this%n_types, this%n_types))
!  allocate(this%n(this%n_types, this%n_types))
!  allocate(this%nc(this%n_types, this%n_types))
!  allocate(this%dc(this%n_types, this%n_types))
!  allocate(this%m(this%n_types, this%n_types))
!  allocate(this%mc(this%n_types, this%n_types))
!  call rewind_xmlfile(fxml)
!
!  this%H_coeff = 0.0_dp
!  this%Vrep    = 0.0_dp
!  this%rc  = 1.0_dp
!  this%r0  = 0.0_dp
!  this%n   = 0.0_dp
!  this%nc  = 0.0_dp
!  this%dc  = 1.0_dp
!  this%m   = 0.0_dp
!  this%mc  = 0.0_dp
!
!  do tt=1, (this%n_types*(this%n_types+1))/2
!    call get_node(fxml, path='//Bowler_params/per_pair_data', attributes=attributes, status=status)
!    if (status /= 0) call system_abort("Couldn't find per_pair_data for " // trim(str(ti)) // ' ' // trim(str(tj)))
!
!    call QUIP_FoX_get_value(attributes, "Z1", value, status);
!    if (status == 0) read (value, *) Zi
!    call QUIP_FoX_get_value(attributes, "Z2", value, status);
!    if (status == 0) read (value, *) Zj
!
!    ti = get_type(this%type_of_atomic_num,Zi)
!    tj = get_type(this%type_of_atomic_num,Zj)
!
!    call QUIP_FoX_get_value(attributes, "rc", value, status);
!    if (status == 0) read (value, *) this%rc(ti,tj)
!    call QUIP_FoX_get_value(attributes, "r0", value, status);
!    if (status == 0) read (value, *) this%r0(ti,tj)
!    call QUIP_FoX_get_value(attributes, "n", value, status);
!    if (status == 0) read (value, *) this%n(ti,tj)
!    call QUIP_FoX_get_value(attributes, "nc", value, status);
!    if (status == 0) read (value, *) this%nc(ti,tj)
!    call QUIP_FoX_get_value(attributes, "dc", value, status);
!    if (status == 0) read (value, *) this%dc(ti,tj)
!    call QUIP_FoX_get_value(attributes, "m", value, status);
!    if (status == 0) read (value, *) this%m(ti,tj)
!    call QUIP_FoX_get_value(attributes, "mc", value, status);
!    if (status == 0) read (value, *) this%mc(ti,tj)
!
!    call QUIP_FoX_get_value(attributes, "H_sss", value, status);
!    if (status == 0) read (value, *) this%H_coeff(SK_SSS,ti,tj)
!    call QUIP_FoX_get_value(attributes, "H_sps", value, status);
!    if (status == 0) read (value, *) this%H_coeff(SK_SPS,ti,tj)
!    call QUIP_FoX_get_value(attributes, "H_pps", value, status);
!    if (status == 0) read (value, *) this%H_coeff(SK_PPS,ti,tj)
!    call QUIP_FoX_get_value(attributes, "H_ppp", value, status);
!    if (status == 0) read (value, *) this%H_coeff(SK_PPP,ti,tj)
!    call QUIP_FoX_get_value(attributes, "Vrep", value, status);
!    if (status == 0) read (value, *) this%Vrep(ti,tj)
!
!    if (ti /= tj) then
!      this%rc(tj,ti) = this%rc(ti,tj)
!      this%r0(tj,ti) = this%r0(ti,tj)
!      this%n(tj,ti) = this%n(ti,tj)
!      this%nc(tj,ti) = this%nc(ti,tj)
!      this%dc(tj,ti) = this%dc(ti,tj)
!      this%m(tj,ti) = this%m(ti,tj)
!      this%mc(tj,ti) = this%mc(ti,tj)
!      this%H_coeff(SK_SSS:SK_PPP,tj,ti) = this%H_coeff(SK_SSS:SK_PPP,ti,tj)
!      this%Vrep(tj,ti) = this%Vrep(ti,tj)
!
!      call QUIP_FoX_get_value(attributes, "H_pss", value, status);
!      if (status == 0) read (value, *) this%H_coeff(SK_SPS,tj,ti)
!    endif
!  end do
!
!  call TBModel_Bowler_fix_tails(this)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine TBModel_Bowler_fix_tails(this)
  type(TBModel_Bowler), intent(inout) :: this

  integer ti, tj
  real(dp) :: x(2), y(2), yd1, yd2

  x(1) = this%tailx0
  x(2) = this%cutoff
  y(2) = 0.0_dp
  yd2 = 0.0_dp

  allocate(this%H_tail_spline(this%n_types, this%n_types))
  allocate(this%Vrep_tail_spline(this%n_types, this%n_types))

  do ti=1, this%n_types
  do tj=1, this%n_types
    y(1) = TBModel_Bowler_dist_scaling(x(1), this%r0(ti,tj), this%n(ti,tj), this%rc(ti,tj), this%nc(ti,tj))
    yd1 = TBModel_Bowler_dist_scaling_deriv(x(1), this%r0(ti,tj), this%n(ti,tj), this%rc(ti,tj), this%nc(ti,tj))
    call initialise(this%H_tail_spline(ti,tj), x, y, yd1, yd2)

    y(1) = TBModel_Bowler_dist_scaling(x(1), this%r0(ti,tj), this%m(ti,tj), this%dc(ti,tj), this%mc(ti,tj))
    yd1 = TBModel_Bowler_dist_scaling_deriv(x(1), this%r0(ti,tj), this%m(ti,tj), this%dc(ti,tj), this%mc(ti,tj))
    call initialise(this%Vrep_tail_spline(ti,tj), x, y, yd1, yd2)
  end do
  end do

end subroutine TBModel_Bowler_fix_tails

subroutine TBModel_Bowler_Print(this,file)
  type(TBModel_Bowler),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  integer:: ti, tj

  if(this%n_types == 0) call System_Abort('TBModel_Bowler_Print: TBModel_Bowler structure not initialised')

  call Print('TBModel_Bowler Structure:', file=file)

  !! print bowler params
  call Print ("TBModel_Bowler tailx0 " // this%tailx0 // " cutoff " // this%cutoff, file=file)
  do ti=1, this%n_types
    call Print ("TBModel_Bowler type " // ti // " Z " // this%atomic_num(ti) //  " n_orbs " // this%n_orbs(ti) // &
      " n_elecs " // this%n_elecs(ti), file=file)
    if (this%n_orb_sets(ti) == 1) then
      call Print ("TBModel_Bowler E " // this%E(1,ti), file=file)
    else
      call Print ("TBModel_Bowler E " // this%E(1:2,ti), file=file)
    endif
  end do

  call verbosity_push_decrement()
  do ti=1, this%n_types
  do tj=1, this%n_types
    call Print ("TBModel_Bowler interaction " // &
      ti // " " //  tj // " Z " //  this%atomic_num(ti) // " " // this%atomic_num(tj), file=file)

    call Print ("TBModel_Bowler r0 " // this%r0(ti,tj), file=file)

    call Print ("TBModel_Bowler SK " // this%H_coeff(:,ti,tj), file=file)
    call Print ("TBModel_Bowler H scaling rc n nc " // &
       this%rc(ti,tj) // " " //  this%n(ti,tj) //  " " // this%nc(ti,tj), file=file)
    call Print (this%H_tail_spline(ti,tj), file=file)

    call Print ("TBModel_Bowler Vrep " // this%Vrep(ti,tj), file=file)
    call Print ("TBModel_Bowler Vrep scaling dc m mc " // &
       this%dc(ti,tj) // " " //  this%m(ti,tj) //  this%mc(ti,tj), file=file)
    call Print (this%Vrep_tail_spline(ti,tj), file=file)

  end do
  end do
  call verbosity_pop()

end subroutine TBModel_Bowler_Print

function TBModel_Bowler_n_orbs_of_Z(this, Z)
  type(TBModel_Bowler), intent(in) :: this
  integer, intent(in) :: Z
  integer TBModel_Bowler_n_orbs_of_Z

  TBModel_Bowler_n_orbs_of_Z = this%n_orbs(get_type(this%type_of_atomic_num,Z))
end function TBModel_Bowler_n_orbs_of_Z

function TBModel_Bowler_n_orb_sets_of_Z(this, Z)
  type(TBModel_Bowler), intent(in) :: this
  integer, intent(in) :: Z
  integer TBModel_Bowler_n_orb_sets_of_Z

  TBModel_Bowler_n_orb_sets_of_Z = this%n_orb_sets(get_type(this%type_of_atomic_num,Z))
end function TBModel_Bowler_n_orb_sets_of_Z

function TBModel_Bowler_n_orbs_of_orb_set_of_Z(this, Z, i_set)
  type(TBModel_Bowler), intent(in) :: this
  integer, intent(in) :: Z, i_set
  integer TBModel_Bowler_n_orbs_of_orb_set_of_Z

  TBModel_Bowler_n_orbs_of_orb_set_of_Z = N_ORBS_OF_SET(this%orb_set_type(i_set,get_type(this%type_of_atomic_num,Z)))

end function TBModel_Bowler_n_orbs_of_orb_set_of_Z

function TBModel_Bowler_orb_type_of_orb_set_of_Z(this, Z, i_set)
  type(TBModel_Bowler), intent(in) :: this
  integer, intent(in) :: Z, i_set
  integer TBModel_Bowler_orb_type_of_orb_set_of_Z

  TBModel_Bowler_orb_type_of_orb_set_of_Z = this%orb_set_type(i_set,get_type(this%type_of_atomic_num,Z))

end function TBModel_Bowler_orb_type_of_orb_set_of_Z

function TBModel_Bowler_n_elecs_of_Z(this, Z)
  type(TBModel_Bowler), intent(in) :: this
  integer, intent(in) :: Z
  real(dp) TBModel_Bowler_n_elecs_of_Z

  TBModel_Bowler_n_elecs_of_Z = this%n_elecs(get_type(this%type_of_atomic_num,Z))
end function TBModel_Bowler_n_elecs_of_Z

subroutine TBModel_Bowler_get_HS_blocks(this, at, at_i, at_j, dv_hat, dv_mag, b_H, b_S, i_mag)
  type(TBModel_Bowler), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_i, at_j
  real(dp), intent(in) :: dv_hat(3), dv_mag
  real(dp), intent(out) :: b_H(:,:)
  real(dp), intent(out) :: b_S(:,:)
  integer, intent(in), optional :: i_mag

  integer ti, tj, is, js, i_set, j_set
  integer i, j
  real(dp) SK_frad_H(N_SK)
  real(dp) dv_hat_sq(3)

  ti = get_type(this%type_of_atomic_num,at%Z(at_i))
  tj = get_type(this%type_of_atomic_num,at%Z(at_j))

  b_S = 0.0_dp

  if (dv_mag .feq. 0.0_dp) then
    b_H = 0.0_dp

    i = 1
    do i_set = 1, this%n_orb_sets(ti)
    do is=1, N_ORBS_OF_SET(this%orb_set_type(i_set,ti))
      b_H(i,i) = onsite_function(this, ti, this%orb_set_type(i_set,ti))
      b_S(i,i) = 1.0_dp
    i = i + 1
    end do
    end do

  else

    dv_hat_sq = dv_hat**2

    i = 1
    do i_set=1, this%n_orb_sets(ti)
    do is=1, N_ORBS_OF_SET(this%orb_set_type(i_set,ti))
      j = 1
      do j_set=1, this%n_orb_sets(tj)
      call radial_functions(this, ti, tj, dv_mag, this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	is, js, SK_frad_H)
      do js=1, N_ORBS_OF_SET(this%orb_set_type(j_set,tj))

	b_H(i,j) = angular_function(dv_hat, dv_hat_sq, this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	  is, js, SK_frad_H)
	j = j + 1
      end do
      end do
      i = i + 1
    end do
    end do
  endif

end subroutine TBModel_Bowler_get_HS_blocks

subroutine TBModel_Bowler_get_dHS_masks(this, at, at_ind, d_mask, od_mask)
  type(TBModel_Bowler), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_ind
  logical, intent(out), optional :: d_mask(:), od_mask(:)

  if (at_ind < 0) then
    if (present(d_mask)) d_mask = .true.
    if (present(od_mask)) od_mask = .true.
  else
    if (present(d_mask)) d_mask = .false.
    if (present(od_mask)) then
      od_mask = .false.
      od_mask(at_ind) = .true.
    endif
  endif
end subroutine TBModel_Bowler_get_dHS_masks

function TBModel_Bowler_get_dHS_blocks(this, at, at_i, at_j, dv_hat, dv_mag, at_ind, b_dH, b_dS, i_mag)
  type(TBModel_Bowler), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_i, at_j
  real(dp), intent(in) :: dv_hat(3), dv_mag
  integer, intent(in) :: at_ind
  real(dp), intent(out) :: b_dH(:,:,:)
  real(dp), intent(out) :: b_dS(:,:,:)
  integer, intent(in), optional :: i_mag
  logical :: TBModel_Bowler_get_dHS_blocks

  integer ti, tj, is, js, i_set, j_set
  integer i, j
  real(dp) SK_frad_H(N_SK)
  real(dp) SK_dfrad_H(N_SK)
  real(dp) dv_hat_sq(3)
  real(dp) virial_outerprod_fac

  if ((at_ind > 0 .and. at_i /= at_ind .and. at_j /= at_ind) .or. &
      (at_ind > 0 .and. at_i == at_j) .or. &
      (dv_mag > this%cutoff) .or. &
      (dv_mag .feq. 0.0_dp)) then
    TBModel_Bowler_get_dHS_blocks = .false.
    return
  endif

  ti = get_type(this%type_of_atomic_num,at%Z(at_i))
  tj = get_type(this%type_of_atomic_num,at%Z(at_j))

  b_dS = 0.0_dp

  dv_hat_sq = dv_hat**2


  i = 1
  do i_set=1, this%n_orb_sets(ti)
  do is=1, N_ORBS_OF_SET(this%orb_set_type(i_set,ti))
    j = 1
    do j_set=1, this%n_orb_sets(tj)
    call radial_functions(this, ti, tj, dv_mag, this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
      is, js, SK_frad_H)
    call dradial_functions(this, ti, tj, dv_mag, this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
      is, js, SK_dfrad_H)
    do js=1, N_ORBS_OF_SET(this%orb_set_type(j_set,tj))

      if (at_ind > 0) then
	virial_outerprod_fac = 1.0_dp
      else
	virial_outerprod_fac = -dv_mag*dv_hat(-at_ind)
      end if
      b_dH(i,j,:) = dangular_function(dv_mag, dv_hat, dv_hat_sq, &
	this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	is, js, SK_frad_H, SK_dfrad_H)*virial_outerprod_fac
      j = j + 1
    end do
    end do
    i = i + 1
  end do
  end do

  if (at_ind == at_j) b_dH = -b_dH

  TBModel_Bowler_get_dHS_blocks = .true.
  return

end function TBModel_Bowler_get_dHS_blocks

subroutine radial_functions(this, ti, tj, dv_mag, orb_set_type_i, orb_set_type_j, is, js, f_H)
  type(TBModel_Bowler), intent(in) :: this
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: dv_mag
  integer, intent(in) :: orb_set_type_i, orb_set_type_j
  integer, intent(in) :: is, js
  real(dp), intent(out) :: f_H(N_SK)

  if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_S) then
    f_H(SK_SSS) = calc_H_coeff(this, SK_SSS, dv_mag, ti, tj)
  else if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_P) then
    f_H(SK_SPS) = calc_H_coeff(this, SK_SPS, dv_mag, ti, tj)
  else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_S) then
    f_H(SK_SPS) = -calc_H_coeff( this, SK_SPS, dv_mag, tj, ti)
  else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_P) then
    f_H(SK_PPS) = calc_H_coeff(this, SK_PPS, dv_mag, ti, tj)
    f_H(SK_PPP) = calc_H_coeff(this, SK_PPP, dv_mag, ti, tj)
  else
    call system_abort("TBModel_Bowler radial_functions got invalide orb_set_type "//orb_set_type_i//" or "//orb_set_type_j)
  endif
end subroutine radial_functions

subroutine dradial_functions(this, ti, tj, dv_mag, orb_set_type_i, orb_set_type_j, is, js, f_dH)
  type(TBModel_Bowler), intent(in) :: this
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: dv_mag
  integer, intent(in) :: orb_set_type_i, orb_set_type_j
  integer, intent(in) :: is, js
  real(dp), intent(out) :: f_dH(N_SK)

  if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_S) then
    f_dH(SK_SSS) = calc_H_coeff_deriv(this, SK_SSS, dv_mag, ti, tj)
  else if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_P) then
    f_dH(SK_SPS) = calc_H_coeff_deriv(this, SK_SPS, dv_mag, ti, tj)
  else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_S) then
    f_dH(SK_SPS) = -calc_H_coeff_deriv( this, SK_SPS, dv_mag, tj, ti)
  else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_P) then
    f_dH(SK_PPS) = calc_H_coeff_deriv(this, SK_PPS, dv_mag, ti, tj)
    f_dH(SK_PPP) = calc_H_coeff_deriv(this, SK_PPP, dv_mag, ti, tj)
  else
    call system_abort("TBModel_Bowler dradial_functions got invalide orb_set_type "//orb_set_type_i//" or "//orb_set_type_j)
  endif
end subroutine dradial_functions

function TBModel_Bowler_calc_H_coeff(this, sk_ind, dist, ti, tj)
  type(TBModel_Bowler), intent(in) :: this
  integer, intent(in) :: sk_ind
  real(dp), intent(in) :: dist
  integer, intent(in) :: ti, tj
  real(dp) :: TBModel_Bowler_calc_H_coeff

  TBModel_Bowler_calc_H_coeff = this%H_coeff(sk_ind, ti, tj)* &
    TBModel_Bowler_H_dist_func(this, dist, ti, tj)

end function TBModel_Bowler_calc_H_coeff

function TBModel_Bowler_calc_H_coeff_deriv(this, sk_ind, dist, ti, tj)
  type(TBModel_Bowler), intent(in) :: this
  integer, intent(in) :: sk_ind
  real(dp), intent(in) :: dist
  integer, intent(in) :: ti, tj
  real(dp) :: TBModel_Bowler_calc_H_coeff_deriv

  TBModel_Bowler_calc_H_coeff_deriv = this%H_coeff(sk_ind, ti, tj)* &
    TBModel_Bowler_H_dist_func_deriv(this, dist, ti, tj)

end function TBModel_Bowler_calc_H_coeff_deriv

function TBModel_Bowler_H_dist_func(this, dist, ti, tj)
  type(TBModel_Bowler), intent(in) :: this
  real(dp), intent(in) :: dist
  integer, intent(in) :: ti, tj
  real(dp) :: TBModel_Bowler_H_dist_func

  if (dist <= this%cutoff) then
    if (dist <= this%tailx0) then
      TBModel_Bowler_H_dist_func = TBModel_Bowler_dist_scaling(dist, &
	this%r0(ti,tj), this%n(ti,tj), this%rc(ti, tj), this%nc(ti, tj))
    else
      TBModel_Bowler_H_dist_func = spline_value(this%H_tail_spline(ti,tj), dist)
    endif
  else
    TBModel_Bowler_H_dist_func = 0.0_dp
  endif
end function TBModel_Bowler_H_dist_func

function TBModel_Bowler_H_dist_func_deriv(this, dist, ti, tj)
  type(TBModel_Bowler), intent(in) :: this
  real(dp), intent(in) :: dist
  integer, intent(in) :: ti, tj
  real(dp) :: TBModel_Bowler_H_dist_func_deriv

  if (dist <= this%cutoff) then
    if (dist <= this%tailx0) then
      TBModel_Bowler_H_dist_func_deriv = TBModel_Bowler_dist_scaling_deriv(dist, &
	this%r0(ti,tj), this%n(ti,tj), this%rc(ti, tj), this%nc(ti, tj))
    else
      TBModel_Bowler_H_dist_func_deriv = spline_deriv(this%H_tail_spline(ti,tj), dist)
    endif
  else
    TBModel_Bowler_H_dist_func_deriv = 0.0_dp
  endif
end function TBModel_Bowler_H_dist_func_deriv

function TBModel_Bowler_Vrep_dist_func(this, dist, ti, tj)
  type(TBModel_Bowler), intent(in) :: this
  real(dp), intent(in) :: dist
  integer, intent(in) :: ti, tj
  real(dp) :: TBModel_Bowler_Vrep_dist_func

  if (dist <= this%cutoff) then
    if (dist <= this%tailx0) then
      TBModel_Bowler_Vrep_dist_func = TBModel_Bowler_dist_scaling(dist, &
	this%r0(ti,tj), this%m(ti,tj), this%dc(ti, tj), this%mc(ti, tj))
    else
      TBModel_Bowler_Vrep_dist_func = spline_value(this%Vrep_tail_spline(ti,tj), dist)
    endif
  else
    TBModel_Bowler_Vrep_dist_func = 0.0_dp
  endif
end function TBModel_Bowler_Vrep_dist_func

function TBModel_Bowler_Vrep_dist_func_deriv(this, dist, ti, tj)
  type(TBModel_Bowler), intent(in) :: this
  real(dp), intent(in) :: dist
  integer, intent(in) :: ti, tj
  real(dp) :: TBModel_Bowler_Vrep_dist_func_deriv

  if (dist <= this%cutoff) then
    if (dist <= this%tailx0) then
      TBModel_Bowler_Vrep_dist_func_deriv = TBModel_Bowler_dist_scaling_deriv(dist, &
	this%r0(ti,tj), this%m(ti,tj), this%dc(ti, tj), this%mc(ti, tj))
    else
      TBModel_Bowler_Vrep_dist_func_deriv = spline_deriv(this%Vrep_tail_spline(ti,tj), dist)
    endif
  else
    TBModel_Bowler_Vrep_dist_func_deriv = 0.0_dp
  endif
end function TBModel_Bowler_Vrep_dist_func_deriv

function TBModel_Bowler_dist_scaling(r, r0, n, rc, nc)
  real(dp), intent(in) :: r, r0, n, rc, nc
  real(dp) :: TBModel_Bowler_dist_scaling

  TBModel_Bowler_dist_scaling = (r0/r)**n * exp(n* ((r0/rc)**nc - (r/rc)**nc))
end function TBModel_Bowler_dist_scaling

function TBModel_Bowler_dist_scaling_deriv(r, r0, n, rc, nc)
  real(dp), intent(in) :: r, r0, n, rc, nc
  real(dp) :: TBModel_Bowler_dist_scaling_deriv

  real(dp) tmp1, tmp2

  tmp1 = (r0/r)**n
  tmp2 = exp(n*((r0/rc)**nc - (r/rc)**nc))

  TBModel_Bowler_dist_scaling_deriv = n*tmp1/(-r)*tmp2 + tmp1*tmp2*(-n)*nc*((r/rc)**nc)/r 
end function TBModel_Bowler_dist_scaling_deriv

function onsite_function(this, ti, orb_set_type)
  type(TBModel_Bowler), intent(in) :: this
  integer, intent(in) :: ti, orb_set_type
  real(dp) :: onsite_function

  onsite_function = this%E(orb_set_type, ti)

end function onsite_function

function TBModel_Bowler_get_local_rep_E(this, at, i)
  type(TBModel_Bowler), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: TBModel_Bowler_get_local_rep_E

  real(dp) E, dist
  integer ji, j, ti, tj

  TBModel_Bowler_get_local_rep_E = 0.0_dp

  ti = get_type(this%type_of_atomic_num, at%Z(i))
  do ji=1, n_neighbours(at, i)
    j = neighbour(at, i, ji, dist)
    if (dist .feq. 0.0_dp) cycle
    tj = get_type(this%type_of_atomic_num, at%Z(j))

    E = this%Vrep(ti,tj) * TBModel_Bowler_Vrep_dist_func(this, dist, ti, tj)
    TBModel_Bowler_get_local_rep_E = TBModel_Bowler_get_local_rep_E + E/2.0_dp
  end do

end function TBModel_Bowler_get_local_rep_E

function TBModel_Bowler_get_local_rep_E_force(this, at, i) result(force)
  type(TBModel_Bowler), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: force(3,at%N)

  real(dp) dE_dr, dist, dv_hat(3)
  integer ji, j, ti, tj

  force = 0.0_dp

  ti = get_type(this%type_of_atomic_num, at%Z(i))
  do ji=1, n_neighbours(at, i)
    j = neighbour(at, i, ji, dist, cosines = dv_hat)
    if (dist .feq. 0.0_dp) cycle
    tj = get_type(this%type_of_atomic_num, at%Z(j))

    dE_dr = this%Vrep(ti,tj) * TBModel_Bowler_Vrep_dist_func_deriv(this, dist, ti, tj)
    force(:,i) = force(:,i) + dE_dr*dv_hat(:)/2.0_dp
    force(:,j) = force(:,j) - dE_dr*dv_hat(:)/2.0_dp
  end do

end function TBModel_Bowler_get_local_rep_E_force

function TBModel_Bowler_get_local_rep_E_virial(this, at, i) result(virial)
  type(TBModel_Bowler), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: virial(3,3)

  real(dp) dE_dr, dist, dv_hat(3)
  integer ji, j, ti, tj

  virial = 0.0_dp

  ti = get_type(this%type_of_atomic_num, at%Z(i))
  do ji=1, n_neighbours(at, i)
    j = neighbour(at, i, ji, dist, cosines = dv_hat)
    if (dist .feq. 0.0_dp) cycle
    tj = get_type(this%type_of_atomic_num, at%Z(j))

    dE_dr = this%Vrep(ti,tj) * TBModel_Bowler_Vrep_dist_func_deriv(this, dist, ti, tj)
    virial = virial - dE_dr*(dv_hat .outer. dv_hat) * dist / 2.0_dp
  end do

end function TBModel_Bowler_get_local_rep_E_virial

end module TBModel_Bowler_module
