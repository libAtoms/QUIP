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
!X TBModel_NRL_TB module
!X
!% Calculate energies using NRL-TB tight-binding model
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module TBModel_NRL_TB_module

use system_module, only : dp, print, inoutput, optional_default, system_abort, verbosity_push_decrement, verbosity_pop
use periodictable_module
use units_module
use linearalgebra_module
use dictionary_module
use extendable_str_module
use paramreader_module
use atoms_module

use TB_Common_module
use QUIP_Common_module
use TBModel_NRL_TB_defs_module

implicit none
private

! character(30) :: nrltb_default_file = "tightbind.parms.NRL_TB.xml"

include 'TBModel_interface.h'

public :: TBModel_NRL_TB
type TBModel_NRL_TB
  integer :: n_types = 0, n_mag = 0
  logical :: is_orthogonal = .true., is_magnetic = .false., has_pair_repulsion = .false., &
	     overlap_zero_limit = .true., force_harrison_signs = .false.
  character(len=STRING_LENGTH) label

  real(dp) :: cutoff = 0.0_dp

  integer, allocatable :: type_of_atomic_num(:)
  integer, allocatable :: n_orbs(:), n_elecs(:), n_orb_sets(:), orb_set_type(:,:)
  integer, allocatable :: atomic_num(:)

  real(dp), allocatable :: atomic_mass(:)
  real(dp), allocatable :: r_cut(:,:), screen_l(:,:)
  real(dp), allocatable :: pair_rep_inner(:,:), pair_rep_outer(:,:)

  real(dp), allocatable :: lambda_sq(:,:)
  real(dp), allocatable :: abcd(:,:,:,:,:)

  real(dp), allocatable :: H_coeff(:,:,:,:,:), S_coeff(:,:,:,:,:)

  logical :: has_default_fermi_e = .false., has_default_fermi_T = .true., has_default_band_width = .false., has_default_k_density=.false.
  real(dp) :: default_fermi_e, default_fermi_T = 0.001_dp, default_band_width, default_k_density

end type TBModel_NRL_TB

integer :: parse_cur_type_i, parse_cur_type_j
type(extendable_str), private, save :: parse_cur_data
logical, private :: parse_in_tbm, parse_matched_label
type (TBModel_NRL_TB), pointer :: parse_tbm

interface Initialise
  module procedure TBModel_NRL_TB_Initialise_str
end interface Initialise

interface Finalise
  module procedure TBModel_NRL_TB_Finalise
end interface Finalise

interface Print
  module procedure TBModel_NRL_TB_Print
end interface Print

interface n_orbs_of_Z
  module procedure TBModel_NRL_TB_n_orbs_of_Z
end interface n_orbs_of_Z

interface n_orb_sets_of_Z
  module procedure TBModel_NRL_TB_n_orb_sets_of_Z
end interface n_orb_sets_of_Z

interface n_orbs_of_orb_set_of_Z
  module procedure TBModel_NRL_TB_n_orbs_of_orb_set_of_Z
end interface n_orbs_of_orb_set_of_Z

interface orb_type_of_orb_set_of_Z
  module procedure TBModel_NRL_TB_orb_type_of_orb_set_of_Z
end interface orb_type_of_orb_set_of_Z

interface n_elecs_of_Z
  module procedure TBModel_NRL_TB_n_elecs_of_Z
end interface n_elecs_of_Z

interface get_HS_blocks
  module procedure TBModel_NRL_TB_get_HS_blocks
end interface get_HS_blocks

interface get_dHS_masks
  module procedure TBModel_NRL_TB_get_dHS_masks
end interface get_dHS_masks

interface get_dHS_blocks
  module procedure TBModel_NRL_TB_get_dHS_blocks
end interface get_dHS_blocks

interface get_local_rep_E
  module procedure TBModel_NRL_TB_get_local_rep_E
end interface get_local_rep_E

interface get_local_rep_E_force
  module procedure TBModel_NRL_TB_get_local_rep_E_force
end interface get_local_rep_E_force

interface get_local_rep_E_virial
  module procedure TBModel_NRL_TB_get_local_rep_E_virial
end interface get_local_rep_E_virial

contains

subroutine TBModel_NRL_TB_Initialise_str(this, args_str, param_str)
  type(TBModel_NRL_TB), intent(inout) :: this
  character(len=*), intent(in) :: args_str
  character(len=*), intent(in) :: param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='TBModel_NRL_TB_Initialise_str args_str')) then
    call system_abort("TBModel_NRL_TB_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call TBModel_NRL_TB_read_params_xml(this, param_str)
end subroutine TBModel_NRL_TB_Initialise_str

subroutine TBModel_NRL_TB_Finalise(this)
  type(TBModel_NRL_TB), intent(inout) :: this

  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%n_orbs)) deallocate(this%n_orbs)
  if (allocated(this%n_elecs)) deallocate(this%n_elecs)
  if (allocated(this%n_orb_sets)) deallocate(this%n_orb_sets)
  if (allocated(this%orb_set_type)) deallocate(this%orb_set_type)

  if (allocated(this%atomic_mass)) deallocate(this%atomic_mass)
  if (allocated(this%r_cut)) deallocate(this%r_cut)
  if (allocated(this%screen_l)) deallocate(this%screen_l)
  if (allocated(this%pair_rep_inner)) deallocate(this%pair_rep_inner)
  if (allocated(this%pair_rep_outer)) deallocate(this%pair_rep_outer)

  if (allocated(this%lambda_sq)) deallocate(this%lambda_sq)
  if (allocated(this%abcd)) deallocate(this%abcd)
  if (allocated(this%H_coeff)) deallocate(this%H_coeff)
  if (allocated(this%S_coeff)) deallocate(this%S_coeff)

  this%n_types = 0
  this%n_mag = 0
  this%label = ''
end subroutine TBModel_NRL_TB_Finalise

subroutine TBM_characters_handler(in)
  character(len=*), intent(in) :: in

  if (parse_in_tbm) then
    call concat(parse_cur_data, in, keep_lf=.false.)
  endif
end subroutine TBM_characters_handler

subroutine TBM_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer status
  character(len=1024) :: value

  integer ii, ti, tj

  call zero(parse_cur_data)

  if (name == 'NRL_TB_params') then ! new NRL_TB stanza

    if (parse_in_tbm) &
      call system_abort("IPModel_startElement_handler entered NRL_TB_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

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
      call Initialise(parse_cur_data)
      if (parse_tbm%n_types /= 0) then
        call finalise(parse_tbm)
      endif
    endif

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

  elseif (parse_in_tbm .and. name == 'header') then

    call QUIP_FoX_get_value(attributes, "is_orthogonal", value, status);
    if (status == 0) read (value, *) parse_tbm%is_orthogonal
    call QUIP_FoX_get_value(attributes, "is_magnetic", value, status)
    if (status == 0) read (value, *) parse_tbm%is_magnetic
    if (parse_tbm%is_magnetic) then
      parse_tbm%n_mag = 2
    else
      parse_tbm%n_mag = 1
    endif
    call QUIP_FoX_get_value(attributes, "has_pair_repulsion", value, status)
    if (status == 0) read (value, *) parse_tbm%has_pair_repulsion
    call QUIP_FoX_get_value(attributes, "overlap_zero_limit", value, status)
    if (status == 0) read (value, *) parse_tbm%overlap_zero_limit
    call QUIP_FoX_get_value(attributes, "force_harrison_signs", value, status)
    if (status == 0) read (value, *) parse_tbm%force_harrison_signs

  elseif (parse_in_tbm .and. name == 'n_types') then

    call QUIP_FoX_get_value(attributes, "v", value, status);
    if (status == 0) read (value, *) parse_tbm%n_types

    allocate(parse_tbm%atomic_num(parse_tbm%n_types))
    parse_tbm%atomic_num = 0
    allocate(parse_tbm%atomic_mass(parse_tbm%n_types))
    allocate(parse_tbm%n_orbs(parse_tbm%n_types))
    allocate(parse_tbm%n_elecs(parse_tbm%n_types))
    allocate(parse_tbm%n_orb_sets(parse_tbm%n_types))
    allocate(parse_tbm%orb_set_type(max_n_orb_sets,parse_tbm%n_types))
    allocate(parse_tbm%r_cut(parse_tbm%n_types,parse_tbm%n_types))
    allocate(parse_tbm%screen_l(parse_tbm%n_types,parse_tbm%n_types))
    if (parse_tbm%has_pair_repulsion) then
      allocate(parse_tbm%pair_rep_inner(parse_tbm%n_types,parse_tbm%n_types))
      allocate(parse_tbm%pair_rep_outer(parse_tbm%n_types,parse_tbm%n_types))
    endif
    allocate(parse_tbm%lambda_sq(parse_tbm%n_types,parse_tbm%n_mag))
    allocate(parse_tbm%abcd(4,max_n_orb_sets,parse_tbm%n_types,parse_tbm%n_types,parse_tbm%n_mag))
    allocate(parse_tbm%H_coeff(4,N_SK,parse_tbm%n_types,parse_tbm%n_types,parse_tbm%n_mag))
    allocate(parse_tbm%S_coeff(4,N_SK,parse_tbm%n_types,parse_tbm%n_types,parse_tbm%n_mag))

  elseif (parse_in_tbm .and. name == 'per_type_data') then

    if (parse_tbm%n_types == 0) &
      call system_abort("NRL_TB_params got to per_type_data before finding n_types")

    call QUIP_FoX_get_value(attributes, "type", value, status);
    if (status == 0) read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_mass", value, status);
    if (status == 0) read (value, *) parse_tbm%atomic_mass(ti)
    parse_tbm%atomic_num(ti) = 0
    call QUIP_FoX_get_value(attributes, "atomic_num", value, status);
    if (status == 0) read (value, *) parse_tbm%atomic_num(ti)
    if (parse_tbm%atomic_num(ti) == 0) then
      parse_tbm%atomic_num(ti) = Atomic_Number_from_Mass(parse_tbm%atomic_mass(ti))
    endif
    call QUIP_FoX_get_value(attributes, "n_orbs", value, status);
    if (status == 0) read (value, *) parse_tbm%n_orbs(ti)
    call QUIP_FoX_get_value(attributes, "n_elecs", value, status);
    if (status == 0) read (value, *) parse_tbm%n_elecs(ti)
    call QUIP_FoX_get_value(attributes, "n_orb_sets", value, status);
    if (status == 0) read (value, *) parse_tbm%n_orb_sets(ti)
    if (parse_tbm%is_magnetic) then
      call QUIP_FoX_get_value(attributes, "lambda_sq_up", value, status);
      if (status == 0) read (value, *) parse_tbm%lambda_sq(ti,1)
      call QUIP_FoX_get_value(attributes, "lambda_sq_down", value, status);
      if (status == 0) read (value, *) parse_tbm%lambda_sq(ti,2)
    else
      call QUIP_FoX_get_value(attributes, "lambda_sq", value, status);
      if (status == 0) read (value, *) parse_tbm%lambda_sq(ti,:)
    endif

    if (allocated(parse_tbm%type_of_atomic_num)) deallocate(parse_tbm%type_of_atomic_num)
    allocate(parse_tbm%type_of_atomic_num(maxval(parse_tbm%atomic_num(:))))
    parse_tbm%type_of_atomic_num(:) = 0
    do ii=1, parse_tbm%n_types
      if (parse_tbm%atomic_num(ii) > 0) &
	parse_tbm%type_of_atomic_num(parse_tbm%atomic_num(ii)) = ii
    end do

    parse_cur_type_i = ti
    parse_cur_type_j = ti

  elseif (parse_in_tbm .and. name == 'orb_set_type') then

    call zero(parse_cur_data)

  elseif (parse_in_tbm .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "type1", value, status);
    if (status == 0) read (value, *) ti
    call QUIP_FoX_get_value(attributes, "type2", value, status);
    if (status == 0) read (value, *) tj
    call QUIP_FoX_get_value(attributes, "r_cut", value, status);
    if (status == 0) read (value, *) parse_tbm%r_cut(ti,tj)
    call QUIP_FoX_get_value(attributes, "screen_l", value, status);
    if (status == 0) read (value, *) parse_tbm%screen_l(ti,tj)
    if (parse_tbm%has_pair_repulsion) then
      call QUIP_FoX_get_value(attributes, "pair_rep_inner", value, status);
      if (status == 0) read (value, *) parse_tbm%pair_rep_inner(ti,tj)
      call QUIP_FoX_get_value(attributes, "pair_rep_outer", value, status);
      if (status == 0) read (value, *) parse_tbm%pair_rep_outer(ti,tj)
    endif
    parse_cur_type_i = ti
    parse_cur_type_j = tj

  elseif (parse_in_tbm .and. name == 'abcd') then

    call zero(parse_cur_data)

  elseif (parse_in_tbm .and. name == 'H_coeff') then

    call zero(parse_cur_data)

  elseif (parse_in_tbm .and. name == 'S_coeff') then

    call zero(parse_cur_data)

  endif

end subroutine TBM_startElement_handler

subroutine TBM_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  character(len=10240) :: val
  real(dp) :: abcd_v(32), coeff_v(160)
  integer k

  if (parse_in_tbm) then
    if (name == 'NRL_TB_params') then
      parse_in_tbm = .false.
    elseif (name == 'orb_set_type') then
      val = string(parse_cur_data)
      read (val, *) parse_tbm%orb_set_type(1:parse_tbm%n_orb_sets(parse_cur_type_i),parse_cur_type_i)
    elseif (name == 'abcd') then
      val = string(parse_cur_data)
      if (parse_tbm%is_magnetic) then
	if (parse_tbm%n_orb_sets(parse_cur_type_i) == 4) then ! f orbitals
	   read (val, *) abcd_v(1:32)
	   parse_tbm%abcd(1:4,1,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(1:4)
	   parse_tbm%abcd(1:4,1,parse_cur_type_i,parse_cur_type_j,2) = abcd_v(5:8)
	   parse_tbm%abcd(1:4,2,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(9:12)
	   parse_tbm%abcd(1:4,2,parse_cur_type_i,parse_cur_type_j,2) = abcd_v(13:16)
	   parse_tbm%abcd(1:4,3,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(17:20)
	   parse_tbm%abcd(1:4,3,parse_cur_type_i,parse_cur_type_j,2) = abcd_v(21:24)
	   parse_tbm%abcd(1:4,4,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(25:28)
	   parse_tbm%abcd(1:4,4,parse_cur_type_i,parse_cur_type_j,2) = abcd_v(29:32)
	else ! no f orbitals
	   read (val, *) abcd_v(1:24)
	   parse_tbm%abcd(1:4,1,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(1:4)
	   parse_tbm%abcd(1:4,1,parse_cur_type_i,parse_cur_type_j,2) = abcd_v(5:8)
	   parse_tbm%abcd(1:4,2,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(9:12)
	   parse_tbm%abcd(1:4,2,parse_cur_type_i,parse_cur_type_j,2) = abcd_v(13:16)
	   parse_tbm%abcd(1:4,3,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(17:20)
	   parse_tbm%abcd(1:4,3,parse_cur_type_i,parse_cur_type_j,2) = abcd_v(21:24)
	endif
      else
	if (parse_tbm%n_orb_sets(parse_cur_type_i) == 4) then ! f orbitals
	   read (val, *) abcd_v(1:16)
	   parse_tbm%abcd(1:4,1,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(1:4)
	   parse_tbm%abcd(1:4,2,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(5:8)
	   parse_tbm%abcd(1:4,3,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(9:12)
	   parse_tbm%abcd(1:4,4,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(13:16)
	else ! no f orbitals
	   read (val, *) abcd_v(1:12)
	   parse_tbm%abcd(1:4,1,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(1:4)
	   parse_tbm%abcd(1:4,2,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(5:8)
	   parse_tbm%abcd(1:4,3,parse_cur_type_i,parse_cur_type_j,1) = abcd_v(9:12)
	endif
      endif
    elseif (name == 'H_coeff') then
      val = string(parse_cur_data)
      if (parse_tbm%is_magnetic) then
	if (parse_tbm%n_orb_sets(parse_cur_type_i) == 4) then ! f orbitals
	  read (val, *) coeff_v(1:160)
	else ! no f orbitals
	  read (val, *) coeff_v(1:80)
	  coeff_v(81:160) = 0.0_dp
	endif
	do k=1, 20
	  parse_tbm%H_coeff(1:4,k,parse_cur_type_i,parse_cur_type_j,1) = coeff_v(8*(k-1)+1:8*(k-1)+4)
	  parse_tbm%H_coeff(1:4,k,parse_cur_type_i,parse_cur_type_j,2) = coeff_v(8*(k-1)+5:8*(k-1)+8)
	end do
      else
	if (parse_tbm%n_orb_sets(parse_cur_type_i) == 4) then ! f orbitals
	  read (val, *) coeff_v(1:80)
	else ! no f orbitals
	  read (val, *) coeff_v(1:40)
	  coeff_v(41:80) = 0.0_dp
	endif
	do k=1, 20
	  parse_tbm%H_coeff(1:4,k,parse_cur_type_i,parse_cur_type_j,1) = coeff_v(4*(k-1)+1:4*(k-1)+4)
	end do
      endif
    elseif (name == 'S_coeff') then
      val = string(parse_cur_data)
      if (parse_tbm%is_magnetic) then
	if (parse_tbm%n_orb_sets(parse_cur_type_i) == 4) then ! f orbitals
	  read (val, *) coeff_v(1:160)
	else ! no f orbitals
	  read (val, *) coeff_v(1:80)
	  coeff_v(81:160) = 0.0_dp
	endif
	do k=1, 20
	  parse_tbm%S_coeff(1:4,k,parse_cur_type_i,parse_cur_type_j,1) = coeff_v(8*(k-1)+1:8*(k-1)+4)
	  parse_tbm%S_coeff(1:4,k,parse_cur_type_i,parse_cur_type_j,2) = coeff_v(8*(k-1)+5:8*(k-1)+8)
	end do
      else
	if (parse_tbm%n_orb_sets(parse_cur_type_i) == 4) then ! f orbitals
	  read (val, *) coeff_v(1:80)
	else ! no f orbitals
	  read (val, *) coeff_v(1:40)
	  coeff_v(41:80) = 0.0_dp
	endif
	do k=1, 20
	  parse_tbm%S_coeff(1:4,k,parse_cur_type_i,parse_cur_type_j,1) = coeff_v(4*(k-1)+1:4*(k-1)+4)
	end do
      endif
    endif
  endif

end subroutine TBM_endElement_handler

subroutine TBModel_NRL_TB_read_params_xml(this, param_str)
  type(TBModel_NRL_TB), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type(xml_t) :: fxml

  if (len(trim(param_str)) <= 0) return

  parse_in_tbm = .false.
  parse_matched_label = .false.
  parse_tbm => this
  call Initialise(parse_cur_data)

  call open_xml_string(fxml, param_str)

  call parse(fxml, &
    characters_handler = TBM_characters_handler, &
    startElement_handler = TBM_startElement_handler, &
    endElement_handler = TBM_endElement_handler)

  call close_xml_t(fxml)

  if (this%n_types == 0) call system_abort("Called TBModel_NRL_TB_read_params_xml but no types are defined")

  this%cutoff = maxval(this%r_cut)

  call TBModel_NRL_TB_change_units(this)

  call finalise(parse_cur_data)

end subroutine TBModel_NRL_TB_read_params_xml

subroutine TBModel_NRL_TB_change_units(this)
  type(TBModel_NRL_TB), intent(inout) :: this

  integer :: i, j

  this%r_cut = this%r_cut * BOHR
  this%cutoff = this%cutoff * BOHR
  this%screen_l = this%screen_l * BOHR

  if (this%has_pair_repulsion) then
    this%pair_rep_inner = this%pair_rep_inner * BOHR
    this%pair_rep_outer = this%pair_rep_outer * BOHR
  endif

  this%lambda_sq = this%lambda_sq / BOHR
  this%abcd = this%abcd * RYDBERG

  this%H_coeff(MC_E,:,:,:,:) = this%H_coeff(MC_E,:,:,:,:) * RYDBERG
  this%H_coeff(MC_F,:,:,:,:) = this%H_coeff(MC_F,:,:,:,:) * RYDBERG / BOHR
  this%H_coeff(MC_FB,:,:,:,:) = this%H_coeff(MC_FB,:,:,:,:) * RYDBERG / BOHR**2
  this%H_coeff(MC_G_SQ,:,:,:,:) = this%H_coeff(MC_G_SQ,:,:,:,:) / BOHR

  do i=1, this%n_types
  do j=1, this%n_types
    if (i == j .and. this%overlap_zero_limit) then
      this%S_coeff(MC_E,:,i,j,:) = this%S_coeff(MC_E,:,i,j,:) / BOHR
      this%S_coeff(MC_F,:,i,j,:) = this%S_coeff(MC_F,:,i,j,:) / BOHR**2
      this%S_coeff(MC_FB,:,i,j,:) = this%S_coeff(MC_FB,:,i,j,:) / BOHR**3
      this%S_coeff(MC_G_SQ,:,i,j,:) = this%S_coeff(MC_G_SQ,:,i,j,:) / BOHR
    else
      this%S_coeff(MC_F,:,i,j,:) = this%S_coeff(MC_F,:,i,j,:) / BOHR
      this%S_coeff(MC_FB,:,i,j,:) = this%S_coeff(MC_FB,:,i,j,:) / BOHR**2
      this%S_coeff(MC_G_SQ,:,i,j,:) = this%S_coeff(MC_G_SQ,:,i,j,:) / BOHR
    endif
  end do
  end do

end subroutine TBModel_NRL_TB_change_units

subroutine TBModel_NRL_TB_Print(this,file)
  type(TBModel_NRL_TB),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  integer::i, j, ii, i_mag, n_sk_use, max_spdf

  if(this%n_types == 0) call System_Abort('TBModel_NRL_TB_Print: TBModel_NRL_TB structure not initialised')

  call Print('TBModel_NRL_TB n_types:' // this%n_types, file=file)

  call Print ('TBModel_NRL_TB: is_orthogonal ' // this%is_orthogonal // &
   ' is_magnetic '// this%is_magnetic// ' has_pair_repulsion '// this%has_pair_repulsion//  &
   ' overlap_zero_limit '// this%overlap_zero_limit// ' force_harrison_signs '// this%force_harrison_signs, file=file)

  if (any(this%orb_set_type == ORB_F)) then
     n_sk_use = 20
     max_spdf = SPDF_F
  else
     n_sk_use = 10
     max_spdf = SPDF_D
  endif

  do i=1, this%n_types
    call Print ('TBModel_NRL_TB: atomic num mass n_orbs n_elecs ' // this%atomic_num(i) //  " " // this%atomic_mass(i) //  &
      " " // this%n_orbs(i) //  " " // this%n_elecs(i), file=file)
    call verbosity_push_decrement()
    do j=1, this%n_types
      call Print ("TBModel_NRL_TB: types " //  i // " " //  j, file=file)
      call Print ('TBModel_NRL_TB: r_cut screen_l ' // this%r_cut(i,j) // " " //  this%screen_l(i,j), file=file)
      if (this%has_pair_repulsion) then
	call Print ('TBModel_NRL_TB: pair_rep inner outer ' // this%pair_rep_inner(i,j) // " " // this%pair_rep_outer(i,j), file=file)
      endif
      do i_mag=1, this%n_mag
	if (this%is_magnetic .and. i_mag == 1) call print("TBModel_NRL_TB: spin up")
	if (this%is_magnetic .and. i_mag == 2) call print("TBModel_NRL_TB: spin down")
	if (i == j) then
	  call Print ('TBModel_NRL_TB: lambda_sq ' // this%lambda_sq(i,i_mag), file=file)
	endif
	call Print ('TBModel_NRL_TB: a_s b_s c_s d_s ' // this%abcd(ABCD_A:ABCD_D,SPDF_S,i,j,i_mag), file=file)
	call Print ('TBModel_NRL_TB: a_p b_p c_p d_p ' // this%abcd(ABCD_A:ABCD_D,SPDF_P,i,j,i_mag), file=file)
	call Print ('TBModel_NRL_TB: a_d b_d c_d d_d ' // this%abcd(ABCD_A:ABCD_D,SPDF_D,i,j,i_mag), file=file)
	if (max_spdf >= SPDF_F) &
	   call Print ('TBModel_NRL_TB: a_f b_f c_f d_f ' // this%abcd(ABCD_A:ABCD_D,SPDF_F,i,j,i_mag), file=file)

	do ii=1, n_sk_use
	  call Print ('TBModel_NRL_TB: H e f fb g_sq ' // ii // " " //  this%H_coeff(MC_E:MC_G_SQ,ii,i,j,i_mag), file=file)
	end do

	do ii=1, n_sk_use
	  call Print ('TBModel_NRL_TB: S e f fb g_sq ' // ii // " " //  this%S_coeff(MC_E:MC_G_SQ,ii,i,j,i_mag), file=file)
	end do
      end do ! i_mag

    end do
    call verbosity_pop()
  end do

end subroutine TBModel_NRL_TB_Print

function TBModel_NRL_TB_n_orbs_of_Z(this, Z)
  type(TBModel_NRL_TB), intent(in) :: this
  integer, intent(in) :: Z
  integer TBModel_NRL_TB_n_orbs_of_Z

  TBModel_NRL_TB_n_orbs_of_Z = this%n_orbs(get_type(this%type_of_atomic_num,Z))
end function TBModel_NRL_TB_n_orbs_of_Z

function TBModel_NRL_TB_n_orb_sets_of_Z(this, Z)
  type(TBModel_NRL_TB), intent(in) :: this
  integer, intent(in) :: Z
  integer TBModel_NRL_TB_n_orb_sets_of_Z

  TBModel_NRL_TB_n_orb_sets_of_Z = this%n_orb_sets(get_type(this%type_of_atomic_num,Z))
end function TBModel_NRL_TB_n_orb_sets_of_Z

function TBModel_NRL_TB_n_orbs_of_orb_set_of_Z(this, Z, i_set)
  type(TBModel_NRL_TB), intent(in) :: this
  integer, intent(in) :: Z, i_set
  integer TBModel_NRL_TB_n_orbs_of_orb_set_of_Z

  TBModel_NRL_TB_n_orbs_of_orb_set_of_Z = N_ORBS_OF_SET(this%orb_set_type(i_set,get_type(this%type_of_atomic_num,Z)))

end function TBModel_NRL_TB_n_orbs_of_orb_set_of_Z

function TBModel_NRL_TB_orb_type_of_orb_set_of_Z(this, Z, i_set)
  type(TBModel_NRL_TB), intent(in) :: this
  integer, intent(in) :: Z, i_set
  integer TBModel_NRL_TB_orb_type_of_orb_set_of_Z

  TBModel_NRL_TB_orb_type_of_orb_set_of_Z = this%orb_set_type(i_set,get_type(this%type_of_atomic_num,Z))

end function TBModel_NRL_TB_orb_type_of_orb_set_of_Z

function TBModel_NRL_TB_n_elecs_of_Z(this, Z)
  type(TBModel_NRL_TB), intent(in) :: this
  integer, intent(in) :: Z
  real(dp) TBModel_NRL_TB_n_elecs_of_Z

  TBModel_NRL_TB_n_elecs_of_Z = this%n_elecs(get_type(this%type_of_atomic_num,Z))
end function TBModel_NRL_TB_n_elecs_of_Z

subroutine TBModel_NRL_TB_get_HS_blocks(this, at, at_i, at_j, dv_hat, dv_mag, b_H, b_S, i_mag)
  type(TBModel_NRL_TB), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_i, at_j
  real(dp), intent(in) :: dv_hat(3), dv_mag
  real(dp), intent(out) :: b_H(:,:), b_S(:,:)
  integer, intent(in), optional :: i_mag

  integer ti, tj, is, js, i_set, j_set
  integer i, j
  real(dp) SK_frad_H(N_SK), SK_frad_S(N_SK)
  real(dp) dv_hat_sq(3)

  ti = get_type(this%type_of_atomic_num,at%Z(at_i))
  tj = get_type(this%type_of_atomic_num,at%Z(at_j))

  b_H = 0.0_dp
  b_S = 0.0_dp

  if (dv_mag .feq. 0.0_dp) then
    i = 1
    do i_set = 1, this%n_orb_sets(ti)
    do is=1, N_ORBS_OF_SET(this%orb_set_type(i_set,ti))
      b_H(i,i) = onsite_function(this, at, at_i, this%orb_set_type(i_set,ti), i_mag)
      b_S(i,i) = 1.0D0
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
      SK_frad_H = 0.0_dp
      SK_frad_S = 0.0_dp
      call radial_functions(this, ti, tj, dv_mag, this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	SK_frad_H, SK_frad_S, i_mag)
      do js=1, N_ORBS_OF_SET(this%orb_set_type(j_set,tj))

	b_H(i,j) = angular_function(dv_hat, dv_hat_sq, this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	  is, js, SK_frad_H)

	if (this%is_orthogonal) then
	  b_S(i,j) = 0.0_dp
	else
	  b_S(i,j) = angular_function(dv_hat, dv_hat_sq, this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	    is, js, SK_frad_S)
	endif

	j = j + 1
      end do
      end do
      i = i + 1
    end do
    end do
  endif ! onsite

end subroutine TBModel_NRL_TB_get_HS_blocks

subroutine TBModel_NRL_TB_get_dHS_masks(this, at, at_ind, d_mask, od_mask)
  type(TBModel_NRL_TB), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_ind
  logical, intent(out), optional :: d_mask(:), od_mask(:)

  integer ji, j

  if (present(d_mask)) then
    if (at_ind < 0) then
      d_mask = .true.
    else
      d_mask = .false.
      d_mask(at_ind) = .true.
      do ji=1, n_neighbours(at, at_ind)
	j = neighbour(at, at_ind, ji)
	d_mask(j) = .true.
      end do
    endif
  endif

  if (present(od_mask)) then
    if (at_ind < 0) then
      od_mask = .true.
    else
      od_mask = .false.
      od_mask(at_ind) = .true.
    endif
  endif

end subroutine TBModel_NRL_TB_get_dHS_masks


function TBModel_NRL_TB_get_dHS_blocks(this, at, at_i, at_j, dv_hat, dv_mag, at_ind, b_dH, b_dS, i_mag)
  type(TBModel_NRL_TB), intent(inout) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_i, at_j
  real(dp), intent(in) :: dv_hat(3), dv_mag
  integer, intent(in) :: at_ind
  real(dp), intent(out) :: b_dH(:,:,:), b_dS(:,:,:)
  integer, intent(in), optional :: i_mag
  logical TBModel_NRL_TB_get_dHS_blocks

  integer ti, tj, is, js, i_set, j_set
  integer i, j
  real(dp) SK_frad_H(N_SK), SK_frad_S(N_SK), SK_dfrad_H(N_SK), SK_dfrad_S(N_SK)
  real(dp) dv_hat_sq(3)
  real(dp) virial_outerprod_fac

  if (at_i /= at_j) then ! off-diagonal block
    ! for force derivatives, only matching rows contribute
    if (at_ind > 0 .and. at_ind /= at_i .and. at_ind /= at_j) then
      TBModel_NRL_TB_get_dHS_blocks = .false.
      return
    endif
  else ! diagonal block
    ! only atoms within cutoff contribute
    if (dv_mag > this%cutoff) then
      TBModel_NRL_TB_get_dHS_blocks = .false.
      return
    endif
  endif

  ! force derivative off-diagonal masquerading as diagonal don't contribute
  if (at_ind > 0 .and. at_i == at_j .and. (dv_mag .fne. 0.0_dp)) then
    b_dH = 0.0_dp
    b_dS = 0.0_dp
    TBModel_NRL_TB_get_dHS_blocks = .true.
    return
  endif

  ti = get_type(this%type_of_atomic_num,at%Z(at_i))
  tj = get_type(this%type_of_atomic_num,at%Z(at_j))

  b_dH = 0.0_dp
  b_dS = 0.0_dp
  
  if (dv_mag .feq. 0.0_dp) then ! on-site
    i = 1
    do i_set = 1, this%n_orb_sets(ti)
    do is=1, N_ORBS_OF_SET(this%orb_set_type(i_set,ti))
      b_dH(i,i,:) = donsite_function(this, at, at_i, this%orb_set_type(i_set,ti), at_ind, i_mag)
      b_dS(i,i,:) = 0.0D0
    i = i + 1
    end do
    end do

  else ! hopping (even if on diagonal matrix block)

    dv_hat_sq = dv_hat**2

    i = 1
    do i_set=1, this%n_orb_sets(ti)
    do is=1, N_ORBS_OF_SET(this%orb_set_type(i_set,ti))
      j = 1
      do j_set=1, this%n_orb_sets(tj)
      call radial_functions(this, ti, tj, dv_mag, &
	this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	SK_frad_H, SK_frad_S, i_mag)
      call dradial_functions(this, ti, tj, dv_mag, &
	this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	SK_dfrad_H, SK_dfrad_S, i_mag)
      do js=1, N_ORBS_OF_SET(this%orb_set_type(j_set,tj))

	if (at_ind > 0) then
	  virial_outerprod_fac = 1.0_dp
	else
	  virial_outerprod_fac = -dv_mag*dv_hat(-at_ind)
	endif
	b_dH(i,j,:) = dangular_function(dv_mag, dv_hat, dv_hat_sq, &
	  this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	  is, js, SK_frad_H, SK_dfrad_H)*virial_outerprod_fac

	if (this%is_orthogonal) then
	  b_dS(i,j,:) = 0.0_dp
	else
	  b_dS(i,j,:) = dangular_function(dv_mag, dv_hat, dv_hat_sq, &
	    this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	    is, js, SK_frad_S, SK_dfrad_S)*virial_outerprod_fac
	endif

	j = j + 1
      end do
      end do
      i = i + 1
    end do
    end do

    if (at_ind == at_j) then
      b_dH = -b_dH
      if (.not. this%is_orthogonal) then
	b_dS = -b_dS
      endif
    endif

  endif ! hopping

  TBModel_NRL_TB_get_dHS_blocks = .true.
  return

end function TBModel_NRL_TB_get_dHS_blocks

function TBModel_NRL_TB_get_local_rep_E(this, at, i)
  type(TBModel_NRL_TB), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: TBModel_NRL_TB_get_local_rep_E

  TBModel_NRL_TB_get_local_rep_E = 0.0_dp

end function TBModel_NRL_TB_get_local_rep_E

function TBModel_NRL_TB_get_local_rep_E_force(this, at, i) result(force)
  type(TBModel_NRL_TB), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: force(3,at%N)

  force = 0.0_dp

end function TBModel_NRL_TB_get_local_rep_E_force

function TBModel_NRL_TB_get_local_rep_E_virial(this, at, i) result(virial)
  type(TBModel_NRL_TB), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: virial(3,3)

  virial = 0.0_dp

end function TBModel_NRL_TB_get_local_rep_E_virial

subroutine radial_functions(this, ti, tj, dv_mag, orb_set_type_i, orb_set_type_j, f_H, f_S, i_mag)
  type(TBModel_NRL_TB), intent(in) :: this
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: dv_mag
  integer, intent(in) :: orb_set_type_i, orb_set_type_j
  real(dp), intent(out) :: f_H(N_SK), f_S(N_SK)
  integer, intent(in), optional :: i_mag

  integer :: u_i_mag

  u_i_mag = optional_default(1, i_mag)

  if (this%n_mag == 1) u_i_mag = 1
  if (u_i_mag > this%n_mag) call system_abort("NRL_TB_Model radial_functions called for spin " // u_i_mag  // " max spin only " // this%n_mag)

  if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_S) then
    f_H(SK_SSS) = calc_SK_coeff_H(this, SK_SSS,ti,tj, dv_mag,u_i_mag)
    f_S(SK_SSS) = calc_SK_coeff_S_zero_limit(this, SK_SSS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j, u_i_mag)
    if (this%force_harrison_signs) f_H(SK_SSS) = harrison_sign_sss*abs(f_H(SK_SSS))
    if (this%force_harrison_signs) f_S(SK_SSS) = harrison_sign_sss*abs(f_S(SK_SSS))
  else if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_P) then
    f_H(SK_SPS) = calc_SK_coeff_H(this, SK_SPS,ti,tj, dv_mag,u_i_mag)
    f_S(SK_SPS) = calc_SK_coeff_S_zero_limit(this, SK_SPS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_H(SK_SPS) = harrison_sign_sps*abs(f_H(SK_SPS))
    if (this%force_harrison_signs) f_S(SK_SPS) = harrison_sign_sps*abs(f_S(SK_SPS))
  else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_S) then
    f_H(SK_SPS) = -calc_SK_coeff_H(this, SK_SPS,tj,ti, dv_mag,u_i_mag)
    f_S(SK_SPS) = -calc_SK_coeff_S_zero_limit(this, SK_SPS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_H(SK_SPS) = -harrison_sign_sps*abs(f_H(SK_SPS))
    if (this%force_harrison_signs) f_S(SK_SPS) = -harrison_sign_sps*abs(f_S(SK_SPS))
  else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_P) then
    f_H(SK_PPS) = calc_SK_coeff_H(this, SK_PPS,ti,tj, dv_mag,u_i_mag)
    f_H(SK_PPP) = calc_SK_coeff_H(this, SK_PPP,ti,tj, dv_mag,u_i_mag)
    f_S(SK_PPS) = calc_SK_coeff_S_zero_limit(this, SK_PPS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_PPP) = calc_SK_coeff_S_zero_limit(this, SK_PPP,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_H(SK_PPS) = harrison_sign_pps*abs(f_H(SK_PPS))
    if (this%force_harrison_signs) f_S(SK_PPS) = harrison_sign_pps*abs(f_S(SK_PPS))
    if (this%force_harrison_signs) f_H(SK_PPP) = harrison_sign_ppp*abs(f_H(SK_PPP))
    if (this%force_harrison_signs) f_S(SK_PPP) = harrison_sign_ppp*abs(f_S(SK_PPP))
  else if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_D) then
    f_H(SK_SDS) = calc_SK_coeff_H(this, SK_SDS,ti,tj, dv_mag,u_i_mag)
    f_S(SK_SDS) = calc_SK_coeff_S_zero_limit(this, SK_SDS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_H(SK_SDS) = harrison_sign_sds*abs(f_H(SK_SDS))
    if (this%force_harrison_signs) f_S(SK_SDS) = harrison_sign_sds*abs(f_S(SK_SDS))
  else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_S) then
    f_H(SK_SDS) = calc_SK_coeff_H(this, SK_SDS,tj,ti, dv_mag,u_i_mag)
    f_S(SK_SDS) = calc_SK_coeff_S_zero_limit(this, SK_SDS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_H(SK_SDS) = harrison_sign_sds*abs(f_H(SK_SDS))
    if (this%force_harrison_signs) f_S(SK_SDS) = harrison_sign_sds*abs(f_S(SK_SDS))
  else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_D) then
    f_H(SK_PDS) = calc_SK_coeff_H(this, SK_PDS,ti,tj, dv_mag,u_i_mag)
    f_H(SK_PDP) = calc_SK_coeff_H(this, SK_PDP,ti,tj, dv_mag,u_i_mag)
    f_S(SK_PDS) = calc_SK_coeff_S_zero_limit(this, SK_PDS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_PDP) = calc_SK_coeff_S_zero_limit(this, SK_PDP,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_H(SK_PDS) = harrison_sign_pds*abs(f_H(SK_PDS))
    if (this%force_harrison_signs) f_S(SK_PDS) = harrison_sign_pds*abs(f_S(SK_PDS))
    if (this%force_harrison_signs) f_H(SK_PDP) = harrison_sign_pdp*abs(f_H(SK_PDP))
    if (this%force_harrison_signs) f_S(SK_PDP) = harrison_sign_pdp*abs(f_S(SK_PDP))
  else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_P) then
    f_H(SK_PDS) = -calc_SK_coeff_H(this, SK_PDS,tj,ti, dv_mag,u_i_mag)
    f_H(SK_PDP) = -calc_SK_coeff_H(this, SK_PDP,tj,ti, dv_mag,u_i_mag)
    f_S(SK_PDS) = -calc_SK_coeff_S_zero_limit(this, SK_PDS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_PDP) = -calc_SK_coeff_S_zero_limit(this, SK_PDP,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_H(SK_PDS) = -harrison_sign_pds*abs(f_H(SK_PDS))
    if (this%force_harrison_signs) f_S(SK_PDS) = -harrison_sign_pds*abs(f_S(SK_PDS))
    if (this%force_harrison_signs) f_H(SK_PDP) = -harrison_sign_pdp*abs(f_H(SK_PDP))
    if (this%force_harrison_signs) f_S(SK_PDP) = -harrison_sign_pdp*abs(f_S(SK_PDP))
  else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_D) then
    f_H(SK_DDS) = calc_SK_coeff_H(this, SK_DDS,tj,ti, dv_mag,u_i_mag)
    f_H(SK_DDP) = calc_SK_coeff_H(this, SK_DDP,tj,ti, dv_mag,u_i_mag)
    f_H(SK_DDD) = calc_SK_coeff_H(this, SK_DDD,tj,ti, dv_mag,u_i_mag)
    f_S(SK_DDS) = calc_SK_coeff_S_zero_limit(this, SK_DDS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_DDP) = calc_SK_coeff_S_zero_limit(this, SK_DDP,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_DDD) = calc_SK_coeff_S_zero_limit(this, SK_DDD,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_H(SK_DDS) = harrison_sign_dds*abs(f_H(SK_DDS))
    if (this%force_harrison_signs) f_S(SK_DDS) = harrison_sign_dds*abs(f_S(SK_DDS))
    if (this%force_harrison_signs) f_H(SK_DDP) = harrison_sign_ddp*abs(f_H(SK_DDP))
    if (this%force_harrison_signs) f_S(SK_DDP) = harrison_sign_ddp*abs(f_S(SK_DDP))

  else if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_F) then
    f_H(SK_SFS) = calc_SK_coeff_H(this, SK_SFS,ti,tj, dv_mag,u_i_mag)
    f_S(SK_SFS) = calc_SK_coeff_S_zero_limit(this, SK_SFS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  else if (orb_set_type_i == ORB_F .and. orb_set_type_j == ORB_S) then
    f_H(SK_SFS) = -calc_SK_coeff_H(this, SK_SFS,tj,ti, dv_mag,u_i_mag)
    f_S(SK_SFS) = -calc_SK_coeff_S_zero_limit(this, SK_SFS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_F) then
    f_H(SK_PFS) = calc_SK_coeff_H(this, SK_PFS,ti,tj, dv_mag,u_i_mag)
    f_H(SK_PFP) = calc_SK_coeff_H(this, SK_PFP,ti,tj, dv_mag,u_i_mag)
    f_S(SK_PFS) = calc_SK_coeff_S_zero_limit(this, SK_PFS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_PFP) = calc_SK_coeff_S_zero_limit(this, SK_PFP,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  else if (orb_set_type_i == ORB_F .and. orb_set_type_j == ORB_P) then
    f_H(SK_PFS) = calc_SK_coeff_H(this, SK_PFS,tj,ti, dv_mag,u_i_mag)
    f_H(SK_PFP) = calc_SK_coeff_H(this, SK_PFP,tj,ti, dv_mag,u_i_mag)
    f_S(SK_PFS) = calc_SK_coeff_S_zero_limit(this, SK_PFS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_PFP) = calc_SK_coeff_S_zero_limit(this, SK_PFP,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_F) then
    f_H(SK_DFS) = calc_SK_coeff_H(this, SK_DFS,ti,tj, dv_mag,u_i_mag)
    f_H(SK_DFP) = calc_SK_coeff_H(this, SK_DFP,ti,tj, dv_mag,u_i_mag)
    f_H(SK_DFD) = calc_SK_coeff_H(this, SK_DFD,ti,tj, dv_mag,u_i_mag)
    f_S(SK_DFS) = calc_SK_coeff_S_zero_limit(this, SK_DFS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_DFP) = calc_SK_coeff_S_zero_limit(this, SK_DFP,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_DFD) = calc_SK_coeff_S_zero_limit(this, SK_DFD,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  else if (orb_set_type_i == ORB_F .and. orb_set_type_j == ORB_D) then
    f_H(SK_DFS) = -calc_SK_coeff_H(this, SK_DFS,tj,ti, dv_mag,u_i_mag)
    f_H(SK_DFP) = -calc_SK_coeff_H(this, SK_DFP,tj,ti, dv_mag,u_i_mag)
    f_H(SK_DFD) = -calc_SK_coeff_H(this, SK_DFD,tj,ti, dv_mag,u_i_mag)
    f_S(SK_DFS) = -calc_SK_coeff_S_zero_limit(this, SK_DFS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_DFP) = -calc_SK_coeff_S_zero_limit(this, SK_DFP,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_DFD) = -calc_SK_coeff_S_zero_limit(this, SK_DFD,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  else if (orb_set_type_i == ORB_F .and. orb_set_type_j == ORB_F) then
    f_H(SK_FFS) = calc_SK_coeff_H(this, SK_FFS,ti,tj, dv_mag,u_i_mag)
    f_H(SK_FFP) = calc_SK_coeff_H(this, SK_FFP,ti,tj, dv_mag,u_i_mag)
    f_H(SK_FFD) = calc_SK_coeff_H(this, SK_FFD,ti,tj, dv_mag,u_i_mag)
    f_H(SK_FFF) = calc_SK_coeff_H(this, SK_FFF,ti,tj, dv_mag,u_i_mag)
    f_S(SK_FFS) = calc_SK_coeff_S_zero_limit(this, SK_FFS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_FFP) = calc_SK_coeff_S_zero_limit(this, SK_FFP,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_FFD) = calc_SK_coeff_S_zero_limit(this, SK_FFD,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_S(SK_FFF) = calc_SK_coeff_S_zero_limit(this, SK_FFF,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  endif

end subroutine radial_functions

subroutine dradial_functions(this, ti, tj, dv_mag, orb_set_type_i, orb_set_type_j, f_dH, f_dS, i_mag)
  type(TBModel_NRL_TB), intent(in) :: this
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: dv_mag
  integer, intent(in) :: orb_set_type_i, orb_set_type_j
  real(dp), intent(out) :: f_dH(N_SK), f_dS(N_SK)
  integer, intent(in), optional :: i_mag

  integer :: u_i_mag

  u_i_mag = optional_default(1, i_mag)

  if (this%n_mag == 1) u_i_mag = 1
  if (u_i_mag > this%n_mag) call system_abort("NRL_TB_Model dradial_functions called for spin " // u_i_mag  // " max spin only " // this%n_mag)

  if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_S) then
    f_dH(SK_SSS) = calc_SK_coeff_H_d(this, SK_SSS,ti,tj, dv_mag,u_i_mag)
    f_dS(SK_SSS) = calc_SK_coeff_S_d_zero_limit(this, SK_SSS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_dH(SK_SSS) = harrison_sign_sss*abs(f_dH(SK_SSS))
    if (this%force_harrison_signs) f_dS(SK_SSS) = harrison_sign_sss*abs(f_dS(SK_SSS))
  else if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_P) then
    f_dH(SK_SPS) = calc_SK_coeff_H_d(this, SK_SPS,ti,tj, dv_mag,u_i_mag)
    f_dS(SK_SPS) = calc_SK_coeff_S_d_zero_limit(this, SK_SPS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_dH(SK_SPS) = harrison_sign_sps*abs(f_dH(SK_SPS))
    if (this%force_harrison_signs) f_dS(SK_SPS) = harrison_sign_sps*abs(f_dS(SK_SPS))
  else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_S) then
    f_dH(SK_SPS) = -calc_SK_coeff_H_d(this, SK_SPS,tj,ti, dv_mag,u_i_mag)
    f_dS(SK_SPS) = -calc_SK_coeff_S_d_zero_limit(this, SK_SPS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_dH(SK_SPS) = -harrison_sign_sps*abs(f_dH(SK_SPS))
    if (this%force_harrison_signs) f_dS(SK_SPS) = -harrison_sign_sps*abs(f_dS(SK_SPS))
  else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_P) then
    f_dH(SK_PPS) = calc_SK_coeff_H_d(this, SK_PPS,ti,tj, dv_mag,u_i_mag)
    f_dH(SK_PPP) = calc_SK_coeff_H_d(this, SK_PPP,ti,tj, dv_mag,u_i_mag)
    f_dS(SK_PPS) = calc_SK_coeff_S_d_zero_limit(this, SK_PPS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_PPP) = calc_SK_coeff_S_d_zero_limit(this, SK_PPP,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_dH(SK_PPS) = harrison_sign_pps*abs(f_dH(SK_PPS))
    if (this%force_harrison_signs) f_dS(SK_PPS) = harrison_sign_pps*abs(f_dS(SK_PPS))
    if (this%force_harrison_signs) f_dH(SK_PPP) = harrison_sign_ppp*abs(f_dH(SK_PPP))
    if (this%force_harrison_signs) f_dS(SK_PPP) = harrison_sign_ppp*abs(f_dS(SK_PPP))
  else if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_D) then
    f_dH(SK_SDS) = calc_SK_coeff_H_d(this, SK_SDS,ti,tj, dv_mag,u_i_mag)
    f_dS(SK_SDS) = calc_SK_coeff_S_d_zero_limit(this, SK_SDS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_dH(SK_SDS) = harrison_sign_sds*abs(f_dH(SK_SDS))
    if (this%force_harrison_signs) f_dS(SK_SDS) = harrison_sign_sds*abs(f_dS(SK_SDS))
  else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_S) then
    f_dH(SK_SDS) = calc_SK_coeff_H_d(this, SK_SDS,tj,ti, dv_mag,u_i_mag)
    f_dS(SK_SDS) = calc_SK_coeff_S_d_zero_limit(this, SK_SDS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_dH(SK_SDS) = harrison_sign_sds*abs(f_dH(SK_SDS))
    if (this%force_harrison_signs) f_dS(SK_SDS) = harrison_sign_sds*abs(f_dS(SK_SDS))
  else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_D) then
    f_dH(SK_PDS) = calc_SK_coeff_H_d(this, SK_PDS,ti,tj, dv_mag,u_i_mag)
    f_dH(SK_PDP) = calc_SK_coeff_H_d(this, SK_PDP,ti,tj, dv_mag,u_i_mag)
    f_dS(SK_PDS) = calc_SK_coeff_S_d_zero_limit(this, SK_PDS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_PDP) = calc_SK_coeff_S_d_zero_limit(this, SK_PDP,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_dH(SK_PDS) = harrison_sign_pds*abs(f_dH(SK_PDS))
    if (this%force_harrison_signs) f_dS(SK_PDS) = harrison_sign_pds*abs(f_dS(SK_PDS))
    if (this%force_harrison_signs) f_dH(SK_PDP) = harrison_sign_pdp*abs(f_dH(SK_PDP))
    if (this%force_harrison_signs) f_dS(SK_PDP) = harrison_sign_pdp*abs(f_dS(SK_PDP))
  else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_P) then
    f_dH(SK_PDS) = -calc_SK_coeff_H_d(this, SK_PDS,tj,ti, dv_mag,u_i_mag)
    f_dH(SK_PDP) = -calc_SK_coeff_H_d(this, SK_PDP,tj,ti, dv_mag,u_i_mag)
    f_dS(SK_PDS) = -calc_SK_coeff_S_d_zero_limit(this, SK_PDS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_PDP) = -calc_SK_coeff_S_d_zero_limit(this, SK_PDP,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_dH(SK_PDS) = -harrison_sign_pds*abs(f_dH(SK_PDS))
    if (this%force_harrison_signs) f_dS(SK_PDS) = -harrison_sign_pds*abs(f_dS(SK_PDS))
    if (this%force_harrison_signs) f_dH(SK_PDP) = -harrison_sign_pdp*abs(f_dH(SK_PDP))
    if (this%force_harrison_signs) f_dS(SK_PDP) = -harrison_sign_pdp*abs(f_dS(SK_PDP))
  else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_D) then
    f_dH(SK_DDS) = calc_SK_coeff_H_d(this, SK_DDS,tj,ti, dv_mag,u_i_mag)
    f_dH(SK_DDP) = calc_SK_coeff_H_d(this, SK_DDP,tj,ti, dv_mag,u_i_mag)
    f_dH(SK_DDD) = calc_SK_coeff_H_d(this, SK_DDD,tj,ti, dv_mag,u_i_mag)
    f_dS(SK_DDS) = calc_SK_coeff_S_d_zero_limit(this, SK_DDS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_DDP) = calc_SK_coeff_S_d_zero_limit(this, SK_DDP,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_DDD) = calc_SK_coeff_S_d_zero_limit(this, SK_DDD,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    if (this%force_harrison_signs) f_dH(SK_DDS) = harrison_sign_dds*abs(f_dH(SK_DDS))
    if (this%force_harrison_signs) f_dS(SK_DDS) = harrison_sign_dds*abs(f_dS(SK_DDS))
    if (this%force_harrison_signs) f_dH(SK_DDP) = harrison_sign_ddp*abs(f_dH(SK_DDP))
    if (this%force_harrison_signs) f_dS(SK_DDP) = harrison_sign_ddp*abs(f_dS(SK_DDP))

  else if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_F) then
    f_dH(SK_SFS) = calc_SK_coeff_H_d(this, SK_SFS,ti,tj, dv_mag,u_i_mag)
    f_dS(SK_SFS) = calc_SK_coeff_S_d_zero_limit(this, SK_SFS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  else if (orb_set_type_i == ORB_F .and. orb_set_type_j == ORB_S) then
    f_dH(SK_SFS) = -calc_SK_coeff_H_d(this, SK_SFS,tj,ti, dv_mag,u_i_mag)
    f_dS(SK_SFS) = -calc_SK_coeff_S_d_zero_limit(this, SK_SFS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_F) then
    f_dH(SK_PFS) = calc_SK_coeff_H_d(this, SK_PFS,ti,tj, dv_mag,u_i_mag)
    f_dH(SK_PFP) = calc_SK_coeff_H_d(this, SK_PFP,ti,tj, dv_mag,u_i_mag)
    f_dS(SK_PFS) = calc_SK_coeff_S_d_zero_limit(this, SK_PFS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_PFP) = calc_SK_coeff_S_d_zero_limit(this, SK_PFP,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  else if (orb_set_type_i == ORB_F .and. orb_set_type_j == ORB_P) then
    f_dH(SK_PFS) = calc_SK_coeff_H_d(this, SK_PFS,tj,ti, dv_mag,u_i_mag)
    f_dH(SK_PFP) = calc_SK_coeff_H_d(this, SK_PFP,tj,ti, dv_mag,u_i_mag)
    f_dS(SK_PFS) = calc_SK_coeff_S_d_zero_limit(this, SK_PFS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_PFP) = calc_SK_coeff_S_d_zero_limit(this, SK_PFP,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_F) then
    f_dH(SK_DFS) = calc_SK_coeff_H_d(this, SK_DFS,ti,tj, dv_mag,u_i_mag)
    f_dH(SK_DFP) = calc_SK_coeff_H_d(this, SK_DFP,ti,tj, dv_mag,u_i_mag)
    f_dH(SK_DFD) = calc_SK_coeff_H_d(this, SK_DFD,ti,tj, dv_mag,u_i_mag)
    f_dS(SK_DFS) = calc_SK_coeff_S_d_zero_limit(this, SK_DFS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_DFP) = calc_SK_coeff_S_d_zero_limit(this, SK_DFP,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_DFD) = calc_SK_coeff_S_d_zero_limit(this, SK_DFD,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  else if (orb_set_type_i == ORB_F .and. orb_set_type_j == ORB_D) then
    f_dH(SK_DFS) = -calc_SK_coeff_H_d(this, SK_DFS,tj,ti, dv_mag,u_i_mag)
    f_dH(SK_DFP) = -calc_SK_coeff_H_d(this, SK_DFP,tj,ti, dv_mag,u_i_mag)
    f_dH(SK_DFD) = -calc_SK_coeff_H_d(this, SK_DFD,tj,ti, dv_mag,u_i_mag)
    f_dS(SK_DFS) = -calc_SK_coeff_S_d_zero_limit(this, SK_DFS,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_DFP) = -calc_SK_coeff_S_d_zero_limit(this, SK_DFP,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_DFD) = -calc_SK_coeff_S_d_zero_limit(this, SK_DFD,tj,ti, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  else if (orb_set_type_i == ORB_F .and. orb_set_type_j == ORB_F) then
    f_dH(SK_FFS) = calc_SK_coeff_H_d(this, SK_FFS,ti,tj, dv_mag,u_i_mag)
    f_dH(SK_FFP) = calc_SK_coeff_H_d(this, SK_FFP,ti,tj, dv_mag,u_i_mag)
    f_dH(SK_FFD) = calc_SK_coeff_H_d(this, SK_FFD,ti,tj, dv_mag,u_i_mag)
    f_dH(SK_FFF) = calc_SK_coeff_H_d(this, SK_FFF,ti,tj, dv_mag,u_i_mag)
    f_dS(SK_FFS) = calc_SK_coeff_S_d_zero_limit(this, SK_FFS,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_FFP) = calc_SK_coeff_S_d_zero_limit(this, SK_FFP,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_FFD) = calc_SK_coeff_S_d_zero_limit(this, SK_FFD,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
    f_dS(SK_FFF) = calc_SK_coeff_S_d_zero_limit(this, SK_FFF,ti,tj, dv_mag, &
      ti == tj, orb_set_type_i == orb_set_type_j,u_i_mag)
  endif
end subroutine  dradial_functions

function calc_SK_coeff_H(this, SK_type, ti, tj, dist, i_mag)
  type(TBModel_NRL_TB), intent(in) :: this
  integer, intent(in) :: SK_type
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: dist
  integer, intent(in), optional :: i_mag
  real(dp) :: calc_SK_coeff_H

  integer :: u_i_mag

  u_i_mag = optional_default(1, i_mag)

  calc_SK_coeff_H = calc_SK_poly(this%H_coeff(1:4,SK_type, ti, tj, u_i_mag),  dist) 
  calc_SK_coeff_H = calc_SK_coeff_H * NRLTB_cutoff_function(this,  dist, ti, tj)
end function calc_SK_coeff_H

function calc_SK_coeff_H_d(this, SK_type, ti, tj, dist, i_mag)
  type(TBModel_NRL_TB), intent(in) :: this
  integer, intent(in) :: SK_type
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: dist
  integer, intent(in), optional :: i_mag
  real(dp) :: calc_SK_coeff_H_d

  real(dp) :: SK, SKd, co, cod
  integer :: u_i_mag

  u_i_mag = optional_default(1, i_mag)

  SKd = calc_SK_poly_deriv(this%H_coeff(1:4,SK_type, ti, tj,u_i_mag),  dist) 
  SK = calc_SK_poly(this%H_coeff(1:4,SK_type, ti, tj,u_i_mag),  dist) 
  co = NRLTB_cutoff_function(this,  dist, ti, tj)
  cod = NRLTB_cutoff_function_d(this,  dist, ti, tj)

  calc_SK_coeff_H_d = SK*cod + SKd*co
end function calc_SK_coeff_H_d

function calc_SK_coeff_S_zero_limit(this, SK_type, ti, tj, dist, ti_eq_tj, orb_ti_eq_orb_tj, i_mag)
  type(TBModel_NRL_TB), intent(in) :: this
  integer, intent(in) :: SK_type
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: dist
  logical, intent(in) :: ti_eq_tj, orb_ti_eq_orb_tj
  integer, intent(in), optional :: i_mag
  real(dp) :: calc_SK_coeff_S_zero_limit

  integer :: u_i_mag

  u_i_mag = optional_default(1, i_mag)

  calc_SK_coeff_S_zero_limit = calc_SK_poly_zero_limit(this%S_coeff(1:4,SK_type, ti, tj, u_i_mag),  dist, &
    this%overlap_zero_limit, ti_eq_tj, orb_ti_eq_orb_tj)
  calc_SK_coeff_S_zero_limit = calc_SK_coeff_S_zero_limit * NRLTB_cutoff_function(this,  dist, ti, tj)
end function calc_SK_coeff_S_zero_limit


function calc_SK_coeff_S_d_zero_limit(this, SK_type, ti, tj, dist, ti_eq_tj, orb_ti_eq_orb_tj, i_mag)
  type(TBModel_NRL_TB), intent(in) :: this
  integer, intent(in) :: SK_type
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: dist
  logical, intent(in) :: ti_eq_tj, orb_ti_eq_orb_tj
  integer, intent(in), optional :: i_mag
  real(dp) :: calc_SK_coeff_S_d_zero_limit

  real(dp) :: SK, SKd, co, cod
  integer :: u_i_mag

  u_i_mag = optional_default(1, i_mag)

  SK = calc_SK_poly_zero_limit(this%S_coeff(1:4,SK_type, ti, tj, u_i_mag),  dist, &
    this%overlap_zero_limit, ti_eq_tj, orb_ti_eq_orb_tj)
  SKd = calc_SK_poly_zero_limit_deriv(this%S_coeff(1:4,SK_type, ti, tj, u_i_mag),  dist, &
    this%overlap_zero_limit, ti_eq_tj, orb_ti_eq_orb_tj)
  co = NRLTB_cutoff_function(this,  dist, ti, tj)
  cod = NRLTB_cutoff_function_d(this,  dist, ti, tj)

  calc_SK_coeff_S_d_zero_limit = SK*cod + SKd*co
end function calc_SK_coeff_S_d_zero_limit

function calc_SK_poly(coeff, dist)
  real(dp), intent(in) :: coeff(4)
  real(dp), intent(in) :: dist
  real(dp) :: calc_SK_poly

  calc_SK_poly = (coeff(MC_E) + dist*(coeff(MC_F) + coeff(MC_FB)*dist))*exp(-coeff(MC_G_SQ)*dist) 
end function calc_SK_poly

function calc_SK_poly_deriv(coeff, dist)
  real(dp), intent(in) :: coeff(4)
  real(dp), intent(in) :: dist
  real(dp) :: calc_SK_poly_deriv

  real(dp) poly, poly_d, expv, expv_d

  expv = exp(-coeff(MC_G_SQ)*dist)
  expv_d = expv*(-coeff(MC_G_SQ))
  poly = (coeff(MC_E) + dist*(coeff(MC_F) + coeff(MC_FB)*dist))
  poly_d = coeff(MC_F) + 2.0_dp*coefF(MC_FB)*dist

  calc_SK_poly_deriv = expv*poly_d + expv_d*poly
end function calc_SK_poly_deriv

function calc_SK_poly_zero_limit(coeff, dist, overlap_zero_limit, ti_eq_tj, orb_ti_eq_orb_tj)
  real(dp), intent(in) :: coeff(4)
  real(dp), intent(in) :: dist
  logical, intent(in) :: overlap_zero_limit, ti_eq_tj, orb_ti_eq_orb_tj
  real(dp) :: calc_SK_poly_zero_limit

  real(dp) zero_limit

  if (overlap_zero_limit .and. ti_eq_tj .and. orb_ti_eq_orb_tj) then
    zero_limit = 1.0_dp
  else
    zero_limit = 0.0_dp
  end if

  if (overlap_zero_limit .and. ti_eq_tj) then
    calc_SK_poly_zero_limit = (zero_limit + dist*(coeff(MC_E) + dist*(coeff(MC_F) + &
      coeff(MC_FB)*dist)))*exp(-COEFF(MC_G_SQ)*dist)
  else
    calc_SK_poly_zero_limit = (coeff(MC_E) + dist*(coeff(MC_F) + coeff(MC_FB)*dist))* &
      exp(-coeff(MC_G_SQ)*dist)
  endif
end function calc_SK_poly_zero_limit

function calc_SK_poly_zero_limit_deriv(coeff, dist, overlap_zero_limit, ti_eq_tj, orb_ti_eq_orb_tj)
  real(dp), intent(in) :: coeff(4)
  real(dp), intent(in) :: dist
  logical, intent(in) :: overlap_zero_limit, ti_eq_tj, orb_ti_eq_orb_tj
  real(dp) :: calc_SK_poly_zero_limit_deriv

  real(dp) zero_limit
  real(dp) poly, poly_d, expv, expv_d

  if (overlap_zero_limit .and. ti_eq_tj .and. orb_ti_eq_orb_tj) then
    zero_limit = 1.0_dp
  else
    zero_limit = 0.0_dp
  end if

  expv = exp(-COEFF(MC_G_SQ)*dist)
  expv_d = expv*(-COEFF(MC_G_SQ))

  if (overlap_zero_limit .and. ti_eq_tj) then
    poly = (zero_limit + dist*(coeff(MC_E) + dist*(coeff(MC_F) + coeff(MC_FB)*dist)))
    poly_d = coeff(MC_E) + dist*(2.0_dp*coeff(MC_F) + 3.0_dp*coeff(MC_FB)*dist)
  else
    poly = (coeff(MC_E) + dist*(coeff(MC_F) + coeff(MC_FB)*dist))
    poly_d = coeff(MC_F) + 2.0_dp*coeff(MC_FB)*dist
  endif

  calc_SK_poly_zero_limit_deriv = expv * poly_d + expv_d * poly
end function calc_SK_poly_zero_limit_deriv

function onsite_function(this, at, at_i, orb_set_type,i_mag)
  type(TBModel_NRL_TB), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_i, orb_set_type
  integer, intent(in), optional :: i_mag
  real(dp) :: onsite_function

  integer :: ji, j, ti, tj
  real(dp) :: dist
  real(dp) :: density(this%n_types)
  real(dp) :: f_cut, expv
  integer :: u_i_mag

  u_i_mag = optional_default(1, i_mag)

  if (this%n_mag == 1) u_i_mag = 1
  if (u_i_mag > this%n_mag) call system_abort("NRL_TB_Model onsite_functions called for spin " // u_i_mag  // " max spin only " // this%n_mag)

  ti = get_type(this%type_of_atomic_num,at%Z(at_i))

  density = 0.0_dp
  do ji=1, n_neighbours(at, at_i)
    j = neighbour(at, at_i, ji, dist)
    tj = get_type(this%type_of_atomic_num,at%Z(j))

    f_cut = NRLTB_cutoff_function(this, dist, ti, tj)
    expv = exp(-this%lambda_sq(tj,u_i_mag)*dist)
    density(tj) = density(tj) + expv*f_cut
  end do

  onsite_function = 0.0_dp
  do tj=1, this%n_types
    onsite_function = onsite_function + onsite_poly(this%abcd(1:4,orb_set_type,tj,ti,u_i_mag), density(tj))
  end do

end function onsite_function

function donsite_function(this, at, at_i, orb_set_type, at_ind, i_mag)
  type(TBModel_NRL_TB), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_i, orb_set_type, at_ind
  integer, intent(in), optional :: i_mag
  real(dp) :: donsite_function(3)

  integer :: ji, j, ti, tj
  real(dp) :: dist
  real(dp) :: density(this%n_types), density_d(this%n_types,3)
  real(dp) :: f_cut, expv, f_cut_d, expv_d, dircos(3)
  real(dp) :: virial_outerprod_fac
  integer :: u_i_mag

  u_i_mag = optional_default(1, i_mag)

  if (this%n_mag == 1) u_i_mag = 1
  if (u_i_mag > this%n_mag) call system_abort("NRL_TB_Model donsite_functions called for spin " // u_i_mag  // " max spin only " // this%n_mag)

  ti = get_type(this%type_of_atomic_num,at%Z(at_i))

  density = 0.0_dp
  density_d = 0.0_dp
  do ji=1, n_neighbours(at, at_i)
    j = neighbour(at, at_i, ji, dist, cosines = dircos)
    tj = get_type(this%type_of_atomic_num,at%Z(j))

    f_cut = NRLTB_cutoff_function(this, dist, ti, tj)
    f_cut_d = NRLTB_cutoff_function_d(this, dist, ti, tj)
    expv = exp(-this%lambda_sq(tj,u_i_mag)*dist)
    expv_d = expv*(-this%lambda_sq(tj,u_i_mag))
    density(tj) = density(tj) + expv*f_cut
    if (at_i == at_ind) then
      density_d(tj,:) = density_d(tj,:) + (expv*f_cut_d + expv_d*f_cut)*(-dircos(:))
    else if (j == at_ind) then
      density_d(tj,:) = density_d(tj,:) + (expv*f_cut_d + expv_d*f_cut)*dircos(:)
    else if (at_ind < 0) then
      virial_outerprod_fac = -dist*dircos(-at_ind)
      density_d(tj,:) = density_d(tj,:) + (expv*f_cut_d + expv_d*f_cut)*(-dircos(:))*virial_outerprod_fac
    endif
  end do

  donsite_function = 0.0_dp
  do tj=1, this%n_types
    donsite_function(1) = donsite_function(1) + donsite_poly(this%abcd(1:4,orb_set_type,tj,ti,u_i_mag), &
						     density(tj), density_d(tj,1))
    donsite_function(2) = donsite_function(2) + donsite_poly(this%abcd(1:4,orb_set_type,tj,ti,u_i_mag), &
						     density(tj), density_d(tj,2))
    donsite_function(3) = donsite_function(3) + donsite_poly(this%abcd(1:4,orb_set_type,tj,ti,u_i_mag), &
						     density(tj), density_d(tj,3))
  end do

end function donsite_function

function onsite_poly(abcd, density)
  real(dp), intent(in) :: abcd(4)
  real(dp), intent(in) :: density
  real(dp) :: onsite_poly

  onsite_poly = abcd(ABCD_A) + abcd(ABCD_B)*density**(2.0_dp/3.0_dp) + &
    abcd(ABCD_C)*density**(4.0_dp/3.0_dp) + abcd(ABCD_D)*density**2

end function onsite_poly

function donsite_poly(abcd, density, density_d)
  real(dp), intent(in) :: abcd(4)
  real(dp), intent(in) :: density, density_d
  real(dp) :: donsite_poly

  if (density .feq. 0.0_dp) then
    donsite_poly = 0.0_dp
  else
    donsite_poly = (abcd(ABCD_B)*(2.0_dp/3.0_dp)*density**(-1.0_dp/3.0_dp) + &
      abcd(ABCD_C)*(4.0_dp/3.0_dp)*density**(1.0_dp/3.0_dp) + abcd(ABCD_D)*2.0_dp*density)*density_d
  endif

end function donsite_poly


function NRLTB_cutoff_function(this, r, ti, tj)
  type(TBModel_NRL_TB), intent(in) :: this
  real(dp), intent(in) :: r
  integer, intent(in) :: ti, tj
  real(dp) :: NRLTB_cutoff_function

  double precision screen_R0
  double precision cutoff_smooth
  double precision expv

  if (r .gt. 10.0_dp**(-4)) then

    screen_R0 = this%r_cut(ti,tj) - 5.0D0*abs(this%screen_l(ti,tj))

    expv = exp((r-screen_R0)/abs(this%screen_l(ti,tj)))

    cutoff_smooth = NRLTB_cutoff_func_smooth(this, r, ti, tj)
    NRLTB_cutoff_function = ( 1.0_dp/(1.0_dp+expv) ) * cutoff_smooth
  else
    NRLTB_cutoff_function = 0.0_dp
  end if

end function NRLTB_cutoff_function

function NRLTB_cutoff_function_d(this, r, ti, tj)
  type(TBModel_NRL_TB), intent(in) :: this
  real(dp), intent(in) :: r
  integer, intent(in) :: ti, tj
  real(dp) :: NRLTB_cutoff_function_d

  double precision screen_R0
  double precision cutoff_smooth, cutoff_smooth_d
  double precision expv, expv_d

  if (r .gt. 10.0_dp**(-4)) then

    screen_R0 = this%r_cut(ti,tj) - 5.0D0*abs(this%screen_l(ti,tj))

    expv = exp((r-screen_R0)/abs(this%screen_l(ti,tj)))
    expv_d = expv * 1.0_dp/abs(this%screen_l(ti,tj))

    cutoff_smooth = NRLTB_cutoff_func_smooth(this, r, ti, tj)
    cutoff_smooth_d = NRLTB_cutoff_func_smooth_d(this, r, ti, tj)

    NRLTB_cutoff_function_d = ( 1.0_dp/(1.0_dp+expv) ) * cutoff_smooth_d - &
      (1.0_dp/(1.0_dp+expv)**2)*expv_d * cutoff_smooth
  else
    NRLTB_cutoff_function_d = 0.0_dp
  end if

end function NRLTB_cutoff_function_d

function NRLTB_cutoff_func_smooth(this, r, ti, tj)
  type(TBModel_NRL_TB), intent(in) :: this
  real(dp), intent(in) :: r
  integer, intent(in) :: ti, tj
  real(dp) :: NRLTB_cutoff_func_smooth

  double precision R_MIN, R_MAX

  double precision PI
  parameter (PI = 3.14159265358979323846264338327950288_dp)


  if (this%screen_l(ti,tj) .lt. 0.0_dp) then
    NRLTB_cutoff_func_smooth = 1.0_dp
    return
  endif

  R_MAX = this%r_cut(ti,tj)
  R_MIN = R_MAX - abs(this%screen_l(ti,tj))

  if (r .lt. R_MIN) then 
    NRLTB_cutoff_func_smooth = 1.0_dp
  else if (r .gt. R_MAX) then
    NRLTB_cutoff_func_smooth = 0.0_dp
  else
    NRLTB_cutoff_func_smooth = 1.0_dp - (1.0_dp - cos( (r-R_MIN)*PI / (R_MAX-R_MIN) ))/2.0_dp
  end if

end function NRLTB_cutoff_func_smooth

function NRLTB_cutoff_func_smooth_d(this, r, ti, tj)
  type(TBModel_NRL_TB), intent(in) :: this
  real(dp), intent(in) :: r
  integer, intent(in) :: ti, tj
  real(dp) :: NRLTB_cutoff_func_smooth_d

  double precision R_MIN, R_MAX

  double precision PI
  parameter (PI = 3.14159265358979323846264338327950288_dp)


  if (this%screen_l(ti,tj) .lt. 0.0_dp) then
    NRLTB_cutoff_func_smooth_d = 1.0_dp
    return
  endif

  R_MAX = this%r_cut(ti,tj)
  R_MIN = R_MAX - abs(this%screen_l(ti,tj))

  if (r .lt. R_MIN) then 
    NRLTB_cutoff_func_smooth_d = 0.0_dp
  else if (r .gt. R_MAX) then
    NRLTB_cutoff_func_smooth_d = 0.0_dp
  else
    NRLTB_cutoff_func_smooth_d = -sin( (r-R_MIN)*PI / (R_MAX-R_MIN) )/2.0_dp * (PI/(R_MAX-R_MIN))
  end if

end function NRLTB_cutoff_func_smooth_d

end module TBModel_NRL_TB_module
