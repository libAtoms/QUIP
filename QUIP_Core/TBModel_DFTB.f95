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
!X TBModel_DFTB module
!X
!% Calculate energies using DFTB tight-binding model
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module TBModel_DFTB_module

use libatoms_module

use TB_Common_module
use QUIP_Common_module

implicit none
private

include 'TBModel_interface.h'

public :: TBModel_DFTB
type TBModel_DFTB
  integer :: n_types = 0
  character(len=FIELD_LENGTH) label

  real(dp) :: cutoff = 0.0_dp
  logical :: is_orthogonal = .false.

  integer, allocatable :: type_of_atomic_num(:)
  integer, allocatable :: n_orbs(:), n_elecs(:), n_orb_sets(:), orb_set_type(:,:)
  integer, allocatable :: atomic_num(:)

  !! DFTB parameters
  real(dp), allocatable:: E(:,:)

  real(dp), allocatable:: SK_cutoff(:,:), Vrep_cutoff(:,:)

  type(spline), allocatable :: H_spline(:,:,:), S_spline(:,:,:), Vrep_spline(:,:)

  logical :: has_default_fermi_E = .false., has_default_fermi_T = .true., has_default_band_width = .false., has_default_k_density = .false.
  real(dp) :: default_fermi_E, default_fermi_T = 0.001_dp, default_band_width, default_k_density

end type TBModel_DFTB

type(extendable_str), private, save :: parse_cur_data
logical, private :: parse_in_tbm, parse_matched_label
logical :: parse_in_H_spline, parse_in_S_spline, parse_in_Vrep_spline
integer  :: parse_cur_type_i, parse_cur_type_j, parse_cur_point_i
type (TBModel_DFTB), pointer :: parse_tbm
real(dp), allocatable :: parse_SK_r(:), parse_H_vals(:,:), parse_S_vals(:,:)
real(dp), allocatable :: parse_Vrep_r(:), parse_Vrep_vals(:)

interface Initialise
  module procedure TBModel_DFTB_Initialise_str
end interface Initialise

interface Finalise
  module procedure TBModel_DFTB_Finalise
end interface Finalise

interface Print
  module procedure TBModel_DFTB_Print
end interface Print

interface n_orbs_of_Z
  module procedure TBModel_DFTB_n_orbs_of_Z
end interface n_orbs_of_Z

interface n_orb_sets_of_Z
  module procedure TBModel_DFTB_n_orb_sets_of_Z
end interface n_orb_sets_of_Z

interface n_orbs_of_orb_set_of_Z
  module procedure TBModel_DFTB_n_orbs_of_orb_set_of_Z
end interface n_orbs_of_orb_set_of_Z

interface orb_type_of_orb_set_of_Z
  module procedure TBModel_DFTB_orb_type_of_orb_set_of_Z
end interface orb_type_of_orb_set_of_Z

interface n_elecs_of_Z
  module procedure TBModel_DFTB_n_elecs_of_Z
end interface n_elecs_of_Z

interface get_HS_blocks
  module procedure TBModel_DFTB_get_HS_blocks
end interface get_HS_blocks

interface get_dHS_masks
  module procedure TBModel_DFTB_get_dHS_masks
end interface get_dHS_masks

interface get_dHS_blocks
  module procedure TBModel_DFTB_get_dHS_blocks
end interface get_dHS_blocks

interface get_local_rep_E
  module procedure TBModel_DFTB_get_local_rep_E
end interface get_local_rep_E

interface get_local_rep_E_force
  module procedure TBModel_DFTB_get_local_rep_E_force
end interface get_local_rep_E_force

interface get_local_rep_E_virial
  module procedure TBModel_DFTB_get_local_rep_E_virial
end interface get_local_rep_E_virial

contains

subroutine TBModel_DFTB_Initialise_str(this, args_str, param_str)
  type(TBModel_DFTB), intent(inout) :: this
  character(len=*) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label)
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='TBModel_DFTB_Initialise_str args_str')) then
    call system_abort("TBModel_DFTB_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call TBModel_DFTB_read_params_xml(this, param_str)
end subroutine TBModel_DFTB_Initialise_str

subroutine TBModel_DFTB_Finalise(this)
  type(TBModel_DFTB), intent(inout) :: this

  integer ti, tj, i

  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%n_orbs)) deallocate(this%n_orbs)
  if (allocated(this%n_elecs)) deallocate(this%n_elecs)
  if (allocated(this%n_orb_sets)) deallocate(this%n_orb_sets)
  if (allocated(this%orb_set_type)) deallocate(this%orb_set_type)
  if (allocated(this%atomic_num)) deallocate(this%atomic_num)

  if (allocated(this%E)) deallocate(this%E)
  if (allocated(this%SK_cutoff)) deallocate(this%SK_cutoff)
  if (allocated(this%Vrep_cutoff)) deallocate(this%Vrep_cutoff)

  if (allocated(this%H_spline)) then
    do ti=1, size(this%H_spline,2)
    do tj=1, size(this%H_spline,3)
      do i=1, size(this%H_spline,1)
	call spline_finalise(this%H_spline(i,ti,tj))
      end do
    end do
    end do
    deallocate(this%H_spline)
  endif
  if (allocated(this%S_spline)) then
    do ti=1, size(this%S_spline,2)
    do tj=1, size(this%S_spline,3)
      do i=1, size(this%S_spline,1)
	call spline_finalise(this%S_spline(i,ti,tj))
      end do
    end do
    end do
    deallocate(this%S_spline)
  endif
  if (allocated(this%Vrep_spline)) then
    do ti=1, size(this%Vrep_spline,1)
    do tj=1, size(this%Vrep_spline,2)
      call spline_finalise(this%Vrep_spline(ti,tj))
    end do
    end do
    deallocate(this%Vrep_spline)
  endif

  this%n_types = 0
  this%label = ''
end subroutine TBModel_DFTB_Finalise

subroutine TBM_characters_handler(in)
  character(len=*), intent(in) :: in

  if (parse_in_tbm) then
    call concat(parse_cur_data, in, keep_lf=.false.)
  endif
end subroutine

subroutine TBM_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer status
  character(len=1024) :: value

  integer max_n_orb_sets, SK_npts, Vrep_npts
  integer ti, tj

  call zero(parse_cur_data)

  if (name == 'DFTB_params') then ! new DFTB stanza
    call print("DFTB_params startElement_handler", PRINT_NERD)

    if (parse_in_tbm) &
      call system_abort("IPModel_startElement_handler entered DFTB_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

    if (parse_matched_label) then
      call print("DFTB_params startElement_handler bailing because we already matched our label", PRINT_NERD)
      return ! we already found an exact match for this label
    endif

    call QUIP_FoX_get_value(attributes, 'label', value, status)
    if (status /= 0) value = ''

    call print("DFTB_params startElement_handler found xml label '"//trim(value)//"'", PRINT_NERD)

    if (len(trim(parse_tbm%label)) > 0) then ! we were passed in a label
      call print("DFTB_params startElement_handler was passed in label '"//trim(parse_tbm%label)//"'", PRINT_NERD)
      if (value == parse_tbm%label) then ! exact match
	call print("DFTB_params startElement_handler got label exact match", PRINT_NERD)
        parse_matched_label = .true.
        parse_in_tbm = .true.
      else ! no match
	call print("DFTB_params startElement_handler got label didn't match", PRINT_NERD)
        parse_in_tbm = .false.
      endif
    else ! no label passed in
      call print("DFTB_params startElement_handler was not passed in a label", PRINT_NERD)
      parse_in_tbm = .true.
      call print("DFTB_params startElement_handler was not passed in a label", PRINT_NERD)
    endif

    call print("DFTB_params startElement_handler parse_in_tbm " // parse_in_tbm, PRINT_NERD)

    if (parse_in_tbm) then
      call initialise(parse_cur_data)
      if (parse_tbm%n_types /= 0) then
	call print("DFTB_params startElement_handler finalising old data, restarting to parse new section", PRINT_NERD)
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

  elseif (parse_in_tbm .and. name == 'n_types') then
    call QUIP_FoX_get_value(attributes, "v", value, status);
    if (status == 0) read (value, *) parse_tbm%n_types

    allocate(parse_tbm%atomic_num(parse_tbm%n_types))
    parse_tbm%atomic_num = 0
    allocate(parse_tbm%n_orbs(parse_tbm%n_types))
    allocate(parse_tbm%n_elecs(parse_tbm%n_types))
    allocate(parse_tbm%n_orb_sets(parse_tbm%n_types))

    allocate(parse_tbm%H_spline(10,parse_tbm%n_types,parse_tbm%n_types))
    allocate(parse_tbm%S_spline(10,parse_tbm%n_types,parse_tbm%n_types))
    allocate(parse_tbm%Vrep_spline(parse_tbm%n_types,parse_tbm%n_types))

    allocate(parse_tbm%SK_cutoff(parse_tbm%n_types,parse_tbm%n_types))
    allocate(parse_tbm%Vrep_cutoff(parse_tbm%n_types,parse_tbm%n_types))

  elseif (parse_in_tbm .and. name == 'cutoff') then
    call QUIP_FoX_get_value(attributes, "v", value, status);
    if (status == 0) read (value, *) parse_tbm%cutoff
    parse_tbm%cutoff = parse_tbm%cutoff*BOHR
  elseif (parse_in_tbm .and. name == 'max_n_orb_sets') then
    call QUIP_FoX_get_value(attributes, "v", value, status);
    if (status == 0) read (value, *) max_n_orb_sets
    allocate(parse_tbm%orb_set_type(max_n_orb_sets,parse_tbm%n_types))
    allocate(parse_tbm%E(max_n_orb_sets,parse_tbm%n_types))
  elseif (parse_in_tbm .and. name == 'per_type_data') then
    call QUIP_FoX_get_value(attributes, "type", value, status);
    if (status == 0) read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status);
    if (status == 0) read (value, *) parse_tbm%atomic_num(ti)
    call QUIP_FoX_get_value(attributes, "n_orbs", value, status);
    if (status == 0) read (value, *) parse_tbm%n_orbs(ti)
    call QUIP_FoX_get_value(attributes, "n_orb_sets", value, status);
    if (status == 0) read (value, *) parse_tbm%n_orb_sets(ti)
    call QUIP_FoX_get_value(attributes, "n_elecs", value, status);
    if (status == 0) read (value, *) parse_tbm%n_elecs(ti)

    parse_cur_type_i = ti
    parse_cur_type_j = ti

    if (allocated(parse_tbm%type_of_atomic_num)) deallocate(parse_tbm%type_of_atomic_num)
    allocate(parse_tbm%type_of_atomic_num(maxval(parse_tbm%atomic_num(:))))
    parse_tbm%type_of_atomic_num(:) = 0
    do ti=1, parse_tbm%n_types
      if (parse_tbm%atomic_num(ti) > 0) &
	parse_tbm%type_of_atomic_num(parse_tbm%atomic_num(ti)) = ti
    end do

  elseif (parse_in_tbm .and. name == 'per_pair_data') then
    call QUIP_FoX_get_value(attributes, "type1", value, status);
    if (status == 0) read (value, *) ti
    call QUIP_FoX_get_value(attributes, "type2", value, status);
    if (status == 0) read (value, *) tj
    call QUIP_FoX_get_value(attributes, "SK_npts", value, status);
    if (status == 0) read (value, *) SK_npts
    call QUIP_FoX_get_value(attributes, "Vrep_npts", value, status);
    if (status == 0) read (value, *) Vrep_npts
    call QUIP_FoX_get_value(attributes, "SK_cutoff", value, status);
    if (status == 0) read (value, *) parse_tbm%SK_cutoff(ti,tj)
    parse_tbm%SK_cutoff(ti,tj) = parse_tbm%SK_cutoff(ti,tj)*BOHR
    call QUIP_FoX_get_value(attributes, "Vrep_cutoff", value, status);
    if (status == 0) read (value, *) parse_tbm%Vrep_cutoff(ti,tj)
    parse_tbm%Vrep_cutoff(ti,tj) = parse_tbm%Vrep_cutoff(ti,tj)*BOHR

    allocate(parse_SK_r(SK_npts))
    allocate(parse_H_vals(10,SK_npts))
    allocate(parse_S_vals(10,SK_npts))
    allocate(parse_Vrep_r(Vrep_npts))
    allocate(parse_Vrep_vals(Vrep_npts))

    parse_cur_type_i = ti
    parse_cur_type_j = tj

  elseif (parse_in_tbm .and. name == 'H_spline') then
    parse_in_H_spline = .true.
    parse_cur_point_i = 1
  elseif (parse_in_tbm .and. name == 'S_spline') then
    parse_in_S_spline = .true.
    parse_cur_point_i = 1
  elseif (parse_in_tbm .and. name == 'Vrep_spline') then
    parse_in_Vrep_spline = .true.
    parse_cur_point_i = 1
  elseif (parse_in_tbm .and. name == 'point') then
    if (parse_in_H_spline .or. parse_in_S_spline) then
      call QUIP_FoX_get_value(attributes, "r", value, status);
      if (status == 0) read (value, *) parse_SK_r(parse_cur_point_i)
      parse_SK_r(parse_cur_point_i) = parse_SK_r(parse_cur_point_i) * BOHR
    elseif (parse_in_Vrep_spline) then
      call QUIP_FoX_get_value(attributes, "r", value, status);
      if (status == 0) read (value, *) parse_Vrep_r(parse_cur_point_i)
      parse_Vrep_r(parse_cur_point_i) = parse_Vrep_r(parse_cur_point_i) * BOHR
    else
      call system_abort('Found point in DFTB_params but not in H_spline, S_spline, or Vrep_spline')
    endif
  endif

end subroutine TBM_startElement_handler

subroutine TBM_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  character(len=10240) :: val
  integer ii

  if (parse_in_tbm) then
    if (name == 'DFTB_params') then
      parse_in_tbm = .false.
    elseif (name == 'orb_set_type') then
      val = string(parse_cur_data)
      read (val, *) parse_tbm%orb_set_type(1:parse_tbm%n_orb_sets(parse_cur_type_i),parse_cur_type_i)
    elseif (name == 'E') then
      val = string(parse_cur_data)
      read (val, *) parse_tbm%E(1:parse_tbm%n_orb_sets(parse_cur_type_i),parse_cur_type_i)
      parse_tbm%E(:,parse_cur_type_i) = parse_tbm%E(:,parse_cur_type_i) * HARTREE
    elseif (name == 'per_pair_data') then
      if (allocated(parse_SK_r)) deallocate(parse_SK_r)
      if (allocated(parse_H_vals)) deallocate(parse_H_vals)
      if (allocated(parse_S_vals)) deallocate(parse_S_vals)
      if (allocated(parse_Vrep_r)) deallocate(parse_Vrep_r)
      if (allocated(parse_Vrep_vals)) deallocate(parse_Vrep_vals)
    elseif (parse_in_tbm .and. name == 'H_spline') then
      do ii=1, 10
	call spline_init(parse_tbm%H_spline(ii,parse_cur_type_i,parse_cur_type_j), parse_SK_r, parse_H_vals(ii,:), 0.0_dp, 0.0_dp)
      end do
      parse_in_H_spline = .false.
    elseif (parse_in_tbm .and. name == 'S_spline') then
      do ii=1, 10
	call spline_init(parse_tbm%S_spline(ii,parse_cur_type_i,parse_cur_type_j), parse_SK_r, parse_S_vals(ii,:), 0.0_dp, 0.0_dp)
      end do
      parse_in_S_spline = .false.
    elseif (parse_in_tbm .and. name == 'Vrep_spline') then
      call spline_init(parse_tbm%Vrep_spline(parse_cur_type_i,parse_cur_type_j), parse_Vrep_r, parse_Vrep_vals, 0.0_dp, 0.0_dp)
      parse_in_Vrep_spline = .false.
    elseif (name == 'point') then
      if (parse_in_H_spline) then
	val = string(parse_cur_data)
	parse_H_vals(1:10,parse_cur_point_i) = 0.0_dp
	read (val, *, end=1000) parse_H_vals(1:10,parse_cur_point_i)
1000    continue
	parse_H_vals(1:10,parse_cur_point_i) = parse_H_vals(1:10,parse_cur_point_i) * HARTREE
      elseif (parse_in_S_spline) then
	val = string(parse_cur_data)
	parse_S_vals(1:10,parse_cur_point_i) = 0.0_dp
	read (val, *, end=1010) parse_S_vals(1:10,parse_cur_point_i)
1010    continue
      elseif (parse_in_Vrep_spline) then
	val = string(parse_cur_data)
	read (val, *) parse_Vrep_vals(parse_cur_point_i)
	parse_Vrep_vals(parse_cur_point_i) = parse_Vrep_vals(parse_cur_point_i) * HARTREE
      else
	call system_abort('Found point in DFTB_params but not in H_spline, S_spline, or Vrep_Spline')
      endif
      parse_cur_point_i = parse_cur_point_i + 1
    endif
  endif

end subroutine TBM_endElement_handler


subroutine TBModel_DFTB_read_params_xml(this, param_str)
  type(TBModel_DFTB), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  integer :: err
  type(xml_t) :: fxml

  if (len(trim(param_str)) <= 0) return

  parse_tbm => this
  parse_in_tbm = .false.
  parse_matched_label = .false.
  parse_in_H_spline = .false.
  parse_in_S_spline = .false.
  parse_in_Vrep_spline = .false.

  err = increase_stack(5*len(param_str))
  if (err /= 0)then
    call print("TBModel_DfTB_Initialise_str Failed to increase stack before calling parse_xml", PRINT_ALWAYS)
  endif

  call open_xml_string(fxml, param_str)

  call parse(fxml, &
    characters_handler = TBM_characters_handler, &
    startElement_handler = TBM_startElement_handler, &
    endElement_handler = TBM_endElement_handler)

  call close_xml_t(fxml)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  type(inoutput), pointer :: in
!
!  integer ti, tj, Zi, Zj, tt, i, j, ii, jj
!  integer max_n_orb_sets
!
!  integer SK_npts, Vrep_npts
!  real(dp), allocatable :: SK_r(:), H_vals(:,:), S_vals(:,:), Vrep_r(:), Vrep_vals(:)
!
!  integer status
!  type(dictionary_t) :: attributes
!  type(xml_t) :: fxml
!  character(len=1024) :: value
!  integer :: ndata
!  character(len=10240) :: pcdata
!  integer cur_i
!
!  if (present(from_io)) then
!    in => from_io
!  else
!    allocate(in)
!    call Initialise(in, trim(dftb_default_file), INPUT)
!  endif
!
!  call open_xmlfile(in%unit, fxml, status)
!! call enable_debug(sax=.true.)
!  if (status /= 0) call system_abort('TBModel_DFTB_read_params_xml Cannot open xml file')
!
!  call rewind_xmlfile(fxml)
!
!  call get_node(fxml, path='//DFTB_params', attributes=attributes, status=status)
!  if (status /= 0) call system_abort('TBModel_DFTB_read_params_xml Cannot find /DFTB_params/header')
!
!  call get_node(fxml, path='//DFTB_params/n_types', attributes=attributes, status=status)
!  call QUIP_FoX_get_value(attributes, "v", value, status);
!  if (status == 0) read (value, *) this%n_types
!
!  call get_node(fxml, path='//DFTB_params/cutoff', attributes=attributes, status=status)
!  call QUIP_FoX_get_value(attributes, "v", value, status);
!  if (status == 0) read (value, *) this%cutoff
!  this%cutoff = this%cutoff*BOHR
!
!  call get_node(fxml, path='//DFTB_params/max_n_orb_sets', attributes=attributes, status=status)
!  call QUIP_FoX_get_value(attributes, "v", value, status);
!  if (status == 0) read (value, *) max_n_orb_sets
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
!  !! read DFTB per-type params
!  do tt=1, this%n_types
!    call get_node(fxml, path='//DFTB_params/per_type_data', attributes=attributes, status=status)
!    if (status /= 0) call system_abort("Couldn't find per_type_data for " // str(ti))
!
!    call QUIP_FoX_get_value(attributes, "type", value, status);
!    if (status == 0) then
!      read (value, *) ti
!    else
!      call system_abort ("failed to find type in per_type_data " // str(tt))
!    endif
!
!    call QUIP_FoX_get_value(attributes, "atomic_num", value, status);
!    if (status == 0) read (value, *) this%atomic_num(ti)
!
!    call QUIP_FoX_get_value(attributes, "n_orbs", value, status);
!    if (status == 0) read (value, *) this%n_orbs(ti)
!
!    call QUIP_FoX_get_value(attributes, "n_orb_sets", value, status);
!    if (status == 0) read (value, *) this%n_orb_sets(ti)
!
!    call QUIP_FoX_get_value(attributes, "n_elecs", value, status);
!    if (status == 0) read (value, *) this%n_elecs(ti)
!
!    call get_node(fxml, path='orb_set_type', attributes=attributes, pcdata=pcdata, status=status)
!    ndata = 0
!    call build_data_array(pcdata, this%orb_set_type(1:this%n_orb_sets(ti),ti), ndata)
!    if (ndata /= this%n_orb_sets(ti)) call system_abort ('Mismatch in amount of data reading orb_set_type '// trim(str(ti)) // &
!	' ' // trim(str(ndata)))
!
!    call get_node(fxml, path='E', attributes=attributes, pcdata=pcdata, status=status)
!    ndata = 0
!    call build_data_array(pcdata, this%E(1:this%n_orb_sets(ti),ti), ndata)
!    if (ndata /= this%n_orb_sets(ti)) call system_abort ('Mismatch in amount of data reading E '// trim(str(ti)) // &
!	' ' // trim(str(ndata)))
!    this%E(:,ti) = this%E(:,ti)*HARTREE
!
!  end do
!
!  allocate(this%type_of_atomic_num(maxval(this%atomic_num(:))))
!  this%type_of_atomic_num(:) = 0
!  do ti=1, this%n_types
!    this%type_of_atomic_num(this%atomic_num(ti)) = ti
!  end do
!
!  allocate(this%H_spline(10,this%n_types,this%n_types))
!  allocate(this%S_spline(10,this%n_types,this%n_types))
!  allocate(this%Vrep_spline(this%n_types,this%n_types))
!
!  allocate(this%SK_cutoff(this%n_types,this%n_types))
!  allocate(this%Vrep_cutoff(this%n_types,this%n_types))
!
!  call rewind_xmlfile(fxml)
!
!  do ii=1, this%n_types
!  do jj=1, this%n_types
!    call get_node(fxml, path='//DFTB_params/per_pair_data', attributes=attributes, status=status)
!    if (status /= 0) then
!      call system_abort("Can't find per_pair_data for " // str(ii) // ' ' // str(jj))
!    endif
!    call QUIP_FoX_get_value(attributes, "type1", value, status);
!    if (status == 0) read (value, *) ti
!    call QUIP_FoX_get_value(attributes, "type2", value, status);
!    if (status == 0) read (value, *) tj
!    call QUIP_FoX_get_value(attributes, "SK_npts", value, status);
!    if (status == 0) read (value, *) SK_npts
!    call QUIP_FoX_get_value(attributes, "Vrep_npts", value, status);
!    if (status == 0) read (value, *) Vrep_npts
!    call QUIP_FoX_get_value(attributes, "SK_cutoff", value, status);
!    if (status == 0) read (value, *) this%SK_cutoff(ti,tj)
!    this%SK_cutoff(ti,tj) = this%SK_cutoff(ti,tj)*BOHR
!    call QUIP_FoX_get_value(attributes, "Vrep_cutoff", value, status);
!    if (status == 0) read (value, *) this%Vrep_cutoff(ti,tj)
!    this%Vrep_cutoff(ti,tj) = this%Vrep_cutoff(ti,tj)*BOHR
!
!    allocate(SK_r(SK_npts))
!    allocate(H_vals(10,SK_npts))
!    allocate(S_vals(10,SK_npts))
!
!    allocate(Vrep_r(Vrep_npts))
!    allocate(Vrep_vals(Vrep_npts))
!
!    call mark_node(fxml, path='//DFTB_params/per_pair_data/H_spline', attributes=attributes, status=status)
!    if (status /= 0) then
!      call system_abort("Can't find H_spline for " // str(ii) // ' ' // str(jj))
!    endif
!    do i=1, SK_npts
!      call get_node(fxml, path='point', attributes=attributes, pcdata=pcdata, status=status)
!      if (status /= 0) then
!	call system_abort("Can't find H_spline/point for " // str(ii) // ' ' // str(jj))
!      endif
!
!      call QUIP_FoX_get_value(attributes, "r", value, status);
!      if (status == 0) read (value, *) SK_r(i)
!
!      ndata = 0
!      call build_data_array(pcdata, H_vals(1:10,i), ndata)
!    end do
!    SK_r = SK_r*BOHR
!    H_vals = H_vals*HARTREE
!    do i=1, 10
!      call spline_init(this%H_spline(i,ti,tj), SK_r, H_vals(i,:), 0.0_dp, 0.0_dp)
!    end do
!
!    call mark_node(fxml, path='//DFTB_params/per_pair_data/S_spline', attributes=attributes, status=status)
!    if (status /= 0) then
!      call system_abort("Can't find S_spline for " // str(ii) // ' ' // str(jj))
!    endif
!    do i=1, SK_npts
!      call get_node(fxml, path='point', attributes=attributes, pcdata=pcdata, status=status)
!      if (status /= 0) then
!	call system_abort("Can't find S_spline/point for " // str(ii) // ' ' // str(jj))
!      endif
!
!      call QUIP_FoX_get_value(attributes, "r", value, status);
!      if (status == 0) read (value, *) SK_r(i)
!
!      ndata = 0
!      call build_data_array(pcdata, S_vals(1:10,i), ndata)
!    end do
!    SK_r = SK_r*BOHR
!    do i=1, 10
!      call spline_init(this%S_spline(i,ti,tj), SK_r, S_vals(i,:), 0.0_dp, 0.0_dp)
!    end do
!
!    call mark_node(fxml, path='//DFTB_params/per_pair_data/Vrep_spline', attributes=attributes, status=status)
!    if (status /= 0) then
!      call system_abort("Can't find Vrep_spline for " // str(ii) // ' ' // str(jj))
!    endif
!    do i=1, Vrep_npts
!      call get_node(fxml, path='point', attributes=attributes, pcdata=pcdata, status=status)
!      if (status /= 0) then
!	call system_abort("Can't find Vrep_spline/point for " // str(ii) // ' ' // str(jj))
!      endif
!
!      call QUIP_FoX_get_value(attributes, "r", value, status);
!      if (status == 0) read (value, *) Vrep_r(i)
!
!      ndata = 0
!      call build_data_array(pcdata, Vrep_vals(i:i), ndata)
!    end do
!    Vrep_r = Vrep_r*BOHR
!    Vrep_vals = Vrep_vals*HARTREE
!    call spline_init(this%Vrep_spline(ti,tj), Vrep_r, Vrep_vals, 0.0_dp, 0.0_dp)
!
!    deallocate(SK_r)
!    deallocate(H_vals)
!    deallocate(S_vals)
!    deallocate(Vrep_r)
!    deallocate(Vrep_vals)
!  end do
!  end do
!
!end subroutine TBModel_DFTB_read_params_xml
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine TBModel_DFTB_Print(this,file)
  type(TBModel_DFTB),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  integer:: ti, tj, i

  if(this%n_types == 0) call System_Abort('TBModel_DFTB_Print: TBModel_DFTB structure not initialised')

  call Print('TBModel_DFTB Structure:', file=file)

  !! print DFTB params
  call Print ("TBModel_DFTB cutoff " // this%cutoff, file=file)
  do ti=1, this%n_types
    call Print ("TBModel_DFTB type " // ti // " Z " // this%atomic_num(ti) // " n_orbs " // this%n_orbs(ti) // &
      " n_elecs " // this%n_elecs(ti), file=file)
    if (this%n_orb_sets(ti) == 1) then
      call Print ("TBModel_DFTB E " // this%E(1,ti), file=file)
    else if (this%n_orb_sets(ti) == 2) then
      call Print ("TBModel_DFTB E " // this%E(1:2,ti), file=file)
    else 
      call Print ("TBModel_DFTB E " //this%E(1:3,ti), file=file)
    endif
  end do

  call verbosity_push_decrement(PRINT_NERD)
  if (current_verbosity() >= PRINT_NORMAL) then
    do ti=1, this%n_types
    do tj=1, this%n_types
      call Print ("TBModel_DFTB interaction " // ti // " " // tj // " Z " // this%atomic_num(ti) //  " " // &
	 this%atomic_num(tj), file=file)

      do i=1, 10
	call Print ("H_spline " // i, file=file)
	call Print (this%H_spline(i,ti,tj), file=file)
      end do
      do i=1, 10
	call Print ("S_spline " // i, file=file)
	call Print (this%S_spline(i,ti,tj), file=file)
      end do
      call Print ("Vrep_spline", file=file)
      call Print (this%Vrep_spline(ti,tj), file=file)

    end do
    end do
  endif
  call verbosity_pop()

end subroutine TBModel_DFTB_Print

function TBModel_DFTB_n_orbs_of_Z(this, Z)
  type(TBModel_DFTB), intent(in) :: this
  integer, intent(in) :: Z
  integer TBModel_DFTB_n_orbs_of_Z

  TBModel_DFTB_n_orbs_of_Z = this%n_orbs(get_type(this%type_of_atomic_num,Z))
end function TBModel_DFTB_n_orbs_of_Z

function TBModel_DFTB_n_orb_sets_of_Z(this, Z)
  type(TBModel_DFTB), intent(in) :: this
  integer, intent(in) :: Z
  integer TBModel_DFTB_n_orb_sets_of_Z

  TBModel_DFTB_n_orb_sets_of_Z = this%n_orb_sets(get_type(this%type_of_atomic_num,Z))
end function TBModel_DFTB_n_orb_sets_of_Z

function TBModel_DFTB_n_orbs_of_orb_set_of_Z(this, Z, i_set)
  type(TBModel_DFTB), intent(in) :: this
  integer, intent(in) :: Z, i_set
  integer TBModel_DFTB_n_orbs_of_orb_set_of_Z

  TBModel_DFTB_n_orbs_of_orb_set_of_Z = N_ORBS_OF_SET(this%orb_set_type(i_set,get_type(this%type_of_atomic_num,Z)))

end function TBModel_DFTB_n_orbs_of_orb_set_of_Z

function TBModel_DFTB_orb_type_of_orb_set_of_Z(this, Z, i_set)
  type(TBModel_DFTB), intent(in) :: this
  integer, intent(in) :: Z, i_set
  integer TBModel_DFTB_orb_type_of_orb_set_of_Z

  TBModel_DFTB_orb_type_of_orb_set_of_Z = this%orb_set_type(i_set,get_type(this%type_of_atomic_num,Z))

end function TBModel_DFTB_orb_type_of_orb_set_of_Z

function TBModel_DFTB_n_elecs_of_Z(this, Z)
  type(TBModel_DFTB), intent(in) :: this
  integer, intent(in) :: Z
  real(dp) TBModel_DFTB_n_elecs_of_Z

  TBModel_DFTB_n_elecs_of_Z = this%n_elecs(get_type(this%type_of_atomic_num,Z))
end function TBModel_DFTB_n_elecs_of_Z

subroutine TBModel_DFTB_get_HS_blocks(this, at, at_i, at_j, dv_hat, dv_mag, b_H, b_S, i_mag)
  type(TBModel_DFTB), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_i, at_j
  real(dp), intent(in) :: dv_hat(3), dv_mag
  real(dp), intent(out) :: b_H(:,:)
  real(dp), intent(out) :: b_S(:,:)
  integer, intent(in), optional :: i_mag

  integer ti, tj, is, js, i_set, j_set
  integer i, j
  real(dp) SK_frad_H(10), SK_frad_S(10)
  real(dp) dv_hat_sq(3)

  ti = get_type(this%type_of_atomic_num,at%Z(at_i))
  tj = get_type(this%type_of_atomic_num,at%Z(at_j))


  if (dv_mag .feq. 0.0_dp) then
    b_H = 0.0_dp
    b_S = 0.0_dp

    i = 1
    do i_set = 1, this%n_orb_sets(ti)
    do is=1, N_ORBS_OF_SET(this%orb_set_type(i_set,ti))
      b_H(i,i) = this%E(i_set,ti)
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
	is, js, SK_frad_H, SK_frad_S)
      do js=1, N_ORBS_OF_SET(this%orb_set_type(j_set,tj))

	b_H(i,j) = angular_function(dv_hat, dv_hat_sq, this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	  is, js, SK_frad_H)
	b_S(i,j) = angular_function(dv_hat, dv_hat_sq, this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	  is, js, SK_frad_S)
	j = j + 1
      end do
      end do
      i = i + 1
    end do
    end do
  endif

end subroutine TBModel_DFTB_get_HS_blocks

subroutine TBModel_DFTB_get_dHS_masks(this, at, at_ind, d_mask, od_mask)
  type(TBModel_DFTB), intent(in) :: this
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

end subroutine

function TBModel_DFTB_get_dHS_blocks(this, at, at_i, at_j, dv_hat, dv_mag, at_ind, b_dH, b_dS, i_mag)
  type(TBModel_DFTB), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_i, at_j
  real(dp), intent(in) :: dv_hat(3), dv_mag
  integer, intent(in) :: at_ind
  real(dp), intent(out) :: b_dH(:,:,:)
  real(dp), intent(out) :: b_dS(:,:,:)
  integer, intent(in), optional :: i_mag
  logical :: TBModel_DFTB_get_dHS_blocks

  integer ti, tj, is, js, i_set, j_set
  integer i, j
  real(dp) SK_frad_H(10), SK_frad_S(10)
  real(dp) SK_dfrad_H(10), SK_dfrad_S(10)
  real(dp) dv_hat_sq(3)
  real(dp) virial_outerprod_fac

  if ((at_ind > 0 .and. at_i /= at_ind .and. at_j /= at_ind) .or. &
      (at_ind > 0 .and. at_i == at_j) .or. &
      (dv_mag > this%cutoff) .or. &
      (dv_mag .feq. 0.0_dp)) then
    TBModel_DFTB_get_dHS_blocks = .false.
    return
  endif

  ti = get_type(this%type_of_atomic_num,at%Z(at_i))
  tj = get_type(this%type_of_atomic_num,at%Z(at_j))

  if (dv_mag > this%SK_cutoff(ti,tj)) then
    TBModel_DFTB_get_dHS_blocks = .false.
    return
  endif

  dv_hat_sq = dv_hat**2

  i = 1
  do i_set=1, this%n_orb_sets(ti)
  do is=1, N_ORBS_OF_SET(this%orb_set_type(i_set,ti))
    j = 1
    do j_set=1, this%n_orb_sets(tj)
    call radial_functions(this, ti, tj, dv_mag, this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
      is, js, SK_frad_H, SK_frad_S)
    call dradial_functions(this, ti, tj, dv_mag, this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
      is, js, SK_dfrad_H, SK_dfrad_S)
    do js=1, N_ORBS_OF_SET(this%orb_set_type(j_set,tj))

      if (at_ind > 0) then
	virial_outerprod_fac = 1.0_dp
      else
	virial_outerprod_fac = -dv_mag*dv_hat(-at_ind)
      end if
      b_dH(i,j,:) = dangular_function(dv_mag, dv_hat, dv_hat_sq, &
	this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	is, js, SK_frad_H, SK_dfrad_H)*virial_outerprod_fac
      b_dS(i,j,:) = dangular_function(dv_mag, dv_hat, dv_hat_sq, &
	this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
	is, js, SK_frad_S, SK_dfrad_S)*virial_outerprod_fac
      j = j + 1
    end do
    end do
    i = i + 1
  end do
  end do

  if (at_ind == at_j) then
    b_dH = -b_dH
    b_dS = -b_dS
  end if

  TBModel_DFTB_get_dHS_blocks = .true.
  return


end function TBModel_DFTB_get_dHS_blocks
 
subroutine radial_functions(this, ti, tj, dv_mag, orb_set_type_i, orb_set_type_j, is, js, f_H, f_S)
   type(TBModel_DFTB), intent(in) :: this
   integer, intent(in) :: ti, tj
   real(dp), intent(in) :: dv_mag
   integer, intent(in) :: orb_set_type_i, orb_set_type_j
   integer, intent(in) :: is, js
   real(dp), intent(out) :: f_H(10), f_S(10)

   if (dv_mag > this%SK_cutoff(ti,tj)) then
     f_H = 0.0_dp
     f_S = 0.0_dp
     return
   endif

   if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_S) then
     f_H(SK_SSS) = spline_value(this%H_spline(1,ti,tj), dv_mag)
     f_S(SK_SSS) = spline_value(this%S_spline(1,ti,tj), dv_mag)
   else if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_P) then
     f_H(SK_SPS) = spline_value(this%H_spline(2,ti,tj), dv_mag)
     f_S(SK_SPS) = spline_value(this%S_spline(2,ti,tj), dv_mag)
   else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_S) then
     f_H(SK_SPS) = -spline_value(this%H_spline(2,tj,ti), dv_mag)
     f_S(SK_SPS) = -spline_value(this%S_spline(2,tj,ti), dv_mag)
   else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_P) then
     f_H(SK_PPS) = spline_value(this%H_spline(3,ti,tj), dv_mag)
     f_S(SK_PPS) = spline_value(this%S_spline(3,ti,tj), dv_mag)
     f_H(SK_PPP) = spline_value(this%H_spline(4,ti,tj), dv_mag)
     f_S(SK_PPP) = spline_value(this%S_spline(4,ti,tj), dv_mag)
   else if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_D) then
     f_H(SK_SDS) = spline_value(this%H_spline(5,ti,tj), dv_mag)
     f_S(SK_SDS) = spline_value(this%S_spline(5,ti,tj), dv_mag)
   else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_S) then
     f_H(SK_SDS) = spline_value(this%H_spline(5,tj,ti), dv_mag)
     f_S(SK_SDS) = spline_value(this%S_spline(5,tj,ti), dv_mag)
   else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_D) then
     f_H(SK_PDS) = spline_value(this%H_spline(6,ti,tj), dv_mag)
     f_S(SK_PDS) = spline_value(this%S_spline(6,ti,tj), dv_mag)
     f_H(SK_PDP) = spline_value(this%H_spline(7,ti,tj), dv_mag)
     f_S(SK_PDP) = spline_value(this%S_spline(7,ti,tj), dv_mag)
   else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_P) then
     f_H(SK_PDS) = -spline_value(this%H_spline(6,tj,ti), dv_mag)
     f_S(SK_PDS) = -spline_value(this%S_spline(6,tj,ti), dv_mag)
     f_H(SK_PDP) =  spline_value(this%H_spline(7,tj,ti), dv_mag)
     f_S(SK_PDP) =  spline_value(this%S_spline(7,tj,ti), dv_mag)
   else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_D) then
     f_H(SK_DDS) = spline_value(this%H_spline(9,ti,tj), dv_mag)
     f_S(SK_DDS) = spline_value(this%S_spline(8,ti,tj), dv_mag)
     f_H(SK_DDP) = spline_value(this%H_spline(9,ti,tj), dv_mag)
     f_S(SK_DDP) = spline_value(this%S_spline(9,ti,tj), dv_mag)
     f_H(SK_DDD) = spline_value(this%H_spline(10,ti,tj), dv_mag)
     f_S(SK_DDD) = spline_value(this%S_spline(10,ti,tj), dv_mag)
   endif
 end subroutine radial_functions

 subroutine dradial_functions(this, ti, tj, dv_mag, orb_set_type_i, orb_set_type_j, is, js, f_dH, f_dS)
   type(TBModel_DFTB), intent(in) :: this
   integer, intent(in) :: ti, tj
   real(dp), intent(in) :: dv_mag
   integer, intent(in) :: orb_set_type_i, orb_set_type_j
   integer, intent(in) :: is, js
   real(dp), intent(out) :: f_dH(10), f_dS(10)

   if (dv_mag > this%SK_cutoff(ti,tj)) then
     f_dH = 0.0_dp
     f_dS = 0.0_dp
     return
   endif

   if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_S) then
     f_dH(SK_SSS) = spline_deriv(this%H_spline(1,ti,tj), dv_mag)
     f_dS(SK_SSS) = spline_deriv(this%S_spline(1,ti,tj), dv_mag)
   else if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_P) then
     f_dH(SK_SPS) = spline_deriv(this%H_spline(2,ti,tj), dv_mag)
     f_dS(SK_SPS) = spline_deriv(this%S_spline(2,ti,tj), dv_mag)
   else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_S) then
     f_dH(SK_SPS) = -spline_deriv(this%H_spline(2,tj,ti), dv_mag)
     f_dS(SK_SPS) = -spline_deriv(this%S_spline(2,tj,ti), dv_mag)
   else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_P) then
     f_dH(SK_PPS) = spline_deriv(this%H_spline(3,ti,tj), dv_mag)
     f_dS(SK_PPS) = spline_deriv(this%S_spline(3,ti,tj), dv_mag)
     f_dH(SK_PPP) = spline_deriv(this%H_spline(4,ti,tj), dv_mag)
     f_dS(SK_PPP) = spline_deriv(this%S_spline(4,ti,tj), dv_mag)
   else if (orb_set_type_i == ORB_S .and. orb_set_type_j == ORB_D) then
     f_dH(SK_SDS) = spline_deriv(this%H_spline(5,ti,tj), dv_mag)
     f_dS(SK_SDS) = spline_deriv(this%S_spline(5,ti,tj), dv_mag)
   else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_S) then
     f_dH(SK_SDS) = spline_deriv(this%H_spline(5,tj,ti), dv_mag)
     f_dS(SK_SDS) = spline_deriv(this%S_spline(5,tj,ti), dv_mag)
   else if (orb_set_type_i == ORB_P .and. orb_set_type_j == ORB_D) then
     f_dH(SK_PDS) = spline_deriv(this%H_spline(6,ti,tj), dv_mag)
     f_dS(SK_PDS) = spline_deriv(this%S_spline(6,ti,tj), dv_mag)
     f_dH(SK_PDP) = spline_deriv(this%H_spline(7,ti,tj), dv_mag)
     f_dS(SK_PDP) = spline_deriv(this%S_spline(7,ti,tj), dv_mag)
   else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_P) then
     f_dH(SK_PDS) = -spline_deriv(this%H_spline(6,tj,ti), dv_mag)
     f_dS(SK_PDS) = -spline_deriv(this%S_spline(6,tj,ti), dv_mag)
     f_dH(SK_PDP) =  spline_deriv(this%H_spline(7,tj,ti), dv_mag)
     f_dS(SK_PDP) =  spline_deriv(this%S_spline(7,tj,ti), dv_mag)
   else if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_D) then
     f_dH(SK_DDS) = spline_deriv(this%H_spline(9,ti,tj), dv_mag)
     f_dS(SK_DDS) = spline_deriv(this%S_spline(8,ti,tj), dv_mag)
     f_dH(SK_DDP) = spline_deriv(this%H_spline(9,ti,tj), dv_mag)
     f_dS(SK_DDP) = spline_deriv(this%S_spline(9,ti,tj), dv_mag)
     f_dH(SK_DDD) = spline_deriv(this%H_spline(10,ti,tj), dv_mag)
     f_dS(SK_DDD) = spline_deriv(this%S_spline(10,ti,tj), dv_mag)
   endif

end subroutine dradial_functions
 
function TBModel_DFTB_get_local_rep_E(this, at, i)
  type(TBModel_DFTB), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: TBModel_DFTB_get_local_rep_E

  real(dp) E, dist
  integer ji, j, ti, tj

  TBModel_DFTB_get_local_rep_E = 0.0_dp

  ti = get_type(this%type_of_atomic_num, at%Z(i))
  do ji=1, atoms_n_neighbours(at, i)
    j = atoms_neighbour(at, i, ji, dist)
    if (dist .feq. 0.0_dp) cycle
    tj = get_type(this%type_of_atomic_num, at%Z(j))

    if (dist < this%Vrep_cutoff(ti,tj)) then
      E = spline_value(this%Vrep_spline(ti,tj), dist)
      TBModel_DFTB_get_local_rep_E = TBModel_DFTB_get_local_rep_E + E/2.0_dp
    endif
  end do

end function TBModel_DFTB_get_local_rep_E

function TBModel_DFTB_get_local_rep_E_force(this, at, i) result(force)
  type(TBModel_DFTB), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: force(3,at%N)

  real(dp) dE_dr, dist, dv_hat(3)
  integer ji, j, ti, tj

  force = 0.0_dp

  ti = get_type(this%type_of_atomic_num, at%Z(i))
  do ji=1, atoms_n_neighbours(at, i)
    j = atoms_neighbour(at, i, ji, dist, cosines = dv_hat)
    if (dist .feq. 0.0_dp) cycle
    tj = get_type(this%type_of_atomic_num, at%Z(j))

    if (dist < this%Vrep_cutoff(ti,tj)) then
      dE_dr = spline_deriv(this%Vrep_spline(ti,tj),dist)
      force(:,i) = force(:,i) + dE_dr*dv_hat(:)/2.0_dp
      force(:,j) = force(:,j) - dE_dr*dv_hat(:)/2.0_dp
    end if
  end do

end function TBModel_DFTB_get_local_rep_E_force

function TBModel_DFTB_get_local_rep_E_virial(this, at, i) result(virial)
  type(TBModel_DFTB), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: virial(3,3)

  real(dp) dE_dr, dist, dv_hat(3)
  integer ji, j, ti, tj

  virial = 0.0_dp

  ti = get_type(this%type_of_atomic_num, at%Z(i))
  do ji=1, atoms_n_neighbours(at, i)
    j = atoms_neighbour(at, i, ji, dist, cosines = dv_hat)
    if (dist .feq. 0.0_dp) cycle
    tj = get_type(this%type_of_atomic_num, at%Z(j))

    if (dist < this%Vrep_cutoff(ti,tj)) then
      dE_dr = spline_deriv(this%Vrep_spline(ti,tj),dist)
      virial = virial - dE_dr*(dv_hat .outer. dv_hat) * dist / 2.0_dp
    end if
  end do

end function TBModel_DFTB_get_local_rep_E_virial

end module TBModel_DFTB_module
