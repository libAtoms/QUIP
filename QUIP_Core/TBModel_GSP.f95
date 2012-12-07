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
!X TBModel_GSP module
!X
!% Calculate energies using Goodwin-Skinner-Pettifor tight-binding model
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module TBModel_GSP_module

use libatoms_module

use TB_Common_module
use QUIP_Common_module

implicit none
private

integer, parameter :: max_n_orb_sets = 3

include 'TBModel_interface.h'

public :: TBModel_GSP
type TBModel_GSP
  integer :: n_types = 0
  character(len=STRING_LENGTH) label

  real(dp) :: cutoff = 0.0_dp
  real(dp) :: cutoff_H = 0.0_dp
  logical :: is_orthogonal = .true.

  integer, allocatable :: type_of_atomic_num(:)
  integer, allocatable :: n_orbs(:), n_orb_sets(:), orb_set_type(:,:)
  real(dp), allocatable :: n_elecs(:)
  integer, allocatable :: atomic_num(:)

!       GSP parameters
  real(dp), allocatable :: H_coeff(:,:,:)
  real(dp), allocatable :: O_coeff(:,:,:)
  real(dp), allocatable :: Atau(:,:,:)
  real(dp), allocatable :: deltatau0(:,:,:)
  real(dp), allocatable :: r0(:,:), rc(:,:), na(:,:), nb(:,:), nc(:,:)

  real(dp), allocatable :: B(:,:), Rcore(:,:), lambda0(:), C(:), nu(:), m(:)
  real(dp), allocatable :: Rtail(:), Rcut(:) 

  real(dp), allocatable :: Ai(:,:,:), Ri(:,:,:)

  real(dp) :: tailx0

  real(dp), allocatable :: E(:,:)

  type(spline), allocatable :: H_tail_spline(:,:)
  type(spline), allocatable :: Vrep_env_emb_tail_spline(:)

  logical :: has_default_fermi_E = .false., has_default_fermi_T = .true., has_default_band_width = .false., has_default_k_density=.false.
  real(dp) :: default_fermi_E, default_fermi_T = 0.001_dp, default_band_width, default_k_density

end type

integer, private :: parse_cur_type
logical, private :: parse_in_tbm, parse_matched_label
type(TBModel_GSP), pointer :: parse_tbm

interface Initialise
  module procedure TBModel_GSP_Initialise_str
end interface Initialise

interface Finalise
  module procedure TBModel_GSP_Finalise
end interface Finalise

interface Print
  module procedure TBModel_GSP_Print
end interface Print

interface n_orbs_of_Z
  module procedure TBModel_GSP_n_orbs_of_Z
end interface n_orbs_of_Z

interface n_orb_sets_of_Z
  module procedure TBModel_GSP_n_orb_sets_of_Z
end interface n_orb_sets_of_Z

interface n_orbs_of_orb_set_of_Z
  module procedure TBModel_GSP_n_orbs_of_orb_set_of_Z
end interface n_orbs_of_orb_set_of_Z

interface orb_type_of_orb_set_of_Z
  module procedure TBModel_GSP_orb_type_of_orb_set_of_Z
end interface orb_type_of_orb_set_of_Z

interface n_elecs_of_Z
  module procedure TBModel_GSP_n_elecs_of_Z
end interface n_elecs_of_Z

interface get_HS_blocks
  module procedure TBModel_GSP_get_HS_blocks
end interface get_HS_blocks

interface get_dHS_masks
  module procedure TBModel_GSP_get_dHS_masks
end interface get_dHS_masks

interface get_dHS_blocks
  module procedure TBModel_GSP_get_dHS_blocks
end interface get_dHS_blocks

interface get_local_rep_E
  module procedure TBModel_GSP_get_local_rep_E
end interface get_local_rep_E

interface get_local_rep_E_force
  module procedure TBModel_GSP_get_local_rep_E_force
end interface get_local_rep_E_force

interface get_local_rep_E_virial
  module procedure TBModel_GSP_get_local_rep_E_virial
end interface get_local_rep_E_virial

interface calc_H_coeff
  module procedure TBModel_GSP_calc_H_coeff
end interface calc_H_coeff

interface calc_H_coeff_deriv
  module procedure TBModel_GSP_calc_H_coeff_deriv
end interface calc_H_coeff_deriv

contains

subroutine TBModel_GSP_Initialise_str(this, args_str, param_str)
  type(TBModel_GSP), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='TBModel_GSP_Initialise_str args_str')) then
  call system_abort("TBModel_GSP_Initialise_str parse problems"//trim(args_str))
  endif
  call finalise(params)

  call TBModel_GSP_read_params_xml(this, param_str)
end subroutine TBModel_GSP_Initialise_str

subroutine TBModel_GSP_Finalise(this)
  type(TBModel_GSP), intent(inout) :: this

  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%n_orbs)) deallocate(this%n_orbs)
  if (allocated(this%n_elecs)) deallocate(this%n_elecs)
  if (allocated(this%n_orb_sets)) deallocate(this%n_orb_sets)
  if (allocated(this%orb_set_type)) deallocate(this%orb_set_type)
  if (allocated(this%atomic_num)) deallocate(this%atomic_num)

  if (allocated(this%H_coeff)) deallocate(this%H_coeff)
  if (allocated(this%O_coeff)) deallocate(this%O_coeff)
  if (allocated(this%Atau)) deallocate(this%Atau)
  if (allocated(this%deltatau0)) deallocate(this%deltatau0)

  if (allocated(this%r0)) deallocate(this%r0)
  if (allocated(this%rc)) deallocate(this%rc)
  if (allocated(this%na)) deallocate(this%na)
  if (allocated(this%nb)) deallocate(this%nb)
  if (allocated(this%nc)) deallocate(this%nc)

  if (allocated(this%B)) deallocate(this%B)
  if (allocated(this%Rcore)) deallocate(this%Rcore)
  if (allocated(this%lambda0)) deallocate(this%lambda0)
  if (allocated(this%C)) deallocate(this%C)
  if (allocated(this%nu)) deallocate(this%nu)
  if (allocated(this%m)) deallocate(this%m)
  if (allocated(this%Rtail)) deallocate(this%Rtail)
  if (allocated(this%Rcut))  deallocate(this%Rcut)

  if (allocated(this%Ai)) deallocate(this%Ai)
  if (allocated(this%Ri)) deallocate(this%Ri)

  if (allocated(this%E)) deallocate(this%E)

  if (allocated(this%H_tail_spline)) deallocate(this%H_tail_spline)
  if (allocated(this%Vrep_env_emb_tail_spline)) deallocate(this%Vrep_env_emb_tail_spline)

  this%n_types = 0
  this%label = ''
end subroutine TBModel_GSP_Finalise

subroutine TBM_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  character(len=1024) :: value
  integer status
  integer ti, tj, Zi, Zj

  if (name == 'GSP_params') then ! new GSP stanza

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
      call QUIP_FoX_get_value(attributes, "cutoff_H", value, status);
      if (status == 0) read (value, *) parse_tbm%cutoff_H

      allocate(parse_tbm%atomic_num(parse_tbm%n_types))
      parse_tbm%atomic_num = 0
      allocate(parse_tbm%n_orbs(parse_tbm%n_types))
      allocate(parse_tbm%n_elecs(parse_tbm%n_types))
      allocate(parse_tbm%n_orb_sets(parse_tbm%n_types))
      allocate(parse_tbm%orb_set_type(max_n_orb_sets,parse_tbm%n_types))
      allocate(parse_tbm%E(max_n_orb_sets,parse_tbm%n_types))
      allocate(parse_tbm%B(parse_tbm%n_types,parse_tbm%n_types))
      allocate(parse_tbm%Rcore(parse_tbm%n_types,parse_tbm%n_types))
      allocate(parse_tbm%lambda0(parse_tbm%n_types))
      allocate(parse_tbm%C(parse_tbm%n_types))
      allocate(parse_tbm%nu(parse_tbm%n_types))
      allocate(parse_tbm%m(parse_tbm%n_types))
      allocate(parse_tbm%Rtail(parse_tbm%n_types))
      allocate(parse_tbm%Rcut(parse_tbm%n_types))
      allocate(parse_tbm%Ai(4,parse_tbm%n_types,parse_tbm%n_types))
      allocate(parse_tbm%Ri(4,parse_tbm%n_types,parse_tbm%n_types))

      allocate(parse_tbm%H_coeff(10, parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%O_coeff(10, parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%Atau(10, parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%deltatau0(10, parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%rc(parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%r0(parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%na(parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%nb(parse_tbm%n_types, parse_tbm%n_types))
      allocate(parse_tbm%nc(parse_tbm%n_types, parse_tbm%n_types))

      parse_tbm%H_coeff = 0.0_dp
      parse_tbm%O_coeff = 0.0_dp
      parse_tbm%Atau    = 0.0_dp
      parse_tbm%deltatau0 = 0.0_dp
      parse_tbm%rc     = 1.0_dp
      parse_tbm%r0     = 0.0_dp
      parse_tbm%na     = 0.0_dp
      parse_tbm%nb     = 1.0_dp
      parse_tbm%nc     = 0.0_dp
      parse_tbm%B      = 1.0_dp
      parse_tbm%Rcore  = 0.0_dp
      parse_tbm%lambda0 = 0.0_dp
      parse_tbm%C      = 0.0_dp
      parse_tbm%nu     = 0.0_dp
      parse_tbm%m      = 0.0_dp
      parse_tbm%Rtail  = 0.0_dp
      parse_tbm%Rcut   = 0.0_dp
      parse_tbm%Ai     = 0.0_dp
      parse_tbm%Ri     = 0.0_dp
    endif

  elseif (parse_in_tbm .and. name == 'per_type_data') then
    if (parse_cur_type > parse_tbm%n_types) &
      call system_abort('Too many types defined in GSP_params')

    call QUIP_FoX_get_value(attributes, "Z", value, status);
    if (status == 0) read (value, *) parse_tbm%atomic_num(parse_cur_type)

    call QUIP_FoX_get_value(attributes, "n_orbs", value, status);
    if (status == 0) read (value, *) parse_tbm%n_orbs(parse_cur_type)
    if (parse_tbm%n_orbs(parse_cur_type) == 5) then
      parse_tbm%n_orb_sets(parse_cur_type) = 1
      parse_tbm%orb_set_type(1,parse_cur_type) = ORB_D
      call QUIP_FoX_get_value(attributes, "E_d", value, status);
      if (status == 0) read (value, *) parse_tbm%E(ORB_D,parse_cur_type)
    else
     call system_abort("TBModel_GSP_read_params_xml can only do d")
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

!         Parameters of the repulsive environment-dependent term
    call QUIP_FoX_get_value(attributes, "lambda0", value, status);
    if (status == 0) read (value, *) parse_tbm%lambda0(parse_cur_type)
    call QUIP_FoX_get_value(attributes, "C", value, status);
    if (status == 0) read (value, *) parse_tbm%C(parse_cur_type)
    call QUIP_FoX_get_value(attributes, "nu", value, status);
    if (status == 0) read (value, *) parse_tbm%nu(parse_cur_type)
    call QUIP_FoX_get_value(attributes, "m", value, status);
    if (status == 0) read (value, *) parse_tbm%m(parse_cur_type)
    call QUIP_FoX_get_value(attributes, "Rtail", value, status);
    if (status == 0) read (value, *) parse_tbm%Rtail(parse_cur_type)
    call QUIP_FoX_get_value(attributes, "Rcut", value, status);
    if (status == 0) read (value, *) parse_tbm%Rcut(parse_cur_type)

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
    call QUIP_FoX_get_value(attributes, "na", value, status);
    if (status == 0) read (value, *) parse_tbm%na(ti,tj)
    call QUIP_FoX_get_value(attributes, "nb", value, status);
    if (status == 0) read (value, *) parse_tbm%nb(ti,tj)
    call QUIP_FoX_get_value(attributes, "nc", value, status);
    if (status == 0) read (value, *) parse_tbm%nc(ti,tj)

!         parameters of the repulsive environment-dependent term
    call QUIP_FoX_get_value(attributes, "B", value, status);
    if (status == 0) read (value, *) parse_tbm%B(ti,tj)
    call QUIP_FoX_get_value(attributes, "Rcore", value, status);
    if (status == 0) read (value, *) parse_tbm%Rcore(ti,tj)

    call QUIP_FoX_get_value(attributes, "A1", value, status);
    if (status == 0) read (value, *) parse_tbm%Ai(1,ti,tj)
    call QUIP_FoX_get_value(attributes, "A2", value, status);
    if (status == 0) read (value, *) parse_tbm%Ai(2,ti,tj)
    call QUIP_FoX_get_value(attributes, "A3", value, status);
    if (status == 0) read (value, *) parse_tbm%Ai(3,ti,tj)
    call QUIP_FoX_get_value(attributes, "A4", value, status);
    if (status == 0) read (value, *) parse_tbm%Ai(4,ti,tj)
    call QUIP_FoX_get_value(attributes, "R1", value, status);
    if (status == 0) read (value, *) parse_tbm%Ri(1,ti,tj)
    call QUIP_FoX_get_value(attributes, "R2", value, status);
    if (status == 0) read (value, *) parse_tbm%Ri(2,ti,tj)
    call QUIP_FoX_get_value(attributes, "R3", value, status);
    if (status == 0) read (value, *) parse_tbm%Ri(3,ti,tj)
    call QUIP_FoX_get_value(attributes, "R4", value, status);
    if (status == 0) read (value, *) parse_tbm%Ri(4,ti,tj)

    call QUIP_FoX_get_value(attributes, "H_dds", value, status);
    if (status == 0) read (value, *) parse_tbm%H_coeff(SK_DDS,ti,tj)
    call QUIP_FoX_get_value(attributes, "H_ddp", value, status);
    if (status == 0) read (value, *) parse_tbm%H_coeff(SK_DDP,ti,tj)
    call QUIP_FoX_get_value(attributes, "H_ddd", value, status);
    if (status == 0) read (value, *) parse_tbm%H_coeff(SK_DDD,ti,tj)

    call QUIP_FoX_get_value(attributes, "O_dds", value, status);
    if (status == 0) read (value, *) parse_tbm%O_coeff(SK_DDS,ti,tj)
    call QUIP_FoX_get_value(attributes, "O_ddp", value, status);
    if (status == 0) read (value, *) parse_tbm%O_coeff(SK_DDP,ti,tj)
    call QUIP_FoX_get_value(attributes, "O_ddd", value, status);
    if (status == 0) read (value, *) parse_tbm%O_coeff(SK_DDD,ti,tj)

    call QUIP_FoX_get_value(attributes, "Atau_dds", value, status);
    if (status == 0) read (value, *) parse_tbm%Atau(SK_DDS,ti,tj)
    call QUIP_FoX_get_value(attributes, "Atau_ddp", value, status);
    if (status == 0) read (value, *) parse_tbm%Atau(SK_DDP,ti,tj)
    call QUIP_FoX_get_value(attributes, "Atau_ddd", value, status);
    if (status == 0) read (value, *) parse_tbm%Atau(SK_DDD,ti,tj)

    call QUIP_FoX_get_value(attributes, "deltatau0_dds", value, status);
    if (status == 0) read (value, *) parse_tbm%deltatau0(SK_DDS,ti,tj)
    call QUIP_FoX_get_value(attributes, "deltatau0_ddp", value, status);
    if (status == 0) read (value, *) parse_tbm%deltatau0(SK_DDP,ti,tj)
    call QUIP_FoX_get_value(attributes, "deltatau0_ddd", value, status);
    if (status == 0) read (value, *) parse_tbm%deltatau0(SK_DDD,ti,tj)

    if (ti /= tj) then
      parse_tbm%rc(tj,ti) = parse_tbm%rc(ti,tj)
      parse_tbm%r0(tj,ti) = parse_tbm%r0(ti,tj)
      parse_tbm%na(tj,ti) = parse_tbm%na(ti,tj)
      parse_tbm%nb(tj,ti) = parse_tbm%nb(ti,tj)
      parse_tbm%nc(tj,ti) = parse_tbm%nc(ti,tj)
      parse_tbm%B(tj,ti) = parse_tbm%B(ti,tj)
      parse_tbm%Rcore(tj,ti) = parse_tbm%Rcore(ti,tj)
      parse_tbm%Ai(1,tj,ti) = parse_tbm%Ai(1,ti,tj)
      parse_tbm%Ai(2,tj,ti) = parse_tbm%Ai(2,ti,tj)
      parse_tbm%Ai(3,tj,ti) = parse_tbm%Ai(3,ti,tj)
      parse_tbm%Ai(4,tj,ti) = parse_tbm%Ai(4,ti,tj)
      parse_tbm%Ri(1,tj,ti) = parse_tbm%Ri(1,ti,tj)
      parse_tbm%Ri(2,tj,ti) = parse_tbm%Ri(2,ti,tj)
      parse_tbm%Ri(3,tj,ti) = parse_tbm%Ri(3,ti,tj)
      parse_tbm%Ri(4,tj,ti) = parse_tbm%Ri(4,ti,tj)
      parse_tbm%H_coeff(SK_DDS,tj,ti) = parse_tbm%H_coeff(SK_DDS,ti,tj)
      parse_tbm%H_coeff(SK_DDP,tj,ti) = parse_tbm%H_coeff(SK_DDP,ti,tj)
      parse_tbm%H_coeff(SK_DDD,tj,ti) = parse_tbm%H_coeff(SK_DDD,ti,tj)
      parse_tbm%O_coeff(SK_DDS,tj,ti) = parse_tbm%O_coeff(SK_DDS,ti,tj)
      parse_tbm%O_coeff(SK_DDP,tj,ti) = parse_tbm%O_coeff(SK_DDP,ti,tj)
      parse_tbm%O_coeff(SK_DDD,tj,ti) = parse_tbm%O_coeff(SK_DDD,ti,tj)
      parse_tbm%Atau(SK_DDS,tj,ti)    = parse_tbm%Atau(SK_DDS,ti,tj)
      parse_tbm%Atau(SK_DDP,tj,ti)    = parse_tbm%Atau(SK_DDP,ti,tj)
      parse_tbm%Atau(SK_DDD,tj,ti)    = parse_tbm%Atau(SK_DDD,ti,tj)
      parse_tbm%deltatau0(SK_DDS,tj,ti) = parse_tbm%deltatau0(SK_DDS,ti,tj)
      parse_tbm%deltatau0(SK_DDP,tj,ti) = parse_tbm%deltatau0(SK_DDP,ti,tj)
      parse_tbm%deltatau0(SK_DDD,tj,ti) = parse_tbm%deltatau0(SK_DDD,ti,tj)
    end if
  endif

end subroutine TBM_startElement_handler

subroutine TBM_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (name == 'GSP_params') then
    parse_in_tbm = .false.
  endif

end subroutine TBM_endElement_handler

subroutine TBModel_GSP_read_params_xml(this, param_str)
  type(TBModel_GSP), intent(inout), target :: this
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

  if (this%n_types <= 0) call system_abort("TBModel_GSP_read_params_xml couldn't find any types defined")
!"      
  call TBModel_GSP_fix_tails(this)

end subroutine

subroutine TBModel_GSP_fix_tails(this)
  type(TBModel_GSP), intent(inout) :: this

  integer ti, tj
  real(dp) :: x(2), y(2), yd1, yd2, x_env(2)

  x(1) = this%tailx0
  x(2) = this%cutoff_H
  y(2) = 0.0_dp
  yd2 = 0.0_dp

  allocate(this%H_tail_spline(this%n_types, this%n_types))
  allocate(this%Vrep_env_emb_tail_spline(this%n_types))

  do ti=1, this%n_types
    do tj=1, this%n_types
       y(1) = TBModel_GSP_dist_scaling(x(1),this%r0(ti,tj), &
         this%rc(ti,tj),this%na(ti,tj),this%nb(ti,tj),this%nc(ti,tj))
       yd1 = TBModel_GSP_dist_scaling_deriv(x(1),this%r0(ti,tj), &
         this%rc(ti,tj),this%na(ti,tj),this%nb(ti,tj),this%nc(ti,tj))
       call spline_init(this%H_tail_spline(ti,tj), x, y, yd1, yd2)
    end do

    x_env(1) = this%Rtail(ti)
    x_env(2) = this%Rcut(ti)
    y(1) = TBModel_GSP_Vrep_env_emb_term(x_env(1),this%C(ti),this%nu(ti)) 
    yd1  = TBModel_GSP_Vrep_env_emb_term_deriv(x_env(1),this%C(ti),this%nu(ti)) 
    call spline_init(this%Vrep_env_emb_tail_spline(ti), x_env, y, yd1, yd2)
  end do

end subroutine

subroutine TBModel_GSP_Print(this,file)
  type(TBModel_GSP),    intent(in)           :: this
  type(Inoutput), intent(inout),optional      :: file

  integer:: ti, tj

  if(this%n_types == 0) call System_Abort('TBModel_GSP_Print: TBModel_GSP structure not initialised')
!'      
  call Print('TBModel_GSP Structure:', file=file)

  !! print GSP params
  call Print ("TBModel_GSP tailx0 " // this%tailx0 // " cutoff " // this%cutoff, file=file)
  call Print ("                    " // " cutoff_H " // this%cutoff_H, file=file)
  do ti=1, this%n_types
    call Print ("TBModel_GSP type " // ti // " Z " // this%atomic_num(ti) //  " n_orbs " // this%n_orbs(ti) // &
      " n_elecs " // this%n_elecs(ti), file=file)
    if (this%n_orb_sets(ti) == 1) then
      call Print ("TBModel_GSP E " // this%E(1,ti), file=file)
    else
      call Print ("TBModel_GSP Es " // this%E(ORB_S,ti), file=file)
      call Print ("TBModel_GSP Ed " // this%E(ORB_D,ti), file=file)
    endif
  end do

  call verbosity_push_decrement()
  do ti=1, this%n_types
    do tj=1, this%n_types
     call Print ("TBModel_GSP interaction " // &
      ti // " " //  tj // " Z " //  this%atomic_num(ti) // " " // this%atomic_num(tj), file=file)

     call Print ("TBModel_GSP r0 " // this%r0(ti,tj), file=file)

     call Print ("TBModel_GSP SK dds " // this%H_coeff(SK_DDS,ti,tj), file=file)
     call Print ("TBModel_GSP SK ddp " // this%H_coeff(SK_DDP,ti,tj), file=file)
     call Print ("TBModel_GSP SK ddd " // this%H_coeff(SK_DDD,ti,tj), file=file)
     call Print ("TBModel_GSP Overlap coefficients dds " // this%O_coeff(SK_DDS,ti,tj), file=file)
     call Print ("TBModel_GSP Overlap coefficients ddp " // this%O_coeff(SK_DDP,ti,tj), file=file)
     call Print ("TBModel_GSP Overlap coefficients ddd " // this%O_coeff(SK_DDD,ti,tj), file=file)
     call Print ("TBModel_GSP Atau dds " // this%Atau(SK_DDS,ti,tj), file=file)
     call Print ("TBModel_GSP Atau ddp " // this%Atau(SK_DDP,ti,tj), file=file)
     call Print ("TBModel_GSP Atau ddd " // this%Atau(SK_DDD,ti,tj), file=file)
     call Print ("TBModel_GSP deltatau0 dds " // this%deltatau0(SK_DDS,ti,tj), file=file)
     call Print ("TBModel_GSP deltatau0 ddp " // this%deltatau0(SK_DDP,ti,tj), file=file)
     call Print ("TBModel_GSP deltatau0 ddd " // this%deltatau0(SK_DDD,ti,tj), file=file)
     call Print ("TBModel_GSP H scaling rc na nb nc " // &
       this%rc(ti,tj) // " " //  this%na(ti,tj) // " " //  this%nb(ti,tj) //  " " // this%nc(ti,tj), file=file)

     call Print (this%H_tail_spline(ti,tj), file=file)
     call Print (this%Vrep_env_emb_tail_spline(ti), file=file)

    call Print ("TBModel_GSP pair potential Vrep A1 A2 A3 A4 " // &
       this%Ai(1,ti,tj) // " " //  this%Ai(2,ti,tj) // " " // &
        this%Ai(3,ti,tj) // " " //  this%Ai(4,ti,tj), file=file)
    call Print ("TBModel_GSP pair potential Vrep R1 R2 R3 R4 " // &
       this%Ri(1,ti,tj) // " " //  this%Ri(2,ti,tj) // " " // &
        this%Ri(3,ti,tj) // " " //  this%Ri(4,ti,tj), file=file)

     call Print ("TBModel_GSP Environmental Vrep : B Rcore " // &
           this%B(ti,tj) // " " //  this%Rcore(ti,tj), file=file)
    end do
      call Print ("TBModel_GSP Environmental Vrep : lambda0 nu " // &
       this%lambda0(ti) // " " //  this%nu(ti) , file=file)
      call Print ("TBModel_GSP Environmental Vrep : C m " // &
       this%c(ti) // " " // this%m(ti), file=file)
      call Print ("TBModel_GSP Environmental Vrep : Rtail Rcut " // &
       this%Rtail(ti) // " " // this%Rcut(ti), file=file)
  end do
  call verbosity_pop()

end subroutine TBModel_GSP_Print

function TBModel_GSP_n_orbs_of_Z(this, Z)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in) :: Z
  integer TBModel_GSP_n_orbs_of_Z

  TBModel_GSP_n_orbs_of_Z = this%n_orbs(get_type(this%type_of_atomic_num,Z))
end function TBModel_GSP_n_orbs_of_Z

function TBModel_GSP_n_orb_sets_of_Z(this, Z)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in) :: Z
  integer TBModel_GSP_n_orb_sets_of_Z

  TBModel_GSP_n_orb_sets_of_Z = this%n_orb_sets(get_type(this%type_of_atomic_num,Z))
end function TBModel_GSP_n_orb_sets_of_Z

function TBModel_GSP_n_orbs_of_orb_set_of_Z(this, Z, i_set)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in) :: Z, i_set
  integer TBModel_GSP_n_orbs_of_orb_set_of_Z

  TBModel_GSP_n_orbs_of_orb_set_of_Z = N_ORBS_OF_SET(this%orb_set_type(i_set,get_type(this%type_of_atomic_num,Z)))

end function TBModel_GSP_n_orbs_of_orb_set_of_Z

function TBModel_GSP_orb_type_of_orb_set_of_Z(this, Z, i_set)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in) :: Z, i_set
  integer TBModel_GSP_orb_type_of_orb_set_of_Z

  TBModel_GSP_orb_type_of_orb_set_of_Z = this%orb_set_type(i_set,get_type(this%type_of_atomic_num,Z))

end function TBModel_GSP_orb_type_of_orb_set_of_Z

function TBModel_GSP_n_elecs_of_Z(this, Z)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in) :: Z
  real(dp) TBModel_GSP_n_elecs_of_Z

  TBModel_GSP_n_elecs_of_Z = this%n_elecs(get_type(this%type_of_atomic_num,Z))
end function TBModel_GSP_n_elecs_of_Z

subroutine TBModel_GSP_get_HS_blocks(this, at, at_i, at_j, dv_hat, dv_mag, b_H, b_S)
  type(TBModel_GSP), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_i, at_j
  real(dp), intent(in) :: dv_hat(3), dv_mag
  real(dp), intent(out) :: b_H(:,:)
  real(dp), intent(out) :: b_S(:,:)

  integer i, j, ti, tj, is, js, i_set, j_set
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
         SK_frad_H)
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

end subroutine TBModel_GSP_get_HS_blocks

subroutine TBModel_GSP_get_dHS_masks(this, at_ind, d_mask, od_mask)
  type(TBModel_GSP), intent(in) :: this
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

function TBModel_GSP_get_dHS_blocks(this, at, at_i, at_j, dv_hat, dv_mag, at_ind, b_dH, b_dS)
  type(TBModel_GSP), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_i, at_j
  real(dp), intent(in) :: dv_hat(3), dv_mag
  integer, intent(in) :: at_ind
  real(dp), intent(out) :: b_dH(:,:,:)
  real(dp), intent(out) :: b_dS(:,:,:)
  logical :: TBModel_GSP_get_dHS_blocks

  integer :: i, j, ti, tj, is, js, i_set, j_set
  real(dp) SK_frad_H(N_SK) 
  real(dp) SK_dfrad_H(N_SK)
  real(dp) dv_hat_sq(3)
  real(dp) virial_outerprod_fac

  if ((at_ind > 0 .and. at_i /= at_ind .and. at_j /= at_ind) .or. &
      (at_ind > 0 .and. at_i == at_j) .or. &
      (dv_mag > this%cutoff_H) .or. &
      (dv_mag .feq. 0.0_dp)) then
    TBModel_GSP_get_dHS_blocks = .false.
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
      SK_frad_H)
    call dradial_functions(this, ti, tj, dv_mag, this%orb_set_type(i_set,ti), this%orb_set_type(j_set, tj), &
      SK_dfrad_H)
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

  TBModel_GSP_get_dHS_blocks = .true.
  return

end function TBModel_GSP_get_dHS_blocks

subroutine radial_functions(this, ti, tj, dv_mag, orb_set_type_i, orb_set_type_j, f_H)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: dv_mag
  integer, intent(in) :: orb_set_type_i, orb_set_type_j
  real(dp), intent(out) :: f_H(N_SK)

  if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_D) then
    f_H(SK_DDS) = calc_H_coeff(this, SK_DDS, dv_mag, ti, tj)
    f_H(SK_DDP) = calc_H_coeff(this, SK_DDP, dv_mag, ti, tj)
    f_H(SK_DDD) = calc_H_coeff(this, SK_DDD, dv_mag, ti, tj)
  endif
end subroutine radial_functions

subroutine dradial_functions(this, ti, tj, dv_mag, orb_set_type_i, orb_set_type_j, f_dH)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: dv_mag
  integer, intent(in) :: orb_set_type_i, orb_set_type_j
  real(dp), intent(out) :: f_dH(N_SK)

  if (orb_set_type_i == ORB_D .and. orb_set_type_j == ORB_D) then
    f_dH(SK_DDS) = calc_H_coeff_deriv(this, SK_DDS, dv_mag, ti, tj)
    f_dH(SK_DDP) = calc_H_coeff_deriv(this, SK_DDP, dv_mag, ti, tj)
    f_dH(SK_DDD) = calc_H_coeff_deriv(this, SK_DDD, dv_mag, ti, tj)
  endif
end subroutine dradial_functions

function TBModel_GSP_calc_H_coeff(this, sk_ind, dist, ti, tj)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in) :: sk_ind
  real(dp), intent(in) :: dist
  integer, intent(in) :: ti, tj
  real(dp) :: TBModel_GSP_calc_H_coeff

  TBModel_GSP_calc_H_coeff = this%H_coeff(sk_ind, ti, tj)* &
    TBModel_GSP_H_dist_func(this, dist, ti, tj)

end function TBModel_GSP_calc_H_coeff

function TBModel_GSP_calc_H_coeff_deriv(this, sk_ind, dist, ti, tj)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in) :: sk_ind
  real(dp), intent(in) :: dist
  integer, intent(in) :: ti, tj
  real(dp) :: TBModel_GSP_calc_H_coeff_deriv

  TBModel_GSP_calc_H_coeff_deriv = this%H_coeff(sk_ind, ti, tj)* &
    TBModel_GSP_H_dist_func_deriv(this, dist, ti, tj)

end function TBModel_GSP_calc_H_coeff_deriv

function TBModel_GSP_calc_O_coeff(this, sk_ind, dist, ti, tj)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in) :: sk_ind
  real(dp), intent(in) :: dist
  integer, intent(in) :: ti, tj
  real(dp) :: TBModel_GSP_calc_O_coeff

  TBModel_GSP_calc_O_coeff = this%O_coeff(sk_ind, ti, tj)* &
    TBModel_GSP_H_dist_func(this, dist, ti, tj)

end function TBModel_GSP_calc_O_coeff

function TBModel_GSP_calc_O_coeff_deriv(this, sk_ind, dist, ti, tj)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in) :: sk_ind
  real(dp), intent(in) :: dist
  integer, intent(in) :: ti, tj
  real(dp) :: TBModel_GSP_calc_O_coeff_deriv

  TBModel_GSP_calc_O_coeff_deriv = this%O_coeff(sk_ind, ti, tj)* &
    TBModel_GSP_H_dist_func_deriv(this, dist, ti, tj)

end function TBModel_GSP_calc_O_coeff_deriv

function TBModel_GSP_H_dist_func(this, dist, ti, tj)
  type(TBModel_GSP), intent(in) :: this
  real(dp), intent(in) :: dist
  integer, intent(in) :: ti, tj
  real(dp) :: TBModel_GSP_H_dist_func

  if (dist <= this%cutoff_H) then
    if (dist <= this%tailx0) then
      TBModel_GSP_H_dist_func = TBModel_GSP_dist_scaling(dist, &
        this%r0(ti,tj),this%rc(ti,tj),this%na(ti,tj),this%nb(ti,tj),this%nc(ti,tj)) 
    else
      TBModel_GSP_H_dist_func = spline_value(this%H_tail_spline(ti,tj), dist)
    endif
  else
    TBModel_GSP_H_dist_func = 0.0_dp
  endif
end function TBModel_GSP_H_dist_func

function TBModel_GSP_H_dist_func_deriv(this, dist, ti, tj)
  type(TBModel_GSP), intent(in) :: this
  real(dp), intent(in) :: dist
  integer, intent(in) :: ti, tj
  real(dp) :: TBModel_GSP_H_dist_func_deriv

  if (dist <= this%cutoff_H) then
    if (dist <= this%tailx0) then
      TBModel_GSP_H_dist_func_deriv = TBModel_GSP_dist_scaling_deriv(dist, &
        this%r0(ti,tj),this%rc(ti,tj),this%na(ti,tj),this%nb(ti,tj),this%nc(ti,tj))
    else
      TBModel_GSP_H_dist_func_deriv = spline_deriv(this%H_tail_spline(ti,tj), dist)
    endif
  else
    TBModel_GSP_H_dist_func_deriv = 0.0_dp
  endif
end function TBModel_GSP_H_dist_func_deriv

! Compute the screening function Sij_tau
function TBModel_GSP_screening(this,sk_ind,at,i,j,dist,ti,tj)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in)            :: sk_ind
  type(Atoms)                    :: at
  real(dp), intent(in)           :: dist
  integer, intent(in)            :: ti, tj
  integer                        :: i,j
  real(dp)                       :: TBModel_GSP_screening

  TBModel_GSP_screening =  &
       (TBModel_GSP_screening_c(this,sk_ind,at,i,j,dist,ti,tj) &
              - TBModel_GSP_screening_mu(this) )/ &
     (1._dp + TBModel_GSP_calc_O_coeff(this,sk_ind,dist,ti,tj)**2.0_dp - &
            2._dp * TBModel_GSP_screening_mu(this) )

end function TBModel_GSP_screening

function TBModel_GSP_screening_c(this,sk_ind,at, i,j,dist,ti,tj) 
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in)            :: sk_ind
  type(Atoms)                    :: at
  integer, intent(in)            :: ti, tj
  real(dp), intent(in)           :: dist
  integer                        :: i,j
  real(dp)                       :: TBModel_GSP_screening_c
  real(dp)                       :: c_temp
  real(dp)                       :: rik, rjk
  real(dp)                       :: g_i, g_j, g_j2, theta_i, theta_j
  integer                        :: k, kk, ki, kj, tk
  logical                        :: lgo
  real(dp)                       :: cterm(15)

  c_temp = 0._dp
  do ki = 1, atoms_n_neighbours(at, i)
     k = atoms_neighbour(at, i, ki, rik)
     if (rik .feq. 0.0_dp) cycle
     tk = get_type(this%type_of_atomic_num, at%Z(k))
     rjk = distance(at, j, k, (/0, 0, 0/))  
   
     theta_i = cosine_neighbour(at,i,j,k)
     theta_j = cosine_neighbour(at,j,i,k)
     g_i = TBModel_GSP_screening_g(sk_ind,theta_i)
     g_j = TBModel_GSP_screening_g(sk_ind,theta_j)
     g_j2 = TBModel_GSP_screening_g(sk_ind,-theta_j)

     call TBModel_GSP_calc_c_terms(this)
     !call TBModel_GSP_calc_c_terms(this,sk_ind,dist,rik,rjk,ti,tj,tk,g_i,g_j,g_j2,cterm)

     c_temp = c_temp + (cterm(1)*cterm(2) + cterm(3)*cterm(4))*cterm(5)*cterm(6) - &
                        cterm(7) *cterm(8) *cterm(9) *cterm(10)*cterm(10) - &
                        cterm(11)*cterm(12)*cterm(13)*cterm(14)*cterm(14) 
  enddo

  do kj =1, atoms_n_neighbours(at, j)
     k =  atoms_neighbour(at, j, kj, rjk)
     if (rjk .feq. 0.0_dp) cycle
     tk = get_type(this%type_of_atomic_num, at%Z(k))
     lgo = .true.
     do  ki=1, atoms_n_neighbours(at, i)
       kk = atoms_neighbour(at, i, ki)  
       if(k.eq.kk) then
          exit
          lgo = .false.
       endif
     enddo
     if(.not.lgo) cycle 
     rik = distance(at, i, k, (/0, 0, 0/))

     theta_i = cosine_neighbour(at,i,j,k)
     theta_j = cosine_neighbour(at,j,i,k)
     g_i = TBModel_GSP_screening_g(sk_ind,theta_i)
     g_j = TBModel_GSP_screening_g(sk_ind,theta_j)
     g_j2 = TBModel_GSP_screening_g(sk_ind,-theta_j)

     !call TBModel_GSP_calc_c_terms(this,sk_ind,dist,rik,rjk,ti,tj,tk,g_i,g_j,g_j2,cterm)
     call TBModel_GSP_calc_c_terms(this)

     c_temp = c_temp + (cterm(1)*cterm(2) + cterm(3)*cterm(4))*cterm(5)*cterm(6) - &
                        cterm(7) *cterm(8) *cterm(9) *cterm(10)*cterm(10) - &
                        cterm(11)*cterm(12)*cterm(13)*cterm(14)*cterm(14) 
  enddo

  TBModel_GSP_screening_c = &
  this%Atau(sk_ind,ti,tj)*(1._dp + this%deltatau0(sk_ind,ti,tj))/4._dp/ &
               TBModel_GSP_calc_H_coeff (this,sk_ind,dist,ti,tj) * c_temp

end function TBModel_GSP_screening_c

function TBModel_GSP_screening_c_deriv(this)
  type(TBModel_GSP), intent(in) :: this
  real(dp)                       :: TBModel_GSP_screening_c_deriv
    
   TBModel_GSP_screening_c_deriv = 0.d0

end function TBModel_GSP_screening_c_deriv

function TBModel_GSP_screening_mu(this) 
  type(TBModel_GSP), intent(in) :: this
  real(dp)                       :: TBModel_GSP_screening_mu

  TBModel_GSP_screening_mu = 0._dp

end function TBModel_GSP_screening_mu

function TBModel_GSP_screening_g(sk_ind,theta)
  integer, intent(in)            :: sk_ind
  real(dp)                       :: theta
  real(dp)                       ::  TBModel_GSP_screening_g

  if(sk_ind.eq.SK_DDS) then  
     TBModel_GSP_screening_g = (1._dp/4._dp) * (1._dp + 3._dp* dcos(2._dp*theta))      
  elseif(sk_ind.eq.SK_DDP) then  
     TBModel_GSP_screening_g = (dsqrt(3._dp)/2._dp)*dsin(2._dp*theta)
  elseif(sk_ind.eq.SK_DDD) then  
     TBModel_GSP_screening_g = (dsqrt(3._dp)/4._dp)*(1._dp - dcos(2._dp*theta))
  endif

end function TBModel_GSP_screening_g

subroutine TBModel_GSP_calc_c_terms(this)  
  type(TBModel_GSP), intent(in) :: this
! integer, intent(in)            :: sk_ind
! integer, intent(in)            :: ti, tj, tk
! real(dp), intent(in)           :: rij, rik, rjk
! real(dp)                       :: g_i, g_j, g_j2
! real(dp)                       :: cterm(15)

!     cterm(1) = TBModel_GSP_calc_H_coeff(this,SK_SDS,rik,ti,tk)
!     cterm(2) = TBModel_GSP_calc_O_coeff(this,SK_SDS,rjk,tk,tj)
!     cterm(3) = TBModel_GSP_calc_O_coeff(this,SK_SDS,rik,ti,tk)
!     cterm(4) = TBModel_GSP_calc_H_coeff(this,SK_SDS,rjk,tk,tj)
!     cterm(5) = g_i 
!     cterm(6) = g_j2
!     cterm(7) = TBModel_GSP_calc_H_coeff(this,SK_SDS,rik,ti,tk)
!     cterm(8) = TBModel_GSP_calc_O_coeff(this,SK_SDS,rik,tk,ti)
!     cterm(9) = TBModel_GSP_calc_O_coeff(this,sk_ind,rij,ti,tj)
!     cterm(10)= g_i
!     cterm(11)= TBModel_GSP_calc_O_coeff(this,sk_ind,rij,ti,tj)
!     cterm(12)= TBModel_GSP_calc_O_coeff(this,SK_SDS,rjk,tj,tk)
!     cterm(13)= TBModel_GSP_calc_H_coeff(this,SK_SDS,rjk,tk,tj)
!     cterm(14)= g_j 
end subroutine TBModel_GSP_calc_c_terms

subroutine TBModel_GSP_calc_c_terms_deriv(this)  
  type(TBModel_GSP), intent(in) :: this
!  integer, intent(in)            :: sk_ind
!  integer, intent(in)            :: ti, tj, tk
!  real(dp), intent(in)           :: rij, rik, rjk
!  real(dp)                       :: dg_i, dg_j, dg_j2
!  real(dp)                       :: cterm_deriv(15)

!       cterm_deriv(1) = TBModel_GSP_calc_H_coeff_deriv(this,SK_SDS,rik,ti,tk)
!       cterm_deriv(2) = TBModel_GSP_calc_O_coeff_deriv(this,SK_SDS,rjk,tk,tj)
!       cterm_deriv(3) = TBModel_GSP_calc_O_coeff_deriv(this,SK_SDS,rik,ti,tk)
!       cterm_deriv(4) = TBModel_GSP_calc_H_coeff_deriv(this,SK_SDS,rjk,tk,tj)
!       cterm_deriv(5) = dg_i
!       cterm_deriv(6) = dg_j2
!       cterm_deriv(7) = TBModel_GSP_calc_H_coeff_deriv(this,SK_SDS,rik,ti,tk)
!       cterm_deriv(8) = TBModel_GSP_calc_O_coeff_deriv(this,SK_SDS,rik,tk,ti)
!       cterm_deriv(9) = TBModel_GSP_calc_O_coeff_deriv(this,sk_ind,rij,ti,tj)
!       cterm_deriv(10)= dg_i
!       cterm_deriv(11)= TBModel_GSP_calc_O_coeff_deriv(this,sk_ind,rij,ti,tj)
!       cterm_deriv(12)= TBModel_GSP_calc_O_coeff_deriv(this,SK_SDS,rjk,tj,tk)
!       cterm_deriv(13)= TBModel_GSP_calc_H_coeff_deriv(this,SK_SDS,rjk,tk,tj)
!       cterm_deriv(14)= dg_j
end subroutine TBModel_GSP_calc_c_terms_deriv

!Compute the repulsive pair potential V(Rij) = sum_k A_k (R_k - Rij)**3
function TBModel_GSP_Vrep_dist_func(this, dist, ti, tj)
  type(TBModel_GSP), intent(in) :: this
  real(dp), intent(in)           :: dist
  integer, intent(in)            :: ti, tj
  integer                        :: i
  real(dp)                       :: TBModel_GSP_Vrep_dist_func

  TBModel_GSP_Vrep_dist_func = 0.0_dp
         
  do i = 1, 4
    if(dist.gt.this%Ri(i,ti,tj) ) cycle
    TBModel_GSP_Vrep_dist_func = TBModel_GSP_Vrep_dist_func + this%Ai(i,ti,tj) * (this%Ri(i,ti,tj) - dist)**3.0_dp
  enddo

end function

! Compute the derivative of the repulsive pair potential V(Rij) = sum_k A_k (R_k - Rij)**3
function TBModel_GSP_Vrep_dist_func_deriv(this, dist, ti, tj)
  type(TBModel_GSP), intent(in) :: this
  real(dp), intent(in)           :: dist
  integer, intent(in)            :: ti, tj
  integer                        :: i
  real(dp)                       :: TBModel_GSP_Vrep_dist_func_deriv

  TBModel_GSP_Vrep_dist_func_deriv = 0.0_dp
  do i = 1, 4
    if(dist.gt.this%Ri(i,ti,tj) ) cycle
    TBModel_GSP_Vrep_dist_func_deriv =   TBModel_GSP_Vrep_dist_func_deriv - 3.0_dp * this%Ai(i,ti,tj) * (this%Ri(i,ti,tj) - dist)**2.0_dp
  enddo
end function TBModel_GSP_Vrep_dist_func_deriv

! Compute the environmental Repulsive pair potential 
function TBModel_GSP_Vrep_env(this,dist,lambda,ti,tj)
  type(TBModel_GSP), intent(in) :: this
  real(dp), intent(in)           :: dist
  real(dp)                       :: lambda
  integer, intent(in)            :: ti, tj
  real(dp)                       :: TBModel_GSP_Vrep_env 

  TBModel_GSP_Vrep_env = this%B(ti,tj)/dist * &
              dexp(-lambda*(dist - 2.0_dp* this%Rcore(ti,tj)))

end function TBModel_GSP_Vrep_env

! Compute the derivative (term ij) of the environmental Repulsive pair potential 
function TBModel_GSP_Vrep_env_deriv_ij(this,dist,ti,tj,lambda)

  type(TBModel_GSP), intent(in) :: this
  real(dp), intent(in)           :: dist
  integer, intent(in)            :: ti, tj
  real(dp), intent(in)           :: lambda
  real(dp)                       :: TBModel_GSP_Vrep_env_deriv_ij

  TBModel_GSP_Vrep_env_deriv_ij = - this%B(ti,tj)/dist * &
              dexp(-lambda*(dist - 2.0_dp* this%Rcore(ti,tj))) * &
              (1.0_dp/dist +  lambda) 

end function TBModel_GSP_Vrep_env_deriv_ij

! Compute the derivative (term wk) of the environmental Repulsive pair potential 
function TBModel_GSP_Vrep_env_deriv_wk(this,dist,dist_k,ti,tj,emb_i,lambda)
  type(TBModel_GSP), intent(in) :: this
  real(dp), intent(in)           :: dist, dist_k
  integer, intent(in)            :: ti, tj
  real(dp)                       :: emb_i 
  real(dp)                       :: lambda
  real(dp)                       :: emb_term_deriv 
  real(dp)                       :: TBModel_GSP_Vrep_env_deriv_wk

  if(dist_k.lt.this%Rtail(ti)) then
    emb_term_deriv = TBModel_GSP_Vrep_env_emb_term_deriv(dist_k,this%C(ti),this%nu(ti)) 
  elseif(dist_k.lt.this%Rcut(ti)) then
    emb_term_deriv = spline_deriv(this%Vrep_env_emb_tail_spline(ti), dist_k)
  else
    emb_term_deriv = 0.0_dp
  endif

  TBModel_GSP_Vrep_env_deriv_wk = - this%B(ti,tj)/dist/this%m(ti)/2.0_dp * &
          (dist - 2.0_dp* this%Rcore(ti,tj)) * &
          dexp(-lambda * (dist - 2.0_dp* this%Rcore(ti,tj)))*&
          emb_i ** ((1-this%m(ti))/this%m(ti)) * &
          emb_term_deriv 

end function TBModel_GSP_Vrep_env_deriv_wk

! Compute the embedded-atom expression =>  Sum C exp(-nu Rik) 
function TBModel_GSP_Vrep_env_emb(this, at, i, ti)
  type(TBModel_GSP), intent(in) :: this
  type(Atoms), intent(in)        :: at 
  integer                        :: i, ti, k 
  integer                        :: ki 
  real(dp)                       :: emb, rik 
  real(dp)                       :: TBModel_GSP_Vrep_env_emb

  emb = 0.0_dp
  do ki=1, atoms_n_neighbours(at, i)
    k = atoms_neighbour(at, i, ki, rik)
    if (rik .feq. 0.0_dp) cycle
    if(rik.lt.this%Rtail(ti)) then  
      emb = emb + TBModel_GSP_Vrep_env_emb_term(rik,this%C(ti),this%nu(ti))  
    elseif(rik.lt.this%Rcut(ti)) then
      emb = emb + spline_value(this%Vrep_env_emb_tail_spline(ti),rik)
    endif 
  enddo
  TBModel_GSP_Vrep_env_emb = emb
end function TBModel_GSP_Vrep_env_emb

! compute C exp(-nu Rij)
function TBModel_GSP_Vrep_env_emb_term(dist,C, nu)
  real(dp), intent(in)           :: C, nu
  real(dp), intent(in)           :: dist 
  real(dp)                       :: TBModel_GSP_Vrep_env_emb_term

  TBModel_GSP_Vrep_env_emb_term = C * dexp(-nu*dist)

end function TBModel_GSP_Vrep_env_emb_term

! compute derivative of C exp(-nu Rij)
function TBModel_GSP_Vrep_env_emb_term_deriv(dist,C, nu)
  real(dp), intent(in)           :: C, nu
  real(dp), intent(in)           :: dist 
  real(dp)                       :: TBModel_GSP_Vrep_env_emb_term_deriv

  TBModel_GSP_Vrep_env_emb_term_deriv = - nu * C * dexp(-nu*dist)

end function TBModel_GSP_Vrep_env_emb_term_deriv

! Compute lambda, necessary for the embedded-atom like repulsive potential
function TBModel_GSP_Vrep_env_lambda(this,ti,tj,emb_i,emb_j)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in)            :: ti, tj
  real(dp)                       :: emb_i, emb_j 
  real(dp)                       :: lambda_i, lambda_j
  real(dp)                       :: TBModel_GSP_Vrep_env_lambda

  lambda_i = this%lambda0(ti) + emb_i**(1/this%m(ti))
  lambda_j = this%lambda0(tj) + emb_j**(1/this%m(tj))
  TBModel_GSP_Vrep_env_lambda = (lambda_i + lambda_j)/2.0_dp 

end function TBModel_GSP_Vrep_env_lambda

! Scaling factor for the hopping parameters 
function TBModel_GSP_dist_scaling(r, r0, rc, na, nb, nc)
  real(dp), intent(in) :: r, r0, rc, na, nb, nc
  real(dp) :: TBModel_GSP_dist_scaling

  TBModel_GSP_dist_scaling = (r0/r)**na * exp(nb* ((r0/rc)**nc - (r/rc)**nc))

end function TBModel_GSP_dist_scaling

! Scaling factor derivative for the hopping parameters 
function TBModel_GSP_dist_scaling_deriv(r, r0, rc, na, nb, nc)
  real(dp), intent(in) :: r, r0, rc, na, nb, nc
  real(dp) :: TBModel_GSP_dist_scaling_deriv

  real(dp) tmp1, tmp2

  tmp1 = (r0/r)**na
  tmp2 = exp(nb*((r0/rc)**nc - (r/rc)**nc))

  TBModel_GSP_dist_scaling_deriv = na*tmp1/(-r)*tmp2 + tmp1*tmp2*(-nb)*nc*((r/rc)**nc)/r 
end function TBModel_GSP_dist_scaling_deriv

function onsite_function(this, ti, orb_set_type)
  type(TBModel_GSP), intent(in) :: this
  integer, intent(in) :: ti, orb_set_type
  real(dp) :: onsite_function

  onsite_function = this%E(orb_set_type, ti)

end function onsite_function

!Repulsive term: Pure two body potential + Repulsive embedded-atom like potential
function TBModel_GSP_get_local_rep_E(this, at, i)
  type(TBModel_GSP), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: emb_i, emb_j, lambda_ij
  real(dp) :: TBModel_GSP_get_local_rep_E

  real(dp) :: E, dist
  integer ji, j, ti, tj

  TBModel_GSP_get_local_rep_E = 0.0_dp

  ti = get_type(this%type_of_atomic_num, at%Z(i))
  emb_i = TBModel_GSP_Vrep_env_emb(this, at, i, ti) 
  do ji=1, atoms_n_neighbours(at, i)
    j = atoms_neighbour(at, i, ji, dist)
    if (dist .feq. 0.0_dp) cycle
    tj = get_type(this%type_of_atomic_num, at%Z(j))
    emb_j = TBModel_GSP_Vrep_env_emb(this, at, j, tj)   
    lambda_ij = TBModel_GSP_Vrep_env_lambda(this,ti,tj,emb_i,emb_j)

    E = TBModel_GSP_Vrep_dist_func(this,dist,ti,tj) + &
        TBModel_GSP_Vrep_env(this,dist,lambda_ij,ti,tj)

    TBModel_GSP_get_local_rep_E = TBModel_GSP_get_local_rep_E + E/2.0_dp
  end do

end function TBModel_GSP_get_local_rep_E

! Forces due to the repulsive term: two body potential + embedded-atom like potential
function TBModel_GSP_get_local_rep_E_force(this, at, i) result(force)
  type(TBModel_GSP), intent(in) :: this
  type(Atoms), intent(in)        :: at
  integer, intent(in)            :: i
  real(dp)                       :: force(3,at%N)

  real(dp) :: dE_dr, dist, dist_k, dv_hat(3)
  real(dp) :: emb_i, emb_j 
  real(dp) :: lambda_ij
  integer  :: ji, j, ti, tj
  integer  :: k, ik, jk

  force = 0.0_dp

  ti = get_type(this%type_of_atomic_num, at%Z(i))
  emb_i = TBModel_GSP_Vrep_env_emb(this, at, i, ti)
  do ji=1, atoms_n_neighbours(at, i)
    j = atoms_neighbour(at, i, ji, dist, cosines = dv_hat)
    if (dist .feq. 0.0_dp) cycle
    tj = get_type(this%type_of_atomic_num, at%Z(j))
    emb_j = TBModel_GSP_Vrep_env_emb(this, at, j, tj)
    lambda_ij = TBModel_GSP_Vrep_env_lambda(this,ti,tj,emb_i,emb_j)

    dE_dr =  TBModel_GSP_Vrep_dist_func_deriv(this, dist, ti, tj) + &
             TBModel_GSP_Vrep_env_deriv_ij(this,dist,ti,tj,lambda_ij) 
    force(:,i) = force(:,i) + dE_dr*dv_hat(:)/2.0_dp
    force(:,j) = force(:,j) - dE_dr*dv_hat(:)/2.0_dp

    do ik = 1, atoms_n_neighbours(at, i)
     k = atoms_neighbour(at, i, ik, dist_k, cosines = dv_hat)
     if (dist_k .feq. 0.0_dp) cycle
     dE_dr = TBModel_GSP_Vrep_env_deriv_wk(this,dist,dist_k,ti,tj,emb_i,lambda_ij)
     force(:,i) = force(:,i) + dE_dr*dv_hat(:)/2.0_dp
     force(:,k) = force(:,k) - dE_dr*dv_hat(:)/2.0_dp
    enddo

    do jk = 1, atoms_n_neighbours(at, j)
     k = atoms_neighbour(at, j, jk, dist_k, cosines = dv_hat)
     if (dist_k .feq. 0.0_dp) cycle
     dE_dr = TBModel_GSP_Vrep_env_deriv_wk(this,dist,dist_k,tj,ti,emb_j,lambda_ij)
     force(:,j) = force(:,j) + dE_dr*dv_hat(:)/2.0_dp
     force(:,k) = force(:,k) - dE_dr*dv_hat(:)/2.0_dp
    enddo
  end do
end function TBModel_GSP_get_local_rep_E_force

function TBModel_GSP_get_local_rep_E_virial(this, at, i) result(virial)
  type(TBModel_GSP), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: virial(3,3)
  real(dp) :: lambda_ij

  real(dp) :: emb_i, emb_j 
  real(dp) :: dE_dr, dist, dist_k, dv_hat(3)
  integer  :: ji, j, ti, tj, k, ik, jk

  virial = 0.0_dp

  ti = get_type(this%type_of_atomic_num, at%Z(i))
  emb_i = TBModel_GSP_Vrep_env_emb(this, at, i, ti)
  do ji=1, atoms_n_neighbours(at, i)
    j = atoms_neighbour(at, i, ji, dist, cosines = dv_hat)
    if (dist .feq. 0.0_dp) cycle
    tj = get_type(this%type_of_atomic_num, at%Z(j))
    emb_j = TBModel_GSP_Vrep_env_emb(this, at, j, tj)
    lambda_ij = TBModel_GSP_Vrep_env_lambda(this,ti,tj,emb_i,emb_j)

    dE_dr =  TBModel_GSP_Vrep_dist_func_deriv(this, dist, ti, tj) + &
             TBModel_GSP_Vrep_env_deriv_ij(this,dist,ti,tj,lambda_ij) 
    virial = virial - dE_dr*(dv_hat .outer. dv_hat) * dist / 2.0_dp

    do ik = 1, atoms_n_neighbours(at, i)
     k = atoms_neighbour(at, i, ik, dist_k, cosines = dv_hat)
     if (dist_k .feq. 0.0_dp) cycle
     dE_dr = TBModel_GSP_Vrep_env_deriv_wk(this,dist,dist_k,ti,tj,emb_i,lambda_ij)
     virial = virial - dE_dr*(dv_hat .outer. dv_hat) * dist_k / 2.0_dp
    enddo

    do jk = 1, atoms_n_neighbours(at, j)
     k = atoms_neighbour(at, j, jk, dist_k, cosines = dv_hat)
     if (dist_k .feq. 0.0_dp) cycle
     dE_dr = TBModel_GSP_Vrep_env_deriv_wk(this,dist,dist_k,tj,ti,emb_j,lambda_ij)
     virial = virial - dE_dr*(dv_hat .outer. dv_hat) * dist_k / 2.0_dp
    enddo
  end do

end function TBModel_GSP_get_local_rep_E_virial

end module TBModel_GSP_module

