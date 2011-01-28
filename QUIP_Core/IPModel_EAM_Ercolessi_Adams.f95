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
!X IPModel_EAM_ErcolAd module 
!X  
!% Module for the Embedded Atom Potential of Liu, Ercolessi, Adams. 
!% (Ref. Liu, Ercolessi, Adams, {\it Modelling Simul. Mater. Sci. Eng.} {\bf 12}, 665-670, (2004)).
!% In this potential no analytic function is assumed for the actractive and repulsive terms.
!% Each function is described as a set of points, whose values are the IP parameters.
!% A cubic spline function is then used for interpolation between those points.
!% The 'IPModel_EAM_ErcolAd' object contains all the parameters (the grid points) 
!% read from an 'EAM_ErcolAd_params' XML stanza.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_EAM_ErcolAd_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

#define PAIR
#define EMBED

include 'IPModel_interface.h'

public :: IPModel_EAM_ErcolAd
type IPModel_EAM_ErcolAd
  integer :: n_types = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)

  real(dp) :: cutoff = 0.0_dp

  real(dp), allocatable :: r_min(:,:), r_cut(:,:)
  real(dp), allocatable :: spline_V(:,:,:,:), spline_rho(:,:,:), spline_F(:,:,:) !% Cubic spline for interpolation between the mesh points.

  type(Spline), allocatable :: V(:,:),rho(:),F(:)
  real(dp), allocatable :: V_F_shift(:)

  character(len=FIELD_LENGTH) :: label

end type IPModel_EAM_ErcolAd

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_EAM_ErcolAd), private, pointer :: parse_ip
integer :: parse_cur_type_i, parse_cur_type_j, parse_cur_point
logical :: parse_in_spline_V, parse_in_spline_rho, parse_in_spline_F

interface Initialise
  module procedure IPModel_EAM_ErcolAd_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_EAM_ErcolAd_Finalise
end interface Finalise

interface Print
  module procedure IPModel_EAM_ErcolAd_Print
end interface Print

interface Calc
  module procedure IPModel_EAM_ErcolAd_Calc
end interface Calc

contains

subroutine IPModel_EAM_ErcolAd_Initialise_str(this, args_str, param_str)
  type(IPModel_EAM_ErcolAd), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_EAM_ErcolAd_Initialise_str args_str')) then
    call system_abort("IPModel_EAM_ErcolAd_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_EAM_ErcolAd_read_params_xml(this, param_str)

end subroutine IPModel_EAM_ErcolAd_Initialise_str

subroutine IPModel_EAM_ErcolAd_Finalise(this)
  type(IPModel_EAM_ErcolAd), intent(inout) :: this

  integer i, j

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if (allocated(this%r_min)) deallocate(this%r_min)
  if (allocated(this%r_cut)) deallocate(this%r_cut)

  if (allocated(this%spline_V)) deallocate(this%spline_V)
  if (allocated(this%spline_rho)) deallocate(this%spline_rho)
  if (allocated(this%spline_F)) deallocate(this%spline_F)

  if (allocated(this%rho)) then
    do i=1,this%n_types
      call Finalise(this%rho(i))
    end do
    deallocate(this%rho)
  endif

  if (allocated(this%F)) then
    do i=1,this%n_types
      call Finalise(this%F(i))
    end do
    deallocate(this%F)
  endif

  if (allocated(this%V)) then
    do i=1,this%n_types
    do j=1,this%n_types
      call finalise(this%V(i,j))
    end do
    end do
    deallocate(this%V)
  endif

  this%n_types = 0
  this%label = ''
end subroutine IPModel_EAM_ErcolAd_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator.  It computes energy, forces and virial.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_EAM_ErcolAd_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_EAM_ErcolAd), intent(inout) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial
  character(len=*), optional      :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  real(dp) :: de, rho_i, drho_i_dri(3), drho_i_drj(3), drho_i_drij_outer_rij(3,3), w_f, V_r, rho_r
  real(dp), pointer :: w_e(:)
  integer :: i, ji, j, ti, tj
  real(dp) :: r_ij_mag, r_ij_hat(3)
  real(dp) :: F_n, dF_n
  real(dp) :: spline_rho_d_val, spline_V_d_val, virial_factor(3,3), virial_i(3,3)

  type(Dictionary) :: params
  logical, dimension(:), pointer :: atom_mask_pointer
  logical :: has_atom_mask_name
  character(FIELD_LENGTH) :: atom_mask_name
  real(dp) :: r_scale, E_scale
  logical :: do_rescale_r, do_rescale_E

  INIT_ERROR(error)

  if (present(e)) then
     e = 0.0_dp
  endif

  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_EAM_ErcolAd_Calc', error)
     local_e = 0.0_dp
  endif

  if (present(f)) then 
     call check_size('Force',f,(/3,at%N/),'IPModel_EAM_ErcolAd_Calc', error)
     f = 0.0_dp
  end if

  if (present(virial)) then
     virial = 0.0_dp
  endif

  if (present(local_virial)) then
     call check_size('Local_virial',local_virial,(/9,at%N/),'IPModel_EAM_ErcolAd_Calc', error)
     local_virial = 0.0_dp
  endif

  if (.not. assign_pointer(at, "weight", w_e)) nullify(w_e)

  atom_mask_pointer => null()
  if(present(args_str)) then
     call initialise(params)
     call param_register(params, 'atom_mask_name', 'NONE',atom_mask_name,has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
     call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")

     if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='IPModel_EAM_ErcolAd_Calc args_str')) &
     call system_abort("IPModel_EAM_ErcolAd_Calc failed to parse args_str='"//trim(args_str)//"'")
     call finalise(params)

     if( has_atom_mask_name ) then
        if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) &
        call system_abort("IPModel_EAM_ErcolAd_Calc did not find "//trim(atom_mask_name)//" propery in the atoms object.")
     else
        atom_mask_pointer => null()
     endif
  else
     do_rescale_r = .false.
     do_rescale_E = .false.
  endif

  if (do_rescale_r) call print('IPModel_Tersoff_Calc: rescaling distances by factor '//r_scale, PRINT_VERBOSE)
  if (do_rescale_E) call print('IPModel_Tersoff_Calc: rescaling energy by factor '//E_scale, PRINT_VERBOSE)

  do i=1, at%N
    if (present(mpi)) then
       if (mpi%active) then
	 if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
       endif
    endif

    if(associated(atom_mask_pointer)) then
       if(.not. atom_mask_pointer(i)) cycle
    endif

    ti = get_type(this%type_of_atomic_num, at%Z(i))

    w_f = 1.0_dp

    rho_i = 0.0_dp
    drho_i_dri = 0.0_dp
    drho_i_drij_outer_rij = 0.0_dp
    do ji=1, atoms_n_neighbours(at, i)
      j = atoms_neighbour(at, i, ji, r_ij_mag, cosines = r_ij_hat)
      if (r_ij_mag .feq. 0.0_dp) cycle

      if (do_rescale_r) r_ij_mag = r_ij_mag*r_scale

      tj = get_type(this%type_of_atomic_num, at%Z(j))

      if (associated(w_e)) w_f = 0.5_dp*(w_e(i)+w_e(j))

      if (r_ij_mag < this%r_min(ti,tj)) cycle
      if (r_ij_mag >= this%r_cut(ti,tj)) cycle

      V_r = eam_spline_V(this, ti, tj, r_ij_mag)
      rho_r = eam_spline_rho(this, tj, r_ij_mag)

      V_r = V_r - 2.0_dp*this%V_F_shift(ti)*rho_r

      rho_i = rho_i + rho_r

      de = 0.5_dp*V_r
#ifdef PAIR
      if (present(e)) e = e + de*w_f
      if (present(local_e)) local_e(i) = local_e(i) + de
#endif

      if (present(f) .or. present(virial) .or. present(local_virial)) then
         spline_rho_d_val = eam_spline_rho_d(this,tj,r_ij_mag)
         spline_V_d_val = eam_spline_V_d(this,ti,tj,r_ij_mag)
         spline_V_d_val = spline_V_d_val - 2.0_dp*this%V_F_shift(ti)*spline_rho_d_val
         if (present(f)) then
            drho_i_dri = drho_i_dri + spline_rho_d_val*r_ij_hat
#ifdef PAIR
            f(:,i) = f(:,i) + 0.5_dp*w_f*spline_V_d_val*r_ij_hat
            f(:,j) = f(:,j) - 0.5_dp*w_f*spline_V_d_val*r_ij_hat
#endif
         endif
         if (present(virial) .or. present(local_virial)) then
            virial_factor = (r_ij_hat .outer. r_ij_hat)*r_ij_mag
            drho_i_drij_outer_rij = drho_i_drij_outer_rij + spline_rho_d_val*virial_factor
            virial_i = 0.5_dp * w_f*spline_V_d_val*virial_factor
         endif
#ifdef PAIR
         if (present(virial) ) virial = virial - virial_i
         if (present(local_virial) ) local_virial(:,i) = local_virial(:,i) - reshape(virial_i,(/9/))
#endif
      endif

    enddo ! ji

    if (associated(w_e)) w_f = w_e(i)

#ifdef EMBED
    if (present(e) .or. present(local_e)) then
      F_n = eam_spline_F(this, ti, rho_i)
      F_n = F_n + this%V_F_shift(ti)*rho_i
      de = F_n
      if (present(e)) e = e + de*w_f
      if (present(local_e)) local_e(i) = local_e(i) + de
    endif

    if (present(f) .or. present(virial) .or. present(local_virial)) then
      dF_n = eam_spline_F_d(this, ti, rho_i) 
      dF_n = dF_n + this%V_F_shift(ti)
      if (present(f)) f(:,i) = f(:,i) + w_f*dF_n*drho_i_dri
      if (present(virial) .or. present(local_virial)) virial_i = w_f*dF_n*drho_i_drij_outer_rij
      if (present(virial))  virial = virial - virial_i
      if (present(local_virial)) local_virial(:,i) = local_virial(:,i) - reshape(virial_i,(/9/))

      if (present(f)) then
	! cross terms for forces
	do ji=1, atoms_n_neighbours(at, i)
	  j = atoms_neighbour(at, i, ji, r_ij_mag, cosines = r_ij_hat)
	  if (r_ij_mag .feq. 0.0_dp) cycle

	  if (do_rescale_r) r_ij_mag = r_ij_mag*r_scale

	  tj = get_type(this%type_of_atomic_num, at%Z(j))

	  drho_i_drj = -eam_spline_rho_d(this, tj, r_ij_mag)*r_ij_hat
	  f(:,j) = f(:,j) + w_f*dF_n*drho_i_drj
	end do
      end if

    endif
#endif

  enddo ! i

  if (present(mpi)) then
     if (present(e)) e = sum(mpi, e)
     if (present(local_e)) call sum_in_place(mpi, local_e)
     if (present(f)) call sum_in_place(mpi, f)
     if (present(virial)) call sum_in_place(mpi, virial)
     if (present(local_virial)) call sum_in_place(mpi, local_virial)
  endif

  if (do_rescale_r) then
     if (present(f)) f = f*r_scale
  end if

  if (do_rescale_E) then
     if (present(e)) e = e*E_scale
     if (present(local_e)) local_e = local_e*E_scale
     if (present(f)) f = f*E_scale
     if (present(virial)) virial=virial*E_scale
     if (present(local_virial)) local_virial=local_virial*E_scale
  end if

end subroutine IPModel_EAM_ErcolAd_Calc

function eam_spline_V(this, ti, tj, r)
  type(IPModel_EAM_ErcolAd), intent(in) :: this
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: r
  real(dp) :: eam_spline_V

  if (r < min_knot(this%V(ti,tj)) .or. r >= max_knot(this%V(ti,tj))) then
    eam_spline_V = 0.0_dp
  else
    eam_spline_V = spline_value(this%V(ti,tj),r)
  endif

end function eam_spline_V

function eam_spline_rho(this, ti, r)
  type(IPModel_EAM_ErcolAd), intent(in) :: this
  integer, intent(in) :: ti
  real(dp), intent(in) :: r
  real(dp) :: eam_spline_rho

  if (r < min_knot(this%rho(ti)) .or. r >= max_knot(this%rho(ti))) then
    eam_spline_rho = 0.0_dp
  else
    eam_spline_rho = spline_value(this%rho(ti),r)
  endif

end function eam_spline_rho

function eam_spline_F(this, ti, rho)
  type(IPModel_EAM_ErcolAd), intent(in) :: this
  integer, intent(in) :: ti
  real(dp), intent(in) :: rho
  real(dp) :: eam_spline_F

  if (rho < min_knot(this%F(ti)) .or. rho >= max_knot(this%F(ti))) then
    eam_spline_F = 0.0_dp
  else
    eam_spline_F = spline_value(this%F(ti),rho)
  endif

end function eam_spline_F

function eam_spline_V_d(this, ti, tj, r)
  type(IPModel_EAM_ErcolAd), intent(in) :: this
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: r
  real(dp) :: eam_spline_V_d

  if (r < min_knot(this%V(ti,tj)) .or. r >= max_knot(this%V(ti,tj))) then
    eam_spline_V_d = 0.0_dp
  else
    eam_spline_V_d = spline_deriv(this%V(ti,tj),r)
  endif

end function eam_spline_V_d

function eam_spline_rho_d(this, ti, r)
  type(IPModel_EAM_ErcolAd), intent(in) :: this
  integer, intent(in) :: ti
  real(dp), intent(in) :: r
  real(dp) :: eam_spline_rho_d

  if (r < min_knot(this%rho(ti)) .or. r >= max_knot(this%rho(ti))) then
    eam_spline_rho_d = 0.0_dp
  else
    eam_spline_rho_d = spline_deriv(this%rho(ti),r)
  endif

end function eam_spline_rho_d

function eam_spline_F_d(this, ti, rho)
  type(IPModel_EAM_ErcolAd), intent(in) :: this
  integer, intent(in) :: ti
  real(dp), intent(in) :: rho
  real(dp) :: eam_spline_F_d

  if (rho < min_knot(this%F(ti)) .or. rho >= max_knot(this%F(ti))) then
    eam_spline_F_d = 0.0_dp
  else
    eam_spline_F_d = spline_deriv(this%F(ti),rho)
  endif

end function eam_spline_F_d

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!X XML param reader functions
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer status
  character(len=1024) :: value

  integer ti, Zi, Zj
  integer n_types, n_spline_V, n_spline_rho, n_spline_F
  real(dp) :: v

  if (name == 'EAM_ErcolAd_params') then ! new Ercolessi Adams stanza

    if (parse_in_ip) &
      call system_abort("IPModel_startElement_handler entered EAM_ErcolAd_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

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

      call QUIP_FoX_get_value(attributes, "n_types", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find n_types")
      read (value, *) parse_ip%n_types
      n_types = parse_ip%n_types

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      allocate(parse_ip%r_min(n_types,n_types))
      allocate(parse_ip%r_cut(n_types,n_types))
      parse_ip%r_min = 0.0_dp
      parse_ip%r_cut = 0.0_dp

      allocate(parse_ip%V_F_shift(n_types))
      parse_ip%V_F_shift = 0.0_dp

      call QUIP_FoX_get_value(attributes, "n_spline_V", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find n_spline_V")
      read (value, *) n_spline_V
      call QUIP_FoX_get_value(attributes, "n_spline_rho", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find n_spline_rho")
      read (value, *) n_spline_rho
      call QUIP_FoX_get_value(attributes, "n_spline_F", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find n_spline_F")
      read (value, *) n_spline_F

      allocate(parse_ip%spline_V(5,n_spline_V,n_types,n_types))
      allocate(parse_ip%spline_rho(5,n_spline_rho,n_types))
      allocate(parse_ip%spline_F(5,n_spline_F,n_types))

      allocate(parse_ip%V(n_types,n_types))
      allocate(parse_ip%rho(n_types))
      allocate(parse_ip%F(n_types))
    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find type")
    read (value, *) ti

    if (ti < 1 .or. ti > parse_ip%n_types) call system_abort("IPModel_EAM_ErcolAd_read_params_xml got" // &
      " per_type_data type out of range " // ti // " n_types " // parse_ip%n_types)

    parse_cur_type_i = ti
    parse_cur_type_j = ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    call QUIP_FoX_get_value(attributes, "V_F_shift", value, status)
    if (status == 0) then
      read (value, *) parse_ip%V_F_shift(parse_cur_type_i)
    endif

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
	parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "atomic_num_i", value, status)
    if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find atomic_num_i")
    read (value, *) Zi
    call QUIP_FoX_get_value(attributes, "atomic_num_j", value, status)
    if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find atomic_num_j")
    read (value, *) Zj

    parse_cur_type_i = get_type(parse_ip%type_of_atomic_num,Zi)
    parse_cur_type_j = get_type(parse_ip%type_of_atomic_num,Zj)

    call QUIP_FoX_get_value(attributes, "r_min", value, status)
    if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find r_min")
    read (value, *) parse_ip%r_min(parse_cur_type_i, parse_cur_type_j)
    call QUIP_FoX_get_value(attributes, "r_cut", value, status)
    if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find r_cut")
    read (value, *) parse_ip%r_cut(parse_cur_type_i, parse_cur_type_j)

  elseif (parse_in_ip .and. name == 'spline_V') then
    parse_in_spline_V = .true.
    parse_cur_point = 1
  elseif (parse_in_ip .and. name == 'spline_rho') then
    parse_in_spline_rho = .true.
    parse_cur_point = 1
  elseif (parse_in_ip .and. name == 'spline_F') then
    parse_in_spline_F = .true.
    parse_cur_point = 1
  elseif (parse_in_ip .and. name == 'point') then

    if (parse_in_spline_V) then

      if (parse_cur_point > size(parse_ip%spline_V,2)) call system_abort ("IPModel_EAM_ErcolAd got too " // &
	"many points " // parse_cur_point // " type " // parse_cur_type_i // " " // parse_cur_type_j // " in_spline_V")

      call QUIP_FoX_get_value(attributes, "r", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find r")
      read (value, *) v
      parse_ip%spline_V(1,parse_cur_point,parse_cur_type_i,parse_cur_type_j) = v

      call QUIP_FoX_get_value(attributes, "y", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find y")
      read (value, *) v
      parse_ip%spline_V(2,parse_cur_point,parse_cur_type_i,parse_cur_type_j) = v

      call QUIP_FoX_get_value(attributes, "b", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find b")
      read (value, *) v
      parse_ip%spline_V(3,parse_cur_point,parse_cur_type_i,parse_cur_type_j) = v

      call QUIP_FoX_get_value(attributes, "c", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find c")
      read (value, *) v
      parse_ip%spline_V(4,parse_cur_point,parse_cur_type_i,parse_cur_type_j) = v

      call QUIP_FoX_get_value(attributes, "d", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find d")
      read (value, *) v
      parse_ip%spline_V(5,parse_cur_point,parse_cur_type_i,parse_cur_type_j) = v

      if (parse_cur_type_i /= parse_cur_type_j) then
	parse_ip%spline_V(1:5,parse_cur_point,parse_cur_type_j,parse_cur_type_i) = &
	  parse_ip%spline_V(1:5,parse_cur_point,parse_cur_type_i,parse_cur_type_j)
      endif

    else if (parse_in_spline_rho) then

      if (parse_cur_point > size(parse_ip%spline_rho,2)) call system_abort ("IPModel_EAM_ErcolAd got " // &
	"too many points " // parse_cur_point // " type " // parse_cur_type_i // " in_spline_rho")

      call QUIP_FoX_get_value(attributes, "r", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find r")
      read (value, *) v
      parse_ip%spline_rho(1,parse_cur_point,parse_cur_type_i) = v

      call QUIP_FoX_get_value(attributes, "y", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find y")
      read (value, *) v
      parse_ip%spline_rho(2,parse_cur_point,parse_cur_type_i) = v

      call QUIP_FoX_get_value(attributes, "b", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find b")
      read (value, *) v
      parse_ip%spline_rho(3,parse_cur_point,parse_cur_type_i) = v

      call QUIP_FoX_get_value(attributes, "c", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find c")
      read (value, *) v
      parse_ip%spline_rho(4,parse_cur_point,parse_cur_type_i) = v

      call QUIP_FoX_get_value(attributes, "d", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find d")
      read (value, *) v
      parse_ip%spline_rho(5,parse_cur_point,parse_cur_type_i) = v

    else if (parse_in_spline_F) then

      if (parse_cur_point > size(parse_ip%spline_F,2)) call system_abort ("IPModel_EAM_ErcolAd got " // &
	"too many points " // parse_cur_point // " type " // parse_cur_type_i // " in_spline_F")

      call QUIP_FoX_get_value(attributes, "r", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find r")
      read (value, *) v
      parse_ip%spline_F(1,parse_cur_point,parse_cur_type_i) = v

      call QUIP_FoX_get_value(attributes, "y", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find y")
      read (value, *) v
      parse_ip%spline_F(2,parse_cur_point,parse_cur_type_i) = v

      call QUIP_FoX_get_value(attributes, "b", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find b")
      read (value, *) v
      parse_ip%spline_F(3,parse_cur_point,parse_cur_type_i) = v

      call QUIP_FoX_get_value(attributes, "c", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find c")
      read (value, *) v
      parse_ip%spline_F(4,parse_cur_point,parse_cur_type_i) = v

      call QUIP_FoX_get_value(attributes, "d", value, status)
      if (status /= 0) call system_abort ("IPModel_EAM_ErcolAd_read_params_xml cannot find d")
      read (value, *) v
      parse_ip%spline_F(5,parse_cur_point,parse_cur_type_i) = v
    endif

    parse_cur_point = parse_cur_point + 1
  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if (name == 'EAM_ErcolAd_params') then
      parse_in_ip = .false.
    elseif (parse_in_ip .and. name == 'spline_V') then
      parse_in_spline_V = .false.
    elseif (parse_in_ip .and. name == 'spline_rho') then
      parse_in_spline_rho = .false.
    elseif (parse_in_ip .and. name == 'spline_F') then
      parse_in_spline_F = .false.
    endif
  endif

  end subroutine IPModel_endElement_handler

subroutine IPModel_EAM_ErcolAd_read_params_xml(this, param_str)
  type(IPModel_EAM_ErcolAd), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type (xml_t) :: fxml
  integer ti, tj

  if (len(trim(param_str)) <= 0) return

  parse_ip => this
  parse_in_ip = .false.
  parse_matched_label = .false.

  call open_xml_string(fxml, param_str)

  call parse(fxml, &
    startElement_handler = IPModel_startElement_handler, &
    endElement_handler = IPModel_endElement_handler)

  call close_xml_t(fxml)

  if (this%n_types == 0) call system_abort("EAM_ErcolAd Tried to parse, but n_types still 0")

  do ti=1, this%n_types
    call initialise(this%rho(ti),this%spline_rho(1,:,ti),this%spline_rho(2,:,ti),huge(1.0_dp),0.0_dp)
    call initialise(this%F(ti),this%spline_F(1,:,ti),this%spline_F(2,:,ti),huge(1.0_dp),huge(1.0_dp))
    do tj=1, this%n_types
      call initialise(this%V(ti,tj),this%spline_V(1,:,ti,tj),this%spline_V(2,:,ti,tj),huge(1.0_dp),0.0_dp)
    end do
  end do

  parse_ip%cutoff = maxval(parse_ip%r_cut)

end subroutine IPModel_EAM_ErcolAd_read_params_xml

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of potential parameters: number of different types, cutoff radius, atomic numbers, ect.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_EAM_ErcolAd_Print (this, file)
  type(IPModel_EAM_ErcolAd), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_EAM_ErcolAd : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print("IPModel_EAM_ErcolAd : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    call Print("IPModel_EAM_ErcolAd : spline_rho ", file=file)
    call Print(this%rho(ti), file=file)
    call Print("IPModel_EAM_ErcolAd : spline_F ", file=file)
    call Print(this%F(ti), file=file)
    do tj=1, this%n_types
      call Print("IPModel_EAM_ErcolAd : pair "// ti // " " // tj // " r_min " // this%r_min(ti,tj) // " r_cut " // this%r_cut(ti,tj), file=file)
      call Print("IPModel_EAM_ErcolAd : pair spline_V ", file=file)
      call Print(this%V(ti,tj), file=file)
    end do
    call verbosity_pop()
  end do
    
end subroutine IPModel_EAM_ErcolAd_Print

end module IPModel_EAM_ErcolAd_module
