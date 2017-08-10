! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   Portions of this code were written by
! H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! H0 X
! H0 X   Copyright 2006-2017.
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
!X IPModel_DispTS
!X
!% Pairwise dispersion correction from Tkatchenko and Scheffler,
!% PRL 102(7) 073005 (2009).
!%
!% Free-atom reference data must be supplied in the XML file; Hirshfeld volumes
!% must be supplied in the Atoms object.
!%
!% Free-atom reference data is available e.g. from Chu and Dalgarno,
!% JCP 121, 4083 (2004) (used in original implementation) or from
!% Gould and Bučko, JCTC 12(8), 3603 (2016) (better coverage)
!%
!% Tail corrections for the energy and virial are available for
!% this potential.  Doing them properly, however, requires
!% an integral over the smooth cutoff function; this involves the
!% the special Ci and Si functions (sine and cosine integrals).
!% Those are not yet implemented here, so instead this potential
!% asks for a precomputed constant that depends only on the cutoff
!% and transition width:
!% \[
!%     I = -\frac{2}{3} \pi (r_{in}^{-3} - 3 \int_{r_{in}}^{r_{out}} r^{-4} S(r) dr)
!% \]
!% where r_{out} is this potential's cutoff, r_{in} is the cutoff
!% minus the cutoff_transition_width, and S(r) is the switching
!% function that goes smoothly from 1 at r_{in} to 0 at r_{out}.
!% In this potential, S(r) is the first half of a cosine (shifted
!% and scaled).
!%
!% To use tail corrections, compute this integral using any
!% method you like and supply it ($I$) as 'tail_correction_const'
!% in the parameter file or on the command line (or use an existing
!% parameter file that already has it).
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_DispTS_module

use error_module
use system_module, only : dp, inoutput, print, verbosity_push_decrement, verbosity_pop, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_DispTS
type IPModel_DispTS
  integer :: n_types = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)
  real(dp), allocatable :: c6_free(:), alpha_free(:), r_vdW_free(:)
  logical :: only_inter_resid = .false.
  logical :: do_tail_corrections = .true.

  real(dp) :: cutoff = 0.0_dp
  real(dp) :: cutoff_transition_width = 0.5_dp

  real(dp) :: tail_correction_const = 0.0_dp

  ! Constants for now...
  real(dp) :: damp_steepness = 20.0_dp
  real(dp) :: damp_scale = 0.94 !PBE

  character(len=STRING_LENGTH) :: label

end type IPModel_DispTS

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_DispTS), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_DispTS_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_DispTS_Finalise
end interface Finalise

interface Print
  module procedure IPModel_DispTS_Print
end interface Print

interface Calc
  module procedure IPModel_DispTS_Calc
end interface Calc

contains

subroutine IPModel_DispTS_Initialise_str(this, args_str, param_str)
  type(IPModel_DispTS), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="Label to identify the potential")
  call param_register(params, 'only_inter_resid', 'F', this%only_inter_resid, &
      help_string="If True, only calculate interactions between atoms with different ResIDs (requires 'resid' property to be present)")
  !call param_register(params, 'do_tail_corrections', this%do_tail_corrections, &
  !    help_string="If True, apply long-range corrections to the energy and virial")
  call param_register(params, 'tail_correction_const', '0.0', this%tail_correction_const, has_value_target=this%do_tail_corrections, &
      help_string="Constant used to calculate tail corrections.  Both the energy and virial corrections are equal to this value divided by the cell volume &
      multiplied by the sum of all pairwise C6 coefficients.  Units: Ang^-3.")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_DispTS_Initialise_str args_str')) then
    call system_abort("IPModel_DispTS_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_DispTS_read_params_xml(this, param_str)

  !  Add initialisation code here

end subroutine IPModel_DispTS_Initialise_str

subroutine IPModel_DispTS_Finalise(this)
  type(IPModel_DispTS), intent(inout) :: this

  ! Add finalisation code here

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%c6_free)) deallocate(this%c6_free)
  if (allocated(this%alpha_free)) deallocate(this%alpha_free)
  if (allocated(this%r_vdW_free)) deallocate(this%r_vdW_free)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_DispTS_Finalise


subroutine IPModel_DispTS_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
    type(IPModel_DispTS), intent(inout):: this
    type(Atoms), intent(inout)      :: at
    real(dp), intent(out), optional :: e, local_e(:)
    real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
    real(dp), intent(out), optional :: virial(3,3)

    character(len=*), optional      :: args_str
    type(Dictionary) :: params
    logical :: has_atom_mask_name
    character(STRING_LENGTH) :: atom_mask_name
    real(dp) :: r_scale, E_scale
    logical :: do_rescale_r, do_rescale_E
    type(MPI_Context), intent(in), optional :: mpi
    integer, intent(out), optional :: error

    integer, pointer, dimension(:)  :: resid
    real(dp), pointer, dimension(:) :: my_hirshfeld_volume
    character(STRING_LENGTH)        :: hirshfeld_vol_name

    integer     :: i, j, d, ji, ti, tj
    real(dp)    :: dr(3), dr_mag, de, de_dr, vi, vj
    real(dp)    :: damp, dfdamp
    real(dp)    :: tail_correction, c6_sum

    INIT_ERROR(error)

    if (present(e)) e = 0.0_dp
    if (present(local_e)) then
       call check_size('Local_E',local_e,(/at%N/),'IPModel_DispTS_Calc', error)
       local_e = 0.0_dp
    endif
    if (present(f)) then
       call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_DispTS_Calc', error)
       f = 0.0_dp
    end if
    if (present(virial)) virial = 0.0_dp
    if (present(local_virial)) then
       call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_DispTS_Calc', error)
       local_virial = 0.0_dp
       RAISE_ERROR("IPModel_DispTS_Calc: local_virial calculation requested but not supported yet.", error)
    endif

    tail_correction = 0.0_dp
    c6_sum = 0.0_dp

    if (this%only_inter_resid) then
       if (.not. assign_pointer(at, "resid", resid)) then
         RAISE_ERROR("IPModel_DispTS_Calc calculation with only_inter_resid=T requires resid field", error)
       endif
    end if

    if (present(args_str)) then
        if (len_trim(args_str) > 0) then
            call initialise(params)
            call param_register(params, 'hirshfeld_vol_name', 'hirshfeld_rel_volume', hirshfeld_vol_name, &
                                help_string='Name of the Atoms property containing relative Hirshfeld volumes $v/v_{free}$')
            call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
            call param_register(params, 'r_scale', '1.0', r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
            call param_register(params, 'E_scale', '1.0', E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")
 
            if (.not. param_read_line(params, args_str, ignore_unknown=.true., task='IPModel_MBD_Calc args_str')) then
                RAISE_ERROR("IPModel_DispTS_Calc failed to parse args_str '"//trim(args_str)//"'", error)
            endif
            call finalise(params)
            call assign_property_pointer(at, trim(hirshfeld_vol_name), my_hirshfeld_volume, error)
            PASS_ERROR_WITH_INFO("IPModel_DispTS_Calc could not find '"//trim(hirshfeld_vol_name)//"' property in the Atoms object", error)
        endif
        if(has_atom_mask_name) then
           RAISE_ERROR('IPModel_LJ_Calc: atom_mask_name found, but not supported', error)
        endif
        if (do_rescale_r .or. do_rescale_E) then
           RAISE_ERROR("IPModel_LJ_Calc: rescaling of potential with r_scale and E_scale not yet implemented!", error)
        end if
    else
        call assign_property_pointer(at, 'hirshfeld_rel_volume', my_hirshfeld_volume, error)
        PASS_ERROR_WITH_INFO("IPModel_DispTS_Calc could not find 'hirshfeld_rel_volume' property in the Atoms object", error)
    endif

    ! Adapted from IPModel_LJ.f95 as a general pair potential
    do i = 1, at%N

      if (present(mpi)) then
         if (mpi%active) then
           if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
         endif
      endif

      do ji = 1, n_neighbours(at, i)
        j = neighbour(at, i, ji, dr_mag, cosines = dr)

        if (dr_mag .feq. 0.0_dp) cycle
        if ((i < j)) cycle

        if (this%only_inter_resid) then
          if (resid(i) == resid(j)) cycle
        endif

        ti = get_type(this%type_of_atomic_num, at%Z(i))
        tj = get_type(this%type_of_atomic_num, at%Z(j))
        vi = my_hirshfeld_volume(i)
        vj = my_hirshfeld_volume(j)
        c6_sum = c6_sum + IPModel_DispTS_pair_c6(this, ti, tj, vi, vj)

        damp = 0.0_dp
        dfdamp = 0.0_dp
        if (present(f) .or. present(virial)) then
            call IPModel_DispTS_fdamp(this, ti, tj, vi, vj, dr_mag, damp, dfdamp)
        else if (present(e) .or. present(local_e)) then
            call IPModel_DispTS_fdamp(this, ti, tj, vi, vj, dr_mag, damp)
        endif

        if (present(e) .or. present(local_e)) then
          de = IPModel_DispTS_pairenergy(this, ti, tj, vi, vj, dr_mag, damp)

          if (present(local_e)) then
            local_e(i) = local_e(i) + 0.5_dp*de
            if(i/=j) local_e(j) = local_e(j) + 0.5_dp*de
          endif
          if (present(e)) then
            if(i==j) then
               e = e + 0.5_dp*de
            else
               e = e + de
            endif
          endif
        endif
        if (present(f) .or. present(virial)) then
          de_dr = IPModel_DispTS_pairenergy_deriv(this, ti, tj, vi, vj, dr_mag, damp, dfdamp)
          if (present(f)) then
            f(:,i) = f(:,i) + de_dr*dr
            if(i/=j) f(:,j) = f(:,j) - de_dr*dr
          endif
          if (present(virial)) then
             if(i==j) then
                virial = virial - 0.5_dp*de_dr*(dr .outer. dr)*dr_mag
             else
                virial = virial - de_dr*(dr .outer. dr)*dr_mag
             endif
          endif
        endif
      end do
   end do

   c6_sum = c6_sum * 2.0_dp ! Convention is to double count this sum

   if (present(mpi)) then
      if (present(e)) e = sum(mpi, e)
      if (present(local_e)) call sum_in_place(mpi, local_e)
      if (present(virial)) call sum_in_place(mpi, virial)
      if (present(f)) call sum_in_place(mpi, f)
   endif

   if (this%do_tail_corrections) then
      tail_correction = c6_sum * this%tail_correction_const / cell_volume(at)
      if (present(e)) e = e + tail_correction
      if (present(virial)) then
         do d = 1, 3
            virial(d, d) = virial(d, d) + tail_correction
         enddo
      endif
   endif

end subroutine IPModel_DispTS_Calc

subroutine IPModel_DispTS_fdamp(this, ti, tj, vi, vj, r, damp, dfdamp)
    type(IPModel_DispTS), intent(in) :: this
    integer, intent(in) :: ti, tj
    real(dp), intent(in) :: r, vi, vj  ! Pair distance and relative Hirshfeld volumes
    real(dp), intent(out), optional :: damp
    real(dp), intent(out), optional :: dfdamp

    real(dp) :: rfree ! Inner cutoff distance

    rfree = this%r_vdW_free(ti)*vi**(1/3.0_dp) + this%r_vdW_free(tj)*vj**(1/3.0_dp)
    if (present(damp)) then
        damp = 1.0_dp / (1.0_dp + exp(-1.0_dp * this%damp_steepness * (r/rfree/this%damp_scale - 1.0_dp)))
    endif
    if (present(dfdamp)) then
        dfdamp = this%damp_steepness / (2.0_dp * this%damp_scale * rfree) &
            / (cosh(this%damp_steepness * (r/rfree/this%damp_scale - 1.0_dp)) + 1.0_dp)
    endif

end subroutine IPModel_DispTS_fdamp


function IPModel_DispTS_pair_c6(this, ti, tj, vi, vj)
    type(IPModel_DispTS), intent(in) :: this
    integer, intent(in) :: ti, tj
    real(dp), intent(in) :: vi, vj  ! relative Hirshfeld volumes
    real(dp) :: IPModel_DispTS_pair_c6

    real(dp) :: c6, c6i, c6j ! VdW R^-6 coefficients
    real(dp) :: alphai, alphaj

    c6i = this%c6_free(ti) * vi**2
    c6j = this%c6_free(tj) * vj**2
    alphai = this%alpha_free(ti) * vi
    alphaj = this%alpha_free(tj) * vj
    c6 = 2*c6i*c6j / (alphai/alphaj * c6j + alphaj/alphai * c6i)
    IPModel_DispTS_pair_c6 = c6

end function IPModel_DispTS_pair_c6

function IPModel_DispTS_pairenergy(this, ti, tj, vi, vj, r, damp)
    type(IPModel_DispTS), intent(in) :: this
    integer, intent(in) :: ti, tj
    real(dp), intent(in) :: r, vi, vj  ! Pair distance and relative Hirshfeld volumes
    real(dp), intent(in) :: damp ! Damping factor (precalculated)
    real(dp) :: IPModel_DispTS_pairenergy

    real(dp) :: c6 ! VdW R^-6 coefficient
    real(dp) :: cut ! Smooth cutoff function

    c6 = IPModel_DispTS_pair_c6(this, ti, tj, vi, vj)
    cut = coordination_function(r, this%cutoff, this%cutoff_transition_width)
    IPModel_DispTS_pairenergy = -1.0_dp * c6 * r**(-6) * damp * cut

end function IPModel_DispTS_pairenergy


function IPModel_DispTS_pairenergy_deriv(this, ti, tj, vi, vj, r, damp, dfdamp)
    type(IPModel_DispTS), intent(in) :: this
    integer, intent(in) :: ti, tj
    real(dp), intent(in) :: r, vi, vj  ! Pair distance and relative Hirshfeld volumes
    real(dp), intent(in) :: damp, dfdamp ! Damping function and derivative
    real(dp) :: IPModel_DispTS_pairenergy_deriv

    real(dp) :: c6 ! VdW R^-6 coefficient
    real(dp) :: cut, dcut ! Cutoff function and derivative

    c6 = IPModel_DispTS_pair_c6(this, ti, tj, vi, vj)
    cut = coordination_function(r, this%cutoff, this%cutoff_transition_width)
    dcut = dcoordination_function(r, this%cutoff, this%cutoff_transition_width)

    IPModel_DispTS_pairenergy_deriv = -1.0_dp * c6 * r**(-6) &
        * (dfdamp*cut + damp*dcut - 6.0_dp/r * damp*cut)

end function IPModel_DispTS_pairenergy_deriv


subroutine IPModel_DispTS_Print(this, file)
  type(IPModel_DispTS), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti

  call Print("IPModel_DispTS : T-S dispersion correction potential", file=file)
  call Print("IPModel_DispTS : n_types = " // this%n_types // " cutoff = " // this%cutoff // &
             " transition width = " // this%cutoff_transition_width, file=file)

  do ti=1, this%n_types
    call Print("IPModel_DispTS : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call Print("IPModel_DispTS : type " // ti // " free-atom: C6 " // this%c6_free(ti) &
               // " polarizability " // this%alpha_free(ti) &
               // " vdW radius " // this%r_vdW_free(ti), file=file)
  end do

end subroutine IPModel_DispTS_Print

subroutine IPModel_DispTS_read_params_xml(this, param_str)
  type(IPModel_DispTS), intent(inout), target :: this
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
    call system_abort("IPModel_DispTS_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_DispTS_read_params_xml

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

  if (name == 'DispTS_params') then ! new DispTS stanza

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
        call system_abort("Can't find n_types in DispTS_params")
      endif

      call QUIP_FoX_get_value(attributes, 'only_inter_resid', value, status)
      if (status == 0) then
         read (value, *), parse_ip%only_inter_resid
      else
         parse_ip%only_inter_resid = .false.
      endif
      call QUIP_FoX_get_value(attributes, 'tail_correction_const', value, status)
      if (status == 0) then
         read (value, *), parse_ip%tail_correction_const
         parse_ip%do_tail_corrections = .true.
      else
         parse_ip%tail_correction_const = 0.0_dp
         parse_ip%do_tail_corrections = .false.
      endif


      allocate(parse_ip%atomic_num(parse_ip%n_types))
      allocate(parse_ip%c6_free(parse_ip%n_types))
      allocate(parse_ip%alpha_free(parse_ip%n_types))
      allocate(parse_ip%r_vdW_free(parse_ip%n_types))
      parse_ip%atomic_num = 0
      parse_ip%c6_free = 0.0_dp
      parse_ip%alpha_free = 0.0_dp
      parse_ip%r_vdW_free = 0.0_dp

      call QUIP_FoX_get_value(attributes, "cutoff", value, status)
      if (status /= 0) call system_abort ("IPModel_DispTS_read_params_xml cannot find cutoff")
      read (value, *) parse_ip%cutoff
      call QUIP_FoX_get_value(attributes, "cutoff_transition_width", value, status)
      if (status == 0) then
          read (value, *) parse_ip%cutoff_transition_width
      endif
    endif


  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_DispTS_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_DispTS_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    call QUIP_FoX_get_value(attributes, "c6_free", value, status)
    if (status /= 0) call system_abort ("IPModel_DispTS_read_params_xml cannot find c6_free")
    read (value, *) parse_ip%c6_free(ti)

    call QUIP_FoX_get_value(attributes, "alpha_free", value, status)
    if (status /= 0) call system_abort ("IPModel_DispTS_read_params_xml cannot find alpha_free")
    read (value, *) parse_ip%alpha_free(ti)

    call QUIP_FoX_get_value(attributes, "r_vdW_free", value, status)
    if (status /= 0) call system_abort ("IPModel_DispTS_read_params_xml cannot find r_vdW_free")
    read (value, *) parse_ip%r_vdW_free(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
        parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do


  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'DispTS_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_DispTS_module
