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
!X TBModel module 
!X
!% General object which handles the possible TB scheme, addressing
!% the calls to the right modules. 
!% The available models are the following
!%   \begin{itemize}
!%    \item    "NRL_TB", whose related object is \texttt{TBModel_NRL_TB}
!%    \item    "Bowler",  whose related object is \texttt{TBModel_Bowler} 
!%    \item    "DFTB" , whose related object is \texttt{TBModel_DFTB}
!%    \item    "GSP" , whose related object is \texttt{TBModel_GSP}
!%   \end{itemize}
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module TBModel_module

use system_module, only : dp, print, inoutput, PRINT_VERBOSE, system_abort, operator(//)
use dictionary_module
use paramreader_module
use atoms_module

use QUIP_Common_module
use TBModel_NRL_TB_module
use TBModel_Bowler_module
use TBModel_DFTB_module
use TBModel_GSP_module

implicit none
private

integer, parameter :: FF_NONE = 0, FF_NRL_TB = 1, FF_Bowler = 2, FF_DFTB = 3, FF_GSP = 4
character(len=6) :: ff_names(0:4) = (/ &
      'NONE  ', &
      'NRL_TB', &
      'Bowler', &
      'DFTB  ', &
      'GSP   '/)

include 'TBModel_interface.h'

public :: TBModel
type TBModel
  real(dp) cutoff

  logical is_orthogonal
  integer :: functional_form = FF_NONE

  type(TBModel_NRL_TB) tbmodel_nrl_tb
  type(TBModel_Bowler) tbmodel_bowler
  type(TBModel_DFTB)   tbmodel_dftb
  type(TBModel_GSP)    tbmodel_gsp

  logical :: has_default_fermi_E, has_default_fermi_T, has_default_band_width, has_default_k_density
  real(dp) :: default_fermi_E, default_fermi_T, default_band_width, default_k_density

end type TBModel

interface Initialise
  module procedure TBModel_Initialise_str
end interface Initialise

interface Finalise
  module procedure TBModel_Finalise
end interface Finalise

interface Print
  module procedure TBModel_Print
end interface Print

interface n_orbs_of_Z
  module procedure TBModel_n_orbs_of_Z
end interface n_orbs_of_Z

interface n_orb_sets_of_Z
  module procedure TBModel_n_orb_sets_of_Z
end interface n_orb_sets_of_Z

interface n_orbs_of_orb_set_of_Z
  module procedure TBModel_n_orbs_of_orb_set_of_Z
end interface n_orbs_of_orb_set_of_Z

interface orb_type_of_orb_set_of_Z
  module procedure TBModel_orb_type_of_orb_set_of_Z
end interface orb_type_of_orb_set_of_Z

interface n_elecs_of_Z
  module procedure TBModel_n_elecs_of_Z
end interface n_elecs_of_Z

interface get_HS_blocks
  module procedure TBModel_get_HS_blocks
end interface get_HS_blocks

interface get_dHS_masks
  module procedure TBModel_get_dHS_masks
end interface get_dHS_masks

interface get_dHS_blocks
  module procedure TBModel_get_dHS_blocks
end interface get_dHS_blocks

interface get_local_rep_E
  module procedure TBModel_get_local_rep_E
end interface get_local_rep_E

interface get_local_rep_E_force
  module procedure TBModel_get_local_rep_E_force
end interface get_local_rep_E_force

interface get_local_rep_E_virial
  module procedure TBModel_get_local_rep_E_virial
end interface get_local_rep_E_virial

public :: has_fermi_E
interface has_fermi_E
  module procedure TBModel_has_fermi_E
end interface has_fermi_E

public :: has_fermi_T
interface has_fermi_T
  module procedure TBModel_has_fermi_T
end interface has_fermi_T

public :: has_band_width
interface has_band_width
  module procedure TBModel_has_band_width
end interface has_band_width

public :: has_k_density
interface has_k_density
  module procedure TBModel_has_k_density
end interface has_k_density

contains

subroutine TBModel_Initialise_str(this, args_str, param_str)
  type(TBModel), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params
  logical :: is_nrl_tb, is_bowler, is_dftb, is_gsp
  real(dp) :: str_fermi_e, str_fermi_T, str_band_width, str_k_density
  logical, target :: has_str_fermi_e, has_str_fermi_T, has_str_band_width, has_str_k_density

  call Initialise(params)
  call param_register(params, 'NRL-TB', 'false', is_nrl_tb, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'Bowler', 'false', is_bowler, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'DFTB', 'false', is_dftb, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'GSP', 'false', is_gsp, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'fermi_e', '0.0', str_fermi_e, has_value_target=has_str_fermi_e, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'fermi_T', '0.0', str_fermi_T, has_value_target=has_str_fermi_T, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'band_width', '0.0', str_band_width, has_value_target=has_str_band_width, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'k_density', '0.0', str_k_density, has_value_target=has_str_k_density, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='TBModel_Initialise_str args_str')) then
    call system_abort("TBModel_Initialise_str failed to parse args_str='"//trim(args_str)//"'")
  endif
  call finalise(params)

  if (count((/is_nrl_tb, is_bowler, is_dftb, is_gsp/)) /= 1) then
    call system_abort("TBModel_Initialise_str found too few or too many TB Model types args_str='"//trim(args_str)//"'")
  endif

  if (is_nrl_tb) then
    this%functional_form = FF_NRL_TB
    call Initialise(this%tbmodel_nrl_tb, args_str, param_str)
    this%cutoff=this%tbmodel_nrl_tb%cutoff
    this%is_orthogonal = this%tbmodel_nrl_tb%is_orthogonal
    this%has_default_fermi_e = this%tbmodel_nrl_tb%has_default_fermi_e
    this%has_default_fermi_T = this%tbmodel_nrl_tb%has_default_fermi_T
    this%has_default_band_width = this%tbmodel_nrl_tb%has_default_band_width
    this%has_default_k_density = this%tbmodel_nrl_tb%has_default_k_density
    this%default_fermi_e = this%tbmodel_nrl_tb%default_fermi_e
    this%default_fermi_T = this%tbmodel_nrl_tb%default_fermi_T
    this%default_band_width = this%tbmodel_nrl_tb%default_band_width
    this%default_k_density = this%tbmodel_nrl_tb%default_k_density
  else if (is_bowler) then
    this%functional_form = FF_Bowler
    call Initialise(this%tbmodel_bowler, args_str, param_str)
    this%cutoff=this%tbmodel_bowler%cutoff
    this%is_orthogonal = this%tbmodel_bowler%is_orthogonal
    this%has_default_fermi_e = this%tbmodel_bowler%has_default_fermi_e
    this%has_default_fermi_T = this%tbmodel_bowler%has_default_fermi_T
    this%has_default_band_width = this%tbmodel_bowler%has_default_band_width
    this%has_default_k_density = this%tbmodel_bowler%has_default_k_density
    this%default_fermi_e = this%tbmodel_bowler%default_fermi_e
    this%default_fermi_T = this%tbmodel_bowler%default_fermi_T
    this%default_band_width = this%tbmodel_bowler%default_band_width
    this%default_k_density = this%tbmodel_bowler%default_k_density
  else  if(is_dftb) then
    this%functional_form = FF_DFTB
    call Initialise(this%tbmodel_dftb, args_str, param_str)
    this%cutoff=this%tbmodel_dftb%cutoff
    this%is_orthogonal = this%tbmodel_dftb%is_orthogonal
    this%has_default_fermi_e = this%tbmodel_dftb%has_default_fermi_e
    this%has_default_fermi_T = this%tbmodel_dftb%has_default_fermi_T
    this%has_default_band_width = this%tbmodel_dftb%has_default_band_width
    this%has_default_k_density = this%tbmodel_dftb%has_default_k_density
    this%default_fermi_e = this%tbmodel_dftb%default_fermi_e
    this%default_fermi_T = this%tbmodel_dftb%default_fermi_T
    this%default_band_width = this%tbmodel_dftb%default_band_width
    this%default_k_density = this%tbmodel_dftb%default_k_density
  else if (is_gsp) then
    this%functional_form = FF_GSP
    call Initialise(this%tbmodel_gsp, args_str, param_str)
    this%cutoff=this%tbmodel_gsp%cutoff
    this%is_orthogonal = this%tbmodel_gsp%is_orthogonal
    this%has_default_fermi_e = this%tbmodel_gsp%has_default_fermi_e
    this%has_default_fermi_T = this%tbmodel_gsp%has_default_fermi_T
    this%has_default_band_width = this%tbmodel_gsp%has_default_band_width
    this%has_default_k_density = this%tbmodel_gsp%has_default_k_density
    this%default_fermi_e = this%tbmodel_gsp%default_fermi_e
    this%default_fermi_T = this%tbmodel_gsp%default_fermi_T
    this%default_band_width = this%tbmodel_gsp%default_band_width
    this%default_k_density = this%tbmodel_gsp%default_k_density
  end if

  if (has_str_fermi_e) then
    this%has_default_fermi_e = .true.
    this%default_fermi_e = str_fermi_e
  endif
  if (has_str_fermi_T) then
    this%has_default_fermi_T = .true.
    this%default_fermi_T = str_fermi_T
  endif
  if (has_str_band_width) then
    this%has_default_band_width = .true.
    this%default_band_width = str_band_width
  endif
  if (has_str_k_density) then
    this%has_default_k_density = .true.
    this%default_k_density = str_k_density
  endif

end subroutine TBModel_Initialise_str

subroutine TBModel_Finalise(this)
  type(TBModel), intent(inout) :: this

  select case(this%functional_form)
    case (FF_NRL_TB)
      call Finalise(this%tbmodel_nrl_tb)
    case (FF_Bowler)
      call Finalise(this%tbmodel_bowler)
    case (FF_DFTB)
      call Finalise(this%tbmodel_dftb)
    case (FF_GSP)
      call Finalise(this%tbmodel_gsp)
    case (FF_NONE)
      return
    case default
      call system_abort ('TBModel_Finalise confused by functional_form' // this%functional_form)
  end select

end subroutine TBModel_Finalise

function TBModel_n_orbs_of_Z(this, Z)
  type(TBModel), intent(in) :: this
  integer, intent(in) :: Z
  integer TBModel_n_orbs_of_Z

  select case(this%functional_form)
    case (FF_NRL_TB)
      TBModel_n_orbs_of_Z = n_orbs_of_Z(this%tbmodel_nrl_tb,Z)
    case (FF_Bowler)
      TBModel_n_orbs_of_Z = n_orbs_of_Z(this%tbmodel_bowler,Z)
    case (FF_DFTB)
      TBModel_n_orbs_of_Z = n_orbs_of_Z(this%tbmodel_dftb,Z)
    case (FF_GSP)
      TBModel_n_orbs_of_Z = n_orbs_of_Z(this%tbmodel_gsp,Z)
    case default
      call system_abort ('TBModel_n_orbs_of_Z confused by functional_form' // this%functional_form)
  end select

end function TBModel_n_orbs_of_Z

function TBModel_n_orb_sets_of_Z(this, Z)
  type(TBModel), intent(in) :: this
  integer, intent(in) :: Z
  integer TBModel_n_orb_sets_of_Z

  select case(this%functional_form)
    case (FF_NRL_TB)
      TBModel_n_orb_sets_of_Z = n_orb_sets_of_Z(this%tbmodel_nrl_tb,Z)
    case (FF_Bowler)
      TBModel_n_orb_sets_of_Z = n_orb_sets_of_Z(this%tbmodel_bowler,Z)
    case (FF_DFTB)
      TBModel_n_orb_sets_of_Z = n_orb_sets_of_Z(this%tbmodel_dftb,Z)
    case (FF_GSP)
      TBModel_n_orb_sets_of_Z = n_orb_sets_of_Z(this%tbmodel_gsp,Z)
    case default
      call system_abort ('TBModel_n_orb_sets_of_Z confused by functional_form' // this%functional_form)
  end select

end function TBModel_n_orb_sets_of_Z

function TBModel_n_orbs_of_orb_set_of_Z(this, Z, i_set)
  type(TBModel), intent(in) :: this
  integer, intent(in) :: Z, i_set
  integer TBModel_n_orbs_of_orb_set_of_Z

  select case(this%functional_form)
    case (FF_NRL_TB)
      TBModel_n_orbs_of_orb_set_of_Z = n_orbs_of_orb_set_of_Z(this%tbmodel_nrl_tb,Z,i_set)
    case (FF_Bowler)
      TBModel_n_orbs_of_orb_set_of_Z = n_orbs_of_orb_set_of_Z(this%tbmodel_bowler,Z,i_set)
    case (FF_DFTB)
      TBModel_n_orbs_of_orb_set_of_Z = n_orbs_of_orb_set_of_Z(this%tbmodel_dftb,Z,i_set)
    case (FF_GSP)
      TBModel_n_orbs_of_orb_set_of_Z = n_orbs_of_orb_set_of_Z(this%tbmodel_gsp,Z,i_set)
    case default
      call system_abort ('TBModel_n_orbs_of_orb_set_of_Z confused by functional_form' // this%functional_form)
  end select

end function TBModel_n_orbs_of_orb_set_of_Z

function TBModel_orb_type_of_orb_set_of_Z(this, Z, i_set)
  type(TBModel), intent(in) :: this
  integer, intent(in) :: Z, i_set
  integer TBModel_orb_type_of_orb_set_of_Z

  select case(this%functional_form)
    case (FF_NRL_TB)
      TBModel_orb_type_of_orb_set_of_Z = orb_type_of_orb_set_of_Z(this%tbmodel_nrl_tb,Z,i_set)
    case (FF_Bowler)
      TBModel_orb_type_of_orb_set_of_Z = orb_type_of_orb_set_of_Z(this%tbmodel_bowler,Z,i_set)
    case (FF_DFTB)
      TBModel_orb_type_of_orb_set_of_Z = orb_type_of_orb_set_of_Z(this%tbmodel_dftb,Z,i_set)
    case (FF_GSP)
      TBModel_orb_type_of_orb_set_of_Z = orb_type_of_orb_set_of_Z(this%tbmodel_gsp,Z,i_set)
    case default
      call system_abort ('TBModel_orb_type_of_orb_set_of_Z confused by functional_form' // this%functional_form)
  end select

end function TBModel_orb_type_of_orb_set_of_Z

function TBModel_n_elecs_of_Z(this, Z)
  type(TBModel), intent(in) :: this
  integer, intent(in) :: Z
  real(dp) TBModel_n_elecs_of_Z

  select case(this%functional_form)
    case (FF_NRL_TB)
      TBModel_n_elecs_of_Z = n_elecs_of_Z(this%tbmodel_nrl_tb,Z)
    case (FF_Bowler)
      TBModel_n_elecs_of_Z = n_elecs_of_Z(this%tbmodel_bowler,Z)
    case (FF_DFTB)
      TBModel_n_elecs_of_Z = n_elecs_of_Z(this%tbmodel_dftb,Z)
    case (FF_GSP)
      TBModel_n_elecs_of_Z = n_elecs_of_Z(this%tbmodel_gsp,Z)
    case default
      call system_abort ('TBModel_n_elecs_of_Z confused by functional_form' // this%functional_form)
  end select

end function TBModel_n_elecs_of_Z

subroutine TBModel_Print(this,file)
  type(TBModel),    intent(in)           :: this
  type(Inoutput), intent(inout),optional :: file

  call print("TBModel: is_orthogonal " // this%is_orthogonal, PRINT_VERBOSE)
  call print("TBModel: functional_form " // ff_names(this%functional_form), PRINT_VERBOSE)
  call print("TBModel: has_default_fermi_e " // this%has_default_fermi_e // " default_fermi_e " // &
    this%default_fermi_e, PRINT_VERBOSE)
  call print("TBModel: has_default_fermi_T " // this%has_default_fermi_T // " default_fermi_T " // &
    this%default_fermi_T, PRINT_VERBOSE)
  call print("TBModel: has_default_band_width " // this%has_default_band_width // " default_band_width " // &
    this%default_band_width, PRINT_VERBOSE)
  call print("TBModel: has_default_k_density " // this%has_default_k_density // " default_k_density " // &
    this%default_k_density, PRINT_VERBOSE)

  select case(this%functional_form)
    case (FF_NRL_TB)
      call Print(this%tbmodel_nrl_tb, file=file)
    case (FF_Bowler)
      call Print(this%tbmodel_bowler, file=file)
    case (FF_DFTB)
      call Print(this%tbmodel_dftb, file=file)
    case (FF_GSP)
      call Print(this%tbmodel_gsp, file=file)
    case default
      call system_abort ('TBModel_Print confused by functional_form' // this%functional_form)
  end select

end subroutine

subroutine TBModel_get_HS_blocks(this, at, i, j, dv_hat, dv_mag, b_H, b_S, i_mag)
  type(TBModel), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i, j
  real(dp), intent(in) :: dv_hat(3), dv_mag
  real(dp), intent(out) :: b_H(:,:), b_S(:,:)
  integer, intent(in), optional :: i_mag

  select case(this%functional_form)
    case (FF_NRL_TB)
      call get_HS_blocks(this%tbmodel_nrl_tb, at, i, j, dv_hat, dv_mag, b_H, b_S, i_mag)
    case (FF_Bowler)
      call get_HS_blocks(this%tbmodel_bowler, at, i, j, dv_hat, dv_mag, b_H, b_S, i_mag)
    case (FF_DFTB)
      call get_HS_blocks(this%tbmodel_dftb, at, i, j, dv_hat, dv_mag, b_H, b_S, i_mag)
    case (FF_GSP)
      call get_HS_blocks(this%tbmodel_gsp, at, i, j, dv_hat, dv_mag, b_H, b_S)
    case default
      call system_abort ('TBModel_get_HS_blocks confused by functional_form' // this%functional_form)
  end select

end subroutine

subroutine TBModel_get_dHS_masks(this, at, at_ind, d_mask, od_mask)
  type(TBModel), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_ind
  logical, intent(out), optional :: d_mask(:), od_mask(:)

  select case(this%functional_form)
    case (FF_NRL_TB)
      call get_dHS_masks(this%tbmodel_nrl_tb, at, at_ind, d_mask, od_mask)
    case (FF_Bowler)
      call get_dHS_masks(this%tbmodel_bowler, at, at_ind, d_mask, od_mask)
    case (FF_DFTB)
      call get_dHS_masks(this%tbmodel_dftb, at, at_ind, d_mask, od_mask)
    case (FF_GSP)
      call get_dHS_masks(this%tbmodel_gsp, at_ind, d_mask, od_mask)
    case default
      call system_abort ('TBModel_get_HS_masks confused by functional_form' // this%functional_form)
  end select

end subroutine

function TBModel_get_dHS_blocks(this, at, i, j, dv_hat, dv_mag, at_ind, b_dH, b_dS, i_mag)
  type(TBModel), intent(inout) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i, j
  real(dp), intent(in) :: dv_hat(3), dv_mag
  integer, intent(in) :: at_ind
  real(dp), intent(out) :: b_dH(:,:,:), b_dS(:,:,:)
  integer, intent(in), optional :: i_mag
  logical TBModel_get_dHS_blocks

  select case(this%functional_form)
    case (FF_NRL_TB)
      TBModel_get_dHS_blocks = get_dHS_blocks(this%tbmodel_nrl_tb, at, i, j, dv_hat, dv_mag, at_ind, b_dH, b_dS, i_mag)
    case (FF_Bowler)
      TBModel_get_dHS_blocks = get_dHS_blocks(this%tbmodel_bowler, at, i, j, dv_hat, dv_mag, at_ind, b_dH, b_dS, i_mag)
    case (FF_DFTB)
      TBModel_get_dHS_blocks = get_dHS_blocks(this%tbmodel_dftb, at, i, j, dv_hat, dv_mag, at_ind, b_dH, b_dS, i_mag)
    case (FF_GSP)
      TBModel_get_dHS_blocks = get_dHS_blocks(this%tbmodel_gsp, at, i, j, dv_hat, dv_mag, at_ind, b_dH, b_dS)
    case default
      call system_abort ('TBModel_get_HS_blocks confused by functional_form' // this%functional_form)
  end select

end function TBModel_get_dHS_blocks

function TBModel_get_local_rep_E(this, at, i)
  type(TBModel), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: TBModel_get_local_rep_E

  select case(this%functional_form)
    case (FF_NRL_TB)
      TBModel_get_local_rep_E = get_local_rep_E(this%tbmodel_nrl_tb, at, i)
    case (FF_Bowler)
      TBModel_get_local_rep_E = get_local_rep_E(this%tbmodel_bowler, at, i)
    case (FF_DFTB)
      TBModel_get_local_rep_E = get_local_rep_E(this%tbmodel_dftb, at, i)
    case (FF_GSP)
      TBModel_get_local_rep_E = get_local_rep_E(this%tbmodel_gsp, at, i)
    case default
      call system_abort ('TBModel_get_local_rep_E confused by functional_form' // this%functional_form)
  end select
end function TBModel_get_local_rep_E

function TBModel_get_local_rep_E_force(this, at, i) result(force)
  type(TBModel), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: force(3,at%N)

  select case(this%functional_form)
    case (FF_NRL_TB)
      force = get_local_rep_E_force(this%tbmodel_nrl_tb, at, i)
    case (FF_Bowler)
      force = get_local_rep_E_force(this%tbmodel_bowler, at, i)
    case (FF_DFTB)
      force = get_local_rep_E_force(this%tbmodel_dftb, at, i)
    case (FF_GSP)
      force = get_local_rep_E_force(this%tbmodel_gsp, at, i)
    case default
      call system_abort ('TBModel_get_local_rep_E_force confused by functional_form' // this%functional_form)
  end select
end function TBModel_get_local_rep_E_force

function TBModel_get_local_rep_E_virial(this, at, i) result(virial)
  type(TBModel), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i
  real(dp) :: virial(3,3)

  select case(this%functional_form)
    case (FF_NRL_TB)
      virial = get_local_rep_E_virial(this%tbmodel_nrl_tb, at, i)
    case (FF_Bowler)
      virial = get_local_rep_E_virial(this%tbmodel_bowler, at, i)
    case (FF_DFTB)
      virial = get_local_rep_E_virial(this%tbmodel_dftb, at, i)
    case (FF_GSP)
      virial = get_local_rep_E_virial(this%tbmodel_gsp, at, i)
    case default
      call system_abort ('TBModel_get_local_rep_E_virial confused by functional_form' // this%functional_form)
  end select
end function TBModel_get_local_rep_E_virial

function TBModel_has_Fermi_E(out_Fermi_E, this, Fermi_E, args_str)
  real(dp), intent(out) :: out_Fermi_E
  type(TBModel), intent(in) :: this
  real(dp), optional, intent(in) :: Fermi_E
  character(len=*), optional, intent(in) :: args_str
  logical :: TBModel_has_Fermi_E

  type(Dictionary) :: params
  logical :: has_str_fermi_e

  TBModel_has_Fermi_E = .false.
  if (present(Fermi_E)) then
    out_Fermi_E = Fermi_E
    TBModel_has_Fermi_E = .true.
  else
    if (this%has_default_fermi_E) then
      out_Fermi_E = this%default_Fermi_E
      TBModel_has_Fermi_E = .true.
    endif
    if (present(args_str)) then
      call initialise(params)
      call param_register(params, 'fermi_e', ''//out_fermi_E, out_fermi_E, has_value_target=has_str_fermi_e, help_string="No help yet.  This source file was $LastChangedBy$")
      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='TBModel_has_Fermi_E args_str')) then
	call system_abort("TBModel_has_fermi_e failed to parse args_str='"//trim(args_str)//"'")
      endif
      call finalise(params)
      if (has_str_fermi_e) TBModel_has_fermi_E = .true.
    endif
  endif

end function TBModel_has_Fermi_E

function TBModel_has_Fermi_T(out_Fermi_T, this, Fermi_T, args_str)
  real(dp), intent(out) :: out_Fermi_T
  type(TBModel), intent(in) :: this
  real(dp), optional, intent(in) :: Fermi_T
  character(len=*), optional, intent(in) :: args_str
  logical :: TBModel_has_Fermi_T

  type(Dictionary) :: params
  logical :: has_str_fermi_T

  TBModel_has_Fermi_T = .false.
  if (present(Fermi_T)) then
    out_Fermi_T = Fermi_T
    TBModel_has_Fermi_T = .true.
  else
    if (this%has_default_fermi_T) then
      out_Fermi_T = this%default_Fermi_T
      TBModel_has_Fermi_T = .true.
    endif
    if (present(args_str)) then
      call initialise(params)
      call param_register(params, 'fermi_T', ''//out_fermi_T, out_fermi_T, has_value_target=has_str_fermi_T, help_string="No help yet.  This source file was $LastChangedBy$")
      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='TBModel_has_Fermi_T args_str')) then
	call system_abort("TBModel_has_fermi_T failed to parse args_str='"//trim(args_str)//"'")
      endif
      call finalise(params)
      if (has_str_fermi_T) TBModel_has_fermi_T = .true.
    endif
  endif

end function TBModel_has_Fermi_T

function TBModel_has_band_width(out_band_width, this, band_width, args_str)
  real(dp), intent(out) :: out_band_width
  type(TBModel), intent(in) :: this
  real(dp), optional, intent(in) :: band_width
  character(len=*), optional, intent(in) :: args_str
  logical :: TBModel_has_band_width

  type(Dictionary) :: params
  logical :: has_str_band_width

  TBModel_has_band_width = .false.
  if (present(band_width)) then
    out_band_width = band_width
    TBModel_has_band_width = .true.
  else
    if (this%has_default_band_width) then
      out_band_width = this%default_band_width
      TBModel_has_band_width = .true.
    endif
    if (present(args_str)) then
      call initialise(params)
      call param_register(params, 'band_width', ''//out_band_width, out_band_width, has_value_target=has_str_band_width, help_string="No help yet.  This source file was $LastChangedBy$")
      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='TBModel_has_band_width args_str')) then
	call system_abort("TBModel_has_band_width failed to parse args_str='"//trim(args_str)//"'")
      endif
      call finalise(params)
      if (has_str_band_width) TBModel_has_band_width = .true.
    endif
  endif

end function TBModel_has_band_width

function TBModel_has_k_density(out_k_density, this, k_density, args_str)
  real(dp), intent(out) :: out_k_density
  type(TBModel), intent(in) :: this
  real(dp), optional, intent(in) :: k_density
  character(len=*), optional, intent(in) :: args_str
  logical :: TBModel_has_k_density

  type(Dictionary) :: params
  logical :: has_str_k_density

  TBModel_has_k_density = .false.
  if (present(k_density)) then
    out_k_density = k_density
    TBModel_has_k_density = .true.
  else
    if (this%has_default_k_density) then
      out_k_density = this%default_k_density
      TBModel_has_k_density = .true.
    endif
    if (present(args_str)) then
      call initialise(params)
      call param_register(params, 'k_density', ''//out_k_density, out_k_density, has_value_target=has_str_k_density, help_string="No help yet.  This source file was $LastChangedBy$")
      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='TBModel_has_k_density args_str')) then
	call system_abort("TBModel_has_k_density failed to parse args_str='"//trim(args_str)//"'")
      endif
      call finalise(params)
      if (has_str_k_density) TBModel_has_k_density = .true.
    endif
  endif

end function TBModel_has_k_density

end module TBModel_module
