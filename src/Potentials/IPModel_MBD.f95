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
!X IPModel_MBD
!X
!% Many-body van der Waals correction from
!% Tkatchenko, DiStasio, Car, & Scheffler
!% Phys. Rev. Lett. 108, 236402 (2012).
!%
!% Requires the MBD code from http://www.fhi-berlin.mpg.de/~tkatchen/MBD/MBD.tar
!% Download and extract into ThirdParty/MBD
!% TODO the Makefile should automatically patch their UTILS.F90 and compile it
!%      in
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_MBD_module

use error_module
use system_module, only : dp, inoutput, print, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module
use MBD_UTILS

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_MBD
type IPModel_MBD
  real(dp) :: cutoff = 0.0_dp ! Cutoff in the quip sense doesn't really apply here
  integer  :: xc = 1 ! PBE
  real(dp) :: mbd_cfdm_dip_cutoff = 100.d0 ! Angstrom
  real(dp) :: mbd_supercell_cutoff= 25.d0  ! Angstrom
  real(dp) :: mbd_scs_dip_cutoff  = 120.0  ! Angstrom
  logical, dimension(3) :: mbd_scs_vacuum_axis = {.false., .false., .false.}

end type IPModel_MBD

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_MBD), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_MBD_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_MBD_Finalise
end interface Finalise

interface Print
  module procedure IPModel_MBD_Print
end interface Print

interface Calc
  module procedure IPModel_MBD_Calc
end interface Calc

contains

subroutine IPModel_MBD_Initialise_str(this, args_str, param_str, error)
  type(IPModel_MBD), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(Dictionary) :: params
  integer, optional, intent(out):: error


  INIT_ERROR(error)
  call Finalise(this)

  call initialise(params)
  ! Some of these should definitely be in the xml file, but use command line for now
  call param_register(params, 'cfdm_dip_cutoff', '100.0', this%mbd_cfdm_dip_cutoff, help_string='MBD dipole field integration cutoff')
  call param_register(params, 'scs_dip_cutoff', '120.0', this%mbd_scs_dip_cutoff, help_string='Periodic SCS integration cutoff - important for low-dim systems')
  call param_register(params, 'supercell_cutoff', '25.0', this%mbd_supercell_cutoff, help_string='Radius used to make periodic supercell - important convergence parameter')
  call param_register(params, 'scs_vacuum_axis', {.false., .false., .false.}, this%mbd_scs_vacuum_axis, help_string='Which directions should be treated as vacuum instead of periodic')
  if(.not. param_read_line(params, args_str, ignore_unknown=.true., task='IPModel_MBD_Initialise args_str')) then
     RAISE_ERROR("IPModel_MBD_Init failed to parse args_str='"//trim(args_str)//"'", error)
  end if
  call finalise(params)

end subroutine IPModel_MBD_Initialise_str

subroutine IPModel_MBD_Finalise(this)
  type(IPModel_MBD), intent(inout) :: this

  ! Add finalisation code here

end subroutine IPModel_MBD_Finalise


subroutine IPModel_MBD_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_MBD), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   ! Add calc() code here


   INIT_ERROR(error)

   if (present(e)) e = energy
   if (present(local_e)) then
       !TODO lammps usually asks for this though - is there another way?
       RAISE_ERROR("IPModel_MBD does not have local energies", error)
   endif
   ! Need finite differences for forces and virials
   if (present(f) .or. present(virial) .or. present(local_virial)) then
       RAISE_ERROR("IPModel_MBD does not yet provide analytical gradients", error)
   endif

end subroutine IPModel_MBD_Calc


subroutine IPModel_MBD_Print(this, file)
  type(IPModel_MBD), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  call Print("IPModel_MBD : Many-body dispersion", file=file)
  call Print("IPModel_MBD : mbd_supercell_cutoff = " // this%mbd_supercell_cutoff, file=file)
  ! And a few more

end subroutine IPModel_MBD_Print


end module IPModel_MBD_module
