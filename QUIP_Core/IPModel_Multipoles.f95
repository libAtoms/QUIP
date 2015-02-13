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
!% Multipolesised interatomic potential: 
!%
!% Energy and Force routines are hardwired
!% Cutoff is hardwired
!% 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_Multipoles_module

use error_module
use system_module, only : dp, inoutput, print, operator(//)
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

public :: IPModel_Multipoles
type IPModel_Multipoles
  real(dp) :: cutoff = 0.0_dp
  real(dp) :: ewald_error
  real(dp) :: smooth_coulomb_cutoff


end type IPModel_Multipoles

type :: Ewald_arrays
  real(dp), dimension(:,:),allocatable :: Q,P,cosk,sink,kvectors,etc
end type Ewald_arrays

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
  call param_register(params, 'kconf', '0.0', this%kconf, help_string='strength of quadratic confinement potential on O atoms. potential is kconf*(rO)^2')
  if(.not. param_read_line(params, args_str, ignore_unknown=.true., task='IPModel_Multipoles_Initialise args_str')) then
     RAISE_ERROR("IPModel_Multipoles_Init failed to parse args_str='"//trim(args_str)//"'", error)
  end if
  call finalise(params)

end subroutine IPModel_Multipoles_Initialise_str

subroutine IPModel_Multipoles_Finalise(this)
  type(IPModel_Multipoles), intent(inout) :: this

  ! Add finalisation code here

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

   call multipole_sites_setup(atoms,this,dummy_atoms,multipoles) ! multipoles includes exclude_list

   if (this%pbc_method == PBC_Method_Ewald) then
     call ewald_setup(dummy_atoms,multipoles,ewald)
   end if

   if (this%polarisation_method /= Polarisation_None) then
     call calc_permanent_e_field(dummy_atoms,multipoles,ewald,e_field)
     call build_polarisation_matrix(dummy_atoms,multipoles,ewald,pol_matrix) ! A^-1-T
     call calc_induced_dipoles(pol_matrix,e_field,induced_dipoles) ! solve linear sys 
     call update_dipoles(multipoles,induced_dipoles)
   end if

   call calc_potentials_fields_grads(dummy_atoms,multipoles,ewald)
   call calc_site_energies_forces(multipoles)
   call sites_to_atoms(atoms,multipoles,e,f)

   ! clean up

end subroutine IPModel_Multipoles_Calc


subroutine IPModel_Multipoles_Print(this, file)
  type(IPModel_Multipoles), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  call Print("IPModel_Multipoles : Multipoles Potential", file=file)
  call Print("IPModel_Multipoles : cutoff = " // this%cutoff, file=file)
  call Print("IPModel_Multipoles : kconf = " // this%kconf, file=file)

end subroutine IPModel_Multipoles_Print


end module IPModel_Multipoles_module
