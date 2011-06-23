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
!X IPModel_Custom
!X
!% Customised interatomic potential: 
!%
!% Energy and Force routines are hardwired
!% Cutoff is hardwired
!% 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_Custom_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_Custom
type IPModel_Custom
  real(dp) :: cutoff = 10.0_dp
  real(dp) :: kconf = 0.0_dp
end type IPModel_Custom

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_Custom), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_Custom_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Custom_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Custom_Print
end interface Print

interface Calc
  module procedure IPModel_Custom_Calc
end interface Calc

contains

subroutine IPModel_Custom_Initialise_str(this, args_str, param_str)
  type(IPModel_Custom), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  !  Add initialisation code here

  this%kconf = 0.05_dp

end subroutine IPModel_Custom_Initialise_str

subroutine IPModel_Custom_Finalise(this)
  type(IPModel_Custom), intent(inout) :: this

  ! Add finalisation code here

end subroutine IPModel_Custom_Finalise


subroutine IPModel_Custom_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_Custom), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   ! Add calc() code here


   real(dp) :: energy, force(3,at%N)
   real(dp) :: rO1, rO2, drO1(3), drO2(3)

   INIT_ERROR(error)

   ! Harmonic confining potential on Os

   rO1 = distance_min_image(at, 1, (/0.0_dp, 0.0_dp, 0.0_dp/))
   rO2 = distance_min_image(at, 4, (/0.0_dp, 0.0_dp, 0.0_dp/))


   energy = this%kConf*rO1**2 + this%kConf*rO2**2 

   !Forces

   drO1 = diff_min_image(at, 1, (/0.0_dp, 0.0_dp, 0.0_dp/))/rO1
   drO2 = diff_min_image(at, 4, (/0.0_dp, 0.0_dp, 0.0_dp/))/rO2

   force(:,1) = 2.0_dp*this%kConf*rO1*drO1 
   force(:,4) = 2.0_dp*this%kConf*rO2*drO2 


   if (present(e)) e = energy
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_Custom_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_Custom_Calc', error)
      f = force
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_Custom_Calc', error)
      local_virial = 0.0_dp
   endif


end subroutine IPModel_Custom_Calc


subroutine IPModel_Custom_Print(this, file)
  type(IPModel_Custom), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  call Print("IPModel_Custom : Custom Potential", file=file)
  call Print("IPModel_Custom : cutoff = " // this%cutoff, file=file)
  call Print("IPModel_Custom : kconf = " // this%kconf, file=file)

end subroutine IPModel_Custom_Print


end module IPModel_Custom_module
