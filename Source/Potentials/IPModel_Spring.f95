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
!X IPModel_Spring
!X
!% Tethers two groups of atoms (their Centers of Geometry) together with a spring: 
!%
!% Energy and Force routines are hardwired
!% Cutoff is hardwired
!% 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_Spring_module

use error_module
use system_module, only : dp, inoutput, print, operator(//), split_string,string_to_int
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

public :: IPModel_Spring
type IPModel_Spring
  real(dp) :: cutoff = 0.0_dp
  real(dp) :: kconf = 0.0_dp
  real(dp) :: right = 0.0_dp
  integer,allocatable,dimension(:) :: spring_indices1
  integer,allocatable,dimension(:) :: spring_indices2
end type IPModel_Spring

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_Spring), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_Spring_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Spring_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Spring_Print
end interface Print

interface Calc
  module procedure IPModel_Spring_Calc
end interface Calc

contains

subroutine IPModel_Spring_Initialise_str(this, args_str, param_str, error)
  type(IPModel_Spring), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(Dictionary) :: params
  integer, optional, intent(out):: error
  character(len=STRING_LENGTH) :: indices1_string, indices2_string
  character(len=STRING_LENGTH), dimension(99) :: indices1_fields, indices2_fields
  integer:: i,n_group1,n_group2

  INIT_ERROR(error)
  call Finalise(this)

  call initialise(params)
  call param_register(params, 'kconf', '0.0', this%kconf, help_string='strength of quadratic confinement potential on atoms. potential is kconf*(rO)^2')
  call param_register(params, 'right', '0.0', this%right, help_string='distance up to which the confinement potential is flat')
  call param_register(params, 'cutoff', '0.0', this%cutoff, help_string='cutoff')
  call param_register(params, 'indices1', PARAM_MANDATORY, indices1_string, help_string="Indices (1-based) of the first group of atoms you wish to tether, format {i1 i2 i3 ...}")
  call param_register(params, 'indices2', PARAM_MANDATORY, indices2_string, help_string="Indices (1-based) of the second group of atoms you wish to tether, format {i1 i2 i3 ...}")

  if(.not. param_read_line(params, args_str, ignore_unknown=.true., task='IPModel_Spring_Initialise args_str')) then
     RAISE_ERROR("IPModel_Spring_Init failed to parse args_str='"//trim(args_str)//"'", error)
  end if
  call finalise(params)

  call split_string(indices1_string,' ','{}',indices1_fields(:),n_group1,matching=.true.)
  call split_string(indices2_string,' ','{}',indices2_fields(:),n_group2,matching=.true.)
  allocate(this%spring_indices1(n_group1))
  allocate(this%spring_indices2(n_group2))

  do i=1,n_group1
    this%spring_indices1(i) = string_to_int(indices1_fields(i))
  end do
  do i=1,n_group2
    this%spring_indices2(i) = string_to_int(indices2_fields(i))
  end do


end subroutine IPModel_Spring_Initialise_str

subroutine IPModel_Spring_Finalise(this)
  type(IPModel_Spring), intent(inout) :: this

  ! Add finalisation code here

end subroutine IPModel_Spring_Finalise


subroutine IPModel_Spring_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_Spring), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error


   real(dp) :: energy, force(3,at%N), r, dr(3), com1(3), com2(3), theforce(3), disp
   integer :: i , n_group1, n_group2, i_group11, i_group1i, i_group21, i_group2i
  
   n_group1=size(this%spring_indices1)
   n_group2=size(this%spring_indices2)
  
   INIT_ERROR(error)

   com1 = 0.0_dp
   i_group11 = this%spring_indices1(1)
   do i=2,n_group1
     i_group1i = this%spring_indices1(i)
     com1 = com1 + diff_min_image(at, i_group11, i_group1i)
   end do
   com1 = com1 / n_group1
   com1 = com1 + at%pos(:,i_group11)

   com2 = 0.0_dp
   i_group21 = this%spring_indices2(1)
   do i=2,n_group2
     i_group2i = this%spring_indices2(i)
     com2 = com2 + diff_min_image(at, i_group21, i_group2i)
   end do
   com2 = com2 / n_group2
   com2 = com2 + at%pos(:,i_group21)

   r  = distance_min_image(at, com1, com2)
   
   energy = 0.0_dp
   force = 0.0_dp
   if (r .fgt. this%right) then
     disp = r - this%right
     dr = diff_min_image(at, com1, com2) / r
     ! Harmonic confining potential on tethered atoms
     ! energy
     energy = energy + this%kConf*disp**2

     ! force
     theforce = ( 2.0_dp*this%kConf*disp ) * dr
     do i=1,n_group1
       force(:,this%spring_indices1(i)) = theforce / n_group1
     end do
     do i=1,n_group2
       force(:,this%spring_indices2(i)) = - theforce / n_group2
     end do
   end if

   if (present(e)) e = energy
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_Spring_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_Spring_Calc', error)
      f = force
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_Spring_Calc', error)
      local_virial = 0.0_dp
   endif

end subroutine IPModel_Spring_Calc


subroutine IPModel_Spring_Print(this, file)
  type(IPModel_Spring), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  call Print("IPModel_Spring : Spring Potential", file=file)
  call Print("IPModel_Spring : cutoff = " // this%cutoff, file=file)
  call Print("IPModel_Spring : kconf = " // this%kconf, file=file)
  call Print("IPModel_Spring : right = " // this%right, file=file)
  call Print("IPModel_Spring : group 1 atoms = " // this%spring_indices1, file=file)
  call Print("IPModel_Spring : group 2 atoms = " // this%spring_indices2, file=file)

end subroutine IPModel_Spring_Print


end module IPModel_Spring_module
