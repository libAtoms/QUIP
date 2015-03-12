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
!X IPModel_Tether
!X
!% Tether a selection of atoms to the origin with a spring: 
!%
!% Energy and Force routines are hardwired
!% Cutoff is hardwired
!% 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_Tether_module

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

public :: IPModel_Tether
type IPModel_Tether
  real(dp) :: cutoff = 0.0_dp
  real(dp) :: r0 = 0.0_dp
  real(dp) :: kconf = 0.0_dp
  integer,allocatable,dimension(:) :: tether_indices
end type IPModel_Tether

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_Tether), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_Tether_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Tether_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Tether_Print
end interface Print

interface Calc
  module procedure IPModel_Tether_Calc
end interface Calc

contains

subroutine IPModel_Tether_Initialise_str(this, args_str, param_str, error)
  type(IPModel_Tether), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(Dictionary) :: params
  integer, optional, intent(out):: error
  character(len=STRING_LENGTH) :: indices_string
  character(len=STRING_LENGTH), dimension(99) :: indices_fields
  integer:: i,n_tethered

  INIT_ERROR(error)
  call Finalise(this)

  call initialise(params)
  call param_register(params, 'kconf', '0.0', this%kconf, help_string='strength of quadratic confinement potential on atoms. potential is kconf*(r - r0)^2')
  call param_register(params, 'r0', '0.0', this%r0, help_string='distance at which quadratic confinement potential on atoms begins. potential is kconf*(r - r0)^2')
  call param_register(params, 'indices', PARAM_MANDATORY, indices_string, help_string="Indices (1-based) of the atoms you wish to tether, format {i1 i2 i3 ...}")

  if(.not. param_read_line(params, args_str, ignore_unknown=.true., task='IPModel_Tether_Initialise args_str')) then
     RAISE_ERROR("IPModel_Tether_Init failed to parse args_str='"//trim(args_str)//"'", error)
  end if
  call finalise(params)

  call split_string(indices_string,' ','{}',indices_fields(:),n_tethered,matching=.true.)
  allocate(this%tether_indices(n_tethered))

  do i=1,n_tethered
    this%tether_indices(i) = string_to_int(indices_fields(i))
  end do


end subroutine IPModel_Tether_Initialise_str

subroutine IPModel_Tether_Finalise(this)
  type(IPModel_Tether), intent(inout) :: this

  ! Add finalisation code here

end subroutine IPModel_Tether_Finalise


subroutine IPModel_Tether_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_Tether), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error


   real(dp) :: energy, force(3,at%N), r, disp, dr(3), origin(3)=0.0_dp
   integer :: i , n_tethered, i_teth
  
   n_tethered=size(this%tether_indices)
  
   INIT_ERROR(error)

   ! Harmonic confining potential on tethered atoms
   energy = 0.0_dp
   do i=1,n_tethered
     i_teth = this%tether_indices(i)

     r  = distance_min_image(at, i_teth , origin)
     disp = r - this%r0

     if(r .fle. this%r0) then
       force(:,i_teth) = 0.0_dp
     else
     ! energy
       energy = energy + this%kConf*disp**2
     ! force
       dr = diff_min_image(at, i_teth, origin)/r
       force(:,i_teth) = ( 2.0_dp*this%kConf*disp ) * dr 
     end if
   end do

   if (present(e)) e = energy
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_Tether_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_Tether_Calc', error)
      f = force
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_Tether_Calc', error)
      local_virial = 0.0_dp
   endif

end subroutine IPModel_Tether_Calc


subroutine IPModel_Tether_Print(this, file)
  type(IPModel_Tether), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  call Print("IPModel_Tether : Tether Potential", file=file)
  call Print("IPModel_Tether : cutoff = " // this%cutoff, file=file)
  call Print("IPModel_Tether : kconf = " // this%kconf, file=file)
  call Print("IPModel_Tether : tethered atoms = " // this%tether_indices, file=file)

end subroutine IPModel_Tether_Print


end module IPModel_Tether_module
