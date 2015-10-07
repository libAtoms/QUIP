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

use error_module
use system_module, only : dp, inoutput, print, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module

use mpi_context_module
use QUIP_Common_module

use units_module

implicit none
private

include 'IPModel_interface.h'

!  n	A_{OO}^n
!  2	−4797.380 894 912
!  3	321 106.944 415 7
!  4	−9 435 745.917 272
!  5	159 805 932.3305
!  6	−1 719 523 144.815
!  7	12 193 219 799.32
!  8	−56 992 670 350.61
!  9	169 359 294 902.4
!  10	−290 395 612 488.0
!  11	218 964 657 251.3


! XXX TODO test potential / forces; compare with as given...

! Cutoff: 5.9 Angstrom
! delta should be 0.11 Angstrom (actually 0.109172373466 ...)


public :: IPModel_Custom
type IPModel_Custom
  real(dp) :: cutoff = 0.0_dp
  real(dp) :: delta = 0.0_dp
  logical :: energyshift = .true., forceshift = .true.
  integer :: nmax = 11
  ! coefficients for n=2 to n=nmax:
  real(dp) :: coeff(10) = (/ -4797.380894912_dp, 321106.9444157_dp, -9435745.917272_dp, 159805932.3305_dp, -1719523144.815_dp, 12193219799.32_dp, -56992670350.61_dp, 169359294902.4_dp, -290395612488.0_dp, 218964657251.3_dp /)
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

subroutine IPModel_Custom_Initialise_str(this, args_str, param_str, error)
  type(IPModel_Custom), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(Dictionary) :: params
  integer, optional, intent(out):: error


  INIT_ERROR(error)
  call Finalise(this)

  call initialise(params)
  call param_register(params, 'cutoff', '6.01', this%cutoff, help_string='total cutoff for interaction calculation')
  call param_register(params, 'delta', '0.11', this%delta, help_string='distance over which to change force linearly to zero (original cutoff + delta should give total cutoff)')
  call param_register(params, 'energyshift', 'true', this%energyshift, help_string='shift energy to zero at (total) cutoff?')
  call param_register(params, 'forceshift', 'false', this%forceshift, help_string='linearly shift force to zero at cutoff? (delta needs to be 0.0 for this)')
  if(.not. param_read_line(params, args_str, ignore_unknown=.true., task='IPModel_Custom_Initialise args_str')) then
     RAISE_ERROR("IPModel_Custom_Init failed to parse args_str='"//trim(args_str)//"'", error)
  end if
  call finalise(params)

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


   real(dp) :: eval, fval
   real(dp) :: orig_cutoff, cutoff_energy, cutoff_force
   real(dp) :: energy, force(3,at%N)

   integer :: i, ji, j !, ti, tj
   real(dp) :: dr(3), dr_mag,  deltar, rinv
   integer :: n

   INIT_ERROR(error)

   if (present(e)) e = 0.0_dp
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_Custom_Calc', error)
      f = 0.0_dp
   end if
 

   energy = 0.0_dp
   force = 0.0_dp

   ! Izvekov&Voth FM spline potential, polynomial approximation
   ! cf. http://jcp.aip.org/resource/1/jcpsa6/v123/i13/p134105_s1

   orig_cutoff = this%cutoff - this%delta ! actual force/energy calculation wants original cutoff

   cutoff_energy = 0.0_dp

   if (this%energyshift .or. this%delta > 0) then
     ! calculate energy at original cutoff
     rinv = BOHR / orig_cutoff
     do n=2, this%nmax
       cutoff_energy = cutoff_energy + this%coeff(n-1) * rinv**(n-1) / (n-1)
     end do
     cutoff_energy = cutoff_energy * HARTREE ! convert to eV
   end if

   cutoff_force = 0.0_dp

   if (this%forceshift .or. this%delta > 0) then
     ! calculate force at original cutoff
     rinv = BOHR / orig_cutoff
     do n=2, this%nmax
       cutoff_force = cutoff_force + this%coeff(n-1) * rinv**n
     end do
     cutoff_force = cutoff_force * HARTREE/BOHR ! convert to eV/Angstrom
   end if


   do i=1, at%N
     !i_is_min_image = is_min_image(at, i)
     do ji=1, n_neighbours(at, i)
       j = neighbour(at, i, ji, dr_mag, cosines=dr)

       if (dr_mag .feq. 0.0_dp) cycle
       if ((i <= j)) cycle
       if (dr_mag > this%cutoff) cycle

       if (dr_mag > orig_cutoff) then
         ! linear extrapolation

         deltar = this%cutoff - dr_mag
         fval = cutoff_force * deltar / this%delta
         ! force is negative derivative of energy, and we get one minus sign from the inner
         ! derivative of deltar, so no minus sign
         eval = 0.5_dp * cutoff_force * deltar**2 / this%delta

       else
         ! tabulated polynomial potential

         rinv = BOHR / dr_mag ! take inverse, and convert to a.u.

         eval = 0.0_dp
         fval = 0.0_dp
         do n=2, this%nmax
           ! coefficients are in a.u.
           fval = fval + this%coeff(n-1) * rinv**n
           eval = eval + this%coeff(n-1) * rinv**(n-1) / (n-1)
         end do
         eval = eval * HARTREE ! convert to eV
         fval = fval * HARTREE/BOHR ! convert to eV/Angstrom

         eval = eval - cutoff_energy ! eval is now zero at orig_cutoff

         if (this%delta .feq. 0.0_dp) then
           ! linear force shift:
           eval = eval + cutoff_force * (dr_mag - orig_cutoff)
           fval = fval - cutoff_force
         else
           ! ensure energy is continuous at orig_cutoff
           eval = eval + 0.5_dp * cutoff_force * this%delta
         end if
       end if

       energy = energy + eval
       force(:,i) = force(:,i) - fval * dr
       force(:,j) = force(:,j) + fval * dr
     end do
   end do


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

  call Print("IPModel_Custom : Izvekov&Voth tabulated polynomial for one-site water", file=file)
  call Print("IPModel_Custom : cutoff = " // this%cutoff, file=file)
  call Print("IPModel_Custom : delta = " // this%delta, file=file)
  call Print("IPModel_Custom : energyshift = " // this%energyshift, file=file)
  call Print("IPModel_Custom : forceshift = " // this%forceshift, file=file)

end subroutine IPModel_Custom_Print


end module IPModel_Custom_module
