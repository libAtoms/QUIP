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
!X  ApproxFermi module 
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module ApproxFermi_module

use System_module
use linearalgebra_module ! for .feq.
use Units_module
use Functions_module

implicit none

private

public :: ApproxFermi
type ApproxFermi
  integer :: n_poles = 0
  complex(dp), allocatable :: a(:), z(:)
  real(dp) :: Fermi_E, band_width
end type ApproxFermi

public :: Initialise
interface Initialise
  module procedure ApproxFermi_Initialise, ApproxFermi_Initialise_lowlevel
end interface Initialise

public :: Finalise
interface Finalise
  module procedure ApproxFermi_Finalise
end interface Finalise

public :: Wipe
interface Wipe
  module procedure ApproxFermi_Wipe
end interface Wipe

public :: Print
interface Print
  module procedure ApproxFermi_Print
end interface Print

public :: approx_f_fermi, approx_f_fermi_deriv

contains

subroutine ApproxFermi_Finalise(this)
  type(ApproxFermi), intent(inout) :: this

  call Wipe(this)

end subroutine

subroutine ApproxFermi_Wipe(this)
  type(ApproxFermi), intent(inout) :: this

  if (allocated(this%a)) deallocate(this%a)
  if (allocated(this%z)) deallocate(this%z)
  this%n_poles = 0

end subroutine

subroutine ApproxFermi_Print(this,out)
  type(ApproxFermi),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: out

  integer::i

  call Print('ApproxFermi: Fermi_E ' // this%Fermi_E // ' band_width ' // this%band_width // ' n_poles ' // this%n_poles, file=out)

  if (this%n_poles > 0) then
    do i=1, size(this%a)
      call Print ('GreensFunc : a z ' //  i // ' ' // this%a(i) // ' ' // this%z(i), PRINT_VERBOSE, out)
    end do
  endif

end subroutine ApproxFermi_Print

subroutine ApproxFermi_Initialise(this, Fermi_E, Fermi_T, band_width)
  type(ApproxFermi), intent(inout) :: this
  real(dp), intent(in) :: Fermi_E, Fermi_T, band_width

  call Finalise(this)

  call calc_poles_from_T(this, Fermi_E, Fermi_T, band_width)

end subroutine

subroutine ApproxFermi_Initialise_lowlevel(this, n_poles, Fermi_E, band_width)
  type(ApproxFermi), intent(inout) :: this
  integer, intent(in) :: n_poles
  real(dp), intent(in) :: Fermi_E, band_width

  call calc_poles(this, Fermi_E, n_poles, band_width)

end subroutine

subroutine calc_poles(this, Fermi_E, n_poles, band_width)
  type(ApproxFermi), intent(inout) :: this
  real(dp), intent(in) :: Fermi_E, band_width
  integer, intent(in) :: n_poles

  if (n_poles == 0 .or. (band_width .feq. 0.0_dp)) then
    call system_abort("Can't calc_poles with n_poles or band_width = 0")
  endif

  this%Fermi_E = Fermi_E
  this%n_poles = n_poles
  this%band_width = band_width

  if (allocated(this%a)) then
    if (size(this%a) /= this%n_poles) deallocate(this%a)
  endif
  if (allocated(this%z)) then
    if (size(this%z) /= this%n_poles) deallocate(this%z)
  endif

  if (.not.allocated(this%z)) allocate(this%z(this%n_poles))
  if (.not.allocated(this%a)) allocate(this%a(this%n_poles))

  call calc_pole_values(this)

end subroutine calc_poles

subroutine calc_poles_from_T(this, Fermi_E, Fermi_T, band_width)
  type(ApproxFermi), intent(inout) :: this
  real(dp), intent(in) :: Fermi_E, Fermi_T, band_width

  if (band_width .feq. 0.0_dp) then
    call system_abort("Can't calc_poles_from_T with band_width = 0")
  endif

  this%n_poles = guess_n_poles(Fermi_T, band_width, this%band_width)

  call calc_poles(this, Fermi_E, this%n_poles, this%band_width)

end subroutine calc_poles_from_T

function guess_n_poles(Fermi_T, band_width_in, band_width)
  real(dp), intent(in) :: Fermi_T, band_width_in
  real(dp), intent(out) :: band_width
  integer :: guess_n_poles

  type(ApproxFermi) t_AF
  real(dp) :: Fermi_slope, AF_slope
  integer np

  Fermi_slope = f_Fermi_deriv(0.0_dp, Fermi_T, 0.0_dp)

  do np=4, 1000, 4
    call Initialise(t_AF, np, 0.0_dp, band_width_in)
    AF_slope = approx_f_Fermi_deriv(t_AF, 0.0_dp)
    if (AF_slope < Fermi_slope) then
      guess_n_poles = np
      band_width = band_width_in*(AF_slope/Fermi_slope)
      return
    endif
    call Finalise(t_AF)
  end do

  call system_abort ("Couldn't find appropriate band_width and n_poles in guess_n_poles")

end function guess_n_poles

elemental function approx_f_Fermi(this, E)
  type(ApproxFermi), intent(in) :: this
  real(dp), intent(in) :: E
  real(dp) :: approx_f_Fermi

  integer i

  approx_f_Fermi = 0.0_dp
  do i=1, this%n_poles
    approx_f_Fermi = approx_f_Fermi + this%a(i)/(E-this%z(i))
  end do
  
end function approx_f_Fermi

elemental function approx_f_Fermi_deriv(this, E)
  type(ApproxFermi), intent(in) :: this
  real(dp), intent(in) :: E
  real(dp) :: approx_f_Fermi_deriv

  integer i

  approx_f_Fermi_deriv = 0.0_dp
  do i=1, this%n_poles
    approx_f_Fermi_deriv = approx_f_Fermi_deriv - this%a(i)/((E-this%z(i))**2)
  end do

end function approx_f_Fermi_deriv

subroutine calc_pole_values(this)
  type(ApproxFermi), intent(inout) :: this

  integer :: n_poles
  real(dp) :: n_poles_d

  complex(dp) :: a, b, c, aa, bb, cc, x
  integer i, j
  real(dp) i_d

  real(dp) :: gamma, x_bot, scale, shift
  ! real(dp) :: f_max

  n_poles = this%n_poles * 2
  n_poles_d = n_poles

  gamma = 3.0_dp-sqrt(8.0_dp)
  ! f_max = 1.0_dp/( ((4.0_dp*((1.0_dp-gamma)**-2) - 1.0_dp)**(n_poles_d/2)) + 1.0_dp )
  x_bot = -4.0_dp*(2.0_dp*n_poles_d-12.0_dp+6.0_dp*gamma)/((1.0_dp+gamma)**2)

  do i=1, this%n_poles/2
    i_d = i-1
    a = (1.0_dp + gamma) / (2.0_dp*n_poles_d)
    b = (1.0_dp - gamma) / (n_poles_d)
    c = cmplx( cos(PI/(n_poles_d/2.0_dp) + i_d*PI/(n_poles_d/4.0_dp)), &
	       sin(PI/(n_poles_d/2.0_dp) + i_d*PI/(n_poles_d/4.0_dp)), dp )

    aa = a*a
    bb = (2.0_dp*a + b*c)
    cc = (1.0_dp-c)

    this%z(i) = conjg ( (-bb - sqrt(bb*bb - 4.0_dp*aa*cc) ) / (2.0_dp * aa) )
    this%z(n_poles/2-i+1) = (-bb + sqrt(bb*bb - 4.0_dp*aa*cc) ) / (2.0_dp * aa)
  end do

  do i=1, this%n_poles/2
    x = this%z(i)
    this%z(i) = this%z(this%n_poles-i+1)
    this%z(this%n_poles-i+1) = x
  end do

  scale = -this%band_width/x_bot
  shift = this%Fermi_E

  do i=1, this%n_poles
    this%z(i) = this%z(i) * scale + shift
  end do

  do i=1, this%n_poles
    x = this%z(i)

    this%a(i) = (scale**(n_poles/2))*((1.0_dp-((x-shift)/scale)*(1.0_dp-gamma)/n_poles_d)**(n_poles/2))

    do j=1, this%n_poles
      if (j == i) then
	  this%a(i) = scale*(((n_poles_d*2.0_dp)/(1.0_dp+gamma))**2)*this%a(i) /  &
			      (x-conjg(this%z(j)))
      else
	  this%a(i) = scale*(((n_poles_d*2.0_dp)/(1.0_dp+gamma))**2)*this%a(i) / &
			       ( (x-this%z(j))*(x-conjg(this%z(j))) )
      endif
    end do

    this%a(i) = this%a(i) * 2.0_dp
  end do

end subroutine

end module ApproxFermi_module
