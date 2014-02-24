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


#include "error.inc"

module multipole_interactions_module

  use system_module          
  use units_module           
  use extendable_str_module  
  use linearalgebra_module   
  use dictionary_module      
  use table_module           
  use periodictable_module   
  use connection_module           
  use atoms_types_module           
  use atoms_module           
  use clusters_module        
  use structures_module      
  use units_module
  use functions_module

  implicit none
  private

public :: Multipole_Moments_Site_Site_Interaction, assignment(=), dipole_propagator, T_rank_zero, T_rank_one, T_rank_two, T_rank_three

public :: Multipole_Interactions_Site
type Multipole_Interactions_Site
  ! alpha is scalar polarisability
  integer :: pos_type, charge_method, dipole_method, atomic_number, d ! d is the number of multipole components on this site, 1 if just a charge, 3 if dipole, 4 if charge+dipole
  real(dp) :: charge = 0.0_dp, potential=0.0_dp, alpha=0.0_dp , e_grad_charge 
  real(dp), dimension(3) :: position, dipole, e_grad_pos, e_field, e_grad_dipole
  real(dp), dimension(3,3) :: quadrupole
  real(dp), dimension(:,:,:), allocatable :: charge_grad_positions, dipole_grad_positions, pos_grad_positions ! are derivatives of the multipole position and components with respect to atomic positions
  logical :: initialised, polarisable
end type Multipole_Interactions_Site

!% Overloaded assigment operator. 
private :: Multipole_Site_Assignment
interface assignment(=)
   module procedure Multipole_Site_Assignment
end interface assignment(=)

contains

subroutine Multipole_Site_Assignment(to,from)

  type(Multipole_Interactions_Site), intent(inout) :: to
  type(Multipole_Interactions_Site), intent(in)    :: from

  call Multipole_Site_Finalise(to)

  to%pos_type = from%pos_type
  to%polarisable = from%polarisable
  to%alpha = from%alpha
  to%charge = from%charge
  to%potential = from%potential
  to%e_field = from%e_field
  to%d=from%d
  to%atomic_number = from%atomic_number
  to%position = from%position
  to%dipole=from%dipole
  to%quadrupole=from%quadrupole
  to%e_grad_pos = from%e_grad_pos
  to%charge_method = from%charge_method
  to%dipole_method = from%dipole_method
  to%e_grad_charge=from%e_grad_charge
  to%e_grad_dipole=from%e_grad_dipole

  if (allocated(from%charge_grad_positions)) then
    allocate(to%charge_grad_positions(size(from%charge_grad_positions,1),size(from%charge_grad_positions,2),size(from%charge_grad_positions,3)))
    to%charge_grad_positions=from%charge_grad_positions
  end if

  if (allocated(from%dipole_grad_positions)) then
    allocate(to%dipole_grad_positions(size(from%dipole_grad_positions,1),size(from%dipole_grad_positions,2),size(from%dipole_grad_positions,3)))
    to%dipole_grad_positions=from%dipole_grad_positions
  end if

  if (allocated(from%pos_grad_positions)) then
    allocate(to%pos_grad_positions(size(from%pos_grad_positions,1),size(from%pos_grad_positions,2),size(from%pos_grad_positions,3)))
    to%pos_grad_positions=from%pos_grad_positions
  end if

  to%initialised=.true.

end subroutine Multipole_Site_Assignment

subroutine Multipole_Site_Finalise(this)

  type(Multipole_Interactions_Site), intent(inout) :: this

  if (allocated(this%charge_grad_positions)) deallocate(this%charge_grad_positions)
  if (allocated(this%dipole_grad_positions)) deallocate(this%dipole_grad_positions)
  if (allocated(this%pos_grad_positions)) deallocate(this%pos_grad_positions)

  this%d=0
  this%charge = 0.0_dp
  this%dipole = 0.0_dp

  this%atomic_number=0
  this%pos_type = 0
  this%charge_method = 0
  this%dipole_method = 0

  this%alpha = 0.0_dp
  this%potential = 0.0_dp
  this%e_field = 0.0_dp
  this%position = 0.0_dp

  this%initialised=.false.

end subroutine Multipole_Site_Finalise

!  calculates interactions between multipole moments on two sites.
!  site%d refers to the total number of multipole components on that site. If d=1 there is only a charge, if d=3 there is only a dipole
!  if d=4 then the site has both a charge and a dipole.
!  if multiple sites are present then the energies and forces from the different interactions are added together
recursive subroutine Multipole_Moments_Site_Site_Interaction(energy,site_one, site_two, do_energy, do_force, do_field, do_pot, &
   erf_kappa, erfc_kappa, yukawa_alpha, cutoff,smoothlength, error,test)
  real(dp), intent(out) :: energy
  type(Multipole_Interactions_Site), intent(inout) :: site_one, site_two
  logical :: do_energy, do_force, do_test, do_field, do_pot
  real(dp), optional :: erf_kappa, erfc_kappa, yukawa_alpha, cutoff, smoothlength
  logical, optional,intent(in) :: test
  integer, optional, intent(out) :: error

  integer :: i
  real(dp) :: e0, e_plus, step=1.0D-8, q
  real(dp), dimension(3) :: pos, dip
  logical, dimension(2) :: has_charge, has_dipole 

  do_test=optional_default(.false.,test)

  if (do_energy)  energy=0.0_dp

  if (do_force) then
    site_one%e_grad_pos = 0.0_dp
    site_one%e_grad_charge = 0.0_dp
    site_one%e_grad_dipole = 0.0_dp

    site_two%e_grad_pos = 0.0_dp
    site_two%e_grad_charge = 0.0_dp
    site_two%e_grad_dipole = 0.0_dp
  end if
  if (do_field) then
    site_one%e_field = 0.0_dp
    site_two%e_field = 0.0_dp
  end if
  if (do_pot) then
    site_one%potential = 0.0_dp
    site_two%potential = 0.0_dp
  end if

  has_charge(1) = (site_one%d .eq. 1 .or. site_one%d .eq. 4)
  has_charge(2) = (site_two%d .eq. 1 .or. site_two%d .eq. 4)

  has_dipole(1) = (site_one%d .eq. 3 .or. site_one%d .eq. 4)
  has_dipole(2) = (site_two%d .eq. 3 .or. site_two%d .eq. 4)

  if (has_charge(1)) then
    if (has_charge(2)) then                      
      call Multipole_Interactions_Charge_Charge(energy,site_one, site_two,do_energy, do_force, do_field, do_pot,erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, &
             yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength)
    end if
    if (has_dipole(2)) then 
      call Multipole_Interactions_Charge_Dipole(energy,site_one, site_two,do_energy, do_force, do_field, do_pot,erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, &
             yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength)
    end if
  end if
  if (has_dipole(1)) then
    if (has_charge(2)) then 
      call Multipole_Interactions_Charge_Dipole(energy,site_two, site_one,do_energy, do_force, do_field, do_pot,erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, &
             yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength)
    end if
    if (has_dipole(2)) then 
      call Multipole_Interactions_Dipole_Dipole(energy,site_one, site_two,do_energy, do_force, do_field, do_pot,erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, &
             yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength)
    end if
  end if
  if (do_field) then ! we will have double counted the field contributions
    if (has_charge(1) .and. has_dipole(1)) then 
      site_one%e_field = 0.5_dp * site_one%e_field 
    end if
    if (has_charge(2) .and. has_dipole(2)) then
      site_two%e_field = 0.5_dp * site_two%e_field 
    end if
  end if
  if (do_pot) then ! we will have double counted the pot contributions
    if (has_charge(1) .and. has_dipole(1)) then 
      site_one%potential = 0.5_dp * site_one%potential
    end if
    if (has_charge(2) .and. has_dipole(2)) then
      site_two%potential = 0.5_dp * site_two%potential 
    end if
  end if


  if (do_test) then
    call print("testing gradients of site-site interactions")
!    call print(" diff vec : "//site_two%position-site_one%position)
    e0=energy
    pos=site_one%position
    do i=1,3
      site_one%position(i) = pos(i) + step
      call Multipole_Moments_Site_Site_Interaction(e_plus,site_one, site_two,.true.,.false.,.false.,.false., erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, &
             yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength,test=.false.)
!      call print("r_"//i//" : "//site_two%position(i)-site_one%position(i)//" gradient : "//site_one%e_grad_pos(i)//" e diff : "//(e_plus-e0)//" fd grad : " // (e_plus-e0)/(step) )
      site_one%position(i) = pos(i)
    end do
    pos=site_two%position
    do i=1,3
      site_two%position(i) = pos(i) + step
      call Multipole_Moments_Site_Site_Interaction(e_plus,site_one, site_two,.true.,.false.,.false.,.false., erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, &
             yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength,test=.false.)
!      call print("r_"//i//" : "//site_two%position(i)-site_one%position(i)//" gradient : "//site_two%e_grad_pos(i)//" e diff : "//(e_plus-e0)//" fd grad : " // (e_plus-e0)/(step) )
      site_two%position(i) = pos(i)
    end do
  end if

end subroutine Multipole_Moments_Site_Site_Interaction


subroutine Multipole_Interactions_Charge_Charge(energy,site_one, site_two, do_energy,do_force, do_field, do_pot, erf_kappa, erfc_kappa, yukawa_alpha,cutoff,smoothlength,error)
  type(Multipole_Interactions_Site), intent(inout) :: site_one, site_two
  logical,intent(in) :: do_energy,do_force, do_field, do_pot
  real(dp),intent(out) :: energy
  real(dp), optional :: erf_kappa, erfc_kappa, yukawa_alpha,cutoff,smoothlength
  integer, intent(out), optional :: error
  real(dp), dimension(3) :: r_one_two, e_grad_pos, T1
  real(dp) :: T0

  ! vector pointing from site one to site two. position of second site should already be set to correct image position
  r_one_two =  site_two%position - site_one%position  

  T0 = T_rank_zero(r_one_two,erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength)
  T1 = T_rank_one(r_one_two,erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength)

  if (do_energy) energy = energy + site_one%charge * T0 * site_two%charge  

  ! force
  if (do_force) then
    e_grad_pos =  site_one%charge * T1 * site_two%charge ! NB this should be negative of the force
    site_one%e_grad_pos    = site_one%e_grad_pos    - e_grad_pos
    site_two%e_grad_pos    = site_two%e_grad_pos    + e_grad_pos

    site_one%e_grad_charge = site_one%e_grad_charge + T0 * site_two%charge
    site_two%e_grad_charge = site_two%e_grad_charge + T0 * site_one%charge
  end if

  if (do_pot) then
    site_one%potential = site_one%potential + T0 * site_two%charge
    site_two%potential = site_two%potential + T0 * site_one%charge
  end if

  if (do_field) then
    site_one%e_field = site_one%e_field + T1 * site_two%charge
    site_two%e_field = site_two%e_field - T1 * site_one%charge
  end if

  return

end subroutine Multipole_Interactions_Charge_Charge

subroutine Multipole_Interactions_Charge_Dipole(energy,site_one, site_two,do_energy, do_force, do_field, do_pot, erf_kappa, erfc_kappa, yukawa_alpha,cutoff,smoothlength,error)
  type(Multipole_Interactions_Site), intent(inout) :: site_one, site_two
  logical,intent(in) :: do_energy, do_force, do_field, do_pot
  real(dp),intent(out) :: energy
  real(dp), optional :: erf_kappa, erfc_kappa, yukawa_alpha,cutoff,smoothlength
  integer, intent(out), optional :: error
  real(dp), dimension(3) :: r_one_two, e_grad_pos, T1
  real(dp), dimension(3,3) :: T2
  real(dp) :: T0 

  ! site one is the charge, site two is the dipole

  ! vector pointing from site one to site two. position of second site should already be set to correct image position
  r_one_two =  site_two%position - site_one%position  

  T0 = T_rank_zero(r_one_two,erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength)
  T1 = T_rank_one(r_one_two,erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength)
  T2 = T_rank_two(r_one_two,erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength)

  if (do_energy) energy = energy + site_one%charge * dot_product(T1, site_two%dipole)  

  ! force
  if (do_force) then
    e_grad_pos = - site_one%charge * matmul(T2,site_two%dipole) !NB grad of energy so minus the force
    site_one%e_grad_pos    = site_one%e_grad_pos     + e_grad_pos
    site_two%e_grad_pos    = site_two%e_grad_pos     - e_grad_pos

    site_one%e_grad_charge = site_one%e_grad_charge  + dot_product(T1, site_two%dipole)
    site_two%e_grad_dipole = site_two%e_grad_dipole  + site_one%charge * T1
  end if

  if (do_pot) then
    site_one%potential = site_one%potential + dot_product(T1, site_two%dipole) 
    site_two%potential = site_two%potential + site_one%charge * T0
  end if

  if (do_field) then
    site_one%e_field = site_one%e_field + matmul(T2,site_two%dipole)
    site_two%e_field = site_two%e_field - T1 *site_one%charge
  end if

  return

end subroutine Multipole_Interactions_Charge_Dipole

subroutine Multipole_Interactions_Dipole_Dipole(energy,site_one, site_two,do_energy, do_force, do_field, do_pot, erf_kappa, erfc_kappa, yukawa_alpha,cutoff,smoothlength,error)
  type(Multipole_Interactions_Site), intent(inout) :: site_one, site_two
  logical,intent(in) :: do_energy, do_force, do_field, do_pot
  real(dp),intent(out) :: energy
  real(dp), optional :: erf_kappa, erfc_kappa, yukawa_alpha,cutoff,smoothlength
  integer, intent(out), optional :: error
  real(dp), dimension(3) :: r_one_two, e_grad_pos, temp1, T1
  real(dp), dimension(3,3) :: T2, temp2
  real(dp), dimension(3,3,3) :: T3
  integer :: i,j,k


  ! vector pointing from site one to site two. position of second site should already be set to correct image position
  r_one_two =  site_two%position - site_one%position  

  T1 = T_rank_one(r_one_two,erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength)
  T2 = T_rank_two(r_one_two,erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength)
  T3 = T_rank_three(r_one_two,erf_kappa=erf_kappa, erfc_kappa=erfc_kappa, yukawa_alpha=yukawa_alpha,cutoff=cutoff,smoothlength=smoothlength)
  temp1 = matmul(T2,site_two%dipole)

  if (do_energy) energy = energy - dot_product(site_one%dipole, temp1) 

  ! force
  if (do_force) then
    do i=1,3
      do j=1,3
        temp1 = T3(i,j,:)
        temp2(i,j) = dot_product(temp1,site_two%dipole)
      end do
    end do
    temp2 = transpose(temp2)
    e_grad_pos =  matmul(temp2,site_one%dipole)
    site_one%e_grad_pos    = site_one%e_grad_pos    + e_grad_pos
    site_two%e_grad_pos    = site_two%e_grad_pos    - e_grad_pos

    site_one%e_grad_dipole = site_one%e_grad_dipole - matmul(T2, site_two%dipole)
    site_two%e_grad_dipole = site_two%e_grad_dipole - matmul(T2, site_one%dipole)
  end if

  if (do_pot) then
    site_one%potential = site_one%potential - dot_product(T1, site_two%dipole) 
    site_two%potential = site_two%potential + dot_product(T1, site_one%dipole) 
  end if

  if (do_field) then
    site_one%e_field = site_one%e_field + matmul(T2,site_two%dipole)
    site_two%e_field = site_two%e_field + matmul(T2,site_one%dipole)
  end if

  return

end subroutine Multipole_Interactions_Dipole_Dipole

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!______________________________________________________________________________________________________
!
! T tensors with damping using notation from J Chem Phys 133, 234101
! note parity of these tensors is (-1)^r where r is the rank of the tensor, minus signs in brackets
! indicate that the sign of 'diff' (the first argument) matters.
!
! Units are as follows : 
! T_rank_zero   :                eV /(electron charge)^2
! T_rank_one    :  Angstrom^-1 * eV /(electron charge)^2
! T_rank_two    :  Angstrom^-2 * eV /(electron charge)^2
! T_rank_three  :  Angstrom^-3 * eV /(electron charge)^2
!
! so that the potential, electric field, and field gradient due to a charge are given by:
! pot = charge * T_rank_zero
! E_field = (-) charge * T_rank_one
! grad(E_field) = charge * T_rank_two
!
! And for a dipole in units of Ang*(electron charge) 
! pot = (-) dot_product(dipole, T_rank_one)
! E_field = matmul(T_rank_two,dipole)
! grad(E_field) = (-) T_rank_three * dipole  (a rank two tensor arising from multiplication of a rank 3 and rank 1 tensor)
!
!______________________________________________________________________________________________________

function T_rank_zero(diff, erf_kappa, erfc_kappa, yukawa_alpha,cutoff,smoothlength,error)
  ! generalised 1/r including damping and screening
  real(dp), dimension(3),intent(in) :: diff
  real(dp), intent(in), optional :: erf_kappa, erfc_kappa, yukawa_alpha,cutoff,smoothlength
  integer, optional :: error

  real(dp) :: r, T_rank_zero, a, f, s0, s0_damp, s0_screen
  logical :: do_screen, do_damp
  integer :: yukawa_count

  r=norm(diff)

  yukawa_count = count((/present(yukawa_alpha),present(cutoff),present(smoothlength)/))
  if (yukawa_count == 3) then
    a = yukawa_alpha
    f = poly_switch(r,cutoff,smoothlength)
  else if (present(yukawa_alpha) .or. present(smoothlength) ) then
    RAISE_ERROR("T_rank_zero was called with some yukawa parameters missing",error)
  end if
  if (present(yukawa_alpha) .and. present(erfc_kappa)) then
    RAISE_ERROR("Simultaneously requesting two types of screening, check erf and erfc haven't been confused",error)
  end if

  do_screen = present(erfc_kappa) .or. present(yukawa_alpha)
  do_damp = present(erf_kappa)

  if (do_damp .and. do_screen) then
    s0 = -1.0_dp
  else if (do_damp .or. do_screen) then 
    s0 = 0.0_dp
  else
    s0 = 1.0_dp
  end if

  s0_damp=0.0_dp
  s0_screen=0.0_dp

  if (do_damp) s0_damp=erf(erf_kappa*r)
  if (do_screen) then
    if (present(erfc_kappa)) then 
      s0_screen= erfc(erfc_kappa*r)
    else
      s0_screen= exp(-a*r)*f
    end if
  end if

  s0 = s0 + s0_damp + s0_screen
  T_rank_zero = (HARTREE*BOHR)*s0/r
end function T_rank_zero

function T_rank_one(diff, erf_kappa, erfc_kappa, yukawa_alpha,cutoff,smoothlength,error)
  ! generalised grad(1/r) including damping and screening
  real(dp), dimension(3),intent(in) :: diff
  real(dp), intent(in), optional :: erf_kappa, erfc_kappa, yukawa_alpha,cutoff,smoothlength
  integer, optional :: error
  real(dp) :: r,r3, a,f,df, s0_damp, s0_screen, s1,s1_damp,s1_screen
  real(dp), dimension(3) :: T_rank_one
  logical :: do_damp, do_screen
  integer :: yukawa_count


  r=norm(diff)
  r3=r*r*r

  yukawa_count = count((/present(yukawa_alpha),present(cutoff),present(smoothlength)/))
  if (yukawa_count == 3) then
    a = yukawa_alpha
    f = poly_switch(r,cutoff,smoothlength)
    df = dpoly_switch(r,cutoff,smoothlength)
  else if (present(yukawa_alpha) .or. present(smoothlength) ) then
    RAISE_ERROR("T_rank_zero was called with some yukawa parameters missing",error)
  end if
  if (present(yukawa_alpha) .and. present(erfc_kappa)) then
    RAISE_ERROR("Simultaneously requesting two types of screening, check erf and erfc haven't been confused",error)
  end if

  do_screen = present(erfc_kappa) .or. present(yukawa_alpha)
  do_damp = present(erf_kappa)

  if (do_damp .and. do_screen) then
    s1 = -1.0_dp
  else if (do_damp .or. do_screen) then 
    s1 = 0.0_dp
  else
    s1 = 1.0_dp
  end if

  s0_damp=0.0_dp
  s0_screen=0.0_dp
  s1_damp=0.0_dp
  s1_screen=0.0_dp

  if (do_damp) then
    s0_damp=erf(erf_kappa*r)
    s1_damp = s0_damp - (2.0_dp*r*erf_kappa/sqrt(PI)) * exp(-(erf_kappa*r)**2)
  end if
  if (do_screen) then
    if (present(erfc_kappa)) then 
      s0_screen= erfc(erfc_kappa*r)
      s1_screen = s0_screen + (2.0_dp*r*erfc_kappa/sqrt(PI)) * exp(-(erfc_kappa*r)**2)
    else
      s0_screen = exp(-a*r)*f
      s1_screen= s0_screen - r*exp(-a*r)*(df - a*f)
    end if
  end if

  s1 = s1+s1_damp+s1_screen
  T_rank_one = -(HARTREE*BOHR)*s1/r3 * diff
end function T_rank_one

function T_rank_two(diff, erf_kappa, erfc_kappa, yukawa_alpha,cutoff,smoothlength,error)
  real(dp), dimension(3),intent(in) :: diff
  real(dp), intent(in), optional :: erf_kappa, erfc_kappa, yukawa_alpha, cutoff,smoothlength
  integer, optional :: error
  real(dp) :: r,r3,r5, a, f,df,d2f, &
    s0_damp=0.0_dp, s0_screen=0.0_dp, s1=0.0_dp,s1_damp=0.0_dp,s1_screen=0.0_dp,s2=0.0_dp,s2_damp=0.0_dp,s2_screen=0.0_dp
  real(dp), dimension(3,3) :: T_rank_two, identity
  logical :: do_damp, do_screen
  integer :: yukawa_count

  r=norm(diff)
  r3=r*r*r
  r5=r3*r*r

  yukawa_count = count((/present(yukawa_alpha),present(cutoff),present(smoothlength)/))
  if (yukawa_count == 3) then
    a = yukawa_alpha
    f = poly_switch(r,cutoff,smoothlength)
    df = dpoly_switch(r,cutoff,smoothlength)
    d2f = d2poly_switch(r,cutoff,smoothlength)
  else if (present(yukawa_alpha) .or. present(smoothlength) ) then
    RAISE_ERROR("T_rank_zero was called with some yukawa parameters missing",error)
  end if
  if (present(yukawa_alpha) .and. present(erfc_kappa)) then
    RAISE_ERROR("Simultaneously requesting two types of screening, check erf and erfc haven't been confused",error)
  end if

  identity=0.0_dp
  call add_identity(identity)

  do_screen = present(erfc_kappa) .or. present(yukawa_alpha)
  do_damp = present(erf_kappa)

  s0_damp=0.0_dp
  s0_screen=0.0_dp
  s1_damp=0.0_dp
  s1_screen=0.0_dp
  s2_damp=0.0_dp
  s2_screen=0.0_dp

  if (do_damp .and. do_screen) then
    s1 = -1.0_dp
    s2 = -1.0_dp
  else if (do_damp .or. do_screen) then 
    s1 = 0.0_dp
    s2 = 0.0_dp
  else
    s1 = 1.0_dp
    s2 = 1.0_dp
  end if

  if (do_damp) then
    s0_damp=erf(erf_kappa*r)
    s1_damp = s0_damp - (2.0_dp*r*erf_kappa/sqrt(PI)) * exp(-(erf_kappa*r)**2)
    s2_damp = s1_damp - (4.0_dp/(3.0_dp*sqrt(PI)))*((r*erf_kappa)**3) * exp(-(erf_kappa*r)**2)
  end if
  if (do_screen) then
    if (present(erfc_kappa)) then 
      s0_screen= erfc(erfc_kappa*r)
      s1_screen = s0_screen + (2.0_dp*r*erfc_kappa/sqrt(PI)) * exp(-(erfc_kappa*r)**2)
      s2_screen = s1_screen + (4.0_dp/(3.0_dp*sqrt(PI)))*((r*erfc_kappa)**3)* exp(-(erfc_kappa*r)**2)
    else
      s0_screen = exp(-a*r)*f
      s1_screen= s0_screen - r*exp(-a*r)*(df - a*f)
      s2_screen = s1_screen + (r*r/3.0_dp)* exp(-a*r) *( (a**2)*f - 2.0_dp*a*df + d2f )
    end if
  end if

  s1 = s1+s1_damp+s1_screen
  s2 = s2+s2_damp+s2_screen

  T_rank_two = (3.0_dp*s2/r5) * (diff .outer. diff) - (s1/r3) * identity
  T_rank_two =   (HARTREE*BOHR)*T_rank_two
end function T_rank_two


function T_rank_three(diff, erf_kappa, erfc_kappa, yukawa_alpha,cutoff,smoothlength,error)
  real(dp), dimension(3),intent(in) :: diff
  real(dp), intent(in), optional :: erf_kappa, erfc_kappa, yukawa_alpha, cutoff,smoothlength
  integer, optional :: error
  real(dp) :: r,r5,r7,c1, c2, a, f, df, d2f, d3f, &
                    s0_damp, s0_screen,s1_damp,s1_screen, &
                                       s2_damp,s2_screen, s2,  &
                                       s3_damp,s3_screen, s3
  real(dp), dimension(3,3,3) :: T_rank_three
  logical :: do_damp, do_screen
  integer :: yukawa_count, i,j,k

  r=norm(diff)
  r5=r**5
  r7=r5*r*r

  yukawa_count = count((/present(yukawa_alpha),present(cutoff),present(smoothlength)/))
  if (yukawa_count == 3) then
    a = yukawa_alpha
    f = poly_switch(r,cutoff,smoothlength)
    df = dpoly_switch(r,cutoff,smoothlength)
    d2f = d2poly_switch(r,cutoff,smoothlength)
    d3f = d3poly_switch(r,cutoff,smoothlength)
  else if (present(yukawa_alpha) .or. present(smoothlength) ) then
    RAISE_ERROR("T_rank_zero was called with some yukawa parameters missing",error)
  end if
  if (present(yukawa_alpha) .and. present(erfc_kappa)) then
    RAISE_ERROR("Simultaneously requesting two types of screening, check erf and erfc haven't been confused",error)
  end if

  do_screen = present(erfc_kappa) .or. present(yukawa_alpha)
  do_damp = present(erf_kappa)

  s0_damp=0.0_dp
  s0_screen=0.0_dp
  s1_damp=0.0_dp
  s1_screen=0.0_dp
  s2_damp=0.0_dp
  s2_screen=0.0_dp
  s3_damp=0.0_dp
  s3_screen=0.0_dp

  if (do_damp .and. do_screen) then
    s2 = -1.0_dp
    s3 = -1.0_dp
  else if (do_damp .or. do_screen) then 
    s2 = 0.0_dp
    s3 = 0.0_dp
  else
    s2 = 1.0_dp
    s3 = 1.0_dp
  end if

  if (do_damp) then
    s0_damp=erf(erf_kappa*r)
    s1_damp = s0_damp - (2.0_dp/sqrt(PI))*(r*erf_kappa) * exp(-(erf_kappa*r)**2)
    s2_damp = s1_damp - (4.0_dp/(3.0_dp*sqrt(PI)))*((r*erf_kappa)**3) * exp(-(erf_kappa*r)**2)
    s3_damp = s2_damp - (8.0_dp/(15.0_dp*sqrt(PI)))*((r*erf_kappa)**5) * exp(-(erf_kappa*r)**2)
  end if
  if (do_screen) then
    if (present(erfc_kappa)) then 
      s0_screen= erfc(erfc_kappa*r)
      s1_screen = s0_screen + (2.0_dp*r*erfc_kappa/sqrt(PI)) * exp(-(erfc_kappa*r)**2)
      s2_screen = s1_screen + (4.0_dp/(3.0_dp*sqrt(PI)))*((r*erfc_kappa)**3)* exp(-(erfc_kappa*r)**2)
      s3_screen = s2_screen + (8.0_dp/(15.0_dp*sqrt(PI)))*((r*erfc_kappa)**5)* exp(-(erfc_kappa*r)**2)
    else
      s0_screen = exp(-a*r)*f
      s1_screen= s0_screen - r*exp(-a*r)*(df - a*f)
      s2_screen = s1_screen + (r*r/3.0_dp)* exp(-a*r) *( (a**2)*f - 2.0_dp*a*df + d2f )
      c1 = (  (a**2)*f*(1.0_dp+a*r) - a*df*(2.0_dp+3.0_dp*a*r) + d2f*(1.0_dp + 3.0_dp*a*r) - d3f*r  )
      s3_screen = s2_screen + (r*r/15.0_dp)* c1 * exp(-a*r)
    end if
  end if

  s2 = s2+s2_damp+s2_screen
  s3 = s3+s3_damp+s3_screen
  c1 = (-s3*15.0_dp/r7)
  c2 = (s2*3.0_dp/r5)
  T_rank_three = 0.0_dp
  do i=1,3
    do j=1,3
      do k=1,3
         T_rank_three(i,j,k) = c1*diff(i)*diff(j)*diff(k)
         if (j .eq. k ) T_rank_three(i,j,k) = T_rank_three(i,j,k) + c2*diff(i)
         if (i .eq. k ) T_rank_three(i,j,k) = T_rank_three(i,j,k) + c2*diff(j)
         if (i .eq. j ) T_rank_three(i,j,k) = T_rank_three(i,j,k) + c2*diff(k)
      end do
    end do
  end do
  T_rank_three =   (HARTREE*BOHR)*T_rank_three

end function T_rank_three


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!______DEPRECATED_____________________________________________________________________________________

! propagator is the T tensor in Stone's notation which generates dipolar electric field 
! at a position 'diff' with respect to the dipole's position.
! the electric field is given by the matrix product : matmul(propagator,dipole)
! diff in Angstroms, dipole in Angstroms * electron charge, field in eV / (Angstrom * electron charge)
function dipole_propagator(diff) result (propagator)
  real(dp), dimension(3),intent(in) :: diff
  real(dp), dimension(3) :: r_hat
  real(dp), dimension(3,3) :: propagator
  real(dp) :: r

  propagator = 0.0_dp
  call add_identity(propagator)
  r=norm(diff)
  r_hat = diff / r

  propagator = (3.0_dp)*(r_hat .outer. r_hat) - propagator
  propagator = (HARTREE*BOHR/(r*r*r))*propagator

end function dipole_propagator

function coordination_function(r,cutoff_in,transition_width)

   real(dp)             :: coordination_function
   real(dp), intent(in) :: r, cutoff_in, transition_width

   if( r > cutoff_in ) then
       coordination_function = 0.0_dp
   elseif( r > (cutoff_in-transition_width) ) then
       coordination_function = 0.5_dp * ( cos(PI*(r-cutoff_in+transition_width)/transition_width) + 1.0_dp )
   else
       coordination_function = 1.0_dp
   endif

endfunction coordination_function

function dcoordination_function(r,cutoff_in,transition_width)

   real(dp)             :: dcoordination_function
   real(dp), intent(in) :: r, cutoff_in,transition_width

   if( r > cutoff_in ) then
       dcoordination_function = 0.0_dp
   elseif( r > (cutoff_in-transition_width) ) then
       dcoordination_function = - 0.5_dp * (PI/ transition_width) * sin(PI*(r-cutoff_in+transition_width)/transition_width)
   else
       dcoordination_function = 0.0_dp
   endif

endfunction dcoordination_function

function d2coordination_function(r,cutoff_in,transition_width)

   real(dp)             :: d2coordination_function
   real(dp), intent(in) :: r, cutoff_in,transition_width

   if( r > cutoff_in ) then
       d2coordination_function = 0.0_dp
   elseif( r > (cutoff_in-transition_width) ) then
       d2coordination_function = - 0.5_dp * ((PI/ transition_width)**2) * cos(PI*(r-cutoff_in+transition_width)/transition_width) 
   else
       d2coordination_function = 0.0_dp
   endif

endfunction d2coordination_function

!NB this third derivative is not continuous! 
function d3coordination_function(r,cutoff_in,transition_width)

   real(dp)             :: d3coordination_function
   real(dp), intent(in) :: r, cutoff_in,transition_width

   if( r > cutoff_in ) then
       d3coordination_function = 0.0_dp
   elseif( r > (cutoff_in-transition_width) ) then
       d3coordination_function = 0.5_dp * ((PI/ transition_width)**3) * sin(PI*(r-cutoff_in+transition_width)/transition_width) 
   else
       d3coordination_function = 0.0_dp
   endif

endfunction d3coordination_function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module multipole_interactions_module
