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

  implicit none
  private

public :: Multipole_Moments_Site_Site_Interaction, assignment(=)

public :: Multipole_Interactions_Site
type Multipole_Interactions_Site
  integer :: pos_type
  real(dp) :: charge = 0.0_dp
  real(dp), dimension(3) :: position, dipole, e_grad_pos
  real(dp), dimension(:), allocatable :: e_grad_moment
  real(dp), dimension(3,3) :: quadrupole
  real(dp), dimension(:,:,:), allocatable :: moment_grad_positions, pos_grad_positions ! are derivatives of the multipole position and components with respect to atomic positions
  integer :: atomic_number, d
  logical :: initialised
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
  to%charge = from%charge
  to%d=from%d
  to%atomic_number = from%atomic_number
  to%position = from%position
  to%dipole=from%dipole
  to%quadrupole=from%quadrupole
  to%e_grad_pos = from%e_grad_pos

  if (allocated(from%e_grad_moment)) then
    allocate(to%e_grad_moment(size(from%e_grad_moment)))
    to%e_grad_moment=from%e_grad_moment
  end if

  if (allocated(from%moment_grad_positions)) then
    allocate(to%moment_grad_positions(size(from%moment_grad_positions,1),size(from%moment_grad_positions,2),size(from%moment_grad_positions,3)))
    to%moment_grad_positions=from%moment_grad_positions
  end if

  if (allocated(from%pos_grad_positions)) then
    allocate(to%pos_grad_positions(size(from%pos_grad_positions,1),size(from%pos_grad_positions,2),size(from%pos_grad_positions,3)))
    to%pos_grad_positions=from%pos_grad_positions
  end if

  to%initialised=.true.

end subroutine Multipole_Site_Assignment

subroutine Multipole_Site_Finalise(this)

  type(Multipole_Interactions_Site), intent(inout) :: this

  if (allocated(this%e_grad_moment)) deallocate(this%e_grad_moment)
  if (allocated(this%moment_grad_positions)) deallocate(this%moment_grad_positions)
  if (allocated(this%pos_grad_positions)) deallocate(this%pos_grad_positions)

  this%atomic_number=0
  this%d=0
  this%charge = 0.0_dp
  this%e_grad_pos = 0.0_dp
  this%pos_type = 0
  this%initialised=.false.

end subroutine Multipole_Site_Finalise


recursive subroutine Multipole_Moments_Site_Site_Interaction(energy,site_one, site_two, do_f, error ,test)
  ! calculate mutual forces and torques due to multipole moments at two sites
  real(dp), intent(out) :: energy
  real(dp) :: e0, e_plus, step=1.0D-8, q
  type(Multipole_Interactions_Site), intent(inout) :: site_one, site_two
  logical :: do_f, do_test
  logical, optional,intent(in) :: test
  integer, optional, intent(out) :: error
  integer :: i
  real(dp), dimension(3) :: pos, dip

  do_test=optional_default(.false.,test)

  if (site_one%d .eq. 1) then
    if (site_two%d .eq. 1) then ! both are charges
      call Multipole_Interactions_Charge_Charge(energy,site_one, site_two, do_f)
    end if
    if (site_two%d .eq. 3) then ! site one is a charge and site two is a dipole
      call Multipole_Interactions_Charge_Dipole(energy,site_one, site_two, do_f)
    end if
  else if (site_one%d .eq. 3) then
    if (site_two%d .eq. 1) then ! site one is a dipole and site two is a charge
      call Multipole_Interactions_Charge_Dipole(energy,site_two, site_one, do_f)
    end if
    if (site_two%d .eq. 3) then ! both are dipoles
      call Multipole_Interactions_Dipole_Dipole(energy,site_one, site_two, do_f)
    end if
  end if

  if (do_test) then
    call print("testing gradients of site-site interactions")
    e0=energy
    pos=site_one%position
    do i=1,3
      site_one%position(i) = pos(i) + step
      call Multipole_Moments_Site_Site_Interaction(e_plus,site_one, site_two,.false.)
      call print("step size : "// step//" gradient : "//site_one%e_grad_pos(i)//" e diff : "//(e_plus-e0)//" fd grad : " // (e_plus-e0)/(step) )
      site_one%position(i) = pos(i)
    end do
    pos=site_two%position
    do i=1,3
      site_two%position(i) = pos(i) + step
      call Multipole_Moments_Site_Site_Interaction(e_plus,site_one, site_two,.false.)
      call print("step size : "// step//" gradient : "//site_two%e_grad_pos(i)//" e diff : "//(e_plus-e0)//" (e-e0) / (grad*step) : " // (e_plus-e0)/(site_two%e_grad_pos(i)*step) )
      site_two%position(i) = pos(i)
    end do
  end if

end subroutine Multipole_Moments_Site_Site_Interaction


subroutine Multipole_Interactions_Charge_Charge(energy,site_one, site_two, do_f)
  type(Multipole_Interactions_Site), intent(inout) :: site_one, site_two
  logical,intent(in) :: do_f
  real(dp),intent(out) :: energy
  real(dp), dimension(3) :: r_one_two, e_grad_pos
  real(dp) :: r_norm

  ! vector pointing from site one to site two. position of second site should already be set according to minimum image convention
  r_one_two = site_two%position - site_one%position
  r_one_two = r_one_two  / BOHR
  r_norm = norm(r_one_two)  

  energy = site_one%charge * site_two%charge / r_norm

  ! force
  if (do_f) then
    e_grad_pos =  (energy / (r_norm*r_norm) ) * r_one_two 
    e_grad_pos = e_grad_pos * (HARTREE / BOHR) !  convert from Hartree/Bohr to eV/Angstrom
    site_one%e_grad_pos = e_grad_pos
    site_two%e_grad_pos = -e_grad_pos
  end if

  energy = HARTREE * energy ! convert to eV

  if (do_f) then
    site_one%e_grad_moment = energy / site_one%charge
    site_two%e_grad_moment = energy / site_two%charge
  end if

  return

end subroutine Multipole_Interactions_Charge_Charge



subroutine Multipole_Interactions_Charge_Dipole(energy,site_one, site_two, do_f)
  type(Multipole_Interactions_Site), intent(inout) :: site_one, site_two
  logical,intent(in) :: do_f
  real(dp),intent(out) :: energy
  real(dp), dimension(3) :: r_one_two, dipole, r_hat
  real(dp) :: r_norm, energy_ha, r_cubed

  ! site one is the charge, site two is the dipole

  ! Convert everything to atomic units: Hartrees, e*BOHR, BOHR
  r_one_two = (site_two%position - site_one%position)*(1.0_dp/BOHR)
  r_norm = norm(r_one_two)  
  r_cubed = r_norm*r_norm*r_norm
  r_hat = r_one_two * (1.0_dp /r_norm)

  dipole = site_two%dipole *(1.0_dp/BOHR) 

  energy_ha = -site_one%charge * dot_product(dipole,r_one_two) / (r_cubed) ! this in atomic units: 
  energy = energy_ha * HARTREE ! convert energy to eV

  if (do_f) then
    site_one%e_grad_moment = energy / site_one%charge ! energy already in eV and charge just in units of e
    site_one%e_grad_pos = ( site_one%charge / (r_cubed) ) * (dipole - (3.0_dp* dot_product(dipole,r_hat) * r_hat) ) ! this is in Hartrees per Bohr
    site_one%e_grad_pos = site_one%e_grad_pos * (HARTREE/BOHR)

    site_two%e_grad_moment = - ( site_one%charge / (r_cubed) )* r_one_two
    site_two%e_grad_moment = site_two%e_grad_moment * (HARTREE/BOHR) 
    site_two%e_grad_pos = - site_one%e_grad_pos
  end if


  return

end subroutine Multipole_Interactions_Charge_Dipole


subroutine Multipole_Interactions_Dipole_Dipole(energy,site_one, site_two, do_f)
  type(Multipole_Interactions_Site), intent(inout) :: site_one, site_two
  logical,intent(in) :: do_f
  real(dp),intent(out) :: energy
  real(dp), dimension(3) :: r_one_two, dip1,dip2, r_hat
  real(dp) :: r_norm, energy_ha, r_cubed, r_power_four, dip1_dot_rhat, dip2_dot_rhat, dip1_dot_dip2

  ! vector pointing from site one to site two
  r_one_two = site_two%position - site_one%position
  r_one_two = r_one_two*(1.0_dp/BOHR)
  r_norm = norm(r_one_two)  
  r_hat = r_one_two * (1.0_dp/r_norm)
  r_cubed = r_norm*r_norm*r_norm
  r_power_four = r_cubed*r_norm

  ! convert dipole moments from e * Angstrom to e * Bohr (atomic units)
  dip1 = site_one%dipole *(1.0_dp/BOHR)
  dip2 = site_two%dipole *(1.0_dp/BOHR)

  ! pre-calculate dot products
  dip1_dot_rhat =  dot_product(dip1,r_hat) 
  dip2_dot_rhat =  dot_product(dip2,r_hat) 
  dip1_dot_dip2 =  dot_product(dip1,dip2) 

  energy_ha = (1.0_dp/r_cubed) * ( dip1_dot_dip2 - (3.0_dp*dip1_dot_rhat*dip2_dot_rhat) )

  energy = energy_ha * HARTREE

  if (do_f) then
    site_one%e_grad_moment = (1.0_dp/r_cubed) * (dip2 - 3.0_dp * dip2_dot_rhat* r_hat )
    site_one%e_grad_moment = site_one%e_grad_moment * (HARTREE/BOHR) 
    site_one%e_grad_pos = (3.0_dp/(r_power_four)) * ( (dip1_dot_dip2 - 5.0_dp*dip1_dot_rhat*dip2_dot_rhat)*r_hat + (dip1_dot_rhat)*dip2 + (dip2_dot_rhat)*dip1  )
    site_one%e_grad_pos = site_one%e_grad_pos * (HARTREE/BOHR)

    site_two%e_grad_moment = (1.0_dp/r_cubed) * (dip1 - 3.0_dp * dip1_dot_rhat* r_hat )
    site_two%e_grad_moment = site_two%e_grad_moment * (HARTREE/BOHR) 
    site_two%e_grad_pos = - site_one%e_grad_pos
  end if

  return

end subroutine Multipole_Interactions_Dipole_Dipole

end module multipole_interactions_module
