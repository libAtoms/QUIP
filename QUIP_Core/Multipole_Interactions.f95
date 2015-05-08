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
  use gamma_module

  implicit none
  private

integer, parameter :: Damping_None = 0
integer, parameter :: Damping_Exp = 1
integer, parameter :: Damping_Erf = 2
integer, parameter :: Damping_Erf_Uniform = 3

integer, parameter :: Screening_None = 0
integer, parameter :: Screening_Yukawa = 1
integer, parameter :: Screening_Erfc_Uniform = 2



public :: Multipole_Moments_Site_Site_Interaction, assignment(=),  T_rank_zero, T_rank_one, T_rank_two, T_rank_three, finalise

public :: Multipole_Calc_Opts
type Multipole_Calc_Opts
  integer :: damping=Damping_None,screening=Screening_None,damp_exp_order
  real(dp) :: erf_kappa_uniform,erfc_kappa_uniform,damp_exp_scale,yukawa_alpha,yukawa_smooth_length
  logical :: do_energy, do_force, do_test, do_field, do_pot
  !______________________________________________________________________________________________________
  ! Only certain dmaping/screening combinations make sense. The interactions can be damped(smeared) at short range and screened at long range
  ! Short range damping can be gaussian or exponential
  !
  ! erf_kappa -     
  !             inverse length parameter for gaussian short-range damping, used in the reciprical part of Ewald
  !             but also in preventing divergence in polarisable models
  ! 
  ! damp_exp_scale, damp_exp_order
  !     
  !             Exponential short-range damping as used in the TTM , FX, and  AMOEBA water potentials, among others.
  !             the damping goes like exp(-a(r/r0)**m) where m=damp_exp_order. Models use m=3 or m=4. a is the damp_exp_scale
  !             (optional, default 1.0), and r0 is related to the damp_rad of each site, see function site_site_params.
  !
  ! yukawa_alpha,smooth_length,cutoff
  !             parameters for yukawa exponential screening of long-range interactions
  !
  ! erfc_kappa -  
  !             inverse length parameter for Ewald long-range screening
  !
  !
  !______________________________________________________________________________________________________

end type Multipole_Calc_Opts

public :: Multipole_Interactions_Site
type Multipole_Interactions_Site
  ! alpha is scalar polarisability
  integer :: pos_type, charge_method, dipole_method, atomic_number, d ! d is the number of multipole components on this site, 1 if just a charge, 3 if dipole, 4 if charge+dipole
  real(dp) :: charge = 0.0_dp, potential=0.0_dp, alpha=0.0_dp , e_grad_charge=0.0_dp , damp_rad=0.0_dp
  real(dp), dimension(3) :: position, dipole, e_grad_pos=(/0.0_dp,0.0_dp,0.0_dp/), e_field=(/0.0_dp,0.0_dp,0.0_dp/), e_grad_dipole=(/0.0_dp,0.0_dp,0.0_dp/)
  ! real(dp), dimension(3,3) :: quadrupole ! no quadrupolar interactions implemented but could be straighforwardly added
  integer, dimension(:), allocatable :: atom_indices
  real(dp), dimension(:,:,:), allocatable :: charge_grad_positions, dipole_grad_positions, pos_grad_positions ! are derivatives of the multipole position and components with respect to atomic positions
  logical :: initialised, polarisable,damped
end type Multipole_Interactions_Site

!% Overloaded assigment operator. 
private :: Multipole_Site_Assignment
interface assignment(=)
   module procedure Multipole_Site_Assignment
end interface assignment(=)

interface finalise
   module procedure Multipole_Site_Finalise
end interface finalise

contains

subroutine Multipole_Site_Assignment(to,from)

  type(Multipole_Interactions_Site), intent(inout) :: to
  type(Multipole_Interactions_Site), intent(in)    :: from

  call Multipole_Site_Finalise(to)

  to%pos_type = from%pos_type
  to%polarisable = from%polarisable
  to%alpha = from%alpha
  to%damped = from%damped
  to%damp_rad = from%damp_rad
  to%charge = from%charge
  to%potential = from%potential
  to%e_field = from%e_field
  to%d=from%d
  to%atomic_number = from%atomic_number
  to%position = from%position
  to%dipole=from%dipole
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
  this%damp_rad = 0.0_dp
  this%potential = 0.0_dp
  this%e_field = 0.0_dp
  this%position = 0.0_dp

  this%initialised=.false.

end subroutine Multipole_Site_Finalise

! calculate pairwise damping coeffs for exponential, gaussian damping
subroutine Site_Site_Params(do_damp,do_screen,damp_rad_pairwise,erf_kappa_pairwise,radius_one,radius_two,calc_opts)
  real (dp), intent(in) :: radius_one,radius_two
  type(Multipole_Calc_Opts) :: calc_opts
  real(dp),intent(out) :: erf_kappa_pairwise,damp_rad_pairwise
  logical,intent(out)::do_damp,do_screen

  do_screen  =  calc_opts%screening /= Screening_None
  do_damp    =  calc_opts%damping /= Damping_None

  if(calc_opts%damping==Damping_Erf_Uniform) then 
    erf_kappa_pairwise = calc_opts%erf_kappa_uniform
  else if (do_damp) then
    damp_rad_pairwise = sqrt(radius_one*radius_two)
    damp_rad_pairwise = damp_rad_pairwise/(calc_opts%damp_exp_scale**(1.0_dp/calc_opts%damp_exp_order))
    erf_kappa_pairwise =  1.0_dp/norm((/radius_one,radius_two/)) !! equations 12,13 in J. Phys.: Condens. Matter 26 (2014) 213202    
  end if

end subroutine Site_Site_Params


!  calculates interactions between multipole moments on two sites.
!  site%d refers to the total number of multipole components on that site. If d=1 there is only a charge, if d=3 there is only a dipole
!  if d=4 then the site has both a charge and a dipole.
!  if multiple sites are present then the energies and forces from the different interactions are added together
recursive subroutine Multipole_Moments_Site_Site_Interaction(energy,site_one,site_two,calc_opts, cutoff,error,test)
  real(dp), intent(out) :: energy
  type(Multipole_Interactions_Site), intent(inout) :: site_one, site_two
  type(Multipole_Calc_Opts) :: calc_opts
  real(dp), optional ::   cutoff
  logical, optional,intent(in) :: test
  integer, optional, intent(out) :: error

  integer :: i
  real(dp) :: e0, e_plus, step=1.0D-8, q
  real(dp), dimension(3) :: pos, dip
  logical, dimension(2) :: has_charge, has_dipole 
  logical :: do_test

  do_test=optional_default(.false.,test)

  if (calc_opts%do_energy)  energy=0.0_dp
!call print("do energy ? "//calc_opts%do_energy)
!call print("site site interaction")
!call print("position 1 "//site_one%position)
!call print("position 2 "//site_two%position)

  if (calc_opts%do_force) then
    site_one%e_grad_pos = 0.0_dp
    site_one%e_grad_charge = 0.0_dp
    site_one%e_grad_dipole = 0.0_dp

    site_two%e_grad_pos = 0.0_dp
    site_two%e_grad_charge = 0.0_dp
    site_two%e_grad_dipole = 0.0_dp
  end if
  if (calc_opts%do_field) then
    site_one%e_field = 0.0_dp
    site_two%e_field = 0.0_dp
  end if
  if (calc_opts%do_pot) then
    site_one%potential = 0.0_dp
    site_two%potential = 0.0_dp
  end if

  has_charge(1) = (site_one%d .eq. 1 .or. site_one%d .eq. 4)
  has_charge(2) = (site_two%d .eq. 1 .or. site_two%d .eq. 4)

  has_dipole(1) = (site_one%d .eq. 3 .or. site_one%d .eq. 4)
  has_dipole(2) = (site_two%d .eq. 3 .or. site_two%d .eq. 4)

  if (has_charge(1)) then
    if (has_charge(2)) then                      
      call Multipole_Interactions_Charge_Charge(energy,site_one, site_two,calc_opts,cutoff=cutoff)
    end if
    if (has_dipole(2)) then 
      call Multipole_Interactions_Charge_Dipole(energy,site_one, site_two,calc_opts,cutoff=cutoff)
    end if
  end if
  if (has_dipole(1)) then
    if (has_charge(2)) then 
      call Multipole_Interactions_Charge_Dipole(energy,site_two, site_one,calc_opts,cutoff=cutoff)
    end if
    if (has_dipole(2)) then 
      call Multipole_Interactions_Dipole_Dipole(energy,site_one, site_two,calc_opts,cutoff=cutoff)
    end if
  end if
  ! if a site has both a charge and a dipole we iwll have double counted the pot and field
  if (site_one%d == 4) then
      site_one%e_field = 0.5_dp * site_one%e_field 
      site_one%potential = 0.5_dp * site_one%potential
  end if
  if (site_two%d == 4) then
      site_two%e_field = 0.5_dp * site_two%e_field 
      site_two%potential = 0.5_dp * site_two%potential
  end if


!  call print("site site interaction energy is "//energy)


  if (do_test) then
    call print("testing gradients of site-site interactions with step size "//step)
!    call print(" diff vec : "//site_two%position-site_one%position)
    e0=energy
    pos=site_one%position
    do i=1,3
      site_one%position(i) = pos(i) + step
      call Multipole_Moments_Site_Site_Interaction(e_plus,site_one, site_two,calc_opts,cutoff=cutoff,test=.false.)
      call print("r_"//i//" : "//(site_two%position(i)-site_one%position(i))//" gradient : "//site_one%e_grad_pos(i)//" e diff : "//(e_plus-e0)//" fd grad : " // (e_plus-e0)/(step) )
      site_one%position(i) = pos(i)
    end do
    pos=site_two%position
    do i=1,3
      site_two%position(i) = pos(i) + step
      call Multipole_Moments_Site_Site_Interaction(e_plus,site_one, site_two,calc_opts,cutoff=cutoff,test=.false.)
      call print("r_"//i//" : "//(site_two%position(i)-site_one%position(i))//" gradient : "//site_two%e_grad_pos(i)//" e diff : "//(e_plus-e0)//" fd grad : " // (e_plus-e0)/(step) )
      site_two%position(i) = pos(i)
    end do
  end if

end subroutine Multipole_Moments_Site_Site_Interaction


subroutine Multipole_Interactions_Charge_Charge(energy,site_one, site_two,calc_opts,cutoff,error)
  type(Multipole_Interactions_Site), intent(inout) :: site_one, site_two
  type(Multipole_Calc_Opts) :: calc_opts
  real(dp),intent(out) :: energy
  real(dp), optional :: cutoff
  integer, intent(out), optional :: error
  real(dp), dimension(3) :: r_one_two, e_grad_pos, T1
  real(dp) :: T0

  ! vector pointing from site one to site two. position of second site should already be set to correct image position
  r_one_two =  site_two%position - site_one%position  


  T0 = T_rank_zero(r_one_two,calc_opts,site_one%damp_rad,site_two%damp_rad,cutoff=cutoff)
  T1 = T_rank_one(r_one_two,calc_opts,site_one%damp_rad,site_two%damp_rad,cutoff=cutoff)
  if (calc_opts%do_energy) energy = energy + site_one%charge * T0 * site_two%charge  

  ! force
  if (calc_opts%do_force) then
    e_grad_pos =  site_one%charge * T1 * site_two%charge ! NB this should be negative of the force
    site_one%e_grad_pos    = site_one%e_grad_pos    - e_grad_pos
    site_two%e_grad_pos    = site_two%e_grad_pos    + e_grad_pos

    site_one%e_grad_charge = site_one%e_grad_charge + T0 * site_two%charge
    site_two%e_grad_charge = site_two%e_grad_charge + T0 * site_one%charge
  end if

  if (calc_opts%do_pot) then
    site_one%potential = site_one%potential + T0 * site_two%charge
    site_two%potential = site_two%potential + T0 * site_one%charge
  end if

  if (calc_opts%do_field) then
    site_one%e_field = site_one%e_field + T1 * site_two%charge
    site_two%e_field = site_two%e_field - T1 * site_one%charge
  end if

  return

end subroutine Multipole_Interactions_Charge_Charge

subroutine Multipole_Interactions_Charge_Dipole(energy,site_one, site_two,calc_opts,cutoff,error)
  type(Multipole_Interactions_Site), intent(inout) :: site_one, site_two
  type(Multipole_Calc_Opts) :: calc_opts
  real(dp),intent(out) :: energy
  real(dp), optional :: cutoff
  integer, intent(out), optional :: error
  real(dp), dimension(3) :: r_one_two, e_grad_pos, T1
  real(dp), dimension(3,3) :: T2
  real(dp) :: T0 

  ! site one is the charge, site two is the dipole

  ! vector pointing from site one to site two. position of second site should already be set to correct image position
  r_one_two =  site_two%position - site_one%position  

  T0 = T_rank_zero(r_one_two,calc_opts,site_one%damp_rad,site_two%damp_rad,cutoff=cutoff)
  T1 = T_rank_one(r_one_two,calc_opts,site_one%damp_rad,site_two%damp_rad,cutoff=cutoff)
  T2 = T_rank_two(r_one_two,calc_opts,site_one%damp_rad,site_two%damp_rad,cutoff=cutoff)

  if (calc_opts%do_energy) energy = energy + site_one%charge * dot_product(T1, site_two%dipole)  

  ! force
  if (calc_opts%do_force) then
    e_grad_pos = - site_one%charge * matmul(T2,site_two%dipole) !NB grad of energy so minus the force
    site_one%e_grad_pos    = site_one%e_grad_pos     + e_grad_pos
    site_two%e_grad_pos    = site_two%e_grad_pos     - e_grad_pos

    site_one%e_grad_charge = site_one%e_grad_charge  + dot_product(T1, site_two%dipole)
    site_two%e_grad_dipole = site_two%e_grad_dipole  + site_one%charge * T1
  end if

  if (calc_opts%do_pot) then
    site_one%potential = site_one%potential + dot_product(T1, site_two%dipole) 
    site_two%potential = site_two%potential + site_one%charge * T0
  end if

  if (calc_opts%do_field) then
    site_one%e_field = site_one%e_field + matmul(T2,site_two%dipole)
    site_two%e_field = site_two%e_field - T1 *site_one%charge
  end if

  return

end subroutine Multipole_Interactions_Charge_Dipole

subroutine Multipole_Interactions_Dipole_Dipole(energy,site_one, site_two,calc_opts,cutoff,error)
  type(Multipole_Interactions_Site), intent(inout) :: site_one, site_two
  type(Multipole_Calc_Opts) :: calc_opts
  real(dp),intent(out) :: energy
  real(dp), optional :: cutoff
  integer, intent(out), optional :: error
  real(dp), dimension(3) :: r_one_two, e_grad_pos, temp1, T1
  real(dp), dimension(3,3) :: T2, temp2
  real(dp), dimension(3,3,3) :: T3
  integer :: i,j,k

  ! vector pointing from site one to site two. position of second site should already be set to correct image position
  r_one_two =  site_two%position - site_one%position  

  T1 = T_rank_one(r_one_two,calc_opts,site_one%damp_rad,site_two%damp_rad,cutoff=cutoff)
  T2 = T_rank_two(r_one_two,calc_opts,site_one%damp_rad,site_two%damp_rad,cutoff=cutoff)
  T3 = T_rank_three(r_one_two,calc_opts,site_one%damp_rad,site_two%damp_rad,cutoff=cutoff)

  temp1 = matmul(T2,site_two%dipole)

  if (calc_opts%do_energy) energy = energy - dot_product(site_one%dipole, temp1) 

  ! force
  if (calc_opts%do_force) then
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

  if (calc_opts%do_pot) then
    site_one%potential = site_one%potential - dot_product(T1, site_two%dipole) 
    site_two%potential = site_two%potential + dot_product(T1, site_one%dipole) 
  end if

  if (calc_opts%do_field) then
    site_one%e_field = site_one%e_field + matmul(T2,site_two%dipole)
    site_two%e_field = site_two%e_field + matmul(T2,site_one%dipole)
  end if

  return

end subroutine Multipole_Interactions_Dipole_Dipole


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!______________________________________________________________________________________________________
!
! T tensors with damping/screening using notation from J Chem Phys 133, 234101
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

function T_rank_zero(diff, calc_opts,radius_one,radius_two,cutoff,error)
  ! generalised 1/r including damping and screening
  real(dp), dimension(3),intent(in) :: diff
  type(Multipole_Calc_Opts) :: calc_opts
  real(dp),intent(in) :: radius_one,radius_two
  real(dp), intent(in), optional ::  cutoff
  integer, optional :: error

  real(dp) :: r, T_rank_zero, a, f, s0, s0_damp, s0_screen,u,um
  real(dp) :: damp_rad_exp,erf_kappa,erfc_kappa
  logical :: do_screen, do_damp
  integer :: m

  call Site_Site_Params(do_screen,do_damp,damp_rad_exp,erf_kappa,radius_one,radius_two,calc_opts)
  r=norm(diff)

  if (do_damp .and. do_screen) then
    s0 = -1.0_dp
  else if (do_damp .or. do_screen) then 
    s0 = 0.0_dp
  else
    s0 = 1.0_dp
  end if

  s0_damp=0.0_dp
  s0_screen=0.0_dp

  if (calc_opts%damping == Damping_Erf .or. calc_opts%damping == Damping_Erf_Uniform) then
    s0_damp=erf(erf_kappa*r)
  end if
  if (calc_opts%damping == Damping_Exp ) then 
    m=calc_opts%damp_exp_order
    u=(r/damp_rad_exp)
    um=u**m
    s0_damp=1.0_dp-exp(-um)+u*gamma_incomplete_upper((1.0_dp-1.0_dp/m),um)
  end if

  if (calc_opts%screening == Screening_Erfc_Uniform) then
    erfc_kappa=calc_opts%erfc_kappa_uniform
    s0_screen= erfc(erfc_kappa*r)
  end if
  if (calc_opts%screening == Screening_Yukawa) then
    a = calc_opts%yukawa_alpha
    f = poly_switch(r,cutoff,calc_opts%yukawa_smooth_length)
    s0_screen= exp(-a*r)*f
  end if
  s0 = s0 + s0_damp + s0_screen
  T_rank_zero = (HARTREE*BOHR)*s0/r

end function T_rank_zero

function T_rank_one(diff, calc_opts,radius_one,radius_two,cutoff,error)
  ! generalised grad(1/r) including damping and screening
  real(dp), dimension(3),intent(in) :: diff
  type(Multipole_Calc_Opts) :: calc_opts
  real(dp),intent(in) :: radius_one,radius_two
  real(dp), intent(in), optional ::  cutoff
  integer, optional :: error
  real(dp) :: damp_rad_exp,erf_kappa,erfc_kappa
  real(dp) :: r,r3, a,f,df, u,um,s0_damp, s0_screen, s1,s1_damp,s1_screen
  real(dp), dimension(3) :: T_rank_one
  logical :: do_damp, do_screen
  integer:: m
  r=norm(diff)
  r3=r*r*r

  call Site_Site_Params(do_screen,do_damp,damp_rad_exp,erf_kappa,radius_one,radius_two,calc_opts)

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

  if (calc_opts%damping == Damping_Erf .or. calc_opts%damping == Damping_Erf_Uniform) then
    s0_damp=erf(erf_kappa*r)
    s1_damp = s0_damp - (2.0_dp*r*erf_kappa/sqrt(PI)) * exp(-(erf_kappa*r)**2)
  end if
  if (calc_opts%damping == Damping_Exp ) then 
    m=calc_opts%damp_exp_order
    u=(r/damp_rad_exp)
    um=u**m
    s1_damp=1.0_dp-exp(-um)
  end if

  if (calc_opts%screening == Screening_Erfc_Uniform) then
    erfc_kappa=calc_opts%erfc_kappa_uniform
    s0_screen= erfc(erfc_kappa*r)
    s1_screen = s0_screen + (2.0_dp*r*erfc_kappa/sqrt(PI)) * exp(-(erfc_kappa*r)**2)
  end if
  if (calc_opts%screening == Screening_Yukawa) then
    a = calc_opts%yukawa_alpha
    f = poly_switch(r,cutoff,calc_opts%yukawa_smooth_length)
    df = dpoly_switch(r,cutoff,calc_opts%yukawa_smooth_length)
    s0_screen = exp(-a*r)*f
    s1_screen= s0_screen - r*exp(-a*r)*(df - a*f)
  end if

  s1 = s1+s1_damp+s1_screen
  T_rank_one = -(HARTREE*BOHR)*s1/r3 * diff
end function T_rank_one

function T_rank_two(diff, calc_opts,radius_one,radius_two,cutoff,error)
  real(dp), dimension(3),intent(in) :: diff
  type(Multipole_Calc_Opts) :: calc_opts
  real(dp),intent(in) :: radius_one,radius_two
  real(dp), intent(in), optional ::  cutoff
  integer, optional :: error
  real(dp) :: damp_rad_exp,erf_kappa,erfc_kappa
  real(dp) :: r,r3,r5, a, f,df,d2f,u,um, &
    s0_damp=0.0_dp, s0_screen=0.0_dp, s1=0.0_dp,s1_damp=0.0_dp,s1_screen=0.0_dp,s2=0.0_dp,s2_damp=0.0_dp,s2_screen=0.0_dp
  real(dp), dimension(3,3) :: T_rank_two, identity
  logical :: do_damp, do_screen
  integer :: m

  r=norm(diff)
  r3=r*r*r
  r5=r3*r*r

  call Site_Site_Params(do_screen,do_damp,damp_rad_exp,erf_kappa,radius_one,radius_two,calc_opts)

  identity=0.0_dp
  call add_identity(identity)

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

  if (calc_opts%damping == Damping_Erf .or. calc_opts%damping == Damping_Erf_Uniform) then
    s0_damp=erf(erf_kappa*r)
    s1_damp = s0_damp - (2.0_dp*r*erf_kappa/sqrt(PI)) * exp(-(erf_kappa*r)**2)
    s2_damp = s1_damp - (4.0_dp/(3.0_dp*sqrt(PI)))*((r*erf_kappa)**3) * exp(-(erf_kappa*r)**2)
  end if
  if (calc_opts%damping == Damping_Exp ) then 
    m=calc_opts%damp_exp_order
    u=(r/damp_rad_exp)
    um=u**m
    s1_damp=1.0_dp-exp(-um)
    s2_damp=1.0_dp-(1.0_dp+(m/3.0_dp)*um)*exp(-um)
  end if
  if (calc_opts%screening == Screening_Erfc_Uniform) then
    erfc_kappa=calc_opts%erfc_kappa_uniform
    s0_screen= erfc(erfc_kappa*r)
    s1_screen = s0_screen + (2.0_dp*r*erfc_kappa/sqrt(PI)) * exp(-(erfc_kappa*r)**2)
    s2_screen = s1_screen + (4.0_dp/(3.0_dp*sqrt(PI)))*((r*erfc_kappa)**3)* exp(-(erfc_kappa*r)**2)
  end if
  if (calc_opts%screening == Screening_Yukawa) then
    a = calc_opts%yukawa_alpha
    f = poly_switch(r,cutoff,calc_opts%yukawa_smooth_length)
    df = dpoly_switch(r,cutoff,calc_opts%yukawa_smooth_length)
    d2f = d2poly_switch(r,cutoff,calc_opts%yukawa_smooth_length)
    s0_screen = exp(-a*r)*f
    s1_screen= s0_screen - r*exp(-a*r)*(df - a*f)
    s2_screen = s1_screen + (r*r/3.0_dp)* exp(-a*r) *( (a**2)*f - 2.0_dp*a*df + d2f )
  end if

  s1 = s1+s1_damp+s1_screen
  s2 = s2+s2_damp+s2_screen

  T_rank_two = (3.0_dp*s2/r5) * (diff .outer. diff) - (s1/r3) * identity
  T_rank_two =   (HARTREE*BOHR)*T_rank_two
end function T_rank_two


function T_rank_three(diff,calc_opts,radius_one,radius_two,cutoff,error)
  real(dp), dimension(3),intent(in) :: diff
  type(Multipole_Calc_Opts) :: calc_opts
  real(dp),intent(in) :: radius_one,radius_two
  real(dp), intent(in), optional ::  cutoff
  integer, optional :: error
  real(dp) :: damp_rad_exp, erf_kappa, erfc_kappa
  real(dp) :: r,r5,r7,c1, c2, a, f, df, d2f, d3f, u,um,&
                    s0_damp, s0_screen,s1_damp,s1_screen, &
                                       s2_damp,s2_screen, s2,  &
                                       s3_damp,s3_screen, s3
  real(dp), dimension(3,3,3) :: T_rank_three
  logical :: do_damp, do_screen
  integer :: i,j,k,m

  r=norm(diff)
  r5=r**5
  r7=r5*r*r

  call Site_Site_Params(do_screen,do_damp,damp_rad_exp,erf_kappa,radius_one,radius_two,calc_opts)

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

  if (calc_opts%damping == Damping_Erf .or. calc_opts%damping == Damping_Erf_Uniform) then
    s0_damp=erf(erf_kappa*r)
    s1_damp = s0_damp - (2.0_dp/sqrt(PI))*(r*erf_kappa) * exp(-(erf_kappa*r)**2)
    s2_damp = s1_damp - (4.0_dp/(3.0_dp*sqrt(PI)))*((r*erf_kappa)**3) * exp(-(erf_kappa*r)**2)
    s3_damp = s2_damp - (8.0_dp/(15.0_dp*sqrt(PI)))*((r*erf_kappa)**5) * exp(-(erf_kappa*r)**2)
  end if
  if (calc_opts%damping == Damping_Exp ) then 
    m=calc_opts%damp_exp_order
    u=(r/damp_rad_exp)
    um=u**m
    s1_damp=1.0_dp-exp(-um)
    s2_damp=1.0_dp-(1.0_dp+(m/3.0_dp)*um)*exp(-um)
    s3_damp=1.0_dp-(1.0_dp+(m*(8.0_dp-m)/15.0)*um+(m**2/15.0_dp)*um*um)*exp(-um)
  end if
  if (calc_opts%screening == Screening_Erfc_Uniform) then
    erfc_kappa=calc_opts%erfc_kappa_uniform
    s0_screen= erfc(erfc_kappa*r)
    s1_screen = s0_screen + (2.0_dp*r*erfc_kappa/sqrt(PI)) * exp(-(erfc_kappa*r)**2)
    s2_screen = s1_screen + (4.0_dp/(3.0_dp*sqrt(PI)))*((r*erfc_kappa)**3)* exp(-(erfc_kappa*r)**2)
    s3_screen = s2_screen + (8.0_dp/(15.0_dp*sqrt(PI)))*((r*erfc_kappa)**5)* exp(-(erfc_kappa*r)**2)
  end if
  if (calc_opts%screening == Screening_Yukawa) then
    a = calc_opts%yukawa_alpha
    f = poly_switch(r,cutoff,calc_opts%yukawa_smooth_length)
    df = dpoly_switch(r,cutoff,calc_opts%yukawa_smooth_length)
    d2f = d2poly_switch(r,cutoff,calc_opts%yukawa_smooth_length)
    d3f = d3poly_switch(r,cutoff,calc_opts%yukawa_smooth_length)

    s0_screen = exp(-a*r)*f
    s1_screen= s0_screen - r*exp(-a*r)*(df - a*f)
    s2_screen = s1_screen + (r*r/3.0_dp)* exp(-a*r) *( (a**2)*f - 2.0_dp*a*df + d2f )
    c1 = (  (a**2)*f*(1.0_dp+a*r) - a*df*(2.0_dp+3.0_dp*a*r) + d2f*(1.0_dp + 3.0_dp*a*r) - d3f*r  )
    s3_screen = s2_screen + (r*r/15.0_dp)* c1 * exp(-a*r)
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


end module multipole_interactions_module
