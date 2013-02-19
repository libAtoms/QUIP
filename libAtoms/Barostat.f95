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
!X  Barostat module
!X
!% This module contains the implementations for all the barostats available
!% in libAtoms: Langevin thermalized Hoover
!%
!% With Q > 0, use adaptive (aka open) Langevin, from Jones \& Leimkuhler 
!% preprint "Adaptive Stochastic Methods for Nonequilibrium Sampling"
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module barostat_module

  use system_module
  use units_module
  use linearalgebra_module
  use atoms_module

  implicit none
  private
  public :: barostat, initialise, finalise, print, update_barostat, set_degrees_of_freedom, barostat_mass
  public :: BAROSTAT_NONE, BAROSTAT_HOOVER_LANGEVIN
  public :: barostat_pre_vel1, barostat_post_vel1_pre_pos, barostat_post_pos_pre_calc, barostat_post_calc_pre_vel2, barostat_post_vel2

  real(dp), parameter :: BAROSTAT_MIN_TEMP = 1.0_dp

  integer, parameter :: BAROSTAT_NONE                 = 0, &
                        BAROSTAT_HOOVER_LANGEVIN      = 1

  !% Langevin barostat - fluctuation/dissipation ---
  !% originally implemented for hydrostatic strain using
  !% EOM from Quigley, D. and Probert, M.I.J., \emph{J. Chem. Phys.}, 
  !% {\bfseries 120} 11432
  !% rediscretized by Noam Bernstein based on ideas from and discussion
  !% with B. Leimkuhler.
  !% modified for arbitrary strain based on (but not exactly using)
  !% formulation in E. Tadmor and R. Miller Modeling Materials: Continuum, Atomistic and Multiscale 
  !% Techniques (Cambridge University Press, 2011). Chap 9.5.  Their 
  !% definition of stress (factors of F, F^T, J, and F^-T) for finite strain is 
  !% optionally used, but certain terms in EOM are excluded, and their discretization 
  !% is entirely ignored.

  type barostat

     integer  :: type  = BAROSTAT_NONE     !% One of the types listed above
     real(dp) :: gamma_epsilon = 0.0_dp   !% Friction coefficient for Langevin part
     real(dp) :: T = -1.0_dp   !% Temperature for Langevin part
     real(dp) :: Ndof  = 0.0_dp   !% The number of degrees of freedom of atoms attached to this barostat
     real(dp) :: stress_ext(3,3) = 0.0_dp       !% External applied stress
     real(dp) :: W_epsilon = 0.0_dp     !% Fictious cell mass in Langevin NPT
     real(dp) :: epsilon_v(3,3) = 0.0_dp !% Velocity of deformation gradient for barostat
     real(dp) :: epsilon_f(3,3) = 0.0_dp !% Force on deformation gradient
     real(dp) :: epsilon_r(3,3) = reshape( (/ 1.0_dp, 0.0_dp, 0.0_dp,    0.0_dp, 1.0_dp, 0.0_dp,    0.0_dp, 0.0_dp, 1.0_dp /), (/ 3, 3 /) ) !% Value of deformation gradient
     logical :: hydrostatic_strain = .true., diagonal_strain = .true., finite_strain_formulation=.false.

     real(dp) :: lattice0(3,3), lattice0_inv(3,3), F_tmp_norm
     logical :: lattice0_initialised = .false.

  end type barostat

  interface initialise
     module procedure barostat_initialise
  end interface initialise

  interface finalise
     module procedure barostat_finalise
  end interface finalise

  interface assignment(=)
     module procedure barostat_assignment
  end interface assignment(=)

  interface print
     module procedure barostat_print
  end interface

  interface update_barostat
     module procedure barostat_update_barostat
  end interface update_barostat

  interface set_degrees_of_freedom
     module procedure set_degrees_of_freedom_int, set_degrees_of_freedom_real
  end interface set_degrees_of_freedom

  interface barostat_mass
     module procedure barostat_mass_int, barostat_mass_real
  end interface barostat_mass

contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X INITIALISE / FINALISE
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine barostat_initialise(this,type,p_ext,stress_ext,hydrostatic_strain,diagonal_strain,finite_strain_formulation,deformation_grad,cell_volume,W_epsilon,Ndof,gamma_epsilon,T,W_epsilon_factor)

    type(barostat),   intent(inout) :: this
    integer,            intent(in)    :: type
    real(dp), optional, intent(in)    :: p_ext
    real(dp), optional, intent(in)    :: stress_ext(3,3)
    logical, optional, intent(in)     :: hydrostatic_strain, diagonal_strain
    logical, optional, intent(in)     :: finite_strain_formulation
    real(dp), optional, intent(in)    :: deformation_grad(3,3)
    real(dp), optional, intent(in)    :: cell_volume
    real(dp), optional, intent(in)    :: W_epsilon
    real(dp), optional, intent(in)    :: Ndof
    real(dp), optional, intent(in)    :: gamma_epsilon
    real(dp), optional, intent(in)    :: T
    real(dp), optional, intent(in)    :: W_epsilon_factor

    real(dp) :: use_T

    if (type /= BAROSTAT_NONE .and. .not.present(p_ext) .and. .not. present(stress_ext)) &
         call system_abort('barostat_initialise: p_ext or stress_ext must be specified when turning on a barostat')

    if (present(p_ext) .and. present(stress_ext)) &
         call system_abort('barostat_initialise: got both p_ext and stress_ext')

    this%type = type
    this%epsilon_v = 0.0_dp
    this%epsilon_f = 0.0_dp
    this%epsilon_f = 0.0_dp
    if (present(deformation_grad)) then
      this%epsilon_r = deformation_grad
    else
      this%epsilon_r = 0.0_dp; call add_identity(this%epsilon_r)
    endif
    this%lattice0_initialised = .false.

    select case(this%type)

    case(BAROSTAT_NONE) 

       this%stress_ext = 0.0_dp
       this%gamma_epsilon = 0.0_dp
       this%W_epsilon = 0.0_dp
       this%Ndof = 0.0_dp
       this%T = -1.0_dp

    case(BAROSTAT_HOOVER_LANGEVIN)

       if (.not.present(p_ext) .and. .not. present(stress_ext)) call system_abort('barostat initialise: p_ext or stress_ext is required for Hoover-Langevin barostat')
       if (present(p_ext) .and. present(stress_ext)) call system_abort('barostat initialise: p_ext and stress_ext conflict for Hoover-Langevin barostat')
       if (.not.present(gamma_epsilon)) call system_abort('barostat initialise: gamma_epsilon is required for Hoover-Langevin barostat')
       if (gamma_epsilon <= 0.0_dp) call system_abort('barostat initialise: gamma_epsilon must be > 0.0 for Hoover-Langevin')
       if (.not.present(Ndof)) call system_abort('barostat initialise: Ndof is required for Hoover-Langevin barostat')
       if (Ndof <= 0.0_dp) call system_abort('barostat initialise: Ndof must be > 0.0 for Hoover-Langevin')

       if (.not.present(W_epsilon) .and. .not. present(cell_volume)) call system_abort('barostat initialise: W_epsilon or cell_volume is required for Hoover-Langevin barostat')
       if (present(W_epsilon)) then
	  if (W_epsilon <= 0.0_dp) call system_abort('barostat initialise: W must be > 0.0 for Hoover-Langevin')
       endif
       if (present(cell_volume)) then
	  if (cell_volume <= 0.0_dp) call system_abort('barostat initialise: cell_volume must be > 0.0 for Hoover-Langevin')
       endif

       if (present(p_ext)) then ! S = - p_ext I
	  this%stress_ext = 0.0_dp
	  call add_identity(this%stress_ext)
	  this%stress_ext = -p_ext * this%stress_ext
       else
	  this%stress_ext = stress_ext
       endif

       if (present(hydrostatic_strain)) this%hydrostatic_strain = hydrostatic_strain
       if (present(diagonal_strain)) this%diagonal_strain = diagonal_strain
       if (present(finite_strain_formulation)) this%finite_strain_formulation = finite_strain_formulation

       this%Ndof = Ndof
       if (present(W_epsilon)) then
	  this%W_epsilon = W_epsilon
       else
	  this%W_epsilon = barostat_mass(maxval(abs(this%stress_ext)), cell_volume, Ndof, gamma_epsilon, T, W_epsilon_factor)
       endif
       use_T = optional_default(-1.0_dp, T)
       if (use_T > 0.0_dp) then
	 this%gamma_epsilon = gamma_epsilon
       else
	 this%gamma_epsilon = 0.0_dp
       endif
       this%T = use_T

    case default
       call system_abort("Invalid barostat type " // type)

    end select

  end subroutine barostat_initialise

  subroutine barostat_finalise(this)

    type(barostat), intent(inout) :: this

    this%type  = BAROSTAT_NONE  
    this%gamma_epsilon = 0.0_dp
    this%stress_ext     = 0.0_dp
    this%hydrostatic_strain     = .true.
    this%diagonal_strain     = .true.
    this%finite_strain_formulation     = .false.
    this%Ndof  = 0.0_dp
    this%W_epsilon = 0.0_dp
    this%epsilon_v = 0.0_dp
    this%epsilon_f = 0.0_dp
    this%epsilon_r = 0.0_dp; call add_identity(this%epsilon_r)
    this%T = -1.0_dp
    this%lattice0_initialised = .false.

  end subroutine barostat_finalise

  subroutine barostat_assignment(to,from)

    type(barostat), intent(out) :: to
    type(barostat), intent(in)  :: from

    to%type  = from%type      
    to%gamma_epsilon = from%gamma_epsilon 
    to%stress_ext     = from%stress_ext
    to%hydrostatic_strain     = from%hydrostatic_strain
    to%diagonal_strain     = from%diagonal_strain
    to%finite_strain_formulation     = from%finite_strain_formulation
    to%W_epsilon = from%W_epsilon
    to%epsilon_v = from%epsilon_v
    to%epsilon_f = from%epsilon_f
    to%epsilon_r = from%epsilon_r
    to%T = from%T

  end subroutine barostat_assignment

  subroutine barostat_update_barostat(this,p_ext,stress_ext,hydrostatic_strain,diagonal_strain,finite_strain_formulation,W_epsilon,T)

    type(barostat),   intent(inout) :: this
    real(dp), optional, intent(in)    :: p_ext
    real(dp), optional, intent(in)    :: stress_ext(3,3)
    logical, optional, intent(in)    :: hydrostatic_strain, diagonal_strain
    logical, optional, intent(in)    :: finite_strain_formulation
    real(dp), optional, intent(in)    :: W_epsilon
    real(dp), optional, intent(in)    :: T

    if( present(p_ext) .and. present(stress_ext)) call system_abort("barostat_update_barostat: got both p_ext and stress_ext, conflict")
    if( present(p_ext) ) then
      this%stress_ext = 0.0_dp
      call add_identity(this%stress_ext)
      this%stress_ext = -p_ext * this%stress_ext  
    endif
    if( present(stress_ext) ) this%stress_ext = stress_ext
    if( present(hydrostatic_strain) ) this%hydrostatic_strain = hydrostatic_strain
    if( present(diagonal_strain) ) this%diagonal_strain = diagonal_strain
    if( present(finite_strain_formulation) ) this%finite_strain_formulation = finite_strain_formulation
    if( present(W_epsilon) ) this%W_epsilon = W_epsilon
    if( present(T) ) this%T = T

  endsubroutine barostat_update_barostat

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X PRINTING
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine barostat_print(this,file)

    type(barostat),           intent(in) :: this
    type(inoutput), optional, intent(in) :: file

    select case(this%type)

       case(BAROSTAT_NONE)
	  call print('Barostat off',file=file)

       case(BAROSTAT_HOOVER_LANGEVIN)
	  call print('Hoover-Langevin, stress_ext = '// &
	  round(this%stress_ext(1,1)/GPA,5)//' '// &
	  round(this%stress_ext(1,2)/GPA,5)//' '// &
	  round(this%stress_ext(1,3)/GPA,5)//' '// &
	  round(this%stress_ext(2,1)/GPA,5)//' '// &
	  round(this%stress_ext(2,2)/GPA,5)//' '// &
	  round(this%stress_ext(2,3)/GPA,5)//' '// &
	  round(this%stress_ext(3,1)/GPA,5)//' '// &
	  round(this%stress_ext(3,2)/GPA,5)//' '// &
	  round(this%stress_ext(3,3)/GPA,5)//' GPa , hydrostatic_strain = '//this%hydrostatic_strain//&
	    ' diagonal_strain = '//this%diagonal_strain//' finite_strain_formulation = '//this%finite_strain_formulation//&
	    ' gamma_epsilon = '//round(this%gamma_epsilon,5)//' fs^-1, '// &
	    ' W_epsilon = '//round(this%W_epsilon,5)//' eV/fs, T = '//round(this%T,2)//' K, Ndof = '// round(this%Ndof,1),file=file)

    end select
    
  end subroutine barostat_print

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X SETTING Ndof
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine set_degrees_of_freedom_int(this,Ndof)

    type(barostat), intent(inout) :: this
    integer,          intent(in)    :: Ndof

    this%Ndof = real(Ndof,dp)

  end subroutine set_degrees_of_freedom_int

  subroutine set_degrees_of_freedom_real(this,Ndof)

    type(barostat), intent(inout) :: this
    real(dp),         intent(in)    :: Ndof

    this%Ndof = Ndof

  end subroutine set_degrees_of_freedom_real

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X CHOOSING NOSE-HOOVER(-LANGEVIN) MASS
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  pure function barostat_mass_int(p_ext, cell_volume, Ndof, gamma_epsilon, T, W_epsilon_factor) result(W_epsilon)

    real(dp), intent(in) :: p_ext, cell_volume
    integer,  intent(in) :: Ndof
    real(dp), intent(in) :: gamma_epsilon
    real(dp), intent(in), optional   :: T
    real(dp), intent(in), optional   :: W_epsilon_factor
    real(dp) :: W_epsilon ! result

    W_epsilon = barostat_mass_real(p_ext, cell_volume, real(Ndof,dp), gamma_epsilon, T, W_epsilon_factor)

  end function barostat_mass_int

  pure function barostat_mass_real(p_ext, cell_volume, Ndof, gamma_epsilon, T, W_epsilon_factor) result(W_epsilon)

    real(dp), intent(in) :: p_ext, cell_volume
    real(dp), intent(in) :: Ndof
    real(dp), intent(in) :: gamma_epsilon
    real(dp), intent(in), optional   :: T
    real(dp), intent(in), optional   :: W_epsilon_factor
    real(dp) :: W_epsilon ! result

    real(dp) :: use_T, use_W_epsilon_factor

    use_T = optional_default(BAROSTAT_MIN_TEMP, T)
    if (use_T <= 0.0_dp) use_T = BAROSTAT_MIN_TEMP

    use_W_epsilon_factor=optional_default(1.0_dp, W_epsilon_factor)

    W_epsilon = use_W_epsilon_factor*max( 9.0_dp*abs(p_ext)*cell_volume/((gamma_epsilon*2.0_dp*PI)**2), &
		     (Ndof+3.0_dp)*BOLTZMANN_K*use_T/((gamma_epsilon*2.0_dp*PI)**2) )

  end function barostat_mass_real

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X BAROSTAT ROUTINES
  !X
  !X These routines interleave the usual velocity Verlet steps and thermostat steps
  !X should modify the velocities accelerations, and positions as required:
  !X
  !X advance_verlet1
  !X   00 (barostat_pre_vel1) *
  !X   10 (thermostat_pre_vel1) *
  !X   20 v(t+dt/2) = v(t) + a(t)dt/2
  !X   30 (thermostat_post_vel1_pre_pos)
  !X   40 r(t+dt) = r(t) + v(t+dt/2)dt
  !X   50 (thermostat_post_pos_pre_calc)
  !X   60 (barostat_post_pos_pre_calc) *
  !X calc F, virial
  !X advance_verlet2
  !X   70 (thermostat_post_calc_pre_vel2) *
  !X   80 v(t+dt) = v(t+dt/2) + a(t+dt)dt/2
  !X   90 (thermostat_post_vel2) *
  !X   100 (barostat_post_vel2) *
  !X
  !X * marks routines that are needed in the refactored barostat/thermostat plan
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine barostat_pre_vel1(this,at,dt,virial)

    type(barostat), intent(inout) :: this
    type(atoms),      intent(inout) :: at
    real(dp),         intent(in)    :: dt
    real(dp), dimension(3,3), optional, intent(in) :: virial

    real(dp) :: vel_decay(3,3), pos_scale(3,3), lattice_p(3,3), F_inv(3,3), F_tmp(3,3)

    select case(this%type)

       !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !X
       !X HOOVER-LANGEVIN
       !X
       !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    case(BAROSTAT_HOOVER_LANGEVIN)
       !TIME_PROPAG_TEX 00 
       !TIME_PROPAG_TEX 00 {\color {blue}
       !TIME_PROPAG_TEX 00 before all Verlet and thermostat steps (barostat\_pre\_vel1)
       !TIME_PROPAG_TEX 00
       !TIME_PROPAG_TEX 00 $$ \epsilon_v = \epsilon_v \exp\left( -(\tau/2) \gamma_\epsilon \epsilon_r^{-1}\right)^T$$
       !TIME_PROPAG_TEX 00 $$ \epsilon_v = \epsilon_v +  (\tau/2) \epsilon_f/W_\epsilon $$
       !TIME_PROPAG_TEX 00 $$ v = \exp \left(  -(\tau/2) (1 + 3/N_d) \epsilon_v \right) v $$
       !TIME_PROPAG_TEX 00 $$ (r,h) = \exp\left( (\tau/2) \epsilon_v \right) (r,h) $$
       !TIME_PROPAG_TEX 00 }
       !TIME_PROPAG_TEX 00

       if (.not. present(virial)) call system_abort("barostat_pre_vel1 needs virial")

       if (.not. this%lattice0_initialised) then
	  this%lattice0_initialised = .true.
	  this%lattice0 = at%lattice
	  call matrix3x3_inverse(this%lattice0, this%lattice0_inv)
	  this%F_tmp_norm = 0.0_dp
       endif

       ! half step epsilon_v drag 
       if (this%finite_strain_formulation) then
	  call matrix3x3_inverse(this%epsilon_r, F_inv)
	  vel_decay = transpose(matrix_exp(-0.5_dp*dt*this%gamma_epsilon*F_inv))
	  vel_decay = 0.5_dp*(vel_decay + transpose(vel_decay))
	  this%epsilon_v = this%epsilon_v .mult. vel_decay
       else
	  this%epsilon_v = this%epsilon_v * exp(-0.5_dp*dt*this%gamma_epsilon)
       endif
       ! half step epsilon_v force (from last step) parts 
       this%epsilon_v = this%epsilon_v + 0.5_dp*dt*this%epsilon_f/this%W_epsilon

       ! half step barostat drag part of v
       vel_decay = matrix_exp(-0.5_dp*dt*((1.0_dp + 3.0_dp/this%Ndof)*this%epsilon_v))
       vel_decay = 0.5_dp*(vel_decay + transpose(vel_decay))
       at%velo = vel_decay .mult. at%velo

       ! half step position affine defomration
       pos_scale = matrix_exp(0.5_dp*dt*this%epsilon_v)
       pos_scale = 0.5_dp*(pos_scale + transpose(pos_scale))
       lattice_p = pos_scale .mult. at%lattice
       ! remove any rotation from new lattice
       F_tmp = lattice_p .mult. this%lattice0_inv
       F_tmp = 0.5_dp*(F_tmp + transpose(F_tmp))
       if (abs(deviation_from_identity(F_tmp) - this%F_tmp_norm) > 5.0e-2) then
	 call print_warning("barostat transformation projecting away rotations is very different from previous one.  Did the lattice change?")
       endif
       this%F_tmp_norm = deviation_from_identity(F_tmp)
       lattice_p = F_tmp .mult. this%lattice0
       call set_lattice(at, lattice_p, scale_positions=.true.)

    end select

  end subroutine barostat_pre_vel1
  
  subroutine barostat_post_vel1_pre_pos(this,at,dt,virial)

    type(barostat), intent(inout) :: this
    type(atoms),      intent(inout) :: at
    real(dp),         intent(in)    :: dt
    real(dp), dimension(3,3), optional, intent(in) :: virial

  end subroutine barostat_post_vel1_pre_pos

  subroutine barostat_post_pos_pre_calc(this,at,dt,virial)

    type(barostat), intent(inout) :: this
    type(atoms),      intent(inout) :: at
    real(dp),         intent(in)    :: dt
    real(dp), dimension(3,3), optional, intent(in) :: virial

    real(dp) :: pos_scale(3,3), lattice_p(3,3), F_tmp(3,3)

    select case(this%type)

      case(BAROSTAT_HOOVER_LANGEVIN)
        !TIME_PROPAG_TEX 60 
        !TIME_PROPAG_TEX 60 {\color {blue}
        !TIME_PROPAG_TEX 60 after Verlet pos step, before force calc (barostat\_post\_pos\_pre\_calc)
        !TIME_PROPAG_TEX 60
        !TIME_PROPAG_TEX 60 $$ (r,h) = \exp\left( (\tau/2) \epsilon_v \right) (r,h) $$
        !TIME_PROPAG_TEX 60 }
        !TIME_PROPAG_TEX 60 
        if (.not. present(virial)) call system_abort("barostat_post_pos_pre_calc needs virial")

        ! half step position affine defomration
        pos_scale = matrix_exp(0.5_dp*dt*this%epsilon_v)
        pos_scale = 0.5_dp*(pos_scale + transpose(pos_scale))
        lattice_p = pos_scale .mult. at%lattice
        ! remove any rotation from new lattice
        F_tmp = lattice_p .mult. this%lattice0_inv
        F_tmp = 0.5_dp*(F_tmp + transpose(F_tmp))
        if (abs(deviation_from_identity(F_tmp) - this%F_tmp_norm) > 5.0e-2) then
	   call print_warning("barostat transformation projecting away rotations is very different from previous one.  Did the lattice change?")
        endif
        this%F_tmp_norm = deviation_from_identity(F_tmp)
        lattice_p = F_tmp .mult. this%lattice0
        call set_lattice(at, lattice_p, scale_positions=.true.)

    end select

  end subroutine barostat_post_pos_pre_calc

  function deviation_from_identity(a)
     real(dp), intent(in) :: a(:, :)
     real(dp) :: deviation_from_identity

     integer :: i, j

     deviation_from_identity = 0.0_dp
     do i=1, size(a, 1)
     do j=1, size(a, 2)
       if (i == j) then
	  deviation_from_identity = deviation_from_identity + (a(i,j)-1.0_dp)**2
       else
	  deviation_from_identity = deviation_from_identity + (a(i,j))**2
       endif
     end do
     end do

     deviation_from_identity = sqrt(deviation_from_identity)

  end function deviation_from_identity

  subroutine barostat_post_calc_pre_vel2(this,at,dt,virial)
    
    type(barostat), intent(inout) :: this
    type(atoms),      intent(inout) :: at
    real(dp),         intent(in)    :: dt
    real(dp), dimension(3,3), optional, intent(in) :: virial

  end subroutine barostat_post_calc_pre_vel2
  
  subroutine barostat_post_vel2(this,at,dt,virial)

    type(barostat), intent(inout) :: this
    type(atoms),      intent(inout) :: at
    real(dp),         intent(in)    :: dt
    real(dp), dimension(3,3), optional, intent(in) :: virial

    integer :: i, j
    real(dp) :: vel_decay(3,3), volume_p, rand_f_cell(3,3), kinetic_virial(3,3), F_det, F_T_inv(3,3), F_inv(3,3), epsilon_f_hydrostatic
    real(dp) :: rand_f_scale

    select case(this%type)

    case(BAROSTAT_HOOVER_LANGEVIN)
       !TIME_PROPAG_TEX 100
       !TIME_PROPAG_TEX 100 {\color {blue}
       !TIME_PROPAG_TEX 100 after all Verlet and thermostat steps (barostat\_post\_vel2)
       !TIME_PROPAG_TEX 100
       !TIME_PROPAG_TEX 100 $$ v = \exp \left(  -(\tau/2) (1 + 3/N_d) \epsilon_v \right) v $$
       !TIME_PROPAG_TEX 100 $$ \epsilon_r = \epsilon_r + \tau \epsilon_v $$
       !TIME_PROPAG_TEX 100 $r.v.$ is symmetric matrix, normal distribution entries, variance 3 on diagonal and 6 off diagonal
       !TIME_PROPAG_TEX 100 \begin{eqnarray*}
       !TIME_PROPAG_TEX 100 \epsilon_f  & = & \left[ 3 \left( \mathrm{vir} + \sum_i m_i v_i \otimes v_i + 
       !TIME_PROPAG_TEX 100                                 3 V \epsilon_r \sigma \epsilon_r^T / \left|\epsilon_r\right|  + 
       !TIME_PROPAG_TEX 100                  3/N_d \sum_i m_i  v_i \otimes v_i  \right) + \right.
       !TIME_PROPAG_TEX 100             \\  &  & \sqrt{2 k_B T \gamma_\epsilon W_\epsilon/\tau} r.v. \bigg] F^{-T}
       !TIME_PROPAG_TEX 100 \end{eqnarray*}
       !TIME_PROPAG_TEX 100 $$ \epsilon_v = \epsilon_v +  (\tau/2) \epsilon_f/W_\epsilon $$
       !TIME_PROPAG_TEX 100 $$ \epsilon_v = \epsilon_v \exp\left( -(\tau/2) \gamma_\epsilon \epsilon_r^{-1}\right)^T$$
       !TIME_PROPAG_TEX 100  }
       !TIME_PROPAG_TEX 100
       if (.not. present(virial)) call system_abort("barostat_pre_vel1 needs virial")

       call system_resync_rng()

       !Decay the velocities for dt/2 again barostat part
       vel_decay = matrix_exp(-0.5_dp*dt*((1.0_dp + 3.0_dp/this%Ndof)*this%epsilon_v))
       vel_decay = 0.5_dp*(vel_decay + transpose(vel_decay))
       at%velo(:,:) = vel_decay .mult. at%velo(:,:)

       this%epsilon_r = this%epsilon_r + dt*this%epsilon_v

       ! calc epsilon force
       volume_p = cell_volume(at)
       rand_f_cell = 0.0_dp
       do i=1, 3
       do j=i, 3
	  if (i == j) then
	     rand_f_scale = 3.0_dp
	  else
	     ! empirical - maybe have something to do with fact that matrix is symmetric so there are
	     ! fewer degrees of freedom than first appears?
	     rand_f_scale = 6.0_dp
	  endif
	  rand_f_cell(i,j) = sqrt(rand_f_scale)*sqrt(2.0_dp*BOLTZMANN_K*this%T*this%gamma_epsilon*this%W_epsilon/dt)*ran_normal()
       end do
       end do

       ! symmetrize rand_f_cell
       rand_f_cell(2,1) = rand_f_cell(1,2)
       rand_f_cell(3,1) = rand_f_cell(1,3)
       rand_f_cell(3,2) = rand_f_cell(2,3)

       kinetic_virial = matmul(at%velo*spread(at%mass,dim=1,ncopies=3),transpose(at%velo))
       kinetic_virial = 0.5_dp*(kinetic_virial + transpose(kinetic_virial))
       if (this%finite_strain_formulation) then
	  F_det = matrix3x3_det(this%epsilon_r)
	  call matrix3x3_inverse(transpose(this%epsilon_r), F_T_inv)

	  this%epsilon_f = (3.0_dp*(virial+kinetic_virial+volume_p*(this%epsilon_r .mult. this%stress_ext .mult. transpose(this%epsilon_r))/F_det + &
				    3.0_dp/this%Ndof*kinetic_virial) + rand_f_cell) .mult. F_T_inv
       else
	  this%epsilon_f = 3.0_dp*(virial+kinetic_virial+volume_p*this%stress_ext + &
				   3.0_dp/this%Ndof*kinetic_virial) + rand_f_cell
       endif

       if (this%diagonal_strain .or. this%hydrostatic_strain) then
	 this%epsilon_f(1,2) = 0.0_dp
	 this%epsilon_f(1,3) = 0.0_dp
	 this%epsilon_f(2,3) = 0.0_dp
	 this%epsilon_f(2,1) = 0.0_dp
	 this%epsilon_f(3,1) = 0.0_dp
	 this%epsilon_f(3,2) = 0.0_dp
	 if (this%hydrostatic_strain) then
	   epsilon_f_hydrostatic = (this%epsilon_f(1,1)+ this%epsilon_f(2,2)+ this%epsilon_f(3,3)) / 3.0_dp
	   this%epsilon_f(1,1) = epsilon_f_hydrostatic
	   this%epsilon_f(2,2) = epsilon_f_hydrostatic
	   this%epsilon_f(3,3) = epsilon_f_hydrostatic
	 endif
       endif

       ! half step with epsilon force
       this%epsilon_v = this%epsilon_v + 0.5_dp*dt*this%epsilon_f/this%W_epsilon
       ! half step with epsilon drag
       if (this%finite_strain_formulation) then
	  call matrix3x3_inverse(this%epsilon_r, F_inv)
	  vel_decay = transpose(matrix_exp(-0.5_dp*dt*this%gamma_epsilon*F_inv))
	  vel_decay = 0.5_dp*(vel_decay + transpose(vel_decay))
	  this%epsilon_v = this%epsilon_v .mult. vel_decay 
       else
	  this%epsilon_v = this%epsilon_v * exp(-0.5_dp*dt*this%gamma_epsilon)
       endif

       call print("BAROSTAT_KE "//barostat_KE(this))

    end select

  end subroutine barostat_post_vel2

  function barostat_KE(this)
    type(barostat), intent(in) :: this
    real(dp) :: barostat_KE

    select case(this%type)
      case(BAROSTAT_HOOVER_LANGEVIN)
	 barostat_KE = 0.5_dp*this%W_epsilon*sum(this%epsilon_v*transpose(this%epsilon_v))
    end select
  end function barostat_KE

end module barostat_module
