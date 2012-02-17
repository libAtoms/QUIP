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
  use atoms_module

  implicit none

  real(dp), parameter :: BAROSTAT_MIN_TEMP = 1.0_dp

  integer, parameter :: BAROSTAT_NONE                 = 0, &
                        BAROSTAT_HOOVER_LANGEVIN      = 1

  !% Langevin barostat - fluctuation/dissipation ---
  !% EOM from Quigley, D. and Probert, M.I.J., \emph{J. Chem. Phys.}, 
  !% {\bfseries 120} 11432
  !% rediscretized by Noam Bernstein

  type barostat

     integer  :: type  = BAROSTAT_NONE     !% One of the types listed above
     real(dp) :: gamma_epsilon = 0.0_dp   !% Friction coefficient for Langevin part
     real(dp) :: T = -1.0_dp   !% Temperature for Langevin part
     real(dp) :: Ndof  = 0.0_dp   !% The number of degrees of freedom of atoms attached to this barostat
     real(dp) :: p_ext = 0.0_dp       !% External pressure
     real(dp) :: W_epsilon = 0.0_dp     !% Fictious cell mass in Langevin NPT
     real(dp) :: epsilon_v = 0.0_dp !% Velocity of barostat variable
     real(dp) :: epsilon_f = 0.0_dp !% Force on barostat variable

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

  subroutine barostat_initialise(this,type,p_ext,cell_volume,W_epsilon,Ndof,gamma_epsilon,T)

    type(barostat),   intent(inout) :: this
    integer,            intent(in)    :: type
    real(dp), optional, intent(in)    :: p_ext
    real(dp), optional, intent(in)    :: cell_volume
    real(dp), optional, intent(in)    :: W_epsilon
    real(dp), optional, intent(in)    :: Ndof
    real(dp), optional, intent(in)    :: gamma_epsilon
    real(dp), optional, intent(in)    :: T

    real(dp) :: use_T

    if (type /= BAROSTAT_NONE .and. .not.present(p_ext)) &
         call system_abort('initialise: p_ext must be specified when turning on a barostat')

    this%type = type
    this%epsilon_v = 0.0_dp
    this%epsilon_f = 0.0_dp

    select case(this%type)

    case(BAROSTAT_NONE) 

       this%p_ext = 0.0_dp
       this%gamma_epsilon = 0.0_dp
       this%W_epsilon = 0.0_dp
       this%Ndof = 0.0_dp
       this%T = -1.0_dp

    case(BAROSTAT_HOOVER_LANGEVIN)

       if (.not.present(p_ext)) call system_abort('barostat initialise: p_ext is required for Hoover-Langevin barostat')
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

       this%p_ext = p_ext
       this%Ndof = Ndof
       if (present(W_epsilon)) then
	  this%W_epsilon = W_epsilon
       else
	  this%W_epsilon = barostat_mass(p_ext, cell_volume, Ndof, gamma_epsilon, T)
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
    this%p_ext     = 0.0_dp
    this%Ndof  = 0.0_dp
    this%W_epsilon = 0.0_dp
    this%epsilon_v = 0.0_dp
    this%epsilon_f = 0.0_dp
    this%T = -1.0_dp

  end subroutine barostat_finalise

  subroutine barostat_assignment(to,from)

    type(barostat), intent(out) :: to
    type(barostat), intent(in)  :: from

    to%type  = from%type      
    to%gamma_epsilon = from%gamma_epsilon 
    to%p_ext     = from%p_ext
    to%W_epsilon = from%W_epsilon
    to%epsilon_v = from%epsilon_v
    to%epsilon_f = from%epsilon_f
    to%T = from%T

  end subroutine barostat_assignment

  subroutine barostat_update_barostat(this,p_ext,W_epsilon,T)

    type(barostat),   intent(inout) :: this
    real(dp), optional, intent(in)    :: p_ext
    real(dp), optional, intent(in)    :: W_epsilon
    real(dp), optional, intent(in)    :: T

    if( present(p_ext) ) this%p_ext = p_ext
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
	  call print('Hoover-Langevin, p_ext = '//round(this%p_ext/GPA,5)//' GPa , gamma_epsilon = '//round(this%gamma_epsilon,5)//' fs^-1, '// &
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

  pure function barostat_mass_int(p_ext, cell_volume, Ndof, gamma_epsilon, T) result(W_epsilon)

    real(dp), intent(in) :: p_ext, cell_volume
    integer,  intent(in) :: Ndof
    real(dp), intent(in) :: gamma_epsilon
    real(dp), intent(in), optional   :: T
    real(dp) :: W_epsilon ! result

    W_epsilon = barostat_mass_real(p_ext, cell_volume, real(Ndof,dp), gamma_epsilon, T)

  end function barostat_mass_int

  pure function barostat_mass_real(p_ext, cell_volume, Ndof, gamma_epsilon, T) result(W_epsilon)

    real(dp), intent(in) :: p_ext, cell_volume
    real(dp), intent(in) :: Ndof
    real(dp), intent(in) :: gamma_epsilon
    real(dp), intent(in), optional   :: T
    real(dp) :: W_epsilon ! result

    real(dp) :: use_T

    use_T = optional_default(BAROSTAT_MIN_TEMP, T)
    if (use_T <= 0.0_dp) use_T = BAROSTAT_MIN_TEMP

    W_epsilon = max( 9.0_dp*abs(p_ext)*cell_volume/((gamma_epsilon*2.0_dp*PI)**2), &
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

    real(dp) :: decay, lattice_p(3,3)

    select case(this%type)

       !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !X
       !X HOOVER-LANGEVIN
       !X
       !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    case(BAROSTAT_HOOVER_LANGEVIN)
       !TIME_PROPAG_TEX 00
       !TIME_PROPAG_TEX 00 before all Verlet and thermostat steps (barostat\_pre\_vel1)
       !TIME_PROPAG_TEX 00
       !TIME_PROPAG_TEX 00 $$ \epsilon_v = \epsilon_v \exp\left( -(\tau/2) \gamma_\epsilon \right)$$
       !TIME_PROPAG_TEX 00 $$ \epsilon_v = \epsilon_v +  (\tau/2) \epsilon_f/W_\epsilon $$
       !TIME_PROPAG_TEX 00 $$ v = v \exp \left(  -(\tau/2) (1 + 3/N_d) \epsilon_v \right) $$
       !TIME_PROPAG_TEX 00 $$ (r,h) = (r,h) \exp\left( (\tau/2) \epsilon_v \right) $$
       !TIME_PROPAG_TEX 00 

       if (.not. present(virial)) call system_abort("barostat_pre_vel1 needs virial")

       ! half step epsilon_v drag 
       this%epsilon_v = this%epsilon_v*exp(-0.5_dp*dt*this%gamma_epsilon)
       ! half step epsilon_v force (from last step) parts 
       this%epsilon_v = this%epsilon_v + 0.5_dp*dt*this%epsilon_f/this%W_epsilon

       ! half step barostat drag part of v
       decay = exp(-0.5_dp*dt*((1.0_dp + 3.0_dp/this%Ndof)*this%epsilon_v))
       at%velo = at%velo*decay

       ! half step position affine defomration
       lattice_p = at%lattice
       call set_lattice(at, lattice_p*exp(0.5_dp*dt*this%epsilon_v), scale_positions=.false.)
       at%pos = at%pos*exp(0.5_dp*dt*this%epsilon_v)


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

    real(dp) :: lattice_p(3,3)

    select case(this%type)

      case(BAROSTAT_HOOVER_LANGEVIN)
       !TIME_PROPAG_TEX 60
       !TIME_PROPAG_TEX 60 after Verlet pos step, before force calc (barostat\_post\_pos\_pre\_calc)
       !TIME_PROPAG_TEX 60
       !TIME_PROPAG_TEX 60 $$ (r,h) = (r,h) \exp\left( (\tau/2) \epsilon_v \right) $$
       !TIME_PROPAG_TEX 60 
       if (.not. present(virial)) call system_abort("barostat_post_pos_pre_calc needs virial")

        ! half step position affine defomration
        lattice_p = at%lattice
        call set_lattice(at, lattice_p*exp(0.5_dp*dt*this%epsilon_v), scale_positions=.false.)
        at%pos = at%pos*exp(0.5_dp*dt*this%epsilon_v)

    end select

  end subroutine barostat_post_pos_pre_calc

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

    real(dp) :: decay, volume_p, rand_f_cell

    select case(this%type)

    case(BAROSTAT_HOOVER_LANGEVIN)
       !TIME_PROPAG_TEX 100
       !TIME_PROPAG_TEX 100 after all Verlet and thermostat steps (barostat\_post\_vel2)
       !TIME_PROPAG_TEX 100
       !TIME_PROPAG_TEX 100 $$ v = v \exp \left(  -(\tau/2) (1 + 3/N_d) \epsilon_v \right) $$
       !TIME_PROPAG_TEX 100 $$ \epsilon_f = \left( Tr[\mathrm{vir}] + \sum_i m_i \left|v_i\right|^2 - 3 V P\right) + 3/N_d \sum_i m_i \left| v_i \right|^2 + \sqrt{2 k_B T \gamma_\epsilon W_\epsilon/\tau} r $$
       !TIME_PROPAG_TEX 100 $$ \epsilon_v = \epsilon_v +  (\tau/2) \epsilon_f/W_\epsilon $$
       !TIME_PROPAG_TEX 100 $$ \epsilon_v = \epsilon_v \exp\left( -(\tau/2) \gamma_\epsilon \right)$$
       !TIME_PROPAG_TEX 100 
       if (.not. present(virial)) call system_abort("barostat_pre_vel1 needs virial")

       call system_resync_rng()

       !Decay the velocities for dt/2 again barostat part
       decay = exp(-0.5_dp*dt*((1.0_dp + 3.0_dp/this%Ndof)*this%epsilon_v))
       at%velo(:,:) = at%velo(:,:)*decay

       ! calc epsilon force
       volume_p = cell_volume(at)
       rand_f_cell = sqrt(2.0_dp*BOLTZMANN_K*this%T*this%gamma_epsilon*this%W_epsilon/dt)*ran_normal()
       this%epsilon_f = (trace(virial)+sum(at%mass*sum(at%velo**2,dim=1))-3.0_dp*volume_p*this%p_ext) + &
	  3.0_dp/this%Ndof*sum(at%mass*sum(at%velo**2,dim=1)) + rand_f_cell

       ! half step with epsilon force
       this%epsilon_v = this%epsilon_v + 0.5_dp*dt*this%epsilon_f/this%W_epsilon
       ! half step with epsilon drag
       this%epsilon_v = this%epsilon_v*exp(-0.5_dp*dt*this%gamma_epsilon)

       call print("BAROSTAT_KE "//(0.5_dp*this%W_epsilon*this%epsilon_v**2))

    end select

  end subroutine barostat_post_vel2

  function barostat_KE(this)
    type(barostat), intent(in) :: this
    real(dp) :: barostat_KE

    select case(this%type)
      case(BAROSTAT_HOOVER_LANGEVIN)
	 barostat_KE = 0.5_dp*this%W_epsilon*this%epsilon_v**2
    end select
  end function barostat_KE

end module barostat_module
