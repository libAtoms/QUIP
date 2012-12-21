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

module elasticity_module

  use libatoms_module

  use Potential_module
  implicit none

  private

  public :: calc_elastic_constants
  interface calc_elastic_constants
     module procedure pot_calc_elastic_constants
  end interface

  public :: youngs_modulus, poisson_ratio, graphene_elastic, einstein_frequencies, elastic_fields


contains

  subroutine pot_calc_elastic_constants(this, at, fd, args_str, c, c0, relax_initial, return_relaxed, relax_tol, relax_method)
    type(Potential), intent(inout) :: this
    type(Atoms), intent(inout) :: at !% Atoms object for which to compute $C_{ij}$
    real(dp), intent(in), optional :: fd !% Finite strain to apply. Default $10^{-2}$.
    character(len=*), intent(in), optional :: args_str !% Optional args_str to pass to 'minim'
    real(dp), intent(out), optional :: c(6,6) !% Elastic constants (with relaxation)
    real(dp), intent(out), optional :: c0(6,6) !% Elastic constants (without relaxation)
    logical, optional :: relax_initial !% Should the initial cell be relaxed?
    logical, optional :: return_relaxed !% If true, overwrite 'at' with relaxed positions and lattice (default false)
    real(dp), optional :: relax_tol !% Relaxation df$^2$ tolerance. Default 1e-8
    character(len=*), optional :: relax_method !% method to pass for minim

    logical my_relax_initial, my_return_relaxed
    integer ii, jj
    type(atoms) :: at_bulk
    type(atoms) :: at_t
    real(dp) :: my_fd, my_relax_tol
    real(dp) :: volume
    real(dp) :: Fp(3,3), Fm(3,3)
    real(dp) :: V0p(3,3), V0m(3,3), Vp(3,3), Vm(3,3)
    integer iter
    character(len=100) :: my_relax_method 

    real(dp) :: e0, V0(3,3)

    my_fd = optional_default(1.0e-2_dp, fd)
    my_relax_tol = optional_default(1e-8_dp, relax_tol)
    my_relax_initial = optional_default(.true., relax_initial)
    my_return_relaxed = optional_default(.false., return_relaxed)
    my_relax_method = trim(optional_default('cg_n', relax_method))

    at_bulk = at
    call set_cutoff(at_bulk, cutoff(this))
    call calc_connect(at_bulk)

    if (current_verbosity() > PRINT_VERBOSE) then
       call verbosity_push_decrement(PRINT_NERD)
       call calc(this, at_bulk, energy=e0, virial=v0, args_str=args_str)
       call verbosity_pop()
       call print("initial config")
       call write(at_bulk, 'stdout')
       call print("e " // e0)
       call print("V")
       call print(V0)
    endif

    if (my_relax_initial) then
       call verbosity_push_decrement(PRINT_VERBOSE)
       iter = minim(this, at_bulk, my_relax_method, my_relax_tol, 1000, 'FAST_LINMIN', do_print=.false., &
            do_pos=.true., do_lat=.true., args_str=args_str)
       call verbosity_pop()

       if (current_verbosity() >= PRINT_VERBOSE) then
          call verbosity_push_decrement(PRINT_VERBOSE)
          call calc(this, at_bulk, energy=e0, virial=v0, args_str=args_str)
          call verbosity_pop()
          call print("relaxed config")
          call write(at_bulk, 'stdout')
          call print("e " // e0)
          call print("V")
          call print(V0)
       endif

       if (my_return_relaxed) at = at_bulk
    endif

    volume = cell_volume(at_bulk)

    call print("volume " // volume, PRINT_VERBOSE)

    do ii=1, 3
       do jj=ii, 3
          call print("doing strain [" // ii // "," // jj//"]", PRINT_VERBOSE)
          Fp = 0.0_dp; call add_identity(Fp)
          Fm = 0.0_dp; call add_identity(Fm)

          Fp(ii,jj) = Fp(ii,jj) + my_fd
          Fm(ii,jj) = Fm(ii,jj) - my_fd

          at_t = at_bulk
          call set_lattice(at_t, Fp .mult. at_t%lattice, scale_positions=.true.)
          call calc_connect(at_t)
          call calc(this, at_t, energy=e0, virial=V0p, args_str=args_str)
          call verbosity_push_decrement(PRINT_VERBOSE)
          call print("plus perturbed config")
          call write(at_t, 'stdout')
          call print("E0p" // e0)
          call print("V0p")
          call print(V0p)
          call verbosity_pop()

          if (present(c)) then
             call verbosity_push_decrement(PRINT_VERBOSE)
             iter = minim(this, at_t, my_relax_method, my_relax_tol, 1000, 'FAST_LINMIN', do_print=.false., &
                  do_pos=.true., do_lat=.false., args_str=args_str)
             call verbosity_pop()

             call calc_connect(at_t)
             call calc(this, at_t, energy=e0, virial=Vp, args_str=args_str)
             call verbosity_push_decrement(PRINT_VERBOSE)
             call print("plus perturbed relaxed config")
             call write(at_t, 'stdout')
             call print("Ep" // e0)
             call print("Vp")
             call print(Vp)
             call verbosity_pop()
          endif

          at_t = at_bulk
          call set_lattice(at_t, Fm .mult. at_t%lattice, scale_positions=.true.)
          call calc_connect(at_t)
          call calc(this, at_t, energy=e0, virial=V0m, args_str=args_str)
          call verbosity_push_decrement(PRINT_VERBOSE)
          call print("minus perturbed config")
          call write(at_t, 'stdout')
          call print("E0m" // e0)
          call print("V0m")
          call print(V0m)
          call verbosity_pop()

          if (present(c)) then
             call verbosity_push_decrement(PRINT_VERBOSE)
             iter = minim(this, at_t, my_relax_method, my_relax_tol, 1000, 'FAST_LINMIN', do_print=.false., &
                  do_pos=.true., do_lat=.false., args_str=args_str)
             call verbosity_pop()

             call calc_connect(at_t)
             call calc(this, at_t, energy=e0, virial=Vm, args_str=args_str)
             call verbosity_push_decrement(PRINT_VERBOSE)
             call print("minus perturbed relaxed config")
             call write(at_t, 'stdout')
             call print("Em" // e0)
             call print("Vm")
             call print(Vm)
             call verbosity_pop()
          endif

          if (present(c0)) then
             c0(strain_index(ii,jj),strain_index(1,1)) = V0m(1,1)-V0p(1,1)
             c0(strain_index(ii,jj),strain_index(1,2)) = V0m(1,2)-V0p(1,2)
             c0(strain_index(ii,jj),strain_index(1,3)) = V0m(1,3)-V0p(1,3)
             c0(strain_index(ii,jj),strain_index(2,2)) = V0m(2,2)-V0p(2,2)
             c0(strain_index(ii,jj),strain_index(2,3)) = V0m(2,3)-V0p(2,3)
             c0(strain_index(ii,jj),strain_index(3,3)) = V0m(3,3)-V0p(3,3)
             c0(strain_index(1,1),strain_index(ii,jj)) = c0(strain_index(ii,jj),strain_index(1,1))
             c0(strain_index(1,2),strain_index(ii,jj)) = c0(strain_index(ii,jj),strain_index(1,2))
             c0(strain_index(1,3),strain_index(ii,jj)) = c0(strain_index(ii,jj),strain_index(1,3))
             c0(strain_index(2,2),strain_index(ii,jj)) = c0(strain_index(ii,jj),strain_index(2,2))
             c0(strain_index(2,3),strain_index(ii,jj)) = c0(strain_index(ii,jj),strain_index(2,3))
             c0(strain_index(3,3),strain_index(ii,jj)) = c0(strain_index(ii,jj),strain_index(3,3))
          endif
          if (present(c)) then
             c(strain_index(ii,jj),strain_index(1,1)) = Vm(1,1)-Vp(1,1)
             c(strain_index(ii,jj),strain_index(1,2)) = Vm(1,2)-Vp(1,2)
             c(strain_index(ii,jj),strain_index(1,3)) = Vm(1,3)-Vp(1,3)
             c(strain_index(ii,jj),strain_index(2,2)) = Vm(2,2)-Vp(2,2)
             c(strain_index(ii,jj),strain_index(2,3)) = Vm(2,3)-Vp(2,3)
             c(strain_index(ii,jj),strain_index(3,3)) = Vm(3,3)-Vp(3,3)
             c(strain_index(1,1),strain_index(ii,jj)) = c(strain_index(ii,jj),strain_index(1,1))
             c(strain_index(1,2),strain_index(ii,jj)) = c(strain_index(ii,jj),strain_index(1,2))
             c(strain_index(1,3),strain_index(ii,jj)) = c(strain_index(ii,jj),strain_index(1,3))
             c(strain_index(2,2),strain_index(ii,jj)) = c(strain_index(ii,jj),strain_index(2,2))
             c(strain_index(2,3),strain_index(ii,jj)) = c(strain_index(ii,jj),strain_index(2,3))
             c(strain_index(3,3),strain_index(ii,jj)) = c(strain_index(ii,jj),strain_index(3,3))
          endif
       end do
    end do

    if (present(c0)) c0 = c0 / (2.0_dp*my_fd*volume)
    if (present(c)) c = c / (2.0_dp*my_fd*volume)

  end subroutine pot_calc_elastic_constants

  function strain_index(i, j)
    integer, intent(in) :: i, j
    integer :: strain_index

    if (i == j) strain_index = i
    if ((i == 2 .and. j == 3) .or. (i == 3 .and. j == 2)) strain_index = 4
    if ((i == 3 .and. j == 1) .or. (i == 1 .and. j == 3)) strain_index = 5
    if ((i == 1 .and. j == 2) .or. (i == 2 .and. j == 1)) strain_index = 6
  end function strain_index


  !% Calculate Youngs modulus $E_l$ from $6\times6$ elastic constants matrix $C_{ij}$
  !% This is the modulus for loading in the $l$ direction.
  !% Formula is from W. Brantley, Calculated elastic constants for stress problems associated 
  !% with semiconductor devices. J. Appl. Phys., 44, 534 (1973).
  function youngs_modulus(c, l) result(E)
    real(dp), intent(in) :: c(6,6)
    real(dp), dimension(3), intent(in) :: l
    real(dp) :: E

    real(dp) :: s(6,6)
    real(dp), dimension(3) :: lhat

    call inverse(c, s)

    ! Normalise directions
    lhat = l/norm(l)

    ! Youngs modulus in direction l, ratio of stress sigma_l 
    ! to strain response epsilon_l
    E = 1.0_dp/(S(1,1)-2.0_dp*(s(1,1)-s(1,2)-0.5_dp*s(4,4))*(lhat(1)*lhat(1)*lhat(2)*lhat(2) + &
         lhat(2)*lhat(2)*lhat(3)*lhat(3) + &
         lhat(1)*lhat(1)*lhat(3)*lhat(3)))

  end function youngs_modulus

  !% Calculate Poisson ratio $v_{lm}$ from $6\times6$ elastic
  !% constant matrix $C_{ij}$. This is the response in $m$ direction 
  !% to pulling in $l$ direction. Result is dimensionless.
  !% Formula is from W. Brantley, Calculated elastic constants for stress problems associated 
  !% with semiconductor devices. J. Appl. Phys., 44, 534 (1973).
  function poisson_ratio(cc, l, m) result(v)
    real(dp), intent(in) :: cc(6,6)
    real(dp), dimension(3), intent(in) :: l, m
    real(dp) :: v

    real(dp) :: s(6,6)
    real(dp), dimension(3) :: lhat, mhat

    call inverse(cc, s)

    ! Normalise directions
    lhat = l/norm(l)
    mhat = m/norm(m)

    ! Poisson ratio v_lm: response in m direction to strain in 
    ! l direction, v_lm = - epsilon_m/epsilon_l
    v = -((s(1,2) + (s(1,1)-s(1,2)-0.5_dp*s(4,4))*(lhat(1)*lhat(1)*mhat(1)*mhat(1) + &
         lhat(2)*lhat(2)*mhat(2)*mhat(2) + &
         lhat(3)*lhat(3)*mhat(3)*mhat(3))) / &
         (s(1,1) - 2.0_dp*(s(1,1)-s(1,2)-0.5_dp*s(4,4))*(lhat(1)*lhat(1)*lhat(2)*lhat(2) + &
         lhat(2)*lhat(2)*lhat(3)*lhat(3) + &
         lhat(1)*lhat(1)*lhat(3)*lhat(3))))
  end function poisson_ratio

  !% Calculate in-plane elastic constants of a graphene sheet with lattice
  !% parameter 'a' using the Potential 'pot'. On exit, 'poisson'
  !% will contain the in plane poisson ratio (dimensionless) and 'young' the
  !% in plane Young's modulus (GPa).
  subroutine Graphene_Elastic(pot, a, poisson, young, args_str, cb)
    type(Potential), intent(inout) :: pot
    real(dp), intent(out) :: a, poisson, young
    character(len=*), intent(in), optional :: args_str !% arg_str for potential_calc
    real(dp), intent(out), optional :: cb

    type(Atoms) :: cube, at, at2, tube
    real(dp) :: a0, v(3,3), eps(3,3)
    integer :: steps, tube_i, i, n, m
    real(dp) :: fix(3,3)
    integer, parameter :: nx = 2, ny = 3
    real(dp) :: tube_r(12), tube_energy(12), graphene_e_per_atom, radius, energy, c(1), chisq
    
    a0 = 1.45_dp
    cube = Graphene_Cubic(a0)

    call Supercell(at, cube, nx, ny, 1)
    call Set_Cutoff(at, cutoff(pot))
    call randomise(at%pos, 0.01_dp)
    call calc_connect(at)

    ! fix lattice in x and y directions
    fix = 1.0_dp
    fix(1,1) = 0.0_dp
    fix(2,2) = 0.0_dp
    call set_param_value(at, "Minim_Lattice_Fix", fix)

    ! Geometry optimise with variable lattice
    steps = minim(pot, at, 'cg', 1e-6_dp, 100, &
         'FAST_LINMIN', do_pos=.true.,do_lat=.true.,do_print=.false.)

    ! Set a to average of x and y lattice constants
    a = 0.5_dp*(at%lattice(1,1)/(3.0_dp*nx) + at%lattice(2,2)/(sqrt(3.0_dp)*ny))

    call calc(pot, at, energy=graphene_e_per_atom, args_str=args_str)
    graphene_e_per_atom = graphene_e_per_atom/at%N

    cube = Graphene_Cubic(a)
    call Supercell(at2, cube, nx, ny, 1)
    call Set_Cutoff(at2, cutoff(pot))
    call calc_connect(at2)

    ! Apply small strain in x direction
    eps = 0.0_dp; call add_identity(eps)
    eps(1,1) = eps(1,1)+0.001_dp
    call set_lattice(at2, eps .mult. at2%lattice, scale_positions=.true.)

    ! Fix lattice in x direction
    fix = 1.0_dp
    fix(1,1) = 0.0_dp
    call set_param_value(at, "Minim_Lattice_Fix", fix)

    ! Geometry optimse, allowing to contract in y direction
    steps = minim(pot, at2, 'cg', 1e-6_dp, 100, &
         'FAST_LINMIN', do_print=.false., do_pos=.true.,do_lat=.true.)

    poisson = -((at2%lattice(2,2) - at%lattice(2,2))/at%lattice(2,2))/ &
         ((at2%lattice(1,1) - at%lattice(1,1))/at%lattice(1,1))

    ! Calculate stress to find Young's modulus
    call Calc(pot, at2, virial=v)

    young = -v(1,1)/(eps(1,1)-1.0_dp)*GPA/cell_volume(at2)

    if (present(cb)) then
       ! Finally, calculate bending modulus by fitting strain energy per atom to
       ! curve 1/2 c_b/r**2 for a series of nanotubes. We use 12 nanotubes (n,m) 
       ! where n runs from 10 to 20 in steps of 2 and m is either 0 or n.
       tube_i = 1
       tube_r = 0.0_dp
       tube_energy = 0.0_dp
       call verbosity_push_decrement(PRINT_VERBOSE)
       do n=10,20,2
          do m=0,n,n
             radius = graphene_tube(tube, a, n, m, 3)
             call set_cutoff(tube, cutoff(pot))
             call calc_connect(tube)

             
             i = minim(pot, tube, 'cg', 1e-5_dp, 1000, 'FAST_LINMIN', do_print=.false., &
                  do_pos=.true., do_lat=.false.)
             
             call calc(pot, tube, energy=energy, args_str=args_str)

             tube_r(tube_i) = tube_radius(tube)
             tube_energy(tube_i) = energy/tube%N - graphene_e_per_atom
             tube_i = tube_i + 1
          end do
       end do
       call verbosity_pop()
       
       call least_squares(tube_r, tube_energy, (/ ( 1.0_dp, i=1,12) /), &
            c, chisq, inverse_square)
       cb = 2.0_dp*c(1)

       call finalise(tube)
    end if

    call Finalise(cube)
    call Finalise(at)
    call Finalise(at2)

  end subroutine Graphene_Elastic

  subroutine inverse_square(x,afunc)
    real(dp) :: x, afunc(:)

    afunc(1) = 1.0_dp/(x*x)
        
  end subroutine inverse_square


  function einstein_frequencies(pot, at, args_str, ii, delta) result(w_e)
    type(Potential), intent(inout) :: pot  !% Potential to use
    type(Atoms), intent(in) :: at            !% Atoms structure - should be equilibrium bulk configuation
    character(len=*), intent(in), optional :: args_str !% arg_str for potential_calc
    integer, optional, intent(in) :: ii       !% The atom to displace (default 1)
    real(dp), optional, intent(in) :: delta  !% How much to displace it (default 1e-4_dp)
    real(dp), dimension(3) :: w_e

    type(Atoms) :: myatoms
    integer :: myi, j
    real(dp) :: mydelta, mass
    real(dp), allocatable, dimension(:,:) :: f

    myatoms = at
    myi = optional_default(1, ii)
    mydelta = optional_default(1e-4_dp, delta)

    if (has_property(at, 'mass')) then
       mass = at%mass(myi)
    else
       mass = ElementMass(at%Z(myi))
    end if

    allocate(f(3,at%N))
 
    call set_cutoff(myatoms, cutoff(pot)+0.5_dp)
    call calc_connect(myatoms)

    do j=1,3
       myatoms%pos = at%pos
       myatoms%pos(j,myi) = myatoms%pos(j,myi) + mydelta
       call calc_connect(myatoms)
       call calc(pot, myatoms, force=f, args_str=args_str)
       w_e(j) = sqrt(-f(j,myi)/(mass*mydelta))*ONESECOND
    end do

    deallocate(f)
    
  end function einstein_frequencies
  
  
    subroutine elastic_fields(at, a, C11, C12, C44, Cij)
    type(Atoms), intent(inout) :: at
    real(dp), intent(in) :: a
    real(dp), intent(in), optional :: C11, C12, C44
    real(dp), intent(in), optional :: Cij(6,6)

    real(dp) :: C(6,6), strain(6), stress(6), b
    real(dp), dimension(3) :: n1,n2,n3, d
    real(dp), dimension(3,3) :: rotXYZ, E123, EEt, V, S, R, SS, Sig, RSigRt, RtE, SigEvecs
    integer :: i, j, m, nn, ngood

    real(dp), pointer, dimension(:) :: S_xx_sub1, S_yy_sub1, S_zz_sub1, S_yz, S_xz, S_xy, &
         Sig_xx, Sig_yy, Sig_zz, Sig_yz, Sig_xz, Sig_xy, SigEval1, SigEval2, SigEval3, &
         von_mises_stress, von_mises_strain, atomic_vol, energy_density

    real(dp), pointer, dimension(:,:) :: SigEvec1, SigEvec2, SigEvec3
    logical :: dum
    real(dp), dimension(4) :: dist
    real(dp), dimension(3,4) :: ndiff


    rotXYZ = 0.0_dp
    rotXYZ(1,2) = 1.0_dp
    rotXYZ(2,3) = 1.0_dp
    rotXYZ(3,1) = 1.0_dp

    ! Elastic constants matrix C_{ij}
    if (present(C11) .and. present(C12) .and. present(C44)) then
       C = 0.0_dp
       C(1,1) = C11; C(2,2) = C11; C(3,3) = C11;
       C(4,4) = C44; C(5,5) = C44; C(6,6) = C44;
       C(1,2) = C12; C(1,3) = C12; C(2,3) = C12;
       C(2,1) = C12; C(3,1) = C12; C(3,2) = C12;
    else if (present(Cij)) then
       C = Cij
    else
       call system_abort('elastic_fields: either C11, C12 and C44 (for cubic cell) or full Cij matrix must be present')
    end if


    ! Create properties and assign pointers

    call add_property(at, 'S_xx_sub1', 0.0_dp)
    call add_property(at, 'S_yy_sub1', 0.0_dp)
    call add_property(at, 'S_zz_sub1', 0.0_dp)
    call add_property(at, 'S_yz', 0.0_dp)
    call add_property(at, 'S_xz', 0.0_dp)
    call add_property(at, 'S_xy', 0.0_dp)

    call add_property(at, 'Sig_xx', 0.0_dp)
    call add_property(at, 'Sig_yy', 0.0_dp)
    call add_property(at, 'Sig_zz', 0.0_dp)
    call add_property(at, 'Sig_yz', 0.0_dp)
    call add_property(at, 'Sig_xz', 0.0_dp)
    call add_property(at, 'Sig_xy', 0.0_dp)

    call add_property(at, 'SigEval1', 0.0_dp)
    call add_property(at, 'SigEval2', 0.0_dp)
    call add_property(at, 'SigEval3', 0.0_dp)
    call add_property(at, 'SigEvec1', 0.0_dp, n_cols=3)
    call add_property(at, 'SigEvec2', 0.0_dp, n_cols=3)
    call add_property(at, 'SigEvec3', 0.0_dp, n_cols=3)

    call add_property(at, 'von_mises_stress', 0.0_dp)
    call add_property(at, 'von_mises_strain', 0.0_dp)

    call add_property(at, 'atomic_vol', 0.0_dp)
    call add_property(at, 'energy_density', 0.0_dp)

    dum = assign_pointer(at, 'S_xx_sub1', S_xx_sub1)
    dum = assign_pointer(at, 'S_yy_sub1', S_yy_sub1)
    dum = assign_pointer(at, 'S_zz_sub1', S_zz_sub1)
    dum = assign_pointer(at, 'S_yz', S_yz)
    dum = assign_pointer(at, 'S_xz', S_xz)
    dum = assign_pointer(at, 'S_xy', S_xy)

    dum = assign_pointer(at, 'Sig_xx', Sig_xx)
    dum = assign_pointer(at, 'Sig_yy', Sig_yy)
    dum = assign_pointer(at, 'Sig_zz', Sig_zz)
    dum = assign_pointer(at, 'Sig_yz', Sig_yz)
    dum = assign_pointer(at, 'Sig_xz', Sig_xz)
    dum = assign_pointer(at, 'Sig_xy', Sig_xy)

    dum = assign_pointer(at, 'SigEval1', SigEval1)
    dum = assign_pointer(at, 'SigEval2', SigEval2)
    dum = assign_pointer(at, 'SigEval3', SigEval3)
    dum = assign_pointer(at, 'SigEvec1', SigEvec1)
    dum = assign_pointer(at, 'SigEvec2', SigEvec2)
    dum = assign_pointer(at, 'SigEvec3', SigEvec3)

    dum = assign_pointer(at, 'von_mises_stress', von_mises_stress)
    dum = assign_pointer(at, 'von_mises_strain', von_mises_strain)

    dum = assign_pointer(at, 'atomic_vol', atomic_vol)
    dum = assign_pointer(at, 'energy_density', energy_density)

    call calc_connect(at)

    ! Loop over all atoms
    ngood = 0
    do i=1,at%N

       call print('Atom '//i)

       ! Count number of nearest neighbours
       nn = 0
       do m=1,n_neighbours(at, i) 
          if (is_nearest_neighbour(at,i,m)) then
             nn = nn + 1
             j = neighbour(at, i, m, distance=dist(nn), diff=ndiff(:,nn))

             if (at%Z(j) == 1) then ! Skip hydrogen neighbours
                nn = nn - 1
                cycle
             end if

             call print(nn//': '//dist(nn)//ndiff(:,nn), PRINT_VERBOSE)
             if (nn > 5) exit
          end if
       end do

       if (nn == 4) then

          ngood = ngood + 1

          ! Find cubic axes from neighbours
          n1 = ndiff(:,2)-ndiff(:,1)
          n2 = ndiff(:,3)-ndiff(:,1)
          n3 = ndiff(:,4)-ndiff(:,1)

          call print('n1 '//norm(n1)//n1, PRINT_VERBOSE)
          call print('n2 '//norm(n2)//n2, PRINT_VERBOSE)
          call print('n3 '//norm(n3)//n3, PRINT_VERBOSE)

          ! Estimate volume of Voronoi cell as rhombic dodecahedron with volume 16/9*sqrt(3)*b**3 * 1/2
          b = (norm(ndiff(:,1)) + norm(ndiff(:,2)) + norm(ndiff(:,3)))/3.0_dp
          atomic_vol(i) = 16.0_dp/9.0_dp*sqrt(3.0_dp)*b**3/2.0_dp
          call print('atomic volume '//atomic_vol(i), PRINT_VERBOSE)

          e123(:,1) = (n1 + n2 - n3)/a
          e123(:,2) = (n2 + n3 - n1)/a
          e123(:,3) = (n3 + n1 - n2)/a

          ! Kill near zero elements
          where (abs(e123) < 1.0e-6_dp) e123 = 0.0_dp

          if (all(e123 < 0.0_dp)) e123 = -e123

          call print('det(e) = '//matrix3x3_det(e123), PRINT_VERBOSE)
          call print(e123, PRINT_VERBOSE)

          if (matrix3x3_det(e123) < 0) then
             e123(:,3) = -e123(:,3)
             call print('reflected axis 3', PRINT_VERBOSE)
             call print(e123, PRINT_VERBOSE)
          end if

          ! Find polar decomposition: e123 = S*R where S is symmetric, 
          ! and R is a rotation
          !
          !  EEt = E*E', EEt = VDV' D diagonal, S = V D^1/2 V', R = S^-1*E

          EEt = e123 .mult. transpose(e123) ! Normal
          call diagonalise(EEt, D, V)

          ! Check positive definite
          if (any(D < 0)) then
             call print(e123, PRINT_VERBOSE)
             call system_abort("EE' is not positive definite")
          end if

          S = V .mult. diag(sqrt(D)) .mult. transpose(V)
          R = V .mult. diag(D ** (-0.5_dp)) .mult. transpose(V) .mult. e123

          call print('S:', PRINT_VERBOSE); call print(S, PRINT_VERBOSE)
          call print('R:', PRINT_VERBOSE); call print(R, PRINT_VERBOSE)

          RtE = transpose(R) .mult. e123

          ! Check for permutations - which way does x point?
          if (RtE(2,1) > RtE(1,1) .and. RtE(2,1) > RtE(3,1)) then
             ! y direction
             R = rotXYZ .mult. R
          else if (RtE(3,1) > RtE(1,1) .and. RtE(3,1) > RtE(2,1)) then
             ! z direction
             R = transpose(rotXYZ) .mult. R
          end if

          SS = transpose(R) .mult. S .mult. R

          call print('R:', PRINT_VERBOSE);    call print(R, PRINT_VERBOSE)
          call print('RtE:', PRINT_VERBOSE);  call print(RtE, PRINT_VERBOSE)
          call print('RtSR:', PRINT_VERBOSE); call print(SS, PRINT_VERBOSE)

          ! Strain(1:6) = (/eps11,eps22,eps33,eps23,eps13,eps12/)
          strain(1) = SS(1,1) - 1.0_dp
          strain(2) = SS(2,2) - 1.0_dp
          strain(3) = SS(3,3) - 1.0_dp
          strain(4) = 2.0_dp*SS(2,3)
          strain(5) = 2.0_dp*SS(1,3)
          strain(6) = 2.0_dp*SS(1,2)
       else
          strain = 0.0_dp
          S = 0.0_dp
          S(1,1) = 1.0_dp; S(2,2) = 1.0_dp; S(3,3) = 1.0_dp
       end if

       stress = C .mult. strain

       ! Now stress(1:6) = (/sig11,sig22,sig33,sig23,sig13,sig12/)

       sig = 0.0_dp
       sig(1,1) = stress(1)
       sig(2,2) = stress(2)
       sig(3,3) = stress(3)

       sig(1,2) = stress(6)
       sig(1,3) = stress(5)
       sig(2,3) = stress(4)

       sig(2,1) = stress(6)
       sig(3,1) = stress(5)
       sig(3,2) = stress(4)

       von_mises_stress(i) = sqrt(0.5*((stress(1) - stress(2))**2.0_dp + &
            (stress(2) - stress(3))**2.0_dp + &
            (stress(3) - stress(1))**2.0_dp))

       von_mises_strain(i) = sqrt(strain(6)**2.0_dp + strain(5)**2.0_dp + &
            strain(4)**2.0_dp + &
            1.0_dp/6.0_dp*((strain(2) - strain(3))**2.0_dp + &
            (strain(1) - strain(3))**2.0_dp + &
            (strain(1) - strain(2))**2.0_dp))

       RSigRt = R .mult. sig .mult. transpose(R)

       call symmetrise(RSigRt)
       call diagonalise(RSigRt,D,SigEvecs)

       ! Return strain in such a way that vector (S_xx_sub1, S_yy_sub1, S_zz_sub1, S_yz, S_xz, S_xy)
       ! yields (sig_xx, sig_yy, sig_zz, sig_yz, sig_xz, sig_xy) when multiplied by C_ij matrix
       ! i.e. we return (e1,e2,e3,e4,e5,e6) according to compressed Voigt notation.
       S_xx_sub1(i) = S(1,1) - 1.0_dp
       S_yy_sub1(i) = S(2,2) - 1.0_dp
       S_zz_sub1(i) = S(3,3) - 1.0_dp
       S_yz(i)      = 2.0_dp*S(2,3)
       S_xz(i)      = 2.0_dp*S(1,3)
       S_xy(i)      = 2.0_dp*S(1,2)

       Sig_xx(i)    = RSigRt(1,1)
       Sig_yy(i)    = RSigRt(2,2)
       Sig_zz(i)    = RSigRt(3,3)
       Sig_yz(i)    = RSigRt(2,3)
       Sig_xz(i)    = RSigRt(1,3)
       Sig_xy(i)    = RSigRt(1,2)

       SigEval1(i)  = D(1)
       SigEval2(i)  = D(2)
       SigEval3(i)  = D(3)

       SigEvec1(:,i) = SigEvecs(:,1)
       SigEvec2(:,i) = SigEvecs(:,2)
       SigEvec3(:,i) = SigEvecs(:,3)

       ! energy density in eV/A**3, from u = 1/2*stress*strain
       energy_density(i) = 0.5_dp*( (/Sig_xx(i), Sig_yy(i), Sig_zz(i), Sig_yz(i), Sig_xz(i), Sig_xy(i) /) .dot. &
            (/S_xx_sub1(i), S_yy_sub1(i), S_zz_sub1(i), S_yz(i), S_xz(i), S_xy(i) /) )/GPA

    end do

    call print('Processed '//ngood//' of '//at%N//' atoms.')

  end subroutine elastic_fields


end module elasticity_module
