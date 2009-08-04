module elasticity_module

  use libatoms_module

  use QUIP_module
  use MetaPotential_module
  implicit none

  private

  public :: calc_elastic_constants
  interface calc_elastic_constants
     module procedure metapot_calc_elastic_constants
  end interface

  public :: youngs_modulus, poisson_ratio, graphene_elastic, einstein_frequencies


contains

  subroutine metapot_calc_elastic_constants(this, at, fd, args_str, c, c0, relax_initial, return_relaxed, relax_tol)
    type(MetaPotential), intent(inout) :: this
    type(Atoms), intent(inout) :: at !% Atoms object for which to compute $C_{ij}$
    real(dp), intent(in), optional :: fd !% Finite strain to apply. Default $10^{-3}$.
    character(len=*), intent(in), optional :: args_str !% Optional args_str to pass to 'minim'
    real(dp), intent(out), optional :: c(6,6) !% Elastic constants (with relaxation)
    real(dp), intent(out), optional :: c0(6,6) !% Elastic constants (without relaxation)
    logical, optional :: relax_initial !% Should the initial cell be relaxed?
    logical, optional :: return_relaxed !% If true, overwrite 'at' with relaxed positions and lattice (default false)
    real(dp), optional :: relax_tol !% Relaxation df^2 tolerance. Default 1e-8

    logical my_relax_initial, my_return_relaxed
    integer ii, jj
    type(atoms) :: at_bulk
    type(atoms) :: at_t
    real(dp) :: my_fd, my_relax_tol
    real(dp) :: volume
    real(dp) :: Fp(3,3), Fm(3,3)
    real(dp) :: V0p(3,3), V0m(3,3), Vp(3,3), Vm(3,3)
    integer iter

    real(dp) :: e0, V0(3,3)

    my_fd = optional_default(1.0e-3_dp, fd)
    my_relax_tol = optional_default(1e-8_dp, relax_tol)
    my_relax_initial = optional_default(.true., relax_initial)
    my_return_relaxed = optional_default(.false., return_relaxed)

    at_bulk = at
    call calc_connect(at_bulk)

    if (current_verbosity() > VERBOSE) then
       call verbosity_push_decrement(NERD)
       call calc(this, at_bulk, e=e0, virial=v0)
       call verbosity_pop()
       call print("initial config")
       call print_xyz(at_bulk, mainlog)
       call print("e " // e0)
       call print("V")
       call print(V0)
    endif

    if (my_relax_initial) then
       call verbosity_push_decrement(VERBOSE)
       iter = minim(this, at_bulk, 'cg', my_relax_tol, 1000, 'FAST_LINMIN', do_print=.false., &
            do_pos=.true., do_lat=.true., args_str=args_str, use_n_minim=.true.)
       call verbosity_pop()

       if (current_verbosity() >= VERBOSE) then
          call verbosity_push_decrement(VERBOSE)
          call calc(this, at_bulk, e=e0, virial=v0)
          call verbosity_pop()
          call print("relaxed config")
          call print_xyz(at_bulk, mainlog)
          call print("e " // e0)
          call print("V")
          call print(V0)
       endif

       if (my_return_relaxed) at = at_bulk
    endif

    volume = cell_volume(at_bulk)

    call print("volume " // volume, VERBOSE)

    do ii=1, 3
       do jj=ii, 3
          call print("doing ii " // ii // " jj " // jj, VERBOSE)
          Fp = 0.0_dp; call add_identity(Fp)
          Fm = 0.0_dp; call add_identity(Fm)

          Fp(ii,jj) = Fp(ii,jj) + my_fd
          Fm(ii,jj) = Fm(ii,jj) - my_fd

          at_t = at_bulk
          call set_lattice(at_t, Fp .mult. at_t%lattice)
          at_t%pos = Fp .mult. at_t%pos
          call calc_connect(at_t)
          call calc(this, at_t, e=e0, virial=V0p)
          call verbosity_push_decrement(VERBOSE)
          call print("plus perturbed config")
          call print_xyz(at_t, mainlog)
          call print("E0p" // e0)
          call print("V0p")
          call print(V0p)
          call verbosity_pop()

          if (present(c)) then
             call verbosity_push_decrement(VERBOSE)
             iter = minim(this, at_t, 'cg', my_relax_tol, 1000, 'FAST_LINMIN', do_print=.false., &
                  do_pos=.true., do_lat=.false., args_str=args_str, use_n_minim=.true.)
             call verbosity_pop()

             call calc_connect(at_t)
             call calc(this, at_t, e=e0, virial=Vp)
             call verbosity_push_decrement(VERBOSE)
             call print("plus perturbed relaxed config")
             call print_xyz(at_t, mainlog)
             call print("Ep" // e0)
             call print("Vp")
             call print(Vp)
             call verbosity_pop()
          endif

          at_t = at_bulk
          call set_lattice(at_t, Fm .mult. at_t%lattice)
          at_t%pos = Fm .mult. at_t%pos
          call calc_connect(at_t)
          call calc(this, at_t, e=e0, virial=V0m)
          call verbosity_push_decrement(VERBOSE)
          call print("minus perturbed config")
          call print_xyz(at_t, mainlog)
          call print("E0m" // e0)
          call print("V0m")
          call print(V0m)
          call verbosity_pop()

          if (present(c)) then
             call verbosity_push_decrement(VERBOSE)
             iter = minim(this, at_t, 'cg', my_relax_tol, 1000, 'FAST_LINMIN', do_print=.false., &
                  do_pos=.true., do_lat=.false., args_str=args_str, use_n_minim=.true.)
             call verbosity_pop()

             call calc_connect(at_t)
             call calc(this, at_t, e=e0, virial=Vm)
             call verbosity_push_decrement(VERBOSE)
             call print("minus perturbed relaxed config")
             call print_xyz(at_t, mainlog)
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

  end subroutine metapot_calc_elastic_constants

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
  function poisson_ratio(c, l, m) result(v)
    real(dp), intent(in) :: c(6,6)
    real(dp), dimension(3), intent(in) :: l, m
    real(dp) :: v

    real(dp) :: s(6,6)
    real(dp), dimension(3) :: lhat, mhat

    call inverse(c, s)

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
  !% parameter 'a' using the MetaPotential 'metapot'. On exit, 'poisson'
  !% will contain the in plane poisson ratio (dimensionless) and 'young' the
  !% in plane Young's modulus (GPa).
  subroutine Graphene_Elastic(metapot, a, poisson, young, cb)
    type(MetaPotential), intent(inout) :: metapot
    real(dp), intent(out) :: a, poisson, young
    real(dp), intent(out), optional :: cb

    type(Atoms) :: cube, at, at2, tube
    real(dp) :: a0, v(3,3), eps(3,3)
    integer :: steps, tube_i, i, n, m
    logical :: fix(3,3)
    integer, parameter :: nx = 2, ny = 3
    real(dp) :: tube_r(12), tube_energy(12), graphene_e_per_atom, radius, energy, c(1), chisq
    
    a0 = 1.45_dp
    cube = Graphene_Cubic(a0)

    call Supercell(at, cube, nx, ny, 1)
    call Atoms_Set_Cutoff(at, cutoff(metapot))
    call randomise(at%pos, 0.01_dp)
    call calc_connect(at)

    fix = .true.
    fix(1,1) = .false.
    fix(2,2) = .false.

    ! Geometry optimise with variable lattice
    steps = minim(metapot, at, 'cg', 1e-6_dp, 100, &
         'FAST_LINMIN', do_pos=.true.,do_lat=.true.,lattice_fix=fix, do_print=.false.)

    ! Set a to average of x and y lattice constants
    a = 0.5_dp*(at%lattice(1,1)/(3.0_dp*nx) + at%lattice(2,2)/(sqrt(3.0_dp)*ny))

    call calc(metapot, at, e=graphene_e_per_atom)
    graphene_e_per_atom = graphene_e_per_atom/at%N

    cube = Graphene_Cubic(a)
    call Supercell(at2, cube, nx, ny, 1)
    call Atoms_Set_Cutoff(at2, cutoff(metapot))
    call calc_connect(at2)

    ! Apply small strain in x direction
    eps = 0.0_dp; call add_identity(eps)
    eps(1,1) = eps(1,1)+0.001_dp
    call set_lattice(at2, eps .mult. at2%lattice)
    at2%pos = eps .mult. at2%pos

    fix = .false.
    fix(1,1) = .true.
    ! Fix lattice in x direction

    ! Geometry optimse, allowing to contract in y direction
    steps = minim(metapot, at2, 'cg', 1e-6_dp, 100, &
         'FAST_LINMIN', do_print=.false., do_pos=.true.,do_lat=.true.,lattice_fix=fix)

    poisson = -((at2%lattice(2,2) - at%lattice(2,2))/at%lattice(2,2))/ &
         ((at2%lattice(1,1) - at%lattice(1,1))/at%lattice(1,1))

    ! Calculate stress to find Young's modulus
    call Calc(metapot, at2, virial=v)

    young = -v(1,1)/(eps(1,1)-1.0_dp)*GPA/cell_volume(at2)

    if (present(cb)) then
       ! Finally, calculate bending modulus by fitting strain energy per atom to
       ! curve 1/2 c_b/r**2 for a series of nanotubes. We use 12 nanotubes (n,m) 
       ! where n runs from 10 to 20 in steps of 2 and m is either 0 or n.
       tube_i = 1
       tube_r = 0.0_dp
       tube_energy = 0.0_dp
       call verbosity_push_decrement(VERBOSE)
       do n=10,20,2
          do m=0,n,n
             radius = graphene_tube(tube, a, n, m, 3)
             call set_cutoff(tube, cutoff(metapot))
             call calc_connect(tube)

             
             i = minim(metapot, tube, 'cg', 1e-5_dp, 1000, 'FAST_LINMIN', do_print=.false., &
                  do_pos=.true., do_lat=.false.)
             
             call calc(metapot, tube, e=energy)

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


  function einstein_frequencies(metapot, at, i, delta) result(w_e)
    type(MetaPotential), intent(inout) :: metapot  !% MetaPotential to use
    type(Atoms), intent(in) :: at            !% Atoms structure - should be equilibrium bulk configuation
    integer, optional, intent(in) :: i       !% The atom to displace (default 1)
    real(dp), optional, intent(in) :: delta  !% How much to displace it (default 1e-4_dp)
    real(dp), dimension(3) :: w_e

    type(Atoms) :: myatoms
    integer :: myi, j
    real(dp) :: mydelta, mass
    real(dp), allocatable, dimension(:,:) :: f

    myatoms = at
    myi = optional_default(1, i)
    mydelta = optional_default(1e-4_dp, delta)

    if (has_property(at, 'mass')) then
       mass = at%mass(myi)
    else
       mass = ElementMass(at%Z(myi))
    end if

    allocate(f(3,at%N))
 
    call set_cutoff(myatoms, cutoff(metapot)+0.5_dp)
    call calc_connect(myatoms)

    do j=1,3
       myatoms%pos = at%pos
       myatoms%pos(j,myi) = myatoms%pos(j,myi) + mydelta
       call calc_connect(myatoms)
       call calc(metapot, myatoms, f=f)
       w_e(j) = sqrt(-f(j,myi)/(mass*mydelta))*ONESECOND
    end do

    deallocate(f)
    
  end function einstein_frequencies
  
  
  

end module elasticity_module
