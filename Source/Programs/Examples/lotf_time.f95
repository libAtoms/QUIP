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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     Learn-on-the-fly (LOTF) hybrid molecular dynamics code
!X
!X     When using this software, the following should be referenced:
!X
!X     Gabor Csanyi, Tristan Albaret, Mike C. Payne and Alessandro De Vita
!X     "Learn on the fly": a hybrid classical and quantum-mechanical
!X         molecular dynamics simulation
!X     Physical Review Letters 93 p. 175503 (2004) >>PDF [626 KB]
!X
!X     Gabor Csanyi, T. Albaret, G. Moras, M. C. Payne, A. De Vita
!X     Multiscale hybrid simulation methods for material systems
!X     J. Phys. Cond. Mat. 17 R691-R703 Topical Review (2005)
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

program lotf_time

  ! 
  ! this is an example program for LOTF which does time-embedding
  ! i.e. fit springs for all atoms, and use these springs to advance the dynamics
  ! for some number of steps in between calling a force model
  !


  use libAtoms_module
  use QUIP_module
  use adjustablepotential_module

  implicit none

  type(Atoms) :: at
  type(Potential) :: pot
  real(dp), allocatable :: f(:,:)
  integer :: i, j, n_cycle, n_extrap
  type(Table) :: fitlist
  type(DynamicalSystem) :: ds, ds_saved
  type(CInOutput) :: movie
  
  
  call system_initialise(PRINT_NORMAL)
  call initialise(pot, 'FilePot command=./castep_driver.sh')
  call initialise(movie, "movie.xyz", action=OUTPUT)

  ! read in initial configuration
  call read(at, "start.xyz")
  call set_cutoff(at, 4.0_dp)
  call calc_connect(at)

  call initialise(ds, at)
  call rescale_velo(ds, 300.0_dp)
  
  ! number of LOTF extrapolation cycles
  n_cycle = 100

  ! number of LOTF extrapolation steps in each cycle
  n_extrap = 10


  allocate(f(3,at%N))


  ! create fitlist - all atoms
  call wipe(fitlist)
  do i=1,ds%atoms%N
     call append(fitlist, (/i, 0, 0, 0/))
  end do


  call print_title('Time embedding')
     
  ! bootstrap the adjustable potential
  call adjustable_potential_init(ds%atoms, fitlist, directionN=ds%atoms%N, &
          method='SVD', nnonly=.true., &
          spring_hops=2, map=.false.)
  call calc(pot, ds%atoms, f=f) 
  call adjustable_potential_optimise(ds%atoms, f, method='SVD')

  ! main loop
  do i = 1,n_cycle

     ! reinitialise LOTF interpolation 'springs'
     call adjustable_potential_init(ds%atoms, fitlist, directionN=ds%atoms%N, &
          method='SVD', nnonly=.true., &
          spring_hops=2, map=.true.)
     
     call print_title('LOTF: Extrapolation')
  
     ds_saved = ds

     ! extrapolation loop
     do j = 1,n_extrap

        ! get force from optimised adjustable potential
        call adjustable_potential_force(ds%atoms, f)

        ! advance the dynamics
        call advance_verlet(ds, 1.0_dp, f)
        call ds_print_status(ds, 'E')

     end do

     call print_title('LOTF: Optimisation')

     ! get new force from force-model (e.g. CASTEP)
     call calc(pot, ds%atoms, f=f)

     call adjustable_potential_optimise(ds%atoms, f, method='SVD')


     call print_title('LOTF: Interpolation')

     ! revert dynamical system
     ds = ds_saved

     ! interpolation loop
     do j = 1,n_extrap

        ! get force from optimised adjustable potential
        call adjustable_potential_force(ds%atoms, f, interp=real(j-1,dp)/real(n_extrap,dp), &
             interp_space=.false., interp_order='linear')

        ! advance the dynamics
        call advance_verlet(ds, 1.0_dp, f)
        call ds_print_status(ds, 'I', instantaneous = .true.)

        ! write movie of configurations
        call write(movie, ds%atoms)
     end do

     ! update connectivity
     call print_title('Connectivity update')
     call calc_connect(ds%atoms)


  end do

end program lotf_time

