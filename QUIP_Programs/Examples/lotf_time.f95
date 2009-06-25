!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     Learn-on-the-fly (LOTF) hybrid molecular dynamics code
!X
!X    
!X     Authors: Gabor Csanyi, Alessio Commisso, Steven Winfield
!X     James Kermode, Gianpietro Moras, Michael Payne, Alessandro De Vita
!X
!X     
!X     Copyright 2005, All Rights Reserved 
!X
!X     This source code is confidential, all distribution is 
!X     prohibited. Making unauthorized copies is also prohibited
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

! $Id: lotf_time.f95,v 1.3 2008-04-09 12:46:59 jrk33 Exp $

! $Log: not supported by cvs2svn $
! Revision 1.2  2008/04/09 11:26:56  jrk33
! use LOTF_module
!
! Revision 1.1.1.1  2008/03/13 17:36:35  gc121
! this module contains the LOTF top level programs
!
! Revision 1.2  2006/10/30 12:48:54  gc121
! synced calling conventions
!
! Revision 1.1  2006/10/06 11:41:36  gc121
! subroutine version of the simple time-embedded LOTF extrapolation loop
!

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
  
  
  call system_initialise(NORMAL)
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

