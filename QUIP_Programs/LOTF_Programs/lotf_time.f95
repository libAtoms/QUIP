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

subroutine lotf_time(pos, vel, force, Z, cell, timestep, thermo_time, n_extrap_steps)

  use libAtoms_module
  use LOTF_module

  implicit none
  
  real(dp), intent(inout):: pos(:,:), vel(:,:),force(:,:)
  real(dp), intent(in)   :: cell(3,3), timestep, thermo_time
  integer, intent(in)    :: Z(:)
  integer, intent(in)    :: n_extrap_steps
  
  type(DynamicalSystem)::ds
  type(Atoms)::at, di
  type(table)::fitforce
  real(dp), allocatable::f0(:,:), f1(:,:), f(:,:), df(:,:)
  real(dp)::de
  integer:: n, i, Nat
  integer, allocatable::fitlist(:)
  logical, save::firstcall=.true.


  Nat = size(pos,2)
 
  if(firstcall) then
     ! initialise program
     call system_initialise(NORMAL,1) ! use SILENT for no output
     firstcall = .false.
  end if

  ! create some atoms
  call atoms_initialise(at, Nat, cell)
  call atoms_set_cutoff(at, 4.0_dp)
  at%pos = pos
  at%Z = Z
  call calc_connect(at)
  call calc_dists(at)


  ! initialise dynamics
  call ds_initialise(ds, at, vel)

  ! allocate some force arrays
  allocate(f0(3,Nat))
  allocate(f1(3,Nat))
  allocate(f(3,Nat))
  allocate(fitlist(Nat))
  allocate(df(3,Nat))

  ! do something!


  do i=1,Nat
     fitlist(i) = i
  end do


  call print('==================== Time embedding ====================')
     
  ! compute default forces
  call hybrid_force(ds%atoms, f0, SW_Force_noopt)
  ! compute target forces
  f1 = force

  ! copy force differences for the embedding zone 
  df = f1-f0
  call wipe(fitforce)
  call append(fitforce, reshape(fitlist, (/1,Nat/)), df)

  write(line, '(i0,a,i0,a)') Nat, " atoms, fitting on ", Nat, " atoms" 
  call print(line)
  call print("Fitting forces:", VERBOSE)   
  call print(fitforce, VERBOSE)
  call print("", VERBOSE)

  call print('====================  Adjustable potential   ====================')

  ! create and optimise adjustable potential
  call adjustable_potential_init(ds%atoms, fitforce, Nat)
  call adjustable_potential_optimise(ds%atoms, real_part(fitforce))
  
  
  call print('====================        Dynamics         ====================')

  do i = 1, n_extrap_steps

     ! get force from optimised adjustable potential
     call adjustable_potential_force(ds%atoms, df)
     call hybrid_force(ds%atoms, f, SW_Force_noopt)
     f = f+df
        
     ! advance the dynamics
     call advance_verlet(ds, 1.0_dp, f)
     call ds_print_status(ds, 'D')
  end do


  pos = ds%atoms%pos
  vel = ds%atoms%velo

  ! need to compute final forces separately
  call adjustable_potential_force(ds%atoms, df)
  call hybrid_force(ds%atoms, f, SW_Force_noopt)
  force = f+df

  call Finalise(at)
  call Finalise(ds)
  call Finalise(fitforce)
  deallocate(f0,f1,f,df, fitlist)

end subroutine lotf_time

