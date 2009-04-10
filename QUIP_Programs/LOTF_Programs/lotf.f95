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

! $Id: lotf.f95,v 1.3 2008-04-22 10:30:51 jrk33 Exp $

! $Log: not supported by cvs2svn $
! Revision 1.2  2008/03/18 14:26:57  jrk33
! at%Z = 14 -> call set_atoms(at, 14)
!
! Revision 1.1.1.1  2008/03/13 17:36:35  gc121
! this module contains the LOTF top level programs
!
! Revision 1.21  2007/11/07 11:33:34  gc121
! fitlist was in the wrong place
!
! Revision 1.20  2007/11/07 10:49:43  gc121
! updated style
!
! Revision 1.19  2007/11/07 10:47:56  gc121
! updated to work with rest of code
!
! Revision 1.18  2007/07/02 17:31:40  gc121
! updated to reflect other changes
!
! Revision 1.17  2006/12/12 15:19:42  gc121
! sync
!
! Revision 1.16  2006/10/30 12:49:28  gc121
! synced calling conventions
!
! Revision 1.15  2006/06/29 17:55:15  gc121
! updated functionames
!
! Revision 1.14  2006/06/21 13:45:50  gc121
! added mandatory classical_force argument to hybrid_force() and overloaded SW routines to enable ignoring unknown atom types on request
!
! Revision 1.13  2006/06/20 17:23:19  gc121
! added new copyright notice to include James, Gian, Mike and Alessandro
!



program lotf
  use libAtoms_module
  use LOTF_module
!  use test_force_model_module
  

  implicit none
  
  type(DynamicalSystem)::ds
  type(Atoms)::at, dia
  type(inoutput)::movie
  type(table)::embedlist, fitlist, embedforce, fitforce
  real(dp), allocatable::f0(:,:), f1(:,:), f(:,:), df(:,:)
  real(dp)::de
  integer:: n

  ! initialise program
  call system_initialise(NORMAL,1)
  !call initialise(movie, "movie.xyz")

  ! create some atoms
  call diamond(dia, 5.44_dp)
  call supercell(at, dia, 3,3,3)
  call set_atoms(at, 14)
  call atoms_set_cutoff(at, 4.0_dp)
  call randomise(at%pos, 0.01_dp)
  call calc_connect(at)
  !call test_force_model(at)

  !call print(at, verbosity=VERBOSE)

  ! create list of embedded atoms
  call append(embedlist, (/1,0,0,0/))
  call BFS_grow(at, embedlist, 2)
  !call print(embedlist)


  ! allocate some force arrays
  allocate(f0(3,at%N))
  allocate(f1(3,at%N))
  allocate(f(3,at%N))

  ! initialise dynamics
  call ds_initialise(ds, at)
  call rescale_velo(ds, 1000.0_dp)
  call zero_momentum(ds)
  !call print(ds)


  ! grow the embed list to include a buffer zone
  fitlist = embedlist
  call BFS_grow(ds%atoms, fitlist, 3)
  call reallocate(df, 3, fitlist%N, zero=.true.)

  ! do something!

  do n=1,10

     ! print movie and recompute connectivity
     if(mod(n, 5) == 1) then
        call calc_connect(ds%atoms)
        !call print_xyz(ds%atoms, movie)
     end if

     call print('==================== Simple hybrid embedding ====================')
     
     ! compute default forces
     call hybrid_force(ds%atoms, f0, SW_Force_noopt)
     ! compute hybrid forces
     call hybrid_force(ds%atoms, f1, embedlist, SW_Force_noopt, SWeps_Force, little_clusters = .false. )

     df = 0.0_dp
     ! copy force differences for the embedding zone 
     df(:,1:embedlist%N) = f1(:,int_part(embedlist,1))-f0(:,int_part(embedlist,1))
     call wipe(fitforce)
     call append(fitforce, intpart_2D=int_part(fitlist), realpart_2D=df)

     call print("Embedding "//embedlist%N//" atoms, fitting on "//fitlist%N//" atoms") 

     call print("Fitting forces:", VERBOSE)   
     call print(fitforce, VERBOSE)
     call print("", VERBOSE)

     call print('====================  Adjustable potential   ====================')

     ! create and optimise adjustable potential
     call adjustable_potential_init(ds%atoms, fitlist, embedlist%N)
     call adjustable_potential_optimise(ds%atoms, real_part(fitforce))

     ! get force from optimised adjustable potential
     call adjustable_potential_force(ds%atoms, df)

     ! add it onto default forces
     f = f0
     f(:,int_part(fitforce,1)) = f(:,int_part(fitforce,1))+df

     call print('====================        Dynamics         ====================')

     ! advance the dynamics
     call advance_verlet(ds, 1.0_dp, f)
     call ds_print_status(ds, 'D')
  end do
  

  call Atoms_Finalise(at)
  call Atoms_Finalise(dia)
  call DS_Finalise(ds)
  !call Free(movie)
  call finalise(embedlist)
  call finalise(fitlist)
  call finalise(embedforce)
  call finalise(fitforce)
  deallocate(f0,f1,f,df)

  call system_finalise()

end program lotf

