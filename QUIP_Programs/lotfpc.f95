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

! $Id: lotfpc.f95,v 1.4 2008-04-09 11:23:36 jrk33 Exp $

! $Log: not supported by cvs2svn $
! Revision 1.3  2008/04/09 11:14:23  jrk33
! added keyword args to hybrid_force calls
!
! Revision 1.2  2008/04/09 10:56:52  jrk33
! lotf.h has gone!
!
! Revision 1.1.1.1  2008/03/13 17:36:35  gc121
! this module contains the LOTF top level programs
!
! Revision 1.28  2007/07/13 14:23:40  jrk33
! assign_pointer now a function
!
! Revision 1.27  2007/07/02 17:49:24  gc121
! updated to reflect other changes
!
! Revision 1.26  2007/07/02 17:31:40  gc121
! updated to reflect other changes
!
! Revision 1.25  2007/04/11 08:52:08  saw44
! Added usage of verbosity stack
!
! Revision 1.24  2007/03/28 17:35:24  saw44
! Cleaned up unused variables
!
! Revision 1.23  2006/12/12 15:19:42  gc121
! sync
!
! Revision 1.22  2006/10/30 20:47:59  gc121
! sync
!
! Revision 1.21  2006/10/30 12:48:20  gc121
! sync
!
! Revision 1.20  2006/10/30 12:30:02  gc121
! synced calling convention changes with adjustable potential
!
! Revision 1.19  2006/10/12 08:31:50  jrk33
! Converted all force calls from functions to subroutines
!
! Revision 1.18  2006/06/21 13:45:50  gc121
! added mandatory classical_force argument to hybrid_force() and overloaded SW routines to enable ignoring unknown atom types on request
!
! Revision 1.17  2006/06/20 17:23:19  gc121
! added new copyright notice to include James, Gian, Mike and Alessandro
!


program lotfpc

  use libAtoms_module
  use LOTF_module

  implicit none

  ! Parameters for the program (we should probably put a nice command line argument
  ! reader in System...)
 
  logical,  parameter :: do_lots_of_little_clusters = .false.
  integer,  parameter :: n_pc_iter = 1          ! number of predictor/corrector iterations in each cycle
  integer,  parameter :: n_interp_steps = 10    ! number of steps to extrapolate/interpolate for
  integer,  parameter :: n_lotf_cycles = 1      ! number of extrap/interp cycles to do
  real(dp), parameter :: dt = 1.0_dp  ! time step
  real(dp), parameter :: init_temp = 2000.0_dp
  integer, pointer ::qm(:)
  
  !local variables
  type(DynamicalSystem)::ds, ds_saved
  type(Atoms)::at, dia
  type(inoutput)::movie, datafile
  type(table)::embedlist, fitlist, fitforce
  real(dp), allocatable::f0(:,:), f1(:,:), f(:,:), df(:,:), f_hyb(:,:), data(:,:)
  integer::i, j, n
  character(5)::prefix

  ! initialise program
  call system_initialise(NORMAL,1)
  !call initialise(movie, "movie.xyz")
  !call initialise(datafile, 'forcediff.dat')

  ! create some atoms
  call diamond(dia, 5.44_dp)
  call supercell(at, dia, 3,3,3)
  at%Z = 14
  call atoms_set_cutoff(at, 4.0_dp)
  call randomise(at%pos, 0.01_dp)
  call calc_connect(at)


  ! allocate some force arrays
  allocate( f0(3,at%N), f1(3,at%N), f(3,at%N), f_hyb(3,at%N) )

  ! initialise dynamics
  call ds_initialise(ds, at)
  call rescale_velo(ds, init_temp)
  call zero_momentum(ds)

  !allocate space for the test data array
  allocate(data(3,n_interp_steps+1))



  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! X
  ! X  bootstrap the adjustable potential
  ! X
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  ! create list of embedded atoms
  call append(embedlist, (/1,0,0,0/))
  call BFS_grow(ds%atoms, embedlist, 2)

  call Add_Property(ds%atoms, 'QM', 0)
  if (.not. assign_pointer(ds%atoms, 'QM', qm)) &
       call system_abort('error assigning pointer to QM property')
  qm(int_part(embedlist,1)) = 1

  ! compute default forces
  call hybrid_force(ds%atoms, f0, SW_Force_noopt) ! with no embedded atoms


  ! compute hybrid forces
  call hybrid_force(ds%atoms, f1, embedlist, classical_force=SW_force_noopt, qm_force=SWeps_force, little_clusters = do_lots_of_little_clusters, terminate = .false., periodic_clusters = (/.true., .true., .true./), buf_grow_hops = 10)  

  ! grow the embed list to include a fit zone
  fitlist = embedlist
  call BFS_grow(ds%atoms, fitlist, 3, min_images_only = .true.)

  ! copy force differences for the embedding zone... 
  call reallocate(df, 3, fitlist%N, zero=.true.)
  df(:,1:embedlist%N) = f1(:,int_part(embedlist,1))-f0(:,int_part(embedlist,1))
  call wipe(fitforce)
  call append(fitforce, intpart_2d=int_part(fitlist), realpart_2d=df) 

  call adjustable_potential_init(ds%atoms, fitlist, embedlist%N)
  call adjustable_potential_optimise(ds%atoms, real_part(fitforce))

  !call wipe(fitlist)
  call wipe(fitforce)

  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! X
  ! X  main LOTF loop
  ! X
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  ! Operation                                         parameters   saved parameters
  !
  !                                                         P1          P0
  ! set up new structure for parameters                   
  !                                                         P1          P1
  ! Extrapolate  using P1
  !                                                         P1          P1
  ! Optimise parameters, get P2                            
  !                                                         P2          P1
  ! Interpolate
  !                                                         P2          P1
  ! calcConnect
  !                                                         P2          P1
  
  do n=1,n_lotf_cycles

     call Print('')
     call print('====================   Quantum Selection     ====================')
     call Print('')

     ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     ! X
     ! X Code could go here to set up new embedlist
     ! X
     ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     write(line,'(i0,a)') embedlist%N,' atoms flagged for quantum treatment'
     call Print(line)

     !call Print('Building fit zone...')
     !fitlist = embedlist
     !call BFS_grow(ds%atoms, fitlist, 3, min_images_only = .true.)
     
     write(line,'(a,i0,a)') 'Fitting on ',fitlist%N,' atoms in total'
     call Print(line)
     
     call adjustable_potential_init(ds%atoms, fitlist, fitlist%N, map=.true., spring_hops = 2, method='SVD')

     call Print('')
     call print('====================     Extrapolation     ====================')
     call Print('')



     !only reallocate df if needed
     call reallocate(df,3,fitlist%N)

     ds_saved = ds

     do i = 1, n_interp_steps

        ! get force from optimised adjustable potential
        call adjustable_potential_force(ds%atoms, df)

        ! add it onto default forces
        call hybrid_force(ds%atoms, f, SW_Force_noopt)
        f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1)) + df

        ! ********** TESTING ***********
        ! Get the 'actual' force
        !call hybrid_force(ds%atoms, f_hyb, embedlist, SW_Force_noopt, SWeps_Force,do_lots_of_little_clusters)
        !data(1,i) = ds%t              ! time
        !data(2,i) = RMS_diff(f_hyb,f,fitlist) ! RMS force differences
        ! ******************************

        ! advance the dynamics
        call advance_verlet(ds, dt, f)
        call ds_print_status(ds, 'E')

     end do


     do j = 1, n_pc_iter
        call Print('')
        call print('====================   Computation of forces   ====================')
        call Print('')
        
        call Print('Computing new forces')
        ! compute default force
        call hybrid_force(ds%atoms, f0, SW_Force_noopt)
        ! compute hybrid forces
        call hybrid_force(ds%atoms, f1, embedlist, classical_force=SW_Force_noopt, qm_force=SWeps_force, little_clusters = do_lots_of_little_clusters, terminate = .false., periodic_clusters = (/.false., .false., .false./), buf_grow_hops = 10, randomise_buffer = .true.)  
        call print('f0:')
        call print(f0)
        call print('f1:')
        call print(f1)

        call Print('')
        call print('==================== Optimisation of parameters ====================')
        call Print('')
        
        ! copy new force differences for the same embedding zone as before...
        df = 0.0_dp
        df(:,1:embedlist%N) = f1(:,int_part(embedlist,1))-f0(:,int_part(embedlist,1))
        call wipe(fitforce)
        call append(fitforce, intpart_2d=int_part(fitlist), realpart_2d=df)
        
        ! now optimise the parameters for these forces at the new positions
        call adjustable_potential_optimise(ds%atoms, real_part(fitforce))
        ! ********** TESTING ***********
        !call adjustable_potential_force(ds%atoms, df)
        !f = f0;
        !f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1)) + df
        !data(1,i) = ds%t              ! time
        !data(2,i) = RMS_diff(f1,f,fitlist) ! RMS force differences
        !data(3,i) = RMS_diff(real_part(fitforce), df)
        ! ******************************


        call Print('')
        call print('====================        Interpolation        ====================')
        call Print('')
        
        ! revert to the saved dynamical system
        ds = ds_saved
        
        do i = 1, n_interp_steps
           
           ! get force differences from interpolated parameters
           call adjustable_potential_force(ds%atoms, df, interp=real(i-1,dp)/real(n_interp_steps,dp))
           !call adjustable_potential_force(ds%atoms, df, interp_space=.true.)
           
           ! add it onto default forces
           call hybrid_force(ds%atoms, f, SW_Force_noopt) 
           f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1))+df
           
           !*********** TESTING *************
           !get the actual force
           !call hybrid_force(ds%atoms, f_hyb, embedlist, classical_force=SW_Force_noopt, qm_force=SWeps_Force,do_lots_of_little_clusters)
           !data(3,i) = RMS_diff(f_hyb,f,fitlist)
           !*********************************
           
           ! advance the dynamics
           call advance_verlet(ds, dt, f)
           write(prefix, '(a1,i1)') 'I', j
           call ds_print_status(ds, prefix)
           
        end do
     end do

     !Write the data to the data file
     !do i = 1, n_interp_steps+1
     !   write(line,'(3(f0.18,tr2))') data(:,i)
     !   call Print(datafile,line)
     !end do

     call Print('')
     call print('==================== Recalculate Connectivity ====================')   
     call Print('')

     ! print movie and recompute connectivity
     call calc_connect(ds%atoms)
     call print_xyz(ds%atoms, movie, properties = 'pos:QM')

  end do
  

  call adjustable_potential_finalise
  call Atoms_Finalise(at)
  call Atoms_Finalise(dia)
  call DS_Finalise(ds)
  call DS_Finalise(ds_saved)
  !call finalise(movie)
  call finalise(datafile)
  call finalise(embedlist)
  call finalise(fitlist)
  call finalise(fitforce)
  deallocate(f0,f1,f,f_hyb,df)

  call system_finalise()

end program lotfpc

