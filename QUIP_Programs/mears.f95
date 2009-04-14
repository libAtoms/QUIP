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

! $Id: mears.f95,v 1.2 2008-04-09 11:14:52 jrk33 Exp $

! $Log: not supported by cvs2svn $
! Revision 1.1.1.1  2008/03/13 17:36:35  gc121
! this module contains the LOTF top level programs
!


program mears


  use libAtoms_module
  use hybrid_force_module
  use AdjustablePotential_module
  use FileQM_module
  use Watanabe_module

  implicit none

  ! Parameters for the program (we should probably put a nice command line argument
  ! reader in System...)
 
  logical,  parameter :: do_lots_of_little_clusters = .false.
  integer,  parameter :: n_interp_steps = 10   ! number of steps to extrapolate/interpolate for
  integer,  parameter :: n_lotf_cycles = 10 ! number of extrap/interp cycles to do
  integer,  parameter :: buffer_size = 3
  real(dp), parameter :: dt = 1.0_dp  ! time step
  real(dp), parameter :: init_temp = 300.0_dp

  !local variables
  type(DynamicalSystem)::ds, ds_saved
  type(Atoms)::at, dia
  type(inoutput)::movie, datafile
  type(table)::embedlist, fitlist, fitforce, opt_table
  real(dp), allocatable::f0(:,:), f1(:,:), f(:,:), df(:,:), f_hyb(:,:), data(:,:)
  integer::i, j, k, n
  integer, pointer ::qm(:)
  logical :: dummy

  ! initialise program
  call system_initialise
  call initialise(movie, "movie.xyz")
  !call initialise(datafile, 'forcediff.dat')



  ! create some atoms

  call read_xyz(at, "start.xyz")
  call atoms_set_cutoff(at, 5.0_dp)
  call randomise(at%pos, 0.001_dp)
  call calc_connect(at)

  ! allocate some force arrays
  allocate( f0(3,at%N), f1(3,at%N), f(3,at%N), f_hyb(3,at%N) )


  ! initialise dynamics
  call ds_initialise(ds, at)
  call rescale_velo(ds, init_temp)
  call zero_momentum(ds)
  ds%thermostat = LANGEVIN_THERM
  ds%sim_temp = init_temp
  ds%thermal_tau = 100.0_dp

  !allocate space for the test data array
  allocate(data(3,n_interp_steps))

  ! set the QM command
  call FileQM_Set_Command('./castep_driver')

  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! X
  ! X  bootstrap the adjustable potential
  ! X
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  ! create list of embedded atoms
  call append(embedlist, (/419,0,0,0/)) ! just picked an atom 
  call BFS_grow(at, embedlist, 1)

  ! add the QM property to the atoms so that we can produce a nice movie later on
  call Add_Property(ds%atoms, 'QM', 0)
  dummy = assign_pointer(ds%atoms, 'QM', qm)
  qm(int_part(embedlist,1)) = 1
  
  ! compute default forces
  call hybrid_force(ds%atoms, f0, Watanabe_force) ! with no embedded atoms
  
  ! compute hybrid forces
  call hybrid_force(ds%atoms, f1, embedlist, classical_force=Watanabe_force, qm_force=FileQM_Force, little_clusters = do_lots_of_little_clusters, buf_grow_hops = buffer_size) 

  ! grow the embed list to include a fit zone
  fitlist = embedlist
  call BFS_grow(ds%atoms, fitlist, 3)

  ! copy force differences for the embedding zone... 
  call reallocate(df, 3, fitlist%N, zero=.true.)
  df(:,1:embedlist%N) = f1(:,int_part(embedlist,1))-f0(:,int_part(embedlist,1))
  call wipe(fitforce)
  call append(fitforce, int_part(fitlist), df) 
    
  ! initialise and optimise the parameters
  call adjustable_potential_init(ds%atoms, fitlist, embedlist%N)  ! no mapping here
  call adjustable_potential_optimise(ds%atoms, real_part(fitforce))

  call wipe(fitlist)
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
     ! X Code should go here to set up embedlist
     ! X
     ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     write(line,'(i0,a)') embedlist%N,' atoms flagged for quantum treatment'
     call Print(line)

     call Print('Building fit zone...')
     fitlist = embedlist
     call BFS_grow(ds%atoms, fitlist, 3)
     
     write(line,'(a,i0,a)') 'Fitting on ',fitlist%N,' atoms in total'
     call Print(line)

     call adjustable_potential_init(ds%atoms, fitlist, embedlist%N, map = .true.)

     call Print('')
     call print('====================     Extrapolation     ====================')
     call Print('')

     ds_saved = ds

     !only reallocate df if needed
     call reallocate(df,3,fitlist%N)

     do i = 1, n_interp_steps

        ! get force from optimised adjustable potential
        call adjustable_potential_force(ds%atoms, df)

        ! add it onto default forces
        call hybrid_force(ds%atoms, f, Watanabe_force)
        f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1)) + df

        ! ********** TESTING ***********
        ! Get the 'actual' force
        !call hybrid_force(ds%atoms, f_hyb, embedlist, Watanabe_force, FileQM_Force,do_lots_of_little_clusters)
        !data(1,i) = ds%t              ! time
        !data(2,i) = RMS_diff(f_hyb,f,fitlist) ! RMS force differences
        ! ******************************

        ! advance the dynamics
        call advance_verlet(ds, dt, f)
        call ds_print_status(ds, 'E')

     end do


     call Print('')
     call print('====================   Computation of forces   ====================')
     call Print('')

     call Print('Computing new forces')
     !compute default force
     call hybrid_force(ds%atoms, f0, Watanabe_force)
     ! compute hybrid forces
     call hybrid_force(ds%atoms, f1, embedlist, Watanabe_force, FileQM_Force, little_clusters = do_lots_of_little_clusters, buf_grow_hops = buffer_size)     

     call Print('')
     call print('==================== Optimisation of parameters ====================')
     call Print('')
     
     ! copy new force differences for the same embedding zone as before...
     df = 0.0_dp
     df(:,1:embedlist%N) = f1(:,int_part(embedlist,1))-f0(:,int_part(embedlist,1))
     call wipe(fitforce)
     call append(fitforce, int_part(fitlist), df)

     ! now optimise the parameters for these forces at the new positions
     call adjustable_potential_optimise(ds%atoms, real_part(fitforce))

     call Print('')
     call print('====================        Interpolation        ====================')
     call Print('')

     ! revert to the save dynamical system
     ds = ds_saved

     do i = 1, n_interp_steps
        
        ! get force differences from interpolated parameters
        call adjustable_potential_force(ds%atoms, df, interp = real(i-1,dp)/real(n_interp_steps,dp))

        ! add it onto default forces
        call hybrid_force(ds%atoms, f, Watanabe_force) 
        f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1))+df

        !*********** TESTING *************
        !get the actual force
        !call hybrid_force(ds%atoms, f_hyb, embedlist, Watanabe_force, FileQM_Force,do_lots_of_little_clusters,buffer_size)
        !data(3,i) = RMS_diff(f_hyb,f,fitlist)
        !*********************************

        ! advance the dynamics
        call advance_verlet(ds, dt, f)
        call ds_print_status(ds, 'I')

     end do

     !Write the data to the data file
     !do i = 1, n_interp_steps
     !   write(line,'(3(f0.18,tr2))') data(:,i)
     !   call Print(datafile,line)
     !end do

     call Print('')
     call print('==================== Recalculate Connectivity ====================')   
     call Print('')

     ! print movie and recompute connectivity
     call calc_connect(ds%atoms)
     call print_xyz(ds%atoms, movie, properties='pos:QM')

  end do
  

  call adjustable_potential_finalise
  call Atoms_Finalise(at)
  call Atoms_Finalise(dia)
  call DS_Finalise(ds)
  call DS_Finalise(ds_saved)
  call finalise(movie)
  !call finalise(datafile)
  call finalise(embedlist)
  call finalise(fitlist)
  call finalise(fitforce)
  deallocate(f0,f1,f,f_hyb,df)

  call system_finalise()

end program mears

