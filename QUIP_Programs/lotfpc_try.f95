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




program lotfpc

  use libAtoms_module
  use LOTF_module

  implicit none

  ! Parameters
  integer :: n_pc_iter, n_interp_steps, n_lotf_cycles, seed, refit, thermostat
  integer :: single_cluster, full_qm, classical, grow_hops, fit_grow_hops, spring_hops
  real(dp)  :: dt, init_temp, sim_temp, nneightol, inject_heat

  integer, pointer ::qm(:)
  
  !local variables
  type(DynamicalSystem)::ds, ds_saved
  type(Atoms)::at, dia
  type(inoutput):: movie, datafile, qmtrajectory
  type(table)::embedlist, fitlist, fitforce, twobody_save, threebody_save
  real(dp), allocatable::f0(:,:), f1(:,:), f(:,:), df(:,:), f_hyb(:,:), data(:,:), qm_pos(:,:,:), fdiff(:,:)
  real(dp) :: energy, dV_dt, Epot
  integer::i, j, n
  logical :: do_periodic(3)

  ! initialise program
  call system_initialise(NORMAL,1)
  !call initialise(movie, "movie.xyz")
  !call initialise(datafile, 'forcediff.dat')

  ! create some atoms
  call diamond(dia, 5.44_dp)
  call supercell(at, dia, 3,3,3)
  at%Z = 14

  !call randomise(at%pos, .1_dp)

  call atoms_set_cutoff(at, 4.0_dp)
  call randomise(at%pos, 0.01_dp)
  call calc_connect(at)

  ! allocate some force arrays
  allocate( f0(3,at%N), f1(3,at%N), f(3,at%N), f_hyb(3,at%N), fdiff(3,at%N) )


  ! initialise dynamics
  call ds_initialise(ds, at)
  call rescale_velo(ds, 2000.0_dp)
  call  zero_momentum(ds)
  dt = 1.0_dp


  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! X
  ! X  bootstrap the adjustable potential
  ! X
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  ! create list of embedded atoms
  call append(embedlist, (/1,0,0,0/))
  call bfs_grow(ds%atoms, embedlist, 2)
!  call discard_non_min_images(embedlist)

  call add_property(ds%atoms, 'QM', 0)
  if (.not. assign_pointer(ds%atoms, 'QM', qm)) &
       call system_abort('QM property missing')
  qm(int_part(embedlist,1)) = 1


  
  ! compute default forces
  call hybrid_force(ds%atoms, f0, SW_force_noopt) 
  call hybrid_force(ds%atoms, f1, embedlist, qm_force=SWeps_force, &
          classical_force=SW_force_noopt, buf_grow_hops=3, &
          little_clusters=.false.,terminate=.true., &
          periodic_clusters=(/.false.,.false., .false./), randomise_buffer=.true.)


  ! grow the embed list to include a fit zone
  fitlist = embedlist
  call bfs_grow(ds%atoms, fitlist, 3, min_images_only=.true.)
  !call discard_non_min_images(fitlist)

  call reallocate(df, 3, fitlist%N, zero=.true.)
  ! copy force differences for the embedding zone... 
  df(:,1:embedlist%N) = f1(:,int_part(embedlist,1))-f0(:,int_part(embedlist,1))

     
  call wipe(fitforce)
  call append(fitforce, int_part(fitlist), df) 

  call adjustable_potential_init(ds%atoms, fitlist, fitlist%N, method='SVD', &
       spring_hops=2) !threebody_springs=.true.)

  call adjustable_potential_optimise(ds%atoms, real_part(fitforce), method='SVD')
!  call adjustable_potential_force(ds%atoms, df, power = dV_dt)

!  f = f0
!  f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1)) + df

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
  
  do n=1, 100000

     call print('')
     call print('====================   Quantum Selection     ====================')
     call print('')

     ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     ! X
     ! X Code should go here to set up embedlist
     ! X
     ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     write(line,'(i0,a)') embedlist%N,' atoms flagged for quantum treatment'
     call print(line)

     !call print('Building fit zone...')
     !fitlist = embedlist
     !call bfs_grow(ds%atoms, fitlist, 3, min_images_only=.true.)
     !call discard_non_min_images(fitlist)

     
     write(line,'(a,i0,a)') 'Fitting on ',fitlist%N,' atoms in total'
     call print(line)
     
     call adjustable_potential_init(ds%atoms, fitlist, fitlist%N, map=.true., &
          method='SVD', spring_hops=2)! threebody_springs=.true.)


     call print('')
     call print('====================     Extrapolation     ====================')
     call print('')

     !only reallocate df if needed
     call reallocate(df,3,fitlist%N)

     ds_saved = ds

     do i = 1, 10

        Epot = SW_energy(ds%atoms)
        call hybrid_force(ds%atoms, f, SW_force_noopt)
        

        ! get force from optimised adjustable potential
           
        call adjustable_potential_force(ds%atoms, df, power=dV_dt, energy=energy)
        Epot = Epot + energy

        f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1)) + df


        ! advance the dynamics
        call advance_verlet(ds, dt, f, dV_dt)


        call ds_print_status(ds, 'E', Epot)

        call print_xyz(ds%atoms,  movie, properties='pos:QM')

     end do


     do j = 1, 1

           call print('')
           call print('====================   Computation of forces   ====================')
           call print('')

           call print('Computing new forces')
           !compute default force

           call hybrid_force(ds%atoms, f0, SW_force_noopt)


           ! compute hybrid forces
           call hybrid_force(ds%atoms, f1, embedlist, qm_force=SWeps_force, &
                classical_force=SW_force_noopt, buf_grow_hops=10, &
                little_clusters=.false.,terminate=.false., &
                periodic_clusters=(/.true.,.true., .true./), randomise_buffer=.true.)

           call print('f0:')
           call print(f0)
           call print('f1:')
           call print(f1)
           call print('')
           call print('==================== Optimisation of parameters ====================')
           call print('')
     

           df = 0.0_dp
           ! copy new force differences for the same embedding zone as before...
           df(:,1:embedlist%N) = f1(:,int_part(embedlist,1))-f0(:,int_part(embedlist,1))
           call wipe(fitforce)
           call append(fitforce, int_part(fitlist), df)

           ! now optimise the parameters for these forces at the new positions
           !     call adjustable_potential_optimise(ds%atoms, real_part(fitforce), &
           !          max_force=max(maxval(abs(f0)), maxval(abs(f1))))


           call verbosity_push(NERD)
           call adjustable_potential_optimise(ds%atoms, real_part(fitforce), method='SVD')
           call verbosity_pop()

           call print('')
           call print('====================        Interpolation        ====================')
           call print('')

           ! revert to the saved dynamical system
           ds = ds_saved

           do i = 1, 10

              call hybrid_force(ds%atoms, f, SW_force_noopt) 
              Epot = SW_energy(ds%atoms)
              
              ! get force differences from interpolated parameters
              call adjustable_potential_force(ds%atoms, df, power=dV_dt, &
                   interp=real(i-1,dp)/real(n_interp_steps,dp), energy=energy)
              Epot = Epot + energy

              ! add it onto default forces
        

              ! advance the dynamics
              call advance_verlet(ds, dt, f, dV_dt)
              call ds_print_status(ds, 'I', Epot)

           end do
        end do




     call print('')
     call print('==================== Recalculate Connectivity ====================')   
     call print('')

     ! print movie and recompute connectivity
     call calc_connect(ds%atoms)
     call print_xyz(ds%atoms, movie, properties='pos:QM')

  end do
  

  call adjustable_potential_finalise
  call atoms_finalise(at)
  call atoms_finalise(dia)
  call ds_finalise(ds)
  call ds_finalise(ds_saved)
  call finalise(embedlist)
  call finalise(fitlist)
  call finalise(fitforce)
  deallocate(f0,f1,f,f_hyb,df)

  call system_finalise()

end program lotfpc

