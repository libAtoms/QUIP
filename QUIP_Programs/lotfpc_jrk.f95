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

! $Id: lotfpc_jrk.f95,v 1.3 2008-04-09 11:26:56 jrk33 Exp $

! $Log: not supported by cvs2svn $
! Revision 1.2  2008/04/09 11:00:46  jrk33
! lotf.h has gone!
!
! Revision 1.1.1.1  2008/03/13 17:36:35  gc121
! this module contains the LOTF top level programs
!
! Revision 1.2  2007/08/31 14:15:12  jrk33
! Add thermostat and inject_heat params
!
! Revision 1.1  2007/08/21 09:16:56  jrk33
! Experimental version of lotfpc
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


!!#define FIT_FORCE

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
  type(Dictionary) :: params

  ! initialise program
  call system_initialise()


  call initialise(params)
  call param_register(params, 'single_cluster', '1', single_cluster)
  call param_register(params, 'full_qm', '0', full_qm)
  call param_register(params, 'grow_hops', '3', grow_hops)
  call param_register(params, 'classical', '0', classical)
  call param_register(params, 'n_pc_iter', '1', n_pc_iter)
  call param_register(params, 'n_interp_steps', '10', n_interp_steps)
  call param_register(params, 'n_lotf_cycles', '100000', n_lotf_cycles)
  call param_register(params, 'dt', '1.0', dt)
  call param_register(params, 'init_temp', '2000.0', init_temp)
  call param_register(params, 'sim_temp', '1000.0', sim_temp)
  call param_register(params, 'seed', '0', seed)
  call param_register(params, 'nneightol', '1.3', nneightol)
  call param_register(params, 'spring_hops', '2', spring_hops)
  call param_register(params, 'fit_grow_hops', '3', fit_grow_hops)
  call param_register(params, 'refit', '0', refit)
  call param_register(params, 'thermostat', '0', thermostat)
  call param_register(params, 'inject_heat', '0.0', inject_heat)

  
  
  if (.not. param_read_args(params)) call system_abort('Error parsing command line')

  call print('Parameters:')
  call print('')
  call param_print(params)
  call print('')

  do_periodic = .false.
!!$  if (.not. single_cluster) then
!!$     do_periodic = .false.
!!$  else
!!$     do_periodic = .true.
!!$  end if

  ! Reseed random number generator if necessary
  if (seed /= 0) call hello_world(seed, common_seed = .true.)

  !call initialise(movie, "movie.xyz")
  !call initialise(datafile, 'posdiff.dat')

  ! create some atoms
  call diamond(dia, 5.44_dp)
  call supercell(at, dia, 3,3,3)
  at%Z = 14

!!$  call graphene_sheet(at, 1.42_dp, 6, 0, 2, 5)
!!$  at%Z = 6

  !call randomise(at%pos, .1_dp)

  call atoms_set_cutoff(at, 5.0_dp)
  call calc_connect(at)

  ! allocate some force arrays
  allocate( f0(3,at%N), f1(3,at%N), f(3,at%N), f_hyb(3,at%N), fdiff(3,at%N) )

!!$  allocate (qm_pos(3,at%N,n_interp_steps+1))

  ! initialise dynamics
  call ds_initialise(ds, at)
  call rescale_velo(ds, init_temp)
  !call zero_momentum(ds)

  ds%thermostat = thermostat
  ds%sim_temp = sim_temp
!!$  ds%Q_nose = suggested_Q_nose(ds)


  ds%atoms%nneightol = nneightol

  ! load fully QM trajectory
!!$  if (full_qm == 0) then
!!$     call initialise(qmtrajectory,'fullqm.xyz')
!!$     do i=1,n_interp_steps+1
!!$        call read_xyz(at,qmtrajectory)
!!$        qm_pos(:,:,i) = at%pos
!!$     end do
!!$     call finalise(qmtrajectory)
!!$  else
!!$     call print_xyz(ds%atoms,movie)
!!$  end if

  !allocate space for the test data array
  allocate(data(9,n_interp_steps+1))

  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! X
  ! X  bootstrap the adjustable potential
  ! X
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  ! create list of embedded atoms
  call append(embedlist, (/1,0,0,0/))
  call bfs_grow(ds%atoms, embedlist, 2)
  call discard_non_min_images(embedlist)

  ! treat all atoms with QM
!!$  do i=1,64
!!$     call append(embedlist, (/i,0,0,0/))
!!$  end do

  call add_property(ds%atoms, 'QM', 0)
  if (.not. assign_pointer(ds%atoms, 'QM', qm)) &
       call system_abort('QM property missing')
  qm(int_part(embedlist,1)) = 1

  if (.not. (full_qm == 1 .or. classical == 1)) then

  
     ! compute default forces
     call hybrid_force(ds%atoms, f0, SW_force_noopt) 

     ! compute hybrid forces
!!$  call hybrid_force(ds%atoms, f1, embedlist, SW_force_noopt2, SW_force_noopt, &
!!$       little_clusters=do_lots_of_little_clusters, terminate=.false.)  

     call hybrid_force(ds%atoms, f1, embedlist, qm_force=SWeps_force, &
          classical_force=SW_force_noopt, buf_grow_hops=grow_hops, &
          little_clusters=(single_cluster == 0),terminate=.true., &
          periodic_clusters=do_periodic, randomise_buffer=.true.)


!!$     f1(3,:) = 0.0_dp

     ! grow the embed list to include a fit zone
     fitlist = embedlist
     call bfs_grow(ds%atoms, fitlist, fit_grow_hops)
     call discard_non_min_images(fitlist)

     if(thermostat) then
        ! Only thermostat the not-fit atoms so we can measure heat generated
        ! by the QM region
        do i=1,ds%atoms%N
           if (is_in_array(fitlist%int(1,1:fitlist%N),i)) then
              ds%atoms%thermostat_mask(i) = 0
           else
              ds%atoms%thermostat_mask(i) = 1
           end if
        end do
     end if

     call print(embedlist)

     call print('Thermostat applied to '//count(ds%atoms%thermostat_mask == 1)//' atoms')
     
     call reallocate(df, 3, fitlist%N, zero=.true.)
#ifdef FIT_FORCE
     ! fit to QM force not force difference
     df(:,1:fitlist%N) = f1(:,int_part(fitlist,1))
#else
     ! copy force differences for the embedding zone... 
     df(:,1:embedlist%N) = f1(:,int_part(embedlist,1))-f0(:,int_part(embedlist,1))
#endif
     
!!$     call print('f1(embed)')
!!$     call print(f1(:,int_part(embedlist,1)))
!!$     call print('f0(embed)')
!!$     call print(f0(:,int_part(embedlist,1)))

     call wipe(fitforce)
     call append(fitforce, int_part(fitlist), df) 

!!$     call print('fitforce=')
!!$     call print(fitforce)
!!$     call print('')

#ifdef FIT_FORCE
     call adjustable_potential_init(ds%atoms, fitlist)
#else
     call adjustable_potential_init(ds%atoms, fitlist, fitlist%N, method='SVD', &
          spring_hops=spring_hops) !threebody_springs=.true.)
#endif

     call print(fitlist)

#ifdef FIT_FORCE
     call adjustable_potential_optimise(ds%atoms, real_part(fitforce))
#else
     call adjustable_potential_optimise(ds%atoms, real_part(fitforce), method='SVD')
#endif
     call adjustable_potential_force(ds%atoms, df, power = dV_dt)

!!$     call print('after optim, df=')
!!$     call print(df)
!!$     call print('')

!!$     call print('maxval(df)='//maxval(df))
!!$     call print('rms(df)='//sqrt(norm2(reshape(df,(/3*fitforce%N/))/(3*fitforce%N))))

     f = f0
     f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1)) + df

!!$     call print('Max diff ='//maxval(abs(f-f1)))
!!$     call print('RMS diff = '//RMS_diff(f1,f,fitlist))
!!$
!!$  call adjustable_potential_print_params()

     call wipe(fitlist)
     call wipe(fitforce)

  end if

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

     if (.not. (full_qm == 1 .or. classical == 1)) then


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

     call print('Building fit zone...')
     fitlist = embedlist
     call bfs_grow(ds%atoms, fitlist, fit_grow_hops)
     call discard_non_min_images(fitlist)

     
     write(line,'(a,i0,a)') 'Fitting on ',fitlist%N,' atoms in total'
     call print(line)
     
#ifdef FIT_FORCE
     call adjustable_potential_init(ds%atoms, fitlist, map=.true.)
#else
     call adjustable_potential_init(ds%atoms, fitlist, fitlist%N, map=.true., &
          method='SVD', spring_hops=spring_hops)! threebody_springs=.true.)
#endif

     end if

     call print('')
     call print('====================     Extrapolation     ====================')
     call print('')

     !only reallocate df if needed
     call reallocate(df,3,fitlist%N)

     ds_saved = ds

     do i = 1, n_interp_steps

!!$        call hybrid_force(ds%atoms, f, SW_force_noopt)

        Epot = SW_energy(ds%atoms)
        call hybrid_force(ds%atoms, f, SW_force_noopt)
        
!!$        call print('classical force=')
!!$        call print(f)

!!$        f(3,:) = 0.0_dp

        if (.not. (full_qm == 1 .or. classical == 1)) then

        ! get force from optimised adjustable potential
!!$        call print('Extrapolated parameters <-- param')

           
        call adjustable_potential_force(ds%atoms, df, power=dV_dt, energy=energy)
        Epot = Epot + energy

#ifdef FIT_FORCE
        ! replace with fit forces in fit zone
        f(:,int_part(fitlist,1)) = df
#else
        ! add it onto default forces
        f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1)) + df
#endif
        end if

        if (full_qm == 1) then
!!$           call hybrid_force(ds%atoms, f, SW_force_noopt2) ! move with the QM force
           call SWeps_force(ds%atoms, f) ! move with the QM force

!!$           call print('TB force=')
!!$           call print(f)

        else if (refit == 1) then

           ! ********** TESTING ***********
           ! Get the 'actual' force

           call hybrid_force(ds%atoms, f0, SW_force_noopt)

           call hybrid_force(ds%atoms, f_hyb, embedlist, qm_force=SWeps_force, &
                classical_force=SW_force_noopt, buf_grow_hops=grow_hops, &
                little_clusters=(single_cluster == 0),terminate=.true., &
                periodic_clusters=do_periodic, randomise_buffer=.true.)

           data(1,i) = ds%t
           data(6,i) = RMS_diff(f_hyb,f,fitlist) ! RMS force differences
           fdiff = abs(f_hyb-f)
           data(7,i) = maxval(fdiff)      ! Max force error


           df = 0.0_dp
#ifdef FIT_FORCE 
           ! here we fit to QM force not force difference
           df(:,1:fitlist%N) = f_hyb(:,int_part(fitlist,1))
#else
           ! copy new force differences for the same embedding zone as before...
           df(:,1:embedlist%N) = f_hyb(:,int_part(embedlist,1))-f0(:,int_part(embedlist,1))
#endif
           call wipe(fitforce)
           call append(fitforce, int_part(fitlist), df)

           ! save the optimised spring constants
           twobody_save = twobody
           !threebody_save = threebody

           ! start from where we were after last fit
           twobody = twobody_old
           !threebody = threebody_old


#ifdef FIT_FORCE
           call adjustable_potential_optimise(ds%atoms, real_part(fitforce))
#else
           call adjustable_potential_optimise(ds%atoms, real_part(fitforce), method='SVD')
#endif

!!$           call Print('Refit parameters (extrap) <-- param')
!!$           call adjustable_potential_force(ds%atoms, df, dV_dt)

           twobody = twobody_save
           !threebody = threebody_save

!!$           call print('hybrid force=')
!!$           call print(f_hyb)

!!$           f_hyb(3,:) = 0.0_dp

!!$           data(8,i) = maxval(fdiff(2,:))
!!$           data(9,i) = maxval(fdiff(3,:))

!!$           data(2,i) = RMS_diff(ds%atoms%pos, qm_pos(:,:,int(ds%t+1)), embedlist)
!!$           data(3,i) = maxval(abs(ds%atoms%pos(:,int_part(embedlist,1))-&
!!$                qm_pos(:,int_part(embedlist,1),int(ds%t+1))))

           
           ! ******************************

        end if



        ! advance the dynamics
        call advance_verlet(ds, dt, f, dV_dt)


        ! inject some heat to QM atoms?
        if (inject_heat /= 0.0_dp) then
           
           do j=1,embedlist%N
              ds%atoms%velo(:,embedlist%int(1,j)) = inject_heat*ds%atoms%velo(:,embedlist%int(1,j))
           end do

        end if

        call ds_print_status(ds, 'E', Epot)

        call print_xyz(ds%atoms,  movie, properties='pos:QM')

     end do


     if (.not. (full_qm == 1 .or. classical == 1)) then

        do j = 1, n_pc_iter

           call print('')
           call print('====================   Computation of forces   ====================')
           call print('')

           call print('Computing new forces')
           !compute default force
!!$     call hybrid_force(ds%atoms, f0, SW_force_noopt)
           call hybrid_force(ds%atoms, f0, SW_force_noopt)

!!$           f0(3,:) = 0.0_dp

           ! compute hybrid forces
!!$     call hybrid_force(ds%atoms, f1, embedlist, SW_force_noopt2, SW_force_noopt, &
!!$          little_clusters=do_lots_of_little_clusters,terminate=.false.)  
           call hybrid_force(ds%atoms, f1, embedlist, qm_force=SWeps_force, &
                classical_force=SW_force_noopt, buf_grow_hops=grow_hops, &
                little_clusters=(single_cluster == 0),terminate=.true., &
                periodic_clusters=do_periodic, randomise_buffer=.true.)

!!$           f1(3,:) = 0.0_dp


           call print('')
           call print('==================== Optimisation of parameters ====================')
           call print('')
     

           df = 0.0_dp
#ifdef FIT_FORCE 
           ! here we fit to QM force not force difference
           df(:,1:fitlist%N) = f1(:,int_part(fitlist,1))
#else
           ! copy new force differences for the same embedding zone as before...
           df(:,1:embedlist%N) = f1(:,int_part(embedlist,1))-f0(:,int_part(embedlist,1))
#endif
           call wipe(fitforce)
           call append(fitforce, int_part(fitlist), df)

           ! now optimise the parameters for these forces at the new positions
           !     call adjustable_potential_optimise(ds%atoms, real_part(fitforce), &
           !          max_force=max(maxval(abs(f0)), maxval(abs(f1))))


#ifdef FIT_FORCE
           call adjustable_potential_optimise(ds%atoms, real_part(fitforce))
#else
           call adjustable_potential_optimise(ds%atoms, real_part(fitforce), method='SVD')
#endif

           ! ********** TESTING ***********
           call Print('Optimised parameters <-- param')
           call adjustable_potential_force(ds%atoms, df, power=dV_dt)
           f = f0

#ifdef FIT_FORCE
           f(:,int_part(fitlist,1)) = df
#else
           f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1)) + df
#endif
        
           data(1,i) = ds%t              ! time
!!$           data(2,i) = RMS_diff(ds%atoms%pos, qm_pos(:,:,int(ds%t+1)), embedlist)
!!$           data(3,i) = maxval(abs(ds%atoms%pos(:,int_part(embedlist,1))-&
!!$             qm_pos(:,int_part(embedlist,1),int(ds%t+1))))
           data(6,i) = RMS_diff(f1,f,fitlist) ! RMS force differences

           fdiff = abs(f1-f)
           data(7,i) = maxval(fdiff)      ! Max force error

           data(8,i) = RMS_diff(real_part(fitforce), df)  ! RMS error in df
           data(9,i) = maxval(abs(real_part(fitforce)-df)) ! Max error in df
           ! ******************************


           call print('')
           call print('====================        Interpolation        ====================')
           call print('')

           ! revert to the saved dynamical system
           ds = ds_saved

           do i = 1, n_interp_steps

              call hybrid_force(ds%atoms, f, SW_force_noopt) 
              Epot = SW_energy(ds%atoms)
              
              ! get force differences from interpolated parameters
!!$              call Print('Interpolated parameters')
              call adjustable_potential_force(ds%atoms, df, power=dV_dt, &
                   interp=real(i-1,dp)/real(n_interp_steps,dp), energy=energy)
              Epot = Epot + energy

              ! add it onto default forces
        
!!$        call hybrid_force(ds%atoms, f, SW_force_noopt) 

!!$              f(3,:) = 0.0_dp


              if (refit == 1) then
#ifdef FIT_FORCE
                 ! replace with fit forces in fit zone
                 f(:,int_part(fitlist,1)) = df
#else
                 f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1))+df
#endif

                 !*********** TESTING *************
                 
                 call hybrid_force(ds%atoms, f0, SW_force_noopt)
              
                 !get the actual force
                 call hybrid_force(ds%atoms, f_hyb, embedlist, qm_force=SWeps_force, &
                      classical_force=SW_force_noopt, buf_grow_hops=grow_hops, &
                      little_clusters=(single_cluster == 0),terminate=.true., &
                      periodic_clusters=do_periodic, randomise_buffer=.true.)

                 twobody_save = twobody
                 !threebody_save = threebody

                 twobody = twobody_old
                 !threebody = threebody_old

                 df = 0.0_dp
#ifdef FIT_FORCE 
                 ! here we fit to QM force not force difference
                 df(:,1:fitlist%N) = f_hyb(:,int_part(fitlist,1))
#else
                 ! copy new force differences for the same embedding zone as before...
                 df(:,1:embedlist%N) = f_hyb(:,int_part(embedlist,1))-f0(:,int_part(embedlist,1))
#endif
                 call wipe(fitforce)
                 call append(fitforce, int_part(fitlist), df)

#ifdef FIT_FORCE
                 call adjustable_potential_optimise(ds%atoms, real_part(fitforce))
#else
                 call adjustable_potential_optimise(ds%atoms, real_part(fitforce), method='SVD')
#endif

!!$              call Print('Refit parameters (interp) <-- param')
!!$              call adjustable_potential_force(ds%atoms, df, dV_dt)
                 
                 twobody = twobody_save
                 !threebody = threebody_save

                 data(8,i) = RMS_diff(f_hyb,f,fitlist)  ! RMS force error

                 fdiff = abs(f_hyb-f)
                 data(9,i) = maxval(fdiff)      ! Max force error
!!$              data(12,i) = maxval(fdiff(2,:))
!!$              data(13,i) = maxval(fdiff(3,:))
                 
!!$                 data(4,i) = RMS_diff(ds%atoms%pos, qm_pos(:,:,int(ds%t+1)), embedlist)
!!$                 data(5,i) = maxval(abs(ds%atoms%pos(:,int_part(embedlist,1))- &
!!$                      qm_pos(:,int_part(embedlist,1),int(ds%t+1))))
                 
                 !*********************************

              else
                 f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1))+df
              end if

              ! advance the dynamics
              call advance_verlet(ds, dt, f, dV_dt)
              call ds_print_status(ds, 'I', Epot)

           end do
        end do

     end if ! .not. (do_full_qm .or. do_classical)

!!$     data(4,i) = RMS_diff(ds%atoms%pos, qm_pos(:,:,int(ds%t+1)), embedlist)
!!$     data(5,i) = maxval(abs(ds%atoms%pos(:,int_part(embedlist,1))- &
!!$                   qm_pos(:,int_part(embedlist,1),int(ds%t+1))))

     !Write the data to the data file
     do i = 1, n_interp_steps+1
        write(line,'(13(f0.18,tr2))') data(:,i)
        call print(line,file=datafile)
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
  !call finalise(movie)
  !call finalise(datafile)
  call finalise(embedlist)
  call finalise(fitlist)
  call finalise(fitforce)
  deallocate(f0,f1,f,f_hyb,df)

  call system_finalise()

end program lotfpc

