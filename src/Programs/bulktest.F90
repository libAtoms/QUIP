! HJ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HJ X
! HJ X   libAtoms+QUIP: atomistic simulation library
! HJ X
! HJ X   Portions of this code were written by
! HJ X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! HJ X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! HJ X
! HJ X   Copyright 2006-2010.
! HJ X
! HJ X   These portions of the source code are released under the GNU General
! HJ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
! HJ X
! HJ X   If you would like to license the source code under different terms,
! HJ X   please contact Gabor Csanyi, gabor@csanyi.net
! HJ X
! HJ X   Portions of this code were written by Noam Bernstein as part of
! HJ X   his employment for the U.S. Government, and are not subject
! HJ X   to copyright in the USA.
! HJ X
! HJ X
! HJ X   When using this software, please cite the following reference:
! HJ X
! HJ X   http://www.libatoms.org
! HJ X
! HJ X  Additional contributions by
! HJ X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! HJ X
! HJ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

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

program bulktest

  use libAtoms_module
  use LOTF_module

  implicit none

  ! Parameters for the program (we should probably put a nice command line argument
  ! reader in System...)
 
  logical,  parameter :: do_lots_of_little_clusters = .true.
  integer,  parameter :: n_pc_iter = 1        ! number of predictor/corrector iterations in each cycle
  integer,  parameter :: n_interp_steps = 10   ! number of steps to extrapolate/interpolate for
  integer,  parameter :: n_lotf_cycles = 100000 ! number of extrap/interp cycles to do
  real(dp), parameter :: dt = 1.0_dp  ! time step
  real(dp), parameter :: init_temp = 300.0_dp
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
  call system_initialise(PRINT_NORMAL,1)
  call initialise(movie, "movie.xyz")
  !call initialise(datafile, 'forcediff.dat')

  ! create some atoms
  call diamond(dia, 5.44_dp)
  call supercell(at, dia, 5,5,5)
  at%Z = 14
  call set_cutoff(at, 4.0_dp)
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
  ! with no embedded atoms
  call hybrid_force(ds%atoms, f0, classical_force=SW_Force_noopt) 
  
  ! compute hybrid forces
  call hybrid_force(ds%atoms, f1, embedlist, qm_force=SWeps_force, &
       classical_force=SW_force_noopt, buf_grow_hops=3, &
       little_clusters=.true.,terminate=.true., &
       periodic_clusters=(/.false.,.false.,.false./), randomise_buffer=.true.)

  
  ! grow the embed list to include a fit zone
  fitlist = embedlist
  call BFS_grow(ds%atoms, fitlist, 3)

  ! copy force differences for the embedding zone... 
  call reallocate(df, 3, fitlist%N, zero=.true.)
  df(:,1:embedlist%N) = f1(:,int_part(embedlist,1))-f0(:,int_part(embedlist,1))
  call wipe(fitforce)
  call append(fitforce, int_part(fitlist), df) 

  call adjustable_potential_init(ds%atoms, fitlist, embedlist%N)
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
     ! X Code could go here to set up new embedlist
     ! X
     ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     write(line,'(i0,a)') embedlist%N,' atoms flagged for quantum treatment'
     call Print(line)

     call Print('Building fit zone...')
     fitlist = embedlist
     call BFS_grow(ds%atoms, fitlist, 3)
     
     write(line,'(a,i0,a)') 'Fitting on ',fitlist%N,' atoms in total'
     call Print(line)
     
     call adjustable_potential_init(ds%atoms, fitlist, embedlist%N, map=.true.)

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
        call hybrid_force(ds%atoms, f, classical_force=SW_Force_noopt)
        f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1)) + df

        ! ********** TESTING ***********
        ! Get the 'actual' force
        !call hybrid_force(ds%atoms, f_hyb, embedlist, SW_Force_noopt, SWeps_force,do_lots_of_little_clusters)
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
        !compute default force
        call hybrid_force(ds%atoms, f0, classical_force=SW_Force_noopt)
        ! compute hybrid forces
        call hybrid_force(ds%atoms, f1, embedlist, qm_force=SWeps_force, &
             classical_force=SW_force_noopt, buf_grow_hops=3, &
             little_clusters=.true.,terminate=.true., &
             periodic_clusters=(/.false., .false.,.false./),randomise_buffer=.true.)
        
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
           call hybrid_force(ds%atoms, f, classical_force=SW_Force_noopt) 
           f(:,int_part(fitlist,1)) = f(:,int_part(fitlist,1))+df
           
           !*********** TESTING *************
           !get the actual force
           !call hybrid_force(ds%atoms, f_hyb, embedlist, SW_Force_noopt, SWeps_force,do_lots_of_little_clusters)
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
     call print_xyz(ds%atoms, movie, properties='pos:QM:avg_ke')

  end do
  

  call adjustable_potential_finalise
  call Atoms_Finalise(at)
  call Atoms_Finalise(dia)
  call DS_Finalise(ds)
  call DS_Finalise(ds_saved)
  call finalise(movie)
  call finalise(embedlist)
  call finalise(fitlist)
  call finalise(fitforce)
  deallocate(f0,f1,f,f_hyb,df)

  call system_finalise()

end program bulktest
