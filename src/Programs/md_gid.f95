program md_gid

  use libatoms_module
  use potential_module
  implicit none

  type(Dictionary)      :: params

  type(MPI_context) :: mpi_glob
  type(Potential)       :: pot
  type(Atoms) :: at1, at2
  type(DynamicalSystem) :: ds1, ds2

  real(dp) :: p0, p_init, p, temp_init, temp, dt, deltaTemp, deltaBeta, &
     lnP0, lnP, beta0, beta, deltaV, deltaH, f0, f1, tau, tau_cell, cell_oscillation_time, &
     bulk_modulus_estimate1, bulk_modulus_estimate2
  integer :: nTotal, nSteps, nCorrector, nSmooth
  logical :: iso1, iso2, restart
  character(STRING_LENGTH) :: at_file1, at_file2, param_file, init_args
  integer :: i, j
  real(dp), dimension(:,:), allocatable :: vol, kine, pote, enth, tp

  call system_initialise(verbosity=PRINT_NORMAL, enable_timing=.false.)

  call initialise(params)
  call param_register(params, 'at_file1', PARAM_MANDATORY,at_file1,"Phase 1 input configuration")
  call param_register(params, 'at_file2', PARAM_MANDATORY,at_file2,"Phase 2 input configuration")
  call param_register(params, 'param_file', PARAM_MANDATORY,param_file,"XML file")
  call param_register(params, 'init_args', PARAM_MANDATORY,init_args,"Init args")
  call param_register(params, 'p_init', PARAM_MANDATORY,p_init,"Initial pressure (GPa) at a coexistence point")
  call param_register(params, 'temp_init', PARAM_MANDATORY,temp_init,"Initial temperature at a coexistence point")
  call param_register(params, 'iso1', "false",iso1,"Allow only isotropic lattice deformation in phase 1")
  call param_register(params, 'iso2', "false",iso1,"Allow only isotropic lattice deformation in phase 2")
  call param_register(params, 'tau', "100.0",tau,"Thermostat time constant")
  call param_register(params, 'tau_cell', "2000.0",tau_cell,"Barostat time constant")
  call param_register(params, 'cell_oscillation_time', "1000.0",cell_oscillation_time,"Cell vibration time period")
  call param_register(params, 'bulk_modulus_estimate1', "150.0",bulk_modulus_estimate1,"Estimated bulk modulus (GPa) for phase 1")
  call param_register(params, 'bulk_modulus_estimate2', "150.0",bulk_modulus_estimate2,"Estimated bulk modulus (GPa) for phase 2")
  call param_register(params, 'restart', "false",restart,"Do not rescale velocities")
  call param_register(params, 'nsteps', "30",nsteps,"GD integration time steps")
  call param_register(params, 'ncorrector', "6",ncorrector,"GD integration corrector steps")
  call param_register(params, 'ntotal', "15000",ntotal,"MD steps per GD integration iteration")
  call param_register(params, 'nsmooth', "5000",nSmooth,"Time constant (in steps) for averaging")
  call param_register(params, 'dt', "1.0",dt,"MD time step")
  call param_register(params, 'deltatemp', "25",deltatemp,"temperature step in GD integration")

  if (.not. param_read_args(params, check_mandatory = .true.)) then
     call print("Usage: md_gid ")
     call system_abort('Exit: Mandatory argument(s) missing...')
  endif

  call print("==================== Input parameters ====================")
  call param_print(params)
  call print("==========================================================")

  call finalise(params)

  call initialise(mpi_glob)
  call Potential_Filename_Initialise(pot, args_str=trim(init_args), param_filename=trim(param_file), mpi_obj=mpi_glob)

  p = p_init/GPA
  temp = temp_init

  ! System 1
  call read(at1,trim(at_file1),mpi=mpi_glob)
  call set_cutoff(at1,cutoff(pot),1.0_dp)
  call initialise(ds1, at1)

  if(iso1) then
     call add_thermostat(ds1,THERMOSTAT_LANGEVIN_NPT,t=temp,p=p,tau=tau,tau_cell=tau_cell, &
        bulk_modulus_estimate=bulk_modulus_estimate1/GPA,cell_oscillation_time=cell_oscillation_time)
  else
     call add_thermostat(ds1,THERMOSTAT_LANGEVIN_PR,t=temp,p=p,tau=tau,tau_cell=tau_cell, &
        bulk_modulus_estimate=bulk_modulus_estimate1/GPA,cell_oscillation_time=cell_oscillation_time)
  endif
  if( .not. restart) call rescale_velo(ds1,2*temp_init)

  ! System 1
  call read(at2,trim(at_file2),mpi=mpi_glob)
  call set_cutoff(at2,cutoff(pot),1.0_dp)
  call initialise(ds2, at2)

  if(iso2) then
     call add_thermostat(ds2,THERMOSTAT_LANGEVIN_NPT,t=temp,p=p,tau=tau,tau_cell=tau_cell, &
        bulk_modulus_estimate=bulk_modulus_estimate2/GPA,cell_oscillation_time=cell_oscillation_time)
  else
     call add_thermostat(ds2,THERMOSTAT_LANGEVIN_PR,t=temp,p=p,tau=tau,tau_cell=tau_cell, &
        bulk_modulus_estimate=bulk_modulus_estimate2/GPA,cell_oscillation_time=cell_oscillation_time)
  endif
  if( .not. restart) call rescale_velo(ds2,2*temp_init)

  allocate(vol(2,nTotal),kine(2,nTotal),pote(2,nTotal))
  allocate(tp(2,nSteps))

  do j = 1, nSteps
     lnP0 = log(p)
     p0 = p
     beta0=1.0_dp/(BOLTZMANN_K*temp)
     beta=beta0

     vol=0.0_dp
     kine=0.0_dp
     pote=0.0_dp

     call runMD(ds1, ds2, pot, temp,p,vol,kine,pote,nTotal)
     enth = kine + pote + p*vol

     deltaV = expav(vol(1,:),dt*nSmooth) - expav(vol(2,:),dt*nSmooth)
     deltaH = expav(enth(1,:),dt*nSmooth) - expav(enth(2,:),dt*nSmooth)
     
     !f0=-deltaH/(beta*p*deltaV)
     f0=-deltaH/(beta*deltaV)
     print*, "delta in predictor", deltaV, deltaH, f0

     temp = temp + deltaTemp
     beta=1.0_dp/(BOLTZMANN_K*temp)
     deltaBeta=beta-beta0

     p = p0 + f0*deltaBeta
     lnP = log(p)
     !lnP = lnP0+f0*deltaBeta
     !p = exp(lnP)
     print*, "temperature", temp
     print*, "pressure in predictor step", p*GPA, lnP

     do i = 1, nCorrector
        call runMD(ds1, ds2, pot,temp,p,vol,kine,pote,nTotal)
        enth = kine + pote + p*vol

        deltaV = expav(vol(1,:),dt*nSmooth) - expav(vol(2,:),dt*nSmooth)
        deltaH = expav(enth(1,:),dt*nSmooth) - expav(enth(2,:),dt*nSmooth)

        !f1=-deltaH/(beta*p*deltaV)
        f1=-deltaH/(beta*deltaV)
        print*, "delta in corrector", i, deltaV, deltaH, f1
        p = p0 + (f0+f1)*deltaBeta/2.0_dp
        lnP = log(p)
        !lnP = lnP0 + (f0+f1)*deltaBeta/2.0_dp
        !p = exp(lnP)
        print*, "pressure in corrector step", p*GPA, lnP
        tp(:,j) = (/temp,p*GPA/)
     enddo
  enddo

  call print("final coexistence (T,p)")
  call print(tp)

  call finalise(ds1)
  call finalise(ds2)
  call finalise(at1, at2)
  deallocate(vol,kine,pote,tp)
  call system_finalise()

  contains

     subroutine runMD(ds1, ds2, pot, t,p,vol,kine,pote,nTotal)
        type(DynamicalSystem), intent(inout) :: ds1, ds2
        type(Potential), intent(inout) :: pot

        real(dp), intent(in) :: t, p
        real(dp), dimension(:,:), intent(out) :: vol, kine, pote
        integer, intent(in) :: nTotal

        integer :: i
        real(dp), dimension(3,3) :: virial1, virial2
        real(dp) :: energy1, energy2
        real(dp), dimension(:,:), pointer :: force1, force2
        
        call add_property(ds1%atoms,'force',n_cols=3,value=0.0_dp,ptr2=force1)
        call add_property(ds2%atoms,'force',n_cols=3,value=0.0_dp,ptr2=force2)

        call update_thermostat(ds1,t=t,p=p)
        call update_thermostat(ds2,t=t,p=p)

        call calc(pot,ds1%atoms, virial=virial1)
        call calc(pot,ds2%atoms, virial=virial2)

        do i = 1, nTotal
           call advance_verlet1(ds1,dt,virial=virial1)
           call calc(pot,ds1%atoms, force=force1, energy=energy1, virial=virial1)
           call advance_verlet2(ds1, dt, force1,virial=virial1)

           call advance_verlet1(ds2,dt,virial=virial2)
           call calc(pot,ds2%atoms, force=force2, energy=energy2, virial=virial2)
           call advance_verlet2(ds2, dt, force2,virial=virial2)

           vol(:,i) = (/cell_volume(ds1%atoms),cell_volume(ds2%atoms)/)
           kine(:,i) = (/[kinetic_energy(ds1),kinetic_energy(ds2)]/)
           pote(:,i) = (/energy1,energy2/)
        enddo
        call map_into_cell(ds1%atoms)
        call map_into_cell(ds2%atoms)
     endsubroutine runMD

     function expav(x,time)
        real(dp), dimension(:), intent(in) :: x
        real(dp), intent(in) :: time
        real(dp) :: expav

        real(dp) :: f1, f2
        real(dp) :: ave
        integer :: i

        f1=exp(-1.0_dp/time)
        f2=1.0_dp-f1
        ave=x(1)
     
        do i = 2, size(x)
           ave = f1*ave + f2*x(i)
        enddo

        expav = ave
     endfunction expav

endprogram md_gid
