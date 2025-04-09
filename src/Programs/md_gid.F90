program md_gid

  use libatoms_module
  use potential_module
  implicit none

  type(Dictionary)      :: params

  type(MPI_context) :: mpi_glob
  type(Potential)       :: pot
  type(Atoms) :: at1, at2
  type(DynamicalSystem) :: ds1, ds2
  type(CInoutput) :: savexyz

  real(dp) :: p0, p_init, p, temp_init, temp, dt, deltaTemp, deltaBeta, &
     lnP0, lnP, beta0, beta, deltaV, deltaH, f0, f1, tau, tau_cell, cell_oscillation_time, &
     bulk_modulus_estimate1, bulk_modulus_estimate2
  integer :: nTotal, nSteps, nCorrector, nSmooth, nPrint
  logical :: iso1, iso2, restart, has_trajectory_file
  character(STRING_LENGTH) :: at_file1, at_file2, param_file, init_args
  integer :: i, j
  real(dp), dimension(:,:), allocatable :: vol, kine, pote, enth, tp
  character(STRING_LENGTH) :: trajectory_file

  call system_initialise(verbosity=PRINT_NORMAL, enable_timing=.false.)

  call initialise(params)
  call param_register(params, 'at_file1', PARAM_MANDATORY,at_file1,"Phase 1 input configuration")
  call param_register(params, 'at_file2', PARAM_MANDATORY,at_file2,"Phase 2 input configuration")
  call param_register(params, 'param_file', PARAM_MANDATORY,param_file,"XML file")
  call param_register(params, 'init_args', PARAM_MANDATORY,init_args,"Init args")
  call param_register(params, 'p_init', PARAM_MANDATORY,p_init,"Initial pressure (GPa) at a coexistence point")
  call param_register(params, 'temp_init', PARAM_MANDATORY,temp_init,"Initial temperature at a coexistence point")
  call param_register(params, 'iso1', "false",iso1,"Allow only isotropic lattice deformation in phase 1")
  call param_register(params, 'iso2', "false",iso2,"Allow only isotropic lattice deformation in phase 2")
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
  call param_register(params, 'nprint', "100",nPrint,"Status print")
  call param_register(params, 'dt', "1.0",dt,"MD time step")
  call param_register(params, 'deltatemp', "25",deltatemp,"temperature step in GD integration")
  call param_register(params, 'trajectory_file', "",trajectory_file,"file to save trajectory along integration",has_value_target=has_trajectory_file)

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

  if(has_trajectory_file) call initialise(savexyz,trim(trajectory_file),action=OUTPUT,append=.true.,mpi=mpi_glob)

  p = p_init/EV_A3_IN_GPA
  temp = temp_init

  ! System 1
  call read(at1,trim(at_file1),mpi=mpi_glob)
  call set_cutoff(at1,cutoff(pot),1.0_dp)
  call initialise(ds1, at1)

  if(iso1) then
     call add_thermostat(ds1,THERMOSTAT_LANGEVIN_NPT,t=temp,p=p,tau=tau,tau_cell=tau_cell, &
        bulk_modulus_estimate=bulk_modulus_estimate1/EV_A3_IN_GPA,cell_oscillation_time=cell_oscillation_time)
  else
     call add_thermostat(ds1,THERMOSTAT_LANGEVIN_PR,t=temp,p=p,tau=tau,tau_cell=tau_cell, &
        bulk_modulus_estimate=bulk_modulus_estimate1/EV_A3_IN_GPA,cell_oscillation_time=cell_oscillation_time)
  endif
  if( .not. restart) call rescale_velo(ds1,2*temp_init)

  ! System 1
  call read(at2,trim(at_file2),mpi=mpi_glob)
  call set_cutoff(at2,cutoff(pot),1.0_dp)
  call initialise(ds2, at2)

  if(iso2) then
     call add_thermostat(ds2,THERMOSTAT_LANGEVIN_NPT,t=temp,p=p,tau=tau,tau_cell=tau_cell, &
        bulk_modulus_estimate=bulk_modulus_estimate2/EV_A3_IN_GPA,cell_oscillation_time=cell_oscillation_time)
  else
     call add_thermostat(ds2,THERMOSTAT_LANGEVIN_PR,t=temp,p=p,tau=tau,tau_cell=tau_cell, &
        bulk_modulus_estimate=bulk_modulus_estimate2/EV_A3_IN_GPA,cell_oscillation_time=cell_oscillation_time)
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

     !deltaV = expav(vol(1,:),dt*nSmooth) - expav(vol(2,:),dt*nSmooth)
     !deltaH = expav(enth(1,:),dt*nSmooth) - expav(enth(2,:),dt*nSmooth)
     deltaV = runav(vol(1,:),nSmooth) / ds1%atoms%n - runav(vol(2,:),nSmooth) / ds2%atoms%n
     deltaH = runav(enth(1,:),nSmooth) / ds1%atoms%n - runav(enth(2,:),nSmooth) / ds2%atoms%n
     
     !f0=-deltaH/(beta*p*deltaV)
     f0=-deltaH/(beta*deltaV)
     call print("delta in predictor "//deltaV//"  "//deltaH//"  "//f0)

     temp = temp + deltaTemp
     beta=1.0_dp/(BOLTZMANN_K*temp)
     deltaBeta=beta-beta0

     p = p0 + f0*deltaBeta
     lnP = log(p)
     !lnP = lnP0+f0*deltaBeta
     !p = exp(lnP)
     call print("temperature  "//temp)
     call print("pressure in predictor step   "//p*EV_A3_IN_GPA//"   "//lnP)

     do i = 1, nCorrector
        call runMD(ds1, ds2, pot,temp,p,vol,kine,pote,nTotal)
        enth = kine + pote + p*vol

        !deltaV = expav(vol(1,:),dt*nSmooth) - expav(vol(2,:),dt*nSmooth)
        !deltaH = expav(enth(1,:),dt*nSmooth) - expav(enth(2,:),dt*nSmooth)
        deltaV = runav(vol(1,:),nSmooth) / ds1%atoms%n - runav(vol(2,:),nSmooth) / ds2%atoms%n
        deltaH = runav(enth(1,:),nSmooth) / ds1%atoms%n - runav(enth(2,:),nSmooth) / ds2%atoms%n

        !f1=-deltaH/(beta*p*deltaV)
        f1=-deltaH/(beta*deltaV)
        call print("delta in predictor "//i//"   "//deltaV//"  "//deltaH//"  "//f0)
        p = p0 + (f0+f1)*deltaBeta/2.0_dp
        lnP = log(p)
        !lnP = lnP0 + (f0+f1)*deltaBeta/2.0_dp
        !p = exp(lnP)
        call print("pressure in corrector step   "//p*EV_A3_IN_GPA//"   "//lnP)
        tp(:,j) = (/temp,p*EV_A3_IN_GPA/)
     enddo
     if(has_trajectory_file) then
        call set_value(ds1%atoms%params, 'temperature', "" // temp)
        call set_value(ds2%atoms%params, 'temperature', "" // temp)
        call set_value(ds1%atoms%params, 'pressure', "" // p)
        call set_value(ds2%atoms%params, 'pressure', "" // p)
        call write(ds1%atoms,savexyz,properties="species:pos:velo:travel:force")
        call write(ds2%atoms,savexyz,properties="species:pos:velo:travel:force")
     endif
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
           if(mod(i,nPrint)==0) then
              call ds_print_status(ds1,label="SYS1",epot=energy1,instantaneous=.true.,mpi_obj=mpi_glob)
              call ds_print_status(ds2,label="SYS2",epot=energy2,instantaneous=.true.,mpi_obj=mpi_glob)
           endif
        enddo
        call map_into_cell(ds1%atoms)
        call map_into_cell(ds2%atoms)

        call set_value(ds1%atoms%params, 'Energy', "" // energy1)
        call set_value(ds2%atoms%params, 'Energy', "" // energy2)
        call set_value(ds1%atoms%params, 'Virial', "" // reshape(virial1,(/9/)) )
        call set_value(ds2%atoms%params, 'Virial', "" // reshape(virial2,(/9/)) )
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

     function runav(x,time)
        real(dp), dimension(:), intent(in) :: x
        integer, intent(in) :: time
        real(dp) :: runav

        integer :: n

        n = min(size(x,1), time)
        runav = sum(x(size(x,1)-n+1:)) / n

     endfunction runav

endprogram md_gid
