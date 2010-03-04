module md_module
use libatoms_module
use quip_module
implicit none
private

  public :: md_params
  type md_params
    character(len=FIELD_LENGTH) :: atoms_in_file, params_in_file, trajectory_out_file
    integer :: N_steps
    real(dp) :: dt,  thermalise_wait_time, damping_time
    real(dp) :: T_initial, T, T_final, freq_DT, langevin_tau, p_ext
    real(dp) :: velocity_rescaling_time 
    real(dp) :: cutoff_buffer 
    integer :: velocity_rescaling_freq
    logical :: calc_virial, const_T, const_P
    character(len=FIELD_LENGTH) :: metapot_init_args, metapot_calc_args
    integer :: summary_interval, params_print_interval, at_print_interval, pot_print_interval
    integer :: rng_seed
    logical :: damping, rescale_velocity 
    logical :: annealing, zero_momentum, zero_angular_momentum
    logical :: calc_dipole_moment, quiet_calc, do_timing
    integer :: advance_md_substeps
  end type md_params

public :: get_params, print_params, print_usage, do_prints, initialise_md_thermostat, update_thermostat

contains

subroutine get_params(params, mpi_glob)
  type(md_params), intent(inout) :: params
  type(MPI_context), intent(in) :: mpi_glob

  type(Dictionary) :: md_params_dict
  type(Extendable_Str) :: es
  logical :: md_params_exist

  call initialise(md_params_dict)
  params%atoms_in_file=''
  params%trajectory_out_file=''
  params%params_in_file=''
  params%metapot_init_args=''
  params%metapot_calc_args=''
  call param_register(md_params_dict, 'atoms_in_file', 'stdin', params%atoms_in_file)
  call param_register(md_params_dict, 'trajectory_out_file', 'traj.xyz', params%trajectory_out_file)
  call param_register(md_params_dict, 'params_in_file', 'quip_params.xml', params%params_in_file)
  call param_register(md_params_dict, 'rng_seed', '-1', params%rng_seed)
  call param_register(md_params_dict, 'N_steps', '1', params%N_steps)
  call param_register(md_params_dict, 'dt', '1.0', params%dt)
  call param_register(md_params_dict, 'T', '0.0', params%T_initial)
  call param_register(md_params_dict, 'T_final', '-100.0', params%T_final)
  call param_register(md_params_dict, 'freq_DT', '10.0', params%freq_DT)
  call param_register(md_params_dict, 'p_ext', '0.0', params%p_ext)
  call param_register(md_params_dict, 'thermalise_wait_time', '10.0', params%thermalise_wait_time)
  call param_register(md_params_dict, 'velocity_rescaling_freq', '2', params%velocity_rescaling_freq)
  call param_register(md_params_dict, 'velocity_rescaling_time', '0.0', params%velocity_rescaling_time)
  call param_register(md_params_dict, 'const_T', 'F', params%const_T)
  call param_register(md_params_dict, 'const_P', 'F', params%const_P)
  call param_register(md_params_dict, 'annealing', 'F', params%annealing)
  call param_register(md_params_dict, 'damping', 'F', params%damping)
  call param_register(md_params_dict, 'damping_time', '10.0', params%damping_time)
  call param_register(md_params_dict, 'rescale_velocity', 'F', params%rescale_velocity)
  call param_register(md_params_dict, 'langevin_tau', '100.0', params%langevin_tau)
  call param_register(md_params_dict, 'calc_virial', 'F', params%calc_virial)
  call param_register(md_params_dict, 'metapot_init_args', PARAM_MANDATORY, params%metapot_init_args)
  call param_register(md_params_dict, 'cutoff_buffer', '0.5', params%cutoff_buffer)
  call param_register(md_params_dict, 'summary_interval', '1', params%summary_interval)
  call param_register(md_params_dict, 'params_print_interval', '-1', params%params_print_interval)
  call param_register(md_params_dict, 'at_print_interval', '100', params%at_print_interval)
  call param_register(md_params_dict, 'pot_print_interval', '-1', params%pot_print_interval)
  call param_register(md_params_dict, 'zero_momentum', 'F', params%zero_momentum)
  call param_register(md_params_dict, 'zero_angular_momentum', 'F', params%zero_angular_momentum)
  call param_register(md_params_dict, 'metapot_calc_args', '', params%metapot_calc_args)
  call param_register(md_params_dict, 'dipole_moment', 'F', params%calc_dipole_moment)
  call param_register(md_params_dict, 'quiet_calc', 'T', params%quiet_calc)
  call param_register(md_params_dict, 'do_timing', 'F', params%do_timing)
  call param_register(md_params_dict, 'advance_md_substeps', '-1', params%advance_md_substeps)

  inquire(file='md_params', exist=md_params_exist)
  if (md_params_exist) then
    call initialise(es)
    call read(es, 'md_params', convert_to_string=.true., mpi_comm=mpi_glob%communicator)
    if (.not. param_read_line(md_params_dict, string(es))) then
      call print_usage()
      call system_abort("Error reading params from md_params file")
    endif
    call finalise(es)
  end if

  if (.not. param_read_args(md_params_dict)) then
    call print_usage()
    call system_abort("Error reading params from command line")
  endif
  call finalise(md_params_dict)

  if (len(trim(params%metapot_init_args)) == 0) then
     call print_usage()
     call system_abort("get_params got empty metapot_init_args")
  end if

  if (params%N_steps < 1) call system_abort("get_params got N_steps " // params%N_steps // " < 1")
  if (params%dt <= 0.0_dp) call system_abort("get_params got dt " // params%dt // " <= 0.0")
  if (params%const_T) then
    if (params%langevin_tau <= 0.0_dp) call system_abort("get_params got const_T, but langevin_tau " // params%langevin_tau // " <= 0.0")
  endif

  params%T = params%T_initial

end subroutine get_params

subroutine print_params(params)
  type(md_params), intent(in) :: params

  call print("md_params%metapot_init_args='" // trim(params%metapot_init_args) // "'")
  call print("md_params%metapot_calc_args='" // trim(params%metapot_calc_args) // "'")
  call print("md_params%atoms_in_file='" // trim(params%atoms_in_file) // "'")
  call print("md_params%params_in_file='" // trim(params%params_in_file) // "'")
  call print("md_params%trajectory_out_file='" // trim(params%trajectory_out_file) // "'")
  call print("md_params%cutoff_buffer='" // params%cutoff_buffer // "'")
  call print("md_params%rng_seed=" // params%rng_seed)
  call print("md_params%N_steps=" // params%N_steps)
  call print("md_params%dt=" // params%dt)
  call print("md_params%T=" // params%T_initial)
  if(params%T_final.lt.0.0) then
   call print("md_params%T_final=" //  params%T_initial)
  else
   call print("md_params%T_final=" // params%T_final)
   call print("md_params%freq_DT=" // params%freq_DT)
  endif
  call print("md_params%const_T=" // params%const_T)
  if (params%const_T) then
    call print("md_params%langevin_tau=" // params%langevin_tau)
    call print("md_params%rescale_velocity=" // params%rescale_velocity)
  endif
  call print("md_params%annealing=" // params%annealing)
  call print("md_params%damping=" // params%damping)
  if(params%damping) then
    call print("md_params%damping_time=" // params%damping_time)
  endif
  call print("md_params%const_P=" // params%const_P)
  call print("md_params%p_ext=" // params%p_ext)
  call print("md_params%thermalise_wait_time=" // params%thermalise_wait_time)
  call print("md_params%calc_virial=" // params%calc_virial)
  call print("md_params%summary_interval=" // params%summary_interval)
  call print("md_params%params_print_interval=" // params%params_print_interval)
  call print("md_params%at_print_interval=" // params%at_print_interval)
  call print("md_params%pot_print_interval=" // params%pot_print_interval)
  call print("md_params%zero_momentum=" // params%zero_momentum)
  call print("md_params%zero_angular_momentum=" // params%zero_angular_momentum)
  call print("md_params%calc_dipole_moment=" // params%calc_dipole_moment)
  call print("md_params%quiet_calc=" // params%quiet_calc)
  call print("md_params%do_timing=" // params%do_timing)
  call print("md_params%advance_md_substeps=" // params%advance_md_substeps)
end subroutine print_params

subroutine print_usage()
  call Print('Usage: md <command line arguments>', ERROR)
  call Print('available parameters (in md_params file or on command line):', ERROR)
  call Print('  metapot_init_args="args" [metapot_calc_args="args"]', ERROR)
  call Print('  [atoms_in_file=file(stdin)] [params_in_file=file(quip_params.xml)]', ERROR)
  call Print('  [trajectory_out_file=file(traj.xyz)] [cutoff_buffer=(0.5)] [rng_seed=n(none)]', ERROR)
  call Print('  [N_steps=n(1)] [dt=dt(1.0)] [const_T=logical(F)] [T=T(0.0)] [langevin_tau=tau(100.0)]', ERROR)
  call Print('  [calc_virial=logical(F)]', ERROR)
  call Print('  [summary_interval=n(1)] [params_print_interval=n(-1)] [at_print_interval=n(100)] [pot_print_interval=n(-1)]', ERROR)
  call Print('  [zero_momentum=T/F(F)] [zero_angular_momentum=T/F(F)] [dipole_moment=T/F(F)]', ERROR)
  call Print('  [quiet_calc=T/F(T)] [do_timing=T/F(F)] [advance_md_substeps=N(-1)', ERROR)
end subroutine print_usage

subroutine do_prints(params, ds, e, metapot, traj_out, i_step, override_intervals)
  type(md_params), intent(in) :: params
  type(DynamicalSystem), intent(inout) :: ds
  real(dp), intent(in) :: e
  type(metapotential), intent(in) :: metapot
  type(Cinoutput), optional, intent(inout) :: traj_out
  integer, intent(in) :: i_step
  logical, optional :: override_intervals

  logical my_override_intervals

  my_override_intervals = optional_default(.false., override_intervals)

  if (params%summary_interval > 0) then
    if (my_override_intervals .or. mod(i_step,params%summary_interval) == 0) call print_summary(params, ds, e)
  endif
  if (params%params_print_interval > 0) then
    if (my_override_intervals .or. mod(i_step,params%params_print_interval) == 0) call print_atoms_params(params, ds%atoms)
  endif
  if (params%pot_print_interval > 0) then
    if (my_override_intervals .or. mod(i_step,params%pot_print_interval) == 0) call print_pot(params, metapot)
  endif
  if (params%at_print_interval > 0) then
      if (my_override_intervals .or. mod(i_step,params%at_print_interval) == 0) &
	call print_at(params, ds, e, metapot, traj_out, real_format='f18.10')
  endif
end subroutine

subroutine print_summary(params, ds, e)
  type(md_params), intent(in) :: params
  type(DynamicalSystem), intent(in) :: ds
  real(dp), intent(in) :: e

  real(dp) :: mu(3)

  call ds_print_status(ds, "STAT", e)
  if (get_value(ds%atoms%params, 'Dipole_Moment', mu)) then
    call print("MU " // ds%t // " " // mu)
  endif
end subroutine print_summary

subroutine print_atoms_params(params, at)
  type(md_params), intent(in) :: params
  type(Atoms), intent(in) :: at
  call print("PA " // write_string(at%params, real_format='f18.10'))
end subroutine print_atoms_params

subroutine print_pot(params, metapot)
  type(md_params), intent(in) :: params
  type(metapotential), intent(in) :: metapot

  mainlog%prefix="MD_POT"
  if (metapot%is_simple) then
    if (associated(metapot%pot%tb)) then
      call print(metapot%pot%tb)
    endif
  endif
  mainlog%prefix=""
end subroutine print_pot

subroutine print_at(params, ds, e, metapot, out, real_format)
  type(md_params), intent(in) :: params
  type(DynamicalSystem), intent(inout) :: ds
  real(dp), intent(in) :: e
  type(metapotential), intent(in) :: metapot
  type(Cinoutput), intent(inout) :: out
  character(len=*), intent(in), optional :: real_format

  call write(out, ds%atoms)
  !call print_xyz(ds%atoms, out, all_properties=.true., comment="t="//ds%t//" e="//(kinetic_energy(ds)+e), real_format=real_format)
end subroutine print_at

subroutine initialise_md_thermostat(ds, params, do_rescale)
  type(md_params), intent(inout) :: params
  type(DynamicalSystem), intent(inout) :: ds
  logical, optional :: do_rescale
  integer           :: N_iter

  if(params%T_final.gt.0.0_dp) then
    params%annealing = .true.
    N_iter = abs(params%T_initial - params%T_final)/params%freq_DT 
    params%N_steps = N_iter*nint(params%thermalise_wait_time/params%dt)  
    call print('Number of steps : ' //  params%N_steps)
    params%freq_DT = sign(params%freq_DT, params%T_final - params%T_initial)
    call print('Annealing from ' //  params%T_initial // "K to " // params%T_final // "K")
    call print('Temperature increment ' // params%freq_DT)
  endif

  if (params%damping) then
    call print('MD damping')
    call enable_damping(ds, params%damping_time)
  elseif (params%const_T.and..not.params%const_P) then
    call print('Running NVT at T = ' // params%T // " K")
    call add_thermostat(ds, LANGEVIN, params%T, tau=params%langevin_tau)
    ds%atoms%thermostat_region = 1
    if(present(do_rescale).and..not.do_rescale) then
      call print('Velocities not rescaled') 
    elseif(params%rescale_velocity) then
      call rescale_velo(ds, params%T, mass_weighted=.true., zero_L=.true.)
    endif
  elseif (params%const_T.and.params%const_P) then
    call print('Running NPT at T = '// params%T // " K and external p = " // params%p_ext )
    call add_thermostat(ds,  LANGEVIN_NPT, params%T, tau=params%langevin_tau, p=params%p_ext)
    ds%atoms%thermostat_region = 1
    if (.not.params%calc_virial) then
       params%calc_virial = .true.
       call print('Set true calc_virial option')
    endif
    if(present(do_rescale).and..not.do_rescale) then
      call print('Velocities not rescaled') 
    elseif(params%rescale_velocity) then
      call rescale_velo(ds, params%T, mass_weighted=.true., zero_L=.true.)
    endif
  elseif (params%const_P.and..not.params%const_T) then
    call system_abort('Const_P and not Const_T')
  else
    call print('Running NVE')
    if (params%T > 0.0_dp) then
      call print("Rescaling initial velocities to T=" // params%T)
      call rescale_velo(ds, params%T, mass_weighted=.true., zero_L=.true.)
    endif
  endif

end subroutine initialise_md_thermostat

subroutine update_thermostat(ds, params, do_rescale)
  type(md_params), intent(inout) :: params
  type(DynamicalSystem), intent(inout) :: ds
  logical, optional :: do_rescale

  if (params%const_T.and..not.params%const_P) then
    call print('Running NVT at T = ' // params%T // " K")
  else
    call print('Running NPT at T = '// params%T // " K and external p = " // params%p_ext )
  endif
  ds%thermostat%T  = params%T 
  if(present(do_rescale).and..not.do_rescale) then
    call print('Velocities not rescaled')
  elseif(params%rescale_velocity) then
    call rescale_velo(ds, params%T, mass_weighted=.true., zero_L=.true.)
  endif
end subroutine update_thermostat

end module md_module

program md
use libatoms_module
use quip_module
use md_module
use libatoms_misc_utils_module

implicit none
  type (Potential) :: pot
  type (MetaPotential) :: metapot
  type(MPI_context) :: mpi_glob
  type(extendable_str) :: es
  type(Cinoutput) :: traj_out
  type(Atoms) :: at_in
  type(DynamicalSystem) :: ds

  real(dp) :: E, virial(3,3)
  real(dp), pointer :: forces(:,:)
  real(dp) :: cutoff_buffer, max_moved, last_state_change_time
  real(dp), pointer :: local_dn(:)
  real(dp) :: mu(3)
  real(dp) :: toll_time = 0.001_dp

  integer i, i_step
  type(md_params) :: params

  call system_initialise()

  call initialise(mpi_glob)

  call get_params(params, mpi_glob)

  if (params%do_timing) call enable_timing()
  call system_timer("md_prep")

  call print_params(params)

  if (params%rng_seed >= 0) call system_reseed_rng(params%rng_seed)

  call read_xyz(at_in, params%atoms_in_file, mpi_comm=mpi_glob%communicator)

  call potential_initialise_filename(pot, params%metapot_init_args, params%params_in_file, mpi_obj=mpi_glob)
  call initialise(metapot, "Simple", pot, mpi_obj=mpi_glob)

  call print(metapot)

  call initialise(ds, at_in)
  call finalise(at_in)

  call initialise(traj_out, params%trajectory_out_file, OUTPUT)

  call initialise_md_thermostat(ds, params) 

  if (params%zero_momentum) call zero_momentum(ds)
  if (params%zero_angular_momentum) call zero_angular_momentum(ds%atoms)

  ! add properties
  call add_property(ds%atoms, 'forces', 0.0_dp, n_cols=3)
  if (params%calc_dipole_moment) then
    call add_property(ds%atoms, 'local_dn', 0.0_dp, 1)
  endif

  ! make pointers to properties
  if (.not. assign_pointer(ds%atoms, 'forces', forces)) &
    call system_abort('Impossible failure to assign_ptr for forces')
  if (params%calc_dipole_moment) then
    if (.not. assign_pointer(ds%atoms, 'local_dn', local_dn)) &
      call system_abort('failure to assign pointer for local_dn')
  endif

  cutoff_buffer=params%cutoff_buffer
  call set_cutoff(ds%atoms, cutoff(metapot)+cutoff_buffer)

  ! start with p(t), v(t)
  ! calculate f(t)
  call calc_connect(ds%atoms)
  if (params%quiet_calc) call verbosity_push_decrement()
  if (params%calc_virial) then
    call calc(metapot, ds%atoms, e=E, f=forces, virial=virial, args_str=params%metapot_calc_args)
  else
    call calc(metapot, ds%atoms, e=E, f=forces, args_str=params%metapot_calc_args)
  endif
  if (params%quiet_calc) call verbosity_pop()
  call set_value(ds%atoms%params, 'time', ds%t)
  call do_prints(params, ds, e, metapot, traj_out, 0, override_intervals = .true.)

  ! calculate a(t) from f(t)
  forall(i = 1:ds%N) ds%atoms%acc(:,i) = forces(:,i) / ElementMass(ds%atoms%Z(i))

  call calc_connect(ds%atoms)
  max_moved = 0.0_dp

  last_state_change_time = ds%t

  call system_timer("md_prep")

  ! on entry, we have p(t), v(t), a(t), like advance verlet 1 wants
  call system_timer("md_loop")
  do i_step=1, params%N_steps

    if (params%annealing.and.ds%t - last_state_change_time >= params%thermalise_wait_time) then
        call print('Rescaling velocities at time '//ds%t//' from '//params%T// &
                  ' K  to '//(params%T + params%freq_DT)//' K.')
        last_state_change_time = ds%t
        params%T = params%T + params%freq_DT
        call update_thermostat(ds, params, do_rescale=(ds%cur_temp.gt.params%T))
    endif

    call advance_md(ds, params, metapot, forces, virial, E)

    ! now we have p(t+dt), v(t+dt), a(t+dt)

    if (params%calc_dipole_moment) then
      mu = dipole_moment(ds%atoms%pos, local_dn)
      call set_value(ds%atoms%params, 'Dipole_Moment', mu)
    endif

    call system_timer("md/print")
    call set_value(ds%atoms%params, 'time', ds%t)
    call do_prints(params, ds, e, metapot, traj_out, i_step)

    call system_timer("md/print")

  end do
  call system_timer("md_loop")

  call do_prints(params, ds, e, metapot, traj_out, params%N_steps, override_intervals = .true.)

  call system_finalise()

contains

  subroutine advance_md(ds, params, metapot, forces, virial, E)
    type(DynamicalSystem), intent(inout) :: ds
    type(md_params), intent(in) :: params
    type(MetaPotential), intent(inout) :: metapot
    real(dp), intent(inout) :: forces(:,:), virial(3,3)
    real(dp), intent(inout) :: E

    integer :: i_substep
    real(dp), pointer :: new_pos(:,:), new_velo(:,:), new_forces(:,:)
    real(dp) :: new_E
    logical :: has_new_pos, has_new_velo, has_new_forces, has_new_E

    if (params%advance_md_substeps > 0) then
      call calc_connect(ds%atoms)
      if (params%quiet_calc) call verbosity_push_decrement()
      if (params%calc_virial) then
	call calc(metapot, ds%atoms, e=E, f=forces, virial=virial, args_str=params%metapot_calc_args // &
	  'do_md md_time_step='//params%dt // ' md_n_steps='//params%advance_md_substeps)
      else
	call calc(metapot, ds%atoms, e=E, f=forces, args_str=params%metapot_calc_args // &
	  'do_md md_time_step='//params%dt // ' md_n_steps='//params%advance_md_substeps)
      endif
      if (params%quiet_calc) call verbosity_pop()
      has_new_pos = assign_pointer(ds%atoms, 'new_pos', new_pos)
      has_new_velo = assign_pointer(ds%atoms, 'new_velo', new_velo)
      has_new_forces = assign_pointer(ds%atoms, 'new_forces', new_forces)
      has_new_E = get_value(ds%atoms%params, 'New_Energy', new_E)
      if (count ( (/ has_new_pos, has_new_velo, has_new_forces /) ) > 0) then
	if (count ( (/ has_new_pos, has_new_velo, has_new_forces /) ) /= 3) then
	    call system_abort("advance_md tried to do md within driver, got only some of " // &
	     "has_new_pos=" // has_new_pos // " has_new_velo="//has_new_velo//" has_new_forces="//has_new_forces)
	  endif
	ds%atoms%pos = new_pos
	ds%atoms%velo = new_velo
	forces = new_forces
	if (has_new_E) E = new_E
	ds%t = ds%t + params%dt*params%advance_md_substeps
	ds%nSteps = ds%nSteps + params%advance_md_substeps
	call calc_connect(ds%atoms)
	max_moved = 0.0_dp
      else ! failed to find new_pos
	do i_substep=1, params%advance_md_substeps
	  call advance_md_one(ds, params, metapot, forces, virial, E)
	end do
      endif
    else
      call advance_md_one(ds, params, metapot, forces, virial, E)
    endif

  end subroutine advance_md


  subroutine advance_md_one(ds, params, metapot, forces, virial, E)
    type(DynamicalSystem), intent(inout) :: ds
    type(md_params), intent(in) :: params
    type(MetaPotential), intent(inout) :: metapot
    real(dp), intent(inout) :: forces(:,:), virial(3,3)
    real(dp), intent(inout) :: E

    ! first Verlet half-step
    call system_timer("md/advance_verlet1")
    if(params%const_P) then
       call advance_verlet1(ds, params%dt, forces, virial=virial)
    else
       call advance_verlet1(ds, params%dt, forces)
    endif
    call system_timer("md/advance_verlet1")
    ! now we have p(t+dt), v(t+dt/2), a(t)

    max_moved = max_moved + params%dt*maxval(abs(ds%atoms%velo))*sqrt(3.0_dp)
    call system_timer("md/calc_connect")
    if (max_moved > 0.9_dp*cutoff_buffer) then
      call calc_connect(ds%atoms)
      max_moved = 0.0_dp
    else
      call calc_dists(ds%atoms)
    endif
    call system_timer("md/calc_connect")

    ! calc f(t+dt)
    call system_timer("md/calc")
    if (params%quiet_calc) call verbosity_push_decrement()
    if (params%calc_virial) then
      call calc(metapot, ds%atoms, e=E, f=forces, virial=virial, args_str=params%metapot_calc_args)
    else
      call calc(metapot, ds%atoms, e=E, f=forces, args_str=params%metapot_calc_args)
    endif
    if (params%quiet_calc) call verbosity_pop()
    call system_timer("md/calc")
    ! now we have a(t+dt)

    ! second Verlet half-step
    call system_timer("md/advance_verlet2")
    if(params%const_P) then
       call advance_verlet2(ds, params%dt, forces, virial=virial)
    else
       call advance_verlet2(ds, params%dt, forces)
    endif
    call system_timer("md/advance_verlet2")

  end subroutine advance_md_one

end program
