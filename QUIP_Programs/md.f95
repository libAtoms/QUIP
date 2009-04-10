module md_module
use libatoms_module
use quip_module
implicit none
private

  public :: md_params
  type md_params
    character(len=FIELD_LENGTH) :: atoms_in_file, params_in_file, trajectory_out_file
    integer :: N_steps
    real(dp) :: dt
    real(dp) :: T, langevin_tau
    logical :: calc_virial, const_T
    character(len=FIELD_LENGTH) :: metapot_init_args, metapot_calc_args
    integer :: summary_interval, at_print_interval, pot_print_interval
    integer :: rng_seed
    logical :: zero_momentum, zero_angular_momentum
    logical :: calc_dipole_moment, quiet_calc, do_timing
  end type md_params

public :: get_params, print_params, print_usage, do_prints

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
  call param_register(md_params_dict, 'T', '0.0', params%T)
  call param_register(md_params_dict, 'const_T', 'F', params%const_T)
  call param_register(md_params_dict, 'langevin_tau', '100.0', params%langevin_tau)
  call param_register(md_params_dict, 'calc_virial', 'F', params%calc_virial)
  call param_register(md_params_dict, 'metapot_init_args', PARAM_MANDATORY, params%metapot_init_args)
  call param_register(md_params_dict, 'summary_interval', '1', params%summary_interval)
  call param_register(md_params_dict, 'at_print_interval', '100', params%at_print_interval)
  call param_register(md_params_dict, 'pot_print_interval', '-1', params%pot_print_interval)
  call param_register(md_params_dict, 'zero_momentum', 'F', params%zero_momentum)
  call param_register(md_params_dict, 'zero_angular_momentum', 'F', params%zero_angular_momentum)
  call param_register(md_params_dict, 'metapot_calc_args', '', params%metapot_calc_args)
  call param_register(md_params_dict, 'dipole_moment', 'F', params%calc_dipole_moment)
  call param_register(md_params_dict, 'quiet_calc', 'T', params%quiet_calc)
  call param_register(md_params_dict, 'do_timing', 'T', params%do_timing)

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

end subroutine get_params

subroutine print_params(params)
  type(md_params), intent(in) :: params

  call print("md_params%metapot_init_args='" // trim(params%metapot_init_args) // "'")
  call print("md_params%metapot_calc_args='" // trim(params%metapot_calc_args) // "'")
  call print("md_params%atoms_in_file='" // trim(params%atoms_in_file) // "'")
  call print("md_params%params_in_file='" // trim(params%params_in_file) // "'")
  call print("md_params%trajectory_out_file='" // trim(params%trajectory_out_file) // "'")
  call print("md_params%rng_seed=" // params%rng_seed)
  call print("md_params%N_steps=" // params%N_steps)
  call print("md_params%dt=" // params%dt)
  call print("md_params%T=" // params%T)
  call print("md_params%const_T=" // params%const_T)
  if (params%const_T) then
    call print("md_params%langevin_tau=" // params%langevin_tau)
  endif
  call print("md_params%calc_virial=" // params%calc_virial)
  call print("md_params%summary_interval=" // params%summary_interval)
  call print("md_params%at_print_interval=" // params%at_print_interval)
  call print("md_params%pot_print_interval=" // params%pot_print_interval)
  call print("md_params%zero_momentum=" // params%zero_momentum)
  call print("md_params%zero_angular_momentum=" // params%zero_angular_momentum)
  call print("md_params%calc_dipole_moment=" // params%calc_dipole_moment)
  call print("md_params%quiet_calc=" // params%quiet_calc)
  call print("md_params%do_timing=" // params%do_timing)
end subroutine print_params

subroutine print_usage()
  call Print('Usage: md <command line arguments>', ERROR)
  call Print('available parameters (in md_params file or on command line):', ERROR)
  call Print('  metapot_init_args="args" [metapot_calc_args="args"]', ERROR)
  call Print('  [atoms_in_file=file(stdin)] [params_in_file=file(quip_params.xml)]', ERROR)
  call Print('  [trajectory_out_file=file(traj.xyz)] [rng_seed=n(none)]', ERROR)
  call Print('  [N_steps=n(1)] [dt=dt(1.0)] [const_T=logical(F)] [T=T(0.0)] [langevin_tau=tau(100.0)]', ERROR)
  call Print('  [calc_virial=logical(F)]', ERROR)
  call Print('  [summary_interval=n(1)] [at_print_interval=n(100)] [pot_print_interval=n(-1)]', ERROR)
  call Print('  [zero_momentum=T/F(F)] [zero_angular_momentum=T/F(F)] [dipole_moment=T/F(F)]', ERROR)
  call Print('  [quiet_calc=T/F(T)] [do_timing=T/F(F)]', ERROR)
end subroutine print_usage

subroutine do_prints(params, ds, e, metapot, traj_out, i_step, override_intervals)
  type(md_params), intent(in) :: params
  type(DynamicalSystem), intent(inout) :: ds
  real(dp), intent(in) :: e
  type(metapotential), intent(in) :: metapot
  type(inoutput), optional, intent(inout) :: traj_out
  integer, intent(in) :: i_step
  logical, optional :: override_intervals

  logical my_override_intervals

  my_override_intervals = optional_default(.false., override_intervals)

  if (params%summary_interval > 0) then
    if (my_override_intervals .or. mod(i_step,params%summary_interval) == 0) call print_summary(params, ds, e)
  endif
  if (params%pot_print_interval > 0) then
    if (my_override_intervals .or. mod(i_step,params%pot_print_interval) == 0) call print_pot(params, metapot)
  endif
  if (params%at_print_interval > 0) then
    if (.not. present(traj_out)) then
      if (my_override_intervals .or. mod(i_step,params%at_print_interval) == 0) &
	call print_at(params, ds, e, metapot, mainlog, real_format='f18.10')
    else
      if (my_override_intervals .or. mod(i_step,params%at_print_interval) == 0) &
	call print_at(params, ds, e, metapot, traj_out, real_format='f18.10')
    endif
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
  type(inoutput), intent(inout) :: out
  character(len=*), intent(in), optional :: real_format

  call print_xyz(ds%atoms, out, all_properties=.true., comment="t="//ds%t//" e="//(kinetic_energy(ds)+e), real_format=real_format)
end subroutine print_at

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
  type(inoutput) :: traj_out
  type(Atoms) :: at_in
  type(DynamicalSystem) :: ds

  real(dp) :: E, virial(3,3)
  real(dp), pointer :: forces(:,:)
  real(dp) :: cutoff_buffer, max_moved
  real(dp), pointer :: local_dn(:)
  real(dp) :: mu(3)

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

  if (params%const_T) then
    call add_thermostat(ds, LANGEVIN, params%T, tau=params%langevin_tau)
    ds%atoms%thermostat_region = 1
  else
    if (params%T > 0.0_dp) then
      call print("Rescaling initial velocities to T=" // params%T)
      call rescale_velo(ds, params%T, mass_weighted=.true., zero_L=.true.)
    endif
  endif

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

  cutoff_buffer=0.5_dp
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
  call do_prints(params, ds, e, metapot, traj_out, 0, override_intervals = .true.)

  ! calculate a(t) from f(t)
  forall(i = 1:ds%N) ds%atoms%acc(:,i) = forces(:,i) / ElementMass(ds%atoms%Z(i))

  call calc_connect(ds%atoms)
  max_moved = 0.0_dp

  call system_timer("md_prep")

  ! on entry, we have p(t), v(t), a(t), like advance verlet 1 wants
  call system_timer("md_loop")
  do i_step=1, params%N_steps

    ! first Verlet half-step
    call system_timer("md/advance_verlet1")
    call advance_verlet1(ds, params%dt, forces)
    call system_timer("md/advance_verlet1")
    ! now we have p(t+dt), v(t+dt/2), a(t)

    max_moved = max_moved + maxval(abs(ds%atoms%velo))*sqrt(3.0_dp)
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
    call advance_verlet2(ds, params%dt, forces)
    call system_timer("md/advance_verlet2")

    ! now we have p(t+dt), v(t+dt), a(t+dt)

    if (params%calc_dipole_moment) then
      mu = dipole_moment(ds%atoms%pos, local_dn)
      call set_value(ds%atoms%params, 'Dipole_Moment', mu)
    endif

    call system_timer("md/print")
    call do_prints(params, ds, e, metapot, traj_out, i_step)

    call system_timer("md/print")
  end do
  call system_timer("md_loop")

  call do_prints(params, ds, e, metapot, traj_out, params%N_steps, override_intervals = .true.)

  call system_finalise()

end program
