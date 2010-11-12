! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   Portions of this code were written by
! H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! H0 X
! H0 X   Copyright 2006-2010.
! H0 X
! H0 X   These portions of the source code are released under the GNU General
! H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
! H0 X
! H0 X   If you would like to license the source code under different terms,
! H0 X   please contact Gabor Csanyi, gabor@csanyi.net
! H0 X
! H0 X   Portions of this code were written by Noam Bernstein as part of
! H0 X   his employment for the U.S. Government, and are not subject
! H0 X   to copyright in the USA.
! H0 X
! H0 X
! H0 X   When using this software, please cite the following reference:
! H0 X
! H0 X   http://www.libatoms.org
! H0 X
! H0 X  Additional contributions by
! H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! H0 X
! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"
module md_module
use libatoms_module
use potential_module
implicit none
private

  public :: md_params
  type md_params
    character(len=FIELD_LENGTH) :: atoms_in_file, params_in_file, trajectory_out_file
    integer :: N_steps
    real(dp) :: max_time
    real(dp) :: dt,  thermalise_wait_time, damping_time
    real(dp) :: T_initial, T, T_final, freq_DT, langevin_tau, p_ext
    real(dp) :: velocity_rescaling_time 
    real(dp) :: cutoff_buffer 
    integer :: velocity_rescaling_freq
    logical :: calc_virial, calc_energy, const_T, const_P
    character(len=FIELD_LENGTH) :: pot_init_args, pot_calc_args, first_pot_calc_args
    integer :: summary_interval, params_print_interval, at_print_interval, pot_print_interval
    character(len=FIELD_LENGTH), allocatable :: print_property_list(:)
    integer :: rng_seed
    logical :: damping, rescale_velocity 
    logical :: annealing, zero_momentum, zero_angular_momentum
    logical :: quiet_calc, do_timing
    integer :: advance_md_substeps
    logical :: v_dep_quants_extra_calc
    logical :: continuation
  end type md_params

public :: get_params, print_params, print_usage, do_prints, initialise_md_thermostat, update_thermostat

contains

subroutine get_params(params, mpi_glob)
  type(md_params), intent(inout) :: params
  type(MPI_context), intent(in) :: mpi_glob

  type(Dictionary) :: md_params_dict
  type(Extendable_Str) :: es
  logical :: md_params_exist
  character(len=FIELD_LENGTH) :: print_property_list_str
  integer :: n_print_property_list
  logical :: has_N_steps

  call initialise(md_params_dict)
  call param_register(md_params_dict, 'atoms_in_file', 'stdin', params%atoms_in_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'trajectory_out_file', 'traj.xyz', params%trajectory_out_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'params_in_file', 'quip_params.xml', params%params_in_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'rng_seed', '-1', params%rng_seed, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'N_steps', '1', params%N_steps, has_value_target=has_N_steps, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'max_time', '-1.0', params%max_time, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'dt', '1.0', params%dt, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'T', '0.0', params%T_initial, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'T_final', '-100.0', params%T_final, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'freq_DT', '10.0', params%freq_DT, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'p_ext', '0.0', params%p_ext, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'thermalise_wait_time', '10.0', params%thermalise_wait_time, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'velocity_rescaling_freq', '2', params%velocity_rescaling_freq, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'velocity_rescaling_time', '0.0', params%velocity_rescaling_time, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'const_T', 'F', params%const_T, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'const_P', 'F', params%const_P, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'annealing', 'F', params%annealing, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'damping', 'F', params%damping, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'damping_time', '10.0', params%damping_time, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'rescale_velocity', 'F', params%rescale_velocity, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'langevin_tau', '100.0', params%langevin_tau, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'calc_virial', 'F', params%calc_virial, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'calc_energy', 'T', params%calc_energy, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'pot_init_args', PARAM_MANDATORY, params%pot_init_args, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'cutoff_buffer', '0.5', params%cutoff_buffer, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'summary_interval', '1', params%summary_interval, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'params_print_interval', '-1', params%params_print_interval, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'at_print_interval', '100', params%at_print_interval, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'print_property_list', '', print_property_list_str, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'pot_print_interval', '-1', params%pot_print_interval, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'zero_momentum', 'F', params%zero_momentum, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'zero_angular_momentum', 'F', params%zero_angular_momentum, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'pot_calc_args', '', params%pot_calc_args, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'first_pot_calc_args', '', params%first_pot_calc_args, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'quiet_calc', 'T', params%quiet_calc, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'do_timing', 'F', params%do_timing, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'advance_md_substeps', '-1', params%advance_md_substeps, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'v_dep_quants_extra_calc', 'F', params%v_dep_quants_extra_calc, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(md_params_dict, 'continuation', 'F', params%continuation, help_string="No help yet.  This source file was $LastChangedBy$")

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

  if (len(trim(params%pot_init_args)) == 0) then
     call print_usage()
     call system_abort("get_params got empty pot_init_args")
  end if

  if (len_trim(params%first_pot_calc_args) == 0) params%first_pot_calc_args = params%pot_calc_args


  if (has_N_steps .and. params%N_steps > 0 .and. params%max_time > 0.0_dp) then
    call system_abort("get_params got both N_steps="//params%N_steps//" > 0 and max_time="//params%max_time//" > 0.0")
  endif
  if (params%max_time <= 0.0_dp .and. params%N_steps < 1) call system_abort("get_params got max_time="//params%max_time//" <= 0.0 and N_steps=" // params%N_steps // " < 1")
  if (params%dt <= 0.0_dp) call system_abort("get_params got dt " // params%dt // " <= 0.0")
  if (params%const_T) then
    if (params%langevin_tau <= 0.0_dp) call system_abort("get_params got const_T, but langevin_tau " // params%langevin_tau // " <= 0.0")
  endif

  if (len_trim(print_property_list_str) > 0) then
    allocate(params%print_property_list(500))
    call split_string_simple(trim(print_property_list_str), params%print_property_list, n_print_property_list, ':')
    deallocate(params%print_property_list)
    allocate(params%print_property_list(n_print_property_list))
    call split_string_simple(trim(print_property_list_str), params%print_property_list, n_print_property_list, ':')
  else
    if (allocated(params%print_property_list)) deallocate(params%print_property_list)
  endif

  params%T = params%T_initial

end subroutine get_params

subroutine print_params(params)
  type(md_params), intent(in) :: params

  integer :: i

  call print("md_params%pot_init_args='" // trim(params%pot_init_args) // "'")
  call print("md_params%pot_calc_args='" // trim(params%pot_calc_args) // "'")
  call print("md_params%first_pot_calc_args='" // trim(params%first_pot_calc_args) // "'")
  call print("md_params%atoms_in_file='" // trim(params%atoms_in_file) // "'")
  call print("md_params%params_in_file='" // trim(params%params_in_file) // "'")
  call print("md_params%trajectory_out_file='" // trim(params%trajectory_out_file) // "'")
  call print("md_params%cutoff_buffer='" // params%cutoff_buffer // "'")
  call print("md_params%rng_seed=" // params%rng_seed)
  call print("md_params%N_steps=" // params%N_steps)
  call print("md_params%max_time=" // params%max_time)
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
  call print("md_params%calc_energy=" // params%calc_energy)
  call print("md_params%summary_interval=" // params%summary_interval)
  call print("md_params%params_print_interval=" // params%params_print_interval)
  call print("md_params%at_print_interval=" // params%at_print_interval)
  call print("md_params%print_property_list=", nocr=.true.)
  if (allocated(params%print_property_list)) then
    if (size(params%print_property_list) > 0) then
      call print(trim(params%print_property_list(1)), nocr=.true.)
      do i=2, size(params%print_property_list)
	call print(":"//trim(params%print_property_list(i)), nocr=.true.)
      end do
    endif
  endif
  call print("", nocr=.false.)
  call print("md_params%pot_print_interval=" // params%pot_print_interval)
  call print("md_params%zero_momentum=" // params%zero_momentum)
  call print("md_params%zero_angular_momentum=" // params%zero_angular_momentum)
  call print("md_params%quiet_calc=" // params%quiet_calc)
  call print("md_params%do_timing=" // params%do_timing)
  call print("md_params%advance_md_substeps=" // params%advance_md_substeps)
  call print("md_params%v_dep_quants_extra_calc=" // params%v_dep_quants_extra_calc)
  call print("md_params%continuation=" // params%continuation)
end subroutine print_params

subroutine print_usage()
  call Print('Usage: md <command line arguments>', PRINT_ALWAYS)
  call Print('available parameters (in md_params file or on command line):', PRINT_ALWAYS)
  call Print('  pot_init_args="args" [pot_calc_args="args"] [first_pot_calc_args="args"]', PRINT_ALWAYS)
  call Print('  [atoms_in_file=file(stdin)] [params_in_file=file(quip_params.xml)]', PRINT_ALWAYS)
  call Print('  [trajectory_out_file=file(traj.xyz)] [cutoff_buffer=(0.5)] [rng_seed=n(none)]', PRINT_ALWAYS)
  call Print('  [N_steps=n(1)] [max_time=t(-1.0)] [dt=dt(1.0)] [const_T=logical(F)] [T=T(0.0)] [langevin_tau=tau(100.0)]', PRINT_ALWAYS)
  call Print('  [calc_virial=logical(F)]', PRINT_ALWAYS)
  call Print('  [summary_interval=n(1)] [params_print_interval=n(-1)] [at_print_interval=n(100)]', PRINT_ALWAYS)
  call Print('  [print_property_list=prop1:prop2:...()] [pot_print_interval=n(-1)]', PRINT_ALWAYS)
  call Print('  [zero_momentum=T/F(F)] [zero_angular_momentum=T/F(F)]', PRINT_ALWAYS)
  call Print('  [quiet_calc=T/F(T)] [do_timing=T/F(F)] [advance_md_substeps=N(-1)] [v_dep_quants_extra_calc=T/F(F)]', PRINT_ALWAYS)
  call Print('  [continuation=T/F(F)]', PRINT_ALWAYS)
end subroutine print_usage

subroutine do_prints(params, ds, e, pot, restraint_stuff, restraint_stuff_timeavg, traj_out, i_step, override_intervals)
  type(md_params), intent(in) :: params
  type(DynamicalSystem), intent(inout) :: ds
  real(dp), intent(in) :: e
  type(potential), intent(in) :: pot
  real(dp), allocatable, intent(in) :: restraint_stuff(:,:), restraint_stuff_timeavg(:,:)
  type(Cinoutput), optional, intent(inout) :: traj_out
  integer, intent(in) :: i_step
  logical, optional :: override_intervals

  logical my_override_intervals

  my_override_intervals = optional_default(.false., override_intervals)

  if (params%summary_interval > 0) then
    if (my_override_intervals .or. mod(i_step,params%summary_interval) == 0) call print_summary(params, ds, e)
  endif
  if (allocated(restraint_stuff)) then
     if (params%summary_interval > 0) then
	if (my_override_intervals .or. mod(i_step, params%summary_interval) == 0) &
	   call print_restraint_stuff(params, ds, restraint_stuff, restraint_stuff_timeavg)
     endif
  endif

  if (params%params_print_interval > 0) then
    if (my_override_intervals .or. mod(i_step,params%params_print_interval) == 0) call print_atoms_params(params, ds%atoms)
  endif

  if (params%pot_print_interval > 0) then
    if (my_override_intervals .or. mod(i_step,params%pot_print_interval) == 0) call print_pot(params, pot)
  endif

  if (params%at_print_interval > 0) then
      if (my_override_intervals .or. mod(i_step,params%at_print_interval) == 0) &
	call print_at(params, ds, e, pot, traj_out)
  endif

end subroutine

subroutine print_restraint_stuff(params, ds, restraint_stuff, restraint_stuff_timeavg)
  type(md_params), intent(in) :: params
  type(DynamicalSystem), intent(in) :: ds
  real(dp), intent(in) :: restraint_stuff(:,:), restraint_stuff_timeavg(:,:)

  logical, save :: firstcall = .true.

  if (firstcall) then
    call print("#RI t  target_v C E dE/dcoll dE/dk ....")
    call print("#R t  target_v C E dE/dcoll dE/dk ....")
    firstcall = .false.
  endif
  call print("RI " // ds%t // " " // reshape( restraint_stuff, (/ size(restraint_stuff) /) ))
  call print("R " // ds%t // " " // reshape( restraint_stuff_timeavg, (/ size(restraint_stuff_timeavg) /) ))
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

subroutine print_pot(params, pot)
  type(md_params), intent(in) :: params
  type(potential), intent(in) :: pot

  mainlog%prefix="MD_POT"
#ifdef HAVE_TB
  if (pot%is_simple) then
    if (associated(pot%simple%tb)) then
      call print(pot%simple%tb)
    endif
  endif
#endif
  mainlog%prefix=""
end subroutine print_pot

subroutine print_at(params, ds, e, pot, out)
  type(md_params), intent(in) :: params
  type(DynamicalSystem), intent(inout) :: ds
  real(dp), intent(in) :: e
  type(Potential), intent(in) :: pot
  type(CInOutput), intent(inout) :: out

  if (allocated(params%print_property_list)) then
    call write(out, ds%atoms, properties_array=params%print_property_list, real_format='%18.10f')
  else
    call write(out, ds%atoms, real_format='%18.10f')
  endif
  !call print_xyz(ds%atoms, out, all_properties=.true., comment="t="//ds%t//" e="//(kinetic_energy(ds)+e), real_format='f18.10')
end subroutine print_at

subroutine initialise_md_thermostat(ds, params, do_rescale)
  type(md_params), intent(inout) :: params
  type(DynamicalSystem), intent(inout) :: ds
  logical, optional :: do_rescale
  integer           :: N_iter

  if (params%T_final.gt.0.0_dp) then
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
use potential_module
use md_module
use libatoms_misc_utils_module
use restraints_constraints_xml_module

implicit none
  type (Potential) :: pot
  type(MPI_context) :: mpi_glob
  type(Cinoutput) :: traj_out, atoms_in_cio
  type(Atoms) :: at_in
  type(DynamicalSystem) :: ds

  real(dp) :: E, virial(3,3)
  real(dp), pointer :: force_p(:,:)
  real(dp) :: cutoff_buffer, max_moved, last_state_change_time

  integer :: i, i_step, initial_i_step
  type(md_params) :: params
  integer :: error = ERROR_NONE
  type(extendable_str) :: params_es

  logical :: store_constraint_force
  real(dp), allocatable :: restraint_stuff(:,:), restraint_stuff_timeavg(:,:)
  character(STRING_LENGTH) :: extra_calc_args
  integer :: l_error

  call system_initialise()

  call initialise(mpi_glob)

  call get_params(params, mpi_glob)

  if (params%do_timing) call enable_timing()
  call system_timer("md_prep")

  call print_params(params)

  if (params%rng_seed >= 0) call system_reseed_rng(params%rng_seed)

  call initialise(atoms_in_cio, params%atoms_in_file, INPUT, mpi=mpi_glob)
  call read(at_in, atoms_in_cio, error=error)
  HANDLE_ERROR(error)

  call read(params_es, params%params_in_file, convert_to_string=.true., mpi_comm=mpi_glob%communicator)

  call initialise(pot, args_str=params%pot_init_args, param_str=string(params_es), mpi_obj=mpi_glob)

  call print(pot)

  call initialise(ds, at_in)
  call finalise(at_in)

  ! set some initial values, in particular if this run is a continuation
  initial_i_step = 1
  if (params%continuation) then
    call get_param_value(ds%atoms, 'time', ds%t, l_error)
    call get_param_value(ds%atoms, 'i_step', initial_i_step, l_error)
    CLEAR_ERROR(error)
  endif

  call init_restraints_constraints(ds, string(params_es))
  store_constraint_force = has_property(ds%atoms, "constraint_force")
  if (ds%Nrestraints > 0) then
     allocate(restraint_stuff(5,ds%Nrestraints))
     allocate(restraint_stuff_timeavg(5,ds%Nrestraints))
  endif

  call finalise(params_es)

  call initialise(traj_out, params%trajectory_out_file, OUTPUT)

  call initialise_md_thermostat(ds, params) 

  if (params%zero_momentum) call zero_momentum(ds)
  if (params%zero_angular_momentum) call zero_angular_momentum(ds%atoms)

  ! add properties
  call add_property(ds%atoms, 'force', 0.0_dp, n_cols=3, ptr2=force_p)

  cutoff_buffer=params%cutoff_buffer
  call set_cutoff(ds%atoms, cutoff(pot)+cutoff_buffer)

  ! start with p(t), v(t)
  ! calculate f(t)
  call calc_connect(ds%atoms)
  if (params%quiet_calc) call verbosity_push_decrement()

  extra_calc_args="force"
  if (params%calc_energy) extra_calc_args = trim(extra_calc_args) // " energy"
  if (params%calc_virial) extra_calc_args = trim(extra_calc_args) // " virial"
  call calc(pot, ds%atoms, args_str=trim(params%first_pot_calc_args)//" "//trim(extra_calc_args))
  if (params%calc_energy) call get_param_value(ds%atoms, "energy", E)
  if (params%calc_virial) call get_param_value(ds%atoms, "virial", virial)
  if (params%quiet_calc) call verbosity_pop()
  call set_value(ds%atoms%params, 'time', ds%t)

  ! calculate a(t) from f(t)
  forall(i = 1:ds%N) ds%atoms%acc(:,i) = force_p(:,i) / ElementMass(ds%atoms%Z(i))

  if (ds%Nrestraints > 0) then
     call get_restraint_stuff(ds, restraint_stuff)
     restraint_stuff_timeavg = restraint_stuff
  end if
  call do_prints(params, ds, e, pot, restraint_stuff, restraint_stuff_timeavg, traj_out, 0, override_intervals = .true.)

  call calc_connect(ds%atoms)
  max_moved = 0.0_dp

  last_state_change_time = ds%t

  call system_timer("md_prep")

  ! on entry, we have p(t), v(t), a(t), like advance verlet 1 wants
  call system_timer("md_loop")
  i_step = initial_i_step
  do while ((params%N_steps > 0 .and. i_step <= params%N_steps) .or. (params%max_time > 0.0_dp .and. ds%t <= params%max_time))

    if (params%annealing.and.ds%t - last_state_change_time >= params%thermalise_wait_time) then
        call print('Rescaling velocities at time '//ds%t//' from '//params%T// &
                  ' K  to '//(params%T + params%freq_DT)//' K.')
        last_state_change_time = ds%t
        params%T = params%T + params%freq_DT
        call update_thermostat(ds, params, do_rescale=(ds%cur_temp.gt.params%T))
    endif

    call advance_md(ds, params, pot, store_constraint_force)
    if (ds%Nrestraints > 0) then
       call get_restraint_stuff(ds, restraint_stuff)
       call update_exponential_average(restraint_stuff_timeavg, params%dt/ds%avg_time, restraint_stuff)
    end if
    if (params%calc_energy) call get_param_value(ds%atoms, "energy", E)

    ! now we have p(t+dt), v(t+dt), a(t+dt)

    call system_timer("md/print")
    call set_value(ds%atoms%params, 'time', ds%t)
    call set_value(ds%atoms%params, 'i_step', i_step)
    call do_prints(params, ds, e, pot, restraint_stuff, restraint_stuff_timeavg, traj_out, i_step)

    call system_timer("md/print")

    i_step = i_step + 1
  end do
  call system_timer("md_loop")

  call do_prints(params, ds, e, pot, restraint_stuff, restraint_stuff_timeavg, traj_out, params%N_steps, override_intervals = .true.)

  call system_finalise()

contains

  subroutine get_restraint_stuff(ds, restraint_stuff)
    type(DynamicalSystem), intent(in) :: ds
    real(dp), intent(out) :: restraint_stuff(:,:)

    integer i_r

    do i_r = 1, ds%Nrestraints
       restraint_stuff(1,i_r) =  ds%restraint(i_r)%target_v
       restraint_stuff(2,i_r) =  ds%restraint(i_r)%C
       restraint_stuff(3,i_r) =  ds%restraint(i_r)%E
       restraint_stuff(4,i_r) = -ds%restraint(i_r)%dE_dcoll
       restraint_stuff(5,i_r) = -ds%restraint(i_r)%dE_dk
    end do
  end subroutine

  subroutine advance_md(ds, params, pot, store_constraint_force)
    type(DynamicalSystem), intent(inout) :: ds
    type(md_params), intent(in) :: params
    type(Potential), intent(inout) :: pot
    logical, intent(in) :: store_constraint_force

    integer :: i_substep
    real(dp), pointer :: new_pos(:,:), new_velo(:,:), new_forces(:,:), force_p(:,:)
    real(dp) :: new_E
    logical :: has_new_pos, has_new_velo, has_new_forces, has_new_E
    character(STRING_LENGTH) :: extra_calc_args

    if (params%advance_md_substeps > 0) then
      call calc_connect(ds%atoms)
      if (params%quiet_calc) call verbosity_push_decrement()

      extra_calc_args="force"
      if (params%calc_energy) extra_calc_args = trim(extra_calc_args) // " energy"
      if (params%calc_virial) extra_calc_args = trim(extra_calc_args) // " virial"
      call calc(pot, ds%atoms, args_str=trim(params%pot_calc_args)//" "//trim(extra_calc_args))

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
	if (.not. assign_pointer(ds%atoms, "force", force_p)) then
	  call system_abort("advance_md failed to assign pointer for force")
	endif
	force_p = new_forces
	if (has_new_E) call set_param_value(ds%atoms, "energy", new_E)
	ds%t = ds%t + params%dt*params%advance_md_substeps
	ds%nSteps = ds%nSteps + params%advance_md_substeps
	call calc_connect(ds%atoms)
	max_moved = 0.0_dp
      else ! failed to find new_pos
	do i_substep=1, params%advance_md_substeps
	  call advance_md_one(ds, params, pot, store_constraint_force)
	end do
      endif
    else
      call advance_md_one(ds, params, pot, store_constraint_force)
    endif
  end subroutine advance_md

  subroutine advance_md_one(ds, params, pot, store_constraint_force)
    type(DynamicalSystem), intent(inout) :: ds
    type(md_params), intent(in) :: params
    type(Potential), intent(inout) :: pot
    logical, intent(in) :: store_constraint_force

    real(dp), pointer :: force_p(:,:)
    real(dp) :: E, virial(3,3)
    character(STRING_LENGTH) :: extra_calc_args

    ! start with have p(t), v(t), a(t)

    ! first Verlet half-step
    call system_timer("md/advance_verlet1")
    if(params%const_P) then
       call advance_verlet1(ds, params%dt, virial=virial, store_constraint_force=store_constraint_force)
    else
       call advance_verlet1(ds, params%dt, store_constraint_force=store_constraint_force)
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

    extra_calc_args="force"
    if (params%calc_energy) extra_calc_args = trim(extra_calc_args) // " energy"
    if (params%calc_virial) extra_calc_args = trim(extra_calc_args) // " virial"
    call calc(pot, ds%atoms, args_str=trim(params%pot_calc_args)//" "//trim(extra_calc_args))
    if (params%calc_energy) call get_param_value(ds%atoms, "energy", E)
    if (params%calc_virial) call get_param_value(ds%atoms, "virial", virial)
    if (.not. assign_pointer(ds%atoms, "force", force_p)) call system_abort("md failed to get force")

    if (params%quiet_calc) call verbosity_pop()
    call system_timer("md/calc")
    ! now we have a(t+dt)

    ! second Verlet half-step
    call system_timer("md/advance_verlet2")
    if(params%const_P) then
       call advance_verlet2(ds, params%dt, force_p, virial=virial, E=E, store_constraint_force=store_constraint_force)
    else
       call advance_verlet2(ds, params%dt, force_p, E=E, store_constraint_force=store_constraint_force)
    endif
    call system_timer("md/advance_verlet2")

    ! now we have p(t+dt), v(t+dt), a(t+dt)

    ! call calc again if needed for v dep. forces
    if (params%v_dep_quants_extra_calc) then
      if (params%quiet_calc) call verbosity_push_decrement()
      call calc(pot, ds%atoms, args_str=trim(params%pot_calc_args)//" "//trim(extra_calc_args))
      if (params%quiet_calc) call verbosity_pop()
    end if
    call system_timer("md/calc")

  end subroutine advance_md_one

end program
