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

program mc
  use libatoms_module
  use potential_module
  implicit none
  
  type(Potential) :: pot
  type(MPI_context) :: mpi_glob
  type(Cinoutput) :: atoms_out_file, atoms_in_file
  type(Atoms) :: at
  type(Dictionary) :: params
  type(extendable_str) :: params_es
  character(len=FIELD_LENGTH) :: atoms_in_filename, atoms_out_filename, params_in_filename
  character(len=FIELD_LENGTH) :: pot_init_args, pot_calc_args, first_pot_calc_args
  character(len=FIELD_LENGTH), allocatable :: atoms_print_property_list(:)
  character(len=FIELD_LENGTH) :: atoms_print_property_list_str
  integer :: N_steps, summary_interval, atoms_print_interval, rng_seed, atoms_print_property_list_n, i, rejections
  integer :: error = ERROR_NONE
  real(dp) :: E, Enew, T, proposal_amplitude
  logical :: do_timing, quiet_calc, fix_Roo
  real(dp), allocatable :: oldpos(:,:)

  call system_initialise()
  call initialise(mpi_glob)
  
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! get program parameters
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  call initialise(params)
  call param_register(params, 'atoms_in_file', 'stdin', atoms_in_filename, help_string="Initial atomic data filename")
  call param_register(params, 'atoms_out_file', 'mc.xyz', atoms_out_filename, help_string="Output filename for MC samples")
  call param_register(params, 'params_in_file', 'quip_params.xml', params_in_filename, help_string="QUIP XML parameter filename")
  call param_register(params, 'rng_seed', '-1', rng_seed, help_string="Random seed")
  call param_register(params, 'N_steps', '1', N_steps,  help_string="Number of MC steps to perform")
  call param_register(params, 'T', '-1.0', T, "Simulation temperature (Kelvin)")
  call param_register(params, 'pot_init_args', PARAM_MANDATORY, pot_init_args, help_string="string to initialise potential")
  call param_register(params, 'summary_interval', '1', summary_interval, help_string="how often to print summary line")
  call param_register(params, 'atoms_print_interval', '100', atoms_print_interval, help_string="how often to print atomic config to log file")
  call param_register(params, 'print_property_list', '', atoms_print_property_list_str, help_string="list of properties to print for atoms")
  call param_register(params, 'pot_calc_args', '', pot_calc_args, help_string="args string for potential calc")
  call param_register(params, 'first_pot_calc_args', '', first_pot_calc_args, help_string="args string for first potential calc")
  call param_register(params, 'quiet_calc', 'T', quiet_calc, help_string="do calc() quietly")
  call param_register(params, 'do_timing', 'F', do_timing, help_string="if true, do timing")
  call param_register(params, 'use_fortran_random', 'F', system_use_fortran_random, help_string="if true, use fortran builtin random_number() routine")
  call param_register(params, 'proposal_amplitude', '1.0', proposal_amplitude, help_string="multiplies the proposal step size")
  call param_register(params, 'fix_Roo', 'F', fix_Roo, help_string="if true, does not change the initial Roo distance. defaults to false")

  ! parse params
  if (.not. param_read_args(params)) then
    call param_print_help(params)
    call system_abort("Error reading params from command line")
  endif

  if (len(trim(pot_init_args)) == 0) then
     call param_print_help(params)
     call system_abort("get_params() got empty pot_init_args")
  end if
  call finalise(params)

  if (len_trim(first_pot_calc_args) == 0) first_pot_calc_args = pot_calc_args

  if (len_trim(atoms_print_property_list_str) > 0) then
    allocate(atoms_print_property_list(500))
    call split_string_simple(trim(atoms_print_property_list_str), atoms_print_property_list, atoms_print_property_list_n, ':')
    deallocate(atoms_print_property_list)
    allocate(atoms_print_property_list(atoms_print_property_list_n))
    call split_string_simple(trim(atoms_print_property_list_str), atoms_print_property_list, atoms_print_property_list_n, ':')
  else
    if (allocated(atoms_print_property_list)) deallocate(atoms_print_property_list)
  endif

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Initialise simulation
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  ! Check temperature
  if (T <= 0.0_dp) then
     call system_abort("Cannot run at negative temperature: "//T)
  end if

  ! Turn on timing
  if (do_timing) call enable_timing()

  ! set random seed if asked for
  if (rng_seed >= 0) call system_reseed_rng(rng_seed)

  ! Read in atoms
  call initialise(atoms_in_file, atoms_in_filename, INPUT, mpi=mpi_glob)
  call read(at, atoms_in_file, error=error)
  HANDLE_ERROR(error)

  ! Read in QUIP parameters and initialise potential
  if (len_trim(params_in_filename) > 0) then
     call read(params_es, params_in_filename, convert_to_string=.true., mpi_comm=mpi_glob%communicator)
  else
     call initialise(params_es)
  endif
  call initialise(pot, args_str=pot_init_args, param_str=string(params_es), mpi_obj=mpi_glob)
  call print(pot)

  call finalise(params_es)

  call initialise(atoms_out_file, atoms_out_filename, OUTPUT, mpi=mpi_glob)

  call set_cutoff(at, cutoff(pot))


  ! first calc
  call calc_connect(at)
  if (quiet_calc) call verbosity_push_decrement()
  call calc(pot, at, args_str=trim(first_pot_calc_args)//" energy")
  call get_param_value(at, "energy", E)
  if (quiet_calc) call verbosity_pop()

  ! print header
  if(summary_interval > 0) then
     call print("   Iteration      Energy      Acceptance ratio")
  end if


  i = 0
  rejections = 0
  allocate(oldpos(3,at%N))
  call system_timer("mc_loop")
  do while(i < N_steps)
     ! propose new move
     oldpos = at%pos
     call propose_move(at)
     ! calculate new energy
     call calc_connect(at)

     if(mod(i, atoms_print_interval) == 0) then
!        if(fix_Roo) then
           call calc(pot, at, args_str=trim(pot_calc_args)//" energy force")
!        end if
        if (allocated(atoms_print_property_list)) then
           call write(atoms_out_file, at, properties_array=atoms_print_property_list, real_format='%18.10f')
        else
           call write(atoms_out_file, at, real_format='%18.10f')
        endif

     end if


     if (quiet_calc) call verbosity_push_decrement()
     call calc(pot, at, args_str=trim(pot_calc_args)//" energy")
     call get_param_value(at, "energy", Enew)
     if (quiet_calc) call verbosity_pop()

     ! Metropolis test
     if( Enew > E ) then
        if ( ran_uniform() > exp(-(Enew-E)/(BOLTZMANN_K*T)) ) then
           !reject the move
           at%pos = oldpos
           rejections = rejections+1
        else
           E = Enew ! accept
        end if
     else
        E = Enew ! accept
     end if

     ! print summary
     if(summary_interval > 0) then
        if(mod(i, summary_interval) == 0) then
           call print("MC "//i//" "//E//" "//(1.0_dp - real(rejections,dp)/real(summary_interval,dp)))
           rejections = 0
        end if
     end if


     i = i+1
  end do
  call system_timer("mc_loop")

  deallocate(oldpos)

  call system_finalise()

contains

  subroutine propose_move(at)
    type(Atoms), intent(inout) :: at
    type(Quaternion) :: quat
    real(dp) :: dr(3), u1, u2, u3, axis(3)
    real(dp), parameter :: drOHamp = 0.05_dp, drOOamp = 0.5_dp, dHOHamp = 1.0_dp

    ! Atoms are in the order OHHOHH

    ! first choose between monomer rotations vs. internal changes coupled with OO distance changes

    if(ran_uniform() < 0.5) then
       
       ! change OO distance
       if(.not. fix_Roo) then
          dr = (at%pos(:,4)-at%pos(:,1))/norm(at%pos(:,4)-at%pos(:,1))*(ran_uniform()-0.5_dp)*drOOamp*proposal_amplitude
          at%pos(:,4) = at%pos(:,4)+dr
          at%pos(:,5) = at%pos(:,5)+dr
          at%pos(:,6) = at%pos(:,6)+dr
       end if
       
       ! change OH distances
       at%pos(:,2) = at%pos(:,2) + (at%pos(:,2)-at%pos(:,1))/norm(at%pos(:,2)-at%pos(:,1))*(ran_uniform()-0.5_dp)*drOHamp*proposal_amplitude
       at%pos(:,3) = at%pos(:,3) + (at%pos(:,3)-at%pos(:,1))/norm(at%pos(:,3)-at%pos(:,1))*(ran_uniform()-0.5_dp)*drOHamp*proposal_amplitude
       at%pos(:,5) = at%pos(:,5) + (at%pos(:,5)-at%pos(:,4))/norm(at%pos(:,5)-at%pos(:,4))*(ran_uniform()-0.5_dp)*drOHamp*proposal_amplitude
       at%pos(:,6) = at%pos(:,6) + (at%pos(:,6)-at%pos(:,4))/norm(at%pos(:,6)-at%pos(:,4))*(ran_uniform()-0.5_dp)*drOHamp*proposal_amplitude
       
       ! change OH angles
       axis = (at%pos(:,2)-at%pos(:,1)) .cross. (at%pos(:,3)-at%pos(:,1))
       quat = rotation(axis, dHOHamp*proposal_amplitude/180.0_dp*PI*(ran_uniform()-0.5_dp))
       at%pos(:,2) = at%pos(:,1)+ rotate(at%pos(:,2)-at%pos(:,1), quat)
       
       axis = (at%pos(:,5)-at%pos(:,4)) .cross. (at%pos(:,6)-at%pos(:,4))
       quat = rotation(axis, dHOHamp*proposal_amplitude/180.0_dp*PI*(ran_uniform()-0.5_dp))
       at%pos(:,5) = at%pos(:,4)+ rotate(at%pos(:,5)-at%pos(:,4), quat)
       
    else
       ! uniform random rotation. formula from:
       ! K. Shoemake. Uniform random rotations. In D. Kirk, editor, Graphics Gems III, pages 124-132. Academic, New York, 1992.
       ! via  http://planning.cs.uiuc.edu/node198.html

       ! rotate first water monomer around its oxygen
       u1 = ran_uniform()
       u2 = ran_uniform()
       u3 = ran_uniform()
       call initialise(quat, sqrt(1.0_dp-u1)*sin(2.0_dp*PI*u2), sqrt(1.0_dp-u1)*cos(2.0_dp*PI*u2), sqrt(u1)*sin(2.0_dp*PI*u3), sqrt(u1)*cos(2.0_dp*PI*u3)) ! quaternion representing random rotation
       at%pos(:,2) = at%pos(:,1)+rotate(at%pos(:,2)-at%pos(:,1),quat)
       at%pos(:,3) = at%pos(:,1)+rotate(at%pos(:,3)-at%pos(:,1),quat)
       
       ! rotate second water monomer around its oxygen
       u1 = ran_uniform()
       u2 = ran_uniform()
       u3 = ran_uniform()
       call initialise(quat, sqrt(1.0_dp-u1)*sin(2.0_dp*PI*u2), sqrt(1.0_dp-u1)*cos(2.0_dp*PI*u2), sqrt(u1)*sin(2.0_dp*PI*u3), sqrt(u1)*cos(2.0_dp*PI*u3)) ! quaternion representing random rotation
       at%pos(:,5) = at%pos(:,4)+rotate(at%pos(:,5)-at%pos(:,4),quat)
       at%pos(:,6) = at%pos(:,4)+rotate(at%pos(:,6)-at%pos(:,4),quat)
    end if
  end subroutine propose_move

end program mc
