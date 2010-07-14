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

program test_potential

use libAtoms_module
use QUIP_module
use libatoms_misc_utils_module
use elasticity_module
use phonons_module


implicit none

  type(Potential) pot1, pot2
  type(MetaPotential) metapot
  type(MPI_context) mpi_glob
  type(Atoms) at

  type(Dictionary) :: cli_params

  character(len=100) verbosity
  logical :: do_hybrid
  character(len=FIELD_LENGTH) init_args, calc_args, at_file, param_file, init_args_pot1, init_args_pot2

  real(dp) :: E0
  real(dp), allocatable :: local_E_p(:), local_E_m(:), local_N_p(:), local_N_m(:)
  real(dp), pointer :: local_E_fd(:), local_N(:), local_N_fd(:)
  integer :: fd_index
  real(dp) :: fd_vec(3), p0(3)

  call system_initialise()

  init_args = ''
  init_args_pot1 = ''
  init_args_pot2 = ''
  calc_args = ''
  call initialise(cli_params)
  call param_register(cli_params, 'fd_index', PARAM_MANDATORY, fd_index)
  call param_register(cli_params, 'fd_vec', PARAM_MANDATORY, fd_vec)
  call param_register(cli_params, 'init_args', PARAM_MANDATORY, init_args)
  call param_register(cli_params, 'at_file', 'stdin', at_file)
  call param_register(cli_params, 'param_file', 'quip_params.xml', param_file)
  call param_register(cli_params, 'calc_args', '', calc_args)
  call param_register(cli_params, 'verbosity', 'NORMAL', verbosity)
  call param_register(cli_params, 'hybrid', 'F', do_hybrid)
  call param_register(cli_params, 'init_args_pot1', '', init_args_pot1)
  call param_register(cli_params, 'init_args_pot2', '', init_args_pot2)

  call print("n_args " // cmd_arg_count())

  if (.not. param_read_args(cli_params, do_check = .true.)) then
    call print("Usage: eval fd_index=i fd_vec='x y z'", PRINT_ALWAYS)
    call print("  init_args='str' [at_file=file(stdin)] [param_file=file(quip_parms.xml)]",PRINT_ALWAYS)
    call print("  [calc_args='str'] [verbosity=VERBOSITY(PRINT_NORMAL)]", PRINT_ALWAYS)
    call print("  [hybrid] [init_args_pot1] [init_args_pot2]", PRINT_ALWAYS)
    call system_abort("Confused by CLI arguments")
  end if
  call finalise(cli_params)

  call print ("Using init args " // trim(init_args))
  call print ("Using calc args " // trim(calc_args))
  if (do_hybrid) then
    call print("Hybrid using init args pot1 " // trim(init_args_pot1))
    call print("Hybrid using init args pot2 " // trim(init_args_pot2))
  endif

  call Initialise(mpi_glob)

  call read_xyz(at, at_file, mpi_comm=mpi_glob%communicator)
  if (do_hybrid) then
    call Potential_Initialise_filename(pot1, init_args_pot1, param_file, mpi_obj=mpi_glob)
    call Potential_Initialise_filename(pot2, init_args_pot2, param_file, mpi_obj=mpi_glob)
    call Initialise(metapot, init_args, pot1, pot2, mpi_obj=mpi_glob)
  else
    call Potential_Initialise_filename(pot1, init_args, param_file, mpi_obj=mpi_glob)
    call Initialise(metapot, "Simple", pot1, mpi_obj=mpi_glob)
  endif

  select case(verbosity)
    case ("NORMAL")
      call verbosity_push(PRINT_NORMAL)
    case ("VERBOSE")
      call verbosity_push(PRINT_VERBOSE)
    case ("NERD")
      call verbosity_push(PRINT_NERD)
    case ("ANAL")
      call verbosity_push(PRINT_ANAL)
    case default
      call system_abort("confused by verbosity " // trim(verbosity))
  end select

  call set_cutoff(at, cutoff(metapot)+0.5_dp)
  call calc_connect(at)


  call calc(metapot, at, e = E0, args_str = calc_args)

  allocate(local_E_p(at%N))
  allocate(local_E_m(at%N))
  allocate(local_N_p(at%N))
  allocate(local_N_m(at%N))

  call add_property(at, 'local_E_fd', 0.0_dp)
  call add_property(at, 'local_N', 0.0_dp)
  call add_property(at, 'local_N_fd', 0.0_dp)

  p0 = at%pos(:,fd_index)

  at%pos(:,fd_index) = p0 + fd_vec
  call calc_connect(at)
  call calc(metapot, at, local_e = local_E_p, args_str = trim(calc_args) // " do_at_local_N")
  if (.not. assign_pointer(at, "local_N", local_N)) &
    call system_abort("Impossible: Failed to assign pointer for local_N")
  local_N_p = local_N

  at%pos(:,fd_index) = p0 - fd_vec
  call calc_connect(at)
  call calc(metapot, at, local_e = local_E_m, args_str = trim(calc_args) // " do_at_local_N")
  if (.not. assign_pointer(at, "local_N", local_N)) &
    call system_abort("Impossible: Failed to assign pointer for local_N")
  local_N_m = local_N

  if (.not. assign_pointer(at, 'local_E_fd', local_E_fd)) &
    call system_abort("Impossible failure to assign pointer for local_E_fd")
  if (.not. assign_pointer(at, 'local_N_fd', local_N_fd)) &
    call system_abort("Impossible failure to assign pointer for local_N_fd")

  local_E_fd = (local_E_p-local_E_m)/(2.0_dp*norm(fd_vec))
  local_N_fd = (local_N_p-local_N_m)/(2.0_dp*norm(fd_vec))

  mainlog%prefix="LOCAL_E_FD"
  call print_xyz(at, mainlog, all_properties=.true.)
  mainlog%prefix=""

  call system_finalise()
  stop

end program
