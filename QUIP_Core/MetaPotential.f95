!% MetaPotential module to abstract hybrid schemes such as LOTF and the
!% local energy scheme into a Potential like object to facilitate
!% activities such as geometry optimisation.
!%
!% From the outside, a MetaPotential looks rather like a standard
!% QUIP Potential, with interfaces for 'Initialise', 'Finalise',
!% 'Print', 'Cutoff' and 'Calc'.
!%
!% It can also be used to geometry optimise an Atoms structure,
!% using the 'Minim' interface.
!%
!% MetaPotential_LOTF also has to live in this module so that it can call
!% MetaPotential_minim() without any circular depenendancies.

#ifdef HAVE_LOCAL_E_MIX
#define HAVE_HYBRID
#endif
#ifdef HAVE_ONIOM
#define HAVE_HYBRID
#endif

module MetaPotential_module

  use libAtoms_module
  use QUIP_Common_module
  use MPI_context_module
  use Potential_module

  use TB_module

#ifdef HAVE_LOTF
  use adjustablepotential_module, only: adjustable_potential_init, adjustable_potential_optimise, &
       adjustable_potential_force, adjustable_potential_finalise
!  use adjustablepotential_sw_module, only: adjustable_potential_sw_init, adjustable_potential_sw_optimise, &
!       adjustable_potential_sw_force, adjustable_potential_sw_finalise
#endif

  implicit none

  private


  !*************************************************************************
  !*
  !*  MetaPotential header
  !*
  !*************************************************************************

  integer :: n_iter = 0

  public :: MetaPotential
  type MetaPotential
     type(MPI_context) :: mpi
     
     logical :: is_simple
     type(Potential), pointer :: pot => null()


     logical :: is_sum
     type(MetaPotential_Sum), pointer :: sum => null()

     logical :: is_forcemixing
     type(MetaPotential_FM), pointer :: forcemixing => null()

#ifdef HAVE_LOCAL_E_MIX
     logical :: is_local_e_mix
     type(MetaPotential_Local_E_Mix), pointer :: local_e_mix => null()
#endif

#ifdef HAVE_ONIOM
     logical :: is_oniom
     type(MetaPotential_ONIOM), pointer :: oniom => null()
#endif
  end type MetaPotential

  public :: MetaPotential_Sum
  type MetaPotential_Sum
     type(MPI_context) :: mpi

     type(Potential), pointer :: pot1 => null() 
     type(Potential), pointer :: pot2 => null() 

     logical  :: subtract_pot1
     logical  :: subtract_pot2
  end type MetaPotential_Sum

  public :: Initialise
  interface Initialise
     module procedure MetaPotential_Initialise, MetaPotential_Sum_Initialise
  end interface

  public :: Finalise
  interface Finalise
     module procedure MetaPotential_Finalise, MetaPotential_Sum_Finalise
  end interface

  public :: Print
  interface Print
     module procedure MetaPotential_Print, MetaPotential_Sum_Print
  end interface

  public :: Cutoff
  interface Cutoff
     module procedure MetaPotential_Cutoff, MetaPotential_Sum_Cutoff
  end interface

  public :: Calc
  interface Calc
     module procedure MetaPotential_Calc, MetaPotential_Sum_Calc
  end interface

  public :: Minim
  interface Minim
     module procedure MetaPotential_Minim
  end interface

  public :: test_gradient
  interface test_gradient
     module procedure metapot_test_gradient
  end interface

  public :: n_test_gradient
  interface n_test_gradient
     module procedure metapot_n_test_gradient
  end interface

  public :: metapotential_minimise
  type metapotential_minimise
     real(dp) :: minim_pos_lat_preconditioner = 1.0_dp

     real(dp) :: minim_save_lat(3,3)
     character(len=1024) :: minim_args_str
     logical :: minim_do_pos, minim_do_lat
     integer :: minim_n_eval_e, minim_n_eval_f, minim_n_eval_ef
     real(dp) :: pos_lat_preconditioner_factor
     type(Atoms), pointer :: minim_at => null()
     real(dp), pointer :: last_connect_x(:) => null()
     type(MetaPotential), pointer :: minim_metapot => null()
     type(Inoutput), pointer :: minim_inoutput_movie => null()
     type(CInoutput), pointer :: minim_cinoutput_movie => null()
     logical, dimension(3,3) :: lattice_fix
  end type metapotential_minimise


#include "MetaPotential_ForceMixing_header.f95"

#ifdef HAVE_LOCAL_E_MIX
#include "MetaPotential_Local_E_Mix_header.f95"
#endif
#ifdef HAVE_ONIOM
#include "MetaPotential_ONIOM_header.f95"
#endif

  contains

  !*************************************************************************
  !*
  !*  MetaPotential routines
  !*
  !*************************************************************************

subroutine metapotential_initialise(this, args_str, pot, pot2, bulk_scale, mpi_obj)
  type(MetaPotential), intent(inout) :: this
  character(len=*), intent(in) :: args_str !% Valid arguments are 'Simple', 'Sum', 'ForceMixing', 'Local_E_Mix' and 'ONIOM'
  type(Potential), intent(in), target :: pot !% Potential upon which this MetaPotential is based
  type(Potential), optional, intent(in), target :: pot2 !% Optional second potential
  type(Atoms), optional, intent(inout) :: bulk_scale !% optional bulk structure for calculating space and E rescaling
  type(MPI_Context), intent(in), optional :: mpi_obj

  type(Dictionary) :: params


  call finalise(this)
  call initialise(params)
  call param_register(params, 'Simple', 'false', this%is_simple)
  call param_register(params, 'Sum', 'false', this%is_sum)
  call param_register(params, 'ForceMixing', 'false', this%is_forcemixing)
#ifdef HAVE_LOCAL_E_MIX
  call param_register(params, 'Local_E_Mix', 'false', this%is_local_e_mix)
#endif HAVE_LOCAL_E_MIX
#ifdef HAVE_ONIOM
  call param_register(params, 'ONIOM', 'false', this%is_oniom)
#endif HAVE_ONIOM

  if(.not. param_read_line(params, args_str, ignore_unknown=.true.)) then
    call system_abort("MetaPotential_initialise failed to parse args_str='"//trim(args_str)//"'")
  end if
  call finalise(params)

  if (present(mpi_obj)) this%mpi = mpi_obj

  if (count( (/this%is_simple, this%is_sum, this%is_forcemixing &
#ifdef HAVE_LOCAL_E_MIX
          , this%is_local_e_mix &
#endif HAVE_LOCAL_E_MIX
#ifdef HAVE_ONIOM
          , this%is_oniom &
#endif HAVE_ONIOM
          /) ) /= 1) then
  call system_abort("MetaPotential_initialise found too few or two many MetaPotential types args_str='"&
      //trim(args_str)//"'")
  end if

  if (this%is_simple) then
    if (present(bulk_scale)) call print("MetaPotential_initialise Simple ignoring bulk_scale passed in", ERROR)

    this%pot => pot

  else if (this%is_sum) then
    if (present(bulk_scale)) call print("MetaPotential_initialise Sum ignoring bulk_scale passed in", ERROR)

    if (.not. present(pot2)) &
      call system_abort('MetaPotential_initialise: two potentials needs for sum metapotential')

    allocate(this%sum)
    call initialise(this%sum, args_str, pot, pot2, mpi_obj)

  else if (this%is_forcemixing) then

    allocate(this%forcemixing)
    if(present(pot2)) then
       call initialise(this%forcemixing, args_str, pot, pot2, bulk_scale, mpi_obj)
    else
       ! if only one pot is given, assign it to QM. this is useful for time-embedding LOTF
       call initialise(this%forcemixing, args_str, qmpot=pot, reference_bulk=bulk_scale, mpi=mpi_obj)
    endif

#ifdef HAVE_LOCAL_E_MIX
  else if (this%is_local_e_mix) then
    if (.not. present(pot2)) &
      call system_abort('MetaPotential_initialise: two potentials needs for local_e_mix metapotential')

    allocate(this%local_e_mix)
    call initialise(this%local_e_mix, args_str, pot, pot2, bulk_scale, mpi_obj)
#endif HAVE_LOCAL_E_MIX
#ifdef HAVE_ONIOM
  else if (this%is_oniom) then
    if (.not. present(pot2)) &
      call system_abort('MetaPotential_initialise: two potentials needs for oniom metapotential')

    allocate(this%oniom)
    call initialise(this%oniom, args_str, pot, pot2, bulk_scale, mpi_obj)
#endif HAVE_ONIOM
  end if

  end subroutine metapotential_initialise

  subroutine metapotential_finalise(this)
    type(MetaPotential), intent(inout) :: this

    if (this%is_simple) then
       nullify(this%pot)
    else if (this%is_sum) then
       call finalise(this%sum)
       deallocate(this%sum)
    else if (this%is_forcemixing) then
       call finalise(this%forcemixing)
       deallocate(this%forcemixing)
#ifdef HAVE_LOCAL_E_MIX
    else if (this%is_local_e_mix) then
       call finalise(this%local_e_mix)
       deallocate(this%local_e_mix)
#endif HAVE_LOCAL_E_MIX
#ifdef HAVE_ONIOM
    else if (this%is_oniom) then
       call finalise(this%oniom)
       deallocate(this%oniom)
#endif HAVE_ONIOM
    end if

    this%is_simple = .false.
    this%is_sum = .false.
    this%is_forcemixing = .false.
#ifdef HAVE_LOCAL_E_MIX
    this%is_local_e_mix = .false.
#endif
#ifdef HAVE_ONIOM
    this%is_oniom = .false.
#endif

  end subroutine metapotential_finalise


  subroutine metapotential_calc(this, at, e, local_e, f, df, virial, args_str, err)
    type(MetaPotential), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    real(dp), intent(out), optional :: e
    real(dp), intent(out), optional :: local_e(:)
    real(dp), intent(out), optional :: f(:,:)
    real(dp), intent(out), optional :: df(:,:)
    real(dp), intent(out), optional :: virial(3,3)
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: err

    if (present(err)) err = 0

    if (this%is_simple) then
       call Calc(this%pot, at, e, local_e, f, df, virial, args_str, err)
    else if (this%is_sum) then
       call Calc(this%sum, at, e, local_e, f, df, virial, args_str, err)
    else if (this%is_forcemixing) then
       call Calc(this%forcemixing, at, e, local_e, f, virial, args_str, err)
#ifdef HAVE_LOCAL_E_MIX
    else if (this%is_local_e_mix) then
      call Calc(this%local_e_mix, at, e, local_e, f, virial, args_str, err)
#endif HAVE_local_e_mix
#ifdef HAVE_ONIOM
    else if (this%is_oniom) then
      call Calc(this%oniom, at, e, local_e, f, virial, args_str, err)
#endif HAVE_ONIOM
    else
      call system_abort('MetaPotential_Calc: meta potential is not initialised')
    endif

  end subroutine

  subroutine metapotential_print(this, file)
    type(MetaPotential), intent(inout) :: this
    type(Inoutput), intent(inout), optional :: file

    if (this%is_simple) then
       call Print('MetaPotential containing potential')
       call Print(this%pot, file=file)
    else if (this%is_sum) then
       call Print(this%sum, file=file)
    else if (this%is_forcemixing) then
       call Print(this%forcemixing, file=file)
#ifdef HAVE_LOCAL_E_MIX
    else if (this%is_local_e_mix) then
       call Print(this%local_e_mix, file=file)
#endif HAVE_LOCAL_E_MIX
#ifdef HAVE_ONIOM
    else if (this%is_oniom) then
       call Print(this%oniom, file=file)
#endif HAVE_oniom
    else
       call system_abort('MetaPotential_Print: metapotential is not initialised')
    end if
  end subroutine metapotential_print


  function metapotential_cutoff(this)
    type(MetaPotential), intent(in) :: this
    real(dp) :: metapotential_cutoff

    if (this%is_simple) then
       metapotential_cutoff = cutoff(this%pot)
    else if (this%is_sum) then
       metapotential_cutoff = cutoff(this%sum)
    else if (this%is_forcemixing) then
       metapotential_cutoff = cutoff(this%forcemixing)
#ifdef HAVE_LOCAL_E_MIX
    else if (this%is_local_e_mix) then
       metapotential_cutoff = cutoff(this%local_e_mix)
#endif HAVE_local_e_mix
#ifdef HAVE_ONIOM
    else if (this%is_oniom) then
       metapotential_cutoff = cutoff(this%oniom)
#endif HAVE_ONIOM
    else
       call system_abort('MetaPotential_Cutoff: metapotential is not initialised')
    end if
  end function metapotential_cutoff



  !*************************************************************************
  !*
  !*  Minimisation routines
  !*
  !*************************************************************************


  function metapotential_minim(this, at, method, convergence_tol, max_steps, linminroutine, do_print, print_inoutput, print_cinoutput, &
       do_pos, do_lat, args_str, eps_guess, use_n_minim, use_fire, lattice_fix, hook_print_interval)
    type(Atoms), intent(inout), target :: at !% starting configuration
    type(Metapotential), intent(inout), target :: this !% metapotential to evaluate energy/forces with
    character(*), intent(in)    :: method !% passed to minim()
    real(dp),     intent(in)    :: convergence_tol !% Minimisation is treated as converged once $|\mathbf{\nabla}f|^2 <$
                                                    !% 'convergence_tol'. 
    integer,      intent(in)    :: max_steps  !% Maximum number of steps
    character(*), intent(in),optional    :: linminroutine !% Name of the line minisation routine to use, passed to base minim()
    logical, optional :: do_print !% if true, print configurations using minim's hook()
    type(inoutput), intent(inout), optional, target :: print_inoutput !% inoutput object to print configs to, needed if do_print is true
    type(cinoutput), intent(inout), optional, target :: print_cinoutput !% cinoutput object to print configs to, needed if do_print is true
    logical, optional :: do_pos, do_lat !% do relaxation w.r.t. positions and/or lattice (is neither is included, do both)
    character(len=*), intent(in), optional :: args_str !% arguments to pass to calc()
    real(dp), intent(in), optional :: eps_guess !% eps_guess argument to pass to minim
    logical, intent(in), optional :: use_n_minim !% if true, use n_minim instead of minim
    logical, intent(in), optional :: use_fire   !% if true, use fire_minim instead of minim
    logical, dimension(3,3), optional :: lattice_fix !% Mask to fix some components of lattice. Defaults to all false.
    integer, intent(in), optional :: hook_print_interval !% how often to print xyz from hook function
    integer::metapotential_minim

    integer n_iter, n_iter_tot
    real(dp), allocatable :: x(:)
    real(dp) :: deform_grad(3,3)
    logical my_do_print
    logical done
    real(dp) :: my_eps_guess
    logical my_use_n_minim, my_use_fire

    real(dp) :: initial_E, final_E, mass
    type(metapotential_minimise) am
    integer :: am_data_size
    character, allocatable :: am_data(:)


    am_data_size = size(transfer(am, am_data))
    allocate(am_data(am_data_size))

    my_use_n_minim = optional_default(.false., use_n_minim)
    my_use_fire = optional_default(.false., use_fire)

    my_eps_guess = optional_default(1.0e-2_dp/at%N, eps_guess)
    if (my_eps_guess .feq. 0.0_dp) my_eps_guess = 1.0e-2_dp/at%N

    call calc_connect(at)
    am%minim_at => at
    am%pos_lat_preconditioner_factor = am%minim_pos_lat_preconditioner*am%minim_at%N

    if (present(args_str)) then
      am%minim_args_str = args_str
    else
      am%minim_args_str = ""
    endif
    am%minim_metapot => this

    if (.not.present(do_pos) .and. .not. present(do_lat)) then
      am%minim_do_pos = .true.
      am%minim_do_lat = .true.
    else
      am%minim_do_pos = optional_default(.false., do_pos)
      am%minim_do_lat = optional_default(.false., do_lat)
    endif

    am%lattice_fix = .false.
    if (present(lattice_fix)) am%lattice_fix = lattice_fix

    my_do_print = optional_default(.false., do_print)

    if (my_do_print .and. .not. present(print_inoutput) .and. .not. present(print_cinoutput)) &
         call system_abort("metapotential_minim: do_print is true, but no print_inoutput or print_cinoutput present")

    if (my_do_print) then
       if (present(print_cinoutput)) then
          am%minim_cinoutput_movie => print_cinoutput
          if (this%is_forcemixing) &
               this%forcemixing%minim_cinoutput_movie => print_cinoutput
#ifdef HAVE_LOCAL_E_MIX
          if (this%is_local_e_mix) &
               this%local_e_mix%minim_cinoutput_movie => print_cinoutput
#endif HAVE_LOCAL_E_MIX
#ifdef HAVE_ONIOM
          if (this%is_oniom) &
               this%oniom%minim_cinoutput_movie => print_cinoutput
#endif HAVE_ONIOM

       else
          am%minim_inoutput_movie => print_inoutput
          if (this%is_forcemixing) &
               this%forcemixing%minim_inoutput_movie => print_inoutput
#ifdef HAVE_LOCAL_E_MIX
          if (this%is_local_e_mix) &
               this%local_e_mix%minim_inoutput_movie => print_inoutput
#endif HAVE_LOCAL_E_MIX
#ifdef HAVE_ONIOM
          if (this%is_oniom) &
               this%oniom%minim_inoutput_movie => print_inoutput
#endif HAVE_ONIOM
       end if
    else
      nullify(am%minim_inoutput_movie)
    endif

    am%minim_n_eval_e = 0
    am%minim_n_eval_f = 0

    allocate(x(9+am%minim_at%N*3))
    allocate(am%last_connect_x(size(x)))
    am%last_connect_x=1.0e38_dp
    deform_grad = 0.0_dp; call add_identity(deform_grad)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor, am%minim_at%travel, am%minim_at%lattice)

    am_data = transfer(am, am_data)
    if (my_use_n_minim) then
       n_iter = n_minim(x, both_func, initial_E, final_E, my_eps_guess, max_steps, convergence_tol, print_hook, &
            hook_print_interval=hook_print_interval, data=am_data)
    else if (my_use_fire) then
       if (has_property(at, 'mass')) then
          mass = at%mass(1)
       else
          mass = ElementMass(at%Z(1))
       end if
       n_iter = fire_minim(x, mass, dummy_energy_func, gradient_func, 1.0_dp, convergence_tol, max_steps, &
            print_hook, hook_print_interval=hook_print_interval, data=am_data)
    else
       n_iter = minim(x, energy_func, gradient_func, method, convergence_tol, max_steps, linminroutine, &
            print_hook, hook_print_interval=hook_print_interval, eps_guess=my_eps_guess, data=am_data)
    endif
    call print("minim relax w.r.t. both n_iter " // n_iter, VERBOSE)
    call unpack_pos_dg(x, am%minim_at, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)
    call calc_connect(am%minim_at)
    n_iter_tot = n_iter
    done = .true.
    deallocate(am%last_connect_x)
    deallocate(x)
    call print("MINIM_N_EVAL E " // am%minim_n_eval_e // " F " // am%minim_n_eval_f // &
      " EF " // am%minim_n_eval_ef, VERBOSE)

    metapotential_minim = n_iter_tot

    deallocate(am_data)

  end function metapotential_minim

  subroutine undo_travel(at)
    type(Atoms), intent(inout) :: at

    if (any(at%travel /= 0)) then
      at%pos = at%pos + (at%lattice .mult. at%travel)
      at%travel = 0
    endif
  end subroutine

  ! test the gradient given a metapotential, string to pass to calc, and atoms structure
  function metapot_test_gradient(metapot, at, do_pos, do_lat, args_str, dir_field)
    type(Atoms), intent(inout), target :: at
    type(Metapotential), target, intent(inout) :: metapot
    logical, intent(in), optional :: do_pos, do_lat
    character(len=*), intent(in), optional :: args_str
    character(len=*), intent(in), optional :: dir_field
    logical metapot_test_gradient

    real(dp) :: dg(3,3)
    real(dp), allocatable :: x(:), dir_x(:)
    real(dp), pointer :: dir_p(:,:)
    integer :: am_data_size
    type(metapotential_minimise) :: am
    character, allocatable :: am_data(:)

    am_data_size = size(transfer(am, am_data))
    allocate(am_data(am_data_size))

    metapot_test_gradient = .true.

    if (present(args_str)) then
      am%minim_args_str = args_str
    else
      am%minim_args_str = ""
    endif
    am%minim_metapot => metapot

    ! neither do_pos not do_lat is specified, so do both by default
    if (.not. present(do_pos) .and. .not. present(do_lat)) then
      am%minim_do_pos= .true.
      am%minim_do_lat= .true.
    else ! something is specified, so do exactly that
      am%minim_do_pos = optional_default(.false., do_pos)
      am%minim_do_lat = optional_default(.false., do_lat)
    endif

    am%pos_lat_preconditioner_factor = am%minim_pos_lat_preconditioner*at%N

    if (present(dir_field)) then
      if (.not. assign_pointer(at, trim(dir_field), dir_p)) &
	call system_abort("metapot_test_gradient got dir_field '"//trim(dir_field)//"' but can't assign pointer to it")
      call print("metapot_test_gradient using '"//trim(dir_field)//"' field for gradient test direction")
    else
      nullify(dir_p)
    endif

    am%minim_at => at
    allocate(x(at%N*3+9))
    allocate(am%last_connect_x(at%N*3+9))
    am%last_connect_x=1.0e38_dp
    dg = 0.0_dp; call add_identity(dg)
    call pack_pos_dg(at%pos, dg, x, am%pos_lat_preconditioner_factor, at%travel, at%lattice)
    if (associated(dir_p)) then
      allocate(dir_x(at%N*3+9))
      dg = 0.0_dp
      call pack_pos_dg(dir_p, dg, dir_x, am%pos_lat_preconditioner_factor)
    endif
    am_data = transfer(am, am_data)
    if (associated(dir_p)) then
      metapot_test_gradient = test_gradient(x, energy_func, gradient_func, dir=dir_x, data=am_data)
    else
      metapot_test_gradient = test_gradient(x, energy_func, gradient_func, data=am_data)
    endif
    deallocate(x)
    if (allocated(dir_x)) deallocate(dir_x)
    deallocate(am%last_connect_x)

    deallocate(am_data)

  end function metapot_test_gradient

  ! test the gradient given a metapotential, string to pass to calc, and atoms structure
  subroutine metapot_n_test_gradient(metapot, at, do_pos, do_lat, args_str, dir_field)
    type(Metapotential), target, intent(inout) :: metapot
    type(Atoms), intent(inout), target :: at
    logical, intent(in), optional :: do_pos, do_lat
    character(len=*), intent(in), optional :: args_str
    character(len=*), intent(in), optional :: dir_field

    real(dp), allocatable :: x(:), dir_x(:)
    real(dp), pointer :: dir_p(:,:)
    real(dp) :: dg(3,3)
    integer :: am_data_size
    type(metapotential_minimise) :: am
    character, allocatable :: am_data(:)

    am_data_size = size(transfer(am, am_data))
    allocate(am_data(am_data_size))

    if (present(args_str)) then
      am%minim_args_str = args_str
    else
      am%minim_args_str = ""
    endif
    am%minim_metapot => metapot

    ! neither do_pos not do_lat is specified, so do both by default
    if (.not. present(do_pos) .and. .not. present(do_lat)) then
      am%minim_do_pos= .true.
      am%minim_do_lat= .true.
    else ! something is specified, so do exactly that
      am%minim_do_pos = optional_default(.false., do_pos)
      am%minim_do_lat = optional_default(.false., do_lat)
    endif

    am%pos_lat_preconditioner_factor = am%minim_pos_lat_preconditioner*at%N

    if (present(dir_field)) then
      if (.not. assign_pointer(at, trim(dir_field), dir_p)) &
	call system_abort("metapot_test_gradient got dir_field '"//trim(dir_field)//"' but can't assign pointer to it")
      call print("metapot_test_gradient using '"//trim(dir_field)//"' field for gradient test direction")
    else
      nullify(dir_p)
    endif

    am%minim_at => at
    allocate(x(at%N*3+9))
    allocate(am%last_connect_x(at%N*3+9))
    am%last_connect_x=1.0e38_dp
    dg = 0.0_dp; call add_identity(dg)
    call pack_pos_dg(at%pos, dg, x, am%pos_lat_preconditioner_factor, at%travel, at%lattice)
    if (associated(dir_p)) then
      allocate(dir_x(at%N*3+9))
      dg = 0.0_dp
      call pack_pos_dg(dir_p, dg, dir_x, am%pos_lat_preconditioner_factor)
    endif
    am_data = transfer(am, am_data)
    if (associated(dir_p)) then
      call n_test_gradient(x, energy_func, gradient_func, data=am_data, dir=dir_x)
    else
      call n_test_gradient(x, energy_func, gradient_func, data=am_data)
    endif
    deallocate(x)
    if (allocated(dir_x)) deallocate(dir_x)
    deallocate(am%last_connect_x)

    deallocate(am_data)

  end subroutine metapot_n_test_gradient

  ! hook for minim to print configuration during position relaxation
  subroutine print_hook(x,dx,E,done,do_print,am_data)
    real(dp), intent(in) :: x(:), dx(:)
    real(dp), intent(in) :: E
    logical, intent(out) :: done
    logical, optional, intent(in) :: do_print
    character, intent(in), optional :: am_data(:)

    real(dp), allocatable :: f(:,:)
    real(dp), pointer :: f_p(:,:)
    integer, pointer, dimension(:) :: hybrid_mark
    real(dp) :: virial(3,3), deform_grad(3,3), deform_grad_inv(3,3)
    type(metapotential_minimise) :: am
    logical :: my_do_print

    if (.not. present(am_data)) call system_abort("metapotential_minimise print_hook must have am_data")
    my_do_print = optional_default(.true., do_print)

    am = transfer(am_data, am)
    ! beware of transfer and pointers !!!
    call atoms_repoint(am%minim_at)

    if (associated(am%minim_inoutput_movie) .or. associated(am%minim_cinoutput_movie)) then
      if (size(x) /= am%minim_at%N*3+9) call system_abort("Called gradient_func() with size mismatch " // &
        size(x) // " /= " // am%minim_at%N // "*3+9")

      call unpack_pos_dg(x, am%minim_at, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
      call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

      allocate(f(3,am%minim_at%N))
      f = reshape(dx(10:), (/3,am%minim_at%N/))
      f = transpose(deform_grad) .mult. f

      virial = reshape(dx(1:9), (/3,3/) )
      call inverse(deform_grad, deform_grad_inv)
      virial = virial .mult. transpose(deform_grad_inv)

      ! Set atom parameters to print energy, maxforce and maxvirial
      call set_value(am%minim_at%params, 'E', E)
      if (am%minim_do_pos) then
         call set_value(am%minim_at%params, 'MaxForce', maxval(norm(f,1)))
         call set_value(am%minim_at%params, 'df2', norm2(reshape(f,(/3*am%minim_at%N/))))


         if (am%minim_metapot%is_forcemixing) then
            if (assign_pointer(am%minim_at, 'hybrid_mark', hybrid_mark)) then 
               if (.not. am%minim_metapot%forcemixing%minimise_mm) then
                  call set_value(am%minim_at%params, 'QM_MaxForce', &
                       maxval(abs(f(:,find(hybrid_mark == HYBRID_ACTIVE_MARK)))))
            
                  call set_value(am%minim_at%params, 'QM_df2', &
                       norm2(reshape(f(:,find(hybrid_mark == HYBRID_ACTIVE_MARK)),&
                       (/3*count(hybrid_mark == HYBRID_ACTIVE_MARK)/))))

               end if
            end if
         end if

      end if

      if (am%minim_do_lat) &
	call set_value(am%minim_at%params, 'MaxVirial', maxval(abs(virial)))

      if (assign_pointer(am%minim_at, "forces", f_p)) f_p = -f
      if (assign_pointer(am%minim_at, "force", f_p)) f_p = -f ! compatibility with crack program

      if (my_do_print) then
         if (associated(am%minim_inoutput_movie)) then
            call print_xyz(am%minim_at, am%minim_inoutput_movie, all_properties=.true., real_format='f12.5')
         end if

         if (associated(am%minim_cinoutput_movie)) then
            if (.not. am%minim_metapot%mpi%active .or. (am%minim_metapot%mpi%active .and. am%minim_metapot%mpi%my_proc == 0)) &
                 call write(am%minim_cinoutput_movie, am%minim_at)
         end if
      end if

      deallocate(f)
      call fix_atoms_deform_grad(deform_grad, am%minim_at,am)
    endif
    done = .false.

  end subroutine print_hook

  function dummy_energy_func(x, am_data)
    real(dp) :: x(:)
    character, optional :: am_data(:)
    real(dp) :: dummy_energy_func

    dummy_energy_func = 0.0_dp

  end function dummy_energy_func


  ! compute energy
  function energy_func(x, am_data)
    real(dp) :: x(:)
    character, optional :: am_data(:)
    real(dp) :: energy_func

    real(dp) :: max_atom_rij_change
    real(dp) :: deform_grad(3,3)
    type(metapotential_minimise)  :: am

    call system_timer("energy_func")

    if (.not. present(am_data)) call system_abort("metapotential_minimise energy_func must have am_data")
    am = transfer(am_data, am)
    call atoms_repoint(am%minim_at)

    am%minim_n_eval_e = am%minim_n_eval_e + 1

    if (size(x) /= am%minim_at%N*3+9) call system_abort("Called energy_func() with size mismatch " // &
      size(x) // " /= " // am%minim_at%N // "*3+9")

    ! Note: TB will return 0 for cutoff(am%minim_metapot), but TB does its own calc_connect, so doesn't matter
    max_atom_rij_change = max_rij_change(am%last_connect_x, x, cutoff(am%minim_metapot), &
      1.0_dp/am%pos_lat_preconditioner_factor)

    call print("energy_func got x " // x, NERD)

    call unpack_pos_dg(x, am%minim_at, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

    ! Safety factor of 1.1, just in case
    ! Note: TB will return 0 for cutoff(am%minim_metapot), but TB does its own calc_connect, so doesn't matter
    if (1.1*max_atom_rij_change >= am%minim_at%cutoff - cutoff(am%minim_metapot)) then
      call print("gradient_func: Do calc_connect, atoms moved " // max_atom_rij_change // &
        "*1.1 >= buffer " // (am%minim_at%cutoff - cutoff(am%minim_metapot)), NERD)
      call calc_connect(am%minim_at)
      am%last_connect_x = x
    else
      call print("gradient_func: Do calc_dists, atoms moved " // max_atom_rij_change // &
        " *1.1 < buffer " // (am%minim_at%cutoff - cutoff(am%minim_metapot)), NERD)
      call calc_dists(am%minim_at)
    end if

    call print("energy_func using am%minim_at", NERD)
    if (current_verbosity() >= NERD) call print_xyz(am%minim_at,mainlog,real_format="f10.6")

    call calc(am%minim_metapot, am%minim_at, e = energy_func, args_str = am%minim_args_str)
    call print ("energy_func got energy " // energy_func, NERD)

    call fix_atoms_deform_grad(deform_grad, am%minim_at, am)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor, am%minim_at%travel, am%minim_at%lattice)

    am_data = transfer(am, am_data)

    call system_timer("energy_func")

  end function energy_func

  ! compute gradient (forces and virial)
  ! x is vectorized version of atomic positions
  ! result is vectorized version of forces
  function gradient_func(x, am_data)
    real(dp) :: x(:)
    character, optional :: am_data(:)
    real(dp) :: gradient_func(size(x))

    real(dp) :: deform_grad(3,3), virial(3,3), deform_grad_inv(3,3)
    real(dp), allocatable :: f(:,:)
    real(dp) :: max_atom_rij_change
    integer :: t_i(2), i
    integer, pointer, dimension(:) :: move_mask, fixed_metapot

    type(metapotential_minimise) :: am

    call system_timer("gradient_func")

    if (.not. present(am_data)) call system_abort("metapotential_minimise gradient_func must have am_data")
    am = transfer(am_data, am)
    ! beware of transfer and pointers !!!
    call atoms_repoint(am%minim_at)

    am%minim_n_eval_e = am%minim_n_eval_e + 1

    am%minim_n_eval_f = am%minim_n_eval_f + 1

    if (size(x) /= am%minim_at%N*3+9) call system_abort("Called gradient_func() with size mismatch " // &
      size(x) // " /= " // am%minim_at%N // "*3+9")

    ! Note: TB will return 0 for cutoff(am%minim_metapot), but TB does its own calc_connect, so doesn't matter
    max_atom_rij_change = max_rij_change(am%last_connect_x, x, cutoff(am%minim_metapot), &
      1.0_dp/am%pos_lat_preconditioner_factor)

    if (current_verbosity() >= NERD) call print("gradient_func got x " // x, NERD)

    call unpack_pos_dg(x, am%minim_at, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

    ! Safety factor of 1.1, just in case
    ! Note: TB will return 0 for cutoff(am%minim_metapot), but TB does its own calc_connect, so doesn't matter
    if (1.1*max_atom_rij_change >= am%minim_at%cutoff - cutoff(am%minim_metapot)) then
      call print("gradient_func: Do calc_connect, atoms moved " // max_atom_rij_change // &
        "*1.1 >= buffer " // (am%minim_at%cutoff - cutoff(am%minim_metapot)), NERD)
      call calc_connect(am%minim_at)
      am%last_connect_x = x
    else
      call print("gradient_func: Do calc_dists, atoms moved " // max_atom_rij_change // &
        " *1.1 < buffer " // (am%minim_at%cutoff - cutoff(am%minim_metapot)), NERD)
      call calc_dists(am%minim_at)
    end if

    call print("gradient_func using am%minim_at", NERD)
    if (current_verbosity() >= NERD) call print_xyz(am%minim_at,mainlog,real_format="f10.6")

    allocate(f(3,am%minim_at%N))
    f = 0.0_dp
    virial = 0.0_dp
    if (am%minim_do_pos .and. am%minim_do_lat) then
      call calc(am%minim_metapot, am%minim_at, f = f, virial = virial, args_str = am%minim_args_str)
    else if (am%minim_do_pos) then
      call calc(am%minim_metapot, am%minim_at, f = f, args_str = am%minim_args_str)
    else
      call calc(am%minim_metapot, am%minim_at, virial = virial, args_str = am%minim_args_str)
    endif

    ! zero forces if fixed by metapotential
    if (am%minim_do_pos .and. assign_pointer(am%minim_at, "fixed_metapot", fixed_metapot)) then
      do i=1, am%minim_at%N
	if (fixed_metapot(i) /= 0) f(:,i) = 0.0_dp
      end do
    endif

    ! Zero force on any fixed atoms
    if (assign_pointer(am%minim_at, 'move_mask', move_mask)) then
       do i=1,am%minim_at%N
          if (move_mask(i) == 0) f(:,i) = 0.0_dp
       end do
    end if

    if (current_verbosity() >= NERD) then
      call print ("gradient_func got f", NERD)
      call print(f, NERD)
      call print ("gradient_func got virial", NERD)
      call print(virial, NERD)
    endif

    f = transpose(deform_grad) .mult. f

    call inverse(deform_grad, deform_grad_inv)
    virial = virial .mult. transpose(deform_grad_inv)

    ! Zero out components corresponding to fixed lattice elements
    virial = merge(0.0_dp, virial, am%lattice_fix)

    call pack_pos_dg(-f, -virial, gradient_func, 1.0_dp/am%pos_lat_preconditioner_factor)
    if (current_verbosity() >= NERD) then
      call print ("gradient_func packed as", NERD)
      call print (gradient_func, NERD)
    endif

    t_i = maxloc(abs(f))
    call print ("gradient_func got max force component " // f(t_i(1),t_i(2)) // " on atom " // t_i(2) // &
      " max virial component " // maxval(abs(virial)))

    deallocate(f)

    call fix_atoms_deform_grad(deform_grad, am%minim_at,am)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor, am%minim_at%travel, am%minim_at%lattice)
    am_data = transfer(am, am_data)

    call system_timer("gradient_func")

  end function gradient_func

  ! compute energy and gradient (forces and virial)
  ! x is vectorized version of atomic positions
  ! result is vectorized version of forces
  subroutine both_func(x,val,grad,error,am_data)
    real(dp) :: x(:)
    real(dp) :: val
    real(dp) :: grad(:)
    character, optional :: am_data(:)
    integer :: error

    real(dp) :: deform_grad(3,3), virial(3,3), deform_grad_inv(3,3)
    real(dp), allocatable :: f(:,:)
    real(dp) :: max_atom_rij_change
    integer :: t_i(2), i
    integer, pointer, dimension(:) :: move_mask, fixed_metapot

    type(metapotential_minimise) :: am

    call system_timer("both_func")

    if (.not. present(am_data)) call system_abort("metapotential_minimise both_func must have am_data")
    am = transfer(am_data, am)
    ! beware of transfer and pointers !!!
    call atoms_repoint(am%minim_at)

    error = 0

    am%minim_n_eval_ef = am%minim_n_eval_ef + 1

    if (size(x) /= am%minim_at%N*3+9) call system_abort("Called both_func() with size mismatch " // &
      size(x) // " /= " // am%minim_at%N // "*3+9")

    ! Note: TB will return 0 for cutoff(am%minim_metapot), but TB does its own calc_connect, so doesn't matter
    max_atom_rij_change = max_rij_change(am%last_connect_x, x, cutoff(am%minim_metapot), &
      1.0_dp/am%pos_lat_preconditioner_factor)
max_atom_rij_change = 1.038_dp

    call print("both_func got x " // x, NERD)

    call unpack_pos_dg(x, am%minim_at, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

    ! Safety factor of 1.1, just in case
    ! Note: TB will return 0 for cutoff(am%minim_metapot), but TB does its own calc_connect, so doesn't matter
    if (1.1*max_atom_rij_change >= am%minim_at%cutoff - cutoff(am%minim_metapot)) then
      call print("both_func: Do calc_connect, atoms moved " // max_atom_rij_change // "*1.1 >= buffer " // &
        (am%minim_at%cutoff - cutoff(am%minim_metapot)), NERD)
      call calc_connect(am%minim_at)
      am%last_connect_x = x
    else
      call print("both_func: Do calc_dists, atoms moved " // max_atom_rij_change // " *1.1 < buffer " // &
        (am%minim_at%cutoff - cutoff(am%minim_metapot)), NERD)
      call calc_dists(am%minim_at)
    end if

    call print("both_func using am%minim_at", NERD)
    if (current_verbosity() >= NERD) call print_xyz(am%minim_at,mainlog,real_format="f10.6")

    allocate(f(3,am%minim_at%N))
    f = 0.0_dp
    virial = 0.0_dp
    if (am%minim_do_pos .and. am%minim_do_lat) then
      call calc(am%minim_metapot, am%minim_at, e = val, f = f, virial = virial, args_str = am%minim_args_str)
    else if (am%minim_do_pos) then
      call calc(am%minim_metapot, am%minim_at, e = val, f = f, args_str = am%minim_args_str)
    else
      call calc(am%minim_metapot, am%minim_at, e = val, virial = virial, args_str = am%minim_args_str)
    endif
    call print ("both_func got f", NERD)
    call print(f, NERD)
    call print ("both_func got virial", NERD)
    call print(virial, NERD)

    ! zero forces if fixed by metapotential
    if (am%minim_do_pos .and. assign_pointer(am%minim_at, "fixed_metapot", fixed_metapot)) then
      do i=1, am%minim_at%N
	if (fixed_metapot(i) /= 0) f(:,i) = 0.0_dp
      end do
    endif

    ! Zero force on any fixed atoms
    if (assign_pointer(am%minim_at, 'move_mask', move_mask)) then
       do i=1,am%minim_at%N
          if (move_mask(i) == 0) f(:,i) = 0.0_dp
       end do
    end if

    f = transpose(deform_grad) .mult. f

    call inverse(deform_grad, deform_grad_inv)
    virial = virial .mult. transpose(deform_grad_inv)

    call pack_pos_dg(-f, -virial, grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call print ("both_func gradient packed as", NERD)
    call print (grad, NERD)

    t_i = maxloc(abs(f))
    call print ("both_func got max force component " // f(t_i(1),t_i(2)) // " on atom " // t_i(2) // &
      " max virial component " // maxval(abs(virial)))

    deallocate(f)

    call fix_atoms_deform_grad(deform_grad, am%minim_at,am)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor, am%minim_at%travel, am%minim_at%lattice)

    am_data = transfer(am, am_data)

    call system_timer("both_func")

  end subroutine both_func

  function max_rij_change(last_connect_x, x, r_cut, lat_factor)
    real(dp), intent(inout) :: last_connect_x(:)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: r_cut, lat_factor
    real(dp) :: max_rij_change

    real(dp) :: F0(3,3), F0_evals(3), dF(3,3), dF_evals(3)

    call print("max_rij_change, size(last_connect_x) " // size(last_connect_x) // &
      " size(x) " // size(x), NERD)
    if (size(last_connect_x) == size(x)) then
      F0 = lat_factor*reshape(last_connect_x(1:9), (/ 3,3 /))
      dF = lat_factor*reshape(x(1:9), (/ 3,3 /)) - F0
      call diagonalise(matmul(F0,transpose(F0)), F0_evals)
      call diagonalise(matmul(dF,transpose(dF)), dF_evals)
      F0_evals = sqrt(F0_evals)
      dF_evals = sqrt(dF_evals)
      call print("max_rij_change = " // maxval(abs(F0_evals)) //"*2.0_dp*"// &
        maxval(abs(last_connect_x-x))*sqrt(3.0_dp)  // &
        " + " // maxval(abs(dF_evals))//"*"//r_cut, VERBOSE)
      max_rij_change = maxval(abs(F0_evals))*2.0_dp*maxval(abs(last_connect_x-x))*sqrt(3.0_dp) + &
                 maxval(abs(dF_evals))*r_cut
    else
      call system_abort("max_rij_change: size(last_connect_x)="//size(last_connect_x)// &
        " /= size(x)="//size(x))
    endif

  end function max_rij_change

  ! prepare atoms structure given a deformation gradient matrix
  subroutine prep_atoms_deform_grad(deform_grad, at, am)
    real(dp), intent(in) :: deform_grad(3,3)
    type(Atoms), intent(inout) :: at
    type(metapotential_minimise), intent(inout) :: am

    real(dp) :: new_lat(3,3)

    am%minim_save_lat = at%lattice

    at%pos = deform_grad .mult. at%pos
    new_lat = deform_grad .mult. am%minim_save_lat
    call set_lattice(at, new_lat)
  end subroutine prep_atoms_deform_grad

  ! remove effect of deformation gradient applied by prep_atoms_deform_grad
  subroutine fix_atoms_deform_grad(deform_grad, at,am)
    real(dp), intent(in) :: deform_grad(3,3)
    type(Atoms), intent(inout) :: at
    type(metapotential_minimise), intent(inout) :: am

    real(dp) :: deform_grad_inv(3,3)

    call inverse(deform_grad, deform_grad_inv)
    at%pos = deform_grad_inv .mult. at%pos
    call set_lattice(at, am%minim_save_lat)

  end subroutine fix_atoms_deform_grad

  ! unpack a deformation grad and a pos array from 1-D array into a atoms%pos 2-D array
  subroutine unpack_pos_dg(xx, at, dg, lat_factor)
    real(dp), intent(in) :: xx(:)
    type(Atoms), intent(inout) :: at
    real(dp) :: dg(3,3)
    real(dp), intent(in) :: lat_factor

    if (3*at%N+9 /= size(xx)) call system_abort("Called unpack_pos with mismatching sizes xx " // size(xx) // " at " // at%N)

    dg = lat_factor*reshape(xx(1:9), (/ 3,3 /) )
    at%pos = reshape(xx(10:), (/ 3, at%N /) )

    at%pos = at%pos - (at%lattice .mult. at%travel)
  end subroutine

  ! pack a 3xN 2-D array into a 1-D array
  subroutine pack_pos_dg(x2d, dg2d, x, lat_factor, travel, lattice)
    real(dp), intent(in) :: x2d(:,:)
    real(dp), intent(in)  :: dg2d(3,3)
    real(dp), intent(out) :: x(:)
    real(dp), intent(in)  :: lat_factor
    integer, intent(in), optional :: travel(:,:)
    real(dp), intent(in), optional :: lattice(3,3)

    if ((present(travel) .and. .not. present(lattice)) .or. &
	(.not. present(travel) .and. present(lattice))) then
      call system_abort("Called pack with travel and no lattice or vice versa")
    endif

    if (size(x2d,1) /= 3) call system_abort("Called pack with mismatching size(x2d,1) " // &
      size(x2d,1) // " != 3")
    
    if (size(dg2d,1) /= 3 .or. size(dg2d,2) /= 3) &
      call system_abort("Called pack with mismatching size(dg2d,1) " // shape(dg2d) // " != 3x3")
    if (size(x) /= (size(x2d)+size(dg2d))) &
      call system_abort("Called pack with mismatching size x " // size(x) // " x2d " // size(x2d) // &
      " dg2d " // size(dg2d))

    x(1:9) = lat_factor*reshape(dg2d, (/ 9 /) )

    if (present(travel)) then
      x(10:) = reshape(x2d + (lattice .mult. travel), (/ size(x2d) /))
    else
      x(10:) = reshape(x2d, (/ size(x2d) /) )
    endif

  end subroutine

  !*************************************************************************
  !*
  !*  MetaPotential_Sum routines
  !*
  !*************************************************************************

  subroutine MetaPotential_Sum_Initialise(this, args_str, pot1, pot2, mpi)
    type(MetaPotential_Sum), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    type(Potential), intent(in), target :: pot1, pot2
    type(MPI_Context), intent(in), optional :: mpi

    type(Dictionary) :: params

    call finalise(this)

    call initialise(params)
    call param_register(params, 'subtract_pot1', 'F', this%subtract_pot1)
    call param_register(params, 'subtract_pot2', 'F', this%subtract_pot2)

    if (.not. param_read_line(params, args_str, ignore_unknown=.true.)) then
       call system_abort('MetaPotential_Sum_Initialise failed to parse args_str="'//trim(args_str)//'"')
    end if

    call finalise(params)

    this%pot1 => pot1
    this%pot2 => pot2

    if (present(mpi)) this%mpi = mpi

  end subroutine MetaPotential_Sum_Initialise

  subroutine MetaPotential_Sum_Finalise(this)
    type(MetaPotential_Sum), intent(inout) :: this
    
    nullify(this%pot1)
    nullify(this%pot2)

    this%subtract_pot1 = .false.
    this%subtract_pot2 = .false.

  end subroutine MetaPotential_Sum_Finalise

  subroutine MetaPotential_Sum_Print(this, file)
    type(MetaPotential_Sum), intent(inout) :: this
    type(Inoutput), intent(inout), optional :: file

    call print('MetaPotential_Sum:', file=file)
    call print(' subtract_pot1=' // this%subtract_pot1, file=file)
    call print(' subtract_pot2=' // this%subtract_pot2, file=file)
    call print('', file=file)
    if (associated(this%pot1)) then
       call print('Potential 1:', file=file)
       call print(this%pot1, file=file)
       call print('', file=file)
    else
       call print('Potential 1 not initialised', file=file)
       call print('', file=file)
    end if
    if (associated(this%pot2)) then
       call print('Potential 2:', file=file)
       call Print(this%pot2, file=file)
       call print('', file=file)
    else
       call print('Potential 2 not initialised', file=file)
       call print('', file=file)
    end if

  end subroutine MetaPotential_Sum_Print

  subroutine MetaPotential_Sum_Calc(this, at, e, local_e, f, df, virial, args_str, err)
    type(MetaPotential_Sum), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    real(dp), intent(out), optional :: e
    real(dp), intent(out), optional :: local_e(:)
    real(dp), intent(out), optional :: f(:,:)
    real(dp), intent(out), optional :: df(:,:)
    real(dp), intent(out), optional :: virial(3,3)
    character(*), intent(in), optional :: args_str
    integer, intent(out), optional :: err

    real(dp) :: my_e_1, my_e_2
    real(dp), allocatable :: my_local_e_1(:), my_local_e_2(:)
    real(dp), allocatable :: my_f_1(:,:), my_f_2(:,:)
    real(dp), allocatable :: my_df_1(:,:), my_df_2(:,:)
    real(dp) :: my_virial_1(3,3), my_virial_2(3,3)
    integer :: my_err_1, my_err_2

    if (present(local_e)) allocate(my_local_e_1(size(local_e)), my_local_e_2(size(local_e)))
    if (present(f)) allocate(my_f_1(size(f, 1),size(f, 2)), my_f_2(size(f, 1),size(f, 2)))
    if (present(df)) allocate(my_df_1(size(df, 1),size(df, 2)), my_df_2(size(df, 1),size(df, 2)))

    call calc(this%pot1, at, e, local_e, f, df, virial, args_str, err)

    if (present(err)) then
      if (err /= 0) then
        return
      end if
    end if

    if (present(e)) my_e_1 = e
    if (present(local_e)) my_local_e_1 = local_e
    if (present(f)) my_f_1 = f
    if (present(df)) my_df_1 = df
    if (present(virial)) my_virial_1 = virial
    if (present(err)) my_err_1 = err

    if (this%subtract_pot1) then
      if (present(e)) my_e_1 = - my_e_1
      if (present(local_e)) my_local_e_1 = - my_local_e_1
      if (present(f)) my_f_1 = - my_f_1
      if (present(df)) my_df_1 = - my_df_1
      if (present(virial)) my_virial_1 = - my_virial_1
    end if

    call calc(this%pot2, at, e, local_e, f, df, virial, args_str, err)

    if (present(err)) then
      if (err /= 0) then
        return
      end if
    end if

    if (present(e)) my_e_2 = e
    if (present(local_e)) my_local_e_2 = local_e
    if (present(f)) my_f_2 = f
    if (present(df)) my_df_2 = df
    if (present(virial)) my_virial_2 = virial
    if (present(err)) my_err_2 = err

    if (this%subtract_pot2) then
      if (present(e)) my_e_2 = - my_e_2
      if (present(local_e)) my_local_e_2 = - my_local_e_2
      if (present(f)) my_f_2 = - my_f_2
      if (present(df)) my_df_2 = - my_df_2
      if (present(virial)) my_virial_2 = - my_virial_2
    end if

    if (present(e)) e = my_e_1 + my_e_2
    if (present(local_e)) local_e = my_local_e_1 + my_local_e_2
    if (present(f)) f = my_f_1 + my_f_2
    if (present(df)) df = my_df_1 + my_df_2
    if (present(virial)) virial = my_virial_1 + my_virial_2

    if (present(local_e)) deallocate(my_local_e_1, my_local_e_2)
    if (present(f)) deallocate(my_f_1, my_f_2)
    if (present(df)) deallocate(my_df_1, my_df_2)
  
  end subroutine MetaPotential_Sum_Calc

  function MetaPotential_Sum_Cutoff(this)
    type(MetaPotential_Sum), intent(in) :: this
    real(dp) :: metapotential_sum_cutoff

    if(associated(this%pot1) .and. associated(this%pot2)) then
       metapotential_sum_cutoff = max(cutoff(this%pot1), cutoff(this%pot2))
    else
       metapotential_sum_cutoff = 0.0_dp
    endif

  end function MetaPotential_Sum_Cutoff


#include "MetaPotential_ForceMixing_routines.f95"

#ifdef HAVE_LOCAL_E_MIX
#include "MetaPotential_Local_E_Mix_routines.f95"
#endif
#ifdef HAVE_ONIOM
#include "MetaPotential_ONIOM_routines.f95"
#endif

#include "MetaPotential_Hybrid_utils.f95"

end module MetaPotential_module
