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

module DomainDecomposition_module
  use atoms_module
  use error_module
  use linearalgebra_module
  use MPI_context_module
  use table_module

  implicit none

  character(TABLE_STRING_LENGTH), parameter :: Z_STR = "z"
  character(TABLE_STRING_LENGTH), parameter :: MASS_STR = "mass"
  character(TABLE_STRING_LENGTH), parameter :: SPECIES_STR = "species"
  character(TABLE_STRING_LENGTH), parameter :: TRAVEL_STR = "travel"
  character(TABLE_STRING_LENGTH), parameter :: POS_STR = "pos"
  character(TABLE_STRING_LENGTH), parameter :: VELO_STR = "velo"
  character(TABLE_STRING_LENGTH), parameter :: LOCAL_TO_GLOBAL_STR = &
       "local_to_global"
  character(TABLE_STRING_LENGTH), parameter :: GHOSTS_R_STR = "ghosts_r"
  character(TABLE_STRING_LENGTH), parameter :: GHOSTS_L_STR = "ghosts_l"

  type DomainDecomposition

     integer                    :: Ntotal               !% Number of total particles in this simulation
     integer, pointer           :: local_to_global(:)   !% Local index to global index
     integer, allocatable       :: global_to_local(:)   !% Global index to local index

     integer                    :: decomposition(3)    = (/ 2, 2, 2 /)   !% Type of decomposition

     real(DP)                   :: requested_border    = 0.0_DP
     real(DP)                   :: border(3)           = 0.0_DP
     real(DP)                   :: verlet_shell        = 0.0_DP

     logical                    :: communicate_forces  = .false.

     real(DP)                   :: lower(3)               !% Lower domain boundary, in fraction of the total cell
     real(DP)                   :: upper(3)               !% Upper domain boundary, in fraction of the total cell
     real(DP)                   :: lower_with_border(3)   !% Lower domain boundary, including border
     real(DP)                   :: upper_with_border(3)   !% Upper domain boundary, including border

     type(MPI_context)          :: mpi                  !% MPI communicator
     logical                    :: periodic(3)          !% Periodicity for domain decomposition
     integer                    :: l(3), r(3)           !% Ranks of neighboring domains in x-, y- and z-direction
     real(DP)                   :: off_l(3), off_r(3)   !% Distance vector to neighboring domaind

     type(Table)                :: state_properties   !% Fields to communicate if for particles
     type(Table)                :: ghost_properties   !% Fields to communicate for ghosts
     type(Table)                :: backc_properties   !% Fields for back-communication after force computations     

     logical, allocatable       :: state_mask(:)
     logical, allocatable       :: ghost_mask(:)
     logical, allocatable       :: backc_mask(:)

     integer                    :: state_buffer_size
     integer                    :: ghost_buffer_size
     integer                    :: backc_buffer_size

     character(1), allocatable  :: send_l(:)   !% buffer for sending to the left
     character(1), allocatable  :: send_r(:)   !% buffer for sending to the right
     character(1), allocatable  :: recv_l(:)   !% buffer for receiving from the left
     character(1), allocatable  :: recv_r(:)   !% buffer for receiving from the right

     integer                    :: n_ghosts_r(3), n_ghosts_l(3)  !% length of the ghost particle lists

     integer, pointer           :: ghosts_r(:)  => NULL()   !% particles send to the right (where they become ghosts)
     integer, pointer           :: ghosts_l(:)  => NULL()   !% particles send to the left (where they become ghosts)

     integer                    :: n_send_p_tot   !% Statistics: Number of total particles send
     integer                    :: n_recv_p_tot   !% Statistics: Number of total particles received
     integer                    :: n_send_g_tot   !% Statistics: Number of total ghosts send
     integer                    :: n_recv_g_tot   !% Statistics: Number of total ghosts received
     integer                    :: nit_p          !% Statistics: Number of particle send events
     integer                    :: nit_g          !% Statistics: Number of ghost send events
     
  endtype DomainDecomposition

  interface initialise
     module procedure domain_decomposition_initialise
  endinterface

  interface finalise
     module procedure domain_decomposition_finalise
  endinterface

  interface allocate
     module procedure domain_decomposition_allocate
  endinterface

  interface enable
     module procedure domain_decomposition_enable
  endinterface enable

  interface is_in_domain
     module procedure domain_decomposition_is_in_domain
  endinterface is_in_domain

  interface set_border
     module procedure domain_decomposition_set_border
  endinterface

  interface communicate_particles
     module procedure domain_decomposition_communicate_particles
  endinterface

  interface communicate_ghosts
     module procedure domain_decomposition_communicate_ghosts
  endinterface

  interface communicate_forces
     module procedure domain_decomposition_communicate_forces
  endinterface


  !
  ! Supplement to the Dictionary class required by the domain
  ! decomposition module.
  ! Move it do Dictionary.f95?
  !

  interface keys_to_mask
     module procedure dictionary_keys_to_mask
  endinterface keys_to_mask

  interface pack_buffer
     module procedure dictionary_pack_buffer
  endinterface

  interface unpack_buffer
     module procedure dictionary_unpack_buffer
  endinterface

contains

  !% Initialize the domain decomposition module
  subroutine domain_decomposition_initialise(this, at, mpi, &
       decomposition, verlet_shell, error)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    type(MPI_context), intent(in)             :: mpi
    integer, intent(in), optional             :: decomposition(3)
    real(DP), intent(in), optional            :: verlet_shell
    integer, intent(out), optional :: error

    ! ---

    integer   :: i, d
    real(DP)  :: l(3)
    logical   :: periods_for_mpi(3)

    ! ---

    INIT_ERROR(error)

    call print("- domain_decomposition_initialise -", PRINT_VERBOSE)

    if (present(decomposition)) then
       this%decomposition  = decomposition
    endif

    if (present(verlet_shell)) then
       this%verlet_shell   = verlet_shell
    endif

    call print("Parallelization using domain decomposition (via MPI) over " // &
         this%decomposition // "domains.")

    call initialise(this%mpi, &
         context  = mpi, &
         dims     = this%decomposition, &
         error   = error)
    PASS_ERROR(error)

    if (this%decomposition(1)*this%decomposition(2)*this%decomposition(3) /= &
         this%mpi%n_procs) then
       RAISE_ERROR("Decomposition geometry requires " // this%decomposition(1)*this%decomposition(2)*this%decomposition(3) // " processes, however, MPI returns " // this%mpi%n_procs // " processes.", error)
    endif

! For now this needs to be true
!    this%periodic          = at%periodic
    this%periodic          = .true.
    this%requested_border  = 0.0_DP

    call print("coords             = ( " // this%mpi%my_coords // ")", PRINT_VERBOSE)

    do d = 1, 3
       call cart_shift( &
            this%mpi, d-1, 1, this%l(d), this%r(d), error)
       PASS_ERROR(error)
    enddo

    if ((at%lattice(:, 1) .dot. at%lattice(:, 2)) > 1e-6 .or. &
        (at%lattice(:, 2) .dot. at%lattice(:, 3)) > 1e-6 .or. &
        (at%lattice(:, 3) .dot. at%lattice(:, 1)) > 1e-6) then
       RAISE_ERROR("Only orthorhombic lattices are support for now.", error)
    endif

    l        = (/ at%lattice(1, 1), at%lattice(2, 2), at%lattice(3, 3) /)

    this%lower  = (1.0_DP*this%mpi%my_coords)     / this%decomposition
    this%upper  = (1.0_DP*(this%mpi%my_coords+1)) / this%decomposition

    call print("lower              = ( " // this%lower // " )", PRINT_VERBOSE)
    call print("upper              = ( " // this%upper // " )", PRINT_VERBOSE)

    this%off_r  = 0.0_DP
    this%off_l  = 0.0_DP

    do d = 1, 3
!       if (at%periodic(d) .and. this%decomposition(d) > 1) then
       if (this%decomposition(d) > 1) then
          ! No periodicity in binning because we explicitly copy the atoms
          ! from the other processors
!          this%locally_periodic(d) = .false.
       else
          this%periodic(d) = .false.
       endif

       if (this%mpi%my_coords(d) == 0) then
          this%off_l(d) = -l(d)
       else if (this%mpi%my_coords(d) == this%decomposition(d)-1) then
          this%off_r(d) = l(d)
       endif
    enddo

!    call print("periodic (global)  = ( " // at%periodic // " )", PRINT_VERBOSE)
    call print("periodic (par.)    = ( " // this%periodic // " )", PRINT_VERBOSE)
!    call print("periodic (local)   = ( " // this%locally_periodic // " )", PRINT_VERBOSE)

    call print("off_l              = ( " // this%off_l // " )", PRINT_VERBOSE)
    call print("off_r              = ( " // this%off_r // " )", PRINT_VERBOSE)

    this%n_send_p_tot  = 0
    this%n_recv_p_tot  = 0
    this%n_send_g_tot  = 0
    this%n_recv_g_tot  = 0
    this%nit_p         = 0
    this%nit_g         = 0


    ! Additional properties

    call add_property(at, LOCAL_TO_GLOBAL_STR, 0)
    call add_property(at, GHOSTS_R_STR, 0)
    call add_property(at, GHOSTS_L_STR, 0)

    call append(this%state_properties, Z_STR)
    call append(this%state_properties, SPECIES_STR)
    call append(this%state_properties, MASS_STR)
    call append(this%state_properties, TRAVEL_STR)
    call append(this%state_properties, POS_STR)
    call append(this%state_properties, VELO_STR)
    call append(this%state_properties, LOCAL_TO_GLOBAL_STR)

    call append(this%ghost_properties, Z_STR)
    call append(this%ghost_properties, SPECIES_STR)
    call append(this%ghost_properties, MASS_STR)
    call append(this%ghost_properties, TRAVEL_STR)
    call append(this%ghost_properties, POS_STR)
    call append(this%ghost_properties, VELO_STR)
    call append(this%ghost_properties, LOCAL_TO_GLOBAL_STR)

  endsubroutine domain_decomposition_initialise


  !% Allocate buffers
  subroutine domain_decomposition_allocate(this, at, error)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    integer, intent(out), optional :: error

    ! ---

    integer  :: i, s
    integer, pointer :: ptr(:)

    ! ---

    INIT_ERROR(error)

    call print("- domain_decomposition_allocate -", PRINT_VERBOSE)

    allocate(this%global_to_local(at%N))

!    call log_memory_start("domain_decomposition_allocate")

    this%state_buffer_size = 0
    this%ghost_buffer_size = 0
    this%backc_buffer_size = 0

    if (this%state_properties%N > 0) then
       allocate(this%state_mask(at%properties%N))
       call keys_to_mask( &
            at%properties, this%state_properties, this%state_mask, &
            s = this%state_buffer_size, error = error)
       PASS_ERROR(error)
    endif
    if (this%ghost_properties%N > 0) then
       allocate(this%ghost_mask(at%properties%N))
       call keys_to_mask( &
            at%properties, this%ghost_properties, this%ghost_mask, &
            s = this%ghost_buffer_size, error = error)
       PASS_ERROR(error)
    endif
    if (this%backc_properties%N > 0) then
       allocate(this%backc_mask(at%properties%N))
       call keys_to_mask( &
            at%properties, this%backc_properties, this%backc_mask, &
            s = this%backc_buffer_size, error = error)
       PASS_ERROR(error)
    endif

    call print("size(state_mask) = " // size(this%state_mask))
    call print("state_mask = " // this%state_mask)

    ! Additional rank information
    !this%ghost_buffer_size  = this%ghost_buffer_size + 1

    call print("     state_buffer_size  = " // this%state_buffer_size, PRINT_VERBOSE)
    call print("     ghost_buffer_size  = " // this%ghost_buffer_size, PRINT_VERBOSE)
    call print("     backc_buffer_size  = " // this%backc_buffer_size, PRINT_VERBOSE)

    s = max(this%state_buffer_size, this%ghost_buffer_size) * at%N

    allocate(this%send_l(s))
    allocate(this%send_r(s))
    allocate(this%recv_l(s))
    allocate(this%recv_r(s))

    if (.not. assign_pointer(at, LOCAL_TO_GLOBAL_STR, this%local_to_global)) then
       RAISE_ERROR("Could not assign pointer for property '" // LOCAL_TO_GLOBAL_STR // "'.", error)
    endif
    if (.not. assign_pointer(at, GHOSTS_R_STR, this%ghosts_r)) then
       RAISE_ERROR("Could not assign pointer for property '" // GHOSTS_R_STR // "'.", error)
    endif
    if (.not. assign_pointer(at, GHOSTS_L_STR, this%ghosts_l)) then
       RAISE_ERROR("Could not assign pointer for property '" // GHOSTS_L_STR // "'.", error)
    endif

!    call log_memory_estimate(this%send_l)
!    call log_memory_estimate(this%send_r)
!    call log_memory_estimate(this%recv_l)
!    call log_memory_estimate(this%recv_r)
!
!    call log_memory_estimate(this%ghosts_r)
!    call log_memory_estimate(this%ghosts_l)
!
!    call log_memory_stop("domain_decomposition_allocate")

  endsubroutine domain_decomposition_allocate


  !% Enable domain decomposition, after this call every process
  !% remains only with the atoms in its local domain.
  subroutine domain_decomposition_enable(this, at, error)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    integer, intent(out), optional :: error

    ! ---

    integer :: i, j

    ! ---

    INIT_ERROR(error)

    j = 0
    do i = 1, at%N
       !% Assign global indices
       this%local_to_global(i) = i
       this%global_to_local(i) = -1

       !% If this is in the domain, keep it
       if (is_in_domain(this, at, at%pos(:, i))) then
          j = j+1
          if (i /= j) then
             call copy_entry(at, i, j)
          endif
          this%global_to_local(i) = j
       endif
    enddo
    at%N = j
    at%Ndomain = j

    this%Ntotal = sum(this%mpi, at%N, error)
    PASS_ERROR(error)

    !
    ! Copy ghost particles (for first integration step)
    !

    call communicate_ghosts(this, at, .true.)
    if (this%communicate_forces) then
       call communicate_forces(this, at)
    endif

  endsubroutine domain_decomposition_enable


  !% Is the position in the current domain?
  function domain_decomposition_is_in_domain(this, at, r) result(id)
    implicit none

    type(DomainDecomposition), intent(in)  :: this
    type(Atoms), intent(in)                :: at
    real(DP), intent(in)                   :: r(3)

    logical                                :: id

    ! ---

    real(DP)  :: s(3)

    ! ---

    s = at%g .mult. r
    id = all(s >= this%lower) .and. all(s < this%upper)

  endfunction domain_decomposition_is_in_domain


  !% Set the communication border
  subroutine domain_decomposition_set_border(this, at, &
       border, verlet_shell, error)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    real(DP), intent(in)                      :: border
    real(DP), intent(in), optional            :: verlet_shell
    integer, intent(out), optional :: error

    ! ---

    integer  :: d

    ! ---

    INIT_ERROR(error)

    call print("- domain_decomposition_set_border -", PRINT_VERBOSE)

    if (present(verlet_shell)) then
       this%verlet_shell  = verlet_shell
    endif

    this%requested_border  = max(this%requested_border, border)
    this%border            = this%requested_border + this%verlet_shell

    call print("requested_border  = " // this%requested_border, PRINT_VERBOSE)
    call print("verlet_shell      = " // this%verlet_shell, PRINT_VERBOSE)
    call print("border            = " // this%border, PRINT_VERBOSE)

    if (any((at%lattice .mult. (this%upper - this%lower)) < 2*this%border)) then
       RAISE_ERROR("Domain smaller than twice the border. This does not work (yet).", error)
    endif

    do d = 1, 3
       if (this%periodic(d) .or. this%mpi%my_coords(d) /= 0) then
          this%lower_with_border(d)  = this%lower(d) - (at%g(d, :) .dot. this%border)
       else
          this%lower_with_border(d)  = this%lower(d)
       endif
       if (this%periodic(d) .or. this%mpi%my_coords(d) /= this%decomposition(d)-1) then
          this%upper_with_border(d)  = this%upper(d) + (at%g(d, :) .dot. this%border)
       else
          this%upper_with_border(d)  = this%upper(d)
       endif
    enddo

    call print("lower_with_border  = ( " // this%lower_with_border // " )", PRINT_VERBOSE)
    call print("upper_with_border  = ( " // this%upper_with_border // " )", PRINT_VERBOSE)

  endsubroutine domain_decomposition_set_border

  
  !% Free memory, clean up
  subroutine domain_decomposition_finalise(this)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this

    ! ---

    call print("- domain_decomposition_finalise -", PRINT_VERBOSE)

    deallocate(this%global_to_local)

    deallocate(this%send_l)
    deallocate(this%send_r)
    deallocate(this%recv_l)
    deallocate(this%recv_r)

    deallocate(this%ghosts_l)
    deallocate(this%ghosts_r)

    call print("Average number of particles sent/received per iteration:", PRINT_VERBOSE)
    call print("Particles send  = " // (1.0_DP*this%n_send_p_tot)/this%nit_p, PRINT_VERBOSE)
    call print("Particles recv  = " // (1.0_DP*this%n_recv_p_tot)/this%nit_p, PRINT_VERBOSE)
    call print("Ghosts send     = " // (1.0_DP*this%n_send_g_tot)/this%nit_g, PRINT_VERBOSE)
    call print("Ghosts recv     = " // (1.0_DP*this%n_recv_g_tot)/this%nit_g, PRINT_VERBOSE)

    call finalise(this%state_properties)
    call finalise(this%ghost_properties)
    call finalise(this%backc_properties)

    call finalise(this%mpi)

  endsubroutine domain_decomposition_finalise


  !% Add particle data to the send buffer
  subroutine copy_to_send_buffer(this, at, i, n, buffer)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    integer, intent(in)                       :: i
    integer, intent(inout)                    :: n
    character(1), intent(inout)               :: buffer(:)

    ! ---

! WAIT FOR DICTIONARY IN ATOMS
!    call pack_buffer(at%data, mask, i, n, buffer)

    this%global_to_local(this%local_to_global(i)) = 0 ! This one is gone

!    write (ilog, '(5X,A,3I5,6F20.10)')  "send: ", at%index(i), i, at%Z(i), POS3(p, i), PNC3(p, i)

  endsubroutine copy_to_send_buffer


  !% Copy particle data from the receive buffer
  subroutine copy_from_recv_buffer(this, at, n, buffer, off)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    integer, intent(in)                       :: n
    character(1), intent(in)                  :: buffer(:)
    real(DP), intent(in)                      :: off(3)

    ! ---

    integer  :: i

    ! ---

!    do i = 1, n
    i = 0
    do while (i < n)
       at%Ndomain = at%Ndomain+1

! WAIT FOR DICTIONARY IN ATOMS
!       call unpack_buffer(at%data, mask, i, buffer, at%Ndomain)

       this%global_to_local(this%local_to_global(at%Ndomain))  = at%Ndomain

!       write (ilog, '(5X,A,3I5,9F20.10)')  "recv1: ", at%index(at%Ndomain), this%global_to_local(at%index(at%Ndomain)), at%Z(at%Ndomain), POS3(p, at%Ndomain), PNC3(p, at%Ndomain), off(:)

!       POS3(p, at%Ndomain)  = POS3(p, at%Ndomain) + off(:)
!       PNC3(at, at%Ndomain)  = PNC3(at, at%Ndomain) + off

!       write (ilog, '(5X,A,3I5,9F20.10)')  "recv2: ", at%index(at%Ndomain), this%global_to_local(at%index(at%Ndomain)), at%Z(at%Ndomain), POS3(p, at%Ndomain), PNC3(p, at%Ndomain), off(:)
    enddo

  endsubroutine copy_from_recv_buffer


  !% Communicate particles which left the domains
  !% to the neighboring domains (former order routine)
  subroutine domain_decomposition_communicate_particles(this, at, error)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    integer, intent(out), optional :: error

    ! ---

    ! 
    ! General and auxiliary variables      
    !

    integer   :: i, d
    integer   :: oldnatloc

    ! 
    ! Structure variables, mpi and system structure
    !

    integer   :: n_send_l, n_send_r, n_recv_l, n_recv_r

    real(DP)  :: off_l(3), off_r(3)

    ! ---

    INIT_ERROR(error)

    call system_timer("domain_decomposition_communicate_particles")

    this%nit_p = this%nit_p + 1

    do i = at%Ndomain+1, at%N
       this%global_to_local(this%local_to_global(i)) = 0
    enddo

    !
    ! Loop over dimensions and distribute particle in the
    ! respective direction
    !

    do d = 1, 3

       if (this%decomposition(d) > 1) then

          oldnatloc  = at%Ndomain
          at%Ndomain   = 0

          n_send_r   = 0
          n_send_l   = 0

          do i = 1, oldnatloc

             if ((at%g(d, :) .dot. at%pos(:, i)) >= this%upper(d)) then
                ! Send to the right

                call copy_to_send_buffer(this, at, i, n_send_r, this%send_r)
!                n_send_r = n_send_r + 1

             else if ((at%g(d, :) .dot. at%pos(:, i)) < this%lower(d)) then
                ! Send to the left

                call copy_to_send_buffer(this, at, i, n_send_l, this%send_l)
!                n_send_l = n_send_l + 1

             else
                ! Keep on this processor and reorder

                at%Ndomain = at%Ndomain+1

                if (at%Ndomain /= i) then
                   call copy_entry(at, at%Ndomain, i)
                endif

             endif

          enddo

          !write (ilog, *)  d, "r: ", n_send_r
          !write (ilog, *)  d, "l: ", n_send_l

          this%n_send_p_tot = this%n_send_p_tot + n_send_r + n_send_l

!         call mpi_sendrecv( &
!               this%send_r, n_send_r, MPI_DOUBLE_PRECISION, this%r(d), 0, &
!               this%recv_l, this%state_buffer_size*at%maxnatloc, MPI_DOUBLE_PRECISION, this%l(d), 0, &
!               this%comm, status, i)
!          PASS_MPI_ERROR(i, error)
!
!          call mpi_get_count(status, MPI_DOUBLE_PRECISION, n_recv_l, i)
!          PASS_MPI_ERROR(i, error)

          call sendrecv(this%mpi, &
               this%send_r(1:n_send_r), this%r(d), 0, &
               this%recv_l, this%l(d), 0, &
               n_recv_l, error)
          PASS_ERROR(error)

!          call mpi_sendrecv( &
!               this%send_l, n_send_l, MPI_DOUBLE_PRECISION, this%l(d), 1, &
!               this%recv_r, this%state_buffer_size*at%maxnatloc, MPI_DOUBLE_PRECISION, this%r(d), 1, &
!               this%comm, status, i)
!          PASS_MPI_ERROR(i, error)
!
!          call mpi_get_count(status, MPI_DOUBLE_PRECISION, n_recv_r, i)
!          PASS_MPI_ERROR(i, error)

          call sendrecv(this%mpi, &
               this%send_l(1:n_send_l), this%l(d), 1, &
               this%recv_r, this%r(d), 1, &
               n_recv_r, error)
          PASS_ERROR(error)

          this%n_recv_p_tot = this%n_recv_p_tot + n_recv_r/this%state_buffer_size + n_recv_l/this%state_buffer_size

          off_l    = 0.0_DP
          ! This will be done by inbox
          off_l(d) = this%off_l(d)

          off_r   = 0.0_DP
          ! This will be done by inbox
          off_r(d) = this%off_r(d)

          call copy_from_recv_buffer(this, at, n_recv_l, this%recv_l, off_l)
          call copy_from_recv_buffer(this, at, n_recv_r, this%recv_r, off_r)

       endif

    enddo

    at%N = at%Ndomain

    call system_timer("domain_decomposition_communicate_particles")

  endsubroutine domain_decomposition_communicate_particles


  !% Copy particle data to the (ghost) send buffer
  subroutine copy_to_send_ghosts(this, at, i, n, buffer)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    integer, intent(in)                       :: i
    integer, intent(inout)                    :: n
    character(1), intent(inout)               :: buffer(:)

    ! ---

! WAIT FOR DICTIONARY IN ATOMS
!    call pack_buffer(at%data, mask, i, n, buffer)

    buffer(n+1) = 'h'
    n = n+1

!    if (i > at%Ndomain) then
!       buffer(n)  = at%from_rank(i)
!    else
!       buffer(n)  = this%mpi%my_proc
!    endif

!    call print("g-send: " // this%local_to_global(i) // ", " // at%Z(i) // ", " // at%pos(:, i))

  endsubroutine copy_to_send_ghosts


  !% Copy particle data from the (ghost) receive buffer
  subroutine copy_from_recv_ghosts(this, at, n, buffer, off)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    integer, intent(in)                       :: n
    character(1), intent(in)                  :: buffer(:)
    real(DP), intent(in)                      :: off(3)

    ! ---

    integer  :: i

    ! ---

    i = 0
    do while (i < n)
       at%N = at%N+1

! WAIT FOR DICTIONARY IN ATOMS
!       call unpack_buffer(at%data, mask, i, buffer, at%N)

       i  = i + 1
!       at%from_rank(at%N)  = buffer(i)

       this%global_to_local(this%local_to_global(at%N)) = at%N

!       PNC3(at, at%N)  = PNC3(at, at%N) + off
!       POS3(at, at%N)  = cyclic_in_cell(at, PNC3(at, at%N))

#ifdef DEBUG
       if (.not. ( all(POS3(at, at%N) >= 0.0_DP) .and. all(POS3(at, at%N) < (/ at%lattice(1, 1), at%lattice(2, 2), at%lattice(3, 3) /)) )) then
          call particles_dump_info(at, at%N)
          EXIT_ON_ERROR("Particle outside of the simulation domain.", i)
       endif
#endif

!       if (at%index(at%N) == 914635) then
!          write (ilog, '(5X,A,2I5,9F20.10)')  "g-recv: ", at%index(at%N), at%el(at%N), POS3(p, at%N), PNC3(p, at%N), off(:)
!       endif
    enddo

  endsubroutine copy_from_recv_ghosts


  !% Communicate ghost particles to neighboring domains
  !% (former prebinning routine). *bwidth* is the width of
  !% the border.
  subroutine domain_decomposition_communicate_ghosts(this, at, new_list, error)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    logical, intent(in)                       :: new_list
    integer, intent(out), optional :: error

    ! ---

    real(DP)  :: upper(3), lower(3)
    integer   :: i, d, list_off_r, list_off_l
    integer   :: n_send_r, n_send_l, n_recv_r, n_recv_l

    real(DP)  :: off_l(3), off_r(3)

    ! ---

    INIT_ERROR(error)

    call system_timer("domain_decomposition_communicate_ghosts")

    this%nit_g = this%nit_g + 1

    do d = 1, 3

       if (this%periodic(d) .or. this%mpi%my_coords(d) /= 0) then
          lower(d)  = this%lower(d) + (at%g(d, :) .dot. this%border)
       else
          lower(d)  = this%lower(d)
       endif

       if (this%periodic(d) .or. this%mpi%my_coords(d) /= this%decomposition(d)-1) then
          upper(d)  = this%upper(d) - (at%g(d, :) .dot. this%border)
       else
          upper(d)  = this%upper(d)
       endif

    enddo

    do i = at%Ndomain+1, at%N
       this%global_to_local(this%local_to_global(i)) = 0
    enddo

    at%N  = at%Ndomain

    !
    ! Loop over dimensions and distribute particle in the
    ! respective direction
    !

    list_off_r = 0
    list_off_l = 0
    do d = 1, 3

       if (this%decomposition(d) > 1) then

          n_send_r  = 0
          n_send_l  = 0

          if (new_list) then

             this%n_ghosts_r(d)  = 0
             this%n_ghosts_l(d)  = 0

             do i = 1, at%N
!                if (PNC(at, i, d) >= upper(d)) then
                if ((at%g(d, :) .dot. at%pos(:, i)) >= upper(d)) then
                   call copy_to_send_ghosts(this, at, i, n_send_r, this%send_r)

                   this%n_ghosts_r(d) = this%n_ghosts_r(d)+1
                   this%ghosts_r(list_off_r+this%n_ghosts_r(d)) = this%local_to_global(i)

!                else if (PNC(at, i, d) < lower(d)) then
                else if ((at%g(d, :) .dot. at%pos(:, i)) < lower(d)) then
                   call copy_to_send_ghosts(this, at, i, n_send_l, this%send_l)

                   this%n_ghosts_l(d) = this%n_ghosts_l(d)+1
                   this%ghosts_l(list_off_l+this%n_ghosts_l(d)) = this%local_to_global(i)

                endif
             enddo

          else

             do i = 1, this%n_ghosts_r(d)
                call copy_to_send_ghosts( &
                     this, at, &
                     this%global_to_local(this%ghosts_r(list_off_r+i)), &
                     n_send_r, this%send_r)
             enddo

             do i = 1, this%n_ghosts_l(d)
                call copy_to_send_ghosts( &
                     this, at, &
                     this%global_to_local(this%ghosts_l(list_off_l+i)), &
                     n_send_l, this%send_l)
             enddo

          endif

          this%n_send_g_tot = this%n_send_g_tot + this%n_ghosts_r(d) + this%n_ghosts_l(d)

!          call mpi_sendrecv( &
!               this%send_r, n_send_r, MPI_DOUBLE_PRECISION, this%r(d), 0, &
!               this%recv_l, this%ghost_buffer_size*at%maxnatloc, MPI_DOUBLE_PRECISION, this%l(d), 0, &
!               this%comm, status, i)
!          PASS_MPI_ERROR(i, error)
!
!          call mpi_get_count(status, MPI_DOUBLE_PRECISION, n_recv_l, i)
!          PASS_MPI_ERROR(i, error)

          call sendrecv(this%mpi, &
               this%send_r(1:n_send_r), this%r(d), 0, &
               this%recv_l, this%l(d), 0, &
               n_recv_l, error)
          PASS_ERROR(error)

!          call mpi_sendrecv( &
!               this%send_l, n_send_l, MPI_DOUBLE_PRECISION, this%l(d), 1, &
!               this%recv_r, this%ghost_buffer_size*at%maxnatloc, MPI_DOUBLE_PRECISION, this%r(d), 1, &
!               this%comm, status, i)
!          PASS_MPI_ERROR(i, error)
!
!          call mpi_get_count(status, MPI_DOUBLE_PRECISION, n_recv_r, i)
!          PASS_MPI_ERROR(i, error)

          call sendrecv(this%mpi, &
               this%send_l(1:n_send_l), this%l(d), 1, &
               this%recv_r, this%r(d), 1, &
               n_recv_r, error)
          PASS_ERROR(error)

          call print("Received ghost buffers. Sizes: l = " // n_recv_l // ", r = " // n_recv_r)

          this%n_recv_g_tot = this%n_recv_g_tot + n_recv_r/this%ghost_buffer_size + n_recv_l/this%ghost_buffer_size

          off_l    = 0.0_DP
          off_l(d) = this%off_l(d)

          off_r    = 0.0_DP
          off_r(d) = this%off_r(d)

          call copy_from_recv_ghosts(this, at, n_recv_l, this%recv_l, off_l)
          call copy_from_recv_ghosts(this, at, n_recv_r, this%recv_r, off_r)

          list_off_r = list_off_r + this%n_ghosts_r(d)
          list_off_l = list_off_l + this%n_ghosts_l(d)

       endif

    enddo

    call system_timer("domain_decomposition_communicate_ghosts")

  endsubroutine domain_decomposition_communicate_ghosts


  !% Copy forces to the (ghost) send buffer
  subroutine copy_forces_to_send_ghosts(this, at, i, n, buffer)
    implicit none

    type(DomainDecomposition), intent(in)  :: this
    type(Atoms), intent(inout)             :: at
    integer, intent(in)                    :: i
    integer, intent(in)                    :: n
    character(1), intent(inout)            :: buffer(:)

    ! ---

    integer  :: m

    ! ---

    m = n
! WAIT FOR DICTIONARY IN ATOMS
!    call pack_buffer(at%data, mask, i, m, buffer)

    !write (ilog, '(A,I5,3F20.10)')  "Fsend: ", at%index(i), FOR3(p, i)

  endsubroutine copy_forces_to_send_ghosts
  

  !% Copy particle data from the (ghost) receive buffer
  subroutine copy_forces_from_recv_ghosts(this, at, cur, n, buffer)
    implicit none

    type(DomainDecomposition), intent(in)  :: this
    type(Atoms), intent(inout)             :: at
    integer, intent(inout)                 :: cur
    integer, intent(in)                    :: n
    character(1), intent(in)               :: buffer(:)

    ! ---

    integer  :: i

    ! ---

    i = 0
    do while (i < n)
       cur = cur+1

! WAIT FOR DICTIONARY IN ATOMS
!       call unpack_buffer(at%data, mask, i, buffer, cur)
    enddo

  endsubroutine copy_forces_from_recv_ghosts


  !% Communicate forces of ghost particles back.
  !% This is needed for rigid object (i.e., water) or to reduce the
  !% border size in BOPs.
  subroutine domain_decomposition_communicate_forces(this, at, error)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    integer, intent(out), optional :: error

    ! ---

    integer   :: i, d, list_off_r, list_off_l, n_recv_r, n_recv_l, cur

    ! ---

    INIT_ERROR(error)

    call system_timer("domain_decomposition_communicate_forces")

    !
    ! Loop over dimensions and distribute particle in the
    ! respective direction
    !

    list_off_r = 0
    list_off_l = 0
    cur = at%Ndomain
    do d = 1, 3

       if (this%decomposition(d) > 1) then

          do i = 1, this%n_ghosts_r(d)
             call copy_forces_to_send_ghosts(this, at, this%global_to_local(this%ghosts_r(list_off_r+i)), (i-1)*this%backc_buffer_size, this%send_r)
          enddo

          do i = 1, this%n_ghosts_l(d)
             call copy_forces_to_send_ghosts(this, at, this%global_to_local(this%ghosts_l(list_off_l+i)), (i-1)*this%backc_buffer_size, this%send_l)
          enddo

!          call mpi_sendrecv( &
!               this%send_r, this%backc_buffer_size*this%n_ghosts_r(d), MPI_DOUBLE_PRECISION, this%r(d), 0, &
!               this%recv_l, this%backc_buffer_size*at%maxnatloc, MPI_DOUBLE_PRECISION, this%l(d), 0, &
!               this%comm, status, i)
!          PASS_MPI_ERROR(i, error)
!
!          call mpi_get_count(status, MPI_DOUBLE_PRECISION, n_recv_l, i)
!          PASS_MPI_ERROR(i, error)

          call sendrecv(this%mpi, &
               this%send_r(1:this%backc_buffer_size*this%n_ghosts_r(d)), &
               this%r(d), 0, &
               this%recv_l, this%l(d), 0, &
               n_recv_l, error)
          PASS_ERROR(error)

!          call mpi_sendrecv( &
!               this%send_l, this%backc_buffer_size*this%n_ghosts_l(d), MPI_DOUBLE_PRECISION, this%l(d), 1, &
!               this%recv_r, this%backc_buffer_size*at%maxnatloc, MPI_DOUBLE_PRECISION, this%r(d), 1, &
!               this%comm, status, i)
!          PASS_MPI_ERROR(i, error)
!
!          call mpi_get_count(status, MPI_DOUBLE_PRECISION, n_recv_r, i)
!          PASS_MPI_ERROR(i, error)

          call sendrecv(this%mpi, &
               this%send_l(1:this%backc_buffer_size*this%n_ghosts_l(d)), &
               this%l(d), 1, &
               this%recv_r, this%r(d), 1, &
               n_recv_r, error)
          PASS_ERROR(error)

          call copy_forces_from_recv_ghosts(this, at, cur, n_recv_l, this%recv_l)
          call copy_forces_from_recv_ghosts(this, at, cur, n_recv_r, this%recv_r)

          list_off_r = list_off_r + this%n_ghosts_r(d)
          list_off_l = list_off_l + this%n_ghosts_l(d)

       endif

    enddo

    call system_timer("domain_decomposition_communicate_forces")

  endsubroutine domain_decomposition_communicate_forces


  !
  ! Supplement to the Dictionary object
  !

  !% Convert a list of keys to a mask
  subroutine dictionary_keys_to_mask(this, keys, mask, s, error)
    implicit none

    type(Dictionary), intent(in)      :: this          !% Dictionary object
    type(Table), intent(in)           :: keys          !% List of keys
    logical, intent(out)              :: mask(this%N)  !% True if in keys
    integer, intent(out), optional    :: s             !% Size of the buffer
    integer, intent(out), optional :: error

    ! ---

    integer :: i, entry_i

    ! ---

    INIT_ERROR(error)

    mask = .false.
    s = 0
    do i = 1, keys%N
       entry_i = lookup_entry_i(this, keys%str(1, i))
       if (entry_i == -1) then
          RAISE_ERROR("Could not find key '" // trim(keys%str(1, i)) // "'.", error)
       endif

       select case(this%entries(i)%type)

       case (T_REAL_A)
          s = s + sizeof(this%entries(i)%r_a(1))

       case (T_INTEGER_A)
          s = s + sizeof(this%entries(i)%i_a(1))

       case (T_LOGICAL_A)
          s = s + sizeof(this%entries(i)%l_a(1))

       case (T_REAL_A2)
          s = s + sizeof(this%entries(i)%r_a2(:, 1))

       case default
          RAISE_ERROR("Don't know how to handle entry type " // this%entries(i)%type, error)

       endselect

       mask(entry_i) = .true.
    enddo

  endsubroutine dictionary_keys_to_mask


  !% Pack buffer, copy all entries with for which tag is true to the buffer
  subroutine dictionary_pack_buffer(this, mask, &
       data_i, buffer_i, buffer, error)
    implicit none

    type(Dictionary), intent(in)      :: this
    logical, intent(in)               :: mask(:)
    integer, intent(in)               :: data_i
    integer, intent(inout)            :: buffer_i
    character(1), intent(inout)       :: buffer(:)
    integer, intent(out), optional :: error

    ! ---

    integer  :: i

    ! ---

    INIT_ERROR(error)

    do i = 1, this%N
       if (mask(i)) then

          select case(this%entries(i)%type)

          case (T_REAL_A)
             buffer(buffer_i+1:)  = &
                  transfer(this%entries(i)%r_a(data_i), buffer)
             buffer_i          = buffer_i + sizeof(this%entries(i)%r_a(1))

          case (T_INTEGER_A)
             buffer(buffer_i+1:)  = &
                  transfer(this%entries(i)%i_a(data_i), buffer)
             buffer_i          = buffer_i + sizeof(this%entries(i)%i_a(1))

          case (T_LOGICAL_A)
             buffer(buffer_i+1:)  = &
                  transfer(this%entries(i)%l_a(data_i), buffer)
             buffer_i          = buffer_i + sizeof(this%entries(i)%l_a(1))

          case (T_REAL_A2)
             buffer(buffer_i+1:) = &
                  transfer(this%entries(i)%r_a2(:, data_i), buffer)
             buffer_i = buffer_i + sizeof(this%entries(i)%r_a2(:, 1))

          case default
             RAISE_ERROR("Don't know how to handle entry type " // this%entries(i)%type, error)

          endselect

       endif
    enddo

  endsubroutine dictionary_pack_buffer


  !% Unpack buffer, copy all entries with for which mask is true from the buffer
  subroutine dictionary_unpack_buffer(this, mask, &
       buffer_i, buffer, data_i, error)
    implicit none

    type(Dictionary), intent(inout)   :: this
    logical, intent(in)               :: mask(:)
    integer, intent(inout)            :: buffer_i
    character(1), intent(in)          :: buffer(:)
    integer, intent(in)               :: data_i
    integer, intent(out), optional :: error

    ! ---

    integer  :: i, ndims

    ! ---

    INIT_ERROR(error)

    do i = 1, this%N
       if (mask(i)) then

          select case(this%entries(i)%type)

          case (T_REAL_A)
             this%entries(i)%r_a(data_i) = &
                  transfer(buffer(buffer_i+1:), 1.0_DP)
             buffer_i = buffer_i + sizeof(this%entries(i)%r_a(1))

          case (T_INTEGER_A)
             this%entries(i)%i_a(data_i) = &
                  transfer(buffer(buffer_i+1:), 1)
             buffer_i = buffer_i + sizeof(this%entries(i)%i_a(1))

          case (T_LOGICAL_A)
             this%entries(i)%l_a(data_i) = &
                  transfer(buffer(buffer_i+1:), .true.)
             buffer_i = buffer_i + sizeof(this%entries(i)%l_a(1))

          case (T_REAL_A2)
             ndims = size(this%entries(i)%r_a2, 1)
             this%entries(i)%r_a2(:, data_i) = &
                  transfer(buffer(buffer_i+1:), 1.0_DP, ndims)
             buffer_i = buffer_i + sizeof(this%entries(i)%r_a2(:, 1))

          case default
             RAISE_ERROR("Don't know how to handle entry type " // this%entries(i)%type, error)
             
          endselect

       endif
    enddo

  endsubroutine dictionary_unpack_buffer

endmodule DomainDecomposition_module
