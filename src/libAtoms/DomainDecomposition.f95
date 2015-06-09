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
  use error_module
  use system_module
  use extendable_str_module
  use MPI_context_module
  use linearalgebra_module
  use table_module
  use atoms_types_module
  use dictionary_module

  implicit none

  private

  public :: DomainDecomposition

  public :: initialise
  interface initialise
     module procedure domaindecomposition_initialise
  endinterface

  public :: finalise
  interface finalise
     module procedure domaindecomposition_finalise
  endinterface

  interface update_sendrecv_masks
     module procedure domaindecomposition_update_sendrecv_masks
  endinterface

  public :: enable
  interface enable
     module procedure domaindecomposition_enable
  endinterface enable

  public :: disable
  interface disable
     module procedure domaindecomposition_disable
  endinterface disable

  public :: allocate
  interface allocate
     module procedure domaindecomposition_allocate
  endinterface allocate

  interface is_in_domain
     module procedure domaindecomposition_is_in_domain
  endinterface is_in_domain

  interface sdompos
     module procedure sdompos1, sdompos3
  endinterface

  public :: set_comm_property
  interface set_comm_property
     module procedure domaindecomposition_set_comm_property
  endinterface set_comm_property

  public :: set_border
  interface set_border
     module procedure domaindecomposition_set_border
  endinterface

  public :: comm_atoms
  interface comm_atoms
     module procedure domaindecomposition_comm_atoms
  endinterface

  public :: comm_ghosts
  interface comm_ghosts
     module procedure domaindecomposition_comm_ghosts
  endinterface

  public :: comm_reverse
  interface comm_reverse
     module procedure domaindecomposition_comm_reverse
  endinterface

  public :: comm_atoms_to_all
  interface comm_atoms_to_all
     module procedure domaindecomposition_comm_atoms_to_all
  endinterface

  public :: deepcopy
  interface deepcopy
     module procedure domaindecomposition_deepcopy
  endinterface deepcopy


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
  subroutine domaindecomposition_initialise(this, &
       decomposition, verlet_shell, mode, error)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    integer,         optional, intent(in)     :: decomposition(3)
    real(dp),        optional, intent(in)     :: verlet_shell
    integer,         optional, intent(in)     :: mode
    integer,         optional, intent(out)    :: error

    ! ---

    INIT_ERROR(error)

    call print("DomainDecomposition : initialise", PRINT_VERBOSE)

    this%decomposition = (/ 1, 1, 1 /)
    if (present(decomposition)) then
       this%decomposition  = decomposition
    endif

    this%verlet_shell = 0.0_DP
    if (present(verlet_shell)) then
       this%verlet_shell   = verlet_shell
    endif

    this%mode = DD_WRAP_TO_CELL
    if (present(mode)) then
       if (mode /= DD_WRAP_TO_CELL .and. mode /= DD_WRAP_TO_DOMAIN) then
          RAISE_ERROR("Unknown domain decomposition mode '" // mode // "'.", error)
       endif

       this%mode = mode
    else
       this%mode = DD_WRAP_TO_CELL
    endif

    this%decomposed        = .false.
    this%requested_border  = 0.0_DP

    this%n_send_p_tot  = 0
    this%n_recv_p_tot  = 0
    this%n_send_g_tot  = 0
    this%n_recv_g_tot  = 0
    this%nit_p         = 0
    this%nit_g         = 0

    ! Initialise property lists
    call initialise(this%atoms_properties)
    call initialise(this%ghost_properties)
    call initialise(this%reverse_properties)

  endsubroutine domaindecomposition_initialise


  !% Update send and receive masks, only necessary after properties have
  !% been added to or removed from the Atoms object
  subroutine domaindecomposition_update_sendrecv_masks(this, at, error)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(in)                   :: at
    integer, intent(out), optional            :: error

    ! ---

    integer  :: s

    ! ---

    INIT_ERROR(error)

    ! XXX FIXME This should become unecessary if DynamicalSystem only
    ! pointer to an atoms object
    if (.not. assign_pointer(at%properties, "local_to_global", &
                             this%local_to_global)) then
       RAISE_ERROR("Could not find 'local_to_global'", error)
    endif

    if (allocated(this%atoms_mask)) then
       if (size(this%atoms_mask) < at%properties%N) then
          deallocate(this%atoms_mask)
          this%atoms_buffer_size = 0
       endif
    endif

    if (allocated(this%ghost_mask)) then
       if (size(this%ghost_mask) < at%properties%N) then
          deallocate(this%ghost_mask)
          this%ghost_buffer_size = 0
       endif
    endif

    if (allocated(this%reverse_mask)) then
       if (size(this%reverse_mask) < at%properties%N) then
          deallocate(this%reverse_mask)
          this%reverse_buffer_size = 0
       endif
    endif

    if (.not. allocated(this%atoms_mask) .and. this%atoms_properties%N > 0) then
       allocate(this%atoms_mask(at%properties%N))
       call keys_to_mask( &
            at%properties, this%atoms_properties, this%atoms_mask, &
            s = this%atoms_buffer_size, error = error)
       PASS_ERROR(error)

       call print("DomainDecomposition : atoms_buffer_size  = " // &
            this%atoms_buffer_size, PRINT_VERBOSE)
    endif
    if (.not. allocated(this%ghost_mask) .and. this%ghost_properties%N > 0) then
       allocate(this%ghost_mask(at%properties%N))
       call keys_to_mask( &
            at%properties, this%ghost_properties, this%ghost_mask, &
            s = this%ghost_buffer_size, error = error)
       PASS_ERROR(error)

       call print("DomainDecomposition : ghost_buffer_size  = " // &
            this%ghost_buffer_size, PRINT_VERBOSE)
    endif
    if (.not. allocated(this%reverse_mask) .and. this%reverse_properties%N > 0) then
       allocate(this%reverse_mask(at%properties%N))
       call keys_to_mask( &
            at%properties, this%reverse_properties, this%reverse_mask, &
            s = this%reverse_buffer_size, error = error)
       PASS_ERROR(error)

       call print("DomainDecomposition : reverse_buffer_size  = " // &
            this%reverse_buffer_size, PRINT_VERBOSE)
    endif

    s = max(this%atoms_buffer_size, this%ghost_buffer_size) * at%Nbuffer

    if (allocated(this%send_l)) then
       if (s > size(this%send_l, 1)) then
          deallocate(this%send_l)
          deallocate(this%send_r)
          deallocate(this%recv_l)
          deallocate(this%recv_r)
       endif
    endif
       
    if (.not. allocated(this%send_l)) then
       allocate(this%send_l(s))
       allocate(this%send_r(s))
       allocate(this%recv_l(s))
       allocate(this%recv_r(s))
    endif

  endsubroutine domaindecomposition_update_sendrecv_masks


  !% Enable domain decomposition, after this call every process
  !% remains only with the atoms in its local domain.
  subroutine domaindecomposition_enable(this, at, &
       mpi, decomposition, verlet_shell, mode, error)
    implicit none

    type(DomainDecomposition),   intent(inout)  :: this
    type(Atoms),                 intent(inout)  :: at
    type(MPI_context), optional, intent(in)     :: mpi
    integer,           optional, intent(in)     :: decomposition(3)
    real(dp),          optional, intent(in)     :: verlet_shell
    integer,           optional, intent(in)     :: mode
    integer,           optional, intent(out)    :: error

    ! ---

    integer   :: i, j

    ! ---

    INIT_ERROR(error)

    call print("DomainDecomposition : enable", PRINT_VERBOSE)

    if (present(decomposition)) then
       this%decomposition  = decomposition
    endif

    if (present(verlet_shell)) then
       this%verlet_shell   = verlet_shell
    endif

    if (present(mode)) then
       if (mode /= DD_WRAP_TO_CELL .and. mode /= DD_WRAP_TO_DOMAIN) then
          RAISE_ERROR("Unknown domain decomposition mode '" // mode // "'.", error)
       endif

       this%mode = mode
    else
       this%mode = DD_WRAP_TO_CELL
    endif

    call print("DomainDecomposition : Parallelisation using domain decomposition over " // &
         this%decomposition // "domains.", PRINT_NORMAL)

    call initialise(this%mpi, &
         context  = mpi, &
         dims     = this%decomposition, &
         error    = error)
    PASS_ERROR(error)

    call print("DomainDecomposition : " // this%mpi%n_procs // " MPI processes available.", PRINT_NORMAL)

    if (this%decomposition(1)*this%decomposition(2)*this%decomposition(3) /= &
         this%mpi%n_procs) then
       RAISE_ERROR("Decomposition geometry requires " // this%decomposition(1)*this%decomposition(2)*this%decomposition(3) // " processes, however, MPI returns " // this%mpi%n_procs // " processes.", error)
    endif

! For now this needs to be true
!    this%periodic          = at%periodic
    this%periodic          = .true.
    this%requested_border  = 0.0_DP

    call print("DomainDecomposition : coords             = ( " // this%mpi%my_coords // ")", PRINT_VERBOSE)

    do i = 1, 3
       call cart_shift( &
            this%mpi, i-1, 1, this%l(i), this%r(i), error)
       PASS_ERROR(error)
    enddo

    this%lower  = (1.0_DP*this%mpi%my_coords)     / this%decomposition
    this%upper  = (1.0_DP*(this%mpi%my_coords+1)) / this%decomposition
    this%center = (this%lower + this%upper)/2

    call print("DomainDecomposition : lower              = ( " // this%lower // " )", PRINT_VERBOSE)
    call print("DomainDecomposition : upper              = ( " // this%upper // " )", PRINT_VERBOSE)
    call print("DomainDecomposition : center             = ( " // this%center // " )", PRINT_VERBOSE)

    this%n_send_p_tot  = 0
    this%n_recv_p_tot  = 0
    this%n_send_g_tot  = 0
    this%n_recv_g_tot  = 0
    this%nit_p         = 0
    this%nit_g         = 0

    ! allocate internal buffers
    call allocate(this, at, error=error)
    PASS_ERROR(error)

    ! at is initialised, this means we need to retain only the atoms that
    ! are in the current domain.
    j = 0
    do i = 1, at%N
       ! If this is in the domain, keep it
       if (is_in_domain(this, at, at%pos(:, i))) then
          j = j+1
          if (i /= j) then
             call copy_entry(at, i, j)
          endif
       endif
    enddo
    at%N = j
    at%Ndomain = j

    this%decomposed = .true.

  endsubroutine domaindecomposition_enable


  !% Disable domain decomposition, after this call every process
  !% retains an identical copy of the system.
  subroutine domaindecomposition_disable(this, at, error)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    integer, intent(out), optional            :: error

    ! ---

    INIT_ERROR(error)

    if (.not. this%decomposed) then
       RAISE_ERROR("DomainDecomposition-object is not in a domain-decomposed state. Cannot disable domain decomposition.", error)
    endif

    call comm_atoms_to_all(this, at, error=error)
    PASS_ERROR(error)

    this%decomposed  = .false.

    this%local_to_global = 0

    if (allocated(this%global_to_local))  deallocate(this%global_to_local)
    if (allocated(this%ghosts_r))         deallocate(this%ghosts_r)
    if (allocated(this%ghosts_l))         deallocate(this%ghosts_l)

  endsubroutine domaindecomposition_disable


  !% Allocate internal buffer
  subroutine domaindecomposition_allocate(this, at, range, error)
    implicit none

    type(DomainDecomposition),           intent(inout)  :: this
    type(Atoms),                         intent(inout)  :: at
    integer,                   optional, intent(in)     :: range(2)
    integer,                   optional, intent(inout)  :: error

    ! ---

    integer   :: i, d, my_range(2)
    real(dp)  :: l(3)

    ! ---

    if ((at%lattice(:, 1) .dot. at%lattice(:, 2)) > 1e-6 .or. &
        (at%lattice(:, 2) .dot. at%lattice(:, 3)) > 1e-6 .or. &
        (at%lattice(:, 3) .dot. at%lattice(:, 1)) > 1e-6) then
       RAISE_ERROR("Only orthorhombic lattices are support for now.", error)
    endif

    l = (/ at%lattice(1, 1), at%lattice(2, 2), at%lattice(3, 3) /)

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
    call print("DomainDecomposition : periodic (par.)    = ( " // this%periodic // " )", PRINT_VERBOSE)
!    call print("periodic (local)   = ( " // this%locally_periodic // " )", PRINT_VERBOSE)

    if (this%mode == DD_WRAP_TO_CELL) then
       this%off_l  = 0.0_DP
       this%off_r  = 0.0_DP
    endif

    call print("DomainDecomposition : off_l              = ( " // this%off_l // " )", PRINT_VERBOSE)
    call print("DomainDecomposition : off_r              = ( " // this%off_r // " )", PRINT_VERBOSE)

    if (present(range)) then
       my_range = range
       if (my_range(2) - my_range(1) + 1 /= at%N) then
          RAISE_ERROR("*range* and number of atoms do not match.", error)
       endif
    else
       my_range = (/ 1, at%N /)
    endif

    ! Additional properties (integers)
    call add_property(at, "local_to_global", 0, &
         ptr=this%local_to_global, error=error)
    PASS_ERROR(error)

    ! Allocate local buffers
    allocate(this%global_to_local(at%Nbuffer))
    allocate(this%ghosts_r(at%Nbuffer))
    allocate(this%ghosts_l(at%Nbuffer))

    do i = 1, at%N
       !% Assign global indices
       this%local_to_global(i) = i+my_range(1)-1
       this%global_to_local(i+my_range(1)-1) = i
    enddo

    this%Ntotal = sum(this%mpi, at%N, error)
    PASS_ERROR(error)

    call print("DomainDecomposition : Ntotal = " // this%Ntotal, PRINT_VERBOSE)

    call set_border(this, at, this%requested_border, error=error)
    PASS_ERROR(error)

  endsubroutine domaindecomposition_allocate


  !% Is the position in the current domain?
  function domaindecomposition_is_in_domain(this, at, r) result(id)
    implicit none

    type(DomainDecomposition), intent(in)  :: this
    type(Atoms), intent(in)                :: at
    real(dp), intent(in)                   :: r(3)

    logical                                :: id

    ! ---

    real(dp)  :: s(3)

    ! ---

    s = at%g .mult. r
    id = all(s >= this%lower) .and. all(s < this%upper)

  endfunction domaindecomposition_is_in_domain


  !% Set which properties to communicate on which events
  !% comm_atoms:   Communicate when atom is moved to different domain.
  !%               Forces, for example, may be excluded since they are updated
  !%               on every time step.
  !% comm_ghosts:  Communicate when atoms becomes a ghost on some domain
  !%               Masses, for example, might be excluded since atoms are
  !%               propagated on the domain they reside in only.
  !% comm_reverse: Communicate back from ghost atoms to the original domain atom
  !%               and accumulate
  subroutine domaindecomposition_set_comm_property(this, propname, &
       comm_atoms, comm_ghosts, comm_reverse)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    character(*),              intent(in)     :: propname
    logical,         optional, intent(in)     :: comm_atoms
    logical,         optional, intent(in)     :: comm_ghosts
    logical,         optional, intent(in)     :: comm_reverse

    ! ---

    if (present(comm_atoms)) then
       if (comm_atoms) then
          call set_value(this%atoms_properties, propname, .true.)
       else
          call remove_value(this%atoms_properties, propname)
       endif
    endif
    if (present(comm_ghosts)) then
       if (comm_ghosts) then
          call set_value(this%ghost_properties, propname, .true.)
       else
          call remove_value(this%ghost_properties, propname)
       endif
    endif
    if (present(comm_reverse)) then
       if (comm_reverse) then
          call set_value(this%reverse_properties, propname, .true.)
       else
          call remove_value(this%reverse_properties, propname)
       endif
    endif

  endsubroutine domaindecomposition_set_comm_property


  !% Set the communication border
  subroutine domaindecomposition_set_border(this, at, &
       border, verlet_shell, error)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(in)                   :: at
    real(dp), intent(in)                      :: border
    real(dp), intent(in), optional            :: verlet_shell
    integer, intent(out), optional :: error

    ! ---

    integer  :: d

    ! ---

    INIT_ERROR(error)

    call print("DomainDecomposition : set_border", PRINT_VERBOSE)

    if (present(verlet_shell)) then
       this%verlet_shell  = verlet_shell
    endif

    this%requested_border  = max(this%requested_border, border)
    this%border            = this%requested_border + this%verlet_shell

    call print("DomainDecomposition : requested_border  = " // this%requested_border, PRINT_VERBOSE)
    call print("DomainDecomposition : verlet_shell      = " // this%verlet_shell, PRINT_VERBOSE)
    call print("DomainDecomposition : border            = " // this%border, PRINT_VERBOSE)

    if (any(((at%lattice .mult. (this%upper - this%lower)) < 2*this%border) .and. this%periodic)) then
       RAISE_ERROR("Domain smaller than twice the border. This does not work (yet). border = " // this%border // ", lattice = " // at%lattice(:, 1) // ", " // at%lattice(:, 2) // ", " // at%lattice(:, 3) // ", lower = " // this%lower // ", upper = " // this%upper, error)
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

    call print("DomainDecomposition : lower_with_border  = ( " // this%lower_with_border // " )", PRINT_VERBOSE)
    call print("DomainDecomposition : upper_with_border  = ( " // this%upper_with_border // " )", PRINT_VERBOSE)

  endsubroutine domaindecomposition_set_border

  
  !% Free memory, clean up
  subroutine domaindecomposition_finalise(this)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this

    ! ---

    call print("DomainDecomposition : finalise", PRINT_VERBOSE)

    if (allocated(this%send_l))  deallocate(this%send_l)
    if (allocated(this%send_r))  deallocate(this%send_r)
    if (allocated(this%recv_l))  deallocate(this%recv_l)
    if (allocated(this%recv_r))  deallocate(this%recv_r)

    if (allocated(this%global_to_local))  deallocate(this%global_to_local)
    if (allocated(this%ghosts_r))         deallocate(this%ghosts_r)
    if (allocated(this%ghosts_l))         deallocate(this%ghosts_l)

    if (allocated(this%atoms_mask))    deallocate(this%atoms_mask)
    if (allocated(this%ghost_mask))    deallocate(this%ghost_mask)
    if (allocated(this%reverse_mask))  deallocate(this%reverse_mask)

    if (this%n_send_p_tot > 0 .and. this%n_recv_p_tot > 0 .and. &
        this%n_send_g_tot > 0 .and. this%n_recv_g_tot > 0) then
       call print("DomainDecomposition : Average number of particles sent/received per iteration:", PRINT_VERBOSE)
       call print("DomainDecomposition : Particles send  = " // (1.0_DP*this%n_send_p_tot)/this%nit_p, PRINT_VERBOSE)
       call print("DomainDecomposition : Particles recv  = " // (1.0_DP*this%n_recv_p_tot)/this%nit_p, PRINT_VERBOSE)
       call print("DomainDecomposition : Ghosts send     = " // (1.0_DP*this%n_send_g_tot)/this%nit_g, PRINT_VERBOSE)
       call print("DomainDecomposition : Ghosts recv     = " // (1.0_DP*this%n_recv_g_tot)/this%nit_g, PRINT_VERBOSE)
    endif

    call finalise(this%atoms_properties)
    call finalise(this%ghost_properties)
    call finalise(this%reverse_properties)

    call finalise(this%mpi)

  endsubroutine domaindecomposition_finalise


  !% Compute the scaled position as the minimum image from the domain center
  !% for a single coordinate only
  function sdompos1(this, at, r, d)
    implicit none

    type(DomainDecomposition), intent(in)  :: this
    type(Atoms), intent(in)                :: at
    real(dp), intent(in)                   :: r(3)
    integer, intent(in)                    :: d

    real(dp)                               :: sdompos1

    ! ---

    real(dp)  :: s

    ! ---

    ! center is in scaled coordinates
    s         = (at%g(d, :) .dot. r) - this%center(d)
    sdompos1  = s - nint(s) + this%center(d)

  endfunction sdompos1


  !% Compute the scaled position as the minimum image from the domain center
  function sdompos3(this, at, r)
    implicit none

    type(DomainDecomposition), intent(in)  :: this
    type(Atoms), intent(in)                :: at
    real(dp), intent(in)                   :: r(3)

    real(dp)                               :: sdompos3(3)

    ! ---

    real(dp)  :: s(3)

    ! ---

    ! center is in scaled coordinates
    s         = (at%g .mult. r) - this%center
    sdompos3  = s - nint(s) + this%center

  endfunction sdompos3


  !% Add particle data to the send buffer
  subroutine pack_atoms_buffer(this, at, i, n, buffer)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    integer, intent(in)                       :: i
    integer, intent(inout)                    :: n
    character(1), intent(inout)               :: buffer(:)

    ! ---

!    write (*, *)  "Packing: i = " // i // ", n = " // n // ", index = " // this%local_to_global(i) // ", pos = " // at%pos(:, i)

    call pack_buffer(at%properties, this%atoms_mask, i, n, buffer)
    this%global_to_local(this%local_to_global(i)) = 0 ! This one is gone

  endsubroutine pack_atoms_buffer


  !% Copy particle data from the receive buffer
  subroutine unpack_atoms_buffer(this, at, n, buffer)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    integer, intent(in)                       :: n
    character(1), intent(in)                  :: buffer(:)

    ! ---

    integer  :: i
!    real(dp), pointer :: pos(:, :)

    ! ---

!    if (.not. assign_pointer(at%properties, "pos", pos)) then
!       call system_abort("...")
!    endif

!    do i = 1, n
    i = 0
    do while (i < n)
       at%Ndomain = at%Ndomain+1

       call unpack_buffer(at%properties, this%atoms_mask, i, buffer, at%Ndomain)
!       write (*, *) "Unpacked: n = " // at%Ndomain // ", index = " // this%local_to_global(at%Ndomain) // ", pos = " // pos(:, at%Ndomain)
       this%global_to_local(this%local_to_global(at%Ndomain))  = at%Ndomain

!       call print("n = " // at%Ndomain // ", pos = " // at%pos(:, at%Ndomain))
    enddo

  endsubroutine unpack_atoms_buffer


  !% Copy particle data from the receive buffer
  subroutine unpack_atoms_buffer_if_in_domain(this, at, n, buffer, error)
    implicit none

    type(DomainDecomposition),           intent(inout)  :: this
    type(Atoms),                         intent(inout)  :: at
    integer,                             intent(in)     :: n
    character(1),                        intent(in)     :: buffer(:)
    integer,                   optional, intent(out)    :: error

    ! ---

    integer            :: i
    real(dp)           :: s(3)
    real(dp), pointer  :: pos(:, :)

    ! ---

    INIT_ERROR(error)

    if (.not. assign_pointer(at%properties, "pos", pos)) then
       RAISE_ERROR("Property 'pos' not found.", error)
    endif

    i = 0
    do while (i < n)
       at%Ndomain = at%Ndomain+1

       call unpack_buffer(at%properties, this%atoms_mask, i, buffer, at%Ndomain)

       s = sdompos(this, at, pos(1:3, at%Ndomain))
       if (all(s >= this%lower) .and. all(s < this%upper)) then
          this%global_to_local(this%local_to_global(at%Ndomain))  = at%Ndomain
       else
          at%Ndomain = at%Ndomain-1
       endif
    enddo

  endsubroutine unpack_atoms_buffer_if_in_domain


  !% Communicate particles which left the domains
  !% to the neighboring domains (former order routine)
  subroutine domaindecomposition_comm_atoms(this, at, error)
    implicit none

    type(DomainDecomposition),           intent(inout)  :: this
    type(Atoms),                         intent(inout)  :: at
    integer,                   optional, intent(out)    :: error

    ! ---

    ! 
    ! General and auxiliary variables      
    !

    integer   :: i, d
    integer   :: last_N

    ! 
    ! Structure variables, mpi and system structure
    !

    integer   :: n_send_l, n_send_r, n_recv_l, n_recv_r

    real(dp)  :: off_l(3), off_r(3), s

    ! ---

    INIT_ERROR(error)

    call print("DomainDecomposition : comm_atoms", PRINT_VERBOSE)

    call system_timer("domaindecomposition_comm_atoms")

    ! Update internal buffers
    call update_sendrecv_masks(this, at)

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

          last_N  = at%Ndomain
          at%Ndomain   = 0

          n_send_r   = 0
          n_send_l   = 0

          do i = 1, last_N

             s = sdompos(this, at, at%pos(:, i), d)
!             call print("d = " // d // ", i = " // i // ", s = " // s // ", lower = " // this%lower(d) // ", upper = " // this%upper(d))
             if (s >= this%upper(d)) then
                ! Send to the right

                call pack_atoms_buffer(this, at, i, n_send_r, this%send_r)

             else if (s < this%lower(d)) then
                ! Send to the left

                call pack_atoms_buffer(this, at, i, n_send_l, this%send_l)

             else
                ! Keep on this processor and reorder

                at%Ndomain = at%Ndomain+1

                if (at%Ndomain /= i) then
                   call copy_entry(at, i, at%Ndomain)
                endif

             endif

          enddo

          this%n_send_p_tot = this%n_send_p_tot + n_send_r + n_send_l

          call sendrecv(this%mpi, &
               this%send_r(1:n_send_r), this%r(d), 0, &
               this%recv_l, this%l(d), 0, &
               n_recv_l, error)
          PASS_ERROR(error)

          call sendrecv(this%mpi, &
               this%send_l(1:n_send_l), this%l(d), 1, &
               this%recv_r, this%r(d), 1, &
               n_recv_r, error)
          PASS_ERROR(error)

          call print("DomainDecomposition : Send state buffers. Sizes: l = " // n_send_l/this%atoms_buffer_size // ", r = " // n_send_r/this%atoms_buffer_size, PRINT_NERD)

          call print("DomainDecomposition : Received state buffers. Sizes: l = " // n_recv_l/this%atoms_buffer_size // ", r = " // n_recv_r/this%atoms_buffer_size, PRINT_NERD)

          this%n_recv_p_tot = this%n_recv_p_tot + n_recv_r/this%atoms_buffer_size + n_recv_l/this%atoms_buffer_size

          off_l    = 0.0_DP
          ! This will be done by inbox
          off_l(d) = this%off_l(d)

          off_r   = 0.0_DP
          ! This will be done by inbox
          off_r(d) = this%off_r(d)

          call unpack_atoms_buffer(this, at, n_recv_l, this%recv_l)
          call unpack_atoms_buffer(this, at, n_recv_r, this%recv_r)

       endif

    enddo

    at%N = at%Ndomain
    if (at%Ndomain /= last_N) then
       call print("DomainDecomposition : Repointing atoms object", PRINT_ANAL)
       call atoms_repoint(at)
    endif

    this%Ntotal = sum(this%mpi, at%N, error)

!    call print("DomainDecomposition : Total number of atoms = " // this%Ntotal, PRINT_VERBOSE)

    call system_timer("domaindecomposition_comm_atoms")

  endsubroutine domaindecomposition_comm_atoms


  !% Communicate all particles from this domain to all other domains
  subroutine domaindecomposition_comm_atoms_to_all(this, at, error)
    implicit none

    type(DomainDecomposition),           intent(inout)  :: this
    type(Atoms),                         intent(inout)  :: at
    integer,                   optional, intent(out)    :: error

    ! ---

    ! 
    ! General and auxiliary variables      
    !

    integer   :: i, j
    integer   :: last_N

    ! 
    ! Structure variables, mpi and system structure
    !

    integer   :: n_send, n_recv

    real(dp)  :: s(3)

    ! ---

    INIT_ERROR(error)

    call print("DomainDecomposition : comm_atoms_to_all", PRINT_VERBOSE)
    call print("DomainDecomposition : decomposed = " // this%decomposed, &
         PRINT_VERBOSE)

    ! Update internal buffers
    call update_sendrecv_masks(this, at)

    !
    ! Loop over dimensions and distribute particle in the
    ! respective direction
    !

    last_N      = at%Ndomain
    at%Ndomain  = 0

    n_send      = 0

    do i = 1, last_N

       ! Send always
       call pack_atoms_buffer(this, at, i, n_send, this%send_r)

       s = sdompos(this, at, at%pos(:, i))
       if (.not. this%decomposed .or. &
           all(s >= this%lower) .and. all(s < this%upper)) then
          ! Keep on this processor and reorder

          at%Ndomain = at%Ndomain+1

          if (at%Ndomain /= i) then
             call copy_entry(at, i, at%Ndomain)
          endif

       endif

    enddo

    do j = 0, this%mpi%n_procs-1

       if (j == this%mpi%my_proc) then
          call bcast(this%mpi, n_send, root=j, error=error)
          PASS_ERROR(error)

!          call print("n_send = " // n_send // ", size(send_r) = " // size(this%send_r), PRINT_ALWAYS)

          call bcast(this%mpi, this%send_r(1:n_send), root=j, error=error)
          PASS_ERROR(error)
       else
          call bcast(this%mpi, n_recv, root=j, error=error)
          PASS_ERROR(error)

!          call print("n_recv = " // n_recv // ", size(recv_r) = " // size(this%recv_r), PRINT_ALWAYS)
!          call print("n = " // (n_recv/this%atoms_buffer_size), PRINT_ALWAYS)

          call bcast(this%mpi, this%recv_r(1:n_recv), root=j, error=error)
          PASS_ERROR(error)

          if (this%decomposed) then
             call unpack_atoms_buffer_if_in_domain(this, at, &
                  n_recv, this%recv_r, &
                  error = error)
             PASS_ERROR(error)
          else
             call unpack_atoms_buffer(this, at, n_recv, this%recv_r)
          endif
       endif

    enddo

    call print("DomainDecomposition : " // &
         "Ndomain = " // at%Ndomain // ", Nbuffer = " // at%Nbuffer, &
         PRINT_VERBOSE)

    if (this%decomposed) then
       call print("DomainDecomposition : total number of atoms = " // &
            sum(this%mpi, at%Ndomain), PRINT_VERBOSE)
    endif

    at%N = at%Ndomain
    if (at%Ndomain /= last_N) then
       call print("DomainDecomposition : Repointing atoms object", PRINT_ANAL)
       call atoms_repoint(at)
    endif

    call print("DomainDecomposition : Sorting atoms object", PRINT_ANAL)

    ! Sort atoms by global index
    call atoms_sort(at, "local_to_global", error=error)
    PASS_ERROR(error)

  endsubroutine domaindecomposition_comm_atoms_to_all


  !% Copy particle data from the (ghost) receive buffer
  subroutine unpack_ghost_buffer(this, at, n, buffer, off)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms),               intent(inout)  :: at
    integer,                   intent(in)     :: n
    character(1),              intent(in)     :: buffer(:)
    real(dp),                  intent(in)     :: off(3)

    ! ---

    integer  :: i

    ! ---

    i = 0
    do while (i < n)
       at%N = at%N+1
       call unpack_buffer(at%properties, this%ghost_mask, i, buffer, at%N)
       at%pos(:, at%N)  = at%pos(:, at%N) + off
       this%global_to_local(this%local_to_global(at%N)) = at%N
    enddo

  endsubroutine unpack_ghost_buffer


  !% Communicate ghost particles to neighboring domains
  !% (former prebinning routine). *bwidth* is the width of
  !% the border.
  subroutine domaindecomposition_comm_ghosts(this, at, new_list, error)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    logical, intent(in)                       :: new_list
    integer, intent(out), optional            :: error

    ! ---

    real(dp)  :: upper(3), lower(3)
    integer   :: i, d, list_off_r, list_off_l, last_N
    integer   :: n_send_r, n_send_l, n_recv_r, n_recv_l

    real(dp)  :: off_l(3), off_r(3), s

    ! ---

    INIT_ERROR(error)

    call print("DomainDecomposition : comm_ghosts", PRINT_NERD)

    call update_sendrecv_masks(this, at)

    call system_timer("domaindecomposition_comm_ghosts")

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

    call print("DomainDecomposition : upper = " // upper, PRINT_NERD)
    call print("DomainDecomposition : lower = " // lower, PRINT_NERD)

    last_N = at%n
    do i = at%Ndomain+1, at%N
       this%global_to_local(this%local_to_global(i)) = 0
    enddo
    at%N = at%Ndomain

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
                s = sdompos(this, at, at%pos(:, i), d)
                if (s >= upper(d)) then
                   call pack_buffer(at%properties, this%ghost_mask, &
                        i, &
                        n_send_r, this%send_r)
                   this%n_ghosts_r(d) = this%n_ghosts_r(d)+1
                   this%ghosts_r(list_off_r+this%n_ghosts_r(d)) = this%local_to_global(i)

                else if (s < lower(d)) then
                   call pack_buffer(at%properties, this%ghost_mask, &
                        i, &
                        n_send_l, this%send_l)

                   this%n_ghosts_l(d) = this%n_ghosts_l(d)+1
                   this%ghosts_l(list_off_l+this%n_ghosts_l(d)) = this%local_to_global(i)

                endif
             enddo

          else

             do i = 1, this%n_ghosts_r(d)
                call pack_buffer(at%properties, this%ghost_mask, &
                     this%global_to_local(this%ghosts_r(list_off_r+i)), &
                     n_send_r, this%send_r)
             enddo

             do i = 1, this%n_ghosts_l(d)
                call pack_buffer(at%properties, this%ghost_mask, &
                     this%global_to_local(this%ghosts_l(list_off_l+i)), &
                     n_send_l, this%send_l)
             enddo

          endif

          this%n_send_g_tot = this%n_send_g_tot + this%n_ghosts_r(d) + this%n_ghosts_l(d)

          call sendrecv(this%mpi, &
               this%send_r(1:n_send_r), this%r(d), 0, &
               this%recv_l, this%l(d), 0, &
               n_recv_l, error)
          PASS_ERROR(error)

          call sendrecv(this%mpi, &
               this%send_l(1:n_send_l), this%l(d), 1, &
               this%recv_r, this%r(d), 1, &
               n_recv_r, error)
          PASS_ERROR(error)

          call print("DomainDecomposition : Send ghost buffers. Sizes: l = " // n_send_l/this%ghost_buffer_size // ", r = " // n_send_r/this%ghost_buffer_size, PRINT_NERD)
          call print("DomainDecomposition : Received ghost buffers. Sizes: l = " // n_recv_l/this%ghost_buffer_size // ", r = " // n_recv_r/this%ghost_buffer_size, PRINT_NERD)

          this%n_recv_g_tot = this%n_recv_g_tot + &
               n_recv_r/this%ghost_buffer_size + &
               n_recv_l/this%ghost_buffer_size

          off_l    = 0.0_DP
          off_l(d) = this%off_l(d)

          off_r    = 0.0_DP
          off_r(d) = this%off_r(d)

          call unpack_ghost_buffer(this, at, n_recv_r, this%recv_r, off_r)
          call unpack_ghost_buffer(this, at, n_recv_l, this%recv_l, off_l)

          list_off_r = list_off_r + this%n_ghosts_r(d)
          list_off_l = list_off_l + this%n_ghosts_l(d)

       endif

    enddo

    if (at%N /= last_N) then
       call print("DomainDecomposition : Repointing atoms object", PRINT_ANAL)
       call atoms_repoint(at)
    endif

    call system_timer("domaindecomposition_comm_ghosts")

  endsubroutine domaindecomposition_comm_ghosts


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
    call pack_buffer(at%properties, this%reverse_mask, i, m, buffer)

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

       call unpack_buffer(at%properties, this%reverse_mask, i, buffer, cur)
    enddo

  endsubroutine copy_forces_from_recv_ghosts


  !% Communicate forces of ghost particles back.
  !% This is needed for rigid object (i.e., water) or to reduce the
  !% border size in BOPs.
  subroutine domaindecomposition_comm_reverse(this, at, error)
    implicit none

    type(DomainDecomposition), intent(inout)  :: this
    type(Atoms), intent(inout)                :: at
    integer, intent(out), optional :: error

    ! ---

    integer   :: i, d, list_off_r, list_off_l, n_recv_r, n_recv_l, cur

    ! ---

    INIT_ERROR(error)

    call print("DomainDecomposition : comm_ghosts", PRINT_NERD)

    call update_sendrecv_masks(this, at)

    call system_timer("domaindecomposition_comm_reverse")

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
             call copy_forces_to_send_ghosts(this, at, this%global_to_local(this%ghosts_r(list_off_r+i)), (i-1)*this%reverse_buffer_size, this%send_r)
          enddo

          do i = 1, this%n_ghosts_l(d)
             call copy_forces_to_send_ghosts(this, at, this%global_to_local(this%ghosts_l(list_off_l+i)), (i-1)*this%reverse_buffer_size, this%send_l)
          enddo

          call sendrecv(this%mpi, &
               this%send_r(1:this%reverse_buffer_size*this%n_ghosts_r(d)), &
               this%r(d), 0, &
               this%recv_l, this%l(d), 0, &
               n_recv_l, error)
          PASS_ERROR(error)

          call sendrecv(this%mpi, &
               this%send_l(1:this%reverse_buffer_size*this%n_ghosts_l(d)), &
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

    call system_timer("domaindecomposition_comm_reverse")

  endsubroutine domaindecomposition_comm_reverse


  !% Create a deep copy
  subroutine domaindecomposition_deepcopy(to, from)
    implicit none

    type(DomainDecomposition), intent(out)  :: to
    type(DomainDecomposition), intent(in)   :: from

    ! ---

    to = from
    
  endsubroutine domaindecomposition_deepcopy

  !
  ! Supplement to the Dictionary object
  !

  !% Convert a list of keys to a mask
  subroutine dictionary_keys_to_mask(this, keys, mask, s, error)
    implicit none

    type(Dictionary), intent(in)      :: this          !% Dictionary object
    type(Dictionary), intent(in)      :: keys          !% List of keys
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
       entry_i = lookup_entry_i(this, string(keys%keys(i)))
       if (entry_i == -1) then
!          RAISE_ERROR("Could not find key '" // keys%keys(i) // "'.", error)
          call print("DomainDecomposition : WARNING - Could not find key '" // keys%keys(i) // "', this property will not be communicated.", PRINT_ANAL)
       else

          select case(this%entries(entry_i)%type)

          case (T_REAL_A)
             s = s + sizeof(this%entries(entry_i)%r_a(1))

          case (T_INTEGER_A)
             s = s + sizeof(this%entries(entry_i)%i_a(1))

          case (T_LOGICAL_A)
             s = s + sizeof(this%entries(entry_i)%l_a(1))

          case (T_REAL_A2)
             s = s + sizeof(this%entries(entry_i)%r_a2(:, 1))

          case (T_INTEGER_A2)
             s = s + sizeof(this%entries(entry_i)%i_a2(:, 1))

          case default
             RAISE_ERROR("Don't know how to handle entry type " // this%entries(entry_i)%type // " (key '" // this%keys(entry_i) // "').", error)

          endselect

          mask(entry_i) = .true.

       endif
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

    integer  :: i, s

    ! ---

    INIT_ERROR(error)

    do i = 1, this%N
       if (mask(i)) then

          select case(this%entries(i)%type)

          case (T_REAL_A)
             s = sizeof(this%entries(i)%r_a(1))
             buffer(buffer_i+1:buffer_i+s) = &
                  transfer(this%entries(i)%r_a(data_i), buffer)
             buffer_i = buffer_i + s

          case (T_INTEGER_A)
             s = sizeof(this%entries(i)%i_a(1))
             buffer(buffer_i+1:buffer_i+s) = &
                  transfer(this%entries(i)%i_a(data_i), buffer)
             buffer_i = buffer_i + s

          case (T_LOGICAL_A)
             s = sizeof(this%entries(i)%l_a(1))
             buffer(buffer_i+1:buffer_i+s) = &
                  transfer(this%entries(i)%l_a(data_i), buffer)
             buffer_i = buffer_i + s

          case (T_REAL_A2)
             s = sizeof(this%entries(i)%r_a2(:, 1))
             buffer(buffer_i+1:buffer_i+s) = &
                  transfer(this%entries(i)%r_a2(:, data_i), buffer)
             buffer_i = buffer_i + s

          case (T_INTEGER_A2)
             s = sizeof(this%entries(i)%i_a2(:, 1))
             buffer(buffer_i+1:buffer_i+s) = &
                  transfer(this%entries(i)%i_a2(:, data_i), buffer)
             buffer_i = buffer_i + s

          case default
             RAISE_ERROR("Don't know how to handle entry type " // this%entries(i)%type // " (key '" // this%keys(i) // "').", error)

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

    integer  :: i, ndims, s

    ! ---

    INIT_ERROR(error)

    do i = 1, this%N
       if (mask(i)) then

          select case(this%entries(i)%type)

          case (T_REAL_A)
             s = sizeof(this%entries(i)%r_a(1))
             this%entries(i)%r_a(data_i) = &
                  transfer(buffer(buffer_i+1:buffer_i+s), 1.0_DP)
             buffer_i = buffer_i + s

          case (T_INTEGER_A)
             s = sizeof(this%entries(i)%i_a(1))
             this%entries(i)%i_a(data_i) = &
                  transfer(buffer(buffer_i+1:buffer_i+s), 1)
             buffer_i = buffer_i + s

          case (T_LOGICAL_A)
             this%entries(i)%l_a(data_i) = &
                  transfer(buffer(buffer_i+1:), .true.)
             buffer_i = buffer_i + sizeof(this%entries(i)%l_a(1))

          case (T_REAL_A2)
             s = sizeof(this%entries(i)%r_a2(:, 1))
             ndims = size(this%entries(i)%r_a2, 1)
             this%entries(i)%r_a2(:, data_i) = &
                  transfer(buffer(buffer_i+1:buffer_i+s), 1.0_DP, ndims)
             buffer_i = buffer_i + s

          case (T_INTEGER_A2)
             s = sizeof(this%entries(i)%i_a2(:, 1))
             ndims = size(this%entries(i)%i_a2, 1)
             this%entries(i)%i_a2(:, data_i) = &
                  transfer(buffer(buffer_i+1:buffer_i+s), 1.0_DP, ndims)
             buffer_i = buffer_i + s

          case default
             RAISE_ERROR("Don't know how to handle entry type " // this%entries(i)%type // " (key '" // this%keys(i) // "').", error)
             
          endselect

       endif
    enddo

  endsubroutine dictionary_unpack_buffer

endmodule DomainDecomposition_module
