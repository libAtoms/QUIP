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

!X
!X  Connection module
!X
!% A Connection object stores the connectivity information (i.e. 
!% list) for a specific Atoms object.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module Connection_module

  use error_module
  use system_module
  use units_module
  use periodictable_module
  use linearalgebra_module
  use dictionary_module
  use table_module
  use Atoms_types_module

  implicit none
  private

  real(dp), parameter :: CONNECT_LATTICE_TOL = 1e-8_dp

  public :: Connection

  public :: initialise
  interface initialise
     module procedure connection_initialise
  end interface initialise

  public :: finalise
  interface finalise
     module procedure connection_finalise
  end interface finalise

  public :: wipe
  interface wipe
     module procedure connection_wipe
  end interface wipe

  !% Print a verbose textual description of a Connection object to the default
  !% logger or to a specificied Inoutput object.
  public :: print
  interface print
     module procedure connection_print
  end interface print

  !% Overloaded assigment operators for Connection objects.
  public :: assignment(=)
  interface assignment(=)
     module procedure connection_assignment
  end interface assignment(=)

  public :: deepcopy
  interface deepcopy
     module procedure connection_assignment
  endinterface deepcopy

  public :: calc_connect
  interface calc_connect
     module procedure connection_calc_connect
  endinterface

  public :: calc_dists
  interface calc_dists
     module procedure connection_calc_dists
  end interface calc_dists

  public :: calc_connect_hysteretic
  interface calc_connect_hysteretic
     module procedure connection_calc_connect_hysteretic
  endinterface

  public :: n_neighbours
  interface n_neighbours
     module procedure connection_n_neighbours, connection_n_neighbours_with_dist
  end interface

  public :: n_neighbours_total
  interface n_neighbours_total
     module procedure connection_n_neighbours_total
  endinterface

  public :: is_min_image
  interface is_min_image
     module procedure connection_is_min_image
  end interface

  public :: neighbour
  interface neighbour
     module procedure connection_neighbour, connection_neighbour_minimal
  endinterface

  public :: neighbour_index
  interface neighbour_index
     module procedure connection_neighbour_index
  endinterface

  public :: neighbour_minimal
  interface neighbour_minimal
     module procedure connection_neighbour_minimal
  endinterface

  public :: add_bond, remove_bond, remove_bonds, cell_of_pos, connection_cells_initialise
  public :: connection_fill, divide_cell, fit_box_in_cell, get_min_max_images
  public :: max_cutoff, partition_atoms, cell_n
  public :: is_in_subregion

contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Connectivity procedures: Initialise and Finalise
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  !% Initialise a Connection object for a given number of atoms 'N'.
  !% If the optional Atoms argument is present then we calculate
  !% the atomic density to initialise the default lengths of the neighbour
  !% list for efficient memory usage.
   subroutine connection_initialise(this, N, Nbuffer, pos, lattice, g,  origin, extent, nn_guess, store_rij, fill)
    type(Connection),   intent(inout) :: this
    integer,            intent(in)    :: N    ! No. of atoms
    integer,            intent(in)    :: Nbuffer    ! Buffer size
    real(dp), optional, intent(in) :: pos(:,:), lattice(3,3), g(3,3)
    real(dp), optional, intent(in) :: origin(3), extent(3,3)
    integer, optional, intent(in) :: nn_guess
    logical, optional, intent(in) :: store_rij
    logical, optional, intent(in) :: fill

    logical :: do_fill

    do_fill = optional_default(.true., fill)

    ! If already initialised, destroy the existing data and start again
    if (this%initialised) call connection_finalise(this)

    if (do_fill) call connection_fill(this, N, Nbuffer, pos, lattice, g, origin, extent, nn_guess, store_rij)

    this%last_connect_cutoff = 0.0_dp
    this%last_connect_lattice(:,:) = 0.0_dp

  end subroutine connection_initialise

  subroutine connection_fill(this, N, Nbuffer, pos, lattice, g, origin, extent, nn_guess, store_rij, error)
    type(Connection),   intent(inout) :: this
    integer,            intent(in)    :: N    ! No. of atoms
    integer,            intent(in)    :: Nbuffer    ! Buffer size
    real(dp), optional, intent(in) :: pos(:,:), lattice(3,3), g(3,3)
    real(dp), optional, intent(in) :: origin(3), extent(3,3)
    integer, optional, intent(in) :: nn_guess
    logical, optional, intent(in) :: store_rij
    integer, intent(out), optional :: error     
    
    integer                           :: i, do_nn_guess
    real(dp)                          :: extent_inv(3,3), subregion_center(3)
    logical :: do_subregion, do_store_rij

    INIT_ERROR(error)
    do_nn_guess = optional_default(5, nn_guess)

    if (present(origin) .and. present(extent)) then
      if (.not.present(lattice) .or. .not.present(g)) then
         RAISE_ERROR("connection_fill got origin and extent, so trying to do subregion, but lattice or g are missing", error)
      end if
      do_subregion = .true.
      call matrix3x3_inverse(extent,extent_inv)
      subregion_center = origin + 0.5_dp*sum(extent,2)
    else
       do_subregion = .false.
    endif

    if (allocated(this%neighbour1)) then
       if (size(this%neighbour1, 1) < Nbuffer) deallocate(this%neighbour1)
    endif
    if (allocated(this%neighbour2)) then
       if (size(this%neighbour2, 1) < Nbuffer) deallocate(this%neighbour2)
    endif
    
    if (.not. allocated(this%neighbour1)) allocate(this%neighbour1(Nbuffer))
    if (.not. allocated(this%neighbour2)) allocate(this%neighbour2(Nbuffer))
    do i=1,Nbuffer
       if (do_subregion .and. i <= N) then
          if (.not. is_in_subregion(pos(:,i), subregion_center, lattice, g, extent_inv)) then
             if (associated(this%neighbour1(i)%t)) then
                call connection_remove_atom(this, i, error)
                PASS_ERROR(error)
                call finalise(this%neighbour1(i)%t)
                call finalise(this%neighbour2(i)%t)
                deallocate(this%neighbour1(i)%t)
                deallocate(this%neighbour2(i)%t)
             endif
             cycle
          endif
       endif
       if (.not. associated(this%neighbour1(i)%t)) then
          allocate(this%neighbour1(i)%t)
	  do_store_rij = optional_default(.false., store_rij)
	  if (do_store_rij) then
	     call allocate(this%neighbour1(i)%t,4,4, 0, 0, max(do_nn_guess, 1))
	  else
	     call allocate(this%neighbour1(i)%t,4,1, 0, 0, max(do_nn_guess, 1))
	  endif
          this%neighbour1(i)%t%increment = max(do_nn_guess/2, 1)

          allocate(this%neighbour2(i)%t)
          call allocate(this%neighbour2(i)%t,2,0, 0, 0, max(do_nn_guess, 1))
          this%neighbour2(i)%t%increment = max(do_nn_guess/2, 1)
       endif
    end do

    this%initialised = .true.

  end subroutine connection_fill

  function is_in_subregion(p, center, lattice, lattice_inv, extent_inv)
    real(dp), intent(in) :: p(3), center(3)
    real(dp), intent(in) :: lattice(3,3), lattice_inv(3,3), extent_inv(3,3)
    logical :: is_in_subregion

    real(dp) :: relative_p(3), extent_lattice_relative_p(3)

    relative_p = p - center
    call map_into_cell(relative_p, lattice, lattice_inv)

    extent_lattice_relative_p = extent_inv .mult. relative_p

    if (any(extent_lattice_relative_p < -0.5_dp) .or. any(extent_lattice_relative_p > 0.5_dp)) then
      is_in_subregion = .false.
    else
      is_in_subregion = .true.
    endif

  end function is_in_subregion


  !% OMIT
  subroutine connection_cells_initialise(this,cellsNa,cellsNb,cellsNc,Natoms)

    type(Connection),  intent(inout) :: this
    integer,           intent(in)    :: cellsNa,cellsNb,cellsNc ! No. cells in a,b,c directions
    integer, optional, intent(in)    :: Natoms                  ! Number of atoms
    integer                          :: i,j,k,av_atoms,stdev_atoms,Ncells

    if (this%cells_initialised) call connection_cells_finalise(this)

    !Set length and increment based on binomial statistics (if possible)
    Ncells = cellsNa * cellsNb * cellsNc
    if (present(Natoms)) then
       av_atoms = Natoms / Ncells
       stdev_atoms = int(sqrt(real(Natoms,dp) * (real(Ncells,dp) - 1.0_dp) / (real(Ncells,dp) * real(Ncells,dp))))
    else
       av_atoms = 100   !defaults if number of atoms is not given
       stdev_atoms = 10
    end if

    allocate(this%cell_heads(cellsNa, cellsNb, cellsNc))
    this%cell_heads = 0

    this%cellsNa = cellsNa
    this%cellsNb = cellsNb
    this%cellsNc = cellsNc

    this%cells_initialised = .true.

  end subroutine connection_cells_initialise


  !% Finalise this connection object
  subroutine connection_finalise(this)

    type(Connection), intent(inout) :: this
    integer                         :: i

    !do nothing if not initialised / already finalised
    if (.not.this%initialised) return

    if (allocated(this%neighbour1)) then
      do i=1,size(this%neighbour1)
	 if (associated(this%neighbour1(i)%t)) then
	   call finalise(this%neighbour1(i)%t)
	   deallocate(this%neighbour1(i)%t)
	 endif
      end do
    endif

    if (allocated(this%neighbour2)) then
      do i=1,size(this%neighbour2)
	 if (associated(this%neighbour2(i)%t)) then
	   call finalise(this%neighbour2(i)%t)
	   deallocate(this%neighbour2(i)%t)
	 endif
      end do
    endif


    if(allocated(this%neighbour1)) deallocate(this%neighbour1)
    if(allocated(this%neighbour2)) deallocate(this%neighbour2)

    if (allocated(this%is_min_image)) deallocate(this%is_min_image)

    if (allocated(this%last_connect_pos)) deallocate(this%last_connect_pos)

    call connection_cells_finalise(this)

    this%initialised = .false.

  end subroutine connection_finalise

  !% Wipe the contents of the connection tables, but keep the allocation
  subroutine connection_wipe(this)

    type(Connection), intent(inout) :: this
    integer                         :: i

    !do nothing if not initialised / already finalised
    if (.not.this%initialised) return

    do i=1,size(this%neighbour1)
       call wipe(this%neighbour1(i)%t)
    end do

    do i=1,size(this%neighbour2)
       call wipe(this%neighbour2(i)%t)
    end do

    call wipe_cells(this)

  end subroutine connection_wipe

  !% OMIT
  subroutine connection_cells_finalise(this)

    type(Connection), intent(inout) :: this
    integer                         :: i,j,k

    if (.not.this%cells_initialised) return

    deallocate(this%cell_heads)
    if (allocated(this%next_atom_in_cell)) then
       deallocate(this%next_atom_in_cell)
    endif

    this%cells_initialised = .false.

  end subroutine connection_cells_finalise

  subroutine connection_assignment(to,from)

    type(Connection), intent(inout) :: to
    type(Connection), intent(in)    :: from
    integer                         :: i,j,k

    call connection_finalise(to)

    if (.not.from%initialised) return


    call connection_initialise(to,size(from%neighbour1),size(from%neighbour1))


    ! Use Append to append the source tables to our (empty) destination tables
    do i=1,size(from%neighbour1)
       call append(to%neighbour1(i)%t,from%neighbour1(i)%t)
    end do

    do i=1,size(from%neighbour2)
       call append(to%neighbour2(i)%t,from%neighbour2(i)%t)
    end do


    !If cell data is present then copy that too
    if (from%cells_initialised) then
       call connection_cells_initialise(to,from%cellsNa,from%cellsNb,from%cellsNc)
       to%cell_heads = from%cell_heads
       if (allocated(to%next_atom_in_cell) .and. allocated(from%next_atom_in_cell)) then
          to%next_atom_in_cell = from%next_atom_in_cell
       endif
    end if

  end subroutine connection_assignment

  !% OMIT
  subroutine wipe_cells(this)

    type(Connection) :: this
    integer          :: i,j,k

    if (this%cells_initialised) then
       this%cell_heads = 0
    end if

  end subroutine wipe_cells

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Connectivity procedures
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% Test if atom $i$ is a neighbour of atom $j$ and update 'this%connect' as necessary.
  !% Called by 'calc_connect'. The 'shift' vector is added to the position of the $j$ atom
  !% to get the correct image position.
  subroutine test_form_bond(this,cutoff, use_uniform_cutoff, Z, pos, lattice, i,j, shift, check_for_dup, error)

    type(Connection), intent(inout) :: this
    real(dp), intent(in) :: cutoff
    logical, intent(in) :: use_uniform_cutoff
    integer, intent(in) :: Z(:)
    real(dp), intent(in) :: pos(:,:), lattice(3,3)
    integer,     intent(in)    :: i,j
    integer,     intent(in)    :: shift(3)
    logical, intent(in) :: check_for_dup
    integer, intent(out), optional :: error

    integer                    :: index, m, k
    real(dp)                   :: d, dd(3)
    real(dp)                   :: use_cutoff

    INIT_ERROR(error)
    if (i > j) return

    if (.not. associated(this%neighbour1(i)%t) .or. .not. associated(this%neighbour1(j)%t)) return

#ifdef DEBUG
    if(current_verbosity() >= PRINT_ANAL) then
       call print('Entering test_form_bond, i = '//i//' j = '//j, PRINT_ANAL)
       call print('use_uniform_cutoff = '//use_uniform_cutoff, PRINT_ANAL)
    end if
#endif

    !Determine what cutoff distance to use
    if (use_uniform_cutoff) then
       use_cutoff = cutoff
    else
       use_cutoff = bond_length(Z(i),Z(j)) * cutoff
    end if

    !d = norm(pos(:,j)+(lattice .mult. shift) - pos(:,i))
    !OPTIM
    dd = pos(:,j) - pos(:,i)
    do m=1,3
       forall(k=1:3) dd(k) = dd(k) + lattice(k,m) * shift(m)
    end do

    if (.not. check_for_dup) then
       do m=1,3
          if (dd(m) > use_cutoff) return
       end do
    end if

    d = sqrt(dd(1)*dd(1) + dd(2)*dd(2) + dd(3)*dd(3))

    if (check_for_dup) then
      index = find(this%neighbour1(i)%t, (/ j, shift /)) 
      if (index /= 0) then ! bond is already in table
#ifdef DEBUG
	if (current_verbosity() >= PRINT_ANAL) call print('test_form_bond had check_for_dup=T, found bond already in table', PRINT_ANAL)
#endif
	this%neighbour1(i)%t%real(1,index) = d
	if (size(this%neighbour1(i)%t%real,1) == 4) then ! store_rij was set
	   this%neighbour1(i)%t%real(2:4,index) = pos(1:3,j) + (lattice .mult. shift) - pos(1:3,i)
	endif
	return
      endif
    endif

#ifdef DEBUG
    if(current_verbosity() >= PRINT_ANAL)  call print('d = '//d, PRINT_ANAL)
#endif

    if (d < use_cutoff) then
       call add_bond(this, pos, lattice, i, j, shift, d, error=error)
       PASS_ERROR(error)
    end if

#ifdef DEBUG
    if(current_verbosity() >= PRINT_ANAL) call print('Leaving test_form_bond', PRINT_ANAL)
#endif

  end subroutine test_form_bond

  !% Test if atom $i$ is no longer neighbour of atom $j$ with shift $s$, and update 'this%connect' as necessary.
  !% Called by 'calc_connect'. The 'shift' vector is added to the position of the $j$ atom
  !% to get the correct image position.
  function test_break_bond(this,cutoff_break, use_uniform_cutoff, Z, pos, lattice, i,j, shift, error)
    type(Connection), intent(inout) :: this
    real(dp), intent(in) :: cutoff_break
    logical, intent(in) :: use_uniform_cutoff
    integer, intent(in) :: Z(:)
    real(dp), intent(in) :: pos(:,:), lattice(3,3)
    integer,     intent(in)    :: i,j
    integer, intent(in)        :: shift(3)
    logical test_break_bond
    integer, intent(out), optional :: error

    real(dp)                   :: d
    real(dp)                   :: cutoff

    INIT_ERROR(error)
    test_break_bond = .false.
    if (i > j) return

    if (.not. associated(this%neighbour1(i)%t) .or. .not. associated(this%neighbour1(j)%t)) return

#ifdef DEBUG
    if(current_verbosity() >= PRINT_ANAL) then
       call print('Entering test_break_bond, i = '//i//' j = '//j, PRINT_ANAL)
       call print('use_uniform_cutoff = '//use_uniform_cutoff // " cutoff_break = "// cutoff_break, PRINT_ANAL)
    end if
#endif

    !Determine what cutoff distance to use
    if (use_uniform_cutoff) then
       cutoff = cutoff_break
    else
       cutoff = bond_length(Z(i),Z(j)) * cutoff_break
    end if

    d = norm(pos(:,j)+(lattice .mult. shift) - pos(:,i))
#ifdef DEBUG
    if(current_verbosity() >= PRINT_ANAL)  call print('d = '//d//' cutoff = '//cutoff//' i = '//i//' j = '//j, PRINT_ANAL)
#endif
    if (d > cutoff) then
#ifdef DEBUG
       if(current_verbosity() >= PRINT_ANAL) call print('removing bond from tables', PRINT_ANAL)
#endif
       call remove_bond(this, i, j, shift, error)
       PASS_ERROR(error)
       test_break_bond = .true.
    end if

#ifdef DEBUG
    if(current_verbosity() >= PRINT_ANAL) call print('Leaving test_break_bond', PRINT_ANAL)
#endif

  end function test_break_bond

  subroutine add_bond(this, pos, lattice, i, j, shift, d, dv, error)
    type(Connection), intent(inout) :: this
    real(dp), intent(in) :: pos(:,:), lattice(3,3)
    integer,     intent(in)    :: i,j
    integer,     intent(in)    :: shift(3)
    real(dp), intent(in) :: d
    real(dp), intent(in), optional :: dv(3)
    integer, intent(out), optional :: error

    real(dp) :: ddv(3)
    integer :: ii, jj, index

    INIT_ERROR(error)
    if (.not.this%initialised) then
       RAISE_ERROR("add_bond called on uninitialized connection", error)
    endif

    if (.not. associated(this%neighbour1(i)%t) .or. .not. associated(this%neighbour1(j)%t)) then
      RAISE_ERROR("tried to add_bond for atoms i " // i // " j " // j // " which have associated(neighbour1()%t "//associated(this%neighbour1(i)%t) // " " // associated(this%neighbour1(j)%t) // " one of which is false", error)
    endif

    if (i > j) then
      ii = j
      jj = i
    else
      ii = i
      jj = j
    endif

    if (present(dv)) then
      ddv = dv*sign(1,j-i)
    else ! dv not present
      if (size(this%neighbour1(i)%t%real,1) == 4) then
	 ddv = pos(:,jj) + (lattice .mult. shift) - pos(:,ii)
      endif
    endif

    ! Add full details to neighbour1 for smaller of i and j
    if (size(this%neighbour1(i)%t%real,1) == 4) then ! store_rij was set
       call append(this%neighbour1(ii)%t, (/jj, shift /), (/ d, ddv /))
    else
       call append(this%neighbour1(ii)%t, (/jj, shift /), (/ d /))
    endif
    if(ii .ne. jj) then		
       index = this%neighbour1(min(ii,jj))%t%N
       ! Put a reference to this in neighbour2 for larger of i and j
       call append(this%neighbour2(jj)%t, (/ ii, index/))
    end if

  end subroutine add_bond


  subroutine remove_bond(this, i, j, shift, error)
    type(Connection), intent(inout) :: this
    integer,     intent(in)    :: i,j
    integer,     intent(in), optional    :: shift(3)
    integer, intent(out), optional :: error     
    

    integer :: ii, jj, iii, jjj, jjjj, r_index, n_removed, my_shift(3)

    INIT_ERROR(error)
    if (.not. associated(this%neighbour1(i)%t) .or. .not. associated(this%neighbour1(j)%t)) then
       RAISE_ERROR("tried to remove_bond for atoms i " // i // " j " // j // " which have associated(neighbour1()%t " //associated(this%neighbour1(i)%t) // " " // associated(this%neighbour1(j)%t) // " one of which is false", error)
    endif

    if (i > j) then
      ii = j
      jj = i
      if (present(shift)) my_shift = -shift
    else
      ii = i
      jj = j
      if (present(shift)) my_shift = shift
    endif
    ! now ii <= jj

    r_index = 1
    n_removed = 0
    do while (r_index /= 0)
      ! remove entry from neighbour1(ii)
      if (present(shift)) then
	r_index = find(this%neighbour1(ii)%t, (/ jj, my_shift /) )
      else
	r_index = find(this%neighbour1(ii)%t, (/ jj, 0, 0, 0 /), (/ .true., .false., .false., .false./) )
      endif
      if (r_index == 0) then
	if (n_removed == 0) then
	  if (present(shift)) then
	    call print("WARNING: remove bond called for i " // i // " j " // j // " shift " // shift // &
		       " couldn't find a bond to remove", PRINT_ALWAYS)
	  else
	    call print("WARNING: remove bond called for i " // i // " j " // j // &
		       " couldn't find a bond to remove", PRINT_ALWAYS)
	  endif
	endif
      else ! r_index /= 0
	n_removed = n_removed + 1

	call delete(this%neighbour1(ii)%t, r_index, keep_order = .true.)
	! remove entry from neighbour2(jj)
	if (ii /= jj) then
	  call delete(this%neighbour2(jj)%t, (/ ii, r_index /), keep_order = .true.)
	endif
	! renumber other neighbour2 entries
	do iii=r_index, this%neighbour1(ii)%t%N
	  ! jjj is another neighbour of ii
	  jjj = this%neighbour1(ii)%t%int(1,iii)
	  if (jjj > ii) then
	    ! jjjj is r_index in jjj's neighbour2 of pointer back to ii's neighbour1
	    jjjj = find(this%neighbour2(jjj)%t, (/ ii, iii+1 /) )
	    if (jjjj /= 0) then
	      ! decrement reference in jjj's neighbour2 table
	      this%neighbour2(jjj)%t%int(:,jjjj) = (/ ii, iii /)
	    else
	      RAISE_ERROR("remove_bond: Couldn't find neighbor to fix neighbour2 of", error)
	    endif
	  endif
	end do
      end if ! r_index == 0
    end do ! while r_index /= 0

  end subroutine remove_bond

  
   !% Remove all bonds listed in the Table `bonds` from connectivity 
   subroutine remove_bonds(this, at, bonds, error)
     type(Connection), intent(inout) :: this
     type(Atoms), intent(in) :: at
     type(Table), intent(in) :: bonds
     integer, intent(out), optional :: error

     logical :: bond_exists
     integer :: i, j, ji
     integer :: shift(3)

     INIT_ERROR(error)

     do i=1,bonds%N
        bond_exists = .false.
        do ji=1, n_neighbours(this,bonds%int(1,i))
           j = neighbour(this, at, bonds%int(1,i), ji, shift=shift)
           if (j == bonds%int(2,i)) then
              bond_exists = .true.
              exit
           endif
        end do
        if (bond_exists) then
           call remove_bond(this, bonds%int(1,i), bonds%int(2,i), shift, error=error)
           PASS_ERROR(error)
        endif
     end do
   end subroutine remove_bonds



  !%  As for 'calc_connect', but perform the connectivity update
  !%  hystertically: atoms must come within 'cutoff' to be considered
  !%  neighbours, and then will remain connect until them move apart
  !%  further than 'cutoff_break'.
  subroutine connection_calc_connect_hysteretic(this, at, origin, extent, own_neighbour, store_is_min_image, error)
    type(Connection), intent(inout)  :: this
    type(Atoms), intent(inout) :: at
    real(dp), optional :: origin(3), extent(3,3)
    logical, optional, intent(in) :: own_neighbour, store_is_min_image
    integer, intent(out), optional :: error

    integer  :: cellsNa,cellsNb,cellsNc
    integer  :: i,j,k,i2,j2,k2,i3,j3,k3,i4,j4,k4,atom1,atom2
!    integer  :: n1, n2
    integer  :: cell_image_Na, cell_image_Nb, cell_image_Nc
    integer  :: min_cell_image_Na
    integer  :: max_cell_image_Na, min_cell_image_Nb, max_cell_image_Nb
    integer  :: min_cell_image_Nc, max_cell_image_Nc
    real(dp)  :: cutoff
    integer :: ji, s_ij(3), nn_guess
    logical my_own_neighbour, my_store_is_min_image
    logical :: change_i, change_j, change_k, broken
    integer, pointer :: map_shift(:,:)

    INIT_ERROR(error)

    this%N = at%N

    my_own_neighbour = optional_default(.false., own_neighbour)
    my_store_is_min_image = optional_default(.true., store_is_min_image)

    if (at%cutoff < 0.0_dp .or. at%cutoff_break < 0.0_dp) then
       RAISE_ERROR('calc_connect: Negative cutoff radius ' // at%cutoff // ' ' // at%cutoff_break, error)
    end if

    if (at%cutoff > at%cutoff_break) then
       RAISE_ERROR('calc_connect: Negative hysteresis cutoff radius formation ' // at%cutoff // ' > breaking ' // at%cutoff_break, error)
    end if

    if ((at%cutoff .feq. 0.0_dp) .or. (at%cutoff_break .feq. 0.0_dp)) then
      call wipe(this)
      return
    endif

    !Calculate the cutoff value we should use in dividing up the simulation cell
    if (at%use_uniform_cutoff) then
       !The cutoff has been specified by the user, so use that value
       cutoff = at%cutoff
    else
       !Otherwise, find the maximum covalent radius of all the atoms in the structure, double it
       !and multiple by cutoff. At makes sure we can deal with the worst case scenario of big atom 
       !bonded to big atom
       cutoff = 0.0_dp
       do i = 1, at%N
          if (ElementCovRad(at%Z(i)) > cutoff) cutoff = ElementCovRad(at%Z(i))
       end do
       cutoff = (2.0_dp * cutoff) * at%cutoff
    end if

    call print("calc_connect: cutoff calc_connect " // cutoff, PRINT_NERD)

    if (present(origin) .and. present(extent)) then
      cellsNa = 1
      cellsNb = 1
      cellsNc = 1
    else
      call divide_cell(at%lattice, cutoff, cellsNa, cellsNb, cellsNc)
    endif

    call print("calc_connect: cells_N[abc] " // cellsNa // " " // cellsNb // " " // cellsNc, PRINT_NERD)

    ! If the lattice has changed, then the cells need de/reallocating
    if ((cellsNa /= this%cellsNa) .or. &
         (cellsNb /= this%cellsNb) .or. &
         (cellsNc /= this%cellsNc)) call connection_cells_finalise(this)

    ! figure out how many unit cell images we will need to loop over in each direction
    call fit_box_in_cell(cutoff, cutoff, cutoff, at%lattice, cell_image_Na, cell_image_Nb, cell_image_Nc)
    ! cell_image_N{a,b,c} apply to a box of side 2*cutoff. Since we loop from -cell_image_N to
    ! +cell_image_N we can reduce them as follows

    cell_image_Na = max(1,(cell_image_Na+1)/2)
    cell_image_Nb = max(1,(cell_image_Nb+1)/2)
    cell_image_Nc = max(1,(cell_image_Nc+1)/2)

    ! Estimate number of neighbours of each atom. Factor of 1/2 assumes
    ! half will go in neighbour1, half in neighbour2.
    nn_guess = int(0.5_dp*4.0_dp/3.0_dp*PI*cutoff**3*at%N/cell_volume(at%lattice)*cell_image_na*cell_image_nb*cell_image_nc)

    call print('calc_connect: image cells '//cell_image_Na//'x'//cell_image_Nb//'x'//cell_image_Nc, PRINT_NERD)

    ! Allocate space for the connection object if needed
    if (present(origin) .and. present(extent)) then
      if (.not.this%initialised) then
	 call connection_initialise(this, at%N, at%Nbuffer, at%pos, at%lattice, at%g, origin, extent, nn_guess)
      else
	 call connection_fill(this, at%N, at%Nbuffer, at%pos, at%lattice, at%g, origin, extent, nn_guess)
      end if
    else
      if (.not.this%initialised) then
	 call connection_initialise(this, at%N, at%Nbuffer, at%pos, at%lattice, at%g, nn_guess=nn_guess)
      else
	 call connection_fill(this, at%N, at%Nbuffer, at%pos, at%lattice, at%g, nn_guess=nn_guess)
      end if
    endif

    if (.not.this%cells_initialised) then
      call connection_cells_initialise(this, cellsNa, cellsNb, cellsNc,at%N)
    endif

    ! Partition the atoms into cells
    call partition_atoms(this, at, error=error)
    PASS_ERROR(error)
    if (.not. assign_pointer(at, 'map_shift', map_shift)) then
       RAISE_ERROR("calc_connect impossibly failed to assign map_shift pointer", error)
    end if

    ! look for bonds that have been broken, and remove them
    do i=1, at%N
      ji = 1
      do
	if (ji > n_neighbours(this, i)) exit
	j = connection_neighbour(this, at, i, ji, shift = s_ij)
        broken = test_break_bond(this, at%cutoff_break, at%use_uniform_cutoff, &
             at%Z, at%pos, at%lattice, i, j, s_ij, error)
        PASS_ERROR(error)
	if (.not. broken) then
	  ji = ji + 1 ! we didn't break at bond, so go to next one
	              ! if we did break a bond, ji now points to a different bond, so don't increment it
	endif
      end do
!      do ji=1, atoms_n_neighbours(at, i, alt_connect=this)
!	j = atoms_neighbour(at, i, ji, shift = s_ij, alt_connect=this)
!	call test_break_bond(this, at%cutoff_break, at%use_uniform_cutoff, &
!	  at%Z, at%pos, at%lattice, i, j, s_ij)
!      end do
    end do

    ! Here is the main loop:
    ! Go through each cell and update the connectivity between atoms in at cell and neighbouring cells
    ! N.B. test_form_bond updates both atoms i and j, so only update if i <= j to avoid doubling processing

    ! defaults for cellsNx = 1
    k3 = 1; k4 = 1; j3 = 1; j4 = 1; i3 = 1; i4 = 1
    ! Loop over all cells
    do k = 1, cellsNc
       change_k = .true.
       do j = 1, cellsNb
	  change_j = .true.
          do i = 1, cellsNa
	     change_i = .true.
	     call get_min_max_images(at%is_periodic, cellsNa, cellsNb, cellsNc, &
	        cell_image_Na, cell_image_Nb, cell_image_Nc, i, j, k, change_i, change_j, change_k, &
	        min_cell_image_Na, max_cell_image_Na, min_cell_image_Nb, max_cell_image_Nb, min_cell_image_Nc, max_cell_image_Nc)
	     change_i = .false.
	     change_j = .false.
	     change_k = .false.

             !Loop over atoms in cell(i,j,k)
             atom1 = this%cell_heads(i, j, k)
             atom1_loop: do while (atom1 > 0)

                ! Loop over neighbouring cells, applying PBC
                do k2 = -cell_image_Nc, +cell_image_Nc

                   ! the stored cell we are in 
                   if(cellsNc > 1) k3 = mod(k+k2-1+cellsNc,cellsNc)+1 

                   ! the shift we need to get to the cell image
                   k4 = (k+k2-k3)/cellsNc

                   do j2 = -cell_image_Nb, +cell_image_Nb
                      ! the stored cell we are in                 
                      if(cellsNb > 1) j3 = mod(j+j2-1+cellsNb,cellsNb)+1 

                      ! the shift we need to get to the cell image
                      j4 = (j+j2-j3)/cellsNb

                      do i2 = -cell_image_Na, +cell_image_Na
                         ! the stored cell we are in                 
                         if(cellsNa > 1) i3 = mod(i+i2-1+cellsNa,cellsNa)+1 

                         ! the shift we need to get to the cell image
                         i4 = (i+i2-i3)/cellsNa

                         ! The cell we are currently testing atom1 against is cell(i3,j3,k3)
                         ! with shift (i4,j4,k4)
                         ! loop over it's atoms and test connectivity if atom1 < atom2

                         atom2 = this%cell_heads(i3,j3,k3)
                         atom2_loop: do while (atom2 > 0)

                            ! omit atom2 < atom1
                            if (atom1 > atom2) then
                               atom2 = this%next_atom_in_cell(atom2)
                               cycle atom2_loop
                            endif
                            ! omit self in the same cell without shift
                            if (.not. my_own_neighbour .and. (atom1 == atom2 .and. & 
                                 (i4==0 .and. j4==0 .and. k4==0) .and. &
                                 (i==i3 .and. j==j3 .and. k==k3))) then
                               atom2 = this%next_atom_in_cell(atom2)
                               cycle atom2_loop
                            endif

                            call test_form_bond(this, at%cutoff, at%use_uniform_cutoff, &
                                 at%Z, at%pos, at%lattice, atom1,atom2, &
				 (/i4-map_shift(1,atom1)+map_shift(1,atom2),j4-map_shift(2,atom1)+map_shift(2,atom2),k4-map_shift(3,atom1)+map_shift(3,atom2)/), &
				 .true., error)
                            PASS_ERROR(error)

                            atom2 = this%next_atom_in_cell(atom2)
                         end do atom2_loop ! atom2

                      end do ! i2
		   end do ! j2
                end do ! k2

                atom1 = this%next_atom_in_cell(atom1)
             end do atom1_loop ! atom1

          end do ! i
       end do ! j
    end do ! k

    if (my_store_is_min_image) then
       if (allocated(this%is_min_image)) deallocate(this%is_min_image)
       allocate(this%is_min_image(at%n))
       do i=1,at%n
          if (associated(this%neighbour1(i)%t)) then
             this%is_min_image(i) = is_min_image(this, i)
          else
             this%is_min_image(i) = .false.
          end if
       end do
    end if

  end subroutine connection_calc_connect_hysteretic

  subroutine connection_remove_atom(this, i, error)
    type(Connection), intent(inout) :: this
    integer, intent(in) :: i
    integer, intent(out), optional :: error     
    
    integer :: ji, j, jj, s_ij(3), n_entries

    INIT_ERROR(error)
    n_entries=this%neighbour1(i)%t%N
    do ji=n_entries, 1, -1
      j = this%neighbour1(i)%t%int(1,ji)
      s_ij = this%neighbour1(i)%t%int(2:4,ji)
      call remove_bond(this, i, j, s_ij, error)
      PASS_ERROR(error)
    end do

    n_entries=this%neighbour2(i)%t%N
    do ji=n_entries, 1, -1
      j = this%neighbour2(i)%t%int(1,ji)
      jj = this%neighbour2(i)%t%int(2,ji)
      s_ij = this%neighbour1(j)%t%int(2:4,jj)
      call remove_bond(this, j, i, s_ij, error)
      PASS_ERROR(error)
    end do

  end subroutine connection_remove_atom


  !% Fast $O(N)$ connectivity calculation routine. It divides the unit
  !% cell into similarly shaped subcells, of sufficient size that
  !% sphere of radius 'cutoff' is contained in a subcell, at least in
  !% the directions in which the unit cell is big enough. For very
  !% small unit cells, there is only one subcell, so the routine is
  !% equivalent to the standard $O(N^2)$ method.  If 'cutoff_skin' is
  !% present, effective cutoff is increased by this amount, and full
  !% recalculation of connectivity is only done when any atom has
  !% moved more than 0.5*cutoff_skin.
  subroutine connection_calc_connect(this, at, own_neighbour, store_is_min_image, skip_zero_zero_bonds, store_n_neighb, cutoff_skin, max_pos_change, did_rebuild, error)
    type(Connection), intent(inout)  :: this
    type(Atoms), intent(inout) :: at
    logical, optional, intent(in) :: own_neighbour, store_is_min_image, skip_zero_zero_bonds, store_n_neighb
    real(dp), intent(in), optional :: cutoff_skin
    real(dp), intent(out), optional :: max_pos_change
    logical, intent(out), optional :: did_rebuild
    integer, intent(out), optional :: error

    integer  :: cellsNa,cellsNb,cellsNc
    integer  :: i,j,k,i2,j2,k2,i3,j3,k3,i4,j4,k4,atom1,atom2
!    integer  :: n1,n2
    integer  :: cell_image_Na, cell_image_Nb, cell_image_Nc, nn_guess, n_occ
    integer  :: min_cell_image_Na, max_cell_image_Na, min_cell_image_Nb
    integer  :: max_cell_image_Nb, min_cell_image_Nc, max_cell_image_Nc
    real(dp) :: cutoff, density, volume_per_cell, pos_change, my_max_pos_change
    logical my_own_neighbour, my_store_is_min_image, my_skip_zero_zero_bonds, my_store_n_neighb, do_fill
    logical :: change_i, change_j, change_k
    integer, pointer :: map_shift(:,:), n_neighb(:)


    INIT_ERROR(error)

    call system_timer('calc_connect')

    this%N = at%N

    my_own_neighbour = optional_default(.false., own_neighbour)
    my_store_is_min_image = optional_default(.true., store_is_min_image)
    my_skip_zero_zero_bonds = optional_default(.false., skip_zero_zero_bonds)
    my_store_n_neighb = optional_default(.true., store_n_neighb)

    if (at%cutoff < 0.0_dp .or. at%cutoff_break < 0.0_dp) then
       RAISE_ERROR('calc_connect: Negative cutoff radius ' // at%cutoff // ' ' // at%cutoff_break, error)
    end if

    if ((at%cutoff .feq. 0.0_dp) .or. (at%cutoff_break .feq. 0.0_dp)) then
      call wipe(this)
      return
    endif
    !Calculate the cutoff value we should use in dividing up the simulation cell
    if (at%use_uniform_cutoff) then
       !The cutoff has been specified by the user, so use that value
       cutoff = at%cutoff
    else
       !Otherwise, find the maximum covalent radius of all the atoms in the structure, double it
       !and multiple by cutoff. At makes sure we can deal with the worst case scenario of big atom 
       !bonded to big atom
       cutoff = 0.0_dp
       do i = 1, at%N
          if (ElementCovRad(at%Z(i)) > cutoff) cutoff = ElementCovRad(at%Z(i))
       end do
       cutoff = (2.0_dp * cutoff) * at%cutoff
    end if

    if (present(cutoff_skin)) then
       if (cutoff_skin .fne. 0.0_dp) then
          call print('calc_connect: increasing cutoff from '//cutoff//' by cutoff_skin='//cutoff_skin ,PRINT_VERBOSE)
          cutoff = cutoff + cutoff_skin

          if (.not. allocated(this%last_connect_pos) .or. &
               (cutoff > this%last_connect_cutoff) .or. &
               (size(this%last_connect_pos, 2) /= at%n) .or. &
               (at%lattice .fne. this%last_connect_lattice)) then
             call print('calc_connect: forcing a rebuild: either first time, atom number mismatch or lattice mismatch', PRINT_VERBOSE)
             call print('calc_connect: maxval(abs(at%lattice - this%last_connect_lattice)) = '//(maxval(abs(at%lattice - this%last_connect_lattice))), PRINT_VERBOSE)
             if (allocated(this%last_connect_pos)) deallocate(this%last_connect_pos)
             allocate(this%last_connect_pos(3, at%n))
             my_max_pos_change = huge(1.0_dp)
          else
             ! FIXME 1. we should also take into account changes in lattice here - for
             !          now we force a reconnect whenever lattice changes.
             !       2. it may be possible to further speed up calculation of delta_pos by 
             !          not calling distance_min_image() every time
             my_max_pos_change = 0.0_dp
             do i=1, at%N
                pos_change = distance_min_image(at, i, this%last_connect_pos(:, i))
                if (pos_change > my_max_pos_change) my_max_pos_change = pos_change
             end do
          end if

          if (present(max_pos_change)) max_pos_change = my_max_pos_change

          if (my_max_pos_change < 0.5_dp*cutoff_skin) then
             call print('calc_connect: max pos change '//my_max_pos_change//' < 0.5*cutoff_skin, doing a calc_dists() only', PRINT_VERBOSE)
             call calc_dists(this, at)
             call system_timer('calc_connect')
             if (present(did_rebuild)) did_rebuild = .false.
             return
          end if

          ! We need to do a full recalculation of connectivity. Store the current pos and lattice. 
          call print('calc_connect: max pos change '//my_max_pos_change//' >= 0.5*cutoff_skin, doing a full rebuild', PRINT_VERBOSE)
          if (present(did_rebuild)) did_rebuild = .true.
          this%last_connect_pos(:,:) = at%pos
          this%last_connect_lattice(:,:) = at%lattice
          this%last_connect_cutoff = cutoff
       end if
    end if

    call print("calc_connect: cutoff calc_connect " // cutoff, PRINT_VERBOSE)

    call divide_cell(at%lattice, cutoff, cellsNa, cellsNb, cellsNc)

    call print("calc_connect: cells_N[abc] " // cellsNa // " " // cellsNb // " " // cellsNc, PRINT_VERBOSE)

    ! If the lattice has changed, then the cells need de/reallocating
    if ((cellsNa /= this%cellsNa) .or. &
         (cellsNb /= this%cellsNb) .or. &
         (cellsNc /= this%cellsNc)) call connection_finalise(this)

    ! figure out how many unit cell images we will need to loop over in each direction
    call fit_box_in_cell(cutoff, cutoff, cutoff, at%lattice, cell_image_Na, cell_image_Nb, cell_image_Nc)
    ! cell_image_N{a,b,c} apply to a box of side 2*cutoff. Since we loop from -cell_image_N to
    ! +cell_image_N we can reduce them as follows

    cell_image_Na = max(1,(cell_image_Na+1)/2)
    cell_image_Nb = max(1,(cell_image_Nb+1)/2)
    cell_image_Nc = max(1,(cell_image_Nc+1)/2)

    call print('calc_connect: image cells '//cell_image_Na//'x'//cell_image_Nb//'x'//cell_image_Nc, PRINT_VERBOSE)

    ! Allocate space for the connection object if needed
    if (.not.this%initialised) then
       call connection_initialise(this, at%N, at%Nbuffer, at%pos, at%lattice, at%g, fill=.false.)
       do_fill = .true.
    else
       ! otherwise just wipe the connection table
       call wipe(this)
       do_fill = .false.
    end if

    if (.not.this%cells_initialised) &
         call connection_cells_initialise(this, cellsNa, cellsNb, cellsNc,at%N)

    ! Partition the atoms into cells
    call partition_atoms(this, at, error=error)
    PASS_ERROR(error)
    if (.not. assign_pointer(at, 'map_shift', map_shift)) then
       RAISE_ERROR("calc_connect impossibly failed to assign map_shift pointer", error)
    end if

    if (do_fill) then
       volume_per_cell = cell_volume(at%lattice)/real(cellsNa*cellsNb*cellsNc,dp)
 
       ! Count the occupied cells so vacuum does not contribute to average number density
       n_occ = 0
       do k=1,cellsNc
          do j=1,cellsNb
             do i=1,cellsNa
                if (this%cell_heads(i, j, k) > 0) &
                     n_occ = n_occ + 1
             end do
          end do
       end do
       density = at%n/(n_occ*volume_per_cell)

       ! Sphere of radius "cutoff", assume roughly half neighbours in neighbour1 and half in neighbour2
       nn_guess = int(4.0_dp/3.0_dp*PI*cutoff**3*density)/2

       call print('calc_connect: occupied cells '//n_occ//'/'//(cellsNa*cellsNb*cellsNc)//' = '//(n_occ/real(cellsNa*cellsNb*cellsNc,dp)), PRINT_VERBOSE)
       call print('calc_connect: estimated number of neighbours per atom = '//nn_guess, PRINT_VERBOSE)

       call connection_fill(this, at%N, at%Nbuffer, at%pos, at%lattice, at%g, nn_guess=nn_guess)
    end if

    ! Here is the main loop:
    ! Go through each cell and update the connectivity between atoms in at cell and neighbouring cells
    ! N.B. test_form_bond updates both atoms i and j, so only update if i <= j to avoid doubling processing

    ! defaults for cellsNx = 1
    k3 = 1; k4 = 1; j3 = 1; j4 = 1; i3 = 1; i4 = 1

    ! loop over all cells i,j,k, and all atoms in each cell n1
    do k = 1, cellsNc
       change_k = .true.
       do j = 1, cellsNb
	  change_j = .true.
	  do i = 1, cellsNa
	     change_i = .true.
	     call get_min_max_images(at%is_periodic, cellsNa, cellsNb, cellsNc, &
	        cell_image_Na, cell_image_Nb, cell_image_Nc, i, j, k, change_i, change_j, change_k, &
	        min_cell_image_Na, max_cell_image_Na, min_cell_image_Nb, max_cell_image_Nb, min_cell_image_Nc, max_cell_image_Nc)
	     change_i = .false.
	     change_j = .false.
	     change_k = .false.
             atom1 = this%cell_heads(i, j, k)
             atom1_loop: do while (atom1 > 0)

		! Loop over neighbouring cells, applying PBC

		do k2 = min_cell_image_Nc, max_cell_image_Nc

		   ! the stored cell we are in 
		   if(cellsNc > 1) k3 = mod(k+k2-1+cellsNc,cellsNc)+1 

		   ! the shift we need to get to the cell image
		   k4 = (k+k2-k3)/cellsNc

		   do j2 = min_cell_image_Nb, max_cell_image_Nb
		      ! the stored cell we are in                 
		      if(cellsNb > 1) j3 = mod(j+j2-1+cellsNb,cellsNb)+1 

		      ! the shift we need to get to the cell image
		      j4 = (j+j2-j3)/cellsNb

		      do i2 = min_cell_image_Na, max_cell_image_Na
			 ! the stored cell we are in                 
			 if(cellsNa > 1) i3 = mod(i+i2-1+cellsNa,cellsNa)+1 

			 ! the shift we need to get to the cell image
			 i4 = (i+i2-i3)/cellsNa

			 ! The cell we are currently testing atom1 against is cell(i3,j3,k3)
			 ! with shift (i4,j4,k4)
			 ! loop over it's atoms and test connectivity if atom1 < atom2

                         atom2 = this%cell_heads(i3, j3, k3)
                         atom2_loop: do while (atom2 > 0)

			    ! omit atom2 < atom1
			    if (atom1 > atom2) then 
                               atom2 = this%next_atom_in_cell(atom2)
                               cycle atom2_loop
                            endif
                            
                            if (my_skip_zero_zero_bonds .and. &
                                 at%z(atom1) == 0 .and. at%z(atom2) == 0) then
                               atom2 = this%next_atom_in_cell(atom2)
                               cycle atom2_loop
                            endif
       
			    ! omit self in the same cell without shift
			    if (.not. my_own_neighbour .and. &
                                 (atom1 == atom2 .and. & 
				 (i4==0 .and. j4==0 .and. k4==0) .and. &
				 (i==i3 .and. j==j3 .and. k==k3))) then
                               atom2 = this%next_atom_in_cell(atom2)
                               cycle atom2_loop
                            endif
			    call test_form_bond(this,at%cutoff, at%use_uniform_cutoff, &
			      at%Z, at%pos, at%lattice, atom1,atom2, &
				 (/i4-map_shift(1,atom1)+map_shift(1,atom2),j4-map_shift(2,atom1)+map_shift(2,atom2),k4-map_shift(3,atom1)+map_shift(3,atom2)/), &
				 .false., error)
                            PASS_ERROR(error)

                            atom2 = this%next_atom_in_cell(atom2)
			 end do atom2_loop ! atom2

		      end do ! i2
		   end do ! j2
		end do ! k2

                atom1 = this%next_atom_in_cell(atom1)
	     end do atom1_loop ! atom1
	  end do ! i
       end do ! j
    end do ! k

    if (my_store_is_min_image) then
       if (allocated(this%is_min_image)) deallocate(this%is_min_image)
       allocate(this%is_min_image(at%n))
       do i=1,at%n
          this%is_min_image(i) = is_min_image(this, i, error=error)
          PASS_ERROR(error)
       end do
    end if

    if (my_store_n_neighb) then
       call add_property(at, 'n_neighb', 0, ptr=n_neighb, overwrite=.true.)
       do i=1,at%n
          n_neighb(i) = n_neighbours(this, i)
       end do
    end if

    call system_timer('calc_connect')

  end subroutine connection_calc_connect


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% The subroutine 'calc_dists' updates the stored distance tables using 
  !% the stored connectivity and shifts. This should be called every time
  !% any atoms are moved (e.g. it is called by 'DynamicalSystem%advance_verlet').
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine connection_calc_dists(this, at, parallel, error)
    type(Connection), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    logical, optional, intent(in) :: parallel
    integer, optional, intent(out) :: error
    integer                    :: i, j, n, index
    integer, dimension(3)      :: shift
    real(dp), dimension(3)     :: j_pos
    logical :: do_parallel
#ifdef _MPI
    integer:: Nelements, mpi_pos, mpi_old_pos
    include "mpif.h"
    real(dp), allocatable :: mpi_send(:), mpi_recv(:)
    integer err
#endif
    INIT_ERROR(error)

    call system_timer('calc_dists')

    ! Flag to specify whether or not to parallelise calculation.
    ! Only actually run in parallel if parallel==.true. AND
    ! _MPI is #defined. Default to serial mode.
    do_parallel = .false.
    if (present(parallel)) do_parallel = parallel

#ifdef _MPI
    if (do_parallel) then
       ! Nelements = sum(this%connect%neighbour1(i)%t%N)
       Nelements = 0
       do i=1,at%N
          Nelements = Nelements + this%neighbour1(i)%t%N
       end do

       allocate(mpi_send(Nelements))
       allocate(mpi_recv(Nelements))
       if (Nelements > 0) then
	 mpi_send = 0.0_dp
	 mpi_recv = 0.0_dp
	end if
       mpi_pos = 1
    end if
#endif

    if (.not.this%initialised) then
         RAISE_ERROR('CalcDists: Connect is not yet initialised', error)
    endif

    if (at%N > 0) then
       if (size(this%neighbour1(1)%t%real,1) == 4) then ! store_rij was set
	    RAISE_ERROR("CalcDists: can't have store_rij set", error)
       endif
    endif

    do i = 1, at%N

#ifdef _MPI
       if (do_parallel) then
          mpi_old_pos = mpi_pos
          mpi_pos = mpi_pos + this%neighbour1(i)%t%N

          ! cycle loop if processor rank does not match
          if(mod(i, mpi_n_procs()) .ne. mpi_id()) cycle
       end if
#endif

       do n = 1, connection_n_neighbours(this, i) 

          j = connection_neighbour_minimal(this, i, n, shift=shift, index=index)

          ! j_pos = at%pos(:,j) + ( at%lattice .mult. shift )
          j_pos(:) = at%pos(:,j) + ( at%lattice(:,1) * shift(1) + at%lattice(:,2) * shift(2) + at%lattice(:,3) * shift(3) )

          if (i <= j) then
             this%neighbour1(i)%t%real(1,index) = norm(j_pos - at%pos(:,i))
          else
             this%neighbour1(j)%t%real(1,index) = norm(j_pos - at%pos(:,i))
          end if

       end do

#ifdef _MPI
       if (do_parallel) then
	  if (mpi_old_pos <= mpi_pos-1) then
	    mpi_send(mpi_old_pos:mpi_pos-1) = &
		 this%neighbour1(i)%t%real(1,1:this%neighbour1(i)%t%N)
	  end if
       end if
#endif

    end do

#ifdef _MPI
    if (do_parallel) then
       ! collect mpi results
       if (Nelements > 0) then
	 call mpi_allreduce(mpi_send, mpi_recv, &
	      size(mpi_send), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
	 call abort_on_mpi_error(err, "Calc_Dists: MPI_ALL_REDUCE()")
       end if

       mpi_pos = 1
       do i=1, at%N
	  if (this%neighbour1(i)%t%N > 0) then
	    this%neighbour1(i)%t%real(1,1:this%neighbour1(i)%t%N) = &
		 mpi_recv(mpi_pos:mpi_pos+this%neighbour1(i)%t%N-1)
	    mpi_pos = mpi_pos + this%neighbour1(i)%t%N
	  endif
       end do

       if (Nelements > 0) then
	 deallocate(mpi_send, mpi_recv)
       end if
    end if
#endif

    call system_timer('calc_dists')

  end subroutine connection_calc_dists


   subroutine get_min_max_images(is_periodic, cellsNa, cellsNb, cellsNc, cell_image_Na, cell_image_Nb, cell_image_Nc, i, j, k, do_i, do_j, do_k, &
	 min_cell_image_Na, max_cell_image_Na, min_cell_image_Nb, max_cell_image_Nb, min_cell_image_Nc, max_cell_image_Nc)
      logical, intent(in) :: is_periodic(3)
      integer, intent(in) :: cellsNa, cellsNb, cellsNc, cell_image_Na, cell_image_Nb, cell_image_Nc, i, j, k
      logical, intent(in) :: do_i, do_j, do_k
      integer, intent(out) :: min_cell_image_Na, max_cell_image_Na, min_cell_image_Nb, max_cell_image_Nb, min_cell_image_Nc, max_cell_image_Nc

      if (do_i) then
	 if (is_periodic(1)) then
	   min_cell_image_Na = -cell_image_Na
	   max_cell_image_Na = cell_image_Na
	 else
	   if (cell_image_Na < i) then 
	      min_cell_image_Na = -cell_image_Na
	   else
	      min_cell_image_Na = -i+1
	   endif
	   if (cell_image_Na < cellsNa-i) then
	      max_cell_image_Na = cell_image_Na
	   else
	      max_cell_image_Na = cellsNa - i
	   endif
	 endif
      endif
      if (do_j) then
	 if (is_periodic(2)) then
	   min_cell_image_Nb = -cell_image_Nb
	   max_cell_image_Nb = cell_image_Nb
	 else
	   if (cell_image_Nb < j) then 
	      min_cell_image_Nb = -cell_image_Nb
	   else
	      min_cell_image_Nb = -j+1
	   endif
	   if (cell_image_Nb < cellsNb-j) then
	      max_cell_image_Nb = cell_image_Nb
	   else
	      max_cell_image_Nb = cellsNb - j
	   endif
	 endif
      endif
      if (do_k) then
	 if (is_periodic(3)) then
	   min_cell_image_Nc = -cell_image_Nc
	   max_cell_image_Nc = cell_image_Nc
	 else
	   if (cell_image_Nc < k) then 
	      min_cell_image_Nc = -cell_image_Nc
	   else
	      min_cell_image_Nc = -k+1
	   endif
	   if (cell_image_Nc < cellsNc-k) then
	      max_cell_image_Nc = cell_image_Nc
	   else
	      max_cell_image_Nc = cellsNc - k
	   endif
	 endif
      endif
      
      if (current_verbosity() >= PRINT_ANAL) then
         call print('get_min_max_images cell_image_na min='//min_cell_image_na//' max='//max_cell_image_na, PRINT_ANAL)
         call print('get_min_max_images cell_image_nb min='//min_cell_image_nb//' max='//max_cell_image_nb, PRINT_ANAL)
         call print('get_min_max_images cell_image_nc min='//min_cell_image_nc//' max='//max_cell_image_nc, PRINT_ANAL)
      end if

   end subroutine get_min_max_images

  !
  !% Spatially partition the atoms into cells. The number of cells in each dimension must already be
  !% set (cellsNa,b,c). Pre-wiping of the cells can be skipped (e.g. if they are already empty).
  !
  subroutine partition_atoms(this, at, dont_wipe, error)

    type(Connection), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    logical, optional, intent(in)    :: dont_wipe
    integer, intent(out), optional :: error

    logical                          :: my_dont_wipe, neighbour1_allocated
    integer                          :: i,j,k,n
    real(dp) :: lat_pos(3)
    integer, pointer :: map_shift(:,:)

    INIT_ERROR(error)
    ! Check inputs
    if (.not.this%cells_initialised) then 
       RAISE_ERROR('Partition_Atoms: Cells have not been initialised', error)
    end if
    my_dont_wipe = .false.
    if (present(dont_wipe)) my_dont_wipe = dont_wipe

    ! Wipe the cells
    if (.not.my_dont_wipe) then
       call wipe_cells(this)
    else
       RAISE_ERROR("Partition_Atoms: dont_wipe=.true. not supported", error)
    endif

!!    ! Make sure all atomic positions are within the cell
!!    call map_into_cell(at)
    if (.not. assign_pointer(at, 'map_shift', map_shift)) then
       call add_property(at, 'map_shift', 0, 3)
       if (.not. assign_pointer(at, 'map_shift', map_shift)) then
	  RAISE_ERROR("partition_atoms impossibly failed to assign map_shift pointer", error)
       end if
    endif

    neighbour1_allocated = allocated(this%neighbour1)

    if (allocated(this%next_atom_in_cell)) then
       if (size(this%next_atom_in_cell) < this%N) then
          deallocate(this%next_atom_in_cell)
       endif
    endif
    if (.not. allocated(this%next_atom_in_cell) .and. this%N > 0) then
       allocate(this%next_atom_in_cell(this%N))
    endif

    do n = 1, at%N
       if (neighbour1_allocated) then
          if (.not. associated(this%neighbour1(n)%t)) cycle ! not in active subregion
       end if

       ! figure out shift to map atom into cell
       lat_pos = at%g .mult. at%pos(:,n)
       map_shift(:,n) = - floor(lat_pos+0.5_dp)
       ! do the mapping
       lat_pos = lat_pos + map_shift(:,n)

       call cell_of_pos(this, lat_pos, i, j, k)
       !Add the atom to this cell
       this%next_atom_in_cell(n) = this%cell_heads(i, j, k)
       this%cell_heads(i, j, k) = n
    end do

  end subroutine partition_atoms

  subroutine cell_of_pos(this, lat_pos, i, j, k)
    type(Connection), intent(in) :: this
    real(dp), intent(in) :: lat_pos(3)
    integer, intent(out) :: i, j, k

    real(dp) :: t(3)

     t = lat_pos
     i = floor(real(this%cellsNa,dp) * (t(1)+0.5_dp)) + 1
     j = floor(real(this%cellsNb,dp) * (t(2)+0.5_dp)) + 1
     k = floor(real(this%cellsNc,dp) * (t(3)+0.5_dp)) + 1

     ! Very small numerical errors in the lattice inverse can lead to bad i,j,k values.
     ! Test for this:
     if (i < 1) then
	i = 1
     else if (i > this%cellsNa) then
	i = this%cellsNa
     end if

     if (j < 1) then
	j = 1
     else if (j > this%cellsNb) then
	j = this%cellsNb
     end if

     if (k < 1) then
	k = 1
     else if (k > this%cellsNc) then
	k = this%cellsNc
     end if

    end subroutine cell_of_pos


    function cell_n(this, i, j, k, error)
      type(Connection), intent(in) :: this
      integer, intent(in) :: i, j, k
      integer :: cell_n
      integer :: a
      integer, intent(out), optional :: error
      
      INIT_ERROR(error)
      if (.not. this%cells_initialised) then
         RAISE_ERROR('cell_n: cells are not initialised', error)
      end if

      cell_n = 0
      a = this%cell_heads(i, j, k)
      do while (a > 0)
         cell_n = cell_n + 1
         a = this%next_atom_in_cell(a)
      enddo

    end function cell_n

!   Don't understand what this is doing, and doesn't seem to be used
!   anywhere - Lars
!    function cell_contents(this, i, j, k, n, error)
!      type(Connection), intent(in) :: this
!      integer, intent(in) :: i, j, k, n
!      integer :: cell_contents
!      integer, intent(out), optional :: error
!
!      INIT_ERROR(error)
!      if (.not. this%cells_initialised) then
!         RAISE_ERROR('cell_n: cells are not initialised', error)
!      end if
!
!      cell_contents = this%cell(i,j,k)%int(1,n)
!
!    end function cell_contents

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !
   ! Geometry procedures
   !
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   !
   !% Given a simulation cell defined by lattice vectors, how many
   !% times can the cell be divided along the lattice vectors into
   !% subcells such that a sphere of radius 'cutoff' with centre in
   !% one subcell does not spill out of the surrounding $3 \times 3$ subcell block?
   !
   subroutine divide_cell(lattice,cutoff,Na,Nb,Nc)

      real(dp), dimension(3,3), intent(in)  :: lattice  !% Box defined by lattice vectors
      real(dp),                 intent(in)  :: cutoff   !% Radius of sphere
      integer,                  intent(out) :: Na,Nb,Nc !% Number of supercells required along $x$, $y$ and $z$
      real(dp)                              :: cellVol
      real(dp), dimension(3)                :: a, b, c

      a = lattice(:,1); b = lattice(:,2); c = lattice(:,3)
      cellVol = abs( scalar_triple_product(a,b,c) )

      Na = max(1,int( cellVol / (cutoff * norm(b .cross. c) ) )) ! round down to nearest int >= 1
      Nb = max(1,int( cellVol / (cutoff * norm(c .cross. a) ) ))
      Nc = max(1,int( cellVol / (cutoff * norm(a .cross. b) ) ))

      call print('divide_cell: '//Na//'x'//Nb//'x'//Nc//' cells.', PRINT_NERD)

   end subroutine divide_cell

   !
   !% Returns the maximum cutoff radius for 'calc_connect', given the lattice if we want to avoid image neghbours
   !
   function max_cutoff(lattice, error)

     real(dp), dimension(3,3), intent(in) :: lattice
     real(dp)                             :: Max_Cutoff
     real(dp), dimension(3)               :: a,b,c
     real(dp)                             :: cellVol, ra, rb, rc
     integer, intent(out), optional :: error

     INIT_ERROR(error)
     a = lattice(:,1); b = lattice(:,2); c = lattice(:,3)
     cellVol = abs( scalar_triple_product(a,b,c) )

     if(cellVol == 0.0_dp) then
        RAISE_ERROR("Max_cutoff(): cell volume is exactly 0.0!", error)
     end if

     ra = cellVol / norm(b .cross. c)
     rb = cellVol / norm(c .cross. a)
     rc = cellVol / norm(a .cross. b)

     Max_Cutoff = 0.5_dp * min(ra,rb,rc)

   end function max_cutoff

  subroutine connection_print(this,file,error)

    type(Connection), intent(in)    :: this
    type(Inoutput),   optional, intent(inout) :: file
    integer, intent(out), optional :: error
    integer                         :: i,j,k,n,a
    integer                         :: ndata(10)

    INIT_ERROR(error)
    if (.not.this%initialised) then
       RAISE_ERROR('Connection_Print: Connection object not initialised', error)
    end if

    call print('Connectivity data:',file=file)
    call print('-------------------------------------------------',file=file)
    call print('|    I    |    J    |    Shift     | Distance   |',file=file)
    call print('-------------------------------------------------',file=file)

    do i = 1, size(this%neighbour1)

       if (.not. associated(this%neighbour1(i)%t)) cycle

       if ((this%neighbour1(i)%t%N + this%neighbour2(i)%t%N) > 0) then
	 write(line,'(a47)')'| Neighbour1 (i <= j)                           |'
	 call print(line, file=file)
       endif

       do j = 1, this%neighbour1(i)%t%N
          if(j == 1) then
             write(line,'(a2,i7,a3,i7,a3,3i4,a3,f10.5,a2)')'| ',    &
                 & i, ' | ', this%neighbour1(i)%t%int(1,j),' | ',this%neighbour1(i)%t%int(2:4,j),&
                 &' | ',this%neighbour1(i)%t%real(1,j),' |'
          else
             write(line,'(a12,i7,a3,3i4,a3,f10.5,a2)')'|         | ',    &
                  this%neighbour1(i)%t%int(1,j),' | ',this%neighbour1(i)%t%int(2:4,j),' | ',this%neighbour1(i)%t%real(1,j),' |'
          end if
          call print(line,file=file)
       end do

       if (this%neighbour2(i)%t%N > 0) then

          write(line,'(a47)')'-----------------------------------------------'
          call print(line, file=file)
          write(line,'(a47)')'| Neighbour2 (i > j)                            |'
          call print(line, file=file)


          do j = 1, this%neighbour2(i)%t%N

             if(j == 1) then
                write(line,'(a2,i7,a3,i7,a3,i7,a)') '| ', i, ' | ', this%neighbour2(i)%t%int(1,j),' | ',&
                       &this%neighbour2(i)%t%int(2,j),'                 |'
             else
                write(line,'(a12,i7,a3,i7,a)') '|         | ', this%neighbour2(i)%t%int(1,j),' | ',&
                       &this%neighbour2(i)%t%int(2,j),'                 |'
             end if
             call print(line,file=file)
          end do

       end if

       if ((this%neighbour1(i)%t%N + this%neighbour2(i)%t%N) > 0) then
	 write(line,'(a47)')'-----------------------------------------------'
	 call print(line,file=file)
       endif

    end do

    write(line,'(a47)')''
    call print(line,file=file)


    !Print the cell lists if they exist
    if (this%cells_initialised .and. current_verbosity() > PRINT_NORMAL) then
       call verbosity_push_decrement()

       write(line,'(a11)')'Cell Lists:'
       call print(line,file=file)
       write(line,'(a70)')'----------------------------------------------------------------------'
       call print(line,file=file)

       ! Write atomic indices 10 per line
       do k = 1, this%cellsNc
          do j = 1, this%cellsNb
             do i = 1, this%cellsNa
                write(line,'(a7,i0,a1,i0,a1,i0,a3)')'Cell ( ',i,' ',j,' ',k,' ):'
                call print(line,file=file)

                a = this%cell_heads(i, j, k)
                n = 0
                do while (a > 0)
                   n = n + 1
                   ndata(n) = a
                   if (n == 10) then
                      write (line, '(10i7)')  ndata
                      call print(line,file=file)
                      n = 0
                   endif
                   a = this%next_atom_in_cell(a)
                enddo
                if (n > 0) then
                   write (line, '(10i7)')  ndata(1:n)
                   call print(line,file=file)
                endif

                write(line,'(a70)')'----------------------------------------------------------------------'
                call print(line,file=file)
             end do
          end do
       end do

       call verbosity_pop()

    end if

  end subroutine connection_print


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! 
  ! Simple query functions
  !
  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  !% Return the total number of neighbour that atom $i$ has.
  function connection_n_neighbours(this, i, error) result(n)
    type(Connection),   intent(in)   :: this
    integer,            intent(in)   :: i
    integer,  optional, intent(out)  :: error     
    
    integer :: n

    INIT_ERROR(error)

    if (.not. this%initialised) then
       RAISE_ERROR('connection_n_neighbours: Connection structure has no connectivity data. Call calc_connect first.', error)
    end if

    if (.not. associated(this%neighbour1(i)%t)) then
      n = 0
      return
    endif

    ! All neighbours
    n = this%neighbour1(i)%t%N + this%neighbour2(i)%t%N

  end function connection_n_neighbours


  !% Return the total number of neighbour, i.e. the number of bonds in the system
  function connection_n_neighbours_total(this, error) result(n)
    type(Connection),   intent(in)   :: this
    integer,  optional, intent(out)  :: error     
    
    integer :: i, n

    INIT_ERROR(error)

    if (.not. this%initialised) then
       RAISE_ERROR('connection_n_neighbours: Connection structure has no connectivity data. Call calc_connect first.', error)
    end if

    ! All neighbours
    n = 0
    do i = 1, this%N
       if (.not. associated(this%neighbour1(i)%t)) cycle
       n = n + this%neighbour1(i)%t%N + this%neighbour2(i)%t%N
    enddo

  end function connection_n_neighbours_total


  !% Return the number of neighbour that atom $i$ has.
  !% If the optional arguments max_dist or max_factor are present 
  !% then only neighbours closer than this cutoff are included.
  function connection_n_neighbours_with_dist(this, at, i, max_dist, max_factor, error) result(n)
    type(Connection),   intent(in)   :: this
    type(Atoms),        intent(in)   :: at
    integer,            intent(in)   :: i
    real(dp), optional, intent(in)   :: max_dist, max_factor
    integer,  optional, intent(out)  :: error     
    
    integer :: n

    integer :: j, m
    real(dp) :: r_ij

    INIT_ERROR(error)

    if (.not. this%initialised) then
       RAISE_ERROR('connection_n_neighbours: Connection structure has no connectivity data. Call calc_connect first.', error)
    end if

    if (.not. associated(this%neighbour1(i)%t)) then
      n = 0
      return
    endif

    if (.not. present(max_dist) .and. .not. present(max_factor)) then
       ! All neighbours
       n = this%neighbour1(i)%t%N + this%neighbour2(i)%t%N
    else if (present(max_dist)) then
       ! Only count neighbours within max_dist distance of i
       n = 0
       do m=1,this%neighbour1(i)%t%N+this%neighbour2(i)%t%N
          j = connection_neighbour(this, at, i, m, distance=r_ij)
          if (r_ij < max_dist) n = n + 1
       end do
    else if (present(max_factor)) then
       ! Only count neighbours within max_factor of i
       n = 0
       do m=1,this%neighbour1(i)%t%N+this%neighbour2(i)%t%N
          j = connection_neighbour(this, at, i, m, distance=r_ij)
          if (r_ij < bond_length(at%Z(i),at%Z(j))*max_factor) n = n + 1
       end do
    else
       RAISE_ERROR('connection_n_neighbours: optional arguments max_dist and max_factor must not both be present', error)
    end if

  end function connection_n_neighbours_with_dist

  function connection_neighbour_index(this, i, n, index, t, is_j, error) result(j)
    type(Connection), intent(in) :: this
    integer :: i, n
    integer,  intent(out) :: index
    !NB gfortran 4.3.4 seg faults if this is intent(out)
    !NB type(Table), pointer, intent(out) :: t
    type(Table), pointer :: t
    !NB
    logical, intent(out) :: is_j
    integer, intent(out), optional :: error     
    integer :: j

    integer :: i_n1n, j_n1n

    INIT_ERROR(error)

    if (this%initialised) then
       i_n1n = n-this%neighbour2(i)%t%N
       if (n <= this%neighbour2(i)%t%N) then
          j = this%neighbour2(i)%t%int(1,n)
          j_n1n = this%neighbour2(i)%t%int(2,n)
          index = j_n1n
	  t => this%neighbour1(j)%t
	  is_j = .true.
       else if (i_n1n <= this%neighbour1(i)%t%N) then
          j = this%neighbour1(i)%t%int(1,i_n1n)
          index = i_n1n
	  t => this%neighbour1(i)%t
	  is_j = .false.
       else
          RAISE_ERROR('connection_neighbour_index: '//n//' out of range for atom '//i//' Should be in range 1 < n <= '//n_neighbours(this, i), error)
       end if
    else
       RAISE_ERROR('connection_neighbour_index: Connect structure not initialized. Call calc_connect first.', error)
    end if

  end function connection_neighbour_index

  function connection_neighbour_minimal(this, i, n, shift, index) result(j)
    type(Connection), intent(in) :: this
    integer :: i, j, n
    integer,  intent(out) :: shift(3)
    integer,  intent(out) :: index

    type(Table), pointer :: t
    logical :: is_j

    j = connection_neighbour_index(this, i, n, index, t, is_j)

    if (is_j) then
      shift = - t%int(2:4,index)
    else
      shift = t%int(2:4,index)
    endif

  end function connection_neighbour_minimal

  !% Return the index of the $n^{\mbox{\small{th}}}$ neighbour of atom $i$. Together with the
  !% previous function, this facilites a loop over the neighbours of atom $i$. Optionally, we
  !% return other geometric information, such as distance, direction cosines and difference vector,
  !% and also an direct index into the neighbour tables. If $i <= j$, this is an index into 'neighbour1(i)',
  !% if $i > j$, it is an index into 'neighbour1(j)'
  !%
  !%>   do n = 1,atoms_n_neighbours(at, i)
  !%>      j = atoms_neighbour(at, i, n, distance, diff, cosines, shift, index)
  !%>
  !%>      ...
  !%>   end do
  !%
  !% if distance $>$ max_dist, return 0, and do not waste time calculating other quantities
  function connection_neighbour(this, at, i, n, distance, diff, cosines, shift, index, max_dist, jn, error) result(j)
    type(Connection),   intent(in)  :: this
    type(Atoms),        intent(in)  :: at
    integer,            intent(in)  :: i, n
    real(dp), optional, intent(out) :: distance
    real(dp), dimension(3), optional, intent(out) :: diff
    real(dp), optional, intent(out) :: cosines(3)
    integer,  optional, intent(out) :: shift(3)
    integer,  optional, intent(out) :: index
    real(dp), optional, intent(in)  :: max_dist
    integer,  optional, intent(out) :: jn
    integer,  optional, intent(out) :: error     

    real(dp)::mydiff(3), norm_mydiff
    integer ::myshift(3)
    integer ::j, i_n1n, j_n1n, i_njn, m

    INIT_ERROR(error)

    if (.not. associated(this%neighbour1(i)%t)) then
      RAISE_ERROR("called atoms_neighbour on atom " // i // " which has no allocated neighbour1 table", error)
    endif

    ! First we give the neighbour2 entries (i > j) then the neighbour1 (i <= j)
    ! This order chosen to give neighbours in approx numerical order but doesn't matter
    if (this%initialised) then
       i_n1n = n-this%neighbour2(i)%t%N
       if (n <= this%neighbour2(i)%t%N) then
          j = this%neighbour2(i)%t%int(1,n)
          j_n1n = this%neighbour2(i)%t%int(2,n)
          if(present(index)) index = j_n1n
       else if (i_n1n <= this%neighbour1(i)%t%N) then
          j = this%neighbour1(i)%t%int(1,i_n1n)
          if(present(index)) index = i_n1n
       else
          RAISE_ERROR('atoms_neighbour: '//n//' out of range for atom '//i//' Should be in range 1 < n <= '//n_neighbours(this, i), error)
       end if
    else
       RAISE_ERROR('atoms_neighbour: Atoms structure has no connectivity data. Call calc_connect first.', error)
    end if

    if(present(jn)) then
       if(i < j) then
          do i_njn = 1, this%neighbour2(j)%t%N
             if( (this%neighbour2(j)%t%int(1,i_njn)==i) .and. &
                 (this%neighbour2(j)%t%int(2,i_njn)==i_n1n) ) jn = i_njn
          enddo
       elseif(i > j) then
          jn = j_n1n + this%neighbour2(j)%t%N
       else
          do i_njn = 1, this%neighbour1(j)%t%N
             if( (this%neighbour1(j)%t%int(1,i_njn) == i) .and. &
                  all(this%neighbour1(j)%t%int(2:4,i_njn) == -this%neighbour1(i)%t%int(2:4,i_n1n))) &
	       jn = i_njn + this%neighbour2(j)%t%N
          enddo
       endif
    endif
          
    ! found neighbour, now check for optional requests
    if(present(distance)) then
       if(i <= j) then 
          distance = this%neighbour1(i)%t%real(1,i_n1n)
       else
          distance = this%neighbour1(j)%t%real(1,j_n1n)
       end if
       if (present(max_dist)) then
	 if (distance > max_dist) then
	   j = 0
	   return
	 endif
       endif
    else
       if (present(max_dist)) then
	 if (i <= j) then
	   if (this%neighbour1(i)%t%real(1,i_n1n) > max_dist) then
	     j = 0
	     return
	   endif
	 else
	   if (this%neighbour1(j)%t%real(1,j_n1n) > max_dist) then
	     j = 0
	     return
	   endif
	 endif
       endif
    end if

    if(present(diff) .or. present(cosines) .or. present(shift)) then
       if (i <= j) then
          myshift = this%neighbour1(i)%t%int(2:4,i_n1n)
       else
          myshift = -this%neighbour1(j)%t%int(2:4,j_n1n)
       end if

       if(present(shift)) shift = myshift

       if(present(diff) .or. present(cosines)) then
          !mydiff = this%pos(:,j) - this%pos(:,i) + (this%lattice .mult. myshift)
	  if (size(this%neighbour1(i)%t%real,1) == 4) then
	     if (i <= j) then
		mydiff = this%neighbour1(i)%t%real(2:4,i_n1n)
	     else
		mydiff = -this%neighbour1(j)%t%real(2:4,j_n1n)
	     endif
	  else
	     mydiff = at%pos(:,j) - at%pos(:,i)
	  endif
          do m=1,3
             ! forall(k=1:3) mydiff(k) = mydiff(k) + this%lattice(k,m) * myshift(m)
             mydiff(1:3) = mydiff(1:3) + at%lattice(1:3,m) * myshift(m)
          end do
          if(present(diff)) diff = mydiff
          if(present(cosines)) then
            norm_mydiff = sqrt(mydiff(1)*mydiff(1) + mydiff(2)*mydiff(2) + mydiff(3)*mydiff(3))
	    if (norm_mydiff > 0.0_dp) then
	      cosines = mydiff / norm_mydiff
	    else
	      cosines = 0.0_dp
	    endif
	  endif
       end if
    end if

  end function connection_neighbour

  function connection_is_min_image(this, i, error) result(is_min_image)
    type(Connection),  intent(in)   :: this
    integer,           intent(in)   :: i
    integer, optional, intent(out)  :: error

    logical :: is_min_image
    integer :: n, m, NN

    INIT_ERROR(error)

    is_min_image = .true.
    ! First we give the neighbour1 (i <= j) then the neighbour2 entries (i > j) 
    if (this%initialised) then

       if (.not. associated(this%neighbour1(i)%t)) then
          RAISE_ERROR('is_min_image: atoms structure has no connectivity data for atom '//i, error)
       end if

       nn = this%neighbour1(i)%t%N
       do n=1,nn
          if (this%neighbour1(i)%t%int(1,n) == i) then
             is_min_image = .false.
             return
          end if
          do m=n+1,nn
             if (this%neighbour1(i)%t%int(1,n) == this%neighbour1(i)%t%int(1,m)) then
                is_min_image = .false.
                return
             end if
          end do
       end do

       nn = this%neighbour2(i)%t%N
       do n=1,nn
          if (this%neighbour2(i)%t%int(1,n) == i) then
             is_min_image = .false.
             return
          end if
          do m=n+1,nn
             if (this%neighbour2(i)%t%int(1,n) == this%neighbour2(i)%t%int(1,m)) then
                is_min_image = .false.
                return
             end if
          end do
       end do

    else
       RAISE_ERROR('is_min_image: Atoms structure has no connectivity data. Call calc_connect first.', error)
    end if

  endfunction connection_is_min_image 


   !
   ! Fit_Box_In_Cell
   !
   !% Given an orthogonal box, oriented along the cartesian axes with lengths '2*rx', '2*ry' and '2*rz'
   !% and centred on the origin, what parameters must we pass to supercell to make a system big 
   !% enough from our original cell defined by lattice for the box to fit inside?
   !%
   !% The required supercell parameters are returned in 'Na', 'Nb', 'Nc' respectively. If, e.g. 'Na = 1'
   !% then we don\'t need to supercell in the a direction. 'Na > 1' means we must supercell by 'Na' times
   !% in the $x$ direction.
   !
   subroutine fit_box_in_cell(rx,ry,rz,lattice,Na,Nb,Nc)

     real(dp),                 intent(in)  :: rx, ry, rz
     real(dp), dimension(3,3), intent(in)  :: lattice
     integer,                  intent(out) :: Na, Nb, Nc
     !local variables
     real(dp), dimension(3,3)              :: g                 !inverse of lattice
     real(dp)                              :: maxa, maxb, maxc  !the current maximum in the a,b,c directions
     real(dp), dimension(3)                :: lattice_coords    !the coordinates of a corner in terms of a,b,c
     integer                               :: i,j,k
     !This subroutine works by calculating the coordinates of the corners of the box
     !in terms of the lattice vectors, taking the maximum in the lattice directions and
     !rounding up to the nearest integer.

     call matrix3x3_inverse(lattice,g)
     maxa = 0.0_dp; maxb=0.0_dp; maxc=0.0_dp

     !For the time being, do all eight corners. This is most likely overkill, but it's still quick.
     do k = -1, 1, 2
        do j = -1, 1, 2
           do i = -1, 1, 2
              lattice_coords = g .mult. (/i*rx,j*ry,k*rz/)
              if (abs(lattice_coords(1)) > maxa) maxa = abs(lattice_coords(1))
              if (abs(lattice_coords(2)) > maxb) maxb = abs(lattice_coords(2))
              if (abs(lattice_coords(3)) > maxc) maxc = abs(lattice_coords(3))
           end do
        end do
     end do

     Na = ceiling(2.0_dp*maxa)
     Nb = ceiling(2.0_dp*maxb)
     Nc = ceiling(2.0_dp*maxc)

   end subroutine fit_box_in_cell

endmodule Connection_module
