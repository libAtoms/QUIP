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
!X  Atoms_types module
!X
!% 'Header' module that contains the data structures for
!%    Atoms
!%    Connection
!%    DomainDecomposition
!% plus methods for manipulation properties.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module Atoms_types_module

  use error_module
  use system_module
  use units_module
  use periodictable_module
  use linearalgebra_module
  use MPI_context_module
  use extendable_str_module
  use dictionary_module
  use table_module

  private

  public :: DEFAULT_NNEIGHTOL
  real(dp), parameter :: DEFAULT_NNEIGHTOL = 1.2_dp    !% Default value for 'atoms%nneightol'

  public :: Table_pointer
  type Table_pointer
    type(Table), pointer :: t => null()
  end type Table_pointer

  public :: Connection
  type Connection

     !% We do not use a minimum image convention, rather, collect all the images of a neigbouring atom
     !% that fall withing the neighbour cutoff. The different images are made distinct in the connection list
     !% by having different 'shift' vectors associated with them. 
     !%
     !% To save storage, the 'neighbour1' table contains all information about the connection
     !% but is only filled in for $i <= j$. 'neighbour2' is just a list of those of $i$'s neighbours 
     !% with $i > j$ together with an index into the 'neighbour1' table of atom $j$.
     !%
     !% In normal use (i.e. outside this module) you don\'t need direct access to the tables 
     !% so you should use the interface functions 'atoms_n_neighbours' and 'atoms_neighbour' which 
     !% hide the distiction between the cases $i <= j$ and $i > j$.
     !%
     !% \begin{tabular}{|c|c|c|c|c|}
     !% \hline
     !% \multicolumn{4}{|c|}{'neighbour1(i)%int'} & 'neighbour1(i)%real' \\
     !% \hline
     !% 1 & 2 & 3 & 4 & 1 \\
     !% \hline
     !% $j$ & 'shift_a' & 'shift_b' & 'shift_c' & $r_{ij}$ \\
     !% \hline
     !% \end{tabular}
     !%
     !% \begin{tabular}{|c|c|}
     !% \hline
     !% \multicolumn{2}{|c|}{'neighbour2(i)%int}' \\
     !% \hline
     !% 1 & 2 \\
     !% \hline
     !% $j$ & $n$ \\
     !% \hline
     !% \end{tabular}
     !%
     !% N.B. If $i$ and $j$ are neighbours with shift 'shift',
     !% then 'norm(atoms%pos(j) - atoms%pos(i) + shift)'
     !% is a minimum. 
     !% Mnemonic: 'shift' is added to $j$ to get closer to $i$.

     ! neighbour1 (i): i <= j for all j in table
     !
     ! int: Nint=4                            real: Nreal=1
     ! -----------------------------------    --------
     ! | 1 |    2    |    3    |    4    |    |   1  |
     ! -----------------------------------    --------
     ! | j | shift_a | shift_b | shift_c |    | r_ij |
     ! -----------------------------------    --------
     !
     !
     ! neighbour2 (i): i > j for all j in table
     !
     ! int: Nint=1   real: Nreal=0
     ! ---------              --
     ! | 1 | 2 |              ||
     ! ---------              --
     ! | j | n |              ||
     ! ---------              --
     !
     !
     ! cell(i,j,k)
     !
     ! int: Nint=1   real: Nreal=0
     ! --------             --
     ! |   1  |             ||
     ! --------             --
     ! | atom |             ||
     ! --------             --

     logical                                    :: initialised = .false.
     logical                                    :: cells_initialised = .false.
     logical                                    :: too_few_cells_warning_issued = .false.

     integer                                    :: cellsNa !% no. of cells in the lattice directions
     integer                                    :: cellsNb, cellsNc

     integer                                    :: N = -1  !% no. of atoms at last calc_connect

     type(table_pointer), allocatable, dimension(:)     :: neighbour1 !% Neighbour information for pairs $i <= j$. 
                                                              !% Contains full details of $j$, $r_{ij}$ and shift.
     type(table_pointer), allocatable, dimension(:)     :: neighbour2 !% Neighbour information for pairs $i > j$.
                                                              !% Simply contains $j$ and a reference to $j$'s
                                                              !% 'neighbour1' table.

     type(table), allocatable, dimension(:,:,:) :: cell    !% For the linear scaling connection calculator

     logical, allocatable, dimension(:) :: is_min_image   !% True if i is a minimum image

  end type Connection


  public :: Atoms
  type Atoms

     ! Self-deallocating object
     logical                               :: own_this = .false.  !% Do I own myself?
     integer                               :: ref_count = 0  !% Reference counter

     logical                               :: fixed_size = .false. !% Can the number of atoms be changed after initialisation?
     logical                               :: domain_decomposed = .false.  !% Is this Atoms object domain decomposed?
     integer                               :: N = 0 !% The number of atoms held (including ghost particles)
     integer                               :: Ndomain = 0 !% The number of atoms held by the local process (excluding ghost particles)
     integer                               :: Nbuffer = 0 !% The number of atoms that can be stored in the buffers of this Atoms object

     logical                               :: use_uniform_cutoff = .false. !% Rather than covalent radii --- 
                                                                           !% default is variable cutoff.
     real(dp)                              :: cutoff = DEFAULT_NNEIGHTOL, cutoff_break = DEFAULT_NNEIGHTOL
     !% if 'use_uniform_cutoff' is true, cutoff is the cutoff distance in \AA{}.
     !% Otherwise, cutoff is a multiplier for 'bond_length(Zi,Zj)'.

     real(dp)                              :: nneightol = DEFAULT_NNEIGHTOL 
                                              !% Count as nearest neighbour if sum of covalent radii
                                              !% times 'this%nneightol' greater than distance between atoms.
                                              !% Used in cluster carving.

     real(dp),              dimension(3,3) :: lattice    !% Lattice vectors, as columns:
     !%\begin{displaymath}
     !%\left(
     !%\begin{array}{ccc}
     !% | & | & | \\
     !% \mathbf{a} & \mathbf{b} & \mathbf{c} \\
     !% | & | & | \\
     !%\end{array}
     !%\right)
     !% = \left(
     !%\begin{array}{ccc}
     !% R_{11} & R_{12} & R_{13} \\
     !% R_{21} & R_{22} & R_{23} \\
     !% R_{31} & R_{32} & R_{33} \\
     !%\end{array}
     !%\right)
     !%\end{displaymath}
     !% i.e. $\mathbf{a} = $ 'lattice(:,1)', $\mathbf{b} = $ 'lattice(:,2)' and
     !% $\mathbf{c} = $ 'lattice(:,3)'.

     !  ( | | | | | | )          ( (1,1) (1,2) (1,3) ) 
     ! (  | | | | | |  )        (                     )
     ! (  |a| |b| |c|  )    =   (  (2,1) (2,2) (2,3)  )
     ! (  | | | | | |  )        (                     )
     !  ( | | | | | | )          ( (3,1) (3,2) (3,3) ) 
     logical :: is_orthorhombic, is_periodic(3)


     real(dp),              dimension(3,3) :: g          !% Inverse lattice (stored for speed)

     type(Dictionary) :: properties !% Per-atom data, stored as Dictionary arrays of shape (this%N) or (n_cols,this%N)
     type(Dictionary) :: params     !% List of key/value pairs read from comment line of XYZ file


     integer,  pointer, dimension(:)   :: Z      => null()  !% Atomic numbers, dimension is actually $(N)$
     character(1), pointer, dimension(:,:) :: species => null() !% Names of elements

     integer, pointer, dimension(:)    :: move_mask => null()  !% Atoms with 'move_mask' set to false are fixed
     integer, pointer, dimension(:)    :: damp_mask => null()  !% Damping is only applied to those atoms with
                                                              !% 'damp_mask' set to 1. 
                                                              !% By default this is set to 1 for all atoms.

     integer, pointer, dimension(:)    :: thermostat_region => null() !% Which thermostat is applied to each atoms. 
                                                              !% By default this is set to 1 for all atoms.

     integer,  pointer, dimension(:,:) :: travel => null()  !% Travel across periodic conditions. Actually $(3,N)$ array.
                                                            !% See 'map_into_cell' below.
     real(dp), pointer, dimension(:,:) :: pos    => null()  !% $(3,N)$ array of atomic positions, in \AA. 
                                                            !% Position of atom $i$ is 'pos(:,i)'
     real(dp), pointer, dimension(:)   :: mass   => null()  !% Atomic masses, dimension is actually $(N)$

     real(dp), pointer, dimension(:,:) :: velo   => null()  !% $(3,N)$ array  of atomic velocities, in \AA/fs.
     real(dp), pointer, dimension(:,:) :: acc    => null()  !% $(3,N)$ array  of accelerations in \AA/fs$^2$
     real(dp), pointer, dimension(:,:) :: avgpos => null()  !% $(3,N)$ array  of time-averaged atomic positions.
     real(dp), pointer, dimension(:,:) :: oldpos => null()  !% $(3,N)$ array  of positions of atoms at previous time step.
     real(dp), pointer, dimension(:) :: avg_ke => null()    !% Time-averaged atomic kinetic energy
     type(Connection)                      :: connect       !% Connectivity object (see above)
     type(Connection)                      :: hysteretic_connect       !% Hysteretic connectivity object (see above)

  end type Atoms

  public :: atoms_repoint

  !% Add a per-atom property to this atoms object, as extra entry with columns of 
  !% integers, reals, logical, or strings in the 'properties' dictionary. For example, 
  !% this interface is used by the DynamicalSystems module to create the 'velo', 'acc', 
  !% etc. properties.
  !% Optionally, a pointer to the new property is returned.
  public :: add_property
  interface add_property
     module procedure atoms_add_property_int, atoms_add_property_int_a
     module procedure atoms_add_property_real, atoms_add_property_real_a
     module procedure atoms_add_property_str, atoms_add_property_str_2da
     module procedure atoms_add_property_str_a
     module procedure atoms_add_property_logical, atoms_add_property_logical_a 
     module procedure atoms_add_property_int_2Da
     module procedure atoms_add_property_real_2Da
  end interface

  !% Add a per-atom property to this atoms object, but point to existing space
  !% rather than allocating new space for it (as add_property does).
  public :: add_property_from_pointer
  interface add_property_from_pointer
     module procedure atoms_add_property_p_int, atoms_add_property_p_int_a
     module procedure atoms_add_property_p_real, atoms_add_property_p_real_a
     module procedure atoms_add_property_p_str
     module procedure atoms_add_property_p_logical
  end interface


  !% This is a convenience interface to assign pointers to custom properties of this
  !% Atoms object. The pointer is simply directed to the relevant array in the
  !% 'this%properties' Dictionary. Adding or removing atoms will invalidate these pointers.
  !% Returns true if successful, false if property doesn't exist or
  !% if property type and type of pointer don't match.
  !% OUTDATED - replace by subroutine assign_property_pointer with error handling
  public :: assign_pointer
  interface assign_pointer
     module procedure atoms_assign_pointer_int1D, atoms_assign_pointer_int2D
     module procedure atoms_assign_pointer_real1D, atoms_assign_pointer_real2D
     module procedure atoms_assign_pointer_str
     module procedure atoms_assign_pointer_logical
  end interface


  !% This is a convenience interface to assign pointers to custom properties of this
  !% Atoms object. The pointer is simply directed to the relevant array in the
  !% 'this%properties' Dictionary. Adding or removing atoms will invalidate these pointers.
  !% Raises error for failure
  public :: assign_property_pointer
  interface assign_property_pointer
     module procedure atoms_assign_prop_ptr_int1D, atoms_assign_prop_ptr_int2D
     module procedure atoms_assign_prop_ptr_real1D, atoms_assign_prop_ptr_real2D
     module procedure atoms_assign_prop_ptr_str
     module procedure atoms_assign_prop_ptr_logical
  end interface


  ! Note: These are Atoms unrelated, and could be easily move somewhere else
  ! (linearalgebra?)
  public :: map_into_cell
  interface map_into_cell
    module procedure vec_map_into_cell
  end interface

  public :: cell_volume
  interface cell_volume
    module procedure lattice_cell_volume
  end interface

  public :: bond_length

contains
  
  !% OMIT
  ! Initialise pointers for convenient access to special columns of this%properties
  subroutine atoms_repoint(this)
    type(Atoms), target, intent(inout) :: this
    integer :: i

    if(this%N == 0) return

    nullify(this%Z, this%travel, this%pos, this%move_mask, this%damp_mask, &
         this%thermostat_region, this%pos, this%velo, this%acc, this%avgpos, &
         this%oldpos)

    ! Loop over all properties looking for those with special names
    ! which we have pointers for
    do i = 1,this%properties%N

       ! If this%N is zero then point at zero length arrays
       if (this%N == 0) then

          select case(trim(lower_case(string(this%properties%keys(i)))))

          ! Integer properties
          case('z')
             this%Z               => this%properties%entries(i)%i_a(:)
          case('travel')
             this%travel          => this%properties%entries(i)%i_a2(:,:)
          case('move_mask')
             this%move_mask       => this%properties%entries(i)%i_a(:)
          case('damp_mask')
             this%damp_mask       => this%properties%entries(i)%i_a(:)
          case('thermostat_region')
             this%thermostat_region => this%properties%entries(i)%i_a(:)

          ! Real properties
          case('mass')
             this%mass            => this%properties%entries(i)%r_a(:)
          case('pos')
             this%pos             => this%properties%entries(i)%r_a2(:,:)
          case('velo')
             this%velo            => this%properties%entries(i)%r_a2(:,:)
          case('acc')
             this%acc             => this%properties%entries(i)%r_a2(:,:)
          case('avgpos')
             this%avgpos          => this%properties%entries(i)%r_a2(:,:)
          case('oldpos')
             this%oldpos          => this%properties%entries(i)%r_a2(:,:)
          case('avg_ke')
             this%avg_ke          => this%properties%entries(i)%r_a(:)

          ! String properties
          case('species')
             this%species         => this%properties%entries(i)%s_a(:,:)

          end select


       else

          select case(trim(lower_case(string(this%properties%keys(i)))))

          ! Integer properties
          case('z')
             this%Z               => this%properties%entries(i)%i_a(1:this%n)
          case('travel')
             this%travel          => this%properties%entries(i)%i_a2(:,1:this%n)
          case('move_mask')
             this%move_mask       => this%properties%entries(i)%i_a(1:this%n)
          case('damp_mask')
             this%damp_mask       => this%properties%entries(i)%i_a(1:this%n)
          case('thermostat_region')
             this%thermostat_region => this%properties%entries(i)%i_a(1:this%n)

          ! Real properties
          case('mass')
             this%mass            => this%properties%entries(i)%r_a(1:this%n)
          case('pos')
             this%pos             => this%properties%entries(i)%r_a2(:,1:this%N)
          case('velo')
             this%velo            => this%properties%entries(i)%r_a2(:,1:this%N)
          case('acc')
             this%acc             => this%properties%entries(i)%r_a2(:,1:this%N)
          case('avgpos')
             this%avgpos          => this%properties%entries(i)%r_a2(:,1:this%N)
          case('oldpos')
             this%oldpos          => this%properties%entries(i)%r_a2(:,1:this%N)
          case('avg_ke')
             this%avg_ke          => this%properties%entries(i)%r_a(1:this%n)

          ! String properties
          case('species')
             this%species         => this%properties%entries(i)%s_a(:,1:this%N)

          end select

       end if

    end do

  end subroutine atoms_repoint

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Add new properties
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_add_property_p_int(this, name, ptr, error)
    type(Atoms), intent(inout), target :: this
    character(len=*), intent(in) :: name
    integer, intent(in), target :: ptr(:)
    integer, intent(out), optional :: error

    integer :: use_n_cols, i

    INIT_ERROR(error)
    if (size(ptr,1) /= this%Nbuffer) then
      RAISE_ERROR("atoms_add_property_p_int_a: incompatible pointer this%Nbuffer="//this%Nbuffer//" pointer shape "//shape(ptr), error)
    endif
    use_n_cols = 1
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (use_n_cols == 1 .and. this%properties%entries(i)%type /= T_INTEGER_A) then
          RAISE_ERROR("atoms_add_property_p_int: incompatible property "//trim(name)//" already present", error)
       end if
    end if
    call set_value_pointer(this%properties, name, ptr)
    call atoms_repoint(this)
  end subroutine atoms_add_property_p_int

  subroutine atoms_add_property_p_int_a(this, name, ptr, error)
    type(Atoms), intent(inout), target :: this
    character(len=*), intent(in) :: name
    integer, intent(in), target :: ptr(:,:)
    integer, intent(out), optional :: error

    integer :: use_n_cols, i

    INIT_ERROR(error)
    if (size(ptr,2) /= this%Nbuffer) then
      RAISE_ERROR("atoms_add_property_p_int_a: incompatible pointer this%Nbuffer="//this%Nbuffer//" pointer shape "//shape(ptr), error)
    endif
    use_n_cols = size(ptr,1)
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (use_n_cols == 1 .and. this%properties%entries(i)%type /= T_INTEGER_A) then
          RAISE_ERROR("atoms_add_property_p_int_a: incompatible property "//trim(name)//" already present", error)
       end if
       if (use_n_cols /= 1 .and. this%properties%entries(i)%type /= T_INTEGER_A2) then
          RAISE_ERROR("atoms_add_property_p_int_a: incompatible property "//trim(name)//" already present", error)
       end if
    end if
    call set_value_pointer(this%properties, name, ptr)
    call atoms_repoint(this)
  end subroutine atoms_add_property_p_int_a

  subroutine atoms_add_property_p_real(this, name, ptr, error)
    type(Atoms), intent(inout), target :: this
    character(len=*), intent(in) :: name
    real(dp), intent(in), target :: ptr(:)
    integer, intent(out), optional :: error

    integer :: use_n_cols, i

    INIT_ERROR(error)
    if (size(ptr,1) /= this%Nbuffer) then
      RAISE_ERROR("atoms_add_property_p_real_a: incompatible pointer this%Nbuffer="//this%Nbuffer//" pointer shape "//shape(ptr), error)
    endif
    use_n_cols = 1
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (use_n_cols == 1 .and. this%properties%entries(i)%type /= T_REAL_A) then
          RAISE_ERROR("atoms_add_property_p_real: incompatible property "//trim(name)//" already present", error)
       end if
    end if
    call set_value_pointer(this%properties, name, ptr)
    call atoms_repoint(this)
  end subroutine atoms_add_property_p_real

  subroutine atoms_add_property_p_real_a(this, name, ptr, error)
    type(Atoms), intent(inout), target :: this
    character(len=*), intent(in) :: name
    real(dp), intent(in), target :: ptr(:,:)
    integer, intent(out), optional :: error

    integer :: use_n_cols, i

    INIT_ERROR(error)
    if (size(ptr,2) /= this%Nbuffer) then
      RAISE_ERROR("atoms_add_property_p_real_a: incompatible pointer this%Nbuffer="//this%Nbuffer//" pointer shape "//shape(ptr), error)
    endif
    use_n_cols = size(ptr,1)
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (use_n_cols == 1 .and. this%properties%entries(i)%type /= T_REAL_A) then
          RAISE_ERROR("atoms_add_property_p_real_a: incompatible property "//trim(name)//" already present", error)
       end if
       if (use_n_cols /= 1 .and. this%properties%entries(i)%type /= T_REAL_A2) then
          RAISE_ERROR("atoms_add_property_p_real_a: incompatible property "//trim(name)//" already present", error)
       end if
    end if
    call set_value_pointer(this%properties, name, ptr)
    call atoms_repoint(this)
  end subroutine atoms_add_property_p_real_a


  subroutine atoms_add_property_p_logical(this, name, ptr, error)
    type(Atoms), intent(inout), target :: this
    character(len=*), intent(in) :: name
    logical, intent(in), target :: ptr(:)
    integer, intent(out), optional :: error

    integer :: use_n_cols, i

    INIT_ERROR(error)
    if (size(ptr,1) /= this%Nbuffer) then
      RAISE_ERROR("atoms_add_property_p_logical: incompatible pointer this%Nbuffer="//this%Nbuffer//" pointer shape "//shape(ptr), error)
    endif
    use_n_cols = 1
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (use_n_cols == 1 .and. this%properties%entries(i)%type /= T_LOGICAL_A) then
          RAISE_ERROR("atoms_add_property_p_logical: incompatible property "//trim(name)//" already present", error)
       end if
    end if
    call set_value_pointer(this%properties, name, ptr)
    call atoms_repoint(this)
  end subroutine atoms_add_property_p_logical

  subroutine atoms_add_property_p_str(this, name, ptr, error)
    type(Atoms), intent(inout), target :: this
    character(len=*), intent(in) :: name
    character(1), intent(in), target :: ptr(:,:)
    integer, intent(out), optional :: error

    integer :: use_n_cols, i

    INIT_ERROR(error)
    if (size(ptr,2) /= this%N) then
      RAISE_ERROR("atoms_add_property_p_str: incompatible pointer this%Nbuffer="//this%Nbuffer//" pointer shape "//shape(ptr), error)
    endif
    use_n_cols = 1
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (use_n_cols == 1 .and. this%properties%entries(i)%type /= T_CHAR_A) then
          RAISE_ERROR("atoms_add_property_p_str: incompatible property "//trim(name)//" already present", error)
       end if
    end if
    call set_value_pointer(this%properties, name, ptr)
    call atoms_repoint(this)
  end subroutine atoms_add_property_p_str

  subroutine atoms_add_property_int(this, name, value, n_cols, ptr, ptr2, overwrite, error)
    type(Atoms), intent(inout), target :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: value
    integer, intent(in), optional :: n_cols
    integer, intent(out), optional, dimension(:), pointer :: ptr
    integer, intent(out), optional, dimension(:,:), pointer :: ptr2
    logical, optional, intent(in) :: overwrite
    integer, intent(out), optional :: error

    integer :: use_n_cols, i

    INIT_ERROR(error)
    use_n_cols = optional_default(1, n_cols)

    ! Check for incompatible property
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (use_n_cols == 1 .and. this%properties%entries(i)%type /= T_INTEGER_A) then
          RAISE_ERROR("atoms_add_property_int: incompatible property "//trim(name)//" already present", error)
       end if

       if (use_n_cols /= 1 .and. this%properties%entries(i)%type /= T_INTEGER_A2) then
          RAISE_ERROR("atoms_add_property_int: incompatible property "//trim(name)//" already present", error)
       end if       
    end if
    
    if (use_n_cols == 1) then
       call add_array(this%properties, name, value, this%Nbuffer, ptr, overwrite)
    else
       call add_array(this%properties, name, value, (/n_cols, this%Nbuffer/), ptr2, overwrite)
    end if

    call atoms_repoint(this)
  end subroutine atoms_add_property_int

  subroutine atoms_add_property_int_a(this, name, value, n_cols, ptr, ptr2, overwrite, error)
    type(Atoms), intent(inout),target :: this
    character(len=*), intent(in) :: name
    integer, intent(in), dimension(:) :: value
    integer, intent(in), optional :: n_cols
    integer, optional, intent(out), dimension(:), pointer :: ptr
    integer, optional, intent(out), dimension(:,:), pointer :: ptr2
    logical, optional, intent(in) :: overwrite
    integer, optional, intent(out) :: error

    integer :: use_n_cols, i
    integer, allocatable, dimension(:,:) :: tmp_value

    INIT_ERROR(error)
    use_n_cols = optional_default(1, n_cols)

    if (size(value) /= this%Nbuffer) then
       RAISE_ERROR('atoms_add_property_int_a: size(value) ('//size(value)//') /= this%Nbuffer ('//this%Nbuffer//')', error)
    end if

    ! Check for incompatible property
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (use_n_cols == 1 .and. this%properties%entries(i)%type /= T_INTEGER_A) then
          RAISE_ERROR("atoms_add_property_int_a: incompatible property "//trim(name)//" already present", error)
       end if

       if (use_n_cols /= 1 .and. this%properties%entries(i)%type /= T_INTEGER_A2) then
          RAISE_ERROR("atoms_add_property_int_a: incompatible property "//trim(name)//" already present", error)
       end if       
    end if

    if (use_n_cols == 1) then
       call add_array(this%properties, name, value, this%Nbuffer, ptr, overwrite)
    else
       allocate(tmp_value(n_cols, size(value)))
       do i=1,n_cols
          tmp_value(i,:) = value
       end do
       call add_array(this%properties, name, tmp_value, (/n_cols, this%Nbuffer/), ptr2, overwrite)
       deallocate(tmp_value)
    end if
   
    call atoms_repoint(this)
  end subroutine atoms_add_property_int_a


  subroutine atoms_add_property_real(this, name, value, n_cols, ptr, ptr2, overwrite, error)
    type(Atoms), intent(inout),target :: this
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    integer, intent(in), optional :: n_cols
    real(dp), optional, dimension(:), pointer, intent(out) :: ptr
    real(dp), optional, dimension(:,:), pointer, intent(out) :: ptr2
    logical, optional, intent(in) :: overwrite
    integer, intent(out), optional :: error

    integer use_n_cols, i

    INIT_ERROR(error)
    use_n_cols = optional_default(1, n_cols)

    ! Check for incompatible property
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (use_n_cols == 1 .and. this%properties%entries(i)%type /= T_REAL_A) then
          RAISE_ERROR("atoms_add_property_real: incompatible property "//trim(name)//" already present", error)
       end if

       if (use_n_cols /= 1 .and. this%properties%entries(i)%type /= T_REAL_A2) then
          RAISE_ERROR("atoms_add_property_real: incompatible property "//trim(name)//" already present", error)
       end if       
    end if
    
    if (use_n_cols == 1) then
       call add_array(this%properties, name, value, this%Nbuffer, ptr, overwrite)
    else
       call add_array(this%properties, name, value, (/n_cols, this%Nbuffer/), ptr2, overwrite)
    end if

    call atoms_repoint(this)
  end subroutine atoms_add_property_real


  subroutine atoms_add_property_real_a(this, name, value, n_cols, ptr, ptr2, overwrite, error)
    type(Atoms), intent(inout),target :: this
    character(len=*), intent(in) :: name
    real(dp), intent(in), dimension(:) :: value
    integer, intent(in), optional :: n_cols
    real(dp), optional, dimension(:), pointer, intent(out) :: ptr
    real(dp), optional, dimension(:,:), pointer, intent(out) :: ptr2
    logical, optional, intent(in) :: overwrite
    integer, intent(out), optional :: error

    integer :: use_n_cols, i
    real(dp), allocatable, dimension(:,:) :: tmp_value    

    INIT_ERROR(error)
    use_n_cols = optional_default(1, n_cols)

    if (size(value) /= this%Nbuffer) then
       RAISE_ERROR('atoms_add_property_real_a: size(value) ('//size(value)//') /= this%Nbuffer ('//this%Nbuffer//')', error)
    end if

    ! Check for incompatible property
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (use_n_cols == 1 .and. this%properties%entries(i)%type /= T_REAL_A) then
          RAISE_ERROR("atoms_add_property_real_a: incompatible property "//trim(name)//" already present", error)
       end if

       if (use_n_cols /= 1 .and. this%properties%entries(i)%type /= T_REAL_A2) then
          RAISE_ERROR("atoms_add_property_real_a: incompatible property "//trim(name)//" already present", error)
       end if       
    end if
    
    if (use_n_cols == 1) then
       call add_array(this%properties, name, value, this%Nbuffer, ptr, overwrite)
    else
       allocate(tmp_value(n_cols, size(value)))
       do i=1,n_cols
          tmp_value(i,:) = value
       end do
       call add_array(this%properties, name, tmp_value, (/n_cols, this%Nbuffer/), ptr2, overwrite)
       deallocate(tmp_value)
    end if

    call atoms_repoint(this)
  end subroutine atoms_add_property_real_a

  subroutine atoms_add_property_int_2Da(this, name, value, ptr, overwrite, error)
    type(Atoms), intent(inout),target :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: value(:,:)
    integer, optional, dimension(:,:), pointer, intent(out) :: ptr
    logical, optional, intent(in) :: overwrite
    integer, intent(out), optional :: error

    integer i

    INIT_ERROR(error)

    if (size(value,2) /= this%Nbuffer) then
       RAISE_ERROR('atoms_add_property_int_2Da: size(value,2) ('//size(value,2)//') /= this%Nbuffer ('//this%Nbuffer//')', error)
    end if

    ! Check for incompatible property
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (size(value,1) == 1 .and. this%properties%entries(i)%type /= T_INTEGER_A) then
          RAISE_ERROR("atoms_add_property_int_2Da: incompatible property "//trim(name)//" already present", error)
       end if

       if (size(value,1) /= 1 .and. this%properties%entries(i)%type /= T_INTEGER_A2) then
          RAISE_ERROR("atoms_add_property_int_2Da: incompatible property "//trim(name)//" already present", error)
       end if       
    end if
    
    if (size(value,2) /= this%Nbuffer) then
       RAISE_ERROR('atoms_add_property_int_2Da: size(value,2)='//size(value,2)//' != this%Nbuffer ='//this%Nbuffer, error)
    end if

    call add_array(this%properties, name, value, (/size(value,1), this%Nbuffer/), ptr, overwrite)

    call atoms_repoint(this)
  end subroutine atoms_add_property_int_2Da


  subroutine atoms_add_property_real_2Da(this, name, value, ptr, overwrite, error)
    type(Atoms), intent(inout),target :: this
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value(:,:)
    real(dp), optional, dimension(:,:), pointer, intent(out) :: ptr
    logical, optional, intent(in) :: overwrite
    integer, intent(out), optional :: error

    integer i

    INIT_ERROR(error)

    if (size(value,2) /= this%Nbuffer) then
       RAISE_ERROR('atoms_add_property_real_2Da: size(value,2) ('//size(value,2)//') /= this%Nbuffer ('//this%Nbuffer//')', error)
    end if

    ! Check for incompatible property
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (size(value,1) == 1 .and. this%properties%entries(i)%type /= T_REAL_A) then
          RAISE_ERROR("atoms_add_property_real_2Da: incompatible property "//trim(name)//" already present", error)
       end if

       if (size(value,1) /= 1 .and. this%properties%entries(i)%type /= T_REAL_A2) then
          RAISE_ERROR("atoms_add_property_real_2Da: incompatible property "//trim(name)//" already present", error)
       end if       
    end if
    
    if (size(value,2) /= this%Nbuffer) then
       RAISE_ERROR('atoms_add_property_real_2Da: size(value,2)='//size(value,2)//' != this%Nbuffer ='//this%Nbuffer, error)
    end if

    call add_array(this%properties, name, value, (/size(value,1), this%Nbuffer/), ptr, overwrite)

    call atoms_repoint(this)
  end subroutine atoms_add_property_real_2Da


  subroutine atoms_add_property_str(this, name, value, ptr, overwrite, error)
    type(Atoms), intent(inout), target :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value
    character(1), intent(out), optional, dimension(:,:), pointer :: ptr
    logical, optional, intent(in) :: overwrite
    integer, intent(out), optional :: error

    integer i

    INIT_ERROR(error)

    ! Check for incompatible property
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (this%properties%entries(i)%type /= T_CHAR_A) then
          RAISE_ERROR("atoms_add_property_str: incompatible property "//trim(name)//" already present", error)
       end if       
    end if
    
    ! temporary hack - string length fixed to TABLE_STRING_LENGTH
    if (len(value) /= TABLE_STRING_LENGTH) then
       RAISE_ERROR("atoms_add_property_str: string properties much have string length TABLE_STRING_LENGTH but got "//len(value), error)
    end if
    call add_array(this%properties, name, value, (/TABLE_STRING_LENGTH, this%Nbuffer/), ptr, overwrite)

    call atoms_repoint(this)
  end subroutine atoms_add_property_str


  subroutine atoms_add_property_str_2da(this, name, value, ptr, overwrite, error)
    type(Atoms), intent(inout),target :: this
    character(len=*), intent(in) :: name
    character(1), dimension(:,:) :: value
    character(1), intent(out), optional, dimension(:,:), pointer :: ptr
    logical, optional, intent(in) :: overwrite
    integer, intent(out), optional :: error

    integer :: i

    INIT_ERROR(error)

    if (size(value,2) /= this%Nbuffer) then
       RAISE_ERROR('atoms_add_property_str_2da: size(value,2) ('//size(value,2)//') /= this%Nbuffer ('//this%Nbuffer//')', error)
    end if

    ! Check for incompatible property
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (this%properties%entries(i)%type /= T_CHAR_A) then
          RAISE_ERROR("atoms_add_property_str: incompatible property "//trim(name)//" already present", error)
       end if       
    end if
    
    ! temporary hack - string length fixed to TABLE_STRING_LENGTH
    if (size(value,1) /= TABLE_STRING_LENGTH) then
       RAISE_ERROR("atoms_add_property_str: string properties much have string length TABLE_STRING_LENGTH", error)
    end if

    call add_array(this%properties, name, value, (/TABLE_STRING_LENGTH, this%Nbuffer/), ptr, overwrite)

    call atoms_repoint(this)
  end subroutine atoms_add_property_str_2da

  subroutine atoms_add_property_str_a(this, name, value, ptr, overwrite, error)
    type(Atoms), intent(inout),target :: this
    character(len=*), intent(in) :: name
    character(len=*), dimension(:) :: value
    character(1), intent(out), optional, dimension(:,:), pointer :: ptr
    logical, optional, intent(in) :: overwrite
    integer, intent(out), optional :: error

    character(1), allocatable, dimension(:,:) :: tmp_value
    integer :: i, j

    INIT_ERROR(error)

    if (size(value) /= this%Nbuffer) then
       RAISE_ERROR('atoms_add_property_str_a: size(value) ('//size(value)//') /= this%Nbuffer ('//this%Nbuffer//')', error)
    end if

    ! Check for incompatible property
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (this%properties%entries(i)%type /= T_CHAR_A) then
          RAISE_ERROR("atoms_add_property_str: incompatible property "//trim(name)//" already present", error)
       end if       
    end if
    
    ! temporary hack - string length fixed to TABLE_STRING_LENGTH
    if (len(value(1)) /= TABLE_STRING_LENGTH) then
       RAISE_ERROR("atoms_add_property_str: string properties much have string length TABLE_STRING_LENGTH", error)
    end if

    allocate(tmp_value(len(value),this%n))
    do i=1,this%Nbuffer
       do j=1,len(value)
          tmp_value(j,i) = value(i)(j:j)
       end do
    end do
    call add_array(this%properties, name, tmp_value, (/TABLE_STRING_LENGTH, this%Nbuffer/), ptr, overwrite)
    deallocate(tmp_value)

    call atoms_repoint(this)
  end subroutine atoms_add_property_str_a

  subroutine atoms_add_property_logical(this, name, value, ptr, overwrite, error)
    type(Atoms), intent(inout), target :: this
    character(len=*), intent(in) :: name
    logical, intent(in) :: value
    logical, intent(out), optional, dimension(:), pointer :: ptr
    logical, optional, intent(in) :: overwrite
    integer, intent(out), optional :: error

    integer i

    INIT_ERROR(error)

    ! Check for incompatible property
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (this%properties%entries(i)%type /= T_LOGICAL_A) then
          RAISE_ERROR("atoms_add_property_logical: incompatible property "//trim(name)//" already present", error)
       end if
    end if
    
    call add_array(this%properties, name, value, this%Nbuffer, ptr, overwrite)
    
    call atoms_repoint(this)
  end subroutine atoms_add_property_logical


  subroutine atoms_add_property_logical_a(this, name, value, ptr, overwrite, error)
    type(Atoms), intent(inout),target :: this
    character(len=*), intent(in) :: name
    logical, dimension(:) :: value
    logical, intent(out), optional, dimension(:), pointer :: ptr
    logical, optional, intent(in) :: overwrite
    integer, intent(out), optional :: error

    integer i

    INIT_ERROR(error)

    if (size(value) /= this%Nbuffer) then
       RAISE_ERROR('atoms_add_property_logical_a: size(value) ('//size(value)//') /= this%Nbuffer ('//this%Nbuffer//')', error)
    end if

    ! Check for incompatible property
    i = lookup_entry_i(this%properties, name)
    if (i /= -1) then
       if (this%properties%entries(i)%type /= T_LOGICAL_A) then
          RAISE_ERROR("atoms_add_property_logical_a: incompatible property "//trim(name)//" already present", error)
       end if
    end if
    
    call add_array(this%properties, name, value, this%Nbuffer, ptr, overwrite)

    call atoms_repoint(this)
  end subroutine atoms_add_property_logical_a



  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Assigning pointers to properties
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function atoms_assign_pointer_int1D(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    integer, pointer :: ptr(:)
    logical :: atoms_assign_pointer_int1D

    atoms_assign_pointer_int1D = assign_pointer(this%properties, name, ptr)
  end function atoms_assign_pointer_int1D

  function atoms_assign_pointer_int2D(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    integer, pointer :: ptr(:,:)
    logical :: atoms_assign_pointer_int2D

    atoms_assign_pointer_int2D = assign_pointer(this%properties, name, ptr)
  end function atoms_assign_pointer_int2D

  function atoms_assign_pointer_real1D(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    real(dp), pointer :: ptr(:)
    logical :: atoms_assign_pointer_real1D

    atoms_assign_pointer_real1D = assign_pointer(this%properties, name, ptr)
  end function atoms_assign_pointer_real1D

  function atoms_assign_pointer_real2D(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    real(dp), pointer :: ptr(:,:)
    logical :: atoms_assign_pointer_real2D

    atoms_assign_pointer_real2D = assign_pointer(this%properties, name, ptr)
  end function atoms_assign_pointer_real2D

  function atoms_assign_pointer_str(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    character(1), pointer :: ptr(:,:)
    logical :: atoms_assign_pointer_str

    atoms_assign_pointer_str = assign_pointer(this%properties, name, ptr)
  end function atoms_assign_pointer_str

  function atoms_assign_pointer_logical(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    logical, pointer :: ptr(:)
    logical :: atoms_assign_pointer_logical

    atoms_assign_pointer_logical = assign_pointer(this%properties, name, ptr)
  end function atoms_assign_pointer_logical

  subroutine atoms_assign_prop_ptr_int1D(this, name, ptr, error)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    integer, pointer :: ptr(:)
    integer, intent(out), optional :: error
    
    logical :: res

    INIT_ERROR(error)

    res = assign_pointer(this%properties, name, ptr)
    if (.not. res) then
      RAISE_ERROR("atoms_assign_prop_ptr_int1d failed to assign pointer to "//trim(name)//" in this%properties", error)
    endif
  end subroutine atoms_assign_prop_ptr_int1D

  subroutine atoms_assign_prop_ptr_int2D(this, name, ptr, error)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    integer, pointer :: ptr(:,:)
    integer, intent(out), optional :: error
    
    logical :: res

    INIT_ERROR(error)

    res = assign_pointer(this%properties, name, ptr)
    if (.not. res) then
      RAISE_ERROR("atoms_assign_prop_ptr_int2d failed to assign pointer to "//trim(name)//" in this%properties", error)
    endif
  end subroutine atoms_assign_prop_ptr_int2D

  subroutine atoms_assign_prop_ptr_real1D(this, name, ptr, error)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    real(dp), pointer :: ptr(:)
    integer, intent(out), optional :: error
    
    logical :: res

    INIT_ERROR(error)

    res = assign_pointer(this%properties, name, ptr)
    if (.not. res) then
      RAISE_ERROR("atoms_assign_prop_ptr_real1d failed to assign pointer to "//trim(name)//" in this%properties", error)
    endif
  end subroutine atoms_assign_prop_ptr_real1D

  subroutine atoms_assign_prop_ptr_real2D(this, name, ptr, error)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    real(dp), pointer :: ptr(:,:)
    integer, intent(out), optional :: error
    
    logical :: res

    INIT_ERROR(error)

    res = assign_pointer(this%properties, name, ptr)
    if (.not. res) then
      RAISE_ERROR("atoms_assign_prop_ptr_real2d failed to assign pointer to "//trim(name)//" in this%properties", error)
    endif
  end subroutine atoms_assign_prop_ptr_real2D

  subroutine atoms_assign_prop_ptr_str(this, name, ptr, error)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    character(1), pointer :: ptr(:,:)
    integer, intent(out), optional :: error
    
    logical :: res

    INIT_ERROR(error)

    res = assign_pointer(this%properties, name, ptr)
    if (.not. res) then
      RAISE_ERROR("atoms_assign_prop_ptr_str failed to assign pointer to "//trim(name)//" in this%properties", error)
    endif
  end subroutine atoms_assign_prop_ptr_str

  subroutine atoms_assign_prop_ptr_logical(this, name, ptr, error)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    logical, pointer :: ptr(:)
    integer, intent(out), optional :: error
    
    logical :: res

    INIT_ERROR(error)

    res = assign_pointer(this%properties, name, ptr)
    if (.not. res) then
      RAISE_ERROR("atoms_assign_prop_ptr_logical failed to assign pointer to "//trim(name)//" in this%properties", error)
    endif
  end subroutine atoms_assign_prop_ptr_logical



  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Atoms indepedent convenince stuff
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !Bond_Length
  !% Returns the sum of the covalent radii of two atoms
  function bond_length(z1,z2)
    integer, intent(in) :: z1,z2
    real(dp)            :: bond_length
    bond_length = ElementCovRad(z1) + ElementCovRad(z2)
  end function bond_length

  subroutine vec_map_into_cell(pos, lattice, g, shift, mapped)
    real(dp), intent(inout) :: pos(3)
    real(dp), intent(in) :: lattice(3,3), g(3,3)
    integer, intent(out), optional :: shift(3)
    logical, intent(out), optional :: mapped

    integer n, k
    real(dp) :: lattice_coord(3)
    logical :: my_mapped

   lattice_coord = g .mult. pos(:)
   my_mapped = .false.
   if (present(shift)) shift = 0

   do n=1,3
      if ((lattice_coord(n) < -0.5_dp) .or. (lattice_coord(n) >= 0.5_dp)) then
	 k = floor(lattice_coord(n)+0.5_dp)
	 lattice_coord(n) = lattice_coord(n) - k
	 if (present(shift)) shift(n) = -k
	 my_mapped = .true.
      end if
   end do

   ! if this atom has been mapped then recalculate its position and shifts for its neighbours
   if (my_mapped) then
      pos(:) = lattice .mult. lattice_coord
   end if
   if (present(mapped)) mapped = my_mapped
  end subroutine vec_map_into_cell


  !% Returns the (unsigned) volume of the simulation cell of lattice
  function lattice_cell_volume(lattice)
    real(dp), intent(in) :: lattice(3,3)
    real(dp)                :: lattice_Cell_Volume
 
    real(dp), dimension(3)  :: a, b, c

    a = lattice(:,1)
    b = lattice(:,2)
    c = lattice(:,3)

    lattice_Cell_Volume = abs(Scalar_Triple_Product(a,b,c))

  end function lattice_cell_volume

endmodule Atoms_types_module
