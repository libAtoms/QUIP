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
!X 'Header' module that contains the data structures for
!X    1. Atoms
!X    2. Connection
!X    3. DomainDecomposition
!X plus methods for manipulating properties.
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

  implicit none
  private

  public :: DEFAULT_NNEIGHTOL
  real(dp), parameter :: DEFAULT_NNEIGHTOL = 1.2_dp    !% Default value for 'atoms%nneightol'

  public :: DD_WRAP_TO_CELL, DD_WRAP_TO_DOMAIN

  integer, parameter  :: DD_WRAP_TO_CELL    = 1   !% All particles, including ghosts, are wrapped into the cell
  integer, parameter  :: DD_WRAP_TO_DOMAIN  = 2   !% Particles are wrapped into the domain, ghost particles are 
                                                  !% located next to the domain.

  public :: Table_pointer
  type Table_pointer
    type(Table), pointer :: t => null()
  end type Table_pointer

  public :: Connection
  type Connection
     !% The Connection type stores the topology of a set of Atoms
     !%
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
     !% :class:`Table` :attr:`neighbour1` (i): $i \le j$ for all $j$ in table, ``intsize=4``, ``realsize=1``
     !%
     !% 'connect%neighbour1(i)%int'
     !%
     !% +----------+----------+----------+----------+
     !% |    1     |    2     |    3     |    4     |
     !% +----------+----------+----------+----------+
     !% |    j     | shift_a  | shift_b  | shift_c  |
     !% +----------+----------+----------+----------+
     !%      
     !% 'connect%neighbour1(i)%real'
     !%
     !% +----------+
     !% |    1     |
     !% +----------+
     !% |  r_ij    |
     !% +----------+
     !%
     !% :class:`Table` :attr:`neighbour2` (i): $i > j$ for all $j$ in table, ``intsize =2``, ``realsize=0``
     !%
     !% 'connect%neighbour2(i)%int'
     !%
     !% +----------+----------+
     !% |    1     |    2     |
     !% +----------+----------+
     !% |    j     |    n     |
     !% +----------+----------+
     !%      
     !% :class:`Table` :attr:`cell` (i,j,k) with ``intsize = 1``, ``realsize = 0``
     !%
     !% 'connect%cell(i,j,k)2%int'
     !%
     !% +----------+
     !% |    1     |
     !% +----------+
     !% |  atom    |
     !% +----------+
     !%
     !% N.B. If $i$ and $j$ are neighbours with shift 'shift', then
     !% ``norm(atoms%pos(j) - atoms%pos(i) + shift)`` is a minimum.
     !% Mnemonic: 'shift' is added to $j$ to get closer to $i$.

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

     integer, allocatable, dimension(:,:,:) :: cell_heads          !% First entry in cell atoms structure
     integer, allocatable, dimension(:)     :: next_atom_in_cell   !% List of atoms, terminated by zero

     logical, allocatable, dimension(:) :: is_min_image   !% True if i is a minimum image

  end type Connection


  public :: DomainDecomposition
  type DomainDecomposition

     integer                    :: Ntotal                = 0        !% Number of total particles in this simulation
     integer, pointer           :: local_to_global(:)   => NULL()   !% Local index to global index
     integer, allocatable       :: global_to_local(:)               !% Global index to local index

     integer                    :: decomposition(3)      = (/ 1, 1, 1 /)   !% Type of decomposition

     integer                    :: mode                  = DD_WRAP_TO_CELL

     logical                    :: decomposed            = .false.   !% True if domain decomposition is active

     real(dp)                   :: requested_border      = 0.0_DP
     real(dp)                   :: border(3)             = 0.0_DP
     real(dp)                   :: verlet_shell          = 0.0_DP

     logical                    :: communicate_forces    = .false.

     real(dp)                   :: lower(3)              = 0.0_DP   !% Lower domain boundary, in fraction of the total cell
     real(dp)                   :: upper(3)              = 0.0_DP   !% Upper domain boundary, in fraction of the total cell
     real(dp)                   :: center(3)             = 0.0_DP   !% Center of the domain
     real(dp)                   :: lower_with_border(3)  = 0.0_DP   !% Lower domain boundary, including border
     real(dp)                   :: upper_with_border(3)  = 0.0_DP   !% Upper domain boundary, including border

     type(MPI_context)          :: mpi                              !% MPI communicator
     logical                    :: periodic(3)           = .true.   !% Periodicity for domain decomposition
     integer                    :: l(3)                  = 0        !% Ranks of left domains in x-, y- and z-direction
     integer                    :: r(3)                  = 0        !% Ranks of right domains in x-, y- and z-direction
     real(dp)                   :: off_l(3)              = 0.0_DP   !% Distance vector to left domain
     real(dp)                   :: off_r(3)              = 0.0_DP   !% Distance vector to right domain

     type(Dictionary)           :: atoms_properties     !% Fields to communicate if for particles
     type(Dictionary)           :: ghost_properties     !% Fields to communicate for ghosts
     type(Dictionary)           :: reverse_properties   !% Back-communication after force computations     

     logical, allocatable       :: atoms_mask(:)
     logical, allocatable       :: ghost_mask(:)
     logical, allocatable       :: reverse_mask(:)

     integer                    :: atoms_buffer_size    = 0
     integer                    :: ghost_buffer_size    = 0
     integer                    :: reverse_buffer_size  = 0

     character(1), allocatable  :: send_l(:)   !% buffer for sending to the left
     character(1), allocatable  :: send_r(:)   !% buffer for sending to the right
     character(1), allocatable  :: recv_l(:)   !% buffer for receiving from the left
     character(1), allocatable  :: recv_r(:)   !% buffer for receiving from the right

     integer                    :: n_ghosts_r(3)  = 0   !% length of the ghost particle lists (right)
     integer                    :: n_ghosts_l(3)  = 0   !% length of the ghost particle lists (left)

     integer, allocatable       :: ghosts_r(:)   !% particles send to the right (where they become ghosts)
     integer, allocatable       :: ghosts_l(:)   !% particles send to the left (where they become ghosts)

     integer                    :: n_send_p_tot   !% Statistics: Number of total particles send
     integer                    :: n_recv_p_tot   !% Statistics: Number of total particles received
     integer                    :: n_send_g_tot   !% Statistics: Number of total ghosts send
     integer                    :: n_recv_g_tot   !% Statistics: Number of total ghosts received
     integer                    :: nit_p          !% Statistics: Number of particle send events
     integer                    :: nit_g          !% Statistics: Number of ghost send events
     
  endtype DomainDecomposition


  public :: Atoms
  type Atoms
     !% Representation of an atomic configuration and its associated properties
     !%
     !% An atoms object contains atomic numbers, all dynamical variables
     !% and connectivity information for all the atoms in the simulation cell. 
     !% It is initialised like this:
     !%> 	  call initialise(MyAtoms,N,lattice)
     !% where 'N' is the number of atoms to allocate space for and 'lattice' is a $3\times3$
     !% matrix of lattice vectors given as column vectors, so that 'lattice(:,i)' is the i-th lattice vector.
     !% 
     !% Atoms also contains a Connection object, which stores distance information about
     !% the atom neghbours after 'calc_connect' has been called. Rather than using a minimum
     !% image convention, all neighbours are stored up to a radius of 'cutoff', including images

     ! Self-deallocating object
     logical                               :: own_this = .false.  !% Do I own myself?
     integer                               :: ref_count = 0  !% Reference counter

     logical                               :: fixed_size = .false. !% Can the number of atoms be changed after initialisation?
     integer                               :: N = 0 !% The number of atoms held (including ghost particles)
     integer                               :: Ndomain = 0 !% The number of atoms held by the local process (excluding ghost particles)
     integer                               :: Nbuffer = 0 !% The number of atoms that can be stored in the buffers of this Atoms object

     logical                               :: use_uniform_cutoff = .false. !% Rather than covalent radii --- 
                                                                           !% default is variable cutoff.
     real(dp)                              :: cutoff = DEFAULT_NNEIGHTOL   !% if 'use_uniform_cutoff' is true, cutoff
                                                                           !% is the cutoff distance in \AA{}.
                                                                           !% Otherwise, cutoff is a multiplier
                                                                           !% for 'bond_length(Zi,Zj)'.
     real(dp)                              ::  cutoff_break = DEFAULT_NNEIGHTOL  !% Cutoff length for bonds to be considered broken with hysteretic connectivity
     real(dp)                              :: nneightol = DEFAULT_NNEIGHTOL 
                                              !% Count as nearest neighbour if sum of covalent radii
                                              !% times 'this%nneightol' greater than distance between atoms.
                                              !% Used in cluster carving.

     real(dp),              dimension(3,3) :: lattice    !% Lattice vectors, as columns:
     !%\begin{displaymath}
     !%\left(
     !%\begin{array}{ccc}
     !% | & | & | \\ \mathbf{a} & \mathbf{b} & \mathbf{c} \\ | & | & | \\ \end{array}
     !%\right)
     !% = \left(
     !%\begin{array}{ccc}
     !% R_{11} & R_{12} & R_{13} \\ R_{21} & R_{22} & R_{23} \\  R_{31} & R_{32} & R_{33} \\ \end{array}
     !%\right)
     !%\end{displaymath}
     !% i.e. $\mathbf{a}$ = 'lattice(:,1)', $\mathbf{b}$ = 'lattice(:,2)' and
     !% $\mathbf{c}$ 'lattice(:,3)'.

     !  ( | | | | | | )          ( (1,1) (1,2) (1,3) ) 
     ! (  | | | | | |  )        (                     )
     ! (  |a| |b| |c|  )    =   (  (2,1) (2,2) (2,3)  )
     ! (  | | | | | |  )        (                     )
     !  ( | | | | | | )          ( (3,1) (3,2) (3,3) ) 
     logical :: is_orthorhombic, is_periodic(3)


     real(dp),              dimension(3,3) :: g          !% Inverse lattice (stored for speed)

     type(Dictionary) :: properties !% Dictionary of atomic properties. A property is an array
                                    !% of shape (`m`,`n`) where `n` is the number of atoms and `m` is
                                    !% either one (for scalar properties) or three (vector
                                    !% properties). Properties can be integer, real, string or logical.
                                    !% String properties have a fixed length of ``TABLE_STRING_LENGTH=10``
                                    !% characters.
                                    !%
                                    !% From Fortran, the following default properties are aliased with
                                    !% arrays within the Atoms type:
                                    !%
                                    !%  * ``Z`` - Atomic numbers, dimension is actually $(N)$
                                    !%  * ``species`` Names of elements
                                    !%  * ``move_mask`` Atoms with 'move_mask' set to zero are fixed
                                    !%  * ``damp_mask`` Damping is only applied to those atoms with 'damp_mask' set to 1. By default this is set to 1 for all atoms.
                                    !%  * ``thermostat_region`` Which thermostat is applied to each atoms. By default this is set to 1 for all atoms.
                                    !%  * ``travel`` Travel across periodic conditions. $(3,N)$ integer array. See meth:`map_into_cell` below.
                                    !%  * ``pos`` $(3,N)$ array of atomic positions, in $\mathrm{\AA}$. Position of atom $i$ is 'pos(:,i)'
                                    !%  * ``mass`` Atomic masses, dimension is $(N)$
                                    !%  * ``velo`` $(3,N)$ array  of atomic velocities, in $\mathrm{AA}$/fs.
                                    !%  * ``acc`` $(3,N)$ array  of accelerations in $\mathrm{AA}$/fs$^2$
                                    !%  * ``avgpos`` $(3,N)$ array  of time-averaged atomic positions.
                                    !%  * ``oldpos`` $(3,N)$ array  of positions of atoms at previous time step.
                                    !%  * ``avg_ke`` Time-averaged atomic kinetic energy
                                    !% 
                                    !% Custom properties are most conveniently accessed by assign a pointer to
                                    !% them with the :meth:`assign_pointer` routines.
                                    !% 
                                    !% From Python, each property is automatically visible as a
                                    !% array attribute of the :class:`Atoms` object,
                                    !% for example the atomic positions are stored in a real vector
                                    !% property called `pos`, and can be accessed as ``at.pos``.
                                    !% 
                                    !% Properties can be added with the :meth:`add_property` method and
                                    !% removed with :meth:`remove_property`.

     type(Dictionary) :: params     !% Dictionary of parameters. Useful for storing data about this
                                    !% Atoms object, for example the temperature, total energy or
                                    !% applied strain. The data stored here is automatically saved to
                                    !% and loaded from XYZ and NetCDF files.

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
     real(dp), pointer, dimension(:,:) :: pos    => null()  !% $(3,N)$ array of atomic positions, in $\mathrm{AA}$. 
                                                            !% Position of atom $i$ is 'pos(:,i)'
     real(dp), pointer, dimension(:)   :: mass   => null()  !% Atomic masses, dimension is actually $(N)$

     real(dp), pointer, dimension(:,:) :: velo   => null()  !% $(3,N)$ array  of atomic velocities, in $\mathrm{AA}$/fs.
     real(dp), pointer, dimension(:,:) :: acc    => null()  !% $(3,N)$ array  of accelerations in $\mathrm{AA}$/fs$^2$
     real(dp), pointer, dimension(:,:) :: avgpos => null()  !% $(3,N)$ array  of time-averaged atomic positions.
     real(dp), pointer, dimension(:,:) :: oldpos => null()  !% $(3,N)$ array  of positions of atoms at previous time step.
     real(dp), pointer, dimension(:)   :: avg_ke => null()    !% Time-averaged atomic kinetic energy

     type(Connection)                  :: connect             !% Connectivity object
     type(Connection)                  :: hysteretic_connect  !% Hysteretic connectivity object

     type(DomainDecomposition)         :: domain              !% Domain decomposition object

  end type Atoms


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

  
  !% Copy an atom to a different index
  public :: copy_entry
  interface copy_entry
     module procedure atoms_copy_entry
  endinterface copy_entry


  ! Note: These are Atoms unrelated, and could be easily move somewhere else
  ! (linearalgebra?)
  public :: map_into_cell
  interface map_into_cell
    module procedure vec_map_into_cell, array_map_into_cell
  end interface

  public :: cell_volume
  interface cell_volume
    module procedure lattice_cell_volume
  end interface

  public :: atoms_repoint, atoms_sort, bond_length

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
             if (this%properties%entries(i)%type /= T_INTEGER_A) &
                  call system_abort('Confused by Atoms property "Z" not of type T_INTEGER_A')
             this%Z               => this%properties%entries(i)%i_a(:)
          case('travel')
             if (this%properties%entries(i)%type /= T_INTEGER_A2) &
                  call system_abort('Confused by Atoms property "travel" not of type T_INTEGER_A2')
             this%travel          => this%properties%entries(i)%i_a2(:,:)
          case('move_mask')
             if (this%properties%entries(i)%type /= T_INTEGER_A) &
                  call system_abort('Confused by Atoms property "move_mask" not of type T_INTEGER_A')
             this%move_mask       => this%properties%entries(i)%i_a(:)
          case('damp_mask')
             if (this%properties%entries(i)%type /= T_INTEGER_A) &
                  call system_abort('Confused by Atoms property "damp_mask" not of type T_INTEGER_A')
             this%damp_mask       => this%properties%entries(i)%i_a(:)
          case('thermostat_region')
             if (this%properties%entries(i)%type /= T_INTEGER_A) &
                  call system_abort('Confused by Atoms property "thermostat_mask" not of type T_INTEGER_A')
             this%thermostat_region => this%properties%entries(i)%i_a(:)

          ! Real properties
          case('mass')
             if (this%properties%entries(i)%type /= T_REAL_A) &
                  call system_abort('Confused by Atoms property "mass" not of type T_REAL_A')
             this%mass            => this%properties%entries(i)%r_a(:)
          case('pos')
             if (this%properties%entries(i)%type /= T_REAL_A2) &
                  call system_abort('Confused by Atoms property "mass" not of type T_REAL_A2')
             this%pos             => this%properties%entries(i)%r_a2(:,:)
          case('velo')
             if (this%properties%entries(i)%type /= T_REAL_A2) &
                  call system_abort('Confused by Atoms property "velo" not of type T_REAL_A2')
             this%velo            => this%properties%entries(i)%r_a2(:,:)
          case('acc')
             if (this%properties%entries(i)%type /= T_REAL_A2) &
                  call system_abort('Confused by Atoms property "acc" not of type T_REAL_A2')
             this%acc             => this%properties%entries(i)%r_a2(:,:)
          case('avgpos')
             if (this%properties%entries(i)%type /= T_REAL_A2) &
                  call system_abort('Confused by Atoms property "avgpos" not of type T_REAL_A2')
             this%avgpos          => this%properties%entries(i)%r_a2(:,:)
          case('oldpos')
             if (this%properties%entries(i)%type /= T_REAL_A2) &
                  call system_abort('Confused by Atoms property "oldpos" not of type T_REAL_A2')
             this%oldpos          => this%properties%entries(i)%r_a2(:,:)
          case('avg_ke')
             if (this%properties%entries(i)%type /= T_REAL_A) &
                  call system_abort('Confused by Atoms property "avg_ke" not of type T_REAL_A')
             this%avg_ke          => this%properties%entries(i)%r_a(:)

          ! String properties
          case('species')
             if (this%properties%entries(i)%type /= T_CHAR_A) &
                  call system_abort('Confused by Atoms property "species" not of type T_CHAR_A')
             this%species         => this%properties%entries(i)%s_a(:,:)

          end select


       else

          select case(trim(lower_case(string(this%properties%keys(i)))))

          ! Integer properties
          case('z')
             if (this%properties%entries(i)%type /= T_INTEGER_A) &
                  call system_abort('Confused by Atoms property "Z" not of type T_INTEGER_A')
             this%Z               => this%properties%entries(i)%i_a(1:this%n)
          case('travel')
             if (this%properties%entries(i)%type /= T_INTEGER_A2) &
                  call system_abort('Confused by Atoms property "travel" not of type T_INTEGER_A2')
             this%travel          => this%properties%entries(i)%i_a2(:,1:this%n)
          case('move_mask')
             if (this%properties%entries(i)%type /= T_INTEGER_A) &
                  call system_abort('Confused by Atoms property "move_mask" not of type T_INTEGER_A')
             this%move_mask       => this%properties%entries(i)%i_a(1:this%n)
          case('damp_mask')
             if (this%properties%entries(i)%type /= T_INTEGER_A) &
                  call system_abort('Confused by Atoms property "damp_mask" not of type T_INTEGER_A')
             this%damp_mask       => this%properties%entries(i)%i_a(1:this%n)
          case('thermostat_region')
             if (this%properties%entries(i)%type /= T_INTEGER_A) &
                  call system_abort('Confused by Atoms property "thermostat_mask" not of type T_INTEGER_A')
             this%thermostat_region => this%properties%entries(i)%i_a(1:this%n)

          ! Real properties
          case('mass')
             if (this%properties%entries(i)%type /= T_REAL_A) &
                  call system_abort('Confused by Atoms property "mass" not of type T_REAL_A')
             this%mass            => this%properties%entries(i)%r_a(1:this%n)
          case('pos')
             if (this%properties%entries(i)%type /= T_REAL_A2) &
                  call system_abort('Confused by Atoms property "mass" not of type T_REAL_A2')
             this%pos             => this%properties%entries(i)%r_a2(:,1:this%N)
          case('velo')
             if (this%properties%entries(i)%type /= T_REAL_A2) &
                  call system_abort('Confused by Atoms property "velo" not of type T_REAL_A2')
             this%velo            => this%properties%entries(i)%r_a2(:,1:this%N)
          case('acc')
             if (this%properties%entries(i)%type /= T_REAL_A2) &
                  call system_abort('Confused by Atoms property "acc" not of type T_REAL_A2')
             this%acc             => this%properties%entries(i)%r_a2(:,1:this%N)
          case('avgpos')
             if (this%properties%entries(i)%type /= T_REAL_A2) &
                  call system_abort('Confused by Atoms property "avgpos" not of type T_REAL_A2')
             this%avgpos          => this%properties%entries(i)%r_a2(:,1:this%N)
          case('oldpos')
             if (this%properties%entries(i)%type /= T_REAL_A2) &
                  call system_abort('Confused by Atoms property "oldpos" not of type T_REAL_A2')
             this%oldpos          => this%properties%entries(i)%r_a2(:,1:this%N)
          case('avg_ke')
             if (this%properties%entries(i)%type /= T_REAL_A) &
                  call system_abort('Confused by Atoms property "avg_ke" not of type T_REAL_A')
             this%avg_ke          => this%properties%entries(i)%r_a(1:this%n)

          ! String properties
          case('species')
             if (this%properties%entries(i)%type /= T_CHAR_A) &
                  call system_abort('Confused by Atoms property "species" not of type T_CHAR_A')
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

    if (present(ptr) .and. present(ptr2)) then
      RAISE_ERROR('atoms_add_property_real_a got both ptr and ptr2', error)
    endif

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
    
    if (use_n_cols == 1 .and. .not. present(ptr2)) then
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

    if (present(ptr) .and. present(ptr2)) then
      RAISE_ERROR('atoms_add_property_real_a got both ptr and ptr2', error)
    endif

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
    
    if (use_n_cols == 1 .and. .not. present(ptr2)) then
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
  ! Copying individual atoms and sorting
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% Move a single atom from one location to another one.
  !% The destination will be overriden.
  subroutine atoms_copy_entry(this, src, dst, swap, error)

    type(Atoms), intent(inout)  :: this
    integer, intent(in)         :: src
    integer, intent(in)         :: dst
    logical, intent(in), optional :: swap
    integer, optional, intent(out) :: error

    integer i

    logical :: my_swap
    integer :: t_i
    integer, allocatable :: t_i_a(:)
    real(dp) :: t_r
    real(dp), allocatable :: t_r_a(:)
    logical :: t_l
    character(len=1), allocatable :: t_c(:)

    ! ---

    INIT_ERROR(error)
    my_swap = optional_default(.false., swap)

    if (src < 1 .or. src > this%N) then
       RAISE_ERROR('atoms_copy_entry: src='//src//' out of range 1 <= src <= '//this%n, error)
    end if
    if (dst < 1 .or. dst > this%N) then
       RAISE_ERROR('atoms_copy_entry: dst='//dst//' out of range 1 <= dst <= '//this%n, error)
    end if

    do i=1,this%properties%N
       select case (this%properties%entries(i)%type)

       case(T_INTEGER_A)
	  if (my_swap) then
	     t_i = this%properties%entries(i)%i_a(dst)
	     this%properties%entries(i)%i_a(dst) = this%properties%entries(i)%i_a(src)
	     this%properties%entries(i)%i_a(src) = t_i
	  else
	     this%properties%entries(i)%i_a(dst) = this%properties%entries(i)%i_a(src)
	  endif

       case(T_REAL_A)
	  if (my_swap) then
	     t_r = this%properties%entries(i)%r_a(dst)
	     this%properties%entries(i)%r_a(dst) = this%properties%entries(i)%r_a(src)
	     this%properties%entries(i)%r_a(src) = t_r
	  else
	     this%properties%entries(i)%r_a(dst) = this%properties%entries(i)%r_a(src)
	  endif

       case(T_LOGICAL_A)
	  if (my_swap) then
	     t_l = this%properties%entries(i)%l_a(dst)
	     this%properties%entries(i)%l_a(dst) = this%properties%entries(i)%l_a(src)
	     this%properties%entries(i)%l_a(src) = t_l
	  else
	     this%properties%entries(i)%l_a(dst) = this%properties%entries(i)%l_a(src)
	  endif

       case(T_INTEGER_A2)
	  if (my_swap) then
	     allocate(t_i_a(size(this%properties%entries(i)%i_a2,1)))
	     t_i_a = this%properties%entries(i)%i_a2(:,dst)
	     this%properties%entries(i)%i_a2(:,dst) = this%properties%entries(i)%i_a2(:,src)
	     this%properties%entries(i)%i_a2(:,src) = t_i_a
	     deallocate(t_i_a)
	  else
	     this%properties%entries(i)%i_a2(:,dst) = this%properties%entries(i)%i_a2(:,src)
	  endif

       case(T_REAL_A2)
	  if (my_swap) then
	     allocate(t_r_a(size(this%properties%entries(i)%r_a2,1)))
	     t_r_a = this%properties%entries(i)%r_a2(:,dst)
	     this%properties%entries(i)%r_a2(:,dst) = this%properties%entries(i)%r_a2(:,src)
	     this%properties%entries(i)%r_a2(:,src) = t_r_a
	     deallocate(t_r_a)
	  else
	     this%properties%entries(i)%r_a2(:,dst) = this%properties%entries(i)%r_a2(:,src)
	  endif

       case(T_CHAR_A)
	  if (my_swap) then
	     allocate(t_c(size(this%properties%entries(i)%s_a,1)))
	     t_c = this%properties%entries(i)%s_a(:,dst)
	     this%properties%entries(i)%s_a(:,dst) = this%properties%entries(i)%s_a(:,src)
	     this%properties%entries(i)%s_a(:,src) = t_c
	     deallocate(t_c)
	  else
	     this%properties%entries(i)%s_a(:,dst) = this%properties%entries(i)%s_a(:,src)
	  endif

       case default
          RAISE_ERROR('atoms_copy_entry: bad property type '//this%properties%entries(i)%type//' key='//this%properties%keys(i), error)

       end select
    end do

  endsubroutine atoms_copy_entry


  !% sort atoms by one or more (max 2 now) integer or real properties
  subroutine atoms_sort(this, prop1, prop2, prop3, error)
    type(Atoms), intent(inout) :: this
    character(len=*), intent(in) :: prop1
    character(len=*), intent(in), optional :: prop2, prop3
    integer, intent(out), optional :: error

    integer, pointer :: i_p1(:) => null(), i_p2(:) => null(), i_p3(:) => null()
    real(dp), pointer :: r_p1(:) => null(), r_p2(:) => null(), r_p3(:) => null()
    integer :: cur_place, i_a, smallest_i_a
    logical :: is_lt

    INIT_ERROR(error)

    if (.not. assign_pointer(this, prop1, i_p1)) then
       if (.not. assign_pointer(this, prop1, r_p1)) then
	  RAISE_ERROR("atoms_sort can't find 1st integer or real property '" // prop1 //"'", error)
       endif
    endif
    if (present(prop2)) then
       if (.not. assign_pointer(this, prop2, i_p2)) then
	  if (.not. assign_pointer(this, prop2, r_p2)) then
	     RAISE_ERROR("atoms_sort can't find 2nd integer or real property '" // prop2 //"'", error)
	  endif
       endif
    endif
    if (present(prop3)) then
       if (.not. assign_pointer(this, prop3, i_p3)) then
	  if (.not. assign_pointer(this, prop3, r_p3)) then
	     RAISE_ERROR("atoms_sort can't find 3rd integer or real property '" // prop3 //"'", error)
	  endif
       endif
    endif

    do cur_place=1, this%N-1
       smallest_i_a = cur_place
       do i_a = cur_place+1, this%N
	  is_lt = arrays_lt(i_a, smallest_i_a, i_p1=i_p1, r_p1=r_p1, i_p2=i_p2, r_p2=r_p2, i_p3=i_p3, r_p3=r_p3, error=error)
	  PASS_ERROR(error)
	  if (is_lt) then
	     smallest_i_a = i_a
	  endif
       end do
       if (smallest_i_a /= cur_place) then
	  call copy_entry(this, cur_place, smallest_i_a, swap=.true., error=error)
	  PASS_ERROR(error)
       endif
    end do
  end subroutine atoms_sort



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

  subroutine array_map_into_cell(pos, lattice, g)
    real(dp), intent(inout) :: pos(:,:)
    real(dp), intent(in) :: lattice(3,3), g(3,3)

    integer :: i

    do i=1, size(pos, 2)
      call map_into_cell(pos(:,i), lattice, g)
    end do
  end subroutine array_map_into_cell

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
