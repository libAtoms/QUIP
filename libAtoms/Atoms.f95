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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X  Atoms module
!X
!% An atoms object contains atomic numbers, all dynamical variables
!% and connectivity information for all the atoms in the simulation cell. 
!% It is initialised like this:
!%> 	  call initialise(MyAtoms,N,lattice)
!% where 'N' is the number of atoms to allocate space for and 'lattice' is a $3\times3$
!% matrix of lattice vectors given as column vectors, so that lattice(:,i) is the i-th lattice vector.
!% 
!% Atoms also contains a Connection object, which stores distance information about
!% the atom neghbours after 'calc_connect' has been called. Rather than using a minimum
!% image convention, all neighbours are stored up to a radius of 'cutoff', including images
!% 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module  atoms_module

  use system_module
  use extendable_str_module
  use error_module
  use linearalgebra_module
  use table_module
  use dictionary_module
  use periodictable_module
  use minimization_module
  use mpi_context_module

  implicit none

  real(dp), parameter :: DEFAULT_NNEIGHTOL = 1.2_dp    !% Default value for 'atoms%nneightol'

  integer,  parameter :: NOT_NEIGHBOUR = 0     !% Returned by 'Find_Neighbour' if unsuccessful

  logical :: printed_stack_warning = .false.   !% OMIT

  private :: ess_max_len

  type Table_pointer
    type(Table), pointer :: t => null()
  end type Table_pointer

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

     type(table_pointer), allocatable, dimension(:)     :: neighbour1 !% Neighbour information for pairs $i <= j$. 
                                                              !% Contains full details of $j$, $r_{ij}$ and shift.
     type(table_pointer), allocatable, dimension(:)     :: neighbour2 !% Neighbour information for pairs $i > j$.
                                                              !% Simply contains $j$ and a reference to $j$'s
                                                              !% 'neighbour1' table.

     type(table), allocatable, dimension(:,:,:) :: cell    !% For the linear scaling connection calculator

     logical, allocatable, dimension(:) :: is_min_image   !% True if i is a minimum image

  end type Connection


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

  private :: atoms_initialise, connection_initialise
  interface initialise
     module procedure atoms_initialise, connection_initialise
  end interface initialise

  !% Initialise type(Atoms), pointer objects. Shallow copies of these will
  !% survive even if the initial declaration goes out of scope. The object will
  !% automatically deallocate upon calling finalise_ptr when the last shallow
  !% copy goes out of scope
  private :: atoms_initialise_ptr
  interface initialise_ptr
     module procedure atoms_initialise_ptr
  endinterface initialise_ptr

  private :: atoms_is_initialised
  interface is_initialised
     module procedure atoms_is_initialised
  end interface is_initialised

  private :: atoms_shallowcopy
  interface shallowcopy
     module procedure atoms_shallowcopy
  end interface shallowcopy

  !% Free up the memory associated with one or more objects.
  private :: atoms_finalise,atoms_finalise_multi, connection_finalise
  interface finalise
     module procedure atoms_finalise, atoms_finalise_multi, connection_finalise
  end interface finalise

  !% Finalise type(Atoms), pointer objects. Shallow copies of these will
  !% The object will when ref_count == 0.
  private :: atoms_finalise_ptr
  interface finalise_ptr
     module procedure atoms_finalise_ptr
  endinterface finalise_ptr

  private :: connection_wipe
  interface wipe
     module procedure connection_wipe
  end interface wipe

  private :: atoms_zero
  interface zero
     module procedure atoms_zero
  end interface

  !% Overloaded assigment operators for Atoms and Connection objects.
  private :: atoms_assignment, connection_assignment
  interface assignment(=)
     module procedure atoms_assignment, connection_assignment
  end interface assignment(=)

  private :: atoms_set_atoms, atoms_set_atoms_singlez
  interface set_atoms
     module procedure atoms_set_atoms, atoms_set_atoms_singlez
  end interface

  !% increase cutoff
  private :: atoms_set_cutoff_minimum
  interface set_cutoff_minimum
     module procedure atoms_set_cutoff_minimum
  end interface

  !% set (a uniform) cutoff
  private :: atoms_set_cutoff
  interface set_cutoff
     module procedure atoms_set_cutoff
  end interface

  !% set cutoff factor
  private :: atoms_set_cutoff_factor
  interface set_cutoff_factor
     module procedure atoms_set_cutoff_factor
  end interface

  !% Add one or more atoms to an Atoms object.
  private :: add_atom_single, add_atom_multiple
  interface add_atoms
     module procedure add_atom_single, add_atom_multiple
  end interface add_atoms

  !% Remove one or more atoms from an Atoms object.
  private :: remove_atom_single, remove_atom_multiple
  interface remove_atoms
     module procedure remove_atom_single, remove_atom_multiple
  end interface remove_atoms

  !% Add a per-atom property to this atoms object, as extra entry with columns of 
  !% integers, reals, logical, or strings in the 'properties' dictionary. For example, 
  !% this interface is used by the DynamicalSystems module to create the 'velo', 'acc', 
  !% etc. properties.
  !% Optionally, a pointer to the new property is returned.
  private :: atoms_add_property_int, atoms_add_property_int_a
  private :: atoms_add_property_real, atoms_add_property_real_a
  private :: atoms_add_property_str, atoms_add_property_str_2da
  private :: atoms_add_property_str_a
  private :: atoms_add_property_logical, atoms_add_property_logical_a 
  private :: atoms_add_property_int_2Da
  private :: atoms_add_property_real_2Da
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
  private :: atoms_add_property_p_int, atoms_add_property_p_int_a
  private :: atoms_add_property_p_real, atoms_add_property_p_real_a
  private :: atoms_add_property_p_str
  private :: atoms_add_property_p_logical
  interface add_property_from_pointer
     module procedure atoms_add_property_p_int, atoms_add_property_p_int_a
     module procedure atoms_add_property_p_real, atoms_add_property_p_real_a
     module procedure atoms_add_property_p_str
     module procedure atoms_add_property_p_logical
  end interface

  !% get a (per-configuration) value from the atoms%params dictionary
  private :: atoms_get_param_value_int, atoms_get_param_value_int_a
  private :: atoms_get_param_value_real, atoms_get_param_value_real_a, atoms_get_param_value_real_a2
  private :: atoms_get_param_value_str, atoms_get_param_value_es
  private :: atoms_get_param_value_logical
  interface get_param_value
     module procedure atoms_get_param_value_int, atoms_get_param_value_int_a
     module procedure atoms_get_param_value_real, atoms_get_param_value_real_a, atoms_get_param_value_real_a2
     module procedure atoms_get_param_value_str, atoms_get_param_value_es
     module procedure atoms_get_param_value_logical
  end interface get_param_value

  !% set a (per-configuration) value from the atoms%params dictionary
  private :: atoms_set_param_value_int, atoms_set_param_value_int_a
  private :: atoms_set_param_value_real, atoms_set_param_value_real_a, atoms_set_param_value_real_a2
  private :: atoms_set_param_value_str
  private :: atoms_set_param_value_logical
  interface set_param_value
     module procedure atoms_set_param_value_int, atoms_set_param_value_int_a
     module procedure atoms_set_param_value_real, atoms_set_param_value_real_a, atoms_set_param_value_real_a2
     module procedure atoms_set_param_value_str
     module procedure atoms_set_param_value_logical
  end interface set_param_value

  !% Convenience function to test if a property is present. No checking
  !% of property type is done.
  private :: atoms_has_property
  interface has_property
     module procedure atoms_has_property
  end interface

  !% remove a property from this atoms object
  private :: atoms_remove_property
  interface remove_property
     module procedure atoms_remove_property
  end interface

  !% This is a convenience interface to assign pointers to custom properties of this
  !% Atoms object. The pointer is simply directed to the relevant array in the
  !% 'this%properties' Dictionary. Adding or removing atoms will invalidate these pointers.
  !% Returns true if successful, false if property doesn't exist or
  !% if property type and type of pointer don't match.
  !% OUTDATED - replace by subroutine assign_property_pointer with error handling
  private :: atoms_assign_pointer_int1D, atoms_assign_pointer_int2D
  private :: atoms_assign_pointer_real1D, atoms_assign_pointer_real2D
  private :: atoms_assign_pointer_str
  private :: atoms_assign_pointer_logical
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
  private :: atoms_assign_prop_ptr_int1D, atoms_assign_prop_ptr_int2D
  private :: atoms_assign_prop_ptr_real1D, atoms_assign_prop_ptr_real2D
  private :: atoms_assign_prop_ptr_str
  private :: atoms_assign_prop_ptr_logical
  interface assign_property_pointer
     module procedure atoms_assign_prop_ptr_int1D, atoms_assign_prop_ptr_int2D
     module procedure atoms_assign_prop_ptr_real1D, atoms_assign_prop_ptr_real2D
     module procedure atoms_assign_prop_ptr_str
     module procedure atoms_assign_prop_ptr_logical
  end interface

  !% This interface calculates the distance between the nearest periodic images of two points (or atoms).
  private :: distance8_atom_atom, distance8_atom_vec, distance8_vec_atom, distance8_vec_vec
  interface distance_min_image
     module procedure distance8_atom_atom, distance8_atom_vec, distance8_vec_atom, distance8_vec_vec
  end interface

  !% This interface calculates the difference vector between the nearest periodic images of two points (or atoms).
  private :: diff_atom_atom, diff_atom_vec, diff_vec_atom, diff_vec_vec
  interface diff_min_image
     module procedure diff_atom_atom, diff_atom_vec, diff_vec_atom, diff_vec_vec
  end interface

  !% Print a verbose textual description of an Atoms or Connection object to the default logger or to
  !% a specificied Inoutput object.
  private :: atoms_print, connection_print
  interface print
     module procedure atoms_print, connection_print
  end interface print

  !% set the lattice of an atoms object - please use this, rather than setting atoms%lattice
  !% directly, because it also set up reciprocal lattice and orthorhombic/periodic logical flags
  private :: atoms_set_lattice
  interface set_lattice
     module procedure atoms_set_lattice
  end interface set_lattice

  !% Select a subset of the atoms in an atoms object, either using a logical 
  !% mask array or a Table, in which case the first 'int' column is taken to 
  !% be a list of the atoms that should be retained.
  private :: atoms_select
  interface select
     module procedure atoms_select
  end interface select

  !% calculate volume of unit cell
  private :: atoms_cell_volume
  private :: lattice_cell_volume
  interface cell_volume
    module procedure atoms_cell_volume
    module procedure lattice_cell_volume
  end interface

  private :: atoms_map_into_cell
  private :: vec_map_into_cell
  interface map_into_cell
    module procedure atoms_map_into_cell
    module procedure vec_map_into_cell
  end interface

  private :: atoms_bcast
  interface bcast
     module procedure atoms_bcast
  end interface

  private :: atoms_copy_properties
  interface copy_properties
     module procedure atoms_copy_properties
  endinterface copy_properties

  private :: atoms_copy_entry
  interface copy_entry
     module procedure atoms_copy_entry
  endinterface copy_entry

  private :: atoms_index_to_z_index
  interface index_to_z_index
     module procedure atoms_index_to_z_index
  end interface index_to_z_index

  private :: atoms_z_index_to_index
  interface z_index_to_index
     module procedure atoms_z_index_to_index
  end interface z_index_to_index

contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Initialisation and Finalisation
  !
  ! N : number of atoms
  ! lattice: 3x3 matrix of lattice vectors as column vectors
  ! data: optional initial values for data table
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% Initialise an Atoms object to store 'N' atoms and specify the initial lattice.
  subroutine atoms_initialise(this,N,lattice,&
       properties,params,fixed_size,Nbuffer,error)

    type(Atoms),                intent(inout) :: this
    integer,                    intent(in)    :: N
    real(dp), dimension(3,3),   intent(in)    :: lattice
    type(Dictionary), optional, intent(in)    :: properties, params
    logical,          optional, intent(in)    :: fixed_size
    integer,          optional, intent(in)    :: Nbuffer
    integer,          optional, intent(out)   :: error

    integer :: stack_size, stack_size_err, i, type

    INIT_ERROR(error)
    call atoms_finalise(this)

    if (present(params)) then
       this%params = params ! deep copy
    else
       call initialise(this%params)
    end if

    this%N = N
    this%Ndomain = N
    this%Nbuffer = N

    if (present(Nbuffer))  this%Nbuffer = Nbuffer

    this%fixed_size = optional_default(.false., fixed_size)

    if (present(properties)) then
       ! check types and sizes of properties passed in
       do i=1,properties%N
          type = properties%entries(i)%type
          if (type == T_INTEGER_A .or. type == T_REAL_A .or. type == T_LOGICAL_A) then
             if (properties%entries(i)%len /= this%Nbuffer) then
                RAISE_ERROR('atoms_initialise: bad array size for array key='//string(properties%keys(i))//' size('//properties%entries(i)%len//') != this%Nbuffer('//this%Nbuffer//')', error)
             end if
          else if (type == T_INTEGER_A2 .or. type == T_REAL_A2 .or. type == T_CHAR_A) then
             if (properties%entries(i)%len2(2) /= this%Nbuffer) then
                RAISE_ERROR('atoms_initialise: bad array size for array key='//properties%keys(i)//' shape='//properties%entries(i)%len2, error)
             end if
          else
             RAISE_ERROR('atoms_initialise: bad property type='//type//' in properties argument', error)
          end if
       end do

       this%properties = properties  ! deep copy of data from properties argument
    else
       ! By default, we just add properties for Z, species name and position
       call initialise(this%properties)
       
       call add_property(this, 'Z', 0, error=error)
       PASS_ERROR(error)
       call add_property(this, 'pos', 0.0_dp, n_cols=3, error=error)
       PASS_ERROR(error)
       call add_property(this, 'species', repeat(' ', TABLE_STRING_LENGTH), &
            error=error)
       PASS_ERROR(error)
    end if

    call atoms_repoint(this)

    stack_size = 3*N*8/1024
    if (stack_size > 4000 .and. .not. printed_stack_warning) then
       call print('Atoms_Initialise: Stack size must be at least '//stack_size//&
            ' kb, or we get seg faults when functions return arrays.')
       printed_stack_warning = .true.
       call print('Atoms_Initialise: trying to increase stack limit')
       stack_size = int(stack_size/1024)*1024
       stack_size_err = increase_stack(stack_size)
       if (stack_size_err /= 0) then
	call print("Atoms_Initialise: error calling increase_stack err = "// stack_size_err)
       endif
    end if

    call atoms_set_lattice(this, lattice, .false.)

    this%ref_count = 1

    call print("atoms_initialise: Initialised, " // &
         "ref_count = " // this%ref_count, PRINT_ANAL)

  end subroutine atoms_initialise


  subroutine atoms_initialise_ptr(this,N,lattice,&
       properties,params,fixed_size,Nbuffer,error)
    type(Atoms),      pointer                :: this
    integer,                    intent(in)   :: N
    real(dp), dimension(3,3),   intent(in)   :: lattice
    type(Dictionary), optional, intent(in)   :: properties, params
    logical,          optional, intent(in)   :: fixed_size
    integer,          optional, intent(in)   :: Nbuffer
    integer,          optional, intent(out)  :: error

    ! ---

    INIT_ERROR(error)

    if (associated(this)) then
       this%own_this = .false.
    else
       allocate(this)
       this%own_this = .true.
    endif

    call initialise(this, N, lattice, properties, params, fixed_size, &
         Nbuffer, error)
    PASS_ERROR(error)

  endsubroutine atoms_initialise_ptr


  !% Is this atoms object initialised?
  logical function atoms_is_initialised(this)
    type(Atoms), intent(in)  :: this
    atoms_is_initialised = this%ref_count > 0
  end function atoms_is_initialised

  !% Shallow copy of this Atoms object
  subroutine atoms_shallowcopy(this, from)
    type(Atoms), pointer, intent(out)  :: this
    type(Atoms), target                :: from

    ! ---

    this => from
    this%ref_count = this%ref_count + 1

    call print("atoms_shallowcopy: Created shallow copy, " // &
         "ref_count = " // this%ref_count, PRINT_ANAL)

  end subroutine atoms_shallowcopy

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


  subroutine atoms_finalise(this)
    type(Atoms), intent(inout) :: this

    this%ref_count = this%ref_count - 1
    ! Note: this has to be == 0, rather than <= 0. Otherwise it will fail if
    ! finalise is called more than once
    if (this%ref_count == 0) then

       call finalise(this%properties)
       call finalise(this%params)

       ! Nullify pointers
       nullify(this%Z, this%travel, this%mass)
       nullify(this%move_mask, this%thermostat_region, this%damp_mask)
       nullify(this%pos, this%velo, this%acc, this%avgpos, this%oldpos, this%avg_ke)

       call connection_finalise(this%connect)
       call connection_finalise(this%hysteretic_connect)

       this%N = 0
       this%Ndomain = 0
       this%Nbuffer = 0

    else
       if (this%ref_count > 0) then
          call print("atoms_finalise: Not yet finalising, " // &
               "ref_count = " // this%ref_count, PRINT_ANAL)
       endif
    endif

    call print("atoms_finalise: ref_count = " // this%ref_count, PRINT_ANAL)

  end subroutine atoms_finalise


  subroutine atoms_finalise_ptr(this)
    type(Atoms), pointer  :: this

    ! ---

    call finalise(this)
    if (this%own_this .and. this%ref_count == 0) then
       deallocate(this)
       this => NULL()
    endif

  endsubroutine atoms_finalise_ptr

  !Quick multiple finalisations
  subroutine atoms_finalise_multi(at1,at2,at3,at4,at5,at6,at7,at8,at9,at10)
    type(Atoms),           intent(inout) :: at1,at2
    type(Atoms), optional, intent(inout) :: at3,at4,at5,at6,at7,at8,at9,at10
    call atoms_finalise(at1)
    call atoms_finalise(at2)
    if (present(at3)) call atoms_finalise(at3)
    if (present(at4)) call atoms_finalise(at4)
    if (present(at5)) call atoms_finalise(at5)
    if (present(at6)) call atoms_finalise(at6)
    if (present(at7)) call atoms_finalise(at7)
    if (present(at8)) call atoms_finalise(at8)
    if (present(at9)) call atoms_finalise(at9)
    if (present(at10)) call atoms_finalise(at10)
  end subroutine atoms_finalise_multi

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Overloaded Assignment
  !
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_assignment(to,from)

    type(Atoms), intent(inout) :: to
    type(Atoms), intent(in)    :: from

    ! We do not fail if from is unitialised, since overloaded operator routines
    ! are outside scope of error handling mechanism.
    if(.not. is_initialised(from)) then
       to%ref_count = 0
       return
    end if

    call atoms_initialise(to, from%N, from%lattice, from%properties, from%params, from%fixed_size, from%Nbuffer)

    to%use_uniform_cutoff = from%use_uniform_cutoff
    to%cutoff      = from%cutoff
    to%cutoff_break      = from%cutoff_break
    to%nneightol   = from%nneightol

    to%connect     = from%connect
    to%hysteretic_connect     = from%hysteretic_connect

    to%domain_decomposed = from%domain_decomposed

  end subroutine atoms_assignment


   pure function ess_max_len(ess)
      type(extendable_str), intent(in) :: ess(:)
      integer :: ess_max_len

      integer :: i

      ess_max_len = 0
      do i=1, size(ess)
	 if (ess(i)%len > ess_max_len) ess_max_len = ess(i)%len
      end do
   end function ess_max_len

  !% Make a copy of the atoms object 'from' without including
  !% connectivity information. Useful for saving the state of a
  !% dynamical simulation without incurring too great a memory
  !% cost. 
  subroutine atoms_copy_without_connect(to, from, properties, properties_array, error)

    type(Atoms), intent(inout) :: to
    type(Atoms), intent(in)    :: from
    character(len=*), optional, intent(in) :: properties
    character(len=*), optional, intent(in) :: properties_array(:)
    integer, intent(out), optional :: error     
    character(len=ess_max_len(from%properties%keys)) :: tmp_properties_array(from%properties%n)
    integer n_properties    

    INIT_ERROR(error)
    if(.not. is_initialised(from)) then
       RAISE_ERROR("atoms_copy_without_connect: 'from' object is not initialised", error)
    end if

    to%N = from%N
    to%Nbuffer = from%Nbuffer
    to%Ndomain = from%Ndomain
    call set_lattice(to, from%lattice, .false.)

    if (present(properties) .or. present(properties_array)) then
       if (present(properties_array)) then
          call subset(from%properties, properties_array, to%properties, error=error)
          PASS_ERROR(error)
       else
          call parse_string(properties, ':', tmp_properties_array, n_properties, error=error)
          PASS_ERROR(error)
          call subset(from%properties, tmp_properties_array(1:n_properties), to%properties, error=error)
          PASS_ERROR(error)          
       end if
    else
       call deepcopy(to%properties, from%properties, error=error)
       PASS_ERROR(error)
    endif
    to%params = from%params
    to%fixed_size = from%fixed_size
    to%domain_decomposed = from%domain_decomposed
    to%use_uniform_cutoff = from%use_uniform_cutoff
    to%cutoff = from%cutoff
    to%cutoff_break = from%cutoff_break
    to%nneightol = from%nneightol

    call atoms_repoint(to)
    to%ref_count = 1

  end subroutine atoms_copy_without_connect

  subroutine atoms_select(to, from, mask, list, orig_index, error)
    type(Atoms), intent(inout) :: to
    type(Atoms), intent(in) :: from
    logical, intent(in), optional :: mask(:)
    integer, intent(in), optional, target :: list(:)
    logical, optional, intent(in) :: orig_index
    integer, intent(out), optional :: error

    integer :: i, n_list
    integer, allocatable, dimension(:), target :: my_list
    integer, pointer, dimension(:) :: orig_index_ptr, use_list
    logical :: do_orig_index

    INIT_ERROR(error)

    do_orig_index = optional_default(.true., orig_index)

    if ((.not. present(list) .and. .not. present(mask)) .or. (present(list) .and. present(mask))) then
       RAISE_ERROR('atoms_select: either list or mask must be present (but not both)', error)
    end if

    if (present(mask)) then
       if (size(mask) /= from%N) then
          RAISE_ERROR("atoms_select: mismatched sizes of from " // from%N // " and mask " // size(mask), error)
       end if
       call atoms_initialise(to, count(mask), from%lattice)
    else
       call atoms_initialise(to, size(list), from%lattice)
    end if

    call finalise(to%params)
    to%params = from%params
    to%fixed_size = from%fixed_size
    to%use_uniform_cutoff = from%use_uniform_cutoff
    to%cutoff = from%cutoff
    to%cutoff_break = from%cutoff_break
    to%nneightol = from%nneightol

    if(present(mask)) then
       allocate(my_list(count(mask)))
       ! build up a list
       N_list = 1
       do i=1, from%N
          if (mask(i)) then
             my_list(n_list) = i
             n_list = n_list + 1
          endif
       end do
       use_list => my_list
    else
       use_list => list
    end if

    do i=1,from%properties%N
       select case (from%properties%entries(i)%type)

       case(T_INTEGER_A)
          call set_value(to%properties, string(from%properties%keys(i)), from%properties%entries(i)%i_a(use_list))

       case(T_REAL_A)
          call set_value(to%properties, string(from%properties%keys(i)), from%properties%entries(i)%r_a(use_list))

       case(T_LOGICAL_A)
          call set_value(to%properties, string(from%properties%keys(i)), from%properties%entries(i)%l_a(use_list))

       case(T_INTEGER_A2)
          call set_value(to%properties, string(from%properties%keys(i)), from%properties%entries(i)%i_a2(:,use_list))

       case(T_REAL_A2)
          call set_value(to%properties, string(from%properties%keys(i)), from%properties%entries(i)%r_a2(:,use_list))

       case(T_CHAR_A)
          call set_value(to%properties, string(from%properties%keys(i)), from%properties%entries(i)%s_a(:,use_list))

       case default
          RAISE_ERROR('atoms_select: bad property type '//from%properties%entries(i)%type//' key='//from%properties%keys(i), error)

       end select
    end do
    call atoms_repoint(to)

    if (do_orig_index) then
       call add_property(to, 'orig_index', 0, ptr=orig_index_ptr, error=error)
       PASS_ERROR(error)
       orig_index_ptr(:) = use_list(:)
    end if
    if (allocated(my_list)) deallocate(my_list)

    call atoms_repoint(to)
  end subroutine atoms_select

  subroutine atoms_remove_property(this, name, error)
    type(Atoms), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(out), optional :: error

    INIT_ERROR(error)

    if (.not. has_property(this, name)) then
       RAISE_ERROR('atoms_remove_property: no such property "'//trim(name)//'" exists', error)
    end if
    call remove_value(this%properties, name)

  end subroutine atoms_remove_property

  function atoms_has_property(this, name)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    logical :: atoms_has_property

    atoms_has_property = has_key(this%properties, name)
    
  end function atoms_has_property

   subroutine atoms_set_param_value_int(this, key, value)
      type(Atoms), intent(inout) :: this
      character(len=*), intent(in) :: key
      integer, intent(in) :: value
      
      call set_value(this%params, key, value)
   end subroutine atoms_set_param_value_int
   subroutine atoms_set_param_value_int_a(this, key, value)
      type(Atoms), intent(inout) :: this
      character(len=*), intent(in) :: key
      integer, intent(in) :: value(:)
      
      call set_value(this%params, key, value)
   end subroutine atoms_set_param_value_int_a
   subroutine atoms_set_param_value_real(this, key, value)
      type(Atoms), intent(inout) :: this
      character(len=*), intent(in) :: key
      real(dp), intent(in) :: value
      
      call set_value(this%params, key, value)
   end subroutine atoms_set_param_value_real
   subroutine atoms_set_param_value_real_a(this, key, value)
      type(Atoms), intent(inout) :: this
      character(len=*), intent(in) :: key
      real(dp), intent(in) :: value(:)
      
      call set_value(this%params, key, value)
   end subroutine atoms_set_param_value_real_a
   subroutine atoms_set_param_value_real_a2(this, key, value)
      type(Atoms), intent(inout) :: this
      character(len=*), intent(in) :: key
      real(dp), intent(in) :: value(:,:)
      
      call set_value(this%params, key, value)
   end subroutine atoms_set_param_value_real_a2
   subroutine atoms_set_param_value_logical(this, key, value)
      type(Atoms), intent(inout) :: this
      character(len=*), intent(in) :: key
      logical, intent(in) :: value
      
      call set_value(this%params, key, value)
   end subroutine atoms_set_param_value_logical
   subroutine atoms_set_param_value_str(this, key, value)
      type(Atoms), intent(inout) :: this
      character(len=*), intent(in) :: key
      character(1), intent(in) :: value(:)
      
      call set_value(this%params, key, value)
   end subroutine atoms_set_param_value_str

   subroutine atoms_get_param_value_int(this, key, value, error)
      type(Atoms), intent(in) :: this
      character(len=*), intent(in) :: key
      integer, intent(out) :: value
      integer, intent(out), optional :: error
      
      INIT_ERROR(error)
      if (.not. get_value(this%params, key, value)) then
	 RAISE_ERROR("atoms_get_param_value failed to get int value for key='"//trim(key)//"' from this%params", error)
      endif
   end subroutine atoms_get_param_value_int
   subroutine atoms_get_param_value_int_a(this, key, value, error)
      type(Atoms), intent(in) :: this
      character(len=*), intent(in) :: key
      integer, intent(out) :: value(:)
      integer, intent(out), optional :: error
      
      INIT_ERROR(error)
      if (.not. get_value(this%params, key, value)) then
	 RAISE_ERROR("atoms_get_param_value failed to get int array value for key='"//trim(key)//"' from this%params", error)
      endif
   end subroutine atoms_get_param_value_int_a
   subroutine atoms_get_param_value_real(this, key, value, error)
      type(Atoms), intent(in) :: this
      character(len=*), intent(in) :: key
      real(dp), intent(out) :: value
      integer, intent(out), optional :: error
      
      INIT_ERROR(error)
      if (.not. get_value(this%params, key, value)) then
	 RAISE_ERROR("atoms_get_param_value failed to get real value for key='"//trim(key)//"' from this%params", error)
      endif
   end subroutine atoms_get_param_value_real
   subroutine atoms_get_param_value_real_a(this, key, value, error)
      type(Atoms), intent(in) :: this
      character(len=*), intent(in) :: key
      real(dp), intent(out) :: value(:)
      integer, intent(out), optional :: error
      
      INIT_ERROR(error)
      if (.not. get_value(this%params, key, value)) then
	 RAISE_ERROR("atoms_get_param_value failed to get real array value for key='"//trim(key)//"' from this%params", error)
      endif
   end subroutine atoms_get_param_value_real_a
   subroutine atoms_get_param_value_real_a2(this, key, value, error)
      type(Atoms), intent(in) :: this
      character(len=*), intent(in) :: key
      real(dp), intent(out) :: value(:,:)
      integer, intent(out), optional :: error
      
      INIT_ERROR(error)
      if (.not. get_value(this%params, key, value)) then
	 RAISE_ERROR("atoms_get_param_value failed to get real array value for key='"//trim(key)//"' from this%params", error)
      endif
   end subroutine atoms_get_param_value_real_a2
   subroutine atoms_get_param_value_str(this, key, value, error)
      type(Atoms), intent(in) :: this
      character(len=*), intent(in) :: key
      character(1), intent(out) :: value(:)
      integer, intent(out), optional :: error
      
      INIT_ERROR(error)
      if (.not. get_value(this%params, key, value)) then
	 RAISE_ERROR("atoms_get_param_value failed to get str value for key='"//trim(key)//"' from this%params", error)
      endif
   end subroutine atoms_get_param_value_str
   subroutine atoms_get_param_value_es(this, key, value, error)
      type(Atoms), intent(in) :: this
      character(len=*), intent(in) :: key
      type(Extendable_Str), intent(inout) :: value
      integer, intent(out), optional :: error
      
      INIT_ERROR(error)
      if (.not. get_value(this%params, key, value)) then
	 RAISE_ERROR("atoms_get_param_value failed to get es value for key='"//trim(key)//"' from this%params", error)
      endif
   end subroutine atoms_get_param_value_es
   subroutine atoms_get_param_value_logical(this, key, value, error)
      type(Atoms), intent(in) :: this
      character(len=*), intent(in) :: key
      logical, intent(out) :: value
      integer, intent(out), optional :: error
      
      INIT_ERROR(error)
      if (.not. get_value(this%params, key, value)) then
	 RAISE_ERROR("atoms_get_param_value failed to get logical value for key='"//trim(key)//"' from this%params", error)
      endif
   end subroutine atoms_get_param_value_logical

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

  end subroutine atoms_add_property_logical_a


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
  !% Set the cutoff (uniform or factor) to at least the requested value
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_set_cutoff_minimum(this, cutoff, cutoff_break)
    type(Atoms),      intent(inout) :: this
    real(dp),         intent(in)    :: cutoff
    real(dp), optional, intent(in)    :: cutoff_break

    if (present(cutoff_break)) then
      this%cutoff_break = max(this%cutoff_break, cutoff_break)
      this%cutoff = max(this%cutoff, cutoff)
    else
      this%cutoff_break = max(this%cutoff_break, cutoff)
      this%cutoff = max(this%cutoff, cutoff)
    endif

  end subroutine atoms_set_cutoff_minimum

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Specify a uniform neighbour cutoff throughout the system.
  !% If zero, revert to default (uniform_cutoff=false, factor=DEFAULT_NNEIGHTOL)
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_set_cutoff(this, cutoff, cutoff_break)
    type(Atoms),      intent(inout) :: this
    real(dp),         intent(in)    :: cutoff
    real(dp), optional, intent(in)    :: cutoff_break

    if (present(cutoff_break)) then
      if (cutoff .feq. 0.0_dp) then
	this%cutoff = DEFAULT_NNEIGHTOL
      else
	this%cutoff = cutoff
      endif
      if (cutoff_break .feq. 0.0_dp) then
	this%cutoff_break = DEFAULT_NNEIGHTOL
      else
	this%cutoff_break = cutoff_break
      endif
    else
      if (cutoff .feq. 0.0_dp) then
	this%cutoff = DEFAULT_NNEIGHTOL
	this%cutoff_break = DEFAULT_NNEIGHTOL
      else
	this%cutoff = cutoff
	this%cutoff_break = cutoff
      endif
    endif

    if (cutoff .feq. 0.0_dp) then
      this%use_uniform_cutoff = .false.
    else
      this%use_uniform_cutoff = .true.
    endif

  end subroutine atoms_set_cutoff

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Specify the neighbour cutoff to be a mulitple of the bond length
  !% of the two atoms' types.
  !% If zero, revert to default (DEFAULT_NNEIGHTOL)
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_set_cutoff_factor(this, factor, factor_break)
    type(Atoms),      intent(inout) :: this
    real(dp),         intent(in)    :: factor
    real(dp), optional, intent(in)    :: factor_break

    if (present(factor_break)) then
      if (factor .feq. 0.0_dp) then
	this%cutoff = DEFAULT_NNEIGHTOL
      else
	this%cutoff = factor
      endif
      if (factor_break .feq. 0.0_dp) then
	this%cutoff_break = DEFAULT_NNEIGHTOL
      else
	this%cutoff_break = factor_break
      endif
    else
      if (factor .feq. 0.0_dp) then
	this%cutoff = DEFAULT_NNEIGHTOL
	this%cutoff_break = DEFAULT_NNEIGHTOL
      else
	this%cutoff = factor
	this%cutoff_break = factor
      endif
    endif
    this%use_uniform_cutoff = .false.

  end subroutine atoms_set_cutoff_factor

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Return the actual cutoff in \AA{} used by this atoms object
  !% used to form 'Z1---Z2' bonds. If 'this%use_uniform_cutoff' is
  !% true, then this is simply 'this%cutoff', otherwise the
  !% cutoff is used multiplied by the 'Z1---Z2' bond-length.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function atoms_cutoff(this, Z1, Z2)
    type(Atoms), intent(in) :: this
    integer, intent(in) :: Z1, Z2
    real(dp) :: atoms_cutoff

    if (this%use_uniform_cutoff) then
       atoms_cutoff = this%cutoff
    else
       atoms_cutoff = this%cutoff*bond_length(Z1, Z2)
    end if

  end function atoms_cutoff

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Return the actual cutoff in \AA{} used by this atoms object
  !% used to break 'Z1---Z2' bonds. If 'this%use_uniform_cutoff' is
  !% true, then this is simply 'this%cutoff', otherwise the
  !% cutoff is used multiplied by the 'Z1---Z2' bond-length.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function atoms_cutoff_break(this, Z1, Z2)
    type(Atoms), intent(in) :: this
    integer, intent(in) :: Z1, Z2
    real(dp) :: atoms_cutoff_break

    if (this%use_uniform_cutoff) then
       atoms_cutoff_break = this%cutoff_break
    else
       atoms_cutoff_break = this%cutoff_break*bond_length(Z1, Z2)
    end if

  end function atoms_cutoff_break


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Zero data in an Atoms structure ---
  !% this doesn\'t finalise it or change it\'s size. We zero 'this%pos',
  !% 'this%Z' and 'this%species'.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_zero(this,indices)

    type(Atoms), intent(inout)    :: this
    integer, optional, intent(in) :: indices !% Optionally only zero the specified indices.

    if(present(indices)) then
       this%pos(:,indices) = 0.0_dp
       this%Z(indices) = 0
       this%species(:,indices) = ' '
    else
       this%pos = 0.0_dp
       this%Z = 0
       this%species = ' '
    end if

  end subroutine atoms_zero


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Change the lattice vectors, keeping the inverse lattice vectors
  !% up to date. Optionally map the existing atoms into the new cell
  !% and recalculate connectivity.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_set_lattice(this,new_lattice,scale_positions,remap,reconnect)

    type(Atoms),              intent(inout) :: this
    real(dp), dimension(3,3), intent(in)    :: new_lattice
    logical, intent(in)                     :: scale_positions
    logical, optional,        intent(in)    :: remap, reconnect
    real(dp), dimension(:,:), allocatable :: frac

    ! only do the allocation if needed
    if(scale_positions) then
       allocate(frac(3,this%N))
       frac = this%g .mult. this%pos
    end if

    this%lattice = new_lattice

    if(scale_positions) then
       this%pos = this%lattice .mult. frac
       deallocate(frac)
    end if

    call matrix3x3_inverse(this%lattice,this%g)

    if (present(reconnect)) then
       if (reconnect) then
          call calc_connect(this)
          return
       end if           
    end if

    if (present(remap)) then
       if (remap) call set_map_shift(this)
    end if

    this%is_orthorhombic = ( (count(this%lattice(:,1) == 0.0_dp) <= 1) .and. &
			     (count(this%lattice(:,2) == 0.0_dp) <= 1) .and. &
			     (count(this%lattice(:,3) == 0.0_dp) <= 1) )
    this%is_periodic(1) = any(this%lattice(:,1) /= 0.0_dp)
    this%is_periodic(2) = any(this%lattice(:,2) /= 0.0_dp)
    this%is_periodic(3) = any(this%lattice(:,3) /= 0.0_dp)

  end subroutine atoms_set_lattice

  subroutine atoms_set_atoms_singlez(this, Z)
    type(Atoms), intent(inout) :: this
    integer, intent(in)        :: Z
    integer, allocatable, dimension(:) :: Zarray

    allocate(Zarray(this%N))
    Zarray = Z
    call atoms_set_atoms(this, Zarray)
    deallocate(Zarray)
  end subroutine atoms_set_atoms_singlez

  !% Set atomic numbers and optionally masses (if mass property is present)
  !% If 'mass' is not specified then 'ElementMass(Z)' is used.
  subroutine atoms_set_atoms(this, Z, mass)
    type(Atoms), intent(inout) :: this
    integer, dimension(:), intent(in) :: Z
    real(dp), optional, dimension(:), intent(in) :: mass
    
    integer i

    this%Z = Z
    if (has_property(this, 'mass')) then
       this%mass = ElementMass(Z)
       ! Optionally override with user specified masses
       if (present(mass)) this%mass = mass
    end if
    if (has_property(this, 'species')) then
       do i=1,this%N
          this%species(:,i) = ' '
          this%species(1:3,i) = s2a(ElementName(Z(i)))
       end do
    end if
    
  end subroutine atoms_set_atoms

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! 
  ! Simple query functions
  !
  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  !% Return the number of neighbour that atom $i$ has.
  !% If the optional arguments max_dist or max_factor are present 
  !% then only neighbours closer than this cutoff are included.
  function atoms_n_neighbours(this, i, max_dist, max_factor, alt_connect, error) result(n)
    type(Atoms), intent(in), target :: this
    integer, intent(in) :: i
    real(dp), optional, intent(in) :: max_dist, max_factor
    type(Connection), optional, intent(in), target :: alt_connect
    integer, intent(out), optional :: error     
    
    integer :: n

    integer :: j, m
    real(dp) :: r_ij
    type(Connection), pointer :: use_connect

    INIT_ERROR(error)
    if (present(alt_connect)) then
      use_connect => alt_connect
    else
      use_connect => this%connect
    endif

    if (.not. use_connect%initialised) then
       RAISE_ERROR('Atoms_N_Neighbours: Atoms structure has no connectivity data. Call calc_connect first.', error)
    end if

    if (.not. associated(use_connect%neighbour1(i)%t)) then
      n = 0
      return
    endif

    if (.not. present(max_dist) .and. .not. present(max_factor)) then
       ! All neighbours
       n = use_connect%neighbour1(i)%t%N + use_connect%neighbour2(i)%t%N
    else if (present(max_dist)) then
       ! Only count neighbours within max_dist distance of i
       n = 0
       do m=1,use_connect%neighbour1(i)%t%N+use_connect%neighbour2(i)%t%N
          j = atoms_neighbour(this, i, m, distance=r_ij)
          if (r_ij < max_dist) n = n + 1
       end do
    else if (present(max_factor)) then
       ! Only count neighbours within max_factor of i
       n = 0
       do m=1,use_connect%neighbour1(i)%t%N+use_connect%neighbour2(i)%t%N
          j = atoms_neighbour(this, i, m, distance=r_ij)
          if (r_ij < bond_length(this%Z(i),this%Z(j))*max_factor) n = n + 1
       end do
    else
       RAISE_ERROR('Atoms_N_Neighbours: optional arguments max_dist and max_factor must not both be present', error)
    end if

  end function atoms_n_neighbours

  function atoms_neighbour_index(this, i, n, index, t, is_j, alt_connect, error) result(j)
    type(Atoms), intent(in), target :: this
    integer :: i, j, n
    integer,  intent(out) :: index
    type(Table), pointer, intent(out) :: t
    logical, intent(out) :: is_j
    type(Connection), optional, intent(in), target :: alt_connect
    integer, intent(out), optional :: error     

    integer :: i_n1n, j_n1n
    type(Connection), pointer :: use_connect

    INIT_ERROR(error)
    if (present(alt_connect)) then
      use_connect => alt_connect
    else
      use_connect => this%connect
    endif

    if (use_connect%initialised) then
       i_n1n = n-use_connect%neighbour2(i)%t%N
       if (n <= use_connect%neighbour2(i)%t%N) then
          j = use_connect%neighbour2(i)%t%int(1,n)
          j_n1n = use_connect%neighbour2(i)%t%int(2,n)
          index = j_n1n
	  t => use_connect%neighbour1(j)%t
	  is_j = .true.
       else if (i_n1n <= use_connect%neighbour1(i)%t%N) then
          j = use_connect%neighbour1(i)%t%int(1,i_n1n)
          index = i_n1n
	  t => use_connect%neighbour1(i)%t
	  is_j = .false.
       else
          RAISE_ERROR('atoms_neighbour: '//n//' out of range for atom '//i//' Should be in range 1 < n <= '//atoms_n_neighbours(this, i), error)
       end if
    else
       RAISE_ERROR('atoms_neighbour_index: Connect structure not initialized. Call calc_connect first.', error)
    end if

  end function atoms_neighbour_index

  function atoms_neighbour_minimal(this, i, n, shift, index, alt_connect) result(j)
    type(Atoms), intent(in), target :: this
    integer ::i, j, n
    integer,  intent(out) :: shift(3)
    integer,  intent(out) :: index
    type(Connection), optional, intent(in), target :: alt_connect

    type(Table), pointer :: t
    logical :: is_j

    j = atoms_neighbour_index(this, i, n, index, t, is_j, alt_connect)

    if (is_j) then
      shift = - t%int(2:4,index)
    else
      shift = t%int(2:4,index)
    endif

  end function atoms_neighbour_minimal

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
  !% if distance > max_dist, return 0, and do not waste time calculating other quantities
  function atoms_neighbour(this, i, n, distance, diff, cosines, shift, index, max_dist, jn, alt_connect, error) result(j)
    type(Atoms), intent(in), target :: this
    integer ::i, j, n
    real(dp), optional, intent(out) :: distance
    real(dp), dimension(3), optional, intent(out) :: diff
    real(dp), optional, intent(out) :: cosines(3)
    integer,  optional, intent(out) :: shift(3)
    integer,  optional, intent(out) :: index
    real(dp), optional, intent(in)  :: max_dist
    integer,  optional, intent(out) :: jn
    type(Connection), optional, intent(in), target :: alt_connect
    integer, intent(out), optional :: error     

    real(dp)::mydiff(3), norm_mydiff
    integer ::myshift(3)
    integer ::i_n1n, j_n1n, i_njn, m
    type(Connection), pointer :: use_connect

    INIT_ERROR(error)
    if (present(alt_connect)) then
      use_connect => alt_connect
    else
      use_connect => this%connect
    endif

    if (.not. associated(use_connect%neighbour1(i)%t)) then
      RAISE_ERROR("called atoms_neighbour on atom " // i // " which has no allocated neighbour1 table", error)
    endif

    ! First we give the neighbour2 entries (i > j) then the neighbour1 (i <= j)
    ! This order chosen to give neighbours in approx numerical order but doesn't matter
    if (use_connect%initialised) then
       i_n1n = n-use_connect%neighbour2(i)%t%N
       if (n <= use_connect%neighbour2(i)%t%N) then
          j = use_connect%neighbour2(i)%t%int(1,n)
          j_n1n = use_connect%neighbour2(i)%t%int(2,n)
          if(present(index)) index = j_n1n
       else if (i_n1n <= use_connect%neighbour1(i)%t%N) then
          j = use_connect%neighbour1(i)%t%int(1,i_n1n)
          if(present(index)) index = i_n1n
       else
          RAISE_ERROR('atoms_neighbour: '//n//' out of range for atom '//i//' Should be in range 1 < n <= '//atoms_n_neighbours(this, i), error)
       end if
    else
       RAISE_ERROR('atoms_neighbour: Atoms structure has no connectivity data. Call calc_connect first.', error)
    end if

    if(present(jn)) then
       if(i < j) then
          do i_njn = 1, use_connect%neighbour2(j)%t%N
             if( (use_connect%neighbour2(j)%t%int(1,i_njn)==i) .and. &
                 (use_connect%neighbour2(j)%t%int(2,i_njn)==i_n1n) ) jn = i_njn
          enddo
       elseif(i > j) then
          jn = j_n1n + use_connect%neighbour2(j)%t%N
       else
          do i_njn = 1, use_connect%neighbour1(j)%t%N
             if( (use_connect%neighbour1(j)%t%int(1,i_njn) == i) .and. &
                  all(use_connect%neighbour1(j)%t%int(2:4,i_njn) == -use_connect%neighbour1(i)%t%int(2:4,i_n1n))) &
	       jn = i_njn + use_connect%neighbour2(j)%t%N
          enddo
       endif
    endif
          
    ! found neighbour, now check for optional requests
    if(present(distance)) then
       if(i <= j) then 
          distance = use_connect%neighbour1(i)%t%real(1,i_n1n)
       else
          distance = use_connect%neighbour1(j)%t%real(1,j_n1n)
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
	   if (use_connect%neighbour1(i)%t%real(1,i_n1n) > max_dist) then
	     j = 0
	     return
	   endif
	 else
	   if (use_connect%neighbour1(j)%t%real(1,j_n1n) > max_dist) then
	     j = 0
	     return
	   endif
	 endif
       endif
    end if

    if(present(diff) .or. present(cosines) .or. present(shift)) then
       if (i <= j) then
          myshift = use_connect%neighbour1(i)%t%int(2:4,i_n1n)
       else
          myshift = -use_connect%neighbour1(j)%t%int(2:4,j_n1n)
       end if

       if(present(shift)) shift = myshift

       if(present(diff) .or. present(cosines)) then
          !mydiff = this%pos(:,j) - this%pos(:,i) + (this%lattice .mult. myshift)
          mydiff = this%pos(:,j) - this%pos(:,i)
          do m=1,3
             ! forall(k=1:3) mydiff(k) = mydiff(k) + this%lattice(k,m) * myshift(m)
             mydiff(1:3) = mydiff(1:3) + this%lattice(1:3,m) * myshift(m)
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

  end function atoms_neighbour

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Adding and Removing atoms: adding/removing single or multiple atoms
  ! For the atoms variable in Dynamical System, this should be called
  ! shift there to avoid inconsistencies with DS's data
  ! 
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine add_atom_single(this, pos, Z, mass, travel, error)

    type(Atoms),       intent(inout)            :: this
    real(dp),          intent(in), dimension(3) :: pos
    integer,           intent(in)               :: Z
    real(dp), optional,  intent(in)             :: mass
    integer, optional, intent(in), dimension(3) :: travel
    integer, optional, intent(out) :: error

    INIT_ERROR(error)

    if(present(travel)) then
       if(present(mass)) then
          call add_atom_multiple(this, pos=reshape(pos, (/3,1/)), Z=(/Z/), mass=(/mass/), travel=reshape(travel, (/3,1/)), error=error)
       else
          call add_atom_multiple(this, pos=reshape(pos, (/3,1/)), Z=(/Z/), travel=reshape(travel, (/3,1/)), error=error)
       end if
    else
       if(present(mass)) then
          call add_atom_multiple(this, pos=reshape(pos, (/3,1/)), Z=(/Z/), mass=(/mass/), error=error)
       else
          call add_atom_multiple(this, pos=reshape(pos, (/3,1/)), Z=(/Z/), error=error)
       end if
    end if
    PASS_ERROR(error)

  end subroutine add_atom_single

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine add_atom_multiple(this, pos, Z, mass,  velo, acc, travel, error)

    type(Atoms),       intent(inout)              :: this
    real(dp), intent(in), dimension(:,:) :: pos
    integer,  intent(in), dimension(:)   :: Z
    real(dp), optional, intent(in), dimension(:)  :: mass
    integer, optional, intent(in), dimension(:,:) :: travel
    real(dp),optional, intent(in), dimension(:,:) :: velo, acc
    integer, intent(out), optional :: error

    integer                                       :: oldN,i
    integer, allocatable, dimension(:) :: tmp_int
    integer, allocatable, dimension(:,:) :: tmp_int2
    real(dp), allocatable, dimension(:) :: tmp_real
    real(dp), allocatable, dimension(:,:) :: tmp_real2
    logical, allocatable, dimension(:) :: tmp_logical
    character, allocatable, dimension(:,:) :: tmp_char

    INIT_ERROR(error)    

    if (this%fixed_size) then
       RAISE_ERROR("add_atom_multiple: Atoms object cannot be resized (this%fixed_size = .true.)", error)
    end if

    oldN = this%N
    this%N = this%N + size(Z)
    this%Ndomain = this%N
    
    !Check the sizes of the input arrays for consistency
    call check_size('Pos',pos,(/3,size(Z)/),'Add_Atom',error)
    PASS_ERROR(error)

    if (present(mass)) then
       call check_size('Mass',mass,size(Z), 'Add_Atom', error)
       PASS_ERROR(error)
    end if

    if (present(travel)) then
       call check_size('Travel',travel,(/3,size(Z)/),'Add_Atom', error)
       PASS_ERROR(error)
    end if

    if (present(velo)) then
       call check_size('Velo', velo, (/3,size(Z)/), 'Add_Atom', error)
       PASS_ERROR(error)
    end if

    if (present(acc)) then
       call check_size('Acc', acc, (/3,size(Z)/), 'Add_Atom', error)
       PASS_ERROR(error)
    end if

    ! Only resize if the actual number of particles is now larger than the
    ! buffer size.
    if (this%N > this%Nbuffer) then

       ! Resize property data arrays, copying old data.
       ! this will break any existing pointers so we call atoms_repoint() immediately after
       ! (note that user-held pointers will stay broken, there is no way to fix this
       ! since Fortran does not allow pointer-to-pointer types)
       do i=1,this%properties%N
          select case (this%properties%entries(i)%type)

          case(T_INTEGER_A)
             allocate(tmp_int(this%n))
             tmp_int(1:oldN) = this%properties%entries(i)%i_a
             tmp_int(oldn+1:this%n) = 0
             call set_value(this%properties, string(this%properties%keys(i)), tmp_int)
             deallocate(tmp_int)

          case(T_REAL_A)
             allocate(tmp_real(this%n))
             tmp_real(1:oldN) = this%properties%entries(i)%r_a
             tmp_real(oldn+1:this%n) = 0.0_dp
             call set_value(this%properties, string(this%properties%keys(i)), tmp_real)
             deallocate(tmp_real)

          case(T_LOGICAL_A)
             allocate(tmp_logical(this%n))
             tmp_logical(1:oldN) = this%properties%entries(i)%l_a
             tmp_logical(oldn+1:this%n) = .false.
             call set_value(this%properties, string(this%properties%keys(i)), tmp_logical)
             deallocate(tmp_logical)

          case(T_INTEGER_A2)
             allocate(tmp_int2(this%properties%entries(i)%len2(1),this%n))
             tmp_int2(:,1:oldN) = this%properties%entries(i)%i_a2
             tmp_int2(:,oldn+1:this%n) = 0
             call set_value(this%properties, string(this%properties%keys(i)), tmp_int2)
             deallocate(tmp_int2)

          case(T_REAL_A2)
             allocate(tmp_real2(this%properties%entries(i)%len2(1),this%n))
             tmp_real2(:,1:oldN) = this%properties%entries(i)%r_a2
             tmp_real2(:,oldn+1:this%n) = 0.0_dp
             call set_value(this%properties, string(this%properties%keys(i)), tmp_real2)
             deallocate(tmp_real2)

          case(T_CHAR_A)
             allocate(tmp_char(this%properties%entries(i)%len2(1),this%n))
             tmp_char(:,1:oldN) = this%properties%entries(i)%s_a
             tmp_char(:,oldn+1:this%n) = ' '
             call set_value(this%properties, string(this%properties%keys(i)), tmp_char)
             deallocate(tmp_char)

          case default
             RAISE_ERROR('atoms_add: bad property type '//this%properties%entries(i)%type//' key='//this%properties%keys(i), error)

          end select
       end do
       call atoms_repoint(this)
    endif

    this%Nbuffer = max(this%N, this%Nbuffer)

    ! First check the integer properties...
    if (.not. has_key(this%properties, 'Z')) then
       RAISE_ERROR('Atoms_Add: this atoms has no Z property', error)
    end if
    this%z(oldN+1:this%N) = Z
    
    ! set species from Z
    if (.not. has_key(this%properties, 'species')) then
       RAISE_ERROR('Atoms_Add: this atoms has no species property', error)
    end if
    do i=1,size(Z)
       this%species(:,oldN+i) = ' '
       this%species(1:3,oldN+i) = s2a(ElementName(Z(i)))
    end do

    if (present(travel)) then
       if (.not. has_key(this%properties, 'travel')) then
          RAISE_ERROR('Atoms_Add: this atoms has no travel property', error)
       end if
       this%travel(:,oldN+1:this%N) = travel
    else
       if (has_key(this%properties, 'travel')) then
          this%travel(:,oldN+1:this%N) = 0
       end if
    end if
    
    ! Set masks to 1 if properties for them exist
    if (has_key(this%properties, 'move_mask')) &
         this%move_mask(oldN+1:this%N) = 1
    if (has_key(this%properties, 'damp_mask')) &
         this%damp_mask(oldN+1:this%N) = 1
    if (has_key(this%properties, 'thermostat_region')) &
         this%thermostat_region(oldN+1:this%N) = 1
    
    ! ... and now the real properties
    if (has_key(this%properties, 'mass')) then
       if (present(mass)) then
          this%mass(oldN+1:this%N) = mass
       else
          this%mass(oldN+1:this%N) = ElementMass(Z)
       end if
    else if (present(mass)) then
       ! mass specified but property doesn't yet exist, so create it...
       call add_property(this, 'mass', ElementMass(this%Z), ptr=this%mass, error=error)
       PASS_ERROR(error)
       ! ... and then override for new atoms
       this%mass(oldN+1:this%N) = mass
    end if

    if (.not. has_key(this%properties, 'pos')) then
       RAISE_ERROR('Atoms_Add: this atoms has no pos property', error)
    end if
    this%pos(:,oldN+1:this%N) = pos
    
    if (present(velo) .and. has_key(this%properties, 'velo')) &
         this%velo(:,oldN+1:this%N) = velo

    if (present(acc) .and. has_key(this%properties, 'acc')) &
         this%acc(:,oldN+1:this%N) = acc

    call finalise(this%connect)

  end subroutine add_atom_multiple


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine remove_atom_single(this, i, error)

    type(Atoms), intent(inout) :: this
    integer,     intent(in)    :: i
    integer, intent(out), optional :: error

    INIT_ERROR(error)
    call remove_atom_multiple(this,(/i/),error)
    PASS_ERROR(error)

  end subroutine remove_atom_single

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  subroutine remove_atom_multiple(this, atom_indices, error)

    type(Atoms), intent(inout)                 :: this
    integer,     intent(in), dimension(:)      :: atom_indices
    integer, intent(out), optional :: error

    integer i, copysrc
    integer, allocatable, dimension(:), target :: new_indices
    integer, pointer, dimension(:) :: include_list
    integer, allocatable, dimension(:) :: tmp_int
    integer, allocatable, dimension(:,:) :: tmp_int2
    real(dp), allocatable, dimension(:) :: tmp_real
    real(dp), allocatable, dimension(:,:) :: tmp_real2
    logical, allocatable, dimension(:) :: tmp_logical
    character, allocatable, dimension(:,:) :: tmp_char
    integer, dimension(size(atom_indices)) :: sorted
    integer, dimension(:), allocatable :: uniqed

    INIT_ERROR(error)

    if (this%fixed_size) then
       RAISE_ERROR("remove_atom_multiple: Atoms object cannot be resized (this%fixed_size = .true.)", error)
    end if

    !Delete the connection data because the atomic indices become mangled
    call connection_finalise(this%connect)

    ! Permute new_indices, following algorithm in table_record_delete_multiple,
    ! so that atom ordering is same as it was under old scheme when properties were
    ! stored as columns in Table this%data.
    ! (we find first atomc to be removed and last atom to not be removed
    !  and swap them. Repeat until all atoms to be removed are at the end.)

    sorted = atom_indices     ! Get our own copy of the  indices so we can sort them
    call sort_array(sorted)
    call uniq(sorted, uniqed) ! remove duplicates from sorted indices

    allocate(new_indices(this%N))
    do i=1,this%N
       new_indices(i) = i
    end do

    copysrc = this%N
    do i=1,size(uniqed)
       do while(is_in_array(uniqed,copysrc))
          copysrc = copysrc - 1
       end do

       if (uniqed(i) > copysrc) exit
       new_indices(uniqed(i)) = new_indices(copysrc)
       copysrc = copysrc - 1
    end do

    ! update N
    this%N = this%N - size(atom_indices)
    this%Ndomain = this%N
    this%Nbuffer = this%N

    include_list => new_indices(1:this%N)
    
    ! Resize property data arrays, copying old data.
    ! this will break any existing pointers so we call atoms_repoint() immediately after
    ! (note that user-held pointers will stay broken, there is no way to fix this
    ! since Fortran does not allow pointer-to-pointer types)
    do i=1,this%properties%N
       select case (this%properties%entries(i)%type)

       case(T_INTEGER_A)
          allocate(tmp_int(this%n))
          tmp_int(:) = this%properties%entries(i)%i_a(include_list)
          call set_value(this%properties, string(this%properties%keys(i)), tmp_int)
          deallocate(tmp_int)

       case(T_REAL_A)
          allocate(tmp_real(this%n))
          tmp_real(:) = this%properties%entries(i)%r_a(include_list)
          call set_value(this%properties, string(this%properties%keys(i)), tmp_real)
          deallocate(tmp_real)

       case(T_LOGICAL_A)
          allocate(tmp_logical(this%n))
          tmp_logical(:) = this%properties%entries(i)%l_a(include_list)
          call set_value(this%properties, string(this%properties%keys(i)), tmp_logical)
          deallocate(tmp_logical)

       case(T_INTEGER_A2)
          allocate(tmp_int2(this%properties%entries(i)%len2(1),this%n))
          tmp_int2(:,:) = this%properties%entries(i)%i_a2(:,include_list)
          call set_value(this%properties, string(this%properties%keys(i)), tmp_int2)
          deallocate(tmp_int2)

       case(T_REAL_A2)
          allocate(tmp_real2(this%properties%entries(i)%len2(1),this%n))
          tmp_real2(:,:) = this%properties%entries(i)%r_a2(:,include_list)
          call set_value(this%properties, string(this%properties%keys(i)), tmp_real2)
          deallocate(tmp_real2)

       case(T_CHAR_A)
          allocate(tmp_char(this%properties%entries(i)%len2(1),this%n))
          tmp_char(:,:) = this%properties%entries(i)%s_a(:,include_list)
          call set_value(this%properties, string(this%properties%keys(i)), tmp_char)
          deallocate(tmp_char)

       case default
          deallocate(include_list)
          RAISE_ERROR('remove_atom_multiple: bad property type '//this%properties%entries(i)%type//' key='//this%properties%keys(i), error)

       end select
    end do
    call atoms_repoint(this)

    deallocate(uniqed, new_indices)

  end subroutine remove_atom_multiple



  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! 
  !% Map atomic fractional positions back into the unit cell
  !% $-0.5 \le t_x,t_y,t_z < 0.5$
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_map_into_cell(this)
    type(Atoms), intent(inout) :: this

    integer      :: i,j, n, m
    integer      :: shift(3)
    logical      :: mapped
    integer, pointer :: map_shift(:,:)


    ! Loop over all atoms
    ! Convert cartesians to lattice co-ords
    ! If outside unit cell:
    !    add lattice vectors to get back into cell
    !    update travel
    !    update shifts

    if (.not. has_property(this, 'travel')) then
       call add_property(this, 'travel', 0, n_cols=3, ptr2=this%travel)
    end if

    do i=1,this%N
      call map_into_cell(this%pos(:,i), this%lattice, this%g, shift, mapped)
      if (mapped) then
	 this%travel(:,i) = this%travel(:,i) - shift
	 if (this%connect%initialised) then
	   do n=1,atoms_n_neighbours(this, i)  ! Loop over all atom i's neighbours

	      j = atoms_neighbour(this, i, n, index=m) ! get neighbour

	      ! now update the data for atom i and the current neighbour
	      if (i < j) then
		 this%connect%neighbour1(i)%t%int(2:4,m) = this%connect%neighbour1(i)%t%int(2:4,m) + shift
	      else if (i > j) then
		 this%connect%neighbour1(j)%t%int(2:4,m) = this%connect%neighbour1(j)%t%int(2:4,m) - shift
              else
                 ! do nothing when i == j
	      end if
	   end do
	end if ! this%connect%initialised
      end if ! mapped
    end do ! i=1..N

    if (assign_pointer(this, 'map_shift', map_shift)) map_shift = 0

  end subroutine atoms_map_into_cell

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


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% The subroutine 'calc_dists' updates the stored distance tables using 
  !% the stored connectivity and shifts. This should be called every time
  !% any atoms are moved (e.g. it is called by 'advance_verlet').
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine calc_dists(this, parallel, error)

    type(Atoms), intent(inout) :: this
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
#endif
    INIT_ERROR(error)

    ! Flag to specify whether or not to parallelise calculation.
    ! Only actually run in parallel if parallel==.true. AND
    ! _MPI is #defined. Default to serial mode.
    do_parallel = .false.
    if (present(parallel)) do_parallel = parallel

#ifdef _MPI
    if (do_parallel) then
       ! Nelements = sum(this%connect%neighbour1(i)%t%N)
       Nelements = 0
       do i=1,this%N
          Nelements = Nelements + this%connect%neighbour1(i)%t%N
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

    if (.not.this%connect%initialised) then
         RAISE_ERROR('CalcDists: Connect is not yet initialised', error)
      endif

    do i = 1, this%N

#ifdef _MPI
       if (do_parallel) then
          mpi_old_pos = mpi_pos
          mpi_pos = mpi_pos + this%connect%neighbour1(i)%t%N

          ! cycle loop if processor rank does not match
          if(mod(i, mpi_n_procs()) .ne. mpi_id()) cycle
       end if
#endif

       do n = 1, atoms_n_neighbours(this, i) 

          j = atoms_neighbour_minimal(this, i, n, shift=shift, index=index)

          ! j_pos = this%pos(:,j) + ( this%lattice .mult. shift )
          j_pos(:) = this%pos(:,j) + ( this%lattice(:,1) * shift(1) + this%lattice(:,2) * shift(2) + this%lattice(:,3) * shift(3) )

          if (i <= j) then
             this%connect%neighbour1(i)%t%real(1,index) = norm(j_pos - this%pos(:,i))
          else
             this%connect%neighbour1(j)%t%real(1,index) = norm(j_pos - this%pos(:,i))
          end if

       end do

#ifdef _MPI
       if (do_parallel) then
	  if (mpi_old_pos <= mpi_pos-1) then
	    mpi_send(mpi_old_pos:mpi_pos-1) = &
		 this%connect%neighbour1(i)%t%real(1,1:this%connect%neighbour1(i)%t%N)
	  end if
       end if
#endif      

    end do

#ifdef _MPI
    if (do_parallel) then
       ! collect mpi results
       if (Nelements > 0) then
	 call mpi_allreduce(mpi_send, mpi_recv, &
	      size(mpi_send), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, PRINT_ALWAYS)
	 call abort_on_mpi_error(PRINT_ALWAYS, "Calc_Dists: MPI_ALL_REDUCE()")
       end if

       mpi_pos = 1
       do i=1, this%N
	  if (this%connect%neighbour1(i)%t%N > 0) then
	    this%connect%neighbour1(i)%t%real(1,1:this%connect%neighbour1(i)%t%N) = &
		 mpi_recv(mpi_pos:mpi_pos+this%connect%neighbour1(i)%t%N-1)
	    mpi_pos = mpi_pos + this%connect%neighbour1(i)%t%N
	  endif
       end do

       if (Nelements > 0) then
	 deallocate(mpi_send, mpi_recv)
       end if
    end if
#endif

  end subroutine calc_dists


  !% Difference vector between atoms $i$ and $j$ if they are separated by a shift of 'shift'
  function diff(this, i, j, shift)
    type(Atoms), intent(in)    :: this
    integer,     intent(in)    :: i,j
    integer,  dimension(3)     :: shift
    real(dp), dimension(3) :: diff

    diff = this%pos(:,j) - this%pos(:,i) + (this%lattice .mult. shift)

  end function diff


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !  diff_min_image interface
  !  
  !  return relative vector from one position to another, adhering to PBC
  !  and minimum image conventions
  ! 
  !  Flavours are: atom-atom, vector-atom, atom-vector, vector-vector
  !
  !  All are accessible using the 'diff_min_image' interface
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  function diff_atom_atom(this, i, j, shift)

    type(Atoms), intent(in)    :: this
    integer,     intent(in)    :: i,j
    integer,  dimension(3), optional     :: shift
    real(dp), dimension(3)     :: diff_atom_atom
    real(dp)                   :: dummy

    integer, dimension(3) :: myshift

    dummy = distance_min_image(this,i,j, shift=myshift)

    diff_atom_atom = this%pos(:,j) - this%pos(:,i) + (this%lattice .mult. myshift)

    if (present(shift)) shift = myshift

  end function diff_atom_atom


  function diff_vec_atom(this, v, j)

    type(Atoms), intent(in)    :: this
    real(dp), dimension(3)     :: v
    integer,     intent(in)    :: j
    real(dp), dimension(3)     :: diff_vec_atom
    integer,  dimension(3)     :: shift
    real(dp)                   :: dummy

    dummy = distance_min_image(this,v,j, shift=shift)

    diff_vec_atom = this%pos(:,j) - v + (this%lattice .mult. shift)

  end function diff_vec_atom

  function diff_atom_vec(this, i, w)

    type(Atoms), intent(in)    :: this
    integer,     intent(in)    :: i
    real(dp), dimension(3)     :: w
    real(dp), dimension(3)     :: diff_atom_vec
    integer,  dimension(3)     :: shift
    real(dp)                   :: dummy

    dummy = distance_min_image(this,i,w, shift=shift)

    diff_atom_vec = w - this%pos(:,i) + (this%lattice .mult. shift)

  end function diff_atom_vec

  function diff_vec_vec(this, v, w)

    type(Atoms), intent(in)    :: this
    real(dp), dimension(3)     :: v, w
    real(dp), dimension(3)     :: diff_vec_vec
    integer,  dimension(3)     :: shift
    real(dp)                   :: dummy

    dummy = distance_min_image(this,v,w, shift=shift)

    diff_vec_vec = w - v + (this%lattice .mult. shift)

  end function diff_vec_vec

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Return the real position of the atom, taking into account the
  !% stored travel across the periodic boundary conditions.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function realpos(this,i)

    type(Atoms), intent(in) :: this
    integer,     intent(in) :: i
    real(dp), dimension(3)  :: realpos

    if (associated(this%travel)) then
       realpos = (this%lattice .mult. this%travel(:,i)) + this%pos(:,i)
    else
       realpos = this%pos(:,i)
    endif

  end function realpos


  !% Return distance between atoms $i$ and $j$ if they are separated by a shift
  !% of 'shift'.
  function distance(this, i, j, shift)
    type(Atoms), intent(in)::this
    integer,     intent(in)::i, j, shift(3)
    real(dp)::distance

    distance = norm(this%pos(:,j)+(this%lattice .mult. shift)-this%pos(:,i))
  end function distance

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! distance_min_image interface
  !
  ! Actual distance computing routines. 
  !
  ! The real work is done in the function that
  ! computes the distance of two general vector positions.  when
  ! atomic indices are specified, they are first converted to vector
  ! positions.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function distance8_atom_atom(this,i,j,shift)

    type(Atoms),                        intent(in)   :: this
    integer,                            intent(in)   :: i,j
    integer,  optional, dimension(3),   intent(out)  :: shift
    real(dp)                                         :: distance8_atom_atom

    distance8_atom_atom = distance8_vec_vec(this,this%pos(:,i),this%pos(:,j),shift)

  end function distance8_atom_atom

  function distance8_atom_vec(this,i,v,shift)

    type(Atoms),                        intent(in)  :: this
    integer,                            intent(in)  :: i
    real(dp),           dimension(3),   intent(in)  :: v
    integer,  optional, dimension(3),   intent(out) :: shift
    real(dp)                                        :: distance8_atom_vec

    distance8_atom_vec = distance8_vec_vec(this,this%pos(:,i),v,shift)

  end function distance8_atom_vec

  function distance8_vec_atom(this,v,j,shift)

    type(Atoms),                        intent(in)  :: this
    real(dp),           dimension(3),   intent(in)  :: v
    integer,                            intent(in)  :: j
    integer,  optional, dimension(3),   intent(out) :: shift
    real(dp)                                        :: distance8_vec_atom

    distance8_vec_atom = distance8_vec_vec(this,v,this%pos(:,j),shift)

  end function distance8_vec_atom

  ! This is the general function

  function distance8_vec_vec(this,v,w,shift)

    type(Atoms),                        intent(in)  :: this
    real(dp),           dimension(3),   intent(in)  :: v,w
    integer,  optional, dimension(3),   intent(out) :: shift
    real(dp)                                        :: distance8_vec_vec, dist2, tmp
    real(dp),           dimension(3)                :: dvw, lattice_coord
    integer,            dimension(3)                :: init_val
    integer                                         :: i,j,k, i_shift(3)

    !get the difference vector and convert to lattice co-ordinates
    !use the precomputed matrix inverse if possible
    dvw = w - v
    call vec_map_into_cell(dvw, this%lattice, this%g, i_shift)
    lattice_coord = this%g .mult. dvw

    init_val = (/0,0,0/)

    !work out which block of 8 cells we are testing
    where (lattice_coord > 0.0_dp) init_val = -1

    dist2 = huge(1.0_dp) ! effectively +ve infinity

    !now loop over the cells and test
    do k=init_val(3), init_val(3)+1
       do j=init_val(2), init_val(2)+1
          do i=init_val(1), init_val(1)+1

             !construct the shifted vector
             tmp = norm2(dvw + this%lattice(:,1)*i +this%lattice(:,2)*j + this%lattice(:,3)*k)
             !test if it is the smallest so far and store the shift if necessary
             if (tmp < dist2) then
                dist2 = tmp
                if (present(shift)) shift = (/i,j,k/) + i_shift
             end if

          end do
       end do
    end do

    distance8_vec_vec = sqrt(dist2)

  end function distance8_vec_vec

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Cosine of the angle j--i--k
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function cosine(this,i,j,k,error)

    type(Atoms), intent(in) :: this
    integer,     intent(in)    :: i,j,k
    real(dp)                   :: cosine
    real(dp), dimension(3)     :: ij, ik
    integer, intent(out), optional :: error     
    
    INIT_ERROR(error)
    if ((i == j) .or. (i == k)) then
       RAISE_ERROR('Cosine: i == j or i == k', error)
    end if

    ij = diff_min_image(this,i,j)
    ik = diff_min_image(this,i,k)
    cosine = (ij .dot. ik) / (norm(ij)*norm(ik))

  end function cosine

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Cosine of the angle n--i--m where {$n,m$} are the {$n$th, $m$th} neighbours of i 
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function cosine_neighbour(this,i,n,m)

    type(Atoms), intent(in) :: this
    integer,     intent(in)    :: i,n,m
    real(dp)                   :: cosine_neighbour
    real(dp), dimension(3)     :: in, im
    integer::j

    if(n == m) then
       cosine_neighbour = 1.0_dp
       return
    end if
    j = atoms_neighbour(this, i, n, diff=in)
    j = atoms_neighbour(this, i, m, diff=im)
    cosine_neighbour = (in .dot. im) / (norm(in)*norm(im))

  end function cosine_neighbour

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Direction cosines
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% Given two atoms $i$ and $j$ and a shift returns the direction 
  !% cosines of the differnece vector from $i$ to $j$.
  function direction_cosines(this,i,j,shift)

    type(Atoms), intent(in) :: this
    integer,     intent(in)    :: i,j, shift(3)
    real(dp), dimension(3)     :: direction_cosines, diffv

    diffv = diff(this,i,j,shift)
    direction_cosines = diffv / norm(diffv)

  end function direction_cosines

  !% Direction cosines of the difference vector from $i$ to $j$
  function direction_cosines_min_image(this,i,j,error)

    type(Atoms), intent(in) :: this
    integer,     intent(in)    :: i,j
    real(dp), dimension(3)     :: direction_cosines_min_image, diffv
    integer, intent(out), optional :: error     
    
    INIT_ERROR(error)
    if (i == j) then
       RAISE_ERROR('Cosines: i == j', error)
    end if

    diffv = diff_min_image(this,i,j)
    direction_cosines_min_image = diffv / norm(diffv)

  end function direction_cosines_min_image

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Connectivity procedures: Initialise and Finalise
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  !% Initialise a Connection object for a given number of atoms 'N'.
  !% If the optional Atoms argument is present then we calculate
  !% the atomic density to initialise the default lengths of the neighbour
  !% list for efficient memory usage.
   subroutine connection_initialise(this, N, Nbuffer, pos, lattice, g,  origin, extent, nn_guess, fill)
    type(Connection),   intent(inout) :: this
    integer,            intent(in)    :: N    ! No. of atoms
    integer,            intent(in)    :: Nbuffer    ! Buffer size
    real(dp), optional, intent(in) :: pos(:,:), lattice(3,3), g(3,3)
    real(dp), optional, intent(in) :: origin(3), extent(3,3)
    integer, optional, intent(in) :: nn_guess
    logical, optional, intent(in) :: fill

    logical :: do_fill

    do_fill = optional_default(.true., fill)

    ! If already initialised, destroy the existing data and start again
    if (this%initialised) call connection_finalise(this)

    if (do_fill) call connection_fill(this, N, Nbuffer, pos, lattice, g, origin, extent, nn_guess)
  end subroutine connection_initialise

  subroutine connection_fill(this, N, Nbuffer, pos, lattice, g, origin, extent, nn_guess, error)
    type(Connection),   intent(inout) :: this
    integer,            intent(in)    :: N    ! No. of atoms
    integer,            intent(in)    :: Nbuffer    ! Buffer size
    real(dp), optional, intent(in) :: pos(:,:), lattice(3,3), g(3,3)
    real(dp), optional, intent(in) :: origin(3), extent(3,3)
    integer, optional, intent(in) :: nn_guess
    integer, intent(out), optional :: error     
    
    integer                           :: i, do_nn_guess
    real(dp)                          :: extent_inv(3,3), subregion_center(3)
    logical :: do_subregion

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
          call allocate(this%neighbour1(i)%t,4,1, 0, 0, max(do_nn_guess, 1))
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

    allocate( this%cell(cellsNa,cellsNb,cellsNc) )

    do k=1,cellsNc
       do j=1,cellsNb
          do i=1,cellsNa
             call allocate(this%cell(i,j,k),1,0,0,0,max(1,(av_atoms+2*stdev_atoms))) !Good for 97.5% of cells
             call set_increment(this%cell(i,j,k),max(1,stdev_atoms))
          end do
       end do
    end do

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

    do k=1,this%cellsNc
       do j=1,this%cellsNb
          do i=1,this%cellsNa
             call finalise(this%cell(i,j,k))
          end do
       end do
    end do

    deallocate(this%cell)

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
       do k = 1, from%cellsNc
          do j = 1, from%cellsNb
             do i = 1, from%cellsNa
                call append(to%cell(i,j,k),from%cell(i,j,k))
             end do
          end do
       end do
    end if

  end subroutine connection_assignment

  !% OMIT
  subroutine wipe_cells(this)

    type(Connection) :: this
    integer          :: i,j,k

    if (this%cells_initialised) then

       do k = 1, this%cellsNc
          do j = 1, this%cellsNb
             do i = 1, this%cellsNa
                call wipe(this%cell(i,j,k))
             end do
          end do
       end do

    end if

  end subroutine wipe_cells

  !% Test if an atom j is one of i's nearest neighbours
  function is_nearest_neighbour_abs_index(this,i,j, alt_connect)
    type(Atoms), intent(in), target :: this
    integer,     intent(in) :: i,j
    type(Connection), intent(in), optional, target :: alt_connect
    logical                 :: is_nearest_neighbour_abs_index

    integer :: ji
    real(dp) :: d

    is_nearest_neighbour_abs_index = .false.
    do ji=1, atoms_n_neighbours(this, i)
      if (atoms_neighbour(this, i, ji, distance=d, alt_connect=alt_connect) == j) then
	if (d < bond_length(this%Z(i),this%Z(j))*this%nneightol) then
	  is_nearest_neighbour_abs_index = .true.
	  return
	endif
      endif
    end do
  end function is_nearest_neighbour_abs_index


  !% Test if an atom's $n$th neighbour is one if its nearest neighbours
  function is_nearest_neighbour(this,i,n, alt_connect)

    type(Atoms), intent(in), target :: this
    integer,     intent(in) :: i,n
    type(Connection), intent(in), optional, target :: alt_connect
    logical                 :: is_nearest_neighbour

    real(dp)                :: d
    integer :: j

    is_nearest_neighbour = .false.

    j = atoms_neighbour(this, i, n, distance=d, alt_connect=alt_connect)
    if (d < (bond_length(this%Z(i),this%Z(j))*this%nneightol)) &
         is_nearest_neighbour = .true.

  end function is_nearest_neighbour

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
	return
      endif
    endif

#ifdef DEBUG
    if(current_verbosity() >= PRINT_ANAL)  call print('d = '//d, PRINT_ANAL)
#endif

    if (d < use_cutoff) then
       call add_bond(this, pos, lattice, i, j, shift, d, error)
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

  subroutine set_bonds(this, pairs, shifts, error)
    type(Atoms), intent(inout) :: this
    integer, intent(in) :: pairs(:,:)
    integer, intent(in) :: shifts(:,:)
    integer, intent(out), optional :: error

    integer i

    INIT_ERROR(error)
    if (.not.this%connect%initialised) then
       call connection_initialise(this%connect, this%N, this%Nbuffer)
    else
       ! otherwise just wipe the connection table
       call wipe(this%connect)
    end if

    if (size(pairs,1) /= 2) then
       RAISE_ERROR("set_bond pairs not a 2xN array", error)
    end if
    if (size(shifts,1) /= 3) then
       RAISE_ERROR("set_bond shifts not a 3xN array", error)
    end if
    if (size(pairs,2) /= size(shifts,2)) then
       RAISE_ERROR("set_bonds called with mismatching pairs and shifts sizes", error)
    end if

    do i=1, size(pairs,2)
      call add_bond(this%connect, this%pos, this%lattice, pairs(1,i), pairs(2,i), shifts(:,i), error=error)
      PASS_ERROR(error)
    end do
  end subroutine set_bonds

  subroutine add_bond(this, pos, lattice, i, j, shift, d, error)
    type(Connection), intent(inout) :: this
    real(dp), intent(in) :: pos(:,:), lattice(3,3)
    integer,     intent(in)    :: i,j
    integer,     intent(in)    :: shift(3)
    real(dp), intent(in), optional :: d
    integer, intent(out), optional :: error

    real(dp) :: dd
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

    if (present(d)) then
      dd = d
    else
      dd = norm(pos(:,j)+(lattice .mult. shift) - pos(:,i))
    endif

    ! Add full details to neighbour1 for smaller of i and j
    call append(this%neighbour1(ii)%t, (/jj, sign(1,jj-ii)*shift /), (/ dd /))
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


  !%  As for 'calc_connect', but perform the connectivity update
  !%  hystertically: atoms must come within 'cutoff' to be considered
  !%  neighbours, and then will remain connect until them move apart
  !%  further than 'cutoff_break'.
  !%
  !%  Typically 'alt_connect' should be set to the
  !%  'hysteretic_connect' attribute. 'origin' and 'extent'
  !%  vectors can be used to restrict the hysteretic region to only
  !%  part of the entire system -- the 'estimate_origin_extent()'
  !%  routine in clusters.f95 can be used to guess suitable values.
  subroutine calc_connect_hysteretic(this, alt_connect, origin, extent, own_neighbour, store_is_min_image, error)
    type(Atoms), intent(inout), target           :: this
    type(Connection), intent(inout), target, optional :: alt_connect
    real(dp), optional :: origin(3), extent(3,3)
    logical, optional, intent(in) :: own_neighbour, store_is_min_image
    integer, intent(out), optional :: error

    integer                              :: cellsNa,cellsNb,cellsNc,i,j,k,i2,j2,k2,i3,j3,k3,i4,j4,k4,n1,n2,atom1,atom2
    integer                              :: cell_image_Na, cell_image_Nb, cell_image_Nc
    integer                              :: min_cell_image_Na, max_cell_image_Na, min_cell_image_Nb, max_cell_image_Nb, &
                                            min_cell_image_Nc, max_cell_image_Nc
    real(dp)                             :: cutoff
    integer :: ji, s_ij(3), nn_guess
    logical my_own_neighbour, my_store_is_min_image
    type(Connection), pointer :: use_connect
    logical :: change_i, change_j, change_k, broken
    integer, pointer :: map_shift(:,:)

    INIT_ERROR(error)
    if (present(alt_connect)) then
      use_connect => alt_connect
    else
      use_connect => this%connect
    endif

    my_own_neighbour = optional_default(.false., own_neighbour)
    my_store_is_min_image = optional_default(.true., store_is_min_image)

    if (this%cutoff < 0.0_dp .or. this%cutoff_break < 0.0_dp) then
       RAISE_ERROR('calc_connect: Negative cutoff radius ' // this%cutoff // ' ' // this%cutoff_break, error)
    end if

    if (this%cutoff > this%cutoff_break) then
       RAISE_ERROR('calc_connect: Negative hysteresis cutoff radius formation ' // this%cutoff // ' > breaking ' // this%cutoff_break, error)
    end if

    if ((this%cutoff .feq. 0.0_dp) .or. (this%cutoff_break .feq. 0.0_dp)) then
      call wipe(use_connect)
      return
    endif

    !Calculate the cutoff value we should use in dividing up the simulation cell
    if (this%use_uniform_cutoff) then
       !The cutoff has been specified by the user, so use that value
       cutoff = this%cutoff
    else
       !Otherwise, find the maximum covalent radius of all the atoms in the structure, double it
       !and multiple by cutoff. This makes sure we can deal with the worst case scenario of big atom 
       !bonded to big atom
       cutoff = 0.0_dp
       do i = 1, this%N
          if (ElementCovRad(this%Z(i)) > cutoff) cutoff = ElementCovRad(this%Z(i))
       end do
       cutoff = (2.0_dp * cutoff) * this%cutoff
    end if

    call print("calc_connect: cutoff calc_connect " // cutoff, PRINT_NERD)

    if (present(origin) .and. present(extent)) then
      cellsNa = 1
      cellsNb = 1
      cellsNc = 1
    else
      call divide_cell(this%lattice, cutoff, cellsNa, cellsNb, cellsNc)
    endif

    call print("calc_connect: cells_N[abc] " // cellsNa // " " // cellsNb // " " // cellsNc, PRINT_NERD)

    ! If the lattice has changed, then the cells need de/reallocating
    if ((cellsNa /= use_connect%cellsNa) .or. &
         (cellsNb /= use_connect%cellsNb) .or. &
         (cellsNc /= use_connect%cellsNc)) call connection_cells_finalise(use_connect)

    ! figure out how many unit cell images we will need to loop over in each direction
    call fit_box_in_cell(cutoff, cutoff, cutoff, this%lattice, cell_image_Na, cell_image_Nb, cell_image_Nc)
    ! cell_image_N{a,b,c} apply to a box of side 2*cutoff. Since we loop from -cell_image_N to
    ! +cell_image_N we can reduce them as follows

    cell_image_Na = max(1,(cell_image_Na+1)/2)
    cell_image_Nb = max(1,(cell_image_Nb+1)/2)
    cell_image_Nc = max(1,(cell_image_Nc+1)/2)

    ! Estimate number of neighbours of each atom. Factor of 1/2 assumes
    ! half will go in neighbour1, half in neighbour2.
    nn_guess = int(0.5_dp*4.0_dp/3.0_dp*PI*cutoff**3*this%N/cell_volume(this%lattice)*cell_image_na*cell_image_nb*cell_image_nc)

    call print('calc_connect: image cells '//cell_image_Na//'x'//cell_image_Nb//'x'//cell_image_Nc, PRINT_NERD)

    ! Allocate space for the connection object if needed
    if (present(origin) .and. present(extent)) then
      if (.not.use_connect%initialised) then
	 call connection_initialise(use_connect, this%N, this%Nbuffer, this%pos, this%lattice, this%g, origin, extent, nn_guess)
      else
	 call connection_fill(use_connect, this%N, this%Nbuffer, this%pos, this%lattice, this%g, origin, extent, nn_guess)
      end if
    else
      if (.not.use_connect%initialised) then
	 call connection_initialise(use_connect, this%N, this%Nbuffer, this%pos, this%lattice, this%g, nn_guess=nn_guess)
      else
	 call connection_fill(use_connect, this%N, this%Nbuffer, this%pos, this%lattice, this%g, nn_guess=nn_guess)
      end if
    endif

    if (.not.use_connect%cells_initialised) then
      call connection_cells_initialise(use_connect, cellsNa, cellsNb, cellsNc,this%N)
    endif

    ! Partition the atoms into cells
    call partition_atoms(use_connect, this, error=error)
    PASS_ERROR(error)
    if (.not. assign_pointer(this, 'map_shift', map_shift)) then
       RAISE_ERROR("calc_connect impossibly failed to assign map_shift pointer", error)
    end if

    ! look for bonds that have been broken, and remove them
    do i=1, this%N
      ji = 1
      do
	if (ji > atoms_n_neighbours(this, i, alt_connect=use_connect)) exit
	j = atoms_neighbour(this, i, ji, shift = s_ij, alt_connect=use_connect)
        broken = test_break_bond(use_connect, this%cutoff_break, this%use_uniform_cutoff, &
             this%Z, this%pos, this%lattice, i, j, s_ij, error)
        PASS_ERROR(error)
	if (.not. broken) then
	  ji = ji + 1 ! we didn't break this bond, so go to next one
	              ! if we did break a bond, ji now points to a different bond, so don't increment it
	endif
      end do
!      do ji=1, atoms_n_neighbours(this, i, alt_connect=use_connect)
!	j = atoms_neighbour(this, i, ji, shift = s_ij, alt_connect=use_connect)
!	call test_break_bond(use_connect, this%cutoff_break, this%use_uniform_cutoff, &
!	  this%Z, this%pos, this%lattice, i, j, s_ij)
!      end do
    end do

    ! Here is the main loop:
    ! Go through each cell and update the connectivity between atoms in this cell and neighbouring cells
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
	     call get_min_max_images(this%is_periodic, cellsNa, cellsNb, cellsNc, &
	        cell_image_Na, cell_image_Nb, cell_image_Nc, i, j, k, change_i, change_j, change_k, &
	        min_cell_image_Na, max_cell_image_Na, min_cell_image_Nb, max_cell_image_Nb, min_cell_image_Nc, max_cell_image_Nc)
	     change_i = .false.
	     change_j = .false.
	     change_k = .false.

             !Loop over atoms in cell(i,j,k)
             do n1 = 1, use_connect%cell(i,j,k)%N

                atom1 = use_connect%cell(i,j,k)%int(1,n1)

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

                         do n2 = 1, use_connect%cell(i3,j3,k3)%N

                            atom2 = use_connect%cell(i3,j3,k3)%int(1,n2)
                            ! omit atom2 < atom1
                            if (atom1 > atom2) cycle
                            ! omit self in the same cell without shift
                            if (.not. my_own_neighbour .and. (atom1 == atom2 .and. & 
                                 (i4==0 .and. j4==0 .and. k4==0) .and. &
                                 (i==i3 .and. j==j3 .and. k==k3))) cycle

                            call test_form_bond(use_connect, this%cutoff, this%use_uniform_cutoff, &
                                 this%Z, this%pos, this%lattice, atom1,atom2, &
				 (/i4-map_shift(1,atom1)+map_shift(1,atom2),j4-map_shift(2,atom1)+map_shift(2,atom2),k4-map_shift(3,atom1)+map_shift(3,atom2)/), &
				 .true., error)
                            PASS_ERROR(error)

                         end do ! n2

                      end do ! i2
		   end do ! j2
                end do ! k2

             end do ! n1

          end do ! i
       end do ! j
    end do ! k

    if (my_store_is_min_image) then
       if (allocated(use_connect%is_min_image)) deallocate(use_connect%is_min_image)
       allocate(use_connect%is_min_image(this%n))
       do i=1,this%n
          if (associated(use_connect%neighbour1(i)%t)) then
             use_connect%is_min_image(i) = is_min_image(this, i, alt_connect=use_connect)
          else
             use_connect%is_min_image(i) = .false.
          end if
       end do
    end if

  end subroutine calc_connect_hysteretic

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

  !% Fast $O(N)$ connectivity calculation routine. It divides the unit cell into similarly shaped subcells,
  !% of sufficient size that sphere of radius 'cutoff' is contained in a subcell, at least in the directions 
  !% in which the unit cell is big enough. For very small unit cells, there is only one subcell, so the routine
  !% is equivalent to the standard $O(N^2)$ method.
  subroutine calc_connect(this, alt_connect, own_neighbour, store_is_min_image, skip_zero_zero_bonds, error)
    type(Atoms), intent(inout), target    :: this
    type(Connection), intent(inout), target, optional :: alt_connect
    logical, optional, intent(in) :: own_neighbour, store_is_min_image, skip_zero_zero_bonds
    integer, intent(out), optional :: error

    integer                              :: cellsNa,cellsNb,cellsNc,i,j,k,i2,j2,k2,i3,j3,k3,i4,j4,k4,n1,n2,atom1,atom2
    integer                              :: cell_image_Na, cell_image_Nb, cell_image_Nc, nn_guess, n_occ
    integer                              :: min_cell_image_Na, max_cell_image_Na, min_cell_image_Nb, max_cell_image_Nb, &
                                            min_cell_image_Nc, max_cell_image_Nc
    real(dp)                             :: cutoff, density, volume_per_cell
    logical my_own_neighbour, my_store_is_min_image, my_skip_zero_zero_bonds, do_fill
    type(Connection), pointer :: use_connect
    logical :: change_i, change_j, change_k
    integer, pointer :: map_shift(:,:)

    INIT_ERROR(error)
    if (present(alt_connect)) then
      use_connect => alt_connect
    else
      use_connect => this%connect
    endif

    my_own_neighbour = optional_default(.false., own_neighbour)
    my_store_is_min_image = optional_default(.true., store_is_min_image)
    my_skip_zero_zero_bonds = optional_default(.false., skip_zero_zero_bonds)

    if (this%cutoff < 0.0_dp .or. this%cutoff_break < 0.0_dp) then
       RAISE_ERROR('calc_connect: Negative cutoff radius ' // this%cutoff // ' ' // this%cutoff_break, error)
    end if

    if ((this%cutoff .feq. 0.0_dp) .or. (this%cutoff_break .feq. 0.0_dp)) then
      call wipe(use_connect)
      return
    endif
    !Calculate the cutoff value we should use in dividing up the simulation cell
    if (this%use_uniform_cutoff) then
       !The cutoff has been specified by the user, so use that value
       cutoff = this%cutoff
    else
       !Otherwise, find the maximum covalent radius of all the atoms in the structure, double it
       !and multiple by cutoff. This makes sure we can deal with the worst case scenario of big atom 
       !bonded to big atom
       cutoff = 0.0_dp
       do i = 1, this%N
          if (ElementCovRad(this%Z(i)) > cutoff) cutoff = ElementCovRad(this%Z(i))
       end do
       cutoff = (2.0_dp * cutoff) * this%cutoff
    end if

    call print("calc_connect: cutoff calc_connect " // cutoff, PRINT_VERBOSE)

    call divide_cell(this%lattice, cutoff, cellsNa, cellsNb, cellsNc)

    call print("calc_connect: cells_N[abc] " // cellsNa // " " // cellsNb // " " // cellsNc, PRINT_VERBOSE)

    ! If the lattice has changed, then the cells need de/reallocating
    if ((cellsNa /= use_connect%cellsNa) .or. &
         (cellsNb /= use_connect%cellsNb) .or. &
         (cellsNc /= use_connect%cellsNc)) call connection_finalise(use_connect)

    ! figure out how many unit cell images we will need to loop over in each direction
    call fit_box_in_cell(cutoff, cutoff, cutoff, this%lattice, cell_image_Na, cell_image_Nb, cell_image_Nc)
    ! cell_image_N{a,b,c} apply to a box of side 2*cutoff. Since we loop from -cell_image_N to
    ! +cell_image_N we can reduce them as follows

    cell_image_Na = max(1,(cell_image_Na+1)/2)
    cell_image_Nb = max(1,(cell_image_Nb+1)/2)
    cell_image_Nc = max(1,(cell_image_Nc+1)/2)

    call print('calc_connect: image cells '//cell_image_Na//'x'//cell_image_Nb//'x'//cell_image_Nc, PRINT_VERBOSE)

    ! Allocate space for the connection object if needed
    if (.not.use_connect%initialised) then
       call connection_initialise(use_connect, this%N, this%Nbuffer, this%pos, this%lattice, this%g, fill=.false.)
       do_fill = .true.
    else
       ! otherwise just wipe the connection table
       call wipe(use_connect)
       do_fill = .false.
    end if

    if (.not.use_connect%cells_initialised) &
         call connection_cells_initialise(use_connect, cellsNa, cellsNb, cellsNc,this%N)

    ! Partition the atoms into cells
    call partition_atoms(use_connect, this, error=error)
    PASS_ERROR(error)
    if (.not. assign_pointer(this, 'map_shift', map_shift)) then
       RAISE_ERROR("calc_connect impossibly failed to assign map_shift pointer", error)
    end if

    if (do_fill) then
       volume_per_cell = cell_volume(this%lattice)/real(cellsNa*cellsNb*cellsNc,dp)
 
       ! Count the occupied cells so vacuum does not contribute to average number density
       n_occ = 0
       do k=1,cellsNc
          do j=1,cellsNb
             do i=1,cellsNa
                if (use_connect%cell(i,j,k)%n /= 0) n_occ = n_occ + 1
             end do
          end do
       end do
       density = this%n/(n_occ*volume_per_cell)

       ! Sphere of radius "cutoff", assume roughly half neighbours in neighbour1 and half in neighbour2
       nn_guess = int(4.0_dp/3.0_dp*PI*cutoff**3*density)/2

       call print('calc_connect: occupied cells '//n_occ//'/'//(cellsNa*cellsNb*cellsNc)//' = '//(n_occ/real(cellsNa*cellsNb*cellsNc,dp)), PRINT_VERBOSE)
       call print('calc_connect: estimated number of neighbours per atom = '//nn_guess, PRINT_VERBOSE)

       call connection_fill(use_connect, this%N, this%Nbuffer, this%pos, this%lattice, this%g, nn_guess=nn_guess)
    end if

    ! Here is the main loop:
    ! Go through each cell and update the connectivity between atoms in this cell and neighbouring cells
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
	     call get_min_max_images(this%is_periodic, cellsNa, cellsNb, cellsNc, &
	        cell_image_Na, cell_image_Nb, cell_image_Nc, i, j, k, change_i, change_j, change_k, &
	        min_cell_image_Na, max_cell_image_Na, min_cell_image_Nb, max_cell_image_Nb, min_cell_image_Nc, max_cell_image_Nc)
	     change_i = .false.
	     change_j = .false.
	     change_k = .false.
	     do n1=1, use_connect%cell(i,j,k)%N

		atom1 = use_connect%cell(i,j,k)%int(1,n1)

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

			 do n2 = 1, use_connect%cell(i3,j3,k3)%N

			    atom2 = use_connect%cell(i3,j3,k3)%int(1,n2)
			    ! omit atom2 < atom1
			    if (atom1 > atom2) cycle
                            
                            if (my_skip_zero_zero_bonds .and. this%z(atom1) == 0 .and. this%z(atom2) == 0) cycle 
       
			    ! omit self in the same cell without shift
			    if (.not. my_own_neighbour .and. (atom1 == atom2 .and. & 
				 (i4==0 .and. j4==0 .and. k4==0) .and. &
				 (i==i3 .and. j==j3 .and. k==k3))) cycle
			    call test_form_bond(use_connect,this%cutoff, this%use_uniform_cutoff, &
			      this%Z, this%pos, this%lattice, atom1,atom2, &
				 (/i4-map_shift(1,atom1)+map_shift(1,atom2),j4-map_shift(2,atom1)+map_shift(2,atom2),k4-map_shift(3,atom1)+map_shift(3,atom2)/), &
				 .false., error)
                            PASS_ERROR(error)

			 end do ! n2

		      end do ! i2
		   end do ! j2
		end do ! k2

	     end do ! n1
	  end do ! i
       end do ! j
    end do ! k

    if (my_store_is_min_image) then
       if (allocated(use_connect%is_min_image)) deallocate(use_connect%is_min_image)
       allocate(use_connect%is_min_image(this%n))
       do i=1,this%n
          use_connect%is_min_image(i) = is_min_image(this, i, error=error)
          PASS_ERROR(error)
       end do
    end if

  end subroutine calc_connect

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
      
      if (current_verbosity() >= PRINT_NERD) then
         call print('get_min_max_images cell_image_na min='//min_cell_image_na//' max='//max_cell_image_na, PRINT_NERD)
         call print('get_min_max_images cell_image_nb min='//min_cell_image_nb//' max='//max_cell_image_nb, PRINT_NERD)
         call print('get_min_max_images cell_image_nc min='//min_cell_image_nc//' max='//max_cell_image_nc, PRINT_NERD)
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
    if (.not.my_dont_wipe) call wipe_cells(this)

!!    ! Make sure all atomic positions are within the cell
!!    call map_into_cell(at)
    if (.not. assign_pointer(at, 'map_shift', map_shift)) then
       call add_property(at, 'map_shift', 0, 3)
       if (.not. assign_pointer(at, 'map_shift', map_shift)) then
	  RAISE_ERROR("partition_atoms impossibly failed to assign map_shift pointer", error)
       end if
    endif

    neighbour1_allocated = allocated(this%neighbour1)

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
       call append( this%cell(i,j,k), (/n/) )
    end do

  end subroutine partition_atoms

   subroutine set_map_shift(this, error)
      type(Atoms), intent(inout) :: this
      integer, intent(out), optional :: error

      integer, pointer :: map_shift(:,:)
      integer n
      real(dp) :: lat_pos(3)

      INIT_ERROR(error)
      if (.not. assign_pointer(this, 'map_shift', map_shift)) then
	 call add_property(this, 'map_shift', 0, 3)
	 if (.not. assign_pointer(this, 'map_shift', map_shift)) then
            RAISE_ERROR("partition_atoms impossibly failed to assign map_shift pointer", error)
         end if
      endif

      call print("this%N = " // this%N // ", size(map_shift, 2) = " // size(map_shift, 2), PRINT_ALWAYS)
      do n = 1, this%N
         lat_pos = this%g .mult. this%pos(:,n)
         map_shift(:,n) = - floor(lat_pos+0.5_dp)
      end do
    end subroutine set_map_shift


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
      integer, intent(out), optional :: error
      
      INIT_ERROR(error)
      if (.not. this%cells_initialised) then
         RAISE_ERROR('cell_n: cells are not initialised', error)
      end if

      cell_n = this%cell(i,j,k)%n

    end function cell_n

    function cell_contents(this, i, j, k, n, error)
      type(Connection), intent(in) :: this
      integer, intent(in) :: i, j, k, n
      integer :: cell_contents
      integer, intent(out), optional :: error

      INIT_ERROR(error)
      if (.not. this%cells_initialised) then
         RAISE_ERROR('cell_n: cells are not initialised', error)
      end if

      cell_contents = this%cell(i,j,k)%int(1,n)

    end function cell_contents

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



   !
   !% Returns the (unsigned) volume of the simulation cell of 'this'
   !
   function atoms_cell_volume(this)
     type(Atoms), intent(in) :: this
     real(dp)                :: atoms_Cell_Volume

     atoms_cell_volume = cell_volume(this%lattice)
   end function atoms_cell_volume

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

   !
   ! Make_Lattice
   !
   !% Make a matrix of lattice vectors from the lengths 'a','b','c'
   !% and the angles 'alpha', 'beta' and 'gamma'.
   !% One length must be supplied. Any missing angle is assumed to be 90 degrees
   !% and any missing length is assumed to be 'a'.
   !% The vectors are created in a right-handed order.
   !
   function make_lattice(a,b,c,alpha,beta,gamma,error) result(lattice)

     real(dp),           intent(in) :: a
     real(dp), optional, intent(in) :: b,c
     real(dp), optional, intent(in) :: alpha,beta,gamma
     integer, intent(out), optional :: error
     real(dp), dimension(3,3)       :: lattice
     real(dp)                       :: my_b, my_c,            &
                                       cos_alpha, cos2_alpha, &
                                       cos_beta,  cos2_beta,  &
                                       cos_gamma, cos2_gamma, &
                                       sin_gamma, sin2_gamma, &
                                       my_alpha, my_beta,my_gamma
     INIT_ERROR(error)
     my_b = a
     my_c = a

     if (present(b)) my_b = b
     if (present(c)) my_c = c

     my_alpha = PI/2.0_dp
     my_beta  = PI/2.0_dp
     my_gamma = PI/2.0_dp

     if (present(alpha)) my_alpha = alpha
     if (present(beta))  my_beta  = beta
     if (present(gamma)) my_gamma = gamma

     if ( (my_alpha <= 0.0_dp) .or. (my_alpha <= 0.0_dp) .or. (my_gamma <= 0.0_dp) ) then
        RAISE_ERROR('Make_Lattice: Negative angles are not permitted', error)
     end if

     if ( (my_alpha + my_beta) < my_gamma ) then
        RAISE_ERROR('Make_Lattice: alpha + beta < gamma', error)
     end if
     if ( (my_beta + my_gamma) < my_alpha ) then
        RAISE_ERROR('Make_Lattice: beta + gamma < alpha', error)
     end if
     if ( (my_gamma + my_alpha) < my_beta ) then
        RAISE_ERROR('Make_Lattice: gamma + alpha < beta', error)
     end if

     cos_alpha = cos(my_alpha); cos2_alpha = cos_alpha*cos_alpha
     cos_beta  = cos(my_beta);  cos2_beta  = cos_beta *cos_beta
     cos_gamma = cos(my_gamma); cos2_gamma = cos_gamma*cos_gamma
     sin_gamma = sin(my_gamma); sin2_gamma = sin_gamma*sin_gamma

     ! a
     lattice(1,1) = a
     lattice(2,1) = 0.0_dp
     lattice(3,1) = 0.0_dp
     ! b
     lattice(1,2) = my_b * cos_gamma
     lattice(2,2) = my_b * sin_gamma
     lattice(3,2) = 0.0_dp
     ! c
     lattice(1,3) = my_c * cos_beta
     lattice(2,3) = my_c * (cos_alpha - cos_beta*cos_gamma) / sin_gamma
     lattice(3,3) = my_c * sqrt(1.0_dp - (cos2_alpha + cos2_beta - 2.0_dp*cos_alpha*cos_beta*cos_gamma)/ sin2_gamma)

   end function make_lattice

   !
   !% Opposite of Make_Lattice.
   !% Given a lattice, return a,b,c,alpha,beta and gamma (if needed)
   !
   subroutine get_lattice_params(lattice,a,b,c,alpha,beta,gamma)

     real(dp), dimension(3,3), intent(in)  :: lattice
     real(dp), optional,       intent(out) :: a,b,c,alpha,beta,gamma
     !local variables
     real(dp), dimension(3)                :: a_1,a_2,a_3

     a_1 = lattice(:,1)
     a_2 = lattice(:,2)
     a_3 = lattice(:,3)

     if (present(a)) a = norm(a_1)
     if (present(b)) b = norm(a_2)
     if (present(c)) c = norm(a_3)

     if (present(alpha)) alpha = acos( (a_2 .dot. a_3) / (norm(a_2)*norm(a_3)) )
     if (present(beta))  beta  = acos( (a_3 .dot. a_1) / (norm(a_3)*norm(a_1)) )
     if (present(gamma)) gamma = acos( (a_1 .dot. a_2) / (norm(a_1)*norm(a_2)) )

   end subroutine get_lattice_params

   !
   ! Centre_Of_Mass
   !
   !% Calculate the centre of mass of an atoms object, using the closest images to the origin atom,
   !% or first atom if this is not specified.
   !% If an index_list is present, just calculate it for that subset of atoms (then the origin atom is
   !% the first in this list unless it is specified separately).
   !%
   !% Note: Because the origin can be specified separately it need not be one of the atoms in the 
   !% calculation.
   function centre_of_mass(at,index_list,origin,error) result(CoM)

     type(atoms),                      intent(in) :: at
     integer,                optional, intent(in) :: origin
     integer,  dimension(:), optional, intent(in) :: index_list
     real(dp), dimension(3)                       :: CoM
     integer, intent(out), optional :: error

     !local variables
     integer                                      :: i, my_origin
     real(dp)                                     :: M_Tot

     INIT_ERROR(error)
     if (.not. has_property(at, 'mass')) then
        RAISE_ERROR('center_of_mass: Atoms has no mass property', error)
     end if

     if (present(origin)) then
        if (origin > at%N .or. origin < 1) then
           RAISE_ERROR('Centre_Of_Mass: Invalid origin atom', error)
        end if
        my_origin = origin
     else
        if (present(index_list)) then
           my_origin = index_list(1)
        else
           my_origin = 1
        end if
     end if

     CoM = 0.0_dp
     M_Tot = 0.0_dp

     if (present(index_list)) then

        do i = 1, size(index_list)
           if (index_list(i) > at%N .or. index_list(i) < 1) then
              RAISE_ERROR('Centre_Of_Mass: Invalid atom in index_list', error)
           end if
           CoM = CoM + at%mass(index_list(i)) * diff_min_image(at,my_origin,index_list(i))
           M_Tot = M_Tot + at%mass(index_list(i))
        end do

     else

        do i = 1, at%N
           CoM = CoM + at%mass(i) * diff_min_image(at,my_origin,i)
           M_Tot = M_Tot + at%mass(i)
        end do

     end if

     CoM = CoM / M_Tot
     CoM = CoM + at%pos(:,my_origin)

   end function centre_of_mass

   !
   ! Directionality ellipsoid:
   !
   !% Given an origin atom and a list of other atoms, give information as to whether the other atoms
   !% are distributed roughly linearly, planar or spherically around the origin atom.
   !%
   !% The most notable use is to check that the splines in adjustable potential will be able to reproduce
   !% a randomly oriented force difference well.
   !%
   !% The information returned is the set of eigenvectors and associated eigenvalues of the directionality
   !% ellipsoid. One large e-value suggests roughly linear clustering, two similar and one small e-values suggest
   !% a planar distribution, while three similar e-values suggests almost spherical distribution (when copies of
   !% the atoms reflected through the origin atom are also considered).
   !% 
   !% To acheive a more spherical distribution, atoms along the e-vector(s) with the smallest e-value(s) should be
   !% added to the index list (See 'CosAngle_To_Line' below).
   !%
   !% The matrix which is diagonalised is an average of the outer products of the unit vectors from the origin
   !% atom to the other atoms.
   !%
   !% An outer product has 1 eigenvector which is the vector it was constructed from with
   !% eigenvalue 1 and the other eigenvectors have eigenvalue 0.
   !%
   !% The eigenvalues of the averaged matrix sum to 1.
   !% 
   subroutine directionality(this,origin,list,evalues,evectors,method,error)

     type(Atoms),                      intent(in)  :: this     !% The input atoms structure
     integer,                          intent(in)  :: origin   !% The origin atom
     type(table),                      intent(in)  :: list     !% Indices and shifts of the other atoms relative to origin
     real(dp), dimension(3),           intent(out) :: evalues  !% Eigenvalues of the directionality matrix
     real(dp), dimension(3,3),         intent(out) :: evectors !% Eigenvectors of the directionality matrix
     integer, optional,                intent(in)  :: method   !% 'METHOD = 1' Directionality ellipsoid method \\
                                                               !% 'METHOD = 2' Singular Value Decomposition method (default)
     integer, intent(out), optional :: error

     !local variables
     integer                                       :: i, j, k, l, n, my_method, lwork, info, jshift(3)
     real(dp), dimension(3)                        :: r_ij,rhat_ij
     real(dp), dimension(3,3)                      :: A, B 
     real(dp), allocatable, dimension(:,:)         :: vectors, u
     real(dp), allocatable, dimension(:)           :: work

     INIT_ERROR(error)
     my_method = 2

     if (present(method)) then
        if (method < 1 .or. method > 2) then
           write(line,'(a,i0,a)')'Directionality: Method = ',method,' does not exist'
           call system_abort(line)
        end if
        my_method = method
     end if

     if (list%intsize /= 4) then
        RAISE_ERROR('Directionality: list must have 4 int columns for indices and shifts', error)
     end if
     if (list%N == 0) then
        RAISE_ERROR('Directionality: list table has no entries', error)
     end if

     i = origin

     select case(my_method)

     case(1)        ! **** Directionality Ellipsoid method ****

        B = 0.0_dp

        !Construct the directionality matrix
        do n = 1, list%N
           j = list%int(1,n)
           jshift = list%int(2:4,n)
           if (j > this%N) then
              RAISE_ERROR('Directionality: Atom '//j//' is out of range ('//this%N//')', error)
           end if
           r_ij = diff(this,i,j,jshift)
           rhat_ij = r_ij / norm(r_ij)
           forall(k=1:3,l=1:3) A(k,l) = rhat_ij(k)*rhat_ij(l)
           B = B + A
        end do

        B = B / real(list%N,dp)

        !Find eigenvalues/eigenvectors of the directionality matrix
        call diagonalise(B,evalues,evectors)

     case(2)       ! **** Singular Value Decomposition method ****

        lwork = 2 * max(3*min(3,list%N)+max(3,list%N),5*min(3,list%N))

        allocate(work(lwork), vectors(list%N,3), u(list%N,3))

        !Fill 'vectors' with unit vectors
        do n = 1, list%N
           j = list%int(1,n)
           jshift = list%int(2:4,n)
           if (j > this%N) then
              RAISE_ERROR('Directionality: Atom '//j//' is out of range ('//this%N//')', error)
           end if
           r_ij = diff(this,i,j,jshift)
           vectors(n,:) = r_ij / norm(r_ij)
        end do   

        call dgesvd('N','A',list%N,3,vectors,list%N,evalues,u,list%N,evectors,3,work,lwork,info)
        !No left singular vectors, all right singular vectors

        if (info/=0) then
           if (info < 0) then
              RAISE_ERROR('Directionality: Problem with argument '//-info//' passed to DGESVD', error)
           else
              RAISE_ERROR('Directionality: DBDSQR (called from DGESVD) did not converge', error)
           end if
        end if

        deallocate(work, vectors, u)

        evectors = transpose(evectors)

     end select

   end subroutine directionality

   !
   ! CosAngle_To_Line
   !
   !% For use with the 'Directionality' routine above.
   !% Given an atom ('atom') and a direction ('dir') return the absolute value of the cosine of the angle
   !% between the the line running through 'atom' in direction 'dir' and the line
   !% between 'atom' and 'test_atom'
   !%
   !% The idea is that this will be called with the origin atom from 'Directionality' as 'atom',
   !% the eigenvector associated with the smallest eigenvalue as 'dir' and a potentially new 
   !% atom to connect to 'atom' with a spline as 'test_atom'.
   !% 
   !% If the result is close to 1 then accept the 'test_atom', and reject if close to zero

   function cosangle_to_line(this,atom,dir,test_atom,error)

     type(atoms),            intent(in) :: this
     integer,                intent(in) :: atom
     real(dp), dimension(3), intent(in) :: dir
     integer,                intent(in) :: test_atom
     real(dp)                           :: CosAngle_To_Line
     integer, intent(out), optional :: error
     !local variables
     real(dp), dimension(3)             :: r_ab

     INIT_ERROR(error)
     !Sanity checks
     if (atom > this%N) then
        RAISE_ERROR('CosAngle_To_Line: Atom '//atom//' out of range ('//this%N//')', error)
     end if

     if (test_atom > this%N) then
        RAISE_ERROR('CosAngle_To_Line: Test atom '//test_atom//' out of range ('//this%N//')', error)
     end if

     if (norm2(dir) .feq. 0.0_dp) then
        RAISE_ERROR('CosAngle_To_Line: A non-zero direction is required', error)
     end if

     r_ab = diff_min_image(this,atom,test_atom)

     CosAngle_To_Line = abs(dir .dot. r_ab) / (norm(dir)*norm(r_ab))

   end function cosangle_to_line

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !
   ! I/O procedures  
   !
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine atoms_print(this,file,error)
      type(Atoms),    intent(inout)             :: this
      type(Inoutput), intent(inout),optional, target :: file
      integer, intent(out), optional :: error

      type(Inoutput), pointer :: my_out

      INIT_ERROR(error)
      if(.not. is_initialised(this)) then
         RAISE_ERROR('Atoms_Print: Atoms structure not initialised', error)
      end if

      if(current_verbosity() <= PRINT_SILENT) return ! for efficiency

      if (present(file)) then
         my_out => file
      else
         my_out => mainlog
      end if

      call print('Atoms Structure: ', PRINT_NORMAL, my_out)
      call print('Number of atoms = '//this%N, PRINT_NORMAL, my_out)

      if(this%use_uniform_cutoff) then
         call print('Bond-formation cutoff radius = '//this%cutoff//' Angstroms', PRINT_NORMAL, my_out)
         call print('Bond-breaking cutoff radius = '//this%cutoff_break//' Angstroms', PRINT_NORMAL, my_out)
      else
         call print('Bond-formation cutoff radius = '//this%cutoff//' *bond_length', PRINT_NORMAL, my_out)
         call print('Bond-breaking cutoff radius = '//this%cutoff_break//' *bond_length', PRINT_NORMAL, my_out)
      end if

      call print('Lattice vectors:', PRINT_NORMAL, my_out)
      call print('a = ('//this%lattice(:,1)//')', PRINT_NORMAL, my_out)
      call print('b = ('//this%lattice(:,2)//')', PRINT_NORMAL, my_out)
      call print('c = ('//this%lattice(:,3)//')', PRINT_NORMAL, my_out)

      call print('Params')
      call print(this%params)

      call print('Properties')
      call print(this%properties)

      if (this%connect%initialised) then
	call verbosity_push_decrement()
	call connection_print(this%connect, my_out)
	call verbosity_pop()
      end if

      call print('',PRINT_NORMAL, my_out)
   end subroutine atoms_print

   function prop_names_string(this, with_types, error)
     type(Atoms), intent(in) :: this
     logical, optional, intent(in) :: with_types
     integer, intent(out), optional :: error
     character(len=2048) :: prop_names_string

     INIT_ERROR(error)
     prop_names_string=dict_prop_names_string(this%properties, with_types)
     PASS_ERROR(error)

   end function prop_names_string

   function dict_prop_names_string(this,with_types,error)
     type(Dictionary), intent(in) :: this
     logical, intent(in), optional :: with_types
     integer, intent(out), optional :: error
     character(len=2048) :: dict_prop_names_string

     character(len=1) :: prop_type
     character(len=1024) :: tmp
     integer :: i, n_cols, type
     logical :: my_with_types

     INIT_ERROR(error)
     my_with_types = optional_default(.false., with_types)

     dict_prop_names_string = ""
     do i=1,this%N
	if (my_with_types) then

           select case(this%entries(i)%type)
           case(T_INTEGER_A)
              prop_type = 'I'
              n_cols = 1
           case(T_REAL_A)
              prop_type = 'R'
              n_cols = 1
           case(T_LOGICAL_A)
              prop_type = 'L'
              n_cols = 1
           case(T_CHAR_A)
              prop_type = 'S'
              n_cols = 1
           case(T_INTEGER_A2)
              prop_type = 'I'
              n_cols = this%entries(i)%len2(1)
           case(T_REAL_A2)
              prop_type = 'R'
              n_cols = this%entries(i)%len2(1)
           case default
              RAISE_ERROR('dict_prop_names_string: bad property type='//type//' in properties argument', error)
           end select

	  write(tmp,'(i0)') n_cols ! Number of columns for this property
	  dict_prop_names_string=trim(dict_prop_names_string)//string(this%keys(i))//':'// &
				 prop_type//':'//trim(tmp)//':'
	else
	  dict_prop_names_string=trim(dict_prop_names_string)//string(this%keys(i))//':'
	endif
     end do
     ! Remove trailing ':'
     dict_prop_names_string = dict_prop_names_string(1:len_trim(dict_prop_names_string)-1)
   end function dict_prop_names_string

  subroutine connection_print(this,file,error)

    type(Connection), intent(in)    :: this
    type(Inoutput),   optional, intent(inout) :: file
    integer, intent(out), optional :: error
    integer                         :: i,j,k,n,m

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
                n = this%cell(i,j,k)%N / 10 ! integer division will round down
                do m = 1, n
                   write(line,'(10i7)') this%cell(i,j,k)%int(1,(10*m-9):(10*m))
                   call print(line,file=file)
                end do
                m = this%cell(i,j,k)%N - 10 * n
                if (m /= 0) then
                   write(line,'(10i7)') this%cell(i,j,k)%int(1,(10*n+1):(10*n+m))
                   call print(line,file=file)
                end if
                write(line,'(a70)')'----------------------------------------------------------------------'
                call print(line,file=file)
             end do
          end do
       end do

       call verbosity_pop()

    end if

  end subroutine connection_print

  !Bond_Length
  !% Returns the sum of the covalent radii of two atoms
  function bond_length(z1,z2)
    integer, intent(in) :: z1,z2
    real(dp)            :: bond_length
    bond_length = ElementCovRad(z1) + ElementCovRad(z2)
  end function bond_length

  !Termination_Bond_Rescale
  !% Calculates the rescale ratio of a Z1--H bond 
  !% generate from a Z1--Z2 bond.
  function termination_bond_rescale(z1,z2)
    integer,  intent(in) :: z1, z2
    real(dp)             :: termination_bond_rescale
    termination_bond_rescale = (ElementCovRad(z1) + ElementCovRad(1)) &
         / (ElementCovRad(z1) + ElementCovRad(z2))
  end function termination_bond_rescale

  !% Parses an atom_mask, which is string consisting of the '@' symbol followed by a comma separated
  !% list of indices or ranges into a table containing all the indices it represents.
  !% E.g. '@1,37-39,54,99-102' is expanded to a table with 1, 37, 38, 39, 54, 99, 100,
  !% 101, 102 as its first integer column. There must be no spaces in the mask.
  subroutine parse_atom_mask(mask_in,atom_indices,error)

    character(*),  intent(in)    :: mask_in
    type(Table),   intent(out)   :: atom_indices
    integer, intent(out), optional :: error

    character(len(mask_in))      :: mask
    character(20), dimension(50) :: fields
    integer                      :: Nfields, i, j, n, a, b, c

    INIT_ERROR(error)
    call allocate(atom_indices,1,0,0,0,100)

    mask = adjustl(mask_in)

    if (mask(1:1)=='@') then  !Found atom mask
       !Remove @ symbol
       mask = mask(2:)
       !Parse into comma separated numbers or ranges
       call parse_string(mask,',',fields,Nfields)
       if (Nfields==0) then
          RAISE_ERROR('parse_atom_mask: Atom mask contained no indices', error)
       end if
       do i = 1, Nfields
          !Is this field a single atom or a contiguous range?
          n = scan(fields(i),'-')
          if (n /= 0) then
             !it's a range. get start number
             a = string_to_int(fields(i)(:n-1))
             !get end number
             b = string_to_int(fields(i)(n+1:))
             if (a > b) then
                c = b
                b = a
                a = c
             end if
             !add all these atoms to the list
             do j = a, b
                if (.not.is_in_array(int_part(atom_indices,1),j)) call append(atom_indices,j)
             end do
          else
             !it's a single atom.
             a = string_to_int(fields(i))
             if (.not.is_in_array(int_part(atom_indices,1),j)) call append(atom_indices,a)
          end if
       end do
    else
       RAISE_ERROR('parse_atom_mask: Invalid atom mask: '//mask_in, error)
    end if

  end subroutine parse_atom_mask

  !% Find atoms which have integer property 'prop' set to
  !% true (i.e. set to 1) and return them in the Table 'list'.
  subroutine property_to_list(at, prop, list, error)
    type(Atoms), intent(in) :: at
    character(*), intent(in) :: prop
    type(Table), intent(inout) :: list
    integer, intent(out), optional :: error
    
    integer :: i
    integer, dimension(:), pointer :: p

    INIT_ERROR(error)
    if (.not. assign_pointer(at, prop, p)) then
       RAISE_ERROR('property_to_list: at does not have property "'//trim(prop)//'".', error)
    end if

    call allocate(list, 1, 0, 0, 0)
    call append(list, pack((/ (i, i=1,size(p)) /), p == 1))

  end subroutine property_to_list

  !% Convert the Table 'list' to a single column integer property in 'at',
  !% with atoms in list marked with a 1 and absent 
  subroutine list_to_property(at, list, prop, error)
    type(Atoms), intent(inout) :: at
    type(Table), intent(in) :: list
    character(*), intent(in) :: prop
    integer, intent(out), optional :: error

    integer, dimension(:), pointer :: p

    INIT_ERROR(error)
    ! Add property if necessary
    call add_property(at, prop, 0)

    if (.not. assign_pointer(at, prop, p)) then
       RAISE_ERROR('list_to_property: at does not have property "'//trim(prop)//'".', error)
    end if

    p = 0 ! Set all entries to zero
    p(int_part(list,1)) = 1 ! Set entries in 'list' to 1
    
  end subroutine list_to_property

  !% Find atoms which have integer property 'prop' with value 'value'
  !% and return them in a table 'list'.
  subroutine list_matching_prop(at,list,name,value,error)

    type(atoms), intent(in)    :: at
    type(table), intent(inout) :: list
    character(*), intent(in)   :: name
    integer,     intent(in)    :: value
    integer, intent(out), optional :: error

    integer                    :: i
    integer, pointer, dimension(:) :: ptr

    INIT_ERROR(error)
    !find property
    if (.not. assign_pointer(at%properties,name,ptr)) then
       RAISE_ERROR('Property "'//name//'" not found', error)
    end if

    call wipe(list)

    do i = 1, at%N
       if (ptr(i)==value) call append(list,(/i/))
    end do

  end subroutine list_matching_prop

  !% Return the complement of a list, i.e. all those atoms not included
  !% in list. Result is in outlist on exit.
  subroutine complement(at, inlist, outlist)
    type(Atoms), intent(in) :: at
    type(Table), intent(in) :: inlist
    type(Table), intent(out) :: outlist

    integer :: i
    integer, allocatable, dimension(:) :: inarray
    
    call table_allocate(outlist,1,0,0,0)
    
    allocate(inarray(inlist%N))
    inarray = int_part(inlist,1)
    
    do i=1,at%N
       if (.not. is_in_array(inarray,i)) &
            call append(outlist,i)
    end do

    deallocate(inarray)
    
  end subroutine complement


  !% Return the difference between list1 and list2 in outlist.
  !% That is, those elements in list1 but not in list2
  subroutine difference(list1, list2, outlist, error)
    type(Table), intent(in) :: list1, list2
    type(Table), intent(out) :: outlist
    integer, intent(out), optional :: error

    integer :: i
    integer, dimension(:), allocatable :: array1, array2

    INIT_ERROR(error)
    if (list1%N <= list2%N) then
       RAISE_ERROR('difference: list1%N ('//(list1%N)//') <= list2%N ('//(list2%N)//').', error)
    end if

    call table_allocate(outlist, 1, 0, 0, 0)

    allocate(array1(list1%N),array2(list2%N))
    array1 = int_part(list1,1)
    array2 = int_part(list2,1)

    do i=1,list1%N
       if (.not. is_in_array(array2, array1(i))) &
            call append(outlist,list1%int(1,i))
    end do

    deallocate(array1,array2)

  end subroutine difference

  !% move atoms around following neighbor list bonds so that all are in the same periodic image
  !%    (that of 'seed', if present)
  !% poorly tested, especially for situations where not all atoms are in one connected clump
  !% probably needs a better subroutine name
  subroutine coalesce_in_one_periodic_image(this, seed, is_periodic, error)
    type(Atoms), intent(inout) :: this
    integer, intent(in), optional :: seed
    logical, optional, intent(in) :: is_periodic(3)
    integer, intent(out), optional :: error

    integer :: i, ji, jji, j, shift(3), delta_shift(3), jj, k, ki
    integer :: n_neighbours, max_n_neighbours
    integer, allocatable :: shifts(:,:,:), cluster_list(:)
    logical, allocatable :: touched(:), is_in_cluster(:)
    integer :: seed_val, last_n_in_cluster
    logical :: dir_mask(3)
    integer, allocatable, dimension(:) :: dir_indices

    INIT_ERROR(error)
    seed_val=optional_default(1, seed)

    if (present(is_periodic)) then
       dir_mask = .not. is_periodic
       allocate(dir_indices(count(dir_mask)))
       dir_indices(:) = pack((/1,2,3/), dir_mask)
    else
       allocate(dir_indices(3))
       dir_mask = .true.
       dir_indices = (/1,2,3/)
    end if

    max_n_neighbours = 0
    do i=1, this%N
      n_neighbours = atoms_n_neighbours(this, i)
      if (n_neighbours > max_n_neighbours) max_n_neighbours = n_neighbours
    end do
    allocate(shifts(3,max_n_neighbours,this%N))
    shifts = 0

    ! find which atoms are in cluster that includes seed
    allocate(is_in_cluster(this%N))
    is_in_cluster = .false.
    last_n_in_cluster = 0
    is_in_cluster(seed_val) = .true.
    do while (count(is_in_cluster) /= last_n_in_cluster)
      last_n_in_cluster = count(is_in_cluster)
      do i=1, this%N
	if (is_in_cluster(i)) then
	  do ji=1, atoms_n_neighbours(this, i)
	    j = atoms_neighbour(this, i, ji, shift=shift)
	    is_in_cluster(j) = .true.
	  end do
	end if
      end do
    end do
    allocate(cluster_list(count(is_in_cluster)))
    ji = 0
    do i=1, this%N
      if (is_in_cluster(i)) then
	ji = ji + 1
	cluster_list(ji) = i
      endif
    end do

    ! initialize shifts
    do i=1, this%N
      do ji=1, atoms_n_neighbours(this, i)
	j = atoms_neighbour(this, i, ji, shift=shift)
	shifts(:,ji,i) = shift
      end do
    end do

    allocate(touched(this%N))
    touched = .false.
    touched(seed_val) = .true.

    do while (any(shifts(dir_indices,:,cluster_list(:)) /= 0))
      ! look for atoms i that have been touched
      do i=1, this%N
	if (.not. is_in_cluster(i)) cycle
	if (touched(i)) then
	  ! look at neighbors of i
	  do ji=1, atoms_n_neighbours(this,i)
	    j = atoms_neighbour(this, i, ji)
	    ! if neighbor has 0 shift, it doesn't need to move
	    if (any(shifts(dir_indices,ji,i) /= 0)) then
	      ! atom j does need to move
	      ! atoms with non-zero shift should not have been touched before
	      if (touched(j)) then
                 RAISE_ERROR("undo_pbcs tried to move atom " // j // " twice", error)
	      endif

	      ! shift atom to zero out shift
              do k=1,3
                 if (dir_mask(k)) this%pos(:,j) = this%pos(:,j) + shifts(k,ji,i) * this%lattice(:,k)
              end do

	      ! fix shifts of j's neighbours
	      delta_shift = merge(shifts(:,ji,i), (/0, 0, 0/), dir_mask)
	      do jji=1, atoms_n_neighbours(this, j)
		! fix shifts from j to its neighbors
		jj = atoms_neighbour(this, j, jji)
		shifts(:,jji,j) = shifts(:,jji,j) + delta_shift(:)
		! fix shifts from j's neighbours to it
		do ki=1, atoms_n_neighbours(this, jj)
		  k = atoms_neighbour(this, jj, ki)
		  if (k == j) then
		    shifts(:,ki,jj) = shifts(:,ki,jj) - delta_shift(:)
		  endif
		end do ! ki
	      end do ! jji
	      ! shifts(:,ji,i) = 0
	      touched(j) = .true.
	    else
	      ! atom j doesn't need to move
	      touched(j) = .true.
	    endif ! any(shift) /= 0
	  end do ! ji
	endif ! touched i
      end do ! i
    end do ! some shift isn't zero

    deallocate(touched)
    deallocate(shifts)
    deallocate(is_in_cluster)
    deallocate(cluster_list)
    deallocate(dir_indices)

  end subroutine coalesce_in_one_periodic_image

  function closest_atom(this, r, cell_image_Na, cell_image_Nb, cell_image_Nc, mask, dist, diff, error)
    type(Atoms), intent(in) :: this
    real(dp), intent(in) :: r(3)
    integer, intent(in) :: cell_image_Na, cell_image_Nb, cell_image_Nc
    logical, intent(in), optional, dimension(:) :: mask
    real(dp), intent(out), optional :: dist, diff(3)
    integer :: closest_atom
    integer, intent(out), optional :: error

    integer :: i, j, k
    integer i2, j2, k2, i3, j3, k3, i4, j4, k4, n2, atom_i
    integer :: cellsNa, cellsNb, cellsNc
    real(dp) :: pos(3), cur_dist, cur_diff(3), min_dist, min_diff(3)
    integer :: min_cell_image_Na, max_cell_image_Na, min_cell_image_Nb, max_cell_image_Nb, min_cell_image_Nc, max_cell_image_Nc

    INIT_ERROR(error)
    if (.not. this%connect%initialised) then
       RAISE_ERROR("closest_atom must have initialised connection object", error)
    end if

    call cell_of_pos(this%connect, this%g .mult. r, i, j, k)

    cellsNa = this%connect%cellsNa
    cellsNb = this%connect%cellsNb
    cellsNc = this%connect%cellsNc

    call get_min_max_images(this%is_periodic, cellsNa, cellsNb, cellsNc,&
      cell_image_Na, cell_image_Nb, cell_image_Nc, i, j, k, .true., .true., .true., &
      min_cell_image_Na, max_cell_image_Na, min_cell_image_Nb, max_cell_image_Nb, min_cell_image_Nc, max_cell_image_Nc)

    k3 = 1; k4 = 1; j3 = 1; j4 = 1; i3 = 1; i4 = 1

    min_dist = 1.0e38_dp
    closest_atom = 0
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

	     do n2 = 1, this%connect%cell(i3,j3,k3)%N
		atom_i = this%connect%cell(i3,j3,k3)%int(1,n2)
                if (present(mask)) then
                   if (.not. mask(atom_i)) cycle
                end if
		pos = this%pos(:,atom_i) + ( this%lattice .mult. (/ i4, j4, k4 /) )
                cur_diff = pos - r
		cur_dist = norm(cur_diff)
		if (cur_dist < min_dist) then
		  min_dist = cur_dist
                  min_diff = cur_diff
		  closest_atom = atom_i
		endif

	     end do

	  end do
       end do
    end do

    if (present(dist)) dist = min_dist
    if (present(diff)) diff = min_diff

  end function closest_atom

  function is_min_image(this, i, alt_connect, error) 
    type(Atoms), target, intent(in) :: this
    integer, intent(in) ::i
    type(Connection), intent(inout), target, optional :: alt_connect
    integer, intent(out), optional :: error

    logical :: is_min_image
    integer :: n, m, NN
    type(Connection), pointer :: use_connect

    INIT_ERROR(error)
    if (present(alt_connect)) then
      use_connect => alt_connect
    else
      use_connect => this%connect
    endif

    is_min_image = .true.
    ! First we give the neighbour1 (i <= j) then the neighbour2 entries (i > j) 
    if (use_connect%initialised) then

       if (.not. associated(use_connect%neighbour1(i)%t)) then
          RAISE_ERROR('is_min_image: atoms structure has no connectivity data for atom '//i, error)
       end if

       nn = use_connect%neighbour1(i)%t%N
       do n=1,nn
          if (use_connect%neighbour1(i)%t%int(1,n) == i) then
             is_min_image = .false.
             return
          end if
          do m=n+1,nn
             if (use_connect%neighbour1(i)%t%int(1,n) == use_connect%neighbour1(i)%t%int(1,m)) then
                is_min_image = .false.
                return
             end if
          end do
       end do

       nn = use_connect%neighbour2(i)%t%N
       do n=1,nn
          if (use_connect%neighbour2(i)%t%int(1,n) == i) then
             is_min_image = .false.
             return
          end if
          do m=n+1,nn
             if (use_connect%neighbour2(i)%t%int(1,n) == use_connect%neighbour2(i)%t%int(1,m)) then
                is_min_image = .false.
                return
             end if
          end do
       end do

    else
       RAISE_ERROR('is_min_image: Atoms structure has no connectivity data. Call calc_connect first.', error)
    end if

  endfunction is_min_image 

  subroutine atoms_bcast(mpi, at, error)
    type(MPI_context), intent(in) :: mpi
    type(Atoms), intent(inout) :: at
    integer, optional, intent(out) :: error

    character, allocatable, dimension(:) :: char_array
    integer, parameter :: SIZEOF_ATOMS = 1760

    INIT_ERROR(error)

#ifdef __GFORTRAN__
    ! Raise an error if sizeof(Atoms) has changed, indicating fields
    ! have been added or removed from definition of derived type.
    if (size(transfer(at, char_array)) /= SIZEOF_ATOMS) then
       RAISE_ERROR('atoms_bcast: size of Atoms object ('//size(transfer(at, char_array))//' /= '//SIZEOF_ATOMS//' - please make sure atoms_bcast() is up to date, sharing all variables if any new were added', error)
    end if
#endif

    if (.not. mpi%active) return

    if (mpi%my_proc == 0) then

       call print('atoms_bcast: bcasting from  proc '//mpi%my_proc, PRINT_VERBOSE)
       call bcast(mpi, at%n)
       call bcast(mpi, at%Ndomain)
       call bcast(mpi, at%Nbuffer)
       call bcast(mpi, at%domain_decomposed)
       call bcast(mpi, at%use_uniform_cutoff)
       call bcast(mpi, at%cutoff)
       call bcast(mpi, at%cutoff_break)
       call bcast(mpi, at%nneightol)
       call bcast(mpi, at%lattice)
       call bcast(mpi, at%is_orthorhombic)
       call bcast(mpi, at%is_periodic)
       call bcast(mpi, at%fixed_size)
       call bcast(mpi, at%properties)
       call bcast(mpi, at%params)
    else
       call print('atoms_bcast: bcasting to  proc '//mpi%my_proc, PRINT_VERBOSE)
       call finalise(at)
       call bcast(mpi, at%n)
       call bcast(mpi, at%Ndomain)
       call bcast(mpi, at%Nbuffer)
       call bcast(mpi, at%domain_decomposed)
       call bcast(mpi, at%use_uniform_cutoff)
       call bcast(mpi, at%cutoff)
       call bcast(mpi, at%cutoff_break)
       call bcast(mpi, at%nneightol)
       call bcast(mpi, at%lattice)
       call bcast(mpi, at%is_orthorhombic)
       call bcast(mpi, at%is_periodic)
       call bcast(mpi, at%fixed_size)
       call bcast(mpi, at%properties)
       call bcast(mpi, at%params)

       call matrix3x3_inverse(at%lattice,at%g)
       call atoms_repoint(at)
       at%ref_count = 1
    end if

  end subroutine atoms_bcast


  !% Find the indices of the atoms within the cone with its point at the atom 'origin',
  !% defined by direction 'dir' and opening angle with cosine 'cos_theta'. The indices
  !% are returned in the Table 'output' which has a single integer column.
  subroutine atoms_filter_cone(this, origin, dir, cos_theta, output)
    type(Atoms), intent(in) :: this
    integer, intent(in) :: origin
    real(dp), intent(in), dimension(3) :: dir
    real(dp), intent(in) :: cos_theta
    type(Table), intent(out) :: output

    real(dp) :: d(3), ndir(3), p
    integer :: i

    ndir = dir/norm(dir)
    p = ndir .dot. this%pos(:,origin)

    call allocate(output, 1,0,0,0)
    do i=1,this%n
       if (i == origin) cycle
       d = this%pos(:,i) - this%pos(:, origin)
       if ((d .dot. ndir) < p .and. (d .dot. ndir)/(norm(d)) > cos_theta) call append(output, i)
    end do

  end subroutine atoms_filter_cone

  !% Copy some properties from one atoms struct to another
  !% The destination will be overriden.
  subroutine atoms_copy_properties(this, from, property_list, case_sensitive, error)
    type(Atoms), intent(inout) :: this
    type(Atoms), intent(in) :: from
    character(len=*), intent(in) :: property_list
    logical, optional, intent(in) :: case_sensitive
    integer, optional, intent(out) :: error

    character(len=len_trim(property_list)) :: property_a(100)
    integer :: property_a_n

    INIT_ERROR(error)

    call split_string_simple(trim(property_list), property_a, property_a_n, ':', error)
    PASS_ERROR(error)
    call subset(from%properties, property_a(1:property_a_n), this%properties, case_sensitive=case_sensitive, out_no_initialise=.true., error=error)
    PASS_ERROR(error)

  end subroutine atoms_copy_properties

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
	  call atoms_copy_entry(this, cur_place, smallest_i_a, swap=.true., error=error)
	  PASS_ERROR(error)
       endif
    end do
  end subroutine atoms_sort

  !% Convert from a single index in range 1..this%N to a CASTEP-style (element, index) pair
  function atoms_index_to_z_index(this, index) result(z_index)
    type(Atoms), intent(in) :: this
    integer, intent(in) :: index
    integer :: z_index
    integer :: j

    z_index = 0
    do j=1,index
       if (this%z(index) == this%z(j)) z_index = z_index + 1
    end do
  
  end function atoms_index_to_z_index

  !% Inverse of atoms_index_to_z_index
  function atoms_z_index_to_index(this, z, z_index, error) result(index)
    type(Atoms), intent(in) :: this
    integer, intent(in) :: z, z_index
    integer, intent(out), optional :: error
    integer :: index
    integer :: nz
    
    INIT_ERROR(error)

    nz = 0
    do index=1,this%N
       if(this%z(index) == z) nz = nz + 1
       if (nz == z_index) return
    end do
    
    RAISE_ERROR('atoms_z_index_to_index: index pair ('//z//','//z_index//') not found', error)

  end function atoms_z_index_to_index

end module atoms_module
