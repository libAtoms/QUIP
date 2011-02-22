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

  use error_module
  use system_module
  use mpi_context_module
  use linearalgebra_module
  use extendable_str_module
  use dictionary_module
  use table_module
  use paramreader_module
  use periodictable_module
  use Atoms_types_module
  use Connection_module
  use DomainDecomposition_module
  use minimization_module
  use Quaternions_module

  implicit none

  integer,  parameter :: NOT_NEIGHBOUR = 0     !% Returned by 'Find_Neighbour' if unsuccessful

  logical :: printed_stack_warning = .false.   !% OMIT

  private :: ess_max_len

  private :: atoms_initialise
  interface initialise
     module procedure atoms_initialise
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

  private :: atoms_is_domain_decomposed
  interface is_domain_decomposed
     module procedure atoms_is_domain_decomposed
  endinterface is_domain_decomposed

  private :: atoms_shallowcopy
  interface shallowcopy
     module procedure atoms_shallowcopy
  end interface shallowcopy

  !% Free up the memory associated with one or more objects.
  private :: atoms_finalise,atoms_finalise_multi
  interface finalise
     module procedure atoms_finalise, atoms_finalise_multi
  end interface finalise

  !% Finalise type(Atoms), pointer objects. Shallow copies of these will
  !% The object will when ref_count == 0.
  private :: atoms_finalise_ptr
  interface finalise_ptr
     module procedure atoms_finalise_ptr
  endinterface finalise_ptr

  private :: atoms_zero
  interface zero
     module procedure atoms_zero
  end interface

  !% Overloaded assigment operators for Atoms objects.
  private :: atoms_assignment
  interface assignment(=)
     module procedure atoms_assignment
  end interface assignment(=)
  interface deepcopy
     module procedure atoms_assignment
  endinterface

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
  private :: add_atom_single, add_atom_multiple, atoms_join
  interface add_atoms
     module procedure add_atom_single, add_atom_multiple, atoms_join
  end interface add_atoms

  !% Remove one or more atoms from an Atoms object.
  private :: remove_atom_single, remove_atom_multiple
  interface remove_atoms
     module procedure remove_atom_single, remove_atom_multiple
     module procedure remove_atom_multiple_mask
  end interface remove_atoms

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

  !% Print a verbose textual description of an Atoms object to the default logger or to
  !% a specificied Inoutput object.
  private :: atoms_print
  interface print
     module procedure atoms_print
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
  interface cell_volume
    module procedure atoms_cell_volume
  end interface

  private :: atoms_map_into_cell
  interface map_into_cell
    module procedure atoms_map_into_cell
  end interface

  private :: atoms_bcast
  interface bcast
     module procedure atoms_bcast
  end interface

  private :: atoms_copy_properties
  interface copy_properties
     module procedure atoms_copy_properties
  endinterface copy_properties

  private :: atoms_transform_basis
  interface transform_basis
     module procedure atoms_transform_basis
  end interface transform_basis

  private :: atoms_rotate
  interface rotate
     module procedure atoms_rotate
  end interface rotate

  private :: atoms_index_to_z_index
  interface index_to_z_index
     module procedure atoms_index_to_z_index
  end interface index_to_z_index

  private :: atoms_z_index_to_index
  interface z_index_to_index
     module procedure atoms_z_index_to_index
  end interface z_index_to_index

  private :: atoms_calc_connect
  interface calc_connect
     module procedure atoms_calc_connect
  endinterface

  private :: atoms_calc_connect_hysteretic
  interface calc_connect_hysteretic
     module procedure atoms_calc_connect_hysteretic
  endinterface

  private :: atoms_is_min_image
  interface is_min_image
     module procedure atoms_is_min_image
  endinterface

  private :: atoms_set_comm_property
  interface set_comm_property
     module procedure atoms_set_comm_property
  endinterface set_comm_property

  private :: atoms_set_Zs
  interface set_Zs
     module procedure atoms_set_Zs
  endinterface

  private :: atoms_sort_by_rindex
  interface sort
     module procedure atoms_sort, atoms_sort_by_rindex
  endinterface

  private :: atoms_shuffle
  interface shuffle
     module procedure atoms_shuffle
  endinterface

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
    !call atoms_finalise(this)

    if (present(params)) then
       this%params = params ! deep copy
    else
       call initialise(this%params)
    end if

    call initialise(this%domain, error=error)
    PASS_ERROR(error)

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
 
    call set_comm_property(this, 'Z', &
         comm_atoms=.true., comm_ghosts=.true.)
    call set_comm_property(this, 'pos', &
         comm_atoms=.true., comm_ghosts=.true.)
    call set_comm_property(this, 'species', &
         comm_atoms=.true.)
 
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


  !% Is this atoms object domain decomposed?
  logical function atoms_is_domain_decomposed(this)
    type(Atoms), intent(in)  :: this
    atoms_is_domain_decomposed = this%domain%decomposed
  end function atoms_is_domain_decomposed


  !% Shallow copy of this Atoms object
  subroutine atoms_shallowcopy(this, from)
    !NB gfortran 4.3.4 seg faults if this is intent(out)
    !NB type(Atoms), pointer, intent(out)  :: this
    type(Atoms), pointer  :: this
    !NB
    type(Atoms), target                :: from

    ! ---

    this => from
    this%ref_count = this%ref_count + 1

    call print("atoms_shallowcopy: Created shallow copy, " // &
         "ref_count = " // this%ref_count, PRINT_ANAL)

  end subroutine atoms_shallowcopy


  subroutine atoms_finalise(this)
    type(Atoms), intent(inout) :: this

    this%ref_count = this%ref_count - 1
    ! Note: this has to be == 0, rather than <= 0. Otherwise it will fail if
    ! finalise is called more than once
    if (this%ref_count == 0) then

       call print("atoms_finalise: ref_count = 0, finalising", PRINT_ANAL)

       call finalise(this%domain)

       call finalise(this%properties)
       call finalise(this%params)

       ! Nullify pointers
       nullify(this%Z, this%travel, this%mass)
       nullify(this%move_mask, this%thermostat_region, this%damp_mask)
       nullify(this%pos, this%velo, this%acc, this%avgpos, this%oldpos, this%avg_ke)

       call finalise(this%connect)
       call finalise(this%hysteretic_connect)

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

    call deepcopy(to%connect, from%connect)
    call deepcopy(to%hysteretic_connect, from%hysteretic_connect)
    call deepcopy(to%domain, from%domain)

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
    ASSERT(is_initialised(from), "atoms_copy_without_connect: 'from' object is not initialised", error)

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
    to%use_uniform_cutoff = from%use_uniform_cutoff
    to%cutoff = from%cutoff
    to%cutoff_break = from%cutoff_break
    to%nneightol = from%nneightol

    call deepcopy(to%domain, from%domain)

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
    type(Atoms), intent(in) :: this
    integer, intent(in) :: i
    real(dp), optional, intent(in) :: max_dist, max_factor
    type(Connection), optional, intent(in), target :: alt_connect
    integer, intent(out), optional :: error     
    integer :: n
    
    INIT_ERROR(error)
    if (present(alt_connect)) then
       n = n_neighbours(alt_connect, this, i, max_dist, max_factor, error)
       PASS_ERROR(error)
    else
       n = n_neighbours(this%connect, this, i, max_dist, max_factor, error)
       PASS_ERROR(error)
    endif

  end function atoms_n_neighbours

  function atoms_neighbour_index(this, i, n, index, t, is_j, alt_connect, error) result(j)
    type(Atoms), intent(in) :: this
    integer, intent(in) :: i, n
    integer,  intent(out) :: index
    type(Table), pointer, intent(out) :: t
    logical, intent(out) :: is_j
    type(Connection), optional, intent(in), target :: alt_connect
    integer, intent(out), optional :: error     
    integer :: j

    INIT_ERROR(error)
    if (present(alt_connect)) then
       j = neighbour_index(alt_connect, &
            i, n, index, t, is_j, error)
       PASS_ERROR(error)
    else
       j = neighbour_index(this%connect, &
            i, n, index, t, is_j, error)
       PASS_ERROR(error)
    endif

  end function atoms_neighbour_index

  function atoms_neighbour_minimal(this, i, n, shift, index, alt_connect) result(j)
    type(Atoms), intent(in), target :: this
    integer ::i, j, n
    integer,  intent(out) :: shift(3)
    integer,  intent(out) :: index
    type(Connection), optional, intent(in), target :: alt_connect

    if (present(alt_connect)) then
       j = neighbour_minimal(alt_connect, i, n, shift, index)
    else
       j = neighbour_minimal(this%connect, i, n, shift, index)
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
    type(Atoms), intent(in)         :: this
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

    INIT_ERROR(error)
    if (present(alt_connect)) then
       j = neighbour(alt_connect, this, &
            i, n, distance, diff, cosines, shift, index, max_dist, jn, error)
       PASS_ERROR(error)
    else
       j = neighbour(this%connect, this, &
            i, n, distance, diff, cosines, shift, index, max_dist, jn, error)
       PASS_ERROR(error)
    endif

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

    logical :: has_mass, has_velo, has_acc, has_travel

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

    ! If we pass a NULL() pointer as an optional argument, the argument
    ! shows up as present but has size == 0.

    has_mass = .false.
    has_travel = .false.
    has_velo = .false.
    has_acc = .false.

    if (present(mass)) then
       if (size(mass) > 0) then
          has_mass = .true.
          call check_size('Mass',mass,size(Z), 'Add_Atom', error)
          PASS_ERROR(error)
       endif
    end if

    if (present(travel)) then
       if (size(travel, 2) > 0) then
          has_travel = .true.
          call check_size('Travel',travel,(/3,size(Z)/),'Add_Atom', error)
          PASS_ERROR(error)
       endif
    end if

    if (present(velo)) then
       if (size(velo, 2) > 0) then
          has_velo = .true.
          call check_size('Velo', velo, (/3,size(Z)/), 'Add_Atom', error)
          PASS_ERROR(error)
       endif
    end if

    if (present(acc)) then
       if (size(acc, 2) > 0) then
          has_acc = .true.
          call check_size('Acc', acc, (/3,size(Z)/), 'Add_Atom', error)
          PASS_ERROR(error)
       endif
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
    endif

    ! We need to repoint even if the buffers are not reallocated. Otherwise
    ! ifort will complain the array sizes are too small if compiled with
    ! -check bounds
    call atoms_repoint(this)

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

    if (has_travel) then
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
       if (has_mass) then
          this%mass(oldN+1:this%N) = mass
       else
          this%mass(oldN+1:this%N) = ElementMass(Z)
       end if
    else if (has_mass) then
       ! mass specified but property doesn't yet exist, so create it...
       call add_property(this, 'mass', ElementMass(this%Z), ptr=this%mass, error=error)
       PASS_ERROR(error)
       call set_comm_property(this, 'mass', &
            comm_atoms=.true.)
       ! ... and then override for new atoms
       this%mass(oldN+1:this%N) = mass
    end if

    if (.not. has_key(this%properties, 'pos')) then
       RAISE_ERROR('Atoms_Add: this atoms has no pos property', error)
    end if
    this%pos(:,oldN+1:this%N) = pos
    
    if (has_velo .and. has_key(this%properties, 'velo')) &
         this%velo(:,oldN+1:this%N) = velo

    if (has_acc .and. has_key(this%properties, 'acc')) &
         this%acc(:,oldN+1:this%N) = acc

    call finalise(this%connect)

  end subroutine add_atom_multiple


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_join(this, from, error)
    type(Atoms),       intent(inout)  :: this
    type(Atoms),       intent(inout)  :: from
    integer, optional, intent(out)    :: error

    ! ---

    INIT_ERROR(error)

    call atoms_repoint(from)

    call add_atom_multiple(this, from%pos, from%Z, from%mass, from%velo, from%acc, from%travel, error)
    PASS_ERROR(error)

  endsubroutine atoms_join


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


  !% Reshuffle the order of the atomic indices to new_indices.
  subroutine atoms_shuffle(this, new_indices, error)
    type(Atoms),                intent(inout)  :: this
    integer, dimension(this%N), intent(in)     :: new_indices
    integer, optional,          intent(out)    :: error

    ! ---

    integer, allocatable, dimension(:)      :: tmp_int
    integer, allocatable, dimension(:,:)    :: tmp_int2
    real(dp), allocatable, dimension(:)     :: tmp_real
    real(dp), allocatable, dimension(:,:)   :: tmp_real2
    logical, allocatable, dimension(:)      :: tmp_logical
    character, allocatable, dimension(:,:)  :: tmp_char

    integer  :: i

    ! ---

    INIT_ERROR(error)

    ! Resize property data arrays, copying old data.
    ! this will break any existing pointers so we call atoms_repoint() 
    ! immediately after (note that user-held pointers will stay broken, there
    ! is no way to fix this since Fortran does not allow pointer-to-pointer
    ! types)
    do i=1,this%properties%N
       select case (this%properties%entries(i)%type)

       case(T_INTEGER_A)
          allocate(tmp_int(this%n))
          tmp_int(:) = this%properties%entries(i)%i_a(new_indices)
          call set_value(this%properties, string(this%properties%keys(i)), tmp_int)
          deallocate(tmp_int)

       case(T_REAL_A)
          allocate(tmp_real(this%n))
          tmp_real(:) = this%properties%entries(i)%r_a(new_indices)
          call set_value(this%properties, string(this%properties%keys(i)), tmp_real)
          deallocate(tmp_real)

       case(T_LOGICAL_A)
          allocate(tmp_logical(this%n))
          tmp_logical(:) = this%properties%entries(i)%l_a(new_indices)
          call set_value(this%properties, string(this%properties%keys(i)), tmp_logical)
          deallocate(tmp_logical)

       case(T_INTEGER_A2)
          allocate(tmp_int2(this%properties%entries(i)%len2(1),this%n))
          tmp_int2(:,:) = this%properties%entries(i)%i_a2(:,new_indices)
          call set_value(this%properties, string(this%properties%keys(i)), tmp_int2)
          deallocate(tmp_int2)

       case(T_REAL_A2)
          allocate(tmp_real2(this%properties%entries(i)%len2(1),this%n))
          tmp_real2(:,:) = this%properties%entries(i)%r_a2(:,new_indices)
          call set_value(this%properties, string(this%properties%keys(i)), tmp_real2)
          deallocate(tmp_real2)

       case(T_CHAR_A)
          allocate(tmp_char(this%properties%entries(i)%len2(1),this%n))
          tmp_char(:,:) = this%properties%entries(i)%s_a(:,new_indices)
          call set_value(this%properties, string(this%properties%keys(i)), tmp_char)
          deallocate(tmp_char)

       case default
          RAISE_ERROR('remove_atom_multiple: bad property type '//this%properties%entries(i)%type//' key='//this%properties%keys(i), error)
       end select
    end do
    call atoms_repoint(this)

  endsubroutine atoms_shuffle


  subroutine remove_atom_multiple(this, atom_indices, error)
    type(Atoms), intent(inout)             :: this
    integer,     intent(in), dimension(:)  :: atom_indices
    integer,     intent(out), optional     :: error

    integer i, copysrc
    integer, allocatable, dimension(:) :: new_indices
    integer, dimension(size(atom_indices)) :: sorted
    integer, dimension(:), allocatable :: uniqed

    INIT_ERROR(error)

    if (this%fixed_size) then
       RAISE_ERROR("remove_atom_multiple: Atoms object cannot be resized (this%fixed_size = .true.)", error)
    end if

    !Delete the connection data because the atomic indices become mangled
    call finalise(this%connect)

    ! Permute new_indices, following algorithm in table_record_delete_multiple,
    ! so that atom ordering is same as it was under old scheme when properties were
    ! stored as columns in Table this%data.
    ! (we find first atomc to be removed and last atom to not be removed
    !  and swap them. Repeat until all atoms to be removed are at the end.)

    sorted = atom_indices     ! Get our own copy of the  indices so we can sort them
    call heap_sort(sorted)
    call uniq(sorted, uniqed) ! remove duplicates from sorted indices

    allocate(new_indices(this%N))
    do i=1,this%N
       new_indices(i) = i
    end do

    copysrc = this%N
    do i = size(uniqed), 1, -1
       if (uniqed(i) < copysrc) then
          new_indices(uniqed(i)) = new_indices(copysrc)
       else if (uniqed(i) > copysrc) then
          RAISE_ERROR("remove_atom_multiple: Fatal internal error: uniqed(i) > copysrc, should not happen", error)
       endif

       copysrc = copysrc - 1
    end do

    ! update N
    this%N = this%N - size(uniqed)
    this%Ndomain = this%N
    this%Nbuffer = this%N

    if (this%N /= copysrc) then
       RAISE_ERROR("remove_atom_multiple: Fatal internal error: this%N /= copysrc, should not happen", error)
    endif

    ! This will reallocate all buffer to the new size (i.e. this%N)
    call shuffle(this, new_indices(1:this%N), error=error)
    PASS_ERROR(error)

    deallocate(uniqed, new_indices)

  end subroutine remove_atom_multiple


  subroutine remove_atom_multiple_mask(this, mask, error)

    type(Atoms),                intent(inout)  :: this
    logical, dimension(this%N), intent(in)     :: mask
    integer, optional,          intent(out)    :: error

    ! ---

    integer               :: n, i
    integer, allocatable  :: atom_indices(:)
    
    ! ---

    INIT_ERROR(error)

    n = 0
    do i = 1, this%N
       if (mask(i))  n = n + 1
    enddo
    if (n > 0) then
       allocate(atom_indices(n))
       n = 0
       do i = 1, this%N
          if (mask(i)) then
             n = n + 1
             atom_indices(n) = i
          endif
       enddo
       call remove_atom_multiple(this, atom_indices, error)
       PASS_ERROR(error)
       deallocate(atom_indices)
    endif

  endsubroutine remove_atom_multiple_mask


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
       call set_comm_property(this, 'travel', &
            comm_atoms=.true.)
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
    call map_into_cell(dvw, this%lattice, this%g, i_shift)
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
             tmp = normsq(dvw + this%lattice(:,1)*i +this%lattice(:,2)*j + this%lattice(:,3)*k)
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

      do n = 1, this%N
         lat_pos = this%g .mult. this%pos(:,n)
         map_shift(:,n) = - floor(lat_pos+0.5_dp)
      end do
    end subroutine set_map_shift

   !
   !% Returns the (unsigned) volume of the simulation cell of 'this'
   !
   function atoms_cell_volume(this)
     type(Atoms), intent(in) :: this
     real(dp)                :: atoms_Cell_Volume

     atoms_cell_volume = cell_volume(this%lattice)
   end function atoms_cell_volume


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
   function centre_of_mass(at,index_list,mask,origin,error) result(CoM)

     type(atoms),                      intent(in)  :: at
     integer,                optional, intent(in)  :: origin
     integer,  dimension(:), optional, intent(in)  :: index_list
     logical,  dimension(:), optional, intent(in)  :: mask
     integer,                optional, intent(out) :: error
     real(dp), dimension(3)                        :: CoM

     !local variables
     integer                                      :: i, my_origin
     real(dp)                                     :: M_Tot

     INIT_ERROR(error)
     if (.not. has_property(at, 'mass')) then
        RAISE_ERROR('center_of_mass: Atoms has no mass property', error)
     end if

     if (present(index_list) .and. present(mask)) then
        RAISE_ERROR('centre_of_mass: Cannot take both index_list and mask arguments', error)
     endif

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

        if (present(mask)) then
           i = 1
           do while (.not. mask(i) .and. i < at%N)
              i = i+1
           enddo
           if (i >= at%N) then
              RAISE_ERROR('centre_of_mass: No atoms specified in mask.', error)
           endif
           my_origin = i
        else
           my_origin = 1
        endif
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

     else if (present(mask)) then

        do i = 1, at%N
           if (mask(i)) then
              CoM = CoM + at%mass(i) * diff_min_image(at,my_origin,i)
              M_Tot = M_Tot + at%mass(i)
           endif
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

     if (normsq(dir) .feq. 0.0_dp) then
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
	call print(this%connect, my_out)
	call verbosity_pop()
      end if

      call print('',PRINT_NORMAL, my_out)
   end subroutine atoms_print

   function prop_names_string(this, with_types, error)
     type(Atoms), intent(in) :: this
     logical, optional, intent(in) :: with_types
     integer, intent(out), optional :: error
     character(len=STRING_LENGTH) :: prop_names_string

     INIT_ERROR(error)
     prop_names_string=dict_prop_names_string(this%properties, with_types)
     PASS_ERROR(error)

   end function prop_names_string

   function dict_prop_names_string(this,with_types,error)
     type(Dictionary), intent(in) :: this
     logical, intent(in), optional :: with_types
     integer, intent(out), optional :: error
     character(len=STRING_LENGTH) :: dict_prop_names_string

     character(len=1) :: prop_type
     character(len=STRING_LENGTH) :: tmp
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

  function atoms_is_min_image(this, i, alt_connect, error) result(is_min_image)
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

  endfunction atoms_is_min_image

  subroutine atoms_bcast(mpi, at, error)
    type(MPI_context), intent(in) :: mpi
    type(Atoms), intent(inout) :: at
    integer, optional, intent(out) :: error

#ifdef __GFORTRAN__
    character, allocatable, dimension(:) :: char_array
    integer, parameter :: SIZEOF_ATOMS = 1776
#endif

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


  !% sort atoms according to an externally provided field
  subroutine atoms_sort_by_rindex(this, sort_index, error)
    type(Atoms),                 intent(inout)  :: this
    real(DP), dimension(this%N), intent(in)     :: sort_index
    integer,  optional,          intent(out)    :: error

    ! ---

    real(DP), allocatable  :: my_sort_index(:)
    integer, allocatable   :: atom_index(:)

    integer  :: i

    ! ---

    INIT_ERROR(error)

    allocate(my_sort_index(this%N), atom_index(this%N))

    do i = 1, this%N
       my_sort_index(i) = sort_index(i)
       atom_index(i) = i
    enddo

    call heap_sort(my_sort_index, i_data=atom_index)

    call shuffle(this, atom_index, error=error)
    PASS_ERROR(error)

    deallocate(my_sort_index, atom_index)

  endsubroutine atoms_sort_by_rindex


  !% Basis transformation of rank 0, 1 and 2 tensors real values in Atoms object.
  !% This routine transforms rank 1 and rank 2 tensors in this%params and
  !% this%properties. Tensors are identified by having the correct type
  !% (real arrays) and shape (i.e. 3, (3, 3), (3, this%N) (9, this%N) for
  !% vector paramters, tensor parameters, vector properties and tensor
  !% properties respectively), and by having a name which is included in
  !% the relevant list. Extra names can be added to the lists with the
  !% rank1 and rank2 arguments.
  subroutine atoms_transform_basis(this, L, rank1, rank2, error)
    type(Atoms), intent(inout) :: this
    real(dp), dimension(3,3) :: L
    character(len=*), optional, intent(in) :: rank1, rank2
    integer, optional, intent(out) :: error

    character(len=*), parameter :: default_rank1 = "pos:velo:acc:avgpos:oldpos:force:efield:dipoles"
    character(len=*), parameter :: default_rank2 = "virial:local_virial"

    character(len=STRING_LENGTH) :: all_rank1, all_rank2
    character(len=C_KEY_LEN) :: fields(100)
    integer i, j, n_fields
    type(Dictionary) :: rank1_d, rank2_d
    real(dp), pointer :: v(:), v2(:,:)

    INIT_ERROR(error)

    if (.not. is_orthogonal(L)) then
       call print('L=')
       call print(L)
       RAISE_ERROR("atoms_transform_basis: transformation matrix L is not orthogonal", error)
    end if

    all_rank1 = default_rank1
    if (present(rank1)) all_rank1 = trim(all_rank1)//':'//rank1
    call split_string_simple(all_rank1, fields, n_fields, ':', error)
    PASS_ERROR(error)
    call initialise(rank1_d)
    do i=1,n_fields
       call set_value(rank1_d, fields(i), .true.)
    end do

    all_rank2 = default_rank2
    if (present(rank2)) all_rank2 = trim(all_rank2)//':'//rank2
    call split_string_simple(all_rank2, fields, n_fields, ':', error)
    PASS_ERROR(error)
    call initialise(rank2_d)
    do i=1,n_fields
       call set_value(rank2_d, fields(i), .true.)
    end do

    ! Special case: transform lattice as 3 separate rank 1 tensors
    call set_lattice(this,(L .mult. this%lattice), scale_positions=.false.)

    ! First we consider entries in at%params
    do i=1, this%params%n
       ! Is it a real 3-vector?
       if (assign_pointer(this%params, string(this%params%keys(i)), v)) then
          if (has_key(rank1_d, string(this%params%keys(i))) .and. size(v) == 3) then
             call print('atoms_transform_basis: transforming '//this%params%keys(i)//' as a rank 1 tensor', PRINT_VERBOSE)
             v = L .mult. v
          end if
       end if

       ! Is it a real 3x3 matrix?
       if (assign_pointer(this%params, string(this%params%keys(i)), v2)) then
          if (has_key(rank2_d, string(this%params%keys(i))) .and. all(shape(v2) == (/3,3/))) then
             call print('atoms_transform_basis: transforming '//this%params%keys(i)//' as a rank 2 tensor', PRINT_VERBOSE)
             v2 = L .mult. v2 .mult. transpose(L)
          end if
       end if
    end do

    ! Now for the atomic properties in at%properties
    do i=1, this%properties%n
       ! Is it a per-atom 3-vector
       if (assign_pointer(this%properties, string(this%properties%keys(i)), v2)) then
          if (has_key(rank1_d, string(this%properties%keys(i))) .and. size(v2,1) == 3) then
             call print('atoms_transform_basis: transforming '//this%properties%keys(i)//' as a rank 1 tensor', PRINT_VERBOSE)
             do j=1,this%n
                v2(:,j) = L .mult. v2(:,j)
             end do
          end if
       end if
       
       ! Is it a per-atom 3x3 matrix, arranged as a 9 column property?
       if (assign_pointer(this%properties, string(this%properties%keys(i)), v2)) then
          if (has_key(rank2_d, string(this%properties%keys(i))) .and. size(v2, 1) == 9) then
             call print('atoms_transform_basis: transforming '//this%properties%keys(i)//' as a rank 2 tensor', PRINT_VERBOSE)
             do j=1,this%n
                v2(:,j) = reshape(L .mult. reshape(v2(:,j),(/3,3/)) .mult. transpose(L), (/9/))
             end do
          end if
       end if
    end do

    call finalise(rank1_d)
    call finalise(rank2_d)

  end subroutine atoms_transform_basis

  !% Rotate this Atoms object, transforming all rank 1 and rank 2 tensors parameters and properties
  subroutine atoms_rotate(this, axis, angle, rank1, rank2)
    type(Atoms), intent(inout) :: this
    real(dp), intent(in) :: axis(3), angle
    character(len=*), optional, intent(in) :: rank1, rank2

    type(Quaternion) :: q
    real(dp) :: R(3,3)

    q = rotation(axis, angle)
    R = rotation_matrix(q)
    call transform_basis(this, R, rank1, rank2)
    
  end subroutine atoms_rotate

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
  subroutine atoms_calc_connect_hysteretic(this, alt_connect, origin, extent, own_neighbour, store_is_min_image, error)
    type(Atoms), intent(inout), target           :: this
    type(Connection), intent(inout), target, optional :: alt_connect
    real(dp), optional :: origin(3), extent(3,3)
    logical, optional, intent(in) :: own_neighbour, store_is_min_image
    integer, intent(out), optional :: error

    INIT_ERROR(error)
    if (present(alt_connect)) then
       call calc_connect_hysteretic(alt_connect, this, &
            origin, extent, own_neighbour, store_is_min_image, error)
       PASS_ERROR(error)
    else
       call calc_connect_hysteretic(this%connect, this, &
            origin, extent, own_neighbour, store_is_min_image, error)
       PASS_ERROR(error)
    endif

  end subroutine atoms_calc_connect_hysteretic


  !% Fast $O(N)$ connectivity calculation routine. It divides the unit cell into similarly shaped subcells,
  !% of sufficient size that sphere of radius 'cutoff' is contained in a subcell, at least in the directions 
  !% in which the unit cell is big enough. For very small unit cells, there is only one subcell, so the routine
  !% is equivalent to the standard $O(N^2)$ method.
  subroutine atoms_calc_connect(this, alt_connect, own_neighbour, store_is_min_image, skip_zero_zero_bonds, error)
    type(Atoms),                intent(inout)  :: this
    type(Connection), optional, intent(inout)  :: alt_connect
    logical,          optional, intent(in)     :: own_neighbour, store_is_min_image, skip_zero_zero_bonds
    integer,          optional, intent(out)    :: error


    INIT_ERROR(error)
    if (present(alt_connect)) then
       call calc_connect(alt_connect, this, &
            own_neighbour, store_is_min_image, skip_zero_zero_bonds, error)
       PASS_ERROR(error)
    else
       call calc_connect(this%connect, this, &
            own_neighbour, store_is_min_image, skip_zero_zero_bonds, error)
    endif

  end subroutine atoms_calc_connect


  !% Set which properties to communicate when
  !% comm_atoms:   Communicate when atom is moved to different domain.
  !%               Forces, for example, may be excluded since they are updated
  !%               on every time step.
  !% comm_ghosts:  Communicate when atom is dublicated as a ghost on a domain.
  !%               Masses, for example, might be excluded since atoms are
  !%               propagated on the domain they reside in only.
  !% comm_reverse: Communicate back from ghost atoms to the original domain atom
  !%               and accumulate
  !% By default, properties are not communicated.
  subroutine atoms_set_comm_property(this, propname, &
       comm_atoms, comm_ghosts, comm_reverse)
    implicit none

    type(Atoms),       intent(inout)  :: this
    character(*),      intent(in)     :: propname
    logical, optional, intent(in)     :: comm_atoms
    logical, optional, intent(in)     :: comm_ghosts
    logical, optional, intent(in)     :: comm_reverse

    ! ---

    call set_comm_property(this%domain, propname, &
         comm_atoms, comm_ghosts, comm_reverse)

  endsubroutine atoms_set_comm_property


  !% set Zs from species
  subroutine atoms_set_Zs(this, error)
   type(Atoms), intent(inout) :: this
   integer, intent(out), optional :: ERROR

   integer i, ii

   INIT_ERROR(error)

   do i=1, this%N
      ii = find_in_array(ElementName, a2s(this%species(:,i)))
      if (ii < 1) then
	 RAISE_ERROR("failed to match name of atom "//i//" '"//a2s(this%species(:,i))//"'", error)
      endif
      this%Z(i) = ii-1
   end do
  end subroutine atoms_set_Zs

end module atoms_module
