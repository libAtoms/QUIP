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
!X  Group module
!X  
!% A Group is a collection of atoms, with an optional list of objects. These
!% objects can be anything, e.g. rigid bodies, constraints, residues... etc.
!% The idea is that different groups of atoms will require their equations of
!% motion integrating in different ways, so currently AdvanceVerlet loops over all
!% groups, checks the type of the group and integrates the equations of motion
!% for the atoms in that group in the appropriate way.
!% 
!% Two examples:
!% 
!% \begin{itemize}
!% \item A methane molecule with all C-H bonds fixed would have all 5 of its atoms in one
!% group, which would also contain (the indices of) four constraints. During
!% AdvanceVerlet the RATTLE algorithm can find the correct constraint objects to
!% enforce the fixed bond lengths
!% 
!% \item A rigid ammonia molecule could be treated as a rigid body, in which case the
!% group would contain the indices of all 4 of the atoms, and the index of 1 rigid body 
!% in the 'objects' array. Then, during AdvanceVerlet, the NO_SQUISH algorithm would
!% propagate these atoms rigidly all at the same time.
!% \end{itemize}
!% 
!% Initially, when a DynamicalSystem is initialised each atom is put into its own
!% group. As constraints and rigid bodies are added the number of groups remains
!% the same but some of them become empty. These empty groups can be used if extra
!% atoms are added to the system, or they can be left without a problem. If more
!% memory is needed then DS_Free_Groups can be called, which deletes all the empty groups
!% and recreates the group_lookup array inside DynamicalSystem. This is most useful
!% if all the constraints/rigid bodies are set up before a simulation and are known
!% to never change.
!% 
!% A few notes:
!% 
!% \begin{itemize}
!% \item Groups aren't necessarily limited to the usage described here
!% 
!% \item The 'atom' and 'object' arrays are self extensible, but are always exactly
!% the right size (unlike a Table). This is based on the assumption that they will
!% not change much (if at all) after they have been initially set up.
!% 
!% \item Each group has a 'type', which defines how its atoms are integrated. Currently
!% the types are 'TYPE_IGNORE' (numerically zero, and atoms in it are not integrated,
!% which could replace the move_mask variable(?)), 'TYPE_ATOM' (a normal atom),
!% 'TYPE_CONSTRAINED' and 'TYPE_RIGID'. The idea is that if other ways of integrating
!% atoms is required (e.g. linear molecules, positional restraints, belly atoms, etc.)
!% then an new type will be defined and the new will be added to the 'select case'
!% constructs in AdvanceVerlet.
!% \end{itemize}
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module group_module

  use system_module
  use linearalgebra_module

  implicit none
  private

  public :: group, print, initialise, finalise, set_type
  public :: group_n_objects, group_nth_object, group_n_atoms, group_nth_atom
  public :: num_free_groups, free_group, tidy_groups, groups_create_lookup
  public :: group_add_atom, group_delete_atom, group_add_object, merge_groups

  type Group

     integer                            :: type   !% a 'TYPE_' label for all the atoms in the group
     integer, allocatable, dimension(:) :: atom   !% a list of atoms in the group in ascending order
     integer, allocatable, dimension(:) :: object !% a list of objects these atoms are contained in, e.g.
                                                  !%   a single rigid body for the whole group,  or 
                                                  !%   the subset of the constraints involving these atoms
  end type Group

  interface print
     module procedure group_print, groups_print
  end interface print

  interface initialise
     module procedure group_initialise
  end interface initialise

  interface finalise
     module procedure group_finalise, groups_finalise
  end interface finalise

  interface assignment(=)
     module procedure group_assign_group, groups_assign_groups
  end interface assignment(=)

  interface set_type
    module procedure set_type_group
  end interface set_type

contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X Initialise
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine group_initialise(this,type,atoms,objects)

    type(Group),                     intent(inout) :: this
    integer,                         intent(in)    :: type
    integer, optional, dimension(:), intent(in)    :: atoms
    integer, optional, dimension(:), intent(in)    :: objects

    this%type = type
    
    if (allocated(this%atom)) deallocate(this%atom)
    if (present(atoms)) then
       if (size(atoms)>0) then
          allocate(this%atom(size(atoms)))
          this%atom = atoms
          call insertion_sort(this%atom)
       end if
    end if

    if (allocated(this%object)) deallocate(this%object)
    if (present(objects)) then
       if (size(objects)>0) then
          allocate(this%object(size(objects)))
          this%object = objects
          call insertion_sort(this%object)
       end if
    end if   

  end subroutine group_initialise

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X Finalise
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine group_finalise(this)

    type(Group), intent(inout) :: this

    this%type = 0

    if (allocated(this%atom)) deallocate(this%atom)
    if (allocated(this%object)) deallocate(this%object)

  end subroutine group_finalise

  !
  !% Finalise and deallocate an array of groups
  !
  subroutine groups_finalise(this)

    type(Group), dimension(:), allocatable, intent(inout) :: this
    integer                                               :: i
    
    if (allocated(this)) then
       do i = 1, size(this)
          call group_finalise(this(i))
       end do
       deallocate(this)
    end if

  end subroutine groups_finalise

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X Assignment
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine group_assign_group(to,from)

    type(Group), intent(inout) :: to
    type(Group), intent(in)    :: from

    call group_finalise(to)

    to%type = from%type
    if (allocated(from%atom)) then
       allocate(to%atom(size(from%atom)))
       to%atom = from%atom
    end if
    if (allocated(from%object)) then
       allocate(to%object(size(from%object)))
       to%object = from%object
    end if

  end subroutine group_assign_group

  subroutine groups_assign_groups(to,from)

    type(Group), dimension(:), intent(inout) :: to
    type(Group), dimension(:), intent(in)    :: from
    integer                                  :: i

    if (size(to) /= size(from)) call system_abort('Groups_Assign_Groups: Target and source sizes differ')

    do i = 1, size(from)
       to(i) = from(i)
    end do

  end subroutine groups_assign_groups

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X Setting the type
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine set_type_group(this,type)

    type(Group), intent(inout) :: this
    integer,     intent(in)    :: type

    this%type = type

  end subroutine set_type_group

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !% Merge: $g1 + g2 \to g1$. If $g2 = g1$ then the result is $g1$
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine merge_groups(g1,g2)

    type(Group), intent(inout) :: g1,g2
    integer                            :: type   !
    integer, allocatable, dimension(:) :: atom   ! Temporary storage
    integer, allocatable, dimension(:) :: object !
    integer                            :: atom_size, object_size, i

    if (g1%type /= g2%type) call system_abort('Merge_Groups: Groups are of different types')

    type = g1%type

    !Allocate the temporary 'atom' array to the max required
    atom_size = 0
    if (allocated(g1%atom)) atom_size = atom_size + size(g1%atom)
    if (allocated(g2%atom)) atom_size = atom_size + size(g2%atom)
    allocate(atom(atom_size))
    !Put each atom into the array only once, 
    atom_size = 0
    if (allocated(g1%atom)) then
       do i = 1, size(g1%atom)
          if (.not. is_in_array(atom(1:atom_size),g1%atom(i))) then
             atom_size = atom_size + 1
             atom(atom_size) = g1%atom(i)
          end if
       end do
    end if
    if (allocated(g2%atom)) then
       do i = 1, size(g2%atom)
          if (.not. is_in_array(atom(1:atom_size),g2%atom(i))) then
             atom_size = atom_size + 1
             atom(atom_size) = g2%atom(i)
          end if
       end do
    end if
    
    !Now do the same for 'object'
    object_size = 0
    if (allocated(g1%object)) object_size = object_size + size(g1%object)
    if (allocated(g2%object)) object_size = object_size + size(g2%object)
    allocate(object(object_size))
    object_size = 0
    if (allocated(g1%object)) then
       do i = 1, size(g1%object)
          if (.not. is_in_array(object(1:object_size),g1%object(i))) then
             object_size = object_size + 1
             object(object_size) = g1%object(i)
          end if
       end do
    end if
    if (allocated(g2%object)) then
       do i = 1, size(g2%object)
          if (.not. is_in_array(object(1:object_size),g2%object(i))) then
             object_size = object_size + 1
             object(object_size) = g2%object(i)
          end if
       end do
    end if

    !Finalise the groups
    call finalise(g1)
    call finalise(g2)

    !Re-initialise g1 with the new data
    call initialise(g1,type,atom(:atom_size),object(:object_size))   

  end subroutine merge_groups

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !% Deleting from a group
  !%
  !% NOTE: It is, for example, the atom with index 'i' (not 'this%atom(i)') which is deleted
  !%       Also, if atom/object is empty after deletion it is deallocated
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine group_delete_atom(this,i)

    type(Group), intent(inout)         :: this
    integer,     intent(in)            :: i
    integer                            :: pos
    integer, allocatable, dimension(:) :: temp_atom

    if (allocated(this%atom)) then
       if (size(this%atom) == 0) call system_abort('Group_Delete_Atom: Group has no atoms')
       !Find the atoms positions in the array
       pos = find_in_array(this%atom,i)
       if (pos == 0) then
          write(line,'(a,i0,a)')'Group_Delete_Atom: Atom ',i,' is not in this group'
          call system_abort(line)
       end if
       !Copy the others into the temporary
       this%atom(pos) = this%atom(size(this%atom))
       if (size(this%atom) > 1) then
          allocate(temp_atom(size(this%atom)-1))
          temp_atom = this%atom(1:size(this%atom)-1)
          call insertion_sort(temp_atom)
          !Copy back the remaining sorted data into the correctly sized atom array
          deallocate(this%atom)
          allocate(this%atom(size(temp_atom)))
          this%atom = temp_atom
          deallocate(temp_atom)
       else
          deallocate(this%atom)
       end if
    else
       call system_abort('Group_Delete_Atom: Group has no atoms')
    end if

  end subroutine group_delete_atom

  subroutine group_delete_object(this,i)

    type(Group), intent(inout)         :: this
    integer,     intent(in)            :: i
    integer                            :: pos
    integer, allocatable, dimension(:) :: temp_object

    if (allocated(this%object)) then
       if (size(this%object) == 0) call system_abort('Group_Delete_Object: Group has no objects')
       !Find the objects positions in the array
       pos = find_in_array(this%object,i)
       if (pos == 0) then
          write(line,'(a,i0,a)')'Group_Delete_Object: Object ',i,' is not in this group'
          call system_abort(line)
       end if
       !Copy the others into the temporary
       this%object(pos) = this%object(size(this%object))
       if (size(this%object) > 1) then
          allocate(temp_object(size(this%object)-1))
          temp_object = this%object(1:size(this%object)-1)
          call insertion_sort(temp_object)
          !Copy back the remaining sorted data into the correctly sized object array
          deallocate(this%object)
          allocate(this%object(size(temp_object)))
          this%object = temp_object
          deallocate(temp_object)
       else
          deallocate(this%object)
       end if
    else
       call system_abort('Group_Delete_Object: Group has no objects')
    end if

  end subroutine group_delete_object

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !% Adding objects to a group
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine group_add_atom(this,atom)

    type(Group),         intent(inout) :: this
    integer,             intent(in)    :: atom
    integer, allocatable, dimension(:) :: temp_atom
    
    if (allocated(this%atom)) then
       if (.not. is_in_array(this%atom,atom)) then
          allocate(temp_atom(size(this%atom)+1))
          temp_atom = (/this%atom,atom/)
          call insertion_sort(temp_atom)
          deallocate(this%atom)
          allocate(this%atom(size(temp_atom)))
          this%atom = temp_atom
          deallocate(temp_atom)
       end if
    else
       allocate(this%atom(1))
       this%atom(1) = atom
    end if

  end subroutine group_add_atom
  
  subroutine group_add_object(this,object)

    type(Group),         intent(inout) :: this
    integer,             intent(in)    :: object
    integer, allocatable, dimension(:) :: temp_object
    
    if (allocated(this%object)) then
       if (.not. is_in_array(this%object,object)) then
          allocate(temp_object(size(this%object)+1))
          temp_object = (/this%object,object/)
          call insertion_sort(temp_object)
          deallocate(this%object)
          allocate(this%object(size(temp_object)))
          this%object = temp_object
          deallocate(temp_object)
       end if
    else
       allocate(this%object(1))
       this%object(1) = object
    end if
    
  end subroutine group_add_object

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !% Inquiring about the size of a group
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  pure function group_n_atoms(this)
    
    type(Group), intent(in) :: this
    integer                 :: group_n_atoms

    if (allocated(this%atom)) then
       group_n_atoms = size(this%atom)
    else
       group_n_atoms = 0
    end if

  end function group_n_atoms

  function group_nth_atom(this,n)
    
    type(Group), intent(in) :: this
    integer,     intent(in) :: n
    integer                 :: group_nth_atom

    if (allocated(this%atom)) then
       if (n > size(this%atom)) then
          write(line,'(2(a,i0),a)')'Group_Nth_Atom: Requested atom (',n, &
                                   ') is greater than size of group (',size(this%atom),')'
          call system_abort(line)
       else
          group_nth_atom = this%atom(n)
       end if
    else
       call system_abort('Group_Nth_Atom: Group contains no atoms')
    end if

  end function group_nth_atom

  pure function group_n_objects(this)
    
    type(Group), intent(in) :: this
    integer                 :: group_n_objects

    if (allocated(this%object)) then
       group_n_objects = size(this%object)
    else
       group_n_objects = 0
    end if

  end function group_n_objects

  function group_nth_object(this,n)
    
    type(Group), intent(in) :: this
    integer,     intent(in) :: n
    integer                 :: group_nth_object

    if (allocated(this%object)) then
       if (n > size(this%object)) then
          write(line,'(2(a,i0),a)')'Group_Nth_Object: Requested object (',n, &
                          ') is greater than number of objects (',size(this%object),')'
          call system_abort(line)
       else
          Group_Nth_Object = this%object(n)
       end if
    else
       call system_abort('Group_Nth_Object: Group contains no objects')
    end if

  end function group_nth_object

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X I/O
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine group_print(this,verbosity,out,index)

    type(Group),       intent(in) :: this
    integer, intent(in), optional :: verbosity
    type(Inoutput), intent(inout), optional :: out
    integer, optional, intent(in) :: index
    integer                       :: i, nbloc, nrest
    logical                       :: have_atoms, have_objects
    character(10)                 :: form

    if (allocated(this%atom)) then
       if (size(this%atom) > 0) then
          have_atoms = .true.
       else
          have_atoms = .false.
       end if
    else
       have_atoms = .false.
    end if

    if (allocated(this%object)) then
       if (size(this%object) > 0) then
          have_objects = .true.
       else
          have_objects = .false.
       end if
    else
       have_objects = .false.
    end if

    call print('=========',verbosity,out)
    if (present(index)) then
       write(line,'(a,i0)')'Group ',index
       call print(line,verbosity,out)
    else
       call print('  Group  ',verbosity,out)
    end if
    call print('=========',verbosity,out)

    call print('',verbosity,out)

    write(line,'(a,i0)')' Type = ',this%type
    call print(line,verbosity,out)
    
    call print('Atoms:',verbosity,out)

    if (have_atoms) then
       nbloc=size(this%atom)/5
       nrest=mod(size(this%atom),5)
       do i = 1,nbloc
          write(line,'(5(1x,i0))') this%atom((i-1)*5+1:(i*5))
          call print(line,verbosity,out)
       end do
       if (nrest /= 0) then
          write(form,'(a,i0,a)') '(',nrest,'(1x,i0))'
          write(line,form) this%atom((i-1)*5+1:(i-1)*5+nrest)
          call print(line,verbosity,out)
       end if
    else
       call print('(none)',verbosity,out)
    end if

    call print('Objects:',verbosity,out)

    if (have_objects) then
       nbloc=size(this%object)/5
       nrest=mod(size(this%object),5)
       do i = 1,nbloc
          write(line,'(5(1x,i0))') this%object((i-1)*5+1:(i*5))
          call print(line,verbosity,out)
       end do
       if (nrest /= 0) then
          write(form,'(a,i0,a)') '(',nrest,'(1x,i0))'
          write(line,form) this%object((i-1)*5+1:(i-1)*5+nrest)
          call print(line,verbosity,out)
       end if
    else
       call print('(none)',verbosity,out)
    end if

    call print('',verbosity,out)

    call print('=========',verbosity,out)

  end subroutine group_print

  !
  !% Printing arrays of groups, with optional first group index
  !
  subroutine groups_print(this,verbosity,file,first)

    type(Group), dimension(:), intent(in) :: this
    integer, intent(in), optional :: verbosity
    type(Inoutput), intent(inout), optional :: file
    integer,     optional,     intent(in) :: first
    integer                               :: my_first, i

    my_first = 1 !default
    if (present(first)) my_first = first

    do i = 1, size(this)
       call print(this(i),verbosity,file,i-1+my_first)
       call print('', verbosity, file)
    end do

  end subroutine groups_print

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X Working with arrays of Group objects
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  
  !
  !% Create an atom lookup table: The supplied lookup table will be filled with the
  !% group number that atom i belongs to, with an error if it belongs to multiple groups
  !% or if an atom is ungrouped (unless ungrouped errors are turned off).
  !%
  !% The first group is 1 by default, but can be overridden with the 'first' argument
  !
  subroutine groups_create_lookup(this,lookup,first,ungrouped_error)

    type(Group), dimension(:), intent(in)    :: this
    integer,     dimension(:), intent(inout) :: lookup
    integer,     optional,     intent(in)    :: first
    logical,     optional,     intent(in)    :: ungrouped_error
    !local variables
    integer                                  :: i, g, n, my_first, ungrouped
    logical                                  :: ok, my_ungrouped_error
    
    my_ungrouped_error = .true.
    if (present(ungrouped_error)) my_ungrouped_error = ungrouped_error

    !Set the index of the first group correctly
    my_first = 1 !default
    if (present(first)) my_first = first

    !Initialise all lookups to a non-existent group number
    ungrouped = my_first - 1
    lookup = ungrouped

     do g = 1, size(this)
        do n = 1, Group_N_Atoms(this(g))
           i = Group_Nth_Atom(this(g),n)
           if (i > size(lookup)) then
              write(line,'(2(a,i0))')'Groups_Create_Lookup: Tried to store lookup for atom ',i, &
                                     ' but lookup table only has length ',size(lookup)
              call system_abort(line)
           end if
           if (lookup(i) == ungrouped) then
              lookup(i) = g - 1 + my_first
           else
              write(line,'(3(a,i0))')'Groups_Create_Lookup: Atom ',i,' is in groups ',lookup(i),' and ',(g - 1 + my_first)
              call system_abort(line)
           end if
        end do
     end do

     !Now check for ungrouped atoms
     if (my_ungrouped_error) then
        ok = .true.
        do i = 1, size(lookup)
           if (lookup(i) == ungrouped) then
              write(line,'(a,i0,a)')'Groups_Create_Lookup: Atom ',i,' is ungrouped'
              call print_warning(line)
              ok = .false.
           end if
        end do
        if (.not.ok) call system_abort('Refresh_Group_Lookups: Ungrouped atoms found')
     end if

  end subroutine groups_create_lookup
  
  !
  !% Tidy up / free memory used by allocatable group array, optionally leaving 'spare' groups empty. 
  !% Lookups will need to be refreshed after this.
  !
  subroutine tidy_groups(this,spare)

    type(Group), dimension(:), allocatable, intent(inout) :: this
    integer,     optional,                  intent(in)    :: spare
    type(Group), dimension(:), allocatable                :: temp_group
    integer                                               :: Ngroups, i, g

    !Count the number of groups in use (have atoms in them)
    Ngroups = 0
    do i = 1, size(this)
       if (allocated(this(i)%atom)) then
          if (size(this(i)%atom) > 0) Ngroups = Ngroups + 1
       end if
    end do

    !Allocate temporary space 
    allocate(temp_group(Ngroups))

    !Copy non-empty groups into temp space
    g = 0
    do i = 1, size(this)
       if (allocated(this(i)%atom)) then
          if (size(this(i)%atom) > 0) then
             !Copy the group
             g = g + 1
             if (allocated(this(i)%object)) then
                call initialise(temp_group(g),this(i)%type,this(i)%atom,this(i)%object)
             else
                call initialise(temp_group(g),this(i)%type,this(i)%atom)
             end if
          end if
       end if
    end do

    !Finalise and deallocate the input groups
    call finalise(this)

    !Reallocate to the correct size
    if (present(spare)) then
       allocate(this(Ngroups + spare))
    else
       allocate(this(Ngroups))
    end if
       
    !Copy over the data
    this(1:Ngroups) = temp_group !should work with or without 'spare'
    
    !Deallocate temporary space
    call finalise(temp_group)

  end subroutine tidy_groups

  !
  !% Given an array of groups, return the index of an empty one or the lowest index minus 1 (i.e. zero by default)
  !% if none are free. If the first index of the array is not 1, then supply it in 'first'
  !
  function free_group(this,first) result(n)

    type(Group), dimension(:), intent(in) :: this
    integer,     optional,     intent(in) :: first
    integer                               :: i, n, my_first
    logical                               :: atom_free, object_free, free_found
    
    my_first = 1
    if (present(first)) my_first = first
    free_found = .false.

    do i = 1, size(this)

       atom_free = .not.allocated(this(i)%atom)
       object_free = .not.allocated(this(i)%object)
       if (atom_free .and. object_free) then          
          n = i - 1 + my_first
          free_found = .true.
          exit
       end if

    end do

    if (.not.free_found) n = my_first - 1

  end function free_group

  !
  !% Count the number of free groups in an array
  !
  function num_free_groups(this) result(n)

    type(Group), dimension(:), intent(in) :: this
    integer                               :: i,n
    logical                               :: atom_free, object_free

    n = 0
    do i = 1, size(this)
       atom_free = .not.allocated(this(i)%atom)
       object_free = .not.allocated(this(i)%object)
       if (atom_free .and. object_free) n = n + 1
    end do

  end function num_free_groups

end module group_module
