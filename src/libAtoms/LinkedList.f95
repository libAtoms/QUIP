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
!X  Linked list module
!X  
!X 
!% The linked list module contains linked list implementations and retrieval routines.
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module linkedlist_module
  use error_module
  use system_module
  implicit none

  type LinkedList_i
     integer, allocatable :: data
     type(LinkedList_i), pointer :: next => null()
  endtype LinkedList_i

  type LinkedList_i1d
     integer, dimension(:), allocatable :: data
     type(LinkedList_i1d), pointer :: next => null()
  endtype LinkedList_i1d

  type LinkedList_i2d
     integer, dimension(:,:), allocatable :: data
     type(LinkedList_i2d), pointer :: next => null()
  endtype LinkedList_i2d

  interface initialise
     module procedure initialise_LinkedList_i, initialise_LinkedList_i1d, initialise_LinkedList_i2d
  endinterface initialise

  interface finalise
     module procedure finalise_LinkedList_i, finalise_LinkedList_i1d, finalise_LinkedList_i2d
  endinterface finalise

  interface insert
     module procedure insert_LinkedList_i, insert_LinkedList_i1d, insert_LinkedList_i2d
  endinterface insert

  interface delete_node
     module procedure delete_node_LinkedList_i, delete_node_LinkedList_i1d, delete_node_LinkedList_i2d
  endinterface delete_node

  interface update
     module procedure update_LinkedList_i, update_LinkedList_i1d, update_LinkedList_i2d
  endinterface update

  interface retrieve_node
     module procedure retrieve_node_LinkedList_i, retrieve_node_LinkedList_i1d, retrieve_node_LinkedList_i2d
  endinterface retrieve_node

  interface retrieve
     module procedure retrieve_LinkedList_i, retrieve_LinkedList_i1d, retrieve_LinkedList_i2d
  endinterface retrieve

  !interface next
  !   module procedure next_LinkedList_i1d, next_LinkedList_i2d
  !endinterface next

  interface last
     module procedure last_LinkedList_i, last_LinkedList_i1d, last_LinkedList_i2d
  endinterface last

  interface append
     module procedure append_LinkedList_i, append_LinkedList_i1d, append_LinkedList_i2d
  endinterface append

  interface size_LinkedList
     module procedure size_LinkedList_i, size_LinkedList_i1d, size_LinkedList_i2d
  endinterface size_LinkedList

  interface is_in_LinkedList
     module procedure is_in_LinkedList_i, is_in_LinkedList_i1d, is_in_LinkedList_i2d
  endinterface is_in_LinkedList

  contains

  ! integers

  subroutine initialise_LinkedList_i(this,data,error)

     type(LinkedList_i), pointer, intent(inout) :: this
     integer, intent(in), optional :: data
     integer, intent(out), optional :: error

     INIT_ERROR(error)

     if(associated(this)) call finalise(this)

     allocate(this)
     if( present(data) ) then
        this%data = data
     endif
     this%next => null()

  endsubroutine initialise_LinkedList_i

  subroutine finalise_LinkedList_i(this,error)
     type(LinkedList_i), pointer, intent(inout) :: this
     integer, intent(out), optional :: error

     type(LinkedList_i), pointer :: current, next

     INIT_ERROR(error)

     if( .not. associated(this) ) return

     current => this
     do while( associated(current) )
        next => current%next
        deallocate(current)
        current => next
     enddo
     this => null()

  endsubroutine finalise_LinkedList_i

  subroutine insert_LinkedList_i(this,data,error)
     type(LinkedList_i), pointer, intent(inout) :: this
     integer, intent(in), optional :: data
     integer, intent(out), optional :: error

     type(LinkedList_i), pointer :: current

     INIT_ERROR(error)

     if( .not. associated(this)) then
        call initialise(this,data,error)
     else
        allocate(current)
        if(present(data)) then
           current%data = data
        endif

        current%next => this%next
        this%next => current
     endif

  endsubroutine insert_LinkedList_i

  subroutine delete_node_LinkedList_i(this,node,error)
     type(LinkedList_i), pointer, intent(inout) :: this
     type(LinkedList_i), pointer, intent(inout) :: node
     integer, intent(out), optional :: error

     type(LinkedList_i), pointer :: delete_node => null(), previous => null()

     INIT_ERROR(error)

     if( .not. associated(this) ) return

     if( associated(this,node) ) then
        ! remove first node
        delete_node => this
        this => this%next
        deallocate(delete_node)
     else
        delete_node => this
        do
           if(associated(delete_node,node)) then
              previous%next => delete_node%next
              deallocate(delete_node)
              exit
           endif
           previous => delete_node
           delete_node => delete_node%next
           if( .not. associated(delete_node) ) exit
        enddo
     endif

  endsubroutine delete_node_LinkedList_i

  function is_in_LinkedList_i(this,data,error) result(found)
     type(LinkedList_i), pointer, intent(in) :: this
     integer, intent(in), optional :: data
     integer, intent(out), optional :: error

     logical :: found

     type(LinkedList_i), pointer :: current

     INIT_ERROR(error)

     found = .false.
     if( associated(this)) then
        current => this
        do 
           if( associated(current) ) then
              if( current%data == data ) then
                 found = .true.
                 exit
              endif
              current => current%next
           else
              exit
           endif
        enddo
     endif

  endfunction is_in_LinkedList_i

  subroutine update_LinkedList_i(this,data,error)
     type(LinkedList_i), pointer, intent(inout) :: this
     integer, intent(in) :: data
     integer, intent(out), optional :: error

     INIT_ERROR(error)
     
     if( .not. associated(this)) then
        RAISE_ERROR("update_LinkedList_i: linked list not initialised yet.",error)
     endif

     this%data = data

  endsubroutine update_LinkedList_i

  function retrieve_node_LinkedList_i(this,error) result(data)
     type(LinkedList_i), pointer, intent(in) :: this
     integer, intent(out), optional :: error

     integer, pointer :: data

     INIT_ERROR(error)
     
     if( .not. associated(this)) then
        RAISE_ERROR("retrieve_node_LinkedList_i1d: linked list not initialised yet.",error)
     endif

     data => this%data

  endfunction retrieve_node_LinkedList_i

  function next_LinkedList_i(this,error) result(next)
     type(LinkedList_i), pointer, intent(in) :: this
     integer, intent(out), optional :: error

     type(LinkedList_i), pointer :: next

     if( .not. associated(this)) then
        RAISE_ERROR("next_LinkedList_i1d: linked list not initialised yet.",error)
     endif

     next => this%next

  endfunction next_LinkedList_i

  function last_LinkedList_i(this,error) result(last)
     type(LinkedList_i), pointer, intent(in) :: this
     integer, intent(out), optional :: error

     type(LinkedList_i), pointer :: last

     if( .not. associated(this)) then
        RAISE_ERROR("next_LinkedList_i1d: linked list not initialised yet.",error)
     endif

     last => this
     do while( associated(last%next) )
        last => last%next
     enddo

  endfunction last_LinkedList_i

  subroutine append_LinkedList_i(this,data,error)
     type(LinkedList_i), pointer, intent(inout) :: this
     integer, intent(in), optional :: data
     integer, intent(out), optional :: error

     type(LinkedList_i), pointer :: current, new

     INIT_ERROR(error)

     if( .not. associated(this)) then
        call initialise(this,data,error)
     else
        allocate(new)
        if(present(data)) then
           new%data = data
        endif
        new%next => null()

        current => this
        do while( associated(current%next) )
           current => current%next
        enddo
        current%next => new
     endif

  endsubroutine append_LinkedList_i

  subroutine retrieve_LinkedList_i(this,data,error)
     type(LinkedList_i), pointer, intent(inout) :: this
     integer, dimension(:), allocatable, intent(out) :: data
     integer, intent(out), optional :: error

     type(LinkedList_i), pointer :: current
     integer :: i, d1

     INIT_ERROR(error)

     if( .not. associated(this)) then
        call reallocate(data,0)
        return
     endif

     d1 = size_LinkedList(this,error)

     call reallocate(data,d1)
     current => this
     i = 0
     do
        if( associated(current) ) then
           i = i + 1
           data(i) = current%data
           current => current%next
        else
           exit
        endif
     enddo

  endsubroutine retrieve_LinkedList_i

  function size_LinkedList_i(this,error) result(n)
     type(LinkedList_i), pointer, intent(in) :: this
     integer, intent(out), optional :: error

     type(LinkedList_i), pointer :: current

     integer :: n

     n = 0
     if( associated(this) ) then
        current => this
        do
           if( associated(current) ) then
              n = n + 1
              current => current%next
           else
              exit
           endif
        enddo
     endif

  endfunction size_LinkedList_i

  ! 1D integer arrays

  subroutine initialise_LinkedList_i1d(this,data,error)

     type(LinkedList_i1d), pointer, intent(inout) :: this
     integer, dimension(:), intent(in), optional :: data
     integer, intent(out), optional :: error

     INIT_ERROR(error)

     if(associated(this)) call finalise(this)

     allocate(this)
     if( present(data) ) then
        allocate( this%data(size(data,1)) )
        this%data = data
     endif
     this%next => null()

  endsubroutine initialise_LinkedList_i1d

  subroutine finalise_LinkedList_i1d(this,error)
     type(LinkedList_i1d), pointer, intent(inout) :: this
     integer, intent(out), optional :: error

     type(LinkedList_i1d), pointer :: current, next

     INIT_ERROR(error)

     if( .not. associated(this) ) return

     current => this
     do while( associated(current) )
        next => current%next
        if( allocated(current%data) ) deallocate(current%data)
        deallocate(current)
        current => next
     enddo
     this => null()

  endsubroutine finalise_LinkedList_i1d

  subroutine insert_LinkedList_i1d(this,data,error)
     type(LinkedList_i1d), pointer, intent(inout) :: this
     integer, dimension(:), intent(in), optional :: data
     integer, intent(out), optional :: error

     type(LinkedList_i1d), pointer :: current

     INIT_ERROR(error)

     if( .not. associated(this)) then
        call initialise(this,data,error)
     else
        allocate(current)
        if(present(data)) then
           allocate( current%data(size(data,1)) )
           current%data = data
        endif

        current%next => this%next
        this%next => current
     endif

  endsubroutine insert_LinkedList_i1d

  subroutine delete_node_LinkedList_i1d(this,node,error)
     type(LinkedList_i1d), pointer, intent(inout) :: this
     type(LinkedList_i1d), pointer, intent(inout) :: node
     integer, intent(out), optional :: error

     type(LinkedList_i1d), pointer :: delete_node => null(), previous => null()

     INIT_ERROR(error)

     if( .not. associated(this) ) return

     if( associated(this,node) ) then
        ! remove first node
        delete_node => this
        this => this%next
        if( allocated(delete_node%data) ) deallocate(delete_node%data)
        deallocate(delete_node)
     else
        delete_node => this
        do
           if(associated(delete_node,node)) then
              previous%next => delete_node%next
              if( allocated(delete_node%data) ) deallocate(delete_node%data)
              deallocate(delete_node)
              exit
           endif
           previous => delete_node
           delete_node => delete_node%next
           if( .not. associated(delete_node) ) exit
        enddo
     endif

  endsubroutine delete_node_LinkedList_i1d

  function is_in_LinkedList_i1d(this,data,error) result(found)
     type(LinkedList_i1d), pointer, intent(in) :: this
     integer, dimension(:), intent(in), optional :: data
     integer, intent(out), optional :: error

     logical :: found

     type(LinkedList_i1d), pointer :: current

     INIT_ERROR(error)

     found = .false.
     if( associated(this)) then
        current => this
        do 
           if( associated(current) ) then
              if( all( shape(current%data) == shape(data) ) ) then
                 if( all( current%data == data ) ) then
                    found = .true.
                    exit
                 endif
              endif
              current => current%next
           else
              exit
           endif
        enddo
     endif

  endfunction is_in_LinkedList_i1d

  subroutine update_LinkedList_i1d(this,data,error)
     type(LinkedList_i1d), pointer, intent(inout) :: this
     integer, dimension(:), intent(in) :: data
     integer, intent(out), optional :: error

     INIT_ERROR(error)
     
     if( .not. associated(this)) then
        RAISE_ERROR("update_LinkedList_i1d: linked list not initialised yet.",error)
     endif

     call reallocate( this%data, size(data,1) )
     this%data = data

  endsubroutine update_LinkedList_i1d

  function retrieve_node_LinkedList_i1d(this,error) result(data)
     type(LinkedList_i1d), pointer, intent(in) :: this
     integer, intent(out), optional :: error

     integer, dimension(:), pointer :: data

     INIT_ERROR(error)
     
     if( .not. associated(this)) then
        RAISE_ERROR("retrieve_node_LinkedList_i1d: linked list not initialised yet.",error)
     endif

     data => this%data

  endfunction retrieve_node_LinkedList_i1d

  function next_LinkedList_i1d(this,error) result(next)
     type(LinkedList_i1d), pointer, intent(in) :: this
     integer, intent(out), optional :: error

     type(LinkedList_i1d), pointer :: next

     if( .not. associated(this)) then
        RAISE_ERROR("next_LinkedList_i1d: linked list not initialised yet.",error)
     endif

     next => this%next

  endfunction next_LinkedList_i1d

  function last_LinkedList_i1d(this,error) result(last)
     type(LinkedList_i1d), pointer, intent(in) :: this
     integer, intent(out), optional :: error

     type(LinkedList_i1d), pointer :: last

     if( .not. associated(this)) then
        RAISE_ERROR("next_LinkedList_i1d: linked list not initialised yet.",error)
     endif

     last => this
     do while( associated(last%next) )
        last => last%next
     enddo

  endfunction last_LinkedList_i1d

  subroutine append_LinkedList_i1d(this,data,error)
     type(LinkedList_i1d), pointer, intent(inout) :: this
     integer, dimension(:), intent(in), optional :: data
     integer, intent(out), optional :: error

     type(LinkedList_i1d), pointer :: current, new

     INIT_ERROR(error)

     if( .not. associated(this)) then
        call initialise(this,data,error)
     else
        allocate(new)
        if(present(data)) then
           allocate( new%data(size(data,1)) )
           new%data = data
        endif
        new%next => null()

        current => this
        do while( associated(current%next) )
           current => current%next
        enddo
        current%next => new
     endif

  endsubroutine append_LinkedList_i1d

  subroutine retrieve_LinkedList_i1d(this,data,error)
     type(LinkedList_i1d), pointer, intent(inout) :: this
     integer, dimension(:,:), allocatable, intent(out) :: data
     integer, intent(out), optional :: error

     type(LinkedList_i1d), pointer :: current
     integer :: i, d1, d2

     INIT_ERROR(error)

     if( .not. associated(this)) then
        call reallocate(data,0,0)
        return
     endif

     d2 = size_LinkedList(this,error)
     if( allocated(this%data) ) then
        d1 = size(this%data,1)
     else
        RAISE_ERROR("retrieve_LinkedList_i1d: no data content in head node.", error)
     endif

     call reallocate(data,d1,d2)
     current => this
     i = 0
     do
        if( associated(current) ) then
           i = i + 1
           if( .not. allocated(current%data) ) then
              RAISE_ERROR("retrieve_LinkedList_i1d: data missing from node "//i, error)
           endif
           if( any( (/d1/) /= shape(current%data) ) ) then
              RAISE_ERROR("retrieve_LinkedList_i1d: data array inconsistency in node "//i, error)
           endif
           data(:,i) = current%data
           current => current%next
        else
           exit
        endif
     enddo

  endsubroutine retrieve_LinkedList_i1d

  function size_LinkedList_i1d(this,error) result(n)
     type(LinkedList_i1d), pointer, intent(in) :: this
     integer, intent(out), optional :: error

     type(LinkedList_i1d), pointer :: current

     integer :: n

     n = 0
     if( associated(this) ) then
        current => this
        do
           if( associated(current) ) then
              n = n + 1
              current => current%next
           else
              exit
           endif
        enddo
     endif

  endfunction size_LinkedList_i1d

  ! 2D integer arrays

  subroutine initialise_LinkedList_i2d(this,data,error)

     type(LinkedList_i2d), pointer, intent(inout) :: this
     integer, dimension(:,:), intent(in), optional :: data
     integer, intent(out), optional :: error

     INIT_ERROR(error)

     if(associated(this)) call finalise(this)

     allocate(this)
     if( present(data) ) then
        allocate( this%data(size(data,1),size(data,2)) )
        this%data = data
     endif
     this%next => null()

  endsubroutine initialise_LinkedList_i2d

  subroutine finalise_LinkedList_i2d(this,error)
     type(LinkedList_i2d), pointer, intent(inout) :: this
     integer, intent(out), optional :: error

     type(LinkedList_i2d), pointer :: current, next

     INIT_ERROR(error)

     if( .not. associated(this) ) return

     current => this
     do while( associated(current) )
        next => current%next
        if( allocated(current%data) ) deallocate(current%data)
        deallocate(current)
        current => next
     enddo

     this => null()

  endsubroutine finalise_LinkedList_i2d

  subroutine insert_LinkedList_i2d(this,data,error)
     type(LinkedList_i2d), pointer, intent(inout) :: this
     integer, dimension(:,:), intent(in), optional :: data
     integer, intent(out), optional :: error

     type(LinkedList_i2d), pointer :: current

     INIT_ERROR(error)

     if( .not. associated(this)) then
        call initialise(this,data,error)
     else
        allocate(current)
        if(present(data)) then
           allocate( current%data(size(data,1), size(data,2)) )
           current%data = data
        endif

        current%next => this%next
        this%next => current
     endif

  endsubroutine insert_LinkedList_i2d

  subroutine delete_node_LinkedList_i2d(this,node,error)
     type(LinkedList_i2d), pointer, intent(inout) :: this
     type(LinkedList_i2d), pointer, intent(inout) :: node
     integer, intent(out), optional :: error

     type(LinkedList_i2d), pointer :: delete_node => null(), previous => null()

     INIT_ERROR(error)

     if( .not. associated(this) ) return

     if( associated(this,node) ) then
        ! remove first node
        delete_node => this
        this => this%next
        if( allocated(delete_node%data) ) deallocate(delete_node%data)
        deallocate(delete_node)
     else
        delete_node => this
        do
           if(associated(delete_node,node)) then
              previous%next => delete_node%next
              if( allocated(delete_node%data) ) deallocate(delete_node%data)
              deallocate(delete_node)
              exit
           endif
           previous => delete_node
           delete_node => delete_node%next
           if( .not. associated(delete_node) ) exit
        enddo
     endif

  endsubroutine delete_node_LinkedList_i2d

  function is_in_LinkedList_i2d(this,data,error) result(found)
     type(LinkedList_i2d), pointer, intent(in) :: this
     integer, dimension(:,:), intent(in), optional :: data
     integer, intent(out), optional :: error

     logical :: found

     type(LinkedList_i2d), pointer :: current

     INIT_ERROR(error)

     found = .false.
     if( associated(this)) then
        current => this
        do 
           if( associated(current) ) then
              if( all( shape(current%data) == shape(data) ) ) then
                 if( all( current%data == data ) ) then
                    found = .true.
                    exit
                 endif
              endif
              current => current%next
           else
              exit
           endif
        enddo
     endif

  endfunction is_in_LinkedList_i2d

  subroutine update_LinkedList_i2d(this,data,error)
     type(LinkedList_i2d), pointer, intent(inout) :: this
     integer, dimension(:,:), intent(in) :: data
     integer, intent(out), optional :: error

     INIT_ERROR(error)
     
     if( .not. associated(this)) then
        RAISE_ERROR("update_LinkedList_i2d: linked list not initialised yet.",error)
     endif

     call reallocate( this%data, size(data,1), size(data,2) )
     this%data = data

  endsubroutine update_LinkedList_i2d

  function retrieve_node_LinkedList_i2d(this,error) result(data)
     type(LinkedList_i2d), pointer, intent(in) :: this
     integer, intent(out), optional :: error

     integer, dimension(:,:), pointer :: data

     INIT_ERROR(error)
     
     if( .not. associated(this)) then
        RAISE_ERROR("retrieve_node_LinkedList_i2d: linked list not initialised yet.",error)
     endif

     data => this%data

  endfunction retrieve_node_LinkedList_i2d

  function next_LinkedList_i2d(this,error) result(next)
     type(LinkedList_i2d), pointer, intent(in) :: this
     integer, intent(out), optional :: error

     type(LinkedList_i2d), pointer :: next

     if( .not. associated(this)) then
        RAISE_ERROR("next_LinkedList_i2d: linked list not initialised yet.",error)
     endif

     next => this%next

  endfunction next_LinkedList_i2d

  function last_LinkedList_i2d(this,error) result(last)
     type(LinkedList_i2d), pointer, intent(in) :: this
     integer, intent(out), optional :: error

     type(LinkedList_i2d), pointer :: last

     if( .not. associated(this)) then
        RAISE_ERROR("next_LinkedList_i1d: linked list not initialised yet.",error)
     endif

     last => this
     do while( associated(last%next) )
        last => last%next
     enddo
  endfunction last_LinkedList_i2d

  subroutine append_LinkedList_i2d(this,data,error)
     type(LinkedList_i2d), pointer, intent(inout) :: this
     integer, dimension(:,:), intent(in), optional :: data
     integer, intent(out), optional :: error

     type(LinkedList_i2d), pointer :: current, new

     INIT_ERROR(error)

     if( .not. associated(this)) then
        call initialise(this,data,error)
     else
        allocate(new)
        if(present(data)) then
           allocate( new%data(size(data,1), size(data,2)) )
           new%data = data
        endif
        new%next => null()

        current => this
        do while( associated(current%next) )
           current => current%next
        enddo
        current%next => new
     endif

  endsubroutine append_LinkedList_i2d

  subroutine retrieve_LinkedList_i2d(this,data,error)
     type(LinkedList_i2d), pointer, intent(inout) :: this
     integer, dimension(:,:,:), allocatable, intent(out) :: data
     integer, intent(out), optional :: error

     type(LinkedList_i2d), pointer :: current
     integer :: i, d1, d2, d3

     INIT_ERROR(error)

     if( .not. associated(this)) then
        call reallocate(data,0,0,0)
        return
     endif

     d3 = size_LinkedList(this,error)
     if( allocated(this%data) ) then
        d1 = size(this%data,1)
        d2 = size(this%data,2)
     else
        RAISE_ERROR("retrieve_LinkedList_i2d: no data content in head node.", error)
     endif

     call reallocate(data,d1,d2,d3)
     current => this
     i = 0
     do
        if( associated(current) ) then
           i = i + 1
           if( .not. allocated(current%data) ) then
              RAISE_ERROR("retrieve_LinkedList_i2d: data missing from node "//i, error)
           endif
           if( any( (/d1,d2/) /= shape(current%data) ) ) then
              RAISE_ERROR("retrieve_LinkedList_i2d: data array inconsistency in node "//i, error)
           endif
           data(:,:,i) = current%data
           current => current%next
        else
           exit
        endif
     enddo

  endsubroutine retrieve_LinkedList_i2d

  function size_LinkedList_i2d(this,error) result(n)
     type(LinkedList_i2d), pointer, intent(in) :: this
     integer, intent(out), optional :: error

     type(LinkedList_i2d), pointer :: current

     integer :: n

     n = 0
     if( associated(this) ) then
        current => this
        do
           if( associated(current) ) then
              n = n + 1
              current => current%next
           else
              exit
           endif
        enddo
     endif

  endfunction size_LinkedList_i2d

endmodule linkedlist_module
