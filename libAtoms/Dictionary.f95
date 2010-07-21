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

!% A Dictionary object contains a list of keys (strings) and corresponding entries.
!% Entries are a variable type, containing one of:
!% 'integer', 'real(dp)', 'complex(dp)', 'logical', 'character(len=value_len)' 
!% or a 1-D array of any of those. 2-D arrays of integers and reals are also supported.

module dictionary_module

  use system_module
  use linearalgebra_module
  use mpi_context_module
  implicit none

  integer, parameter :: &
       T_NONE = 0, &
       T_INTEGER = 1, T_REAL =2, T_COMPLEX = 3, T_LOGICAL=4, &
       T_INTEGER_A = 5, T_REAL_A = 6, T_COMPLEX_A = 7, T_LOGICAL_A = 8, &
       T_CHAR = 9, T_CHAR_A = 10, T_DATA = 11, T_INTEGER_A2 = 12, T_REAL_A2 = 13 !% OMIT


  integer, parameter :: PROPERTY_INT     = 1 !% Property types, used by Atoms and Table
  integer, parameter :: PROPERTY_REAL    = 2
  integer, parameter :: PROPERTY_STR     = 3
  integer, parameter :: PROPERTY_LOGICAL = 4

  integer, parameter :: key_len = 256 !% maximum length of key string, usually not checked for safety

  integer, parameter :: value_len = 1024 !% maximum length of value string


  integer, parameter :: dict_field_length = 1023  !% Maximum field width during parsing
  integer, parameter :: dict_n_fields = 100       !% Maximum number of fields during parsing



  type DictData
     integer, dimension(:), allocatable :: d
  end type DictData

  type DictEntry 
     !% OMIT
     integer :: type = T_NONE
     integer :: len = 1
     integer :: len2(2) = (/0, 0/)

     logical :: own_data = .true. !% True if we own the data and should free it in finalise()

     integer i
     real(dp) r
     complex(dp) c
     logical l
     character(len=value_len) s

     integer, pointer :: i_a(:) => null()
     real(dp), pointer :: r_a(:) => null()
     complex(dp), pointer :: c_a(:) => null()
     logical, pointer :: l_a(:) => null()
     character(len=value_len), pointer :: s_a(:) => null()

     integer, pointer :: i_a2(:,:) => null()
     real(dp), pointer :: r_a2(:,:) => null()

     type(DictData) :: d

  end type DictEntry

  integer, parameter :: n_entry_block = 10 !% OMIT

  type Dictionary
     integer :: N !% number of entries in use
     character(len=key_len), allocatable :: keys(:) !% array of keys
     type(DictEntry), allocatable :: entries(:)    !% array of entries
  end type Dictionary

  !% Initialise a new empty dictionary
  private :: dictionary_initialise
  interface initialise
     module procedure dictionary_initialise
  end interface initialise

  !% Finalise dictionary
  private :: dictionary_finalise, dictentry_finalise
  interface finalise
     module procedure dictionary_finalise, dictentry_finalise
  end interface finalise

  !% Print a DictEntry or a Dictionary
  private :: dictentry_print, dictionary_print
  interface print
     module procedure dictentry_print, dictionary_print
  end interface print

  !% Set a value in a Dictionary
  private :: dictionary_set_value_i, dictionary_set_value_r, dictionary_set_value_c, dictionary_set_value_l
  private :: dictionary_set_value_i_a, dictionary_set_value_r_a, dictionary_set_value_c_a, dictionary_set_value_l_a
  private :: dictionary_set_value_s, dictionary_set_value_s_a
  private :: dictionary_set_value_d
  interface set_value
     module procedure dictionary_set_value_i, dictionary_set_value_r, dictionary_set_value_c, dictionary_set_value_l
     module procedure dictionary_set_value_i_a, dictionary_set_value_r_a, dictionary_set_value_c_a, dictionary_set_value_l_a
     module procedure dictionary_set_value_s, dictionary_set_value_s_a
     module procedure dictionary_set_value_d
     module procedure dictionary_set_value_i_a2
     module procedure dictionary_set_value_r_a2
  end interface set_value

  !% Get a value from a Dictionary
  private :: dictionary_get_value_i, dictionary_get_value_r, dictionary_get_value_c, dictionary_get_value_l
  private :: dictionary_get_value_i_a, dictionary_get_value_r_a, dictionary_get_value_c_a, dictionary_get_value_l_a
  private :: dictionary_get_value_s, dictionary_get_value_s_a
  private :: dictionary_get_value_d
  interface get_value
     module procedure dictionary_get_value_i, dictionary_get_value_r, dictionary_get_value_c, dictionary_get_value_l
     module procedure dictionary_get_value_i_a, dictionary_get_value_r_a, dictionary_get_value_c_a, dictionary_get_value_l_a
     module procedure dictionary_get_value_s, dictionary_get_value_s_a
     module procedure dictionary_get_value_d
     module procedure dictionary_get_value_i_a2
     module procedure dictionary_get_value_r_a2
  end interface get_value

  private :: dictionary_assign_pointer_i, dictionary_assign_pointer_r
  private :: dictionary_assign_pointer_c, dictionary_assign_pointer_l, dictionary_assign_pointer_s
  private :: dictionary_assign_pointer_i2, dictionary_assign_pointer_r2
  interface assign_pointer
     module procedure dictionary_assign_pointer_i
     module procedure dictionary_assign_pointer_r
     module procedure dictionary_assign_pointer_c
     module procedure dictionary_assign_pointer_l
     module procedure dictionary_assign_pointer_s
     module procedure dictionary_assign_pointer_i2
     module procedure dictionary_assign_pointer_r2
  end interface assign_pointer

  private :: dictionary_set_value_pointer_i, dictionary_set_value_pointer_r
  private :: dictionary_set_value_pointer_c, dictionary_set_value_pointer_l, dictionary_set_value_pointer_s
  private :: dictionary_set_value_pointer_i2, dictionary_set_value_pointer_r2
  interface set_value_pointer
     module procedure dictionary_set_value_pointer_i
     module procedure dictionary_set_value_pointer_r
     module procedure dictionary_set_value_pointer_c
     module procedure dictionary_set_value_pointer_l
     module procedure dictionary_set_value_pointer_s
     module procedure dictionary_set_value_pointer_i2
     module procedure dictionary_set_value_pointer_r2
  end interface set_value_pointer

  private :: dictionary_add_array_i, dictionary_add_array_r
  private :: dictionary_add_array_c, dictionary_add_array_l, dictionary_add_array_s
  private :: dictionary_add_array_i2, dictionary_add_array_r2
  interface add_array
     module procedure dictionary_add_array_i
     module procedure dictionary_add_array_r
     module procedure dictionary_add_array_c
     module procedure dictionary_add_array_l
     module procedure dictionary_add_array_s
     module procedure dictionary_add_array_i2
     module procedure dictionary_add_array_r2
  end interface add_array

  !% Remove an entry from a Dictionary
  private :: dictionary_remove_value
  interface remove_value
     module procedure dictionary_remove_value
  end interface remove_value

  !% Write a string representation of this dictionary
  interface write_string
     module procedure dictionary_write_string
  end interface write_string

  !% Read into this dictionary from a string
  interface read_string
     module procedure dictionary_read_string
  end interface read_string

  interface subset
     module procedure dictionary_subset
  end interface subset

  interface swap
     module procedure dictionary_swap
  end interface swap

  interface has_key
     module procedure dictionary_has_key
  end interface has_key

  interface bcast
     module procedure dictionary_bcast
  end interface bcast

  interface assignment(=)
     module procedure dictionary_deepcopy
  end interface assignment(=)

  private add_entry, extend_entries, remove_entry

contains

  ! ****************************************************************************
  ! *
  ! *    DictEntry methods
  ! * 
  ! ****************************************************************************

  subroutine dictentry_print(this, key, verbosity, file)
    type(DictEntry), intent(in) :: this
    character(len=*), intent(in) :: key
    integer, intent(in), optional :: verbosity
    type(Inoutput), intent(inout), optional :: file
    integer :: i

    if (this%type == T_NONE) then
       call print("Dict entry NONE", verbosity, file)
    else if (this%type == T_INTEGER) then
       write (line,'("Dict entry integer ", A,1x,I0)') trim(key), this%i
       call print(line, verbosity,file)
    else if (this%type == T_REAL) then
       write (line, '("Dict entry real ", A,1x,F30.20)') trim(key), this%r
       call print(line, verbosity,file)
    else if (this%type == T_COMPLEX) then
       write (line,'("Dict entry complex ", A,1x,2F30.20)') trim(key), this%c
       call print(line, verbosity,file)
    else if (this%type == T_LOGICAL) then
       write (line, '("Dict entry logical ", A,1x,L2)') trim(key), this%l
       call print(line, verbosity,file)
    else if (this%type == T_CHAR) then
       write (line,'("Dict entry string ", A,1x,A)') trim(key), this%s
       call print(line, verbosity,file)
    else if (this%type == T_INTEGER_A) then
       write (line, '("Dict entry integer array ",A,1x,I0)') trim(key), size(this%i_a)
       call print(line, verbosity,file)
       call print(this%i_a, verbosity,file)
    else if (this%type == T_REAL_A) then
       write (line, '("Dict entry real array ",A,1x,I0)') trim(key), size(this%r_a)
       call print(line, verbosity,file)
       call print(this%r_a, verbosity,file)
    else if (this%type == T_COMPLEX_A) then
       write (line,'("Dict entry complex array ",A,1x,I0)') trim(key), size(this%c_a)
       call print(line, verbosity,file)
       call print(this%c_a,verbosity,file)
    else if (this%type == T_LOGICAL_A) then
       write (line,'("Dict entry logical array ",A,1x,I0)') trim(key), size(this%l_a)
       call print(line, verbosity,file)
       do i=1,size(this%l_a)
          call print(this%l_a(i),verbosity,file)
       end do
    else if (this%type == T_CHAR_A) then
       write (line,'("Dict entry string array ",A,1x,I0)') trim(key), size(this%s_a)
       call print(line, verbosity,file)
       do i=1,size(this%s_a)
          call print(this%s_a(i),verbosity,file)
       end do
    else if (this%type == T_INTEGER_A2) then
       call print("Dict entry integer array "//shape(this%i_a2), verbosity,file)
       call print(this%i_a2, verbosity,file)
    else if (this%type == T_REAL_A2) then
       call print("Dict entry real array "//shape(this%r_a2), verbosity,file)
       call print(this%r_a2, verbosity,file)
    else if (this%type == T_DATA) then
       print '("Dict entry arbitary data ",A,1x,I0)', trim(key), size(this%d%d)
       call print(line, verbosity,file)
    endif
  end subroutine dictentry_print

  subroutine dictionary_print(this, verbosity, file)
    type(Dictionary), intent(in) :: this
    integer, intent(in), optional :: verbosity
    type(Inoutput), intent(inout), optional :: file

    call print("Dictionary, allocated "// size(this%entries) //" n_entries " // this%N, verbosity=verbosity,file=file)
    call print('', verbosity=verbosity,file=file)
#ifdef NO_F2003_NEW_LINE
    call print(write_string(this,entry_sep=char(13)),verbosity=verbosity,file=file)
#else
    call print(write_string(this,entry_sep=new_line(' ')),verbosity=verbosity,file=file)
#endif

  end subroutine dictionary_print


  ! ****************************************************************************
  ! *
  ! *    Dictionary methods
  ! * 
  ! ****************************************************************************

  subroutine dictionary_initialise(this)
    type(Dictionary), intent(inout) :: this

    call finalise(this)

    call extend_entries(this, n_entry_block)
  end subroutine dictionary_initialise

  subroutine dictionary_finalise(this)
    type(Dictionary), intent(inout) :: this

    integer :: i

    ! Don't trust ifort to dellocate properly
    if (allocated(this%entries)) then
       ! only finalise entries up to this%n, not size(entries):
       ! entries beyond this have been used in the past if some have now been removed,
       ! but they won't point to any different data since remove_entry() does a shallow copy
       do i=1,this%n
          call finalise(this%entries(i))
       end do
       deallocate(this%entries)
    end if
    if (allocated(this%keys)) deallocate(this%keys)
    this%N = 0

  end subroutine dictionary_finalise

  subroutine dictentry_finalise(this)
    type(DictEntry), intent(inout) :: this

    if (allocated(this%d%d)) deallocate(this%d%d)

    if (.not. this%own_data) return

    if (associated(this%i_a)) deallocate(this%i_a)
    this%i_a => null()

    if (associated(this%r_a)) deallocate(this%r_a)
    this%r_a => null()

    if (associated(this%c_a)) deallocate(this%c_a)
    this%c_a => null()

    if (associated(this%l_a)) deallocate(this%l_a)
    this%l_a => null()

    if (associated(this%s_a)) deallocate(this%s_a)
    this%s_a => null()

    if (associated(this%i_a2)) deallocate(this%i_a2)
    this%i_a2 => null()

    if (associated(this%r_a2)) deallocate(this%r_a2)
    this%r_a2 => null()

  end subroutine dictentry_finalise

  subroutine dictionary_get_type_and_size(this, key, type, thesize, thesize2)
    type(Dictionary), intent(in) :: this
    character(len=*), intent(in) :: key
    integer, intent(out) :: type, thesize, thesize2(2)
    integer :: entry_i

    entry_i = lookup_entry_i(this, key)

    if (entry_i <= 0) type = -1
    type = this%entries(entry_i)%type
    thesize = this%entries(entry_i)%len
    thesize2 = this%entries(entry_i)%len2

  end subroutine dictionary_get_type_and_size

  ! ****************************************************************************
  ! *
  ! *    set_value() interface
  ! * 
  ! ****************************************************************************

  subroutine dictionary_set_value_i(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value

    type(DictEntry) entry
    integer entry_i

    entry%type = T_INTEGER
    entry%i = value
    entry_i = add_entry(this, key, entry)
    call finalise(entry)
  end subroutine dictionary_set_value_i

  subroutine dictionary_set_value_r(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: value

    type(DictEntry) entry
    integer entry_i

    entry%type = T_REAL
    entry%r = value
    entry_i = add_entry(this, key, entry)
    call finalise(entry)
  end subroutine dictionary_set_value_r

  subroutine dictionary_set_value_c(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    complex(dp), intent(in) :: value

    type(DictEntry) entry
    integer entry_i

    entry%type = T_COMPLEX
    entry%c = value
    entry_i = add_entry(this, key, entry)
    call finalise(entry)
  end subroutine dictionary_set_value_c

  subroutine dictionary_set_value_l(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in) :: value

    type(DictEntry) entry
    integer entry_i

    entry%type = T_LOGICAL
    entry%l = value
    entry_i = add_entry(this, key, entry)
    call finalise(entry)
  end subroutine dictionary_set_value_l

  subroutine dictionary_set_value_s(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value

    type(DictEntry) entry
    integer entry_i

    entry%type = T_CHAR
    entry%s = value
    entry_i = add_entry(this, key, entry)
    call finalise(entry)
  end subroutine dictionary_set_value_s


  subroutine dictionary_set_value_i_a(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value(:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_INTEGER_A
    entry%len = size(value)
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%i_a(size(value)))
    this%entries(entry_i)%i_a = value
    call finalise(entry)
  end subroutine dictionary_set_value_i_a

  subroutine dictionary_set_value_r_a(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: value(:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_REAL_A
    entry%len = size(value)
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%r_a(size(value)))
    this%entries(entry_i)%r_a = value
    call finalise(entry)
  end subroutine dictionary_set_value_r_a

  subroutine dictionary_set_value_i_a2(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value(:,:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_INTEGER_A2
    entry%len2 = shape(value)
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%i_a2(size(value,1),size(value,2)))
    this%entries(entry_i)%i_a2 = value
    call finalise(entry)
  end subroutine dictionary_set_value_i_a2

  subroutine dictionary_set_value_r_a2(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: value(:,:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_REAL_A2
    entry%len2 = shape(value)
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%r_a2(size(value,1),size(value,2)))
    this%entries(entry_i)%r_a2 = value
    call finalise(entry)
  end subroutine dictionary_set_value_r_a2

  subroutine dictionary_set_value_c_a(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    complex(dp), intent(in) :: value(:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_COMPLEX_A
    entry%len = size(value)
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%c_a(size(value)))
    this%entries(entry_i)%c_a = value
    call finalise(entry)
  end subroutine dictionary_set_value_c_a

  subroutine dictionary_set_value_l_a(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in) :: value(:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_LOGICAL_A
    entry%len = size(value)
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%l_a(size(value)))
    this%entries(entry_i)%l_a = value
    call finalise(entry)
  end subroutine dictionary_set_value_l_a

  subroutine dictionary_set_value_s_a(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value(:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_CHAR_A
    entry%len = size(value)
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%s_a(size(value)))
    this%entries(entry_i)%s_a = value
    call finalise(entry)
  end subroutine dictionary_set_value_s_a

  subroutine dictionary_set_value_d(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    type(DictData), intent(in) :: value

    type(DictEntry) entry
    integer entry_i

    entry%type = T_DATA
    entry%d = value
    entry_i = add_entry(this, key, entry)
    call finalise(entry)
  end subroutine dictionary_set_value_d


  ! ****************************************************************************
  ! *
  ! *    get_value() interface
  ! * 
  ! ****************************************************************************

  function dictionary_get_value_i(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    integer, intent(out) :: v
    logical :: dictionary_get_value_i
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_i = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_INTEGER) then
       v = this%entries(entry_i)%i
       dictionary_get_value_i = .true.
    else
       dictionary_get_value_i = .false.
    endif
  end function dictionary_get_value_i

  function dictionary_get_value_r(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    real(dp), intent(out) :: v
    logical :: dictionary_get_value_r
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_r = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_REAL) then
       v = this%entries(entry_i)%r
       dictionary_get_value_r = .true.
    else
       dictionary_get_value_r = .false.
    endif
  end function dictionary_get_value_r


  function dictionary_get_value_c(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    complex(dp), intent(out) :: v
    logical :: dictionary_get_value_c
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_c = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_COMPLEX) then
       v = this%entries(entry_i)%c
       dictionary_get_value_c = .true.
    else
       dictionary_get_value_c = .false.
    endif
  end function dictionary_get_value_c

  function dictionary_get_value_l(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    logical, intent(out) :: v
    logical :: dictionary_get_value_l
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_l = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_LOGICAL) then
       v = this%entries(entry_i)%l
       dictionary_get_value_l = .true.
    else
       dictionary_get_value_l = .false.
    endif
  end function dictionary_get_value_l


  function dictionary_get_value_s(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    character(len=value_len), intent(out) :: v
    logical :: dictionary_get_value_s
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_s = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_CHAR) then
       v = this%entries(entry_i)%s
       dictionary_get_value_s = .true.
    else
       dictionary_get_value_s = .false.
    endif
  end function dictionary_get_value_s



  function dictionary_get_value_i_a(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    integer, intent(out) :: v(:)
    logical :: dictionary_get_value_i_a
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_i_a = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_INTEGER_A) then
       if (size(v) >= size(this%entries(entry_i)%i_a)) then
          v(1:size(this%entries(entry_i)%i_a)) = this%entries(entry_i)%i_a
          dictionary_get_value_i_a = .true.
       else
          dictionary_get_value_i_a = .false.
       endif
    else
       dictionary_get_value_i_a = .false.
    endif
  end function dictionary_get_value_i_a

  function dictionary_get_value_r_a(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    real(dp), intent(out) :: v(:)
    logical :: dictionary_get_value_r_a
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_r_a = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_REAL_A) then
       if (size(v) >= size(this%entries(entry_i)%r_a)) then
          v(1:size(this%entries(entry_i)%r_a)) = this%entries(entry_i)%r_a
          dictionary_get_value_r_a = .true.
       else
          dictionary_get_value_r_a = .false.
       endif
    else
       dictionary_get_value_r_a = .false.
    endif
  end function dictionary_get_value_r_a

  function dictionary_get_value_i_a2(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    integer, intent(out) :: v(:,:)
    logical :: dictionary_get_value_i_a2
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_i_a2 = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_INTEGER_A2) then
       if (size(v,1) >= size(this%entries(entry_i)%i_a2,1) .and. &
            size(v,2) >= size(this%entries(entry_i)%i_a2,2)) then
          v(1:size(this%entries(entry_i)%i_a2,1),1:size(this%entries(entry_i)%i_a2,2)) = this%entries(entry_i)%i_a2
          dictionary_get_value_i_a2 = .true.
       else
          dictionary_get_value_i_a2 = .false.
       endif
    else
       dictionary_get_value_i_a2 = .false.
    endif
  end function dictionary_get_value_i_a2

  function dictionary_get_value_r_a2(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    real(dp), intent(out) :: v(:,:)
    logical :: dictionary_get_value_r_a2
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_r_a2 = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_REAL_A2) then
       if (size(v,1) >= size(this%entries(entry_i)%r_a2,1) .and. &
            size(v,2) >= size(this%entries(entry_i)%r_a2,2) ) then
          v(1:size(this%entries(entry_i)%r_a2,1),1:size(this%entries(entry_i)%r_a2,2)) = this%entries(entry_i)%r_a2
          dictionary_get_value_r_a2 = .true.
       else
          dictionary_get_value_r_a2 = .false.
       endif
    else
       dictionary_get_value_r_a2 = .false.
    endif
  end function dictionary_get_value_r_a2


  function dictionary_get_value_c_a(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    complex(dp), intent(out) :: v(:)
    logical :: dictionary_get_value_c_a
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_c_a = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_COMPLEX_A) then
       if (size(v) >= size(this%entries(entry_i)%c_a)) then
          v(1:size(this%entries(entry_i)%c_a)) = this%entries(entry_i)%c_a
          dictionary_get_value_c_a = .true.
       else
          dictionary_get_value_c_a = .false.
       endif
    else
       dictionary_get_value_c_a = .false.
    endif
  end function dictionary_get_value_c_a

  function dictionary_get_value_l_a(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    logical, intent(out) :: v(:)
    logical :: dictionary_get_value_l_a
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_l_a = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_LOGICAL_A) then
       if (size(v) >= size(this%entries(entry_i)%l_a)) then
          v(1:size(this%entries(entry_i)%l_a)) = this%entries(entry_i)%l_a
          dictionary_get_value_l_a = .true.
       else
          dictionary_get_value_l_a = .false.
       endif
    else
       dictionary_get_value_l_a = .false.
    endif
  end function dictionary_get_value_l_a


  function dictionary_get_value_s_a(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    character(len=value_len), intent(out) :: v(:)
    logical :: dictionary_get_value_s_a
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_s_a = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_CHAR_A) then
       if (size(v) >= size(this%entries(entry_i)%s_a)) then
          v(1:size(this%entries(entry_i)%s_a)) = this%entries(entry_i)%s_a
          dictionary_get_value_s_a = .true.
       else
          dictionary_get_value_s_a = .false.
       endif
    else
       dictionary_get_value_s_a = .false.
    endif
  end function dictionary_get_value_s_a



  function dictionary_get_value_d(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    type(DictData), intent(out) :: v
    logical :: dictionary_get_value_d
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_d = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_DATA) then
       v = this%entries(entry_i)%d
       dictionary_get_value_d = .true.
    else
       dictionary_get_value_d = .false.
    endif
  end function dictionary_get_value_d


  ! ****************************************************************************
  ! *
  ! *    assign_pointer() interface
  ! * 
  ! ****************************************************************************

  function dictionary_assign_pointer_i(this, key, v, case_sensitive)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    integer, intent(out), pointer, dimension(:) :: v
    logical :: dictionary_assign_pointer_i
    logical, optional :: case_sensitive

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)

    if (entry_i <= 0) then
       dictionary_assign_pointer_i = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_INTEGER_A) then
       v => this%entries(entry_i)%i_a
       dictionary_assign_pointer_i = .true.
    else
       dictionary_assign_pointer_i = .false.
    endif
  end function dictionary_assign_pointer_i

  function dictionary_assign_pointer_r(this, key, v, case_sensitive)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    real(dp), intent(out), pointer, dimension(:) :: v
    logical :: dictionary_assign_pointer_r
    logical, optional :: case_sensitive

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    
    if (entry_i <= 0) then
       dictionary_assign_pointer_r = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_REAL_A) then
       v => this%entries(entry_i)%r_a
       dictionary_assign_pointer_r = .true.
    else
       dictionary_assign_pointer_r = .false.
    endif
  end function dictionary_assign_pointer_r

  function dictionary_assign_pointer_c(this, key, v, case_sensitive)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    complex(dp), intent(out), pointer, dimension(:) :: v
    logical :: dictionary_assign_pointer_c
    logical, optional :: case_sensitive

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)

    if (entry_i <= 0) then
       dictionary_assign_pointer_c = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_COMPLEX_A) then
       v => this%entries(entry_i)%c_a
       dictionary_assign_pointer_c = .true.
    else
       dictionary_assign_pointer_c = .false.
    endif
  end function dictionary_assign_pointer_c

  function dictionary_assign_pointer_l(this, key, v, case_sensitive)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    logical, intent(out), pointer, dimension(:) :: v
    logical :: dictionary_assign_pointer_l
    logical, optional :: case_sensitive

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)

    if (entry_i <= 0) then
       dictionary_assign_pointer_l = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_LOGICAL_A) then
       v => this%entries(entry_i)%l_a
       dictionary_assign_pointer_l = .true.
    else
       dictionary_assign_pointer_l = .false.
    endif
  end function dictionary_assign_pointer_l

  function dictionary_assign_pointer_s(this, key, v, case_sensitive)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    character(len=*), intent(out), pointer, dimension(:) :: v
    logical :: dictionary_assign_pointer_s
    logical, optional :: case_sensitive

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)

    if (entry_i <= 0) then
       dictionary_assign_pointer_s = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_CHAR_A) then
       v => this%entries(entry_i)%s_a
       dictionary_assign_pointer_s = .true.
    else
       dictionary_assign_pointer_s = .false.
    endif
  end function dictionary_assign_pointer_s

  function dictionary_assign_pointer_i2(this, key, v, case_sensitive)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    integer, intent(out), pointer, dimension(:,:) :: v
    logical :: dictionary_assign_pointer_i2
    logical, optional :: case_sensitive

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)

    if (entry_i <= 0) then
       dictionary_assign_pointer_i2 = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_INTEGER_A2) then
       v => this%entries(entry_i)%i_a2
       dictionary_assign_pointer_i2 = .true.
    else
       dictionary_assign_pointer_i2 = .false.
    endif
  end function dictionary_assign_pointer_i2

  function dictionary_assign_pointer_r2(this, key, v, case_sensitive)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    real(dp), intent(out), pointer, dimension(:,:) :: v
    logical :: dictionary_assign_pointer_r2
    logical, optional :: case_sensitive

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    
    if (entry_i <= 0) then
       dictionary_assign_pointer_r2 = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_REAL_A2) then
       v => this%entries(entry_i)%r_a2
       dictionary_assign_pointer_r2 = .true.
    else
       dictionary_assign_pointer_r2 = .false.
    endif
  end function dictionary_assign_pointer_r2


  ! ****************************************************************************
  ! *
  ! *    set_value_pointer() interface
  ! * 
  ! ****************************************************************************

  subroutine dictionary_set_value_pointer_i(this, key, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in), pointer :: ptr(:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_INTEGER_A
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry_i = add_entry(this, key, entry)
    entry%len = size(ptr)
    this%entries(entry_i)%i_a => ptr
    call finalise(entry)

  end subroutine dictionary_set_value_pointer_i

  subroutine dictionary_set_value_pointer_r(this, key, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), intent(in), pointer :: ptr(:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_REAL_A
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry_i = add_entry(this, key, entry)
    entry%len = size(ptr)
    this%entries(entry_i)%r_a => ptr
    call finalise(entry)

  end subroutine dictionary_set_value_pointer_r

  subroutine dictionary_set_value_pointer_c(this, key, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    complex(dp), intent(in), pointer :: ptr(:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_COMPLEX_A
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry_i = add_entry(this, key, entry)
    entry%len = size(ptr)
    this%entries(entry_i)%c_a => ptr
    call finalise(entry)

  end subroutine dictionary_set_value_pointer_c

  subroutine dictionary_set_value_pointer_l(this, key, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in), pointer :: ptr(:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_LOGICAL_A
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry_i = add_entry(this, key, entry)
    entry%len = size(ptr)
    this%entries(entry_i)%l_a => ptr
    call finalise(entry)

  end subroutine dictionary_set_value_pointer_l

  subroutine dictionary_set_value_pointer_s(this, key, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(in), pointer :: ptr(:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_CHAR_A
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry_i = add_entry(this, key, entry)
    entry%len = size(ptr)
    this%entries(entry_i)%s_a => ptr
    call finalise(entry)

  end subroutine dictionary_set_value_pointer_s

  subroutine dictionary_set_value_pointer_i2(this, key, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in), pointer :: ptr(:,:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_INTEGER_A2
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry_i = add_entry(this, key, entry)
    entry%len2 = shape(ptr)
    this%entries(entry_i)%i_a2 => ptr
    call finalise(entry)

  end subroutine dictionary_set_value_pointer_i2

  subroutine dictionary_set_value_pointer_r2(this, key, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), intent(in), pointer :: ptr(:,:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_REAL_A2
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry_i = add_entry(this, key, entry)
    entry%len2 = shape(ptr)
    this%entries(entry_i)%r_a2 => ptr
    call finalise(entry)

  end subroutine dictionary_set_value_pointer_r2


  ! ****************************************************************************
  ! *
  ! *    add_array() interface
  ! * 
  ! ****************************************************************************

  subroutine dictionary_add_array_i(this, key, value, len, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value
    integer, intent(in) :: len
    integer, pointer, optional, intent(out) :: ptr(:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_INTEGER_A
    entry%len = len
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%i_a(len))
    this%entries(entry_i)%i_a(:) = value
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%i_a

  end subroutine dictionary_add_array_i

  subroutine dictionary_add_array_r(this, key, value, len, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: value
    integer, intent(in) :: len
    real(dp), intent(out), pointer, optional :: ptr(:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_REAL_A
    entry%len = len
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%r_a(len))
    this%entries(entry_i)%r_a(:) = value
    call finalise(entry)

    if (present(ptr)) ptr => this%entries(entry_i)%r_a

  end subroutine dictionary_add_array_r

  subroutine dictionary_add_array_c(this, key, value, len, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    complex(dp), intent(in) :: value
    integer, intent(in) :: len
    complex(dp), pointer, optional, intent(out) :: ptr(:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_COMPLEX_A
    entry%len = len
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%c_a(len))
    this%entries(entry_i)%c_a(:) = value
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%c_a

  end subroutine dictionary_add_array_c

  subroutine dictionary_add_array_l(this, key, value, len, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in) :: value
    integer, intent(in) :: len
    logical, pointer, optional, intent(out) :: ptr(:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_LOGICAL_A
    entry%len = len
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%l_a(len))
    this%entries(entry_i)%l_a(:) = value
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%l_a

  end subroutine dictionary_add_array_l

  subroutine dictionary_add_array_s(this, key, value, len, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value
    integer, intent(in) :: len
    character(len=*), pointer, optional, intent(out) :: ptr(:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_CHAR_A
    entry%len = len
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%s_a(len))
    this%entries(entry_i)%s_a(:) = value
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%s_a

  end subroutine dictionary_add_array_s

  subroutine dictionary_add_array_i2(this, key, value, len2, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value
    integer, intent(in) :: len2(2)
    integer, pointer, optional, intent(out) :: ptr(:,:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_INTEGER_A2
    entry%len2 = len2
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%i_a2(len2(1),len2(2)))
    this%entries(entry_i)%i_a2(:,:) = value
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%i_a2

  end subroutine dictionary_add_array_i2

  subroutine dictionary_add_array_r2(this, key, value, len2, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: value
    integer, intent(in) :: len2(2)
    real(dp), intent(out), pointer, optional :: ptr(:,:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_REAL_A2
    entry%len2 = len2
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) allocate(this%entries(entry_i)%r_a2(len2(1),len2(2)))
    this%entries(entry_i)%r_a2(:,:) = value
    call finalise(entry)

    if (present(ptr)) ptr => this%entries(entry_i)%r_a2

  end subroutine dictionary_add_array_r2


  ! ****************************************************************************
  ! *
  ! *    Miscellaneous routines
  ! * 
  ! ****************************************************************************

  subroutine dictionary_remove_value(this, key)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key

    integer entry_i

    entry_i = lookup_entry_i(this, key)
    if (entry_i > 0) then
       call remove_entry(this, entry_i)
    endif
  end subroutine dictionary_remove_value

  subroutine dictionary_read_string(this, str, append)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: str
    logical, optional, intent(in) :: append !% If true, append to dictionary (default false)

    logical :: do_append
    character(len=dict_field_length) :: field
    integer equal_pos
    character(len=dict_field_length), dimension(dict_n_fields) :: final_fields
    character(len=dict_field_length) :: key, value
    integer :: i, num_pairs

    do_append = optional_default(.false., append)

    if (.not. do_append) then
       call initialise(this)
    end if

    call split_string(str, ' ,', '""'//"''"//'{}', final_fields, num_pairs, matching=.true.)

    ! Set the new entries
    do i=1,num_pairs
       field = final_fields(i)
       equal_pos = index(trim(field),'=')
       if (equal_pos == 0) then
          key=field
          value=''
       else if (equal_pos == 1) then
          call system_abort("Malformed field '"//trim(field)//"'")
       else
          key = field(1:equal_pos-1)
          value = field(equal_pos+1:len(trim(field)))
       endif
       call print("dictionary_read_string key='"//trim(key)//"' value='"//trim(value)//"'", PRINT_NERD)
       if (len_trim(value) > value_len) then
          call system_abort("dictionary_read_string: value "//trim(value)//" too long")
          return
       end if

       if (.not. dictionary_parse_value(this,key,value)) &
            call system_abort('dictionary_read_string: error parsing '//trim(key)//'='//trim(value))

    end do

  end subroutine dictionary_read_string


  !% Parse a string corresponding to a dictionary value and
  !% attempt to work out what type it should be, then set it
  function dictionary_parse_value(this, key, strvalue, char_a_sep) result(status)
    type(Dictionary), intent(inout) :: this
    character(*), intent(in) :: key
    character(*), intent(in) :: strvalue
    character(1), optional, intent(in) :: char_a_sep
    logical :: status

    character(len=dict_field_length), dimension(dict_n_fields) :: fields
    character(len=len(strvalue)) :: datastr, shapestr, myvalue
    integer :: num_fields, i, j
    real(dp) :: r
    logical :: l, err, all_int, all_real, all_logical
    type(DictData) :: data
    integer, allocatable, dimension(:) :: i_a
    real(dp), allocatable, dimension(:) :: r_a
    logical, allocatable, dimension(:) :: l_a
    integer, allocatable, dimension(:,:) :: i_a2
    real(dp), allocatable, dimension(:,:) :: r_a2
    character(1) :: my_char_a_sep
    logical got_2d
    integer shape_2d(2)
    integer :: parse_status

    my_char_a_sep = optional_default(',',char_a_sep)

    ! First, let's check if value is arbitrary data, represented
    ! as 'DATA"i1 i2 i3 i4"' with i1,i2, etc. integers
    if (strvalue(1:4) == 'DATA') then
       datastr = strvalue(5:len(trim(strvalue)))
       call parse_string(datastr, ' ',fields, num_fields)

       allocate(data%d(num_fields))
       do j=1,num_fields
          data%d(j) = string_to_int(fields(j),err)
          if (err) then
             status = .false.
             return
          end if
       end do

       call set_value(this,key,data)
       deallocate(data%d)
       status = .true.
       return
    end if

    ! 2D arrays are represented with shape in ()s at beginning of value string
    got_2d = .false.
    if (strvalue(1:1) == '(') then
       got_2d = .true.
       shapestr = strvalue(2:index(strvalue,')')-1)
       read (shapestr, *) shape_2d
       myvalue = strvalue(index(strvalue,')')+1:len_trim(strvalue))
    else
       myvalue = strvalue
    end if

    ! Otherwise, start by splitting value into fields
    call parse_string(myvalue, ' ', fields, num_fields, status=parse_status)

    if (parse_status /= 0) then
       ! parsing failed, treat as string
       call set_value(this,key,myvalue)
       status = .true.
       return
    else if (num_fields == 0) then
       ! Nothing there, assume it's logical and true
       l = .true.
       call set_value(this, key, l)
       status = .true.
       return

    else if (num_fields == 1) then

       ! Scalar data, try to guess type

       do i=1, len_trim(fields(1))
          if (fields(1)(i:i) == '/') then
             ! slash is technically a delimiter, but we don't like that, so we'll call it a string manually (otherswise string_to_real will be confused)
             call set_value(this,key,fields(1))
             status = .true.
             return
          endif
       end do

       ! Is it an integer?
       i = string_to_int(fields(1),err)
       if (.not. err) then
          call set_value(this, key, i)
          status = .true.
          return
       end if

       ! Is it a logical?
       if (trim(fields(1)) == 'T' .or. trim(fields(1)) == 't' .or. &
            trim(fields(1)) == 'F' .or. trim(fields(1)) == 'f') then
          l = string_to_logical(fields(1))
          call set_value(this, key, l)
          status = .true.
          return
       end if

       ! How about a real?
       r = string_to_real(fields(1),err)
       if (.not. err) then
          call set_value(this, key, r)
          status = .true.
          return
       end if

       ! Add complex scalar here...

       ! Does it contain one or more array separator characters?
       if (scan(fields(1),my_char_a_sep) /= 0) then
          call parse_string(myvalue, my_char_a_sep, fields, num_fields)
          call set_value(this,key,fields(1:num_fields))
          status = .true.
          return
       else
          ! All else has failed, treat it as a single string
          call set_value(this,key,fields(1))
          status = .true.
          return
       end if

    else
       ! Array data, again have to guess type

       allocate(i_a(num_fields),r_a(num_fields),l_a(num_fields))

       ! Are all fields integers?
       all_int = .true.
       do j=1,num_fields
          i_a(j) = string_to_int(fields(j),err)

          if (err) then
             all_int = .false.
             exit ! Found something that's not an integer
          end if
       end do

       if (all_int) then
          if (got_2d) then
             allocate(i_a2(shape_2d(1), shape_2d(2)))
             i_a2 = reshape(i_a, shape_2d)
             call set_value(this,key,i_a2)
          else
             call set_value(this,key,i_a)
          end if
          status = .true.
          deallocate(i_a,r_a,l_a)
          if (allocated(i_a2)) deallocate(i_a2)
          return
       end if

       ! Are all fields logicals?
       all_logical = .true.
       do j=1,num_fields
          if (.not. (trim(fields(j)) == 'T' .or. trim(fields(j)) == 't' .or. &
               trim(fields(j)) == 'F' .or. trim(fields(j)) == 'f')) then 
             all_logical = .false.
             exit
          end if
          l_a(j) = string_to_logical(fields(j),err)
       end do

       if (all_logical) then
          call set_value(this,key,l_a)
          status = .true.
          deallocate(i_a,r_a,l_a)
          return
       end if

       ! Are all fields real numbers?
       all_real = .true.
       do j=1,num_fields
          r_a(j) = string_to_real(fields(j),err)

          if (err) then
             all_real = .false.
             exit
          end if
       end do

       if (all_real) then
          if (got_2d) then
             allocate(r_a2(shape_2d(1), shape_2d(2)))
             r_a2 = reshape(r_a, shape_2d)
             call set_value(this,key,r_a2)
          else
             call set_value(this,key,r_a)
          end if
          status = .true.
          deallocate(i_a,r_a,l_a)
          if (allocated(r_a2)) deallocate(r_a2)
          return
       end if

       ! Add complex array here...

       ! We're left with strings. Does it contain array seperator?
       if (scan(myvalue,my_char_a_sep) /= 0) then
          call parse_string(myvalue, my_char_a_sep, fields, num_fields)
          call set_value(this,key,fields(1:num_fields))
          status = .true.
          return
       else
          ! Fall back option: treat entire myvalue as single string
          call set_value(this,key,myvalue)
          status =.true.
          return
       end if
    end if

  end function dictionary_parse_value


  function dictionary_write_string(this, real_format, entry_sep, char_a_sep) result(str)
    type(Dictionary), intent(in) :: this
    character(len=*), optional, intent(in) :: real_format !% Output format for reals, default is 'f9.3'
    character(1), optional, intent(in) :: entry_sep !% Entry seperator, default is single space
    character(1), optional, intent(in) :: char_a_sep !% Output separator for character arrays, default is ','
    character(len=2048) :: str ! Beware, fixed length string

    integer :: i,j
    character(1) :: my_char_a_sep, my_entry_sep
    character(255) :: my_real_format


    my_real_format = optional_default('f9.3',real_format)
    my_char_a_sep = optional_default(',',char_a_sep)
    my_entry_sep = optional_default(' ',entry_sep)

    str = ''

    do i=1,this%N

       if (i == 1) then
          str = trim(this%keys(i))//'='
       else
          str = trim(str)//my_entry_sep//trim(this%keys(i))//'='
       end if
       select case(this%entries(i)%type)

       case(T_INTEGER)
          str = trim(str)//this%entries(i)%i

       case(T_REAL)
          write (line,'('//my_real_format//')') this%entries(i)%r
          str = trim(str)//trim(adjustl(line))

       case(T_COMPLEX)
          str = trim(str)//this%entries(i)%c

       case(T_LOGICAL)
          str = trim(str)//this%entries(i)%l

       case(T_CHAR)
          if (index(trim(this%entries(i)%s), ' ') == 0) then
             str = trim(str)//trim(this%entries(i)%s)
          else
             str = trim(str)//'"'//trim(this%entries(i)%s)//'"'
          end if

       case(T_INTEGER_A)
          str = trim(str)//'"'//this%entries(i)%i_a//'"'

       case(T_REAL_A)
          write (line, '('//size(this%entries(i)%r_a)//trim(my_real_format)//')') this%entries(i)%r_a
          str = trim(str)//'"'//trim(adjustl(line))//'"'

       case(T_COMPLEX_A)
          str = trim(str)//'"'//this%entries(i)%c_a//'"'

       case(T_LOGICAL_A)
          str = trim(str)//'"'//this%entries(i)%l_a//'"'

       case(T_CHAR_A)
          str = trim(str)//'"'
          do j=1,size(this%entries(i)%s_a)
             str = trim(str)//trim(this%entries(i)%s_a(j))//my_char_a_sep
          end do
          str = trim(str(1:len(trim(str))-1))//'"'

       case(T_INTEGER_A2)
          str = trim(str)//'"('//shape(this%entries(i)%i_a2)//') '//reshape(this%entries(i)%i_a2,(/size(this%entries(i)%i_a2)/))//'"'

       case(T_REAL_A2)
          str = trim(str)//'"('//shape(this%entries(i)%r_a2)//') '//reshape(this%entries(i)%r_a2,(/size(this%entries(i)%r_a2)/))//'"'

       case(T_DATA)
          str = trim(str)//'DATA"'//this%entries(i)%d%d//'"'
       end select
    end do


  end function  dictionary_write_string


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !% OMIT
  subroutine remove_entry(this, entry_i)
    type(Dictionary), intent(inout) :: this
    integer, intent(in) :: entry_i

    if (entry_i <= 0 .or. entry_i > this%N) &
         call system_abort ("Called remove_entry with invalid entry number")

    call finalise(this%entries(entry_i))

    if (entry_i < this%N) then
       this%keys(entry_i:this%N-1) = this%keys(entry_i+1:this%N)
       this%entries(entry_i:this%N-1) = this%entries(entry_i+1:this%N)
    endif

    this%N = this%N - 1
  end subroutine remove_entry

  !% OMIT
  function add_entry(this, key, entry, array_alloc) result(entry_i)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    type(DictEntry), intent(in) :: entry
    logical, intent(out), optional :: array_alloc

    integer entry_i

    if (this%N >= size(this%entries)) then
       call extend_entries(this, n_entry_block)
    endif

    entry_i = lookup_entry_i(this, key)
    if (present(array_alloc)) array_alloc = .true.

    if (entry_i > 0) then
       ! we only free entry if new type/size/shape is incompatible with existing entry
       ! or if new entry will not be allocated by us (i.e. entry%own_data == .false.)
       if (.not. entry%own_data .or. &
            this%entries(entry_i)%type == T_DATA .or. &
            this%entries(entry_i)%type /= entry%type .or. &
            this%entries(entry_i)%len /= entry%len .or. &
            any(this%entries(entry_i)%len2 /= entry%len2)) then
          call finalise(this%entries(entry_i))
          this%entries(entry_i) = entry
       else
          if (present(array_alloc)) array_alloc = .false.
          ! we're keeping allocation, so we overwrite only scalar components of entry
          this%entries(entry_i)%type = entry%type
          this%entries(entry_i)%len = entry%len
          this%entries(entry_i)%len2 = entry%len2
          this%entries(entry_i)%i = entry%i
          this%entries(entry_i)%r = entry%r
          this%entries(entry_i)%c = entry%c
          this%entries(entry_i)%l = entry%l
          this%entries(entry_i)%s = entry%s
       end if
     else
       this%N = this%N + 1
       entry_i = this%N
       this%keys(entry_i) = key
       this%entries(entry_i) = entry
    endif
  end function add_entry


  !% OMIT
  subroutine extend_entries(this, n)
    type(Dictionary), intent(inout) :: this
    integer, intent(in) :: n

    type(DictEntry), allocatable :: t_entries(:)
    character(len=key_len), allocatable :: t_keys(:)
    integer old_size

    if (allocated(this%entries)) then
       allocate(t_entries(size(this%entries)))
       allocate(t_keys(size(this%entries)))
       t_entries = this%entries
       t_keys = this%keys

       old_size = size(this%entries)

       deallocate(this%entries)
       deallocate(this%keys)

       allocate(this%entries(old_size + n))
       allocate(this%keys(old_size + n))

       this%entries(1:this%N) = t_entries(1:this%N)
       this%keys(1:this%N) = t_keys(1:this%N)

       if (allocated(t_entries)) deallocate(t_entries)
       deallocate (t_keys)
    else
       allocate(this%entries(n))
       allocate(this%keys(n))
    endif

  end subroutine extend_entries

  !% Convert a word to lower case
  function lower_case(word)
    character(*) , intent(in) :: word
    character(len(word)) :: lower_case

    integer :: i,ic,nlen

    nlen = len(word)
    do i=1,nlen
       ic = ichar(word(i:i))
       if (ic >= 65 .and. ic <= 90) then
          lower_case(i:i) = char(ic+32)
       else
          lower_case(i:i) = char(ic)
       end if
    end do
  end function lower_case

  !% OMIT
  function lookup_entry_i(this, key, case_sensitive)
    type(Dictionary), intent(in) :: this
    character(len=*) :: key
    logical, optional :: case_sensitive
    integer :: lookup_entry_i

    logical :: do_case_sensitive
    integer i

    do_case_sensitive = optional_default(.false., case_sensitive)

    lookup_entry_i = -1
    do i=1, this%N
       if(do_case_sensitive) then
          if (trim(this%keys(i)) == trim(key)) then
             lookup_entry_i = i
             return
          endif
       else
          if (trim(lower_case(this%keys(i))) == trim(lower_case(key))) then
             lookup_entry_i = i
             return
          endif
       end if
    end do
  end function lookup_entry_i

  !% Return true if 'key' is in Dictionary or false if not
  function dictionary_has_key(this, key, case_sensitive)
    type(Dictionary), intent(in) :: this
    character(len=*) :: key
    logical :: dictionary_has_key
    logical, optional :: case_sensitive

    dictionary_has_key = lookup_entry_i(this, key, case_sensitive) /= -1

  end function dictionary_has_key

  !% Return a dictionary that is a subset of this,
  !% with only the keys in the arrays 'keys' present.
  function dictionary_subset(this, keys, case_sensitive)
    type(Dictionary), intent(in) :: this
    character(len=*), dimension(:) :: keys
    type(Dictionary) :: dictionary_subset
    logical, optional :: case_sensitive

    integer :: i, entry_i

    call initialise(dictionary_subset)

    do i=1,size(keys)
       entry_i = lookup_entry_i(this, keys(i), case_sensitive)
       if (entry_i <= 0) &
            call system_abort('dictionary_subset: key '//trim(keys(i))//' not in dictionary')
       entry_i = add_entry(dictionary_subset, keys(i), this%entries(entry_i))
    end do

  end function dictionary_subset


  !% Swap the positions of two entries in the dictionary. Arrays are not moved in memory.
  subroutine dictionary_swap(this, key1, key2, case_sensitive)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key1, key2
    logical, optional :: case_sensitive

    integer :: i1,i2
    type(DictEntry) :: tmp_entry
    character(len=key_len) :: tmp_key

    i1 = lookup_entry_i(this,key1,case_sensitive)
    i2 = lookup_entry_i(this,key2,case_sensitive)

    if (i1 <= 0) &
         call system_abort('dictionary_swap: key '//key1//' not in dictionary')
    if (i2 <= 0) &
         call system_abort('dictionary_swap: key '//key2//' not in dictionary')

    tmp_entry = this%entries(i2)
    tmp_key   = this%keys(i2)
    this%entries(i2) = this%entries(i1)
    this%keys(i2)    = this%keys(i1)
    this%entries(i1) = tmp_entry
    this%keys(i1)    = tmp_key

  end subroutine dictionary_swap

  subroutine dictionary_bcast(mpi, dict)
    type(MPI_context), intent(in) :: mpi
    type(Dictionary), intent(inout) :: dict

    integer :: i, size_tmp, shape_tmp(2)

    if (.not. mpi%active) return

    if (mpi%my_proc == 0) then
       call bcast(mpi, dict%n)
       do i=1,dict%n
          call bcast(mpi, dict%entries(i)%type)
          call bcast(mpi, dict%keys(i))

          if (dict%entries(i)%type == T_NONE) then
             ! nothing to do
          else if (dict%entries(i)%type == T_INTEGER) then
             call bcast(mpi, dict%entries(i)%i)
          else if (dict%entries(i)%type == T_REAL) then
             call bcast(mpi, dict%entries(i)%r)
          else if (dict%entries(i)%type == T_COMPLEX) then
             call bcast(mpi, dict%entries(i)%c)
          else if (dict%entries(i)%type == T_LOGICAL) then
             call bcast(mpi, dict%entries(i)%l)
          else if (dict%entries(i)%type == T_CHAR) then
             call bcast(mpi, dict%entries(i)%s)
          else if (dict%entries(i)%type == T_INTEGER_A) then
             size_tmp = size(dict%entries(i)%i_a)
             call bcast(mpi, size_tmp)
             call bcast(mpi, dict%entries(i)%i_a)
          else if (dict%entries(i)%type == T_REAL_A) then
             size_tmp = size(dict%entries(i)%r_a)
             call bcast(mpi, size_tmp)
             call bcast(mpi, dict%entries(i)%r_a)
          else if (dict%entries(i)%type == T_COMPLEX_A) then
             size_tmp = size(dict%entries(i)%c_a)
             call bcast(mpi, size_tmp)
             call bcast(mpi, dict%entries(i)%c_a)
          else if (dict%entries(i)%type == T_LOGICAL_A) then
             size_tmp = size(dict%entries(i)%l_a)
             call bcast(mpi, size_tmp)
             call bcast(mpi, dict%entries(i)%l_a)
          else if (dict%entries(i)%type == T_CHAR_A) then
             size_tmp = size(dict%entries(i)%s_a)
             call bcast(mpi, size_tmp)
             call bcast(mpi, dict%entries(i)%s_a)
          else if (dict%entries(i)%type == T_INTEGER_A2) then
             shape_tmp = shape(dict%entries(i)%i_a2)
             call bcast(mpi, shape_tmp)
             call bcast(mpi, dict%entries(i)%i_a2)
          else if (dict%entries(i)%type == T_REAL_A2) then
             shape_tmp = shape(dict%entries(i)%r_a2)
             call bcast(mpi, shape_tmp)
             call bcast(mpi, dict%entries(i)%r_a2)
          else if (dict%entries(i)%type == T_DATA) then
             size_tmp = size(dict%entries(i)%d%d)
             call bcast(mpi, size_tmp)
             call bcast(mpi, dict%entries(i)%d%d)
          end if
       end do
    else
       call finalise(dict)
       call bcast(mpi, dict%n)
       allocate(dict%keys(dict%n))
       allocate(dict%entries(dict%n))
       do i=1,dict%n
          call bcast(mpi, dict%entries(i)%type)
          call bcast(mpi, dict%keys(i))

          if (dict%entries(i)%type == T_NONE) then
             ! nothing to do
          else if (dict%entries(i)%type == T_INTEGER) then
             call bcast(mpi, dict%entries(i)%i)
          else if (dict%entries(i)%type == T_REAL) then
             call bcast(mpi, dict%entries(i)%r)
          else if (dict%entries(i)%type == T_COMPLEX) then
             call bcast(mpi, dict%entries(i)%c)
          else if (dict%entries(i)%type == T_LOGICAL) then
             call bcast(mpi, dict%entries(i)%l)
          else if (dict%entries(i)%type == T_CHAR) then
             call bcast(mpi, dict%entries(i)%s)
          else if (dict%entries(i)%type == T_INTEGER_A) then
             size_tmp = size(dict%entries(i)%i_a)
             call bcast(mpi, size_tmp)
             allocate(dict%entries(i)%i_a(size_tmp))
             call bcast(mpi, dict%entries(i)%i_a)
          else if (dict%entries(i)%type == T_REAL_A) then
             size_tmp = size(dict%entries(i)%r_a)
             call bcast(mpi, size_tmp)
             allocate(dict%entries(i)%r_a(size_tmp))
             call bcast(mpi, dict%entries(i)%r_a)
          else if (dict%entries(i)%type == T_COMPLEX_A) then
             size_tmp = size(dict%entries(i)%c_a)
             call bcast(mpi, size_tmp)
             allocate(dict%entries(i)%c_a(size_tmp))
             call bcast(mpi, dict%entries(i)%c_a)
          else if (dict%entries(i)%type == T_LOGICAL_A) then
             size_tmp = size(dict%entries(i)%l_a)
             call bcast(mpi, size_tmp)
             allocate(dict%entries(i)%l_a(size_tmp))
             call bcast(mpi, dict%entries(i)%l_a)
          else if (dict%entries(i)%type == T_CHAR_A) then
             size_tmp = size(dict%entries(i)%s_a)
             call bcast(mpi, size_tmp)
             allocate(dict%entries(i)%s_a(size_tmp))
             call bcast(mpi, dict%entries(i)%s_a)
          else if (dict%entries(i)%type == T_INTEGER_A2) then
             shape_tmp = shape(dict%entries(i)%i_a2)
             call bcast(mpi, shape_tmp)
             allocate(dict%entries(i)%i_a2(shape_tmp(1),shape_tmp(2)))
             call bcast(mpi, dict%entries(i)%i_a2)
          else if (dict%entries(i)%type == T_REAL_A2) then
             shape_tmp = shape(dict%entries(i)%r_a2)
             call bcast(mpi, shape_tmp)
             allocate(dict%entries(i)%r_a2(shape_tmp(1),shape_tmp(2)))
             call bcast(mpi, dict%entries(i)%r_a2)
          else if (dict%entries(i)%type == T_DATA) then
             size_tmp = size(dict%entries(i)%d%d)
             call bcast(mpi, size_tmp)
             allocate(dict%entries(i)%d%d(size_tmp))
             call bcast(mpi, dict%entries(i)%d%d)
          end if
       end do
    end if

  end subroutine dictionary_bcast

  !% Make a deep copy of 'from' in 'this', allocating new memory for array components
  subroutine dictionary_deepcopy(this, from)
    type(Dictionary), intent(inout) :: this
    type(Dictionary), intent(in) :: from
    
    integer :: i

    call initialise(this)
    do i=1,from%n

       select case(from%entries(i)%type)
       case(T_INTEGER)
          call set_value(this, from%keys(i), from%entries(i)%i)

       case(T_REAL)
          call set_value(this, from%keys(i), from%entries(i)%r)

       case(T_COMPLEX)
          call set_value(this, from%keys(i), from%entries(i)%c)

       case(T_LOGICAL)
          call set_value(this, from%keys(i), from%entries(i)%l)

       case(T_CHAR)
          call set_value(this, from%keys(i), from%entries(i)%s)

       case(T_INTEGER_A)
          call set_value(this, from%keys(i), from%entries(i)%i_a)

       case(T_REAL_A)
          call set_value(this, from%keys(i), from%entries(i)%r_a)

       case(T_COMPLEX_A)
          call set_value(this, from%keys(i), from%entries(i)%c_a)

       case(T_LOGICAL_A)
          call set_value(this, from%keys(i), from%entries(i)%l_a)

       case(T_CHAR_A)
          call set_value(this, from%keys(i), from%entries(i)%s_a)

       case(T_INTEGER_A2)
          call set_value(this, from%keys(i), from%entries(i)%i_a2)

       case(T_REAL_A2)
          call set_value(this, from%keys(i), from%entries(i)%r_a2)

       case(T_DATA)
          call set_value(this, from%keys(i), from%entries(i)%d)
       end select
    end do
    
  end subroutine dictionary_deepcopy

#ifdef HAVE_QUIPPY
  subroutine dictionary_get_array(this, key, nd, dtype, dshape, dloc)
    type(Dictionary), intent(in) :: this
    character(len=*), intent(in) :: key
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer(NUMPY_PTR_SIZE), intent(out) :: dloc
    
    integer entry_i

    nd = 0
    dtype = 0
    dshape(:) = 0
    dloc = 0

    entry_i = lookup_entry_i(this, key)
    if (entry_i == 0) return

    select case(this%entries(entry_i)%type)
       case(T_INTEGER_A)
          nd = 1
          dtype = NUMPY_INTEGER
          dshape(1) = size(this%entries(entry_i)%i_a)
          dloc = loc(this%entries(entry_i)%i_a)

       case(T_REAL_A)
          nd = 1
          dtype = NUMPY_REAL_DP
          dshape(1) = size(this%entries(entry_i)%r_a)
          dloc = loc(this%entries(entry_i)%r_a)

       case(T_COMPLEX_A)
          nd = 1
          dtype = NUMPY_COMPLEX
          dshape(1) = size(this%entries(entry_i)%c_a)
          dloc = loc(this%entries(entry_i)%c_a)

       case(T_LOGICAL_A)
          nd = 1
          dtype = NUMPY_LOGICAL
          dshape(1) = size(this%entries(entry_i)%l_a)
          dloc = loc(this%entries(entry_i)%l_a)

       case(T_CHAR_A)
          nd = 2
          dtype = NUMPY_CHAR
          dshape(1) = len(this%entries(entry_i)%s_a(1))
          dshape(2) = size(this%entries(entry_i)%s_a)
          dloc = loc(this%entries(entry_i)%s_a)

       case(T_INTEGER_A2)
          nd = 2
          dtype = NUMPY_INTEGER
          dshape(1) = size(this%entries(entry_i)%i_a2, 1)
          dshape(2) = size(this%entries(entry_i)%i_a2, 2)
          dloc = loc(this%entries(entry_i)%i_a2)

       case(T_REAL_A2)
          nd = 2
          dtype = NUMPY_REAL_DP
          dshape(1) = size(this%entries(entry_i)%r_a2, 1)
          dshape(2) = size(this%entries(entry_i)%r_a2, 2)
          dloc = loc(this%entries(entry_i)%r_a2)
    end select

  end subroutine dictionary_get_array
#endif

end module dictionary_module
