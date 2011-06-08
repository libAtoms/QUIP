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
!% 'integer', 'real(dp)', 'complex(dp)', 'logical', extendable_str
!% or a 1-D array of any of those. 2-D arrays of integers and reals are also supported.

#include "error.inc"

module dictionary_module

  use system_module
  use system_module, only: string_to_real
  use linearalgebra_module
  use mpi_context_module
  use extendable_str_module

  implicit none

  private
  
  public :: T_NONE, T_INTEGER, T_REAL, T_COMPLEX, T_LOGICAL, &
       T_INTEGER_A, T_REAL_A, T_COMPLEX_A, T_LOGICAL_A, &
       T_CHAR, T_CHAR_A, T_DATA, T_INTEGER_A2, T_REAL_A2, T_DICT
  integer, parameter :: &
       T_NONE = 0, &
       T_INTEGER = 1, T_REAL =2, T_COMPLEX = 3, T_LOGICAL=4, &
       T_INTEGER_A = 5, T_REAL_A = 6, T_COMPLEX_A = 7, T_LOGICAL_A = 8, &
       T_CHAR = 9, T_CHAR_A = 10, T_DATA = 11, T_INTEGER_A2 = 12, T_REAL_A2 = 13, T_DICT = 14 !% OMIT

  ! Maintained for backwards compatibility with old NetCDF files using type attribute
  public :: PROPERTY_INT, PROPERTY_REAL, PROPERTY_STR, PROPERTY_LOGICAL
  integer, parameter :: &
       PROPERTY_INT = 1, PROPERTY_REAL = 2, PROPERTY_STR = 3, PROPERTY_LOGICAL = 4

  public :: C_KEY_LEN, DICT_FIELD_LENGTH, DICT_N_FIELDS
  integer, parameter :: C_KEY_LEN = 256
  integer, parameter :: DICT_FIELD_LENGTH = 1024  !% Maximum field width during parsing
  integer, parameter :: DICT_STRING_LENGTH = 20480  !% Maximum string width during parsing, should be greater than the paramreader's string_length
  integer, parameter :: DICT_N_FIELDS = 100       !% Maximum number of fields during parsing

  public :: dictdata
  type DictData
     integer, dimension(:), allocatable :: d
  end type DictData

  public :: dictentry
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
     type(extendable_str) s

     integer, pointer :: i_a(:) => null()
     real(dp), pointer :: r_a(:) => null()
     complex(dp), pointer :: c_a(:) => null()
     logical, pointer :: l_a(:) => null()
     character(len=1), pointer :: s_a(:,:) => null()

     integer, pointer :: i_a2(:,:) => null()
     real(dp), pointer :: r_a2(:,:) => null()

     type(DictData) :: d

  end type DictEntry

  integer, parameter :: n_entry_block = 10 !% OMIT

  public Dictionary
  type Dictionary
     integer :: N !% number of entries in use
     type(extendable_str), allocatable :: keys(:) !% array of keys
     type(DictEntry), allocatable :: entries(:)    !% array of entries
     integer :: cache_invalid !% non-zero on exit from set_value(), set_value_pointer(), add_array(), remove_entry() if any array memory locations changed
  end type Dictionary

  public c_dictionary_ptr_type
  type c_dictionary_ptr_type
     type(Dictionary), pointer :: p
  end type c_dictionary_ptr_type
  
  !% Initialise a new empty dictionary
  public initialise
  interface initialise
     module procedure dictionary_initialise
  end interface initialise

  !% Finalise dictionary
  public finalise
  interface finalise
     module procedure dictionary_finalise, dictentry_finalise
  end interface finalise

  !% Print a DictEntry or a Dictionary
  public print
  interface print
     module procedure dictentry_print, dictionary_print
  end interface print

  public print_keys
  interface print_keys
     module procedure dictionary_print_keys
  end interface print_keys

  !% Set a value in a Dictionary
  public set_value
  interface set_value
     module procedure dictionary_set_value_none
     module procedure dictionary_set_value_i, dictionary_set_value_r, dictionary_set_value_c, dictionary_set_value_l
     module procedure dictionary_set_value_i_a, dictionary_set_value_r_a, dictionary_set_value_c_a, dictionary_set_value_l_a
     module procedure dictionary_set_value_s, dictionary_set_value_s_es, dictionary_set_value_s_a2
     module procedure dictionary_set_value_s_a
     module procedure dictionary_set_value_d
     module procedure dictionary_set_value_i_a2
     module procedure dictionary_set_value_r_a2
     module procedure dictionary_set_value_dict
  end interface set_value

  public set_value_pointer
  interface set_value_pointer
     module procedure dictionary_set_value_pointer_i
     module procedure dictionary_set_value_pointer_r
     module procedure dictionary_set_value_pointer_c
     module procedure dictionary_set_value_pointer_l
     module procedure dictionary_set_value_pointer_s
     module procedure dictionary_set_value_pointer_i2
     module procedure dictionary_set_value_pointer_r2
  end interface set_value_pointer

  !% Get a value from a Dictionary
  public get_value
  interface get_value
     module procedure dictionary_get_value_i, dictionary_get_value_r, dictionary_get_value_c, dictionary_get_value_l
     module procedure dictionary_get_value_i_a, dictionary_get_value_r_a, dictionary_get_value_c_a, dictionary_get_value_l_a
     module procedure dictionary_get_value_s, dictionary_get_value_s_es, dictionary_get_value_s_a, dictionary_get_value_s_a2
     module procedure dictionary_get_value_d
     module procedure dictionary_get_value_i_a2
     module procedure dictionary_get_value_r_a2
     module procedure dictionary_get_value_dict
  end interface get_value

  public assign_pointer
  interface assign_pointer
     module procedure dictionary_assign_pointer_r0
     module procedure dictionary_assign_pointer_i
     module procedure dictionary_assign_pointer_r
     module procedure dictionary_assign_pointer_c
     module procedure dictionary_assign_pointer_l
     module procedure dictionary_assign_pointer_s
     module procedure dictionary_assign_pointer_i2
     module procedure dictionary_assign_pointer_r2
  end interface assign_pointer

  public add_array
  interface add_array
     module procedure dictionary_add_array_i
     module procedure dictionary_add_array_r
     module procedure dictionary_add_array_c
     module procedure dictionary_add_array_l
     module procedure dictionary_add_array_s
     module procedure dictionary_add_array_i2
     module procedure dictionary_add_array_r2
     module procedure dictionary_add_array_i_a
     module procedure dictionary_add_array_r_a
     module procedure dictionary_add_array_c_a
     module procedure dictionary_add_array_l_a
     module procedure dictionary_add_array_s_a
     module procedure dictionary_add_array_i2_a
     module procedure dictionary_add_array_r2_a
  end interface add_array

  !% Remove an entry from a Dictionary
  public remove_value
  interface remove_value
     module procedure dictionary_remove_value
  end interface remove_value

  !% Write a string representation of this dictionary
  public write_string
  interface write_string
     module procedure dictionary_write_string
  end interface write_string

  !% Read into this dictionary from a string
  public read_string
  interface read_string
     module procedure dictionary_read_string
  end interface read_string

  public subset
  interface subset
     module procedure dictionary_subset
     module procedure dictionary_subset_es
  end interface subset

  public swap
  interface swap
     module procedure dictionary_swap
  end interface swap

  public has_key
  interface has_key
     module procedure dictionary_has_key
  end interface has_key

  public bcast
  interface bcast
     module procedure dictionary_bcast
  end interface bcast

  public expand_string
  interface expand_string
     module procedure dictionary_expand_string
  end interface expand_string

  public deepcopy
  interface deepcopy
     module procedure dictionary_deepcopy
  endinterface deepcopy

  public assignment(=)
  interface assignment(=)
     module procedure dictionary_deepcopy_no_error
#ifdef POINTER_COMPONENT_MANUAL_COPY
     module procedure dictentry_assign
     module procedure dictdata_assign
#endif
  end interface assignment(=)

  public :: dictionary_get_key, dictionary_get_type_and_size, dictionary_get_array, lookup_entry_i, lower_case

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
    integer :: i, j

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
       write (line,'("Dict entry string ", A,1x,A)') trim(key), string(this%s)
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
       do i=1,size(this%s_a,1)
          do j=1,size(this%s_a,2)
             call print(this%s_a(i,j),verbosity,file)
          end do
       end do
    else if (this%type == T_INTEGER_A2) then
       call print("Dict entry integer array "//shape(this%i_a2), verbosity,file)
       call print(this%i_a2, verbosity,file)
    else if (this%type == T_REAL_A2) then
       call print("Dict entry real array "//shape(this%r_a2), verbosity,file)
       call print(this%r_a2, verbosity,file)
    else if (this%type == T_DATA .or. this%type == T_DICT) then
       print '("Dict entry arbitary data ",A,1x,I0)', trim(key), size(this%d%d)
       call print(line, verbosity,file)
    endif
  end subroutine dictentry_print

  subroutine dictionary_print_keys(this, verbosity, file)
    type(Dictionary), intent(in) :: this
    integer, intent(in), optional :: verbosity
    type(Inoutput), intent(inout), optional :: file

    integer :: i
    type(extendable_str) :: es

    call initialise(es)
    do i=1,this%N
      if (i < this%N) then
	call concat(es, string(this%keys(i))//":")
      else
	call concat(es, string(this%keys(i)))
      endif
    end do
    call print(es, verbosity=verbosity, file=file)
    call finalise(es)
  end subroutine dictionary_print_keys

  subroutine dictionary_print(this, verbosity, file)
    type(Dictionary), intent(in) :: this
    integer, intent(in), optional :: verbosity
    type(Inoutput), intent(inout), optional :: file
    type(extendable_str) :: str

    call print("Dictionary, allocated "// size(this%entries) //" n_entries " // this%N, verbosity=verbosity,file=file)
    call print('', verbosity=verbosity,file=file)
    call print(write_string(this,entry_sep=quip_new_line),verbosity=verbosity,file=file)
    call finalise(str)

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
    this%cache_invalid = 1
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
!          call print('finalising entry i='//i//' key='//this%keys(i)//' type='//this%entries(i)%type)
          call finalise(this%entries(i))
       end do
       deallocate(this%entries)
    end if
    if (allocated(this%keys)) then
       do i=1,size(this%keys)
          call finalise(this%keys(i))
       end do
       deallocate(this%keys)
    end if
    this%N = 0
    this%cache_invalid = 1

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

    call finalise(this%s)

  end subroutine dictentry_finalise

  subroutine dictionary_get_key(this, i, key, error)
    type(Dictionary), intent(in) :: this
    integer, intent(in) :: i
    character(len=C_KEY_LEN), intent(out) :: key
    integer, intent(out), optional :: error

    INIT_ERROR(error)

    if (i < 1 .or. 1 > this%n) then
       RAISE_ERROR('dictionary_get_key: index '//i//' out of range/', error)
    end if
    key = string(this%keys(i))

  end subroutine dictionary_get_key

  subroutine dictionary_get_type_and_size(this, key, type, thesize, thesize2, error)
    type(Dictionary), intent(in) :: this
    character(len=*), intent(in) :: key
    integer, intent(out) :: type, thesize, thesize2(2)
    integer :: entry_i
    integer, intent(out), optional :: error

    INIT_ERROR(error)
   
    entry_i = lookup_entry_i(this, key)

    if (entry_i <= 0) then
       RAISE_ERROR('dictionary_get_type_and_size: key "'//key//'" not found', error)
    end if

    type = this%entries(entry_i)%type
    thesize = this%entries(entry_i)%len
    thesize2 = this%entries(entry_i)%len2

  end subroutine dictionary_get_type_and_size

  ! ****************************************************************************
  ! *
  ! *    set_value() interface
  ! * 
  ! ****************************************************************************

  subroutine dictionary_set_value_none(this, key)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key

    type(DictEntry) entry
    integer entry_i

    entry%type = T_NONE
    entry_i = add_entry(this, key, entry)
    call finalise(entry)
    
  end subroutine dictionary_set_value_none

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
    call initialise(entry%s)
    call concat(entry%s, value)
    entry_i = add_entry(this, key, entry)
    call finalise(entry)

  end subroutine dictionary_set_value_s

  subroutine dictionary_set_value_s_es(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    type(extendable_str), intent(in) :: value

    type(DictEntry) entry
    integer entry_i

    entry%type = T_CHAR
    call initialise(entry%s, value)
    entry_i = add_entry(this, key, entry)
    call finalise(entry)

  end subroutine dictionary_set_value_s_es

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
    if (do_alloc) then
       allocate(this%entries(entry_i)%i_a(size(value)))
       this%cache_invalid = 1
    end if
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
    if (do_alloc) then
       allocate(this%entries(entry_i)%r_a(size(value)))
       this%cache_invalid = 1
    end if
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
    entry%len = 0
    entry%len2 = shape(value)
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%i_a2(size(value,1),size(value,2)))
       this%cache_invalid = 1
    end if
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
    entry%len = 0
    entry%len2 = shape(value)
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%r_a2(size(value,1),size(value,2)))
       this%cache_invalid = 1
    end if
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
    if (do_alloc) then
       allocate(this%entries(entry_i)%c_a(size(value)))
       this%cache_invalid = 1
    end if
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
    if (do_alloc) then
       allocate(this%entries(entry_i)%l_a(size(value)))
       this%cache_invalid = 1
    end if
    this%entries(entry_i)%l_a = value
    call finalise(entry)
  end subroutine dictionary_set_value_l_a

  subroutine dictionary_set_value_s_a(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value(:)

    type(DictEntry) entry
    integer entry_i, i, j
    logical do_alloc

    entry%type = T_CHAR_A
    entry%len = 0
    entry%len2 = (/len(value(1)), size(value) /)
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%s_a(entry%len2(1), entry%len2(2)))
       this%cache_invalid = 1
    end if
    do i=1,entry%len2(1)
       do j=1,entry%len2(2)
          this%entries(entry_i)%s_a(i,j) = value(j)(i:i)
       end do
    end do
    call finalise(entry)
  end subroutine dictionary_set_value_s_a

  subroutine dictionary_set_value_s_a2(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    character, intent(in) :: value(:,:)

    type(DictEntry) entry
    integer entry_i
    logical do_alloc

    entry%type = T_CHAR_A
    entry%len = 0
    entry%len2 = shape(value)
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%s_a(entry%len2(1), entry%len2(2)))
       this%cache_invalid = 1
    end if
    this%entries(entry_i)%s_a = value
    call finalise(entry)
  end subroutine dictionary_set_value_s_a2

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
    this%cache_invalid = 1

  end subroutine dictionary_set_value_d


  subroutine dictionary_set_value_dict(this, key, value)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    type(Dictionary), intent(in) :: value

    type(DictEntry) entry
    integer entry_i

    entry%type = T_DICT
    allocate(entry%d%d(size(transfer(value,entry%d%d))))
    entry%d%d = transfer(value, entry%d%d)
    entry_i = add_entry(this, key, entry)
    call finalise(entry)
    this%cache_invalid = 1

  end subroutine dictionary_set_value_dict
  

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
    character(len=*), intent(out) :: v
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

    v="" ! fill with blanks first
    if (this%entries(entry_i)%type == T_CHAR) then
       v(1:this%entries(entry_i)%s%len) = string(this%entries(entry_i)%s)
       dictionary_get_value_s = .true.
    else
       dictionary_get_value_s = .false.
    endif
  end function dictionary_get_value_s

  function dictionary_get_value_s_es(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    type(extendable_str), intent(out) :: v
    logical :: dictionary_get_value_s_es
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_s_es = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_CHAR) then
       call initialise(v, this%entries(entry_i)%s)
       dictionary_get_value_s_es = .true.
    else
       dictionary_get_value_s_es = .false.
    endif

  end function dictionary_get_value_s_es


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
    character(len=*), intent(out) :: v(:)
    logical :: dictionary_get_value_s_a
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i, j, k

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_s_a = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_CHAR_A) then
       if (len(v(1)) >= size(this%entries(entry_i)%s_a,1) .and. size(v) >= size(this%entries(entry_i)%s_a,2)) then
          do k=1,this%entries(entry_i)%len2(1)
             do j=1,this%entries(entry_i)%len2(2)
                v(j)(k:k) = this%entries(entry_i)%s_a(k,j)
             end do
          end do
          dictionary_get_value_s_a = .true.
       else
          dictionary_get_value_s_a = .false.
       endif
    else
       dictionary_get_value_s_a = .false.
    endif
  end function dictionary_get_value_s_a

  function dictionary_get_value_s_a2(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    character, intent(out) :: v(:,:)
    logical :: dictionary_get_value_s_a2
    logical, optional :: case_sensitive
    integer, optional :: i

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_s_a2 = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_CHAR_A) then
       if (size(v,1) >= size(this%entries(entry_i)%s_a,1) .and. size(v,2) >= size(this%entries(entry_i)%s_a,2)) then
          v(1:size(this%entries(entry_i)%s_a,1),1:size(this%entries(entry_i)%s_a,2)) = this%entries(entry_i)%s_a
          dictionary_get_value_s_a2 = .true.
       else
          dictionary_get_value_s_a2 = .false.
       endif
    else
       dictionary_get_value_s_a2 = .false.
    endif

  end function dictionary_get_value_s_a2


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

  function dictionary_get_value_dict(this, key, v, case_sensitive, i)
    type(Dictionary), intent(in) :: this
    character(len=*) key
    type(Dictionary), intent(out) :: v
    logical :: dictionary_get_value_dict
    logical, optional :: case_sensitive
    integer, optional :: i

    type(dictionary) :: tmp_dict
    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    if (present(i)) i = entry_i

    if (entry_i <= 0) then
       dictionary_get_value_dict = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_DICT) then
       ! bug: overloaded assignment operator is invoked
       call initialise(tmp_dict)
       tmp_dict = transfer(this%entries(entry_i)%d%d, v)
       v = tmp_dict
       call finalise(tmp_dict)
       dictionary_get_value_dict = .true.
    else
       dictionary_get_value_dict = .false.
    endif
  end function dictionary_get_value_dict



  ! ****************************************************************************
  ! *
  ! *    assign_pointer() interface
  ! * 
  ! ****************************************************************************

  function dictionary_assign_pointer_r0(this, key, v, case_sensitive)
    type(Dictionary), intent(in), target :: this
    character(len=*) key
    real(dp), intent(out), pointer :: v
    logical :: dictionary_assign_pointer_r0
    logical, optional :: case_sensitive

    integer entry_i

    entry_i = lookup_entry_i(this, key, case_sensitive)
    
    if (entry_i <= 0) then
       dictionary_assign_pointer_r0 = .false.
       return
    endif

    if (this%entries(entry_i)%type == T_REAL) then
       v => this%entries(entry_i)%r
       dictionary_assign_pointer_r0 = .true.
    else
       dictionary_assign_pointer_r0 = .false.
    endif
  end function dictionary_assign_pointer_r0

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
    character(1), intent(out), pointer, dimension(:,:) :: v
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
    integer, intent(in), target :: ptr(:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_INTEGER_A
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry%len = size(ptr)
    entry_i = add_entry(this, key, entry)
    this%entries(entry_i)%i_a => ptr
    call finalise(entry)
    this%cache_invalid = 1

  end subroutine dictionary_set_value_pointer_i

  subroutine dictionary_set_value_pointer_r(this, key, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), intent(in), target :: ptr(:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_REAL_A
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry%len = size(ptr)
    entry_i = add_entry(this, key, entry)
    this%entries(entry_i)%r_a => ptr
    call finalise(entry)
    this%cache_invalid = 1

  end subroutine dictionary_set_value_pointer_r

  subroutine dictionary_set_value_pointer_c(this, key, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    complex(dp), intent(in), target :: ptr(:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_COMPLEX_A
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry%len = size(ptr)
    entry_i = add_entry(this, key, entry)
    this%entries(entry_i)%c_a => ptr
    call finalise(entry)
    this%cache_invalid = 1

  end subroutine dictionary_set_value_pointer_c

  subroutine dictionary_set_value_pointer_l(this, key, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in), target :: ptr(:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_LOGICAL_A
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry%len = size(ptr)
    entry_i = add_entry(this, key, entry)
    this%entries(entry_i)%l_a => ptr
    call finalise(entry)
    this%cache_invalid = 1

  end subroutine dictionary_set_value_pointer_l

  subroutine dictionary_set_value_pointer_s(this, key, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(1), intent(in), target :: ptr(:,:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_CHAR_A
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry%len = 0
    entry%len2 = size(ptr)
    entry_i = add_entry(this, key, entry)
    this%entries(entry_i)%s_a => ptr
    call finalise(entry)
    this%cache_invalid = 1

  end subroutine dictionary_set_value_pointer_s

  subroutine dictionary_set_value_pointer_i2(this, key, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in), target :: ptr(:,:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_INTEGER_A2
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry%len = 0
    entry%len2 = shape(ptr)
    entry_i = add_entry(this, key, entry)
    this%entries(entry_i)%i_a2 => ptr
    call finalise(entry)
    this%cache_invalid = 1

  end subroutine dictionary_set_value_pointer_i2

  subroutine dictionary_set_value_pointer_r2(this, key, ptr)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), intent(in), target :: ptr(:,:)

    type(DictEntry) entry
    integer entry_i

    entry%type = T_REAL_A2
    entry%own_data = .false. ! force any possible previous entry with same key to be finalised
    entry%len = 0
    entry%len2 = shape(ptr)
    entry_i = add_entry(this, key, entry)
    this%entries(entry_i)%r_a2 => ptr
    call finalise(entry)
    this%cache_invalid = 1

  end subroutine dictionary_set_value_pointer_r2


  ! ****************************************************************************
  ! *
  ! *    add_array() interface
  ! * 
  ! ****************************************************************************

  subroutine dictionary_add_array_i(this, key, value, len, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value
    integer, intent(in) :: len
    integer, pointer, optional, intent(out) :: ptr(:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_INTEGER_A
    entry%len = len
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%i_a(len))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key)//'" with value '//value, PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%i_a(:) = value
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%i_a

  end subroutine dictionary_add_array_i

  subroutine dictionary_add_array_r(this, key, value, len, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: value
    integer, intent(in) :: len
    real(dp), intent(out), pointer, optional :: ptr(:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_REAL_A
    entry%len = len
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%r_a(len))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key)//'" with value '//value, PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%r_a(:) = value
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%r_a

  end subroutine dictionary_add_array_r

  subroutine dictionary_add_array_c(this, key, value, len, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    complex(dp), intent(in) :: value
    integer, intent(in) :: len
    complex(dp), pointer, optional, intent(out) :: ptr(:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_COMPLEX_A
    entry%len = len
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%c_a(len))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key)//'" with value '//value, PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%c_a(:) = value
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%c_a

  end subroutine dictionary_add_array_c

  subroutine dictionary_add_array_l(this, key, value, len, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in) :: value
    integer, intent(in) :: len
    logical, pointer, optional, intent(out) :: ptr(:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_LOGICAL_A
    entry%len = len
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%l_a(len))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key)//'" with value '//value, PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%l_a(:) = value
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%l_a

  end subroutine dictionary_add_array_l

  subroutine dictionary_add_array_s(this, key, value, len2, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(1), intent(in) :: value
    integer, intent(in) :: len2(2)
    character(1), pointer, optional, intent(out) :: ptr(:,:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_CHAR_A
    entry%len = 0
    entry%len2 = len2
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%s_a(len2(1),len2(2)))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key)//'" with value '//value, PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%s_a(:,:) = value
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%s_a

  end subroutine dictionary_add_array_s

  subroutine dictionary_add_array_i2(this, key, value, len2, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value
    integer, intent(in) :: len2(2)
    integer, pointer, optional, intent(out) :: ptr(:,:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_INTEGER_A2
    entry%len = 0
    entry%len2 = len2
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%i_a2(len2(1),len2(2)))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key)//'" with value '//value, PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%i_a2(:,:) = value
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%i_a2

  end subroutine dictionary_add_array_i2

  subroutine dictionary_add_array_r2(this, key, value, len2, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: value
    integer, intent(in) :: len2(2)
    real(dp), intent(out), pointer, optional :: ptr(:,:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_REAL_A2
    entry%len = 0
    entry%len2 = len2
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%r_a2(len2(1),len2(2)))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key)//'" with value '//value, PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%r_a2(:,:) = value
    call finalise(entry)

    if (present(ptr)) ptr => this%entries(entry_i)%r_a2

  end subroutine dictionary_add_array_r2


  subroutine dictionary_add_array_i_a(this, key, value, len, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value(:)
    integer, intent(in) :: len
    integer, pointer, optional, intent(out) :: ptr(:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_INTEGER_A
    entry%len = len
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%i_a(len))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key), PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%i_a(:) = value
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%i_a

  end subroutine dictionary_add_array_i_a

  subroutine dictionary_add_array_r_a(this, key, value, len, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: value(:)
    integer, intent(in) :: len
    real(dp), intent(out), pointer, optional :: ptr(:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_REAL_A
    entry%len = len
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%r_a(len))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key), PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%r_a(:) = value(:)
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%r_a

  end subroutine dictionary_add_array_r_a

  subroutine dictionary_add_array_c_a(this, key, value, len, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    complex(dp), intent(in) :: value(:)
    integer, intent(in) :: len
    complex(dp), pointer, optional, intent(out) :: ptr(:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_COMPLEX_A
    entry%len = len
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then 
       allocate(this%entries(entry_i)%c_a(len))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key), PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%c_a(:) = value(:)
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%c_a

  end subroutine dictionary_add_array_c_a

  subroutine dictionary_add_array_l_a(this, key, value, len, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in) :: value(:)
    integer, intent(in) :: len
    logical, pointer, optional, intent(out) :: ptr(:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_LOGICAL_A
    entry%len = len
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%l_a(len))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key), PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%l_a(:) = value(:)
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%l_a

  end subroutine dictionary_add_array_l_a

  subroutine dictionary_add_array_s_a(this, key, value, len2, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(1), intent(in) :: value(:,:)
    integer, intent(in) :: len2(2)
    character(1), pointer, optional, intent(out) :: ptr(:,:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_CHAR_A
    entry%len = 0
    entry%len2 = len2
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%s_a(len2(1),len2(2)))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key), PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%s_a(:,:) = value(:,:)
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%s_a

  end subroutine dictionary_add_array_s_a

  subroutine dictionary_add_array_i2_a(this, key, value, len2, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value(:,:)
    integer, intent(in) :: len2(2)
    integer, pointer, optional, intent(out) :: ptr(:,:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_INTEGER_A2
    entry%len = 0
    entry%len2 = len2
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then
       allocate(this%entries(entry_i)%i_a2(len2(1),len2(2)))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key), PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%i_a2(:,:) = value(:,:)
    call finalise(entry)
    if (present(ptr)) ptr => this%entries(entry_i)%i_a2

  end subroutine dictionary_add_array_i2_a

  subroutine dictionary_add_array_r2_a(this, key, value, len2, ptr, overwrite)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: value(:,:)
    integer, intent(in) :: len2(2)
    real(dp), intent(out), pointer, optional :: ptr(:,:)
    logical, optional, intent(in) :: overwrite

    type(DictEntry) entry
    integer entry_i
    logical do_alloc, do_overwrite

    do_overwrite = optional_default(.false., overwrite)
    entry%type = T_REAL_A2
    entry%len = 0
    entry%len2 = len2
    entry_i = add_entry(this, key, entry, do_alloc)
    if (do_alloc) then 
       allocate(this%entries(entry_i)%r_a2(len2(1),len2(2)))
       this%cache_invalid = 1
    end if
    if (do_overwrite .and. .not. do_alloc) call print('WARNING: overwriting array "'//trim(key), PRINT_VERBOSE)
    if (do_alloc .or. do_overwrite) this%entries(entry_i)%r_a2(:,:) = value(:,:)
    call finalise(entry)

    if (present(ptr)) ptr => this%entries(entry_i)%r_a2

  end subroutine dictionary_add_array_r2_a

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
    if (entry_i > 0) call remove_entry(this, entry_i)

  end subroutine dictionary_remove_value

  subroutine dictionary_read_string(this, str, append, error)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: str
    logical, optional, intent(in) :: append !% If true, append to dictionary (default false)
    integer, intent(out), optional :: error

    logical :: do_append
    character(len=DICT_FIELD_LENGTH) :: field
    integer equal_pos
    character(len=DICT_FIELD_LENGTH), dimension(dict_n_fields) :: final_fields
    character(len=C_KEY_LEN) :: key
    character(len=DICT_FIELD_LENGTH) :: value
    integer :: i, num_pairs

    INIT_ERROR(error);
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
          RAISE_ERROR("dictionary_read_string: Malformed field '"//trim(field)//"'", error)
       else
          key = field(1:equal_pos-1)
          value = field(equal_pos+1:len(trim(field)))
       endif
       call print("dictionary_read_string: key='"//trim(key)//"' value='"//trim(value)//"'", PRINT_NERD)

       if (.not. dictionary_parse_value(this,key,value)) then
          RAISE_ERROR('dictionary_read_string: error parsing '//trim(key)//'='//trim(value), error)
       end if

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

    character(len=DICT_FIELD_LENGTH), dimension(dict_n_fields) :: fields
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
    call parse_string(myvalue, ' ', fields, num_fields, error=parse_status)

    if (parse_status /= 0) then
       ! parsing failed, treat as string
       call set_value(this,key,myvalue)
       status = .true.
       return
    else if (num_fields == 0) then
       ! Nothing there, so type is T_NONE
       call set_value(this, key)
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

  function dictionary_write_string(this, real_format, entry_sep, char_a_sep, error)
    type(Dictionary), intent(in) :: this
    character(len=*), optional, intent(in) :: real_format !% Output format for reals, default is 'f9.3'
    character(1), optional, intent(in) :: entry_sep !% Entry seperator, default is single space
    character(1), optional, intent(in) :: char_a_sep !% Output separator for character arrays, default is ','
    type(extendable_str) :: str
    character(len=DICT_STRING_LENGTH) :: dictionary_write_string
    integer, intent(out), optional :: error

    integer :: i, j, k
    character(1) :: my_char_a_sep, my_entry_sep
    character(255) :: my_real_format, tmp_string

    INIT_ERROR(error)
    call initialise(str)

    my_real_format = optional_default('f9.3',real_format)
    my_char_a_sep = optional_default(',',char_a_sep)
    my_entry_sep = optional_default(' ',entry_sep)

    do i=1,this%N

       if (i == 1) then
          call concat(str, string(this%keys(i)))
       else
          call concat(str, my_entry_sep//string(this%keys(i)))
       end if

       if (this%entries(i)%type /= T_NONE) call concat(str, '=')
       
       select case(this%entries(i)%type)

       case(T_NONE)
          ! no value to write

       case(T_INTEGER)
          call concat(str, ''//this%entries(i)%i)

       case(T_REAL)
          write (line,'('//my_real_format//')') this%entries(i)%r
          call concat(str, trim(adjustl(line)))

       case(T_COMPLEX)
          call concat(str, ''//this%entries(i)%c)

       case(T_LOGICAL)
          call concat(str, ''//this%entries(i)%l)

       case(T_CHAR)
          if (index(string(this%entries(i)%s), ' ') == 0) then
             call concat(str, string(this%entries(i)%s))
          else
             call concat(str, '"'//string(this%entries(i)%s)//'"')
          end if

       case(T_INTEGER_A)
          call concat(str, '"'//this%entries(i)%i_a//'"')

       case(T_REAL_A)
          write (line, '('//size(this%entries(i)%r_a)//trim(my_real_format)//')') this%entries(i)%r_a
          call concat(str, '"'//trim(adjustl(line))//'"')

       case(T_COMPLEX_A)
          call concat(str, '"'//this%entries(i)%c_a//'"')

       case(T_LOGICAL_A)
          call concat(str, '"'//this%entries(i)%l_a//'"')

       case(T_CHAR_A)
          call concat(str,'"')
          do j=1,size(this%entries(i)%s_a,2)
             tmp_string = ''
             do k=1,size(this%entries(i)%s_a,1)
                tmp_string(k:k) = this%entries(i)%s_a(k,j)
             end do
             call concat(str, trim(tmp_string)//my_char_a_sep)
          end do
          call concat(str,'"')

       case(T_INTEGER_A2)
          call concat(str, '"('//shape(this%entries(i)%i_a2)//') ')
          call concat(str, ' '//reshape(this%entries(i)%i_a2,(/size(this%entries(i)%i_a2)/))//'"')

       case(T_REAL_A2)
          call concat(str, '"('//shape(this%entries(i)%r_a2)//') ')
          call concat(str, reshape(this%entries(i)%r_a2,(/size(this%entries(i)%r_a2)/))//'"')

       case(T_DATA)
          call concat(str, 'DATA"'//this%entries(i)%d%d//'"')

       case(T_DICT)
          call concat(str, 'DATA"'//this%entries(i)%d%d//'"')
       end select
    end do

    if (str%len > len(dictionary_write_string)) then
       RAISE_ERROR('dictionary_write_string: string too long (' // str%len // ' > ' // len(dictionary_write_string) // ')', error)
    end if

    dictionary_write_string = ''//str
    call finalise(str)

  end function dictionary_write_string


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !% OMIT
  subroutine remove_entry(this, entry_i, error)
    type(Dictionary), intent(inout) :: this
    integer, intent(in) :: entry_i
    integer, intent(out), optional :: error
    integer i

    INIT_ERROR(error)

    if (entry_i <= 0 .or. entry_i > this%N) then
       RAISE_ERROR("remove_entry: Called remove_entry with invalid entry number", error)
    end if

    if (entry_i < this%N) then
       do i=entry_i,this%n-1
          call dictionary_swap(this, string(this%keys(i)), string(this%keys(i+1)))
       end do
    endif

    call finalise(this%keys(this%n))

    this%N = this%N - 1
    this%cache_invalid = 1

  end subroutine remove_entry

  !% OMIT
  function add_entry(this, key, entry, array_alloc) result(entry_i)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    type(DictEntry), intent(in) :: entry
    logical, intent(out), optional :: array_alloc

    integer entry_i

    ! must test separately because size(this%entries) is undefined if allocated(this%entries) is false
    if (.not. allocated(this%entries)) call extend_entries(this, n_entry_block)
    if (this%N >= size(this%entries)) call extend_entries(this, n_entry_block)

    entry_i = lookup_entry_i(this, key)
    if (present(array_alloc)) array_alloc = .true.

    if (entry_i > 0) then
       ! we only free entry if new type/size/shape is incompatible with existing entry
       ! or if new entry will not be allocated by us (i.e. entry%own_data == .false.)
       if (.not. entry%own_data .or. &
            this%entries(entry_i)%type == T_DATA .or. &
            this%entries(entry_i)%type == T_DICT .or. &
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
       call initialise(this%keys(entry_i))
       call concat(this%keys(entry_i),  key)
       this%entries(entry_i) = entry
    endif
  end function add_entry


  !% OMIT
  subroutine extend_entries(this, n)
    type(Dictionary), intent(inout) :: this
    integer, intent(in) :: n

    type(DictEntry), allocatable :: t_entries(:)
    type(Extendable_str), allocatable :: t_keys(:)
    integer old_size
#ifdef ALLOCATABLE_COMPONENT_MANUAL_COPY
    integer :: i
#endif


    if (allocated(this%entries)) then
       allocate(t_entries(size(this%entries)))
       allocate(t_keys(size(this%entries)))
#ifdef ALLOCATABLE_COMPONENT_MANUAL_COPY
       allocate(t_keys(size(this%keys)))
       do i=1,size(this%entries)
	  t_entries(i) = this%entries(i)
	  t_keys(i) = this%keys(i)
       end do
#else
       t_entries = this%entries
       t_keys = this%keys
#endif

       old_size = size(this%entries)

       deallocate(this%entries)
       deallocate(this%keys)

       allocate(this%entries(old_size + n))
       allocate(this%keys(old_size + n))

#ifdef ALLOCATABLE_COMPONENT_MANUAL_COPY
       do i=1, size(t_entries)
	  this%entries(i) = t_entries(i)
	  this%keys(i) = t_keys(i)
       end do
#else
       this%entries(1:this%N) = t_entries(1:this%N)
       this%keys(1:this%N) = t_keys(1:this%N)
#endif

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
          if (string(this%keys(i)) == trim(key)) then
             lookup_entry_i = i
             return
          endif
       else
          if (lower_case(string(this%keys(i))) == trim(lower_case(key))) then
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

  subroutine dictionary_subset(this, keys, out, case_sensitive, out_no_initialise, error)
    type(Dictionary), intent(in) :: this
    character(len=*), dimension(:) :: keys
    type(Dictionary), intent(inout) :: out
    logical, intent(in), optional :: case_sensitive, out_no_initialise
    integer, intent(out), optional :: error

    type(extendable_str), dimension(:), allocatable :: tmp_keys
    integer i
    
    INIT_ERROR(error)
    
    allocate(tmp_keys(size(keys)))
    do i=1,size(keys)
       call initialise(tmp_keys(i))
       call concat(tmp_keys(i), keys(i))
    end do

    call dictionary_subset_es(this, tmp_keys, out, case_sensitive, out_no_initialise, error)
    PASS_ERROR(error)

    do i=1,size(tmp_keys)
       call finalise(tmp_keys(i))
    end do
    deallocate(tmp_keys)

  end subroutine dictionary_subset


  !% Return a dictionary that is a subset of this in 'out'
  !% with only the keys in the arrays 'keys' present.
  subroutine dictionary_subset_es(this, keys, out, case_sensitive, out_no_initialise, error)
    type(Dictionary), intent(in) :: this
    type(extendable_str), dimension(:), intent(in) :: keys
    type(Dictionary), intent(inout) :: out
    logical, intent(in), optional :: case_sensitive, out_no_initialise
    integer, intent(out), optional :: error

    type(Dictionary) :: tmp_dict
    logical :: my_out_no_initialise
    integer :: i, j, io
    
    INIT_ERROR(error)


    my_out_no_initialise = optional_default(.false., out_no_initialise)
    if (.not. my_out_no_initialise) call initialise(out)

    do j=1,size(keys)

       i = lookup_entry_i(this, string(keys(j)), case_sensitive)
       if (i == -1) then
          RAISE_ERROR('dictionary_subset_es: key '//string(keys(j))//' not in dictionary', error)
       end if
       if (my_out_no_initialise) then 
	  io = lookup_entry_i(out, string(keys(j)), case_sensitive)
	  if (io >= 1) then
	     if (this%entries(i)%type /= out%entries(io)%type) then
		RAISE_ERROR('entry type for key '//string(keys(i))//' does not match in this, out', error)
	     endif
	     select case(this%entries(i)%type)
	     case(T_INTEGER_A)
		if (any(shape(this%entries(i)%i_a) /= shape(out%entries(io)%i_a))) then
		   RAISE_ERROR('entry size for key '//string(keys(i))//' does not match in this, out', error)
		endif
	     case(T_REAL_A)
		if (any(shape(this%entries(i)%r_a) /= shape(out%entries(io)%r_a))) then
		   RAISE_ERROR('entry size for key '//string(keys(i))//' does not match in this, out', error)
		endif
	     case(T_COMPLEX_A)
		if (any(shape(this%entries(i)%c_a) /= shape(out%entries(io)%c_a))) then
		   RAISE_ERROR('entry size for key '//string(keys(i))//' does not match in this, out', error)
		endif
	     case(T_LOGICAL_A)
		if (any(shape(this%entries(i)%l_a) /= shape(out%entries(io)%l_a))) then
		   RAISE_ERROR('entry size for key '//string(keys(i))//' does not match in this, out', error)
		endif
	     case(T_CHAR_A)
		if (any(shape(this%entries(i)%s_a) /= shape(out%entries(io)%s_a))) then
		   RAISE_ERROR('entry size for key '//string(keys(i))//' does not match in this, out', error)
		endif
	     case(T_INTEGER_A2)
		if (any(shape(this%entries(i)%i_a2) /= shape(out%entries(io)%i_a2))) then
		   RAISE_ERROR('entry size for key '//string(keys(i))//' does not match in this, out', error)
		endif
	     case(T_REAL_A2)
		if (any(shape(this%entries(i)%r_a2) /= shape(out%entries(io)%r_a2))) then
		   RAISE_ERROR('entry size for key '//string(keys(i))//' does not match in this, out', error)
		endif
	     case(T_DATA)
		if (any(shape(this%entries(i)%d) /= shape(out%entries(io)%d))) then
		   RAISE_ERROR('entry size for key '//string(keys(i))//' does not match in this, out', error)
		endif
	     end select
	     call remove_entry(out, io)
	  endif ! io >= 1
       endif ! my_out_no_initialise

       select case(this%entries(i)%type)
       case(T_INTEGER)
          call set_value(out, string(this%keys(i)), this%entries(i)%i)

       case(T_REAL)
          call set_value(out, string(this%keys(i)), this%entries(i)%r)

       case(T_COMPLEX)
          call set_value(out, string(this%keys(i)), this%entries(i)%c)

       case(T_LOGICAL)
          call set_value(out, string(this%keys(i)), this%entries(i)%l)

       case(T_CHAR)
          call set_value(out, string(this%keys(i)), this%entries(i)%s)

       case(T_INTEGER_A)
          call set_value(out, string(this%keys(i)), this%entries(i)%i_a)

       case(T_REAL_A)
          call set_value(out, string(this%keys(i)), this%entries(i)%r_a)

       case(T_COMPLEX_A)
          call set_value(out, string(this%keys(i)), this%entries(i)%c_a)

       case(T_LOGICAL_A)
          call set_value(out, string(this%keys(i)), this%entries(i)%l_a)

       case(T_CHAR_A)
          call set_value(out, string(this%keys(i)), this%entries(i)%s_a)

       case(T_INTEGER_A2)
          call set_value(out, string(this%keys(i)), this%entries(i)%i_a2)

       case(T_REAL_A2)
          call set_value(out, string(this%keys(i)), this%entries(i)%r_a2)

       case(T_DATA)
          call set_value(out, string(this%keys(i)), this%entries(i)%d)

       case(T_DICT)
          if (.not. get_value(this, string(this%keys(i)), tmp_dict)) then
             RAISE_ERROR('dictionary_subset_es: cannot get_value() as Dictionary type.', error)
          end if
          call set_value(out, string(this%keys(i)), tmp_dict)
          call finalise(tmp_dict)

       end select
    end do

  end subroutine dictionary_subset_es


  !% Swap the positions of two entries in the dictionary. Arrays are not moved in memory.
  subroutine dictionary_swap(this, key1, key2, case_sensitive, error)
    type(Dictionary), intent(inout) :: this
    character(len=*), intent(in) :: key1, key2
    logical, optional :: case_sensitive
    integer, intent(out), optional :: error

    integer :: i1,i2
    type(DictEntry) :: tmp_entry
    type(Extendable_str) :: tmp_key

    INIT_ERROR(error);

    i1 = lookup_entry_i(this,key1,case_sensitive)
    i2 = lookup_entry_i(this,key2,case_sensitive)

    if (i1 <= 0) then
       RAISE_ERROR('dictionary_swap: key '//key1//' not in dictionary', error)
    end if
    if (i2 <= 0) then
       RAISE_ERROR('dictionary_swap: key '//key2//' not in dictionary', error)
    end if

    tmp_entry = this%entries(i2)
    tmp_key   = this%keys(i2)
    this%entries(i2) = this%entries(i1)
    this%keys(i2)    = this%keys(i1)
    this%entries(i1) = tmp_entry
    this%keys(i1)    = tmp_key

    ! Avoid memory leaks by freeing Extanable_str memory
    ! (other entry types are shallow copied).
    call finalise(tmp_key)
    if (tmp_entry%type == T_CHAR) call finalise(tmp_entry%s)

  end subroutine dictionary_swap

  subroutine dictionary_bcast(mpi, dict, error)
    type(MPI_context), intent(in) :: mpi
    type(Dictionary), intent(inout) :: dict
    integer, intent(out), optional :: error

    integer :: i, size_tmp, shape_tmp(2)
    character, allocatable, dimension(:) :: char_array
    type(DictEntry) :: entry

    INIT_ERROR(error)

    if (.not. mpi%active) return

    if (mpi%my_proc == 0) then
       call bcast(mpi, dict%n)
       do i=1,dict%n
          call bcast(mpi, dict%entries(i)%type)
          call bcast(dict%keys(i), mpi%communicator)

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
             call bcast(dict%entries(i)%s, mpi%communicator)
          else if (dict%entries(i)%type == T_INTEGER_A) then
             !size_tmp = size(dict%entries(i)%i_a)
             call bcast(mpi, dict%entries(i)%len)
             call bcast(mpi, dict%entries(i)%i_a)
          else if (dict%entries(i)%type == T_REAL_A) then
             !size_tmp = size(dict%entries(i)%r_a)
             call bcast(mpi, dict%entries(i)%len)
             call bcast(mpi, dict%entries(i)%r_a)
          else if (dict%entries(i)%type == T_COMPLEX_A) then
             !size_tmp = size(dict%entries(i)%c_a)
             call bcast(mpi, dict%entries(i)%len)
             call bcast(mpi, dict%entries(i)%c_a)
          else if (dict%entries(i)%type == T_LOGICAL_A) then
             !size_tmp = size(dict%entries(i)%l_a)
             call bcast(mpi, dict%entries(i)%len)
             call bcast(mpi, dict%entries(i)%l_a)
          else if (dict%entries(i)%type == T_CHAR_A) then
             !shape_tmp = shape(dict%entries(i)%s_a)
             call bcast(mpi, dict%entries(i)%len2)
             call bcast(mpi, dict%entries(i)%s_a)
          else if (dict%entries(i)%type == T_INTEGER_A2) then
             !shape_tmp = shape(dict%entries(i)%i_a2)
             call bcast(mpi, dict%entries(i)%len2)
             call bcast(mpi, dict%entries(i)%i_a2)
          else if (dict%entries(i)%type == T_REAL_A2) then
             !shape_tmp = shape(dict%entries(i)%r_a2)
             call bcast(mpi, dict%entries(i)%len2)
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
          call bcast(dict%keys(i), mpi%communicator)

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
             call bcast(dict%entries(i)%s, mpi%communicator)
          else if (dict%entries(i)%type == T_INTEGER_A) then
             !size_tmp = size(dict%entries(i)%i_a)
             call bcast(mpi, size_tmp)
             dict%entries(i)%len = size_tmp
             dict%entries(i)%len2 = (/0, 0/)
             allocate(dict%entries(i)%i_a(size_tmp))
             call bcast(mpi, dict%entries(i)%i_a)
          else if (dict%entries(i)%type == T_REAL_A) then
             !size_tmp = size(dict%entries(i)%r_a)
             call bcast(mpi, size_tmp)
             dict%entries(i)%len = size_tmp
             dict%entries(i)%len2 = (/0, 0/)
             allocate(dict%entries(i)%r_a(size_tmp))
             call bcast(mpi, dict%entries(i)%r_a)
          else if (dict%entries(i)%type == T_COMPLEX_A) then
             !size_tmp = size(dict%entries(i)%c_a)
             call bcast(mpi, size_tmp)
             dict%entries(i)%len = size_tmp
             dict%entries(i)%len2 = (/0, 0/)
             allocate(dict%entries(i)%c_a(size_tmp))
             call bcast(mpi, dict%entries(i)%c_a)
          else if (dict%entries(i)%type == T_LOGICAL_A) then
             !size_tmp = size(dict%entries(i)%l_a)
             call bcast(mpi, size_tmp)
             dict%entries(i)%len = size_tmp
             dict%entries(i)%len2 = (/0, 0/)
             allocate(dict%entries(i)%l_a(size_tmp))
             call bcast(mpi, dict%entries(i)%l_a)
          else if (dict%entries(i)%type == T_CHAR_A) then
             !shape_tmp = shape(dict%entries(i)%s_a)
             call bcast(mpi, shape_tmp)
             dict%entries(i)%len = 0
             dict%entries(i)%len2 = shape_tmp
             allocate(dict%entries(i)%s_a(shape_tmp(1),shape_tmp(2)))
             call bcast(mpi, dict%entries(i)%s_a)
          else if (dict%entries(i)%type == T_INTEGER_A2) then
             !shape_tmp = shape(dict%entries(i)%i_a2)
             call bcast(mpi, shape_tmp)
             dict%entries(i)%len = 0
             dict%entries(i)%len2 = shape_tmp
             allocate(dict%entries(i)%i_a2(shape_tmp(1),shape_tmp(2)))
             call bcast(mpi, dict%entries(i)%i_a2)
          else if (dict%entries(i)%type == T_REAL_A2) then
             !shape_tmp = shape(dict%entries(i)%r_a2)
             call bcast(mpi, shape_tmp)
             dict%entries(i)%len = 0
             dict%entries(i)%len2 = shape_tmp
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

  !% Copy extendable string from 'in' to 'out', expanding variables formatted
  !% as '$KEY' or '${KEY}' using values in this Dictionary. An error
  !% is raised a key is not found, or if 
  subroutine dictionary_expand_string(this, in, out, error)
    type(Dictionary), intent(in) :: this
    type(Extendable_str), intent(inout) :: in
    type(Extendable_str), intent(out) :: out
    integer, intent(out), optional :: error

    integer s, e, save_cur
    character(len=C_KEY_LEN) :: key
    character(len=DICT_FIELD_LENGTH) :: val
    character(len=63), parameter :: valid_chars =  'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'
    logical got_brace

    INIT_ERROR(error)

    save_cur = in%cur
    in%cur = 1
    call initialise(out)
    do
       s = index(in, '$')
       if (s == 0) then
          call concat(out, substr(in, in%cur, in%len))
          exit
       end if
       if (s > len(in)-1) then
          RAISE_ERROR('dictionary_expand_string: $ found at end of input string "'//in//'"', error)
       end if

       if (in%s(s+1) == '{') then
          got_brace = .true.
          e = index(in, '}')
          if (e == 0) then
             RAISE_ERROR('dictionary_expand_string: unmatched { in input string "'//in//'"', error)
          end if
          key = substr(in, s+2, e-1, error)
          PASS_ERROR(error);
       else
          got_brace = .false.
          e = s+1
          do while (all(verify(in%s(e:e), valid_chars) == 0))
             e = e + 1
             if (e > len(in)) exit
          end do
          e = e - 1
          key = substr(in, s+1, e, error)
          PASS_ERROR(error);
       end if

       if (.not. get_value(this, key, val)) then
          RAISE_ERROR('dictionary_expand_string: unknown key "'//trim(key)//'" or value is not of type T_CHAR', error)
       end if

       call concat(out, substr(in, in%cur, s-1, error))
       PASS_ERROR(error);
       call concat(out, val)
       in%cur = e+1
    end do

    in%cur = save_cur

  end subroutine dictionary_expand_string


  !% Make a deep copy of 'from' in 'this', allocating new memory for array components
  subroutine dictionary_deepcopy(this, from, error)
    type(Dictionary), intent(inout) :: this
    type(Dictionary), intent(in) :: from
    integer, optional, intent(out) :: error
    
    INIT_ERROR(error)
    call subset(from, from%keys(1:from%n), this, error=error)
    PASS_ERROR(error)
    
  end subroutine dictionary_deepcopy


  !% Make a deep copy of 'from' in 'this', allocating new memory for array components; no error passing, for the assignment operator
  subroutine dictionary_deepcopy_no_error(this, from)
    type(Dictionary), intent(inout) :: this
    type(Dictionary), intent(in) :: from
    
    call deepcopy(this, from)
    
  end subroutine dictionary_deepcopy_no_error

  subroutine dictionary_get_array(this, key, nd, dtype, dshape, dloc)
    use iso_c_binding, only: c_intptr_t
    type(Dictionary), intent(in) :: this
    character(len=*), intent(in) :: key
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer(c_intptr_t), intent(out) :: dloc
    
    integer entry_i

    nd = 0
    dtype = 0
    dshape(:) = 0
    dloc = 0

    entry_i = lookup_entry_i(this, key)
    if (entry_i == 0) return

    dtype = this%entries(entry_i)%type

    select case(this%entries(entry_i)%type)
       case(T_INTEGER_A)
          nd = 1
          dshape(1) = size(this%entries(entry_i)%i_a)
          dloc = loc(this%entries(entry_i)%i_a)

       case(T_REAL_A)
          nd = 1
          dshape(1) = size(this%entries(entry_i)%r_a)
          dloc = loc(this%entries(entry_i)%r_a)

       case(T_COMPLEX_A)
          nd = 1
          dshape(1) = size(this%entries(entry_i)%c_a)
          dloc = loc(this%entries(entry_i)%c_a)

       case(T_LOGICAL_A)
          nd = 1
          dshape(1) = size(this%entries(entry_i)%l_a)
          dloc = loc(this%entries(entry_i)%l_a)

       case(T_CHAR_A)
          nd = 2
          dshape(1) = size(this%entries(entry_i)%s_a,1)
          dshape(2) = size(this%entries(entry_i)%s_a,2)
          dloc = loc(this%entries(entry_i)%s_a)

       case(T_INTEGER_A2)
          nd = 2
          dshape(1) = size(this%entries(entry_i)%i_a2, 1)
          dshape(2) = size(this%entries(entry_i)%i_a2, 2)
          dloc = loc(this%entries(entry_i)%i_a2)

       case(T_REAL_A2)
          nd = 2
          dshape(1) = size(this%entries(entry_i)%r_a2, 1)
          dshape(2) = size(this%entries(entry_i)%r_a2, 2)
          dloc = loc(this%entries(entry_i)%r_a2)
    end select

  end subroutine dictionary_get_array

     
#ifdef POINTER_COMPONENT_MANUAL_COPY
subroutine dictentry_assign(to, from)
   type(DictEntry), intent(inout) :: to
   type(DictEntry), intent(in) :: from

   to%type = from%type
   to%len = from%len
   to%len2 = from%len2

   to%own_data = from%own_data
   to%i = from%i
   to%r = from%r
   to%c = from%c

   to%s = from%s

   to%i_a => from%i_a
   to%r_a => from%r_a
   to%c_a => from%c_a
   to%l_a => from%l_a
   to%s_a => from%s_a

   to%i_a2 => from%i_a2
   to%r_a2 => from%r_a2

   to%d = from%d
end subroutine dictentry_assign

subroutine dictdata_assign(to, from)
   type(DictData), intent(inout) :: to
   type(DictData), intent(in) :: from

   if (allocated(to%d)) deallocate(to%d)
   if (allocated(from%d)) then
      allocate(to%d(size(from%d)))
      to%d = from%d
   endif
end subroutine dictdata_assign
#endif

end module dictionary_module
