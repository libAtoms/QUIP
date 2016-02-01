#include "error.inc"

subroutine lattice_abc_to_xyz(cell_lengths, cell_angles, lattice)
  use system_module
  use atoms_module
  implicit none
  real(dp), intent(in) :: cell_lengths(3), cell_angles(3)
  real(dp), intent(out) :: lattice(3,3)

  lattice = make_lattice(cell_lengths(1), cell_lengths(2), cell_lengths(3), &
       cell_angles(1), cell_angles(2), cell_angles(3))

end subroutine lattice_abc_to_xyz

subroutine lattice_xyz_to_abc(lattice, cell_lengths, cell_angles)
  use system_module
  use atoms_module
  implicit none
  real(dp), intent(in) :: lattice(3,3)
  real(dp), intent(out) :: cell_lengths(3), cell_angles(3)

  call get_lattice_params(lattice, cell_lengths(1), cell_lengths(2), cell_lengths(3), &
       cell_angles(1), cell_angles(2), cell_angles(3))

end subroutine lattice_xyz_to_abc

subroutine c_system_initialise(verbosity)
  use system_module
  integer verbosity

  call system_initialise(verbosity)

end subroutine c_system_initialise

function quippy_running()
  use system_module, only: get_quippy_running
  logical quippy_running
  quippy_running = get_quippy_running()
end function quippy_running

! Error handling routines callable from C

subroutine c_push_error_with_info(doc, fn, line, kind)
  use error_module
  character(*), intent(in)       :: doc
  character(*), intent(in)       :: fn
  integer, intent(in)            :: line
  integer, intent(in), optional  :: kind

  call push_error_with_info(doc, fn, line, kind)

end subroutine c_push_error_with_info

subroutine c_push_error(fn, line, kind)
  use error_module

  character(*), intent(in)       :: fn
  integer, intent(in)            :: line
  integer, intent(in), optional  :: kind

  call push_error(fn, line, kind)

end subroutine c_push_error

subroutine c_error_abort(error)
  use error_module
  integer, intent(inout), optional :: error

  call error_abort(error)
end subroutine c_error_abort

subroutine c_error_clear_stack()
  use error_module
  call error_clear_stack
end subroutine c_error_clear_stack


! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
! Dictionary routines callable from C
!
! Routines used by CInOutput to manipulate Fortran Dictionary objects from C.
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine c_dictionary_initialise(this)
  use dictionary_module
  integer, intent(inout) :: this(12)
  type(c_dictionary_ptr_type) :: this_ptr

  allocate(this_ptr%p)
  call initialise(this=this_ptr%p)
  this = transfer(this_ptr, this)

end subroutine c_dictionary_initialise

subroutine c_dictionary_finalise(this)
  use dictionary_module
  integer, intent(in) :: this(12)
  type(c_dictionary_ptr_type) :: this_ptr

  this_ptr = transfer(this, this_ptr)
  call finalise(this=this_ptr%p)
  deallocate(this_ptr%p)

end subroutine c_dictionary_finalise

subroutine c_dictionary_get_n(this, n)
  use dictionary_module
  integer, intent(in) :: this(12)
  integer, intent(out) :: n
  type(c_dictionary_ptr_type) :: this_ptr

  this_ptr = transfer(this, this_ptr)
  n = this_ptr%p%n

end subroutine c_dictionary_get_n

subroutine c_dictionary_get_key(this, i, key, length, error)
  use error_module
  use dictionary_module
  integer, intent(in) :: this(12), i
  character(len=*), intent(out) :: key
  integer, intent(out), optional :: length, error
  type(c_dictionary_ptr_type) :: this_ptr

  INIT_ERROR(error)
  length = 0
  this_ptr = transfer(this, this_ptr)
  call dictionary_get_key(this_ptr%p, i, key, error)
  PASS_ERROR(error)
  length = len_trim(key)

end subroutine c_dictionary_get_key

subroutine c_dictionary_query_key(this, key, type, dshape, dloc, error)

  use error_module
  use dictionary_module
  use iso_c_binding, only: c_intptr_t

  use extendable_str_module
  use system_module

#ifdef __GFORTRAN__
  INTERFACE
     subroutine c_dictionary_query_index(this, entry_i, key, type, dshape, dloc, error)
       use iso_c_binding, only: c_intptr_t
       integer, intent(in) :: this(12)
       integer, intent(in) :: entry_i
       character(len=*), intent(out) :: key
       integer, intent(out) :: type
       integer, dimension(2), intent(out) :: dshape
       integer(c_intptr_t), intent(out) :: dloc
       integer, intent(out), optional :: error
     end subroutine c_dictionary_query_index
  END INTERFACE
#endif

  integer, intent(in) :: this(12)
  character(len=*), intent(inout) :: key
  integer, intent(out) :: type
  integer, dimension(2), intent(out) :: dshape
  integer(c_intptr_t), intent(out) :: dloc
  integer, intent(out), optional :: error

  integer entry_i
  type(c_dictionary_ptr_type) :: this_ptr

  INIT_ERROR(error)
  this_ptr = transfer(this, this_ptr)

  entry_i = lookup_entry_i(this_ptr%p, key)
  if (entry_i <= 0) then
     RAISE_ERROR('c_dictionary_query_key: key '//trim(key)//' not found', error)
  end if
  call c_dictionary_query_index(this, entry_i, key, type, dshape, dloc, error)
  PASS_ERROR(error)

end subroutine c_dictionary_query_key

subroutine c_dictionary_query_index(this, entry_i, key, type, dshape, dloc, error)
  use error_module
  use extendable_str_module
  use system_module
  use dictionary_module
  use iso_c_binding, only: c_intptr_t
  integer, intent(in) :: this(12)
  integer, intent(in) :: entry_i
  character(len=*), intent(out) :: key
  integer, intent(out) :: type
  integer, dimension(2), intent(out) :: dshape
  integer(c_intptr_t), intent(out) :: dloc
  integer, intent(out), optional :: error

  type(c_dictionary_ptr_type) :: this_ptr

  INIT_ERROR(error)
  this_ptr = transfer(this, this_ptr)

  type = 0
  dshape(:) = 0
  dloc = 0

  if (entry_i <= 0 .or. entry_i > this_ptr%p%N) then
     RAISE_ERROR('c_dictionary_query_index: entry_i = '//entry_i//' outside range 1 < entry_i <= '//this_ptr%p%N, error)
  end if

  type = this_ptr%p%entries(entry_i)%type
  key  = string(this_ptr%p%keys(entry_i))

  select case(this_ptr%p%entries(entry_i)%type)
  case(T_NONE)
     dloc = 0
     return

  case(T_INTEGER)
     dloc = loc(this_ptr%p%entries(entry_i)%i)

  case(T_REAL)
     dloc = loc(this_ptr%p%entries(entry_i)%r)

  case(T_COMPLEX)
     dloc = loc(this_ptr%p%entries(entry_i)%c)

  case(T_CHAR)
     dshape(1) = this_ptr%p%entries(entry_i)%s%len
     dloc = loc(this_ptr%p%entries(entry_i)%s%s)

  case(T_LOGICAL)
     dloc = loc(this_ptr%p%entries(entry_i)%l)

  case(T_INTEGER_A)
     dshape(1) = size(this_ptr%p%entries(entry_i)%i_a)
     dloc = loc(this_ptr%p%entries(entry_i)%i_a)

  case(T_REAL_A)
     dshape(1) = size(this_ptr%p%entries(entry_i)%r_a)
     dloc = loc(this_ptr%p%entries(entry_i)%r_a)

  case(T_COMPLEX_A)
     dshape(1) = size(this_ptr%p%entries(entry_i)%c_a)
     dloc = loc(this_ptr%p%entries(entry_i)%c_a)

  case(T_LOGICAL_A)
     dshape(1) = size(this_ptr%p%entries(entry_i)%l_a)
     dloc = loc(this_ptr%p%entries(entry_i)%l_a)

  case(T_CHAR_A)
     dshape(1) = size(this_ptr%p%entries(entry_i)%s_a,1)
     dshape(2) = size(this_ptr%p%entries(entry_i)%s_a,2)
     dloc = loc(this_ptr%p%entries(entry_i)%s_a)

  case(T_INTEGER_A2)
     dshape(1) = size(this_ptr%p%entries(entry_i)%i_a2, 1)
     dshape(2) = size(this_ptr%p%entries(entry_i)%i_a2, 2)
     dloc = loc(this_ptr%p%entries(entry_i)%i_a2)

  case(T_REAL_A2)
     dshape(1) = size(this_ptr%p%entries(entry_i)%r_a2, 1)
     dshape(2) = size(this_ptr%p%entries(entry_i)%r_a2, 2)
     dloc = loc(this_ptr%p%entries(entry_i)%r_a2)

  case(T_DATA)
     ! not supported
     dloc = 0
     RAISE_ERROR('c_dictionary_query_index: data type T_DATA not supported.', error)

  end select

end subroutine c_dictionary_query_index

subroutine c_dictionary_add_key(this, key, type, dshape, loc, error)
  use system_module
  use dictionary_module
  use iso_c_binding, only: c_intptr_t

#ifdef __GFORTRAN__
  INTERFACE
     subroutine c_dictionary_query_key(this, key, type, dshape, dloc, error)

       use error_module
       use dictionary_module
       use iso_c_binding, only: c_intptr_t

       use extendable_str_module
       use system_module

       integer, intent(in) :: this(12)
       character(len=*), intent(inout) :: key
       integer, intent(out) :: type
       integer, dimension(2), intent(out) :: dshape
       integer(c_intptr_t), intent(out) :: dloc
       integer, intent(out), optional :: error

     end subroutine c_dictionary_query_key
  END INTERFACE
#endif

  integer, intent(in) :: this(12)
  character(len=*), intent(inout) :: key
  integer, intent(in) :: type
  integer, intent(in) :: dshape(2)
  integer(c_intptr_t), intent(out) :: loc
  integer, intent(out), optional :: error

  type(c_dictionary_ptr_type) this_ptr
  integer mytype, mydshape(2), tmp_error

  INIT_ERROR(error)
  mytype = -1
  mydshape(:) = -1

  ! check if a compatiable entry already exists
  tmp_error = 0
  call c_dictionary_query_key(this, key, mytype, mydshape, loc, tmp_error)
  CLEAR_ERROR(tmp_error)
  if (tmp_error == 0 .and. type == mytype .and. all(dshape == mydshape)) then
     return
  end if

  ! otherwise we need to add a new entry
  this_ptr = transfer(this, this_ptr)
  select case(type)
  case(T_NONE)
     call set_value(this_ptr%p, key)

  case(T_INTEGER)
     call set_value(this_ptr%p, key, 0)

  case(T_REAL)
     call set_value(this_ptr%p, key, 0.0_dp)

  case(T_COMPLEX)
     call set_value(this_ptr%p, key, (0.0_dp, 0.0_dp))

  case(T_CHAR)
     call set_value(this_ptr%p, key, repeat('X', dshape(1)))

  case(T_LOGICAL)
     call set_value(this_ptr%p, key, .false.)

  case(T_INTEGER_A)
     call add_array(this_ptr%p, key, 0, dshape(1))

  case(T_REAL_A)
     call add_array(this_ptr%p, key, 0.0_dp, dshape(1))

  case(T_COMPLEX_A)
     call add_array(this_ptr%p, key, 0.0_dp, dshape(1))

  case(T_LOGICAL_A)
     call add_array(this_ptr%p, key, .false., dshape(1))

  case(T_CHAR_A)
     call add_array(this_ptr%p, key, ' ', dshape)

  case(T_INTEGER_A2)
     call add_array(this_ptr%p, key, 0, dshape)

  case(T_REAL_A2)
     call add_array(this_ptr%p, key, 0.0_dp, dshape)

  case default
     RAISE_ERROR('c_dictionary_add_key: unknown data type '//type, error)

  end select

  ! Finally, update 'loc' pointer
  call c_dictionary_query_key(this, key, mytype, mydshape, loc, error)
  PASS_ERROR(error)

end subroutine c_dictionary_add_key

subroutine c_extendable_str_concat(this, str, keep_lf, add_lf_if_missing)
  use extendable_str_module
  type c_extendable_str_ptr_type
     type(Extendable_str), pointer :: p
  end type c_extendable_str_ptr_type
  integer, intent(in) :: this(12)
  character(*), intent(in) :: str
  integer, intent(in) :: keep_lf, add_lf_if_missing

  type(c_extendable_str_ptr_type) :: this_ptr
  logical do_keep_lf, do_add_lf_if_missing

  do_keep_lf = .false.
  if (keep_lf == 1) do_keep_lf = .true.
  do_add_lf_if_missing = .false.
  if (add_lf_if_missing == 1) do_add_lf_if_missing = .true.

  this_ptr = transfer(this, this_ptr)
  call concat(this_ptr%p, str, do_keep_lf, do_add_lf_if_missing)

end subroutine c_extendable_str_concat

