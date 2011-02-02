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
!X  Extendable string module
!X  
!%  Defines strings that can extend as required, using a character array
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module extendable_str_module

use system_module

implicit none
private

integer, parameter :: EXTENDABLE_STRING_LENGTH_INCREMENT = 10240 !% Increment of extendable string
integer, parameter :: extendable_string_reading_buffer = 1024 !% Read ext_string by chunks of this size

public :: extendable_str
type extendable_str
  character(len=1), allocatable :: s(:)
  integer :: len = 0
  integer :: increment = EXTENDABLE_STRING_LENGTH_INCREMENT
  integer :: cur = 0
end type extendable_str


public :: concat
interface concat
  module procedure extendable_str_concat
end interface concat

public :: string
interface string
  module procedure extendable_str_string
end interface string

public :: substr
interface substr
   module procedure extendable_str_substr
end interface substr

public :: read
interface read
  module procedure extendable_str_read_unit, extendable_str_read_file
end interface read

public :: initialise
interface initialise
  module procedure extendable_str_initialise
end interface initialise

public :: zero
interface zero
  module procedure extendable_str_zero
end interface zero

public :: finalise
interface finalise
  module procedure extendable_str_finalise
end interface finalise

public :: print
interface print
  module procedure extendable_str_print
end interface print

public :: len
interface len
   module procedure extendable_str_len
endinterface

public :: read_line
interface read_line
  module procedure extendable_str_read_line
end interface read_line

public :: parse_line
interface parse_line
  module procedure extendable_str_parse_line
end interface parse_line

public :: index
interface index
  module procedure extendable_str_index
end interface index

public :: bcast
interface bcast
   module procedure extendable_str_bcast
end interface bcast

public :: operator(//)
interface operator(//)
   module procedure extendable_str_cat_string, string_cat_extendable_str
   module procedure extendable_str_cat_extendable_str
   module procedure string_cat_extendable_str_array
end interface operator(//)

public :: assignment(=)
interface assignment(=)
   module procedure extendable_str_assign_string
   module procedure string_assign_extendable_str
#ifdef ALLOCATABLE_COMPONENT_MANUAL_COPY
   module procedure extendable_str_assign_extendable_str
#endif
endinterface

contains

subroutine extendable_str_initialise(this, copy_from)
  type(extendable_str), intent(inout) :: this
  type(extendable_str), intent(in), optional :: copy_from

  call finalise(this)

  if (present(copy_from)) then
     this%len = copy_from%len
     this%increment = copy_from%increment
     this%cur = copy_from%cur
     allocate(this%s(this%len))
     this%s = copy_from%s
  endif
end subroutine extendable_str_initialise

subroutine extendable_str_finalise(this)
  type(extendable_str), intent(inout) :: this

  if (allocated(this%s)) deallocate(this%s)
  this%len = 0
end subroutine extendable_str_finalise

subroutine extendable_str_zero(this)
  type(extendable_str), intent(inout) :: this

  this%len = 0
end subroutine extendable_str_zero

subroutine extendable_str_print(this,verbosity,file)
  use iso_c_binding
  implicit none
  type(extendable_str),    intent(in)           :: this
  integer,        intent(in), optional :: verbosity
  type(Inoutput), intent(inout),optional:: file

  if (allocated(this%s)) then
     call print ("extendable_str, len " // this%len // " size(s) " // size(this%s), verbosity, file)
     call print(this%s(1:this%len))
! OLD IMPLEMENTATION
!    do i=1, this%len
!      call print (this%s(i), verbosity, out)
!    end do
  else
     call print("extendable_str, len " // this%len // " not allocated", verbosity, file)
  endif
end subroutine extendable_str_print

pure function extendable_str_len(this)
  type(extendable_str), intent(in) :: this
  integer :: extendable_str_len
  extendable_str_len = this%len
end function extendable_str_len

subroutine extendable_str_concat(this, str, keep_lf, add_lf_if_missing)
  type(extendable_str), intent(inout) :: this
  character(len=*), intent(in) :: str
  logical, intent(in), optional :: keep_lf, add_lf_if_missing

  character, allocatable :: t(:)
  integer str_len
  integer add_len, new_len
  integer i
  logical my_keep_lf, my_add_lf_if_missing
  logical :: add_lf

  my_keep_lf = optional_default(.true., keep_lf)
  my_add_lf_if_missing = optional_default(.false., add_lf_if_missing)

  str_len = len_trim(str)
  add_lf = .false.
  if (str_len > 0) then
     if (my_add_lf_if_missing .and. str(str_len:str_len) /= quip_new_line) add_lf = .true.
  else
     add_lf = my_add_lf_if_missing
  endif
  add_len = str_len
  if (add_lf) add_len = add_len + 1
  if (add_len > 0) then
    if (allocated(this%s)) then ! already allocated contents
      new_len = this%len + add_len
      if (new_len > size(this%s)) then ! new length too big for current allocated array
	 if (this%increment > 0) then
           new_len = size(this%s)
           do while(new_len < this%len + add_len)
             new_len = new_len + this%increment
             this%increment = this%increment*2
           enddo
	 endif
	 if (this%len > 0) then ! need to save old data
	   allocate(t(size(this%s)))
	   t = this%s
	   deallocate(this%s)
	   allocate(this%s(new_len))
	   this%s(1:this%len) = t(1:this%len)
	   deallocate(t)
	 else ! don't need to save old data
	   if (this%increment > 0) then
             new_len = this%increment
             this%increment = this%increment*2
             do while (new_len < this%len + add_len)
               new_len = new_len + this%increment
               this%increment = this%increment*2
             enddo
           endif
	   deallocate(this%s)
	   allocate(this%s(new_len))
	 endif
      endif
    else ! array not already allocated
      allocate(this%s(add_len))
      this%len = 0
    endif

    do i=1, str_len
       if (.not. my_keep_lf .and. str(i:i) == quip_new_line) cycle
       this%s(this%len+1) = str(i:i)
       this%len = this%len + 1
    end do

    if (add_lf) then
       this%s(this%len+1) = quip_new_line
       this%len = this%len + 1
    endif

 endif

 if (this%len > 0 .and. this%cur <= 0) this%cur = 1

end subroutine extendable_str_concat

function extendable_str_string(this)
  type(extendable_str), intent(in) :: this
  character(len=this%len) :: extendable_str_string
  integer i

  do i=1, this%len
    extendable_str_string(i:i) = this%s(i)
  end do

end function extendable_str_string

function extendable_str_substr(this, start, end, error)
  type(extendable_str), intent(in) :: this
  integer, intent(in) :: start, end
  integer, optional, intent(out) :: error
  character(len=(end-start+1)) :: extendable_str_substr
  integer i, j

  INIT_ERROR(error)
  
  if (start < 1) then
     RAISE_ERROR('extendable_str_substr: start('//start//') < 1', error)
  end if
  if (end > this%len) then
     RAISE_ERROR('extendable_str_substr: end('//end//') > len('//this%len//')', error)
  end if

  j = 1
  do i=start, end
    extendable_str_substr(j:j) = this%s(i)
    j = j + 1
  end do

end function extendable_str_substr


subroutine extendable_str_read_file(this, file, convert_to_string, mpi_comm, keep_lf)
  type(extendable_str), intent(inout) :: this
  character(len=*), intent(in) :: file
  logical, intent(in), optional :: convert_to_string
  integer, intent(in), optional :: mpi_comm
  logical, intent(in), optional :: keep_lf

  type(inoutput) :: in
  logical do_read

  do_read = .true.
  if (present(mpi_comm) .and. mpi_id() /= 0) do_read = .false.

  if (do_read) call initialise(in, trim(file), INPUT)
  call read(this, in%unit, convert_to_string, mpi_comm, keep_lf)
  if (do_read) call finalise(in)

end subroutine extendable_str_read_file

#ifdef PGI_IOSTAT_FUNCS
function is_iostat_end(stat)
  integer, intent(in) :: stat
  logical :: is_iostat_end

  if (stat == -1) then
    is_iostat_end = .true.
  else
    is_iostat_end = .false.
  endif
end function is_iostat_end

function is_iostat_eor(stat)
  integer, intent(in) :: stat
  logical :: is_iostat_eor

  if (stat == -2) then
    is_iostat_eor = .true.
  else
    is_iostat_eor = .false.
  endif
end function is_iostat_eor

#endif

subroutine extendable_str_read_unit(this, unit, convert_to_string, mpi_comm, keep_lf)
  type(extendable_str), intent(inout) :: this
  integer, intent(in) :: unit
  logical, intent(in), optional :: convert_to_string
  integer, intent(in), optional :: mpi_comm
  logical, intent(in), optional :: keep_lf

  character(len=EXTENDABLE_STRING_READING_BUFFER) :: line
  integer n_read
  integer stat
  logical last_was_incomplete
  logical done
  logical my_convert_to_string
  integer stack_size, stack_size_err
  logical do_read
  logical my_keep_lf

  this%len = 0
  my_convert_to_string = optional_default(.false., convert_to_string)
  my_keep_lf = optional_default(.false., keep_lf)

  do_read = .true.
  if (present(mpi_comm) .and. mpi_id() /= 0) do_read = .false.

  if (do_read) then
    done = .false.
    last_was_incomplete = .true.
    do while (.not. done)
      read (unit=unit, fmt='(A)', iostat=stat, advance='no', size=n_read) line
      if (.not. is_iostat_end(stat)) then
	if (n_read == 0) cycle
	if (.not.last_was_incomplete .and. .not. my_keep_lf) then
	  call concat(this, " " // trim(line))
	else
	  call concat(this, trim(line))
	endif
	if (is_iostat_eor(stat) .and. my_keep_lf) then
	  call concat(this, quip_new_line)
	endif
	last_was_incomplete = (stat == 0)
      else
	done = .true. ! EOF
      endif
    end do
  endif

  call extendable_str_bcast(this, mpi_comm)

  if (my_convert_to_string) then
    stack_size = floor(this%len/1024.0_dp) + 10
    stack_size_err = increase_stack(stack_size)
    if (stack_size_err /= 0) then
      call print("extendable_str_read_unit: error calling c_increase_stack size = " // stack_size // &
	" err = "// stack_size_err)
    endif
  endif

  if (this%len > 0) this%cur = 1

end subroutine extendable_str_read_unit

subroutine extendable_str_bcast(this, mpi_comm)
  type(extendable_str), intent(inout) :: this
  integer, intent(in), optional :: mpi_comm

#ifdef _MPI
  include 'mpif.h'
  integer err, size_this_s

  if (present(mpi_comm)) then
    if (mpi_id() == 0) then
      call mpi_bcast(size(this%s), 1, MPI_INTEGER, 0, mpi_comm, err)
      call mpi_bcast(this%len, 1, MPI_INTEGER, 0, mpi_comm, err)
      call mpi_bcast(this%increment, 1, MPI_INTEGER, 0, mpi_comm, err)
      call mpi_bcast(this%s, this%len, MPI_CHARACTER, 0, mpi_comm, err)
    else
      call finalise(this)
      call mpi_bcast(size_this_s, 1, MPI_INTEGER, 0, mpi_comm, err)
      call mpi_bcast(this%len, 1, MPI_INTEGER, 0, mpi_comm, err)
      call mpi_bcast(this%increment, 1, MPI_INTEGER, 0, mpi_comm, err)
      allocate(this%s(size_this_s))
      call mpi_bcast(this%s, this%len, MPI_CHARACTER, 0, mpi_comm, err)
    endif
  endif
#endif

  return
end subroutine extendable_str_bcast

pure function extendable_str_index(this, substr)
  type(extendable_str), intent(in) :: this
  character(len=*), intent(in) :: substr
  integer extendable_str_index

  logical found_it
  integer i, j

  if (this%cur <= 0 .or. size(this%s) <= 0) return

  extendable_str_index = 0

  do i=this%cur, this%len-len(substr)+1
    found_it = .true.
    do j=0, len(substr)-1
      if (this%s(i+j) /= substr(j+1:j+1)) found_it = .false.
      if (.not. found_it) exit
    end do
    if (found_it) then
      extendable_str_index = i
      return
    endif
  end do

end function extendable_str_index

function extendable_str_read_line(this, status)
  type(extendable_str), intent(inout) :: this
  integer, intent(out), optional :: status
  character(len=max(1,index(this,quip_new_line)-this%cur)) :: extendable_str_read_line

  integer line_len
  integer i

  line_len = index(this,quip_new_line)-this%cur

  if (this%cur <= this%len) then
    do i=1, line_len
      extendable_str_read_line(i:i) = this%s(this%cur+i-1)
    end do
    this%cur = this%cur + max(line_len+1,0)
  endif

  if (present(status)) then
    if (line_len < 1) then
      status = -1
    else
      status = 0
    endif
  endif

end function extendable_str_read_line

subroutine extendable_str_parse_line(this, delimiters, fields, num_fields, status)
  type(extendable_str), intent(inout) :: this
  character(*),               intent(in)    :: delimiters
  character(*), dimension(:), intent(inout) :: fields
  integer,                    intent(out)   :: num_fields
  integer, optional,          intent(out)   :: status

  integer my_status
  character(len=EXTENDABLE_STRING_LENGTH_INCREMENT) :: local_line

  local_line = read_line(this, my_status)
  if (present(status)) status = my_status
  if (my_status == 0) then
    call parse_string(local_line, delimiters, fields, num_fields)
  endif
end subroutine extendable_str_parse_line

function extendable_str_cat_string(this, str)
  type(extendable_str), intent(in)  :: this
  character(*), intent(in)          :: str

  type(extendable_str)              :: extendable_str_cat_string

  ! ---

  call initialise(extendable_str_cat_string, this)
  call concat(extendable_str_cat_string, str)
end function extendable_str_cat_string

function string_cat_extendable_str(str, this)
  character(*), intent(in)          :: str
  type(extendable_str), intent(in)  :: this

  character(len(str)+this%len)      :: string_cat_extendable_str

  ! ---

  string_cat_extendable_str = str
  string_cat_extendable_str(max(1,len(str)+1):) = string(this)
end function string_cat_extendable_str

pure function sumlen(this)
  type(extendable_str), intent(in)  :: this(:)
  integer :: sumlen, i
  sumlen = 0
  do i = lbound(this, 1), ubound(this, 1)
     sumlen = sumlen + len(this(i))
  enddo
endfunction sumlen

function string_cat_extendable_str_array(str, this)
  character(*), intent(in)          :: str
  type(extendable_str), intent(in)  :: this(:)
  character(len(str)+sumlen(this)+3*size(this)) :: string_cat_extendable_str_array
  integer :: i, c

  string_cat_extendable_str_array = str
  c = max(1,len(str)+1)
  do i = lbound(this, 1), ubound(this, 1)
     call print("i = " // i // ", c = " // c // ", len(str) = " // len(this(i)) // ", str = " // this(i))
     string_cat_extendable_str_array(c:c+len(this(i))+2) = "'" // string(this(i)) // "'"
     c = c+len(this(i))+3
  enddo
end function string_cat_extendable_str_array

function extendable_str_cat_extendable_str(this, str)
  type(extendable_str), intent(in)  :: this
  type(extendable_str), intent(in)  :: str

  type(extendable_str)              :: extendable_str_cat_extendable_str

  ! ---

  call initialise(extendable_str_cat_extendable_str, this)
  call concat(extendable_str_cat_extendable_str, string(str))
end function extendable_str_cat_extendable_str

subroutine extendable_str_assign_string(to, from)
  type(extendable_str), intent(out)  :: to
  character(*), intent(in)           :: from

  call initialise(to)
  call concat(to, from)
end subroutine extendable_str_assign_string

#ifdef ALLOCATABLE_COMPONENT_MANUAL_COPY
subroutine extendable_str_assign_extendable_str(to, from)
  type(extendable_str), intent(out)  :: to
  type(extendable_str), intent(in)   :: from

  call initialise(to)
  call concat(to, string(from))
end subroutine extendable_str_assign_extendable_str
#endif

subroutine string_assign_extendable_str(to, from)
  type(extendable_str), intent(in)  :: from
  character(from%len), intent(out)  :: to

  to = string(from)
end subroutine string_assign_extendable_str


end module extendable_str_module
