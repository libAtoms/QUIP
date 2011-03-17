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
!X Exception handling
!X
!% This modules keeps track an error stack. When an error occurs a list
!% of files/lines that allow to trace back the position in the code where
!% the error occured first is constructed.
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module error_module
  implicit none

  private

  ! ---

  integer, parameter  :: ERROR_STACK_SIZE  = 100

  integer, parameter  :: ERROR_DOC_LENGTH = 1000
  integer, parameter  :: ERROR_FN_LENGTH  = 100

  ! ---

  !% Error kinds
  integer, parameter  :: ERROR_NONE                 = 0
  integer, parameter  :: ERROR_UNSPECIFIED          = -1
  integer, parameter  :: ERROR_IO                   = -2
  integer, parameter  :: ERROR_IO_EOF               = -3
  integer, parameter  :: ERROR_MPI                  = -4
  integer, parameter  :: ERROR_MINIM_NOT_CONVERGED  = -5

  !% Strings
  integer, parameter                      :: ERROR_STR_LENGTH = 20
  character(ERROR_STR_LENGTH), parameter  :: ERROR_STR_UNSPECIFIED = &
       "unspecified"
  character(ERROR_STR_LENGTH), parameter  :: ERROR_STR_IO = "IO"
  character(ERROR_STR_LENGTH), parameter  :: ERROR_STR_IO_EOF = "IO EOF"
  character(ERROR_STR_LENGTH), parameter  :: ERROR_STR_MPI = "MPI"
  character(ERROR_STR_LENGTH), parameter  :: ERROR_STR_MINIM_NOT_CONVERGED = &
       "MINIM_NOT_CONVERGED"
  character(ERROR_STR_LENGTH), parameter  :: ERROR_STRINGS(5) = &
       (/ ERROR_STR_UNSPECIFIED, ERROR_STR_IO, ERROR_STR_IO_EOF, &
          ERROR_STR_MPI, ERROR_STR_MINIM_NOT_CONVERGED /)

  public :: ERROR_NONE, ERROR_UNSPECIFIED, ERROR_IO, ERROR_IO_EOF, ERROR_MPI
  public :: ERROR_MINIM_NOT_CONVERGED

  ! ---

  type ErrorDescriptor
     integer                      :: kind     !% kind of error
     logical                      :: has_doc  !% does the documentation string exist?
     character(ERROR_DOC_LENGTH)  :: doc      !% documentation string
     character(ERROR_FN_LENGTH)   :: fn       !% file name where the error occured
     integer                      :: line     !% code line where the error occured
  endtype ErrorDescriptor

  ! ---

  save
  integer                :: error_stack_position  = 0       !% If this is zero, no error has occured
  integer                :: error_mpi_myid = 0              !% MPI rank of process. Initialised by system_initialise()
  type(ErrorDescriptor)  :: error_stack(ERROR_STACK_SIZE)   !% Error stack

  ! ---

  public :: system_abort
  interface system_abort
     module procedure error_abort_with_message
  endinterface

  interface quippy_running
     function quippy_running()
       logical :: quippy_running
     end function quippy_running
  end interface quippy_running

  interface quippy_error_abort
     subroutine quippy_error_abort(message)
       character(*), intent(in) :: message
     end subroutine quippy_error_abort
  end interface quippy_error_abort

  public :: error_abort
  interface error_abort
     module procedure error_abort_from_stack
  endinterface

  public :: error_clear_stack

  public :: push_error, push_error_with_info
  public :: get_error_string_and_clear, clear_error

  ! ---

  public  :: error_unit, error_mpi_myid
  integer :: error_unit = -1

contains

  !% Push a new error callback to the stack
  subroutine push_error(fn, line, kind)
    implicit none

    character(*), intent(in)       :: fn
    integer, intent(in)            :: line
    integer, intent(in), optional  :: kind

    ! ---

    !$omp critical

    error_stack_position  = error_stack_position + 1

    if (error_stack_position > ERROR_STACK_SIZE) then
       call system_abort("Fatal error: Error stack size too small.")
    endif

    if (present(kind)) then
       error_stack(error_stack_position)%kind = kind
    else
       error_stack(error_stack_position)%kind = ERROR_UNSPECIFIED
    endif
    error_stack(error_stack_position)%fn    = fn
    error_stack(error_stack_position)%line  = line

    !$omp end critical

  endsubroutine push_error

  subroutine error_clear_stack()
    error_stack_position = 0
  end subroutine error_clear_stack


  !% Push a new information string onto the error stack
  subroutine push_error_with_info(doc, fn, line, kind)
    implicit none

    character(*), intent(in)       :: doc
    character(*), intent(in)       :: fn
    integer, intent(in)            :: line
    integer, intent(in), optional  :: kind

    ! ---

    !$omp critical

    error_stack_position  = error_stack_position + 1

    if (error_stack_position > ERROR_STACK_SIZE) then
       call system_abort("Fatal error: Error stack size too small.")
    endif

    if (present(kind)) then
       error_stack(error_stack_position)%kind = kind
    else
       error_stack(error_stack_position)%kind = ERROR_UNSPECIFIED
    endif
    error_stack(error_stack_position)%fn = fn
    error_stack(error_stack_position)%line = line
    error_stack(error_stack_position)%has_doc = .true.
    error_stack(error_stack_position)%doc = doc

    !$omp end critical

  endsubroutine push_error_with_info


  !% This error has been handled, clear it.
  subroutine clear_error(error)
    implicit none

    integer, intent(inout) :: error

    ! ---

    error = ERROR_NONE

    error_stack_position = 0
    
  endsubroutine clear_error


  !% Construct a string describing the error.
  function get_error_string_and_clear(error) result(str)
    use iso_c_binding

    implicit none

    integer, intent(inout), optional  :: error

    character(ERROR_DOC_LENGTH)       :: str

    ! ---

    integer        :: i
    character(10)  :: linestr

    ! ---

    if (present(error)) then
       if (-error < lbound(ERROR_STRINGS, 1) .or. &
           -error > ubound(ERROR_STRINGS, 1)) then
          call system_abort("Fatal: error descriptor out of bounds. Did you initialise the error variable?")
       endif

       str = "Traceback (most recent call last - error kind " &
            // trim(ERROR_STRINGS(-error)) // "):"
    else
       str = "Traceback (most recent call last)"
    endif
    do i = error_stack_position, 1, -1

       write (linestr, '(I10)')  error_stack(i)%line

       if (error_stack(i)%has_doc) then
          
          str = trim(str) // C_NEW_LINE // &
               '  File "' // &
               trim(error_stack(i)%fn) // &
               '", line ' // &
               trim(adjustl(linestr)) // &
               C_NEW_LINE // &
               "    " // &
               trim(error_stack(i)%doc)

       else

          str = trim(str) // C_NEW_LINE // &
               '  File "' // &
               trim(error_stack(i)%fn) // &
               '", line ' // &
               trim(adjustl(linestr))

       endif

    enddo

    error_stack_position = 0

    if (present(error)) then
       error  = ERROR_NONE
    endif

  endfunction get_error_string_and_clear


  !% Quit with an error message. Calls 'MPI_Abort' for MPI programs.
  subroutine error_abort_with_message(message)
    character(*),      intent(in) :: message
#ifdef IFORT_TRACEBACK_ON_ABORT
    integer :: j
#endif /* IFORT_TRACEBACK_ON_ABORT */
#ifdef SIGNAL_ON_ABORT
    integer :: status
    integer, parameter :: SIGUSR1 = 30
#endif
#ifdef _MPI
    integer::PRINT_ALWAYS
    include "mpif.h"
#endif

    if (quippy_running()) call quippy_error_abort(message)

#ifdef _MPI
    write(unit=error_unit, fmt='(a,i0," ",a)') 'SYSTEM ABORT: proc=',error_mpi_myid,error_linebreak_string(trim(message),80)
#else
    write(unit=error_unit, fmt='(a," ",a)') 'SYSTEM ABORT:', error_linebreak_string(trim(message),80)
#endif
    call flush(error_unit)

#ifdef _MPI
    call MPI_Abort(MPI_COMM_WORLD, 1, PRINT_ALWAYS)
#endif

#ifdef IFORT_TRACEBACK_ON_ABORT
    ! Cause an integer divide by zero error to persuade
    ! ifort to issue a traceback
    j = 1/0
#endif

#ifdef DUMP_CORE_ON_ABORT
    call fabort()
#else
#ifdef SIGNAL_ON_ABORT
    ! send ourselves a USR1 signal rather than aborting
    call kill(getpid(), SIGUSR1, status)
#else
    stop
#endif
#endif
  end subroutine error_abort_with_message


  !% Stop program execution since this error is not handled properly
  subroutine error_abort_from_stack(error)
    implicit none

    integer, intent(inout), optional :: error

    ! ---

    ! This is for compatibility with quippy, change to error_abort
    call system_abort(get_error_string_and_clear(error))

  endsubroutine error_abort_from_stack

  pure function error_linebreak_string_length(str, line_len) result(length)
    character(len=*), intent(in) :: str
    integer, intent(in) :: line_len
    integer :: length

    length = len_trim(str)+2*len_trim(str)/line_len+3

  end function error_linebreak_string_length

  function error_linebreak_string(str, line_len) result(lb_str)
    character(len=*), intent(in) :: str
    integer, intent(in) :: line_len

    character(len=error_linebreak_string_length(str, line_len)) :: lb_str

    logical :: word_break
    integer :: copy_len, last_space
    character(len=len(lb_str)) :: tmp_str
    character :: quip_new_line

#ifdef NO_F2003_NEW_LINE
    quip_new_line = char(13)
#else
    quip_new_line = new_line(' ')
#endif

    lb_str=""
    tmp_str=trim(str)
    do while (len_trim(tmp_str) > 0)
      copy_len = min(len_trim(tmp_str),line_len)
      if (tmp_str(copy_len:copy_len) /= " ") then
	 last_space=scan(tmp_str(1:copy_len), " ", .true.)
	 if ( last_space > 0 .and. (len_trim(tmp_str(1:copy_len)) - last_space) < 4) then
	    copy_len=last_space
	 endif
      endif

      if (len_trim(lb_str) > 0) then ! we already have some text, add newline before concatenating next line
	lb_str = trim(lb_str)//quip_new_line//trim(tmp_str(1:copy_len))
      else ! just concatenate next line
	lb_str = trim(tmp_str(1:copy_len))
      endif
      ! if we broke in mid word, add "-"
      word_break = .true.
      if (tmp_str(copy_len:copy_len) == " ") then ! we broke right after a space, so no wordbreak
	word_break = .false.
      else ! we broke after a character
	if (copy_len < len_trim(tmp_str)) then ! there's another character after this one, check if it's a space
	  if (tmp_str(copy_len+1:copy_len+1) == " ") then
	    word_break = .false.
	  endif
	else ! we broke after the last character
	  word_break = .false.
	endif
      endif
      if (word_break) lb_str = trim(lb_str)//"-"
      tmp_str(1:copy_len) = ""
      tmp_str=adjustl(tmp_str)
    end do

   end function error_linebreak_string


endmodule error_module
