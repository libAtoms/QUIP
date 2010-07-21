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

  ! ---

  integer, parameter  :: ERROR_STACK_SIZE  = 100

  integer, parameter  :: ERROR_DOC_LENGTH = 1000
  integer, parameter  :: ERROR_FN_LENGTH  = 100

  ! ---

  !% Error kinds
  integer, parameter  :: ERROR_NONE         = 0
  integer, parameter  :: ERROR_UNSPECIFIED  = -1
  integer, parameter  :: ERROR_IO           = -2
  integer, parameter  :: ERROR_MPI          = -3

  public :: ERROR_NONE, ERROR_UNSPECIFIED, ERROR_IO, ERROR_MPI

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
  type(ErrorDescriptor)  :: error_stack(ERROR_STACK_SIZE)   !% Error stack

  ! ---

  public :: system_abort
#ifdef HAVE_QUIPPY
  interface system_abort
     subroutine system_abort(message)
       character(*), intent(in) :: message
     end subroutine system_abort
  end interface system_abort

  public :: error_abort
  interface error_abort
     module procedure error_abort_from_stack
  endinterface
#else
  interface system_abort
     module procedure error_abort_with_message
  endinterface

  public :: error_abort
  interface error_abort
     module procedure error_abort_with_message, error_abort_from_stack
  endinterface
#endif

  public :: push_error, push_error_with_info
  public :: get_error_string_and_clear, clear_error

  ! ---

  public  :: error_unit
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
       stop "Fatal error: Error stack size too small."
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

    str = "Traceback (most recent call last):"
    do i = error_stack_position, 1, -1

       write (linestr, *)  error_stack(i)%line

       if (error_stack(i)%has_doc) then
       
          str = str // C_NEW_LINE // &
               '  File "' // &
               trim(error_stack(i)%fn) // &
               '", line ' // &
               trim(linestr) // &
               C_NEW_LINE // &
               "    " // &
               trim(error_stack(i)%doc)

       else

          str = str // C_NEW_LINE // &
               '  File "' // &
               trim(error_stack(i)%fn) // &
               '", line ' // &
               trim(linestr)

       endif

    enddo

    error_stack_position = 0

    if (present(error)) then
       error  = ERROR_NONE
    endif

  endfunction get_error_string_and_clear


#ifndef HAVE_QUIPPY
  !% Quit with an error message. Calls 'MPI_Abort' for MPI programs.
  subroutine error_abort_with_message(message)
    character(*),      intent(in) :: message
#ifdef IFORT_TRACEBACK_ON_ABORT
    integer :: j
#endif IFORT_TRACEBACK_ON_ABORT
#ifdef SIGNAL_ON_ABORT
    integer :: status
    integer, parameter :: SIGUSR1 = 30
#endif
#ifdef _MPI
    integer::PRINT_ALWAYS
    include "mpif.h"
#endif

#ifdef _MPI
    write(unit=error_unit, fmt='(a,i0," ",a)') 'SYSTEM ABORT: proc=',mpi_id(),trim(message)
#else
    write(unit=error_unit, fmt='(a," ",a)') 'SYSTEM ABORT:', trim(message)
#endif

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
#endif

  !% Stop program execution since this error is not handled properly
  subroutine error_abort_from_stack(error)
    implicit none

    integer, intent(inout), optional :: error

    ! ---

    ! This is for compatibility with quippy, change to error_abort
    call system_abort(get_error_string_and_clear(error))

  endsubroutine error_abort_from_stack

endmodule error_module
