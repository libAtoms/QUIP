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
  use extendable_str_module
  use system_module

  implicit none
  private

  ! ---

  integer, parameter  :: ERROR_STACK_SIZE  = 100

  ! ---

  integer, parameter  :: ERROR_NONE        = 0
  integer, parameter  :: ERROR_OCCURED     = -1
  integer, parameter  :: ERROR_HAS_INFO    = -2

  public :: ERROR_NONE, ERROR_OCCURED, ERROR_HAS_INFO

  ! ---

  type ErrorDescriptor
     integer               :: kind  !% kind of error
     type(extendable_str)  :: doc   !% documentation string
     type(extendable_str)  :: fn    !% file name where the error occured
     integer               :: line  !% code line where the error occured
  endtype ErrorDescriptor

  ! ---

  save
  integer                :: error_stack_position  = 0       !% If this is zero, no error has occured
  type(ErrorDescriptor)  :: error_stack(ERROR_STACK_SIZE)   !% Error stack

  ! ---

  public :: push_error, push_error_with_info, get_error_string_and_clear
  public :: abort_on_error

contains

  !% Push a new error callback to the stack
  subroutine push_error(fn, line)
    implicit none

    character(*), intent(in)  :: fn
    integer, intent(in)       :: line

    ! ---

    !$omp critical

    error_stack_position  = error_stack_position + 1

    if (error_stack_position > ERROR_STACK_SIZE) then
       stop "Fatal error: Error stack size too small."
    endif

    error_stack(error_stack_position)%kind = ERROR_OCCURED
    call initialise(error_stack(error_stack_position)%fn)
    call concat(error_stack(error_stack_position)%fn, fn)
    error_stack(error_stack_position)%line = line

    !$omp end critical

  endsubroutine push_error


  !% Push a new information string onto the error stack
  subroutine push_error_with_info(doc, fn, line)
    implicit none

    character(*), intent(in)  :: doc
    character(*), intent(in)  :: fn
    integer, intent(in)       :: line

    ! ---

    !$omp critical

    error_stack_position  = error_stack_position + 1

    if (error_stack_position > ERROR_STACK_SIZE) then
       call system_abort("Fatal error: Error stack size too small.")
    endif

    error_stack(error_stack_position)%kind = ERROR_HAS_INFO
    call initialise(error_stack(error_stack_position)%fn)
    call concat(error_stack(error_stack_position)%fn, fn)
    error_stack(error_stack_position)%line = line
    call initialise(error_stack(error_stack_position)%doc)
    call concat(error_stack(error_stack_position)%doc, doc)

    !$omp end critical

  endsubroutine push_error_with_info


  !% Construct a string describing the error.
  function get_error_string_and_clear(ierror) result(str)
    implicit none

    integer, intent(inout), optional  :: ierror

    type(extendable_str)              :: str

    ! ---

    integer               :: i

    ! ---

    call initialise(str)
    call concat(str, "Internal error. Traceback (most recent call last):")
    do i = error_stack_position, 1, -1

       if (error_stack(i)%kind == ERROR_HAS_INFO) then
       
          call concat(str, "\n" // &
               '  File "' // &
               string(error_stack(i)%fn) // &
               '", line ' // &
               error_stack(i)%line // &
               "\n" // &
               "    " // &
               string(error_stack(i)%doc))

       else

          call concat(str, "\n" // &
               '  File "' // &
               string(error_stack(i)%fn) // &
               '", line ' // &
               error_stack(i)%line)

       endif

    enddo

    error_stack_position = 0

    if (present(ierror)) then
       ierror  = ERROR_NONE
    endif

  endfunction get_error_string_and_clear


  !% Stop program execution since this error is not handled properly
  subroutine abort_on_error(ierror)
    implicit none

    integer, intent(inout), optional :: ierror

    ! ---

    call system_abort(string(get_error_string_and_clear(ierror)))

  endsubroutine abort_on_error

endmodule error_module
