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
!X Error handling, see error.f95 for the functions called in these macros.
!X
!X Error passing works as follows:
!X - *error* needs to be intent(out) and optional
!X - all functions that receive *error* as an argument must call INIT_ERROR(error)
!X - RAISE_ERROR is used whenever an error occurs. If *error* is not present,
!X   the program execution will be terminated immediately. If *error* is
!X   present it will be set to some value not equal ERROR_NONE and the execution
!X   of the subroutine will be stopped.
!X - PASS_ERROR is used after a function or subroutine that returns error, i.e.
!X     call sub(..., error=error)
!X     PASS_ERROR(error)
!X   If no error occurs (i.e. error==ERROR_NONE), execution will proceed as
!X   usual. If an error occured, the current function will be terminated after
!X   the location of the error is passed to the error module.
!X   If the calling routine handles the error itself, rather than passing
!X   it up with PASS_ERROR(), CLEAR_ERROR() should be used to clear the error
!X   info stack
!X - PASS_ERROR_WITH_INFO is like PASS_ERROR, just an additional string can be
!X   provided describing the error, or parameters.
!X - HANDLE_ERROR will print the error history and stop execution of the program
!X   after an error occured.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#define INIT_ERROR(error)  if (present(error)) then ; error = ERROR_NONE ; endif
#define ASSERT(condition, message, error)   if (.not. (condition)) then ; RAISE_ERROR(message, error) ; endif

#define RAISE_ERROR(message, error)   if (.true.) then ; call push_error_with_info(message, __FILE__, __LINE__) ; if (present(error)) then ; error = ERROR_UNSPECIFIED ; return ; else ; call error_abort(error) ; endif ; endif

#define RAISE_ERROR_WITH_KIND(kind, message, error)   if (.true.) then ; call push_error_with_info(message, __FILE__, __LINE__, kind) ; if (present(error)) then ; error = kind ; return ; else ; call error_abort(error) ; endif ; endif

#define PASS_ERROR(error)  if (present(error)) then ; if (error /= ERROR_NONE) then ; call push_error(__FILE__, __LINE__) ; return ; endif ; endif

#define PASS_ERROR_WITH_INFO(message, error)  if (present(error)) then ; if (error /= ERROR_NONE) then ; call push_error_with_info(message, __FILE__, __LINE__) ; return ; endif ; endif

#define HANDLE_ERROR(error)  if (error /= ERROR_NONE) then ; call push_error(__FILE__, __LINE__) ; call error_abort(error) ; endif

#define CLEAR_ERROR(error) call error_clear_stack()


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X MPI errors
!X
!% MPI error string are obtained using mpi_error_string and then pushed
!% onto the error stack.
!%
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#define PASS_MPI_ERROR(mperror, error)  if (mperror /= MPI_SUCCESS) then ; call push_MPI_error(mperror, __FILE__, __LINE__) ; if (present(error)) then ; error = ERROR_MPI ; return ; else ; call error_abort(error) ; endif ; endif

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X MPI BCAST errors
!X
!X Extension of error handling macros to MPI cases where processes
!X peform different tasks. If an error occurs on one process it will
!X be broadcast to all others before the error is propagated
!X upwards. Replace RAISE_ERROR with BCAST_RAISE_ERROR and PASS_ERROR
!X with BCAST_PASS_ERROR. Additionally, BCAST_CHECK_ERROR must be
!X called on the processes in which no error has occured.  See
!X CInOutput read() for an example usage of these macros.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#define BCAST_ASSERT(condition, message, error, mpi) if (.not. (condition)) then ; BCAST_RAISE_ERROR(message, error, mpi) ; endif

#define BCAST_RAISE_ERROR(message, error, mpi) if (present(error)) call bcast(mpi, error); RAISE_ERROR(message, error)

#define BCAST_RAISE_ERROR_WITH_KIND(kind, message, error, mpi) if (present(error)) then; error = kind; call bcast(mpi, error); end if; RAISE_ERROR_WITH_KIND(kind, message, error)

#define BCAST_PASS_ERROR(error, mpi) if (present(error)) then; if (error /= ERROR_NONE) call bcast(mpi, error); endif; PASS_ERROR(error)

#define BCAST_CHECK_ERROR(error, mpi) if (present(error)) then; call bcast(mpi, error); if (error /= ERROR_NONE) then; RAISE_ERROR_WITH_KIND(error, "An error occured on another MPI process", error); endif; endif


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X Delayed errors - for OpenMP loops
!X
!X A subroutine currently in an OpenMP section cannot be quit using
!X the *return* statement. Hence, the error flag is set using
!X RAISE_DELAYED_ERROR and TRACE_DELAYED_ERROR. After the OpenMP section
!X has finished, INVOKE_DELAYED_ERROR will raise the error and exit
!X the current subroutine if an error occured in the OpenMP section.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#define RAISE_DELAYED_ERROR(message, error_loc)   if (error_loc == ERROR_NONE) then ; call push_error_with_info(message, __FILE__, __LINE__) ; error_loc = ERROR_UNSPECIFIED ; endif

#define TRACE_DELAYED_ERROR(error_loc)  if (error_loc /= ERROR_NONE) then ; call push_error(__FILE__, __LINE__) ; endif

#define TRACE_DELAYED_ERROR_WITH_INFO(message, error_loc)  if (error_loc /= ERROR_NONE) then ; call push_error_with_info(message, __FILE__, __LINE__) ; endif

#define INVOKE_DELAYED_ERROR(error_loc, error)  if (error_loc /= ERROR_NONE) then ; call push_error(__FILE__, __LINE__) ; if (present(error)) then ; error = error_loc ; else ; call error_abort(error) ; endif ; endif

