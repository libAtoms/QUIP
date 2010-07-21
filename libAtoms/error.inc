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
!X Error handling, see error.f90 for the functions called in these macros.
!X
!% Error passing works as follows:
!% - *error* needs to be intent(inout) and optional
!% - RAISE_ERROR is used whenever an error occurs. If *error* is not present,
!%   the program execution will be terminated immediately. If *error* is
!%   present it will be set to some value not equal ERROR_NONE and the execution
!%   of the subroutine will be stopped.
!% - PASS_ERROR is used after a function that returns error, i.e.
!%     call func(..., error=error)
!%     PASS_ERROR
!%   If no error occurs (i.e. error==ERROR_NONE), execution will proceed as
!%   usual. If an error occured, the current function will be terminated after
!%   the location of the error is passed to the error module.
!% - PASS_ERROR_WITH_INFO is like PASS_ERROR, just an additional string can be
!%   provided describing the error, or parameters.
!% - HANDLE_ERROR will print the error history and stop execution of the program
!%   after an error occured.
!%
!% On the top level, the error variable has to be initialized with ERROR_NONE.
!%
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


#define RAISE_ERROR(info, error)   call push_error_with_info(info, __FILE__, __LINE__) ; if (present(error)) then ; error = ERROR_OCCURED ; return ; else ; call error_abort ; endif

#define PASS_ERROR(error)  if (present(error)) then ; if (error /= ERROR_NONE) then ; call push_error(__FILE__, __LINE__) ; return ; endif ; endif

#define PASS_ERROR_WITH_INFO(info, error)  if (present(error)) then ; if (error /= ERROR_NONE) then ; call push_error_with_info(info, __FILE__, __LINE__) ; return ; endif ; endif

#define HANDLE_ERROR(error)  if (error /= ERROR_NONE) then ; call push_error(__FILE__, __LINE__) ; call error_abort ; endif



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X MPI errors
!X
!% MPI error string are obtained using mpi_error_string and then pushed
!% onto the error stack.
!%
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#define PASS_MPI_ERROR(mperror, error)  if (mperror /= MPI_SUCCESS) then ; call push_MPI_error(mperror, __FILE__, __LINE__) ; if (present(error)) then ; error = ERROR_OCCURED ; return ; else ; call error_abort ; endif ; endif



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X Delayed errors - for OpenMP loops
!X
!% A subroutine currently in an OpenMP section cannot be quit using
!% the *return* statement. Hence, the error flag is set using
!% RAISE_DELAYED_ERROR and TRACE_DELAYED_ERROR. After the OpenMP section
!% has finished, INVOKE_DELAYED_ERROR will raise the error and exit
!% the current subroutine if an error occured in the OpenMP section.
!%
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#define RAISE_DELAYED_ERROR(info, error_loc)   if (error_loc == ERROR_NONE) then ; call push_error_with_info(info, __FILE__, __LINE__) ; error_loc = ERROR_OCCURED ; endif

#define TRACE_DELAYED_ERROR(error_loc)  if (error_loc /= ERROR_NONE) then ; call push_error(__FILE__, __LINE__) ; endif

#define TRACE_DELAYED_ERROR_WITH_INFO(info, error_loc)  if (error_loc /= ERROR_NONE) then ; call push_error_with_info(info, __FILE__, __LINE__) ; endif

#define INVOKE_DELAYED_ERROR(error_loc, error)  if (error_loc /= ERROR_NONE) then ; call push_error(__FILE__, __LINE__) ; if (present(error)) then ; error = error_loc ; else ; call error_abort ; endif ; endif

