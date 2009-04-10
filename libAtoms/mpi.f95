!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     libAtoms: atomistic simulation library
!X     
!X     Copyright 2006-2007.
!X
!X     Authors: Gabor Csanyi, Steven Winfield, James Kermode
!X     Contributors: Noam Bernstein, Alessio Comisso
!X
!X     The source code is released under the GNU General Public License,
!X     version 2, http://www.gnu.org/copyleft/gpl.html
!X
!X     If you would like to license the source code under different terms,
!X     please contact Gabor Csanyi, gabor@csanyi.net
!X
!X     When using this software, please cite the following reference:
!X
!X     http://www.libatoms.org
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X  mpi subroutines
!X  
!X  
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! $Id: mpi.f95,v 1.18 2008-07-14 10:23:36 jrk33 Exp $

! $Log: not supported by cvs2svn $
! Revision 1.17  2008/04/23 13:18:19  nb326
! Move mpif.h from header into individual routines, to work around weird Intel compiler but with fox_sax symbols being exported from modules that should make them private
!
! Revision 1.16  2007/10/11 18:01:06  nb326
! Fix #ifdef typo
!
! Revision 1.15  2007/10/11 16:29:38  nb326
! Now only low level routines (don't use system_module), higher level stuff got moved to System.f95
!
! Revision 1.14  2007/10/08 15:20:35  saw44
! removed parallel_print and parallel_warning
!
! Revision 1.13  2007/04/20 11:29:04  jrk33
! Changed size to nproc in get_mpi_size_rank since size is an intrinsic and ifort complains if used as a dummy argument name
!
! Revision 1.12  2007/04/18 01:32:21  gc121
! updated to reflect changes in printing and other naming conventions
!
! Revision 1.11  2007/04/17 17:11:18  jrk33
! Standardised subroutine and function references and printing argument order.
!
! Revision 1.10  2007/04/17 09:57:19  gc121
! put copyright statement in each file
!
! Revision 1.9  2007/03/28 17:48:53  saw44
! Changed no-mpi Decode_MPI_Error to report routine where spurious mpi calls are coming from
!
! Revision 1.8  2007/03/13 13:49:47  jrk33
! All output now goes to logger, apart from System_Abort which does go to stderr.
!
! Revision 1.7  2007/03/12 17:08:27  jrk33
! Reformatted documentaion; IN/OUT/INOUT -> lowercase
!
! Revision 1.6  2007/03/01 14:37:12  jrk33
! Fixed type
!
! Revision 1.5  2007/03/01 13:51:46  jrk33
! Documentation comments reformatted and edited throughout. Anything starting "!(no space)%"
!  is picked up by the documentation generation script
!
! Revision 1.4  2007/02/28 15:49:03  saw44
! Added Parallel_Warning routine
!
! Revision 1.3  2007/01/17 18:31:02  saw44
! Corrected variables in decode_error. Added Parallel_Print routine.
!
! Revision 1.2  2007/01/03 10:38:11  nb326
! MPI is a C++ namespace, replace ifdef MPI with ifdef _MPI
!
! Revision 1.1.1.1  2006/12/04 11:11:30  gc121
! Imported sources
!
! Revision 1.4  2006/06/20 17:23:19  gc121
! added new copyright notice to include James, Gian, Mike and Alessandro
!
! Revision 1.3  2006/02/27 17:00:46  gc121
! now compiles without the -DMPI flag. just gives an error if called
!
! Revision 1.2  2006/02/27 11:11:34  gc121
! added mpi query subroutines and cvs magic
!

module mpi_module
implicit none
private

public :: get_mpi_size_rank, decode_mpi_error

contains

!% Return the mpi size and rank for the communicator 'comm'.
!% this routine aborts of _MPI is not defined
subroutine get_mpi_size_rank(comm, nproc, rank)

  integer, intent(in)  :: comm  !% MPI communicator
  integer, intent(out) :: nproc  !% Total number of processes
  integer, intent(out) :: rank  !% Rank of this process

#ifdef _MPI
include 'mpif.h'
#endif

#ifdef _MPI

  integer::error_code

  call MPI_COMM_SIZE(comm, nproc, error_code)
  if (error_code .ne. MPI_SUCCESS) then
    rank=-1
    nproc=-1
    return
  endif
  call MPI_COMM_RANK(comm, rank, error_code)
  if (error_code .ne. MPI_SUCCESS) then
    rank=-1
    nproc=-1
    return
  endif
#else
  rank = 0
  nproc = 1
#endif
end subroutine get_mpi_size_rank


!% Decode MPI error, abort on failure.
function decode_mpi_error(error_code, routine_name, error_string)
  integer, intent(in) :: error_code
  character(len=*), intent(in) :: routine_name !% Name of MPI routine in which error occured
  character(len=*), intent(out) :: error_string
  logical :: decode_mpi_error

#ifdef _MPI
include 'mpif.h'
#endif

#ifdef _MPI

  character(MPI_MAX_ERROR_STRING)::raw_error_string
  integer::error_string_length, my_error_code

  if(error_code .ne. MPI_SUCCESS) then
     call  MPI_ERROR_STRING(error_code, raw_error_string, error_string_length, my_error_code)
     if(my_error_code .ne. MPI_SUCCESS) then
	write (error_string, '(A,A,I0,A)'), trim(routine_name), " returned with error code = ", error_code, &
	  ", which could not be parsed"
	decode_mpi_error = .true.
     else
	write (error_string, '(A, A, A, A)'), trim(routine_name), " had error '" , trim(raw_error_string), "'"
	decode_mpi_error = .true.
     endif
  else
#endif
    error_string=""
    decode_mpi_error = .false.
#ifdef _MPI
  endif
#endif

end function decode_mpi_error

end module mpi_module
