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
!X ScaLAPACK module 
!X
!% Module wrapping ScaLAPACK routines
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module ScaLAPACK_module

use System_module
use MPI_context_module

implicit none
private

#ifdef SCALAPACK
integer, external :: indxl2g, numroc
#endif

integer dlen_
parameter (dlen_ = 50)

public :: ScaLAPACK
type ScaLAPACK
  logical :: active = .false.
  type(MPI_context) MPI_obj
  integer :: blacs_context
  integer :: n_proc_rows = 1, n_proc_cols = 1, my_proc_row = 0, my_proc_col = 0
end type ScaLAPACK

public :: Matrix_ScaLAPACK_Info
type Matrix_ScaLAPACK_Info
  logical :: active = .false.
  type(ScaLAPACK) :: ScaLAPACK_obj
  integer :: N_R = 0, N_C = 0, NB_R = 0, NB_C = 0
  integer desc(dlen_)
  integer :: l_N_R = 0, l_N_C = 0
end type Matrix_ScaLAPACK_Info

public :: Initialise
interface Initialise
  module procedure ScaLAPACK_Initialise
  module procedure Matrix_ScaLAPACK_Info_Initialise
end interface Initialise

public :: Finalise
interface Finalise
  module procedure ScaLAPACK_Finalise
  module procedure Matrix_ScaLAPACK_Info_Finalise
end interface Finalise

public :: Wipe
interface Wipe
  module procedure Matrix_ScaLAPACK_Info_Wipe
end interface Wipe

public :: Print
interface Print
  module procedure ScaLAPACK_Print, Matrix_ScaLAPACK_Info_Print
  module procedure ScaLAPACK_Matrix_d_print
  module procedure ScaLAPACK_Matrix_z_print
end interface Print

public :: init_matrix_desc
interface init_matrix_desc
  module procedure ScaLAPACK_init_matrix_desc
end interface init_matrix_desc

public :: coords_local_to_global
interface coords_local_to_global
  module procedure Matrix_ScaLAPACK_Info_coords_local_to_global
end interface coords_local_to_global

public :: coords_global_to_local
interface coords_global_to_local
  module procedure Matrix_ScaLAPACK_Info_coords_global_to_local
end interface coords_global_to_local

public :: diagonalise
interface diagonalise
  module procedure ScaLAPACK_diagonalise_r, ScaLAPACK_diagonalise_c
  module procedure ScaLAPACK_diagonalise_gen_r, ScaLAPACK_diagonalise_gen_c
end interface

public :: inverse
interface inverse
  module procedure ScaLAPACK_inverse_r, ScaLAPACK_inverse_c
end interface

public :: add_identity
interface add_identity
  module procedure ScaLAPACK_add_identity_r
end interface

public :: matrix_product_sub
interface matrix_product_sub
  module procedure ScaLAPACK_matrix_product_sub_ddd, ScaLAPACK_matrix_product_sub_zzz
end interface

public :: matrix_product_vect_asdiagonal_sub
interface matrix_product_vect_asdiagonal_sub
  module procedure ScaLAPACK_matrix_product_vect_asdiagonal_sub_ddd
  module procedure ScaLAPACK_matrix_product_vect_asdiagonal_sub_zzd
  module procedure ScaLAPACK_matrix_product_vect_asdiagonal_sub_zzz
end interface matrix_product_vect_asdiagonal_sub

public :: Re_diag
interface Re_diag
  module procedure ScaLAPACK_Re_diagZ, ScaLAPACK_Re_diagD
end interface Re_diag

public :: diag_spinor
interface diag_spinor
  module procedure ScaLAPACK_diag_spinorZ, ScaLAPACK_diag_spinorD
end interface diag_spinor

contains

subroutine Matrix_ScaLAPACK_Info_Initialise(this, N_R, N_C, NB_R, NB_C, scalapack_obj)
  type(Matrix_ScaLAPACK_Info), intent(inout) :: this
  integer, intent(in) :: N_R, NB_R
  integer, intent(in), optional :: N_C, NB_C
  type(ScaLAPACK), intent(in), optional :: scalapack_obj

#ifdef SCALAPACK

  call Finalise(this)

  this%N_R = N_R
  this%NB_R = NB_R

  if (present(N_C)) then
    this%N_C = N_C
  else
    this%N_C = N_R
  endif
  if (present(NB_C)) then
    this%NB_C = NB_C
  else
    this%NB_C = NB_R
  endif

  if (present(scalapack_obj)) then
    if (scalapack_obj%active) then
      this%ScaLAPACK_obj = scalapack_obj
      call init_matrix_desc(this%ScaLAPACK_obj, this%N_R, this%N_C, this%NB_R, this%NB_C, &
	this%desc, this%l_N_R, this%l_N_C)
    endif
  else
    this%l_N_R = this%N_R
    this%l_N_C = this%N_C
  endif

  this%active = this%ScaLAPACK_obj%active
#endif
end subroutine Matrix_ScaLAPACK_Info_Initialise

subroutine ScaLAPACK_Initialise(this, MPI_obj)
  type(ScaLAPACK), intent(inout) :: this
  type(MPI_context), intent(in), optional :: MPI_obj

#ifdef SCALAPACK
  call Finalise(this)

  this%active = .false.
  if (present(MPI_obj)) then
    if (MPI_obj%active .and. MPI_obj%n_procs > 1) then
      this%active = .true.

      this%MPI_obj = MPI_obj
      this%blacs_context = MPI_obj%communicator

      call calc_n_proc_rows_cols(this%MPI_obj%n_procs, this%n_proc_rows, this%n_proc_cols)
      call print("ScaLAPACK_Initialise using proc grid " // this%n_proc_rows // " x " // this%n_proc_cols, PRINT_VERBOSE)

      call blacs_gridinit (this%blacs_context, 'R', this%n_proc_rows, this%n_proc_cols)
      call blacs_gridinfo (this%blacs_context, this%n_proc_rows, this%n_proc_cols, this%my_proc_row, this%my_proc_col)
    endif
  endif

#endif
end subroutine ScaLAPACK_Initialise

subroutine Matrix_ScaLAPACK_Info_Finalise(this)
  type(Matrix_ScaLAPACK_Info), intent(inout) :: this

#ifdef SCALAPACK
  this%active = .false.
#endif
end subroutine

subroutine ScaLAPACK_Finalise(this)
  type(ScaLAPACK), intent(inout) :: this

#ifdef SCALAPACK
  if (this%active) then
    call blacs_gridexit(this%blacs_context)
  endif
  this%active = .false.
#endif
end subroutine

subroutine Matrix_ScaLAPACK_Info_Wipe(this)
  type(Matrix_ScaLAPACK_Info), intent(inout) :: this

#ifdef SCALAPACK
  this%active = .false.
#endif
end subroutine Matrix_ScaLAPACK_Info_Wipe


subroutine Matrix_ScaLAPACK_Info_coords_global_to_local(this, i, j, l_i, l_j, l_row_p, l_col_p)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this
  integer, intent(in) :: i, j
  integer, intent(out) :: l_i, l_j
  integer, intent(out), optional :: l_row_p, l_col_p

#ifdef SCALAPACK
  integer :: l_row_p_use, l_col_p_use

  if (this%active) then
    call infog2l(i, j, this%desc, this%ScaLAPACK_obj%n_proc_rows, this%ScaLAPACK_obj%n_proc_cols, &
      this%ScaLAPACK_obj%my_proc_row, this%ScaLAPACK_obj%my_proc_col, &
      l_i, l_j, l_row_p_use, l_col_p_use)

    if (l_row_p_use /= this%ScaLAPACK_obj%my_proc_row .or. l_col_p_use /= this%ScaLAPACK_obj%my_proc_col) then
      l_i = -1
      l_j = -1
    endif

    if (present(l_row_p)) l_row_p = l_row_p_use
    if (present(l_col_p)) l_col_p = l_col_p_use
  else
#endif
    l_i = i
    l_j = j
    if (present(l_row_p)) l_row_p = 0
    if (present(l_col_p)) l_col_p = 0
#ifdef SCALAPACK
  endif
#endif
end subroutine Matrix_ScaLAPACK_Info_coords_global_to_local

subroutine Matrix_ScaLAPACK_Info_coords_local_to_global(this, l_i, l_j, i, j)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this
  integer, intent(in) :: l_i, l_j
  integer, intent(out) :: i, j

#ifdef SCALAPACK

  if (this%active) then
    i = indxl2g(l_i, this%NB_R, this%ScaLAPACK_obj%my_proc_row, 0, this%ScaLAPACK_obj%n_proc_rows)
    j = indxl2g(l_j, this%NB_C, this%ScaLAPACK_obj%my_proc_col, 0, this%ScaLAPACK_obj%n_proc_cols)
  else
#endif
    i = l_i
    j = l_j
#ifdef SCALAPACK
  endif
#endif

end subroutine Matrix_ScaLAPACK_Info_coords_local_to_global


subroutine ScaLAPACK_init_matrix_desc(this, N_R, N_C, NB_R, NB_C, desc, l_N_R, l_N_C)
  type(ScaLAPACK), intent(in) :: this
  integer, intent(in) :: N_R, NB_R
  integer, intent(in), optional :: N_C, NB_C
  integer, intent(out) :: desc(dlen_)
  integer, intent(out) :: l_N_R, l_N_C

#ifdef SCALAPACK
  integer err
  integer use_N_C, use_NB_C
  integer lld
#endif

  desc = 0
  l_N_R = 0
  l_N_C = 0

#ifdef SCALAPACK
  if (this%active) then
    if (present(N_C)) then
      use_N_C = N_C
    else
      use_N_C = N_R
    endif

    if (present(NB_C)) then
      use_NB_C = NB_C
    else
      use_NB_C = NB_R
    endif

    l_N_R = numroc(N_R, NB_R, this%my_proc_row, 0, this%n_proc_rows)
    l_N_C = numroc(use_N_C, use_NB_C, this%my_proc_col, 0, this%n_proc_cols)

    lld = l_N_R
    if (l_N_r < 1) lld = 1

    call descinit (desc, N_R, use_N_C, NB_R, use_NB_C, 0, 0, this%blacs_context, lld, err)
  endif
#endif
end subroutine ScaLAPACK_init_matrix_desc



subroutine calc_n_proc_rows_cols(n_procs, n_proc_rows, n_proc_cols)
  integer, intent(in) :: n_procs
  integer, intent(out) :: n_proc_rows, n_proc_cols

#ifdef SCALAPACK
  integer n_proc_rows_t, n_proc_cols_t


! long and narrow, or short and wide
!  do n_proc_cols_t=int(sqrt(dble(n_procs))), 1, -1
!    n_proc_rows_t = n_procs/n_proc_cols_t
  do n_proc_rows_t=int(sqrt(dble(n_procs))), 1, -1
    n_proc_cols_t = n_procs/n_proc_rows_t

    if (n_proc_rows_t*n_proc_cols_t .eq. n_procs) then
	n_proc_rows = n_proc_rows_t
	n_proc_cols = n_proc_cols_t
	exit
    endif
  end do
#else
  n_proc_rows = 0
  n_proc_cols = 0
#endif
end subroutine calc_n_proc_rows_cols

subroutine ScaLAPACK_Print(this,file)
  type(ScaLAPACK),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

#ifdef SCALAPACK

  call Print("ScaLAPACK : ", file=file)

  call Print ('ScaLAPACK : active ' // this%active, file=file)
  if (this%active) then
    call Print ('ScaLAPACK : n_proc rows cols ' // this%n_proc_rows // " " // this%n_proc_cols, file=file)
    call Print ('ScaLAPACK : my row col ' //  this%my_proc_row // " " // this%my_proc_col, file=file)
  endif
#endif
end subroutine ScaLAPACK_Print

subroutine Matrix_ScaLAPACK_Info_Print(this,file)
  type(Matrix_ScaLAPACK_Info),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file
#ifdef SCALAPACK

  call Print('Matrix_ScaLAPACK_Info : ', file=file)

  call Print ('Matrix_ScaLAPACK_Info : active ' // this%active, file=file)
  if (this%active) then
    call Print(this%ScaLAPACK_obj, file=file)
    call Print ('Matrix_ScaLAPACK_Info : N_R N_C ' // this%N_R // " " // this%N_C, file=file)
    call Print ('Matrix_ScaLAPACK_Info : NB_R NB_C ' // this%NB_R // " " // this%NB_C, file=file)
    call Print ('Matrix_ScaLAPACK_Info : l_N_R l_N_C ' // this%l_N_R // " " // this%l_N_C, file=file)
  endif
#endif
end subroutine Matrix_ScaLAPACK_Info_Print

subroutine ScaLAPACK_inverse_r(this, data, inv_scalapack_info, inv_data, positive_in)
  type(Matrix_ScaLAPACK_Info), intent(in), target :: this
  real(dp), intent(inout), target :: data(:,:)
  type(Matrix_ScaLAPACK_Info), intent(in), target, optional :: inv_scalapack_info
  real(dp), intent(out), target, optional :: inv_data(:,:)
  logical, intent(in), optional :: positive_in

#ifdef SCALAPACK
  real(dp), pointer :: u_inv_data(:,:)
  integer, pointer :: u_inv_desc(:)
  integer, allocatable :: ipiv(:)
  integer :: info
  integer :: lwork, liwork
  real(dp), allocatable :: work(:)
  integer, allocatable :: iwork(:)

  if (present(inv_data)) then
    u_inv_data => inv_data
    u_inv_data = data
  else
    u_inv_data => data
  endif
  if (present(inv_scalapack_info)) then
    if (.not. present(inv_data)) call system_abort("ScaLAPACK_inverse_r called with inv_scalapack_info but w/o inv_data")
    u_inv_desc => inv_scalapack_info%desc
  else
    u_inv_desc => this%desc
  endif

  allocate(ipiv(numroc(this%N_C, this%NB_C, this%scalapack_obj%my_proc_row, 0, this%scalapack_obj%n_proc_rows)+this%NB_C))
  call pdgetrf(this%N_R, this%N_C, u_inv_data, 1, 1, u_inv_desc, ipiv, info)
  if (info /= 0) then
    call system_abort("ScaLAPACK_inverse_r got pdgetrf info " // info)
  endif

  allocate(work(1))
  allocate(iwork(1))
  lwork = -1
  liwork = -1
  call pdgetri(this%N_R, u_inv_data, 1, 1, u_inv_desc, &
    ipiv, work, lwork, iwork, liwork, info)
  lwork = work(1)
  liwork = iwork(1)
  deallocate(work)
  deallocate(iwork)

  allocate(work(lwork))
  allocate(iwork(liwork))
  call pdgetri(this%N_R, u_inv_data, 1, 1, u_inv_desc, &
    ipiv, work, lwork, iwork, liwork, info)
  if (info /= 0) then
    call system_abort("ScaLAPACK_inverse_r got pdgetrf info " // info)
  endif

  deallocate(work)
  deallocate(iwork)
  deallocate(ipiv)
#endif

end subroutine ScaLAPACK_inverse_r

subroutine ScaLAPACK_inverse_c(this, data, inv_scalapack_info, inv_data, positive_in)
  type(Matrix_ScaLAPACK_Info), intent(in), target :: this
  complex(dp), intent(inout), target :: data(:,:)
  type(Matrix_ScaLAPACK_Info), intent(in), target, optional :: inv_scalapack_info
  complex(dp), intent(out), target, optional :: inv_data(:,:)
  logical, intent(in), optional :: positive_in

#ifdef SCALAPACK
  complex(dp), pointer :: u_inv_data(:,:)
  integer, pointer :: u_inv_desc(:)
  integer, allocatable :: ipiv(:)
  integer :: info

  integer :: lwork, liwork
  complex(dp), allocatable :: work(:)
  integer, allocatable :: iwork(:)

! call system_timer("ScaLAPACK_inverse_c")
! call system_timer("ScaLAPACK_inverse_c/prep")
  if (present(inv_data)) then
    u_inv_data => inv_data
    u_inv_data = data
  else
    u_inv_data => data
  endif
  if (present(inv_scalapack_info)) then
    if (.not. present(inv_data)) call system_abort("ScaLAPACK_inverse_r called with inv_scalapack_info but w/o inv_data")
    u_inv_desc => inv_scalapack_info%desc
  else
    u_inv_desc => this%desc
  endif

  allocate(ipiv(numroc(this%N_C, this%NB_C, this%scalapack_obj%my_proc_row, 0, this%scalapack_obj%n_proc_rows)+this%NB_C))
! call system_timer("ScaLAPACK_inverse_c/prep")
! call system_timer("ScaLAPACK_inverse_c/pzgetrf")
  call pzgetrf(this%N_R, this%N_C, u_inv_data, 1, 1, u_inv_desc, ipiv, info)
! call system_timer("ScaLAPACK_inverse_c/pzgetrf")
! call system_timer("ScaLAPACK_inverse_c/pzgetri_prep")
  if (info /= 0) then
    call system_abort("ScaLAPACK_inverse_r got pzgetrf info " // info)
  endif

  allocate(work(1))
  allocate(iwork(1))
  lwork = -1
  liwork = -1
  call pzgetri(this%N_R, u_inv_data, 1, 1, u_inv_desc, &
    ipiv, work, lwork, iwork, liwork, info)
  lwork = work(1)
  liwork = iwork(1)
  deallocate(work)
  deallocate(iwork)

  allocate(work(lwork))
  allocate(iwork(liwork))
! call system_timer("ScaLAPACK_inverse_c/pzgetri_prep")
! call system_timer("ScaLAPACK_inverse_c/pzgetri")
  call pzgetri(this%N_R, u_inv_data, 1, 1, u_inv_desc, &
    ipiv, work, lwork, iwork, liwork, info)
! call system_timer("ScaLAPACK_inverse_c/pzgetri")
! call system_timer("ScaLAPACK_inverse_c/post")
  if (info /= 0) then
    call system_abort("ScaLAPACK_inverse_r got pzgetrf info " // info)
  endif

  deallocate(work)
  deallocate(iwork)
  deallocate(ipiv)
! call system_timer("ScaLAPACK_inverse_c/post")
! call system_timer("ScaLAPACK_inverse_c")
#endif
end subroutine ScaLAPACK_inverse_c

subroutine ScaLAPACK_diagonalise_r(this, data, evals, evecs_ScaLAPACK_obj, evecs_data, error)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this
  real(dp), intent(in) :: data(:,:)
  real(dp), intent(out) :: evals(:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: evecs_ScaLAPACK_obj
  real(dp), intent(out) :: evecs_data(:,:)
  integer, intent(out), optional :: error

  integer liwork, lwork
  real(dp), allocatable :: data_copy(:,:)
  integer, allocatable :: iwork(:)
  real(dp), allocatable :: work(:)
  integer, allocatable :: icluster(:), ifail(:)
  real(dp), allocatable :: gap(:)

  integer info
  integer M, NZ, i
  real(dp) :: orfac = 1.0e-4_dp

  INIT_ERROR(error)

  evals = 0.0_dp
  evecs_data = 0.0_dp

#ifdef SCALAPACK
  allocate(icluster(2*this%ScaLAPACK_obj%n_proc_rows*this%ScaLAPACK_obj%n_proc_cols))
  allocate(gap(this%ScaLAPACK_obj%n_proc_rows*this%ScaLAPACK_obj%n_proc_cols))
  allocate(ifail(this%N_R))

  allocate(data_copy(this%l_N_R, this%l_N_C))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_r data_copy", size(data_copy)*REAL_SIZE)
  data_copy = data

  lwork = -1
  liwork = -1
  allocate(work(1))
  allocate(iwork(1))
  call pdsyevx('V', 'A', 'U', this%N_R, &
    data_copy, 1, 1, this%desc, &
    0.0_dp, 0.0_dp, 1, 1, -1.0_dp, M, NZ, evals, &
    orfac, evecs_data, 1, 1, evecs_ScaLAPACK_obj%desc, work, lwork, &
    iwork, liwork, ifail, icluster, gap, info)
  lwork=work(1)*4
  liwork=iwork(1)
  deallocate(work)
  deallocate(iwork)
  allocate(work(lwork))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_r work", size(work)*REAL_SIZE)
  allocate(iwork(liwork))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_r iwork", size(iwork)*INTEGER_SIZE)

  call pdsyevx('V', 'A', 'U', this%N_R, &
    data_copy, 1, 1, this%desc, &
    0.0_dp, 0.0_dp, 1, 1, -1.0_dp, M, NZ, evals, &
    orfac, evecs_data, 1, 1, evecs_ScaLAPACK_obj%desc, work, lwork, &
    iwork, liwork, ifail, icluster, gap, info)

  if (info .ne. 0) then
      call print("ScaLAPACK_diagonalise_r got info " // info, PRINT_ALWAYS)
      if (mod(info,2) .ne. 0) then
	  call print("   eigenvectors failed to converge", PRINT_ALWAYS)
      else if (mod(info/2,2) .ne. 0) then
          call print("   eigenvectors failed to orthogonalize", PRINT_ALWAYS)
	  do i=1, size(icluster), 2
	      if (icluster(i) /= 0) then
	          call print (" eval cluster " // icluster(i) // " " // icluster(i+1) // &
		    "("// (icluster(i+1)-icluster(i)+1) // ")", PRINT_ALWAYS)
	      else
	          exit
	      endif
	  end do
          do i=1, this%N_R
            call print(" eigenvalue " // i // " " // evals(i), PRINT_ANAL)
          end do
      else if (mod(info/4,2) .ne. 0) then
	  call print("   not enough space for all eigenvectors in range", PRINT_ALWAYS)
      else if (mod(info/8,2) .ne. 0) then
	  call print("   failed to compute eigenvalues", PRINT_ALWAYS)
      else if (mod(info/16,2) .ne. 0) then
	  call print("   S was not positive definite "//ifail(1), PRINT_ALWAYS)
      endif
  endif

  deallocate (icluster, ifail, gap)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_r work", size(work)*REAL_SIZE)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_r iwork", size(iwork)*INTEGER_SIZE)
  deallocate (iwork, work)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_r data_copy", size(data_copy)*REAL_SIZE)
  deallocate (data_copy)

  if (M /= this%N_R .or. NZ /= this%N_R) then ! bad enough error to report
    RAISE_ERROR("ScaLAPACK_diagonalise_r: Failed to diagonalise info="//info, error)
  endif
#endif

end subroutine ScaLAPACK_diagonalise_r

subroutine ScaLAPACK_diagonalise_c(this, data, evals, evecs_ScaLAPACK_obj, evecs_data, error)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this
  complex(dp), intent(in) :: data(:,:)
  real(dp), intent(out) :: evals(:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: evecs_ScaLAPACK_obj
  complex(dp), intent(out) :: evecs_data(:,:)
  integer, intent(out), optional :: error

  integer liwork, lrwork, lcwork
  integer, allocatable :: iwork(:)
  real(dp), allocatable :: rwork(:)
  complex(dp), allocatable :: cwork(:)
  integer, allocatable :: icluster(:), ifail(:)
  real(dp), allocatable :: gap(:)

  complex(dp), allocatable :: data_copy(:,:)

  integer info
  integer M, NZ, i
  real(dp) :: orfac = 1.0e-4_dp

  INIT_ERROR(error)

  evals = 0.0_dp
  evecs_data = 0.0_dp

#ifdef SCALAPACK
  allocate(icluster(2*this%ScaLAPACK_obj%n_proc_rows*this%ScaLAPACK_obj%n_proc_cols))
  allocate(gap(this%ScaLAPACK_obj%n_proc_rows*this%ScaLAPACK_obj%n_proc_cols))
  allocate(ifail(this%N_R))

  allocate(data_copy(this%l_N_R, this%l_N_C))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_c data_copy", size(data_copy)*COMPLEX_SIZE)
  data_copy = data

  lrwork = -1
  lcwork = -1
  liwork = -1
  allocate(rwork(1))
  allocate(cwork(1))
  allocate(iwork(1))
  call pzheevx('V', 'A', 'U', this%N_R, &
    data_copy, 1, 1, this%desc, &
    0.0_dp, 0.0_dp, 1, 1, -1.0_dp, M, NZ, evals, &
    orfac, evecs_data, 1, 1, evecs_ScaLAPACK_obj%desc, cwork, lcwork, &
    rwork, lrwork, iwork, liwork, ifail, icluster, gap, info)
  lrwork=rwork(1)*4
  lcwork=cwork(1)
  liwork=iwork(1)
  deallocate(rwork)
  deallocate(cwork)
  deallocate(iwork)
  allocate(rwork(lrwork))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_c rwork", size(rwork)*REAL_SIZE)
  allocate(cwork(lcwork))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_c cwork", size(cwork)*COMPLEX_SIZE)
  allocate(iwork(liwork))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_c iwork", size(iwork)*INTEGER_SIZE)

  call pzheevx('V', 'A', 'U', this%N_R, &
    data_copy, 1, 1, this%desc, &
    0.0_dp, 0.0_dp, 1, 1, -1.0_dp, M, NZ, evals, &
    orfac, evecs_data, 1, 1, evecs_ScaLAPACK_obj%desc, cwork, lcwork, &
    rwork, lrwork, iwork, liwork, ifail, icluster, gap, info)

  if (info .ne. 0) then
      call print("ScaLAPACK_diagonlise_c got info " // info, PRINT_ALWAYS)
      if (mod(info,2) .ne. 0) then
	  call print("   eigenvectors failed to converge", PRINT_ALWAYS)
      else if (mod(info/2,2) .ne. 0) then
          call print("   eigenvectors failed to orthogonalize", PRINT_ALWAYS)
	  do i=1, size(icluster), 2
	      if (icluster(i) /= 0) then
	          call print (" eval cluster " // icluster(i) // " " // icluster(i+1) // &
		    "("// (icluster(i+1)-icluster(i)+1) // ")", PRINT_ALWAYS)
	      else
	          exit
	      endif
	  end do
          do i=1, this%N_R
            call print(" eigenvalue " // i // " " // evals(i), PRINT_ANAL)
          end do
      else if (mod(info/4,2) .ne. 0) then
	  call print("   not enough space for all eigenvectors in range", PRINT_ALWAYS)
      else if (mod(info/8,2) .ne. 0) then
	  call print("   failed to compute eigenvalues", PRINT_ALWAYS)
      else if (mod(info/16,2) .ne. 0) then
	  call print("   S was not positive definite "//ifail(1), PRINT_ALWAYS)
      endif
  endif

  deallocate (icluster, ifail, gap)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_c iwork", size(iwork)*INTEGER_SIZE)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_c rwork", size(rwork)*REAL_SIZE)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_c cwork", size(cwork)*COMPLEX_SIZE)
  deallocate (iwork, rwork, cwork)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_c data_copy", size(data_copy)*COMPLEX_SIZE)
  deallocate (data_copy)

  if (M /= this%N_R .or. NZ /= this%N_R) then ! bad enough error to report
    RAISE_ERROR("ScaLAPACK_diagonalise_r: Failed to diagonalise info="//info, error)
  endif
#endif

end subroutine ScaLAPACK_diagonalise_c

subroutine ScaLAPACK_diagonalise_gen_r(this, data, overlap_ScaLAPACK_obj, overlap_data, &
  evals, evecs_ScaLAPACK_obj, evecs_data, error)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this
  real(dp), intent(in) :: data(:,:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: overlap_ScaLAPACK_obj
  real(dp), intent(in) :: overlap_data(:,:)
  real(dp), intent(out) :: evals(:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: evecs_ScaLAPACK_obj
  real(dp), intent(out) :: evecs_data(:,:)
  integer, intent(out), optional :: error

  integer :: liwork, lwork
  integer, allocatable :: iwork(:)
  real(dp), allocatable :: work(:)
  integer, allocatable :: icluster(:), ifail(:)
  real(dp), allocatable :: gap(:)

  real(dp), allocatable :: data_copy(:,:), overlap_data_copy(:,:)

  integer info
  integer M, NZ, i
  real(dp) :: orfac = 1.0e-4_dp

  INIT_ERROR(error)

  evals = 0.0_dp
  evecs_data = 0.0_dp

#ifdef SCALAPACK
  allocate(icluster(2*this%ScaLAPACK_obj%n_proc_rows*this%ScaLAPACK_obj%n_proc_cols))
  allocate(gap(this%ScaLAPACK_obj%n_proc_rows*this%ScaLAPACK_obj%n_proc_cols))
  allocate(ifail(this%N_R))

  allocate(data_copy(this%l_N_R, this%l_N_C))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_gen_r data_copy", size(data_copy)*REAL_SIZE)
  allocate(overlap_data_copy(overlap_ScaLAPACK_obj%l_N_R, overlap_ScaLAPACK_obj%l_N_C))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_gen_r overlap_data_copy", size(overlap_data_copy)*REAL_SIZE)
  data_copy = data
  overlap_data_copy = overlap_data

  lwork = -1
  liwork = -1
  allocate(work(1))
  allocate(iwork(1))
  call pdsygvx(1, 'V', 'A', 'U', this%N_R, &
    data_copy, 1, 1, this%desc, &
    overlap_data_copy, 1, 1, overlap_ScaLAPACK_obj%desc, &
    0.0_dp, 0.0_dp, 1, 1, -1.0_dp, M, NZ, evals, &
    orfac, evecs_data, 1, 1, evecs_ScaLAPACK_obj%desc, work, lwork, &
    iwork, liwork, ifail, icluster, gap, info)
  lwork = work(1)*4
  liwork = iwork(1)
  deallocate(work)
  deallocate(iwork)

  allocate(work(lwork))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_gen_r work", size(work)*REAL_SIZE)
  allocate(iwork(liwork))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_gen_r iwork", size(iwork)*INTEGER_SIZE)

! call print("calling pdsygv : NR " // this%N_R // " shape(data_copy) " // shape(data_copy) // &
! " shape(overlap_data_copy) " // shape(overlap_data_copy) // " shape(evals) " // shape(evals) // &
! " shape(evecs_data) "  // shape(evecs_data) // " shape(work) " // shape(work) // " lwork " // lwork // &
! " shape(iwork) " // shape(iwork) // " liwork " // liwork //  " shape(ifail) " // shape(ifail) // &
! " shape(icluster) " // shape(icluster), PRINT_ALWAYS)
  call pdsygvx(1, 'V', 'A', 'U', this%N_R, &
    data_copy, 1, 1, this%desc, &
    overlap_data_copy, 1, 1, overlap_ScaLAPACK_obj%desc, &
    0.0_dp, 0.0_dp, 1, 1, -1.0_dp, M, NZ, evals, &
    orfac, evecs_data, 1, 1, evecs_ScaLAPACK_obj%desc, work, lwork, &
    iwork, liwork, ifail, icluster, gap, info)
! call print("done pdsygvx", PRINT_ALWAYS)

  if (info .ne. 0) then
      call print("ScaLAPACK_diagonlise_gen_r got info " // info, PRINT_ALWAYS)
      if (mod(info,2) .ne. 0) then
	  call print("   eigenvectors failed to converge", PRINT_ALWAYS)
      else if (mod(info/2,2) .ne. 0) then
	  call print("   eigenvectors failed to orthogonalize", PRINT_ALWAYS)
	  do i=1, size(icluster), 2
	      if (icluster(i) /= 0) then
	          call print (" eval cluster " // icluster(i) // " " // icluster(i+1) // &
		    "("// (icluster(i+1)-icluster(i)+1) // ")", PRINT_ALWAYS)
	      else
	          exit
	      endif
	  end do
	  do i=1, this%N_R
	      call print(" eigenvalue " // i // " " // evals(i), PRINT_ANAL)
	  end do
      else if (mod(info/4,2) .ne. 0) then
	  call print("   not enough space for all eigenvectors in range", PRINT_ALWAYS)
      else if (mod(info/8,2) .ne. 0) then
	  call print("   failed to compute eigenvalues", PRINT_ALWAYS)
      else if (mod(info/16,2) .ne. 0) then
	  call print("   S was not positive definite "//ifail(1), PRINT_ALWAYS)
      endif
  endif

  call DEALLOC_TRACE("ScaLAPACK_diagonalise_gen_r iwork", size(iwork)*INTEGER_SIZE)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_gen_r work", size(work)*REAL_SIZE)
  deallocate (iwork, work)

  deallocate (icluster, ifail, gap)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_gen_r data_copy", size(data_copy)*REAL_SIZE)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_gen_r overlap_data_copy", size(overlap_data_copy)*REAL_SIZE)
  deallocate (data_copy, overlap_data_copy)

  if (M /= this%N_R .or. NZ /= this%N_R) then ! bad enough error to report
    RAISE_ERROR("ScaLAPACK_diagonalise_r: Failed to diagonalise info="//info, error)
  endif
#endif

end subroutine ScaLAPACK_diagonalise_gen_r

subroutine ScaLAPACK_diagonalise_gen_c(this, data, overlap_ScaLAPACK_obj, overlap_data, &
  evals, evecs_ScaLAPACK_obj, evecs_data, error)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this
  complex(dp), intent(in) :: data(:,:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: overlap_ScaLAPACK_obj
  complex(dp), intent(in) :: overlap_data(:,:)
  real(dp), intent(out) :: evals(:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: evecs_ScaLAPACK_obj
  complex(dp), intent(out) :: evecs_data(:,:)
  integer, intent(out), optional :: error

  integer liwork, lrwork, lcwork
  integer, allocatable :: iwork(:)
  real(dp), allocatable :: rwork(:)
  complex(dp), allocatable :: cwork(:)
  integer, allocatable :: icluster(:), ifail(:)
  real(dp), allocatable :: gap(:)

  complex(dp), allocatable :: data_copy(:,:), overlap_data_copy(:,:)

  integer info
  integer M, NZ, i
  real(dp) :: orfac = 1.0e-4_dp

  INIT_ERROR(error)

  evals = 0.0_dp
  evecs_data = 0.0_dp

#ifdef SCALAPACK
  allocate(icluster(2*this%ScaLAPACK_obj%n_proc_rows*this%ScaLAPACK_obj%n_proc_cols))
  allocate(gap(this%ScaLAPACK_obj%n_proc_rows*this%ScaLAPACK_obj%n_proc_cols))
  allocate(ifail(this%N_R))

  allocate(data_copy(this%l_N_R, this%l_N_C))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_gen_r data_copy", size(data_copy)*COMPLEX_SIZE)
  allocate(overlap_data_copy(this%l_N_R, this%l_N_C))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_gen_r overlap_data_copy", size(overlap_data_copy)*COMPLEX_SIZE)
  data_copy = data
  overlap_data_copy = overlap_data

  lrwork = -1
  lcwork = -1
  liwork = -1
  allocate(rwork(1))
  allocate(cwork(1))
  allocate(iwork(1))
  call pzhegvx(1, 'V', 'A', 'U', this%N_R, &
    data_copy, 1, 1, this%desc, &
    overlap_data_copy, 1, 1, overlap_ScaLAPACK_obj%desc, &
    0.0_dp, 0.0_dp, 1, 1, -1.0_dp, M, NZ, evals, &
    orfac, evecs_data, 1, 1, evecs_ScaLAPACK_obj%desc, cwork, lcwork, &
    rwork, lrwork, iwork, liwork, ifail, icluster, gap, info)
  lrwork=rwork(1)*4
  lcwork=cwork(1)
  liwork=iwork(1)
  deallocate(rwork)
  deallocate(cwork)
  deallocate(iwork)
  allocate(rwork(lrwork))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_gen_r rwork", size(rwork)*REAL_SIZE)
  allocate(cwork(lcwork))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_gen_r cwork", size(cwork)*COMPLEX_SIZE)
  allocate(iwork(liwork))
  call ALLOC_TRACE("ScaLAPACK_diagonalise_gen_r iwork", size(iwork)*INTEGER_SIZE)

  call pzhegvx(1, 'V', 'A', 'U', this%N_R, &
    data_copy, 1, 1, this%desc, &
    overlap_data_copy, 1, 1, overlap_ScaLAPACK_obj%desc, &
    0.0_dp, 0.0_dp, 1, 1, -1.0_dp, M, NZ, evals, &
    orfac, evecs_data, 1, 1, evecs_ScaLAPACK_obj%desc, cwork, lcwork, &
    rwork, lrwork, iwork, liwork, ifail, icluster, gap, info)

  if (info .ne. 0) then
      call print("ScaLAPACK_diagonlise_gen_c got info " // info, PRINT_ALWAYS)
      if (mod(info,2) .ne. 0) then
	  call print("   eigenvectors failed to converge", PRINT_ALWAYS)
      else if (mod(info/2,2) .ne. 0) then
          call print("   eigenvectors failed to orthogonalize", PRINT_ALWAYS)
	  do i=1, size(icluster), 2
	      if (icluster(i) /= 0) then
	          call print (" eval cluster " // icluster(i) // " " // icluster(i+1) // &
		    "("// (icluster(i+1)-icluster(i)+1) // ")", PRINT_ALWAYS)
	      else
	          exit
	      endif
	  end do
          do i=1, this%N_R
            call print(" eigenvalue " // i // " " // evals(i), PRINT_ANAL)
          end do
      else if (mod(info/4,2) .ne. 0) then
	  call print("   not enough space for all eigenvectors in range", PRINT_ALWAYS)
      else if (mod(info/8,2) .ne. 0) then
	  call print("   failed to compute eigenvalues", PRINT_ALWAYS)
      else if (mod(info/16,2) .ne. 0) then
	  call print("   S was not positive definite "//ifail(1), PRINT_ALWAYS)
      endif
  endif

  deallocate (icluster, ifail, gap)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_gen_r iwork", size(iwork)*INTEGER_SIZE)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_gen_r rwork", size(rwork)*REAL_SIZE)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_gen_r cwork", size(cwork)*COMPLEX_SIZE)
  deallocate (iwork, rwork, cwork)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_gen_r data_copy", size(data_copy)*COMPLEX_SIZE)
  call DEALLOC_TRACE("ScaLAPACK_diagonalise_gen_r overlap_data_copy", size(overlap_data_copy)*COMPLEX_SIZE)
  deallocate (data_copy, overlap_data_copy)

  if (M /= this%N_R .or. NZ /= this%N_R) then ! bad enough error to report
    RAISE_ERROR("ScaLAPACK_diagonalise_r: Failed to diagonalise info="//info, error)
  endif
#endif

end subroutine ScaLAPACK_diagonalise_gen_c

subroutine ScaLAPACK_add_identity_r(this, data)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this
  real(dp), intent(inout) :: data(:,:)

#ifdef SCALAPACK
  integer g_i
  integer l_i, l_j, p_i, p_j


  do g_i=1, this%N_R
    call coords_global_to_local(this, g_i, g_i, l_i, l_j, p_i, p_j)
    if (p_i == this%ScaLAPACK_obj%my_proc_row .and. &
        p_j == this%ScaLAPACK_obj%my_proc_col) &
      data(l_i, l_j) = data(l_i, l_j) + 1.0_dp
  end do
#endif
end subroutine ScaLAPACK_add_identity_r

subroutine ScaLAPACK_Matrix_d_print(this, data, file)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this
  real(dp), intent(in) :: data(:,:)
  type(Inoutput), intent(inout), optional :: file

#ifdef SCALAPACK
include 'mpif.h'
  integer l_i, l_j, g_i, g_j
  character(len=200), allocatable :: lines(:)

  if (this%l_N_R*this%l_N_c > 0) then
    allocate(lines(this%l_N_R*this%l_N_C))

    do l_i=1, this%l_N_R
    do l_j=1, this%l_N_C
      call coords_local_to_global(this, l_i, l_j, g_i, g_j)
      lines((l_i-1)*this%l_N_C + l_j) = &
	"ScaLAPACK local_matrix li,j " // l_i // " " // l_j // " gi,j " // g_i // " " // g_j &
	// " " // data(l_i,l_j)
    end do
    end do
  else
    allocate(lines(1))
    lines(1) = ""
  endif

  call mpi_print(this%ScaLAPACK_obj%mpi_obj, lines)

  deallocate(lines)
#endif
end subroutine ScaLAPACK_Matrix_d_print

subroutine ScaLAPACK_Matrix_z_print(this, data, file)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this
  complex(dp), intent(in) :: data(:,:)
  type(Inoutput), intent(inout), optional :: file

#ifdef SCALAPACK
  integer l_i, l_j, g_i, g_j
  character(len=200), allocatable :: lines(:)

  if (this%l_N_R*this%l_N_C > 0) then
    allocate(lines(this%l_N_R*this%l_N_C))

    do l_i=1, this%l_N_R
    do l_j=1, this%l_N_C
      call coords_local_to_global(this, l_i, l_j, g_i, g_j)
      lines((l_i-1)*this%l_N_C + l_j) = &
	"ScaLAPACK local_matrix li,j " // l_i // " " // l_j // " gi,j " // g_i // " " // g_j &
	// " " // data(l_i,l_j)
    end do
    end do
  else
    allocate(lines(1))
    lines(1) = ""
  endif

  call mpi_print(this%ScaLAPACK_obj%mpi_obj, lines)

  deallocate(lines)
#endif
end subroutine ScaLAPACK_Matrix_z_print

subroutine ScaLAPACK_matrix_product_sub_ddd(c_scalapack, c_data, a_scalapack, a_data, b_scalapack, b_data, &
  a_transpose, b_transpose, a_conjugate, b_conjugate)
  type(Matrix_ScaLAPACK_Info), intent(in) :: c_scalapack
  real(dp), intent(inout) :: c_data(:,:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: a_scalapack
  real(dp), intent(in) :: a_data(:,:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: b_scalapack
  real(dp), intent(in) :: b_data(:,:)
  logical, intent(in), optional :: a_transpose, b_transpose, a_conjugate, b_conjugate

#ifdef SCALAPACK 

  logical a_transp, b_transp, a_conjg, b_conjg
  character a_op, b_op

  a_transp = .false.
  b_transp = .false.
  a_conjg = .false.
  b_conjg = .false.
  if (present(a_transpose)) a_transp = a_transpose
  if (present(b_transpose)) b_transp = b_transpose
  if (present(a_conjugate)) a_conjg = a_conjugate
  if (present(b_conjugate)) b_conjg = b_conjugate

  if (a_transp .or. a_conjg) then
    a_op = 'T'
  else
    a_op = 'N'
  endif
  if (b_transp .or. b_conjg) then
    b_op = 'T'
  else
    b_op = 'N'
  endif

  call pdgemm(a_op, b_op, c_scalapack%N_R, c_scalapack%N_C, a_scalapack%N_C, &
    1.0_dp, a_data, 1, 1, a_scalapack%desc, b_data, 1, 1, b_scalapack%desc, &
    0.0_dp, c_data, 1, 1, c_scalapack%desc)

#endif

end subroutine ScaLAPACK_matrix_product_sub_ddd

subroutine ScaLAPACK_matrix_product_sub_zzz(c_scalapack, c_data, a_scalapack, a_data, b_scalapack, b_data, &
  a_transpose, b_transpose, a_conjugate, b_conjugate)
  type(Matrix_ScaLAPACK_Info), intent(in) :: c_scalapack
  complex(dp), intent(inout) :: c_data(:,:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: a_scalapack
  complex(dp), intent(in) :: a_data(:,:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: b_scalapack
  complex(dp), intent(in) :: b_data(:,:)
  logical, intent(in), optional :: a_transpose, b_transpose, a_conjugate, b_conjugate

#ifdef SCALAPACK
  logical a_transp, b_transp, a_conjg, b_conjg
  character a_op, b_op

  a_transp = .false.
  b_transp = .false.
  a_conjg = .false.
  b_conjg = .false.
  if (present(a_transpose)) a_transp = a_transpose
  if (present(b_transpose)) b_transp = b_transpose
  if (present(a_conjugate)) a_conjg = a_conjugate
  if (present(b_conjugate)) b_conjg = b_conjugate

  if (a_transp .and. a_conjg) call system_abort("ScaLAPACK_matrix_product_sub_zzz called with a_transp and a_conjg")
  if (b_transp .and. b_conjg) call system_abort("ScaLAPACK_matrix_product_sub_zzz called with b_transp and b_conjg")

  if (a_transp) then
    a_op = 'T'
  else if (a_conjg) then
    a_op = 'C'
  else
    a_op = 'N'
  endif
  if (b_transp) then
    b_op = 'T'
  else if (b_conjg) then
    b_op = 'C'
  else
    b_op = 'N'
  endif

  call pzgemm(a_op, b_op, c_scalapack%N_R, c_scalapack%N_C, a_scalapack%N_C, &
    cmplx(1.0_dp,0.0_dp,dp), a_data, 1, 1, a_scalapack%desc, b_data, 1, 1, b_scalapack%desc, &
    cmplx(0.0_dp,0.0_dp,dp), c_data, 1, 1, c_scalapack%desc)
#endif

end subroutine ScaLAPACK_matrix_product_sub_zzz

subroutine ScaLAPACK_matrix_product_vect_asdiagonal_sub_ddd(this_scalapack, this_data, a_scalapack, a_data, diag)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this_scalapack
  real(dp), intent(out) :: this_data(:,:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: a_scalapack
  real(dp), intent(in) :: a_data(:,:)
  real(dp), intent(in) :: diag(:) 

#ifdef SCALAPACK
  integer l_i, l_j, g_i, g_j

  do l_i=1, this_scalapack%l_N_R
  do l_j=1, this_scalapack%l_N_C
    call coords_local_to_global(this_scalapack, l_i, l_j, g_i, g_j)
    this_data(l_i,l_j) = a_data(l_i,l_j) * diag(g_j)
  end do
  end do
#else
  this_data = 0.0_dp
#endif
end subroutine ScaLAPACK_matrix_product_vect_asdiagonal_sub_ddd

subroutine ScaLAPACK_matrix_product_vect_asdiagonal_sub_zzd(this_scalapack, this_data, a_scalapack, a_data, diag)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this_scalapack
  complex(dp), intent(out) :: this_data(:,:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: a_scalapack
  complex(dp), intent(in) :: a_data(:,:)
  real(dp), intent(in) :: diag(:) 

#ifdef SCALAPACK
  integer l_i, l_j, g_i, g_j

  do l_i=1, this_scalapack%l_N_R
  do l_j=1, this_scalapack%l_N_C
    call coords_local_to_global(this_scalapack, l_i, l_j, g_i, g_j)
    this_data(l_i,l_j) = a_data(l_i,l_j) * diag(g_j)
  end do
  end do
#else
  this_data = 0.0_dp
#endif
end subroutine ScaLAPACK_matrix_product_vect_asdiagonal_sub_zzd

subroutine ScaLAPACK_matrix_product_vect_asdiagonal_sub_zzz(this_scalapack, this_data, a_scalapack, a_data, diag)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this_scalapack
  complex(dp), intent(out) :: this_data(:,:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: a_scalapack
  complex(dp), intent(in) :: a_data(:,:)
  complex(dp), intent(in) :: diag(:) 

#ifdef SCALAPACK
  integer l_i, l_j, g_i, g_j

  do l_i=1, this_scalapack%l_N_R
  do l_j=1, this_scalapack%l_N_C
    call coords_local_to_global(this_scalapack, l_i, l_j, g_i, g_j)
    this_data(l_i,l_j) = a_data(l_i,l_j) * diag(g_j)
  end do
  end do
#else
  this_data = 0.0_dp
#endif
end subroutine ScaLAPACK_matrix_product_vect_asdiagonal_sub_zzz

function ScaLAPACK_Re_diagZ(this_scalapack, this_data)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this_scalapack
  complex(dp), intent(in) :: this_data(:,:)
  real(dp) :: ScaLAPACK_Re_diagZ(this_scalapack%N_R)

  integer :: g_i, l_i, l_j, p_i, p_j
  ScaLAPACK_Re_diagZ = 0.0_dp

#ifdef SCALAPACK
  do g_i=1, this_scalapack%N_R
    call coords_global_to_local(this_scalapack, g_i, g_i, l_i, l_j, p_i, p_j)
    if (p_i == this_scalapack%ScaLAPACK_obj%my_proc_row .and. &
        p_j == this_scalapack%ScaLAPACK_obj%my_proc_col) &
      ScaLAPACK_Re_diagZ(g_i) = real(this_data(l_i, l_j))
  end do
#endif

end function ScaLAPACK_Re_diagZ

function ScaLAPACK_Re_diagD(this_scalapack, this_data)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this_scalapack
  real(dp), intent(in) :: this_data(:,:)
  real(dp) :: ScaLAPACK_Re_diagD(this_scalapack%N_R)

  integer :: g_i, l_i, l_j, p_i, p_j

  ScaLAPACK_Re_diagD = 0.0_dp

#ifdef SCALAPACK
  do g_i=1, this_scalapack%N_R
    call coords_global_to_local(this_scalapack, g_i, g_i, l_i, l_j, p_i, p_j)
    if (p_i == this_scalapack%ScaLAPACK_obj%my_proc_row .and. &
        p_j == this_scalapack%ScaLAPACK_obj%my_proc_col) &
      ScaLAPACK_Re_diagD(g_i) = this_data(l_i, l_j)
  end do
#endif

end function ScaLAPACK_Re_diagD

function ScaLAPACK_diag_spinorZ(this_scalapack, this_data)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this_scalapack
  complex(dp), intent(in) :: this_data(:,:)
  complex(dp) :: ScaLAPACK_diag_spinorZ(2,2,this_scalapack%N_R/2)

  integer :: g_i, l_i, l_j, p_i, p_j
  ScaLAPACK_diag_spinorZ = 0.0_dp

#ifdef SCALAPACK
  do g_i=1, this_scalapack%N_R,2
    call coords_global_to_local(this_scalapack, g_i, g_i, l_i, l_j, p_i, p_j)
    if (p_i == this_scalapack%ScaLAPACK_obj%my_proc_row .and. &
        p_j == this_scalapack%ScaLAPACK_obj%my_proc_col) &
      ScaLAPACK_diag_spinorZ(1,1,(g_i-1)/2+1) = this_data(l_i, l_j)
    call coords_global_to_local(this_scalapack, g_i+1, g_i, l_i, l_j, p_i, p_j)
    if (p_i == this_scalapack%ScaLAPACK_obj%my_proc_row .and. &
        p_j == this_scalapack%ScaLAPACK_obj%my_proc_col) &
      ScaLAPACK_diag_spinorZ(2,1,(g_i-1)/2+1) = this_data(l_i, l_j)
    call coords_global_to_local(this_scalapack, g_i, g_i+1, l_i, l_j, p_i, p_j)
    if (p_i == this_scalapack%ScaLAPACK_obj%my_proc_row .and. &
        p_j == this_scalapack%ScaLAPACK_obj%my_proc_col) &
      ScaLAPACK_diag_spinorZ(1,2,(g_i-1)/2+1) = this_data(l_i, l_j)
    call coords_global_to_local(this_scalapack, g_i+1, g_i+1, l_i, l_j, p_i, p_j)
    if (p_i == this_scalapack%ScaLAPACK_obj%my_proc_row .and. &
        p_j == this_scalapack%ScaLAPACK_obj%my_proc_col) &
      ScaLAPACK_diag_spinorZ(2,2,(g_i-1)/2+1) = this_data(l_i, l_j)
  end do
#endif

end function ScaLAPACK_diag_spinorZ

function ScaLAPACK_diag_spinorD(this_scalapack, this_data)
  type(Matrix_ScaLAPACK_Info), intent(in) :: this_scalapack
  real(dp), intent(in) :: this_data(:,:)
  complex(dp) :: ScaLAPACK_diag_spinorD(2,2,this_scalapack%N_R/2)

  integer :: g_i, l_i, l_j, p_i, p_j

  ScaLAPACK_diag_spinorD = 0.0_dp

#ifdef SCALAPACK
  do g_i=1, this_scalapack%N_R, 2
    call coords_global_to_local(this_scalapack, g_i, g_i, l_i, l_j, p_i, p_j)
    if (p_i == this_scalapack%ScaLAPACK_obj%my_proc_row .and. &
        p_j == this_scalapack%ScaLAPACK_obj%my_proc_col) &
      ScaLAPACK_diag_spinorD(1,1,(g_i-1)/2+1) = this_data(l_i, l_j)
    call coords_global_to_local(this_scalapack, g_i+1, g_i, l_i, l_j, p_i, p_j)
    if (p_i == this_scalapack%ScaLAPACK_obj%my_proc_row .and. &
        p_j == this_scalapack%ScaLAPACK_obj%my_proc_col) &
      ScaLAPACK_diag_spinorD(2,1,(g_i-1)/2+1) = this_data(l_i, l_j)
    call coords_global_to_local(this_scalapack, g_i, g_i+1, l_i, l_j, p_i, p_j)
    if (p_i == this_scalapack%ScaLAPACK_obj%my_proc_row .and. &
        p_j == this_scalapack%ScaLAPACK_obj%my_proc_col) &
      ScaLAPACK_diag_spinorD(1,2,(g_i-1)/2+1) = this_data(l_i, l_j)
    call coords_global_to_local(this_scalapack, g_i+1, g_i+1, l_i, l_j, p_i, p_j)
    if (p_i == this_scalapack%ScaLAPACK_obj%my_proc_row .and. &
        p_j == this_scalapack%ScaLAPACK_obj%my_proc_col) &
      ScaLAPACK_diag_spinorD(2,2,(g_i-1)/2+1) = this_data(l_i, l_j)
  end do
#endif

end function ScaLAPACK_diag_spinorD


end module ScaLAPACK_module
