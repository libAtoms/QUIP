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
!X IP module
!X
!% Contain objects for handling TB matrices (\texttt{type TBMatrix}) and vectors
!% (\texttt{type TBVector}).
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module TBMatrix_module

use error_module
use system_module, only : dp, current_verbosity, inoutput, PRINT_NORMAL, operator(//), optional_default
use mpi_context_module
use linearalgebra_module
use atoms_module

use ScaLAPACK_module
use Matrix_module
use RS_SparseMatrix_module
use TBModel_module

implicit none
private

public :: TBVector
type TBVector
  integer :: N = 0
  integer :: n_vectors = 0
  logical :: is_complex = .false.  !% Logical variable to define if vector elements are complex

  real(dp), allocatable :: data_d(:,:)
  complex(dp), allocatable :: data_z(:,:)
end type TBVector

public :: TBMatrix
type TBMatrix
  integer :: N = 0
  integer :: n_matrices = 0
  logical :: is_complex = .false.   !% Logical variable to define if matrix elements are complex
  logical :: is_sparse = .false.    !% Logical variable to define if such matrix is sparse

  type(MatrixD), allocatable :: data_d(:)
  type(MatrixZ), allocatable :: data_z(:)
  type(RS_SparseMatrixD), allocatable :: sdata_d(:)
  type(RS_SparseMatrixZ), allocatable :: sdata_z(:)
end type TBMatrix

public :: Initialise
interface Initialise
  module procedure TBMatrix_Initialise, TBMatrix_Initialise_tbm, TBVector_Initialise
  module procedure TBMatrix_Initialise_sp
end interface Initialise

public :: Finalise
interface Finalise
  module procedure TBMatrix_Finalise, TBVector_Finalise
end interface Finalise

public :: Wipe
interface Wipe
  module procedure TBMatrix_Wipe, TBVector_Wipe
end interface Wipe

public :: Zero
interface Zero
  module procedure TBMatrix_Zero, TBVector_Zero
end interface Zero

public :: Print
interface Print
  module procedure TBMatrix_Print, TBVector_Print
end interface Print

public :: copy
interface copy
  module procedure TBMatrix_copy_d
  module procedure TBMatrix_copy_z
end interface copy

public :: add_block
interface add_block
  module procedure TBMatrix_add_block_d, TBMatrix_add_block_z
end interface add_block

public :: diagonalise
interface diagonalise
  module procedure TBMatrix_diagonalise, TBmatrix_diagonalise_gen
end interface diagonalise

public :: multDiag
interface multDiag
  module procedure TBMatrix_multDiag, TBMatrix_multDiag_d, TBMatrix_multDiag_z
end interface multDiag

public :: multDiagRL
interface multDiagRL
  module procedure TBMatrix_multDiagRL_d
end interface multDiagRL

public :: matrix_product_sub
interface matrix_product_sub
  module procedure TBMatrix_matrix_product_sub, TBMatrix_matrix_product_sub_r2
end interface matrix_product_sub

public operator(.DOT.)
interface operator(.DOT.)
  module procedure TBVector_dotproduct
end interface

public :: partial_TraceMult
interface partial_TraceMult
  module procedure TBMatrix_partial_TraceMult
end interface partial_TraceMult

public :: partial_TraceMult_spinor
interface partial_TraceMult_spinor
  module procedure TBMatrix_partial_TraceMult_spinor
end interface partial_TraceMult_spinor

public :: TraceMult
interface TraceMult
  module procedure TBMatrix_TraceMult
end interface TraceMult

public :: Re_diag
interface Re_diag
  module procedure TBMatrix_Re_diag
end interface Re_diag

public :: diag_spinor
interface diag_spinor
  module procedure TBMatrix_diag_spinor
end interface diag_spinor

public :: scaled_sum
interface scaled_sum
  module procedure TBMatrix_scaled_sum
end interface scaled_sum

public :: scaled_accum
interface scaled_accum
  module procedure TBMatrix_scaled_accum
end interface scaled_accum

public :: inverse
interface inverse
  module procedure TBMatrix_inverse
end interface inverse

public :: sum_in_place
interface sum_in_place
  module procedure TBMatrix_sum_in_place
end interface sum_in_place

public :: accum_scaled_elem_product
interface accum_scaled_elem_product
  module procedure TBMatrix_accum_scaled_elem_product
end interface accum_scaled_elem_product

public :: sum_matrices
interface sum_matrices
  module procedure TBMatrix_sum_matrices_d
end interface sum_matrices

public :: transpose_sub
interface transpose_sub
  module procedure TBMatrix_transpose_sub
end interface transpose_sub

contains

subroutine TBMatrix_sum_in_place(this, mpi)
  type(TBMatrix), intent(inout) :: this
  type(MPI_context) :: mpi

  integer i

  do i=1, this%n_matrices
    if (this%is_complex) then
      if (this%is_sparse) then
	call sum_in_place(mpi, this%sdata_z(i)%data)
      else
	call sum_in_place(mpi, this%data_z(i)%data)
      endif
    else
      if (this%is_sparse) then
	call sum_in_place(mpi, this%sdata_d(i)%data)
      else
	call sum_in_place(mpi, this%data_d(i)%data)
      endif
    endif
  end do

end subroutine TBMatrix_sum_in_place

subroutine TBMatrix_Initialise(this, N, n_matrices, is_complex, scalapack_obj)
  type(TBMatrix), intent(inout) :: this
  integer, intent(in) :: N
  integer, intent(in), optional :: n_matrices
  logical, intent(in), optional :: is_complex
  type(ScaLAPACK), intent(in), optional :: scalapack_obj

  integer i

  call Finalise(this)

  this%N = N
  if (present(n_matrices)) then
    this%n_matrices = n_matrices
  else
    this%n_matrices = 1
  endif
  if (present(is_complex)) then
    this%is_complex = is_complex
  endif

  if (this%is_complex) then
    allocate(this%data_z(this%n_matrices))
    do i=1, this%n_matrices
      call Initialise(this%data_z(i), this%N, scalapack_obj=scalapack_obj)
    end do
  else
    allocate(this%data_d(this%n_matrices))
    do i=1, this%n_matrices
      call Initialise(this%data_d(i), this%N, scalapack_obj=scalapack_obj)
    end do
  endif

end subroutine TBMatrix_Initialise

subroutine TBMatrix_Initialise_tbm(this, from)
  type(TBMatrix), intent(inout) :: this
  type(TBMatrix), intent(in) :: from

  integer i

  call Finalise(this)

  this%N = from%N
  this%n_matrices = from%n_matrices
  this%is_complex = from%is_complex

  if (this%is_complex) then
    if (this%n_matrices > 0) allocate(this%data_z(this%n_matrices))
    do i=1, this%n_matrices
      call Initialise(this%data_z(i), from%data_z(i))
    end do
  else
    if (this%n_matrices > 0) allocate(this%data_d(this%n_matrices))
    do i=1, this%n_matrices
      call Initialise(this%data_d(i), from%data_d(i))
    end do
  endif

end subroutine TBMatrix_Initialise_tbm

subroutine TBMatrix_Initialise_sp(this, at, first_orb_of_atom, n_matrices, is_complex, mpi_obj)
  type(TBMatrix), intent(inout) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: first_orb_of_atom(:)
  integer, intent(in), optional :: n_matrices
  logical, intent(in), optional :: is_complex
  type(MPI_Context), intent(in), optional :: mpi_obj

  integer i

  call Finalise(this)

  this%is_sparse = .true.

  this%N = first_orb_of_atom(at%N+1)-1
  if (present(n_matrices)) then
    this%n_matrices = n_matrices
  else
    this%n_matrices = 1
  endif
  if (present(is_complex)) then
    this%is_complex = is_complex
  endif

  if (this%is_complex) then
    if (this%n_matrices > 0) allocate(this%sdata_z(this%n_matrices))
    do i=1, this%n_matrices
      call Initialise(this%sdata_z(i), at, first_orb_of_atom, mpi_obj=mpi_obj)
    end do
  else
    if (this%n_matrices > 0) allocate(this%sdata_d(this%n_matrices))
    do i=1, this%n_matrices
      call Initialise(this%sdata_d(i), at, first_orb_of_atom, mpi_obj=mpi_obj)
    end do
  endif

end subroutine TBMatrix_Initialise_sp

subroutine TBMatrix_Finalise(this)
  type(TBMatrix), intent(inout) :: this

  call Wipe(this)

end subroutine TBMatrix_Finalise

subroutine TBMatrix_Zero(this, d_mask, od_mask)
  type(TBMatrix), intent(inout) :: this
  logical, intent(in), optional :: d_mask(:), od_mask(:)

  integer i

  if (allocated(this%data_d)) then
    do i=1, size(this%data_d)
      call Zero(this%data_d(i), d_mask, od_mask)
    end do
  end if
  if (allocated(this%data_z)) then
    do i=1, size(this%data_z)
      call Zero(this%data_z(i), d_mask, od_mask)
    end do
  end if
  if (allocated(this%sdata_d)) then
    do i=1, size(this%sdata_d)
      call Zero(this%sdata_d(i), d_mask, od_mask)
    end do
  end if
  if (allocated(this%sdata_z)) then
    do i=1, size(this%sdata_z)
      call Zero(this%sdata_z(i), d_mask, od_mask)
    end do
  end if
end subroutine TBMatrix_Zero

subroutine TBMatrix_Wipe(this)
  type(TBMatrix), intent(inout) :: this

  integer :: i

  if (allocated(this%data_d)) then
    do i=1, size(this%data_d)
      call Finalise(this%data_d(i))
    end do
    deallocate(this%data_d)
  endif
  if (allocated(this%data_z)) then
    do i=1, size(this%data_z)
      call Finalise(this%data_z(i))
    end do
    deallocate(this%data_z)
  end if
  if (allocated(this%sdata_d)) then
    do i=1, size(this%sdata_d)
      call Finalise(this%sdata_d(i))
    end do
    deallocate(this%sdata_d)
  endif
  if (allocated(this%sdata_z)) then
    do i=1, size(this%sdata_z)
      call Finalise(this%sdata_z(i))
    end do
    deallocate(this%sdata_z)
  endif

  this%N = 0
  this%n_matrices = 0
end subroutine TBMatrix_Wipe

subroutine TBMatrix_Print(this,file)
  type(TBMatrix),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  integer::i

  if (current_verbosity() < PRINT_NORMAL) return

  call Print('TBMatrix : ', file=file)

  call Print ('TBMatrix : N n_matrices ' // this%N // " " // this%n_matrices, file=file)
  call Print ('TBMatrix : is_complex ' // this%is_complex, file=file)

  if (allocated(this%data_d)) then
    do i=1, size(this%data_d)
      call Print ('TBMatrix : data_d ' // i, file=file)
      call Print(this%data_d(i), file=file)
    end do
  endif

  if (allocated(this%data_z)) then
    do i=1, size(this%data_z)
      call Print ('TBMatrix : data_z ' // i, file=file)
      call Print(this%data_z(i), file=file)
    end do
  endif

  if (allocated(this%sdata_d)) then
    do i=1, size(this%sdata_d)
      call Print ('TBMatrix : sdata_d ' // i, file=file)
      call Print(this%sdata_d(i), file=file)
    end do
  endif

  if (allocated(this%sdata_z)) then
    do i=1, size(this%sdata_z)
      call Print ('TBMatrix : sdata_z ' // i, file=file)
      call Print(this%sdata_z(i), file=file)
    end do
  endif

end subroutine  TBMatrix_Print

subroutine TBMatrix_copy_d(this, data_d, index)
  type(TBMatrix),    intent(in)           :: this
  real(dp), intent(inout), dimension(:,:) :: data_d
  integer, optional, intent(in) :: index
  integer my_index

  my_index = optional_default(1, index)

  call Print('TBMatrix : ')
  call Print ('TBMatrix : N n_matrices ' // this%N // " " // this%n_matrices)
  call Print ('TBMatrix : is_complex ' // this%is_complex)

  if (allocated(this%data_d)) then
    if (my_index > size(this%data_d)) call system_abort("index > size(data_d)")
    data_d(:,:) = this%data_d(my_index)%data
  else if (allocated(this%sdata_d)) then
    if (my_index > size(this%sdata_d)) call system_abort("index > size(sdata_d)")
    call copy(this%sdata_d(my_index), data_d)
  endif
end subroutine TBmatrix_copy_d

subroutine TBMatrix_copy_z(this, data_z, index)
  type(TBMatrix),    intent(in)           :: this
  complex(dp), intent(inout), dimension(:,:) :: data_z
  integer, optional, intent(in) :: index
  integer my_index

  my_index = optional_default(1, index)

  call Print('TBMatrix : ')
  call Print ('TBMatrix : N n_matrices ' // this%N // " " // this%n_matrices)
  call Print ('TBMatrix : is_complex ' // this%is_complex)

  if (allocated(this%data_z)) then
    if (my_index > size(this%data_z)) call system_abort("index > size(data_z)")
    if (any(shape(data_z) /= shape(this%data_z(my_index)%data))) call system_abort("data_z size mismatch")
    data_z(:,:) = this%data_z(my_index)%data
  else if (allocated(this%sdata_z)) then
    call copy(this%sdata_z(my_index), data_z)
  endif
end subroutine TBMatrix_copy_z

subroutine TBVector_Initialise(this, N, n_vectors, is_complex)
  type(TBVector), intent(inout) :: this
  integer, intent(in) :: N
  integer, intent(in), optional :: n_vectors
  logical, intent(in), optional :: is_complex

  call Finalise(this)

  this%N = N
  if (present(n_vectors)) then
    this%n_vectors = n_vectors
  else
    this%n_vectors = 1
  endif
  if (present(is_complex)) then
    this%is_complex = is_complex
  endif

  if (this%is_complex) then
    allocate(this%data_z(this%N,this%n_vectors))
  else
    allocate(this%data_d(this%N,this%n_vectors))
  endif

end subroutine TBVector_Initialise

subroutine TBVector_Finalise(this)
  type(TBVector), intent(inout) :: this

  call Wipe(this)

end subroutine TBVector_Finalise

subroutine TBVector_Wipe(this)
  type(TBVector), intent(inout) :: this

  if (allocated(this%data_d)) deallocate(this%data_d)
  if (allocated(this%data_z)) deallocate(this%data_z)

  this%N = 0
  this%n_vectors = 0
end subroutine TBVector_Wipe

subroutine TBVector_Zero(this)
  type(TBVector), intent(inout) :: this

  if (allocated(this%data_d)) this%data_d = 0.0_dp
  if (allocated(this%data_z)) this%data_z = 0.0_dp

end subroutine TBVector_Zero

subroutine TBVector_Print(this,file)
  type(TBVector),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  integer::i

  if (current_verbosity() < PRINT_NORMAL) return

  call Print('TBVector : ', file=file)

  call Print ('TBVector : N n_vectors ' // this%N // ' ' // this%n_vectors, file=file)
  call Print ('TBVector : is_complex ' // this%is_complex, file=file)

  if (allocated(this%data_d)) then
    do i=1, size(this%data_d,2)
      call Print ('TBVector : data_d ' // i, file=file)
      call Print(this%data_d(:,i), file=file)
    end do
  endif

  if (allocated(this%data_z)) then
    do i=1, size(this%data_z,2)
      call Print ('TBVector : data_z ' // i, file=file)
      call Print(this%data_z(:,i), file=file)
    end do
  endif

end subroutine TBVector_Print

subroutine TBMatrix_add_block_d(this, im, block, block_nr, block_nc, first_row, first_col, at_row, at_col)
  type(TBMatrix), intent(inout) :: this
  integer, intent(in) :: im
  real(dp), intent(in) :: block(:,:)
  integer, intent(in) :: block_nr, block_nc
  integer, intent(in) :: first_row, first_col
  integer, intent(in), optional :: at_row, at_col


  if (this%is_complex) then
    if (this%is_sparse) then
      call add_block(this%sdata_z(im), cmplx(block, 0.0_dp, dp), block_nr, block_nc, at_row, at_col)
    else
      call add_block(this%data_z(im), cmplx(block, 0.0_dp, dp), block_nr, block_nc, first_row, first_col)
    endif
  else
    if (this%is_sparse) then
      call add_block(this%sdata_d(im), block, block_nr, block_nc, at_row, at_col)
    else
      call add_block(this%data_d(im), block, block_nr, block_nc, first_row, first_col)
    endif
  endif

end subroutine TBMatrix_add_block_d

subroutine TBMatrix_add_block_z(this, im, block, block_nr, block_nc, first_row, first_col, at_row, at_col)
  type(TBMatrix), intent(inout) :: this
  integer, intent(in) :: im
  complex(dp), intent(in) :: block(:,:)
  integer, intent(in) :: block_nr, block_nc
  integer, intent(in) :: first_row, first_col
  integer, intent(in), optional :: at_row, at_col

  if (this%is_complex) then
    if (this%is_sparse) then
      call add_block(this%sdata_z(im), block, block_nr, block_nc, at_row, at_col)
    else
      call add_block(this%data_z(im), block, block_nr, block_nc, first_row, first_col)
    endif
  else
    call system_abort ("tried to add complex block to real TBMatrix")
  endif

end subroutine TBMatrix_add_block_z

subroutine TBMatrix_diagonalise_gen(this, overlap, evals, evecs, error)
  type(TBMatrix), intent(in) :: this
  type(TBMatrix), intent(in) :: overlap
  type(TBVector), intent(inout) :: evals
  type(TBMatrix), intent(inout), optional :: evecs
  integer, intent(out), optional :: error

  integer i

  INIT_ERROR(error)

  if (this%is_sparse) then
    RAISE_ERROR("can't diagonalize_gen sparse matrix", error)
  endif

  if (this%is_complex) then
    do i=1, this%n_matrices
      if (present(evecs)) then
	call diagonalise(this%data_z(i), overlap%data_z(i), evals%data_d(:,i), evecs%data_z(i), error = error)
      else
	call diagonalise(this%data_z(i), overlap%data_z(i), evals%data_d(:,i), error = error)
      endif
    end do
  else
    do i=1, this%n_matrices
      if (present(evecs)) then
	call diagonalise(this%data_d(i), overlap%data_d(i), evals%data_d(:,i), evecs%data_d(i), error = error)
      else
	call diagonalise(this%data_d(i), overlap%data_d(i), evals%data_d(:,i), error = error)
      endif
    end do
  endif
  PASS_ERROR(error)
end subroutine TBMatrix_diagonalise_gen

subroutine TBMatrix_diagonalise(this, evals, evecs, error)
  type(TBMatrix), intent(in) :: this
  type(TBVector), intent(inout) :: evals
  type(TBMatrix), intent(inout), optional :: evecs
  integer, intent(out), optional :: error

  integer i

  INIT_ERROR(error)

  if (this%is_sparse) then
    RAISE_ERROR("can't diagonalize sparse matrix", error)
  endif

  if (this%is_complex) then
    do i=1, this%n_matrices
      if (present(evecs)) then
	call diagonalise(this%data_z(i), evals%data_d(:,i), evecs%data_z(i), error = error)
      else
	call diagonalise(this%data_z(i), evals%data_d(:,i), error = error)
      endif
    end do
  else
    do i=1, this%n_matrices
      if (present(evecs)) then
	call diagonalise(this%data_d(i), evals%data_d(:,i), evecs%data_d(i), error = error)
      else
	call diagonalise(this%data_d(i), evals%data_d(:,i), error = error)
      endif
    end do
  endif
  PASS_ERROR(error)
end subroutine TBMatrix_diagonalise

function TBVector_dotproduct(this, that)
  type(TBVector), intent(in) :: this, that
  real(dp) :: TBVector_dotproduct

  if (this%is_complex .or. that%is_complex) then
    call system_abort ("TBVector_dotproduct can't handle complex vectors")
  endif

  TBVector_dotproduct = this%data_d .dot. this%data_d
end function TBVector_dotproduct

!subroutine accum_tbmatrix_col_outer_product(this, mati, i, TBV_f)
!  type(TBMatrix), intent(inout) :: this
!  type(TBMatrix), intent(in) :: mati
!  integer, intent(in) :: i
!  type(TBVector), intent(in) :: TBV_f
!
!  integer ik
!
!  if (mati%is_sparse) then
!    call System_abort("can't accum_tbmatrix_col_outer_product of sparse matrix")
!  endif
!
!  do ik=1, mati%n_matrices
!    if (TBV_f%data_d(i,ik) .fne. 0.0_dp) then
!      if (mati%is_complex) then
!	if (this%is_complex) then
!	  call accum_col_outer_product(this%data_z(ik), mati%data_z(ik), i, TBV_f%data_d(i,ik))
!	else
!	  call system_abort("tried to add outer product of columns of real TBMatrix to complex TBMatrix")
!	endif
!      else
!	if (.not. this%is_complex) then
!	  call accum_col_outer_product(this%data_d(ik), mati%data_d(ik), i, TBV_f%data_d(i,ik))
!	else
!	  call system_abort("tried to add outer product of columns of complex TBMatrix to real TBMatrix")
!	endif
!      endif
!    endif
!  end do
!end subroutine

function TBMatrix_TraceMult(a, b, w, a_H, b_H, diag_mask, offdiag_mask)
  type(TBMatrix), intent(in) :: a, b
  real(dp), intent(in), optional, target :: w(:)
  logical, intent(in), optional :: a_H, b_H
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: TBMatrix_TraceMult(a%n_matrices)

  integer im

  do im=1, a%n_matrices
    ! dense-sparse
    if (.not. a%is_complex .and. .not. a%is_sparse .and. .not. b%is_complex .and. b%is_sparse) then
      TBMatrix_TraceMult(im) = TraceMult(a%data_d(im), b%sdata_d(im), w, a_H, b_H, diag_mask, offdiag_mask)
    else if (a%is_complex .and. .not. a%is_sparse .and. b%is_complex .and. b%is_sparse) then
      TBMatrix_TraceMult(im) = TraceMult(a%data_z(im), b%sdata_z(im), w, a_H, b_H, diag_mask, offdiag_mask)
    ! sparse-dense
    else if (.not. a%is_complex .and. a%is_sparse .and. .not. b%is_complex .and. .not. b%is_sparse) then
      TBMatrix_TraceMult(im) = TraceMult(a%sdata_d(im), b%data_d(im), w, a_H, b_H, diag_mask, offdiag_mask)
    else if (.not. a%is_complex .and. a%is_sparse .and. b%is_complex .and. .not. b%is_sparse) then
      TBMatrix_TraceMult(im) = TraceMult(a%sdata_d(im), b%data_z(im), w, a_H, b_H, diag_mask, offdiag_mask)
    else if (a%is_complex .and. a%is_sparse .and. b%is_complex .and. .not. b%is_sparse) then
      TBMatrix_TraceMult(im) = TraceMult(a%sdata_z(im), b%data_z(im), w, a_H, b_H, diag_mask, offdiag_mask)
    else
      call system_abort ("No TBMatrix_TraceMult implemented yet for " // &
          " a%c " // a%is_complex // " a%s " // a%is_sparse // &
          " b%c " // b%is_complex // " b%s " // b%is_sparse )
    endif
  end do

end function TBMatrix_TraceMult

function TBMatrix_partial_TraceMult(a, b, a_H, b_H, diag_mask, offdiag_mask)
  type(TBMatrix), intent(in) :: a, b
  logical, intent(in), optional :: a_H, b_H
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: TBMatrix_partial_TraceMult(a%N, a%n_matrices)

  integer im

  do im=1, a%n_matrices
    if (.not. a%is_complex .and. .not. a%is_sparse .and. .not. b%is_complex .and. .not. b%is_sparse) then
      if (present(diag_mask) .or. present(offdiag_mask)) call system_abort("Can't do diag_mask for TraceMult of 2 dense matrices")
      TBMatrix_partial_TraceMult(:,im) = partial_TraceMult(a%data_d(im), b%data_d(im), a_H, b_H)
    else if (.not. a%is_complex .and. .not. a%is_sparse .and. .not. b%is_complex .and. b%is_sparse) then
      TBMatrix_partial_TraceMult(:,im) = partial_TraceMult(a%data_d(im), b%sdata_d(im), a_H, b_H, diag_mask,offdiag_mask)
    else if (.not. a%is_complex .and. a%is_sparse .and. .not. b%is_complex .and. .not. b%is_sparse) then
      TBMatrix_partial_TraceMult(:,im) = partial_TraceMult(a%sdata_d(im), b%data_d(im), a_H, b_H, diag_mask,offdiag_mask)
    else if (a%is_complex .and. .not. a%is_sparse .and. b%is_complex .and. .not. b%is_sparse) then
      if (present(diag_mask) .or. present(offdiag_mask)) call system_abort("Can't do diag_mask for TraceMult of 2 dense matrices")
      TBMatrix_partial_TraceMult(:,im) = partial_TraceMult(a%data_z(im), b%data_z(im), a_H, b_H)
    else if (a%is_complex .and. .not. a%is_sparse .and. b%is_complex .and. b%is_sparse) then
      TBMatrix_partial_TraceMult(:,im) = partial_TraceMult(a%data_z(im), b%sdata_z(im), a_H, b_H, diag_mask,offdiag_mask)
    else if (a%is_complex .and. a%is_sparse .and. b%is_complex .and. .not. b%is_sparse) then
      TBMatrix_partial_TraceMult(:,im) = partial_TraceMult(a%sdata_z(im), b%data_z(im), a_H, b_H, diag_mask,offdiag_mask)
    else
      call system_abort ("No TBMatrix_partial_TraceMult implemented yet for " // &
        " a%c " // a%is_complex // " a%s " // a%is_sparse // &
        " b%c " // b%is_complex // " b%s " // b%is_sparse )
    endif
  end do

end function TBMatrix_partial_TraceMult

function TBMatrix_partial_TraceMult_spinor(a, b, a_H, b_H, diag_mask, offdiag_mask)
  type(TBMatrix), intent(in) :: a, b
  logical, intent(in), optional :: a_H, b_H
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: TBMatrix_partial_TraceMult_spinor(2,2,a%N/2, a%n_matrices)

  integer im

  do im=1, a%n_matrices
    if (a%is_complex .and. .not. a%is_sparse .and. b%is_complex .and. .not. b%is_sparse) then
      if (present(diag_mask) .or. present(offdiag_mask)) call system_abort("Can't do diag_mask for TraceMult_spinor of 2 dense matrices")
      TBMatrix_partial_TraceMult_spinor(:,:,:,im) = partial_TraceMult_spinor(a%data_z(im), b%data_z(im), a_H, b_H)
    else
      call system_abort ("No TBMatrix_partial_TraceMult_spinor implemented yet for " // &
        " a%c " // a%is_complex // " a%s " // a%is_sparse // &
        " b%c " // b%is_complex // " b%s " // b%is_sparse )
    endif
  end do

end function TBMatrix_partial_TraceMult_spinor


function TBMatrix_Re_diag(a)
  type(TBMatrix), intent(in) :: a
  real(dp) :: TBMatrix_Re_diag(a%N, a%n_matrices)

  integer im

  do im=1, a%n_matrices
    if (a%is_complex) then
      TBMatrix_Re_diag(:,im) = Re_diag(a%data_z(im))
    else
      TBMatrix_Re_diag(:,im) = Re_diag(a%data_d(im))
    endif
  end do

end function TBMatrix_Re_diag

function TBMatrix_diag_spinor(a)
  type(TBMatrix), intent(in) :: a
  complex(dp) :: TBMatrix_diag_spinor(2,2,a%N/2, a%n_matrices)

  integer im

  do im=1, a%n_matrices
    if (a%is_complex) then
      TBMatrix_diag_spinor(:,:,:,im) = diag_spinor(a%data_z(im))
    else
      TBMatrix_diag_spinor(:,:,:,im) = diag_spinor(a%data_d(im))
    endif
  end do

end function TBMatrix_diag_spinor

subroutine TBMatrix_scaled_sum(this, f1_z, m1, f2, m2)
  type(TBMatrix), intent(inout) :: this
  complex(dp), intent(in) :: f1_z
  type(TBMatrix), intent(in) :: m1
  real(dp), intent(in) :: f2
  type(TBMatrix), intent(in) :: m2

  complex(dp) f2_z

  integer im

  if (this%is_sparse) call system_abort("No TBMatrix_scaled_sum for sparse matrices")

  f2_z = f2

  do im=1, this%n_matrices
    if (this%is_complex) then
      if (m1%is_complex) then
	if (m2%is_complex) then
	  call scaled_sum(this%data_z(im), f1_z, m1%data_z(im), f2_z, m2%data_z(im))
	else
	  ! call scaled_sum(this%data_z(im), f1_z, m1%data_z(im), f2_z, m2%data_d(im))
	  call system_abort ("Can't do scaled sum for Z Z D")
	endif
      else
	if (m2%is_complex) then
	  ! call scaled_sum(this%data_z(im), f1_z, m1%data_d(im), f2_z, m2%data_z(im))
	  call system_abort ("Can't do scaled sum for Z D Z")
	else
	  call scaled_sum(this%data_z(im), f1_z, m1%data_d(im), f2_z, m2%data_d(im))
	endif
      endif
    else
      if (m1%is_complex) then
	if (m2%is_complex) then
	  ! call scaled_sum(this%data_d(im), f1_z, m1%data_z(im), f2_z, m2%data_z(im))
	  call system_abort ("Can't do scaled sum for D Z Z")
	else
	  ! call scaled_sum(this%data_d(im), f1_z, m1%data_z(im), f2_z, m2%data_d(im))
	  call system_abort ("Can't do scaled sum for D Z D")
	endif
      else
	if (m2%is_complex) then
	  ! call scaled_sum(this%data_d(im), f1_z, m1%data_d(im), f2_z, m2%data_z(im))
	  call system_abort ("Can't do scaled sum for D D Z")
	else
	  ! call scaled_sum(this%data_d(im), f1_z, m1%data_d(im), f2_z, m2%data_d(im))
	  call system_abort ("Can't do scaled sum for D D D")
	endif
      endif
    endif
  end do

end subroutine TBMatrix_scaled_sum

subroutine TBMatrix_scaled_accum(this, f1_z, m1)
  type(TBMatrix), intent(inout) :: this
  complex(dp), intent(in) :: f1_z
  type(TBMatrix), intent(in) :: m1

  integer im

  if (this%is_sparse) call system_abort("No TBMatrix_accum for sparse matrices")

  do im=1, this%n_matrices
    if (this%is_complex) then
      if (m1%is_complex) then
	call scaled_accum(this%data_z(im), f1_z, m1%data_z(im))
      else
	! call scaled_accum(this%data_z(im), f1_z, m1%data_d(im))
	call system_abort ("Can't do scaled accum for Z D")
      endif
    else
      if (m1%is_complex) then
	call scaled_accum(this%data_d(im), f1_z, m1%data_z(im))
      else
	! call scaled_accum(this%data_d(im), f1_z, m1%data_d(im))
	call system_abort ("Can't do scaled accum for D D")
      endif
    endif
  end do

end subroutine TBMatrix_scaled_accum

subroutine TBMatrix_inverse(this, inv, positive)
  type(TBMatrix), intent(inout) :: this
  type(TBMatrix), intent(out), optional :: inv
  logical, intent(in), optional :: positive

  integer im

  if (this%is_sparse) call system_abort("No TBMatrix_inverse for sparse matrices")

  do im=1, this%n_matrices
    if (this%is_complex) then
      if (present(inv)) then
	if (.not. (inv%is_complex)) then
	    call system_abort("Called TBMatrix_inverse with complex matrix but real inverse")
	endif
	call inverse(this%data_z(im), inv%data_z(im), positive)
      else
	call inverse(this%data_z(im), positive=positive)
      endif
    else
      if (present(inv)) then
	if (inv%is_complex) then
	    call system_abort("Called TBMatrix_inverse with real matrix but complex inverse")
	endif
	call inverse(this%data_d(im), inv%data_d(im), positive)
      else
	call inverse(this%data_d(im), positive=positive)
      endif
    endif
  end do
end subroutine TBMatrix_inverse

subroutine TBMatrix_multDiag(this, A, diag)
  type(TBMatrix), intent(inout) :: this
  type(TBMatrix), intent(in) :: A
  type(TBVector), intent(in) :: diag

  integer im

  if (this%N /= diag%N) then
    call system_abort("Called TBMatrix_multDiag with mismatched sizes")
  endif

  if (this%is_sparse) call system_abort("No TBMatrix_multDiag for sparse matrices")
  if (A%is_sparse) call system_abort("No TBMatrix_multDiag for sparse matrices")

  if ((this%is_complex .and. .not. A%is_complex) .or. &
      (.not. this%is_complex .and.  A%is_complex)) call system_abort ("TBMatrix_multDiag with mismatched types")

  do im=1, this%n_matrices
    if (this%is_complex) then
      if (diag%is_complex) then
	call multDiag(this%data_z(im), A%data_z(im), diag%data_z(:,im))
      else
	call multDiag(this%data_z(im), A%data_z(im), diag%data_d(:,im))
      endif
    else
      if (diag%is_complex) then
	call system_abort("Can't TBMatrix_multDiag of a real matrix time complex diag")
      else
	call multDiag(this%data_d(im), A%data_d(im), diag%data_d(:,im))
      endif
    endif
  end do
end subroutine TBMatrix_multDiag

subroutine TBMatrix_multDiag_d(this, A, diag)
  type(TBMatrix), intent(inout) :: this
  type(TBMatrix), intent(in) :: A
  real(dp), intent(in) :: diag(:)

  integer im

  if (this%N /= size(diag)) then
    call system_abort("Called TBMatrix_multDiag_d with mismatched sizes")
  endif

  if (this%is_sparse) call system_abort("No TBMatrix_multDiag_d for sparse matrices")
  if (A%is_sparse) call system_abort("No TBMatrix_multDiag_d for sparse matrices")

  if ((this%is_complex .and. .not. A%is_complex) .or. &
      (.not. this%is_complex .and.  A%is_complex)) call system_abort ("TBMatrix_multDiag_d with mismatched types")

  do im=1, this%n_matrices
    if (this%is_complex) then
      call multDiag(this%data_z(im), A%data_z(im), diag)
    else
      call multDiag(this%data_d(im), A%data_d(im), diag)
    endif
  end do
end subroutine TBMatrix_multDiag_d


subroutine TBMatrix_multDiag_z(this, A, diag)
  type(TBMatrix), intent(inout) :: this
  type(TBMatrix), intent(in) :: A
  complex(dp), intent(in) :: diag(:)

  integer im

  if (this%N /= size(diag)) then
    call system_abort("Called TBMatrix_multDiag_z with mismatched sizes")
  endif

  if (this%is_sparse) call system_abort("No TBMatrix_multDiag_z for sparse matrices")
  if (A%is_sparse) call system_abort("No TBMatrix_multDiag_z for sparse matrices")

  if ((this%is_complex .and. .not. A%is_complex) .or. &
      (.not. this%is_complex .and.  A%is_complex)) call system_abort ("TBMatrix_multDiag_z with mismatched types")

  do im=1, this%n_matrices
    if (this%is_complex) then
      call multDiag(this%data_z(im), A%data_z(im), diag)
    else
      call system_abort ("TBMatrix_multDiag_z Can't multiply a real matrix by a complex diag")
    endif
  end do
end subroutine TBMatrix_multDiag_z

subroutine TBMatrix_multDiagRL_d(this, A, diag)
  type(TBMatrix), intent(inout) :: this
  type(TBMatrix), intent(in) :: A
  real(dp) :: diag(:)

  integer im

  if (this%N /= size(diag)) &
    call system_abort("Called TBMatrix_multDiagRL_d with mismatched sizes")

  if ((A%is_sparse .and. .not. this%is_sparse) .or. &
      (.not. A%is_sparse .and. this%is_sparse)) &
      call system_abort ("Called TBMatrix_multDiagRL_d with mismatched sparsity")

  if ((this%is_complex .and. .not. A%is_complex) .or. &
      (.not. this%is_complex .and.  A%is_complex)) call system_abort ("TBMatrix_multDiagRL_d with mismatched types")

  do im=1, this%n_matrices
    if (this%is_sparse) then
      if (this%is_complex) then
	call multDiagRL(this%sdata_z(im), A%sdata_z(im), diag)
      else
	call multDiagRL(this%sdata_d(im), A%sdata_d(im), diag)
      endif
    else
      if (this%is_complex) then
	call multDiagRL(this%data_z(im), A%data_z(im), diag)
      else
	call multDiagRL(this%data_d(im), A%data_d(im), diag)
      endif
    endif
  end do
end subroutine TBMatrix_multDiagRL_d

subroutine TBMatrix_matrix_product_sub(C, A, B, A_transpose, A_conjugate, B_transpose, B_conjugate, &
				       diag_mask, offdiag_mask)
  type(TBMatrix), intent(inout) :: C
  type(TBMatrix), intent(in) :: A, B
  logical, intent(in), optional :: A_transpose, A_conjugate, B_transpose, B_conjugate
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)

  logical a_transp, a_conjg, b_transp, b_conjg
  integer im

  if (C%N /= A%N .or. C%N /= B%N) call system_abort ("TBMatrix_matrix_product_sub with C vs A, B size mismatch" &
    // C%N // " " // A%N // " " // B%N )

  if (C%is_sparse .or. A%is_sparse) &
    call system_abort("Can't do TBMatrix_matrix_product_sub C=A*B for sparse C or A matrices")

  a_transp = .false.
  b_transp = .false.
  a_conjg = .false.
  b_conjg = .false.
  if (present(a_transpose)) a_transp = a_transpose
  if (present(b_transpose)) b_transp = b_transpose
  if (present(a_conjugate)) a_conjg = a_conjugate
  if (present(b_conjugate)) b_conjg = b_conjugate

  if (a_transp .and. a_conjg) then
    call system_abort("TBMatrix_matrix_product_sub called with a_transp and a_conj both true")
  endif
  if (b_transp .and. b_conjg) then
    call system_abort("TBMatrix_matrix_product_sub called with b_transp and b_conj both true")
  endif

  do im=1, C%n_matrices
    if (C%is_complex) then
      if (A%is_complex) then
	if (B%is_complex) then
	  if (B%is_sparse) then
	    call matrix_product_sub(C%data_z(im), A%data_z(im), B%sdata_z(im), a_transpose, a_conjugate, b_transpose, b_conjugate, &
				    diag_mask, offdiag_mask)
	  else
	    if (present(diag_mask) .or. present(offdiag_mask)) &
	      call system_abort("TBMatrix_matrix_product_sub can't use masks when B isn't sparse")
	    call matrix_product_sub(C%data_z(im), A%data_z(im), B%data_z(im), a_transpose, a_conjugate, b_transpose, b_conjugate)
	  endif
	else  ! B is real
	  if (B%is_sparse) then
	    call matrix_product_sub(C%data_z(im), A%data_z(im), B%sdata_d(im), a_transpose, a_conjugate, b_transp .or. b_conjg, &
				    diag_mask, offdiag_mask)
	  else
	    if (present(diag_mask) .or. present(offdiag_mask)) &
	      call system_abort("TBMatrix_matrix_product_sub can't use masks when B isn't sparse")
	    call matrix_product_sub(C%data_z(im), A%data_z(im), B%data_d(im), a_transpose, a_conjugate, b_transp .or. b_conjg)
	  endif
	endif
      else ! A is real
	if (B%is_complex) then
	  if (B%is_sparse) then
	    call matrix_product_sub(C%data_z(im), A%data_d(im), B%sdata_z(im), a_transp .or. a_conjg, b_transpose, b_conjugate, &
				    diag_mask, offdiag_mask)
	  else
	    if (present(diag_mask) .or. present(offdiag_mask)) &
	      call system_abort("TBMatrix_matrix_product_sub can't use masks when B isn't sparse")
	    call matrix_product_sub(C%data_z(im), A%data_d(im), B%data_z(im), a_transp .or. a_conjg, b_transpose, b_conjugate)
	  endif
	else
	  call system_abort ("No TBMatrix_matrix_product_sub for C = R R")
	endif
      endif
    else ! C is real
      if (A%is_complex) then
	call system_abort("No TBMatrix_matrix_product for R = C * ?")
      else ! A is real
	if (B%is_complex) then
	  call system_abort("No TBMatrix_matrix_product for R = ?* C")
	else ! B is real
	  if (B%is_sparse) then
	    call matrix_product_sub(C%data_d(im), A%data_d(im), B%sdata_d(im), a_transp .or. a_conjg, b_transp .or. b_conjg, &
				    diag_mask, offdiag_mask)
	  else
	    if (present(diag_mask) .or. present(offdiag_mask)) &
	      call system_abort("TBMatrix_matrix_product_sub can't use masks when B isn't sparse")
	    call matrix_product_sub(C%data_d(im), A%data_d(im), B%data_d(im), a_transp .or. a_conjg, b_transp .or. b_conjg)
	  endif
      endif
      endif
    endif
  end do
end subroutine TBMatrix_matrix_product_sub

subroutine TBMatrix_matrix_product_sub_r2(C, A, B, a_transpose, a_conjugate, b_transpose)
  type(TBMatrix), intent(inout) :: C
  type(TBMatrix), intent(in) :: A
  real(dp), intent(in) :: B(:,:)
  logical, intent(in), optional :: a_transpose, a_conjugate, b_transpose

  logical :: a_transp, a_conjg

  integer im

  a_transp = .false.
  a_conjg = .false.
  if (present(a_transpose)) a_transp = a_transpose
  if (present(a_conjugate)) a_conjg = a_conjugate

  if (a_transp .and. a_conjg) call system_abort("Called TBMatrix_matrix_product_sub_r2 with a_transp and a_conjg")

  if (C%N /= A%N .or. C%N /= size(B,2)) call system_abort ("TBMatrix_matrix_product with C vs A, B size mismatch")

  if (C%is_sparse .or. A%is_sparse) &
    call system_abort("Can't do TBMatrix_matrix_product_r2 for sparse matrices")

  do im=1, C%n_matrices
    if (C%is_complex) then
      call system_abort ("No TBMatrix_matrix_product for C = ? ?")
    else ! C is real
      if (A%is_complex) then
	call system_abort ("No TBMatrix_matrix_product for R = C ?")
      else
	call matrix_product_sub(C%data_d(im), A%data_d(im), B, a_transp .or. a_conjg, b_transpose)
      endif
    endif
  end do
end subroutine TBMatrix_matrix_product_sub_r2

subroutine TBMatrix_accum_scaled_elem_product(A, B, s, C)
  type(TBMatrix), intent(in) :: A, B
  complex(dp), intent(in) :: s
  type(TBMatrix), intent(inout) :: C

  integer im

  if (A%N /= B%N .or. A%N /= C%N) &
    call system_abort ("TBMatrix_accum_scaled_elem_product called with size mismatch")
  if (A%n_matrices /= B%n_matrices .or. A%n_matrices /= C%n_matrices) &
    call system_abort ("TBMatrix_accum_scaled_elem_product called with n_matrices mismatch")

  if (A%is_sparse .or. B%is_sparse .or. C%is_sparse) then
    call system_abort ("TBMatrix_accum_scaled_elem_product called with sparse matrix")
  endif

  do im=1, A%n_matrices
    if (C%is_complex) then
      call system_abort("No TBMatrix_accum_scaled_elem_product for complex C")
    else
      if (A%is_complex) then
	if (B%is_complex) then
	  C%data_d(im)%data = C%data_d(im)%data + real(s * A%data_z(im)%data * B%data_z(im)%data)
	else
	  call system_abort("No TBMatrix_accum_scaled_elem_product for real B")
	endif
      else
	call system_abort("No TBMatrix_accum_scaled_elem_product for real A")
      endif

    endif
  end do
end subroutine TBMatrix_accum_scaled_elem_product

subroutine TBMatrix_sum_matrices_d(this, weights, m)
  type(TBMatrix), intent(in) :: this
  real(dp), intent(in) :: weights(:)
  type(MatrixD), intent(inout) :: m

  integer im

  if (this%N /= m%N .or. this%N /= m%M) call system_abort("TBMatrix_sum_matrices_d called with size mismatch")
  if (this%n_matrices /= size(weights)) call system_abort("TBMatrix_sum_matrices_d called with n_matrices mismatch")

  if (this%is_sparse) call system_abort("Can't do TBMatrix_sum_matrices_d on a sparse TBMatrix")

  m%data = 0.0_dp

  do im=1, this%n_matrices
    if (this%is_complex) then
      m%data = m%data + weights(im)*this%data_z(im)%data
    else
      m%data = m%data + weights(im)*this%data_d(im)%data
    endif
  end do
end subroutine TBMatrix_sum_matrices_d

subroutine TBMatrix_transpose_sub(this, m)
  type(TBMatrix), intent(inout) :: this
  type(TBMatrix), intent(in) :: m

  integer im

  if (this%N /= m%N) call system_abort("TBMatrix_transpose_sub called with size mismatch")
  if (this%n_matrices /= m%n_matrices) call system_abort("TBMatrix_transpose_sub called with n_matrices mismatch")

  if (this%is_sparse .or. m%is_sparse) call system_abort("Can't do TBMatrix_transpose_sub on a sparse TBMatrix")

  do im=1, this%n_matrices
    if (this%is_complex) then
      if (m%is_complex) then
        call transpose_sub(this%data_z(im), m%data_z(im))
      else
	call system_abort("Can't TBMatrix_transpose_sub from real matrix into complex")
        ! call transpose_sub(this%data_z(im), m%data_d(im))
      endif
    else
      if (m%is_complex) then
        ! call transpose_sub(this%data_d(im), m%data_z(im))
	call system_abort("Can't TBMatrix_transpose_sub from complex matrix into real")
      else
        call transpose_sub(this%data_d(im), m%data_d(im))
      endif
    endif
  end do

end subroutine TBMatrix_transpose_sub

end module TBMatrix_module
