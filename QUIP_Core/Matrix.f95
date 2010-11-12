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
!X Matrix module 
!X
!% This module handles the math of real and complex matrices.
!% The objects 'MatrixD' and 'MatrixZ' correspond to real and complex matrices, 
!% respectively.
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module matrix_module

use System_module
use linearalgebra_module

use MPI_context_module
use ScaLAPACK_module

implicit none
private

public :: MatrixD
type MatrixD
  integer :: N = 0, M = 0, l_N = 0, l_M = 0
  real(dp), allocatable:: data(:,:)
  type(Matrix_ScaLAPACK_Info) ScaLAPACK_Info_obj
end type MatrixD

public :: MatrixZ
type MatrixZ
  integer :: N = 0, M = 0, l_N = 0, l_M = 0
  complex(dp), allocatable:: data(:,:)
  type(Matrix_ScaLAPACK_Info) ScaLAPACK_Info_obj
end type MatrixZ

public :: Initialise
interface Initialise
  module procedure MatrixD_Initialise, MatrixZ_Initialise
  module procedure MatrixD_Initialise_mat, MatrixZ_Initialise_mat
end interface Initialise

public :: Finalise
interface Finalise
  module procedure MatrixD_Finalise, MatrixZ_Finalise
end interface Finalise

!% Nullify all the matrix elements.
public :: Wipe
interface Wipe
  module procedure MatrixD_Wipe, MatrixZ_Wipe
end interface Wipe

!% Set at zero some elements of real and complex matrices.
!% Adding the optional logical values (\texttt{d_mask} and \texttt{od_mask}), the diagonal and 
!% off-diagonal elements are selected. 
public :: Zero
interface Zero
  module procedure MatrixD_Zero, MatrixZ_Zero
end interface Zero

public :: Print
interface Print
  module procedure MatrixD_Print, MatrixZ_Print
end interface Print

!% Add a block to a given matrix. The block dimensions have to be specified.
public :: add_block
interface add_block
  module procedure MatrixD_add_block, MatrixZ_add_block
end interface add_block

!% Diagonalise a real or complex matrix.
public :: diagonalise
interface diagonalise
  module procedure MatrixD_diagonalise, MatrixD_diagonalise_gen, MatrixZ_diagonalise, MatrixZ_diagonalise_gen
end interface diagonalise

! public :: accum_col_outer_product
! interface accum_col_outer_product
!   module procedure MatrixDD_accum_col_outer_product, MatrixZZ_accum_col_outer_product
! end interface accum_col_outer_product

!% Compute the trace of the matrix obtained multiplying two input matrices \texttt{a} and \texttt{b}.
public :: TraceMult
interface TraceMult
  module procedure Matrix_TraceMult_DD, Matrix_TraceMult_DZ, Matrix_TraceMult_ZZ, Matrix_TraceMult_ZD
end interface TraceMult

public :: partial_TraceMult
interface partial_TraceMult
  module procedure Matrix_partial_TraceMult_DD, Matrix_partial_TraceMult_DZ, Matrix_partial_TraceMult_ZZ, Matrix_partial_TraceMult_ZD
end interface partial_TraceMult

public :: partial_TraceMult_spinor
interface partial_TraceMult_spinor
  module procedure Matrix_partial_TraceMult_spinor_ZZ
end interface partial_TraceMult_spinor

!% Make real the diagonal elements of a given matrix  \texttt{a} .
public :: Re_diag
interface Re_diag
  module procedure Matrix_Re_diagD, Matrix_Re_diagZ
end interface Re_diag

!% Get diagonal 2x2 blocks a given matrix \texttt{a} .
public :: diag_spinor
interface diag_spinor
  module procedure Matrix_diag_spinorD, Matrix_diag_spinorZ
end interface diag_spinor

!% Compute the sum of two matrices (complex or real) \texttt{m$_1$} and \texttt{m$_2$} weighted by  
!% complex factors  \texttt{f1} and  \texttt{f2}. The result is a complex matrix.
public :: scaled_sum
interface scaled_sum
  module procedure Matrix_scaled_sum_ZZZ, Matrix_scaled_sum_ZDD
end interface scaled_sum

!% Add to a given matrix the elements of an other matrix  \texttt{m$_1$} weighted by a 
!% complex factor \texttt{f1}.
public :: scaled_accum
interface scaled_accum
  module procedure Matrix_scaled_accum_DZ, Matrix_scaled_accum_ZZ
end interface scaled_accum

!% Compute the inverse of a given matrix. 
public :: inverse
interface inverse
  module procedure MatrixD_inverse, MatrixZ_inverse
end interface inverse

public :: multDiag
interface multDiag
  module procedure MatrixD_multDiag_d, MatrixZ_multDiag_d, MatrixZ_multDiag_z
end interface multDiag

public :: multDiagRL
interface multDiagRL
  module procedure MatrixD_multDiagRL_d, MatrixZ_multDiagRL_d
end interface multDiagRL

public :: matrix_product_sub
interface matrix_product_sub
  module procedure MatrixZ_matrix_product_sub_zz, MatrixZ_matrix_product_sub_dz, MatrixZ_matrix_product_sub_zd
  module procedure MatrixD_matrix_product_sub_dr2, MatrixD_matrix_product_sub_dd
end interface matrix_product_sub

!% Add the identity matrix.
public :: add_identity
interface add_identity
  module procedure MatrixD_add_identity
end interface add_identity

!% Scale the matrix elements for a real factor  \texttt{scale}.
public :: scale
interface scale
  module procedure MatrixD_scale
end interface scale

public :: transpose_sub
interface transpose_sub
  module procedure MatrixD_transpose_sub, MatrixZ_transpose_sub
end interface transpose_sub

contains

subroutine MatrixD_Initialise(this, N, M, NB, MB, scalapack_obj)
  type(MatrixD), intent(out) :: this
  integer, intent(in), optional :: N, M, NB, MB
  type(ScaLAPACK), intent(in), optional :: scalapack_obj

  call Finalise(this)

  call matrixany_initialise(N, M, NB, MB, scalapack_obj, this%N, this%M, this%l_N, this%l_M, &
			    this%ScaLAPACK_Info_obj)
  if (this%l_N > 0 .and. this%l_M > 0) then
    allocate(this%data(this%l_N, this%l_M))
    call ALLOC_TRACE("MatrixD_Initialise "//this%l_N//" "//this%l_M, size(this%data)*REAL_SIZE)
  endif

end subroutine MatrixD_Initialise

subroutine MatrixD_Initialise_mat(this, from)
  type(MatrixD), intent(out) :: this
  type(MatrixD), intent(in) :: from

  call Finalise(this)

  this%N = from%N
  this%M = from%M
  this%l_N = from%l_N
  this%l_M = from%l_M
  this%ScaLAPACK_Info_obj = from%ScaLAPACK_Info_obj

  if (this%l_N > 0 .and. this%l_M > 0) then
    allocate(this%data(this%l_N, this%l_M))
    call ALLOC_TRACE("MatrixD_Initialise_mat "//this%l_N//" "//this%l_M,size(this%data)*REAL_SIZE)
  endif

end subroutine MatrixD_Initialise_mat

subroutine MatrixZ_Initialise(this, N, M, NB, MB, scalapack_obj)
  type(MatrixZ), intent(out) :: this
  integer, intent(in), optional :: N, M, NB, MB
  type(ScaLAPACK), intent(in), optional :: scalapack_obj

  call Finalise(this)

  call matrixany_initialise(N, M, NB, MB, scalapack_obj, this%N, this%M, this%l_N, this%l_M, &
			    this%ScaLAPACK_Info_obj)

  if (this%l_N > 0 .and. this%l_M > 0) then
    allocate(this%data(this%l_N, this%l_M))
    call ALLOC_TRACE("MatrixZ_Initialise "//this%l_N//" "//this%l_M, size(this%data)*COMPLEX_SIZE)
  endif

end subroutine MatrixZ_Initialise

subroutine MatrixZ_Initialise_mat(this, from)
  type(MatrixZ), intent(out) :: this
  type(MatrixZ), intent(in) :: from

  call Finalise(this)

  this%N = from%N
  this%M = from%M
  this%l_N = from%l_N
  this%l_M = from%l_M
  this%ScaLAPACK_Info_obj = from%ScaLAPACK_Info_obj

  if (this%l_N > 0 .and. this%l_M > 0) then
    allocate(this%data(this%l_N, this%l_M))
    call ALLOC_TRACE("MatrixZ_Initialise_mat "//this%l_N//" "//this%l_M, size(this%data)*COMPLEX_SIZE)
  endif

end subroutine MatrixZ_Initialise_mat

subroutine matrixany_initialise(N, M, NB, MB, scalapack_obj, this_N, this_M, this_l_N, this_l_M, &
			        this_ScaLAPACK_Info_obj)
  integer, intent(in), optional :: N, M, NB, MB
  type(ScaLAPACK), intent(in), target, optional :: scalapack_obj
  integer, intent(out) :: this_N, this_M, this_l_N, this_l_M
  type(Matrix_ScaLAPACK_Info), intent(out) :: this_ScaLAPACK_Info_obj

  integer use_NB, use_MB

  this_N = 0
  this_M = 0
  this_l_N = 0
  this_l_M = 0

  if (present(N)) then
    this_N = N
    if (present(M)) then
      this_M = M
    else
      this_M = N
    endif

    if (present(NB)) then
      use_NB = NB
    else
      use_NB = 36
    endif
    if (present(MB)) then
      use_MB = MB
    else
      use_MB = 36
    endif

    if (present(scalapack_obj)) then
      call Initialise(this_ScaLAPACK_Info_obj, this_N, this_M, use_NB, use_MB, scalapack_obj)
    else
      call Initialise(this_ScaLAPACK_Info_obj, N, M, use_NB, use_MB)
    endif
    if (this_ScaLAPACK_Info_obj%active) then
      this_l_N = this_ScaLAPACK_Info_obj%l_N_R
      this_l_M = this_ScaLAPACK_Info_obj%l_N_C
    else
      this_l_N = this_N
      this_l_M = this_M
    endif
  endif
end subroutine matrixany_initialise

subroutine MatrixD_Finalise(this)
  type(MatrixD), intent(inout) :: this

  call Wipe(this)
  call Finalise(this%ScaLAPACK_Info_obj)

end subroutine MatrixD_Finalise

subroutine MatrixD_Wipe(this)
  type(MatrixD), intent(inout) :: this

  call Wipe(this%ScaLAPACK_Info_obj)

  if (allocated(this%data)) then
    call DEALLOC_TRACE("MatrixD_Wipe ", size(this%data)*REAL_SIZE)
    deallocate(this%data)
  endif

  this%N = 0
  this%M = 0
  this%l_N = 0
  this%l_M = 0

end subroutine MatrixD_Wipe

subroutine MatrixD_Zero(this, d_mask, od_mask)
  type(MatrixD), intent(inout) :: this
  logical, intent(in), optional :: d_mask(:), od_mask(:)

  integer i

  if (allocated(this%data)) then
    if (present(d_mask) .or. present(od_mask)) then
      if (present(d_mask)) then
	do i=1, min(size(d_mask), this%N, this%M)
	  if (d_mask(i)) this%data(i,i) = 0.0_dp
	end do
      end if
      if (present(od_mask)) then
	do i=1, size(od_mask)
	  if (od_mask(i)) then
	    if (i <= this%M) this%data(:,i) = 0.0_dp
	    if (i <= this%N) this%data(i,:) = 0.0_dp
	  end if
	end do
      end if
    else
      this%data = 0.0_dp
    endif
  endif

end subroutine

subroutine MatrixZ_Finalise(this)
  type(MatrixZ), intent(inout) :: this

  call MatrixZ_Wipe(this)
  call Finalise(this%ScaLAPACK_Info_obj)

end subroutine MatrixZ_Finalise

subroutine MatrixZ_Wipe(this)
  type(MatrixZ), intent(inout) :: this

  call Wipe(this%ScaLAPACK_Info_obj)

  if (allocated(this%data)) then
    call DEALLOC_TRACE("MatrixZ_Wipe ", size(this%data)*COMPLEX_SIZE)
    deallocate(this%data)
  endif

  this%N = 0
  this%M = 0
  this%l_N = 0
  this%l_M = 0

end subroutine MatrixZ_Wipe

subroutine MatrixZ_Zero(this, d_mask, od_mask)
  type(MatrixZ), intent(inout) :: this
  logical, intent(in), optional :: d_mask(:), od_mask(:)

  integer i

  if (allocated(this%data)) then
    if (present(d_mask) .or. present(od_mask)) then
      if (present(d_mask)) then
	do i=1, min(size(d_mask), this%N, this%M)
	  if (d_mask(i)) this%data(i,i) = 0.0_dp
	end do
      end if
      if (present(od_mask)) then
	do i=1, size(od_mask)
	  if (od_mask(i)) then
	    if (i <= this%M) this%data(:,i) = 0.0_dp
	    if (i <= this%N) this%data(i,:) = 0.0_dp
	  end if
	end do
      end if
    else
      this%data = 0.0_dp
    endif
  endif

end subroutine

subroutine MatrixD_Print(this,file)
  type(MatrixD),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  if (current_verbosity() < PRINT_NORMAL) return

  call Print ("MatrixD : ", file=file)

  call Print ("N M " // this%N // " " // this%M // " l_N l_M "  // this%l_N // " " // this%l_M, file=file)
  call Print (this%ScaLAPACK_Info_obj, file=file)
  call Print ("MatrixD data:", file=file)
  if (this%ScaLAPACK_Info_obj%active)  then
    call Print(this%ScaLAPACK_Info_obj, this%data, file=file)
  else
    if (allocated(this%data)) then
      call Print(this%data, file=file)
    endif
  endif

end subroutine MatrixD_Print

subroutine MatrixZ_Print(this,file)
  type(MatrixZ),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  if (current_verbosity() < PRINT_NORMAL) return

  call Print ("MatrixZ : ", file=file)

  call Print ("N M " // this%N // " " // this%M // " l_N l_M "  // this%l_N // " " // this%l_M, file=file)
  call Print (this%ScaLAPACK_Info_obj, file=file)
  call Print ("MatrixZ data:", file=file)
  if (allocated(this%data)) then
    if (this%ScaLAPACK_Info_obj%active)  then
      call Print(this%ScaLAPACK_Info_obj, this%data, file=file)
    else
      call Print(this%data, file=file)
    endif
  endif

end subroutine MatrixZ_Print

subroutine MatrixD_add_block(this, block, block_nr, block_nc, first_row, first_col)
  type(MatrixD), intent(inout) :: this
  real(dp), intent(in) :: block(:,:)
  integer, intent(in) :: block_nr, block_nc   !% Number of rows and columns of the block
  integer, intent(in) :: first_row, first_col  

  integer i, j, io, jo, li, lj

  if (this%ScaLAPACK_Info_obj%active) then
    do jo=1, block_nc
    j = first_col + jo - 1
    do io=1, block_nr
      i = first_row + io - 1
      call coords_global_to_local(this%ScaLAPACK_Info_obj, i, j, li, lj)
      if (li > 0 .and. lj > 0) then
	this%data(li,lj) = this%data(li,lj) + block(io,jo)
      endif
    end do
    end do
  else
    this%data(first_row:first_row+block_nr-1,first_col:first_col+block_nc-1) = &
      this%data(first_row:first_row+block_nr-1,first_col:first_col+block_nc-1) + &
      block(1:block_nr,1:block_nc)
  endif

end subroutine

subroutine MatrixZ_add_block(this, block, block_nr, block_nc, first_row, first_col)
  type(MatrixZ), intent(inout) :: this
  complex(dp), intent(in) :: block(:,:)
  integer, intent(in) :: block_nr, block_nc
  integer, intent(in) :: first_row, first_col

  integer i, j, io, jo, li, lj

  if (this%ScaLAPACK_Info_obj%active) then
    do jo=1, block_nc
    j = first_col + jo - 1
    do io=1, block_nr
      i = first_row + io - 1
      call coords_global_to_local(this%ScaLAPACK_Info_obj, i, j, li, lj)
      if (li > 0 .and. lj > 0) then
	this%data(li,lj) = this%data(li,lj) + block(io,jo)
      endif
    end do
    end do
  else
    this%data(first_row:first_row+block_nr-1,first_col:first_col+block_nc-1) = &
      this%data(first_row:first_row+block_nr-1,first_col:first_col+block_nc-1) + &
      block(1:block_nr,1:block_nc)
  endif
end subroutine

subroutine MatrixD_diagonalise(this, evals, evecs, error)
  type(MatrixD), intent(in), target :: this
  real(dp), intent(inout) :: evals(:)    !% Eigenvalues
  type(MatrixD), intent(inout), target, optional :: evecs  !% Eigenvectors
  integer, intent(out), optional :: error

  real(dp), pointer :: u_evecs(:,:)
  type(Matrix_ScaLAPACK_Info), pointer :: evecs_ScaLAPACK_Info

  INIT_ERROR(error)

  if (present(evecs)) then
    u_evecs => evecs%data
    evecs_ScaLAPACK_Info => evecs%ScaLAPACK_Info_obj
  else
    allocate(u_evecs(this%l_N,this%l_M))
    call ALLOC_TRACE("MatrixD_diagonalise evecs ", size(u_evecs)*REAL_SIZE)
    evecs_ScaLAPACK_Info => this%ScaLAPACK_Info_obj
  endif

  if (this%ScaLAPACK_Info_obj%active) then
    call diagonalise(this%ScaLAPACK_Info_obj, this%data, evals, evecs_ScaLAPACK_Info, u_evecs, error = error)
  else
    call diagonalise(this%data, evals, u_evecs, error = error)
  endif

  if (.not.present(evecs)) then
    call DEALLOC_TRACE("MatrixD_diagonalise evecs ",-size(u_evecs)*REAL_SIZE)
    deallocate(u_evecs)
  endif

  PASS_ERROR(error)

end subroutine MatrixD_diagonalise

subroutine MatrixD_diagonalise_gen(this, overlap, evals, evecs, error)
  type(MatrixD), intent(in), target :: this
  type(MatrixD), intent(in) :: overlap
  real(dp), intent(inout) :: evals(:)
  type(MatrixD), intent(inout), target, optional :: evecs
  integer, intent(out), optional :: error

  real(dp), pointer :: u_evecs(:,:)
  type(Matrix_ScaLAPACK_Info), pointer :: evecs_ScaLAPACK_Info

  INIT_ERROR(error)

  if (present(evecs)) then
    u_evecs => evecs%data
    evecs_ScaLAPACK_Info => evecs%ScaLAPACK_Info_obj
  else
    allocate(u_evecs(this%l_N,this%l_M))
    call ALLOC_TRACE("MatrixD_diagonalise_gen evecs", size(u_evecs)*REAL_SIZE)
    evecs_ScaLAPACK_Info => this%ScaLAPACK_Info_obj
  endif

  if (this%ScaLAPACK_Info_obj%active) then
    call diagonalise(this%ScaLAPACK_Info_obj, this%data, overlap%ScaLAPACK_Info_obj, &
      overlap%data, evals, evecs_ScaLAPACK_Info, u_evecs, error)
  else
    call diagonalise(this%data, overlap%data, evals, u_evecs, error = error)
  endif

  if (.not.present(evecs)) then
    call DEALLOC_TRACE("MatrixD_diagonalise_gen u_evecs", size(u_evecs)*REAL_SIZE)
    deallocate(u_evecs)
  endif

  PASS_ERROR(error)

end subroutine MatrixD_diagonalise_gen

subroutine MatrixZ_diagonalise(this, evals, evecs, error)
  type(MatrixZ), intent(in), target :: this
  real(dp), intent(inout) :: evals(:)
  type(MatrixZ), intent(inout), target, optional :: evecs
  integer, intent(out), optional :: error

  complex(dp), pointer :: u_evecs(:,:)
  type(Matrix_ScaLAPACK_Info), pointer :: evecs_ScaLAPACK_Info

  INIT_ERROR(error)

  if (present(evecs)) then
    u_evecs => evecs%data
    evecs_ScaLAPACK_Info => evecs%ScaLAPACK_Info_obj
  else
    allocate(u_evecs(this%l_N,this%l_M))
    call ALLOC_TRACE("MatrixZ_diagonalise u_evecs", size(u_evecs)*COMPLEX_SIZE)
    evecs_ScaLAPACK_Info => this%ScaLAPACK_Info_obj
  endif

  if (this%ScaLAPACK_Info_obj%active) then
    call diagonalise(this%ScaLAPACK_Info_obj, this%data, evals, evecs_ScaLAPACK_Info, u_evecs, error)
  else
    call diagonalise(this%data, evals, u_evecs, error = error)
  endif

  if (.not.present(evecs)) then
    call DEALLOC_TRACE("MatrixZ_diagonalise u_evecs", size(u_evecs)*COMPLEX_SIZE)
    deallocate(u_evecs)
  endif

  PASS_ERROR(error)

end subroutine MatrixZ_diagonalise

subroutine MatrixZ_diagonalise_gen(this, overlap, evals, evecs, error)
  type(MatrixZ), intent(in), target :: this
  type(MatrixZ), intent(in) :: overlap
  real(dp), intent(inout) :: evals(:)
  type(MatrixZ), intent(inout), target, optional :: evecs
  integer, intent(out), optional :: error

  complex(dp), pointer :: u_evecs(:,:)
  type(Matrix_ScaLAPACK_Info), pointer :: evecs_ScaLAPACK_Info

  INIT_ERROR(error)

  if (present(evecs)) then
    u_evecs => evecs%data
    evecs_ScaLAPACK_Info => evecs%ScaLAPACK_Info_obj
  else
    allocate(u_evecs(this%l_N,this%l_M))
    call ALLOC_TRACE("MatrixZ_diagonalise_gen u_evecs", size(u_evecs)*COMPLEX_SIZE)
    evecs_ScaLAPACK_Info => this%ScaLAPACK_Info_obj
  endif

  if (this%ScaLAPACK_Info_obj%active) then
    call diagonalise(this%ScaLAPACK_Info_obj, this%data, overlap%ScaLAPACK_Info_obj, &
      overlap%data, evals, evecs_ScaLAPACK_Info, u_evecs, error)
  else
    call diagonalise(this%data, overlap%data, evals, u_evecs, error = error)
  endif

  if (.not.present(evecs)) then
    call DEALLOC_TRACE("MatrixZ_diagonalise_gen u_evecs", size(u_evecs)*COMPLEX_SIZE)
    deallocate(u_evecs)
  endif

  PASS_ERROR(error)

end subroutine MatrixZ_diagonalise_gen

!subroutine MatrixDD_accum_col_outer_product(this, mati, i, f)
!  type(MatrixD), intent(inout) :: this
!  type(MatrixD), intent(in) :: mati
!  integer, intent(in) :: i
!  real(dp), intent(in) :: f
!
!  if (this%ScaLAPACK_Info_obj%active) then
!    call system_abort ("No support for ScaLAPACK yet in MatrixD_accum_col_outer_product")
!  else
!    this%data = this%data + f*(mati%data(:,i) .outer. mati%data(:,i))
!  endif
!end subroutine MatrixDD_accum_col_outer_product
!
!subroutine MatrixZZ_accum_col_outer_product(this, mati, i, f)
!  type(MatrixZ), intent(inout) :: this
!  type(MatrixZ), intent(in) :: mati
!  integer, intent(in) :: i
!  real(dp), intent(in) :: f
!
!  if (this%ScaLAPACK_Info_obj%active) then
!    call system_abort ("No support for ScaLAPACK yet in MatrixZ_accum_col_outer_product")
!  else
!    this%data = this%data + f*(mati%data(:,i) .outer. mati%data(:,i))
!  endif
!end subroutine MatrixZZ_accum_col_outer_product

function Matrix_partial_TraceMult_DD(a, b, a_T, b_T)
  type(MatrixD), intent(in) :: a, b
  logical, intent(in), optional :: a_T, b_T
  real(dp) :: Matrix_partial_TraceMult_DD(a%N)

  logical u_a_T, u_b_T
  integer i, j, g_i, g_j

  if ((a%ScaLAPACK_Info_obj%active .and. .not. b%ScaLAPACK_Info_obj%active) .or. &
      (.not. a%ScaLAPACK_Info_obj%active .and. b%ScaLAPACK_Info_obj%active)) then
    call system_abort("Can't do partial tracemult for mixed dense regular and ScaLAPACK matrices")
  endif

  u_a_T = optional_default(.false., a_T)
  u_b_T = optional_default(.false., b_T)

  if (a%ScaLAPACK_Info_obj%active .or. b%ScaLAPACK_Info_obj%active) then
    if (u_a_T .and. .not. u_b_T) then
      Matrix_partial_TraceMult_DD = 0.0_dp
      do j=1, a%l_M
	call coords_local_to_global(a%ScaLAPACK_Info_obj, 1, j, g_i, g_j)
	Matrix_partial_TraceMult_DD(g_j) = sum(a%data(:,j)*b%data(:,j))
      end do
    else if (.not. u_a_T .and. u_b_T) then
      Matrix_partial_TraceMult_DD = 0.0_dp
      do i=1, a%l_N
	call coords_local_to_global(a%ScaLAPACK_Info_obj, i, 1, g_i, g_j)
	Matrix_partial_TraceMult_DD(g_i) = sum(a%data(i,:)*b%data(i,:))
      end do
    else
      call system_abort("Can only do partial_TraceMult for ScaLAPACK matrices if one is transposed and the other not")
    endif
  else
    if (.not. u_a_T .and. .not. u_b_T) then
      do i=1, a%N
	Matrix_partial_TraceMult_DD(i) = sum(a%data(i,:)*b%data(:,i))
      end do
    else if (u_a_T .and. u_b_T) then
      do i=1, a%N
	Matrix_partial_TraceMult_DD(i) = sum(a%data(:,i)*b%data(i,:))
      end do
    else if (.not. u_a_T .and. u_b_T) then
      Matrix_partial_TraceMult_DD(:) = 0.0_dp
      do i=1, a%N
	Matrix_partial_TraceMult_DD(:) = Matrix_partial_TraceMult_DD(:) + a%data(:,i)*b%data(:,i)
      end do
    else
      do i=1, a%N
	Matrix_partial_TraceMult_DD(i) = sum(a%data(:,i)*b%data(:,i))
      end do
    endif
  endif

end function Matrix_partial_TraceMult_DD

function Matrix_partial_TraceMult_DZ(a, b, a_T, b_H)
  type(MatrixD), intent(in) :: a
  type(MatrixZ), intent(in) :: b
  logical, intent(in), optional :: a_T, b_H
  complex(dp) :: Matrix_partial_TraceMult_DZ(a%N)

  logical u_a_T, u_b_H
  integer i, j, g_i, g_j

  if ((a%ScaLAPACK_Info_obj%active .and. .not. b%ScaLAPACK_Info_obj%active) .or. &
      (.not. a%ScaLAPACK_Info_obj%active .and. b%ScaLAPACK_Info_obj%active)) then
    call system_abort("Can't do partial tracemult for mixed dense regular and ScaLAPACK matrices")
  endif

  u_a_T = optional_default(.false., a_T)
  u_b_H = optional_default(.false., b_H)

  if (a%ScaLAPACK_Info_obj%active .or. b%ScaLAPACK_Info_obj%active) then
    if (u_a_T .and. .not. u_b_H) then
      Matrix_partial_TraceMult_DZ = 0.0_dp
      do j=1, a%l_M
	call coords_local_to_global(a%ScaLAPACK_Info_obj, 1, j, g_i, g_j)
	Matrix_partial_TraceMult_DZ(g_j) = sum(a%data(:,j)*b%data(:,j))
      end do
    else if (.not. u_a_T .and. u_b_H) then
      Matrix_partial_TraceMult_DZ = 0.0_dp
      do i=1, a%l_N
	call coords_local_to_global(a%ScaLAPACK_Info_obj, i, 1, g_i, g_j)
	Matrix_partial_TraceMult_DZ(g_i) = sum(a%data(i,:)*conjg(b%data(i,:)))
      end do
    else
      call system_abort("Can only do partial_TraceMult for ScaLAPACK matrices if one is transposed and the other not")
    endif
  else ! no scalapack
    if (.not. u_a_T .and. .not. u_b_H) then
      do i=1, a%N
	Matrix_partial_TraceMult_DZ(i) = sum(a%data(i,:)*b%data(:,i))
      end do
    else if (u_a_T .and. u_b_H) then
      do i=1, a%N
	Matrix_partial_TraceMult_DZ(i) = sum(a%data(:,i)*conjg(b%data(i,:)))
      end do
    else if (.not. u_a_T .and. u_b_H) then
      Matrix_partial_TraceMult_DZ(:) = 0.0_dp
      do i=1, a%N
	Matrix_partial_TraceMult_DZ(:) = Matrix_partial_TraceMult_DZ(:) + a%data(:,i)*conjg(b%data(:,i))
      end do
    else
      do i=1, a%N
	Matrix_partial_TraceMult_DZ(i) = sum(a%data(:,i)*b%data(:,i))
      end do
    endif
  endif

end function Matrix_partial_TraceMult_DZ

function Matrix_partial_TraceMult_ZD(a, b, a_H, b_T)
  type(MatrixZ), intent(in) :: a
  type(MatrixD), intent(in) :: b
  logical, intent(in), optional :: a_H, b_T
  complex(dp) :: Matrix_partial_TraceMult_ZD(a%N)

  logical u_a_H, u_b_T
  integer i, j, g_i, g_j

  if ((a%ScaLAPACK_Info_obj%active .and. .not. b%ScaLAPACK_Info_obj%active) .or. &
      (.not. a%ScaLAPACK_Info_obj%active .and. b%ScaLAPACK_Info_obj%active)) then
    call system_abort("Can't do partial tracemult for mixed dense regular and ScaLAPACK matrices")
  endif

  u_a_H = optional_default(.false., a_H)
  u_b_T = optional_default(.false., b_T)

  if (a%ScaLAPACK_Info_obj%active .or. b%ScaLAPACK_Info_obj%active) then
    if (u_a_H .and. .not. u_b_T) then
      Matrix_partial_TraceMult_ZD = 0.0_dp
      do j=1, a%l_M
	call coords_local_to_global(a%ScaLAPACK_Info_obj, 1, j, g_i, g_j)
	Matrix_partial_TraceMult_ZD(g_j) = sum(conjg(a%data(:,j))*b%data(:,j))
      end do
    else if (.not. u_a_H .and. u_b_T) then
      Matrix_partial_TraceMult_ZD = 0.0_dp
      do i=1, a%l_N
	call coords_local_to_global(a%ScaLAPACK_Info_obj, i, 1, g_i, g_j)
	Matrix_partial_TraceMult_ZD(g_i) = sum(a%data(i,:)*b%data(i,:))
      end do
    else
      call system_abort("Can only do partial_TraceMult for ScaLAPACK matrices if one is transposed and the other not")
    endif
  else
    if (.not. u_a_H .and. .not. u_b_T) then
      do i=1, a%N
	Matrix_partial_TraceMult_ZD(i) = sum(a%data(i,:)*b%data(:,i))
      end do
    else if (u_a_H .and. u_b_T) then
      do i=1, a%N
	Matrix_partial_TraceMult_ZD(i) = sum(conjg(a%data(:,i))*b%data(i,:))
      end do
    else if (.not. u_a_H .and. u_b_T) then
      Matrix_partial_TraceMult_ZD(:) = 0.0_dp
      do i=1, a%N
	Matrix_partial_TraceMult_ZD(:) = Matrix_partial_TraceMult_ZD(:) + a%data(:,i)*b%data(:,i)
      end do
    else
      do i=1, a%N
	Matrix_partial_TraceMult_ZD(i) = sum(conjg(a%data(:,i))*b%data(:,i))
      end do
    endif
  endif

end function Matrix_partial_TraceMult_ZD

function Matrix_partial_TraceMult_ZZ(a, b, a_H, b_H)
  type(MatrixZ), intent(in) :: a
  type(MatrixZ), intent(in) :: b
  logical, intent(in), optional :: a_H, b_H
  complex(dp) :: Matrix_partial_TraceMult_ZZ(a%N)

  logical u_a_H, u_b_H
  integer i, j, g_i, g_j

  if ((a%ScaLAPACK_Info_obj%active .and. .not. b%ScaLAPACK_Info_obj%active) .or. &
      (.not. a%ScaLAPACK_Info_obj%active .and. b%ScaLAPACK_Info_obj%active)) then
    call system_abort("Can't do partial tracemult for mixed dense regular and ScaLAPACK matrices")
  endif

  u_a_H = optional_default(.false., a_H)
  u_b_H = optional_default(.false., b_H)

  if (a%ScaLAPACK_Info_obj%active .or. b%ScaLAPACK_Info_obj%active) then
    if (u_a_H .and. .not. u_b_H) then
      Matrix_partial_TraceMult_ZZ = 0.0_dp
      do j=1, a%l_M
	call coords_local_to_global(a%ScaLAPACK_Info_obj, 1, j, g_i, g_j)
	Matrix_partial_TraceMult_ZZ(g_j) = sum(conjg(a%data(:,j))*b%data(:,j))
      end do
    else if (.not. u_a_H .and. u_b_H) then
      Matrix_partial_TraceMult_ZZ = 0.0_dp
      do i=1, a%l_N
	call coords_local_to_global(a%ScaLAPACK_Info_obj, i, 1, g_i, g_j)
	Matrix_partial_TraceMult_ZZ(g_i) = sum(a%data(i,:)*conjg(b%data(i,:)))
      end do
    else
      call system_abort("Can only do partial_TraceMult for ScaLAPACK matrices if one is transposed and the other not")
    endif
  else
    if (.not. u_a_H .and. .not. u_b_H) then
      do i=1, a%N
	Matrix_partial_TraceMult_ZZ(i) = sum(a%data(i,:)*b%data(:,i))
      end do
    else if (u_a_H .and. u_b_H) then
      do i=1, a%N
	Matrix_partial_TraceMult_ZZ(i) = sum(conjg(a%data(:,i))*conjg(b%data(i,:)))
      end do
    else if (.not. u_a_H .and. u_b_H) then
      Matrix_partial_TraceMult_ZZ(:) = 0.0_dp
      do i=1, a%N
	Matrix_partial_TraceMult_ZZ(:) = Matrix_partial_TraceMult_ZZ(:) + a%data(:,i)*conjg(b%data(:,i))
      end do
    else
      do i=1, a%N
	Matrix_partial_TraceMult_ZZ(i) = sum(conjg(a%data(:,i))*b%data(:,i))
      end do
    endif
  endif

end function Matrix_partial_TraceMult_ZZ

function Matrix_partial_TraceMult_spinor_ZZ(a, b, a_H, b_H)
  type(MatrixZ), intent(in) :: a
  type(MatrixZ), intent(in) :: b
  logical, intent(in), optional :: a_H, b_H
  complex(dp) :: Matrix_partial_TraceMult_spinor_ZZ(2,2,a%N/2)

  logical u_a_H, u_b_H
  integer i, j, g_i, g_j

  if ((a%ScaLAPACK_Info_obj%active .and. .not. b%ScaLAPACK_Info_obj%active) .or. &
      (.not. a%ScaLAPACK_Info_obj%active .and. b%ScaLAPACK_Info_obj%active)) then
    call system_abort("Can't do partial tracemult for mixed dense regular and ScaLAPACK matrices")
  endif

  u_a_H = optional_default(.false., a_H)
  u_b_H = optional_default(.false., b_H)

  if (a%ScaLAPACK_Info_obj%active .or. b%ScaLAPACK_Info_obj%active) then
    if (u_a_H .and. .not. u_b_H) then
      Matrix_partial_TraceMult_spinor_ZZ = 0.0_dp
      do j=1, a%l_M
	call coords_local_to_global(a%ScaLAPACK_Info_obj, 1, j, g_i, g_j)
	Matrix_partial_TraceMult_spinor_ZZ(1,1,(g_j-1)/2+1) = sum(conjg(a%data(1:a%l_N:2,j))*b%data(1:a%l_N:2,j))
	Matrix_partial_TraceMult_spinor_ZZ(1,2,(g_j-1)/2+1) = sum(conjg(a%data(2:a%l_N:2,j))*b%data(2:a%l_N:2,j))
	Matrix_partial_TraceMult_spinor_ZZ(2,1,(g_j-1)/2+1) = sum(conjg(a%data(1:a%l_N:2,j+1))*b%data(1:a%l_N:2,j+1))
	Matrix_partial_TraceMult_spinor_ZZ(2,2,(g_j-1)/2+1) = sum(conjg(a%data(2:a%l_N:2,j+1))*b%data(2:a%l_N:2,j+1))
      end do
    else if (.not. u_a_H .and. u_b_H) then
      Matrix_partial_TraceMult_spinor_ZZ = 0.0_dp
      do i=1, a%l_N, 2
	call coords_local_to_global(a%ScaLAPACK_Info_obj, i, 1, g_i, g_j)
	Matrix_partial_TraceMult_spinor_ZZ(1,1,(g_i-1)/2+1) = sum(a%data(i,1:a%l_M:2)*conjg(b%data(i,1:b%l_M:2)))
	Matrix_partial_TraceMult_spinor_ZZ(1,2,(g_i-1)/2+1) = sum(a%data(i,2:a%l_M:2)*conjg(b%data(i,2:b%l_M:2)))
	Matrix_partial_TraceMult_spinor_ZZ(2,1,(g_i-1)/2+1) = sum(a%data(i+1,1:a%l_M:2)*conjg(b%data(i+1,1:b%l_M:2)))
	Matrix_partial_TraceMult_spinor_ZZ(2,2,(g_i-1)/2+1) = sum(a%data(i+1,2:a%l_M:2)*conjg(b%data(i+1,2:b%l_M:2)))
      end do
    else
      call system_abort("Can only do partial_TraceMult_spinor for ScaLAPACK matrices if one is transposed and the other not")
    endif
  else
    if (.not. u_a_H .and. .not. u_b_H) then
      do i=1, a%N, 2
	Matrix_partial_TraceMult_spinor_ZZ(1,1,(i-1)/2+1) = sum(a%data(i,1:a%l_M:2)*b%data(1:b%l_N:2,i))
	Matrix_partial_TraceMult_spinor_ZZ(1,2,(i-1)/2+1) = sum(a%data(i,2:a%l_M:2)*b%data(2:b%l_N:2,i))
	Matrix_partial_TraceMult_spinor_ZZ(2,1,(i-1)/2+1) = sum(a%data(i+1,1:a%l_M:2)*b%data(1:b%l_N:2,i+1))
	Matrix_partial_TraceMult_spinor_ZZ(2,2,(i-1)/2+1) = sum(a%data(i+1,2:a%l_M:2)*b%data(2:b%l_N:2,i+1))
      end do
    else if (u_a_H .and. u_b_H) then
      do i=1, a%N, 2
	Matrix_partial_TraceMult_spinor_ZZ(1,1,(i-1)/2+1) = sum(conjg(a%data(1:a%l_N:2,i))*conjg(b%data(i,1:b%l_M:2)))
	Matrix_partial_TraceMult_spinor_ZZ(1,2,(i-1)/2+1) = sum(conjg(a%data(2:a%l_N:2,i))*conjg(b%data(i,2:b%l_M:2)))
	Matrix_partial_TraceMult_spinor_ZZ(2,1,(i-1)/2+1) = sum(conjg(a%data(1:a%l_N:2,i+1))*conjg(b%data(i+1,1:b%l_M:2)))
	Matrix_partial_TraceMult_spinor_ZZ(2,2,(i-1)/2+1) = sum(conjg(a%data(2:a%l_N:2,i+1))*conjg(b%data(i+1,2:b%l_M:2)))
      end do
    else if (.not. u_a_H .and. u_b_H) then
      Matrix_partial_TraceMult_spinor_ZZ(:,:,:) = 0.0_dp
      do i=1, a%N, 2
	Matrix_partial_TraceMult_spinor_ZZ(1,1,:) = Matrix_partial_TraceMult_spinor_ZZ(1,1,:) + a%data(1:a%l_N:2,i)*conjg(b%data(1:b%l_N:2,i))
	Matrix_partial_TraceMult_spinor_ZZ(1,2,:) = Matrix_partial_TraceMult_spinor_ZZ(1,2,:) + a%data(1:a%l_N:2,i+1)*conjg(b%data(1:b%l_N:2,i+1))
	Matrix_partial_TraceMult_spinor_ZZ(2,1,:) = Matrix_partial_TraceMult_spinor_ZZ(2,1,:) + a%data(2:a%l_N:2,i)*conjg(b%data(2:b%l_N:2,i))
	Matrix_partial_TraceMult_spinor_ZZ(2,2,:) = Matrix_partial_TraceMult_spinor_ZZ(2,2,:) + a%data(2:a%l_N:2,i+1)*conjg(b%data(2:b%l_N:2,i+1))
      end do
    else
      do i=1, a%N, 2
	Matrix_partial_TraceMult_spinor_ZZ(1,1,(i-1)/2+1) = sum(conjg(a%data(1:a%l_N:2,i))*b%data(1:b%l_N:2,i))
	Matrix_partial_TraceMult_spinor_ZZ(1,2,(i-1)/2+1) = sum(conjg(a%data(2:a%l_N:2,i))*b%data(2:b%l_N:2,i))
	Matrix_partial_TraceMult_spinor_ZZ(2,1,(i-1)/2+1) = sum(conjg(a%data(1:a%l_N:2,i+1))*b%data(1:b%l_N:2,i+1))
	Matrix_partial_TraceMult_spinor_ZZ(2,2,(i-1)/2+1) = sum(conjg(a%data(2:a%l_N:2,i+1))*b%data(2:b%l_N:2,i+1))
      end do
    endif
  endif

end function Matrix_partial_TraceMult_spinor_ZZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Matrix_TraceMult_DD(a, b, w, a_T, b_T)
  type(MatrixD), intent(in) :: a, b
  real(dp), intent(in), optional :: w(:)
  logical, intent(in), optional :: a_T, b_T
  real(dp) :: Matrix_TraceMult_DD

  if (present(w)) then
    Matrix_TraceMult_DD = sum(partial_TraceMult(a, b, a_T, b_T)*w)
  else
    Matrix_TraceMult_DD = sum(partial_TraceMult(a, b, a_T, b_T))
  endif

end function Matrix_TraceMult_DD

function Matrix_TraceMult_DZ(a, b, w, a_T, b_H)
  type(MatrixD), intent(in) :: a
  type(MatrixZ), intent(in) :: b
  real(dp), intent(in), optional :: w(:)
  logical, intent(in), optional :: a_T, b_H
  complex(dp) :: Matrix_TraceMult_DZ

  if (present(w)) then
    Matrix_TraceMult_DZ = sum(partial_TraceMult(a, b, a_T, b_H)*w)
  else
    Matrix_TraceMult_DZ = sum(partial_TraceMult(a, b, a_T, b_H))
  endif

end function Matrix_TraceMult_DZ

function Matrix_TraceMult_ZD(a, b, w, a_H, b_T)
  type(MatrixZ), intent(in) :: a
  type(MatrixD), intent(in) :: b
  real(dp), intent(in), optional :: w(:)
  logical, intent(in), optional :: a_H, b_T
  complex(dp) :: Matrix_TraceMult_ZD

  if (present(w)) then
    Matrix_TraceMult_ZD = sum(partial_TraceMult(a, b, a_H, b_T)*w)
  else
    Matrix_TraceMult_ZD = sum(partial_TraceMult(a, b, a_H, b_T))
  endif

end function Matrix_TraceMult_ZD

function Matrix_TraceMult_ZZ(a, b, w, a_H, b_H)
  type(MatrixZ), intent(in) :: a
  type(MatrixZ), intent(in) :: b
  real(dp), intent(in), optional :: w(:)
  logical, intent(in), optional :: a_H, b_H
  complex(dp) :: Matrix_TraceMult_ZZ

  if (present(w)) then
    Matrix_TraceMult_ZZ = sum(partial_TraceMult(a, b, a_H, b_H)*w)
  else
    Matrix_TraceMult_ZZ = sum(partial_TraceMult(a, b, a_H, b_H))
  endif

end function Matrix_TraceMult_ZZ

function Matrix_Re_diagZ(a)
  type(MatrixZ), intent(in) :: a
  real(dp) :: Matrix_Re_diagZ(a%N)

  integer i

  if (a%ScaLAPACK_Info_obj%active) then
    Matrix_Re_diagZ = Re_diag(a%ScaLAPACK_Info_obj, a%data)
  else
    do i=1, a%N
      Matrix_Re_diagZ(i) = real(a%data(i,i))
    end do
  endif
end function Matrix_Re_diagZ

function Matrix_Re_diagD(a)
  type(MatrixD), intent(in) :: a
  real(dp) :: Matrix_Re_diagD(a%N)

  integer i

  if (a%ScaLAPACK_Info_obj%active) then
    Matrix_Re_diagD = Re_diag(a%ScaLAPACK_Info_obj, a%data)
  else
    do i=1, a%N
      Matrix_Re_diagD(i) = a%data(i,i)
    end do
  endif
end function Matrix_Re_diagD

function Matrix_diag_spinorZ(a)
  type(MatrixZ), intent(in) :: a
  complex(dp) :: Matrix_diag_spinorZ(2,2,a%N/2)

  integer i

  if (a%ScaLAPACK_Info_obj%active) then
    Matrix_diag_spinorZ = diag_spinor(a%ScaLAPACK_Info_obj, a%data)
  else
    Matrix_diag_spinorZ = 0.0_dp
    do i=1, a%N,2
      Matrix_diag_spinorZ(:,:,(i-1)/2+1) = a%data(i:i+1,i:i+1)
    end do
  endif
end function Matrix_diag_spinorZ

function Matrix_diag_spinorD(a)
  type(MatrixD), intent(in) :: a
  complex(dp) :: Matrix_diag_spinorD(2,2,a%N/2)

  integer i

  if (a%ScaLAPACK_Info_obj%active) then
    Matrix_diag_spinorD = diag_spinor(a%ScaLAPACK_Info_obj, a%data)
  else
    Matrix_diag_spinorD = 0.0_dp
    do i=1, a%N, 2
      Matrix_diag_spinorD(:,:,(i-1)/2+1) = a%data(i:i+1,i:i+1)
    end do
  endif
end function Matrix_diag_spinorD

subroutine Matrix_scaled_sum_ZDD(this, f1, m1, f2, m2)
  type(MatrixZ), intent(inout) :: this
  complex(dp), intent(in) :: f1
  type(MatrixD), intent(in) :: m1
  complex(dp), intent(in) :: f2
  type(MatrixD), intent(in) :: m2

  integer i

  if ( (m1%ScaLAPACK_Info_obj%active .and. .not. m2%ScaLAPACK_Info_obj%active) .or. &
       (.not. m1%ScaLAPACK_Info_obj%active .and. m2%ScaLAPACK_Info_obj%active) ) then
    call system_abort("Can't do scaled_sum for mixed ScaLAPCAK non-ScaLAPACK matrices")
  endif

  do i=1, size(this%data,2)
    this%data(:,i) = f1*m1%data(:,i) + f2*m2%data(:,i)
  end do

end subroutine Matrix_scaled_sum_ZDD

subroutine Matrix_scaled_sum_ZZZ(this, f1, m1, f2, m2)
  type(MatrixZ), intent(inout) :: this
  complex(dp), intent(in) :: f1
  type(MatrixZ), intent(in) :: m1
  complex(dp), intent(in) :: f2
  type(MatrixZ), intent(in) :: m2

  integer i

  if ( (m1%ScaLAPACK_Info_obj%active .and. .not. m2%ScaLAPACK_Info_obj%active) .or. &
       (.not. m1%ScaLAPACK_Info_obj%active .and. m2%ScaLAPACK_Info_obj%active) ) then
    call system_abort("Can't do scaled_sum for mixed ScaLAPCAK non-ScaLAPACK matrices")
  endif

  do i=1, size(this%data,2)
    this%data(:,i) = f1*m1%data(:,i) + f2*m2%data(:,i)
  end do

end subroutine Matrix_scaled_sum_ZZZ

subroutine Matrix_scaled_accum_DZ(this, f1, m1)
  type(MatrixD), intent(inout) :: this
  complex(dp), intent(in) :: f1
  type(MatrixZ), intent(in) :: m1

  integer i

  if ( (this%ScaLAPACK_Info_obj%active .and. .not. m1%ScaLAPACK_Info_obj%active) .or. &
       (.not. this%ScaLAPACK_Info_obj%active .and. m1%ScaLAPACK_Info_obj%active) ) then
    call system_abort("Can't do scaled_accum for mixed ScaLAPCAK non-ScaLAPACK matrices")
  endif

  do i=1, size(this%data,2)
    this%data(:,i) = this%data(:,i) + f1*m1%data(:,i)
  end do

end subroutine Matrix_scaled_accum_DZ

subroutine Matrix_scaled_accum_ZZ(this, f1, m1)
  type(MatrixZ), intent(inout) :: this
  complex(dp), intent(in) :: f1
  type(MatrixZ), intent(in) :: m1

  integer i

  if ( (this%ScaLAPACK_Info_obj%active .and. .not. m1%ScaLAPACK_Info_obj%active) .or. &
       (.not. this%ScaLAPACK_Info_obj%active .and. m1%ScaLAPACK_Info_obj%active) ) then
    call system_abort("Can't do scaled_accum for mixed ScaLAPCAK non-ScaLAPACK matrices")
  endif

  do i=1, size(this%data,2)
    this%data(:,i) = this%data(:,i) + f1*m1%data(:,i)
  end do

end subroutine Matrix_scaled_accum_ZZ


subroutine MatrixD_inverse(this, inv, positive)
  type(MatrixD), intent(inout) :: this
  type(MatrixD), intent(out), optional :: inv
  logical, intent(in), optional :: positive

  if (this%ScaLAPACK_Info_obj%active) then
    if (present(inv)) then
      call inverse(this%ScaLAPACK_Info_obj, this%data, inv%ScaLAPACK_Info_obj, inv%data, positive)
    else
      call inverse(this%ScaLAPACK_Info_obj, this%data, positive_in = positive)
    endif
  else
    if (present(inv)) then
      call inverse(this%data, inv%data, positive)
    else
      call inverse(this%data, positive_in = positive)
    endif
  endif
end subroutine MatrixD_inverse

subroutine MatrixZ_inverse(this, inv, positive)
  type(MatrixZ), intent(inout) :: this
  type(MatrixZ), intent(out), optional :: inv
  logical, intent(in), optional :: positive

  if (this%ScaLAPACK_Info_obj%active) then
    if (present(inv)) then
      call inverse(this%ScaLAPACK_Info_obj, this%data, inv%ScaLAPACK_Info_obj, inv%data, positive)
    else
      call inverse(this%ScaLAPACK_Info_obj, this%data, positive_in = positive)
    endif
  else
    if (present(inv)) then
      call inverse(this%data, inv%data, positive)
    else
      call inverse(this%data, positive_in = positive)
    endif
  endif

end subroutine MatrixZ_inverse

subroutine MatrixD_multDiag_d(this, A, diag)
  type(MatrixD), intent(inout) :: this
  type(MatrixD), intent(in) :: A
  real(dp), intent(in) :: diag(:)

  if (this%M /= size(diag) .or. A%M /= size(diag)) &
    call system_abort("Called MatrixD_multDiag_d with mismatched sizes")

  if (this%ScaLAPACK_Info_obj%active .and. A%ScaLAPACK_Info_obj%active) then
    call matrix_product_vect_asdiagonal_sub(this%ScaLAPACK_Info_obj, this%data, A%ScaLAPACK_Info_obj, A%data, diag)
  else if (.not.this%ScaLAPACK_Info_obj%active .and. .not.A%ScaLAPACK_Info_obj%active) then
    call matrix_product_vect_asdiagonal_sub(this%data, A%data, diag)
  else
    call system_abort("Called MatrixD_multDiag_d with mix of ScaLAPACK and non-ScaLAPACK objects")
  endif
end subroutine MatrixD_multDiag_d

subroutine MatrixZ_multDiag_d(this, A, diag)
  type(MatrixZ), intent(inout) :: this
  type(MatrixZ), intent(in) :: A
  real(dp), intent(in) :: diag(:)

  if (this%M /= size(diag) .or. A%M /= size(diag)) &
    call system_abort("Called MatrixZ_multDiag_d with mismatched sizes")

  if (this%ScaLAPACK_Info_obj%active .and. A%ScaLAPACK_Info_obj%active) then
    call matrix_product_vect_asdiagonal_sub(this%ScaLAPACK_Info_obj, this%data, A%ScaLAPACK_Info_obj, A%data, diag)
  else if (.not.this%ScaLAPACK_Info_obj%active .and. .not.A%ScaLAPACK_Info_obj%active) then
    call matrix_product_vect_asdiagonal_sub(this%data, A%data, diag)
  else
    call system_abort("Called MatrixZ_multDiag_d with mix of ScaLAPACK and non-ScaLAPACK objects")
  endif
end subroutine MatrixZ_multDiag_d

subroutine MatrixZ_multDiag_z(this, A, diag)
  type(MatrixZ), intent(inout) :: this
  type(MatrixZ), intent(in) :: A
  complex(dp), intent(in) :: diag(:)

  if (this%M /= size(diag) .or. A%M /= size(diag)) &
    call system_abort("Called MatrixZ_multDiag_z with mismatched sizes")

  if (this%ScaLAPACK_Info_obj%active .and. A%ScaLAPACK_Info_obj%active) then
    call matrix_product_vect_asdiagonal_sub(this%ScaLAPACK_Info_obj, this%data, A%ScaLAPACK_Info_obj, A%data, diag)
  else if (.not.this%ScaLAPACK_Info_obj%active .and. .not.A%ScaLAPACK_Info_obj%active) then
    call matrix_product_vect_asdiagonal_sub(this%data, A%data, diag)
  else
    call system_abort("Called MatrixZ_multDiag_z with mix of ScaLAPACK and non-ScaLAPACK objects")
  endif
end subroutine MatrixZ_multDiag_z

subroutine MatrixD_multDiagRL_d(this, A, diag)
  type(MatrixD), intent(inout) :: this
  type(MatrixD), intent(in) :: A
  real(dp), intent(in) :: diag(:)

  if (this%M /= size(diag) .or. A%M /= size(diag)) &
    call system_abort("Called MatrixD_multDiagRL_d with mismatched sizes")

  call matrix_product_vect_asdiagonal_RL_sub(this%data, A%data, diag)
end subroutine MatrixD_multDiagRL_d

subroutine MatrixZ_multDiagRL_d(this, A, diag)
  type(MatrixZ), intent(inout) :: this
  type(MatrixZ), intent(in) :: A
  real(dp), intent(in) :: diag(:)

  if (this%M /= size(diag) .or. A%M /= size(diag)) &
    call system_abort("Called MatrixZ_multDiagRL_d with mismatched sizes")

  call matrix_product_vect_asdiagonal_RL_sub(this%data, A%data, diag)
end subroutine MatrixZ_multDiagRL_d

subroutine MatrixZ_matrix_product_sub_zz(C, A, B, a_transpose, a_conjugate, b_transpose, b_conjugate)
  type(MatrixZ), intent(inout) :: C
  type(MatrixZ), intent(in) :: A, B
  logical, intent(in), optional :: a_transpose, a_conjugate, b_transpose, b_conjugate

  if (a%ScaLAPACK_Info_obj%active .and. &
      b%ScaLAPACK_Info_obj%active .and. &
      c%ScaLAPACK_Info_obj%active) then
    call matrix_product_sub(C%ScaLAPACK_Info_obj, C%data, A%ScaLAPACK_Info_obj, A%data, B%ScaLAPACK_Info_obj, B%data, &
      a_transpose, a_conjugate, b_transpose, b_conjugate)
  else if (.not.a%ScaLAPACK_Info_obj%active .and. &
	   .not.b%ScaLAPACK_Info_obj%active .and. &
	   .not.c%ScaLAPACK_Info_obj%active) then
    call matrix_product_sub(C%data, A%data, B%data, a_transpose, a_conjugate, b_transpose, b_conjugate)
  else
    call system_abort("Called MatrixZ_matric_product_sub_zz with a mix of ScaLAPACK and non-ScaLAPACK matrices")
  endif
end subroutine MatrixZ_matrix_product_sub_zz

subroutine MatrixZ_matrix_product_sub_dz(C, A, B, a_transpose, b_transpose, b_conjugate)
  type(MatrixZ), intent(inout) :: C
  type(MatrixD), intent(in) :: A
  type(MatrixZ), intent(in) :: B
  logical, intent(in), optional :: a_transpose, b_transpose, b_conjugate

  logical a_transp, b_transp, b_conjg
  real(dp), allocatable :: tC(:,:), tB(:,:)

  allocate(tC(C%l_N, C%l_M))
  call ALLOC_TRACE("MatrixZ_product_sub_dz tC", size(tC)*REAL_SIZE)
  allocate(tB(B%l_N, B%l_M))
  call ALLOC_TRACE("MatrixZ_product_sub_dz tB", size(tB)*REAL_SIZE)

  a_transp = optional_default(.false., a_transpose)
  b_transp = optional_default(.false., b_transpose)
  b_conjg  = optional_default(.false., b_conjugate)

  if (a%ScaLAPACK_Info_obj%active .and. &
      b%ScaLAPACK_Info_obj%active .and. &
      c%ScaLAPACK_Info_obj%active) then
    tB = dble(B%data)
    call matrix_product_sub(C%ScaLAPACK_Info_obj, tC, A%ScaLAPACK_Info_obj, A%data, B%ScaLAPACK_Info_obj, tB, &
      a_transp, b_transp .or. b_conjg)
    C%data = tC
    if (b_conjg) then
      tB = -aimag(B%data)
    else
      tB = aimag(B%data)
    endif
    call matrix_product_sub(C%ScaLAPACK_Info_obj, tC, A%ScaLAPACK_Info_obj, A%data, B%ScaLAPACK_Info_obj, tB, &
      a_transp, b_transp .or. b_conjg)
    c%data = c%data + cmplx(0.0_dp, tC, dp)
  else if (.not.a%ScaLAPACK_Info_obj%active .and. &
	   .not.b%ScaLAPACK_Info_obj%active .and. &
	   .not.c%ScaLAPACK_Info_obj%active) then
    tB = dble(B%data)
    call matrix_product_sub(tC, A%data, tB, a_transp, b_transp .or. b_conjg)
    C%data = tC
    if (b_conjg) then
      tB = -aimag(B%data)
    else
      tB = aimag(B%data)
    endif
    call matrix_product_sub(tC, A%data, tB, a_transp, b_transp .or. b_conjg)
    c%data = c%data + cmplx(0.0_dp, tC, dp)
  else
    call system_abort("Called MatrixZ_matric_product_sub_dz with a mix of ScaLAPACK and non-ScaLAPACK matrices")
  endif

  call DEALLOC_TRACE("MatrixZ_product_sub_dz tC", size(tC)*REAL_SIZE)
  deallocate(tC)
  call DEALLOC_TRACE("MatrixZ_product_sub_dz tB", size(tB)*REAL_SIZE)
  deallocate(tB)
end subroutine MatrixZ_matrix_product_sub_dz

subroutine MatrixZ_matrix_product_sub_zd(C, A, B, a_transpose, a_conjugate, b_transpose)
  type(MatrixZ), intent(inout) :: C
  type(MatrixZ), intent(in) :: A
  type(MatrixD), intent(in) :: B
  logical, intent(in), optional :: a_transpose, a_conjugate, b_transpose

  logical a_transp, a_conjg, b_transp
  real(dp), allocatable :: tC(:,:), tA(:,:)

  allocate(tC(C%l_N, C%l_M))
  call ALLOC_TRACE("MatrixZ_product_sub_zd tC", size(tC)*REAL_SIZE)
  allocate(tA(A%l_N, A%l_M))
  call ALLOC_TRACE("MatrixZ_product_sub_zd tA", size(tA)*REAL_SIZE)

  a_transp = optional_default(.false., a_transpose)
  a_conjg  = optional_default(.false., a_conjugate)
  b_transp = optional_default(.false., b_transpose)

  if (a%ScaLAPACK_Info_obj%active .and. &
      b%ScaLAPACK_Info_obj%active .and. &
      c%ScaLAPACK_Info_obj%active) then
    tA = dble(A%data)
    call matrix_product_sub(C%ScaLAPACK_Info_obj, tC, A%ScaLAPACK_Info_obj, tA, B%ScaLAPACK_Info_obj, B%data, &
      a_transp .or. a_conjg, b_transp)
    C%data = tC
    if (a_conjg) then
      tA = -aimag(A%data)
    else
      tA = aimag(A%data)
    endif
    call matrix_product_sub(C%ScaLAPACK_Info_obj, tC, A%ScaLAPACK_Info_obj, tA, B%ScaLAPACK_Info_obj, B%data, &
      a_transp .or. a_conjg, b_transp)
    c%data = c%data + cmplx(0.0_dp, tC, dp)
  else if (.not.a%ScaLAPACK_Info_obj%active .and. &
	   .not.b%ScaLAPACK_Info_obj%active .and. &
	   .not.c%ScaLAPACK_Info_obj%active) then
    tA = dble(A%data)
    call matrix_product_sub(tC, tA, B%data, a_transp .or. a_conjg, b_transp)
    C%data = tC
    if (a_conjg) then
      tA = -aimag(A%data)
    else
      tA = aimag(A%data)
    endif
    call matrix_product_sub(tC, tA, B%data, a_transp .or. a_conjg, b_transp)
    c%data = c%data + cmplx(0.0_dp, tC, dp)
  else
    call system_abort("Called MatrixZ_matric_product_sub_zd with a mix of ScaLAPACK and non-ScaLAPACK matrices")
  endif

  call DEALLOC_TRACE("MatrixZ_product_sub_zd tC", size(tC)*REAL_SIZE)
  deallocate(tC)
  call DEALLOC_TRACE("MatrixZ_product_sub_zd tA", size(tA)*REAL_SIZE)
  deallocate(tA)
end subroutine MatrixZ_matrix_product_sub_zd

subroutine MatrixD_matrix_product_sub_dd(C, A, B, a_transpose, b_transpose)
  type(MatrixD), intent(inout) :: C
  type(MatrixD), intent(in) :: A, B
  logical, intent(in), optional :: a_transpose, b_transpose

  if (a%ScaLAPACK_Info_obj%active .and. &
      b%ScaLAPACK_Info_obj%active .and. &
      c%ScaLAPACK_Info_obj%active) then
    call matrix_product_sub(C%ScaLAPACK_Info_obj, C%data, A%ScaLAPACK_Info_obj, A%data, B%ScaLAPACK_Info_obj, B%data, &
      a_transpose, b_transpose)
  else if (.not.a%ScaLAPACK_Info_obj%active .and. &
	   .not.b%ScaLAPACK_Info_obj%active .and. &
	   .not.c%ScaLAPACK_Info_obj%active) then
    call matrix_product_sub(C%data, A%data, B%data, a_transpose, b_transpose)
  else
    call system_abort("Called MatrixD_matric_product_sub_dd with a mix of ScaLAPACK and non-ScaLAPACK matrices")
  endif
end subroutine MatrixD_matrix_product_sub_dd

subroutine MatrixD_matrix_product_sub_dr2(C, A, B, a_transpose, b_transpose)
  type(MatrixD), intent(inout) :: C
  type(MatrixD), intent(in) :: A
  real(dp), intent(in) :: B(:,:)
  logical, intent(in), optional :: a_transpose, b_transpose

  if (.not.a%ScaLAPACK_Info_obj%active .and.  .not.c%ScaLAPACK_Info_obj%active) then
    call matrix_product_sub(C%data, A%data, B, a_transpose, b_transpose)
  else
    call system_abort("Called MatrixZ_matric_product_sub_dr2 with a mix of ScaLAPACK and non-ScaLAPACK matrices")
  endif

end subroutine MatrixD_matrix_product_sub_dr2

subroutine MatrixD_add_identity(A)
  type(MatrixD), intent(inout) :: A

  if (A%ScaLAPACK_Info_obj%active) then
    call add_identity(A%ScaLAPACK_Info_obj, A%data)
  else
    call add_identity(A%data)
  endif
end subroutine

subroutine MatrixD_scale(A, scale)
  type(MatrixD), intent(inout) :: A
  real(dp), intent(in) :: scale

  A%data = A%data * scale
end subroutine

subroutine MatrixD_transpose_sub(this, m)
  type(MatrixD), intent(inout) :: this
  type(MatrixD), intent(in) :: m

  if (this%N /= m%M .or. this%M /= m%N) call system_abort("Called MatrixD_transpose_sub with mismatched sizes "// &
    "this " // this%N // " " // this%M // " m " // m%N // " " //m%M)


  if (this%ScaLAPACK_Info_obj%active .and. m%ScaLAPACK_Info_obj%active) then
    call system_abort("MatrixD_transpose_sub not yet implemented for ScaLAPACK matrices")
  else
    this%data = transpose(m%data)
  endif
end subroutine MatrixD_transpose_sub

subroutine MatrixZ_transpose_sub(this, m)
  type(MatrixZ), intent(inout) :: this
  type(MatrixZ), intent(in) :: m

  if (this%N /= m%M .or. this%M /= m%N) call system_abort("Called MatrixZ_transpose_sub with mismatched sizes "// &
    "this " // this%N // " " // this%M // " m " // m%N // " " //m%M)


  if (this%ScaLAPACK_Info_obj%active .and. m%ScaLAPACK_Info_obj%active) then
    call system_abort("MatrixZ_transpose_sub not yet implemented for ScaLAPACK matrices")
  else
    this%data = transpose(m%data)
  endif
end subroutine MatrixZ_transpose_sub

end module matrix_module
