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
!X RS_SparseMatrix module
!X
!% Module for sparce matrix math.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module RS_SparseMatrix_module

use System_module
use MPI_Context_module
use linearalgebra_module
use Atoms_module
use ScaLAPACK_module
use Matrix_module
implicit none
private

public :: RS_SparseMatrixL
type RS_SparseMatrixL
  integer :: N = 0, N_dense_rows = 0
  integer :: n_blocks = 0, data_size = 0
  integer, allocatable :: block_size(:), dense_row_of_row(:)
  integer, allocatable :: row_indices(:), data_ptrs(:)
  integer, allocatable :: col(:)

  integer :: cur_row = 1, cur_col_offset = 0
end type RS_SparseMatrixL

public :: RS_SparseMatrixD
type RS_SparseMatrixD
  type(RS_SparseMatrixL) :: l
  real(dp), allocatable :: data(:)
end type RS_SParseMatrixD

public :: RS_SparseMatrixZ
type RS_SparseMatrixZ
  type(RS_SparseMatrixL) :: l
  complex(dp), allocatable :: data(:)
end type RS_SParseMatrixZ

public :: Initialise
interface Initialise
  module procedure RS_SparseMatrixL_Initialise_at, RS_SparseMatrixL_Initialise_prod
  module procedure RS_SparseMatrixD_Initialise_at, RS_SparseMatrixD_Initialise_prod
  module procedure RS_SparseMatrixZ_Initialise_at, RS_SparseMatrixZ_Initialise_prod
end interface Initialise


public :: Finalise
interface Finalise
  module procedure RS_SparseMatrixL_Finalise
  module procedure RS_SparseMatrixD_Finalise
  module procedure RS_SparseMatrixZ_Finalise
end interface Finalise

public :: Wipe
interface Wipe
  module procedure RS_SparseMatrixL_Wipe
  module procedure RS_SparseMatrixD_Wipe
  module procedure RS_SparseMatrixZ_Wipe
end interface Wipe

public :: Zero
interface Zero
  module procedure RS_SparseMatrixD_Zero
  module procedure RS_SparseMatrixZ_Zero
end interface Zero

public :: Print
interface Print
  module procedure RS_SparseMatrixD_Print
  module procedure RS_SparseMatrixZ_Print
end interface Print

public :: Print_simple
interface Print_simple
  module procedure RS_SparseMatrixD_Print_simple
!  module procedure RS_SparseMatrixZ_Print
end interface Print_simple

public :: copy
interface copy
  module procedure copy_dd_dsp
  module procedure copy_zd_zsp
  module procedure copy_dsp_dd
  module procedure copy_zsp_zd
end interface

public :: matrix_product_sub
interface matrix_product_sub
  module procedure matrix_product_sub_dsp_dsp_dsp, matrix_product_sub_dden_dden_dsp
  module procedure matrix_product_sub_zsp_zsp_zsp
  module procedure matrix_product_sub_dden_dden_zsp, matrix_product_sub_zden_dden_zsp, matrix_product_sub_zden_zden_zsp
  module procedure matrix_product_sub_zden_zden_dsp
end interface matrix_product_sub

public :: add_block
interface add_block
  module procedure RS_SparseMatrixD_add_block
  module procedure RS_SparseMatrixZ_add_block
end interface

public :: partial_TraceMult
interface partial_TraceMult
  module procedure RS_SparseMatrix_partial_TraceMult_dden_dsp
  module procedure RS_SparseMatrix_partial_TraceMult_dden_zsp
  module procedure RS_SparseMatrix_partial_TraceMult_zden_dsp
  module procedure RS_SparseMatrix_partial_TraceMult_zden_zsp
  module procedure RS_SparseMatrix_partial_TraceMult_dsp_zden
  module procedure RS_SparseMatrix_partial_TraceMult_dsp_dden
  module procedure RS_SparseMatrix_partial_TraceMult_zsp_zden
end interface partial_TraceMult

public :: TraceMult
interface TraceMult
  module procedure RS_SparseMatrix_TraceMult_dden_dsp
  module procedure RS_SparseMatrix_TraceMult_dden_zsp
  module procedure RS_SparseMatrix_TraceMult_zden_dsp
  module procedure RS_SparseMatrix_TraceMult_zden_zsp
  module procedure RS_SparseMatrix_TraceMult_dsp_zden
  module procedure RS_SparseMatrix_TraceMult_dsp_dden
  module procedure RS_SparseMatrix_TraceMult_zsp_zden
end interface TraceMult

public :: multDiagRL
interface multDiagRL
  module procedure RS_SparseMatrixD_multDiagRL_d, RS_SparseMatrixZ_multDiagRL_d
end interface multDiagRL

public :: check_sparse
interface check_sparse
  module procedure check_sparse_layout
end interface check_sparse

interface get_dense_block
  module procedure get_dense_blockD, get_dense_blockZ
end interface get_dense_block

contains

subroutine RS_SparseMatrixL_Initialise_at(this, at, first_orb_of_atom, cutoff, mpi_obj)
  type(RS_SparseMatrixL), intent(inout) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: first_orb_of_atom(:)
  real(dp), intent(in), optional :: cutoff
  type(MPI_context), intent(in), optional :: mpi_obj

  integer :: i, ji, j
  integer :: cur_data_pos, cur_row_offset
  integer, allocatable :: neighbour_sorted(:)
  real(dp), allocatable :: neighbour_distance_sorted(:)
  integer :: n_uniq

  real(dp) :: my_cutoff

  if (present(mpi_obj)) then
    call System_abort("No support for parallel sparse matrix")
  end if

  my_cutoff = HUGE(1.0_dp)
  if (present(cutoff)) my_cutoff = cutoff

  call Wipe(this)

  this%N = at%N
  if (size(first_orb_of_atom) /= at%N+1) call system_abort("RS_SparseMatrixL_Initialise_at size(first_orb_of_atom) = " // size(first_orb_of_atom) // "  != at%N+1 = " // (at%N+1))
  this%N_dense_rows = first_orb_of_atom(at%N+1)-1

  allocate(this%block_size(at%N))
  call ALLOC_TRACE("RS_L_init block_size", size(this%block_size)*INTEGER_SIZE)
  allocate(this%dense_row_of_row(at%N))
  call ALLOC_TRACE("RS_L_init dense_row_of_row", size(this%dense_row_of_row)*INTEGER_SIZE)
  this%block_size(:) = first_orb_of_atom(2:at%N+1)-first_orb_of_atom(1:at%N)
  this%dense_row_of_row(:) = first_orb_of_atom(1:at%N)

  allocate(this%row_indices(at%N+1))
  call ALLOC_TRACE("RS_L_init row_indices", size(this%row_indices)*INTEGER_SIZE)

  this%n_blocks = 0
  this%data_size = 0

  this%row_indices = 0
  do i=1, at%N

    allocate(neighbour_sorted(n_neighbours(at, i)))
    allocate(neighbour_distance_sorted(n_neighbours(at, i)))
    do ji=1, n_neighbours(at, i)
      neighbour_sorted(ji) = neighbour(at, i, ji, neighbour_distance_sorted(ji))
    end do
    call sort_array(neighbour_sorted, neighbour_distance_sorted)

    n_uniq = uniq_minval(neighbour_sorted, neighbour_distance_sorted)

    if(n_uniq == 0) call system_abort("atoms with no neighbours, not even itself!?!?!?")

    do ji=1, n_uniq
      j = neighbour_sorted(ji)
      if (ji > 1) then
	if (j == neighbour_sorted(ji-1)) cycle
      end if

      ! insert data for regular element
      if (neighbour_distance_sorted(ji) <= my_cutoff) then
	this%row_indices(i+1) = this%row_indices(i+1) + 1
	this%n_blocks = this%n_blocks + 1
	this%data_size = this%data_size + this%block_size(i)*this%block_size(j)
      end if

    end do

    deallocate(neighbour_sorted)
    deallocate(neighbour_distance_sorted)

  end do

  this%row_indices(1) = 1
  do i=2, at%N+1
    this%row_indices(i) = this%row_indices(i) + this%row_indices(i-1)
  end do

  if (this%n_blocks == 0) then
    allocate(this%col(1))
    allocate(this%data_ptrs(1))
  else
    allocate(this%col(this%n_blocks))
    allocate(this%data_ptrs(this%n_blocks))
  end if
  call ALLOC_TRACE("RS_L_init col", size(this%col)*INTEGER_SIZE)
  call ALLOC_TRACE("RS_L_init data_ptrs", size(this%data_ptrs)*INTEGER_SIZE)

  cur_data_pos = 1
  do i=1, at%N

    allocate(neighbour_sorted(n_neighbours(at, i)))
    allocate(neighbour_distance_sorted(n_neighbours(at, i)))
    do ji=1, n_neighbours(at, i)
      neighbour_sorted(ji) = neighbour(at, i, ji, neighbour_distance_sorted(ji))
    end do
    call sort_array(neighbour_sorted, neighbour_distance_sorted)

    n_uniq = uniq_minval(neighbour_sorted, neighbour_distance_sorted)

    if(n_uniq == 0) call system_abort("atoms with no neighbours, not even itself!?!?!?")

    cur_row_offset = 0

    do ji=1, n_uniq
      j = neighbour_sorted(ji)
      if (ji > 1) then
	if (j == neighbour_sorted(ji-1)) cycle
      end if

      ! insert data for regular element
      if (neighbour_distance_sorted(ji) <= my_cutoff) then
	this%col(this%row_indices(i)+cur_row_offset) = j
	this%data_ptrs(this%row_indices(i)+cur_row_offset) = cur_data_pos
	cur_row_offset = cur_row_offset + 1
	cur_data_pos = cur_data_pos + this%block_size(i)*this%block_size(j)
      end if

    end do

    deallocate(neighbour_sorted)
    deallocate(neighbour_distance_sorted)

  end do
end subroutine RS_SparseMatrixL_Initialise_at

function uniq_minval(i_data, r_data)
  integer, intent(inout) :: i_data(:)
  real(dp), intent(inout) :: r_data(:)
  integer :: uniq_minval

  real(dp) :: min_r_data
  integer :: data_len, i, ii, imin, imax

  data_len = size(i_data)
  if (data_len == 0) then
    uniq_minval=0
    return
  end if

  ii = 1
  do i=1, data_len
    if (i == 1) then
      imin = i
      cycle
    end if
    if (i_data(i) /= i_data(i-1)) then
      imax = i-1
      min_r_data = minval(r_data(imin:imax))
      i_data(ii) = i_data(imin)
      r_data(ii) = min_r_data
      ii = ii + 1
      imin = i
    end if
  end do
  min_r_data = minval(r_data(imin:data_len))
  i_data(ii) = i_data(imin)
  r_data(ii) = min_r_data

  uniq_minval = ii
end function uniq_minval

subroutine RS_SparseMatrixD_Initialise_at(this, at, first_orb_of_atom, cutoff, mpi_obj)
  type(RS_SparseMatrixD), intent(inout) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: first_orb_of_atom(:)
  real(dp), intent(in), optional :: cutoff
  type(MPI_context), intent(in), optional :: mpi_obj

  call Wipe(this)

  call Initialise(this%l, at, first_orb_of_atom, cutoff, mpi_obj)

  allocate(this%data(this%l%data_size))
  call ALLOC_TRACE("RS_D_init data", size(this%data)*REAL_SIZE)
end subroutine RS_SparseMatrixD_Initialise_at

subroutine RS_SparseMatrixZ_Initialise_at(this, at, first_orb_of_atom, cutoff, mpi_obj)
  type(RS_SparseMatrixZ), intent(inout) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: first_orb_of_atom(:)
  real(dp), intent(in), optional :: cutoff
  type(MPI_context), intent(in), optional :: mpi_obj

  call Wipe(this)

  call Initialise(this%l, at, first_orb_of_atom, cutoff, mpi_obj)

  allocate(this%data(this%l%data_size))
  call ALLOC_TRACE("RS_Z_init data", size(this%data)*COMPLEX_SIZE)
end subroutine RS_SparseMatrixZ_Initialise_at

subroutine RS_SparseMatrixL_Initialise_prod(this, al, bl, mpi_obj)
  type(RS_SparseMatrixL), intent(inout) :: this
  type(RS_SparseMatrixL), intent(in) :: al, bl
  type(MPI_context), intent(in), optional :: mpi_obj

  integer is, kk, ks, jj, js
  integer cur_entry, p, pd
  logical, allocatable :: col_used(:)
  integer, allocatable :: col_list(:)

  if (present(mpi_obj)) then
    call System_abort("No support for mpi_obj in sparse matrix")
  end if

  this%N = al%N
  this%N_dense_rows = al%N_dense_rows
  allocate(this%block_size(this%N))
  call ALLOC_TRACE("RS_L_init_product block_size", size(this%block_size)*INTEGER_SIZE)
  allocate(this%dense_row_of_row(this%N))
  call ALLOC_TRACE("RS_L_init_product dense_row_of_row", size(this%dense_row_of_row)*INTEGER_SIZE)
  this%block_size = al%block_size
  this%dense_row_of_row = al%dense_row_of_row

  allocate(this%row_indices(this%N+1))
  call ALLOC_TRACE("RS_L_init_product row_indices", size(this%row_indices)*INTEGER_SIZE)
  this%row_indices = 0
  this%n_blocks = 0
  this%data_size = 0

  allocate(col_used(al%N))
  allocate(col_list(al%N))
  col_used = .false.

  do is=1, al%N
    do kk=al%row_indices(is), al%row_indices(is+1)-1
      ks = al%col(kk)
      cur_entry = 0
      do jj=bl%row_indices(ks), bl%row_indices(ks+1)-1
	js = bl%col(jj)
	if (.not.col_used(js)) then
	  cur_entry = cur_entry + 1
	  col_list(cur_entry) = js
	  col_used(js) = .true.
	  this%n_blocks = this%n_blocks + 1
	  this%data_size = this%data_size + al%block_size(is)*bl%block_size(js)
	end if
      end do
    end do
    do jj=1, cur_entry
      col_used(col_list(jj)) = .false.
    end do
  end do

  allocate(this%col(this%n_blocks))
  allocate(this%data_ptrs(this%n_blocks))
  call ALLOC_TRACE("RS_L_init_product col", size(this%col)*INTEGER_SIZE)
  call ALLOC_TRACE("RS_L_init_product data_ptrs", size(this%data_ptrs)*INTEGER_SIZE)

  p = 1
  pd = 1
  do is=1, al%N
    this%row_indices(is) = p
    do kk=al%row_indices(is), al%row_indices(is+1)-1
      ks = al%col(kk)
      cur_entry = 0
      do jj=bl%row_indices(ks), bl%row_indices(ks+1)-1
	js = bl%col(jj)
	if (.not.col_used(js)) then
	  cur_entry = cur_entry + 1
	  col_list(cur_entry) = js
	  col_used(js) = .true.

	  this%col(p) = js
	  this%data_ptrs(p) = pd
	  p = p + 1
	  pd = pd + al%block_size(is)*bl%block_size(js)
	end if
      end do
    end do
    do jj=1, cur_entry
      col_used(col_list(jj)) = .false.
    end do
  end do
  this%row_indices(al%N+1) = p

  deallocate(col_used)
  deallocate(col_list)
end subroutine RS_SparseMatrixL_Initialise_prod

subroutine RS_SparseMatrixD_Initialise_prod(this, al, bl, mpi_obj)
  type(RS_SparseMatrixD), intent(inout) :: this
  type(RS_SparseMatrixL), intent(in) :: al, bl
  type(MPI_context), intent(in), optional :: mpi_obj

  call Wipe(this)

  call Initialise(this%l, al, bl)

  allocate(this%data(this%l%data_size))
  call ALLOC_TRACE("RS_D_init_product data", size(this%data)*REAL_SIZE)

end subroutine RS_SparseMatrixD_Initialise_prod

subroutine RS_SparseMatrixZ_Initialise_prod(this, al, bl, mpi_obj)
  type(RS_SparseMatrixZ), intent(inout) :: this
  type(RS_SparseMatrixL), intent(in) :: al, bl
  type(MPI_context), intent(in), optional :: mpi_obj

  call Wipe(this)

  call Initialise(this%l, al, bl)

  allocate(this%data(this%l%data_size))
  call ALLOC_TRACE("RS_Z_init_product data", size(this%data)*COMPLEX_SIZE)

end subroutine RS_SparseMatrixZ_Initialise_prod

subroutine RS_SparseMatrixL_Finalise(this)
  type(RS_SparseMatrixL), intent(inout) :: this

  call Wipe(this)
end subroutine RS_SparseMatrixL_Finalise

subroutine RS_SparseMatrixD_Finalise(this)
  type(RS_SparseMatrixD), intent(inout) :: this

  call Wipe(this)
end subroutine RS_SparseMatrixD_Finalise

subroutine RS_SparseMatrixZ_Finalise(this)
  type(RS_SparseMatrixZ), intent(inout) :: this

  call Wipe(this)
end subroutine RS_SparseMatrixZ_Finalise

subroutine RS_SparseMatrixL_Wipe(this)
  type(RS_SparseMatrixL), intent(inout) :: this

  call DEALLOC_TRACE("RS_L_wipe block_size", size(this%block_size)*INTEGER_SIZE)
  if (allocated(this%block_size)) deallocate(this%block_size)
  call DEALLOC_TRACE("RS_L_wipe dense_row_of_row", size(this%dense_row_of_row)*INTEGER_SIZE)
  if (allocated(this%dense_row_of_row)) deallocate(this%dense_row_of_row)
  call DEALLOC_TRACE("RS_L_wipe row_indices", size(this%row_indices)*INTEGER_SIZE)
  if (allocated(this%row_indices)) deallocate(this%row_indices)
  call DEALLOC_TRACE("RS_L_wipe data_ptrs", size(this%data_ptrs)*INTEGER_SIZE)
  if (allocated(this%data_ptrs)) deallocate(this%data_ptrs)
  call DEALLOC_TRACE("RS_L_wipe col", size(this%col)*INTEGER_SIZE)
  if (allocated(this%col)) deallocate(this%col)

  this%N = 0
  this%N_dense_rows = 0
  this%n_blocks = 0
  this%data_size = 0
end subroutine RS_SparseMatrixL_Wipe

subroutine RS_SparseMatrixD_Wipe(this)
  type(RS_SparseMatrixD), intent(inout) :: this

  call Wipe(this%l)
  call DEALLOC_TRACE("RS_D_wipe data", size(this%data)*REAL_SIZE)
  if (allocated(this%data)) deallocate(this%data)

end subroutine RS_SparseMatrixD_Wipe

subroutine RS_SparseMatrixZ_Wipe(this)
  type(RS_SparseMatrixZ), intent(inout) :: this

  call Wipe(this%l)
  call DEALLOC_TRACE("RS_Z_wipe data", size(this%data)*COMPLEX_SIZE)
  if (allocated(this%data)) deallocate(this%data)

end subroutine RS_SparseMatrixZ_Wipe

subroutine RS_SparseMatrixD_Zero(this, d_mask, od_mask)
  type(RS_SparseMatrixD),    intent(inout)           :: this
  logical, intent(in), optional :: d_mask(:), od_mask(:)

  integer i, ji, j, block_nr, block_nc

  if (present(d_mask) .or. present(od_mask)) then
    do i=1, this%l%N
      block_nr = this%l%block_size(i)
      do ji=this%l%row_indices(i), this%l%row_indices(i+1)-1
	j = this%l%col(ji)
	block_nc = this%l%block_size(j)
	if (present(d_mask) .and. i == j) then
	  if (d_mask(i)) this%data(this%l%data_ptrs(ji):this%l%data_ptrs(ji)+block_nr*block_nc-1) = 0.0_dp
	end if
	if (present(od_mask) .and. i /= j) then
	  if (od_mask(i) .or. od_mask(j)) this%data(this%l%data_ptrs(ji):this%l%data_ptrs(ji)+block_nr*block_nc-1) = 0.0_dp
	end if
      end do
    end do
  else
    this%data = 0.0_dp
  end if

end subroutine RS_SparseMatrixD_Zero

subroutine RS_SparseMatrixZ_Zero(this, d_mask, od_mask)
  type(RS_SparseMatrixZ),    intent(inout)           :: this
  logical, intent(in), optional :: d_mask(:), od_mask(:)

  integer i, ji, j, block_nr, block_nc

  if (present(d_mask) .or. present(od_mask)) then
    do i=1, this%l%N
      block_nr = this%l%block_size(i)
      do ji=this%l%row_indices(i), this%l%row_indices(i+1)-1
	j = this%l%col(ji)
	block_nc = this%l%block_size(j)
	if (present(d_mask) .and. i == j) then
	  if (d_mask(i)) this%data(this%l%data_ptrs(ji):this%l%data_ptrs(ji)+block_nr*block_nc-1) = 0.0_dp
	end if
	if (present(od_mask) .and. i /= j) then
	  if (od_mask(i) .or. od_mask(j)) this%data(this%l%data_ptrs(ji):this%l%data_ptrs(ji)+block_nr*block_nc-1) = 0.0_dp
	end if
      end do
    end do
  else
    this%data = 0.0_dp
  end if

end subroutine RS_SparseMatrixZ_Zero

subroutine RS_SparseMatrixD_Print(this,file)
  type(RS_SparseMatrixD),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  integer :: i, j
  integer :: block_nr, block_nc

  call Print ('RS_SparseMatrixD : N ' // this%l%N, file=file)

  call Print ('RS_SparseMatrixD : n_blocks data_size ' // this%l%n_blocks // " " // this%l%data_size, file=file)

  do i=1, this%l%N
    call Print ('RS_SparseMatrixD : row, block_size dense_row ' // i // " " // this%l%block_size(i) // " " // &
					 this%l%dense_row_of_row(i), file=file)
    do j=this%l%row_indices(i), this%l%row_indices(i+1)-1
      call Print ("RS_SparseMatrixD entry j "// j // " col " // this%l%col(j) //  " data_ptr " // this%l%data_ptrs(j), &
	file=file)
      block_nr = this%l%block_size(i)
      block_nc = this%l%block_size(this%l%col(j))
      call Print(reshape( &
	this%data(this%l%data_ptrs(j):this%l%data_ptrs(j)+block_nr*block_nc-1),(/block_nr,block_nc/) ), &
	file=file)
    end do
  end do

end subroutine RS_SparseMatrixD_Print

subroutine RS_SparseMatrixZ_Print(this,file)
  type(RS_SparseMatrixZ),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  integer :: i, j
  integer :: block_nr, block_nc

  call Print('RS_SparseMatrixZ : ', file=file)

  call Print ('RS_SparseMatrixZ : N ' //  this%l%N, file=file)

  call Print ('RS_SparseMatrixZ : n_blocks data_size ' // this%l%n_blocks // this%l%data_size, file=file)

  do i=1, this%l%N
    call Print ('RS_SparseMatrixZ : row, block_size dense_row ' // i // " " //  this%l%block_size(i) //  " " // &
					 this%l%dense_row_of_row(i), file=file)
    do j=this%l%row_indices(i), this%l%row_indices(i+1)-1
      call Print ("RS_SparseMatrixZ entry j " // j // " col " // this%l%col(j) //  " data_ptr " // this%l%data_ptrs(j), file=file)
      block_nr = this%l%block_size(i)
      block_nc = this%l%block_size(this%l%col(j))
      call Print(reshape( &
	this%data(this%l%data_ptrs(j):this%l%data_ptrs(j)+block_nr*block_nc-1),(/block_nr,block_nc/) ), file=file)
    end do
  end do

end subroutine RS_SparseMatrixZ_Print

subroutine RS_SparseMatrixD_Print_simple(this,file)
  type(RS_SparseMatrixD),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  integer :: i, tj, j, ii, jj
  integer :: block_nr, block_nc
  integer :: n_entries

  n_entries = 0
  do i=1, this%l%N
    block_nr = this%l%block_size(i)
    do j=this%l%row_indices(i), this%l%row_indices(i+1)-1
      block_nc = this%l%block_size(this%l%col(j))
      n_entries = n_entries + block_nr*block_nc
    end do
  end do

  call print(this%l%N // " " // n_entries, file=file)

  do i=1, this%l%N
    block_nr = this%l%block_size(i)
    do tj=this%l%row_indices(i), this%l%row_indices(i+1)-1
      j = this%l%col(tj)
      do ii=1, block_nr
      do jj=1, block_nc
	call print((this%l%dense_row_of_row(i)+ii-1) // " " // (this%l%dense_row_of_row(j)+jj-1) // " " // &
	  this%data(this%l%data_ptrs(tj)+ii-1+(jj-1)*block_nr), file=file)
      end do
      end do
    end do
  end do
end subroutine RS_SparseMatrixD_Print_simple

subroutine copy_dsp_dd(this, dense)
  type(RS_SparseMatrixD), intent(in) :: this
  real(dp), intent(inout), dimension(:,:) :: dense

  integer :: i, tj, j, ii, jj
  integer :: block_nr, block_nc
  integer :: n_entries

  if (any(shape(dense) /= (/this%l%N, this%l%N/))) call system_abort("size of dense matrix wrong: expecting "//this%l%N)

  n_entries = 0
  do i=1, this%l%N
    block_nr = this%l%block_size(i)
    do j=this%l%row_indices(i), this%l%row_indices(i+1)-1
      block_nc = this%l%block_size(this%l%col(j))
      n_entries = n_entries + block_nr*block_nc
    end do
  end do

  do i=1, this%l%N
    block_nr = this%l%block_size(i)
    do tj=this%l%row_indices(i), this%l%row_indices(i+1)-1
      j = this%l%col(tj)
      do ii=1, block_nr
        do jj=1, block_nc
  	       dense(this%l%dense_row_of_row(i)+ii-1, &
                 this%l%dense_row_of_row(j)+jj-1) = &
  	             this%data(this%l%data_ptrs(tj)+ii-1+(jj-1)*block_nr)
        end do
      end do
    end do
  end do
end subroutine copy_dsp_dd

subroutine copy_zsp_zd(this, dense)
  type(RS_SparseMatrixZ), intent(in) :: this
  complex(dp), intent(inout), dimension(:,:) :: dense

  integer :: i, tj, j, ii, jj
  integer :: block_nr, block_nc
  integer :: n_entries

  if (any(shape(dense) /= (/this%l%N, this%l%N/))) call system_abort("size of dense matrix wrong: expecting "//this%l%N)

  n_entries = 0
  do i=1, this%l%N
    block_nr = this%l%block_size(i)
    do j=this%l%row_indices(i), this%l%row_indices(i+1)-1
      block_nc = this%l%block_size(this%l%col(j))
      n_entries = n_entries + block_nr*block_nc
    end do
  end do

  do i=1, this%l%N
    block_nr = this%l%block_size(i)
    do tj=this%l%row_indices(i), this%l%row_indices(i+1)-1
      j = this%l%col(tj)
      do ii=1, block_nr
        do jj=1, block_nc
  	       dense(this%l%dense_row_of_row(i)+ii-1, &
                 this%l%dense_row_of_row(j)+jj-1) = &
  	             this%data(this%l%data_ptrs(tj)+ii-1+(jj-1)*block_nr)
        end do
      end do
    end do
  end do
end subroutine copy_zsp_zd

subroutine copy_dd_dsp(a, b, d_mask, od_mask)
  type(MatrixD), intent(inout) :: a
  type(RS_SparseMatrixD), intent(in) :: b
  logical, optional :: d_mask(:), od_mask(:)

  integer is, jj, js, id, jd

  a%data = 0.0_dp
  do is=1, b%l%N
    id = b%l%dense_row_of_row(is)
    do jj=b%l%row_indices(is), b%l%row_indices(is+1)-1
      js = b%l%col(jj)
      if (present(d_mask)) then
	if ((is == js) .and. .not. d_mask(is)) cycle
      endif
      if (present(od_mask)) then
	if ((is /= js) .and. .not. od_mask(is) .and. .not. od_mask(js)) cycle
      endif
      jd = b%l%dense_row_of_row(js)
      a%data(id:id+b%l%block_size(is)-1,jd:jd+b%l%block_size(js)-1) = &
	reshape(b%data(b%l%data_ptrs(jj):b%l%data_ptrs(jj)+b%l%block_size(is)*b%l%block_size(js)-1), &
	        (/b%l%block_size(is), b%l%block_size(js)/) )
    end do
  end do
end subroutine copy_dd_dsp

subroutine copy_zd_zsp(a, b, d_mask, od_mask)
  type(MatrixZ), intent(inout) :: a
  type(RS_SparseMatrixZ), intent(in) :: b
  logical, optional :: d_mask(:), od_mask(:)

  integer is, jj, js, id, jd

  a%data = 0.0_dp
  do is=1, b%l%N
    id = b%l%dense_row_of_row(is)
    do jj=b%l%row_indices(is), b%l%row_indices(is+1)-1
      js = b%l%col(jj)
      if (present(d_mask)) then
	if ((is == js) .and. .not. d_mask(is)) cycle
      endif
      if (present(od_mask)) then
	if ((is /= js) .and. .not. od_mask(is) .and. .not. od_mask(js)) cycle
      endif
      jd = b%l%dense_row_of_row(js)
      a%data(id:id+b%l%block_size(is)-1,jd:jd+b%l%block_size(js)-1) = &
	reshape(b%data(b%l%data_ptrs(jj):b%l%data_ptrs(jj)+b%l%block_size(is)*b%l%block_size(js)-1), &
	        (/b%l%block_size(is), b%l%block_size(js)/) )
    end do
  end do
end subroutine copy_zd_zsp

subroutine matrix_product_sub_dsp_dsp_dsp(c,a,b, a_transpose, b_transpose)
  type(RS_SparseMatrixD), intent(inout) :: c
  type(RS_SparseMatrixD), intent(in) :: a
  type(RS_SparseMatrixD), intent(in) :: b
  logical, intent(in), optional :: a_transpose, b_transpose

  integer is, kk, ks, jj_b, jj_c
  integer block_ni, block_nj, block_nk
  integer max_block_size
  real(dp), allocatable :: a_block(:,:)

  c%data = 0.0_dp

  if (present(a_transpose)) then
    if (a_transpose) call system_abort("Can't do sparse = sparse * sparse with transposed matrices")
  end if
  if (present(b_transpose)) then
    if (b_transpose) call system_abort("Can't do sparse = sparse * sparse with transposed matrices")
  end if

  max_block_size = maxval(a%l%block_size)
  allocate(a_block(max_block_size,max_block_size))

  do is=1, a%l%N
    block_ni = a%l%block_size(is)
    do kk=a%l%row_indices(is), a%l%row_indices(is+1)-1
      ks = a%l%col(kk)
      block_nk = a%l%block_size(ks)
      a_block(1:block_ni,1:block_nk) = reshape(a%data(a%l%data_ptrs(kk):a%l%data_ptrs(kk)+block_ni*block_nk-1), &
					       (/block_ni,block_nk/))

      jj_c = c%l%row_indices(is)
      jj_b = b%l%row_indices(ks)
      do while (jj_c < c%l%row_indices(is+1) .and. jj_b < b%l%row_indices(ks+1))
	if (b%l%col(jj_b) == c%l%col(jj_c)) then
	  block_nj = b%l%block_size(b%l%col(jj_b))
	  c%data(c%l%data_ptrs(jj_c):c%l%data_ptrs(jj_c)+block_ni*block_nj-1) = &
	    c%data(c%l%data_ptrs(jj_c):c%l%data_ptrs(jj_c)+block_ni*block_nj-1) + reshape( &
	      matmul( a_block(1:block_ni,1:block_nk), &
		      reshape(b%data(b%l%data_ptrs(jj_b):b%l%data_ptrs(jj_b)+block_nk*block_nj-1), &
			      (/block_nk,block_nj/)) ), &
	    (/block_ni*block_nj/) )
	  jj_b = jj_b + 1
	  jj_c = jj_c + 1
	else if (b%l%col(jj_b) < c%l%col(jj_c)) then
	  jj_b = jj_b + 1
	else
	  jj_c = jj_c + 1
	end if
      end do

    end do
  end do

  deallocate(a_block)

end subroutine matrix_product_sub_dsp_dsp_dsp

subroutine matrix_product_sub_zsp_zsp_zsp(c,a,b, a_transpose, a_conjugate, b_transpose, b_conjugate)
  type(RS_SparseMatrixZ), intent(inout) :: c
  type(RS_SparseMatrixZ), intent(in) :: a
  type(RS_SparseMatrixZ), intent(in) :: b
  logical, intent(in), optional :: a_transpose, a_conjugate, b_transpose, b_conjugate

  integer is, kk, ks, jj_b, jj_c
  integer block_ni, block_nj, block_nk
  integer max_block_size
  complex(dp), allocatable :: a_block(:,:)

  c%data = 0.0_dp

  if (present(a_transpose)) then
    if (a_transpose) call system_abort("Can't do sparse = sparse * sparse with transposed matrices")
  end if
  if (present(b_transpose)) then
    if (b_transpose) call system_abort("Can't do sparse = sparse * sparse with transposed matrices")
  end if
  if (present(a_conjugate)) then
    if (a_conjugate) call system_abort("Can't do sparse = sparse * sparse with conjugate matrices")
  end if
  if (present(b_conjugate)) then
    if (b_conjugate) call system_abort("Can't do sparse = sparse * sparse with conjugate matrices")
  end if

  max_block_size = maxval(a%l%block_size)
  allocate(a_block(max_block_size,max_block_size))

  do is=1, a%l%N
    block_ni = a%l%block_size(is)
    do kk=a%l%row_indices(is), a%l%row_indices(is+1)-1
      ks = a%l%col(kk)
      block_nk = a%l%block_size(ks)
      a_block(1:block_ni,1:block_nk) = reshape(a%data(a%l%data_ptrs(kk):a%l%data_ptrs(kk)+block_ni*block_nk-1), &
					       (/block_ni,block_nk/))

      jj_c = c%l%row_indices(is)
      jj_b = b%l%row_indices(ks)
      do while (jj_c < c%l%row_indices(is+1) .and. jj_b < b%l%row_indices(ks+1))
	if (b%l%col(jj_b) == c%l%col(jj_c)) then
	  block_nj = b%l%block_size(b%l%col(jj_b))
	  c%data(c%l%data_ptrs(jj_c):c%l%data_ptrs(jj_c)+block_ni*block_nj-1) = &
	    c%data(c%l%data_ptrs(jj_c):c%l%data_ptrs(jj_c)+block_ni*block_nj-1) + reshape( &
	      matmul( a_block(1:block_ni,1:block_nk), &
		      reshape(b%data(b%l%data_ptrs(jj_b):b%l%data_ptrs(jj_b)+block_nk*block_nj-1), &
			      (/block_nk,block_nj/)) ), &
	    (/block_ni*block_nj/) )
	  jj_b = jj_b + 1
	  jj_c = jj_c + 1
	else if (b%l%col(jj_b) < c%l%col(jj_c)) then
	  jj_b = jj_b + 1
	else
	  jj_c = jj_c + 1
	end if
      end do

    end do
  end do

  deallocate(a_block)

end subroutine matrix_product_sub_zsp_zsp_zsp

subroutine matrix_product_sub_dden_dden_dsp(c,a,b,a_transpose, b_transpose, diag_mask, offdiag_mask)
  type(MatrixD), intent(inout) :: c
  type(MatrixD), intent(in) :: a
  type(RS_SparseMatrixD), intent(in) :: b
  logical, intent(in), optional :: a_transpose, b_transpose
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)

  integer id, is, kd, ks, jj, jd, js
  integer block_ni, block_nj, block_nk

  c%data = 0.0_dp

  if (present(a_transpose)) then
    if (a_transpose) call system_abort("Can't do sparse = sparse * sparse with transposed matrices")
  end if
  if (present(b_transpose)) then
    if (b_transpose) call system_abort("Can't do sparse = sparse * sparse with transposed matrices")
  end if

  id = 1
  is = 1
  do while (id <= a%N)
    block_ni = b%l%block_size(is)
    kd = 1
    ks = 1
    do while (kd <= a%M)
      block_nk = b%l%block_size(ks)
      do jj=b%l%row_indices(ks), b%l%row_indices(ks+1)-1
	js = b%l%col(jj)
	jd = b%l%dense_row_of_row(js)
	block_nj = b%l%block_size(js)

	if (present(diag_mask)) then
	  if (js == ks .and. .not. diag_mask(ks)) cycle
	end if
	if (present(offdiag_mask)) then
	  if (js /= ks .and. .not. offdiag_mask(ks) .and. .not. offdiag_mask(js)) cycle
	end if

	c%data(id:id+block_ni-1,jd:jd+block_nj-1) = &
	  c%data(id:id+block_ni-1,jd:jd+block_nj-1) + &
	  matmul( a%data(id:id+block_ni-1,kd:kd+block_nk-1), &
		  reshape(b%data(b%l%data_ptrs(jj):b%l%data_ptrs(jj)+block_nk*block_nj-1),  &
			  (/block_nk,block_nj/)) )
      end do
      kd = kd + b%l%block_size(ks)
      ks = ks + 1
    end do
    id = id + b%l%block_size(is)
    is = is + 1
  end do
end subroutine

subroutine matrix_product_sub_dden_dden_zsp(c,a,b,a_transpose, b_transpose, diag_mask, offdiag_mask)
  type(MatrixD), intent(inout) :: c
  type(MatrixD), intent(in) :: a
  type(RS_SparseMatrixZ), intent(in) :: b
  logical, intent(in), optional :: a_transpose, b_transpose
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)

  integer id, is, kd, ks, jj, jd, js
  integer block_ni, block_nj, block_nk

  c%data = 0.0_dp

  if (present(a_transpose)) then
    if (a_transpose) call system_abort("Can't do sparse = sparse * sparse with transposed matrices")
  end if
  if (present(b_transpose)) then
    if (b_transpose) call system_abort("Can't do sparse = sparse * sparse with transposed matrices")
  end if

  id = 1
  is = 1
  do while (id <= a%N)
    block_ni = b%l%block_size(is)
    kd = 1
    ks = 1
    do while (kd <= a%M)
      block_nk = b%l%block_size(ks)
      do jj=b%l%row_indices(ks), b%l%row_indices(ks+1)-1
	js = b%l%col(jj)
	jd = b%l%dense_row_of_row(js)
	block_nj = b%l%block_size(js)

        if (present(diag_mask)) then
          if (js == ks .and. .not. diag_mask(ks)) cycle
        end if
        if (present(offdiag_mask)) then
          if (js /= ks .and. .not. offdiag_mask(ks) .and. .not. offdiag_mask(js)) cycle
        end if

	c%data(id:id+block_ni-1,jd:jd+block_nj-1) = &
	  c%data(id:id+block_ni-1,jd:jd+block_nj-1) + &
	  matmul( a%data(id:id+block_ni-1,kd:kd+block_nk-1), &
		  reshape(real(b%data(b%l%data_ptrs(jj):b%l%data_ptrs(jj)+block_nk*block_nj-1)),  &
			  (/block_nk,block_nj/)) )
      end do
      kd = kd + b%l%block_size(ks)
      ks = ks + 1
    end do
    id = id + b%l%block_size(is)
    is = is + 1
  end do
end subroutine

subroutine matrix_product_sub_zden_dden_zsp(c,a,b,a_transpose, b_transpose, b_conjugate, diag_mask, offdiag_mask)
  type(MatrixZ), intent(inout) :: c
  type(MatrixD), intent(in) :: a
  type(RS_SparseMatrixZ), intent(in) :: b
  logical, intent(in), optional :: a_transpose, b_transpose, b_conjugate
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)

  integer id, is, kd, ks, jj, jd, js
  integer block_ni, block_nj, block_nk

  c%data = 0.0_dp

  if (present(a_transpose)) then
    if (a_transpose) call system_abort("Can't do sparse = sparse * sparse with transposed matrices")
  end if
  if (present(b_transpose)) then
    if (b_transpose) call system_abort("Can't do sparse = sparse * sparse with transposed matrices")
  end if
  if (present(b_conjugate)) then
    if (b_conjugate) call system_abort("Can't do sparse = sparse * sparse with conjugate matrices")
  end if

  id = 1
  is = 1
  do while (id <= a%N)
    block_ni = b%l%block_size(is)
    kd = 1
    ks = 1
    do while (kd <= a%M)
      block_nk = b%l%block_size(ks)
      do jj=b%l%row_indices(ks), b%l%row_indices(ks+1)-1
	js = b%l%col(jj)
	jd = b%l%dense_row_of_row(js)
	block_nj = b%l%block_size(js)

        if (present(diag_mask)) then
          if (js == ks .and. .not. diag_mask(ks)) cycle
        end if
        if (present(offdiag_mask)) then
          if (js /= ks .and. .not. offdiag_mask(ks) .and. .not. offdiag_mask(js)) cycle
        end if

	c%data(id:id+block_ni-1,jd:jd+block_nj-1) = &
	  c%data(id:id+block_ni-1,jd:jd+block_nj-1) + &
	  matmul( a%data(id:id+block_ni-1,kd:kd+block_nk-1), &
		  reshape(b%data(b%l%data_ptrs(jj):b%l%data_ptrs(jj)+block_nk*block_nj-1),  &
			  (/block_nk,block_nj/)) )
      end do
      kd = kd + b%l%block_size(ks)
      ks = ks + 1
    end do
    id = id + b%l%block_size(is)
    is = is + 1
  end do
end subroutine

subroutine matrix_product_sub_zden_zden_zsp(c,a,b,a_transpose, a_conjugate, b_transpose, b_conjugate, diag_mask, offdiag_mask)
  type(MatrixZ), intent(inout) :: c
  type(MatrixZ), intent(in) :: a
  type(RS_SparseMatrixZ), intent(in) :: b
  logical, intent(in), optional :: a_transpose, a_conjugate, b_transpose, b_conjugate
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)

  integer id, is, kd, ks, jj, jd, js
  integer block_ni, block_nj, block_nk

  c%data = 0.0_dp

  if (present(a_transpose)) then
    if (a_transpose) call system_abort("Can't do sparse = sparse * sparse with transposed matrices")
  end if
  if (present(b_transpose)) then
    if (b_transpose) call system_abort("Can't do sparse = sparse * sparse with transposed matrices")
  end if
  if (present(a_conjugate)) then
    if (a_conjugate) call system_abort("Can't do sparse = sparse * sparse with conjugate matrices")
  end if
  if (present(b_conjugate)) then
    if (b_conjugate) call system_abort("Can't do sparse = sparse * sparse with conjugate matrices")
  end if

  id = 1
  is = 1
  do while (id <= a%N)
    block_ni = b%l%block_size(is)
    kd = 1
    ks = 1
    do while (kd <= a%M)
      block_nk = b%l%block_size(ks)
      do jj=b%l%row_indices(ks), b%l%row_indices(ks+1)-1
	js = b%l%col(jj)
	jd = b%l%dense_row_of_row(js)
	block_nj = b%l%block_size(js)

        if (present(diag_mask)) then
          if (js == ks .and. .not. diag_mask(ks)) cycle
        end if
        if (present(offdiag_mask)) then
          if (js /= ks .and. .not. offdiag_mask(ks) .and. .not. offdiag_mask(js)) cycle
        end if

	c%data(id:id+block_ni-1,jd:jd+block_nj-1) = &
	  c%data(id:id+block_ni-1,jd:jd+block_nj-1) + &
	  matmul( a%data(id:id+block_ni-1,kd:kd+block_nk-1), &
		  reshape(b%data(b%l%data_ptrs(jj):b%l%data_ptrs(jj)+block_nk*block_nj-1),  &
			  (/block_nk,block_nj/)) )
      end do
      kd = kd + b%l%block_size(ks)
      ks = ks + 1
    end do
    id = id + b%l%block_size(is)
    is = is + 1
  end do
end subroutine

subroutine matrix_product_sub_zden_zden_dsp(c,a,b,a_transpose, a_conjugate, b_transpose, diag_mask, offdiag_mask)
  type(MatrixZ), intent(inout) :: c
  type(MatrixZ), intent(in) :: a
  type(RS_SparseMatrixD), intent(in) :: b
  logical, intent(in), optional :: a_transpose, a_conjugate, b_transpose
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)

  integer id, is, kd, ks, jj, jd, js
  integer block_ni, block_nj, block_nk

  if (present(a_transpose)) then
    if (a_transpose) call system_abort("Can't do dense = dense * sparse with transposed matrices")
  end if
  if (present(b_transpose)) then
    if (b_transpose) call system_abort("Can't do dense = dense * sparse with transposed matrices")
  end if
  if (present(a_conjugate)) then
    if (a_conjugate) call system_abort("Can't do dense = dense * sparse with conjugate matrices")
  end if

  c%data = 0.0_dp

  id = 1
  is = 1
  do while (id <= a%N)
    block_ni = b%l%block_size(is)
    kd = 1
    ks = 1
    do while (kd <= a%M)
      block_nk = b%l%block_size(ks)
      do jj=b%l%row_indices(ks), b%l%row_indices(ks+1)-1
	js = b%l%col(jj)
	jd = b%l%dense_row_of_row(js)
	block_nj = b%l%block_size(js)

        if (present(diag_mask)) then
          if (js == ks .and. .not. diag_mask(ks)) cycle
        end if
        if (present(offdiag_mask)) then
          if (js /= ks .and. .not. offdiag_mask(ks) .and. .not. offdiag_mask(js)) cycle
        end if

	c%data(id:id+block_ni-1,jd:jd+block_nj-1) = &
	  c%data(id:id+block_ni-1,jd:jd+block_nj-1) + &
	  matmul( a%data(id:id+block_ni-1,kd:kd+block_nk-1), &
		  reshape(b%data(b%l%data_ptrs(jj):b%l%data_ptrs(jj)+block_nk*block_nj-1),  &
			  (/block_nk,block_nj/)) )
      end do
      kd = kd + b%l%block_size(ks)
      ks = ks + 1
    end do
    id = id + b%l%block_size(is)
    is = is + 1
  end do
end subroutine

subroutine RS_SparseMatrixD_add_block(this, block, block_nr, block_nc, at_row, at_col)
  type(RS_SparseMatrixD), intent(inout) :: this
  real(dp), intent(in) :: block(:,:)
  integer, intent(in) :: block_nr, block_nc
  integer, intent(in), optional :: at_row, at_col

  integer n_cols
  integer cur_col, col_offset

  if (.not. present(at_row) .or. .not. present(at_col)) then
    call System_abort("need at_row and at_col for add_block to RS_SparseMatrixD")
  end if

  if (at_row > this%l%N .or. at_col > this%l%N) then
      call system_abort("RS_SparseMatrixD_add_block tried to add block outside of matrix bounds "//at_row//","//at_col//" "//this%l%N)
  end if

  if (this%l%cur_row /= at_row) then
    this%l%cur_row = at_row
    this%l%cur_col_offset = 0
  end if

  cur_col = this%l%col(this%l%row_indices(this%l%cur_row)+this%l%cur_col_offset)
  if (cur_col == at_col) then
    call add_block_d(this, this%l%cur_row, this%l%cur_col_offset, block_nr, block_nc, block)
    return
  else if (cur_col > at_col) then
    do col_offset = this%l%cur_col_offset-1, 0, -1
      cur_col = this%l%col(this%l%row_indices(this%l%cur_row)+col_offset)
      if (cur_col == at_col) then
	this%l%cur_col_offset = col_offset
	call add_block_d(this, this%l%cur_row, this%l%cur_col_offset, block_nr, block_nc, block)
	return
      end if
    end do
  else
    n_cols = this%l%row_indices(this%l%cur_row+1)-this%l%row_indices(this%l%cur_row)
    do col_offset = this%l%cur_col_offset+1, n_cols-1
      cur_col = this%l%col(this%l%row_indices(this%l%cur_row)+col_offset)
      if (cur_col == at_col) then
	this%l%cur_col_offset = col_offset
	call add_block_d(this, this%l%cur_row, this%l%cur_col_offset, block_nr, block_nc, block)
	return
      end if
    end do
  end if

end subroutine RS_SparseMatrixD_add_block

subroutine RS_SparseMatrixZ_add_block(this, block, block_nr, block_nc, at_row, at_col)
  type(RS_SparseMatrixZ), intent(inout) :: this
  complex(dp), intent(in) :: block(:,:)
  integer, intent(in) :: block_nr, block_nc
  integer, intent(in) :: at_row, at_col

  integer n_cols
  integer cur_col, col_offset

  if (at_row > this%l%N .or. at_col > this%l%N) then
      call system_abort("RS_SparseMatrixZ_add_block_tried to add block outside of matrix bounds "//at_row//","//at_col//" "//this%l%N)
  end if

  if (this%l%cur_row /= at_row) then
    this%l%cur_row = at_row
    this%l%cur_col_offset = 0
  end if

  cur_col = this%l%col(this%l%row_indices(this%l%cur_row)+this%l%cur_col_offset)
  if (cur_col == at_col) then
    call add_block_z(this, this%l%cur_row, this%l%cur_col_offset, block_nr, block_nc, block)
    return
  else if (cur_col > at_col) then
    do col_offset = this%l%cur_col_offset-1, 0, -1
      cur_col = this%l%col(this%l%row_indices(this%l%cur_row)+col_offset)
      if (cur_col == at_col) then
	this%l%cur_col_offset = col_offset
	call add_block_z(this, this%l%cur_row, this%l%cur_col_offset, block_nr, block_nc, block)
	return
      end if
    end do
  else
    n_cols = this%l%row_indices(this%l%cur_row+1)-this%l%row_indices(this%l%cur_row)
    do col_offset = this%l%cur_col_offset+1, n_cols-1
      cur_col = this%l%col(this%l%row_indices(this%l%cur_row)+col_offset)
      if (cur_col == at_col) then
	this%l%cur_col_offset = col_offset
	call add_block_z(this, this%l%cur_row, this%l%cur_col_offset, block_nr, block_nc, block)
	return
      end if
    end do
  end if

end subroutine RS_SparseMatrixZ_add_block

subroutine add_block_d(this, cur_row, cur_col_offset, in_nr, in_nc, block)
  type(RS_SparseMatrixD), intent(inout) :: this
  integer, intent(in) :: cur_row, cur_col_offset
  integer, intent(in) :: in_nr, in_nc
  real(dp), intent(in) :: block(:,:)

  integer ptr_i, this_nr, this_nc

  ptr_i = this%l%row_indices(cur_row)+cur_col_offset
  this_nr = this%l%block_size(cur_row)
  this_nc = this%l%block_size(this%l%col(ptr_i))
  if (in_nr /= this_nr .or. in_nc /= this_nc) then
    call System_abort("add_block_d tried to add block of wrong shape in_nr,nc "//in_nr// " "//in_nc//" sp_nr,nc "//this_nr//" "//this_nc)
  end if
  this%data(this%l%data_ptrs(ptr_i):this%l%data_ptrs(ptr_i)+this_nr*this_nc-1) = &
    this%data(this%l%data_ptrs(ptr_i):this%l%data_ptrs(ptr_i)+this_nr*this_nc-1) + &
    reshape( block(1:in_nr,1:in_nc), (/this_nr*this_nc/) )
end subroutine add_block_d

subroutine add_block_z(this, cur_row, cur_col_offset, in_nr, in_nc, block)
  type(RS_SparseMatrixZ), intent(inout) :: this
  integer, intent(in) :: cur_row, cur_col_offset
  integer, intent(in) :: in_nr, in_nc
  complex(dp), intent(in) :: block(:,:)

  integer ptr_i, this_nr, this_nc

  ptr_i = this%l%row_indices(cur_row)+cur_col_offset
  this_nr = this%l%block_size(cur_row)
  this_nc = this%l%block_size(this%l%col(ptr_i))
  if (in_nr /= this_nr .or. in_nc /= this_nc) then
    call System_abort("add_block_z tried to add block of wrong shape in_nr,nc "//in_nr// " "//in_nc//" sp_nr,nc "//this_nr//" "//this_nc)
  end if
  this%data(this%l%data_ptrs(ptr_i):this%l%data_ptrs(ptr_i)+this_nr*this_nc-1) = &
    this%data(this%l%data_ptrs(ptr_i):this%l%data_ptrs(ptr_i)+this_nr*this_nc-1) + &
    reshape( block(1:in_nr,1:in_nc), (/this_nr*this_nc/) )
end subroutine add_block_z

function RS_SparseMatrix_partial_TraceMult_dden_dsp(a, b, a_T, b_T, diag_mask, offdiag_mask) result(v)
  type(MatrixD), intent(in) :: a
  type(RS_SparseMatrixD), intent(in) ::  b
  logical, intent(in), optional :: a_T, b_T
  logical, intent(in), optional, target :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: v(a%N)

  logical u_a_T, u_b_T
  integer i, ji, j, tt
  integer block_ni, block_nj, max_block_size
  real(dp), allocatable :: a_block(:,:), b_block(:,:)
  logical, pointer :: d_mask(:), od_mask(:)
  logical use_a_ij_transpose
  logical have_any

  max_block_size = maxval(b%l%block_size)
  allocate(a_block(max_block_size,max_block_size))
  allocate(b_block(max_block_size,max_block_size))
  if (present(diag_mask)) then
    d_mask => diag_mask
  else
    allocate(d_mask(b%l%N))
    d_mask = .true.
  end if
  if (present(offdiag_mask)) then
    od_mask => offdiag_mask
  else
    allocate(od_mask(b%l%N))
    od_mask = .true.
  end if

  u_a_T = optional_default(.false., a_T)
  u_b_T = optional_default(.false., b_T)

  use_a_ij_transpose = .false.
  if ((u_a_T .and. u_b_T) .or. (.not. u_a_T .and. .not. u_b_T)) use_a_ij_transpose = .true.

  v = 0.0_dp
  do i=1, b%l%N
    block_ni = b%l%block_size(i)
    do ji=b%l%row_indices(i), b%l%row_indices(i+1)-1
      j = b%l%col(ji)
      if ((i == j .and. .not. d_mask(i)) .or. ((i /= j) .and. .not. od_mask(i) .and. .not. od_mask(j))) cycle
      block_nj = b%l%block_size(j)
      have_any = get_dense_block(a%data, a%scalapack_info_obj, use_a_ij_transpose, i, j, b%l, block_ni, block_nj, a_block)

      if (have_any) then
	b_block(1:block_ni,1:block_nj) = reshape(b%data(b%l%data_ptrs(ji):b%l%data_ptrs(ji)+block_ni*block_nj-1),(/block_ni,block_nj/))

	if (u_b_T) then
	  do tt=1, block_nj
	    v(b%l%dense_row_of_row(i):b%l%dense_row_of_row(i)+block_ni-1) = &
	      v(b%l%dense_row_of_row(i):b%l%dense_row_of_row(i)+block_ni-1) + &
	      a_block(1:block_ni,tt)*b_block(1:block_ni,tt)
	  end do
	else
	  do tt=1, block_ni
	    v(b%l%dense_row_of_row(j):b%l%dense_row_of_row(j)+block_nj-1) = &
	      v(b%l%dense_row_of_row(j):b%l%dense_row_of_row(j)+block_nj-1) + &
	      a_block(tt,1:block_nj)*b_block(tt,1:block_nj)
	  end do
	end if
      end if

    end do
  end do

  deallocate(a_block)
  deallocate(b_block)
  if (.not.present(diag_mask)) deallocate(d_mask)
  if (.not.present(offdiag_mask)) deallocate(od_mask)

end function RS_SparseMatrix_partial_TraceMult_dden_dsp

function RS_SparseMatrix_partial_TraceMult_dden_zsp(a, b, a_T, b_H, diag_mask, offdiag_mask) result(v)
  type(MatrixD), intent(in) :: a
  type(RS_SparseMatrixZ), intent(in) ::  b
  logical, intent(in), optional :: a_T, b_H
  logical, intent(in), optional, target :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: v(a%N)

  logical u_a_T, u_b_H
  integer i, jj, j, tt
  integer block_ni, block_nj, max_block_size
  real(dp), allocatable :: a_block(:,:)
  complex(dp), allocatable :: b_block(:,:)
  logical, pointer :: d_mask(:), od_mask(:)
  logical :: use_a_ij_transpose
  logical have_any

  max_block_size = maxval(b%l%block_size)
  allocate(a_block(max_block_size,max_block_size))
  allocate(b_block(max_block_size,max_block_size))
  if (present(diag_mask)) then
    d_mask => diag_mask
  else
    allocate(d_mask(b%l%N))
    d_mask = .true.
  end if
  if (present(offdiag_mask)) then
    od_mask => offdiag_mask
  else
    allocate(od_mask(b%l%N))
    od_mask = .true.
  end if

  u_a_T = optional_default(.false., a_T)
  u_b_H = optional_default(.false., b_H)

  use_a_ij_transpose = .false.
  if ((u_a_T .and. u_b_H) .or. (.not. u_a_T .and. .not. u_b_H)) use_a_ij_transpose = .true.

  v = 0.0_dp
  do i=1, b%l%N
    block_ni = b%l%block_size(i)
    do jj=b%l%row_indices(i), b%l%row_indices(i+1)-1
      j = b%l%col(jj)
      if ((i == j .and. .not. d_mask(i)) .or. ((i /= j) .and. .not. od_mask(i) .and. .not. od_mask(j))) cycle
      block_nj = b%l%block_size(j)
      have_any = get_dense_block(a%data, a%scalapack_info_obj, use_a_ij_transpose, i, j, b%l, block_ni, block_nj, a_block)

      if (have_any) then
	b_block(1:block_ni,1:block_nj) = reshape(b%data(b%l%data_ptrs(jj):b%l%data_ptrs(jj)+block_ni*block_nj-1),(/block_ni,block_nj/))
	if (u_b_H) b_block = conjg(b_block)

	if (u_b_H) then
	  do tt=1, block_nj
	    v(b%l%dense_row_of_row(i):b%l%dense_row_of_row(i)+block_ni-1) = &
	      v(b%l%dense_row_of_row(i):b%l%dense_row_of_row(i)+block_ni-1) + &
	      a_block(1:block_ni,tt)*b_block(1:block_ni,tt)
	  end do
	else
	  do tt=1, block_ni
	    v(b%l%dense_row_of_row(j):b%l%dense_row_of_row(j)+block_nj-1) = &
	      v(b%l%dense_row_of_row(j):b%l%dense_row_of_row(j)+block_nj-1) + &
	      a_block(tt,1:block_nj)*b_block(tt,1:block_nj)
	  end do
	end if
      end if

    end do
  end do

  deallocate(a_block)
  deallocate(b_block)
  if (.not.present(diag_mask)) deallocate(d_mask)
  if (.not.present(offdiag_mask)) deallocate(od_mask)

end function RS_SparseMatrix_partial_TraceMult_dden_zsp

function RS_SparseMatrix_partial_TraceMult_zden_dsp(a, b, a_H, b_T, diag_mask, offdiag_mask) result(v)
  type(MatrixZ), intent(in) :: a
  type(RS_SparseMatrixD), intent(in) ::  b
  logical, intent(in), optional :: a_H, b_T
  logical, intent(in), optional, target :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: v(a%N)

  logical u_a_H, u_b_T
  integer i, jj, j, tt
  integer block_ni, block_nj, max_block_size
  complex(dp), allocatable :: a_block(:,:)
  real(dp), allocatable :: b_block(:,:)
  logical, pointer :: d_mask(:), od_mask(:)
  logical :: use_a_ij_transpose
  logical have_any

  max_block_size = maxval(b%l%block_size)
  allocate(a_block(max_block_size,max_block_size))
  allocate(b_block(max_block_size,max_block_size))
  if (present(diag_mask)) then
    d_mask => diag_mask
  else
    allocate(d_mask(b%l%N))
    d_mask = .true.
  end if
  if (present(offdiag_mask)) then
    od_mask => offdiag_mask
  else
    allocate(od_mask(b%l%N))
    od_mask = .true.
  end if

  u_a_H = optional_default(.false., a_H)
  u_b_T = optional_default(.false., b_T)

  use_a_ij_transpose = .false.
  if ((u_a_H .and. u_b_T) .or. (.not. u_a_H .and. .not. u_b_T)) use_a_ij_transpose = .true.

  v = 0.0_dp
  do i=1, b%l%N
    block_ni = b%l%block_size(i)
    do jj=b%l%row_indices(i), b%l%row_indices(i+1)-1
      j = b%l%col(jj)
      if ((i == j .and. .not. d_mask(i)) .or. ((i /= j) .and. .not. od_mask(i) .and. .not. od_mask(j))) cycle
      block_nj = b%l%block_size(j)
      have_any = get_dense_block(a%data, a%scalapack_info_obj, use_a_ij_transpose, i, j, b%l, block_ni, block_nj, a_block)
      if (u_a_H) a_block = conjg(a_block)

      if (have_any) then
	b_block(1:block_ni,1:block_nj) = reshape(b%data(b%l%data_ptrs(jj):b%l%data_ptrs(jj)+block_ni*block_nj-1),(/block_ni,block_nj/))

	if (u_b_T) then
	  do tt=1, block_nj
	    v(b%l%dense_row_of_row(i):b%l%dense_row_of_row(i)+block_ni-1) = &
	      v(b%l%dense_row_of_row(i):b%l%dense_row_of_row(i)+block_ni-1) + &
	      a_block(1:block_ni,tt)*b_block(1:block_ni,tt)
	  end do
	else
	  do tt=1, block_ni
	    v(b%l%dense_row_of_row(j):b%l%dense_row_of_row(j)+block_nj-1) = &
	      v(b%l%dense_row_of_row(j):b%l%dense_row_of_row(j)+block_nj-1) + &
	      a_block(tt,1:block_nj)*b_block(tt,1:block_nj)
	  end do
	end if
      end if

    end do
  end do

  deallocate(a_block)
  deallocate(b_block)
  if (.not.present(diag_mask)) deallocate(d_mask)
  if (.not.present(offdiag_mask)) deallocate(od_mask)

end function RS_SparseMatrix_partial_TraceMult_zden_dsp

function RS_SparseMatrix_partial_TraceMult_zden_zsp(a, b, a_H, b_H, diag_mask, offdiag_mask) result(v)
  type(MatrixZ), intent(in) :: a
  type(RS_SparseMatrixZ), intent(in) ::  b
  logical, intent(in), optional :: a_H, b_H
  logical, intent(in), optional, target :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: v(a%N)

  logical u_a_H, u_b_H
  integer i, jj, j, tt
  integer block_ni, block_nj, max_block_size
  complex(dp), allocatable :: a_block(:,:), b_block(:,:)
  logical, pointer :: d_mask(:), od_mask(:)
  logical :: use_a_ij_transpose
  logical have_any

  max_block_size = maxval(b%l%block_size)
  allocate(a_block(max_block_size,max_block_size))
  allocate(b_block(max_block_size,max_block_size))
  if (present(diag_mask)) then
    d_mask => diag_mask
  else
    allocate(d_mask(b%l%N))
    d_mask = .true.
  end if
  if (present(offdiag_mask)) then
    od_mask => offdiag_mask
  else
    allocate(od_mask(b%l%N))
    od_mask = .true.
  end if

  u_a_H = optional_default(.false., a_H)
  u_b_H = optional_default(.false., b_H)

  use_a_ij_transpose = .false.
  if ((u_a_H .and. u_b_H) .or. (.not. u_a_H .and. .not. u_b_H)) use_a_ij_transpose = .true.

  v = 0.0_dp
  do i=1, b%l%N
    block_ni = b%l%block_size(i)
    do jj=b%l%row_indices(i), b%l%row_indices(i+1)-1
      j = b%l%col(jj)
      if ((i == j .and. .not. d_mask(i)) .or. ((i /= j) .and. .not. od_mask(i) .and. .not. od_mask(j))) cycle
      block_nj = b%l%block_size(j)
      have_any = get_dense_block(a%data, a%scalapack_info_obj, use_a_ij_transpose, i, j, b%l, block_ni, block_nj, a_block)
      if (u_a_H) a_block = conjg(a_block)

      if (have_any) then
	b_block(1:block_ni,1:block_nj) = reshape(b%data(b%l%data_ptrs(jj):b%l%data_ptrs(jj)+block_ni*block_nj-1),(/block_ni,block_nj/))
	if (u_b_H) b_block = conjg(b_block)

	if (u_b_H) then
	  do tt=1, block_nj
	    v(b%l%dense_row_of_row(i):b%l%dense_row_of_row(i)+block_ni-1) = &
	      v(b%l%dense_row_of_row(i):b%l%dense_row_of_row(i)+block_ni-1) + &
	      a_block(1:block_ni,tt)*b_block(1:block_ni,tt)
	  end do
	else
	  do tt=1, block_ni
	    v(b%l%dense_row_of_row(j):b%l%dense_row_of_row(j)+block_nj-1) = &
	      v(b%l%dense_row_of_row(j):b%l%dense_row_of_row(j)+block_nj-1) + &
	      a_block(tt,1:block_nj)*b_block(tt,1:block_nj)
	  end do
	end if
      end if

    end do
  end do

  deallocate(a_block)
  deallocate(b_block)
  if (.not.present(diag_mask)) deallocate(d_mask)
  if (.not.present(offdiag_mask)) deallocate(od_mask)

end function RS_SparseMatrix_partial_TraceMult_zden_zsp

function RS_SparseMatrix_partial_TraceMult_dsp_zden(a, b, a_T, b_H, diag_mask, offdiag_mask) result(v)
  type(RS_SparseMatrixD), intent(in) ::  a
  type(MatrixZ), intent(in) :: b
  logical, intent(in), optional :: a_T, b_H
  logical, intent(in), optional, target :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: v(b%N)

  logical u_a_T, u_b_H
  integer i, jj, j, tt
  integer block_ni, block_nj, max_block_size
  real(dp), allocatable :: a_block(:,:)
  complex(dp), allocatable :: b_block(:,:)
  logical, pointer :: d_mask(:), od_mask(:)
  logical use_b_ij_transpose
  logical have_any

  max_block_size = maxval(a%l%block_size)
  allocate(b_block(max_block_size,max_block_size))
  allocate(a_block(max_block_size,max_block_size))
  if (present(diag_mask)) then
    d_mask => diag_mask
  else
    allocate(d_mask(a%l%N))
    d_mask = .true.
  end if
  if (present(offdiag_mask)) then
    od_mask => offdiag_mask
  else
    allocate(od_mask(a%l%N))
    od_mask = .true.
  end if

  u_a_T = optional_default(.false., a_T)
  u_b_H = optional_default(.false., b_H)

  use_b_ij_transpose = .false.
  if ((u_a_T .and. u_b_H) .or. (.not. u_a_T .and. .not. u_b_H)) use_b_ij_transpose = .true.

  v = 0.0_dp
  do i=1, a%l%N
    block_ni = a%l%block_size(i)
    do jj=a%l%row_indices(i), a%l%row_indices(i+1)-1
      j = a%l%col(jj)
      if ((i == j .and. .not. d_mask(i)) .or. ((i /= j) .and. .not. od_mask(i) .and. .not. od_mask(j))) cycle
      block_nj = a%l%block_size(j)
      have_any = get_dense_block(b%data, b%scalapack_info_obj, use_b_ij_transpose, i, j, a%l, block_ni, block_nj, b_block)
      if (u_b_H) b_block = conjg(b_block)

      if (have_any) then
	a_block(1:block_ni,1:block_nj) = reshape(a%data(a%l%data_ptrs(jj):a%l%data_ptrs(jj)+block_ni*block_nj-1),(/block_ni,block_nj/))

	if (u_a_T) then
	  do tt=1, block_ni
	    v(a%l%dense_row_of_row(j):a%l%dense_row_of_row(j)+block_nj-1) = &
	      v(a%l%dense_row_of_row(j):a%l%dense_row_of_row(j)+block_nj-1) + &
	      a_block(tt,1:block_nj)*b_block(tt,1:block_nj)
	  end do
	else
	  do tt=1, block_nj
	    v(a%l%dense_row_of_row(i):a%l%dense_row_of_row(i)+block_ni-1) = &
	      v(a%l%dense_row_of_row(i):a%l%dense_row_of_row(i)+block_ni-1) + &
	      a_block(1:block_ni,tt)*b_block(1:block_ni,tt)
	  end do
	end if
      end if

    end do
  end do

  deallocate(b_block)
  deallocate(a_block)
  if (.not.present(diag_mask)) deallocate(d_mask)
  if (.not.present(offdiag_mask)) deallocate(od_mask)

end function RS_SparseMatrix_partial_TraceMult_dsp_zden

function RS_SparseMatrix_partial_TraceMult_zsp_zden(a, b, a_H, b_H, diag_mask, offdiag_mask) result(v)
  type(RS_SparseMatrixZ), intent(in) ::  a
  type(MatrixZ), intent(in) :: b
  logical, intent(in), optional :: a_H, b_H
  logical, intent(in), optional, target :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: v(b%N)

  logical u_a_H, u_b_H
  integer i, jj, j, tt
  integer block_ni, block_nj, max_block_size
  complex(dp), allocatable :: a_block(:,:)
  complex(dp), allocatable :: b_block(:,:)
  logical, pointer :: d_mask(:), od_mask(:)
  logical use_b_ij_transpose
  logical have_any

  max_block_size = maxval(a%l%block_size)
  allocate(b_block(max_block_size,max_block_size))
  allocate(a_block(max_block_size,max_block_size))
  if (present(diag_mask)) then
    d_mask => diag_mask
  else
    allocate(d_mask(a%l%N))
    d_mask = .true.
  end if
  if (present(offdiag_mask)) then
    od_mask => offdiag_mask
  else
    allocate(od_mask(a%l%N))
    od_mask = .true.
  end if

  u_a_H = optional_default(.false., a_H)
  u_b_H = optional_default(.false., b_H)

  use_b_ij_transpose = .false.
  if ((u_a_H .and. u_b_H) .or. (.not. u_a_H .and. .not. u_b_H)) use_b_ij_transpose = .true.

  v = 0.0_dp
  do i=1, a%l%N
    block_ni = a%l%block_size(i)
    do jj=a%l%row_indices(i), a%l%row_indices(i+1)-1
      j = a%l%col(jj)
      if ((i == j .and. .not. d_mask(i)) .or. ((i /= j) .and. .not. od_mask(i) .and. .not. od_mask(j))) cycle
      block_nj = a%l%block_size(j)
      have_any = get_dense_block(b%data, b%scalapack_info_obj, use_b_ij_transpose, i, j, a%l, block_ni, block_nj, b_block)
      if (u_b_H) b_block = conjg(b_block)

      if (have_any) then
	a_block(1:block_ni,1:block_nj) = reshape(a%data(a%l%data_ptrs(jj):a%l%data_ptrs(jj)+block_ni*block_nj-1),(/block_ni,block_nj/))
	if (u_a_H) a_block = conjg(a_block)

	if (u_a_H) then
	  do tt=1, block_ni
	    v(a%l%dense_row_of_row(j):a%l%dense_row_of_row(j)+block_nj-1) = &
	      v(a%l%dense_row_of_row(j):a%l%dense_row_of_row(j)+block_nj-1) + &
	      a_block(tt,1:block_nj)*b_block(tt,1:block_nj)
	  end do
	else
	  do tt=1, block_nj
	    v(a%l%dense_row_of_row(i):a%l%dense_row_of_row(i)+block_ni-1) = &
	      v(a%l%dense_row_of_row(i):a%l%dense_row_of_row(i)+block_ni-1) + &
	      a_block(1:block_ni,tt)*b_block(1:block_ni,tt)
	  end do
	end if
      end if

    end do
  end do

  deallocate(b_block)
  deallocate(a_block)
  if (.not.present(diag_mask)) deallocate(d_mask)
  if (.not.present(offdiag_mask)) deallocate(od_mask)

end function RS_SparseMatrix_partial_TraceMult_zsp_zden

function RS_SparseMatrix_partial_TraceMult_dsp_dden(a, b, a_T, b_T, diag_mask, offdiag_mask) result(v)
  type(RS_SparseMatrixD), intent(in) ::  a
  type(MatrixD), intent(in) :: b
  logical, intent(in), optional :: a_T, b_T
  logical, intent(in), optional, target :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: v(b%N)

  logical u_a_T, u_b_T
  integer i, jj, j, tt
  integer block_ni, block_nj, max_block_size
  real(dp), allocatable :: a_block(:,:)
  real(dp), allocatable :: b_block(:,:)
  logical, pointer :: d_mask(:), od_mask(:)
  logical use_b_ij_transpose
  logical have_any

  max_block_size = maxval(a%l%block_size)
  allocate(b_block(max_block_size,max_block_size))
  allocate(a_block(max_block_size,max_block_size))
  if (present(diag_mask)) then
    d_mask => diag_mask
  else
    allocate(d_mask(a%l%N))
    d_mask = .true.
  end if
  if (present(offdiag_mask)) then
    od_mask => offdiag_mask
  else
    allocate(od_mask(a%l%N))
    od_mask = .true.
  end if

  u_a_T = optional_default(.false., a_T)
  u_b_T = optional_default(.false., b_T)

  use_b_ij_transpose = .false.
  if ((u_a_T .and. u_b_T) .or. (.not. u_a_T .and. .not. u_b_T)) use_b_ij_transpose = .true.

  v = 0.0_dp
  do i=1, a%l%N
    block_ni = a%l%block_size(i)
    do jj=a%l%row_indices(i), a%l%row_indices(i+1)-1
      j = a%l%col(jj)
      if ((i == j .and. .not. d_mask(i)) .or. ((i /= j) .and. .not. od_mask(i) .and. .not. od_mask(j))) cycle
      block_nj = a%l%block_size(j)
      have_any = get_dense_block(b%data, b%scalapack_info_obj, use_b_ij_transpose, i, j, a%l, block_ni, block_nj, b_block)

      if (have_any) then
	a_block(1:block_ni,1:block_nj) = reshape(a%data(a%l%data_ptrs(jj):a%l%data_ptrs(jj)+block_ni*block_nj-1),(/block_ni,block_nj/))

	if (u_a_T) then
	  do tt=1, block_ni
	    v(a%l%dense_row_of_row(j):a%l%dense_row_of_row(j)+block_nj-1) = &
	      v(a%l%dense_row_of_row(j):a%l%dense_row_of_row(j)+block_nj-1) + &
	      a_block(tt,1:block_nj)*b_block(tt,1:block_nj)
	  end do
	else
	  do tt=1, block_nj
	    v(a%l%dense_row_of_row(i):a%l%dense_row_of_row(i)+block_ni-1) = &
	      v(a%l%dense_row_of_row(i):a%l%dense_row_of_row(i)+block_ni-1) + &
	      a_block(1:block_ni,tt)*b_block(1:block_ni,tt)
	  end do
	end if
      end if

    end do
  end do

  deallocate(b_block)
  deallocate(a_block)
  if (.not.present(diag_mask)) deallocate(d_mask)
  if (.not.present(offdiag_mask)) deallocate(od_mask)

end function RS_SparseMatrix_partial_TraceMult_dsp_dden

function get_dense_blockD(a_data, a_scalapack_info_obj, use_a_ij_transpose, i, j, b_l, block_ni, block_nj, a_block)
  real(dp), intent(in) :: a_data(:,:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: a_scalapack_info_obj
  logical, intent(in) :: use_a_ij_transpose
  integer, intent(in) :: i, j
  type(RS_SparseMatrixL), intent(in) :: b_l
  integer, intent(in) :: block_ni, block_nj
  real(dp), intent(out) :: a_block(:,:)
  logical :: get_dense_blockD

  integer io, jo, g_i, g_j, l_i, l_j

  if (a_scalapack_info_obj%active) then
    get_dense_blockD = .false.
    do io = 0, block_ni-1
    do jo = 0, block_nj-1
      g_i = b_l%dense_row_of_row(i) + io
      g_j = b_l%dense_row_of_row(j) + jo
      if (use_a_ij_transpose) then
	call coords_global_to_local(a_scalapack_info_obj, g_j, g_i, l_i, l_j)
	if (l_i > 0 .and. l_j > 0) then
	  a_block(io+1,jo+1) = a_data(l_i,l_j)
	  get_dense_blockD = .true.
	else
	  a_block(io+1,jo+1) = 0.0_dp
	end if
      else
	call coords_global_to_local(a_scalapack_info_obj, g_i, g_j, l_i, l_j)
	if (l_i > 0 .and. l_j > 0) then
	  a_block(io+1,jo+1) = a_data(l_i,l_j)
	  get_dense_blockD = .true.
	else
	  a_block(io+1,jo+1) = 0.0_dp
	end if
      end if
    end do
    end do
  else
    get_dense_blockD = .true.
    if (use_a_ij_transpose) then
      a_block(1:block_ni,1:block_nj) = transpose(a_data(b_l%dense_row_of_row(j):b_l%dense_row_of_row(j)+block_nj-1, &
				                        b_l%dense_row_of_row(i):b_l%dense_row_of_row(i)+block_ni-1))
    else
      a_block(1:block_ni,1:block_nj) = a_data(b_l%dense_row_of_row(i):b_l%dense_row_of_row(i)+block_ni-1, &
		                              b_l%dense_row_of_row(j):b_l%dense_row_of_row(j)+block_nj-1)
    end if
  end if

end function get_dense_blockD

function get_dense_blockZ(a_data, a_scalapack_info_obj, use_a_ij_transpose, i, j, b_l, block_ni, block_nj, a_block)
  complex(dp), intent(in) :: a_data(:,:)
  type(Matrix_ScaLAPACK_Info), intent(in) :: a_scalapack_info_obj
  logical, intent(in) :: use_a_ij_transpose
  integer, intent(in) :: i, j
  type(RS_SparseMatrixL), intent(in) :: b_l
  integer, intent(in) :: block_ni, block_nj
  complex(dp), intent(out) :: a_block(:,:)
  logical :: get_dense_blockZ

  integer io, jo, g_i, g_j, l_i, l_j

  if (a_scalapack_info_obj%active) then
    get_dense_blockZ = .false.
    do io = 0, block_ni-1
    do jo = 0, block_nj-1
      g_i = b_l%dense_row_of_row(i) + io
      g_j = b_l%dense_row_of_row(j) + jo
      if (use_a_ij_transpose) then
	call coords_global_to_local(a_scalapack_info_obj, g_j, g_i, l_i, l_j)
	if (l_i > 0 .and. l_j > 0) then
	  a_block(io+1,jo+1) = a_data(l_i,l_j)
	  get_dense_blockZ = .true.
	else
	  a_block(io+1,jo+1) = 0.0_dp
	end if
      else
	call coords_global_to_local(a_scalapack_info_obj, g_i, g_j, l_i, l_j)
	if (l_i > 0 .and. l_j > 0) then
	  a_block(io+1,jo+1) = a_data(l_i,l_j)
	  get_dense_blockZ = .true.
	else
	  a_block(io+1,jo+1) = 0.0_dp
	end if
      end if
    end do
    end do
  else
    get_dense_blockZ = .true.
    if (use_a_ij_transpose) then
      a_block(1:block_ni,1:block_nj) = transpose(a_data(b_l%dense_row_of_row(j):b_l%dense_row_of_row(j)+block_nj-1, &
				                        b_l%dense_row_of_row(i):b_l%dense_row_of_row(i)+block_ni-1))
    else
      a_block(1:block_ni,1:block_nj) = a_data(b_l%dense_row_of_row(i):b_l%dense_row_of_row(i)+block_ni-1, &
		                              b_l%dense_row_of_row(j):b_l%dense_row_of_row(j)+block_nj-1)
    end if
  end if

end function get_dense_blockZ

function RS_SparseMatrix_TraceMult_dden_dsp(a, b, w, a_T, b_T, diag_mask, offdiag_mask)
  type(MatrixD), intent(in) :: a
  type(RS_SparseMatrixD), intent(in) ::  b
  real(dp), intent(in), optional :: w(:)
  logical, intent(in), optional :: a_T, b_T
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)
  real(dp) :: RS_SparseMatrix_TraceMult_dden_dsp

  if (present(w)) then
    RS_SparseMatrix_TraceMult_dden_dsp = sum(partial_TraceMult(a, b, a_T, b_T, diag_mask, offdiag_mask)*w)
  else
    RS_SparseMatrix_TraceMult_dden_dsp = sum(partial_TraceMult(a, b, a_T, b_T, diag_mask, offdiag_mask))
  end if

end function RS_SparseMatrix_TraceMult_dden_dsp

function RS_SparseMatrix_TraceMult_dden_zsp(a, b, w, a_T, b_H, diag_mask, offdiag_mask)
  type(MatrixD), intent(in) :: a
  type(RS_SparseMatrixZ), intent(in) ::  b
  real(dp), intent(in), optional :: w(:)
  logical, intent(in), optional :: a_T, b_H
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: RS_SparseMatrix_TraceMult_dden_zsp

  if (present(w)) then
    RS_SparseMatrix_TraceMult_dden_zsp = sum(partial_TraceMult(a, b, a_T, b_H, diag_mask, offdiag_mask)*w)
  else
    RS_SparseMatrix_TraceMult_dden_zsp = sum(partial_TraceMult(a, b, a_T, b_H, diag_mask, offdiag_mask))
  end if

end function RS_SparseMatrix_TraceMult_dden_zsp

function RS_SparseMatrix_TraceMult_zden_dsp(a, b, w, a_H, b_T, diag_mask, offdiag_mask)
  type(MatrixZ), intent(in) :: a
  type(RS_SparseMatrixD), intent(in) ::  b
  real(dp), intent(in), optional :: w(:)
  logical, intent(in), optional :: a_H, b_T
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: RS_SparseMatrix_TraceMult_zden_dsp

  if (present(w)) then
    RS_SparseMatrix_TraceMult_zden_dsp = sum(partial_TraceMult(a, b, a_H, b_T, diag_mask, offdiag_mask)*w)
  else
    RS_SparseMatrix_TraceMult_zden_dsp = sum(partial_TraceMult(a, b, a_H, b_T, diag_mask, offdiag_mask))
  end if

end function RS_SparseMatrix_TraceMult_zden_dsp

function RS_SparseMatrix_TraceMult_dsp_zden(a, b, w, a_T, b_H, diag_mask, offdiag_mask)
  type(RS_SparseMatrixD), intent(in) ::  a
  type(MatrixZ), intent(in) :: b
  real(dp), intent(in), optional :: w(:)
  logical, intent(in), optional :: a_T, b_H
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: RS_SparseMatrix_TraceMult_dsp_zden

  if (present(w)) then
    RS_SparseMatrix_TraceMult_dsp_zden = sum(partial_TraceMult(a, b, a_T, b_H, diag_mask, offdiag_mask)*w)
  else
    RS_SparseMatrix_TraceMult_dsp_zden = sum(partial_TraceMult(a, b, a_T, b_H, diag_mask, offdiag_mask))
  end if

end function RS_SparseMatrix_TraceMult_dsp_zden

function RS_SparseMatrix_TraceMult_dsp_dden(a, b, w, a_T, b_T, diag_mask, offdiag_mask)
  type(RS_SparseMatrixD), intent(in) ::  a
  type(MatrixD), intent(in) :: b
  real(dp), intent(in), optional :: w(:)
  logical, intent(in), optional :: a_T, b_T
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: RS_SparseMatrix_TraceMult_dsp_dden

  if (present(w)) then
    RS_SparseMatrix_TraceMult_dsp_dden = sum(partial_TraceMult(a, b, a_T, b_T, diag_mask, offdiag_mask)*w)
  else
    RS_SparseMatrix_TraceMult_dsp_dden = sum(partial_TraceMult(a, b, a_T, b_T, diag_mask, offdiag_mask))
  end if

end function RS_SparseMatrix_TraceMult_dsp_dden

function RS_SparseMatrix_TraceMult_zsp_zden(a, b, w, a_H, b_H, diag_mask, offdiag_mask)
  type(RS_SparseMatrixZ), intent(in) ::  a
  type(MatrixZ), intent(in) :: b
  real(dp), intent(in), optional :: w(:)
  logical, intent(in), optional :: a_H, b_H
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: RS_SparseMatrix_TraceMult_zsp_zden

  if (present(w)) then
    RS_SparseMatrix_TraceMult_zsp_zden = sum(partial_TraceMult(a, b, a_H, b_H, diag_mask, offdiag_mask)*w)
  else
    RS_SparseMatrix_TraceMult_zsp_zden = sum(partial_TraceMult(a, b, a_H, b_H, diag_mask, offdiag_mask))
  end if

end function RS_SparseMatrix_TraceMult_zsp_zden

function RS_SparseMatrix_TraceMult_zden_zsp(a, b, w, a_H, b_H, diag_mask, offdiag_mask)
  type(MatrixZ), intent(in) :: a
  type(RS_SparseMatrixZ), intent(in) ::  b
  real(dp), intent(in), optional :: w(:)
  logical, intent(in), optional :: a_H, b_H
  logical, intent(in), optional :: diag_mask(:), offdiag_mask(:)
  complex(dp) :: RS_SparseMatrix_TraceMult_zden_zsp

  if (present(w)) then
    RS_SparseMatrix_TraceMult_zden_zsp = sum(partial_TraceMult(a, b, a_H, b_H, diag_mask, offdiag_mask)*w)
  else
    RS_SparseMatrix_TraceMult_zden_zsp = sum(partial_TraceMult(a, b, a_H, b_H, diag_mask, offdiag_mask))
  end if

end function RS_SparseMatrix_TraceMult_zden_zsp

subroutine RS_SparseMatrixD_multDiagRL_d(this, A, diag)
  type(RS_SparseMatrixD), intent(inout) :: this
  type(RS_SparseMatrixD), intent(in) :: A
  real(dp), intent(in) :: diag(:)

  integer i, jj, j, block_ni, block_nj, io, jo
  real(dp) fac

  do i=1, A%l%N
    block_ni = A%l%block_size(i)
    do jj=A%l%row_indices(i), A%l%row_indices(i+1)-1
      j = A%l%col(jj)
      block_nj = A%l%block_size(j)
      do jo=0, block_nj-1
      do io=0, block_ni-1
	fac = 0.5_dp*(diag(this%l%dense_row_of_row(i)+io)+diag(this%l%dense_row_of_row(j)+jo))
	this%data(this%l%data_ptrs(jj)+io+block_ni*jo) = fac*A%data(this%l%data_ptrs(jj)+io+block_ni*jo)
      end do
      end do
    end do
  end do

end subroutine RS_SparseMatrixD_multDiagRL_d

subroutine RS_SparseMatrixZ_multDiagRL_d(this, A, diag)
  type(RS_SparseMatrixZ), intent(inout) :: this
  type(RS_SparseMatrixZ), intent(in) :: A
  real(dp), intent(in) :: diag(:)

  integer i, jj, j, block_ni, block_nj, io, jo
  real(dp) fac

  do i=1, A%l%N
    block_ni = A%l%block_size(i)
    do jj=A%l%row_indices(i), A%l%row_indices(i+1)-1
      j = A%l%col(jj)
      block_nj = A%l%block_size(j)
      do jo=0, block_nj-1
      do io=0, block_ni-1
	fac = 0.5_dp*(diag(this%l%dense_row_of_row(i)+io)+diag(this%l%dense_row_of_row(j)+jo))
	this%data(this%l%data_ptrs(jj)+io+block_ni*jo) = fac*A%data(this%l%data_ptrs(jj)+io+block_ni*jo)
      end do
      end do
    end do
  end do

end subroutine RS_SparseMatrixZ_multDiagRL_d

subroutine check_sparse_layout(l)
  type(RS_SparseMatrixL), intent(in) :: l
  integer i, ji, j

  do i=1, l%N
    do ji=l%row_indices(i), l%row_indices(i+1)-1
      j = l%col(ji)

      call print("check_sparse " // i // " " // j)

      if (ji > l%row_indices(i)) then
	if (j <= l%col(ji-1)) then
	  call print("ERROR row " // i // " col " // j // " (ji="//ji//") out of order")
	  call print(l%col(l%row_indices(i):l%row_indices(i+1)-1))
	end if
      end if

      if (.not. any(l%col(l%row_indices(j):l%row_indices(j+1)-1) == i)) then
	call print ("ERROR row "//i//" col "//j//" has no transpose match")
	call print(l%col(l%row_indices(j):l%row_indices(j+1)-1))
      end if

    end do
  end do
end subroutine

end module RS_SparseMatrix_module
