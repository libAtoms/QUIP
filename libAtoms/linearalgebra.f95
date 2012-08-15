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
!X  Linear algebra module
!X  
!%  This is a general purpose linear algebra module, extending
!%  the Fortran intrinsics
!%  This module defines dot products between matrices and vectors, 
!%  wrappers to \textsc{lapack} for matrix diagonalisation and inversion,
!%  as well as array searching, sorting and averaging.
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module linearalgebra_module
  use error_module
  use system_module
  use units_module
  implicit none
  SAVE

  logical :: use_intrinsic_blas = .false. 
  !% If set to true, use internal routines instead of \textsc{blas} calls for matrix
  !% multiplication. Can be changed at runtime. The default is true.
  integer, parameter :: NOT_FACTORISED = 0
  integer, parameter :: CHOLESKY       = 1
  integer, parameter :: QR             = 2

  type LA_Matrix
     real(qp), dimension(:,:), allocatable :: matrix, factor
     real(qp), dimension(:), allocatable :: s, tau
     integer :: n, m
     logical :: initialised = .false.
     logical :: equilibrated = .false.
     integer :: factorised = NOT_FACTORISED
  endtype LA_Matrix

  interface Initialise
    module procedure LA_Matrix_Initialise
  endinterface Initialise

  interface assignment(=)
     module procedure LA_Matrix_Initialise
  end interface assignment(=)

  interface Finalise
    module procedure LA_Matrix_Finalise
  endinterface Finalise

  interface Matrix_Solve
    module procedure LA_Matrix_Solve_Vector, LA_Matrix_Solve_Matrix
  endinterface Matrix_Solve

  interface Matrix_QR_Solve
    module procedure LA_Matrix_QR_Solve_Vector, LA_Matrix_QR_Solve_Matrix
  endinterface Matrix_QR_Solve

  interface find
     module procedure find_indices
  end interface

  interface sign
     module procedure int_sign, real_sign
  end interface

  private :: k_means_clustering_pick_1d
  interface k_means_clustering_pick
     module procedure k_means_clustering_pick_1d
  end interface

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X Interfaces to mathematical operations
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% Construct a diagonal matrix from a vector, or extract the diagonal elements
  !% of a matrix and return them as a vector.
  ! matrix = diag(vector) , vector = diag(matrix)
  private :: matrix_diagonal_r, matrix_diagonal_c, vector_as_diag_matrix_r, vector_as_diag_matrix_c
  interface diag
    module procedure matrix_diagonal_r, matrix_diagonal_c, vector_as_diag_matrix_r, vector_as_diag_matrix_c
  end interface

  !% Heavily overloaded interface to print matrices and vectors of reals, complex numbers or integers.
  private :: matrix_print_mainlog, matrix_print, vector_print_mainlog, vector_print
  private :: matrix_z_print_mainlog, matrix_z_print, vector_z_print_mainlog, vector_z_print
  private :: int_matrix_print_mainlog, int_matrix_print
  private :: integer_array_print, integer_array_print_mainlog, logical_array_print, logical_array_print_mainlog
  interface print
     module procedure matrix_print_mainlog, matrix_print, vector_print_mainlog, vector_print
     module procedure matrix_z_print_mainlog, matrix_z_print, vector_z_print_mainlog, vector_z_print
     module procedure integer_array_print, integer_array_print_mainlog, logical_array_print_mainlog
     module procedure int_matrix_print_mainlog, int_matrix_print, logical_array_print
  end interface

  private :: matrix_z_print_mathematica, matrix_print_mathematica
  interface print_mathematica
     module procedure matrix_z_print_mathematica, matrix_print_mathematica
  end interface print_mathematica

  ! matrix*matrix , matrix*vector
  !% Overloaded multiplication for matrix $\times$ matrix and matrix $\times$ vector.
  !% Use as 'C = A .mult. B' for matrices and 'C = A .mult. v' for a matrix times a vector.
  private ::  matrix_product_vect, matrix_product_ddd, matrix_product_int_vect, matrix_product_zzz
  private :: matrix_product_int_mat
  interface operator(.mult.)
     module procedure  matrix_product_vect, matrix_product_ddd, matrix_product_int_vect, matrix_product_zzz
     module procedure matrix_product_int_mat
  end interface
  !% Interface for routines that want function to be called 'matrix_product()'
  interface matrix_product
    module procedure matrix_product_ddd
  end interface

  !% Overloaded multiplication for matrix $\times$ matrix in subroutine form,
  !% with no return value allocated on the stack
  interface matrix_product_sub
    module procedure matrix_product_sub_ddd, matrix_product_sub_zzz, matrix_vector_product_sub_ddd
  end interface

  ! matrix*diag(vector)
  !% Matrix product with the diagonal matrix constructed from a vector. For example,
  !% 'B = A .multd. d' is equivalent to 'B = A .mult. diag(d)'
  private :: matrix_product_vect_asdiagonal_dd, matrix_product_vect_asdiagonal_zz
  interface operator(.multd.)
     module procedure matrix_product_vect_asdiagonal_dd, matrix_product_vect_asdiagonal_zz
  end interface

  !% Matrix product with the diagonal matrix constructed from a vector in subroutine form,
  !% with no return value allocated on the stack
#ifdef HAVE_QP
  private :: matrix_product_vect_asdiagonal_sub_qqq
#endif
  private :: matrix_product_vect_asdiagonal_sub_ddd
  private :: matrix_product_vect_asdiagonal_sub_zzd
  private :: matrix_product_vect_asdiagonal_sub_zdz
  private :: matrix_product_vect_asdiagonal_sub_zzz
  private :: vect_asdiagonal_product_matrix_sub_ddd
  private :: vect_asdiagonal_product_matrix_sub_zzd
  private :: vect_asdiagonal_product_matrix_sub_zdz
  private :: vect_asdiagonal_product_matrix_sub_zzz
  interface matrix_product_vect_asdiagonal_sub
#ifdef HAVE_QP  
    module procedure matrix_product_vect_asdiagonal_sub_qqq
#endif
    module procedure matrix_product_vect_asdiagonal_sub_ddd
    module procedure matrix_product_vect_asdiagonal_sub_zzd
    module procedure matrix_product_vect_asdiagonal_sub_zdz
    module procedure matrix_product_vect_asdiagonal_sub_zzz
    module procedure vect_asdiagonal_product_matrix_sub_ddd
    module procedure vect_asdiagonal_product_matrix_sub_zzd
    module procedure vect_asdiagonal_product_matrix_sub_zdz
    module procedure vect_asdiagonal_product_matrix_sub_zzz
  end interface matrix_product_vect_asdiagonal_sub

  private :: matrix_product_vect_asdiagonal_RL_sub_ddd
  private :: matrix_product_vect_asdiagonal_RL_sub_zzd
  interface matrix_product_vect_asdiagonal_RL_sub
    module procedure matrix_product_vect_asdiagonal_RL_sub_ddd
    module procedure matrix_product_vect_asdiagonal_RL_sub_zzd
  end interface matrix_product_vect_asdiagonal_RL_sub

  ! matrix*diag(vector)*matrix.t
  !% Matrix product of matrix and a vector in the form:
  !%> matrix * diag(vector) * transpose(matrix)
  private :: matrix_cfct
  interface matrix_mvmt
     module procedure matrix_cfct
  end interface
  
  ! vector.vector , matrix.matrix
  !% Inner product of two vectors or two matrices. For two vectors $\mathbf{v}$ and 
  !% $\mathbf{w}$ this is simply:
  !%\begin{displaymath}
  !% d = \sum_{i=1}^N v_i w_i
  !%\end{displaymath}
  !% For $N \times M$ matrices $A$ and $B$ it is defined similarily:
  !%\begin{displaymath}
  !% d = \sum_{i=1}^N \sum_{j=1}^M A_{ij} B_{ij}
  !%\end{displaymath}
  private :: matrix_dotproduct_matrix,vector_dotpr
  interface operator(.dot.)
     module procedure matrix_dotproduct_matrix,vector_dotpr
  end interface

  !% Floating point equality testing. Returns false if
  !% '(abs(x-y) > NUMERICAL_ZERO * abs(x))', and true otherwise.
  private :: real_feq,complex_feq,matrix_feq,vector_feq  
  interface operator(.feq.)
     module procedure real_feq,complex_feq,matrix_feq,vector_feq
  end interface

  !% Floating point inequality testing.
  private :: real_fne,complex_fne,matrix_fne,vector_fne
  interface operator(.fne.)
     module procedure real_fne,complex_fne,matrix_fne,vector_fne
  end interface

  !% Overloaded interfaces to \textsc{lapack} matrix diagonlisation
  !% functions for real and complex matrices. Always calls
  !% \textsc{lapack}, regardless of the 'use_intrinsic_blas' setting.
  !% Both diagonalisation and solution of the generalised eigenproblem
  !% are supported, but only for real positive definite symmetric or
  !% complex hermitian postive definite matrices.
  private :: matrix_diagonalise, matrix_diagonalise_generalised
  private :: matrix_z_diagonalise, matrix_z_diagonalise_generalised
  interface diagonalise
     module procedure matrix_diagonalise, matrix_diagonalise_generalised
     module procedure matrix_z_diagonalise, matrix_z_diagonalise_generalised
  end interface

  interface nonsymmetric_diagonalise
     module procedure matrix_nonsymmetric_diagonalise
  end interface

  !% Calculate the inverse of a matrix in-place. Uses \textsc{lapack} to compute the inverse.
  interface inverse
    module procedure matrix_inverse, matrix_z_inverse
  end interface inverse

  !% Return the outer product of two vectors. Usage is 'x .outer. y'.
  private :: outer, z_outer_zz
  interface operator(.outer.)
     module procedure outer, z_outer_zz
#ifdef HAVE_QP  
    module procedure outer_qq
#endif
  end interface

  !% Interface to return the real outer product. Usage is 'x .realouter. y'.
  private :: d_outer_zz
  interface operator(.realouter.)
     module procedure outer, d_outer_zz
  end interface

  !% Cross product between two 3-vectors. Usage is 'x .cross. y'.
  private :: cross_product
  interface operator(.cross.)
     module procedure cross_product
  end interface

  !% Test for matrix symmetry (with floating point equality test '.feq.' as described above).
  private :: matrix_is_symmetric, matrix_z_is_symmetric, int_matrix_is_symmetric
  interface is_symmetric
     module procedure matrix_is_symmetric, matrix_z_is_symmetric, int_matrix_is_symmetric
  end interface is_symmetric

  !% Test for matrix hermiticity (with floating point equals test).
  private :: matrix_z_is_hermitian
  interface is_hermitian
     module procedure matrix_z_is_hermitian
  end interface

  !% Test if matrix is square 
  private :: matrix_square, matrix_z_square, int_matrix_square
  interface is_square
     module procedure matrix_square, matrix_z_square, int_matrix_square
  end interface

  !% Test is matrix is unitary i.e. if $M M^{-1} = I$
  private :: matrix_is_orthogonal
  interface is_orthogonal
     module procedure matrix_is_orthogonal
  end interface is_orthogonal

  !% Symmetrise a matrix: $$A \to \frac{A + A^T}{2}$$
  private :: matrix_symmetrise
  interface symmetrise
     module procedure matrix_symmetrise
  end interface

  !% Return the trace of a matrix.
  private :: matrix_trace
  interface trace
     module procedure matrix_trace
#ifdef HAVE_QP  
    module procedure matrix_trace_q
#endif
  end interface

  private :: matrix_trace_mult
  interface trace_mult
     module procedure matrix_trace_mult
  end interface

  !% Adds the identity to a matrix
  private :: matrix_add_identity_r, matrix_add_identity_c
  interface add_identity
    module procedure matrix_add_identity_r, matrix_add_identity_c
  end interface add_identity

  !% Adds $x\mathbf{I}$ to a matrix
  private :: matrix_add_xidentity_r, matrix_add_xidentity_c
  interface add_xidentity
    module procedure matrix_add_xidentity_r, matrix_add_xidentity_c
  end interface add_xidentity

  !% Return the euclidean norm of a vector or of an array.
  !% For a single vector 'x', 'norm(x)' is equal to 'sqrt(x .dot. x)'
  !% A two-dimensional array is treated as a list of vectors in either
  !% Fortran ('dir=1') or C ('dir=2') style-ordering.
  !% The result is then a one-dimensional array of the norms of each vector. 
  private :: vector_norm, array_norm
  interface norm
     module procedure vector_norm, array_norm
  end interface

  !% Euclidean norm$^2$ of a vector or a of a list of vectors. Result is equal
  !% to 'x .dot. x' for a single vector 'x'.
  private :: vector_normsq, array_normsq
  interface normsq
     module procedure vector_normsq, array_normsq
#ifdef HAVE_QP  
    module procedure vector_normsq_q
#endif
  end interface

  !% Randomise the elements of an array. Uniformly distributed random quantities in the range 
  !% $(-\frac{a}{2},\frac{a}{2})$ are added to each element of the vector or matrix.
  private :: vector_randomise, matrix_randomise, matrix_randomise_vweight, matrix_z_randomise, vector_z_randomise
  interface randomise
     module procedure vector_randomise, matrix_randomise, matrix_randomise_vweight, matrix_z_randomise, vector_z_randomise
  end interface

  !% Search an array by element or by row.
  private :: find_in_array_element_i, find_in_array_element_s, find_in_array_row
  interface find_in_array
     module procedure find_in_array_element_i, find_in_array_element_s, find_in_array_row
  end interface

  !% Root-mean-square difference calculation for components of two vectors or arrays.
  private :: rms_diff1, rms_diff2
  interface rms_diff
     module procedure rms_diff1, rms_diff2
  end interface rms_diff

  
  !% Returns a vector, contining a histogram of frequencies.
  private :: vector_histogram
  interface histogram
     module procedure vector_histogram
  end interface histogram

  private :: sort_array_i, sort_array_r
  interface sort_array
    module procedure sort_array_i, sort_array_r
  end interface sort_array

  private :: heap_sort_i, heap_sort_r
  interface heap_sort
    module procedure heap_sort_i, heap_sort_r
  end interface heap_sort

  private :: insertion_sort_i, insertion_sort_r
  interface insertion_sort
    module procedure insertion_sort_i, insertion_sort_r
  end interface insertion_sort

  private :: binary_search_i, binary_search_r
  interface binary_search
     module procedure binary_search_i, binary_search_r
  endinterface

  ! in addition, these are the intrinsics defined on matrices and vectors:
  !
  ! =    : intrinsic assignment, the left part needs to be allocated
  ! *    : intrinsic product                  : elementwise product
  ! +    : intrinsic sum                      : elementwise sum
  ! -    : intrinsic difference               : elementwise difference
  ! /    : intrinsic division                 : elementwise division
  !
  ! cos,sin,exp,log,sqrt,**,abs .......       : elementwise
  ! matmul : intrinsic row*column product
  ! transpose : intrinsic transposition
  !
  ! in the following you can specify the direction of matrix, as optional argument,
  !    which the operation will be performed along.
  !
  ! maxval: return maximum in the matrix 
  ! minval: return minimum in the matrix
  ! sum  : return sum of matrix or vector elements  
  ! product  : return product of matrix or vector elements
  ! reshape(matrix,newshape) : to change shape
  !       (the order of elements in memory will not change)
  !

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X Interface to array size checking routines
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  ! Usage: call check_size(name,array,sizes,calling routine)
  ! e.g. check_size('Velocity',velo,(/3,ds%N/),'DS_Add_Atom')

  !% Overloaded interface to assert that the size of an array is correct. 
  !% If the test fails 'system_abort' is called with an appropriate error message
  !% constructed from the 'arrayname' and 'caller' arguments.
  private :: check_size_int_dim1, check_size_int_dim1_s, check_size_int_dim2  
  private :: check_size_real_dim1,check_size_real_dim1_s, check_size_real_dim2
#ifdef HAVE_QP
  private :: check_size_quad_dim1,check_size_quad_dim1_s, check_size_quad_dim2
#endif
  private :: check_size_complex_dim1,check_size_complex_dim1_s, check_size_complex_dim2
  private :: check_size_log_dim1, check_size_log_dim1_s, check_size_log_dim2  
  interface check_size
     module procedure check_size_int_dim1, check_size_int_dim1_s, check_size_int_dim2  
     module procedure check_size_real_dim1, check_size_real_dim1_s, check_size_real_dim2
#ifdef HAVE_QP
     module procedure check_size_quad_dim1, check_size_quad_dim1_s, check_size_quad_dim2
#endif
     module procedure check_size_complex_dim1, check_size_complex_dim1_s, check_size_complex_dim2
     module procedure check_size_log_dim1, check_size_log_dim1_s, check_size_log_dim2  
  end interface check_size

  !% Update a measure of a recent average by decaying its current value and adding on a new sample
  private :: update_exponential_average_s, update_exponential_average_v, update_exponential_average_d2
  interface update_exponential_average
     module procedure update_exponential_average_s, update_exponential_average_v, update_exponential_average_d2
  end interface update_exponential_average

  private :: matrix_exp_d
  interface matrix_exp
     module procedure matrix_exp_d
  end interface matrix_exp

CONTAINS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X
!X random crap
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function int_sign(n) 
    integer::n, int_sign
    int_sign = sign(1, n)
  end function int_sign

  function real_sign(n) 
    real(dp)::n, real_sign
    real_sign = sign(1.0_dp, n)
  end function real_sign

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X
!X MATRIX stuff
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% Kronecker delta function $\delta_{ij}$.
  function delta(i,j)
    integer, intent(in) :: i,j
    integer             :: delta
    if (i == j) then
       delta = 1
    else
       delta = 0
    end if
  end function delta

  ! diag(v)
  !
  ! return a matrix with the vector as its diagonal 
  function   vector_as_diag_matrix_r(vect) result(matrix)
    real(dp),intent(in), dimension(:) :: vect
    real(dp), dimension(size(vect),size(vect)) :: matrix
    integer::i

    matrix = 0.0_dp

    do i=1,size(vect)
       matrix(i,i)=vect(i)
    end do

  end function vector_as_diag_matrix_r

  ! diag(v)
  !
  ! return a matrix with the vector as its diagonal 
  function   vector_as_diag_matrix_c(vect) result(matrix)
    complex(dp),intent(in), dimension(:) :: vect
    complex(dp), dimension(size(vect),size(vect)) :: matrix
    integer::i

    matrix = 0.0_dp

    do i=1,size(vect)
       matrix(i,i)=vect(i)
    end do

  end function vector_as_diag_matrix_c

  ! diag(matrix)
  !
  ! returns the diagonal of a matrix as a vector
  function matrix_diagonal_r(matrix) result (vect)
    real(dp),intent(in), dimension(:,:) ::matrix 
    real(dp), dimension(size(matrix,1)) ::vect
    integer::I
    
    if (.NOT.is_square(matrix)) call system_abort('Matrix_diagonal: matrix not squared')

    do i=1,size(matrix,1)
       vect(i)=matrix(i,i)
    end do

  end function matrix_diagonal_r

  ! diag(matrix)
  !
  ! returns the diagonal of a matrix as a vector
  function matrix_diagonal_c(matrix) result (vect)
    complex(dp),intent(in), dimension(:,:) ::matrix 
    complex(dp), dimension(size(matrix,1)) ::vect
    integer::I
    
    if (.NOT.is_square(matrix)) call system_abort('Matrix_diagonal: matrix not squared')

    do i=1,size(matrix,1)
       vect(i)=matrix(i,i)
    end do

  end function matrix_diagonal_c


  ! m(:,:) .mult. v(:)
  !
  ! matrix times a vector
  function matrix_product_vect(matrix,vect) result (prodvect)
    real(dp),intent(in), dimension(:)   :: vect
    real(dp),intent(in), dimension(:,:) :: matrix
    real(dp), dimension(size(matrix,1)) ::prodvect
    integer::N,M,i,j

    N=size(matrix,1)
    M=size(matrix,2)

    if (M /= size(vect)) call check_size('Vector',vect,M,'Matrix_Product_Vect')
    prodvect = 0.0_dp
    do j=1,M
       forall (i=1:N) prodvect(i) = prodvect(i) + matrix(i,j)*vect(j)
    end do

  end function matrix_product_vect

  ! m(:,:) .mult. a(:)
  !
  ! product of real matrix and integer vector
  function matrix_product_int_vect(matrix,intvect) result(prodvect)

     real(dp), intent(in), dimension(:,:) :: matrix
     integer,  intent(in), dimension(:)   :: intvect
     real(dp)                             :: prodvect(size(matrix,1))
     integer                              :: M,N,i,j

     N=size(matrix,1)
     M=size(matrix,2)
     if(M /= size(intvect)) &
          call check_size('Integer Vector',intvect,M,'Matrix_product_int_vect')
     prodvect = 0.0_dp
     do j=1,M
        forall(i=1:N) prodvect(i) = prodvect(i)+matrix(i,j)*intvect(j)
     end do

   end function matrix_product_int_vect

  function matrix_product_int_mat(matrix,intmat) result(prodmat)

     real(dp), intent(in), dimension(:,:) :: matrix
     integer,  intent(in), dimension(:,:)   :: intmat
     real(dp)                             :: prodmat(size(matrix,1),size(intmat,2))
     integer                              :: M,N,L,i

     N=size(matrix,1)
     M=size(matrix,2)
     L=size(intmat,2)
     if (M /= size(intmat,1)) &
      call check_size('Integer Matrix',intmat,(/M,L/),'Matrix_product_int_mat')
     prodmat = 0.0_dp
     do i=1,L
	prodmat(:,i) = matrix_product_int_vect(matrix,intmat(:,i))
     end do

   end function matrix_product_int_mat

   ! subroutine form of m(:,:) .multd. v(:) for real = real * real
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
#ifdef HAVE_QP
   subroutine matrix_product_vect_asdiagonal_sub_qqq(lhs, matrix, vect) 
    real(qp), dimension(:,:), intent(out) :: lhs
    real(qp), dimension(:,:), intent(in) :: matrix
    real(qp), dimension(:), intent(in) :: vect

    integer  :: i
     
    do i = 1, size(vect)
       lhs(:,i) = vect(i) * matrix(:,i)
    enddo
     
   endsubroutine matrix_product_vect_asdiagonal_sub_qqq
#endif

   subroutine matrix_product_vect_asdiagonal_sub_ddd(lhs, matrix, vect) 
    real(dp), dimension(:,:), intent(out) :: lhs
    real(dp), dimension(:,:), intent(in) :: matrix
    real(dp), dimension(:), intent(in) :: vect

    integer  :: i
     
    do i = 1, size(vect)
       lhs(:,i) = vect(i) * matrix(:,i)
    enddo
     
   endsubroutine matrix_product_vect_asdiagonal_sub_ddd

   ! subroutine form of m(:,:) .multd. v(:) for complex = real * complex
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine matrix_product_vect_asdiagonal_sub_zdz(lhs, matrix, vect) 
    complex(dp), dimension(:,:), intent(out) :: lhs
    real(dp), dimension(:,:), intent(in) :: matrix
    complex(dp), dimension(:), intent(in) :: vect

    integer  :: i
     
    do i = 1, size(vect)
       lhs(:,i) = vect(i) * matrix(:,i)
    enddo
     
   endsubroutine matrix_product_vect_asdiagonal_sub_zdz

   ! subroutine form of m(:,:) .multd. v(:) for complex = complex * real
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine matrix_product_vect_asdiagonal_sub_zzd(lhs, matrix, vect) 
    complex(dp), dimension(:,:), intent(out) :: lhs
    complex(dp), dimension(:,:), intent(in) :: matrix
    real(dp), dimension(:), intent(in) :: vect

    integer :: i
     
    do i = 1, size(vect)
       lhs(:,i) = vect(i) * matrix(:,i)
    enddo
     
   endsubroutine matrix_product_vect_asdiagonal_sub_zzd

   ! subroutine form of m(:,:) .multd. v(:) for complex = complex * complex
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine matrix_product_vect_asdiagonal_sub_zzz(lhs, matrix, vect) 
    complex(dp), dimension(:,:), intent(out) :: lhs
    complex(dp), dimension(:,:), intent(in) :: matrix
    complex(dp), dimension(:), intent(in) :: vect

    integer :: i
     
    do i = 1, size(vect)
       lhs(:,i)=vect(i)*matrix(:,i)
    end do
     
   endsubroutine matrix_product_vect_asdiagonal_sub_zzz

   ! subroutine form of v(:) .multd. m(:,:) for real = real * real
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine vect_asdiagonal_product_matrix_sub_ddd(lhs, vectL, matrix) 
     real(dp), dimension(:,:), intent(out) :: lhs
     real(dp), dimension(:), intent(in) :: vectL
     real(dp), dimension(:,:), intent(in) :: matrix

     integer :: i
     
     do i = 1, size(matrix,2)
       lhs(:,i)=vectL(:)*matrix(:,i)
     end do
     
   end subroutine vect_asdiagonal_product_matrix_sub_ddd

   ! subroutine form of v(:) .multd. m(:,:) for complex = complex * real
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine vect_asdiagonal_product_matrix_sub_zzd(lhs, vectL, matrix) 
     complex(dp), dimension(:,:), intent(out) :: lhs
     complex(dp), dimension(:), intent(in) :: vectL
     real(dp), dimension(:,:), intent(in) :: matrix

     integer :: i
     
     do i = 1, size(matrix,2)
       lhs(:,i)=vectL(:)*matrix(:,i)
     end do
     
   end subroutine vect_asdiagonal_product_matrix_sub_zzd

   ! subroutine form of v(:) .multd. m(:,:) for complex = real * complex
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine vect_asdiagonal_product_matrix_sub_zdz(lhs, vectL, matrix) 
     complex(dp), dimension(:,:), intent(out) :: lhs
     real(dp), dimension(:), intent(in) :: vectL
     complex(dp), dimension(:,:), intent(in) :: matrix

     integer :: i
     
     do i = 1, size(matrix,2)
       lhs(:,i)=vectL(:)*matrix(:,i)
     end do
     
   end subroutine vect_asdiagonal_product_matrix_sub_zdz

   ! subroutine form of v(:) .multd. m(:,:) for complex = complex * complex
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine vect_asdiagonal_product_matrix_sub_zzz(lhs, vectL, matrix) 
     complex(dp), dimension(:,:), intent(out) :: lhs
     complex(dp), dimension(:), intent(in) :: vectL
     complex(dp), dimension(:,:), intent(in) :: matrix

     integer :: i
     
     do i = 1, size(matrix,2)
       lhs(:,i)=vectL(:)*matrix(:,i)
     end do
     
   end subroutine vect_asdiagonal_product_matrix_sub_zzz

   ! m(:,:) .multd. v(:)
   !
   ! returns product of matrix and vector as diagonal of another matrix
   function matrix_product_vect_asdiagonal_dd(matrix,vect) result (prodmatrix)
     real(dp),intent(in), dimension(:) :: vect
     real(dp),dimension(size(vect),size(vect))::prodmatrix
     real(dp),intent(in),dimension(:,:)::matrix

     call matrix_product_vect_asdiagonal_sub(prodmatrix, matrix, vect)

   end function matrix_product_vect_asdiagonal_dd

   ! m(:,:) .multd. v(:)
   !
   ! returns product of matrix and vector as diagonal of another matrix
   function matrix_product_vect_asdiagonal_zz(matrix,vect) result (prodmatrix)
     complex(dp),intent(in), dimension(:) :: vect
     complex(dp),dimension(size(vect),size(vect))::prodmatrix
     complex(dp),intent(in),dimension(:,:)::matrix

     call matrix_product_vect_asdiagonal_sub(prodmatrix, matrix, vect)

   end function matrix_product_vect_asdiagonal_zz

   ! multiply a matrix by a vector, interepreted as a diagonal matrix, on 
   ! the left and right hand sides, weighted by 0.5
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine matrix_product_vect_asdiagonal_RL_sub_ddd(lhs, matrix, vect) 
    real(dp), intent(out) :: lhs(:,:)
    real(dp), intent(in) :: matrix(:,:)
    real(dp), intent(in) :: vect(:)

     real(dp)::tmp
     integer::i,j,N
     
     N=size(vect)
     
     do j=1,N
        tmp=vect(j)
        do i=1,N
           lhs(i,j)=0.5_dp*(tmp+vect(i))*matrix(i,j)
        end do
     end do
     
   end subroutine matrix_product_vect_asdiagonal_RL_sub_ddd

   ! multiply a matrix by a vector, interepreted as a diagonal matrix, on 
   ! the left and right hand sides, weighted by 0.5
   !
   ! first argument set to product of matrix and a vector as a diagonal of another matrix (no temp on stack)
   subroutine matrix_product_vect_asdiagonal_RL_sub_zzd(lhs, matrix, vect) 
    complex(dp), intent(out) :: lhs(:,:)
    complex(dp), intent(in) :: matrix(:,:)
    real(dp), intent(in) :: vect(:)

     real(dp)::tmp
     integer::i,j,N
     
     N=size(vect)
     
     do j=1,N
        tmp=vect(j)
        do i=1,N
           lhs(i,j)=0.5_dp*(tmp+vect(i))*matrix(i,j)
        end do
     end do
     
   end subroutine matrix_product_vect_asdiagonal_RL_sub_zzd

   
   ! multiply matrix C with a diagonal matrix that has vect as diagonal, 
   ! and then multiply the result and the tranpose of C. 
   ! NEEDS TO BE OPTIMIZED
   function matrix_cfct(matrix,vect) result (prodmatrix)
     real(dp),intent(in),dimension(:) :: vect
     real(dp),intent(in),dimension(:,:) :: matrix
     real(dp), dimension(size(matrix,1),size(matrix,1)) ::prodmatrix
     integer::i,j,k,N,M
     real(dp)::tmp
     
     N=size(matrix,1)
     M=size(matrix,2) 

     if(M /= size(vect)) &
          call check_size('Vector',vect,M,'Matrix_CFCT')

     do I=1,N
        do j=1,N
           tmp=0.0_dp
           do k=1,size(vect)/2
              tmp=tmp+matrix(I,K)*vect(k)*matrix(j,K)
           end do
           do k=size(vect)/2+1,size(vect)
              if (abs(vect(k)).LT.(real(1.e-10))) cycle
              
              tmp=tmp+matrix(I,K)*vect(k)*matrix(j,K)
           end do
           prodmatrix(I,J)=tmp
        end do
     end do
     
   end function matrix_cfct

   
   ! m1(:,:) .dot. m2(:,:)
   !
   ! scalar product between matrixes 
   function matrix_dotproduct_matrix(matrix1,matrix2) result(prod)
     real(dp),intent(in), dimension(:,:) :: matrix1
     real(dp),intent(in), dimension(:,:) :: matrix2
     real(dp)::prod
     integer::i,j
     
     if(size(matrix1,1) /= size(matrix2,1) .or. size(matrix1,2) /= size(matrix2,2)) & 
          call check_size('Matrix2',matrix2,shape(matrix1),'Matrix_Dotproduct_Matrix')

     prod = 0.0_dp
     do j = 1, size(matrix1,2)
        do i = 1, size(matrix1,1)
           prod = prod + matrix1(i,j)*matrix2(i,j)
        end do
     end do
     
   end function matrix_dotproduct_matrix

   ! subroutine form of m1(:,:) .mult. m2(:,:)
   !
   ! set first argument to matrix product (no temporary on the stack)
   subroutine matrix_product_sub_ddd(lhs, matrix1, matrix2, m1_transpose, m2_transpose, &
    lhs_factor, rhs_factor)
      real(dp), intent(out) :: lhs(:,:)
      real(dp), intent(in) :: matrix1(:,:), matrix2(:,:)
      logical, intent(in), optional :: m1_transpose, m2_transpose
      real(dp), intent(in), optional :: lhs_factor, rhs_factor

     integer::M,N,maxd,K
     character(len=1) :: m1_transp, m2_transp
     integer :: m1_r, m1_c, m2_r, m2_c
     real(dp) :: my_lhs_factor, my_rhs_factor

     m1_transp = 'N'
     m2_transp = 'N'
     if (present(m1_transpose)) then
      if (m1_transpose) m1_transp = 'T'
     endif
     if (present(m2_transpose)) then
      if (m2_transpose) m2_transp = 'T'
     endif

     if (m1_transp == 'N') then
       m1_r = 1; m1_c = 2
     else
       m1_r = 2; m1_c = 1
     endif
     if (m2_transp == 'N') then
       m2_r = 1; m2_c = 2
     else
       m2_r = 2; m2_c = 1
     endif

     my_lhs_factor = optional_default(0.0_dp, lhs_factor)
     my_rhs_factor = optional_default(1.0_dp, rhs_factor)

     call check_size('lhs',lhs,(/size(matrix1,m1_r),size(matrix2,m2_c)/),'Matrix_Product')
     if (m2_transp == 'N') then
       call check_size('Matrix2',matrix2,(/size(matrix1,m1_c),size(matrix2,m2_c)/),'Matrix_Product')
     else
       call check_size('Matrix2',matrix2,(/size(matrix2,m2_c),size(matrix1,m1_c)/),'Matrix_Product')
     endif

     N = size(lhs,1) !# of rows of lhs
     M = size(lhs,2) !# of columns of rhs
     K = size(matrix1,m1_c) !!shared dimension
     maxd=max(N,M,K) 

     if (use_intrinsic_blas) then
        if (m1_transp == 'T') then
          if (m2_transp == 'T') then
            lhs=lhs*my_lhs_factor + my_rhs_factor*matmul(transpose(matrix1),transpose(matrix2))
          else
            lhs=lhs*my_lhs_factor + my_rhs_factor*matmul(transpose(matrix1),matrix2)
          endif
        else
          if (m2_transp == 'T') then
            lhs=lhs*my_lhs_factor + my_rhs_factor*matmul(matrix1,transpose(matrix2))
          else
            lhs=lhs*my_lhs_factor + my_rhs_factor*matmul(matrix1,matrix2)
          endif
        endif
     else
        call DGEMM(m1_transp,m2_transp, N,M,K,my_rhs_factor,matrix1,size(matrix1,1),matrix2,size(matrix2,1),&
	  my_lhs_factor,lhs,size(lhs,1))
     endif

   end subroutine matrix_product_sub_ddd

   ! subroutine form of m1(:,:) .mult. m2(:,:)
   !
   ! set first argument to matrix product (no temporary on the stack)
   subroutine matrix_product_sub_zzz(lhs, matrix1, matrix2, m1_transpose, m1_conjugate, &
    m2_transpose, m2_conjugate, lhs_factor, rhs_factor)
      complex(dp), intent(out) :: lhs(:,:)
      complex(dp), intent(in) :: matrix1(:,:), matrix2(:,:)
      logical, intent(in), optional :: m1_transpose, m1_conjugate, m2_transpose, m2_conjugate
      complex(dp), intent(in), optional :: lhs_factor, rhs_factor

     integer::M,N,maxd,K
     logical :: m1_transp, m2_transp, m1_conjg, m2_conjg
     character(len=1) :: m1_op, m2_op
     integer :: m1_r, m1_c, m2_r, m2_c
     complex(dp) :: my_lhs_factor, my_rhs_factor

     m1_transp = .false.
     m2_transp = .false.
     m1_conjg = .false.
     m2_conjg = .false.
     if (present(m1_transpose)) m1_transp = m1_transpose
     if (present(m2_transpose)) m2_transp = m2_transpose
     if (present(m1_conjugate)) m1_conjg = m1_conjugate
     if (present(m2_conjugate)) m2_conjg = m2_conjugate

     if (m1_conjg .and. m1_transp) call system_abort("Called matrix_product_sub_zzz with m1_transp and m1_conjg true")
     if (m2_conjg .and. m2_transp) call system_abort("Called matrix_product_sub_zzz with m2_transp and m2_conjg true")

     if (m1_transp) then
       m1_op = 'T'
     else if (m1_conjg) then
       m1_op = 'C'
     else
       m1_op = 'N'
     end if

     if (m2_transp) then
       m2_op = 'T'
     else if (m2_conjg) then
       m2_op = 'C'
     else
       m2_op = 'N'
     end if

     if (m1_op == 'N') then
       m1_r = 1; m1_c = 2
     else
       m1_r = 2; m1_c = 1
     endif
     if (m2_op == 'N') then
       m2_r = 1; m2_c = 2
     else
       m2_r = 2; m2_c = 1
     endif

     my_lhs_factor = optional_default(cmplx(0.0_dp, 0.0_dp, dp), lhs_factor)
     my_rhs_factor = optional_default(cmplx(1.0_dp, 0.0_dp, dp), rhs_factor)

     call check_size('lhs',lhs,(/size(matrix1,m1_r),size(matrix2,m2_c)/),'Matrix_Product')
     if (m2_op == 'N') then
       call check_size('Matrix2',matrix2,(/size(matrix1,m1_c),size(matrix2,m2_c)/),'Matrix_Product')
     else
       call check_size('Matrix2',matrix2,(/size(matrix2,m2_c),size(matrix1,m1_c)/),'Matrix_Product')
     endif

     N = size(lhs,1) !# of rows of lhs
     M = size(lhs,2) !# of columns of rhs
     K = size(matrix1,m1_c) !!shared dimension
     maxd=max(N,M,K) 

     if (use_intrinsic_blas) then
	if (m1_transp) then
	  if (m2_transp) then
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(transpose(matrix1),transpose(matrix2))
	  else if (m2_conjg) then
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(transpose(matrix1),conjg(transpose(matrix2)))
	  else
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(transpose(matrix1),matrix2)
	  endif
	else if (m1_conjg) then
	  if (m2_transp) then
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(conjg(transpose(matrix1)),transpose(matrix2))
	  else if (m2_conjg) then
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(conjg(transpose(matrix1)),conjg(transpose(matrix2)))
	  else
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(conjg(transpose(matrix1)),matrix2)
	  endif
	else
	  if (m2_transp) then
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(matrix1,transpose(matrix2))
	  else if (m2_conjg) then
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(matrix1,conjg(transpose(matrix2)))
	  else
	    lhs = lhs*my_lhs_factor + my_rhs_factor*matmul(matrix1,matrix2)
	  endif
	endif
     else
        call ZGEMM(m1_op,m2_op, N,M,K,my_rhs_factor,matrix1,size(matrix1,1),matrix2,size(matrix2,1),my_lhs_factor,lhs,size(lhs,1))
     endif

   end subroutine matrix_product_sub_zzz

   subroutine matrix_vector_product_sub_ddd(lhs, matrix, vector, m_transpose, &
    lhs_factor, rhs_factor)
      real(dp), intent(out) :: lhs(:)
      real(dp), intent(in) :: matrix(:,:), vector(:)
      logical, intent(in), optional :: m_transpose
      real(dp), intent(in), optional :: lhs_factor, rhs_factor

     integer::M,N,maxd
     character(len=1) :: m_transp
     integer :: m_r, m_c
     real(dp) :: my_lhs_factor, my_rhs_factor

     m_transp = 'N'
     if (present(m_transpose)) then
      if (m_transpose) m_transp = 'T'
     endif

     if (m_transp == 'N') then
       m_r = 1; m_c = 2
     else
       m_r = 2; m_c = 1
     endif

     my_lhs_factor = optional_default(0.0_dp, lhs_factor)
     my_rhs_factor = optional_default(1.0_dp, rhs_factor)

     call check_size('lhs',lhs,(/size(matrix,m_r)/),'Matrix_Product')

     N = size(lhs) !# of rows of lhs
     M = size(matrix,m_c) !!shared dimension
     maxd=max(N,M) 

     if (use_intrinsic_blas) then
        if (m_transp == 'T') then
          lhs=lhs*my_lhs_factor + my_rhs_factor*matmul(transpose(matrix),vector)
        else
          lhs=lhs*my_lhs_factor + my_rhs_factor*matmul(matrix,vector)
        endif
     else
        call DGEMV(m_transp,N,M,my_rhs_factor,matrix,size(matrix,1),vector,1,&
	  my_lhs_factor,lhs,1)
     endif

   end subroutine matrix_vector_product_sub_ddd
   ! m1(:,:) .mult. m2(:,:)
   !
   ! matrix multiplication (either lapack or intrinsic)                      
   function matrix_product_ddd(matrix1,matrix2) result (prodmatrix)
     real(dp),intent(in), dimension(:,:) :: matrix1
     real(dp),intent(in), dimension(:,:) :: matrix2
     real(dp), dimension(size(matrix1,1),size(matrix2,2)) ::prodmatrix

     call matrix_product_sub(prodmatrix, matrix1, matrix2)
   end function matrix_product_ddd

   ! m1(:,:) .mult. m2(:,:)
   !
   ! matrix multiplication (either lapack or intrinsic)                      
   function matrix_product_zzz(matrix1,matrix2) result (prodmatrix)
     complex(dp),intent(in), dimension(:,:) :: matrix1
     complex(dp),intent(in), dimension(:,:) :: matrix2
     complex(dp), dimension(size(matrix1,1),size(matrix2,2)) ::prodmatrix

     call matrix_product_sub(prodmatrix, matrix1, matrix2)
   end function matrix_product_zzz

   ! m1(:,:) .mult. m2'(:,:)
   !
   !% Calculate matrix product of 'matrix1' and 'matrix2' after transposing matrix2 
   !% (either using \textsc{lapack} or intrinsic multiplication, depending on the
   !& 'use_intrinsic_blas' setting).
   function matrix_multT(matrix1,matrix2) result (prodmatrix)
     real(dp),intent(in), dimension(:,:) :: matrix1
     real(dp),intent(in), dimension(:,:) :: matrix2
     real(dp), dimension(size(matrix1,1),size(matrix2,2)) ::prodmatrix

     call matrix_product_sub(prodmatrix, matrix1, matrix2, m1_transpose = .false., m2_transpose = .true.)

   end function matrix_multT


   ! a .feq. b
   !
   ! floating point logical comparison
   function real_feq(x,y) result(feq)

     real(dp), intent(in) :: x, y
     logical              :: feq

     feq = (abs(x-y) <= NUMERICAL_ZERO * (abs(x)+abs(y))/2.0_dp) .or. (abs(x-y) <= NUMERICAL_ZERO)
     
   end function real_feq

   ! a .fne. b
   !
   ! floating point logical comparison
   function real_fne(x,y) result(fne)

     real(dp), intent(in) :: x, y
     logical              :: fne
     
     fne = .not. real_feq(x, y) 
     
   end function real_fne

   ! za .feq. zb
   !
   ! complex logical comparison
   function complex_feq(x,y) result(feq)

     complex(dp), intent(in) :: x, y
     logical              :: feq

     if ( (abs(real(x-y)) > NUMERICAL_ZERO * abs(real(x))) .or. &
          (abs(aimag(x-y)) > NUMERICAL_ZERO * abs(aimag(x))) ) then
        feq = .false.
     else
        feq = .true.
     end if
     
   end function complex_feq

   ! za .fne. zb
   !
   ! complex logical comparison
   function complex_fne(x,y) result(fne)

     complex(dp), intent(in) :: x, y
     logical              :: fne
     
     if ( (abs(real(x-y)) > NUMERICAL_ZERO * abs(real(x))) .or. &
          (abs(aimag(x-y)) > NUMERICAL_ZERO * abs(aimag(x))) ) then
        fne = .true.
     else
        fne = .false.
     end if
     
   end function complex_fne
   
   ! m1 .feq. m2
   !
   ! matrix floating point comparison 
   function matrix_feq(matrix1,matrix2) result (feq)
     real(dp),intent(in), dimension(:,:) :: matrix1
     real(dp),intent(in), dimension(:,:) :: matrix2

     integer::i,j
     logical::feq
     
     call check_size('Matrix2',matrix2,shape(matrix1),'Matrix_FEQ')
     
     feq =.true.
     do j=1,size(matrix1,2)
        do i=1,size(matrix1,1)
           if (matrix1(i,j).fne.matrix2(i,j)) then
              feq=.false.
              return
           end if
        end do
     end do
     
   end function matrix_feq
   
   ! m1 .fne. m2
   !
   ! matrix floating point comparison 
   function matrix_fne(matrix1,matrix2) result (fne)
     real(dp),intent(in), dimension(:,:) :: matrix1
     real(dp),intent(in), dimension(:,:) :: matrix2

     integer::i,j
     logical::fne
     
     if(size(matrix1,1) /= size(matrix2,1) .or. size(matrix1,2) /= size(matrix2,2)) & 
          call check_size('Matrix2',matrix2,shape(matrix1),'Matrix_FNE')
     
     fne =.false.
     do j=1,size(matrix1,2)
        do i=1,size(matrix1,1)
           if (matrix1(i,j).FNE.matrix2(i,j)) then
              fne=.true.
              return
           end if
        end do
     end do
     
   end function matrix_fne
   
   ! is_square(matrix)
   !
   ! tells if the matrix is square
   function matrix_square(matrix) result(sq)
     real(dp),intent(in), dimension(:,:) :: matrix
     logical::sq
     
     if (size(matrix,1).EQ.size(matrix,2)) then
        sq=.true.
     else
        sq=.false.
     end if
     
   end function matrix_square

   function int_matrix_square(matrix) result(sq)

     integer, intent(in), dimension(:,:) :: matrix
     logical::sq
     
     sq = size(matrix,1) == size(matrix,2)
     
   end function int_matrix_square

   ! is_square(matrix_z)
   !
   ! tells if the complex matrix is square
   function matrix_z_square(matrix_z) result(sq)
     complex(dp),intent(in), dimension(:,:) ::matrix_z
     logical::sq
     
     if (size(matrix_z,1).EQ.size(matrix_z,2)) then
        sq=.true.
     else
        sq=.false.
     end if
     
   end function matrix_z_square


  ! symmetrise matrix A -> (A + transpose(A))/2
  subroutine matrix_symmetrise(matrix)
    real(dp), intent(inout), dimension(:,:) :: matrix
    integer::i,j,n
    real(dp)::tmp
    
    if (.not.is_square(matrix)) call system_abort('Matrix_Symmetrise: Matrix is not square')

    n=size(matrix,1)
    do i=1,n
       do j=i+1,n
          tmp=0.5_dp*(matrix(i,j)+matrix(j,i))
          matrix(i,j)=tmp
          matrix(j,i)=tmp
       end do
    end do

  end subroutine matrix_symmetrise

  !returns trace of a matrix
  function matrix_trace(matrix) result(tr)
    real(dp),intent(in), dimension(:,:) ::matrix
    real(dp)::tr
    integer::i,N

    N = min(size(matrix,1),size(matrix,2))
   
    tr=0.0_dp
    do i=1,N
       tr=tr+matrix(i,i)
    end do

  end function matrix_trace

#ifdef HAVE_QP
  function matrix_trace_q(matrix) result(tr)
    real(qp),intent(in), dimension(:,:) ::matrix
    real(qp)::tr
    integer::i,N

    N = min(size(matrix,1),size(matrix,2))
   
    tr=0.0_qp
    do i=1,N
       tr=tr+matrix(i,i)
    end do

  end function matrix_trace_q
#endif

  !returns trace of the result matrix
  function matrix_trace_mult(matrixA, matrixB) result(trm)
    real(dp),intent(in), dimension(:,:) ::matrixA, matrixB
    real(dp)::trm
    integer::i, N

    N = size(matrixA, 1)
    if(size(matrixB, 2) /= N) call system_abort("matrix_trace_mult: size(matrixB, 2) /= N")
    if(size(matrixA, 2) /= size(matrixB, 1)) call system_abort("size(matrixA, 2) /= size(matrixB, 1)")
    
    trm=0.0_dp
    do i=1,N
       trm=trm + (matrixA(i,:) .dot. matrixB(:,i))
    end do

  end function matrix_trace_mult

  subroutine matrix_add_identity_r(matrix)
    real(dp), intent(inout) :: matrix(:,:)

    integer i

    if (.not.is_square(matrix)) call system_abort('Matrix_add_identity: Matrix is not square')

    do i=1, size(matrix,1)
      matrix(i,i) = matrix(i,i) + 1.0_dp
    end do
  end subroutine matrix_add_identity_r

  subroutine matrix_add_identity_c(matrix)
    complex(dp), intent(inout) :: matrix(:,:)

    integer i

    if (.not.is_square(matrix)) call system_abort('Matrix_add_identity: Matrix is not square')

    do i=1, size(matrix,1)
      matrix(i,i) = matrix(i,i) + 1.0_dp
    end do
  end subroutine matrix_add_identity_c

  subroutine matrix_add_xidentity_r(matrix,r)
    real(dp), intent(inout) :: matrix(:,:)
    real(dp), intent(in)    :: r

    integer :: i

    if (.not.is_square(matrix)) call system_abort('Matrix_add_xidentity: Matrix is not square')

    do i=1, size(matrix,1)
      matrix(i,i) = matrix(i,i) + r
    end do
  end subroutine matrix_add_xidentity_r

  subroutine matrix_add_xidentity_c(matrix,c)
    complex(dp), intent(inout) :: matrix(:,:)
    complex(dp), intent(in)    :: c

    integer :: i

    if (.not.is_square(matrix)) call system_abort('Matrix_add_xidentity: Matrix is not square')

    do i=1, size(matrix,1)
      matrix(i,i) = matrix(i,i) + c
    end do
  end subroutine matrix_add_xidentity_c

  
  ! the following stuff only with lapack        
  ! diagonalise matrix, only symmetric case
  subroutine matrix_diagonalise(this,evals,evects, error)
    real(dp),intent(in), dimension(:,:) :: this
    real(dp),intent(inout), dimension(:) ::evals
    real(dp),intent(inout), target, optional, dimension(:,:) :: evects
    integer, intent(out), optional :: error

    real(8),allocatable::WORK(:), r8_evals(:)
    real(8), pointer :: r8_evects(:,:)
    integer::N,INFO,LWORK

    INIT_ERROR(error)

    N=size(this,2)
    call check_size('Eigenvalue Vector',evals,N,'Matrix_Diagonalise')

    if (present(evects)) &
      call check_size('Eigenvector Array',evects,shape(this),'Matrix_Diagonalise')

    if (is_symmetric(this)) then

       LWORK=3*N
       allocate(WORK(LWORK))     
       allocate(r8_evals(N))

       if (present(evects) .and. dp == 8) then
	 r8_evects => evects
       else
	 allocate(r8_evects(N,N))
       endif
       r8_evects = this

       if (present(evects)) then
	 call DSYEV('V','U',N,r8_evects,N,r8_evals,WORK,LWORK,INFO)
       else
	 call DSYEV('N','U',N,r8_evects,N,r8_evals,WORK,LWORK,INFO)
	endif

       if (present(evects) .and. dp /= 8) evects = r8_evects
       evals = r8_evals

       deallocate(WORK)
       deallocate(r8_evals)
       if (.not. (present(evects) .and. dp == 8)) deallocate(r8_evects)

       if (INFO /= 0) then
	  mainlog%mpi_all_inoutput_flag=.true.
	  call print ('Matrix_diagonalise: Error in calling DSYEV! (info = '//INFO//')', PRINT_ALWAYS)
	  if (INFO < 0) then
	     call print ('  '//(-INFO)//' argument had illegal value', PRINT_ALWAYS)
	  else if (INFO <= N) then
	     call print ('  '//INFO//' off-diagonal elements of an intermediate tridiagonal form did not converge to zero', PRINT_ALWAYS)
	  endif
	  mainlog%mpi_all_inoutput_flag=.false.
	  RAISE_ERROR ('Matrix_diagonalise: Error in calling DSYEV! (info = '//INFO//')', error)
       endif
    else
       ! Why not print a warning then call more general diagonalise?
       RAISE_ERROR('Matrix_diagonalise: Non symmetric diagonalisation is not permitted',error) 
    end if
 
  end subroutine matrix_diagonalise

  ! the following stuff only with lapack        
  ! diagonalise complex matrix, only hermitian positive definite case
  subroutine matrix_z_diagonalise(this,evals,evects,error)
    complex(dp),intent(in), dimension(:,:) :: this
    real(dp),intent(inout), dimension(:) ::evals
    complex(dp),intent(inout), optional, target, dimension(:,:) :: evects
    integer, intent(out), optional :: error
    integer::N,INFO,LWORK
    integer NB
    integer, external :: ILAENV

    complex(8), pointer :: z8_evects(:,:)
    real(8), allocatable :: r8_evals(:), RWORK(:)
    complex(8), pointer :: WORK(:)

    INIT_ERROR(error)

    N=size(this,2)
    call check_size('Eigenvalue Vector',evals,N,'Matrix_z_Diagonalise')

    if (present(evects)) &
      call check_size('Eigenvector Array',evects,shape(this),'Matrix_z_Diagonalise')

    if (is_hermitian(this)) then

       NB = ILAENV(1, "ZHETRD", "U", N, N, N, N)
       LWORK=(NB+1)*N
       allocate(WORK(LWORK))     
       allocate(RWORK(3*N-2))
       allocate(r8_evals(N))

       if (present(evects) .and. dp == 8) then
	z8_evects => evects
       else
	 allocate(z8_evects(N,N))
       endif
       z8_evects = this

       if (present(evects)) then
	 call ZHEEV('V','U',N,z8_evects,N,r8_evals,WORK,LWORK,RWORK,INFO)
       else
	 call ZHEEV('N','U',N,z8_evects,N,r8_evals,WORK,LWORK,RWORK,INFO)
       endif

       if (present(evects) .and. dp /= 8) evects = z8_evects
       evals = r8_evals

       deallocate(WORK)
       deallocate(RWORK)
       deallocate(r8_evals)
       if (.not.(present(evects) .and. dp == 8)) deallocate(z8_evects)

       if (INFO /= 0) then
	  mainlog%mpi_all_inoutput_flag=.true.
	  call print ('Matrix_z_diagonalise: Error in calling ZHEEV! (info = '//INFO//')', PRINT_ALWAYS)
	  if (INFO < 0) then
	     call print ('  '//(-INFO)//' argument had illegal value', PRINT_ALWAYS)
	  else if (INFO <= N) then
	     call print ('  '//INFO//' off-diagonal elements of an intermediate tridiagonal form did not converge to zero', PRINT_ALWAYS)
	  endif
	  mainlog%mpi_all_inoutput_flag=.false.
	  RAISE_ERROR ('Matrix_z_diagonalise: Error in calling ZHEEV! (info = '//INFO//')', error)
       endif
    else
       ! why not print a warning and call more general diagonalise instead?
       RAISE_ERROR ('Matrix_z_diagonalise: Non hermitian diagonalisation is not permitted', error)
    end if
 
  end subroutine matrix_z_diagonalise

  ! generalised eigenproblem
  ! just works for symmetric systems 
  subroutine  matrix_diagonalise_generalised(this,other,evals,evects,error)
    real(dp),intent(in), dimension(:,:) :: this
    real(dp),intent(in), dimension(:,:) :: other
    real(dp),intent(inout), dimension(:) :: evals
    real(dp),intent(inout), dimension(:,:) :: evects
    integer, intent(out), optional :: error

    real(dp), allocatable :: other_copy(:,:)
    real(dp), allocatable :: WORK(:)
    integer::N,INFO,LWORK
    integer, external :: ILAENV
    integer NB

    INIT_ERROR(error)

    if (dp /= 8) call system_abort("matrix_diagonalise_generalised: no workaround for LAPACK assuming 8 byte double, but dp(="//dp//") /= 8")

    NB = ILAENV(1, "DSYTRD", "U", N, N, N, N)
    N=size(this,1)
    LWORK=(NB+2)*N
    allocate(WORK(LWORK))
    allocate(other_copy(N,N))

    call check_size('Eigenvalue vector',evals,N,'Matrix_Diagonalise_Generalised')
    call check_size('Eigenvector Array',evects,shape(this),'Matrix_Diagonalise_Generalised')

    evects = this
    other_copy = other

    call DSYGV(1,'V','U',N,evects,N,other_copy,N,evals,WORK,LWORK,INFO)

    deallocate(WORK)
    deallocate(other_copy)

    if (INFO /= 0) then
       mainlog%mpi_all_inoutput_flag=.true.
       call print ('Matrix_diagonalise_generalised: Error in calling DSYGV! (info = '//INFO//')', PRINT_ALWAYS)
       if (INFO < 0) then
	  call print ('  '//(-INFO)//' argument had illegal value', PRINT_ALWAYS)
       else if (INFO <= N) then
	  call print ('  '//INFO//' off-diagonal elements of an intermediate tridiagonal form did not converge to zero', PRINT_ALWAYS)
       else
	  call print ('  '//(INFO-N)//' leading minor of B is not positive definite, factorization failed')
       endif
       mainlog%mpi_all_inoutput_flag=.false.
       RAISE_ERROR ('Matrix_diagonalise_generalised: Error in calling DSYGV! (info = '//INFO//')', error)
    endif

  end subroutine matrix_diagonalise_generalised

  ! general eigenproblem
  ! just works for hermitian systems 
  subroutine  matrix_z_diagonalise_generalised(this,other,evals,evects, error )
    complex(dp),intent(in), dimension(:,:) :: this
    complex(dp),intent(in), dimension(:,:) :: other
    real(dp),intent(inout), dimension(:) :: evals
    complex(dp),intent(inout), dimension(:,:) :: evects
    integer, intent(out), optional :: error

    complex(dp), allocatable::WORK(:)
    real(dp), allocatable::RWORK(:)
    integer::N,INFO,LWORK
    integer, external :: ILAENV
    complex(dp), allocatable :: other_copy(:,:)

    integer NB

    INIT_ERROR(error)

    if (dp /= 8) call system_abort("matrix_z_diagonalise_generalised: no workaround for LAPACK assuming 8 byte double, but dp(="//dp//") /= 8")
  
    N=size(this,2)
    NB = ILAENV(1, "ZHETRD", "U", N, N, N, N)
    LWORK=(NB+1)*N
    allocate(WORK(LWORK))
    allocate(RWORK(3*N-2))
    allocate(other_copy(N,N))

    call check_size('Eigenvalue vector',evals,N,'Matrix_z_Diagonalise_Generalised')
    call check_size('Eigenvector Array',evects,shape(this),'Matrix_z_Diagonalise_Generalised')

    evects = this
    other_copy = other

    call ZHEGV(1,'V','U',N,evects,N,other_copy,N,evals,WORK,LWORK,RWORK,INFO)

    deallocate(WORK)
    deallocate(RWORK)
    deallocate(other_copy)

    if (INFO /= 0) then
       mainlog%mpi_all_inoutput_flag=.true.
       call print ('Matrix_z_diagonalise_generalised: Error in calling ZHEGV! (info = '//INFO//')', PRINT_ALWAYS)
       if (INFO < 0) then
	  call print ('  '//(-INFO)//' argument had illegal value', PRINT_ALWAYS)
       else if (INFO <= N) then
	  call print ('  '//INFO//' off-diagonal elements of an intermediate tridiagonal form did not converge to zero', PRINT_ALWAYS)
       else
	  call print ('  '//(INFO-N)//' leading minor of B is not positive definite, factorization failed')
       endif
       mainlog%mpi_all_inoutput_flag=.false.
       RAISE_ERROR ('Matrix_z_diagonalise_generalised: Error in calling ZHEGV! (info = '//INFO//')', error)
    endif

  end subroutine matrix_z_diagonalise_generalised

  subroutine matrix_nonsymmetric_diagonalise(this, eval, l_evects, r_evects, error)
     real(dp), dimension(:,:), intent(in) :: this
     complex(dp), dimension(:), intent(out) :: eval
     complex(dp), dimension(:,:), intent(out), optional, target :: l_evects, r_evects
     integer, optional, intent(out) :: error

     real(dp), dimension(:), allocatable :: wr, wi, work
     real(dp), dimension(:,:), allocatable :: a
     real(dp), dimension(:,:), allocatable :: vl, vr
     integer :: n, lwork, info, i
     character(len=1) :: cl, cr
     integer :: ldvl, ldvr

     INIT_ERROR(error)

     n = size(this,1)
     allocate(a(n,n), wr(n), wi(n))
     if (present(r_evects)) then
	allocate(vr(n,n))
	cr = 'V'
	ldvr = n
     else
	allocate(vr(1,1))
	cr ='N'
	ldvr = 1
     endif
     if (present(l_evects)) then
        allocate(vl(n,n))
	cl='V'
	ldvl = n
     else
	allocate(vl(1,1))
	cl='N'
	ldvl = 1
     endif
     a = this

     allocate(work(1))
     lwork = -1
     call dgeev(cl, cr, n, a, n, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
     if (info /= 0) then
        RAISE_ERROR("matrix_nonsymmetric_diagonalise first dgeev call (determine lwork) failed with INFO="//info, error)
     endif
     lwork = ceiling(work(1))
     deallocate(work)
     allocate(work(lwork))
     call dgeev(cl, cr, n, a, n, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
     if (info /= 0) then
        RAISE_ERROR("matrix_nonsymmetric_diagonalise second dgeev call (actual calculation) failed with INFO="//info, error)
     endif

     eval = cmplx(wr,wi)
     if (present(l_evects) .or. present(r_evects)) then
	i = 1
	do while (i <= n)
	 if (wi(i) == 0) then
	    if (present(l_evects)) l_evects(:,i) = vl(:,i)
	    if (present(r_evects)) r_evects(:,i) = vr(:,i)
	 else
	    if (present(l_evects)) l_evects(:,i) = cmplx(vl(:,i), vl(:,i+1))
	    if (present(l_evects)) l_evects(:,i+1) = cmplx(vl(:,i), -vl(:,i+1))
	    if (present(r_evects)) r_evects(:,i) = cmplx(vr(:,i), vr(:,i+1))
	    if (present(r_evects)) r_evects(:,i+1) = cmplx(vr(:,i), -vr(:,i+1))
	    i = i + 1
	 endif
	 i = i + 1
      end do
     endif

     deallocate( a, wr, wi, work )
     if (.not. present(l_evects)) deallocate(vl)
     if (.not. present(r_evects)) deallocate(vr)


  endsubroutine matrix_nonsymmetric_diagonalise

  subroutine test_eigensys(A, B, evals, evecs)
    real(dp), intent(in) :: A(:,:), B(:,:)
    real(dp), intent(in) :: evals(:)
    real(dp), intent(in) :: evecs(:,:)

    real(dp), allocatable :: t1(:,:), t2(:,:)
    integer n, i, j

    n = size(evals)

    allocate(t1(n,n))
    allocate(t2(n,n))

    t1 = matmul(A, evecs)
    t2 = matmul(B, evecs)

    print '("test eigensys")'

    do i=1, N
    do j=1, N
      if (abs(t2(j,i)) > 1e-10_dp) then
	print '(I0," ",4F30.20)', i, evals(i), t1(j,i), t2(j,i), t1(j,i)/t2(j,i)
      else
	print '(I0," ",3F30.20, "-")', i, evals(i), t1(j,i), t2(j,i)
      endif
    end do
    end do

  end subroutine test_eigensys


  ! calculate inverse, if the matrix is symmetric uses a quicker way, otherways a more general
  ! if inverse is missing, do in-place
  subroutine matrix_inverse(matrix,inverse,positive_in)
    real(dp),intent(in), dimension(:,:), target::matrix
    real(dp),intent(out), dimension(:,:), optional, target::inverse
    logical,optional::positive_in

    integer::N,INFO,i,j,LWORK
    integer,allocatable::IPIV(:)
    real(dp),allocatable::WORK(:)
    logical::positive
    !IMPROVE real(dp), allocatable :: u_lu(:,:), b(:,:), ferr(:), berr(:)
    !IMPROVE integer, allocatable :: iwork(:)

    real(dp), pointer :: u_inverse(:,:)

    if (present(inverse)) then
      u_inverse => inverse
    else
      u_inverse => matrix
    endif
 
    if (present(positive_in)) then
       positive = positive_in
    else
       positive = .true.
    end if

    call check_size('Inverse',u_inverse,shape(matrix),'Matrix_Inverse')   

    N = size(matrix,2)

     if (present(inverse)) u_inverse = matrix

    ! try the symmetric, positive definite case

    if( matrix_is_symmetric(matrix) .and. positive)  then   

       ! cholesky

       ! this only works for positive definite
       call DPOTRF('U', N, u_inverse, N, INFO)

       if(INFO > 0) then
          write (line,'(a)') 'Matrix_Inverse: Matrix is not positive definite, switching to general case inversion!'
          call print(line)

          ! we branch to the general case below
       else
          if (INFO < 0) then
             write(line,'(a,i0,a)')'Matrix_Inverse: Error in calling DPOTRF (',info,')'
             call system_abort(line)
          end if
          ! now do the inverse
          ! again, this call is only for positive definite, but now we should be
          call DPOTRI('U', N, u_inverse, N, INFO)
          if (INFO /= 0) then
             write(line,'(a,i0,a)')'Matrix_Inverse: Error in calling DPOTRI (',info,')'
             call system_abort(line)
          end if
          ! filling the lower part of symmetric matrix   
          do i=1,N
             do j=i+1,N
                u_inverse(j,i) = u_inverse(i,j)
             enddo
          enddo
          return
       end if
    end if

    ! do the general case

    LWORK = 3*N
    allocate(WORK(LWORK))
    allocate(IPIV(N))

    ! LU factorization

    call DGETRF(N, N, u_inverse, N, IPIV, info)
    if (INFO /= 0) call system_abort('Error in calling DGETRF (info = '//info//')')

    !IMPROVE allocate(u_lu(N,N))
    !IMPROVE u_lu = u_inverse

    !inverse

    call DGETRI(N, u_inverse, N, IPIV, WORK, LWORK, info)
    if (INFO /= 0) call system_abort('Error in calling DGETRI (info = '//info//')') 

    !IMPROVE deallocate(WORK)
    !IMPROVE allocate(b(N,N),ferr(N),berr(N),WORK(3*N),iwork(N))
    !IMPROVE b = 0.0_dp
    !IMPROVE call add_identity(b)
    !IMPROVE call dgerfs('N',N,N,matrix,N,u_lu,N,ipiv,b,N,u_inverse,N,ferr,berr,WORK,iwork,info)

    deallocate(WORK,IPIV)

  end subroutine matrix_inverse

  ! calculate inverse, if the matrix is symmetric uses a quicker way, otherways a more general
  ! if inverse is missing, do in-place
  subroutine matrix_z_inverse(matrix,inverse,positive_in)
    complex(dp),intent(in), dimension(:,:), target::matrix
    complex(dp),intent(out), dimension(:,:), optional, target::inverse
    logical,optional::positive_in

    integer::N,INFO,i,j,LWORK
    integer,allocatable::IPIV(:)
    complex(dp),allocatable::WORK(:)
    logical::positive

    complex(dp), pointer :: u_inverse(:,:)

    if (present(inverse)) then
      u_inverse => inverse
      u_inverse = matrix
    else
      u_inverse => matrix
    endif
 
    if (present(positive_in)) then
       positive = positive_in
    else
       positive = .true.
    end if

    call check_size('Inverse',u_inverse,shape(matrix),'Matrix_Inverse')   

    N = size(matrix,2)

    if (matrix_z_is_hermitian(matrix)) then


      if (positive) then

	! cholesky

	! this only works for hermitian positive definite
	call ZPOTRF('U', N, u_inverse, N, INFO)

	if(INFO > 0) then
	  write (line,'(a)') 'Matrix_Inverse: Matrix is not positive definite, switching to general case inversion!'
	  call print(line)
	  ! we branch to the general case below
	else
	  if (INFO < 0) then
	    write(line,'(a,i0,a)')'Matrix_Inverse: Error in calling ZPOTRF (',info,')'
	    call system_abort(line)
	  end if
	  ! now do the inverse
	  ! again, this call is only for positive definite, but now we should be
	  call ZPOTRI('U', N, u_inverse, N, INFO)
	  if (INFO /= 0) then
	    write(line,'(a,i0,a)')'Matrix_Inverse: Error in calling ZPOTRI (',info,')'
	    call system_abort(line)
	  end if
	  ! filling the lower part of hermitian matrix   
	  do i=1,N
	    do j=i+1,N
	      u_inverse(j,i) = conjg(u_inverse(i,j))
	    enddo
	  enddo
	  return
	end if ! zpotrf INFO > 0

      endif ! positive

      lwork = 16*N
      allocate(work(lwork))
      allocate(ipiv(N))

      ! Bunch-Kaufman factorization
      call ZHETRF('U', N, u_inverse, N, ipiv, work, lwork, info)
      if (INFO /= 0) call system_abort('Error in calling ZHETRF (info = '//info//')') 

      ! inverse
      call ZHETRI('U', N, u_inverse, N, ipiv, work, info)
      if (INFO /= 0) call system_abort('Error in calling ZHETRI (info = '//info//')') 

      ! filling the lower part of hermitian matrix   
      do i=1,N
	do j=i+1,N
	  u_inverse(j,i) = conjg(u_inverse(i,j))
	enddo
      enddo

      deallocate(work,ipiv)
      return
    else if (matrix_z_is_symmetric(matrix)) then
      lwork = 16*N
      allocate(work(lwork))
      allocate(ipiv(N))

      ! LU factorization
      call ZSYTRF('U', N, u_inverse, N, ipiv, work, lwork, info)
      if (INFO /= 0) call system_abort('Error in calling ZSYTRF (info = '//info//')') 

      ! inverse
      call ZSYTRI('U', N, u_inverse, N, ipiv, work, info)
      if (INFO /= 0) call system_abort('Error in calling ZSYTRI (info = '//info//')') 

      ! filling the lower part of symmetric matrix   
      do i=1,N
	do j=i+1,N
	  u_inverse(j,i) = u_inverse(i,j)
	enddo
      enddo

      deallocate(work, ipiv)
      return
    end if ! matrix is hermitian

    ! do the general case

    LWORK = 16*N
    allocate(WORK(LWORK))
    allocate(IPIV(N))

    ! LU factorization
    call ZGETRF(N, N, u_inverse, N, IPIV, info)
    if (INFO /= 0) call system_abort('Error in calling ZGETRF (info = '//info//')')

    !inverse
    call ZGETRI(N, u_inverse, N, IPIV, WORK, LWORK, info)
    if (INFO /= 0) call system_abort('Error in calling ZGETRI (info = '//info//')') 

    deallocate(WORK,IPIV)
    return

  end subroutine matrix_z_inverse

  ! Cholesky factorisation of a symmetric matrix.
  ! Various checks (size, symmetricity etc.) might be useful later.
  subroutine LA_Matrix_Initialise(this,matrix)

     type(LA_Matrix), intent(inout) :: this
     real(qp), dimension(:,:), intent(in) :: matrix

     if(this%initialised) call finalise(this)

     this%n = size(matrix,1)
     this%m = size(matrix,2)
     allocate(this%matrix(this%n,this%m), this%factor(this%n,this%m), this%s(this%n),this%tau(this%m) )

     this%matrix = matrix
     this%initialised = .true.

  endsubroutine LA_Matrix_Initialise

  subroutine LA_Matrix_Update(this,matrix)

     type(LA_Matrix), intent(inout) :: this
     real(qp), dimension(:,:), intent(in) :: matrix

     integer :: factorised

     factorised = this%factorised

     if(this%initialised) then
        if( all(shape(matrix) == (/this%n,this%m/)) ) then
           this%matrix = matrix
        else
           call initialise(this,matrix)
        endif
     else
        call initialise(this,matrix)
     endif

     select case(factorised)
     case(CHOLESKY)
        call LA_Matrix_Factorise(this)
     case(QR)
        call LA_Matrix_QR_Factorise(this)
     endselect

  endsubroutine LA_Matrix_Update

  subroutine LA_Matrix_Finalise(this)

     type(LA_Matrix), intent(inout) :: this

     if(.not. this%initialised) return

     this%n = 0
     this%m = 0
     if(allocated(this%matrix) ) deallocate(this%matrix)
     if(allocated(this%factor) ) deallocate(this%factor)
     if(allocated(this%s) ) deallocate(this%s)
     if(allocated(this%tau) ) deallocate(this%tau)
     this%initialised = .false.
     this%equilibrated = .false.
     this%factorised = NOT_FACTORISED

  endsubroutine LA_Matrix_Finalise

  subroutine LA_Matrix_Factorise(this,factor,error)

     type(LA_Matrix), intent(inout) :: this
     real(qp), dimension(:,:), intent(out), optional :: factor
     integer, optional, intent(out) :: error

     integer :: i, j, info
     real(dp) :: scond, amax

     INIT_ERROR(error)

     if(.not. this%initialised) then
        RAISE_ERROR('LA_Matrix_Factorise: object not initialised',error)
     endif

     if( this%n /= this%m ) then
        RAISE_ERROR('LA_Matrix_Factorise: matrix not square',error)
     endif

     this%s = 1.0_qp

     do i = 1, this%n
        this%s(i) = 1.0_qp / sqrt(this%matrix(i,i))
     enddo
     scond = maxval(this%s) / minval(this%s)
     amax = maxval(this%matrix)

     this%equilibrated = ( scond < 0.1_qp )

     if( this%equilibrated ) then
        do i = 1, this%n
           this%factor(:,i) = this%matrix(:,i)*this%s(:)*this%s(i)
        enddo
     else
        this%factor = this%matrix
     endif

#if defined(HAVE_QP) || defined(ALBERT_LAPACK)
     call qpotrf(this%factor,info)
#else
     call dpotrf('L', this%n, this%factor, this%n, info)
     do i = 2, this%n
        do j = 1, i
           this%factor(j,i) = this%factor(i,j)
        enddo
     enddo
#endif

     if( info /= 0 ) then
        RAISE_ERROR('LA_Matrix_Factorise: cannot factorise, error: '//info, error)
     endif

     if( present(factor) ) then
        if( this%equilibrated ) then
           factor = 0.0_qp
           do i = 1, this%n
              do j = 1, this%n
                 factor(j,i) = this%factor(j,i) / this%s(i)
              enddo
           enddo
        else
           factor = this%factor
        endif
     endif

     this%factorised = CHOLESKY
        
  endsubroutine LA_Matrix_Factorise

  subroutine qpotrf(this,info)
     real(qp), dimension(:,:), intent(inout) :: this
     integer, intent(out), optional :: info
     real(qp), dimension(:), allocatable :: v, w

     integer :: i, j, k, p, mu, n
#ifdef _OPENMP
     integer, external :: omp_get_num_threads, omp_get_thread_num
#endif

     n = size(this,1)

     if(present(info)) info = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     another implementation
!     this(:,1) = this(:,1) / sqrt(this(1,1))
!     do j = 2, n
!        this(j:n,j) = this(j:n,j) - matmul( this(j:n,1:j-1), this(j,1:j-1) )
!        if(this(j,j)<0.0_qp) then
!           if(present(info)) then
!              info = j
!              return
!           else
!              call system_abort('my_potrf: trouble')
!           endif
!        endif
!        this(j:n,j) = this(j:n,j)/sqrt(this(j,j))
!     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     p = 1
     mu = 1

!$omp parallel private(mu,v,w) shared(info)
!$ p=omp_get_num_threads()
!$ mu=omp_get_thread_num()+1
     allocate(v(n),w(n))
     do k = 1, n
        if( mu == 1 ) then
           v(k:n) = this(k:n,k)
           if(v(k)<0.0_qp) then
              if(present(info)) then
                 info = k
!                 return
              else
                 call system_abort('my_potrf: trouble')
              endif
           endif
           v(k:n) = v(k:n) / sqrt(v(k))
           this(k:n,k) = v(k:n)
        endif
!$omp barrier        
        v(k+1:n) = this(k+1:n,k)
        do j = k+mu,n,p
           w(j:n) = this(j:n,j)
           w(j:n) = w(j:n) - v(j)*v(j:n)
           this(j:n,j) = w(j:n)
        enddo
!$omp barrier        
     enddo
     deallocate(v,w)
!$omp end parallel
           
     do i = 2, n
        do j = 1, i
           this(j,i) = this(i,j)
        enddo
     enddo

  endsubroutine qpotrf

  subroutine LA_Matrix_Inverse(this,inverse,error)

     type(LA_Matrix), intent(inout) :: this
     real(qp), dimension(:,:), intent(out) :: inverse
     integer, optional, intent(out) :: error

     integer :: i, j, info

     INIT_ERROR(error)

     if( this%factorised == NOT_FACTORISED ) then
        call LA_Matrix_Factorise(this, error=error)
     elseif( this%factorised /= CHOLESKY ) then
        RAISE_ERROR('LA_Matrix_Inverse: matrix not Cholesky-factorised',error)
     endif

     call check_size('inverse',inverse,shape(this%factor),'LA_Matrix_Inverse',error=error)
     PASS_ERROR(error)

     inverse = this%factor

#if defined(HAVE_QP) || defined(ALBERT_LAPACK)
     call qpotri(inverse,info)
#else
     call dpotri('U', this%n, inverse, this%n, info)
#endif

     if( info /= 0 ) then
        RAISE_ERROR('LA_Matrix_Inverse: cannot invert, error: '//info, error)
     endif

     do i = 1, this%n
        do j = i+1, this%n
           inverse(j,i) = inverse(i,j)
        enddo
     enddo
     
     if(this%equilibrated) then
        do i = 1, this%n
           inverse(:,i) = inverse(:,i)*this%s(:)*this%s(i)
        enddo
     endif

  endsubroutine LA_Matrix_Inverse

  subroutine qpotri(this,info)

     real(qp), dimension(:,:), intent(inout) :: this
     integer, optional, intent(out) :: info

     real(qp), dimension(:,:), allocatable :: one
     integer :: i, n

     n = size(this,1)
     allocate(one(n,n))
     one = 0.0_qp
     do i = 1, n
        one(i,i) = 1.0_qp
     enddo
     call qpotrs(this,one)
     this = one
     deallocate(one)

     if(present(info)) info = 0
     
  endsubroutine qpotri

  subroutine LA_Matrix_Solve_Vector(this,b,x,refine,error)

     type(LA_Matrix), intent(inout) :: this
     real(qp), dimension(:), intent(in) :: b
     real(qp), dimension(:), intent(out) :: x
     logical, intent(in), optional :: refine
     integer, intent(out), optional :: error

     real(qp), dimension(:,:), allocatable :: my_x

     INIT_ERROR(error)

     if( (size(b) /= this%n) .or. (size(x) /= this%n) ) then
        RAISE_ERROR('LA_Matrix_Solve_Vector: length of b or x is not n',error)
     endif

     allocate(my_x(this%n,1))
     call LA_Matrix_Solve_Matrix(this,reshape(b,(/this%n,1/)),my_x,refine,error)
     x = my_x(:,1)
     deallocate(my_x)

  endsubroutine LA_Matrix_Solve_Vector

  subroutine LA_Matrix_Solve_Matrix(this,b,x,refine,error)

     type(LA_Matrix), intent(inout) :: this
     real(qp), dimension(:,:), intent(in) :: b
     real(qp), dimension(:,:), intent(out) :: x
     logical, intent(in), optional :: refine
     integer, intent(out), optional :: error

     real(qp), dimension(:,:), allocatable :: my_b, my_x
     integer :: i, m, info

     logical :: my_refine

     INIT_ERROR(error)

     if( size(b,1) /= this%n ) then
        RAISE_ERROR('LA_Matrix_Solve_Matrix: first dimension of b is not n',error)
     endif
     
     my_refine = optional_default(.false.,refine)

     if( this%factorised == NOT_FACTORISED ) then
        call LA_Matrix_Factorise(this,error=error)
     elseif( this%factorised /= CHOLESKY ) then
        call system_abort('LA_Matrix_Solve_Matrix: matrix not Cholesky-factorised')
     endif

     m = size(b,2)
     allocate(my_b(this%n,m),my_x(this%n,m)) !, work(3*this%n),iwork(this%n),ferr(m), berr(m))

     if( this%equilibrated ) then
        do i = 1, m
           my_x(:,i) = b(:,i)*this%s
        enddo
        my_b = my_x
     else
        my_x = b
        my_b = my_x
     endif

#if defined(HAVE_QP) || defined(ALBERT_LAPACK)
     call qpotrs( this%factor, my_x, info )
#else
     call dpotrs( 'U', this%n, m, this%factor, this%n, my_x, this%n, info )
#endif

     if( this%equilibrated ) then
        do i = 1, m
           x(:,i) = my_x(:,i)*this%s
        enddo
     else
        x = my_x
     endif

     if( info /= 0 ) then
        RAISE_ERROR('LA_Matrix_Solve_Matrix: cannot solve, error: '//info, error)
     endif
     deallocate(my_x,my_b)

  endsubroutine LA_Matrix_Solve_Matrix

  subroutine qpotrs(factor,x, info)
    real(qp), dimension(:,:), intent(in) ::factor
    real(qp), dimension(:,:), intent(inout) :: x
    integer, intent(out), optional :: info

    integer :: n, m, i, j

    n = size(factor,1)
    m = size(x,2)
!$omp parallel do    
    do i = 1, m
       do j = 1, n-1
          x(j,i) = x(j,i)/factor(j,j)
          x(j+1:n,i) = x(j+1:n,i) - x(j,i)*factor(j+1:n,j)
       enddo
       x(n,i) = x(n,i)/factor(n,n)

       do j = n, 2, -1
          x(j,i) = x(j,i)/factor(j,j)
          x(1:j-1,i) = x(1:j-1,i) - x(j,i)*factor(1:j-1,j)
       enddo
       x(1,i) = x(1,i)/factor(1,1)
    enddo
    
    if( present(info) ) info = 0

  endsubroutine qpotrs

  function LA_Matrix_LogDet(this,error)

     type(LA_Matrix), intent(inout) :: this         
     integer, intent(out), optional :: error
     real(qp) :: LA_Matrix_LogDet
     integer :: i

     INIT_ERROR(error)

     if( this%factorised == NOT_FACTORISED ) then
        call LA_Matrix_Factorise(this,error=error)
     endif

     LA_Matrix_LogDet = 0.0_qp
     do i = 1, this%m
        LA_Matrix_LogDet = LA_Matrix_LogDet + log(abs(this%factor(i,i)))
     enddo

     if(this%equilibrated) LA_Matrix_LogDet = LA_Matrix_LogDet - sum( log(this%s) )

     if( this%factorised == CHOLESKY ) then
        LA_Matrix_LogDet = LA_Matrix_LogDet*2
     elseif( this%factorised == QR ) then
        LA_Matrix_LogDet = LA_Matrix_LogDet*1
     else
        RAISE_ERROR('LA_Matrix_LogDet: matrix not Cholesky-factorised or QR-factorised',error)
     endif

  endfunction LA_Matrix_LogDet

  function LA_Matrix_Det(this,error)

     type(LA_Matrix), intent(inout) :: this         
     integer, intent(out), optional :: error
     real(qp) :: LA_Matrix_Det
     integer :: i

     INIT_ERROR(error)

     if( this%factorised == NOT_FACTORISED ) then
        call LA_Matrix_Factorise(this,error=error)
     endif

     LA_Matrix_Det = 1.0_qp
     do i = 1, this%n
        LA_Matrix_Det = LA_Matrix_Det * this%factor(i,i)
     enddo

     if(this%equilibrated) LA_Matrix_Det = LA_Matrix_Det / product( this%s )

     if( this%factorised == CHOLESKY ) then
        LA_Matrix_Det = LA_Matrix_Det**2
     elseif( this%factorised == QR ) then
        LA_Matrix_Det = LA_Matrix_Det*1
     else
        RAISE_ERROR('LA_Matrix_LogDet: matrix not Cholesky-factorised or QR-factorised',error)
     endif

  endfunction LA_Matrix_Det

  subroutine LA_Matrix_QR_Factorise(this,q,r,error)
     type(LA_Matrix), intent(inout) :: this         
     real(qp), dimension(:,:), intent(out), optional :: q, r
     integer, intent(out), optional :: error

     integer :: info
#ifndef HAVE_QP
     real(dp), dimension(:), allocatable :: work
     integer :: lwork
#endif

     INIT_ERROR(error)

     this%factor = this%matrix

#if defined(HAVE_QP) || defined(ALBERT_LAPACK)
     call qgeqrf(this%factor,this%tau,info)
#else

     allocate(work(1))
     lwork = -1
     call dgeqrf(this%n, this%m, this%factor, this%n, this%tau, work, lwork, info)
     lwork = nint(work(1))
     deallocate(work)

     allocate(work(lwork))
     call dgeqrf(this%n, this%m, this%factor, this%n, this%tau, work, lwork, info)
     deallocate(work)
#endif

     if( info /= 0 ) then
        RAISE_ERROR('LA_Matrix_QR_Factorise: '//(-info)//'-th parameter had an illegal value.',error)
     endif

     this%factorised = QR

     if( present(q) .or. present(r) ) call LA_Matrix_GetQR(this,q,r,error)

  endsubroutine LA_Matrix_QR_Factorise

  subroutine LA_Matrix_GetQR(this,q,r,error)
     type(LA_Matrix), intent(inout) :: this         
     real(qp), dimension(:,:), intent(out), optional :: q, r
     integer, intent(out), optional :: error

     integer :: j, info
#ifndef HAVE_QP
     real(dp), dimension(:), allocatable :: work
     integer :: lwork
#endif

     INIT_ERROR(error)

     if( this%factorised /= QR ) then
        RAISE_ERROR('LA_Matrix_GetQR: not QR-factorised, call LA_Matrix_QR_Factorise first.',error)
     endif

     if(present(q)) then
        call check_size('q', q, (/this%n,this%m/),'LA_Matrix_GetQR',error)
        q = this%factor

#if defined(HAVE_QP) || defined(ALBERT_LAPACK)
        call qorgqr(this%factor,this%tau,q,info)
#else
        allocate(work(1))
        lwork = -1
        call dorgqr(this%n, this%m, this%m, q, this%n, this%tau, work, lwork, info)
        lwork = nint(work(1))
        deallocate(work)

        allocate(work(lwork))
        call dorgqr(this%n, this%m, this%m, q, this%n, this%tau, work, lwork, info)
        deallocate(work)
#endif
     endif

     if(present(r)) then
        call check_size('r', r, (/this%m,this%m/),'LA_Matrix_GetQR',error)
        r = this%factor(1:this%m,1:this%m)
        do j = 1, this%m - 1
           r(j+1:this%m,j) = 0.0_dp
        enddo
     endif

     if( info /= 0 ) then
        RAISE_ERROR('LA_Matrix_QR_Factorise: '//(info)//'-th parameter had an illegal value.',error)
     endif

  endsubroutine LA_Matrix_GetQR

  subroutine LA_Matrix_QR_Solve_Matrix(this,matrix,result,error)
     type(LA_Matrix), intent(inout) :: this
     real(qp), dimension(:,:), intent(in) :: matrix
     real(qp), dimension(:,:), intent(out) :: result
     integer, intent(out), optional :: error

     real(qp), dimension(:,:), allocatable :: my_result
     integer :: info, i, j, n, o
#ifndef HAVE_QP
     real(dp), dimension(:), allocatable :: work
     integer :: lwork
#endif

     INIT_ERROR(error)

     if(this%factorised == NOT_FACTORISED) then
        call LA_Matrix_QR_Factorise(this,error=error)
     elseif(this%factorised /= QR) then
        RAISE_ERROR('LA_Matrix_QR_Solve_Matrix: matrix not QR-factorised',error)
     endif

     n = size(matrix,1)
     o = size(matrix,2)
     call check_size('result', result, (/this%m,o/),'LA_Matrix_QR_Solve_Matrix',error)

     if( n /= this%n ) then
        RAISE_ERROR('LA_Matrix_QR_Solve_Matrix: dimensions of Q and matrix do not match.',error)
     endif

     allocate(my_result(n,o))
     my_result = matrix
#if defined(HAVE_QP) || defined(ALBERT_LAPACK)
     call qormqr(this%factor, this%tau, my_result, info)
#else
     lwork = -1
     allocate(work(1))
     call dormqr('L', 'T', this%n, o, this%m, this%factor, this%n, this%tau, my_result, this%n, work, lwork, info)
     lwork = nint(work(1))
     deallocate(work)

     allocate(work(lwork))
     call dormqr('L', 'T', this%n, o, this%m, this%factor, this%n, this%tau, my_result, this%n, work, lwork, info)
     deallocate(work)
#endif

     if( info /= 0 ) then
        RAISE_ERROR('LA_Matrix_QR_QR_Solve_Matrix: '//(-info)//'-th parameter had an illegal value.',error)
     endif

     do i = 1, o
        do j = this%m, 2, -1
           my_result(j,i) = my_result(j,i)/this%factor(j,j)
           my_result(1:j-1,i) = my_result(1:j-1,i) - my_result(j,i)*this%factor(1:j-1,j)
        enddo
        my_result(1,i) = my_result(1,i) / this%factor(1,1)
     enddo

     result = my_result(1:this%m,:)
     deallocate(my_result)

  endsubroutine LA_Matrix_QR_Solve_Matrix

  subroutine LA_Matrix_QR_Solve_Vector(this,vector,result,error)
     type(LA_Matrix), intent(inout) :: this
     real(qp), dimension(:), intent(in) :: vector
     real(qp), dimension(:), intent(out) :: result
     integer, intent(out), optional :: error

     real(qp), dimension(:,:), allocatable :: my_result
     integer :: n, m

     INIT_ERROR(error)

     n = size(vector)
     m = size(result)

     allocate(my_result(m,1))

     call LA_Matrix_QR_Solve_Matrix(this,reshape(vector,(/n,1/)),my_result,error=error)
     result = my_result(:,1)

     deallocate(my_result)

  endsubroutine LA_Matrix_QR_Solve_Vector

  subroutine LA_Matrix_QR_Inverse(this,inverse,error)
    type(LA_Matrix), intent(inout) :: this
    real(dp), dimension(:,:), intent(out) :: inverse
    integer, optional, intent(out) :: error

    real(dp), allocatable, dimension(:,:) :: Q, R, R_inv, Q_T
    integer i,j,k

    INIT_ERROR(error)

    if (this%n /= this%m) then
       RAISE_ERROR('LA_Matrix_QR_inverse: cannot invert non-square matrix',error)
    endif

    call check_size('inverse',inverse,shape(this%factor),'LA_Matrix_QR_Inverse',error=error)
    PASS_ERROR(error)
   
    allocate(Q(this%n, this%n), R(this%n, this%n), R_inv(this%n, this%n), Q_T(this%n, this%n))

    if( this%factorised == NOT_FACTORISED ) then
       call LA_Matrix_QR_Factorise(this, error=error)
    elseif( this%factorised /= QR ) then
       RAISE_ERROR('LA_Matrix_QR_Inverse: matrix not QR-factorised',error)
    endif

    call LA_Matrix_GetQR(this, Q, R)

    ! inversion of upper triangular matrix R
    R_inv(:,:) = 0.0_dp
    do j=1,this%n
       do i=1,j-1
          do k=1,j-1
             R_inv(i,j) = R_inv(i,j) + R_inv(i,k)*R(k,j)
          end do
       end do
       do k=1,j-1
          R_inv(k,j) = -R_inv(k,j)/R(j,j)
       end do
       R_inv(j,j) = 1.0_dp/R(j,j)
    end do

    ! A^-1 = (QR)^-1 = R^-1 Q^-1 = R^-1 Q^T
    Q_T = transpose(Q)
    inverse = matmul(R_inv, Q_T)

    deallocate(Q, R, R_inv, Q_T)

  end subroutine LA_Matrix_QR_inverse

  subroutine house_qp(x,v,beta)
     real(qp), dimension(:), intent(in) :: x
     real(qp), dimension(size(x)), intent(out) :: v
     real(qp), intent(out) :: beta

     integer :: n
     real(qp) :: sigma, mu

     n = size(x)
     sigma = dot_product(x(2:n),x(2:n))
     v = (/ 1.0_qp, x(2:n) /)

     if( sigma == 0.0_qp ) then
        beta = 0.0_qp
     else
        mu = sqrt(x(1)**2+sigma)
        if( x(1) <= 0.0_qp ) then
           v(1) = x(1) - mu
        else
           v(1) = -sigma / (x(1)+mu)
        endif
        beta = 2.0_qp * v(1)**2 / (sigma + v(1)**2)
        v = v / v(1)
     endif
  
  endsubroutine house_qp

  subroutine qgeqrf(a,beta,info)
     real(qp), dimension(:,:), intent(inout) :: a
     real(qp), dimension(size(a,2)), intent(out) :: beta
     integer, intent(out), optional :: info

     integer :: m, n, i, j
     real(qp), dimension(:), allocatable :: v, vTa, beta_v

     m = size(a,1)
     n = size(a,2)

     allocate(v(m), vTa(n), beta_v(m))

     do j = 1, n
        call house_qp(a(j:m,j),v(j:m),beta(j))
        vTa(j:n) = matmul(v(j:m),a(j:m,j:n))
        beta_v(j:m) = beta(j)*v(j:m)

        do i = j, n
           a(j:m,i) = a(j:m,i) - beta_v(j:m)*vTa(i)
        enddo
        if (j < m) a(j+1:m,j) = v(j+1:m) 
     enddo

     if( present(info) ) info = 0
     deallocate(v, vTa, beta_v)
     
  endsubroutine qgeqrf
  
  subroutine qormqr(a, beta, c, info)

     real(qp), dimension(:,:), intent(in) :: a
     real(qp), dimension(size(a,2)), intent(in) :: beta
     real(qp), dimension(:,:), intent(inout) :: c
     integer, intent(out), optional :: info

     integer :: m, n, q, i, j
     real(qp), dimension(:), allocatable :: v, vTc, beta_v

     m = size(a,1)
     n = size(a,2)
     q = size(c,2)

     allocate(v(m), vTc(q), beta_v(m))

     do j = 1, n
        v(j:m) = (/ 1.0_qp, a(j+1:m,j) /)
        vTc = matmul(v(j:m),c(j:m,:))
        beta_v(j:m) = beta(j)*v(j:m)

        do i = 1, q
           c(j:m,i) = c(j:m,i) - beta_v(j:m)*vTc(i)
        enddo
     enddo

     if( present(info) ) info = 0
     deallocate(v, vTc, beta_v)

  endsubroutine qormqr
  
  subroutine qorgqr(a,beta,q,info)

     real(qp), dimension(:,:), intent(in) :: a
     real(qp), dimension(size(a,2)), intent(in) :: beta
     real(qp), dimension(:,:), intent(out) :: q
     integer, intent(out), optional :: info

     integer :: m, n, i, j
     real(qp), dimension(:), allocatable :: v, vTq, beta_v

     m = size(a,1)
     n = size(a,2)

     allocate(v(m), vTq(n), beta_v(m))

     q = 0.0_qp
     do j = 1, n
        q(j,j) = 1.0_qp
     enddo

     do j = n, 1, -1
        v(j:m) = (/1.0_qp, a(j+1:m,j) /)
        vTq(j:n) = matmul(v(j:m), q(j:m,j:n))
        beta_v(j:m) = beta(j)*v(j:m)

        do i = j, n
           q(j:m,i) = q(j:m,i) - beta_v(j:m)*vTq(i)
        enddo
     enddo

     if( present(info) ) info = 0
     deallocate(v, vTq, beta_v)

  endsubroutine qorgqr
  
  subroutine Matrix_CholFactorise(A,A_factor)

     real(dp), dimension(:,:), intent(inout), target :: A
     real(dp), dimension(:,:), intent(out), target, optional :: A_factor

     real(dp), dimension(:,:), pointer :: my_A_factor
     integer :: i, j, n, info

     if( .not. is_square(A) ) call system_abort('Matrix_CholFactorise: A is not square')
     if( present(A_factor) ) call check_size('A_factor',A_factor,shape(A), 'Matrix_CholFactorise')

     n = size(A,1)

     if( present(A_factor) ) then
        my_A_factor => A_factor
        my_A_factor = A
     else
        my_A_factor => A
     endif

     call dpotrf('U', n, my_A_factor, n, info)

     if( info /= 0 ) call system_abort('Matrix_CholFactorise: cannot factorise, PRINT_ALWAYS: '//info)

     do i = 1, n
        do j = i+1, n
           my_A_factor(j,i) = 0.0_dp
        enddo
     enddo

     my_A_factor => null()

  endsubroutine Matrix_CholFactorise

  ! Solve system of linear equations with an already Cholesky-factorised matrix.
  subroutine Matrix_BackSubstitute(U,x,b)

     real(dp), dimension(:,:), intent(in) :: U
     real(dp), dimension(:), intent(out)  :: x
     real(dp), dimension(:), intent(in)   :: b

     integer :: n, info

     n = size(U,1)

     x = b

     call dpotrs('U', n, 1, U, n, x, n, info)

     if( info /= 0 ) call system_abort('Matrix_BackSubstitute: cannot solve system')

  endsubroutine Matrix_BackSubstitute

  ! Solve system of linear equations with an upper triagonal matrix.
  subroutine Matrix_Solve_Upper_Triangular(U,x,b)

     real(dp), dimension(:,:), intent(in) :: U
     real(dp), dimension(:), intent(out)  :: x
     real(dp), dimension(:), intent(in)   :: b

     integer :: n, info

     n = size(U,1)

     x = b

     call dtrtrs('U','T','N', n, 1, U, n, x, n, info)

     if( info /= 0 ) call system_abort('Matrix_Solve_Upper_Triangular: cannot solve system')

  endsubroutine Matrix_Solve_Upper_Triangular

  ! Invert an already Cholesky-factorised matrix.
  subroutine Matrix_Factorised_Inverse(U,inverse_A)

     real(dp), dimension(:,:), target, intent(inout)         :: U
     real(dp), dimension(:,:), target, intent(out), optional :: inverse_A

     real(dp), dimension(:,:), pointer :: inverse
     integer :: i, j, n, info

     n = size(U,1)

     if( present(inverse_A) ) then
        inverse => inverse_A
        inverse = U
     else
        inverse => U
     endif

     call dpotri('U', n, inverse, n, info)

     if( info /= 0 ) call system_abort('Matrix_Factorised_Inverse: cannot invert')

     do i = 1, n
        do j = i+1, n
           inverse(j,i) = inverse(i,j)
        enddo
     enddo
     inverse => null()

  endsubroutine Matrix_Factorised_Inverse

  !% Test if this matrix is diagonal
  function is_diagonal(matrix)
    real(dp),intent(in), dimension(:,:) :: matrix
    logical :: is_diagonal

    integer :: i,j

    if (.not. is_square(matrix)) call system_abort('Matrix_diagonal: matrix not squared')
     
    do i=1,size(matrix,1)
       do j=1,size(matrix,2)
          if (i == j) cycle
          if (matrix(i,j) .fne. 0.0_dp) then
             is_diagonal = .false.
             return
          end if
       end do
    end do
    
    is_diagonal = .true.

  end function is_diagonal

  ! matrix_is_symmetric(matrix)
  !
  ! tells if tha matrix is symmetric
  function matrix_is_symmetric(matrix) result(symm)
    real(dp),intent(in), dimension(:,:):: matrix
    logical::symm
    integer::i,j,N

    real(dp) :: maxv
    
    if (.not.Is_Square(matrix)) &
         call system_abort('Matrix_Is_Symmetric: Matrix is not square')

    maxv = maxval(abs(matrix))
    N=size(matrix,1)
    symm=.true.
    do i=1,N
       do j=i+1,N
          if (abs(matrix(i,j)-matrix(j,i)) > NUMERICAL_ZERO*maxv) then
             symm=.false.
             return
          end if
       enddo
    enddo

  end function matrix_is_symmetric

  ! matrix_is_symmetric(matrix)
  !
  ! tells if the matrix is symmetric
  function int_matrix_is_symmetric(matrix) result(symm)

    integer, intent(in), dimension(:,:):: matrix
    logical::symm
    integer::i,j,N
    
    if (.not.is_square(matrix)) &
         call system_abort('is_symmetric: Matrix is not square')

    N=size(matrix,1)
    symm=.true.
    do i=1,N
       do j=i+1,N
          if (matrix(j,i) /= matrix(i,j)) then
             symm=.false.
             return
          end if
       enddo
    enddo

  end function int_matrix_is_symmetric

  ! matrix_z_is_symmetric(matrix)
  !
  ! tells if the matrix is symmetric
  function matrix_z_is_symmetric(matrix) result(symm)
    complex(dp),intent(in), dimension(:,:):: matrix
    logical::symm
    integer::i,j,N

    real(dp) :: maxv
    
    if (.not.Is_Square(matrix)) &
         call system_abort('Matrix_Is_Symmetric: Matrix is not square')

    maxv = maxval(abs(matrix))

    N=size(matrix,1)
    symm=.true.
    do i=1,N
       do j=i+1,N
          if (abs(matrix(j,i)-matrix(i,j)) > maxv*NUMERICAL_ZERO ) then
             symm=.false.
             return
          end if
       enddo
    enddo

  end function matrix_z_is_symmetric

  ! matrix_z_is_hermitian(matrix)
  !
  ! tells if tha matrix is hermitian
  function matrix_z_is_hermitian(matrix) result(herm)
    complex(dp),intent(in), dimension(:,:):: matrix

    logical::herm
    integer::i,j,N

    real(dp) :: maxv
    
    if (.not.Is_Square(matrix)) &
         call system_abort('Matrix_Is_Symmetric: Matrix is not square')

    maxv = maxval(abs(matrix))
    N=size(matrix,1)
    herm=.true.
    do i=1,N
       do j=i,N
          if (abs(matrix(j,i)-conjg(matrix(i,j))) > maxv*NUMERICAL_ZERO) then
             herm=.false.
             return
          end if
       enddo
    enddo

  end function matrix_z_is_hermitian

  function matrix_is_orthogonal(matrix)
    real(dp),intent(in), dimension(:,:):: matrix
    logical :: matrix_is_orthogonal
    real(dp) :: t(size(matrix,1), size(matrix,2)), identity(size(matrix,1),size(matrix,2))
    
    if (.not.Is_Square(matrix)) &
         call system_abort('Matrix_Is_Symmetric: Matrix is not square')

    t = matrix .mult. transpose(matrix)
    identity = 0.0_dp
    call add_identity(identity)
    matrix_is_orthogonal = maxval(abs(t - identity)) < maxval(matrix)*NUMERICAL_ZERO
    
  end function matrix_is_orthogonal

  ! print(matrix)
  ! printing

  subroutine matrix_print_mainlog(this, verbosity)
    real(dp),intent(in),dimension(:,:) :: this
    integer, optional::verbosity
    call matrix_print(this,verbosity,mainlog)
  end subroutine matrix_print_mainlog

  subroutine matrix_z_print_mainlog(this, verbosity)
    complex(dp),intent(in),dimension(:,:) :: this
    integer, optional::verbosity
    call matrix_z_print(this,verbosity,mainlog)
  end subroutine matrix_z_print_mainlog

  subroutine int_matrix_print_mainlog(this, verbosity)
    integer,intent(in),dimension(:,:) :: this
    integer, optional::verbosity
    call int_matrix_print(this,verbosity,mainlog)
  end subroutine int_matrix_print_mainlog

  subroutine int_matrix_print(this,verbosity,file)

    integer, dimension(:,:), intent(in) :: this
    integer, optional,       intent(in) :: verbosity
    type(inoutput),          intent(in) :: file

    integer :: i,N,M,w
    character(20) :: format
    
    w = int(log10(real(maxval(abs(this))+1,dp)))+3 ! leave at least one space between entries, including room for neg. sign
    N = size(this,1)
    M = size(this,2)
    write(format,'(2(a,i0),a)')'(',M,'i',w,')'

    do i = 1, N
       write(line,format) this(i,:)
       call print(line,verbosity,file)
    end do

  end subroutine int_matrix_print

  subroutine matrix_print(this,verbosity,file, always)
    type(inoutput), intent(in)                 :: file
    real(dp),       intent(in),dimension(:,:)  :: this
    integer, optional, intent(in)              :: verbosity
    logical, optional, intent(in)              :: always

    integer                                    :: i, n, w, j
    logical                                    :: t
    character(20)                              :: format
    integer                                    :: max_size = 5
    integer                                    :: absolute_max_width = 50
    logical :: do_always

    do_always = optional_default(.false., always)

    if (do_always) then
      do i=1, size(this,2)
	do j=1, size(this,1)
	  call print(i // " " // j // " " // this(i,j), verbosity, file)
	end do
      end do
    else

      if (size(this,2) > max_size .and. size(this,1) <= max_size) then
	 w = size(this,1)
	 n = size(this,2)
	 t = .true.
      else
	 w = size(this,2)
	 n = size(this,1)
	 t = .false.
      end if

      if (w > absolute_max_width) then
	 call print('Matrix_print: matrix is too large to print', verbosity, file)
	 return
      end if

      if (t) then
	write (line, '(a)') 'Matrix_Print: printing matrix transpose'
	call print(line,verbosity,file)
      endif

      write(format,'(a,i0,a)')'(',w,'(1x,f18.10))'

      do i=1,n
	 if (t) then
	    write(line,format) this(:,i)
	 else
	    write(line,format) this(i,:)
	 end if
	 call print(line,verbosity,file)
      end do

    end if

  end subroutine matrix_print

  subroutine matrix_z_print(this,verbosity,file)
    type(inoutput), intent(in)                  :: file
    complex(dp),    intent(in), dimension(:,:)  :: this
    integer, optional                           :: verbosity
    integer                                     :: i, n, w
    logical                                     :: t
    character(200)                              :: format
    integer                                     :: max_size = 6

    if (size(this,2) > max_size .and. size(this,1) <= max_size) then
       w = size(this,1)
       n = size(this,2)
       t = .true.
    else
       w = size(this,2)
       n = size(this,1)
       t = .false.
    end if


!    if (w > 60) then
!       write(line,'(a,i0,a,i0,a)')'Matrix_z_Print: Both dimensions (',w,',',n,') are too large to print'
!       call system_abort(line)
!    end if

    if (t) then
      write(line, '(a)') 'Matrix_z_Print: printing transpose'
      call print(line,verbosity,file)
    endif

    write(format,'(a,i0,a)')'(',w,'(x,f12.6,"+I*",f12.6))'

    do i=1,n
       if (t) then
          write(line,format) this(:,i)
       else
          write(line,format) this(i,:)
       end if
       call print(line,verbosity,file)
    end do

  end subroutine matrix_z_print

  subroutine matrix_print_mathematica(this, verbosity, file)
    real(dp),    intent(in), dimension(:,:)  :: this
    integer, optional                           :: verbosity
    type(inoutput), intent(in), optional        :: file

    integer i, j

    call print("M = { ", verbosity, file)
    do i=1, size(this,1)
      call print ("{", verbosity, file)
      do j=1, size(this,1)
	if (j == size(this,1)) then
	  call print(this(i,j), verbosity, file)
	else
	  call print(this(i,j)//", ", verbosity, file)
	endif
      end do
      if (i == size(this,1)) then
	call print ("}", verbosity, file)
      else
	call print ("},", verbosity, file)
      endif
    end do
    call print ("};", verbosity, file)
  end subroutine matrix_print_mathematica

  subroutine matrix_z_print_mathematica(this, verbosity, file)
    complex(dp),    intent(in), dimension(:,:)  :: this
    integer, optional                           :: verbosity
    type(inoutput), intent(in), optional        :: file

    integer i, j

    call print("M = { ", verbosity, file)
    do i=1, size(this,1)
      call print ("{", verbosity, file)
      do j=1, size(this,1)
	if (j == size(this,1)) then
	  call print(real(this(i,j), dp) // " + I*"//aimag(this(i,j)), verbosity, file)
	else
	  call print(real(this(i,j), dp) // " + I*"//aimag(this(i,j))//", ", verbosity, file)
	endif
      end do
      if (i == size(this,1)) then
	call print ("}", verbosity, file)
      else
	call print ("},", verbosity, file)
      endif
    end do
    call print ("};", verbosity, file)
  end subroutine matrix_z_print_mathematica

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X
!X Polynomial fitting to boundary conditions
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !
  !% Fit a cubic polynomial to a function given its values ($y_0$ and $y_1$) and 
  !% its derivative ($y^\prime_0$ and $y^\prime_1$) at two points ($x_0$ and $x_1$)
  !
  !% \begin{displaymath}
  !% y = ax^3 + bx^2 + cx + d
  !% \end{displaymath}
  !
  subroutine fit_cubic(x0,y0,y0p,x1,y1,y1p,coeffs)

    real(dp), intent(in)     :: x0,y0,y0p !% $x_0$, $y_0$, $y^\prime_0$
    real(dp), intent(in)     :: x1,y1,y1p !% $x_1$, $y_1$, $y^\prime_1$
    real(dp), intent(out)    :: coeffs(4) !% '(/a,b,c,d/)'
    real(dp), dimension(4,4) :: A, invA
    real(dp)                 :: x02, x03, x12, x13 !x0^2, x0^3 ...
    
    x02 = x0*x0; x03 = x02*x0
    x12 = x1*x1; x13 = x12*x1

    A = reshape( (/    x03,    x13, 3.0_dp*x02, 3.0_dp*x12,           &
                       x02,    x12,  2.0_dp*x0,  2.0_dp*x1,           &
                        x0,     x1,     1.0_dp,     1.0_dp,           &
                    1.0_dp, 1.0_dp,     0.0_dp,     0.0_dp /),(/4,4/) )

    call matrix_inverse(A,invA)

    coeffs = invA .mult. (/y0,y1,y0p,y1p/)

  end subroutine fit_cubic

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X
!X VECTOR stuff
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  ! normsq()
  ! returns (X.dot.X)
  pure function vector_normsq(vector) result(normsq) 

    real(dp), intent(in), dimension(:) :: vector
    real(dp)             :: normsq
   
    normsq = dot_product(vector,vector)

  end function vector_normsq

  pure function vector_normsq_q(vector) result(normsq) 

    real(qp), intent(in), dimension(:) :: vector
    real(qp)             :: normsq
   
    normsq = dot_product(vector,vector)

  end function vector_normsq_q

  ! norm()
  ! returns SQRT((X.dot.X))
  pure function vector_norm(vector) result(norm)

    real(dp), intent(in),dimension(:) :: vector
    real(dp)             :: norm

    norm = sqrt(dot_product(vector,vector))

  end function vector_norm


  ! x .dot. y
  pure function vector_dotpr(this,other) result(dotpr)

    real(dp), intent(in), dimension(:) :: this,other
    real(dp)             :: dotpr
   
    dotpr = dot_product(this,other)
 
  end function vector_dotpr


  !% Return the (smallest) angle between two vectors, in radians.
  !% This is calculated as
  !%\begin{displaymath}
  !%\arccos\left(\frac{\mathbf{a}\cdot\mathbf{b}}{|\mathbf{a}| |\mathbf{b}|}\right)
  !%\end{displaymath}
  pure function angle(a,b)

    real(dp), dimension(:), intent(in) :: a,b
    real(dp)                           :: angle
    real(dp)                           :: arg

    arg = (a .dot. b)/(norm(a)*norm(b))
    if (arg > 1.0_dp) then
       arg = 1.0_dp
    else if (arg < -1.0_dp) then
       arg = -1.0_dp
    end if

    angle = acos(arg)

  end function angle

  ! x .outer. y 
  pure function outer(vector1,vector2) result(outr)
    real(dp),intent(in), dimension(:) ::vector1,vector2
    real(dp), dimension(size(vector1),size(vector2)) ::outr
    integer::i,j
     
    do j=1,size(vector2)
       do i=1,size(vector1)
          outr(i,j)=vector1(i)*vector2(j)
       end do
    end do
  
  end function outer

  pure function outer_qq(vector1,vector2) result(outr)
    real(qp),intent(in), dimension(:) ::vector1,vector2
    real(qp), dimension(size(vector1),size(vector2)) ::outr
    integer::i,j
     
    do j=1,size(vector2)
       do i=1,size(vector1)
          outr(i,j)=vector1(i)*vector2(j)
       end do
    end do
  
  end function outer_qq

  ! x .outer. y 
  pure function d_outer_zz(vector1,vector2) result(outr)
    complex(dp),intent(in), dimension(:) ::vector1,vector2
    real(dp), dimension(size(vector1),size(vector2)) ::outr
    integer::i,j
     
    do j=1,size(vector2)
       do i=1,size(vector1)
          outr(i,j)=vector1(i)*conjg(vector2(j))
       end do
    end do
  
  end function d_outer_zz

  ! x .outer. y 
  pure function z_outer_zz(vector1,vector2) result(outr)
    complex(dp),intent(in), dimension(:) ::vector1,vector2
    complex(dp), dimension(size(vector1),size(vector2)) :: outr
    integer::i,j
     
    do j=1,size(vector2)
       do i=1,size(vector1)
          outr(i,j)=vector1(i)*conjg(vector2(j))
       end do
    end do
  
  end function z_outer_zz


  !% Return the square root of 'r' if it is positive, or zero otherwise. 
  elemental function sqrt_cut(r) 
    real(dp),intent(in)::r
    real(dp)::sqrt_cut
    sqrt_cut = sqrt(max(0.0_dp,r))
  end function sqrt_cut


  ! v1 .feq. v2
  !
  !floating point logical comparison
  function vector_feq(vector1,vector2) 
    real(dp),intent(in), dimension(:) ::vector1,vector2
    logical::vector_feq
    integer::i
   
    if(size(vector1) /= size(vector2)) &
         call check_size('Vector2',vector2,size(vector1),'Vector_FEQ')
  
    vector_feq=.true.
    do i=1,size(vector1)
       if (vector1(i).fne.vector2(i)) then 
          vector_feq=.false.
          return
       end if
    end do

  end function vector_feq

  ! v1 .fne. v2
  !
  function vector_fne(vector1,vector2) 
    real(dp),intent(in), dimension(:) :: vector1,vector2
    logical::vector_fne
    integer::i
   
    if(size(vector1) /= size(vector2)) &
         call check_size('Vector2',vector2,size(vector1),'Vector_FNE')

    vector_fne=.false.
    do i=1,size(vector1)
       if (vector1(i).fne.vector2(i)) then 
          vector_fne=.true.
          return
       end if
    end do
  
  end function vector_fne


  !printing  
  subroutine vector_print_mainlog(this, verbosity)
    real(dp),intent(in), dimension(:) :: this
    integer, optional::verbosity
    call vector_print(this,verbosity,mainlog)
  end subroutine vector_print_mainlog

  subroutine vector_print(this,verbosity,file)
    real(dp),intent(in), dimension(:) :: this
    type(inoutput)                    :: file
    integer, optional                 :: verbosity
    integer                           :: nbloc,i,nrest
   
    nbloc=size(this)/5
    nrest=mod(size(this),5)

    do i=1,nbloc
       call print(""//this((i-1)*5+1:(i-1)*5+5), verbosity, file)
    end do
    if (nrest /= 0) then
       call print(""//this((i-1)*5+1:(i-1)*5+nrest),verbosity,file)
    end if

  end subroutine vector_print

  subroutine vector_z_print_mainlog(this,verbosity)
    complex(dp),intent(in), dimension(:) ::this
    integer, optional::verbosity
    call vector_z_print(this,verbosity,mainlog)
  end subroutine vector_z_print_mainlog

  subroutine vector_z_print(this,verbosity,file)
    complex(dp),intent(in), dimension(:) :: this
    type(inoutput)                       :: file
    integer, optional::verbosity
    integer:: nbloc,i,nrest
   
    nbloc=size(this)/5
    nrest=mod(size(this),5)

    do i=1,nbloc
       call print(""//this((i-1)*5+1:(i-1)*5+5), verbosity, file)
    end do
    if (nrest /= 0) then
       call print(""//this((i-1)*5+1:(i-1)*5+nrest), verbosity, file)
    end if

  end subroutine vector_z_print

  subroutine integer_array_print_mainlog(this, verbosity)
    integer,intent(in), dimension(:) ::this
    integer, optional::verbosity
    call integer_array_print(this,verbosity,mainlog)
  end subroutine integer_array_print_mainlog

  subroutine integer_array_print(this,verbosity,file)
    integer,intent(in), dimension(:) ::this
    type(inoutput)::file
    integer, optional::verbosity
    integer:: nbloc,i,nrest
    character(20)::form

    nbloc=size(this)/5
    nrest=mod(size(this),5)

    do i=1,nbloc
       write(line,'(5i12)')this((i-1)*5+1:(i-1)*5+5)
       call print(line,verbosity,file)
    end do
    if (nrest /= 0) then
       write(form,'(a,i0,a)') '(',nrest,'i12)'
       write(line,form)this((i-1)*5+1:(i-1)*5+nrest)
       call print(line,verbosity,file)
    end if

  end subroutine integer_array_print

  subroutine logical_array_print_mainlog(this, verbosity)
    logical,intent(in), dimension(:) ::this
    integer, optional::verbosity
    call logical_array_print(this,verbosity,mainlog)
  end subroutine logical_array_print_mainlog

  subroutine logical_array_print(this,verbosity,file)
    logical,intent(in), dimension(:) ::this
    type(inoutput)::file
    integer, optional::verbosity
    integer:: nbloc,i,nrest
    character(20)::form

    nbloc=size(this)/5
    nrest=mod(size(this),5)

    do i=1,nbloc
       write(line,'(5l2)')this((i-1)*5+1:(i-1)*5+5)
       call print(line,verbosity,file)
    end do
    if (nrest /= 0) then
       write(form,'(a,i0,a)') '(',nrest,'l2)'
       write(line,form)this((i-1)*5+1:(i-1)*5+nrest)
       call print(line,verbosity,file)
    end if

  end subroutine logical_array_print


  ! add to vector elements, random quantities in the range(-a/2,a/2)     
  subroutine vector_randomise(v,a)
    real(dp), intent(inout), dimension(:)  :: v
    real(dp), intent(in)    :: a

    integer                 :: i
    do i=1,size(v)
       v(i)=v(i)+(ran_uniform()-0.5_dp)*a
    end do
  end subroutine vector_randomise

  ! add to complex vector elements, random quantities in the range(-a/2,a/2)     
  subroutine vector_z_randomise(v,a)
    complex(dp), intent(inout), dimension(:)  :: v
    real(dp), intent(in)    :: a

    integer                 :: i
    do i=1,size(v)
       v(i)=v(i)+cmplx((ran_uniform()-0.5_dp)*a/sqrt(2.0_dp), (ran_uniform()-0.5_dp)*a/sqrt(2.0_dp))
    end do
  end subroutine vector_z_randomise

  !
  !% Return a three vector with normally distributed components
  !
  function ran_normal3()
    real(dp), dimension(3) :: ran_normal3
    ran_normal3(1) = ran_normal()
    ran_normal3(2) = ran_normal()
    ran_normal3(3) = ran_normal()
  end function ran_normal3


  ! returns a vector, contining a histogram of frequencies
  function vector_histogram(vector,min_x,max_x,Nbin,weight_vector, drop_outside)

    real(dp), dimension(:), intent(in) :: vector
    real(dp),               intent(in) :: min_x, max_x
    integer,                intent(in) :: Nbin
    real(dp), dimension(:), intent(in), optional :: weight_vector
    logical, optional, intent(in) :: drop_outside
    real(dp), dimension(Nbin)          :: vector_histogram
    !local variables
    real(dp)                           :: binsize
    integer                            :: i, bin
    logical :: do_drop_outside
  
    if(max_x <= min_x) then
       call system_abort('Vector_Histogram: max_x < min_x')
    end if

    do_drop_outside = optional_default(.false., drop_outside)

    binsize=(max_x-min_x)/(real(Nbin,dp))
    vector_histogram = 0.0_dp

    do i=1,size(vector)

       bin = ceiling((vector(i)-min_x)/binsize)
       if (do_drop_outside) then
	 if (bin < 1 .or. bin > Nbin) cycle
       else
	 if (bin < 1) bin = 1
	 if (bin > Nbin) bin = Nbin
       endif
       if (present(weight_vector)) then
	 vector_histogram(bin) = vector_histogram(bin) + weight_vector(i)
       else
	 vector_histogram(bin) = vector_histogram(bin) + 1.0_dp
       endif

    end do

  end function vector_histogram

  ! Root-Mean-Square difference between elements of two vectors

  !% For two vectors of $N$ dimensions, this is calculated as
  !%\begin{displaymath}
  !%\sum_{i=1}^{N} \left(\mathbf{v_1}_i - \mathbf{v_2}_i\right)^2
  !%\end{displaymath}
  function rms_diff1(vector1, vector2)
    real(dp), dimension(:), intent(in) :: vector1, vector2
    real(dp)                           :: rms_diff1, d
    integer                            :: i
    if(size(vector1) /= size(vector2)) &
         call check_size('Vector 2',vector2,shape(vector1),'RMS_diff')
    rms_diff1 = 0.0_dp
    do i = 1, size(vector1)
       d = vector1(i) - vector2(i)
       rms_diff1 = rms_diff1 + d * d
    end do
    rms_diff1 = rms_diff1 / real(size(vector1),dp)
    rms_diff1 = sqrt(rms_diff1)
  end function rms_diff1


  !% OMIT
  function array_dotprod(this,other, dir) result(c)
    real(dp),intent(in), dimension(:,:) :: this,other
    integer, intent(in)::dir
    real(dp)::c(size(this,3-dir))
    integer::i


    call check_size('this, other',this,shape(other),'Array_Dotprod')

    if(dir==1) then
       do i=1,size(this,2)
          c(i)=sqrt(dot_product(this(:,i),other(:,i)))
       end do
    else if(dir==2) then
       do i=1,size(this,1)
          c(i)=sqrt(dot_product(this(i,:),other(i,:)))
       end do
    else
       call system_abort('array_dotprod: dir must be 1 or 2')
    end if
    
  end function array_dotprod


  ! normsq of each vector
  function array_normsq(this, dir)
    real(dp), dimension(:,:) ::this
    integer, intent(in)::dir
    real(dp)::array_normsq(size(this,3-dir))
    integer::i

    if(dir==1) then
       do i=1,size(this,2)
          array_normsq(i)=dot_product(this(:,i),this(:,i))
       end do
    else if(dir==2) then
       do i=1,size(this,1)
          array_normsq(i)=dot_product(this(i,:),this(i,:))
       end do
    else
       call system_abort('array_normsq: dir must be 1 or 2')
    end if 
  end function array_normsq

  ! norm for each vector
  function array_norm(this, dir) result(sqvalue)
    real(dp),intent(in), dimension(:,:) :: this
    integer, intent(in)::dir
    real(dp)::sqvalue(size(this,3-dir))
    integer::i

    if(dir==1) then
       do i=1,size(this,2)
          sqvalue(i)=sqrt(dot_product(this(:,i),this(:,i)))
       end do
    else if(dir==2) then
       do i=1,size(this,1)
          sqvalue(i)=sqrt(dot_product(this(i,:),this(i,:)))
       end do
    else
       call system_abort('array_norm: dir must be 1 or 2')
    end if
  end function array_norm


  !% Return a list of scalar triple products for each triplet of vectors from 'a', 'b' and 'c', i.e.
  !%> array3_triple(i) = a(:,i) .dot. (b(:,i) .cross. c(:,i))'.
  function array3_triple(a,b,c)
    real(dp),intent(in), dimension(:,:) ::a, b, c
    real(dp)::array3_triple(size(a,2))
    integer::n,m,i
    n=size(a,1)
    m=size(a,2)

    call check_size('B',b,shape(a),'Array3_Triple')
    call check_size('C',c,shape(a),'Array3_Triple')

    do i = 1,m
       array3_triple(i) = scalar_triple_product(a(:,i),b(:,i),c(:,i))
    end do

  end function array3_triple


  ! add random values to the elements of an array in the range(-a/2,a/2)
  subroutine matrix_randomise(m,a)
    real(dp), dimension(:,:), intent(inout) ::m
    real(dp)::a
    !local
    integer::i,j
   
    do i=1,size(m,1)
       do j=1,size(m,2)
          m(i,j)=m(i,j)+(ran_uniform()-0.5_dp)*a
       end do
    end do
  
  end subroutine matrix_randomise

  subroutine matrix_z_randomise(m,a)
    complex(dp), dimension(:,:), intent(inout) ::m
    real(dp)::a
    !local
    integer::i,j
   
    do i=1,size(m,1)
       do j=1,size(m,2)
          m(i,j)=m(i,j)+cmplx((ran_uniform()-0.5_dp)*a/sqrt(2.0_dp), (ran_uniform()-0.5_dp)*a/sqrt(2.0_dp))
       end do
    end do
  
  end subroutine matrix_z_randomise

  ! add to matrix elements, random quantities in the range(-a/2,a/2)     
  subroutine matrix_randomise_vweight(m,a)
    real(dp), intent(inout), dimension(:,:)  :: m
    real(dp), intent(in)    :: a(:)

    integer                 :: i,j

    if (size(m,2) /= size(a)) &
      call system_abort("matrix_randomise_vweight incompatible sizes : m " // shape(m) // " a " // shape(a))
    do i=1,size(m,1)
      do j=1,size(m,2)
        m(i,j)=m(i,j)+(ran_uniform()-0.5_dp)*a(j)
      end do
    end do
  end subroutine matrix_randomise_vweight

  !% For two arrays of dimension $N \times M$, this is calculated as
  !%\begin{displaymath}
  !%\sum_{j=1}^{M} \sum_{i=1}^{N} \left(\mathbf{a_1}_{ij} - \mathbf{a_2}_{ij}\right)^2
  !%\end{displaymath}
  function rms_diff2(array1, array2)
    real(dp), dimension(:,:), intent(in) :: array1, array2
    real(dp)                             :: rms_diff2, d
    integer                              :: i,j
    call check_size('Array 2',array2,shape(array1),'rms_diff')
    rms_diff2 = 0.0_dp
    do j = 1, size(array1,2)
       do i = 1, size(array1,1)
          d = array1(i,j) - array2(i,j)
          rms_diff2 = rms_diff2 + d * d
       end do
    end do
    rms_diff2 = rms_diff2 / real(size(array1),dp)
    rms_diff2 = sqrt(rms_diff2)
  end function rms_diff2


  ! ---------------------------------
  !
  ! 3-vector procedures
  !

  pure function cross_product(x,y) ! x ^ y

    real(dp), dimension(3), intent(in):: x,y
    real(dp), dimension(3)            :: cross_product

    cross_product(1) = + ( x(2)*y(3) - x(3)*y(2) )
    cross_product(2) = - ( x(1)*y(3) - x(3)*y(1) )
    cross_product(3) = + ( x(1)*y(2) - x(2)*y(1) )

  end function cross_product

  !% Return the scalar triple product $\mathbf{x} \cdot \mathbf{y} \times \mathbf{z}$
  !% of the 3-vectors 'x', 'y' and 'z'.
  pure function scalar_triple_product(x,y,z)  ! [x,y,z]

    real(dp), dimension(3), intent(in) :: x,y,z
    real(dp)                           :: scalar_triple_product

    scalar_triple_product = + x(1) * ( y(2)*z(3) - y(3)*z(2) )    &
         - x(2) * ( y(1)*z(3) - y(3)*z(1) )    &
         + x(3) * ( y(1)*z(2) - y(2)*z(1) )

  end function scalar_triple_product

  !% Return the vector triple product $\mathbf{x} \times (\mathbf{y} \times \mathbf{z})$
  !% of the 3-vectors 'x', 'y' and 'z'.
  pure function vector_triple_product(x,y,z) ! x ^ (y ^ z) = y(x.z) - z(x.y)

    real(dp), dimension(3), intent(in) :: x,y,z
    real(dp), dimension(3)             :: vector_triple_product

    vector_triple_product = y * (x.dot.z) - z * (x.dot.y)

  end function vector_triple_product

  !
  !% Return a unit vector in the direction given by the angles 'theta' and 'phi', i.e.
  !% convert the spherical polar coordinates point $(1,\theta,\phi)$ to cartesians.
  !
  pure function unit_vector(theta,phi)

    real(dp),   intent(in) :: theta, phi
    real(dp), dimension(3) :: unit_vector

    unit_vector(1) = sin(theta)*cos(phi)
    unit_vector(2) = sin(theta)*sin(phi)
    unit_vector(3) = cos(theta)

  end function unit_vector

  !
  !% Returns a random unit vector which is
  !% uniformly distributed over the unit sphere
  !% [See Knop, CACM 13 326 (1970)].
  !
  function random_unit_vector() result(v)
    real(dp), dimension(3) :: v
    real(dp)               :: x,y,z,s

    !Find a point in the unit circle
    s = 1.1_dp
    do while(s > 1.0_dp)
       x = 2.0_dp * ran_uniform() - 1.0_dp
       y = 2.0_dp * ran_uniform() - 1.0_dp
       s = x*x + y*y
    end do

    !z must be uniform in [-1,1]. s is uniform in [0,1]
    z = 2.0_dp * s - 1.0_dp
    
    !x and y are constrained by normalisation
    s = sqrt((1-z*z)/s)
    x = x*s
    y = y*s

    v = (/x,y,z/)

  end function random_unit_vector
  

  ! INTERPOLATION

  !% Linearly interpolate between the points $(x_0,y_0)$ and $(x_1,y_1)$.
  !% Returns the interpolated $y$ at position $x$, where $x_0 \le x \le x_1$.
  pure function linear_interpolate(x0,y0,x1,y1,x) result(y)
    real(dp), intent(in) :: x0,y0,x1,y1,x
    real(dp)             :: y
    y = y0 + (y1-y0)*(x-x0)/(x1-x0)
  end function linear_interpolate

  !% Perform a cubic interpolation between the points $(x_0,y_0)$ and $(x_1,y_1)$,
  !% using the cubic function with zero first derivative at $x_0$ and $x_1$.
  !% Returns the interpolated $y$ at position $x$, where $x_0 \le x \le x_1$.
  pure function cubic_interpolate(x0,y0,x1,y1,x) result(y)
    real(dp), intent(in) :: x0,y0,x1,y1,x
    real(dp)             :: u,v,w,y
    u = (x-x0)/(x1-x0) ! relative co-ordinate
    v = u * u
    w = 3*v - 2*v*u    ! w = 3u^2-2u^3
    y = y0 + (y1-y0)*w
  end function cubic_interpolate



  !%OMIT
  subroutine matrix_test()
    real(dp),allocatable ::a(:,:),b(:,:),c(:,:),vect1(:),vect2(:),d(:,:)
    real(dp)::tmp
    write(line,*) 'A is a 2x2 matrix:'; call print(line)
    allocate(a(2,2))
    !b=a assigment dosnt work !make subroutine
    a=3
    call print(a)
    write(line,*) 'Vect1 is a vector of 3 elements:'; call print(line)
    allocate(vect1(3),vect2(3))
    vect1=3
    call print(vect1)

    allocate(b(3,3))
    write(line,*) 'If we assign it to the matrix b:'; call print(line)
    b=diag(vect1)
    call print(b)
    write(line,*) 'Now the product of b * vect1:'; call print(line)
    !works, but you have to initialise the result
    vect2=b.mult.vect1
    call print(vect2)

    !works

    write(line,*) 'The result can be directly assigned to b that becomes:'; call print(line)
    b=b.multd.vect1
    call print(b)

    write(line,*) 'The scalar product b.dot. b is :'; call print(line)
    !works
    tmp=b.dot.b
    write(line,*) tmp; call print(line)
    call print(line)

    ! works
    write(line,*) 'Now b=:'; call print(line)
    b=1
    b(1,2)=5
    b(3,1)=0
    b(2,2)=2
    call print(b)

    write(line,*) 'The rowcol product b*b is'; call print(line)
    call print(b.mult.b)
    write(line,*) 'For comparison: The matmul says that'; call print(line)
    allocate(d(size(b,1),size(b,2)))
    d=matmul(b,b)
    call print(d)
    write(line,*) 'Now b=:'; call print(line)
    ! transposition works
    b(1,2)=-1
    call print(b)
    write(line,*) 'Transposing:'; call print(line)
    b=transpose(b)
    call print(b)

    !works
    write(line,*) 'And Symmetrise:'; call print(line)
    call matrix_symmetrise(b)
    call print(b)

    allocate(c(3,3))

    c=b
#ifndef useintrinsicblas
    write(line,*) 'Eigenproblem of b?:'; call print(line)
    ! works
    call matrix_diagonalise(b,vect1,c)
    write(line,*) 'Eigenvectors:'; call print(line)
    call print(c)
    write(line,*) 'Eigenvalues:'; call print(line)
    call print(vect1)


    ! general diagonalization
    write(line,*) 'Generalise eigenvaue problem:'; call print(line)
    call diagonalise(b,b,vect1,c)
    write(line,*) 'Eigenvectors:'; call print(line)
    call print(c)
    write(line,*) 'Eigenvalues:'; call print(line)
    call print(vect1)


    ! check result of diagonalizations with matlab

    write(line,*) 'Now b=: (simmetric?)'; call print(line)
    b(1,:)=1
    b(2,:)=2
    b(3,1)=1
    b(3,2)=3
    b(3,3)=4
    call matrix_symmetrise(b)
    b(3,:)=0
    b(:,3)=0
    b(2,2)=4
    b(3,3)=6

    call print(b)

    write(line,*) 'Its inverse is(symmetric positive definite case)'; call print(line)
    call matrix_inverse(b,c)!works
    call print(c)
    write(line,*) 'Their rowcol product is:'; call print(line)
    call print(matrix_product(c,b))
    write(line,*) 'Is it the unity matrix? thats good'; call print(line)

    b(2,1)=-2.0
    b(2,1)=3.5
    b(3,2)=7
    write(line,*) 'Now b=: (general case)'; call print(line)
    call print(b)

    write(line,*) 'Its inverse is(general case)'; call print(line)
    call matrix_inverse(b,c)! works
    call print(c)
    write(line,*) 'Their rowcol product is:'; call print(line)
    call print(matrix_product(c,b))
    write(line,*) 'Is it the unity matrix? thats good'; call print(line)
#endif


    write(line,*)'The following should be F and T ', b.FEQ.c, b.FEQ.b
    call print(line)

    write(line,*) 'Putting diag(b) into a vector'; call print(line)
    vect1=diag(b)
    call print(vect1)

    write(line,*)'Checking CFCT, b is (34,48,48,110)? :'
    call print(line)

    deallocate(b)
    allocate(b(2,5))
    b(1,1:3)=3.0
    b(2,1:2)=1.0
    b(2,3:5)=6.0
    b(1,4:5)=2.0

    call print(b)

    deallocate(vect1)
    allocate(vect1(5))
    vect1=1.0_dp

    call print(matrix_cfct(b,vect1))
    
     write(line,*)'Checking write_file_direct, b:'
    call print(line)
     call print(b)
!     call matrix_file_write_direct(b,'pippo.dat')
     b=0.0_dp
 !    call matrix_file_read_direct(b,'pippo.dat')
     write(line,*)'We write it to a file and read back. b:'
    call print(line)
      call print(b)
  end subroutine matrix_test

  !% OMIT
  subroutine vector_test()
    real(dp),allocatable::vector(:),vector1(:),vector2(:)

    write(line,*)'CHECKING VECTORS';call print(line)
    write(line,*)'';call print(line)
    write(line,*)'a ten elements vector of 0s';call print(line)
    allocate(vector(10))
    allocate(vector2(10))
    vector=0.0_dp
    call print(vector)
    vector=1
    write(line,*)'a ten elements vector of 1s';call print(line)
    call print(vector)
    write(line,*)'norm And normsq:',normsq(vector),vector_norm(vector)
    call print(line)
    write(line,*)'Assign it to another vector';call print(line)
    allocate(vector1(10))
    vector1=vector
    call print(vector1)
    write(line,*)'Assign a constant:7';call print(line)
    vector1=7.0_dp
    call print(vector1)

    ! think about overloading
    write(line,*)'Dotproduct vector .DOT. vector',dot_product(vector1,vector1)
    call print(line)
    write(line,*)'The outer product is'; call print(line)
    call print(outer(vector,vector1))
    !    call print(line)
    vector2=sqrt_cut(vector1)
    write(line,*)'sqrt_cut:';call print(line)
    call print(vector2)
    write(line,*)'Vector_reandomise:'; call print(line) 
    call vector_randomise(vector1,6.0_dp)
    call print(vector1)
    if(vector1.FEQ.vector1) then 
       write(line,*)'FEQ OK' 
       call print(line)
    end if
    if(vector1.FNE.vector) then  
       write(line,*)'FNE OK' 
       call print(line)
    end if
    write(line,*)'Checking print:OK if last 2 are different'
    call print(line)

    deallocate (vector)
    allocate(vector(17))
    vector=3.57
    vector(16)=5.0
    vector(17)=-4.1
    call print(vector)

    write(line,*)'Checking vector histogram'
    call print(line)

    write(line,*) vector_histogram(vector,0.0_dp,5.0_dp,5)
    call print(line)
  end subroutine vector_test

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X  Array size and shape checking routines
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  ! The one dimensional checks seem unnecessarily complicated, but they can be
  ! easily generalised to higher dimensions (an extra : in intarray declaration
  ! for example)

  subroutine check_size_int_dim1(arrayname,intarray,n,caller,error)

    character(*),           intent(in) :: arrayname ! The name of the array to be checked
    integer, dimension(:),  intent(in) :: intarray  ! The array to be tested
    integer, dimension(:),  intent(in) :: n         ! The size that intarray should be
    character(*),           intent(in) :: caller    ! The name of the calling routine
    integer, intent(out), optional :: error

    integer, dimension(:), allocatable :: actual_size ! container for the actual size of intarray
    logical                            :: failed      ! Keeps track of any failures 
    integer                            :: i           ! dimension iterator

    INIT_ERROR(error)

    failed = .false.
    allocate( actual_size( size(shape(intarray)) ) )
    actual_size = shape(intarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
       RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_int_dim1

  !overloaded subroutine which allows a scalar 'n' in the one dimensional case
  subroutine check_size_int_dim1_s(arrayname,intarray,n,caller,error)
    character(*),           intent(in) :: arrayname 
    integer, dimension(:),  intent(in) :: intarray  
    integer,                intent(in) :: n         
    character(*),           intent(in) :: caller
    integer, intent(out), optional :: error
    
    INIT_ERROR(error)
    call check_size(arrayname,intarray,(/n/),caller,error)
    PASS_ERROR(error)

  end subroutine check_size_int_dim1_s

  subroutine check_size_int_dim2(arrayname,intarray,n,caller,error)

    character(*),             intent(in) :: arrayname 
    integer, dimension(:,:),  intent(in) :: intarray  
    integer, dimension(:),    intent(in) :: n         
    character(*),             intent(in) :: caller
    integer, intent(out), optional :: error

    integer, dimension(:), allocatable :: actual_size 
    logical                            :: failed       
    integer                            :: i          

    INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(intarray)) ) )
    actual_size = shape(intarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
       RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_int_dim2

  subroutine check_size_real_dim1(arrayname,realarray,n,caller,error)

    character(*),            intent(in) :: arrayname
    real(dp), dimension(:),  intent(in) :: realarray 
    integer,  dimension(:),  intent(in) :: n         
    character(*),            intent(in) :: caller
    integer, intent(out), optional :: error     

    integer, dimension(:), allocatable :: actual_size 
    logical                            :: failed      
    integer                            :: i           

    INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(realarray)) ) )
    actual_size = shape(realarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
       RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_real_dim1

  !overloaded subroutine which allows a scalar 'n' in the one dimensional case
  subroutine check_size_real_dim1_s(arrayname,realarray,n,caller,error)
    character(*),           intent(in) :: arrayname 
    real(dp), dimension(:), intent(in) :: realarray  
    integer,                intent(in) :: n         
    character(*),           intent(in) :: caller
    integer, intent(out), optional :: error     
    
    INIT_ERROR(error)
    call check_size(arrayname,realarray,(/n/),caller,error)
    PASS_ERROR(error)

  end subroutine check_size_real_dim1_s

  subroutine check_size_real_dim2(arrayname,realarray,n,caller, error)

    character(*),              intent(in) :: arrayname 
    real(dp), dimension(:,:),  intent(in) :: realarray 
    integer,  dimension(:),    intent(in) :: n        
    character(*),              intent(in) :: caller
    integer, intent(out), optional :: error  

    integer, dimension(:), allocatable :: actual_size
    logical                            :: failed      
    integer                            :: i          

    INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(realarray)) ) )
    actual_size = shape(realarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
       RAISE_ERROR(trim(caller) //': Size checking failed. Expected: ' // n // ', got: ' // actual_size, error)
    end if

  end subroutine check_size_real_dim2

#ifdef HAVE_QP
  subroutine check_size_quad_dim1(arrayname,quadarray,n,caller,error)

    character(*),            intent(in) :: arrayname
    real(qp), dimension(:),  intent(in) :: quadarray 
    integer,  dimension(:),  intent(in) :: n         
    character(*),            intent(in) :: caller
    integer, intent(out), optional :: error

    integer, dimension(:), allocatable :: actual_size 
    logical                            :: failed      
    integer                            :: i           
    
    INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(quadarray)) ) )
    actual_size = shape(quadarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
       RAISE_ERROR(trim(caller)//': Size checking failed',error)
    end if

  end subroutine check_size_quad_dim1

  !overloaded subroutine which allows a scalar 'n' in the one dimensional case
  subroutine check_size_quad_dim1_s(arrayname,quadarray,n,caller,error)
    character(*),           intent(in) :: arrayname 
    real(qp), dimension(:), intent(in) :: quadarray  
    integer,                intent(in) :: n         
    character(*),           intent(in) :: caller
    integer, intent(out), optional :: error     

    INIT_ERROR(error)
    call check_size(arrayname,quadarray,(/n/),caller,error)
    PASS_ERROR(error)

  end subroutine check_size_quad_dim1_s

  subroutine check_size_quad_dim2(arrayname,quadarray,n,caller,error)

    character(*),              intent(in) :: arrayname 
    real(qp), dimension(:,:),  intent(in) :: quadarray 
    integer,  dimension(:),    intent(in) :: n        
    character(*),              intent(in) :: caller   
    integer, intent(out), optional :: error     

    integer, dimension(:), allocatable :: actual_size
    logical                            :: failed      
    integer                            :: i          

    INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(quadarray)) ) )
    actual_size = shape(quadarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
       RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_quad_dim2
#endif

  subroutine check_size_complex_dim1(arrayname,realarray,n,caller,error)

    character(*),            intent(in) :: arrayname
    complex(dp), dimension(:),  intent(in) :: realarray 
    integer,  dimension(:),  intent(in) :: n         
    character(*),            intent(in) :: caller
    integer, intent(out), optional :: error     
   
    integer, dimension(:), allocatable :: actual_size 
    logical                            :: failed      
    integer                            :: i           

    INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(realarray)) ) )
    actual_size = shape(realarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
       RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_complex_dim1

  !overloaded subroutine which allows a scalar 'n' in the one dimensional case
  subroutine check_size_complex_dim1_s(arrayname,realarray,n,caller,error)
    character(*),           intent(in) :: arrayname 
    complex(dp), dimension(:), intent(in) :: realarray  
    integer,                intent(in) :: n         
    character(*),           intent(in) :: caller
    integer, intent(out), optional :: error     
    
    INIT_ERROR(error)
    call check_size(arrayname,realarray,(/n/),caller)
    PASS_ERROR(error)

  end subroutine check_size_complex_dim1_s

  subroutine check_size_complex_dim2(arrayname,realarray,n,caller, error)

    character(*),              intent(in) :: arrayname 
    complex(dp), dimension(:,:),  intent(in) :: realarray 
    integer,  dimension(:),    intent(in) :: n        
    character(*),              intent(in) :: caller 
    integer, intent(out), optional :: error     

    integer, dimension(:), allocatable :: actual_size
    logical                            :: failed      
    integer                            :: i          

    INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(realarray)) ) )
    actual_size = shape(realarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
       RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_complex_dim2

  subroutine check_size_log_dim1(arrayname,logarray,n,caller,error)

    character(*),           intent(in) :: arrayname 
    logical, dimension(:),  intent(in) :: logarray  
    integer, dimension(:),  intent(in) :: n         
    character(*),           intent(in) :: caller
    integer, intent(out), optional :: error     
    
    integer, dimension(:), allocatable :: actual_size
    logical                            :: failed     
    integer                            :: i          

    INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(logarray)) ) )
    actual_size = shape(logarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
       RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_log_dim1

  !overloaded subroutine which allows a scalar 'n' in the one dimensional case
  subroutine check_size_log_dim1_s(arrayname,logarray,n,caller,error)
    character(*),           intent(in) :: arrayname 
    logical, dimension(:),  intent(in) :: logarray  
    integer,                intent(in) :: n         
    character(*),           intent(in) :: caller
    integer, intent(out), optional :: error     
    
    INIT_ERROR(error)
    call check_size(arrayname,logarray,(/n/),caller)
    PASS_ERROR(error)

  end subroutine check_size_log_dim1_s

  subroutine check_size_log_dim2(arrayname,logarray,n,caller,error)

    character(*),             intent(in) :: arrayname 
    logical, dimension(:,:),  intent(in) :: logarray  
    integer, dimension(:),    intent(in) :: n         
    character(*),             intent(in) :: caller    
    integer, intent(out), optional :: error     
    
    integer, dimension(:), allocatable :: actual_size
    logical                            :: failed     
    integer                            :: i          

    INIT_ERROR(error)
    failed = .false.
    allocate( actual_size( size(shape(logarray)) ) )
    actual_size = shape(logarray)

    if (size(actual_size) /= size(n)) then
       write(line,'(a,i0,a,i0,a)') caller//': '//arrayname//' is ',size(actual_size), &
            ' dimensional and not ',size(n),' dimensional as expected'
       call print(line)
       failed = .true.
    else
       do i = 1, size(actual_size)
          if (actual_size(i) /= n(i)) then
             write(line,'(3(a,i0),a)') caller//': The size of dimension ',i,' of '//arrayname//' is ', &
                  actual_size(i),' and not ',n(i),' as expected'
             call print(line)
             failed = .true.
          end if
       end do
    end if

    if (failed) then
       RAISE_ERROR(trim(caller)//': Size checking failed', error)
    end if

  end subroutine check_size_log_dim2  


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X Array Searching/Sorting/Averaging
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   !% Collect the indices corresponding to .true. elements of a mask
   !% can be used to index other arrays, e.g. array(find(mask))
   function find_indices(mask)
     logical::mask(:)
     integer::find_indices(count(mask))
     integer::i,n
     n = 1
     do i=1,size(mask)
        if(mask(i)) then
           find_indices(n) = i
           n = n+1
        end if
     end do
   end function find_indices
   

   !% Test if a given integer 'val' is in an integer array 'array'.
   pure function is_in_array(this, val)

     
     integer, intent(in) :: val
     integer, intent(in), dimension(:) :: this
     logical             :: Is_In_Array

     if(find_in_array(this, val) > 0) then
        is_in_array = .true.
     else
        is_in_array = .false.
     end if

   end function is_in_array

   pure function find_in_array_element_i(this, val) result(n)
     integer, intent(in) :: val
     integer, intent(in), dimension(:) :: this
     integer            ::i, n

     ! Optimised to avoid overhead of function call to find_in_array_row (jrk33)

     do i = 1,size(this)
        if (this(i) == val) then
           n = i
           return
        end if
     end do

     ! not found
     n = 0

   end function find_in_array_element_i

   pure function find_in_array_element_s(this, val) result(n)
     character(len=*), intent(in) :: val
     character(len=*), intent(in), dimension(:) :: this
     integer            ::i, n

     ! Optimised to avoid overhead of function call to find_in_array_row (jrk33)

     do i = 1,size(this)
        if (this(i) == val) then
           n = i
           return
        end if
     end do

     ! not found
     n = 0

   end function find_in_array_element_s

   function find_in_array_row(this, val, mask) result(n)
     integer, intent(in), dimension(:,:) :: this
     integer, intent(in), dimension(:) :: val
     logical, optional, intent(in), dimension(:) :: mask ! if this exists, we only compare elements where this is 1
     integer             :: i, j, n
     logical :: match

     if(present(mask)) then
        if (size(val) /= size(mask)) &
             call system_abort('Table_Find_Row: mask and row size mismatch')
     end if

     if (present(mask)) then
        do i = 1,size(this,2)
           match = .true.
           do j=1,size(this,1)
              if (.not. mask(j)) cycle
              if (this(j,i) /= val(j)) then
                 match = .false.
                 exit
              end if
           end do
           if (.not. match) cycle
           n = i
           return
        end do
     else
        do i = 1,size(this,2)
           match = .true.
           do j=1,size(this,1)
              if (this(j,i) /= val(j)) then
                 match = .false.
                 exit
              end if
           end do
           if (.not. match) cycle
           n = i
           return
        end do
     end if

     ! not found
     n = 0
   end function find_in_array_row


   !% Return a copy of the integer array 'array' with duplicate
   !% elements removed. 'uniq' assumes the input array is sorted
   !% so that duplicate entries appear consecutively.
   subroutine uniq(array, uniqed)
     integer, intent(in), dimension(:) :: array
     integer, dimension(:), allocatable, intent(out) :: uniqed

     integer, dimension(size(array)) :: seen
     integer :: i, n_seen

     n_seen = 1
     seen(1) = array(1)

     do i=2,size(array)
        if (array(i) == array(i-1)) cycle
        n_seen = n_seen + 1
        seen(n_seen) = array(i)
     end do

     allocate(uniqed(n_seen))
     uniqed = seen(1:n_seen)

   end subroutine uniq


   !% Sort an array of integers into ascending order (slow: scales as N$^2$).
   !% r_data is an accompanying array of reals on which the same reordering is performed
   subroutine sort_array_i(array, r_data)

      integer, dimension(:), intent(inout) :: array
      real(dp), dimension(:), intent(inout), optional :: r_data
      integer                              :: i,j,minpos
      integer :: min, tmp
      real(dp) :: r_tmp


      do i = 1, (size(array) - 1)

         min = huge(0)

         do j = i,size(array)

            if (array(j) < min) then
               min = array(j)
               minpos = j
            end if

         end do

         tmp = array(i)
         array(i) = array(minpos)
         array(minpos) = tmp
	 if (present(r_data)) then
	   r_tmp = r_data(i)
	   r_data(i) = r_data(minpos)
	   r_data(minpos) = r_tmp
	 endif

      end do

   end subroutine sort_array_i


   !% Sort an array of integers into ascending order (slow: scales as N$^2$).
   !% i_data is an accompanying array of integers on which the same reordering is performed
   subroutine sort_array_r(array, i_data,r_data)

      real(dp), dimension(:), intent(inout) :: array
      integer, dimension(:), intent(inout), optional :: i_data
      real(dp), dimension(:), intent(inout), optional :: r_data
      integer                              :: i,j, minpos
      real(dp) :: tmp, min
      integer :: i_tmp
      real(dp) :: r_tmp


      do i = 1, (size(array) - 1)

         min = huge(0.0_dp)

         do j = i,size(array)

            if (array(j) < min) then
               min = array(j)
               minpos = j
            end if

         end do

         tmp = array(i)
         array(i) = array(minpos)
         array(minpos) = tmp
	 if (present(i_data)) then
	   i_tmp = i_data(i)
	   i_data(i) = i_data(minpos)
	   i_data(minpos) = i_tmp
	 endif
	 if (present(r_data)) then
	   r_tmp = r_data(i)
	   r_data(i) = r_data(minpos)
	   r_data(minpos) = r_tmp
	 endif

      end do

   end subroutine sort_array_r


   !% Sort an array of integers into ascending order.
   !% The function uses heapsort, which always scales as N log N.
   !% r_data is an accompanying array of reals on which the same reordering is
   !% performed
   !% (Initial implementation by Andreas Wonisch)
   subroutine heap_sort_i(array, r_data)
     integer,            dimension(:), intent(inout)  :: array
     real(dp), optional, dimension(:), intent(inout)  :: r_data

     ! ---

     integer   :: N, i, j, root, tmpi
     real(DP)  :: tmpr

     ! ---

     N = size(array)

     do i = N/2, 1, -1

        j = i
        call siftdown(j, N)

     enddo

     do i = N, 2, -1
        
        ! Swap
        tmpi       = array(1)
        array(1)   = array(i)
        array(i)   = tmpi
        if (present(r_data)) then
           tmpr       = r_data(1)
           r_data(1)  = r_data(i)
           r_data(i)  = tmpr
        endif

        root = 1
        j = i -1
        call siftdown(root, j)
    
     enddo

   contains

     subroutine siftdown(root, bottom)
       integer, intent(inout)  :: root
       integer, intent(in)     :: bottom

       ! ---

       logical  :: done
       integer  :: maxchild

       ! ---

       done = .false.

       do while ((root*2 <= bottom) .and. .not. done)

          if (root*2 == bottom) then
             maxchild = root * 2
          else if (array(root*2) > array(root*2+1)) then
             maxchild = root * 2
          else
             maxchild = root*2 + 1
          endif

          if (array(root) < array(maxchild)) then

             ! Swap
             tmpi              = array(root)
             array(root)       = array(maxchild)
             array(maxchild)   = tmpi
             if (present(r_data)) then
                tmpr              = r_data(root)
                r_data(root)      = r_data(maxchild)
                r_data(maxchild)  = tmpr
             endif

             root = maxchild
          else
             done = .true.
          endif

       enddo

     endsubroutine siftdown

   endsubroutine heap_sort_i


   !% Sort an array of integers into ascending order.
   !% The function uses heapsort, which always scales as N log N.
   !% i_data and r_data are a accompanying arrays of integers and reals
   !% on which the same reordering is performed
   !% (Initial implementation by Andreas Wonisch)
   subroutine heap_sort_r(array, i_data, r_data)
     real(dp),           dimension(:), intent(inout)  :: array
     integer,  optional, dimension(:), intent(inout)  :: i_data
     real(dp), optional, dimension(:), intent(inout)  :: r_data

     ! ---

     integer   :: N, i, j, root, tmpi
     real(DP)  :: tmpr

     ! ---

     N = size(array)

     do i = N/2, 1, -1

        j = i
        call siftdown(j, N)

     enddo

     do i = N, 2, -1
        
        ! Swap
        tmpr       = array(1)
        array(1)   = array(i)
        array(i)   = tmpr
        if (present(i_data)) then
           tmpi       = i_data(1)
           i_data(1)  = i_data(i)
           i_data(i)  = tmpi
        endif
        if (present(r_data)) then
           tmpr       = r_data(1)
           r_data(1)  = r_data(i)
           r_data(i)  = tmpr
        endif

        root = 1
        j = i -1
        call siftdown(root, j)
    
     enddo

   contains

     subroutine siftdown(root, bottom)
       integer, intent(inout)  :: root
       integer, intent(in)     :: bottom

       ! ---

       logical  :: done
       integer  :: maxchild

       ! ---

       done = .false.

       do while ((root*2 <= bottom) .and. .not. done)

          if (root*2 == bottom) then
             maxchild = root * 2
          else if (array(root*2) > array(root*2+1)) then
             maxchild = root * 2
          else
             maxchild = root*2 + 1
          endif

          if (array(root) < array(maxchild)) then

             ! Swap
             tmpr             = array(root)
             array(root)      = array(maxchild)
             array(maxchild)  = tmpr
             if (present(i_data)) then
                tmpi              = i_data(root)
                i_data(root)      = i_data(maxchild)
                i_data(maxchild)  = tmpi
             endif
             if (present(r_data)) then
                tmpr              = r_data(root)
                r_data(root)      = r_data(maxchild)
                r_data(maxchild)  = tmpr
             endif

             root = maxchild
          else
             done = .true.
          endif

       enddo

     endsubroutine siftdown

   endsubroutine heap_sort_r


   !% Do an in place insertion sort on 'this', in ascending order.
   !% If 'idx' is present  then on exit it will contain the list 
   !% of indices into 'this' permuted in the same way as the entries 
   !% have been.
   subroutine insertion_sort_i(this, idx)
     integer, intent(inout), dimension(:) :: this
     integer, intent(out), dimension(size(this)), optional :: idx

     integer :: i, vi, j, v
     vi = 1

     if(present(idx)) idx = (/ (i, i=1,size(this)) /)
    
     do i = 2, size(this)
        if (this(i) >= this(i-1)) cycle
        v = this(i)
        if (present(idx)) vi = idx(i)
       
        j = i-1
        do while (j >= 1) 
           if (v > this(j)) exit
           this(j+1) = this(j)
           this(j) = v
           if (present(idx)) then
              idx(j+1) = idx(j)
              idx(j) = vi
           end if
           j = j - 1          
        end do
     end do

   end subroutine insertion_sort_i

   subroutine insertion_sort_r(this, idx)
     real(dp), intent(inout), dimension(:) :: this
     integer, intent(out), dimension(size(this)), optional :: idx

     integer :: i, vi, j
     real(dp) :: v
     vi = 1

     if(present(idx)) idx = (/ (i, i=1,size(this)) /)
    
     do i = 2, size(this)
        if (this(i) >= this(i-1)) cycle
        v = this(i)
        if (present(idx)) vi = idx(i)
       
        j = i-1
        do while (j >= 1) 
           if (v > this(j)) exit
           this(j+1) = this(j)
           this(j) = v
           if (present(idx)) then
              idx(j+1) = idx(j)
              idx(j) = vi
           end if
           j = j - 1          
        end do
     end do

   end subroutine insertion_sort_r

   !% Do binary search and return 'index' of element containing 'value', or zero if not found.
   !% 'array' must be sorted into ascending order beforehand.
   !% If the array subscripts don't start at 1, then pass the actual index of the first element as 'first',
   !% then, if the element isn't found, the value returned is 'first' minus 1.
   !% To restrict the search to a subsection of the array supply 'low' and/or 'high'.

   function binary_search_i(array,value,first,low,high,error) result(index)
     
     integer, dimension(:), intent(in)   :: array
     integer,               intent(in)   :: value
     integer, optional,     intent(in)   :: first,low,high
     integer, optional,     intent(out)  :: error
     integer                             :: index
     !local variables
     integer                             :: ilow, ihigh, count, max
     logical                             :: done
     
     INIT_ERROR(error)
     max = size(array)
     ilow = 1; ihigh = max
     if (present(low)) then
        ilow = low
        if (present(first)) ilow = ilow - first + 1
        ASSERT(ilow >= 1, 'binary_search: Bad lower index supplied', error)
     end if
     if (present(high)) then
        ihigh = high
        if (present(first)) ihigh = ihigh - first + 1
        ASSERT(ihigh <= max, 'binary_search: Bad higher index supplied', error)
     end if
     ASSERT(ilow <= ihigh, 'binary_search: Lower index > Higher index!', error)

     index = 0; count = 0
     done = .false.
     
     ASSERT(max >= 1, 'binary_search: Array must have at lease one element!', error)

     if ( (array(ilow) > value) .or. (array(ihigh) < value) ) done = .true.
     if (array(ilow) == value) then
        index = ilow
        done = .true.
     end if
     if (array(ihigh) == value) then
        index = ihigh
        done = .true.
     end if
     
     do while(.not.done)
        count = count + 1
        index = (ihigh + ilow) / 2
        if (index == ilow) then              ! value is not present. exit
           index = 0
           done = .true.
        else if (array(index) == value) then ! value found
           done = .true.
        else if (array(index) < value) then  ! value at this index is too low. shift lower bound
           ilow = index
        else                                 ! value at this index is too high. shift upper bound
           ihigh = index
        end if
        ASSERT(count < max, 'binary_search: Counter hit maximum. Is the array sorted properly?', error)
     end do
     
     if (present(first)) index = index - 1 + first
     
   end function binary_search_i


   !% Do binary search and return 'index' of element being smaller or equal
   !% to 'value'. 'array' must be sorted into ascending order beforehand.
   !% If the array subscripts don't start at 1, then pass the actual index of
   !% the first element as 'first'
   function binary_search_r(array,value,first,low,high,error) result(index)
     
     real(DP), dimension(:), intent(in)   :: array
     real(DP),               intent(in)   :: value
     integer, optional,      intent(in)   :: first,low,high
     integer, optional,      intent(out)  :: error
     integer                              :: index
     !local variables
     integer                              :: ilow, ihigh, count, max
     logical                              :: done

     INIT_ERROR(error)
     max = size(array)
     ilow = 1; ihigh = max
     if (present(low)) then
        ilow = low
        if (present(first)) ilow = ilow - first + 1
        ASSERT(ilow >= 1, 'binary_search: Bad lower index supplied', error)
     end if
     if (present(high)) then
        ihigh = high
        if (present(first)) ihigh = ihigh - first + 1
        ASSERT(ihigh <= max, 'binary_search: Bad higher index supplied', error)
     end if
     ASSERT(ilow <= ihigh, 'binary_search: Lower index > Higher index!', error)

     index = 0; count = 0
     done = .false.
     
     ASSERT(max >= 1, 'binary_search: Array must have at lease one element!', error)

     if (array(ilow) > value) then
        index = ilow-1
        done = .true.
     end if
     if (array(ihigh) <= value) then
        index = ihigh
        done = .true.
     end if
     
     do while(.not.done)
        count = count + 1
        index = (ihigh + ilow) / 2
        if (array(index) <= value .and. array(index+1) > value) then
           ! value found
           done = .true.
        else if (array(index) < value) then
           ! value at this index is too low. shift lower bound
           ilow = index
        else
           ! value at this index is too high. shift upper bound
           ihigh = index
        endif
        ASSERT(count < max, 'binary_search: Counter hit maximum. Is the array sorted properly?', error)
     end do
     
     if (present(first)) index = index - 1 + first
     
   end function binary_search_r


   !% Subtract the average value from each element of a real array.
   subroutine zero_sum(array)

      real(dp), dimension(:,:), intent(inout) :: array
      real(dp), dimension(size(array,1))      :: array_avg
      integer                                 :: i

      array_avg = average_array(array)

      ! Update each component
      do i = 1, size(array,2)
         array(:,i) = array(:,i) - array_avg
      end do

   end subroutine zero_sum

   !% Calculate the average of a two dimenstional real array across its second dimension.
   function average_array(array)

      real(dp), dimension(:,:), intent(in) :: array
      real(dp), dimension(size(array,1))   :: average_array
      integer                              :: i

      average_array = 0.0_dp
      do i = 1, size(array,2)
         average_array = average_array + array(:,i)
      end do

      average_array = average_array / real(size(array,2),dp)

   end function average_array



  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
  !
  !% Linear Least Squares fit using \textsc{lapack} to do singular value decomposition.
  !%
  !% Fit data 'y' with errors 'sig' to model using parameters 'a', 
  !% and calculate $\chi^2$ of the fit. The user defined subroutine 'funcs'
  !% should return the model parameters 'a' for the point 'x' in the array 'afunc'.
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  

  subroutine least_squares(x, y, sig, a, chisq, funcs)

    real(dp), intent(in), dimension(:)  :: x, y, sig
    real(dp), intent(out), dimension(:) :: a
    real(dp), intent(out) :: chisq

    real(dp) :: ain(size(x),size(a)), u(size(x),size(x)), v(size(a),size(a)),&
         w(size(a)), sigmainv(size(a), size(x))
    integer i,j,info
    real(dp), parameter :: TOL = 1e-13_dp
    real(dp) :: wmax, thresh, sum, b(size(x)), afunc(size(a))
    real(dp), allocatable, dimension(:) :: work
    integer ::  ndata, ma, lwork
    interface 
       subroutine funcs(x,afunc)
         use system_module
         real(dp) :: x, afunc(:)
       end subroutine funcs
    end interface

    ndata = size(x)
    ma = size(a)

    ! Call user fit function and scale by sig
    do i=1,ndata
       call funcs(x(i), afunc)
       do j=1,ma
          ain(i,j) = afunc(j)/sig(i)
       end do
       b(i) = y(i)/sig(i)
    end do

    lwork = 2 * max(3*min(ma,ndata)+max(ma,ndata),5*min(ma,ndata))
    allocate(work(lwork))
    call DGESVD('A','A',ndata,ma,ain,ndata,w,u,ndata,v,ma,work,lwork,info)
        
    if (info/=0) then
       if (info < 0) then
          write(line,'(a,i0)')'least_squares: Problem with argument ',-info
          call system_abort(line)
       else
          call system_abort('least_squares: DBDSQR (called from DGESVD) did not converge')
       end if
    end if
    deallocate(work)

    ! Zero singular elements
    wmax = maxval(w)
    thresh = TOL*wmax

    ! Invert singular value matrix
    sigmainv = 0.0_dp
    do i=1,ma
       if (w(i) < thresh) then
          sigmainv(i,i) = 0.0_dp
       else
          sigmainv(i,i) = 1.0_dp/w(i)
       end if
    end do

    a = transpose(v) .mult. sigmainv .mult. transpose(u) .mult. b

    ! Calculate chi^2 for the fitted parameters
    chisq = 0.0_dp
    do i=1,ndata
       call funcs(x(i), afunc)
       sum = 0.0_dp
       do j=1,ma
          sum = sum + a(j)*afunc(j)
       end do
       chisq = chisq + ((y(i)-sum)/sig(i))**2
    end do

  end subroutine least_squares

  !% This routine updates the average and variance of a quantity that is sampled periodically
  subroutine update_running_average_and_variance(average,variance,x,N)

    real(dp), intent(inout) :: average, variance
    real(dp), intent(in)    :: x !% The new sample value
    integer,  intent(in)    :: N !% The sample number

    real(dp)                :: f1, f2, old_av, old_var
    
    f1 = real(N-1,dp) / real(N,dp)
    f2 = 1.0_dp / real(N,dp)
    old_av = average
    old_var = variance

    !Update average
    average = f1 * old_av + f2 * x

    !Update variance
    variance = f1 * (old_var + old_av*old_av) + f2 * x*x - average*average

  end subroutine update_running_average_and_variance

  subroutine update_exponential_average_s(average,decay,x)

    real(dp), intent(inout) :: average
    real(dp), intent(in)    :: decay, x

    real(dp) :: aa(1), xa(1)

    !Pack the scalar data into a 1-component vector and call
    !the vector version of this routine
    
    aa(1) = average
    xa(1) = x
    call update_exponential_average(aa,decay,xa)

    average = aa(1)

  end subroutine update_exponential_average_s

  subroutine update_exponential_average_v(average,decay,x)

    real(dp), intent(inout) :: average(:)
    real(dp), intent(in)    :: decay, x(size(average))

    real(dp) :: f1, f2

    f1 = exp(-decay)
    f2 = 1.0_dp - f1

    average = f1*average + f2*x

  end subroutine update_exponential_average_v

  subroutine update_exponential_average_d2(average,decay,x)

    real(dp), intent(inout) :: average(:,:)
    real(dp), intent(in)    :: decay, x(size(average,1),size(average,2))

    real(dp) :: f1, f2

    f1 = exp(-decay)
    f2 = 1.0_dp - f1

    average = f1*average + f2*x

  end subroutine update_exponential_average_d2

  !
  ! If we need to bin some data in the range [a,b) with bin-width d, which bin (starting at 1) should
  ! the value x go into? Also check that x lies in the correct range.
  !
  function bin(a,b,d,x)

    real(dp), intent(in) :: a,b,d,x
    integer              :: bin

    if (x < a .or. x >= b) call system_abort('bin: data value '//round(x,3)//' does not lie in range [' &
         //round(a,3)//','//round(b,3)//')')

    bin = int((x-a)/d) + 1

  end function bin

  !
  ! Return the coordinate of the centre of the ith bin if we have n bins spanning the range [a,b) 
  !
  function bin_centre(a,b,n,i) result(c)

    real(dp), intent(in) :: a, b
    integer,  intent(in) :: n, i
    real(dp)             :: c

    c = a + (b-a)*(real(i,dp)-0.5_dp)/real(n,dp)

  end function bin_centre

  pure function int_array_ge(array1,array2) result(ge)

    integer, intent(in) :: array1(:)
    integer, intent(in) :: array2(size(array1))
    logical             :: ge
    integer             :: i

    ge = .true.
    i = 1
    do while(i <= size(array1))
       if (array1(i) < array2(i)) then
          ge =.false.
          return
       else if (array1(i) > array2(i)) then
          return
       end if
       i = i + 1
    end do

  end function int_array_ge

  pure function int_array_gt(array1,array2) result(gt)

    integer, intent(in) :: array1(:)
    integer, intent(in) :: array2(size(array1))
    logical             :: gt
    integer             :: i

    gt = .true.
    i = 1
    do while(i <= size(array1))
       if (array1(i) < array2(i)) then
          gt =.false.
          return
       else if (array1(i) > array2(i)) then
          return
       end if
       i = i + 1
    end do
    gt = .false.

  end function int_array_gt

  pure function int_array_lt(array1,array2) result(lt)

    integer, intent(in) :: array1(:)
    integer, intent(in) :: array2(size(array1))
    logical             :: lt

    lt = .not.int_array_ge(array1,array2)

  end function int_array_lt

  !% compare contents of 2 cells in up to N=2 arrays (int or real), return true if 
  !% contents in first cell is less than 2nd, with 1st array taking precendence.  
  function arrays_lt(i, j, i_p1, r_p1, i_p2, r_p2, i_p3, r_p3, error)
     integer, intent(in) :: i, j
     integer, intent(in), pointer, optional :: i_p1(:)
     real(dp), intent(in), pointer, optional :: r_p1(:)
     integer, intent(in), pointer, optional :: i_p2(:)
     real(dp), intent(in), pointer, optional :: r_p2(:)
     integer, intent(in), pointer, optional :: i_p3(:)
     real(dp), intent(in), pointer, optional :: r_p3(:)
     integer, intent(out), optional :: error
     logical :: arrays_lt

     integer :: n_p1, n_p2, n_p3
     logical :: p1_is_i = .false., p1_is_r = .false.
     logical :: p2_is_i = .false., p2_is_r = .false., have_p2 = .false.
     logical :: p3_is_i = .false., p3_is_r = .false., have_p3 = .false.

     INIT_ERROR(error)

     n_p1 = 0
     if (present(i_p1)) then
       if (associated(i_p1)) then
	  n_p1 = n_p1 + 1
	  p1_is_i = .true.
       endif
     end if
     if (present(r_p1)) then
       if (associated(r_p1)) then
	  n_p1 = n_p1 + 1
	  p1_is_r = .true.
       endif
     end if
     if (n_p1 /= 1) then
        RAISE_ERROR("arrays_lt got too many or too few present and associated p1 pointers",error)
     endif

     n_p2 = 0
     if (present(i_p2)) then
       if (associated(i_p2)) then
	  n_p2 = n_p2 + 1
	  p2_is_i = .true.
       endif
     end if
     if (present(r_p2)) then
       if (associated(r_p2)) then
	  n_p2 = n_p2 + 1
	  p2_is_r = .true.
       endif
     endif
     if (n_p2 > 1) then
        RAISE_ERROR("arrays_lt got too many present and associated p2 pointers",error)
     endif
     have_p2 = (n_p2 == 1)

     n_p3 = 0
     if (present(i_p3)) then
       if (associated(i_p3)) then
	  n_p3 = n_p3 + 1
	  p3_is_i = .true.
       endif
     end if
     if (present(r_p3)) then
       if (associated(r_p3)) then
	  n_p3 = n_p3 + 1
	  p3_is_r = .true.
       endif
     endif
     if (n_p3 > 1) then
        RAISE_ERROR("arrays_lt got too many present and associated p3 pointers",error)
     endif
     have_p3 = (n_p3 == 1)

     ! always have some *_p1, test *_p1
     if (p1_is_i) then
	if (i_p1(i) < i_p1(j)) then
	   arrays_lt = .true.
	   return
	else if (i_p1(i) > i_p1(j)) then
	   arrays_lt = .false.
	   return
	endif
     else ! p1_is_r
	if (r_p1(i) < r_p1(j)) then
	   arrays_lt = .true.
	   return
	else if (r_p1(i) > r_p1(j)) then
	   arrays_lt = .false.
	   return
	endif
     endif

     ! *_p1 must be equal, test *_p2
     if (have_p2) then ! must have some present and associated *_p2 pointer, test it
	if (p2_is_i) then
	   if (i_p2(i) < i_p2(j)) then
	      arrays_lt = .true.
	      return
	   else if (i_p2(i) > i_p2(j)) then
	      arrays_lt = .false.
	      return
	   endif
	else ! p2_is_r
	   if (r_p2(i) < r_p2(j)) then
	      arrays_lt = .true.
	      return
	   else if (r_p2(i) > r_p2(j)) then
	      arrays_lt = .false.
	      return
	   endif
	endif
     endif

     ! *_p2 must be equal, test *_p3
     if (have_p3) then ! must have some present and associated *_p3 pointer, test it
	if (p3_is_i) then
	   if (i_p3(i) < i_p3(j)) then
	      arrays_lt = .true.
	      return
	   else if (i_p3(i) > i_p3(j)) then
	      arrays_lt = .false.
	      return
	   endif
	else ! p3_is_r
	   if (r_p3(i) < r_p3(j)) then
	      arrays_lt = .true.
	      return
	   else if (r_p3(i) > r_p3(j)) then
	      arrays_lt = .false.
	      return
	   endif
	endif
     endif

     ! *_p3 must be equal
     arrays_lt = .false.
  end function arrays_lt

  subroutine polar_decomposition(m, S, R)
    real(dp), intent(in) :: m(3,3)
    real(dp), intent(out) :: S(3,3), R(3,3)

    real(dp) :: EEt(3,3), D(3), V(3,3)

    ! Find polar decomposition: m = S*R where S is symmetric, R is a rotation
    !  EEt = E*E', EEt = VDV' D diagonal, S = V D^1/2 V', R = S^-1*E

    EEt = m .mult. transpose(m) ! Normal
    call diagonalise(EEt, D, V)

    ! Check positive definite
    if (any(D < 0)) then
      call print(m, PRINT_VERBOSE)
      call system_abort("polar decomposition 'm' is not positive definite")
    end if

    S = V .mult. diag(sqrt(D)) .mult. transpose(V)
    R = V .mult. diag(D ** (-0.5_dp)) .mult. transpose(V) .mult. m

    call print('S:', PRINT_ANAL); call print(S, PRINT_ANAL)
    call print('R:', PRINT_ANAL); call print(R, PRINT_ANAL)

  end subroutine polar_decomposition

    !% Returns the volume of intersection of two spheres, radius $r_1$ and $r_2$, whose
  !% centres are a distance $d$ apart
  function sphere_intersection_vol(r1,r2,d) result(V)

    real(dp), intent(in) :: r1, r2, d
    real(dp)             :: V
    real(dp)             :: t1, t2, t3

    if (r1<0.0_dp .or. r2<0.0-dp .or. d<0.0_dp .or. (r1 + r2)>d) then
       V = 0.0_dp
       return
    else if (abs(r1-r2)>d) then
       t1 = min(r1,r2)
       V = 4.0_dp*PI*t1*t1*t1/3.0_dp
       return
    end if

    t1 = r1 + r2 - d
    t2 = r1 + r2 + d
    t3 = r1 - r2

    V = PI * t1 * t1 * (d * t2 - 3.0_dp * t3 * t3) / (12.0_dp * d)

  end function sphere_intersection_vol

  function permutation_symbol() result(eps)
    real(dp) :: eps(3,3,3)

    eps = 0.0_dp

    eps(1,2,3) =  1.0_dp
    eps(1,3,2) = -1.0_dp
    eps(2,3,1) =  1.0_dp
    eps(2,1,3) = -1.0_dp
    eps(3,1,2) =  1.0_dp
    eps(3,2,1) = -1.0_dp
  end function permutation_symbol

   ! Matrix3x3_Inverse
   !
   !% Calculate $3\times3$ matrix inverse of 'lattice' and store result in 'g'.
   !% Avoids overhead of calling \textsc{lapack} for simple case of $3\times3$ matrix.
   subroutine matrix3x3_inverse(matrix, g)
     real(dp), intent(in)  :: matrix(3,3)
     real(dp), intent(out) :: g(3,3)

     real(dp) :: stp

     stp = scalar_triple_product(matrix(:,1), matrix(:,2), matrix(:,3))

     g(1,:) = (matrix(:,2) .cross. matrix(:,3))/stp
     g(2,:) = (matrix(:,3) .cross. matrix(:,1))/stp
     g(3,:) = (matrix(:,1) .cross. matrix(:,2))/stp

   end subroutine matrix3x3_inverse

   !% Calulates determinant of $3\times3$ matrix
   function matrix3x3_det(m) result(det)
     real(dp), intent(in) :: m(3,3)
     real(dp) :: det

     det = m(1,1)*(m(2,2)*m(3,3) - m(3,2)*m(2,3)) &
          - m(2,1)*(m(1,2)*m(3,3) - m(3,2)*m(1,3)) &
          + m(3,1)*(m(1,2)*m(2,3) - m(2,2)*m(1,3))

   end function matrix3x3_det

   function pbc_aware_centre(p, lattice, g) result(c)
     real(dp), intent(in) :: p(:,:), lattice(3,3), g(3,3)
     real(dp) :: c(3)

     integer i
     complex(dp) :: c_z(3)

     c_z = 0.0_dp
     do i=1, size(p,2)
       c_z = c_z + exp(cmplx(0.0_dp,1.0_dp)*2.0_dp*PI*(g .mult. p(:,i)))
     end do
     c_z = c_z / real(size(p,2),dp)
     c_z = c_z/abs(c_z)
     c = lattice .mult. real(log(c_z)/(cmplx(0.0_dp,1.0_dp)*2.0_dp*PI),dp)

   end function pbc_aware_centre


   !% K-means clustering. Algorithm is as described in Numerical Recipes
   !% (Third Edition, Section 16.1.2, pp. 848-849).
   subroutine kmeans(data, K, means, assign, err, initialisation)
     real(dp), dimension(:, :), intent(in) :: data
     integer, intent(in) :: K
     real(dp), intent(inout), dimension(:,:) :: means
     integer, dimension(:), intent(inout) :: assign
     real(dp), intent(out), optional :: err
     character(*), intent(in), optional :: initialisation

     character(30) :: do_initialisation
     integer :: N, M, nn, mm, kk, kmin
     logical :: change
     real(dp) :: dmin, d

     N = size(data, 1)
     M = size(data, 2)

     do_initialisation = optional_default('none', initialisation)
     if (trim(do_initialisation) /= 'none' .and. &
         trim(do_initialisation) /= 'random_means' .and. &
         trim(do_initialisation) /= 'random_partition') &
         call system_abort('kmeans: initialisation argument should be one of "none" (default), "random_means", "random_partition"')

     if (trim(do_initialisation) == 'random_means') then
        ! Initialise centroids randomly within range of input data
        do kk=1,K
           do mm=1,M
              means(kk,mm) = ran_uniform()*(maxval(data(:,mm))-minval(data(:,mm))) + minval(data(:,mm))
           end do
        end do
     end if
     
     assign(:) = 0

     if (trim(do_initialisation) == 'random_partition') then
        ! Assign each centroid randomly to one of the K centroids
        do nn=1,N
           assign(nn) = mod(ran(), K)+1
        end do

        ! Compute the initial means from the random partition
        means(:,:) = 0.0_dp
        do nn=1,N
           do mm=1,M
              means(assign(nn),mm) = means(assign(nn),mm) + data(nn,mm)
           end do
        end do
        do kk=1,K
           if (any(assign == kk)) then
              means(kk,:) = means(kk,:) / count(assign == kk)
           end if
        end do
     end if

     change = .true.

     do while (change)
        change = .false.

        ! E-step -- assign each point to the component whose mean it is closest to
        do nn=1,N
           dmin = huge(1.0_dp)
           do kk=1,K
              d = (data(nn,:) - means(kk,:)) .dot. (data(nn,:) - means(kk,:))
              if (d < dmin) then
                 dmin = d
                 kmin = kk
              end if
           end do
           
           if (assign(nn) /= kmin) then
              assign(nn) = kmin
              change = .true.
           end if
        end do

        ! M-step -- update each mean to the average of the data points assigned to it
        means(:,:) = 0.0_dp
        do nn=1,N
           do mm=1,M
              means(assign(nn),mm) = means(assign(nn),mm) + data(nn,mm)
           end do
        end do
        do kk=1,K
           if (any(assign == kk)) then
              means(kk,:) = means(kk,:) / count(assign == kk)
           end if
        end do
     end do

     if (present(err)) then
        ! Compute residual error
        err = 0.0_dp
        do nn=1,N
           err = err + (data(nn,:) - means(assign(nn),:)) .dot. (data(nn,:) - means(assign(nn),:))
        end do
     end if

   end subroutine kmeans

   subroutine symmetric_linear_solve(M, a, M_inv_a)
    real(dp), intent(in) :: M(:,:), a(:)
    real(dp), intent(out):: M_inv_a(:)

    integer :: N
    integer :: lwork, info
    real(dp), allocatable :: work(:)
    integer, allocatable :: ipiv(:)

    N = size(a)

    allocate(ipiv(N))
    M_inv_a = a

    allocate(work(1))
    lwork = -1
    call dsysv('U', N, 1, M, N, ipiv, M_inv_a, N, work, lwork, info)
    if (info /= 0) &
      call system_abort("symmetric_linear_solve failed in lwork computation dsysv " // info)
    lwork = work(1)
    deallocate(work)

    allocate(work(lwork))
    call dsysv('U', N, 1, M, N, ipiv, M_inv_a, N, work, lwork, info)
    if (info /= 0) &
      call system_abort("symmetric_linear_solve failed in dsysv " // info)

    deallocate(work)
    deallocate(ipiv)
   end subroutine symmetric_linear_solve

  !% round n to the nearest number greater or equal to n
  !% with prime factors of only 2, 3, 5, ... max_prime_factor
  subroutine round_prime_factors(n, max_prime_factor, error)
    integer, intent(inout) :: n !% The number to be rounded
    integer, optional, intent(in) :: max_prime_factor
    integer, optional, intent(out) :: error

    integer :: factor, test, i, my_max_prime_factor, num_prime_factors

    INIT_ERROR(error)

    ! Check that n is not zero, otherwise we get stuck in an infinite loop
    ASSERT(n /= 0, 'n=0 has no prime factors', error)

    my_max_prime_factor = optional_default(5, max_prime_factor)
    select case (my_max_prime_factor)
       case ( 2) ; num_prime_factors = 1
       case ( 3) ; num_prime_factors = 2
       case ( 5) ; num_prime_factors = 3
       case ( 7) ; num_prime_factors = 4
       case (11) ; num_prime_factors = 5
       case (13) ; num_prime_factors = 6
       case (17) ; num_prime_factors = 7
       case (19) ; num_prime_factors = 8
       case default
          RAISE_ERROR('unsupported max_prime_factor '//my_max_prime_factor//' in prime_factors()', error)
    end select

    do
       test = n
       do i=1,num_prime_factors
          factor = 1
          select case (i)
             case (1) ; factor = 2
             case (2) ; factor = 3
             case (3) ; factor = 5
             case (4) ; factor = 7
             case (5) ; factor = 11
             case (6) ; factor = 13
             case (7) ; factor = 17
             case (8) ; factor = 19
             case default
                RAISE_ERROR('unsupported prime factor '//i//' in prime_factors()', error)
          end select

          do while(mod(test,factor).eq.0)
             test = test / factor
          enddo
       end do

       if (test == 1) return
       n = n + 1
    end do

  end subroutine round_prime_factors

  !% Calculates integral numerically using the trapezoid formula from (x,y) data
  !% pairs.
  function TrapezoidIntegral(x,y)
  
     real(dp) :: TrapezoidIntegral
     real(dp), dimension(:), intent(in) :: x, y
  
     integer :: n
  
     if( size(x) /= size(y) ) call system_abort('TrapezoidIntegral: dimensions of x and y not the same')
     n = size(x)
  
     TrapezoidIntegral = sum((x(2:n)-x(1:n-1))*(y(2:n)+y(1:n-1)))/2.0_dp
  
  endfunction TrapezoidIntegral

  function matrix_exp_d(a)
     real(dp), intent(in) :: a(:,:)
     real(dp) :: matrix_exp_d(size(a,1),size(a,2))

     integer :: n
     complex(dp), allocatable :: rv(:,:)
     complex(dp), allocatable :: ev(:), rv_inv(:,:)

     if (size(a,1) /= size(a,2)) call system_abort("matrix_exp_d needs a square matrix, not "//size(a,1)//"x"//size(a,2))

     n = size(a,1)
     allocate(rv(n,n), rv_inv(n,n), ev(n))

     call nonsymmetric_diagonalise(a, ev, r_evects=rv)
     call inverse(rv, rv_inv)

     matrix_exp_d = ((rv .multd. exp(ev)) .mult. rv_inv)

     deallocate(rv, rv_inv, ev)
   end function matrix_exp_d

   subroutine k_means_clustering_pick_1d(x, cluster_indices)
      real(dp), intent(in) :: x(:)
      integer, intent(out) :: cluster_indices(:)

      integer :: k
      real(dp), allocatable :: k_means(:)
      integer :: n, rv, i, j, closest_i, closest_j
      real(dp) :: closest_r, r, t_sum
      integer, allocatable :: a(:), prev_a(:)

! print *, "k_means_clustering_pick_1d ", k, size(cluster_indices)
! print *, "x ", x

      k = size(cluster_indices)
      n = size(x)
      allocate(k_means(k))

      ! random initial guess
      cluster_indices = 0
      do i=1, k
	 rv = 1+floor(ran_uniform()*n)
	 do while (any(cluster_indices == rv))
	    rv = 1+floor(ran_uniform()*n)
	 end do
	 cluster_indices(i) = rv
      end do
      k_means = x(cluster_indices)
! print *, "initial k_means ", k_means

      allocate(a(n), prev_a(n))

      ! evaluate initial assignments
      do j=1, n
	 closest_i = 0
	 closest_r = HUGE(1.0_dp)
	 do i=1, k
	    r = sqrt((x(j)-k_means(i))**2)
	    if (r < closest_r) then
	       closest_r = r
	       closest_i = i
	    endif
	 end do
	 a(j) = closest_i
      end do

! print *, "initial assignments ", a

      prev_a = 0
      do while (any(a /= prev_a))
! print *, "start loop"
	 prev_a = a

	 ! update positions
	 do i=1, k
	    t_sum = 0.0_dp
	    do j=1, n
	       if (a(j) == i) t_sum = t_sum + x(j)
	    end do
	    k_means(i) = t_sum / real(count(a == i),dp)
	 end do

! print *, "new k_means ", k_means

	 ! update assignments
	 do j=1, n
	    closest_i = 0
	    closest_r = HUGE(1.0_dp)
	    do i=1, k
	       r = sqrt((x(j)-k_means(i))**2)
	       if (r < closest_r) then
		  closest_r = r
		  closest_i = i
	       endif
	    end do
	    a(j) = closest_i
	 end do
! print *, "new assignments ", a
      end do

! print *, "calculate indices closest to means"
      do i=1, k
	 closest_i = 0
	 closest_r = HUGE(1.0_dp)
	 do j=1, n
	    r = sqrt((x(j)-k_means(i))**2)
	    if (r < closest_r) then
	       closest_r = r
	       closest_j = j
	    endif
	 end do
	 cluster_indices(i) = closest_j
      end do
! do j=1, n
!    closest_i = 0
!    closest_r = HUGE(1.0_dp)
!    do i=1, k
!       r = sqrt((x(j)-x(cluster_indices(i)))**2)
!       if (r < closest_r) then
! 	 closest_r = r
! 	 closest_i = i
!       endif
!    end do
!    a(j) = closest_i
! end do
! print *, "counts"
! do i=1, k
! print *, count(a == i)
! end do

      deallocate(a, prev_a, k_means)

   end subroutine k_means_clustering_pick_1d

end module linearalgebra_module
