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
!X  Sparse module
!X  
!X  Sparse matrix representation and algebra
!X  
!X
!X  Table:
!X  int :    [................]  column indices
!X  real:    [................]  nonzero elements
!X
!X  rows: integer array that holds the starting indices of each row in 
!X  the table
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module sparse_module
  use system_module
  use linearalgebra_module
  use table_module

  implicit none
  private
  SAVE

  public :: sparse, print, print_full, operator(.mult.), assignment(=), sparse_init, sparse_finalise

  type Sparse
     type(table)::table ! values
     integer,allocatable::rows(:) ! indexes
     integer::Nrows ! number of rows
     integer::Ncols ! number of colummns
  end type Sparse

  real(dp), parameter :: SPARSITY_GUESS = 0.01_dp ! A guess at the sparsity of the input matrices
                                                  ! to sparse_assign_matrix

  interface assignment(=)
     module procedure  sparse_assign_matrix,matrix_assign_sparse,&
          sparse_assign_sparse
  end interface assignment(=)

  interface operator(.mult.)
     module procedure sparse_mult_vector, vector_mult_sparse, sparse_mult_sparse
  end interface

  interface print
     module procedure sparse_print
  end interface print

  interface print_full
     module procedure sparse_print_full
  end interface print_full

contains

  ! just init
  subroutine sparse_init(this,Nrows,Ncols,Nmax)
    type(sparse):: this
    integer, optional::Nmax  ! maximum number of nonzero elements
    integer::Nrows,Ncols ! size of the matrix

    call allocate(this%table,1,1,0,0,Nmax)
    call reallocate(this%rows,Nrows+1)
    this%Nrows = Nrows
    this%Ncols = Ncols
    this%rows = 1

  end subroutine sparse_init

  subroutine sparse_finalise(this)
    type(sparse)::this
    call finalise(this%table)
    if(allocated(this%rows))deallocate(this%rows)
  end subroutine sparse_finalise


  !overloading assignment
  subroutine sparse_assign_sparse(this,other)
    type(sparse),intent(OUT):: this
    type(sparse),intent(IN):: other
    this%Nrows = other%Nrows
    this%Ncols = other%Ncols
    this%table = other%table
    call reallocate(this%rows,this%Nrows+1)
    this%rows = other%rows
  end subroutine sparse_assign_sparse

  ! assignment
  subroutine sparse_assign_matrix(this,matrix)
    type(sparse), intent(out) :: this
    real(dp),     intent(in)  :: matrix(:,:)
    integer                   :: Nmax, i, j
    
    Nmax = int( SPARSITY_GUESS * real(size(matrix),dp) )

    call sparse_init(this, size(matrix,1), size(matrix,2), Nmax)

    do i=1,this%Nrows
       this%rows(i) = this%table%N+1

       do j=1,this%Ncols

          if (abs(matrix(i,j)) > NUMERICAL_ZERO ) then
             call append(this%table,realpart=(matrix(i,j)),intpart=j)
          end if

       end do

    end do

    ! last value in rows()
    this%rows(this%Nrows+1) = this%table%N+1

    if(this%table%N == 0) then ! give warning
       call print('sparse_assign_matrix: no nonzero elements!')
    end if

    ! now truncate the length
    call allocate(this%table)

  end subroutine sparse_assign_matrix


  subroutine matrix_assign_sparse(matrix,sp)
    type(sparse),intent(IN):: sp
    real(dp),intent(OUT)::matrix(sp%Nrows,sp%Ncols)
    integer::i,j

    call sparse_check(sp,"matrix_assign_sparse") 
  
    matrix=0.0_dp
    do i = 1,sp%Nrows
       do j = sp%rows(i),sp%rows(i+1)-1
          matrix(i,sp%table%int(1,j)) = sp%table%real(1,j)
       end do
    end do
  end subroutine matrix_assign_sparse

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X                                                       X
  !X Routines for constructing/reading sparse matrices     X
  !X                                                       X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 
  !
  ! Return the (i,j)th element of the sparse matrix
  !
  function sparse_element(this,i,j)

    type(sparse), intent(in) :: this
    integer,      intent(in) :: i, j
    real(dp)                 :: sparse_element
    integer                  :: n

    call sparse_check_bounds(this,i,j,'Sparse_Element')
    
    sparse_element = 0.0_dp

    do n = this%rows(i), this%rows(i+1) - 1
       if (this%table%int(1,n) > j) exit    ! The requested element is 0.0
       if (this%table%int(1,n) == j) then   ! The requested element is stored
          sparse_element = this%table%real(1,n)
          exit
       end if
    end do

  end function sparse_element

  !
  ! Return the number of stored elements
  !
  function sparse_stored_elements(this)
  
    type(sparse), intent(in) :: this
    integer                  :: sparse_stored_elements

    sparse_stored_elements = this%table%N

  end function sparse_stored_elements

  !
  ! Zero all elements of a sparse matrix (keeping storage)
  !
  subroutine sparse_zero(this)
    
    type(sparse), intent(inout) :: this
    
    this%table%real = 0.0_dp


  end subroutine sparse_zero

  !
  ! Copy 'val' into the (i,j)th element of the sparse matrix
  !
  ! The optional logical 'nocheck', when present and true, supresses the deletion of an element 
  ! if we set it to zero. This can be used to insert placeholders for future non-zero values, but note
  ! that these will be classed as 'non-zero' in sparse_non_zero
  !
  subroutine sparse_set_element(this,i,j,val,nocheck)

    type(sparse),      intent(inout) :: this
    integer,           intent(in)    :: i,j
    real(dp),          intent(in)    :: val
    logical, optional, intent(in)    :: nocheck
    !local variables
    integer                          :: n, table_pos, intpart, tmp_intpart
    logical                          :: element_written, my_nocheck
    real(dp)                         :: realpart, tmp_realpart

    if (present(nocheck)) then
       my_nocheck = nocheck
    else
       my_nocheck = .false.
    end if

    ! If we are setting a value to zero and nocheck is false then call delete_element instead
    if (.not.my_nocheck) then
       if (abs(val) < NUMERICAL_ZERO) then
          call sparse_delete_element(this,i,j)
          return
       end if
    end if

    ! Check the element we are setting is a valid one
    call sparse_check_bounds(this,i,j,'Sparse_Set_Element')
    
    !**** SHOULD WE USE A BINARY SEARCH TO FIND THE CORRECT INSERTION POINT? ****

    ! First find the place in table where the element will be inserted, or where it already exists
    table_pos = this%rows(i)
    element_written = .false.
    do n = this%rows(i), this%rows(i+1) - 1

       if (this%table%int(1,table_pos) > j) exit ! We've found the insertion point

       if (this%table%int(1,n) == j) then
          ! The element already exists, so just write the value
          this%table%real(1,n) = val
          element_written = .true.
          exit
       end if

       table_pos = table_pos + 1

    end do

    ! Now, either the element has been written (so we can return) or the insertion point is in table_pos
    if (element_written) return

    ! Append a blank element to 'table', insert the new value and shift the other entries down
    call append(this%table,(/0/),(/0.0_dp/))
    intpart = j
    realpart = val
    do n = table_pos, this%table%N
       !Make a temporary copy of the data that is already there
       tmp_intpart = this%table%int(1,n)
       tmp_realpart = this%table%real(1,n)
       !Copy in the new elements
       this%table%int(1,n) = intpart
       this%table%real(1,n) = realpart
       !Make the tmp_ variables the next to be written back to the table
       intpart = tmp_intpart
       realpart = tmp_realpart
    end do
    
    ! Lastly, all the row pointers for rows i+1 upwards need incrementing
    forall(n=i+1:this%Nrows+1) this%rows(n) = this%rows(n) + 1       

  end subroutine sparse_set_element

  !
  ! Search for the (i,j)th element and if it exists delete it (i.e. set it to zero)
  !
  subroutine sparse_delete_element(this,i,j)

    type(sparse), intent(inout) :: this
    integer,      intent(in)    :: i,j
    integer                     :: n, table_pos

    call sparse_check_bounds(this,i,j,'Sparse_Delete_Element')

    table_pos = 0
    do n = this%rows(i), this%rows(i+1) - 1
       if (this%table%int(1,n) > j) return    ! The element is already 0.0
       if (this%table%int(1,n) == j) then  
          ! We've found the element to be deleted
          table_pos = n
          exit
       end if
    end do

    if (table_pos==0) return !The element didn't exist

    ! Shift all the succeeding table entries back by one
    do n = table_pos, this%table%N - 1
       this%table%int(1,n)  = this%table%int(1,n+1)
       this%table%real(1,n) = this%table%real(1,n+1)
    end do
    
    ! Then delete the last entry
    call delete(this%table,this%table%N)

    ! Decrement the row pointers for rows i+1 upwards by one
    forall(n=i+1:this%Nrows+1) this%rows(n) = this%rows(n) - 1

  end subroutine sparse_delete_element

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X                                                      X
  !X BLAS                                                 X
  !X                                                      X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function sparse_mult_vector(this,vectin) result (vectout)
    real(dp),intent(IN)::vectin(:)
    type(sparse),intent(IN)::this
    real(dp)::vectout(this%Nrows)
    integer::i,j
#ifdef _MPI
    integer::mpi_size, mpi_rank, PRINT_ALWAYS
    include "mpif.h"
    real(dp)::vectout_all(this%Nrows)

    vectout_all = 0.0_dp
    call get_mpi_size_rank(MPI_COMM_WORLD, mpi_size, mpi_rank)
#endif
   
    call sparse_check(this, "sparse_mult_vector")
    if (this%Ncols .ne. size(vectin)) then
       write(line, *) 'sparse_mult_vector: mismatched matrix(',this%Nrows,',',this%Ncols,') and vector(',size(vectin),')!'
       call System_Abort(line)
    end if

    vectout=0.0_dp  
    do i=1,this%Nrows
#ifdef _MPI
       ! cycle loop if processor rank does not match
       if(mod(i, mpi_size) .ne. mpi_rank) cycle       
#endif
       do j= this%rows(i),this%rows(i+1)-1
          vectout(i)=vectout(i)+this%table%real(1,j)*vectin(this%table%int(1,j))
       end do
    end do

#ifdef _MPI
       ! collect mpi results
       call MPI_ALLREDUCE(vectout, vectout_all, &
            size(vectout), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, PRINT_ALWAYS)
       call abort_on_mpi_error(PRINT_ALWAYS, "Sparse .mult. vector: MPI_ALL_REDUCE()")
       vectout = vectout_all
#endif
  end function sparse_mult_vector


  function vector_mult_sparse(vectin,this) result(vectout)
    real(dp),intent(IN)::vectin(:)
    type(sparse),intent(IN)::this
    real(dp)::vectout(this%Ncols)
    integer::N,i,k,M
#ifdef _MPI
    integer::mpi_size, mpi_rank, PRINT_ALWAYS
    include "mpif.h"
    real(dp)::vectout_all(this%Ncols)

    vectout_all = 0.0_dp
    call get_mpi_size_rank(MPI_COMM_WORLD, mpi_size, mpi_rank)
#endif

    N=this%Nrows
    M=this%Ncols
    
    call sparse_check(this, "vector_mult_sparse")

    if (N .NE. size(vectin)) then
       call system_abort("sparse_product_vect:'mismatched vector and matrix!")
    end if  
  
 
    vectout=0.0_dp
    do i=1,N
#ifdef _MPI
       ! cycle loop if processor rank does not match
       if(mod(i, mpi_size) .ne. mpi_rank) cycle       
#endif
       do k= this%rows(i),this%rows(i+1)-1
          vectout(this%table%int(1,k))=vectout(this%table%int(1,k))+this%table%real(1,k)*vectin(i)
       end do
    end do
#ifdef _MPI
       ! collect mpi results
       call MPI_ALLREDUCE(vectout, vectout_all, &
            size(vectout), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, PRINT_ALWAYS)
       call abort_on_mpi_error(PRINT_ALWAYS, "Sparse .mult. vector: MPI_ALL_REDUCE()")
       vectout = vectout_all
#endif
  end function vector_mult_sparse


  !
  ! Tr(sparse*matrix)
  !
  function trace_sparse_mult_matrix(this,matrix) result(tr)
    type(sparse),intent(IN)::this
    real(dp),intent(IN)::matrix(:,:)
    real(dp)::tr
    integer::i,j
   
    call sparse_check(this,"sparce_tracemult_matrix") 
    if(this%Nrows .ne. size(matrix,1)) then 
       call system_abort("trace_sparse_mult_matrix: size of sparse and matrix mismatch ")
    end if

    tr = 0.0_dp

    do i = 1,size(matrix,1)
       do j = this%rows(i),this%rows(i+1)-1
          tr =  tr +this%table%real(1,j) * matrix(i,this%table%int(1,j))
       end do
    end do

  end function trace_sparse_mult_matrix

  !
  ! Sparse_CFCT : (Sparse matrix C) * (diagonal matrix F) * (Transpose of C)
  !
  function sparse_cfct(c,f) result(res)

    type(Sparse),           intent(in) :: c
    real(dp), dimension(:), intent(in) :: f
    type(Sparse)                       :: res
    integer                            :: i,j,k,n
    real(dp)                           :: val

    if (size(f) /= c%Ncols) call System_Abort('Sparse_CFCT: Number of elements in F is not equal to number of rows in C')

    ! Initialise the output matrix
    call Sparse_Init(res,c%Nrows,c%Nrows)
    
    ! Sparse matrix storage stores elements in the same row together, so make the loop over
    ! the columns go faster. Also, the matrix is symmetric so only half the working is needed
    do i = 1, c%Nrows

       do j = 1, i - 1 ! Just copy the pre-worked out values into the lower triangle
          val = sparse_element(res,j,i)
          if (abs(val) > NUMERICAL_ZERO) call sparse_set_element(res,i,j,val)
       end do

       do j = i, c%Nrows ! Actually work out the upper triangle of values        
          val = 0.0_dp
          do n = c%rows(i), c%rows(i+1)-1
             k = c%table%int(1,n)
             val = val + sparse_element(c,i,k) * f(k) * sparse_element(c,j,k)
          end do
          if (abs(val) > NUMERICAL_ZERO) call sparse_set_element(res,i,j,val)
       end do

    end do

  end function sparse_cfct

  function sparse_mult_sparse(a,b) result(res)

    type(sparse), intent(in) :: a,b
    type(sparse)             :: res
    !local variables
    integer                  :: i,j,n
    real(dp)                 :: val

    !Check sizes
    if (a%Ncols /= b%Nrows) call System_Abort('Sparse_Mult_Sparse: Argument 1 rows /= argument 2 columns')

    call Sparse_Init(res,a%Nrows,b%Ncols)

    do i = 1, a%Nrows
       do j = 1, b%Ncols
          val = 0.0_dp
          do n = a%rows(i), a%rows(i+1)-1
             val = val + a%table%real(1,n) * sparse_element(b, a%table%int(1,n), j)
          end do
          if (abs(val) > NUMERICAL_ZERO) call sparse_set_element(res,i,j,val)
       end do
    end do

  end function sparse_mult_sparse

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X                                                      X
  !X Printing                                             X
  !X                                                      X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine sparse_print(this,verbosity,out)
    type(sparse),intent(in)::this
    integer, optional::verbosity
    type(inoutput),intent(inout),optional::out

    call sparse_check(this,"sparse_print")
    write(line,'(a,i0,a,i0)')'Sparse matrix: Nrows = ',this%Nrows,' Ncols = ',this%Ncols
    call print(line,verbosity, out)  
    call print('row indices:', verbosity, out)
    call print(this%rows, verbosity, out)
    call print('  Column    Value', verbosity, out)
    call print(this%table, verbosity, out)
  end subroutine sparse_print

  subroutine sparse_print_full(this,verbosity,out)

    type(sparse), intent(in) :: this
    integer, optional :: verbosity
    type(inoutput),intent(in), optional :: out
    real(dp)                  :: row(this%Ncols)
    integer                   :: i, j
    character(13)             :: format

    call sparse_check(this,"sparse_print_full")
    if (this%Ncols >= 1e6) &
         call system_abort('Sparse_Print_Full: Tried to print a matrix with > 1M columns. Are you crazy?!?')
    write(format,'(a,i0,a)')'(',this%Ncols,'f14.5)'

    do j = 1, this%Nrows
       row = 0.0_dp
       do i = this%rows(j), this%rows(j+1) - 1
          row( this%table%int(1,i) ) = this%table%real(1,i)
       end do
       write(line,format) row
       call Print(line,verbosity,out)
    end do

  end subroutine sparse_print_full

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X                                                       X
  !X Miscellaneous checking routines                       X
  !X                                                       X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !
  ! Check if the element (i,j) is within the bounds of 'this'. caller is the calling routine
  !
  subroutine sparse_check_bounds(this,i,j,caller)

    type(Sparse), intent(in) :: this
    integer,      intent(in) :: i,j
    character(*), intent(in) :: caller

    if (i > this%Nrows .or. j > this%Ncols) then
       write(line,'(4(a,i0),a)') caller//': Element (',i,',',j,') is out of range (',&
                                                                   this%Nrows,',',this%Ncols,')'                     
       call System_Abort(line)
    end if

  end subroutine sparse_check_bounds

  !
  ! Check the input is a valid sparse matrix. caller is the calling routine.
  !
  subroutine sparse_check(this,caller)
   
    type(sparse), intent(in) :: this
    character(*)             :: caller
    integer                  :: m  

    if (.not. allocated(this%rows)) then
       write(line,'(a)')"In "//caller//":"
       call Print(line)
       call System_Abort("Sparse_check: sparse%rows is unallocated")
    end if

    m = size(this%rows) 
    if( (this%rows(m) /= this%table%N+1) .or. (m /= this%Nrows+1) ) then 
       write(line,'(a)')"In "//caller//":"
       call Print(line)
       call System_Abort("Sparse_check: sparse rows and table size mismatch")
    end if

  end subroutine sparse_check

  !
  ! Subroutine to test the sparse module is functioning properly
  ! NOTE: does not include tests for all features
  !
  subroutine sparse_test
    type(sparse)::sparse1,sparse3
    real(dp),allocatable::a(:,:),b(:,:),c(:,:),vect(:),vectout(:)
    allocate(a(5,5),b(5,5),c(5,5),vect(5))
    write(line,*) 'SPARSE TEST';call print(line)


    a=0.0_dp
    a(1,1)=3
    a(1,3)=1
    a(2,2)=4
    a(3,2)=7
    a(3,3)=5
    a(3,4)=9
    a(4,5)=2
    a(5,4)=6
    a(5,5)=5
    write(line,*) 'MATRIX A IS:';call print(line)
    call print(a)  
    c=a

    write(line,*) 'Now we assign it to a sparse';call print(line)
  !  call sparse_assign_matrix(sparse1,a)
    sparse3=a

    write(line,*) 'Now we assign a sparse to a sparse';call print(line)
  !  call sparse_assign_matrix(sparse1,a)
    sparse1=sparse3
  
    call sparse_check(sparse1,"This message should not be seen")
    write(line,*) 'Structure of the sparse';call print(line) 
    write(line,*) 'And its content';call print(line)
    call print(sparse1)
    write(line,*) 'Now we assign it back to a matrix';call print(line)
     call matrix_assign_sparse(b,sparse1)
    b=sparse1
    call print(b)
    if ( a .FNE. b) then
       write(line,*) 'Something is wrong: a .NE. b';call print(line)  
    else
       write(line,*) 'Everything seems alright!';call print(line)
    end if

    write(line,*) 'Now we multiply our sparse with a unity vector';call print(line)
    vect=1
    call print(vect) 
    allocate(vectout(size(vect)))
    vectout = sparse_mult_vector(sparse1,vect)
    write(line,*) 'Result should be 4 4 21 2 11';call print(line)
    call print(vectout) 

     write(line,*) 'Now we multiply the following vector for the sparse';call print(line)
    vect(1)=5
    vect(2)=2
    vect(3)=3
    vect(4)=2 
    vect(5)=1
    
    call print(vect) 
    vectout = vector_mult_sparse(vect,sparse1)
    write(line,*) 'Result should be 15 29 20 33 9';call print(line)
    call print(vectout) 


    write(line,*) 'Lets check if it is working with empty rows ';call print(line)
    a(2,2)=0.0  ! now row 2 has no elements
    call print(a)  
    sparse1=a
    call print(sparse1)
    call matrix_assign_sparse(b,sparse1)
    call print(b)
    if ( a .FNE. b) then
       write(line,*) 'Something is wrong: a .NE. b';call print(line)  
    else
       write(line,*) 'Everything seems alright!';call print(line)
    end if

    write(line,*) 'Now we check an empty matrix';call print(line)
    a=0.0
    call print(a) 
    call sparse_assign_matrix(sparse1,a)
    call print(sparse1)
    call matrix_assign_sparse(b,sparse1)
    call print(b)
    if ( a .FNE. b) then
       write(line,*) 'Something is wrong: a .NE. b';call print(line)  
    else
       write(line,*) 'Everything seems alright!';call print(line)
    end if

    write(line,*) 'Now we check tracemult';call print(line)
    write(line,*) 'The matrix is a=1';call print(line)      
    a=1
    call sparse_assign_matrix(sparse1,c)
    write(line,*) 'The following should be 42:',trace_sparse_mult_matrix(sparse1,a)  
    call print(line) 
    write(line,*) 'Now we check print_full'; call print(line)
    call sparse_print_full(sparse1)

  end subroutine sparse_test

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

end module sparse_module
