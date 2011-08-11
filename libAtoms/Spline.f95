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
!X  Spline module
!X  
!%  Defines derived type for a spline, and associated functions.
!%  Code partially based on Numerical Recipes.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module spline_module
  use linearalgebra_module
  use table_module
  implicit none
  SAVE


  type Spline
     integer                             ::n        !% Number of knot points
     real(dp), allocatable, dimension(:) ::x        !% Knot positions
     real(dp), allocatable, dimension(:) ::y        !% Function values
     real(dp), allocatable, dimension(:) ::y2       !% Second derivative
     real(dp)                            ::yp1,ypn  !% Endpoint derivatives
     ! whether y2 has been initialised or not
     logical                             ::y2_initialised = .false. 
  end type Spline

  interface initialise
     module procedure spline_init
  end interface

  interface finalise
     module procedure spline_finalise
  end interface
  
  interface print
     module procedure spline_print
  end interface

contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X spline_init(this, x, y, yp1, ypn)
  !X
  !% Initialises a spline with given x and y values (and endpoint derivatives)
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  subroutine spline_init(this, x, y, yp1, ypn)
    type(spline), intent(inout)::this
    real(dp), dimension(:) :: x !% Knot points
    real(dp), dimension(:) :: y !%Values of spline at the knot points
    real(dp)::yp1 !% Derivative of the spline at 'x(1)'
    real(dp)::ypn !% Derivative of the spline at 'x(n)'


    call check_size('Y',y,size(x),'Spline_Init')

    ! first deallocate arrays if they have been allocated before
    if(allocated(this%x))  deallocate(this%x)
    if(allocated(this%y))  deallocate(this%y)
    if(allocated(this%y2)) deallocate(this%y2)

    ! allocate new sizes
    this%n = size(x);
    allocate(this%x(this%n))
    allocate(this%y(this%n))
    allocate(this%y2(this%n))

    ! copy data
    this%x = x
    this%y = y
    call sort_array(this%x, r_data=this%y)
    this%yp1 = yp1
    this%ypn = ypn

    ! compute y2
    call spline_y2calc(this)
  end subroutine spline_init

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X spline_finalise(this)
  !X
  !% Deallocates a spline 
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine spline_finalise(this)
    type(spline), intent(inout)::this
    if(allocated(this%x))  deallocate(this%x)
    if(allocated(this%y))  deallocate(this%y)
    if(allocated(this%y2)) deallocate(this%y2)
    this%y2_initialised = .false.
  end subroutine spline_finalise

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X spline_y2calc(this)
  !X
  !% Takes a spline with $x$ and $y$ values, and computes the second derivates
  !% this must be done before interpolation.
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  subroutine spline_y2calc(this)
    type(spline),intent(inout)::this
    real(dp), allocatable, dimension(:) ::u  ! local temporary
    real(dp)::sig,p,un,qn
    integer::i,k, n

    if(.not.allocated(this%x)) call system_abort("spline_y2calc: x not allocated!")
    if(.not.allocated(this%y)) call system_abort("spline_y2calc: y not allocated!")

    n = this%n

    if(.not.allocated(this%y2)) allocate(this%y2(n))

    allocate(u(n))

    if(this%yp1> 0.99e30_dp) then
       ! natural spline with zero second derivative at first point
       this%y2(1) = 0.0_dp
       u(1) = 0.0_dp
    else
       this%y2(1) = -0.5_dp

       u(1) =  (3.0_dp/(this%x(2)-this%x(1)))* &
            ((this%y(2)-this%y(1))/(this%x(2)-this%x(1))-this%yp1)
    end if

    do i=2,n-1
       sig = (this%x(i)-this%x(i-1))/(this%x(i+1)-this%x(i-1))
       p = sig*this%y2(i-1)+2.0_dp
       this%y2(i) = (sig-1.0)/p
       u(i) = (this%y(i+1)-this%y(i))/(this%x(i+1)-this%x(i)) - &
            (this%y(i)-this%y(i-1))/(this%x(i)-this%x(i-1))
       u(i) = (6.0*u(i)/(this%x(i+1)-this%x(i-1))-sig*u(i-1))/p
    end do

    if(this%ypn > 0.99e30_dp)then
       ! natural spline with zero second derivative at last point
       qn = 0.0_dp 
       un = 0.0_dp
    else
       qn = 0.5_dp
       un = (3.0_dp/(this%x(n)-this%x(n-1)))* &
            (this%ypn-(this%y(n)-this%y(n-1))/(this%x(n)-this%x(n-1)))
    end if

    this%y2(n) = (un-qn*u(n-1))/(qn*this%y2(n-1)+1.0_dp)

    do k=n-1,1,-1
       this%y2(k) = this%y2(k)*this%y2(k+1)+u(k)
    end do

    ! set the logical flag
    this%y2_initialised = .true.

    deallocate(u) ! free temporary

    if( this%yp1 > 0.99e30_dp ) then
       this%yp1 = spline_deriv(this,this%x(1))
    endif

    if( this%ypn > 0.99e30_dp ) then
       this%ypn = spline_deriv(this,this%x(n))
    endif

  end subroutine spline_y2calc

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X y = spline_value(this, x)
  !X
  !% Interpolate the spline for a given $x$ point
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  function spline_value(this,x) result(y)
    type(spline)::this
    real(dp)::x,h,a,b,y
    integer::klo,khi,k, n

    if(.NOT.this%y2_initialised) then
       if(allocated(this%x).and.allocated(this%y)) then
          call spline_y2calc(this)
       else
          call system_abort("spline_value: spline has not been initialised")
       end if
    end if

    n = this%n
    if( x < this%x(1) ) then
       y = this%y(1) + (x - this%x(1))*this%yp1
    elseif( x > this%x(n) ) then
       y = this%y(n) + (x - this%x(n))*this%ypn
    else
       klo = 1
       khi = this%n

       do  while(khi-klo > 1)
          k = (khi+klo)/2
          if(this%x(k) > x) then 
             khi = k
          else 
             klo = k
          end if
       end do

       h = this%x(khi)-this%x(klo)
       if(h .EQ. 0.0_dp) then
          call system_abort("spline_interpolate: h=0!!!")
       end if

       a = (this%x(khi)-x)/h
       b = (x-this%x(klo))/h
       y = a*this%y(klo)+b*this%y(khi)+((a*a*a-a)*this%y2(klo)+(b*b*b-b)*this%y2(khi))*(h*h)/6.0
    endif

  end function spline_value

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X dy = spline_deriv(this, x)
  !X
  !% Interpolate the derivative of the spline at the given $x$ point
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  function spline_deriv(this,x) result(dy)
    type(spline)::this
    real(dp)::x,h,a,b,dy
    integer::klo,khi,k,n

    if(.NOT.this%y2_initialised) then
       if(allocated(this%x).and.allocated(this%y)) then
          call spline_y2calc(this)
       else
          call system_abort("spline_deriv: spline has not been initialised")
       end if
    end if

    n = this%n

    if( x < this%x(1) ) then
       dy = this%yp1
    elseif( x > this%x(n) ) then
       dy = this%ypn
    else
       klo = 1
       khi = n

       do  while(khi-klo > 1)
          k = (khi+klo)/2
          if(this%x(k) > x) then 
             khi = k
          else 
             klo = k
          end if
       end do

       h = this%x(khi)-this%x(klo)
       if(h .EQ. 0.0) then
          call system_abort("spline_deriv: h=0!!!")
       end if

       a = (this%x(khi)-x)/h
       b = (x-this%x(klo))/h
       dy = (this%y(khi)-this%y(klo))/h+((3.0_dp*b*b-1.0_dp)*this%y2(khi)-(3.0_dp*a*a-1.0_dp)*this%y2(klo))*h/6.0_dp
    endif

  end function spline_deriv

  function spline_nintegrate(this, x0, x, n_per_knot) result (y_int)
     type(spline), intent(in) :: this
     real(dp), intent(in) :: x0, x
     integer, intent(in), optional :: n_per_knot
     real(dp) :: y_int

     integer :: use_n, total_n, i
     real(dp), allocatable :: int_x(:), int_y(:)

     use_n = optional_default(10, n_per_knot)
     total_n = use_n*abs(x-x0)/minval(abs(this%x(2:this%n)-this%x(1:this%n-1)))

     allocate(int_x(total_n))
     allocate(int_y(total_n))

     do i=1, total_n
        int_x(i) = x0+real(i-1,dp)*(x-x0)/real(total_n-1,dp)
	int_y(i) = spline_value(this, int_x(i))
     end do

     y_int = TrapezoidIntegral(int_x, int_y)
     deallocate(int_x)
     deallocate(int_y)

  end function spline_nintegrate

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X spline_compute_matrices(this)
  !X
  !% The y2 of the spline can also be computed using
  !% linear algebra. Following Numerical Recipes, 
  !%
  !% \begin{displaymath}
  !% A Y'' = BY + C
  !% \end{displaymath}
  !% where $A$, $B$ are matrices, $C$ is a vector
  !%  so
  !% \begin{displaymath}
  !% Y'' = A^{-1} B Y + A^{-1} C
  !% \end{displaymath}
  !%
  !% This subroutine computes these matrices, and stores them in
  !% 'spline%y2_matrix1 = inv(A)*B' and 'spline%y2_matrix0 = inv(A)*C'
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine spline_compute_matrices(this, y2_matrix1, y2_matrix0)
    type(spline), intent(in)::this
    real(dp), dimension(:,:), intent(out)::y2_matrix1
    real(dp), dimension(:), intent(out)::y2_matrix0
    real(dp)::A(this%n,this%n)
    real(dp)::B(this%n,this%n),C(this%n)
    real(dp)::Ainv(this%n,this%n)
    integer::i,n

    n = this%n

    ! Set matrix elements that correspond to
    ! conditions on first derivative at x1 and xn
    ! equation 3.3.5 of Numerical Recipes in c++ second edition pag. 117
    ! evaluated at x=x1 and x=xn
    A=0.0_dp
    B=0.0_dp
    C=0.0_dp
    A(1,1)  = -(1.0_dp/3.0_dp)*(this%x(2)-this%x(1))
    A(1,2)  = -(1.0_dp/6.0_dp)*(this%x(2)-this%x(1))
    A(2,n-1)=  (1.0_dp/6.0_dp)*(this%x(n)-this%x(n-1))
    A(2,n)  =  (1.0_dp/3.0_dp)*(this%x(n)-this%x(n-1))

    B(1,1)  =  1.0_dp/(this%x(2)-this%x(1))
    B(1,2)  = -1.0_dp/(this%x(2)-this%x(1))
    B(2,n-1)=  1.0_dp/(this%x(n)-this%x(n-1))
    B(2,n)  = -1.0_dp/(this%x(n)-this%x(n-1))

    C(1) = this%yp1
    C(n) = this%ypn

    ! Set the matrix elements corresponds to the 
    ! equations 3.3.7 of page 118 in Numerical Recipes

    do I=2,n-1
       A(I+1,I-1)  =(this%x(i)  -this%x(i-1))/6.0_dp
       A(I+1,I)    =(this%x(i+1)-this%x(i-1))/3.0_dp
       A(I+1,I+1)  =(this%x(i+1)-this%x(i)  )/6.0_dp

       B(I+1,I-1)  =  1.0_dp/(this%x(i)  -this%x(i-1))
       B(I+1,I)    = -1.0_dp/(this%x(i+1)-this%x(i))-  &
                      1.0_dp/(this%x(i)  -this%x(i-1))
       B(I+1,I+1)  =  1.0_dp/(this%x(i+1)-this%x(i))
    end do

    call matrix_inverse(A,Ainv)

    y2_matrix1 = Ainv .mult. B

    y2_matrix0 = Ainv .mult. C

  end subroutine spline_compute_matrices


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X spline_print(this)
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine spline_print(this, file)
    type(Spline), intent(in)::this
    type(Inoutput), optional, target, intent(inout) :: file

    real(dp), allocatable :: y2m1(:,:), y2m0(:), y2(:)

    call print("Spline ", file=file)

    call print(this%n, file=file)
    call print(reshape( (/this%x,this%y/), (/this%n,2/)), file=file)

    call verbosity_push_decrement()
    if(this%y2_initialised) then
      call print(this%y2, file=file)
    else

      allocate(y2m1(this%n, this%n))
      allocate(y2m0(this%n))
      allocate(y2(this%n))

      call spline_compute_matrices(this, y2m1, y2m0)
      y2 = (y2m1 .mult. this%y)+y2m0
      call print(y2, file=file)

      deallocate(y2m1)
      deallocate(y2m0)
      deallocate(y2)
    endif

    call verbosity_pop()
  end subroutine spline_print

  function min_knot(this)
    type(Spline), intent(in) :: this
    real(dp) :: min_knot

    min_knot = this%x(1)
  end function min_knot

  function max_knot(this)
    type(Spline), intent(in) :: this
    real(dp) :: max_knot

    max_knot = this%x(size(this%x))
  end function max_knot

end module spline_module



