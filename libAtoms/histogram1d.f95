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
!X Histogram helper functions in 1D
!X
!X Computes periodic and non-periodic histograms within certain bounds.
!X Interpolation can be chosen from
!X   INTERP_LINEAR:  Linear interpolation between neigboring grid
!X                   points
!X   INTERP_LUCY:    Lucy function interpolation between neighboring
!X                   grid points
!X TODO:
!X   INTERP_NONE:    No interpolation, simply sort into bins
!X
!X Provides routine for MPI communication.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

#define ASSERT(x)

module histogram1d_module
  use Error_module
  use System_module
  use Units_module
  use MPI_context_module

  private

  public :: INTERP_LINEAR, INTERP_LUCY
  integer, parameter  :: INTERP_LINEAR  = 0
  integer, parameter  :: INTERP_LUCY    = 1

  public :: Histogram1D
  type Histogram1D

     !
     ! General stuff
     !

     integer            :: interp     ! INTERP_LINEAR, INTERP_LUCY

     integer            :: n = -1     ! number of bins
     real(DP)           :: min_b      ! minimum value
     real(DP)           :: max_b      ! maximun value
     real(DP)           :: db         ! difference
     real(DP)           :: dbin       ! bin size

     logical            :: periodic

     !
     ! For smoothing function
     !

     real(DP)           :: sigma

     real(DP)           :: fac1
!     real(DP)           :: sigma_sq

     !
     ! Values
     !

     real(DP), allocatable  :: x(:)       ! x-values

     real(DP), allocatable  :: h(:)       ! values
     real(DP), allocatable  :: h_sq(:)    ! second moment

     real(DP), allocatable  :: h1(:)      ! histogram with norm 1 (number of values which have been added to each bin)

     integer                :: ng
     real(DP), pointer      :: smoothing_func(:)
   
  endtype Histogram1D

  !
  ! Interface definition
  !

  public :: initialise
  interface initialise
     module procedure histogram1d_initialise, histogram1d_initialise_from_histogram
  endinterface

  public :: finalise
  interface finalise
     module procedure histogram1d_finalise
  endinterface

  public :: clear
  interface clear
     module procedure histogram1d_clear
  endinterface

  public :: set_bounds
  interface set_bounds
     module procedure histogram1d_set_bounds
  endinterface

  public :: write
  interface write
     module procedure histogram1d_write, histogram1d_write_mult
     module procedure histogram1d_write_character_fn, histogram1d_write_mult_character_fn
  endinterface

  public :: add
  interface add
     module procedure histogram1d_add, histogram1d_add_vals, histogram1d_add_vals_norms
     module procedure histogram1d_add_vals_mask, histogram1d_add_vals_norm_mask, histogram1d_add_vals_norms_mask
     module procedure histogram1d_add_histogram
  endinterface

  public :: add_range
  interface add_range
     module procedure histogram1d_add_range, histogram1d_add_range_vals, histogram1d_add_range_vals_norms
  endinterface

  public :: average
  interface average
     module procedure histogram1d_average
  endinterface

  public :: entropy
  interface entropy
     module procedure histogram1d_entropy
  endinterface

  public :: expectation_value
  interface expectation_value
     module procedure histogram1d_expectation_value
  endinterface

  public :: mul
  interface mul
     module procedure histogram1d_mul, histograms_mul, histograms2_mul, histogram1d_mul_vals
  endinterface

  public :: div
  interface div
     module procedure histogram1d_div, histograms_div, histograms2_div, histogram1d_div_vals
  endinterface

  public :: normalize
  interface normalize
     module procedure histogram1d_normalize
  endinterface

  public :: reduce
  interface reduce
     module procedure histogram1d_reduce
  endinterface

  public :: smooth
  interface smooth
     module procedure histogram1d_smooth
  endinterface

  public :: sum_in_place
  interface sum_in_place
     module procedure histogram1d_sum_in_place
  endinterface

! Memory estimation not yet implemented in libAtoms
#if 0
  interface log_memory_estimate
     module procedure log_memory_estimate_histogram, log_memory_estimate_histogram2, log_memory_estimate_histogram3
  endinterface
#endif

contains

  !% Initialize the histogram
  !% Default is linear interpolation, non-periodic
  elemental subroutine histogram1d_initialise(this, n, min_b, max_b, sigma, &
       periodic)
    implicit none

    type(Histogram1D), intent(out)  :: this
    integer, intent(in)             :: n
    real(DP), intent(in)            :: min_b
    real(DP), intent(in)            :: max_b
    real(DP), intent(in), optional  :: sigma
    logical, intent(in), optional   :: periodic

    ! ---

    call finalise(this)

    this%n  = n

    if (present(sigma) .and. sigma > 0.0_DP) then

       this%interp    = INTERP_LUCY

       this%sigma     = sigma

!       this%fac1      = 1.0_DP/(sqrt(2*PI)*this%sigma)
!       this%sigma_sq  = 2*sigma**2
       this%fac1      = 5/(4*this%sigma)

    else

       this%interp    = INTERP_LINEAR

    endif

    if (present(periodic)) then
       this%periodic  = periodic
    else
       this%periodic  = .false.
    endif

    this%smoothing_func => NULL()

    allocate(this%x(n))
    allocate(this%h(n))
    allocate(this%h_sq(n))
    allocate(this%h1(n))

    this%min_b = min_b - 1.0_DP
    this%max_b = max_b - 1.0_DP

    call histogram1d_set_bounds(this, min_b, max_b)

    call histogram1d_clear(this)

  endsubroutine histogram1d_initialise


  !% Initialize the histogram
  elemental subroutine histogram1d_initialise_from_histogram(this, that)
    implicit none

    type(Histogram1D), intent(out)  :: this
    type(Histogram1D), intent(in )  :: that

    ! ---

    call finalise(this)

    this%n         = that%n

    this%interp    = that%interp
    this%sigma     = that%sigma

    if (this%interp == INTERP_LUCY) then
       this%fac1      = 5/(4*this%sigma)
    endif

    this%periodic  = that%periodic

    this%smoothing_func => NULL()

    allocate(this%x(this%n))
    allocate(this%h(this%n))
    allocate(this%h_sq(this%n))
    allocate(this%h1(this%n))

    this%min_b  = that%min_b - 1.0_DP
    this%max_b  = that%max_b - 1.0_DP

    call histogram1d_set_bounds(this, that%min_b, that%max_b)

    call histogram1d_clear(this)
    
  endsubroutine histogram1d_initialise_from_histogram


  !% Clear histogram
  elemental subroutine histogram1d_clear(this)
    implicit none

    type(Histogram1D), intent(inout)  :: this

    ! ---

    this%h1    =  0.0_DP
    this%h     = 0.0_DP
    this%h_sq  = 0.0_DP

  endsubroutine histogram1d_clear


  !% Set histogram bounds
  elemental subroutine histogram1d_set_bounds(this, min_b, max_b)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: min_b
    real(DP), intent(in)              :: max_b

    ! ---

    integer   :: i
    real(DP)  :: x

    ! ---

    if (this%min_b /= min_b .or. this%max_b /= max_b) then

       this%min_b  = min_b
       this%max_b  = max_b
       this%db     = max_b - min_b

       if (this%periodic) then
          this%dbin   = this%db/this%n
       else
          this%dbin   = this%db/(this%n-1)
       endif

       do i = 1, this%n
          this%x(i)  = this%min_b + this%dbin*(i-0.5_DP)
       enddo

       !
       ! Tabulate the smoothing function
       !

       if (this%interp == INTERP_LUCY) then

          this%ng  = floor(this%sigma/this%dbin)
          
          if (.not. associated(this%smoothing_func)) then
             allocate(this%smoothing_func(-this%ng:this%ng))
          else if (ubound(this%smoothing_func, 1) /= this%ng) then
             deallocate(this%smoothing_func)
             allocate(this%smoothing_func(-this%ng:this%ng))
          endif

          do i = -this%ng, this%ng
!          this%smoothing_func(i)  = this%fac1 * exp(-((this%dbin*i)**2)/this%sigma_sq)
             x = abs( ( this%dbin*i )/this%sigma )
             this%smoothing_func(i)  = this%fac1 * (1+3*x)*(1-x)**3
          enddo

       endif

    endif

  endsubroutine histogram1d_set_bounds


  !% Delete a histogram table
  elemental subroutine histogram1d_finalise(this)
    implicit none

    type(Histogram1D), intent(inout)  :: this

    ! ---

    if (allocated(this%x)) then
       deallocate(this%x)
    endif
    if (allocated(this%h)) then
       deallocate(this%h)
    endif
    if (allocated(this%h_sq)) then
       deallocate(this%h_sq)
    endif
    if (allocated(this%h1)) then
       deallocate(this%h1)
    endif

    if (associated(this%smoothing_func)) then
       deallocate(this%smoothing_func)
    endif

  endsubroutine histogram1d_finalise


  !% Output the histogram table to a file (by file unit)
  subroutine histogram1d_write_mult(this, file, xvalues, header, error)
    implicit none

    type(Histogram1D), intent(in)       :: this(:)
    type(InOutput), intent(in)          :: file
    logical, intent(in), optional       :: xvalues
    character(*), intent(in), optional  :: header(lbound(this, 1):ubound(this, 1))
    integer, intent(out), optional      :: error

    ! ---

    integer        :: i, j
    character(80)  :: fmt
    character(20)  :: ext_header(lbound(this, 1):ubound(this, 1))

    integer        :: n
    real(DP)       :: min_b, max_b, dbin, y(lbound(this, 1):ubound(this, 1))
    
    ! ---

    INIT_ERROR(error)

    n      = this(lbound(this, 1))%n
    min_b  = this(lbound(this, 1))%min_b
    max_b  = this(lbound(this, 1))%max_b
    dbin   = this(lbound(this, 1))%dbin

    do i = lbound(this, 1)+1, ubound(this, 1)
       if (this(i)%n /= n) then
          RAISE_ERROR("Number of histogram bins do not match.", error)
       endif

       if (this(i)%min_b /= min_b) then
          RAISE_ERROR("*min_b*s do not match.", error)
       endif

       if (this(i)%max_b /= max_b) then
          RAISE_ERROR("*max_b*s do not match.", error)
       endif
    enddo

    if (present(header)) then
       do i = lbound(this, 1), ubound(this, 1)
          write (fmt, '(I2.2)')  i+2
          ext_header(i)  = trim(fmt) // ":" // trim(header(i))
          ext_header(i)  = adjustr(ext_header(i))
       enddo

       write (fmt, '(A,I4.4,A)')  "(A5,5X,A20,", ubound(this, 1)-lbound(this, 1)+1, "A20)"
       write (file%unit, fmt)  "#01:i", "02:x", ext_header
    endif

    if (present(xvalues) .and. .not. xvalues) then
       write (fmt, '(A,I4.4,A)')  "(A,", n, "ES20.10)"
       write (file%unit, trim(fmt))  "# ", this(lbound(this, 1))%x
       write (fmt, '(A,I4.4,A)')  "(I10,", ubound(this, 1)-lbound(this, 1)+1, "ES20.10)"
    else
       write (fmt, '(A,I4.4,A)')  "(I10,", ubound(this, 1)-lbound(this, 1)+2, "ES20.10)"
    endif

    do i = 1, n
       do j = lbound(this, 1), ubound(this, 1)
          y(j) = this(j)%h(i)
       enddo

       if (present(xvalues) .and. .not. xvalues) then
          write (file%unit, trim(fmt))  i, y
       else
          write (file%unit, trim(fmt))  i, this(lbound(this, 1))%x(i), y
       endif
    enddo

  endsubroutine histogram1d_write_mult


  !% Output the histogram table to a file (by file name)
  subroutine histogram1d_write_mult_character_fn(this, fn, xvalues, header, error)
    implicit none

    type(Histogram1D), intent(in)       :: this(:)
    character(*), intent(in)            :: fn
    logical, intent(in), optional       :: xvalues
    character(*), intent(in), optional  :: header(lbound(this, 1):ubound(this, 1))
    integer, intent(out), optional      :: error

    ! ---

    type(InOutput)  :: file

    ! ---

    INIT_ERROR(error)

    call initialise(file, fn, action=OUTPUT, error=error)
    PASS_ERROR(error)
    call histogram1d_write_mult(this, file, xvalues, header, error=error)
    call finalise(file)
    PASS_ERROR(error)

  endsubroutine histogram1d_write_mult_character_fn


  !% Output the histogram table to a file (by file unit)
  subroutine histogram1d_write(this, file, xvalues, header, error)
    implicit none

    type(Histogram1D), intent(in)       :: this
    type(InOutput), intent(in)          :: file
    logical, intent(in), optional       :: xvalues
    character(*), intent(in), optional  :: header
    integer, intent(out), optional      :: error

    ! ---

    INIT_ERROR(error)

    if (present(header)) then
       call histogram1d_write_mult((/ this /), file, xvalues, (/ header /), error=error)
       PASS_ERROR(error)
    else
       call histogram1d_write_mult((/ this /), file, xvalues, error=error)
       PASS_ERROR(error)
    endif

  endsubroutine histogram1d_write


  !% Output the histogram table to a file (by file name)
  subroutine histogram1d_write_character_fn(this, fn, xvalues, header, error)
    implicit none

    type(Histogram1D), intent(in)       :: this
    character(*), intent(in)            :: fn
    logical, intent(in), optional       :: xvalues
    character(*), intent(in), optional  :: header
    integer, intent(out), optional      :: error

    ! ---

    INIT_ERROR(error)

    if (present(header)) then
       call histogram1d_write_mult_character_fn((/ this /), fn, xvalues, (/ header /), error=error)
       PASS_ERROR(error)
    else
       call histogram1d_write_mult_character_fn((/ this /), fn, xvalues, error=error)
       PASS_ERROR(error)
    endif

  endsubroutine histogram1d_write_character_fn


  !% Add a value to the histogram with linear interpolation
  subroutine histogram1d_add_linear(this, val, norm)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: val
    real(DP), intent(in)              :: norm

    ! ---

    integer   :: i1, i2
    real(DP)  :: d1, d2

    ! ---

    i1  = int(floor((val-this%min_b)/this%dbin))+1
    i2  = i1+1

    d1  = ( (i2-1) - (val-this%min_b)/this%dbin )/this%dbin
    d2  = ( (val-this%min_b)/this%dbin - (i1-1) )/this%dbin

    if ((i1 >= 0 .and. i1 <= this%n) .or. this%periodic) then
       ASSERT(d1 >= 0.0_DP .and. d1*this%dbin <= 1.0_DP)

       if (this%periodic) then
          i1  = modulo(i1-1, this%n)+1
       else
          if (i1 == 0) then
             i1 = 1
          endif
       endif

       this%h1(i1)    = this%h1(i1)   + d1
       this%h(i1)     = this%h(i1)    + d1*norm
       this%h_sq(i1)  = this%h_sq(i1) + d1*norm**2
    endif

    if ((i2 >= 1 .and. i2 <= this%n+1) .or. this%periodic) then
       ASSERT(d2 >= 0.0_DP .and. d2*this%dbin <= 1.0_DP)

       if (this%periodic) then
          i2  = modulo(i2-1, this%n)+1
       else
          if (i2 == this%n+1) then
             i2 = this%n
          endif
       endif

       this%h1(i2)    = this%h1(i2)   + d2
       this%h(i2)     = this%h(i2)    + d2*norm
       this%h_sq(i2)  = this%h_sq(i2) + d2*norm**2
    endif

  endsubroutine histogram1d_add_linear


  !% Add multiple values to the histogram with linear interpolation
  subroutine histogram1d_add_linear_vals(this, vals, norm)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: vals(:)
    real(DP), intent(in)              :: norm

    ! ---

    integer   :: i1(lbound(vals, 1):ubound(vals, 1)), i2(lbound(vals, 1):ubound(vals, 1))
    real(DP)  :: d1(lbound(vals, 1):ubound(vals, 1)), d2(lbound(vals, 1):ubound(vals, 1))

    integer   :: i

    ! ---

    i1  = int(floor((vals-this%min_b)/this%dbin))+1
    i2  = i1+1

    d1  = ( (i2-1) - (vals-this%min_b)/this%dbin )/this%dbin
    d2  = ( (vals-this%min_b)/this%dbin - (i1-1) )/this%dbin

    ASSERT(all(d1 >= 0.0_DP .and. d1*this%dbin <= 1.0_DP))
    ASSERT(all(d2 >= 0.0_DP .and. d2*this%dbin <= 1.0_DP))

    if (this%periodic) then

       i1  = modulo(i1-1, this%n)+1
       i2  = modulo(i2-1, this%n)+1

       do i = lbound(vals, 1), ubound(vals, 1)
          this%h1(i1(i))    = this%h1(i1(i))   + d1(i)
          this%h1(i2(i))    = this%h1(i2(i))   + d2(i)

          this%h(i1(i))     = this%h(i1(i))    + d1(i)*norm
          this%h(i2(i))     = this%h(i2(i))    + d2(i)*norm

          this%h_sq(i1(i))  = this%h_sq(i1(i)) + d1(i)*norm**2
          this%h_sq(i2(i))  = this%h_sq(i2(i)) + d2(i)*norm**2
       enddo

    else

       do i = lbound(vals, 1), ubound(vals, 1)
          if (i1(i) >= 1 .and. i1(i) <= this%n) then
             this%h1(i1(i))    = this%h1(i1(i))   + d1(i)
             this%h(i1(i))     = this%h(i1(i))    + d1(i)*norm
             this%h_sq(i1(i))  = this%h_sq(i1(i)) + d1(i)*norm**2
          else if (i1(i) == 0) then
             this%h1(1)    = this%h1(1)   + d1(i)
             this%h(1)     = this%h(1)    + d1(i)*norm
             this%h_sq(1)  = this%h_sq(1) + d1(i)*norm**2
          endif

          if (i2(i) >= 1 .and. i2(i) <= this%n) then
             this%h1(i2(i))    = this%h1(i2(i))   + d2(i)
             this%h(i2(i))     = this%h(i2(i))    + d2(i)*norm
             this%h_sq(i2(i))  = this%h_sq(i2(i)) + d2(i)*norm**2
          else if (i2(i) == this%n+1) then
             this%h1(this%n)    = this%h1(this%n)   + d2(i)
             this%h(this%n)     = this%h(this%n)    + d2(i)*norm
             this%h_sq(this%n)  = this%h_sq(this%n) + d2(i)*norm**2
          endif
       enddo

    endif

  endsubroutine histogram1d_add_linear_vals


  !% Add multiple values to the histogram with linear interpolation
  subroutine histogram1d_add_linear_vals_norms(this, vals, norms)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: vals(:)
    real(DP), intent(in)              :: norms(:)

    ! ---

    integer   :: i1(lbound(vals, 1):ubound(vals, 1)), i2(lbound(vals, 1):ubound(vals, 1))
    real(DP)  :: d1(lbound(vals, 1):ubound(vals, 1)), d2(lbound(vals, 1):ubound(vals, 1))

    real(DP)  :: d1n(lbound(vals, 1):ubound(vals, 1)), d2n(lbound(vals, 1):ubound(vals, 1))
    real(DP)  :: d1nn(lbound(vals, 1):ubound(vals, 1)), d2nn(lbound(vals, 1):ubound(vals, 1))

    integer   :: i

    ! ---

    i1  = int(floor((vals-this%min_b)/this%dbin))+1
    i2  = i1+1

    d1  = ( (i2-1) - (vals-this%min_b)/this%dbin )/this%dbin
    d2  = ( (vals-this%min_b)/this%dbin - (i1-1) )/this%dbin

    ASSERT(all(d1 >= 0.0_DP .and. d1*this%dbin <= 1.0_DP))
    ASSERT(all(d2 >= 0.0_DP .and. d2*this%dbin <= 1.0_DP))

    d1n   = d1*norms
    d2n   = d2*norms

    d1nn  = d1n*norms
    d2nn  = d2n*norms

    if (this%periodic) then

       i1  = modulo(i1-1, this%n)+1
       i2  = modulo(i2-1, this%n)+1

       do i = lbound(vals, 1), ubound(vals, 1)
          this%h1(i1(i))    = this%h1(i1(i))  + d1(i)
          this%h1(i2(i))    = this%h1(i2(i))  + d2(i)

          this%h(i1(i))     = this%h(i1(i))   + d1n(i)
          this%h(i2(i))     = this%h(i2(i))   + d2n(i)

          this%h_sq(i1(i))  = this%h_sq(i1(i)) + d1nn(i)
          this%h_sq(i2(i))  = this%h_sq(i2(i)) + d2nn(i)
       enddo

    else

       do i = lbound(vals, 1), ubound(vals, 1)
          if (i1(i) >= 1 .and. i1(i) <= this%n) then
             this%h1(i1(i))    = this%h1(i1(i))   + d1(i)
             this%h(i1(i))     = this%h(i1(i))    + d1n(i)
             this%h_sq(i1(i))  = this%h_sq(i1(i)) + d1nn(i)
          else if (i1(i) == 0) then
             this%h1(1)    = this%h1(1)   + d1(i)
             this%h(1)     = this%h(1)    + d1n(i)
             this%h_sq(1)  = this%h_sq(1) + d1nn(i)
          endif

          if (i2(i) >= 1 .and. i2(i) <= this%n) then
             this%h1(i2(i))    = this%h1(i2(i))   + d2(i)
             this%h(i2(i))     = this%h(i2(i))    + d2n(i)
             this%h_sq(i2(i))  = this%h_sq(i2(i)) + d2nn(i)
          else if (i2(i) == this%n+1) then
             this%h1(this%n)    = this%h1(this%n)   + d2(i)
             this%h(this%n)     = this%h(this%n)    + d2n(i)
             this%h_sq(this%n)  = this%h_sq(this%n) + d2nn(i)
          endif
       enddo

    endif

  endsubroutine histogram1d_add_linear_vals_norms


  !% Add multiple values to the histogram with linear interpolation
  !% and an additional mask
  subroutine histogram1d_add_linear_vals_norm_mask(this, vals, norm, mask)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: vals(:)
    real(DP), intent(in)              :: norm
    logical, intent(in)               :: mask(:)

    ! ---

    integer   :: i1(lbound(vals, 1):ubound(vals, 1)), i2(lbound(vals, 1):ubound(vals, 1))
    real(DP)  :: d1(lbound(vals, 1):ubound(vals, 1)), d2(lbound(vals, 1):ubound(vals, 1))

    real(DP)  :: d1n(lbound(vals, 1):ubound(vals, 1)), d2n(lbound(vals, 1):ubound(vals, 1))
    real(DP)  :: d1nn(lbound(vals, 1):ubound(vals, 1)), d2nn(lbound(vals, 1):ubound(vals, 1))

    integer   :: i

    ! ---

    i1  = int(floor((vals-this%min_b)/this%dbin))+1
    i2  = i1+1

    d1  = ( (i2-1) - (vals-this%min_b)/this%dbin )/this%dbin
    d2  = ( (vals-this%min_b)/this%dbin - (i1-1) )/this%dbin

    ASSERT(all(d1 >= 0.0_DP .and. d1*this%dbin <= 1.0_DP))
    ASSERT(all(d2 >= 0.0_DP .and. d2*this%dbin <= 1.0_DP))

    d1n   = d1*norm
    d2n   = d2*norm

    d1nn  = d1n*norm
    d2nn  = d2n*norm

    if (this%periodic) then

       i1  = modulo(i1-1, this%n)+1
       i2  = modulo(i2-1, this%n)+1

       do i = lbound(vals, 1), ubound(vals, 1)
          if (mask(i)) then
             this%h1(i1(i))    = this%h1(i1(i))  + d1(i)
             this%h1(i2(i))    = this%h1(i2(i))  + d2(i)

             this%h(i1(i))     = this%h(i1(i))   + d1n(i)
             this%h(i2(i))     = this%h(i2(i))   + d2n(i)

             this%h_sq(i1(i))  = this%h_sq(i1(i)) + d1nn(i)
             this%h_sq(i2(i))  = this%h_sq(i2(i)) + d2nn(i)
          endif
       enddo

    else

       do i = lbound(vals, 1), ubound(vals, 1)
          if (mask(i)) then
             if (i1(i) >= 1 .and. i1(i) <= this%n) then
                this%h1(i1(i))    = this%h1(i1(i))   + d1(i)
                this%h(i1(i))     = this%h(i1(i))    + d1n(i)
                this%h_sq(i1(i))  = this%h_sq(i1(i)) + d1nn(i)
             else if (i1(i) == 0) then
                this%h1(1)    = this%h1(1)   + d1(i)
                this%h(1)     = this%h(1)    + d1n(i)
                this%h_sq(1)  = this%h_sq(1) + d1nn(i)
             endif

             if (i2(i) >= 1 .and. i2(i) <= this%n) then
                this%h1(i2(i))    = this%h1(i2(i))   + d2(i)
                this%h(i2(i))     = this%h(i2(i))    + d2n(i)
                this%h_sq(i2(i))  = this%h_sq(i2(i)) + d2nn(i)
             else if (i2(i) == this%n+1) then
                this%h1(this%n)    = this%h1(this%n)   + d2(i)
                this%h(this%n)     = this%h(this%n)    + d2n(i)
                this%h_sq(this%n)  = this%h_sq(this%n) + d2nn(i)
             endif
          endif
       enddo

    endif

  endsubroutine histogram1d_add_linear_vals_norm_mask


  !% Add multiple values to the histogram with linear interpolation
  !% and an additional mask
  subroutine histogram1d_add_linear_vals_norms_mask(this, vals, norms, mask)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: vals(:)
    real(DP), intent(in)              :: norms(:)
    logical, intent(in)               :: mask(:)

    ! ---

    integer   :: i1(lbound(vals, 1):ubound(vals, 1)), i2(lbound(vals, 1):ubound(vals, 1))
    real(DP)  :: d1(lbound(vals, 1):ubound(vals, 1)), d2(lbound(vals, 1):ubound(vals, 1))

    real(DP)  :: d1n(lbound(vals, 1):ubound(vals, 1)), d2n(lbound(vals, 1):ubound(vals, 1))
    real(DP)  :: d1nn(lbound(vals, 1):ubound(vals, 1)), d2nn(lbound(vals, 1):ubound(vals, 1))

    integer   :: i

    ! ---

    i1  = int(floor((vals-this%min_b)/this%dbin))+1
    i2  = i1+1

    d1  = ( (i2-1) - (vals-this%min_b)/this%dbin )/this%dbin
    d2  = ( (vals-this%min_b)/this%dbin - (i1-1) )/this%dbin

    ASSERT(all(d1 >= 0.0_DP .and. d1*this%dbin <= 1.0_DP))
    ASSERT(all(d2 >= 0.0_DP .and. d2*this%dbin <= 1.0_DP))

    d1n   = d1*norms
    d2n   = d2*norms

    d1nn  = d1n*norms
    d2nn  = d2n*norms

    if (this%periodic) then

       i1  = modulo(i1-1, this%n)+1
       i2  = modulo(i2-1, this%n)+1

       do i = lbound(vals, 1), ubound(vals, 1)
          if (mask(i)) then
             this%h1(i1(i))    = this%h1(i1(i))  + d1(i)
             this%h1(i2(i))    = this%h1(i2(i))  + d2(i)

             this%h(i1(i))     = this%h(i1(i))   + d1n(i)
             this%h(i2(i))     = this%h(i2(i))   + d2n(i)

             this%h_sq(i1(i))  = this%h_sq(i1(i)) + d1nn(i)
             this%h_sq(i2(i))  = this%h_sq(i2(i)) + d2nn(i)
          endif
       enddo

    else

       do i = lbound(vals, 1), ubound(vals, 1)
          if (mask(i)) then
             if (i1(i) >= 1 .and. i1(i) <= this%n) then
                this%h1(i1(i))    = this%h1(i1(i))   + d1(i)
                this%h(i1(i))     = this%h(i1(i))    + d1n(i)
                this%h_sq(i1(i))  = this%h_sq(i1(i)) + d1nn(i)
             else if (i1(i) == 0) then
                this%h1(1)    = this%h1(1)   + d1(i)
                this%h(1)     = this%h(1)    + d1n(i)
                this%h_sq(1)  = this%h_sq(1) + d1nn(i)
             endif

             if (i2(i) >= 1 .and. i2(i) <= this%n) then
                this%h1(i2(i))    = this%h1(i2(i))   + d2(i)
                this%h(i2(i))     = this%h(i2(i))    + d2n(i)
                this%h_sq(i2(i))  = this%h_sq(i2(i)) + d2nn(i)
             else if (i2(i) == this%n+1) then
                this%h1(this%n)    = this%h1(this%n)   + d2(i)
                this%h(this%n)     = this%h(this%n)    + d2n(i)
                this%h_sq(this%n)  = this%h_sq(this%n) + d2nn(i)
             endif
          endif
       enddo

    endif

  endsubroutine histogram1d_add_linear_vals_norms_mask


  !% Add a value which is broadened by a smoothing function to the histogram
  subroutine histogram1d_add_smoothed(this, val, norm)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: val
    real(DP), intent(in)              :: norm

    ! ---

    real(DP)  :: expvals(this%n)

    integer   :: i1, i2, g1, g2
    real(DP)  :: x2

    ! ---

    x2  = (val-this%min_b)/this%dbin + 0.5_DP
    i1  = floor(x2)
    x2  = x2 - i1

    g1  = -this%ng
    g2  = this%ng-1

    expvals(1:2*this%ng)  = x2*this%smoothing_func(g1:g2) + (1.0_DP-x2)*this%smoothing_func(g1+1:g2+1)

    i1  = i1 - this%ng+1
    i2  = i1 + 2*this%ng - 1

    !
    ! Cut smoothing function if we're at the border
    !

    g1  = 1
    g2  = 2*this%ng

    if (i1 < 1) then
       g1 = g1 + (1 - i1)
       i1 = 1
    endif

    if (i2 > this%n) then
       g2 = g2 - (i2 - this%n)
       i2 = this%n
    endif

    this%h1(i1:i2)    = this%h1(i1:i2)   + expvals(g1:g2)
    this%h(i1:i2)     = this%h(i1:i2)    + expvals(g1:g2)*norm
    this%h_sq(i1:i2)  = this%h_sq(i1:i2) + expvals(g1:g2)*norm**2

  endsubroutine histogram1d_add_smoothed


  !**********************************************************************
  ! Add a value
  !**********************************************************************
  subroutine histogram1d_add(this, val, norm)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: val
    real(DP), intent(in), optional    :: norm

    ! ---

    real(DP)  :: n

    ! ---

    if (present(norm)) then
       n  = norm
    else
       n  = 1.0_DP
    endif

    if (this%interp == INTERP_LINEAR) then
       call histogram1d_add_linear(this, val, n)
    else
       call histogram1d_add_smoothed(this, val, n)
    endif

  endsubroutine histogram1d_add


  !% Add a list of values to the histogram
  subroutine histogram1d_add_vals(this, vals, norm)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: vals(:)
    real(DP), intent(in), optional    :: norm

    ! ---

    integer   :: i
    real(DP)  :: n

    ! ---

    if (present(norm)) then
       n  = norm
    else
       n  = 1.0_DP
    endif

    if (this%interp == INTERP_LINEAR) then

       call histogram1d_add_linear_vals(this, vals,  n)

    else
    
       do i = lbound(vals, 1), ubound(vals, 1)
          call histogram1d_add_smoothed(this, vals(i), n)
       enddo

    endif

  endsubroutine histogram1d_add_vals


  !% Add a list of values to the histogram
  subroutine histogram1d_add_vals_mask(this, vals, mask)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: vals(:)
    logical, intent(in)               :: mask(:)

    ! ---

    integer   :: i

    ! ---

    if (this%interp == INTERP_LINEAR) then

       do i = lbound(vals, 1), ubound(vals, 1)
          if (mask(i)) then
             call histogram1d_add_linear(this, vals(i), 1.0_DP)
          endif
       enddo

    else
    
       do i = lbound(vals, 1), ubound(vals, 1)
          if (mask(i)) then
             call histogram1d_add_smoothed(this, vals(i), 1.0_DP)
          endif
       enddo

    endif

  endsubroutine histogram1d_add_vals_mask


  !% Add a list of values to the histogram
  subroutine histogram1d_add_vals_norms(this, vals, norms)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: vals(:)
    real(DP), intent(in)              :: norms(:)

    ! ---

    integer   :: i

    ! ---

    if (this%interp == INTERP_LINEAR) then

       call histogram1d_add_linear_vals_norms(this, vals, norms)

    else
    
       do i = lbound(vals, 1), ubound(vals, 1)
          call histogram1d_add_smoothed(this, vals(i), norms(i))
       enddo

    endif

  endsubroutine histogram1d_add_vals_norms


  !% Add a list of values to the histogram
  subroutine histogram1d_add_vals_norm_mask(this, vals, norm, mask)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: vals(:)
    real(DP), intent(in)              :: norm
    logical, intent(in)               :: mask(:)

    ! ---

    integer   :: i

    ! ---

    if (this%interp == INTERP_LINEAR) then

       call histogram1d_add_linear_vals_norm_mask(this, vals, norm, mask)

    else
    
       do i = lbound(vals, 1), ubound(vals, 1)
          if (mask(i)) then
             call histogram1d_add_smoothed(this, vals(i), norm)
          endif
       enddo

    endif

  endsubroutine histogram1d_add_vals_norm_mask


  !% Add a list of values to the histogram
  subroutine histogram1d_add_vals_norms_mask(this, vals, norms, mask)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: vals(:)
    real(DP), intent(in)              :: norms(:)
    logical, intent(in)               :: mask(:)

    ! ---

    integer   :: i

    ! ---

    if (this%interp == INTERP_LINEAR) then

       call histogram1d_add_linear_vals_norms_mask(this, vals, norms, mask)

    else
    
       do i = lbound(vals, 1), ubound(vals, 1)
          if (mask(i)) then
             call histogram1d_add_smoothed(this, vals(i), norms(i))
          endif
       enddo

    endif

  endsubroutine histogram1d_add_vals_norms_mask


  !% Add a range of values to the histogram (linear interpolation only)
  subroutine histogram1d_add_range(this, vala, valb, norm, error)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: vala
    real(DP), intent(in)              :: valb
    real(DP), intent(in), optional    :: norm
    integer, intent(out), optional    :: error

    ! ---

    integer   :: i1a, i2a, i1b, i2b, i
    real(DP)  :: d1a, d1b, d2a, d2b
    real(DP)  :: va, vb, fac, n, n_sq

    ! ---

    INIT_ERROR(error)

    if (this%interp /= INTERP_LINEAR) then
       RAISE_ERROR("*histogram1d_add_range* can only be used with linear interpolation.", error)
    endif

    n  = 1.0_DP
    if (present(norm)) then
       n  = norm
    endif
    n_sq  = n**2

    if (vala > valb) then
       va  = valb
       vb  = vala
    else
       va  = vala
       vb  = valb
    endif

    i1a  = int(floor((va-this%min_b)/this%dbin))+1
    i2a  = i1a+1

    d1a  = ( (i2a-1) - (va-this%min_b)/this%dbin )
    d2a  = ( (va-this%min_b)/this%dbin - (i1a-1) )

    i1b  = int(floor((vb-this%min_b)/this%dbin))+1
    i2b  = i1b+1

    d1b  = ( (i2b-1) - (vb-this%min_b)/this%dbin )
    d2b  = ( (vb-this%min_b)/this%dbin - (i1b-1) )

    if (va == vb) then

       this%h1(i1a)    = this%h1(i1a)   + d1a
       this%h(i1a)     = this%h(i1a)    + d1a*n
       this%h_sq(i1a)  = this%h_sq(i1a) + d1a*n_sq

       this%h1(i2a)    = this%h1(i2a)   + d2a
       this%h(i2a)     = this%h(i2a)    + d2a*n
       this%h_sq(i2a)  = this%h_sq(i2a) + d2a*n_sq

    else

       fac  = 1.0_DP/(vb-va)

       if (i1a == i1b) then

          d1a  = fac * ( d2b * ( 1.0_DP - 0.5_DP*d2b ) - d2a * ( 1.0_DP - 0.5_DP*d2a ) )
          d2a  = fac * 0.5_DP * ( d2b**2 - d2a**2 )

          this%h1(i1a)    = this%h1(i1a)   + d1a
          this%h(i1a)     = this%h(i1a)    + d1a*n
          this%h_sq(i1a)  = this%h_sq(i1a) + d1a*n_sq

          this%h1(i2a)    = this%h1(i2a)   + d2a
          this%h(i2a)     = this%h(i2a)    + d2a*n
          this%h_sq(i2a)  = this%h_sq(i2a) + d2a*n_sq

       else

          forall(i = i2a:i1b-1)
             this%h1(i)      = this%h1(i)   + 0.5_DP*fac
             this%h(i)       = this%h(i)    + 0.5_DP*fac*n
             this%h_sq(i)    = this%h_sq(i) + 0.5_DP*fac*n_sq

             this%h1(i+1)    = this%h1(i+1)   + 0.5_DP*fac
             this%h(i+1)     = this%h(i+1)    + 0.5_DP*fac*n
             this%h_sq(i+1)  = this%h_sq(i+1) + 0.5_DP*fac*n_sq
          endforall
       
          d2a  = fac * d1a * ( 1 - 0.5_DP*d1a )
          d1a  = fac * 0.5_DP * d1a**2

          d1b  = fac * d2b * ( 1 - 0.5_DP*d2b )
          d2b  = fac * 0.5_DP * d2b**2

          this%h1(i1a)    = this%h1(i1a)   + d1a
          this%h(i1a)     = this%h(i1a)    + d1a*n
          this%h_sq(i1a)  = this%h_sq(i1a) + d1a*n_sq
       
          this%h1(i1b)    = this%h1(i1b)   + d1b
          this%h(i1b)     = this%h(i1b)    + d1b*n
          this%h_sq(i1b)  = this%h_sq(i1b) + d1b*n_sq
       
          this%h1(i2a)    = this%h1(i2a)   + d2a
          this%h(i2a)     = this%h(i2a)    + d2a*n
          this%h_sq(i2a)  = this%h_sq(i2a) + d2a*n_sq

          this%h1(i2b)    = this%h1(i2b)   + d2b
          this%h(i2b)     = this%h(i2b)    + d2b*n
          this%h_sq(i2b)  = this%h_sq(i2b) + d2b*n_sq

       endif

    endif

  endsubroutine histogram1d_add_range


  !% Add a range of values to the histogram (linear interpolation only)
  subroutine histogram1d_add_range_vals(this, valsa, valsb, norm, mask, error)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: valsa(:)
    real(DP), intent(in)              :: valsb(:)
    real(DP), intent(in), optional    :: norm
    logical, intent(in), optional     :: mask(:)
    integer, intent(out), optional    :: error

    ! ---

    integer   :: i, j
    integer   :: i1a(lbound(valsa, 1):ubound(valsa, 1)), i2a(lbound(valsa, 1):ubound(valsa, 1))
    integer   :: i1b(lbound(valsb, 1):ubound(valsb, 1)), i2b(lbound(valsb, 1):ubound(valsb, 1))
    real(DP)  :: d1a(lbound(valsa, 1):ubound(valsa, 1)), d2a(lbound(valsa, 1):ubound(valsa, 1))
    real(DP)  :: d1b(lbound(valsb, 1):ubound(valsb, 1)), d2b(lbound(valsb, 1):ubound(valsb, 1))
    real(DP)  :: va(lbound(valsa, 1):ubound(valsa, 1)), vb(lbound(valsb, 1):ubound(valsb, 1))
    real(DP)  :: fac(lbound(valsa, 1):ubound(valsa, 1))
    real(DP)  :: n, n_sq
    logical   :: m(lbound(valsa, 1):ubound(valsa, 1))

    ! ---

    INIT_ERROR(error)

    if (this%interp /= INTERP_LINEAR) then
       RAISE_ERROR("*histogram1d_add_range* can only be used with linear interpolation.", error)
    endif

    n  = 1.0_DP
    if (present(norm)) then
       n  = norm
    endif
    n_sq  = n**2

    where (valsa > valsb)
       va  = valsb
       vb  = valsa
    elsewhere
       va  = valsa
       vb  = valsb
    endwhere

    fac  = 1.0_DP
    where (va /= vb)
       fac  = 1.0_DP/(vb-va)
    endwhere

    m  = .true.
    if (present(mask)) then
       m  = mask
    endif

    i1a  = int(floor((va-this%min_b)/this%dbin))+1
    i2a  = i1a+1

    d1a  = ( (i2a-1) - (va-this%min_b)/this%dbin )
    d2a  = ( (va-this%min_b)/this%dbin - (i1a-1) )

    i1b  = int(floor((vb-this%min_b)/this%dbin))+1
    i2b  = i1b+1

    d1b  = ( (i2b-1) - (vb-this%min_b)/this%dbin )
    d2b  = ( (vb-this%min_b)/this%dbin - (i1b-1) )


    do j = lbound(valsa, 1), ubound(valsa, 1)

       if (m(j)) then

          if (va(j) == vb(j)) then

             this%h1(i1a(j))    = this%h1(i1a(j))   + d1a(j)
             this%h(i1a(j))     = this%h(i1a(j))    + d1a(j)*n
             this%h_sq(i1a(j))  = this%h_sq(i1a(j)) + d1a(j)*n_sq

             this%h1(i2a(j))    = this%h1(i2a(j))   + d2a(j)
             this%h(i2a(j))     = this%h(i2a(j))    + d2a(j)*n
             this%h_sq(i2a(j))  = this%h_sq(i2a(j)) + d2a(j)*n_sq

          else

             if (i1a(j) == i1b(j)) then
             
                d1a(j)  = fac(j) * ( d2b(j) * ( 1.0_DP - 0.5_DP*d2b(j) ) - d2a(j) * ( 1.0_DP - 0.5_DP*d2a(j) ) )
                d2a(j)  = fac(j) * 0.5_DP * ( d2b(j)**2 - d2a(j)**2 )

                this%h1(i1a(j))    = this%h1(i1a(j))   + d1a(j)
                this%h(i1a(j))     = this%h(i1a(j))    + d1a(j)*n
                this%h_sq(i1a(j))  = this%h_sq(i1a(j)) + d1a(j)*n_sq

                this%h1(i2a(j))    = this%h1(i2a(j))   + d2a(j)
                this%h(i2a(j))     = this%h(i2a(j))    + d2a(j)*n
                this%h_sq(i2a(j))  = this%h_sq(i2a(j)) + d2a(j)*n_sq

             else

                forall(i = i2a(j):i1b(j)-1)
                   this%h1(i)     = this%h1(i)     + 0.5_DP*fac(j)
                   this%h(i)      = this%h(i)      + 0.5_DP*fac(j)*n
                   this%h_sq(i)   = this%h_sq(i)   + 0.5_DP*fac(j)*n_sq

                   this%h1(i+1)   = this%h1(i+1)   + 0.5_DP*fac(j)
                   this%h(i+1)    = this%h(i+1)    + 0.5_DP*fac(j)*n
                   this%h_sq(i+1) = this%h_sq(i+1) + 0.5_DP*fac(j)*n_sq
                endforall

                d2a(j)  = fac(j) * d1a(j) * ( 1 - 0.5_DP*d1a(j) )
                d1a(j)  = fac(j) * 0.5_DP * d1a(j)**2

                d1b(j)  = fac(j) * d2b(j) * ( 1 - 0.5_DP*d2b(j) )
                d2b(j)  = fac(j) * 0.5_DP * d2b(j)**2

                this%h1(i1a(j))    = this%h1(i1a(j))   + d1a(j)
                this%h(i1a(j))     = this%h(i1a(j))    + d1a(j)*n
                this%h_sq(i1a(j))  = this%h_sq(i1a(j)) + d1a(j)*n_sq
       
                this%h1(i1b(j))    = this%h1(i1b(j))   + d1b(j)
                this%h(i1b(j))     = this%h(i1b(j))    + d1b(j)*n
                this%h_sq(i1b(j))  = this%h_sq(i1b(j)) + d1b(j)*n_sq
       
                this%h1(i2a(j))    = this%h1(i2a(j))   + d2a(j)
                this%h(i2a(j))     = this%h(i2a(j))    + d2a(j)*n
                this%h_sq(i2a(j))  = this%h_sq(i2a(j)) + d2a(j)*n_sq

                this%h1(i2b(j))    = this%h1(i2b(j))   + d2b(j)
                this%h(i2b(j))     = this%h(i2b(j))    + d2b(j)*n
                this%h_sq(i2b(j))  = this%h_sq(i2b(j)) + d2b(j)*n_sq

             endif

          endif

       endif

    enddo

  endsubroutine histogram1d_add_range_vals


  !% Add a range of values to the histogram (linear interpolation only)
  subroutine histogram1d_add_range_vals_norms(this, valsa, valsb, norms, mask, error)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: valsa(:)
    real(DP), intent(in)              :: valsb(:)
    real(DP), intent(in)              :: norms(:)
    logical, intent(in), optional     :: mask(:)
    integer, intent(out), optional    :: error

    ! ---

    integer   :: i, j
    integer   :: i1a(lbound(valsa, 1):ubound(valsa, 1)), i2a(lbound(valsa, 1):ubound(valsa, 1))
    integer   :: i1b(lbound(valsb, 1):ubound(valsb, 1)), i2b(lbound(valsb, 1):ubound(valsb, 1))
    real(DP)  :: d1a(lbound(valsa, 1):ubound(valsa, 1)), d2a(lbound(valsa, 1):ubound(valsa, 1))
    real(DP)  :: d1b(lbound(valsb, 1):ubound(valsb, 1)), d2b(lbound(valsb, 1):ubound(valsb, 1))
    real(DP)  :: va(lbound(valsa, 1):ubound(valsa, 1)), vb(lbound(valsb, 1):ubound(valsb, 1))
    real(DP)  :: fac(lbound(valsa, 1):ubound(valsa, 1))
    logical   :: m(lbound(valsa, 1):ubound(valsa, 1))

    ! ---

    INIT_ERROR(error)

    if (this%interp /= INTERP_LINEAR) then
       RAISE_ERROR("*histogram1d_add_range* can only be used with linear interpolation.", error)
    endif

    where (valsa > valsb)
       va  = valsb
       vb  = valsa
    elsewhere
       va  = valsa
       vb  = valsb
    endwhere

    fac  = 1.0_DP
    where (va /= vb)
       fac  = 1.0_DP/(vb-va)
    endwhere

    m  = .true.
    if (present(mask)) then
       m  = mask
    endif

    i1a  = int(floor((va-this%min_b)/this%dbin))+1
    i2a  = i1a+1

    d1a  = ( (i2a-1) - (va-this%min_b)/this%dbin )
    d2a  = ( (va-this%min_b)/this%dbin - (i1a-1) )

    i1b  = int(floor((vb-this%min_b)/this%dbin))+1
    i2b  = i1b+1

    d1b  = ( (i2b-1) - (vb-this%min_b)/this%dbin )
    d2b  = ( (vb-this%min_b)/this%dbin - (i1b-1) )


    do j = lbound(valsa, 1), ubound(valsa, 1)

       if (m(j)) then

          if (va(j) == vb(j)) then

             this%h1(i1a(j))    = this%h1(i1a(j))   + d1a(j)
             this%h(i1a(j))     = this%h(i1a(j))    + d1a(j)*norms(j)
             this%h_sq(i1a(j))  = this%h_sq(i1a(j)) + d1a(j)*norms(j)*norms(j)

             this%h1(i2a(j))    = this%h1(i2a(j))   + d2a(j)
             this%h(i2a(j))     = this%h(i2a(j))    + d2a(j)*norms(j)
             this%h_sq(i2a(j))  = this%h_sq(i2a(j)) + d2a(j)*norms(j)*norms(j)

          else

             if (i1a(j) == i1b(j)) then

                d1a(j)  = fac(j) * ( d2b(j) * ( 1.0_DP - 0.5_DP*d2b(j) ) - d2a(j) * ( 1.0_DP - 0.5_DP*d2a(j) ) )
                d2a(j)  = fac(j) * 0.5_DP * ( d2b(j)**2 - d2a(j)**2 )

                this%h1(i1a(j))    = this%h1(i1a(j))   + d1a(j)
                this%h(i1a(j))     = this%h(i1a(j))    + d1a(j)*norms(j)
                this%h_sq(i1a(j))  = this%h_sq(i1a(j)) + d1a(j)*norms(j)*norms(j)

                this%h1(i2a(j))    = this%h1(i2a(j))   + d2a(j)
                this%h(i2a(j))     = this%h(i2a(j))    + d2a(j)*norms(j)
                this%h_sq(i2a(j))  = this%h_sq(i2a(j)) + d2a(j)*norms(j)*norms(j)

             else

                forall(i = i2a(j):i1b(j)-1)
                   this%h1(i)      = this%h1(i)     + 0.5_DP*fac(j)
                   this%h(i)       = this%h(i)      + 0.5_DP*fac(j)*norms(j)
                   this%h_sq(i)    = this%h_sq(i)   + 0.5_DP*fac(j)*norms(j)*norms(j)

                   this%h1(i+1)    = this%h1(i+1)   + 0.5_DP*fac(j)
                   this%h(i+1)     = this%h(i+1)    + 0.5_DP*fac(j)*norms(j)
                   this%h_sq(i+1)  = this%h_sq(i+1) + 0.5_DP*fac(j)*norms(j)*norms(j)
                endforall

                d2a(j)  = fac(j) * d1a(j) * ( 1 - 0.5_DP*d1a(j) )
                d1a(j)  = fac(j) * 0.5_DP * d1a(j)**2

                d1b(j)  = fac(j) * d2b(j) * ( 1 - 0.5_DP*d2b(j) )
                d2b(j)  = fac(j) * 0.5_DP * d2b(j)**2

                this%h1(i1a(j))    = this%h1(i1a(j))   + d1a(j)
                this%h(i1a(j))     = this%h(i1a(j))    + d1a(j)*norms(j)
                this%h_sq(i1a(j))  = this%h_sq(i1a(j)) + d1a(j)*norms(j)*norms(j)
       
                this%h1(i1b(j))    = this%h1(i1b(j))   + d1b(j)
                this%h(i1b(j))     = this%h(i1b(j))    + d1b(j)*norms(j)
                this%h_sq(i1b(j))  = this%h_sq(i1b(j)) + d1b(j)*norms(j)*norms(j)
       
                this%h1(i2a(j))    = this%h1(i2a(j))   + d2a(j)
                this%h(i2a(j))     = this%h(i2a(j))    + d2a(j)*norms(j)
                this%h_sq(i2a(j))  = this%h_sq(i2a(j)) + d2a(j)*norms(j)*norms(j)

                this%h1(i2b(j))    = this%h1(i2b(j))   + d2b(j)
                this%h(i2b(j))     = this%h(i2b(j))    + d2b(j)*norms(j)
                this%h_sq(i2b(j))  = this%h_sq(i2b(j)) + d2b(j)*norms(j)*norms(j)

             endif

          endif

       endif

    enddo

  endsubroutine histogram1d_add_range_vals_norms


  !% Add two histograms
  subroutine histogram1d_add_histogram(this, that, fac)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    type(Histogram1D), intent(in)     :: that
    real(DP), intent(in), optional    :: fac

    ! ---
    
    this%h1  = this%h1 + that%h1
    if (present(fac)) then
       this%h   = this%h + fac*that%h
    else
       this%h   = this%h + that%h
    endif

  endsubroutine histogram1d_add_histogram


  !% Compute average value in each bin (does only makes sense if
  !% norm != 1 was used in the add functions)
  subroutine histogram1d_average(this, mpi)
    implicit none

    type(Histogram1D),           intent(inout)  :: this
    type(MPI_context), optional, intent(in)     :: mpi

    ! ---

    if (present(mpi)) then
       call sum_in_place(mpi, this%h1)
       call sum_in_place(mpi, this%h)
       call sum_in_place(mpi, this%h_sq)
    endif

    where (this%h1 == 0.0_DP)
       this%h1 = 1.0_DP
    endwhere

    this%h     = this%h / this%h1
    this%h_sq  = this%h_sq / this%h1

  endsubroutine histogram1d_average


  !% Compute Shannon entropy of this histogram. Needs to be called
  !% after normalize.
  elemental function histogram1d_entropy(this) result(val)
    implicit none

    type(Histogram1D), intent(in)  :: this
    real(DP)                       :: val

    ! ---

    integer   :: i
    real(DP)  :: S

    ! ---

    S  = 0.0_DP
    do i = 1, this%n
       if (this%h(i) > 0.0_DP) then
          S  = S - this%h(i)*log(this%h(i)*this%dbin)
       endif
    enddo

    val  = S*this%dbin

  endfunction histogram1d_entropy


  !% Compute expectation value
  elemental function histogram1d_expectation_value(this) result(val)
    implicit none

    type(Histogram1D), intent(in)  :: this
    real(DP)                       :: val

    ! ---

    real(DP)  :: n

    ! ---

    n = sum(this%h)

    if (n == 0.0_DP) then
       val = 0.0_DP
    else
       val  = sum(this%x * this%h / n)
    endif

  endfunction histogram1d_expectation_value


  !% Multiply the histogram by a value - for normalization
  subroutine histogram1d_mul(this, val)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: val

    ! ---

    this%h  = this%h * val

  endsubroutine histogram1d_mul


  !% Multiply multiple histograms by a value - for normalization
  subroutine histograms_mul(this, val)
    implicit none

    type(Histogram1D), intent(inout)  :: this(:)
    real(DP), intent(in)              :: val

    ! ---

    integer  :: i

    ! ---

    do i = lbound(this, 1), ubound(this, 1)
       this(i)%h  = this(i)%h * val
    enddo

  endsubroutine histograms_mul


  !% Multiply multiple histograms by a value - for normalization
  subroutine histograms2_mul(this, val)
    implicit none

    type(Histogram1D), intent(inout)  :: this(:, :)
    real(DP), intent(in)              :: val

    ! ---

    integer  :: i, j

    ! ---

    do i = lbound(this, 1), ubound(this, 1)
       do j = lbound(this, 2), ubound(this, 2)
          this(i, j)%h  = this(i, j)%h * val
       enddo
    enddo

  endsubroutine histograms2_mul


  !% Multiply the histogram by multiple value - for normalization
  subroutine histogram1d_mul_vals(this, vals)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: vals(this%n)

    ! ---

    this%h  = this%h * vals
    
  endsubroutine histogram1d_mul_vals


  !% Divide histogram by a value
  subroutine histogram1d_div(this, val)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: val

    ! ---

    this%h  = this%h / val

  endsubroutine histogram1d_div


  !% Divide multiple histograms by a value
  subroutine histograms_div(this, val)
    implicit none

    type(Histogram1D), intent(inout)  :: this(:)
    real(DP), intent(in)              :: val

    ! ---

    integer  :: i

    ! ---

    do i = lbound(this, 1), ubound(this, 1)
       this(i)%h  = this(i)%h / val
    enddo

  endsubroutine histograms_div


  !% Divide multiple histograms by a value
  subroutine histograms2_div(this, val)
    implicit none

    type(Histogram1D), intent(inout)  :: this(:, :)
    real(DP), intent(in)              :: val

    ! ---

    integer  :: i, j

    ! ---

    do i = lbound(this, 1), ubound(this, 1)
       do j = lbound(this, 2), ubound(this, 2)
          this(i, j)%h  = this(i, j)%h / val
       enddo
    enddo

  endsubroutine histograms2_div


  !% Divide histogram by multiple values
  subroutine histogram1d_div_vals(this, vals)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: vals(this%n)

    ! ---

    this%h  = this%h / vals
    
  endsubroutine histogram1d_div_vals


  !% Normalize histogram
  elemental subroutine histogram1d_normalize(this)
    implicit none

    type(Histogram1D), intent(inout)  :: this

    ! ---

    real(DP)  :: n

    ! ---

    n = sum(this%h*this%dbin)

    if (n /= 0.0_DP) then
       this%h  = this%h / n
    endif

  endsubroutine histogram1d_normalize


  !% OpenMP reduction - *thpriv* is a threadprivate histogram,
  !% this needs to be called within an *omp parallel* construct.
  subroutine histogram1d_reduce(this, thpriv)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    type(Histogram1D), intent(in)     :: thpriv

    ! ---

    !$omp single
    call initialise(this, thpriv)
    !$omp end single

    !$omp critical
    this%h     = this%h    + thpriv%h
    this%h_sq  = this%h_sq + thpriv%h_sq
    this%h1    = this%h1   + thpriv%h1
    !$omp end critical

  endsubroutine histogram1d_reduce


  !% Smooth histogram by convolution with a Gaussian
  subroutine histogram1d_smooth(this, sigma)
    implicit none

    type(Histogram1D), intent(inout)  :: this
    real(DP), intent(in)              :: sigma

    ! ---

    integer   :: i
    real(DP)  :: fac, sigma_sq, h(this%n)

    ! ---

    h         = 0.0_DP

    fac       = this%dbin/(sqrt(2*PI)*sigma)
    sigma_sq  = sigma**2

    do i = 1, this%n
       h  = h + fac*this%h(i)*exp(-((this%x - this%x(i))**2/(2*sigma_sq)))
    enddo

    this%h  = h

  endsubroutine histogram1d_smooth


  !% Sum histogram from different processors onto root
  subroutine histogram1d_sum_in_place(mpi, this)
    implicit none

    type(MPI_context), intent(in)     :: mpi
    type(Histogram1D), intent(inout)  :: this

    ! ---

    call sum_in_place(mpi, this%h)
    call sum_in_place(mpi, this%h_sq)
    call sum_in_place(mpi, this%h1)

  endsubroutine histogram1d_sum_in_place


#if 0
  !**********************************************************************
  ! Memory estimate logging
  !**********************************************************************
  subroutine log_memory_estimate_histogram(this)
    implicit none

    type(Histogram1D), intent(in)  :: this

    ! ---

    call log_memory_estimate(this%x)
    call log_memory_estimate(this%h)
    call log_memory_estimate(this%h_sq)
    call log_memory_estimate(this%h1)

  endsubroutine log_memory_estimate_histogram


  !**********************************************************************
  ! Memory estimate logging
  !**********************************************************************
  subroutine log_memory_estimate_histogram2(this)
    implicit none

    type(Histogram1D), intent(in)  :: this(:)

    ! ---

    integer  :: i

    ! ---

    do i = lbound(this, 1), ubound(this, 1)
       call log_memory_estimate_histogram(this(i))
    enddo

  endsubroutine log_memory_estimate_histogram2


  !**********************************************************************
  ! Memory estimate logging
  !**********************************************************************
  subroutine log_memory_estimate_histogram3(this)
    implicit none

    type(Histogram1D), intent(in)  :: this(:, :)

    ! ---

    integer  :: i, j

    ! ---

    do j = lbound(this, 2), ubound(this, 2)
       do i = lbound(this, 1), ubound(this, 1)
          call log_memory_estimate_histogram(this(i, j))
       enddo
    enddo

  endsubroutine log_memory_estimate_histogram3
#endif

endmodule histogram1d_module
