!**********************************************************************
! Histogram helper functions
!**********************************************************************

#include "macros.inc"

module histogram2d_module
  use libAtoms_module

  use io
  use logging

  private
  public : histogram2d, init, del, clear, set_bounds, write, add, average, mul, normalize, reduce, log_memory_estimate

  integer, parameter  :: INTERP_LINEAR  = 0
  integer, parameter  :: INTERP_LUCY    = 1

  type histogram2d_t

     !
     ! General stuff
     !

     integer                :: n(2) = -1    ! number of bins
     real(DP)               :: min_b(2)     ! minimum value
     real(DP)               :: max_b(2)     ! minimum value
     real(DP)               :: db(2)        ! difference
     real(DP)               :: dbin(2)      ! bin size

     logical                :: periodic(2)

     !
     ! Values
     !

     real(DP), allocatable  :: x1(:)      ! x-values
     real(DP), allocatable  :: x2(:)      ! x-values

     real(DP), allocatable  :: h(:, :)    ! values
     integer, allocatable   :: h1(:, :)   ! histogram2d with norm 1 (number of values which have been added to each bin)
   
  endtype histogram2d_t

  !
  ! Interface definition
  !

  interface init
     module procedure histogram2d_init, histogram2d_init_from_histogram2d
  endinterface

  interface del
     module procedure histogram2d_del
  endinterface

  interface clear
     module procedure histogram2d_clear
  endinterface

  interface set_bounds
     module procedure histogram2d_set_bounds
  endinterface

  interface write
     module procedure histogram2d_write, histogram2d_write_mult
     module procedure histogram2d_write_character_fn, histogram2d_write_mult_character_fn
  endinterface
,
  interface add
     module procedure histogram2d_add, histogram2d_add_vals, histogram2d_add_vals_norms
     module procedure histogram2d_add_vals_mask, histogram2d_add_vals_norm_mask, histogram2d_add_vals_norms_mask
     module procedure histogram2d_add_histogram2d
  endinterface

  interface average
     module procedure histogram2d_average
  endinterface

  interface mul
     module procedure histogram2d_mul, histogram2ds_mul, histogram2ds2_mul
  endinterface

  interface normalize
     module procedure histogram2d_normalize
  endinterface

  interface reduce
     module procedure histogram2d_reduce
  endinterface

  interface log_memory_estimate
     module procedure log_memory_estimate_histogram2d, log_memory_estimate_histogram2d2, log_memory_estimate_histogram2d3
  endinterface

contains

  !**********************************************************************
  ! Initialize the histogram2d
  !**********************************************************************
  elemental subroutine histogram2d_init(this, n1, n2, min_b1, max_b1, min_b2, max_b2, periodic1, periodic2)
    implicit none

    type(histogram2d_t), intent(out)  :: this
    integer, intent(in)               :: n1
    integer, intent(in)               :: n2
    real(DP), intent(in)              :: min_b1
    real(DP), intent(in)              :: max_b1
    real(DP), intent(in)              :: min_b2
    real(DP), intent(in)              :: max_b2
    logical, intent(in), optional     :: periodic1
    logical, intent(in), optional     :: periodic2

    ! ---

    call del(this)

    this%n  = (/ n1, n2 /)

    this%periodic  = .false.

    if (present(periodic1)) then
       this%periodic(1)  = periodic1
    endif

    if (present(periodic2)) then
       this%periodic(2)  = periodic2
    endif

    allocate(this%x1(this%n(1)))
    allocate(this%x2(this%n(2)))
    allocate(this%h(this%n(1), this%n(2)))
    allocate(this%h1(this%n(1), this%n(2)))

    this%min_b = (/ min_b1, min_b2 /) - 1.0_DP
    this%max_b = (/ max_b1, min_b2 /) - 1.0_DP

    call histogram2d_set_bounds(this, min_b1, max_b1, min_b2, max_b2)

    call histogram2d_clear(this)

  endsubroutine histogram2d_init


  !**********************************************************************
  ! Initialize the histogram2d
  !**********************************************************************
  elemental subroutine histogram2d_init_from_histogram2d(this, that)
    implicit none

    type(histogram2d_t), intent(out)  :: this
    type(histogram2d_t), intent(in)   :: that

    ! ---

    call del(this)

    this%n         = that%n
    this%periodic  = that%periodic

    allocate(this%x1(this%n(1)))
    allocate(this%x2(this%n(2)))
    allocate(this%h(this%n(1), this%n(2)))
    allocate(this%h1(this%n(1), this%n(2)))

    this%min_b = that%min_b - 1.0_DP
    this%max_b = that%max_b - 1.0_DP

    call histogram2d_set_bounds(this, that%min_b(1), that%max_b(1), that%min_b(2), that%max_b(2))

    call histogram2d_clear(this)

  endsubroutine histogram2d_init_from_histogram2d


  !**********************************************************************
  ! Clear histogram2d
  !**********************************************************************
  elemental subroutine histogram2d_clear(this)
    implicit none

    type(histogram2d_t), intent(inout)  :: this

    ! ---

    this%h1    = 0
    this%h     = 0.0_DP

  endsubroutine histogram2d_clear


  !**********************************************************************
  ! Initialize the histogram2d
  !**********************************************************************
  elemental subroutine histogram2d_set_bounds(this, min_b1, max_b1, min_b2, max_b2)
    implicit none

    type(histogram2d_t), intent(inout)  :: this
    real(DP), intent(in)                :: min_b1
    real(DP), intent(in)                :: max_b1
    real(DP), intent(in)                :: min_b2
    real(DP), intent(in)                :: max_b2

    ! ---

    integer   :: i

    ! ---

    if (any(this%min_b /= (/ min_b1, min_b2 /)) .or. any(this%max_b /= (/ max_b1, max_b2 /))) then

       this%min_b  = (/ min_b1, min_b2 /)
       this%max_b  = (/ max_b1, max_b2 /)
       this%db     = this%max_b - this%min_b

       this%dbin   = this%db/this%n

       do i = 1, this%n(1)
          this%x1(i)  = this%min_b(1) + this%dbin(1)*(i-0.5_DP)
          this%x2(i)  = this%min_b(2) + this%dbin(2)*(i-0.5_DP)
       enddo

    endif

  endsubroutine histogram2d_set_bounds


  !**********************************************************************
  ! Delete a histogram2d table
  !**********************************************************************
  elemental subroutine histogram2d_del(this)
    implicit none

    type(histogram2d_t), intent(inout)  :: this

    ! ---

    if (allocated(this%x1)) then
       deallocate(this%x1)
    endif
    if (allocated(this%x2)) then
       deallocate(this%x2)
    endif
    if (allocated(this%h)) then
       deallocate(this%h)
    endif
    if (allocated(this%h1)) then
       deallocate(this%h1)
    endif

  endsubroutine histogram2d_del


  !**********************************************************************
  ! Output the histogram2d table to a file (by file unit)
  !**********************************************************************
  subroutine histogram2d_write_mult(this, un, ierror)
    implicit none

    type(histogram2d_t), intent(in)   :: this(:)
    integer, intent(in)               :: un
    integer, intent(inout), optional  :: ierror

    ! ---

    integer        :: i, j, k
    character(80)  :: fmt

    integer        :: n(2)
    real(DP)       :: min_b(2), max_b(2), dbin(2), y(lbound(this, 1):ubound(this, 1))
    
    ! ---

    n      = this(lbound(this, 1))%n
    min_b  = this(lbound(this, 1))%min_b
    max_b  = this(lbound(this, 1))%max_b
    dbin   = this(lbound(this, 1))%dbin

    do i = lbound(this, 1)+1, ubound(this, 1)
       if (any(this(i)%n /= n)) then
          RAISE_ERROR("Number of histogram2d bins do not match.", ierror)
       endif

       if (any(this(i)%min_b /= min_b)) then
          RAISE_ERROR("*min_b*s do not match.", ierror)
       endif

       if (any(this(i)%max_b /= max_b)) then
          RAISE_ERROR("*max_b*s do not match.", ierror)
       endif
    enddo

    write (fmt, '(A,I4.4,A)')  "(", ubound(this, 1)-lbound(this, 1)+3, "ES20.10)"

    do i = 1, n(1)
       do j = 1, n(2)
          do k = lbound(this, 1), ubound(this, 1)
             y(k) = this(k)%h(i, j)
          enddo

          write (un, trim(fmt))  this(lbound(this, 1))%x1(i), this(lbound(this, 1))%x2(j), y
       enddo

       write (un, *)
    enddo

  endsubroutine histogram2d_write_mult


  !**********************************************************************
  ! Output the histogram2d table to a file (by file name)
  !**********************************************************************
  subroutine histogram2d_write_mult_character_fn(this, fn, ierror)
    implicit none

    type(histogram2d_t), intent(in)   :: this(:)
    character(*), intent(in)          :: fn
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: un

    ! ---

    un = fopen(fn, F_WRITE)
    call histogram2d_write_mult(this, un, ierror=ierror)
    call fclose(un)

    PASS_ERROR(ierror)

  endsubroutine histogram2d_write_mult_character_fn


  !**********************************************************************
  ! Output the histogram2d table to a file (by file unit)
  !**********************************************************************
  subroutine histogram2d_write(this, un, ierror)
    implicit none

    type(histogram2d_t), intent(in)   :: this
    integer, intent(in)               :: un
    integer, intent(inout), optional  :: ierror

    ! ---

    call histogram2d_write_mult((/ this /), un, ierror=ierror)
    PASS_ERROR(ierror)

  endsubroutine histogram2d_write


  !**********************************************************************
  ! Output the histogram2d table to a file (by file name)
  !**********************************************************************
  subroutine histogram2d_write_character_fn(this, fn, ierror)
    implicit none

    type(histogram2d_t), intent(in)   :: this
    character(*), intent(in)          :: fn
    integer, intent(inout), optional  :: ierror

    ! ---

    call histogram2d_write_mult_character_fn((/ this /), fn, ierror)
    PASS_ERROR(ierror)

  endsubroutine histogram2d_write_character_fn


  !**********************************************************************
  ! Add a value to the histogram2d with linear interpolation
  !**********************************************************************
  subroutine histogram2d_add_linear(this, val1, val2, norm)
    implicit none

    type(histogram2d_t), intent(inout)  :: this
    real(DP), intent(in)                :: val1
    real(DP), intent(in)                :: val2
    real(DP), intent(in)                :: norm

    ! ---

    integer   :: i1, i2

    ! ---

    i1  = int(floor((val1-this%min_b(1))/this%dbin(1)))+1
    i2  = int(floor((val2-this%min_b(2))/this%dbin(2)))+1

    if (this%periodic(1)) then
       i1  = modulo(i1-1, this%n(1))+1
    endif
    if (this%periodic(2)) then
       i2  = modulo(i2-1, this%n(2))+1
    endif

    if (i1 >= 1 .and. i1 <= this%n(1) .and. i2 >= 1 .and. i2 <= this%n(2)) then
       this%h1(i1, i2)    = this%h1(i1, i2) + 1
       this%h(i1, i2)     = this%h(i1, i2)  + norm
    endif

  endsubroutine histogram2d_add_linear


  !**********************************************************************
  ! Add a value to the histogram2d with linear interpolation
  !**********************************************************************
  subroutine histogram2d_add_linear_vals(this, vals1, vals2, norm)
    implicit none

    type(histogram2d_t), intent(inout)  :: this
    real(DP), intent(in)                :: vals1(:)
    real(DP), intent(in)                :: vals2(:)
    real(DP), intent(in)                :: norm

    ! ---

    integer   :: i1(lbound(vals1, 1):ubound(vals1, 1))
    integer   :: i2(lbound(vals2, 1):ubound(vals2, 1))

    integer   :: i

    ! ---

    i1  = int(floor((vals1-this%min_b(1))/this%dbin(1)))+1
    i2  = int(floor((vals2-this%min_b(2))/this%dbin(2)))+1

    if (this%periodic(1)) then
       i1  = modulo(i1-1, this%n(1))+1
    endif
    if (this%periodic(2)) then
       i2  = modulo(i2-1, this%n(2))+1
    endif
    
    do i = lbound(vals1, 1), ubound(vals1, 1)
       if (i1(i) >= 1 .and. i1(i) <= this%n(1) .and. i2(i) >= 1 .and. i2(i) <= this%n(2)) then
          this%h1(i1(i), i2(i))    = this%h1(i1(i), i2(i)) + 1
          this%h(i1(i), i2(i))     = this%h(i1(i), i2(i))  + norm
       endif
    enddo

  endsubroutine histogram2d_add_linear_vals


  !**********************************************************************
  ! Add a value
  !**********************************************************************
  subroutine histogram2d_add(this, val1, val2, norm)
    implicit none

    type(histogram2d_t), intent(inout)  :: this
    real(DP), intent(in)                :: val1
    real(DP), intent(in)                :: val2
    real(DP), intent(in), optional      :: norm

    ! ---

    real(DP)  :: n

    ! ---

    if (present(norm)) then
       n  = norm
    else
       n  = 1.0_DP
    endif

    call histogram2d_add_linear(this, val1, val2, n)

  endsubroutine histogram2d_add


  !**********************************************************************
  ! Add a list of values to the histogram2d
  !**********************************************************************
  subroutine histogram2d_add_vals(this, vals1, vals2, norm)
    implicit none

    type(histogram2d_t), intent(inout)  :: this
    real(DP), intent(in)                :: vals1(:)
    real(DP), intent(in)                :: vals2(:)
    real(DP), intent(in), optional      :: norm

    ! ---

    real(DP)  :: n

    ! ---

    if (present(norm)) then
       n  = norm
    else
       n  = 1.0_DP
    endif

    call histogram2d_add_linear_vals(this, vals1, vals2, n)

!    call histogram2d_add_linear_vals(this, vals,  n)

  endsubroutine histogram2d_add_vals


  !**********************************************************************
  ! Add a list of values to the histogram2d
  !**********************************************************************
  subroutine histogram2d_add_vals_mask(this, vals1, vals2, mask)
    implicit none

    type(histogram2d_t), intent(inout)  :: this
    real(DP), intent(in)                :: vals1(:)
    real(DP), intent(in)                :: vals2(:)
    logical, intent(in)                 :: mask(:)

    ! ---

    integer   :: i

    ! ---

    do i = lbound(vals1, 1), ubound(vals1, 1)
       if (mask(i)) then
          call histogram2d_add_linear(this, vals1(i), vals2(i), 1.0_DP)
       endif
    enddo

  endsubroutine histogram2d_add_vals_mask


  !**********************************************************************
  ! Add a list of values to the histogram2d
  !**********************************************************************
  subroutine histogram2d_add_vals_norms(this, vals1, vals2, norms)
    implicit none

    type(histogram2d_t), intent(inout)  :: this
    real(DP), intent(in)                :: vals1(:)
    real(DP), intent(in)                :: vals2(:)
    real(DP), intent(in)                :: norms(:)

    ! ---

    integer  :: i

    ! ---

    do i = lbound(vals1, 1), ubound(vals1, 1)
       call histogram2d_add_linear(this, vals1(i), vals2(i), norms(i))
    enddo

!    call histogram2d_add_linear_vals_norms(this, vals, norms)

  endsubroutine histogram2d_add_vals_norms


  !**********************************************************************
  ! Add a list of values to the histogram2d
  !**********************************************************************
  subroutine histogram2d_add_vals_norm_mask(this, vals1, vals2, norm, mask)
    implicit none

    type(histogram2d_t), intent(inout)  :: this
    real(DP), intent(in)                :: vals1(:)
    real(DP), intent(in)                :: vals2(:)
    real(DP), intent(in)                :: norm
    logical, intent(in)                 :: mask(:)

    ! ---

    integer  :: i

    ! ---

    do i = lbound(vals1, 1), ubound(vals1, 1)
       if (mask(i)) then
          call histogram2d_add_linear(this, vals1(i), vals2(i), norm)
       endif
    enddo

!    call histogram2d_add_linear_vals_norm_mask(this, vals, norm, mask)

  endsubroutine histogram2d_add_vals_norm_mask


  !**********************************************************************
  ! Add a list of values to the histogram2d
  !**********************************************************************
  subroutine histogram2d_add_vals_norms_mask(this, vals1, vals2, norms, mask)
    implicit none

    type(histogram2d_t), intent(inout)  :: this
    real(DP), intent(in)                :: vals1(:)
    real(DP), intent(in)                :: vals2(:)
    real(DP), intent(in)                :: norms(:)
    logical, intent(in)                 :: mask(:)

    ! ---

    integer  :: i

    ! ---

    do i = lbound(vals1, 1), ubound(vals1, 1)
       if (mask(i)) then
          call histogram2d_add_linear(this, vals1(i), vals2(i), norms(i))
       endif
    enddo

!    call histogram2d_add_linear_vals_norms_mask(this, vals, norms, mask)

  endsubroutine histogram2d_add_vals_norms_mask


  !**********************************************************************
  ! Add a list of values to the histogram2d
  !**********************************************************************
  subroutine histogram2d_add_histogram2d(this, that, fac)
    implicit none

    type(histogram2d_t), intent(inout)  :: this
    type(histogram2d_t), intent(in)     :: that
    real(DP), intent(in), optional      :: fac

    ! ---
    
    this%h1  = this%h1 + that%h1
    if (present(fac)) then
       this%h   = this%h + fac*that%h
    else
       this%h   = this%h + that%h
    endif

  endsubroutine histogram2d_add_histogram2d


  !**********************************************************************
  ! Compute average value in each bin (does only makes sense if
  ! norm != 1 was used in the add functions)
  !**********************************************************************
  elemental subroutine histogram2d_average(this)
    implicit none

    type(histogram2d_t), intent(inout)  :: this

    ! ---

    where (this%h1 == 0)
       this%h1 = 1
    endwhere

    this%h     = this%h / this%h1

  endsubroutine histogram2d_average


  !**********************************************************************
  ! Multiply the histogram2d by a value - for normalization
  !**********************************************************************
  subroutine histogram2d_mul(this, val)
    implicit none

    type(histogram2d_t), intent(inout)  :: this
    real(DP), intent(in)              :: val

    ! ---

    this%h  = this%h * val

  endsubroutine histogram2d_mul


  !**********************************************************************
  ! Multiply the histogram2d by a value - for normalization
  !**********************************************************************
  subroutine histogram2ds_mul(this, val)
    implicit none

    type(histogram2d_t), intent(inout)  :: this(:)
    real(DP), intent(in)                :: val

    ! ---

    integer  :: i

    ! ---

    do i = lbound(this, 1), ubound(this, 1)
       this(i)%h  = this(i)%h * val
    enddo

  endsubroutine histogram2ds_mul


  !**********************************************************************
  ! Multiply the histogram2d by a value - for normalization
  !**********************************************************************
  subroutine histogram2ds2_mul(this, val)
    implicit none

    type(histogram2d_t), intent(inout)  :: this(:, :)
    real(DP), intent(in)                :: val

    ! ---

    integer  :: i, j

    ! ---

    do i = lbound(this, 1), ubound(this, 1)
       do j = lbound(this, 2), ubound(this, 2)
          this(i, j)%h  = this(i, j)%h * val
       enddo
    enddo

  endsubroutine histogram2ds2_mul


  !**********************************************************************
  ! Normalize histogram2d
  !**********************************************************************
  elemental subroutine histogram2d_normalize(this)
    implicit none

    type(histogram2d_t), intent(inout)  :: this

    ! ---

    real(DP)  :: n

    ! ---

    n = sum(this%h*this%dbin(1)*this%dbin(2))

    if (n /= 0.0_DP) then
       this%h  = this%h / n
    endif

  endsubroutine histogram2d_normalize


  !**********************************************************************
  ! OpenMP reduction - *thpriv* is a threadprivate histogram,
  ! this needs to be called within an *omp parallel* construct.
  !**********************************************************************
  elemental subroutine histogram2d_reduce(this, thpriv)
    implicit none

    type(histogram2d_t), intent(inout)  :: this
    type(histogram2d_t), intent(in)     :: thpriv

    ! ---

    !$omp single
    call init(this, thpriv)
    !$omp end single

    !$omp critical
    this%h   = this%h  + thpriv%h
    this%h1  = this%h1 + thpriv%h1
    !$omp end critical

  endsubroutine histogram2d_reduce


  !**********************************************************************
  ! Memory estimate logging
  !**********************************************************************
  subroutine log_memory_estimate_histogram2d(this)
    implicit none

    type(histogram2d_t), intent(in)  :: this

    ! ---

    call log_memory_estimate(this%x1)
    call log_memory_estimate(this%x2)
    call log_memory_estimate(this%h)
    call log_memory_estimate(this%h1)

  endsubroutine log_memory_estimate_histogram2d


  !**********************************************************************
  ! Memory estimate logging
  !**********************************************************************
  subroutine log_memory_estimate_histogram2d2(this)
    implicit none

    type(histogram2d_t), intent(in)  :: this(:)

    ! ---

    integer  :: i

    ! ---

    do i = lbound(this, 1), ubound(this, 1)
       call log_memory_estimate_histogram2d(this(i))
    enddo

  endsubroutine log_memory_estimate_histogram2d2


  !**********************************************************************
  ! Memory estimate logging
  !**********************************************************************
  subroutine log_memory_estimate_histogram2d3(this)
    implicit none

    type(histogram2d_t), intent(in)  :: this(:, :)

    ! ---

    integer  :: i, j

    ! ---

    do j = lbound(this, 2), ubound(this, 2)
       do i = lbound(this, 1), ubound(this, 1)
          call log_memory_estimate_histogram2d(this(i, j))
       enddo
    enddo

  endsubroutine log_memory_estimate_histogram2d3

endmodule histogram2d_module
