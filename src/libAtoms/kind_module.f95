module kind_module

  implicit none

  integer, parameter, public :: isp = selected_int_kind(9)
  integer, parameter, public :: idp = selected_int_kind(18)

#ifdef QUAD_PRECISION
  integer, parameter, public :: dp = 16
#elsif TEN_DIGIT_PRECISION
  integer, parameter, public :: dp = selected_real_kind(10)
#else
  integer, parameter, public :: dp = 8
#endif

#ifdef HAVE_QP
  integer, parameter, public :: qp = 16
#elsif TEN_DIGIT_PRECISION
  integer, parameter, public :: qp = selected_real_kind(10)
#else
  integer, parameter, public :: qp = 8
#endif

end module kind_module
