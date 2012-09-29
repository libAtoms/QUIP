! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   libAtoms+QUIP: atomistic simulation library
! HND X
! HND X   Portions of this code were written by
! HND X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! HND X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! HND X
! HND X   Copyright 2006-2010.
! HND X
! HND X   Not for distribution
! HND X
! HND X   Portions of this code were written by Noam Bernstein as part of
! HND X   his employment for the U.S. Government, and are not subject
! HND X   to copyright in the USA.
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   http://www.libatoms.org
! HND X
! HND X  Additional contributions by
! HND X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module descriptors_module

   use libatoms_module

   implicit none

   private

   real(dp), dimension(:,:,:,:,:,:), allocatable, save :: cg_array
   integer :: cg_j1_max=0, cg_m1_max=0, cg_j2_max=0, cg_m2_max=0, cg_j_max=0, cg_m_max=0 
   logical :: cg_initialised = .false.

   real(dp), parameter :: factorial_table(0:16) = (/&
       1.0_dp, &
       1.0_dp, &
       2.0_dp, &
       6.0_dp, &
       24.0_dp, &
       120.0_dp, &
       720.0_dp, &
       5040.0_dp, &
       40320.0_dp, &
       362880.0_dp, &
       3628800.0_dp, &
       39916800.0_dp, &
       479001600.0_dp, &
       6227020800.0_dp, &
       87178291200.0_dp, &
       1307674368000.0_dp, &
       20922789888000.0_dp/)

   integer, parameter, public :: DT_NONE            =  0
   integer, parameter, public :: DT_BISPECTRUM_SO4  =  1
   integer, parameter, public :: DT_BISPECTRUM_SO3  =  2
   integer, parameter, public :: DT_BEHLER          =  3
   integer, parameter, public :: DT_DISTANCE_2B     =  4
   integer, parameter, public :: DT_COORDINATION    =  5
   integer, parameter, public :: DT_ANGLE_3B        =  6
   integer, parameter, public :: DT_CO_ANGLE_3B     =  7
   integer, parameter, public :: DT_CO_DISTANCE_2B  =  8
   integer, parameter, public :: DT_COSNX           =  9
   integer, parameter, public :: DT_TRIHIS          = 10
   integer, parameter, public :: DT_WATER_MONOMER   = 11
   integer, parameter, public :: DT_WATER_DIMER     = 12
   integer, parameter, public :: DT_A2_DIMER        = 13
   integer, parameter, public :: DT_AB_DIMER        = 14
   integer, parameter, public :: DT_BOND_REAL_SPACE = 15
   integer, parameter, public :: DT_ATOM_REAL_SPACE = 16
   integer, parameter, public :: DT_POWER_SO3       = 17
   integer, parameter, public :: DT_POWER_SO4       = 18
   integer, parameter, public :: DT_SOAP            = 19

   integer, parameter :: NP_WATER_DIMER    = 8
   integer, parameter :: NP_A2_DIMER       = 8
   integer, parameter :: NP_AB_DIMER       = 2

   type descriptor_data_mono
      real(dp), dimension(:), allocatable :: data
      real(dp), dimension(:,:,:), allocatable :: grad_data
      integer, dimension(:), allocatable :: ii
      real(dp), dimension(:,:), allocatable :: pos
      logical :: has_data
      logical, dimension(:), allocatable :: has_grad_data

      real(dp) :: covariance_cutoff = 1.0_dp
      real(dp), dimension(:,:), allocatable :: grad_covariance_cutoff
   endtype descriptor_data_mono

   type cplx_2d
      complex(dp), dimension(:,:), allocatable :: mm
   endtype cplx_2d

   type real_2d
      real(dp), dimension(:,:), allocatable :: mm
   endtype real_2d

   type cplx_3d
      complex(dp), dimension(:,:,:), allocatable :: mm
   endtype cplx_3d

   type RadialFunction_type
      integer :: n_max
      real(dp) :: cutoff, min_cutoff
      real(dp), dimension(:,:), allocatable :: RadialTransform
      real(dp), dimension(:), allocatable :: NormFunction

      logical :: initialised = .false.
   endtype RadialFunction_type

   type fourier_SO4_type
      real(dp) :: cutoff
      real(dp) :: z0_ratio
      real(dp) :: z0
      integer :: j_max, Z
      integer, dimension(:), allocatable :: species_Z
      real(dp), dimension(:), allocatable :: w

      logical :: initialised = .false.
   endtype fourier_SO4_type

   type bispectrum_SO4
      real(dp), pointer :: cutoff
      integer, pointer :: j_max, Z
      real(dp), pointer :: z0_ratio
      real(dp), pointer :: z0

      integer, dimension(:), pointer :: species_Z
      real(dp), dimension(:), pointer :: w
      
      type(fourier_SO4_type) :: fourier_SO4

      logical :: initialised = .false.

   endtype bispectrum_SO4

   type bispectrum_SO3

      integer :: l_max, n_max, Z
      real(dp) :: cutoff, min_cutoff

      type(RadialFunction_type) :: radial

      integer, dimension(:), allocatable :: species_Z
      real(dp), dimension(:), allocatable :: w

      logical :: initialised = .false.

   endtype bispectrum_SO3

   type behler_g2
      real(dp) :: eta
      real(dp) :: rs
      real(dp) :: rc
   endtype behler_g2

   type behler_g3
      real(dp) :: eta
      real(dp) :: lambda
      real(dp) :: zeta
      real(dp) :: rc
   endtype behler_g3

   type behler

      real(dp) :: cutoff = 0.0_dp
      logical :: initialised = .false.

      integer :: n_g2, n_g3
      type(behler_g2), dimension(:), allocatable :: g2
      type(behler_g3), dimension(:), allocatable :: g3

   endtype behler

   type distance_2b
      real(dp) :: cutoff
      integer :: Z1, Z2

      logical :: initialised = .false.

   endtype distance_2b

   type coordination
      real(dp) :: cutoff
      real(dp) :: transition_width
      integer :: Z

      logical :: initialised = .false.

   endtype coordination

   type angle_3b
      real(dp) :: cutoff
      integer :: Z, Z1, Z2

      logical :: initialised = .false.

   endtype angle_3b

   type co_angle_3b
      real(dp) :: cutoff
      real(dp) :: coordination_cutoff
      real(dp) :: coordination_transition_width
      integer :: Z, Z1, Z2

      logical :: initialised = .false.

   endtype co_angle_3b

   type co_distance_2b
      real(dp) :: cutoff
      real(dp) :: transition_width
      real(dp) :: coordination_cutoff
      real(dp) :: coordination_transition_width
      integer :: Z1, Z2

      logical :: initialised = .false.

   endtype co_distance_2b

   type cosnx

      integer :: l_max, n_max, Z
      real(dp) :: cutoff, min_cutoff

      type(RadialFunction_type) :: radial

      integer, dimension(:), allocatable :: species_Z
      real(dp), dimension(:), allocatable :: w

      logical :: initialised = .false.

   endtype cosnx

   type trihis
      real(dp) :: cutoff
      integer :: n_gauss

      real(dp), dimension(:,:), allocatable :: gauss_centre
      real(dp), dimension(:,:), allocatable :: gauss_width

      logical :: initialised = .false.

   endtype trihis

   type water_monomer
      real(dp) :: cutoff

      logical :: initialised = .false.

   endtype water_monomer

   type water_dimer
      real(dp) :: cutoff, cutoff_transition_width
      real(dp) :: monomer_cutoff
      logical :: OHH_ordercheck

      logical :: initialised = .false.

   endtype water_dimer

   type A2_dimer
      real(dp) :: cutoff
      real(dp) :: monomer_cutoff
      integer :: atomic_number

      logical :: initialised = .false.

   endtype A2_dimer

   type AB_dimer
      real(dp) :: cutoff
      real(dp) :: monomer_cutoff
      integer :: atomic_number1, atomic_number2

      logical :: initialised = .false.

   endtype AB_dimer

   type bond_real_space
      real(dp) :: bond_cutoff
      real(dp) :: bond_transition_width
      real(dp) :: cutoff
      real(dp) :: transition_width

      logical :: initialised = .false.

   endtype bond_real_space

   type atom_real_space
      real(dp) :: cutoff
      real(dp) :: cutoff_transition_width
      integer :: l_max
      real(dp) :: alpha
      real(dp) :: zeta

      logical :: initialised = .false.

   endtype atom_real_space

   type power_so3
      integer :: l_max, n_max, Z
      real(dp) :: cutoff, min_cutoff

      type(RadialFunction_type) :: radial

      integer, dimension(:), allocatable :: species_Z
      real(dp), dimension(:), allocatable :: w

      logical :: initialised = .false.
   endtype power_so3

   type power_SO4
      real(dp), pointer :: cutoff
      integer, pointer :: j_max, Z
      real(dp), pointer :: z0_ratio
      real(dp), pointer :: z0

      integer, dimension(:), pointer :: species_Z
      real(dp), dimension(:), pointer :: w
      
      type(fourier_SO4_type) :: fourier_SO4

      logical :: initialised = .false.

   endtype power_SO4

   type soap
      real(dp) :: cutoff
      real(dp) :: cutoff_transition_width
      integer :: l_max
      real(dp) :: alpha, sigma

      integer :: n_max
      real(dp), dimension(:), allocatable :: r_basis
      real(dp), dimension(:,:), allocatable :: transform_basis,cholesky_overlap_basis

      logical :: initialised = .false.
   endtype soap

   type descriptor
      integer :: descriptor_type = DT_NONE

      type(bispectrum_SO4)  :: descriptor_bispectrum_SO4
      type(bispectrum_SO3)  :: descriptor_bispectrum_SO3
      type(behler)          :: descriptor_behler
      type(distance_2b)     :: descriptor_distance_2b
      type(coordination)    :: descriptor_coordination
      type(angle_3b)        :: descriptor_angle_3b
      type(co_angle_3b)     :: descriptor_co_angle_3b
      type(co_distance_2b)  :: descriptor_co_distance_2b
      type(cosnx)           :: descriptor_cosnx
      type(trihis)          :: descriptor_trihis
      type(water_monomer)   :: descriptor_water_monomer
      type(water_dimer)     :: descriptor_water_dimer
      type(A2_dimer)        :: descriptor_A2_dimer
      type(AB_dimer)        :: descriptor_AB_dimer
      type(bond_real_space) :: descriptor_bond_real_space
      type(atom_real_space) :: descriptor_atom_real_space
      type(power_so3)       :: descriptor_power_so3
      type(power_SO4)       :: descriptor_power_SO4
      type(soap)            :: descriptor_soap

   endtype

   type descriptor_data
      type(descriptor_data_mono), dimension(:), allocatable :: x
      logical :: atom_based = .false.
   endtype descriptor_data

   type cplx_1d
      complex(dp), dimension(:), allocatable :: m
   endtype cplx_1d

   type real_1d
      real(dp), dimension(:), allocatable :: m
   endtype real_1d

   type spherical_harmonics_type
      type(cplx_1d), dimension(:), allocatable :: spherical_harmonics
      type(cplx_2d), dimension(:), allocatable :: grad_spherical_harmonics
      real(dp) :: r
      real(dp), dimension(3) :: u
   endtype spherical_harmonics_type

   type neighbour_type
      type(spherical_harmonics_type), dimension(:), allocatable :: neighbour
   endtype neighbour_type

   type grad_spherical_harmonics_overlap_type
      type(cplx_3d), dimension(:), allocatable :: grad_integral
   endtype grad_spherical_harmonics_overlap_type

   public :: neighbour_type, real_space_fourier_coefficients, real_space_covariance_coefficient, coordination_function, dcoordination_function
   public :: SphericalYCartesian

   interface initialise
      module procedure descriptor_initialise, RadialFunction_initialise, fourier_so4_initialise, &
      bispectrum_SO4_initialise, bispectrum_SO3_initialise, behler_initialise, distance_2b_initialise, &
      coordination_initialise, angle_3b_initialise, co_angle_3b_initialise, co_distance_2b_initialise, cosnx_initialise, trihis_initialise, &
      water_monomer_initialise, water_dimer_initialise, A2_dimer_initialise, AB_dimer_initialise, &
      bond_real_space_initialise, atom_real_space_initialise, power_so3_initialise, power_SO4_initialise, soap_initialise
   endinterface initialise
   public :: initialise

   interface finalise
      module procedure descriptor_finalise, descriptor_data_finalise, RadialFunction_finalise, fourier_so4_finalise, cplx_2d_array1_finalise, cplx_3d_array2_finalise, &
      bispectrum_SO4_finalise, bispectrum_SO3_finalise, behler_finalise, distance_2b_finalise, coordination_finalise, angle_3b_finalise, co_angle_3b_finalise, &
      co_distance_2b_finalise, cosnx_finalise, trihis_finalise, water_monomer_finalise, water_dimer_finalise, &
      A2_dimer_finalise, AB_dimer_finalise, bond_real_space_finalise, atom_real_space_finalise, power_so3_finalise, power_SO4_finalise, soap_finalise
   endinterface finalise
   public :: finalise

   interface calc
      module procedure descriptor_calc, bispectrum_SO4_calc, bispectrum_SO3_calc, behler_calc, distance_2b_calc, coordination_calc, angle_3b_calc, co_angle_3b_calc, &
      co_distance_2b_calc, cosnx_calc, trihis_calc, water_monomer_calc, water_dimer_calc, A2_dimer_calc, AB_dimer_calc, bond_real_space_calc, atom_real_space_calc, &
      power_so3_calc, power_SO4_calc, soap_calc
   endinterface calc
   public :: calc

   interface cutoff
      module procedure descriptor_cutoff, bispectrum_SO4_cutoff, bispectrum_SO3_cutoff, behler_cutoff, distance_2b_cutoff, coordination_cutoff, angle_3b_cutoff, co_angle_3b_cutoff, &
      co_distance_2b_cutoff, cosnx_cutoff, trihis_cutoff, water_monomer_cutoff, water_dimer_cutoff, A2_dimer_cutoff, AB_dimer_cutoff, bond_real_space_cutoff, atom_real_space_cutoff, &
      power_so3_cutoff, power_SO4_cutoff, soap_cutoff
   endinterface cutoff
   public :: cutoff

   interface descriptor_sizes
      module procedure descriptor_sizes, bispectrum_SO4_sizes, bispectrum_SO3_sizes, behler_sizes, distance_2b_sizes, coordination_sizes, angle_3b_sizes, co_angle_3b_sizes, &
      co_distance_2b_sizes, cosnx_sizes, trihis_sizes, water_monomer_sizes, water_dimer_sizes, A2_dimer_sizes, AB_dimer_sizes, bond_real_space_sizes, atom_real_space_sizes, &
      power_so3_sizes, power_SO4_sizes, soap_sizes
   endinterface descriptor_sizes
   public :: descriptor_sizes

   public :: descriptor, descriptor_data, descriptor_dimensions, descriptor_n_permutations, descriptor_permutations, descriptor_str_add_species
   public :: cg_initialise, real_space_covariance
   public :: cplx_1d, cplx_2d

   contains

   function get_descriptor_type(args_str,error)
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      integer :: get_descriptor_type

      type(Dictionary) :: params
      logical :: is_bispectrum_so4, is_bispectrum_so3, is_behler, is_distance_2b, is_coordination, is_angle_3b, &
         is_co_angle_3b, is_co_distance_2b, is_cosnx, is_trihis, is_water_monomer, is_water_dimer, is_A2_dimer, &
         is_AB_dimer, is_bond_real_space, is_atom_real_space, is_power_so3, is_power_so4, is_soap

      INIT_ERROR(error)

      call initialise(params)
      call param_register(params, 'bispectrum_so4', 'false', is_bispectrum_so4, help_string="Type of descriptor is bispectrum_so4.")
      call param_register(params, 'bispectrum_so3', 'false', is_bispectrum_so3, help_string="Type of descriptor is bispectrum_so3.")
      call param_register(params, 'behler', 'false', is_behler, help_string="Type of descriptor is behler.")
      call param_register(params, 'distance_2b', 'false', is_distance_2b, help_string="Type of descriptor is distance_2b.")
      call param_register(params, 'coordination', 'false', is_coordination, help_string="Type of descriptor is coordination.")
      call param_register(params, 'angle_3b', 'false', is_angle_3b, help_string="Type of descriptor is angle_3b.")
      call param_register(params, 'co_angle_3b', 'false', is_co_angle_3b, help_string="Type of descriptor is co_angle_3b.")
      call param_register(params, 'co_distance_2b', 'false', is_co_distance_2b, help_string="Type of descriptor is co_distance_2b.")
      call param_register(params, 'cosnx', 'false', is_cosnx, help_string="Type of descriptor is cosnx.")
      call param_register(params, 'trihis', 'false', is_trihis, help_string="Type of descriptor is trihis.")
      call param_register(params, 'water_monomer', 'false', is_water_monomer, help_string="Type of descriptor is water_monomer.")
      call param_register(params, 'water_dimer', 'false', is_water_dimer, help_string="Type of descriptor is water_dimer.")
      call param_register(params, 'A2_dimer', 'false', is_A2_dimer, help_string="Type of descriptor is A2_dimer.")
      call param_register(params, 'AB_dimer', 'false', is_AB_dimer, help_string="Type of descriptor is AB_dimer.")
      call param_register(params, 'bond_real_space', 'false', is_bond_real_space, help_string="Type of descriptor is bond_real_space.")
      call param_register(params, 'atom_real_space', 'false', is_atom_real_space, help_string="Type of descriptor is atom_real_space.")
      call param_register(params, 'power_so3', 'false', is_power_so3, help_string="Type of descriptor is power_so3.")
      call param_register(params, 'power_so4', 'false', is_power_so4, help_string="Type of descriptor is power_so4.")
      call param_register(params, 'soap', 'false', is_soap, help_string="Type of descriptor is soap.")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='descriptor_initialise args_str')) then
         RAISE_ERROR("descriptor_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      if (count( (/is_bispectrum_so4, is_bispectrum_so3, is_behler, is_distance_2b, is_coordination, is_angle_3b, is_co_angle_3b, is_co_distance_2b, &
      is_cosnx, is_trihis, is_water_monomer, is_water_dimer, is_A2_dimer, is_AB_dimer, is_bond_real_space, is_atom_real_space, is_power_so3, is_power_so4, is_soap/) ) /= 1) then
         RAISE_ERROR("descriptor_initialise found too few or too many IP Model types args_str='"//trim(args_str)//"'", error)
      endif

      get_descriptor_type = DT_NONE

      if( is_bispectrum_so4 ) then
         get_descriptor_type = DT_BISPECTRUM_SO4
      elseif( is_bispectrum_so3 ) then
         get_descriptor_type = DT_BISPECTRUM_SO3
      elseif( is_behler ) then
         get_descriptor_type = DT_BEHLER
      elseif( is_distance_2b ) then
         get_descriptor_type = DT_DISTANCE_2B
      elseif( is_coordination ) then
         get_descriptor_type = DT_COORDINATION
      elseif( is_angle_3b ) then
         get_descriptor_type = DT_ANGLE_3B
      elseif( is_co_angle_3b ) then
         get_descriptor_type = DT_CO_ANGLE_3B
      elseif( is_co_distance_2b ) then
         get_descriptor_type = DT_CO_DISTANCE_2B
      elseif( is_cosnx ) then
         get_descriptor_type = DT_COSNX
      elseif( is_trihis ) then
         get_descriptor_type = DT_TRIHIS
      elseif( is_water_monomer ) then
         get_descriptor_type = DT_WATER_MONOMER
      elseif( is_water_dimer ) then
         get_descriptor_type = DT_WATER_DIMER
      elseif( is_A2_dimer ) then
         get_descriptor_type = DT_A2_DIMER
      elseif( is_AB_dimer ) then
         get_descriptor_type = DT_AB_DIMER
      elseif( is_bond_real_space ) then
         get_descriptor_type = DT_BOND_REAL_SPACE
      elseif( is_atom_real_space ) then
         get_descriptor_type = DT_ATOM_REAL_SPACE
      elseif( is_power_so3 ) then
         get_descriptor_type = DT_POWER_SO3
      elseif( is_power_so4 ) then
         get_descriptor_type = DT_POWER_SO4
      elseif( is_soap ) then
         get_descriptor_type = DT_SOAP
      endif

   endfunction get_descriptor_type

   subroutine descriptor_initialise(this,args_str,error)
      type(descriptor), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      call finalise(this)

      this%descriptor_type = get_descriptor_type(args_str,error)

      select case(this%descriptor_type)
      case(DT_BISPECTRUM_SO4)
         call initialise(this%descriptor_bispectrum_SO4,args_str,error)
      case(DT_BISPECTRUM_SO3)
         call initialise(this%descriptor_bispectrum_SO3,args_str,error)
      case(DT_BEHLER)
         call initialise(this%descriptor_behler,args_str,error)
      case(DT_DISTANCE_2B)
         call initialise(this%descriptor_distance_2b,args_str,error)
      case(DT_COORDINATION)
         call initialise(this%descriptor_coordination,args_str,error)
      case(DT_ANGLE_3B)
         call initialise(this%descriptor_angle_3b,args_str,error)
      case(DT_CO_ANGLE_3B)
         call initialise(this%descriptor_co_angle_3b,args_str,error)
      case(DT_CO_DISTANCE_2B)
         call initialise(this%descriptor_co_distance_2b,args_str,error)
      case(DT_COSNX)
         call initialise(this%descriptor_cosnx,args_str,error)
      case(DT_TRIHIS)
         call initialise(this%descriptor_trihis,args_str,error)
      case(DT_WATER_MONOMER)
         call initialise(this%descriptor_water_monomer,args_str,error)
      case(DT_WATER_DIMER)
         call initialise(this%descriptor_water_dimer,args_str,error)
      case(DT_A2_DIMER)
         call initialise(this%descriptor_A2_dimer,args_str,error)
      case(DT_AB_DIMER)
         call initialise(this%descriptor_AB_dimer,args_str,error)
      case(DT_BOND_REAL_SPACE)
         call initialise(this%descriptor_bond_real_space,args_str,error)
      case(DT_ATOM_REAL_SPACE)
         call initialise(this%descriptor_atom_real_space,args_str,error)
      case(DT_POWER_SO3)
         call initialise(this%descriptor_power_so3,args_str,error)
      case(DT_POWER_SO4)
         call initialise(this%descriptor_power_so4,args_str,error)
      case(DT_SOAP)
         call initialise(this%descriptor_soap,args_str,error)
      endselect

   endsubroutine descriptor_initialise

   subroutine descriptor_finalise(this,error)
      type(descriptor), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            call finalise(this%descriptor_bispectrum_SO4,error)
         case(DT_BISPECTRUM_SO3)
            call finalise(this%descriptor_bispectrum_SO3,error)
         case(DT_BEHLER)
            call finalise(this%descriptor_behler,error)
         case(DT_DISTANCE_2b)
            call finalise(this%descriptor_distance_2b,error)
         case(DT_COORDINATION)
            call finalise(this%descriptor_coordination,error)
         case(DT_ANGLE_3B)
            call finalise(this%descriptor_angle_3b,error)
         case(DT_CO_ANGLE_3B)
            call finalise(this%descriptor_co_angle_3b,error)
         case(DT_CO_DISTANCE_2b)
            call finalise(this%descriptor_co_distance_2b,error)
         case(DT_COSNX)
            call finalise(this%descriptor_cosnx,error)
         case(DT_TRIHIS)
            call finalise(this%descriptor_trihis,error)
         case(DT_WATER_MONOMER)
            call finalise(this%descriptor_water_monomer,error)
         case(DT_WATER_DIMER)
            call finalise(this%descriptor_water_dimer,error)
         case(DT_A2_dimer)
            call finalise(this%descriptor_A2_dimer,error)
         case(DT_AB_dimer)
            call finalise(this%descriptor_AB_dimer,error)
         case(DT_BOND_REAL_SPACE)
            call finalise(this%descriptor_bond_real_space,error)
         case(DT_ATOM_REAL_SPACE)
            call finalise(this%descriptor_atom_real_space,error)
         case(DT_POWER_SO3)
            call finalise(this%descriptor_power_so3,error)
         case(DT_POWER_SO4)
            call finalise(this%descriptor_power_so4,error)
         case(DT_SOAP)
            call finalise(this%descriptor_soap,error)
      endselect

      this%descriptor_type = DT_NONE

   endsubroutine descriptor_finalise
   
   subroutine descriptor_data_finalise(this,error)
      type(descriptor_data), intent(inout) :: this
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(allocated(this%x)) then
         do i = 1, size(this%x)
            if(allocated(this%x(i)%data)) deallocate(this%x(i)%data)
            if(allocated(this%x(i)%grad_data)) deallocate(this%x(i)%grad_data)
            if(allocated(this%x(i)%ii)) deallocate(this%x(i)%ii)
            if(allocated(this%x(i)%pos)) deallocate(this%x(i)%pos)
            if(allocated(this%x(i)%has_grad_data)) deallocate(this%x(i)%has_grad_data)
            if(allocated(this%x(i)%grad_covariance_cutoff)) deallocate(this%x(i)%grad_covariance_cutoff)
         enddo
      endif

   endsubroutine descriptor_data_finalise

   subroutine RadialFunction_initialise(this,n_max,cutoff, min_cutoff,error)
      type(RadialFunction_type), intent(inout) :: this
      integer, intent(in) :: n_max
      real(dp), intent(in) :: cutoff, min_cutoff
      integer, optional, intent(out) :: error

      real(dp), dimension(:,:), allocatable :: S, vS
      real(dp), dimension(:), allocatable :: eS
      integer :: i, j

      INIT_ERROR(error)

      call finalise(this)

      this%n_max = n_max
      this%cutoff = cutoff
      this%min_cutoff = min_cutoff

      allocate(this%RadialTransform(this%n_max,this%n_max),this%NormFunction(this%n_max))
      allocate(S(this%n_max,this%n_max), vS(this%n_max,this%n_max), eS(this%n_max))

      do i = 1, this%n_max
         this%NormFunction(i) = sqrt(this%cutoff**(2.0_dp*i+5.0_dp)/(2.0_dp*i+5.0_dp))
         do j = 1, this%n_max
            S(j,i) = sqrt((2.0_dp*i+5)*(2.0_dp*j+5))/(i+j+5.0_dp)
         enddo
      enddo

      call diagonalise(S,eS,vS)
      this%RadialTransform = matmul(matmul(vS,diag(1.0_dp/sqrt(eS))),transpose(vS))

      if(allocated(S)) deallocate(S)
      if(allocated(vS)) deallocate(vS)
      if(allocated(eS)) deallocate(eS)

      this%initialised = .true.

   endsubroutine RadialFunction_initialise

   subroutine RadialFunction_finalise(this,error)
      type(RadialFunction_type), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%min_cutoff = 0.0_dp
      this%n_max = 0

      if(allocated(this%RadialTransform)) deallocate(this%RadialTransform)
      if(allocated(this%NormFunction)) deallocate(this%NormFunction)

      this%initialised = .false.

   endsubroutine RadialFunction_finalise

   subroutine cplx_2d_array1_finalise(this)
      type(cplx_2d), dimension(:), allocatable, intent(inout) :: this
      integer :: j

      if(allocated(this)) then
         do j = lbound(this,1), ubound(this,1)
            if(allocated(this(j)%mm)) deallocate(this(j)%mm)
         enddo
         deallocate(this)
      endif
   endsubroutine cplx_2d_array1_finalise

   subroutine cplx_3d_array2_finalise(this)
      type(cplx_3d), dimension(:,:), allocatable, intent(inout) :: this
      integer :: i, j

      if(allocated(this)) then
         do j = lbound(this,2), ubound(this,2)
            do i = lbound(this,1), ubound(this,1)
               if(allocated(this(i,j)%mm)) deallocate(this(i,j)%mm)
            enddo
         enddo
         deallocate(this)
      endif

   endsubroutine cplx_3d_array2_finalise

   subroutine fourier_SO4_calc(this,at,i,U,dU,error)
      type(fourier_SO4_type), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(in) :: i
      type(cplx_2d), dimension(:), allocatable, intent(inout) :: U
      type(cplx_3d), dimension(:,:), allocatable, intent(inout), optional :: dU
      integer, optional, intent(out) :: error

      complex(dp), dimension(:,:), allocatable :: Uc, Up
      complex(dp), dimension(:,:,:), allocatable :: dUc, dUp
      complex(dp) :: z0_pls_Iz, z0_min_Iz, x_pls_Iy, x_min_Iy
      complex(dp), dimension(3) :: dz0_pls_Iz, dz0_min_Iz, dx_pls_Iy, dx_min_Iy
      real(dp), dimension(3) :: diff, u_ij, dfcut, dz0, dr0
      real(dp) :: r0, r, fcut, z0, theta0
      integer :: n, n_i, ji, j, m1, m2
      integer, dimension(116) :: species_map

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR('fourier_SO4_calc: object not initialised',error)
      endif

      species_map = 0
      do j = 1, size(this%species_Z)
         if(this%species_Z(j) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(j)) = j
         endif
      enddo


      if(allocated(U)) then
         if(lbound(U,1) /= 0 .or. ubound(U,1) /= this%j_max) call finalise(U)
      endif

      if(.not.allocated(U)) then
         allocate( U(0:this%j_max) )
         do j = 0, this%j_max
            allocate( U(j)%mm(-j:j,-j:j) )
            U(j)%mm = CPLX_ZERO
         enddo
      endif

      do j = 0, this%j_max
         U(j)%mm = CPLX_ZERO
         do m1 = -j, j, 2
            U(j)%mm(m1,m1) = CPLX_ONE
         enddo
      enddo

      allocate( Uc(-this%j_max:this%j_max, -this%j_max:this%j_max), &
        Up(-this%j_max:this%j_max, -this%j_max:this%j_max) )

      Uc = CPLX_ZERO
      Up = CPLX_ZERO

      if(present(dU)) then
         if(allocated(dU)) call finalise(dU)

         ! dU is not allocated, allocate and zero it
         allocate( dU(0:this%j_max,0:atoms_n_neighbours(at,i,max_dist=this%cutoff)) )
         do j = 0, this%j_max
            allocate( dU(j,0)%mm(3,-j:j,-j:j) )
            dU(j,0)%mm = CPLX_ZERO
         enddo

         allocate( dUc(3,-this%j_max:this%j_max, -this%j_max:this%j_max), &
            dUp(3,-this%j_max:this%j_max, -this%j_max:this%j_max) )
         dUc = CPLX_ZERO
         dUp = CPLX_ZERO
      endif

      n_i = 0
      do n = 1, atoms_n_neighbours(at,i)
         ji = atoms_neighbour(at, i, n, distance=r, diff=diff, cosines=u_ij)
         if( r > this%cutoff ) cycle

         n_i = n_i + 1

         theta0 = r / this%z0
         z0 = r / tan( theta0 )
         r0 = sin( theta0 ) / r

         z0_pls_Iz = ( z0 + CPLX_IMAG*diff(3) ) * r0
         z0_min_Iz = ( z0 - CPLX_IMAG*diff(3) ) * r0
         x_pls_Iy = ( diff(1) + CPLX_IMAG*diff(2) ) * r0
         x_min_Iy = ( diff(1) - CPLX_IMAG*diff(2) ) * r0

         fcut = cutoff_function(r,this%cutoff) * this%w(species_map(at%Z(ji)))

         U(0)%mm(0,0) = U(0)%mm(0,0) + fcut
         Up(0:0,0:0) = CPLX_ONE

         if(present(dU)) then

            dfcut = -dcutoff_function(r,this%cutoff)*u_ij * this%w(species_map(at%Z(ji)))
            dz0 = ( 1.0_dp / tan( theta0 ) - theta0 / sin(theta0)**2 ) * u_ij
            dr0 = ( cos( theta0 ) / (r*this%z0) - r0 / r ) * u_ij

            dz0_pls_Iz = ( z0 + CPLX_IMAG*diff(3) )*dr0 + dz0*r0
            dz0_pls_Iz(3) = dz0_pls_Iz(3) + CPLX_IMAG*r0
               
            dz0_min_Iz = ( z0 - CPLX_IMAG*diff(3) )*dr0 + dz0*r0
            dz0_min_Iz(3) = dz0_min_Iz(3) - CPLX_IMAG*r0
               
            dx_pls_Iy = ( diff(1) + CPLX_IMAG*diff(2) )*dr0
            dx_pls_Iy(1) = dx_pls_Iy(1) + r0
            dx_pls_Iy(2) = dx_pls_Iy(2) + CPLX_IMAG*r0
               
            dx_min_Iy = ( diff(1) - CPLX_IMAG*diff(2) )*dr0
            dx_min_Iy(1) = dx_min_Iy(1) + r0
            dx_min_Iy(2) = dx_min_Iy(2) - CPLX_IMAG*r0

            dUc = CPLX_ZERO
            dUp = CPLX_ZERO

            dU(0,0)%mm(:,0,0) = dU(0,0)%mm(:,0,0) + dfcut*CPLX_ONE

            allocate( dU(0,n_i)%mm(3,-0:0,-0:0) )

            dU(0,n_i)%mm(:,0,0) = - dfcut*CPLX_ONE
         endif

         do j = 1, this%j_max
            Uc(-j:j,-j:j) = CPLX_ZERO
            if(present(dU)) then
               dUc(:,-j:j,-j:j) = CPLX_ZERO
               allocate( dU(j,n_i)%mm(3,-j:j,-j:j) )
               dU(j,n_i)%mm = CPLX_ZERO
            endif

            do m1 = -j, j-2, 2
               do m2 = -j, j, 2
                  if( (j-m2) /= 0 ) then
                     Uc(m2,m1) = Uc(m2,m1) + &
                     sqrt( real(j-m2,dp)/real(j-m1,dp) ) * z0_pls_Iz * Up(m2+1,m1+1)

                     if(present(dU)) dUc(:,m2,m1) = dUc(:,m2,m1) + &
                        sqrt( real(j-m2,dp)/real(j-m1,dp) ) * &
                        ( dz0_pls_Iz * Up(m2+1,m1+1) + z0_pls_Iz * dUp(:,m2+1,m1+1) )
                  endif

                  if( (j+m2) /= 0 ) then
                     Uc(m2,m1) = Uc(m2,m1) - &
                     CPLX_IMAG * sqrt( real(j+m2,dp)/real(j-m1,dp) ) * x_min_Iy * Up(m2-1,m1+1)

                     if(present(dU)) dUc(:,m2,m1) = dUc(:,m2,m1) - &
                        CPLX_IMAG * sqrt( real(j+m2,dp)/real(j-m1,dp) ) * &
                        ( dx_min_Iy * Up(m2-1,m1+1) + x_min_Iy * dUp(:,m2-1,m1+1) )

                  endif
               enddo
            enddo

            m1 = j
            do m2 = -j, j, 2
               if( (j+m2) /= 0 ) then
                  Uc(m2,m1) = Uc(m2,m1) + &
                  sqrt( real(j+m2,dp)/real(j+m1,dp) ) * z0_min_Iz * Up(m2-1,m1-1)

                  if(present(dU)) dUc(:,m2,m1) = dUc(:,m2,m1) + &
                     sqrt( real(j+m2,dp)/real(j+m1,dp) ) * &
                     ( dz0_min_Iz * Up(m2-1,m1-1) + z0_min_Iz * dUp(:,m2-1,m1-1) )
               endif

               if( (j-m2) /= 0 ) then
                  Uc(m2,m1) = Uc(m2,m1) - &
                  CPLX_IMAG * sqrt( real(j-m2,dp)/real(j+m1,dp) ) * x_pls_Iy * Up(m2+1,m1-1)

                  if(present(dU)) dUc(:,m2,m1) = dUc(:,m2,m1) - &
                     CPLX_IMAG * sqrt( real(j-m2,dp)/real(j+m1,dp) ) * &
                     ( dx_pls_Iy * Up(m2+1,m1-1) + x_pls_Iy * dUp(:,m2+1,m1-1) )
               endif
            enddo

            U(j)%mm = U(j)%mm + Uc(-j:j,-j:j) * fcut
            Up(-j:j,-j:j) = Uc(-j:j,-j:j)
            if(present(dU)) then
               dUp(:,-j:j,-j:j) = dUc(:,-j:j,-j:j)
               dU(j,0)%mm = dU(j,0)%mm - dUc(:,-j:j,-j:j) * fcut
               dU(j,n_i)%mm = dU(j,n_i)%mm + dUc(:,-j:j,-j:j) * fcut
               do m1 = -j, j, 2
                  do m2 = -j, j, 2
                     dU(j,0)%mm(:,m2,m1) = dU(j,0)%mm(:,m2,m1) &
                        + Uc(m2,m1) * dfcut
                     dU(j,n_i)%mm(:,m2,m1) = dU(j,n_i)%mm(:,m2,m1) &
                        - Uc(m2,m1) * dfcut
                  enddo
               enddo
            endif

         enddo ! j
      enddo ! n

      if(allocated(Up)) deallocate(Up)
      if(allocated(Uc)) deallocate(Uc)
      if(allocated(dUp)) deallocate(dUp)
      if(allocated(dUc)) deallocate(dUc)

   endsubroutine fourier_SO4_calc

   subroutine fourier_so4_initialise(this,args_str,error)
      type(fourier_SO4_type), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: n_species

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '2.75', this%cutoff, help_string="Cutoff for SO4 bispectrum")
      call param_register(params, 'z0_ratio', '0.0', this%z0_ratio, help_string="Ratio of radius of 4D projection sphere times PI and the cutoff.")
      call param_register(params, 'j_max', '4', this%j_max, help_string="Max of expansion of bispectrum, i.e. resulution")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'n_species', '1', n_species, help_string="Number of species for the descriptor")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='fourier_so4_initialise args_str')) then
         RAISE_ERROR("fourier_so4_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(n_species), this%w(n_species))

      call initialise(params)
      if( n_species == 1 ) then
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
         call param_register(params, 'w', '1.0', this%w(1), help_string="Weight associated to each atomic type")
      else
         call param_register(params, 'species_Z', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species")
         call param_register(params, 'w', PARAM_MANDATORY, this%w, help_string="Weight associated to each atomic type")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='fourier_so4_initialise args_str')) then
         RAISE_ERROR("fourier_so4_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%z0 = max(1.0_dp,this%z0_ratio) * this%cutoff/(PI-0.02_dp)
      
      this%initialised = .true.


   endsubroutine fourier_so4_initialise

   subroutine fourier_so4_finalise(this,error)
      type(fourier_so4_type), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%j_max = 0
      this%z0_ratio = 0.0_dp
      this%z0 = 0.0_dp
      this%Z = 0

      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%w)) deallocate(this%w)

      this%initialised = .false.

   endsubroutine fourier_so4_finalise

   subroutine bispectrum_so4_initialise(this,args_str,error)
      type(bispectrum_so4), intent(inout), target :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      call finalise(this)

      call initialise(this%fourier_SO4,args_str,error)

      this%cutoff => this%fourier_SO4%cutoff
      this%z0_ratio => this%fourier_SO4%z0_ratio
      this%z0 => this%fourier_SO4%z0
      this%j_max => this%fourier_SO4%j_max
      this%Z => this%fourier_SO4%Z
      this%cutoff => this%fourier_SO4%cutoff
      this%species_Z => this%fourier_SO4%species_Z
      this%w => this%fourier_SO4%w
      
      this%initialised = .true.

   endsubroutine bispectrum_so4_initialise

   subroutine bispectrum_so4_finalise(this,error)
      type(bispectrum_so4), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      call finalise(this%fourier_SO4,error)

      this%cutoff => null()
      this%z0_ratio => null()
      this%z0 => null()
      this%j_max => null()
      this%Z => null()
      this%cutoff => null()
      this%species_Z => null()
      this%w => null()

      this%initialised = .false.

   endsubroutine bispectrum_so4_finalise

   subroutine bispectrum_so3_initialise(this,args_str,error)
      type(bispectrum_so3), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: n_species

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)

      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for bispectrum_so3-type descriptors")
      call param_register(params, 'min_cutoff', '0.00', this%min_cutoff, help_string="Cutoff for minimal distances in bispectrum_so3-type descriptors")
      call param_register(params, 'l_max', '4', this%l_max, help_string="L_max for bispectrum_so3-type descriptors")
      call param_register(params, 'n_max', '4', this%n_max, help_string="N_max for bispectrum_so3-type descriptors")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'n_species', '1', n_species, help_string="Number of species for the descriptor")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='bispectrum_so3_initialise args_str')) then
         RAISE_ERROR("bispectrum_so3_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(n_species), this%w(n_species))

      call initialise(params)
      if( n_species == 1 ) then
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
         call param_register(params, 'w', '1.0', this%w(1), help_string="Weight associated to each atomic type")
      else
         call param_register(params, 'species_Z', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species")
         call param_register(params, 'w', PARAM_MANDATORY, this%w, help_string="Weight associated to each atomic type")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='bispectrum_so3_initialise args_str')) then
         RAISE_ERROR("bispectrum_so3_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call initialise(this%Radial,this%n_max,this%cutoff,this%min_cutoff,error)

      this%initialised = .true.

      call print('Dimensions: '//bispectrum_so3_dimensions(this,error))

   endsubroutine bispectrum_so3_initialise

   subroutine bispectrum_so3_finalise(this,error)
      type(bispectrum_so3), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%min_cutoff = 0.0_dp
      this%l_max = 0
      this%n_max = 0
      this%Z = 0

      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%w)) deallocate(this%w)

      call finalise(this%Radial)

      this%initialised = .false.

   endsubroutine bispectrum_so3_finalise

   subroutine behler_initialise(this,args_str,error)
      type(behler), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      !call param_register(params, 'behler_cutoff', '2.75', this%cutoff, help_string="Cutoff for Behler-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='behler_initialise args_str')) then
         RAISE_ERROR("behler_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif

      call finalise(params)

      this%n_g2 = 8
      this%n_g3 = 43

      allocate(this%g2(this%n_g2), this%g3(this%n_g3))

      this%g2(1)%eta = 0.001_dp / BOHR**2; this%g2(1)%rs = 0.000_dp * BOHR; this%g2(1)%rc = 11.338_dp * BOHR
      this%g2(2)%eta = 0.010_dp / BOHR**2; this%g2(2)%rs = 0.000_dp * BOHR; this%g2(2)%rc = 11.338_dp * BOHR
      this%g2(3)%eta = 0.020_dp / BOHR**2; this%g2(3)%rs = 0.000_dp * BOHR; this%g2(3)%rc = 11.338_dp * BOHR
      this%g2(4)%eta = 0.035_dp / BOHR**2; this%g2(4)%rs = 0.000_dp * BOHR; this%g2(4)%rc = 11.338_dp * BOHR
      this%g2(5)%eta = 0.060_dp / BOHR**2; this%g2(5)%rs = 0.000_dp * BOHR; this%g2(5)%rc = 11.338_dp * BOHR
      this%g2(6)%eta = 0.100_dp / BOHR**2; this%g2(6)%rs = 0.000_dp * BOHR; this%g2(6)%rc = 11.338_dp * BOHR
      this%g2(7)%eta = 0.200_dp / BOHR**2; this%g2(7)%rs = 0.000_dp * BOHR; this%g2(7)%rc = 11.338_dp * BOHR
      this%g2(8)%eta = 0.400_dp / BOHR**2; this%g2(8)%rs = 0.000_dp * BOHR; this%g2(8)%rc = 11.338_dp * BOHR

      this%g3( 1)%eta = 0.0001_dp / BOHR**2; this%g3( 1)%lambda = -1.000_dp; this%g3( 1)%zeta =  1.000_dp; this%g3( 1)%rc = 11.338_dp * BOHR
      this%g3( 2)%eta = 0.0001_dp / BOHR**2; this%g3( 2)%lambda =  1.000_dp; this%g3( 2)%zeta =  1.000_dp; this%g3( 2)%rc = 11.338_dp * BOHR
      this%g3( 3)%eta = 0.0001_dp / BOHR**2; this%g3( 3)%lambda = -1.000_dp; this%g3( 3)%zeta =  2.000_dp; this%g3( 3)%rc = 11.338_dp * BOHR
      this%g3( 4)%eta = 0.0001_dp / BOHR**2; this%g3( 4)%lambda =  1.000_dp; this%g3( 4)%zeta =  2.000_dp; this%g3( 4)%rc = 11.338_dp * BOHR
      this%g3( 5)%eta = 0.0030_dp / BOHR**2; this%g3( 5)%lambda = -1.000_dp; this%g3( 5)%zeta =  1.000_dp; this%g3( 5)%rc = 11.338_dp * BOHR
      this%g3( 6)%eta = 0.0030_dp / BOHR**2; this%g3( 6)%lambda =  1.000_dp; this%g3( 6)%zeta =  1.000_dp; this%g3( 6)%rc = 11.338_dp * BOHR
      this%g3( 7)%eta = 0.0030_dp / BOHR**2; this%g3( 7)%lambda = -1.000_dp; this%g3( 7)%zeta =  2.000_dp; this%g3( 7)%rc = 11.338_dp * BOHR
      this%g3( 8)%eta = 0.0030_dp / BOHR**2; this%g3( 8)%lambda =  1.000_dp; this%g3( 8)%zeta =  2.000_dp; this%g3( 8)%rc = 11.338_dp * BOHR
      this%g3( 9)%eta = 0.0080_dp / BOHR**2; this%g3( 9)%lambda = -1.000_dp; this%g3( 9)%zeta =  1.000_dp; this%g3( 9)%rc = 11.338_dp * BOHR
      this%g3(10)%eta = 0.0080_dp / BOHR**2; this%g3(10)%lambda =  1.000_dp; this%g3(10)%zeta =  1.000_dp; this%g3(10)%rc = 11.338_dp * BOHR
      this%g3(11)%eta = 0.0080_dp / BOHR**2; this%g3(11)%lambda = -1.000_dp; this%g3(11)%zeta =  2.000_dp; this%g3(11)%rc = 11.338_dp * BOHR
      this%g3(12)%eta = 0.0080_dp / BOHR**2; this%g3(12)%lambda =  1.000_dp; this%g3(12)%zeta =  2.000_dp; this%g3(12)%rc = 11.338_dp * BOHR
      this%g3(13)%eta = 0.0150_dp / BOHR**2; this%g3(13)%lambda = -1.000_dp; this%g3(13)%zeta =  1.000_dp; this%g3(13)%rc = 11.338_dp * BOHR
      this%g3(14)%eta = 0.0150_dp / BOHR**2; this%g3(14)%lambda =  1.000_dp; this%g3(14)%zeta =  1.000_dp; this%g3(14)%rc = 11.338_dp * BOHR
      this%g3(15)%eta = 0.0150_dp / BOHR**2; this%g3(15)%lambda = -1.000_dp; this%g3(15)%zeta =  2.000_dp; this%g3(15)%rc = 11.338_dp * BOHR
      this%g3(16)%eta = 0.0150_dp / BOHR**2; this%g3(16)%lambda =  1.000_dp; this%g3(16)%zeta =  2.000_dp; this%g3(16)%rc = 11.338_dp * BOHR
      this%g3(17)%eta = 0.0150_dp / BOHR**2; this%g3(17)%lambda = -1.000_dp; this%g3(17)%zeta =  4.000_dp; this%g3(17)%rc = 11.338_dp * BOHR
      this%g3(18)%eta = 0.0150_dp / BOHR**2; this%g3(18)%lambda =  1.000_dp; this%g3(18)%zeta =  4.000_dp; this%g3(18)%rc = 11.338_dp * BOHR
      this%g3(19)%eta = 0.0150_dp / BOHR**2; this%g3(19)%lambda = -1.000_dp; this%g3(19)%zeta = 16.000_dp; this%g3(19)%rc = 11.338_dp * BOHR
      this%g3(20)%eta = 0.0150_dp / BOHR**2; this%g3(20)%lambda =  1.000_dp; this%g3(20)%zeta = 16.000_dp; this%g3(20)%rc = 11.338_dp * BOHR
      this%g3(21)%eta = 0.0250_dp / BOHR**2; this%g3(21)%lambda = -1.000_dp; this%g3(21)%zeta =  1.000_dp; this%g3(21)%rc = 11.338_dp * BOHR
      this%g3(22)%eta = 0.0250_dp / BOHR**2; this%g3(22)%lambda =  1.000_dp; this%g3(22)%zeta =  1.000_dp; this%g3(22)%rc = 11.338_dp * BOHR
      this%g3(23)%eta = 0.0250_dp / BOHR**2; this%g3(23)%lambda = -1.000_dp; this%g3(23)%zeta =  2.000_dp; this%g3(23)%rc = 11.338_dp * BOHR
      this%g3(24)%eta = 0.0250_dp / BOHR**2; this%g3(24)%lambda =  1.000_dp; this%g3(24)%zeta =  2.000_dp; this%g3(24)%rc = 11.338_dp * BOHR
      this%g3(25)%eta = 0.0250_dp / BOHR**2; this%g3(25)%lambda = -1.000_dp; this%g3(25)%zeta =  4.000_dp; this%g3(25)%rc = 11.338_dp * BOHR
      this%g3(26)%eta = 0.0250_dp / BOHR**2; this%g3(26)%lambda =  1.000_dp; this%g3(26)%zeta =  4.000_dp; this%g3(26)%rc = 11.338_dp * BOHR
      this%g3(27)%eta = 0.0250_dp / BOHR**2; this%g3(27)%lambda = -1.000_dp; this%g3(27)%zeta = 16.000_dp; this%g3(27)%rc = 11.338_dp * BOHR
      this%g3(28)%eta = 0.0250_dp / BOHR**2; this%g3(28)%lambda =  1.000_dp; this%g3(28)%zeta = 16.000_dp; this%g3(28)%rc = 11.338_dp * BOHR
      this%g3(29)%eta = 0.0450_dp / BOHR**2; this%g3(29)%lambda = -1.000_dp; this%g3(29)%zeta =  1.000_dp; this%g3(29)%rc = 11.338_dp * BOHR
      this%g3(30)%eta = 0.0450_dp / BOHR**2; this%g3(30)%lambda =  1.000_dp; this%g3(30)%zeta =  1.000_dp; this%g3(30)%rc = 11.338_dp * BOHR
      this%g3(31)%eta = 0.0450_dp / BOHR**2; this%g3(31)%lambda = -1.000_dp; this%g3(31)%zeta =  2.000_dp; this%g3(31)%rc = 11.338_dp * BOHR
      this%g3(32)%eta = 0.0450_dp / BOHR**2; this%g3(32)%lambda =  1.000_dp; this%g3(32)%zeta =  2.000_dp; this%g3(32)%rc = 11.338_dp * BOHR
      this%g3(33)%eta = 0.0450_dp / BOHR**2; this%g3(33)%lambda = -1.000_dp; this%g3(33)%zeta =  4.000_dp; this%g3(33)%rc = 11.338_dp * BOHR
      this%g3(34)%eta = 0.0450_dp / BOHR**2; this%g3(34)%lambda =  1.000_dp; this%g3(34)%zeta =  4.000_dp; this%g3(34)%rc = 11.338_dp * BOHR
      this%g3(35)%eta = 0.0450_dp / BOHR**2; this%g3(35)%lambda = -1.000_dp; this%g3(35)%zeta = 16.000_dp; this%g3(35)%rc = 11.338_dp * BOHR
      this%g3(36)%eta = 0.0450_dp / BOHR**2; this%g3(36)%lambda =  1.000_dp; this%g3(36)%zeta = 16.000_dp; this%g3(36)%rc = 11.338_dp * BOHR
      this%g3(37)%eta = 0.0800_dp / BOHR**2; this%g3(37)%lambda = -1.000_dp; this%g3(37)%zeta =  1.000_dp; this%g3(37)%rc = 11.338_dp * BOHR
      this%g3(38)%eta = 0.0800_dp / BOHR**2; this%g3(38)%lambda =  1.000_dp; this%g3(38)%zeta =  1.000_dp; this%g3(38)%rc = 11.338_dp * BOHR
      this%g3(39)%eta = 0.0800_dp / BOHR**2; this%g3(39)%lambda = -1.000_dp; this%g3(39)%zeta =  2.000_dp; this%g3(39)%rc = 11.338_dp * BOHR
      this%g3(40)%eta = 0.0800_dp / BOHR**2; this%g3(40)%lambda =  1.000_dp; this%g3(40)%zeta =  2.000_dp; this%g3(40)%rc = 11.338_dp * BOHR
      this%g3(41)%eta = 0.0800_dp / BOHR**2; this%g3(41)%lambda = -1.000_dp; this%g3(41)%zeta =  4.000_dp; this%g3(41)%rc = 11.338_dp * BOHR
      this%g3(42)%eta = 0.0800_dp / BOHR**2; this%g3(42)%lambda =  1.000_dp; this%g3(42)%zeta =  4.000_dp; this%g3(42)%rc = 11.338_dp * BOHR
      this%g3(43)%eta = 0.0800_dp / BOHR**2; this%g3(43)%lambda =  1.000_dp; this%g3(43)%zeta = 16.000_dp; this%g3(43)%rc = 11.338_dp * BOHR

      this%cutoff = 11.338_dp * BOHR

      this%initialised = .true.

   endsubroutine behler_initialise

   subroutine behler_finalise(this,error)
      type(behler), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      this%cutoff = 0.0_dp
      this%n_g2 = 0
      this%n_g3 = 0

      if(allocated(this%g2)) deallocate(this%g2)
      if(allocated(this%g3)) deallocate(this%g3)

      this%initialised = .false.

   endsubroutine behler_finalise

   subroutine distance_2b_initialise(this,args_str,error)
      type(distance_2b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for distance_2b-type descriptors")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atom type #1 in bond")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atom type #2 in bond")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='distance_2b_initialise args_str')) then
         RAISE_ERROR("distance_2b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine distance_2b_initialise

   subroutine distance_2b_finalise(this,error)
      type(distance_2b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%Z1 = 0
      this%Z2 = 0

      this%initialised = .false.

   endsubroutine distance_2b_finalise

   subroutine coordination_initialise(this,args_str,error)
      type(coordination), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for coordination-type descriptors")
      call param_register(params, 'transition_width', '0.20', this%transition_width, help_string="Width of transition region from 1 to 0")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='coordination_initialise args_str')) then
         RAISE_ERROR("coordination_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine coordination_initialise

   subroutine coordination_finalise(this,error)
      type(coordination), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%transition_width = 0.0_dp
      this%Z = 0

      this%initialised = .false.

   endsubroutine coordination_finalise

   subroutine angle_3b_initialise(this,args_str,error)
      type(angle_3b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for angle_3b-type descriptors")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atomic number of neighbour #1")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atomic number of neighbour #2")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='angle_3b_initialise args_str')) then
         RAISE_ERROR("angle_3b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine angle_3b_initialise

   subroutine angle_3b_finalise(this,error)
      type(angle_3b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%Z = 0
      this%Z1 = 0
      this%Z2 = 0

      this%initialised = .false.

   endsubroutine angle_3b_finalise

   subroutine co_angle_3b_initialise(this,args_str,error)
      type(co_angle_3b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for co_angle_3b-type descriptors")
      call param_register(params, 'coordination_cutoff', '0.00', this%coordination_cutoff, help_string="Cutoff for coordination function in co_angle_3b-type descriptors")
      call param_register(params, 'coordination_transition_width', '0.00', this%coordination_transition_width, help_string="Transition width for co_angle_3b-type descriptors")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atomic number of neighbour #1")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atomic number of neighbour #2")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='co_angle_3b_initialise args_str')) then
         RAISE_ERROR("co_angle_3b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine co_angle_3b_initialise

   subroutine co_angle_3b_finalise(this,error)
      type(co_angle_3b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%coordination_cutoff = 0.0_dp
      this%coordination_transition_width = 0.0_dp
      this%Z = 0
      this%Z1 = 0
      this%Z2 = 0

      this%initialised = .false.

   endsubroutine co_angle_3b_finalise

   subroutine co_distance_2b_initialise(this,args_str,error)
      type(co_distance_2b), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for co_distance_2b-type descriptors")
      call param_register(params, 'transition_width', '0.50', this%transition_width, help_string="Transition width of cutoff for co_distance_2b-type descriptors")
      call param_register(params, 'coordination_cutoff', '0.00', this%coordination_cutoff, help_string="Cutoff for coordination function in co_distance_2b-type descriptors")
      call param_register(params, 'coordination_transition_width', '0.00', this%coordination_transition_width, help_string="Transition width for co_distance_2b-type descriptors")
      call param_register(params, 'Z1', '0', this%Z1, help_string="Atom type #1 in bond")
      call param_register(params, 'Z2', '0', this%Z2, help_string="Atom type #2 in bond")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='co_distance_2b_initialise args_str')) then
         RAISE_ERROR("co_distance_2b_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine co_distance_2b_initialise

   subroutine co_distance_2b_finalise(this,error)
      type(co_distance_2b), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%coordination_cutoff = 0.0_dp
      this%coordination_transition_width = 0.0_dp
      this%Z1 = 0
      this%Z2 = 0

      this%initialised = .false.

   endsubroutine co_distance_2b_finalise

   subroutine cosnx_initialise(this,args_str,error)
      type(cosnx), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: n_species

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for cosnx-type descriptors")
      call param_register(params, 'min_cutoff', '0.00', this%min_cutoff, help_string="Cutoff for minimal distances in cosnx-type descriptors")
      call param_register(params, 'l_max', '4', this%l_max, help_string="L_max for cosnx-type descriptors")
      call param_register(params, 'n_max', '4', this%n_max, help_string="N_max for cosnx-type descriptors")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'n_species', '1', n_species, help_string="Number of species for the descriptor")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='cosnx_initialise args_str')) then
         RAISE_ERROR("cosnx_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(n_species), this%w(n_species))

      call initialise(params)
      if( n_species == 1 ) then
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
         call param_register(params, 'w', '1.0', this%w(1), help_string="Weight associated to each atomic type")
      else
         call param_register(params, 'species_Z', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species")
         call param_register(params, 'w', PARAM_MANDATORY, this%w, help_string="Weight associated to each atomic type")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='cosnx_initialise args_str')) then
         RAISE_ERROR("cosnx_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call initialise(this%Radial,this%n_max,this%cutoff,this%min_cutoff,error)

      this%initialised = .true.

   endsubroutine cosnx_initialise

   subroutine cosnx_finalise(this,error)
      type(cosnx), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%min_cutoff = 0.0_dp
      this%l_max = 0
      this%n_max = 0

      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%w)) deallocate(this%w)

      call finalise(this%Radial)

      this%initialised = .false.

   endsubroutine cosnx_finalise

   subroutine trihis_initialise(this,args_str,error)
      type(trihis), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      real(dp), dimension(:), allocatable :: gauss_centre1D, gauss_width1D

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for trihis-type descriptors")
      call param_register(params, 'n_gauss', '0', this%n_gauss, help_string="Number of Gaussians for trihis-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='trihis_initialise args_str')) then
         RAISE_ERROR("trihis_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(gauss_centre1D(3*this%n_gauss),gauss_width1D(3*this%n_gauss))
      allocate(this%gauss_centre(3,this%n_gauss),this%gauss_width(3,this%n_gauss))

      call initialise(params)
      call param_register(params, 'trihis_gauss_centre', PARAM_MANDATORY, gauss_centre1D, help_string="Number of Gaussians for trihis-type descriptors")
      call param_register(params, 'trihis_gauss_width', PARAM_MANDATORY, gauss_width1D, help_string="Number of Gaussians for trihis-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='trihis_initialise args_str')) then
         RAISE_ERROR("trihis_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%gauss_centre = reshape(gauss_centre1D,(/3,this%n_gauss/))
      this%gauss_width = reshape(gauss_width1D,(/3,this%n_gauss/))

      deallocate(gauss_centre1D,gauss_width1D)

      this%initialised = .true.

   endsubroutine trihis_initialise

   subroutine trihis_finalise(this,error)
      type(trihis), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%n_gauss = 0

      if(allocated(this%gauss_centre)) deallocate(this%gauss_centre)
      if(allocated(this%gauss_width)) deallocate(this%gauss_width)

      this%initialised = .false.

   endsubroutine trihis_finalise

   subroutine water_monomer_initialise(this,args_str,error)
      type(water_monomer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for water_monomer-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='water_monomer_initialise args_str')) then
         RAISE_ERROR("water_monomer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine water_monomer_initialise

   subroutine water_monomer_finalise(this,error)
      type(water_monomer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp

      this%initialised = .false.

   endsubroutine water_monomer_finalise

   subroutine water_dimer_initialise(this,args_str,error)
      type(water_dimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for water_dimer-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, help_string="Width of smooth cutoff region for water_dimer-type descriptors")
      call param_register(params, 'monomer_cutoff', '1.50', this%monomer_cutoff, help_string="Monomer cutoff for water_dimer-type descriptors")
      call param_register(params, 'OHH_ordercheck', 'T', this%OHH_ordercheck, help_string="T: find water molecules. F: use default order OHH")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='water_dimer_initialise args_str')) then
         RAISE_ERROR("water_dimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine water_dimer_initialise

   subroutine water_dimer_finalise(this,error)
      type(water_dimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%monomer_cutoff = 0.0_dp
      this%OHH_ordercheck = .true.

      this%initialised = .false.

   endsubroutine water_dimer_finalise

   subroutine A2_dimer_initialise(this,args_str,error)
      type(A2_dimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for A2_dimer-type descriptors")
      call param_register(params, 'monomer_cutoff', '1.50', this%monomer_cutoff, help_string="Monomer cutoff for A2_dimer-type descriptors")
      call param_register(params, 'atomic_number', '1', this%atomic_number, help_string="Atomic number in A2_dimer-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='A2_dimer_initialise args_str')) then
         RAISE_ERROR("A2_dimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine A2_dimer_initialise

   subroutine A2_dimer_finalise(this,error)
      type(A2_dimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%monomer_cutoff = 0.0_dp
      this%atomic_number = 0

      this%initialised = .false.

   endsubroutine A2_dimer_finalise

   subroutine AB_dimer_initialise(this,args_str,error)
      type(AB_dimer), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for AB_dimer-type descriptors")
      call param_register(params, 'monomer_cutoff', '1.50', this%monomer_cutoff, help_string="Monomer cutoff for AB_dimer-type descriptors")
      call param_register(params, 'atomic_number1', '1', this%atomic_number1, help_string="Atomic number of atom 1 in AB_dimer-type descriptors")
      call param_register(params, 'atomic_number2', '9', this%atomic_number2, help_string="Atomic number of atom 2 in AB_dimer-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='AB_dimer_initialise args_str')) then
         RAISE_ERROR("AB_dimer_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      if( this%atomic_number1 == this%atomic_number2 ) then
         RAISE_ERROR("AB_dimer_initialise: AB_dimer_atomic_number1 = AB_dimer_atomic_number2 = "//this%atomic_number1//" which would require addtional permutational symmetries. Use A2_dimer descriptor instead.",error)
      endif

      this%initialised = .true.

   endsubroutine AB_dimer_initialise

   subroutine AB_dimer_finalise(this,error)
      type(AB_dimer), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%monomer_cutoff = 0.0_dp
      this%atomic_number1 = 0
      this%atomic_number2 = 0

      this%initialised = .false.

   endsubroutine AB_dimer_finalise

   subroutine bond_real_space_initialise(this,args_str,error)
      type(bond_real_space), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'bond_cutoff', '0.00', this%bond_cutoff, help_string="Bond cutoff for bond_real_space-type descriptors")
      call param_register(params, 'bond_transition_width', '0.00', this%bond_transition_width, help_string="Bond transition width for bond_real_space-type descriptors")
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Space cutoff for bond_real_space-type descriptors")
      call param_register(params, 'transition_width', '0.00', this%transition_width, help_string="Space transition width for bond_real_space-type descriptors")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='bond_real_space_initialise args_str')) then
         RAISE_ERROR("bond_real_space_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine bond_real_space_initialise

   subroutine bond_real_space_finalise(this,error)
      type(bond_real_space), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%bond_cutoff = 0.0_dp
      this%bond_transition_width = 0.0_dp
      this%cutoff = 0.0_dp
      this%transition_width = 0.0_dp

      this%initialised = .false.

   endsubroutine bond_real_space_finalise

   subroutine atom_real_space_initialise(this,args_str,error)
      type(atom_real_space), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Space cutoff for atom_real_space-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.00', this%cutoff_transition_width, help_string="Space transition width for atom_real_space-type descriptors")
      call param_register(params, 'l_max', '0', this%l_max, help_string="Cutoff for spherical harmonics expansion")
      call param_register(params, 'alpha', '1.0', this%alpha, help_string="Width of atomic Gaussians")
      call param_register(params, 'zeta', '1.0', this%zeta, help_string="Exponent of covariance function")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='atom_real_space_initialise args_str')) then
         RAISE_ERROR("atom_real_space_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%initialised = .true.

   endsubroutine atom_real_space_initialise

   subroutine atom_real_space_finalise(this,error)
      type(atom_real_space), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%l_max = 0
      this%alpha = 0.0_dp
      this%zeta = 0.0_dp

      this%initialised = .false.

   endsubroutine atom_real_space_finalise

   subroutine power_so3_initialise(this,args_str,error)
      type(power_so3), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: n_species

      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '0.00', this%cutoff, help_string="Cutoff for power_so3-type descriptors")
      call param_register(params, 'min_cutoff', '0.00', this%min_cutoff, help_string="Cutoff for minimal distances in power_so3-type descriptors")
      call param_register(params, 'l_max', '4', this%l_max, help_string="L_max for power_so3-type descriptors")
      call param_register(params, 'n_max', '4', this%n_max, help_string="N_max for power_so3-type descriptors")
      call param_register(params, 'Z', '0', this%Z, help_string="Atomic number of central atom")
      call param_register(params, 'n_species', '1', n_species, help_string="Number of species for the descriptor")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='power_so3_initialise args_str')) then
         RAISE_ERROR("power_so3_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      allocate(this%species_Z(n_species), this%w(n_species))

      call initialise(params)
      if( n_species == 1 ) then
         call param_register(params, 'species_Z', '0', this%species_Z(1), help_string="Atomic number of species")
         call param_register(params, 'w', '1.0', this%w(1), help_string="Weight associated to each atomic type")
      else
         call param_register(params, 'species_Z', PARAM_MANDATORY, this%species_Z, help_string="Atomic number of species")
         call param_register(params, 'w', PARAM_MANDATORY, this%w, help_string="Weight associated to each atomic type")
      endif

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='power_so3_initialise args_str')) then
         RAISE_ERROR("power_so3_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      call initialise(this%Radial,this%n_max,this%cutoff,this%min_cutoff,error)

      this%initialised = .true.

   endsubroutine power_so3_initialise

   subroutine power_so3_finalise(this,error)
      type(power_so3), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%min_cutoff = 0.0_dp
      this%l_max = 0
      this%n_max = 0
      this%Z = 0

      if(allocated(this%species_Z)) deallocate(this%species_Z)
      if(allocated(this%w)) deallocate(this%w)

      call finalise(this%Radial)

      this%initialised = .false.

   endsubroutine power_so3_finalise

   subroutine power_so4_initialise(this,args_str,error)
      type(power_so4), intent(inout), target :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      call finalise(this)

      call initialise(this%fourier_SO4,args_str,error)

      this%cutoff => this%fourier_SO4%cutoff
      this%z0_ratio => this%fourier_SO4%z0_ratio
      this%z0 => this%fourier_SO4%z0
      this%j_max => this%fourier_SO4%j_max
      this%Z => this%fourier_SO4%Z
      this%cutoff => this%fourier_SO4%cutoff
      this%species_Z => this%fourier_SO4%species_Z
      this%w => this%fourier_SO4%w
      
      this%initialised = .true.

   endsubroutine power_so4_initialise

   subroutine power_so4_finalise(this,error)
      type(power_so4), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      call finalise(this%fourier_SO4,error)

      this%cutoff => null()
      this%z0_ratio => null()
      this%z0 => null()
      this%j_max => null()
      this%Z => null()
      this%cutoff => null()
      this%species_Z => null()
      this%w => null()

      this%initialised = .false.

   endsubroutine power_so4_finalise

   subroutine soap_initialise(this,args_str,error)
      type(soap), intent(inout) :: this
      character(len=*), intent(in) :: args_str
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      real(dp) :: alpha_basis, spacing_basis, cutoff_basis, basis_error_exponent
      real(dp), dimension(:,:), allocatable :: covariance_basis, overlap_basis, cholesky_overlap_basis
      integer :: i, j

      type(LA_Matrix) :: LA_covariance_basis, LA_overlap_basis


      INIT_ERROR(error)

      call finalise(this)

      call initialise(params)
      call param_register(params, 'cutoff', '6.00', this%cutoff, help_string="Cutoff for soap-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.50', this%cutoff_transition_width, help_string="Cutoff transition width for soap-type descriptors")
      call param_register(params, 'l_max', '8', this%l_max, help_string="L_max (spherical harmonics basis band limit) for soap-type descriptors")
      call param_register(params, 'n_max', '40', this%n_max, help_string="N_max (number of radial basis functions) for soap-type descriptors")
      call param_register(params, 'sigma', '0.50', this%sigma, help_string="Width of atomic Gaussians for soap-type descriptors")
      call param_register(params, 'basis_error_exponent', '10.0', basis_error_exponent, help_string="10^(-basis_error_exponent) is the max difference between the target and the expanded function")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='soap_initialise args_str')) then
         RAISE_ERROR("soap_initialise failed to parse args_str='"//trim(args_str)//"'", error)
      endif
      call finalise(params)

      this%alpha = 0.5_dp / this%sigma**2

      alpha_basis = this%alpha
      cutoff_basis = this%cutoff + this%sigma * sqrt(2.0_dp * basis_error_exponent * log(10.0_dp))
      spacing_basis = cutoff_basis / this%n_max

      allocate(this%r_basis(this%n_max), this%transform_basis(this%n_max,this%n_max), &
         covariance_basis(this%n_max,this%n_max), overlap_basis(this%n_max,this%n_max), this%cholesky_overlap_basis(this%n_max,this%n_max))

      !this%r_basis(this%n_max) = cutoff_basis
      !do i = this%n_max-1, 1, -1
      !   this%r_basis(i)  = this%r_basis(i+1) - spacing_basis
      !enddo

      this%r_basis(1) = 0.0_dp
      do i = 2, this%n_max
         this%r_basis(i)  = this%r_basis(i-1) + spacing_basis
      enddo

      do i = 1, this%n_max
         do j = 1, this%n_max
            covariance_basis(j,i) = exp(-alpha_basis * (this%r_basis(i) - this%r_basis(j))**2)
            !overlap_basis(j,i) = exp(-0.5_dp * alpha_basis* (this%r_basis(i) - this%r_basis(j))**2) * ( 1.0_dp + erf( sqrt(alpha_basis/2.0_dp) * (this%r_basis(i) + this%r_basis(j)) ) )
            !print*, 'A', exp( -alpha_basis*(this%r_basis(i)**2+this%r_basis(j)**2) )
            !print*, 'B', sqrt(2.0_dp) * alpha_basis**1.5_dp * (this%r_basis(i) + this%r_basis(j))
            !print*, 'C', alpha_basis*exp(0.5_dp * alpha_basis * (this%r_basis(i) + this%r_basis(j))**2)*sqrt(PI)*(1.0_dp + alpha_basis*(this%r_basis(i) + this%r_basis(j))**2 )
            !print*, 'D', ( 1.0_dp + erf( sqrt(alpha_basis/2.0_dp) * (this%r_basis(i) + this%r_basis(j)) ) )
            !overlap_basis(j,i) = exp( -alpha_basis*(this%r_basis(i)**2+this%r_basis(j)**2) ) * &
            !   ( sqrt(2.0_dp) * alpha_basis**1.5_dp * (this%r_basis(i) + this%r_basis(j)) + &
            !   alpha_basis*exp(0.5_dp * alpha_basis * (this%r_basis(i) + this%r_basis(j))**2)*sqrt(PI)*(1.0_dp + alpha_basis*(this%r_basis(i) + this%r_basis(j))**2 ) * &
            !   ( 1.0_dp + erf( sqrt(alpha_basis/2.0_dp) * (this%r_basis(i) + this%r_basis(j)) ) ) )

            overlap_basis(j,i) = ( exp( -alpha_basis*(this%r_basis(i)**2+this%r_basis(j)**2) ) * &
               sqrt(2.0_dp) * alpha_basis**1.5_dp * (this%r_basis(i) + this%r_basis(j)) + &
               alpha_basis*exp(-0.5_dp * alpha_basis * (this%r_basis(i) - this%r_basis(j))**2)*sqrt(PI)*(1.0_dp + alpha_basis*(this%r_basis(i) + this%r_basis(j))**2 ) * &
               ( 1.0_dp + erf( sqrt(alpha_basis/2.0_dp) * (this%r_basis(i) + this%r_basis(j)) ) ) )
         enddo
      enddo

      !overlap_basis = overlap_basis * sqrt(pi / ( 8.0_dp * alpha_basis ) )
      overlap_basis = overlap_basis / sqrt(128.0_dp * alpha_basis**5)

      call initialise(LA_covariance_basis,covariance_basis)
      call initialise(LA_overlap_basis,overlap_basis)
      call LA_Matrix_Factorise(LA_overlap_basis, this%cholesky_overlap_basis)
      do i = 1, this%n_max
         do j = 1, i-1 !i + 1, this%n_max
            this%cholesky_overlap_basis(j,i) = 0.0_dp
         enddo
      enddo

      call Matrix_Solve(LA_covariance_basis,this%cholesky_overlap_basis,this%transform_basis)

      call finalise(LA_covariance_basis)
      call finalise(LA_overlap_basis)

      if(allocated(covariance_basis)) deallocate(covariance_basis)
      if(allocated(overlap_basis)) deallocate(overlap_basis)
      if(allocated(cholesky_overlap_basis)) deallocate(cholesky_overlap_basis)

      this%initialised = .true.

   endsubroutine soap_initialise

   subroutine soap_finalise(this,error)
      type(soap), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return
      this%cutoff = 0.0_dp
      this%cutoff_transition_width = 0.0_dp
      this%l_max = 0
      this%alpha = 0.0_dp

      this%n_max = 0

      if(allocated(this%r_basis)) deallocate(this%r_basis)
      if(allocated(this%transform_basis)) deallocate(this%transform_basis)
      if(allocated(this%cholesky_overlap_basis)) deallocate(this%cholesky_overlap_basis)

      this%initialised = .false.

   endsubroutine soap_finalise

   subroutine descriptor_str_add_species(this,species,descriptor_str,error)
      character(len=*), intent(in) :: this
      integer, dimension(:), intent(in) :: species
      character(len=STRING_LENGTH), dimension(:), allocatable, intent(out) :: descriptor_str
      integer, optional, intent(out) :: error

      integer :: my_descriptor_type, i, j, k, l, n_species
      real(dp), dimension(:), allocatable :: w

      INIT_ERROR(error)

      if(allocated(descriptor_str)) deallocate(descriptor_str)

      my_descriptor_type = get_descriptor_type(this,error)
      n_species = size(species)

      select case(my_descriptor_type)
      case(DT_BISPECTRUM_SO4,DT_BISPECTRUM_SO3,DT_BEHLER,DT_COSNX,DT_POWER_SO3,DT_POWER_SO4)
         allocate(w(n_species))
         allocate(descriptor_str(n_species))

         if( n_species == 1 ) then
            w = 1.0_dp
         else
            w = real( (/ (i, i=0, n_species-1) /), kind=dp ) / (n_species-1) * 0.5_dp + 0.5_dp
         endif

         do i = 1, n_species
            descriptor_str(i) = trim(this)//" n_species="//n_species//" Z="//species(i)//" species_Z={"//species//"} w={"//w//"}"
         enddo

         deallocate(w)
      case(DT_DISTANCE_2B,DT_CO_DISTANCE_2B)
         allocate(descriptor_str(n_species * (n_species+1) / 2))

         l = 0
         do i = 1, n_species
            do j = i, n_species
               l = l + 1
               descriptor_str(l) = trim(this)//" Z1="//species(i)//" Z2="//species(j)
            enddo
         enddo

      case(DT_COORDINATION)
         allocate(descriptor_str(n_species))
         do i = 1, n_species
            descriptor_str(i) = trim(this)//" Z="//species(i)
         enddo
      case(DT_ANGLE_3B,DT_CO_ANGLE_3B)
         allocate(descriptor_str(n_species * n_species * (n_species+1) / 2))
         l = 0
         do i = 1, n_species
            do j = 1, n_species
               do k = j, n_species
                  l = l + 1
                  descriptor_str(l) = trim(this)//" Z="//species(i)//" Z1="//species(j)//" Z2="//species(k)
               enddo
            enddo
         enddo
      case(DT_WATER_MONOMER,DT_WATER_DIMER,DT_A2_DIMER,DT_AB_DIMER,DT_TRIHIS,DT_BOND_REAL_SPACE,DT_ATOM_REAL_SPACE,DT_SOAP)
         allocate(descriptor_str(1))
         descriptor_str(1) = trim(this)
      case default
         RAISE_ERROR("descriptor_str_add_species: unknown descriptor type "//my_descriptor_type,error)
      endselect

   endsubroutine descriptor_str_add_species

   subroutine descriptor_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(descriptor), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            call calc(this%descriptor_bispectrum_SO4,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_BISPECTRUM_SO3)
            call calc(this%descriptor_bispectrum_SO3,at,descriptor_out,do_descriptor,do_grad_descriptor,error=error)
         case(DT_BEHLER)
            call calc(this%descriptor_behler,at,descriptor_out,do_descriptor,do_grad_descriptor,error=error)
         case(DT_DISTANCE_2b)
            call calc(this%descriptor_distance_2b,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_COORDINATION)
            call calc(this%descriptor_coordination,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_ANGLE_3B)
            call calc(this%descriptor_angle_3b,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_CO_ANGLE_3B)
            call calc(this%descriptor_co_angle_3b,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_CO_DISTANCE_2b)
            call calc(this%descriptor_co_distance_2b,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_COSNX)
            call calc(this%descriptor_cosnx,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_TRIHIS)
            call calc(this%descriptor_trihis,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_WATER_MONOMER)
            call calc(this%descriptor_water_monomer,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_WATER_DIMER)
            call calc(this%descriptor_water_dimer,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_A2_DIMER)
            call calc(this%descriptor_A2_dimer,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_AB_DIMER)
            call calc(this%descriptor_AB_dimer,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_BOND_REAL_SPACE)
            call calc(this%descriptor_bond_real_space,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_ATOM_REAL_SPACE)
            call calc(this%descriptor_atom_real_space,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_POWER_SO3)
            call calc(this%descriptor_power_so3,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_POWER_SO4)
            call calc(this%descriptor_power_so4,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case(DT_SOAP)
            call calc(this%descriptor_soap,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
         case default
            RAISE_ERROR("descriptor_calc: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endsubroutine descriptor_calc

   subroutine bispectrum_SO4_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(bispectrum_SO4), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      type(cplx_2d), dimension(:), allocatable :: U
      type(cplx_3d), dimension(:,:), allocatable :: dU

      complex(dp) :: sub
      complex(dp), dimension(3) :: dsub
      real(dp), dimension(3) :: diff, u_ij
      real(dp) :: r, tmp_cg
      integer :: i, n, n_i, ji, jn, j, m1, m2, j1, j2, m11, m12, m21, m22, i_desc, i_bisp, d, n_descriptors, n_cross, n_neighbours
      integer, dimension(3) :: shift
      integer, dimension(116) :: species_map
      logical :: my_do_descriptor, my_do_grad_descriptor

      INIT_ERROR(error)

      call system_timer('bispectrum_SO4_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO4_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      if( .not. cg_initialised ) then
         call cg_initialise(this%j_max,2)
      elseif( this%j_max > cg_j_max ) then
         call cg_finalise()
         call cg_initialise(this%j_max,2)
      endif

      call finalise(descriptor_out)

      d = bispectrum_SO4_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)
      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle

         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif

         if(my_do_grad_descriptor) then
            n_neighbours = atoms_n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif

      enddo

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         if(my_do_descriptor) descriptor_out%x(i_desc)%has_data = .true.
         if(my_do_grad_descriptor) then
            ! dU is not allocated, allocate and zero it
            allocate( dU(0:this%j_max,0:atoms_n_neighbours(at,i,max_dist=this%cutoff)) )
            do j = 0, this%j_max
               allocate( dU(j,0)%mm(3,-j:j,-j:j) )
               dU(j,0)%mm = CPLX_ZERO
            enddo

            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, atoms_n_neighbours(at,i)
            ji = atoms_neighbour(at, i, n, jn=jn, distance=r, diff=diff, cosines=u_ij,shift=shift)
            if( r > this%cutoff ) cycle

            n_i = n_i + 1

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = ji
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,ji) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif
         enddo

         if(my_do_grad_descriptor) then
            call fourier_SO4_calc(this%fourier_SO4,at,i,U,dU,error=error)
         else
            call fourier_SO4_calc(this%fourier_SO4,at,i,U,error=error)
         endif

         if(my_do_descriptor) then

            i_bisp = 0
            do j1 = 0, this%j_max
               j2 = j1
               !do j2 = 0, this%j_max
                  do j = abs(j1-j2), min(this%j_max,j1+j2)
                     if( mod(j1+j2+j,2) == 1 ) cycle
                     
                     i_bisp = i_bisp + 1

                     !do m1 = -j, j, 2
                     !   do m2 = -j, j, 2
                     !      sub = CPLX_ZERO
                     !      do m11 = max(-j1-m1,-j1), min(j1-m1,j1), 2
                     !         do m21 = max(-j2-m2,-j2), min(j2-m2,j2), 2
                     !            sub = sub + cg_array(j1,m11,j,m1,j1,m11+m1) &
                     !            * cg_array(j2,m21,j,m2,j2,m21+m2) &
                     !            * U(j1)%mm(m11,m11+m1) * U(j2)%mm(m21,m21+m2)
                     !         enddo
                     !      enddo
                     !      descriptor_out%x(i_desc)%data(i_bisp) = descriptor_out%x(i_desc)%data(i_bisp) + sub*conjg(U(j)%mm(-m2,m1))*(-1)**(m2/2)
                     !   enddo
                     !enddo

                     do m1 = -j, j, 2
                        do m2 = -j, j, 2
                           sub = CPLX_ZERO
                           do m11 = max(-j1,m1-j2), min(j1,m1+j2), 2
                              do m12 = max(-j1,m2-j2), min(j1,m2+j2), 2
                                 sub = sub + cg_array(j1,m11,j2,m1-m11,j,m1) &
                                 * cg_array(j1,m12,j2,m2-m12,j,m2) &
                                 * U(j1)%mm(m11,m12) * U(j2)%mm(m1-m11,m2-m12)
                              enddo
                           enddo
                           descriptor_out%x(i_desc)%data(i_bisp) = descriptor_out%x(i_desc)%data(i_bisp) + sub*conjg(U(j)%mm(m1,m2))
                        enddo
                     enddo

                  enddo
               !enddo
            enddo
         endif

         if(my_do_grad_descriptor) then
            n_i = 0
            do n = 0, atoms_n_neighbours(at,i)
               if( n>0 ) then
                  ji = atoms_neighbour(at, i, n, distance=r)
                  if( r > this%cutoff ) cycle
                  n_i = n_i + 1
               endif
               i_bisp = 0
               do j1 = 0, this%j_max
                  j2 = j1
                  !do j2 = 0, this%j_max
                     do j = abs(j1-j2), min(this%j_max,j1+j2)
                        if( mod(j1+j2+j,2) == 1 ) cycle

                        i_bisp = i_bisp + 1

                        !do m1 = -j, j, 2
                        !   do m2 = -j, j, 2
                        !      sub = CPLX_ZERO
                        !      dsub = CPLX_ZERO

                        !      do m11 = max(-j1-m1,-j1), min(j1-m1,j1), 2
                        !         do m21 = max(-j2-m2,-j2), min(j2-m2,j2), 2
                        !            tmp_cg =  cg_array(j1,m11,j,m1,j1,m11+m1) &
                        !              * cg_array(j2,m21,j,m2,j2,m21+m2)

                        !            sub = sub + tmp_cg &
                        !            * U(j1)%mm(m11,m1+m11) * U(j2)%mm(m21,m2+m21)
                        !            dsub = dsub + tmp_cg &
                        !            * ( dU(j1,n_i)%mm(:,m11,m1+m11) * U(j2)%mm(m21,m2+m21) + &
                        !            U(j1)%mm(m11,m1+m11) * dU(j2,n_i)%mm(:,m21,m2+m21) )
                        !         enddo
                        !      enddo
                        !      descriptor_out%x(i_desc)%grad_data(i_bisp,:,n_i) = &
                        !      descriptor_out%x(i_desc)%grad_data(i_bisp,:,n_i) + &
                        !      ( dsub*conjg(U(j)%mm(-m2,m1)) + sub*conjg(dU(j,n_i)%mm(:,-m2,m1)) )*(-1)**(m2/2)
                        !   enddo
                        !enddo
                        do m1 = -j, j, 2
                           do m2 = -j, j, 2
                              sub = CPLX_ZERO
                              dsub = CPLX_ZERO
                              do m11 = max(-j1,m1-j2), min(j1,m1+j2), 2
                                 do m12 = max(-j1,m2-j2), min(j1,m2+j2), 2
                                    
                                    tmp_cg =  cg_array(j1,m11,j2,m1-m11,j,m1) &
                                    * cg_array(j1,m12,j2,m2-m12,j,m2)

                                    sub = sub + tmp_cg &
                                    * U(j1)%mm(m11,m12) * U(j2)%mm(m1-m11,m2-m12)
                                    dsub = dsub + tmp_cg &
                                    * ( dU(j1,n_i)%mm(:,m11,m12) * U(j2)%mm(m1-m11,m2-m12) + &
                                    U(j1)%mm(m11,m12) * dU(j2,n_i)%mm(:,m1-m11,m2-m12) )
                                 enddo
                              enddo
                              descriptor_out%x(i_desc)%grad_data(i_bisp,:,n_i) = &
                              descriptor_out%x(i_desc)%grad_data(i_bisp,:,n_i) + &
                              dsub*conjg(U(j)%mm(m1,m2)) + sub*conjg(dU(j,n_i)%mm(:,m1,m2))
                           enddo
                        enddo

                     enddo
                  !enddo
               enddo 
            enddo
         endif

         call finalise(dU)
      enddo ! i

      ! clear U from the memory
      call finalise(U)

      call system_timer('bispectrum_SO4_calc')

   endsubroutine bispectrum_SO4_calc

   subroutine bispectrum_so3_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(bispectrum_so3), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      type(cplx_1d), dimension(:), allocatable :: SphericalY_ij
      type(cplx_1d), dimension(:,:), allocatable :: fourier_so3

      type(cplx_2d), dimension(:), allocatable :: dSphericalY_ij
      type(cplx_2d), dimension(:,:,:), allocatable :: dfourier_so3

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, n, a, l, m, l1, l2, m1, i_desc, i_pow, n_neighbours, n_i, n_descriptors, n_cross
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij, d_ij
      real(dp), dimension(:), allocatable :: Rad_ij
      real(dp), dimension(:,:), allocatable :: dRad_ij

      complex(dp) :: sub, dsub(3)

      integer, dimension(116) :: species_map

      INIT_ERROR(error)

      call system_timer('bispectrum_so3_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_so3_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      if( .not. cg_initialised ) then
         call cg_initialise(this%l_max)
      elseif( this%l_max > cg_j_max ) then
         call cg_finalise()
         call cg_initialise(this%l_max)
      endif

      call finalise(descriptor_out)

      d = bispectrum_so3_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            n_neighbours = atoms_n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(fourier_so3(0:this%l_max,this%n_max),SphericalY_ij(0:this%l_max),Rad_ij(this%n_max))
      do a = 1, this%n_max
         do l = 0, this%l_max
            allocate(fourier_so3(l,a)%m(-l:l))
            fourier_so3(l,a)%m(:) = CPLX_ZERO
         enddo
      enddo
      do l = 0, this%l_max
         allocate(SphericalY_ij(l)%m(-l:l))
      enddo

      if(my_do_grad_descriptor) then
         allocate( dRad_ij(3,this%n_max), dSphericalY_ij(0:this%l_max) )
         do l = 0, this%l_max
            allocate(dSphericalY_ij(l)%mm(3,-l:l))
         enddo
      endif

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         do a = 1, this%n_max
            do l = 0, this%l_max
               fourier_so3(l,a)%m(:) = CPLX_ZERO
            enddo
         enddo

         if(my_do_descriptor) descriptor_out%x(i_desc)%has_data = .true.
         if(my_do_grad_descriptor) then
            allocate( dfourier_so3(0:this%l_max,this%n_max,0:atoms_n_neighbours(at,i,max_dist=this%cutoff)) )
            do n = 0, atoms_n_neighbours(at,i,max_dist=this%cutoff)
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     allocate(dfourier_so3(l,a,n)%mm(3,-l:l))
                     dfourier_so3(l,a,n)%mm(:,:) = CPLX_ZERO
                  enddo
               enddo
            enddo
         endif

         n_i = 0
         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij > this%cutoff ) cycle

            n_i = n_i + 1
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif

            do a = 1, this%n_max
               Rad_ij(a) = RadialFunction(this%Radial, r_ij, a)
               if(my_do_grad_descriptor) dRad_ij(:,a) = GradRadialFunction(this%Radial, r_ij, a) * u_ij
            enddo

            do l = 0, this%l_max
               do m = -l, l
                  SphericalY_ij(l)%m(m) = SphericalYCartesian(l,m,d_ij)
                  if(my_do_grad_descriptor) dSphericalY_ij(l)%mm(:,m) = GradSphericalYCartesian(l,m,d_ij)
               enddo
            enddo

            do a = 1, this%n_max
               do l = 0, this%l_max
                  do m = -l, l
                     fourier_so3(l,a)%m(m) = fourier_so3(l,a)%m(m) + Rad_ij(a)*SphericalY_ij(l)%m(m)
                     if(my_do_grad_descriptor) then
                        dfourier_so3(l,a,n_i)%mm(:,m) = dfourier_so3(l,a,n_i)%mm(:,m) + &
                        dRad_ij(:,a) * SphericalY_ij(l)%m(m) + Rad_ij(a)*dSphericalY_ij(l)%mm(:,m)
                     endif
                  enddo
               enddo
            enddo

         enddo ! n

         if(my_do_descriptor) then
            i_pow = 0
            do a = 1, this%n_max
               do l1 = 0, this%l_max
                  l2 = l1
                  !do l2 = 0, this%l_max
                     do l = abs(l1-l2), min(this%l_max,l1+l2)
                        if( mod(l1,2)==1 .and. mod(l2,2)==1 .and. mod(l,2)==1 ) cycle
                        i_pow = i_pow + 1

                        do m = -l, l
                           sub = CPLX_ZERO
                           do m1 = max(-l1,m-l2),min(l1,m+l2)
                              sub = sub + cg_array(l1,m1,l2,m-m1,l,m) * conjg(fourier_so3(l1,a)%m(m1)) * conjg(fourier_so3(l2,a)%m(m-m1))
                           enddo

                           descriptor_out%x(i_desc)%data(i_pow) = descriptor_out%x(i_desc)%data(i_pow) + fourier_so3(l,a)%m(m) * sub
                        enddo

                     enddo
                  !enddo
               enddo
            enddo
         endif

         if(my_do_grad_descriptor) then
            do n = 1, atoms_n_neighbours(at,i,max_dist=this%cutoff)
               i_pow = 0
               do a = 1, this%n_max
                  do l1 = 0, this%l_max
                     l2 = l1
                     !do l2 = 0, this%l_max
                        do l = abs(l1-l2), min(this%l_max,l1+l2)
                           if( mod(l1,2)==1 .and. mod(l2,2)==1 .and. mod(l,2)==1 ) cycle
                           i_pow = i_pow + 1

                           do m = -l, l
                              sub = CPLX_ZERO
                              dsub = CPLX_ZERO
                              do m1 = max(-l1,m-l2),min(l1,m+l2)
                                 dsub = dsub + cg_array(l1,m1,l2,m-m1,l,m) * &
                                 ( conjg(dfourier_so3(l1,a,n)%mm(:,m1)) * conjg(fourier_so3(l2,a)%m(m-m1)) + &
                                   conjg(fourier_so3(l1,a)%m(m1)) * conjg(dfourier_so3(l2,a,n)%mm(:,m-m1)) )
                                 sub = sub + cg_array(l1,m1,l2,m-m1,l,m) * conjg(fourier_so3(l1,a)%m(m1)) * conjg(fourier_so3(l2,a)%m(m-m1))
                              enddo

                              descriptor_out%x(i_desc)%grad_data(i_pow,:,n) = descriptor_out%x(i_desc)%grad_data(i_pow,:,n) + &
                              fourier_so3(l,a)%m(m) * dsub + dfourier_so3(l,a,n)%mm(:,m) * sub
                           enddo
                        enddo
                     !enddo
                  enddo
               enddo
               descriptor_out%x(i_desc)%grad_data(:,:,0) = descriptor_out%x(i_desc)%grad_data(:,:,0) - descriptor_out%x(i_desc)%grad_data(:,:,n)
            enddo
         endif

         if(allocated(dfourier_so3)) then
            do n = lbound(dfourier_so3,3), ubound(dfourier_so3,3)
               do a = lbound(dfourier_so3,2), ubound(dfourier_so3,2)
                  do l = lbound(dfourier_so3,1), ubound(dfourier_so3,1)
                     deallocate(dfourier_so3(l,a,n)%mm)
                  enddo
               enddo
            enddo
            deallocate(dfourier_so3)
         endif

      enddo ! i

      if(allocated(Rad_ij)) deallocate(Rad_ij)
      if(allocated(dRad_ij)) deallocate(dRad_ij)

      if(allocated(fourier_so3)) then
         do a = lbound(fourier_so3,2), ubound(fourier_so3,2)
            do l = lbound(fourier_so3,1), ubound(fourier_so3,1)
               deallocate(fourier_so3(l,a)%m)
            enddo
         enddo
         deallocate(fourier_so3)
      endif

      if(allocated(SphericalY_ij)) then
         do l = lbound(SphericalY_ij,1), ubound(SphericalY_ij,1)
            deallocate(SphericalY_ij(l)%m)
         enddo
         deallocate(SphericalY_ij)
      endif

      if(allocated(dSphericalY_ij)) then
         do l = lbound(dSphericalY_ij,1), ubound(dSphericalY_ij,1)
            deallocate(dSphericalY_ij(l)%mm)
         enddo
         deallocate(dSphericalY_ij)
      endif

      call system_timer('bispectrum_so3_calc')

   endsubroutine bispectrum_so3_calc

   subroutine behler_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(behler), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, k, n, m, a, b, i_desc, n_neighbours, n_i, m_i, n_descriptors, n_cross
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, Ang, dAng, Rad, dRad_ij, dRad_ik, dRad_jk, f_cut_ij, f_cut_ik, f_cut_jk, df_cut_ij, df_cut_ik, df_cut_jk, g2, dg2
      real(dp), dimension(3) :: u_ij, u_ik, u_jk, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik

      INIT_ERROR(error)

      call system_timer('behler_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("behler_calc: descriptor object not initialised", error)
      endif

      if( at%cutoff < this%cutoff ) then
         RAISE_ERROR("behler_calc: cutoff of atoms object ("//at%cutoff//") less than cutoff of descriptor ("//this%cutoff//")", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      d = behler_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            n_neighbours = atoms_n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N
         i_desc = i_desc + 1

         if(my_do_descriptor) descriptor_out%x(i_desc)%has_data = .true.
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij > this%cutoff ) cycle

            n_i = n_i + 1

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif

            f_cut_ij = cutoff_function(r_ij,this%cutoff)
            if(my_do_grad_descriptor) df_cut_ij = dcutoff_function(r_ij,this%cutoff)

            do a = 1, this%n_g2
               g2 = exp(-this%g2(a)%eta * (r_ij-this%g2(a)%rs)**2)
               if(my_do_descriptor) descriptor_out%x(i_desc)%data(a) = descriptor_out%x(i_desc)%data(a) + g2 * f_cut_ij
               if(my_do_grad_descriptor) then
                  dg2 = -2.0_dp * this%g2(a)%eta * (r_ij-this%g2(a)%rs) * g2
                  descriptor_out%x(i_desc)%grad_data(a,:,n_i) = ( dg2 * f_cut_ij + g2 * df_cut_ij ) * u_ij
                  descriptor_out%x(i_desc)%grad_data(a,:,0) = descriptor_out%x(i_desc)%grad_data(a,:,0) - descriptor_out%x(i_desc)%grad_data(a,:,n_i)
               endif
            enddo


            m_i = 0
            do m = 1, atoms_n_neighbours(at,i)
               k = atoms_neighbour(at, i, m, distance = r_ik, cosines=u_ik, diff=d_ik)
               if( r_ik > this%cutoff ) cycle

               m_i = m_i + 1

               d_jk = d_ik - d_ij
               r_jk = norm(d_jk)
               if( r_jk .feq. 0.0_dp ) cycle

               u_jk = d_jk / r_jk

               f_cut_ik = cutoff_function(r_ik,this%cutoff)
               f_cut_jk = cutoff_function(r_jk,this%cutoff)

               cos_ijk = dot_product(u_ij,u_ik)

               if(my_do_grad_descriptor) then
                  df_cut_ik = dcutoff_function(r_ik,this%cutoff)
                  df_cut_jk = dcutoff_function(r_jk,this%cutoff)

                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik
               endif

               do b = 1, this%n_g3
                  a = b + this%n_g2

                  Ang = (1.0_dp + this%g3(b)%lambda * cos_ijk)**this%g3(b)%zeta
                  Rad = exp( -this%g3(b)%eta * (r_ij**2 + r_ik**2 + r_jk**2) )

                  if(my_do_descriptor) descriptor_out%x(i_desc)%data(a) = descriptor_out%x(i_desc)%data(a) + 0.5_dp * Ang * Rad * f_cut_ij * f_cut_ik * f_cut_jk
                  if(my_do_grad_descriptor) then
                     dAng = this%g3(b)%zeta * (1.0_dp + this%g3(b)%lambda * cos_ijk)**(this%g3(b)%zeta -1.0_dp) * this%g3(b)%lambda
                     dRad_ij = -this%g3(b)%eta * 2.0_dp * r_ij * Rad
                     dRad_ik = -this%g3(b)%eta * 2.0_dp * r_ik * Rad
                     dRad_jk = -this%g3(b)%eta * 2.0_dp * r_jk * Rad

                     descriptor_out%x(i_desc)%grad_data(a,:,n_i) = descriptor_out%x(i_desc)%grad_data(a,:,n_i) + 0.5_dp * &
                     ( ( dAng * dcosijk_ij * Rad + Ang * ( dRad_ij * u_ij - dRad_jk * u_jk ) ) * f_cut_ij * f_cut_ik * f_cut_jk + &
                     Ang * Rad * f_cut_ik * ( df_cut_ij * u_ij * f_cut_jk - f_cut_ij * df_cut_jk * u_jk ) )

                     descriptor_out%x(i_desc)%grad_data(a,:,m_i) = descriptor_out%x(i_desc)%grad_data(a,:,m_i) + 0.5_dp * &
                     ( ( dAng * dcosijk_ik * Rad + Ang * ( dRad_ik * u_ik + dRad_jk * u_jk ) ) * f_cut_ij * f_cut_ik * f_cut_jk + &
                     Ang * Rad * f_cut_ij * ( df_cut_ik * u_ik * f_cut_jk + f_cut_ik * df_cut_jk * u_jk ) ) 

                     descriptor_out%x(i_desc)%grad_data(a,:,0) = descriptor_out%x(i_desc)%grad_data(a,:,0) - 0.5_dp * &
                     ( ( dAng * (dcosijk_ij+dcosijk_ik) * Rad + Ang * (dRad_ij * u_ij + dRad_ik * u_ik) ) * f_cut_ij * f_cut_ik * f_cut_jk + &
                     Ang * Rad * f_cut_jk * ( df_cut_ij * u_ij * f_cut_ik + f_cut_ij * df_cut_ik * u_ik ) )
                  endif


               enddo

            enddo
         enddo

         do b = 1, this%n_g3
            a = b + this%n_g2

            if(my_do_descriptor) descriptor_out%x(i_desc)%data(a) = descriptor_out%x(i_desc)%data(a) * 2.0_dp**(1.0_dp-this%g3(b)%zeta) 
            if(my_do_grad_descriptor) descriptor_out%x(i_desc)%grad_data(a,:,:) = descriptor_out%x(i_desc)%grad_data(a,:,:) * 2.0_dp**(1.0_dp-this%g3(b)%zeta)
         enddo
      enddo

      call system_timer('behler_calc')

   endsubroutine behler_calc

   subroutine distance_2b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      logical :: my_do_descriptor, my_do_grad_descriptor, Zi1, Zi2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, n
      integer, dimension(3) :: shift
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij

      INIT_ERROR(error)

      call system_timer('distance_2b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("distance_2b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      d = distance_2b_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%has_data = .false.
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,0:1))
            allocate(descriptor_out%x(i)%ii(0:1))
            allocate(descriptor_out%x(i)%pos(3,0:1))
            allocate(descriptor_out%x(i)%has_grad_data(0:1))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,0:1))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)
            if( r_ij > this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type


            i_desc = i_desc + 1
            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%data(1) = r_ij
               descriptor_out%x(i_desc)%has_data = .true.

               descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(r_ij,this%cutoff,0.5_dp)
            endif
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(0) = i
               descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
               descriptor_out%x(i_desc)%has_grad_data(0) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,0) = -u_ij(:)
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = -dcoordination_function(r_ij,this%cutoff,0.5_dp)*u_ij

               descriptor_out%x(i_desc)%ii(1) = j
               descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(1) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,1) = u_ij(:)
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = -descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0)

            endif
         enddo
      enddo

      call system_timer('distance_2b_calc')

   endsubroutine distance_2b_calc

   subroutine coordination_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(coordination), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, n, i_n, n_neighbours, i_desc, n_descriptors, n_cross
      integer, dimension(3) :: shift
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij, df_cut

      INIT_ERROR(error)

      call system_timer('coordination_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("coordination_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      d = coordination_dimensions(this,error)

      call descriptor_sizes(this,at,n_descriptors,n_cross,error)
      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle

         i_desc = i_desc + 1
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            descriptor_out%x(i_desc)%has_data = .false.

            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            n_neighbours = atoms_n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         if(my_do_descriptor) descriptor_out%x(i_desc)%has_data = .true.
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         i_n = 0
         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)

            if( r_ij > this%cutoff ) cycle
            i_n = i_n + 1

            if(my_do_descriptor) &
               descriptor_out%x(i_desc)%data(1) = descriptor_out%x(i_desc)%data(1) + coordination_function(r_ij,this%cutoff,this%transition_width)

            if(my_do_grad_descriptor) then
               df_cut = dcoordination_function(r_ij,this%cutoff,this%transition_width) * u_ij

               descriptor_out%x(i_desc)%grad_data(1,:,0) = descriptor_out%x(i_desc)%grad_data(1,:,0) - df_cut

               descriptor_out%x(i_desc)%ii(i_n) = j
               descriptor_out%x(i_desc)%pos(:,i_n) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(i_n) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,i_n) = df_cut
            endif
         enddo
      enddo

      call system_timer('coordination_calc')

   endsubroutine coordination_calc

   subroutine angle_3b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(angle_3b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      logical :: my_do_descriptor, my_do_grad_descriptor, Zk1, Zk2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, k, n, m
      integer, dimension(3) :: shift_ij, shift_ik
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, fc_j, fc_k, dfc_j, dfc_k
      real(dp), dimension(3) :: u_ij, u_ik, u_jk, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik

      INIT_ERROR(error)

      call system_timer('angle_3b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("angle_3b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      d = angle_3b_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%has_data = .false.
         endif

         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,0:2))
            allocate(descriptor_out%x(i)%ii(0:2))
            allocate(descriptor_out%x(i)%pos(3,0:2))
            allocate(descriptor_out%x(i)%has_grad_data(0:2))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,0:2))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle
         
         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij, cosines = u_ij, diff = d_ij, shift=shift_ij)

            if( r_ij > this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            fc_j = coordination_function(r_ij,this%cutoff,0.5_dp)
            dfc_j = dcoordination_function(r_ij,this%cutoff,0.5_dp)

            do m = 1, atoms_n_neighbours(at,i)

               if( n == m ) cycle

               k = atoms_neighbour(at, i, m, distance = r_ik, cosines = u_ik, diff = d_ik, shift=shift_ik)
               if( r_ik > this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)

               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               d_jk = d_ij - d_ik
               r_jk = norm(d_jk)
               u_jk = d_jk / r_jk

               fc_k = coordination_function(r_ik,this%cutoff,0.5_dp)
               dfc_k = dcoordination_function(r_ik,this%cutoff,0.5_dp)

               cos_ijk = dot_product(d_ij,d_ik)/(r_ij*r_ik)

               i_desc = i_desc + 1

               if(my_do_descriptor) then
                  descriptor_out%x(i_desc)%data(1) = r_ij + r_ik
                  descriptor_out%x(i_desc)%data(2) = (r_ij - r_ik)**2
                  descriptor_out%x(i_desc)%data(3) = r_jk !cos_ijk
                  descriptor_out%x(i_desc)%has_data = .true.

                  descriptor_out%x(i_desc)%covariance_cutoff = fc_j*fc_k
               endif

               if(my_do_grad_descriptor) then
                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik

                  descriptor_out%x(i_desc)%ii(0) = i
                  descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
                  descriptor_out%x(i_desc)%has_grad_data(0) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,0) = - u_ij - u_ik
                  descriptor_out%x(i_desc)%grad_data(2,:,0) = 2.0_dp * (r_ij - r_ik)*(-u_ij + u_ik)
                  descriptor_out%x(i_desc)%grad_data(3,:,0) = 0.0_dp !-dcosijk_ij - dcosijk_ik

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = - dfc_j*fc_k*u_ij - dfc_k*fc_j*u_ik

                  descriptor_out%x(i_desc)%ii(1) = j
                  descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift_ij)
                  descriptor_out%x(i_desc)%has_grad_data(1) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,1) = u_ij
                  descriptor_out%x(i_desc)%grad_data(2,:,1) = 2.0_dp * (r_ij - r_ik)*u_ij
                  descriptor_out%x(i_desc)%grad_data(3,:,1) = u_jk !dcosijk_ij

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = dfc_j*fc_k*u_ij

                  descriptor_out%x(i_desc)%ii(2) = k
                  descriptor_out%x(i_desc)%pos(:,2) = at%pos(:,k) + matmul(at%lattice,shift_ik)
                  descriptor_out%x(i_desc)%has_grad_data(2) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,2) = u_ik
                  descriptor_out%x(i_desc)%grad_data(2,:,2) = 2.0_dp * (r_ij - r_ik)*(-u_ik)
                  descriptor_out%x(i_desc)%grad_data(3,:,2) = -u_jk !dcosijk_ik

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,2) = dfc_k*fc_j*u_ik
               endif
            enddo
         enddo
      enddo

      call system_timer('angle_3b_calc')

   endsubroutine angle_3b_calc

   subroutine co_angle_3b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(co_angle_3b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      type(descriptor) :: my_coordination
      type(descriptor_data) :: descriptor_coordination

      logical :: my_do_descriptor, my_do_grad_descriptor, Zk1, Zk2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, k, n, m, n_neighbours_coordination
      integer, dimension(3) :: shift_ij, shift_ik
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, fc_j, fc_k, dfc_j, dfc_k
      real(dp), dimension(3) :: u_ij, u_ik, u_jk, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik

      INIT_ERROR(error)

      call system_timer('co_angle_3b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("co_angle_3b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      d = co_angle_3b_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0
      do i = 1, at%N
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle
         n_neighbours_coordination = atoms_n_neighbours(at,i,max_dist=this%coordination_cutoff)

         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij)
            if( r_ij > this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            do m = 1, atoms_n_neighbours(at,i)

               if( n == m ) cycle

               k = atoms_neighbour(at, i, m, distance = r_ik)
               if( r_ik > this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)

               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               i_desc = i_desc + 1
               if(my_do_descriptor) then
                  allocate(descriptor_out%x(i_desc)%data(d))
                  descriptor_out%x(i_desc)%data = 0.0_dp
                  descriptor_out%x(i_desc)%has_data = .false.
               endif

               if(my_do_grad_descriptor) then

                  allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:2+n_neighbours_coordination))
                  allocate(descriptor_out%x(i_desc)%ii(0:2+n_neighbours_coordination))
                  allocate(descriptor_out%x(i_desc)%pos(3,0:2+n_neighbours_coordination))
                  allocate(descriptor_out%x(i_desc)%has_grad_data(0:2+n_neighbours_coordination))
                  descriptor_out%x(i_desc)%grad_data = 0.0_dp
                  descriptor_out%x(i_desc)%ii = 0
                  descriptor_out%x(i_desc)%pos = 0.0_dp
                  descriptor_out%x(i_desc)%has_grad_data = .false.

                  allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:2+n_neighbours_coordination))
                  descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
               endif
            enddo
         enddo
      enddo

      call initialise(my_coordination,'coordination cutoff='//this%coordination_cutoff//' coordination_transition_width='//this%coordination_transition_width,error)
      call calc(my_coordination,at,descriptor_coordination,do_descriptor,do_grad_descriptor,error)
      
      i_desc = 0
      do i = 1, at%N
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle

         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij, cosines = u_ij, diff = d_ij, shift=shift_ij)

            if( r_ij > this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            fc_j = coordination_function(r_ij,this%cutoff,0.5_dp)
            dfc_j = dcoordination_function(r_ij,this%cutoff,0.5_dp)

            do m = 1, atoms_n_neighbours(at,i)
               if( n == m ) cycle

               k = atoms_neighbour(at, i, m, distance = r_ik, cosines = u_ik, diff = d_ik, shift=shift_ik)
               if( r_ik > this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)

               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               d_jk = d_ij - d_ik
               r_jk = norm(d_jk)
               u_jk = d_jk / r_jk

               fc_k = coordination_function(r_ik,this%cutoff,0.5_dp)
               dfc_k = dcoordination_function(r_ik,this%cutoff,0.5_dp)

               cos_ijk = dot_product(d_ij,d_ik)/(r_ij*r_ik)

               i_desc = i_desc + 1

               if(my_do_descriptor) then
                  descriptor_out%x(i_desc)%data(1) = r_ij + r_ik
                  descriptor_out%x(i_desc)%data(2) = (r_ij - r_ik)**2
                  descriptor_out%x(i_desc)%data(3) = r_jk !cos_ijk
                  descriptor_out%x(i_desc)%data(4) = descriptor_coordination%x(i)%data(1)
                  descriptor_out%x(i_desc)%has_data = .true.

                  descriptor_out%x(i_desc)%covariance_cutoff = fc_j*fc_k
               endif

               if(my_do_grad_descriptor) then
                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik

                  descriptor_out%x(i_desc)%ii(0) = i
                  descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
                  descriptor_out%x(i_desc)%has_grad_data(0) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,0) = - u_ij - u_ik
                  descriptor_out%x(i_desc)%grad_data(2,:,0) = 2.0_dp * (r_ij - r_ik)*(-u_ij + u_ik)
                  descriptor_out%x(i_desc)%grad_data(3,:,0) = 0.0_dp !-dcosijk_ij - dcosijk_ik
                  descriptor_out%x(i_desc)%grad_data(4,:,0) = descriptor_coordination%x(i)%grad_data(1,:,0)

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = - dfc_j*fc_k*u_ij - dfc_k*fc_j*u_ik

                  descriptor_out%x(i_desc)%ii(1) = j
                  descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift_ij)
                  descriptor_out%x(i_desc)%has_grad_data(1) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,1) = u_ij
                  descriptor_out%x(i_desc)%grad_data(2,:,1) = 2.0_dp * (r_ij - r_ik)*u_ij
                  descriptor_out%x(i_desc)%grad_data(3,:,1) = u_jk !dcosijk_ij

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = dfc_j*fc_k*u_ij

                  descriptor_out%x(i_desc)%ii(2) = k
                  descriptor_out%x(i_desc)%pos(:,2) = at%pos(:,k) + matmul(at%lattice,shift_ik)
                  descriptor_out%x(i_desc)%has_grad_data(2) = .true.
                  descriptor_out%x(i_desc)%grad_data(1,:,2) = u_ik
                  descriptor_out%x(i_desc)%grad_data(2,:,2) = 2.0_dp * (r_ij - r_ik)*(-u_ik)
                  descriptor_out%x(i_desc)%grad_data(3,:,2) = -u_jk !dcosijk_ik

                  descriptor_out%x(i_desc)%grad_covariance_cutoff(:,2) = dfc_k*fc_j*u_ik

                  descriptor_out%x(i_desc)%ii(3:) = descriptor_coordination%x(i)%ii(1:)
                  descriptor_out%x(i_desc)%pos(:,3:) = descriptor_coordination%x(i)%pos(:,1:)
                  descriptor_out%x(i_desc)%has_grad_data(3:) = descriptor_coordination%x(i)%has_grad_data(1:)
                  descriptor_out%x(i_desc)%grad_data(4,:,3:) = descriptor_coordination%x(i)%grad_data(1,:,1:)
               endif
            enddo
         enddo
      enddo

      call finalise(my_coordination)
      call finalise(descriptor_coordination)

      call system_timer('co_angle_3b_calc')

   endsubroutine co_angle_3b_calc

   subroutine co_distance_2b_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(co_distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      type(descriptor) :: my_coordination
      type(descriptor_data) :: descriptor_coordination

      logical :: my_do_descriptor, my_do_grad_descriptor, Zi1, Zi2, Zj1, Zj2
      integer :: d, n_descriptors, n_cross, i_desc, i, j, n, n_neighbours_coordination_i, n_neighbours_coordination_ij
      integer, dimension(3) :: shift
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij

      INIT_ERROR(error)

      call system_timer('co_distance_2b_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("co_distance_2b_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      d = co_distance_2b_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)


      allocate(descriptor_out%x(n_descriptors))
      i_desc = 0

      do i = 1, at%N
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance=r_ij)

            if(r_ij > this%cutoff) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            i_desc = i_desc + 1
            if(my_do_descriptor) then
               allocate(descriptor_out%x(i_desc)%data(d))
               descriptor_out%x(i_desc)%data = 0.0_dp
               descriptor_out%x(i_desc)%has_data = .false.
            endif

            if(my_do_grad_descriptor) then
               n_neighbours_coordination_ij = atoms_n_neighbours(at,i,max_dist=this%coordination_cutoff) + &
               atoms_n_neighbours(at,j,max_dist=this%coordination_cutoff) + 2

               allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%ii(0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%pos(3,0:1+n_neighbours_coordination_ij))
               allocate(descriptor_out%x(i_desc)%has_grad_data(0:1+n_neighbours_coordination_ij))
               descriptor_out%x(i_desc)%grad_data = 0.0_dp
               descriptor_out%x(i_desc)%ii = 0
               descriptor_out%x(i_desc)%pos = 0.0_dp
               descriptor_out%x(i_desc)%has_grad_data = .false.

               allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:1+n_neighbours_coordination_ij))
               descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
            endif
         enddo
      enddo

      call initialise(my_coordination,'coordination cutoff='//this%coordination_cutoff//' coordination_transition_width='//this%coordination_transition_width,error)
      call calc(my_coordination,at,descriptor_coordination,.true.,do_grad_descriptor,error)
      
      i_desc = 0
      do i = 1, at%N
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij, cosines = u_ij, shift=shift)
            if( r_ij > this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            i_desc = i_desc + 1
            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%data(1) = r_ij
               descriptor_out%x(i_desc)%has_data = .true.

               descriptor_out%x(i_desc)%data(2) = descriptor_coordination%x(i)%data(1) + descriptor_coordination%x(j)%data(1)
               descriptor_out%x(i_desc)%data(3) = (descriptor_coordination%x(i)%data(1) - descriptor_coordination%x(j)%data(1))**2

               descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(r_ij,this%cutoff,this%transition_width)
            endif
            if(my_do_grad_descriptor) then
               n_neighbours_coordination_i = atoms_n_neighbours(at,i,max_dist=this%coordination_cutoff)

               descriptor_out%x(i_desc)%ii(0) = i
               descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
               descriptor_out%x(i_desc)%has_grad_data(0) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,0) = -u_ij(:)
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0) = -dcoordination_function(r_ij,this%cutoff,this%transition_width)*u_ij

               descriptor_out%x(i_desc)%ii(1) = j
               descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(1) = .true.
               descriptor_out%x(i_desc)%grad_data(1,:,1) = u_ij(:)
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = -descriptor_out%x(i_desc)%grad_covariance_cutoff(:,0)

               descriptor_out%x(i_desc)%ii(2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%ii(:)
               descriptor_out%x(i_desc)%pos(:,2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%pos(:,:)
               descriptor_out%x(i_desc)%has_grad_data(2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%has_grad_data(:)
               descriptor_out%x(i_desc)%grad_data(2,:,2:n_neighbours_coordination_i+2) = descriptor_coordination%x(i)%grad_data(1,:,:)
               descriptor_out%x(i_desc)%grad_data(3,:,2:n_neighbours_coordination_i+2) = 2.0_dp*(descriptor_coordination%x(i)%data(1) - descriptor_coordination%x(j)%data(1))*&
                  descriptor_coordination%x(i)%grad_data(1,:,:)

               descriptor_out%x(i_desc)%ii(n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%ii(:)
               descriptor_out%x(i_desc)%pos(:,n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%pos(:,:)
               descriptor_out%x(i_desc)%has_grad_data(n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%has_grad_data(:)
               descriptor_out%x(i_desc)%grad_data(2,:,n_neighbours_coordination_i+3:) = descriptor_coordination%x(j)%grad_data(1,:,:)
               descriptor_out%x(i_desc)%grad_data(3,:,n_neighbours_coordination_i+3:) = -2.0_dp*(descriptor_coordination%x(i)%data(1) - descriptor_coordination%x(j)%data(1))*&
                  descriptor_coordination%x(j)%grad_data(1,:,:)

            endif
         enddo
      enddo

      call finalise(my_coordination)
      call finalise(descriptor_coordination)

      call system_timer('co_distance_2b_calc')

   endsubroutine co_distance_2b_calc

   subroutine cosnx_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(cosnx), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, k, n, m, a, b, i_desc, i_cosnx, n_neighbours, n_i, n_descriptors, n_cross
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, T_0_cos_ijk, T_1_cos_ijk, T_n_cos_ijk, U_0_cos_ijk, U_1_cos_ijk, U_n_cos_ijk, Ang
      real(dp), dimension(3) :: u_ij, u_ik, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik, dAng_ij, dAng_ik
      real(dp), dimension(:), allocatable :: Rad_ij, Rad_ik, T_cos_ijk, U_cos_ijk
      real(dp), dimension(:,:), allocatable :: dRad_ij, dRad_ik
      integer, dimension(116) :: species_map

      INIT_ERROR(error)

      call system_timer('cosnx_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("cosnx_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call finalise(descriptor_out)

      d = cosnx_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            n_neighbours = atoms_n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(Rad_ij(this%n_max), Rad_ik(this%n_max))
      allocate(T_cos_ijk(0:this%l_max))
      if(my_do_grad_descriptor) then
         allocate(U_cos_ijk(-1:this%l_max))
         allocate(dRad_ij(3,this%n_max), dRad_ik(3,this%n_max))
      endif

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         if(my_do_descriptor) descriptor_out%x(i_desc)%has_data = .true.
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i) 
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij > this%cutoff ) cycle

            n_i = n_i + 1

            do a = 1, this%n_max
               Rad_ij(a) = RadialFunction(this%Radial, r_ij, a) * this%w(species_map(at%Z(j)))
               if(my_do_grad_descriptor) dRad_ij(:,a) = GradRadialFunction(this%Radial, r_ij, a) * u_ij * this%w(species_map(at%Z(j)))
            enddo

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif

            do m = 1, atoms_n_neighbours(at,i)
               k = atoms_neighbour(at, i, m, distance = r_ik, cosines=u_ik, diff=d_ik)
               if( r_ik > this%cutoff ) cycle

               d_jk = d_ik - d_ij
               r_jk = norm(d_jk)
               if( r_jk .feq. 0.0_dp ) cycle

               cos_ijk = dot_product(u_ij,u_ik)
               if(my_do_grad_descriptor) then
                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik
               endif

               do a = 1, this%n_max
                  Rad_ik(a) = RadialFunction(this%Radial, r_ik, a) * this%w(species_map(at%Z(k)))
                  if(my_do_grad_descriptor) dRad_ik(:,a) = GradRadialFunction(this%Radial, r_ik, a) * u_ik * this%w(species_map(at%Z(k)))
               enddo

               if(this%l_max >= 0) then
                  T_cos_ijk(0) = 1.0_dp
                  T_0_cos_ijk = T_cos_ijk(0)
                  if(my_do_grad_descriptor) then
                     U_cos_ijk(-1) = 0.0_dp
                     U_cos_ijk(0) = 1.0_dp
                     U_0_cos_ijk = U_cos_ijk(0)
                  endif
               endif

               if(this%l_max >= 1) then
                  T_cos_ijk(1) = cos_ijk
                  T_1_cos_ijk = T_cos_ijk(1)
                  if(my_do_grad_descriptor) then
                     U_cos_ijk(1) = 2.0_dp*cos_ijk
                     U_1_cos_ijk = U_cos_ijk(1)
                  endif
               endif

               do b = 2, this%l_max
                  T_n_cos_ijk = 2*cos_ijk*T_1_cos_ijk - T_0_cos_ijk
                  T_0_cos_ijk = T_1_cos_ijk
                  T_1_cos_ijk = T_n_cos_ijk

                  T_cos_ijk(b) = T_n_cos_ijk

                  if(my_do_grad_descriptor) then
                     U_n_cos_ijk = 2*cos_ijk*U_1_cos_ijk - U_0_cos_ijk
                     U_0_cos_ijk = U_1_cos_ijk
                     U_1_cos_ijk = U_n_cos_ijk

                     U_cos_ijk(b) = U_n_cos_ijk
                  endif
               enddo
     
               i_cosnx = 0
               do a = 1, this%n_max
                  do b = 0, this%l_max
                     i_cosnx = i_cosnx + 1

                     Ang = T_cos_ijk(b)

                     if(my_do_descriptor) &
                        descriptor_out%x(i_desc)%data(i_cosnx) = descriptor_out%x(i_desc)%data(i_cosnx) + Rad_ij(a)*Rad_ik(a)*Ang*0.5_dp

                     if(my_do_grad_descriptor) then

                        dAng_ij = b*U_cos_ijk(b-1) * dcosijk_ij
                        dAng_ik = b*U_cos_ijk(b-1) * dcosijk_ik

                        descriptor_out%x(i_desc)%grad_data(i_cosnx,:,0) = descriptor_out%x(i_desc)%grad_data(i_cosnx,:,0) - &
                        ( Rad_ij(a)*Rad_ik(a)*(dAng_ij+dAng_ik) + dRad_ij(:,a)*Rad_ik(a)*Ang + Rad_ij(a)*dRad_ik(:,a)*Ang ) * 0.5_dp

                        descriptor_out%x(i_desc)%grad_data(i_cosnx,:,n_i) = descriptor_out%x(i_desc)%grad_data(i_cosnx,:,n_i) + &
                        (Rad_ij(a)*Rad_ik(a)*dAng_ij + dRad_ij(:,a)*Rad_ik(a)*Ang)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

      if(allocated(Rad_ij)) deallocate(Rad_ij)
      if(allocated(Rad_ik)) deallocate(Rad_ik)
      if(allocated(T_cos_ijk)) deallocate(T_cos_ijk)
      if(allocated(U_cos_ijk)) deallocate(U_cos_ijk)
      if(allocated(dRad_ij)) deallocate(dRad_ij)
      if(allocated(dRad_ik)) deallocate(dRad_ik)

      call system_timer('cosnx_calc')

   endsubroutine cosnx_calc

   subroutine trihis_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(trihis), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, k, n, m, i_desc
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij, r_ik, r_jk, cos_ijk, Sym_Cor_S, Sym_Cor_A, exp_desc
      real(dp), dimension(3) :: u_ij, u_ik, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik, x, exp_arg, dexp_desc
      real(dp), dimension(3,3) :: dx_j, dx_k

      INIT_ERROR(error)

      call system_timer('trihis_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("trihis_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      d = trihis_dimensions(this,error)

      allocate(descriptor_out%x(at%N))
      do i = 1, at%N
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%has_data = .false.
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,0:atoms_n_neighbours(at,i)))
            allocate(descriptor_out%x(i)%ii(0:atoms_n_neighbours(at,i)))
            allocate(descriptor_out%x(i)%pos(3,0:atoms_n_neighbours(at,i)))
            allocate(descriptor_out%x(i)%has_grad_data(0:atoms_n_neighbours(at,i)))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.
         endif
      enddo

      do i = 1, at%N

         if(my_do_descriptor) descriptor_out%x(i)%has_data = .true.
         if(my_do_grad_descriptor) then
            descriptor_out%x(i)%ii(0) = i
            descriptor_out%x(i)%pos(:,0) = at%pos(:,i) 
            descriptor_out%x(i)%has_grad_data(0) = .true.
         endif

         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij > this%cutoff ) cycle

            if(my_do_grad_descriptor) then
               descriptor_out%x(i)%ii(n) = j
               descriptor_out%x(i)%pos(:,n) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i)%has_grad_data(n) = .true.
            endif

            do m = 1, atoms_n_neighbours(at,i)
               k = atoms_neighbour(at, i, m, distance = r_ik, cosines=u_ik, diff=d_ik)
               if( r_ik > this%cutoff ) cycle

               d_jk = d_ik - d_ij
               r_jk = norm(d_jk)
               if( r_jk .feq. 0.0_dp ) cycle

               cos_ijk = dot_product(u_ij,u_ik)
               Sym_Cor_S = r_ij + r_ik
               Sym_Cor_A = (r_ij - r_ik)**2

               x = (/Sym_Cor_S, Sym_Cor_A, cos_ijk/)

               if(my_do_grad_descriptor) then
                  dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                  dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik

                  dx_j(:,1) = u_ij
                  dx_j(:,2) = 2.0_dp*(r_ij - r_ik)*u_ij
                  dx_j(:,3) = dcosijk_ij

                  dx_k(:,1) = u_ik
                  dx_k(:,2) = -2.0_dp*(r_ij - r_ik)*u_ik
                  dx_k(:,3) = dcosijk_ik
               endif

               do i_desc = 1, this%n_gauss

                  exp_arg = (x - this%gauss_centre(:,i_desc))/this%gauss_width(:,i_desc)
                  exp_desc = exp(-0.5_dp*sum(exp_arg**2))

                  if(my_do_descriptor) &
                     descriptor_out%x(i)%data(i_desc) = descriptor_out%x(i)%data(i_desc) + exp_desc

                  if(my_do_grad_descriptor) then
                     dexp_desc = -exp_desc * exp_arg / this%gauss_width(:,i_desc)

                     descriptor_out%x(i)%grad_data(i_desc,:,0) = descriptor_out%x(i)%grad_data(i_desc,:,0) - &
                        matmul(dx_j+dx_k,dexp_desc)
                     descriptor_out%x(i)%grad_data(i_desc,:,n) = descriptor_out%x(i)%grad_data(i_desc,:,n) + &
                        2.0_dp*matmul(dx_j,dexp_desc)
                  endif
               enddo
            enddo
         enddo
      enddo

      call system_timer('trihis_calc')

   endsubroutine trihis_calc

   subroutine water_monomer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(water_monomer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, i, iO, iH1, iH2
      integer, dimension(3) :: shift_1, shift_2
      integer, dimension(:,:), allocatable :: water_monomer_index
      real(dp) :: r1, r2
      real(dp), dimension(3) :: v1, v2, u1, u2

      INIT_ERROR(error)

      call system_timer('water_monomer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("water_monomer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      d = water_monomer_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%has_data = .false.
            descriptor_out%x(i)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,3))
            allocate(descriptor_out%x(i)%ii(3))
            allocate(descriptor_out%x(i)%pos(3,3))
            allocate(descriptor_out%x(i)%has_grad_data(3))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,3))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(water_monomer_index(3,n_descriptors))
      call find_water_monomer(at,water_monomer_index,error=error)

      do i = 1, n_descriptors

         iO = water_monomer_index(1,i)
         iH1 = water_monomer_index(2,i)
         iH2 = water_monomer_index(3,i)


         v1 = diff_min_image(at,iO,iH1,shift=shift_1)
         v2 = diff_min_image(at,iO,iH2,shift=shift_2)
         r1 = sqrt(dot_product(v1,v1))
         r2 = sqrt(dot_product(v2,v2))
         u1 = v1 / r1
         u2 = v2 / r2

         if(my_do_descriptor) then
            descriptor_out%x(i)%has_data = .true.
            descriptor_out%x(i)%data(1) = r1+r2 
            descriptor_out%x(i)%data(2) = (r1-r2)**2
            descriptor_out%x(i)%data(3) = dot_product(v1,v2)
         endif

         if(my_do_grad_descriptor) then
            descriptor_out%x(i)%ii(:) = water_monomer_index(:,i)
            descriptor_out%x(i)%pos(:,1) = at%pos(:,iO)
            descriptor_out%x(i)%pos(:,2) = at%pos(:,iH1) + matmul(at%lattice,shift_1)
            descriptor_out%x(i)%pos(:,3) = at%pos(:,iH2) + matmul(at%lattice,shift_2)
            descriptor_out%x(i)%has_grad_data(:) = .true.

            descriptor_out%x(i)%grad_data(1,:,1) = -u1-u2                  ! 1st descriptor wrt rO
            descriptor_out%x(i)%grad_data(1,:,2) =  u1                     ! 1st descriptor wrt rH1
            descriptor_out%x(i)%grad_data(1,:,3) =  u2                     ! 1st descriptor wrt rH2
            descriptor_out%x(i)%grad_data(2,:,1) =  2.0_dp*(r1-r2)*(u2-u1) ! 2nd descriptor wrt rO
            descriptor_out%x(i)%grad_data(2,:,2) =  2.0_dp*(r1-r2)*u1      ! 2nd descriptor wrt rH1
            descriptor_out%x(i)%grad_data(2,:,3) = -2.0_dp*(r1-r2)*u2      ! 2nd descriptor wrt rH2
            descriptor_out%x(i)%grad_data(3,:,1) =  -v1-v2                 ! 3rd descriptor wrt rO
            descriptor_out%x(i)%grad_data(3,:,2) =  v2                     ! 3rd descriptor wrt rH1
            descriptor_out%x(i)%grad_data(3,:,3) =  v1                     ! 3rd descriptor wrt rH2
         endif

      enddo

      deallocate(water_monomer_index)
      call system_timer('water_monomer_calc')

   endsubroutine water_monomer_calc

   subroutine water_dimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(water_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, n_monomers, i_desc, i, j, n, &
         iAO, iAH1, iAH2, iBO, iBH1, iBH2
      integer, dimension(3) :: shift_AO_BO, shift_AO_AH1, shift_AO_AH2, shift_AO_BH1, shift_AO_BH2, &
         shift_BO_AH1, shift_BO_AH2, shift_BO_BH1, shift_BO_BH2, &
         shift_AH1_AH2, shift_AH1_BH1, shift_AH1_BH2, shift_AH2_BH1, shift_AH2_BH2, shift_BH1_BH2
      real(dp), dimension(3) :: diff_AO_BO, diff_AO_AH1, diff_AO_AH2, diff_AO_BH1, diff_AO_BH2, &
         diff_BO_AH1, diff_BO_AH2, diff_BO_BH1, diff_BO_BH2, &
         diff_AH1_AH2, diff_AH1_BH1, diff_AH1_BH2, diff_AH2_BH1, diff_AH2_BH2, diff_BH1_BH2
      integer, dimension(:,:), allocatable :: water_monomer_index
      real(dp) :: r_AO_BO, r_AO_AH1, r_AO_AH2, r_AO_BH1, r_AO_BH2, r_BO_AH1, r_BO_AH2, r_BO_BH1, r_BO_BH2, &
         r_AH1_AH2, r_AH1_BH1, r_AH1_BH2, r_AH2_BH1, r_AH2_BH2, r_BH1_BH2
      integer, dimension(1) :: j_array

      INIT_ERROR(error)

      call system_timer('water_dimer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("water_dimer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      d = water_dimer_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%has_data = .false.
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,6))
            allocate(descriptor_out%x(i)%ii(6))
            allocate(descriptor_out%x(i)%pos(3,6))
            allocate(descriptor_out%x(i)%has_grad_data(6))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,6))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp

         endif
      enddo

      n_monomers = 0
      do i = 1, at%N
         if(at%Z(i) == 8) n_monomers = n_monomers+1
      enddo

      allocate(water_monomer_index(3,n_monomers))
      call find_water_monomer(at,water_monomer_index,OHH_ordercheck=this%OHH_ordercheck,error=error)

      i_desc = 0
      do i = 1, n_monomers
         iAO = water_monomer_index(1,i)
         iAH1 = water_monomer_index(2,i)
         iAH2 = water_monomer_index(3,i)

         diff_AO_AH1 = diff_min_image(at,iAO,iAH1,shift=shift_AO_AH1)
         diff_AO_AH2 = diff_min_image(at,iAO,iAH2,shift=shift_AO_AH2)
         diff_AH1_AH2 = diff_min_image(at,iAH1,iAH2,shift=shift_AH1_AH2)

         r_AO_AH1 = norm(diff_AO_AH1)
         r_AO_AH2 = norm(diff_AO_AH2)
         r_AH1_AH2 = norm(diff_AH1_AH2)

         do n = 1, atoms_n_neighbours(at,iAO)
            iBO = atoms_neighbour(at,iAO,n,distance=r_AO_BO, diff=diff_AO_BO, shift=shift_AO_BO )
            if(at%Z(iBO) /= 8) cycle
            if( r_AO_BO >= this%cutoff ) cycle
            i_desc = i_desc + 1

            j_array = find(water_monomer_index(1,:) == iBO)
            j = j_array(1)

            iBH1 = water_monomer_index(2,j)
            iBH2 = water_monomer_index(3,j)

            diff_BO_BH1 = diff_min_image(at,iBO,iBH1,shift=shift_BO_BH1)
            diff_BO_BH2 = diff_min_image(at,iBO,iBH2,shift=shift_BO_BH2)
            diff_BH1_BH2 = diff_min_image(at,iBH1,iBH2,shift=shift_BH1_BH2)

            r_BO_BH1 = norm(diff_BO_BH1)
            r_BO_BH2 = norm(diff_BO_BH2)
            r_BH1_BH2 = norm(diff_BH1_BH2)

            diff_AO_BH1 = diff_AO_BO + diff_BO_BH1
            diff_AO_BH2 = diff_AO_BO + diff_BO_BH2
            shift_AO_BH1 = shift_AO_BO + shift_BO_BH1
            shift_AO_BH2 = shift_AO_BO + shift_BO_BH2

            r_AO_BH1 = norm(diff_AO_BH1)
            r_AO_BH2 = norm(diff_AO_BH2)

            diff_BO_AH1 = -diff_AO_BO + diff_AO_AH1
            diff_BO_AH2 = -diff_AO_BO + diff_AO_AH2

            shift_BO_AH1 = -shift_AO_BO + shift_AO_AH1
            shift_BO_AH2 = -shift_AO_BO + shift_AO_AH2

            r_BO_AH1 = norm(diff_BO_AH1)
            r_BO_AH2 = norm(diff_BO_AH2)

            diff_AH1_BH1 = -diff_AO_AH1 + diff_AO_BO + diff_BO_BH1
            diff_AH1_BH2 = -diff_AO_AH1 + diff_AO_BO + diff_BO_BH2
            diff_AH2_BH1 = -diff_AO_AH2 + diff_AO_BO + diff_BO_BH1
            diff_AH2_BH2 = -diff_AO_AH2 + diff_AO_BO + diff_BO_BH2

            shift_AH1_BH1 = -shift_AO_AH1 + shift_AO_BO + shift_BO_BH1
            shift_AH1_BH2 = -shift_AO_AH1 + shift_AO_BO + shift_BO_BH2
            shift_AH2_BH1 = -shift_AO_AH2 + shift_AO_BO + shift_BO_BH1
            shift_AH2_BH2 = -shift_AO_AH2 + shift_AO_BO + shift_BO_BH2

            r_AH1_BH1 = norm(diff_AH1_BH1)
            r_AH1_BH2 = norm(diff_AH1_BH2)
            r_AH2_BH1 = norm(diff_AH2_BH1)
            r_AH2_BH2 = norm(diff_AH2_BH2)

            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%has_data = .true.
               descriptor_out%x(i_desc)%data(:) = (/r_AO_BO, &
                  r_AO_AH1, r_AO_AH2, r_AO_BH1, r_AO_BH2, r_BO_AH1, r_BO_AH2, r_BO_BH1, r_BO_BH2, &
                  r_AH1_AH2, r_AH1_BH1, r_AH1_BH2, r_AH2_BH1, r_AH2_BH2, r_BH1_BH2/)

               descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(r_AO_BO, &
               this%cutoff,this%cutoff_transition_width)
            endif

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(:) = (/ water_monomer_index(:,i),water_monomer_index(:,j) /)
               descriptor_out%x(i_desc)%pos(:,1) = at%pos(:,iAO) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,2) = at%pos(:,iAH1) + matmul(at%lattice,shift_AO_AH1) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,3) = at%pos(:,iAH2) + matmul(at%lattice,shift_AO_AH2) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,4) = at%pos(:,iBO) + matmul(at%lattice,shift_AO_BO) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,5) = at%pos(:,iBH1) + matmul(at%lattice,shift_AO_BH1) ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%pos(:,6) = at%pos(:,iBH2) + matmul(at%lattice,shift_AO_BH2) ! TODO: Have to figure out how to do this.

               descriptor_out%x(i_desc)%has_grad_data(:) = .true.

               descriptor_out%x(i_desc)%grad_data(1,:,1) = -diff_AO_BO / r_AO_BO     ! 1st descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(1,:,4) = -descriptor_out%x(i_desc)%grad_data(1,:,1)        ! 1st descriptor wrt OB

               descriptor_out%x(i_desc)%grad_data(2,:,1) = -diff_AO_AH1 / r_AO_AH1  ! 2nd descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(2,:,2) = -descriptor_out%x(i_desc)%grad_data(2,:,1)        ! 2nd descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(3,:,1) = -diff_AO_AH2 / r_AO_AH2  ! 3rd descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(3,:,3) = -descriptor_out%x(i_desc)%grad_data(3,:,1)        ! 3rd descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(4,:,1) = -diff_AO_BH1 / r_AO_BH1  ! 4th descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(4,:,5) = -descriptor_out%x(i_desc)%grad_data(4,:,1)        ! 4th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(5,:,1) = -diff_AO_BH2 / r_AO_BH2  ! 5th descriptor wrt OA
               descriptor_out%x(i_desc)%grad_data(5,:,6) = -descriptor_out%x(i_desc)%grad_data(5,:,1)        ! 5th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_data(6,:,4) = -diff_BO_AH1 / r_BO_AH1  ! 6th descriptor wrt OB
               descriptor_out%x(i_desc)%grad_data(6,:,2) = -descriptor_out%x(i_desc)%grad_data(6,:,4)        ! 6th descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(7,:,4) = -diff_BO_AH2 / r_BO_AH2  ! 7th descriptor wrt OB
               descriptor_out%x(i_desc)%grad_data(7,:,3) = -descriptor_out%x(i_desc)%grad_data(7,:,4)        ! 7th descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(8,:,4) = -diff_BO_BH1 / r_BO_BH1  ! 8th descriptor wrt OB
               descriptor_out%x(i_desc)%grad_data(8,:,5) = -descriptor_out%x(i_desc)%grad_data(8,:,4)        ! 8th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(9,:,4) = -diff_BO_BH2 / r_BO_BH2  ! 9th descriptor wrt OB
               descriptor_out%x(i_desc)%grad_data(9,:,6) = -descriptor_out%x(i_desc)%grad_data(9,:,4)        ! 9th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_data(10,:,2) = -diff_AH1_AH2 / r_AH1_AH2 ! 10th descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(10,:,3) = -descriptor_out%x(i_desc)%grad_data(10,:,2)         ! 10th descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(11,:,2) = -diff_AH1_BH1 / r_AH1_BH1 ! 11th descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(11,:,5) = -descriptor_out%x(i_desc)%grad_data(11,:,2)         ! 11th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(12,:,2) = -diff_AH1_BH2 / r_AH1_BH2 ! 12th descriptor wrt AH1
               descriptor_out%x(i_desc)%grad_data(12,:,6) = -descriptor_out%x(i_desc)%grad_data(12,:,2)         ! 12th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_data(13,:,3) = -diff_AH2_BH1 / r_AH2_BH1 ! 13th descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(13,:,5) = -descriptor_out%x(i_desc)%grad_data(13,:,3)         ! 13th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(14,:,3) = -diff_AH2_BH2 / r_AH2_BH2 ! 14th descriptor wrt AH2
               descriptor_out%x(i_desc)%grad_data(14,:,6) = -descriptor_out%x(i_desc)%grad_data(14,:,3)         ! 14th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_data(15,:,5) = -diff_BH1_BH2 / r_BH1_BH2 ! 15th descriptor wrt BH1
               descriptor_out%x(i_desc)%grad_data(15,:,6) = -descriptor_out%x(i_desc)%grad_data(15,:,5)         ! 15th descriptor wrt BH2

               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = -dcoordination_function(r_AO_BO,&
               this%cutoff,this%cutoff_transition_width) * diff_AO_BO / r_AO_BO
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,4) = -descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1)
            endif
         enddo
      enddo

      deallocate(water_monomer_index)
      call system_timer('water_dimer_calc')

   endsubroutine water_dimer_calc

   subroutine A2_dimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(A2_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, n_monomers, i_desc, i, j, &
         iA1, iA2, iB1, iB2
      integer, dimension(3) :: shift_A1_A2, shift_A1_B1, shift_A1_B2, shift_A2_B1, shift_A2_B2, shift_B1_B2
      integer, dimension(at%N) :: A2_monomer_index
      real(dp) :: r_A1_A2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2, r_B1_B2

      INIT_ERROR(error)

      call system_timer('A2_dimer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      d = A2_dimer_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%has_data = .false.
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,4))
            allocate(descriptor_out%x(i)%ii(4))
            allocate(descriptor_out%x(i)%pos(3,4))
            allocate(descriptor_out%x(i)%has_grad_data(4))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,4))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      n_monomers = count(at%Z == this%atomic_number) / 2

      call find_A2_monomer(at,this%atomic_number, this%monomer_cutoff, A2_monomer_index,error)

      i_desc = 0
      do i = 1, at%N
         iA1 = i
         iA2 = atoms_neighbour(at,i,A2_monomer_index(i),distance=r_A1_A2,shift=shift_A1_A2)
         if( iA1 > iA2 ) cycle

         do j = i + 1, at%N
            iB1 = j
            iB2 = atoms_neighbour(at,j,A2_monomer_index(j),distance=r_B1_B2,shift=shift_B1_B2)
            if( iB1 > iB2 ) cycle

            r_A1_B1 = distance_min_image(at,iA1,iB1,shift=shift_A1_B1)
            r_A1_B2 = distance_min_image(at,iA1,iB2,shift=shift_A1_B2)

            r_A2_B1 = distance_min_image(at,iA2,iB1,shift=shift_A2_B1)
            r_A2_B2 = distance_min_image(at,iA2,iB2,shift=shift_A2_B2)
            
            if( any( (/r_A1_A2,r_B1_B2,r_A1_B1,r_A1_B2,r_A2_B1,r_A2_B2/) > this%cutoff) ) cycle
            i_desc = i_desc + 1

            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%has_data = .true.
               descriptor_out%x(i_desc)%data(:) = (/ r_A1_A2, r_B1_B2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2/)
            endif

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(:) = (/ iA1, iA2, iB1, iB2 /)
               descriptor_out%x(i_desc)%pos(:,:) = 0.0_dp ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%has_grad_data(:) = .true.

               descriptor_out%x(i_desc)%grad_data(1,:,1) = -diff(at,iA1,iA2,shift=shift_A1_A2) / r_A1_A2      ! 1st descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(1,:,2) = -descriptor_out%x(i_desc)%grad_data(1,:,1)         ! 1st descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(2,:,3) = -diff(at,iB1,iB2,shift=shift_B1_B2) / r_B1_B2      ! 2nd descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(2,:,4) = -descriptor_out%x(i_desc)%grad_data(2,:,3)         ! 2nd descriptor wrt B2

               descriptor_out%x(i_desc)%grad_data(3,:,1) = -diff(at,iA1,iB1,shift=shift_A1_B1) / r_A1_B1      ! 3rd descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(3,:,3) = -descriptor_out%x(i_desc)%grad_data(3,:,1)         ! 3rd descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(4,:,1) = -diff(at,iA1,iB2,shift=shift_A1_B2) / r_A1_B2      ! 4th descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(4,:,4) = -descriptor_out%x(i_desc)%grad_data(4,:,1)         ! 4th descriptor wrt B2

               descriptor_out%x(i_desc)%grad_data(5,:,2) = -diff(at,iA2,iB1,shift=shift_A2_B1) / r_A2_B1      ! 5th descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(5,:,3) = -descriptor_out%x(i_desc)%grad_data(5,:,2)         ! 5th descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(6,:,2) = -diff(at,iA2,iB2,shift=shift_A2_B2) / r_A2_B2      ! 6th descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(6,:,4) = -descriptor_out%x(i_desc)%grad_data(6,:,2)         ! 6th descriptor wrt B2

            endif
         enddo
      enddo

      call system_timer('A2_dimer_calc')

   endsubroutine A2_dimer_calc

   subroutine AB_dimer_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(AB_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, n_descriptors, n_cross, n_monomers, i_desc, i, j, &
         iA1, iA2, iB1, iB2
      integer, dimension(3) :: shift_A1_A2, shift_A1_B1, shift_A1_B2, shift_A2_B1, shift_A2_B2, shift_B1_B2
      integer, dimension(:,:), allocatable :: AB_monomer_index
      real(dp) :: r_A1_A2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2, r_B1_B2

      INIT_ERROR(error)

      call system_timer('AB_dimer_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("AB_dimer_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      d = AB_dimer_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))
      do i = 1, n_descriptors
         if(my_do_descriptor) then
            allocate(descriptor_out%x(i)%data(d))
            descriptor_out%x(i)%data = 0.0_dp
            descriptor_out%x(i)%has_data = .false.
         endif
         if(my_do_grad_descriptor) then
            allocate(descriptor_out%x(i)%grad_data(d,3,4))
            allocate(descriptor_out%x(i)%ii(4))
            allocate(descriptor_out%x(i)%pos(3,4))
            allocate(descriptor_out%x(i)%has_grad_data(4))
            descriptor_out%x(i)%grad_data = 0.0_dp
            descriptor_out%x(i)%ii = 0
            descriptor_out%x(i)%pos = 0.0_dp
            descriptor_out%x(i)%has_grad_data = .false.

            allocate(descriptor_out%x(i)%grad_covariance_cutoff(3,4))
            descriptor_out%x(i)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      if( count(at%Z == this%atomic_number1) == count(at%Z == this%atomic_number2) ) then
         n_monomers = count(at%Z == this%atomic_number1)
      else
         RAISE_ERROR("AB_dimer_calc: number of monomer atoms 1 ("//count(at%Z == this%atomic_number1)//") not equal to number of monomer atoms 2 ("//count(at%Z == this%atomic_number1)//")",error)
      endif

      allocate(AB_monomer_index(2,n_monomers))
      call find_AB_monomer(at,(/this%atomic_number1,this%atomic_number2/), this%monomer_cutoff, AB_monomer_index,error)

      i_desc = 0
      do i = 1, n_monomers
         iA1 = AB_monomer_index(1,i)
         iB1 = AB_monomer_index(2,i)
         do j = i + 1, n_monomers
            iA2 = AB_monomer_index(1,j)
            iB2 = AB_monomer_index(2,j)


            r_A1_B1 = distance_min_image(at,iA1,iB1,shift=shift_A1_B1)
            r_A2_B2 = distance_min_image(at,iA2,iB2,shift=shift_A2_B2)

            r_A1_A2 = distance_min_image(at,iA1,iA2,shift=shift_A1_A2)
            r_B1_B2 = distance_min_image(at,iB1,iB2,shift=shift_B1_B2)

            r_A1_B2 = distance_min_image(at,iA1,iB2,shift=shift_A1_B2)
            r_A2_B1 = distance_min_image(at,iA2,iB1,shift=shift_A2_B1)
            
            if( any( (/r_A1_A2,r_B1_B2,r_A1_B1,r_A1_B2,r_A2_B1,r_A2_B2/) > this%cutoff) ) cycle
            i_desc = i_desc + 1

            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%has_data = .true.
               descriptor_out%x(i_desc)%data(:) = (/ r_A1_B1, r_A2_B2, r_A1_A2, r_B1_B2, r_A1_B2, r_A2_B1 /)
            endif

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(:) = (/ AB_monomer_index(:,i),AB_monomer_index(:,j) /)
               descriptor_out%x(i_desc)%pos(:,:) = 0.0_dp ! TODO: Have to figure out how to do this.
               descriptor_out%x(i_desc)%has_grad_data(:) = .true.

               descriptor_out%x(i_desc)%grad_data(1,:,1) = -diff(at,iA1,iB1,shift=shift_A1_B1) / r_A1_B1      ! 1st descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(1,:,2) = -descriptor_out%x(i_desc)%grad_data(1,:,1)         ! 1st descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(2,:,3) = -diff(at,iA2,iB2,shift=shift_A2_B2) / r_A2_B2      ! 2nd descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(2,:,4) = -descriptor_out%x(i_desc)%grad_data(2,:,3)         ! 2nd descriptor wrt B2

               descriptor_out%x(i_desc)%grad_data(3,:,1) = -diff(at,iA1,iA2,shift=shift_A1_A2) / r_A1_A2      ! 1st descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(3,:,3) = -descriptor_out%x(i_desc)%grad_data(3,:,1)         ! 1st descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(4,:,2) = -diff(at,iB1,iB2,shift=shift_B1_B2) / r_B1_B2      ! 2nd descriptor wrt B1
               descriptor_out%x(i_desc)%grad_data(4,:,4) = -descriptor_out%x(i_desc)%grad_data(4,:,2)         ! 2nd descriptor wrt B2

               descriptor_out%x(i_desc)%grad_data(5,:,1) = -diff(at,iA1,iB2,shift=shift_A1_B2) / r_A1_B2      ! 4th descriptor wrt A1
               descriptor_out%x(i_desc)%grad_data(5,:,4) = -descriptor_out%x(i_desc)%grad_data(5,:,1)         ! 4th descriptor wrt B2
               descriptor_out%x(i_desc)%grad_data(6,:,3) = -diff(at,iA2,iB1,shift=shift_A2_B1) / r_A2_B1      ! 5th descriptor wrt A2
               descriptor_out%x(i_desc)%grad_data(6,:,2) = -descriptor_out%x(i_desc)%grad_data(6,:,3)         ! 5th descriptor wrt B1

            endif
         enddo
      enddo

      deallocate(AB_monomer_index)
      call system_timer('AB_dimer_calc')

   endsubroutine AB_dimer_calc

   subroutine bond_real_space_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(bond_real_space), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      type(atoms) :: at_copy
      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: n_descriptors, n_cross, i_desc, i, j, n, k, m, m_index, l, ij_neighbours
      integer, dimension(3) :: shift_j, shift_k
      real(dp) :: r_ij, r_ijk
      real(dp) :: atom_i(3), atom_j(3), atom_k(3), bond(3), bond_len
      real(dp) :: atom_i_cross_atom_j(3), atom_i_normsq_min_atom_j_normsq
      real(dp), allocatable :: r(:,:), z(:), c(:)
      real(dp), allocatable :: dr(:,:,:,:), dz(:,:,:), dc(:,:,:)
      integer, allocatable :: ii(:)
      real(dp), allocatable :: pos(:,:)
      real(dp) :: r_m_cross_r_l(3)

      INIT_ERROR(error)

      call system_timer('bond_real_space_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("bond_real_space_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0

      do i = 1, at%N
         do n = 1, atoms_n_neighbours(at, i)
            j = atoms_neighbour(at, i, n, shift=shift_j, distance=r_ij, max_dist=this%bond_cutoff)

            if(j == 0) cycle

            i_desc = i_desc + 1

            atom_i = at%pos(:,i)
            atom_j = at%pos(:,j) + matmul(at%lattice, shift_j)

            at_copy = at
            call add_atoms(at_copy, 0.5_dp * (atom_i + atom_j), 1)
            call calc_connect(at_copy)

            ij_neighbours = 0

            do m = 1, atoms_n_neighbours(at_copy, at%N + 1)
               k = atoms_neighbour(at_copy, at%N + 1, m, max_dist=this%cutoff)

               if(k == 0) cycle

               if(at_copy%pos(:,k) .feq. at_copy%pos(:,at%N + 1)) cycle

               ij_neighbours = ij_neighbours + 1
            enddo

            if(my_do_descriptor .or. my_do_grad_descriptor) then
               allocate(r(3,ij_neighbours), z(ij_neighbours), c(ij_neighbours))

               r = 0.0_dp
               z = 0.0_dp
               c = 0.0_dp
               bond = atom_i - atom_j
               bond_len = norm(bond)
               atom_i_cross_atom_j = atom_i .cross. atom_j
               atom_i_normsq_min_atom_j_normsq = normsq(atom_i) - normsq(atom_j)

               if(my_do_grad_descriptor) then
                  allocate(dr(3,ij_neighbours,3,ij_neighbours), dz(ij_neighbours,3,ij_neighbours), dc(ij_neighbours,3,ij_neighbours))
                  allocate(ii(ij_neighbours), pos(3,ij_neighbours))

                  dr = 0.0_dp
                  dz = 0.0_dp
                  dc = 0.0_dp
                  ii = 0
                  pos = 0.0_dp
               endif

               m_index = 2

               do m = 1, atoms_n_neighbours(at_copy, at%N + 1)
                  k = atoms_neighbour(at_copy, at%N + 1, m, shift=shift_k, distance=r_ijk, max_dist=this%cutoff)

                  if(k == 0) cycle

                  if(at_copy%pos(:,k) .feq. at_copy%pos(:,at%N + 1)) cycle

                  atom_k = at_copy%pos(:,k) + matmul(at_copy%lattice, shift_k)

                  if(atom_k .feq. atom_i) then
                     ! r remains zero
                     z(1) = 0.5_dp * bond_len
                     c(1) = coordination_function(r_ijk, this%cutoff, this%transition_width)

                     if(my_do_grad_descriptor) then
                        ! dr remains zero
                        dz(1,:,1) = 0.5_dp * bond / bond_len
                        dz(1,:,2) = - dz(1,:,1)
                        dc(1,:,1) = 0.25_dp * dcoordination_function(r_ijk, this%cutoff, this%transition_width) * bond / r_ijk
                        dc(1,:,2) = - dc(1,:,1)

                        ii(1) = k
                        pos(:,1) = atom_k
                     endif
                  elseif(atom_k .feq. atom_j) then
                     ! r remain zero
                     z(2) = -0.5_dp * bond_len
                     c(2) = coordination_function(r_ijk, this%cutoff, this%transition_width)

                     if(my_do_grad_descriptor) then
                        ! dr remains zero
                        dz(2,:,1) = -0.5_dp * bond / bond_len
                        dz(2,:,2) = - dz(2,:,1)
                        dc(2,:,1) = -0.25_dp * dcoordination_function(r_ijk, this%cutoff, this%transition_width) * bond / r_ijk
                        dc(2,:,2) = - dc(2,:,1)

                        ii(2) = k
                        pos(:,2) = atom_k
                     endif
                  else
                     m_index = m_index + 1

                     r(:,m_index) = ((atom_k .cross. bond) + atom_i_cross_atom_j) / bond_len
                     z(m_index) = ((atom_k .dot. bond) - 0.5_dp * atom_i_normsq_min_atom_j_normsq) / bond_len
                     c(m_index) = coordination_function(r_ijk, this%cutoff, this%transition_width)

                     if(my_do_grad_descriptor) then
                        dr(:,m_index,1,1) = ((/ 0.0_dp, atom_k(3) - atom_j(3), atom_j(2) - atom_k(2) /) / bond_len) - (r(:,m_index) * bond(1) / bond_len**2)
                        dr(:,m_index,2,1) = ((/ atom_j(3) - atom_k(3), 0.0_dp, atom_k(1) - atom_j(1) /) / bond_len) - (r(:,m_index) * bond(2) / bond_len**2)
                        dr(:,m_index,3,1) = ((/ atom_k(2) - atom_j(2), atom_j(1) - atom_k(1), 0.0_dp /) / bond_len) - (r(:,m_index) * bond(3) / bond_len**2)
                        dz(m_index,:,1) = ((atom_k - atom_i) / bond_len) - (z(m_index) * bond / bond_len**2)
                        dc(m_index,:,1) = -0.5_dp * dcoordination_function(r_ijk, this%cutoff, this%transition_width) * (atom_k - at_copy%pos(:,at%N + 1)) / r_ijk

                        dr(:,m_index,1,2) = - dr(:,m_index,1,1) + ((/ 0.0_dp, bond(3), - bond(2) /) / bond_len)
                        dr(:,m_index,2,2) = - dr(:,m_index,2,1) + ((/ - bond(3), 0.0_dp, bond(1) /) / bond_len)
                        dr(:,m_index,3,2) = - dr(:,m_index,3,1) + ((/ bond(2), - bond(1), 0.0_dp /) / bond_len)
                        dz(m_index,:,2) = - dz(m_index,:,1) - (bond / bond_len)
                        dc(m_index,:,2) = dc(m_index,:,1)

                        dr(:,m_index,1,m_index) = (/ 0.0_dp, - bond(3), bond(2) /) / bond_len
                        dr(:,m_index,2,m_index) = (/ bond(3), 0.0_dp, - bond(1) /) / bond_len
                        dr(:,m_index,3,m_index) = (/ - bond(2), bond(1), 0.0_dp /) / bond_len
                        dz(m_index,:,m_index) = bond / bond_len
                        dc(m_index,:,m_index) = -2.0_dp * dc(m_index,:,1)

                        ii(m_index) = k
                        pos(:,m_index) = atom_k
                     endif
                  endif
               enddo
            endif

            if(my_do_descriptor) then
               allocate(descriptor_out%x(i_desc)%data(1 + (1 + 2 * ij_neighbours) * ij_neighbours))

               !descriptor_out%x(i_desc)%data = 0.0_dp

               descriptor_out%x(i_desc)%data(1) = real(ij_neighbours, dp)

               descriptor_out%x(i_desc)%data(2:ij_neighbours + 1) = c

               do m = 1, ij_neighbours
                  if(m <= 2) then
                     descriptor_out%x(i_desc)%data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m)) = z(m)
                  else
                     do l = 3, ij_neighbours
                        if(l == m) then
                           descriptor_out%x(i_desc)%data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1) = normsq(r(:,m))
                           descriptor_out%x(i_desc)%data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l)) = z(m)
                        else
                           descriptor_out%x(i_desc)%data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1) = r(:,m) .dot. r(:,l)
                           descriptor_out%x(i_desc)%data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l)) = ((r(:,m) .cross. r(:,l)) .dot. bond) / bond_len
                        endif
                     enddo
                  endif
               enddo

               descriptor_out%x(i_desc)%covariance_cutoff = coordination_function(r_ij, this%bond_cutoff, this%bond_transition_width)

               descriptor_out%x(i_desc)%has_data = .true.
            endif

            if(my_do_grad_descriptor) then
               allocate(descriptor_out%x(i_desc)%grad_data(1 + (1 + 2 * ij_neighbours) * ij_neighbours,3,ij_neighbours))
               allocate(descriptor_out%x(i_desc)%ii(ij_neighbours))
               allocate(descriptor_out%x(i_desc)%pos(3,ij_neighbours))
               allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,ij_neighbours))
               allocate(descriptor_out%x(i_desc)%has_grad_data(ij_neighbours))

               descriptor_out%x(i_desc)%grad_data = 0.0_dp

               !descriptor_out%x(i_desc)%grad_data(1,:,:) = 0.0_dp

               descriptor_out%x(i_desc)%grad_data(2:ij_neighbours + 1,:,:) = dc

               do m = 1, ij_neighbours
                  if(m <= 2) then
                     descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m),:,1) = dz(m,:,1)
                     descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * m),:,2) = dz(m,:,2)
                  else
                     do l = 3, ij_neighbours
                        if(l == m) then
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,1) = 2.0_dp * (r(:,m) .dot. dr(:,m,1,1))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,1) = 2.0_dp * (r(:,m) .dot. dr(:,m,2,1))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,1) = 2.0_dp * (r(:,m) .dot. dr(:,m,3,1))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),:,1) = dz(m,:,1)

                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,2) = 2.0_dp * (r(:,m) .dot. dr(:,m,1,2))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,2) = 2.0_dp * (r(:,m) .dot. dr(:,m,2,2))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,2) = 2.0_dp * (r(:,m) .dot. dr(:,m,3,2))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),:,2) = dz(m,:,2)

                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,m) = 2.0_dp * (r(:,m) .dot. dr(:,m,1,m))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,m) = 2.0_dp * (r(:,m) .dot. dr(:,m,2,m))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,m) = 2.0_dp * (r(:,m) .dot. dr(:,m,3,m))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),:,m) = dz(m,:,m)
                        else
                           r_m_cross_r_l = r(:,m) .cross. r(:,l)

                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,1) = (dr(:,m,1,1) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,1,1))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,1) = (dr(:,m,2,1) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,2,1))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,1) = (dr(:,m,3,1) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,3,1))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),1,1) = ((((dr(:,m,1,1) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,1,1))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 1.0_dp, 0.0_dp, 0.0_dp /) - (bond * bond(1) / bond_len**2)))) / bond_len
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),2,1) = ((((dr(:,m,2,1) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,2,1))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 0.0_dp, 1.0_dp, 0.0_dp /) - (bond * bond(2) / bond_len**2)))) / bond_len
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),3,1) = ((((dr(:,m,3,1) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,3,1))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 0.0_dp, 0.0_dp, 1.0_dp /) - (bond * bond(3) / bond_len**2)))) / bond_len

                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,2) = (dr(:,m,1,2) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,1,2))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,2) = (dr(:,m,2,2) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,2,2))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,2) = (dr(:,m,3,2) .dot. r(:,l)) + (r(:,m) .dot. dr(:,l,3,2))
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),1,2) = ((((dr(:,m,1,2) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,1,2))) .dot. bond) + (r_m_cross_r_l .dot. ((/ -1.0_dp, 0.0_dp, 0.0_dp /) + (bond * bond(1) / bond_len**2)))) / bond_len 
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),2,2) = ((((dr(:,m,2,2) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,2,2))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 0.0_dp, -1.0_dp, 0.0_dp /) + (bond * bond(2) / bond_len**2)))) / bond_len
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),3,2) = ((((dr(:,m,3,2) .cross. r(:,l)) + (r(:,m) .cross. dr(:,l,3,2))) .dot. bond) + (r_m_cross_r_l .dot. ((/ 0.0_dp, 0.0_dp, -1.0_dp /) + (bond * bond(3) / bond_len**2)))) / bond_len

                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,m) = dr(:,m,1,m) .dot. r(:,l)
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,m) = dr(:,m,2,m) .dot. r(:,l)
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,m) = dr(:,m,3,m) .dot. r(:,l)
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),1,m) = ((dr(:,m,1,m) .cross. r(:,l)) .dot. bond) / bond_len
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),2,m) = ((dr(:,m,2,m) .cross. r(:,l)) .dot. bond) / bond_len
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),3,m) = ((dr(:,m,3,m) .cross. r(:,l)) .dot. bond) / bond_len

                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,1,l) = r(:,m) .dot. dr(:,l,1,l)
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,2,l) = r(:,m) .dot. dr(:,l,2,l)
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l) - 1,3,l) = r(:,m) .dot. dr(:,l,3,l)
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),1,l) = ((r(:,m) .cross. dr(:,l,1,l)) .dot. bond) / bond_len
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),2,l) = ((r(:,m) .cross. dr(:,l,2,l)) .dot. bond) / bond_len
                           descriptor_out%x(i_desc)%grad_data(1 + ij_neighbours + (2 * (m - 1) * ij_neighbours) + (2 * l),3,l) = ((r(:,m) .cross. dr(:,l,3,l)) .dot. bond) / bond_len
                        endif
                     enddo
                  endif
               enddo

               descriptor_out%x(i_desc)%ii = ii
               descriptor_out%x(i_desc)%pos = pos

               descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp

               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1) = dcoordination_function(r_ij, this%bond_cutoff, this%bond_transition_width) * bond / r_ij
               descriptor_out%x(i_desc)%grad_covariance_cutoff(:,2) = - descriptor_out%x(i_desc)%grad_covariance_cutoff(:,1)

               descriptor_out%x(i_desc)%has_grad_data = .true.
            endif

            if(my_do_descriptor .or. my_do_grad_descriptor) then
               deallocate(r, z, c)

               if(my_do_grad_descriptor) then
                  deallocate(dr, dz, dc)
                  deallocate(ii, pos)
               endif
            endif

            call finalise(at_copy)
         enddo
      enddo

      call system_timer('bond_real_space_calc')

   endsubroutine bond_real_space_calc

   subroutine atom_real_space_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(atom_real_space), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, grad_d, n_descriptors, n_cross, descriptor_mould_size, i_desc, i_data, i, j, k, n, l, m, n_neighbours, i_n

      real(dp) :: r
      real(dp), dimension(3) :: diff
      real(dp), dimension(1) :: descriptor_mould
      integer, dimension(3) :: shift

      complex(dp), dimension(:), allocatable :: spherical_harmonics
      complex(dp), dimension(:,:), allocatable :: grad_spherical_harmonics

      INIT_ERROR(error)

      call system_timer('atom_real_space_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("bond_real_space_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))
      descriptor_out%atom_based = .true.

      i_desc = 0
      do i = 1, at%N
         i_desc = i_desc + 1

         n_neighbours = atoms_n_neighbours(at,i,max_dist=this%cutoff)
         d = ( 2 * (this%l_max+1)**2 + 2 ) * n_neighbours

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            grad_d = 2 * (this%l_max+1)**2 + 2

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,1:n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(1:n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,1:n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(1:n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,1:n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(spherical_harmonics(-this%l_max:this%l_max))
      if( my_do_grad_descriptor ) allocate(grad_spherical_harmonics(3,-this%l_max:this%l_max))

      i_desc = 0
      do i = 1, at%N
         i_desc = i_desc + 1
         i_data = 0
         i_n = 0

         if(my_do_descriptor) then
            descriptor_out%x(i_desc)%has_data = .true.
         endif

         if(my_do_grad_descriptor) then
            !descriptor_out%x(i_desc)%ii(0) = i
            !descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            !descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         do n = 1, atoms_n_neighbours(at,i)

            j = atoms_neighbour(at,i,n,distance = r, diff = diff, shift=shift)
            if(r > this%cutoff) cycle
            i_n = i_n + 1

            i_data = i_data + 1
            if(my_do_descriptor) then
               descriptor_out%x(i_desc)%data(i_data) = r
            endif
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(i_n) = j
               descriptor_out%x(i_desc)%pos(:,i_n) = at%pos(:,j) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(i_n) = .true.
               descriptor_out%x(i_desc)%grad_data(i_data,:,i_n) = diff / r
            endif

            i_data = i_data + 1
            if(my_do_descriptor) descriptor_out%x(i_desc)%data(i_data) = real(i_n,dp)
            if(my_do_grad_descriptor) descriptor_out%x(i_desc)%grad_data(i_data,:,i_n) = real(i_n,dp)

            do l = 0, this%l_max
               descriptor_mould_size = size(transfer(spherical_harmonics(-l:l),descriptor_mould))
               
               do m = -l, l
                  if(my_do_descriptor) spherical_harmonics(m) = SphericalYCartesian(l,m,diff)
                  if(my_do_grad_descriptor) grad_spherical_harmonics(:,m) = GradSphericalYCartesian(l,m,diff)
               enddo

               if(my_do_descriptor) then
                  descriptor_out%x(i_desc)%data(i_data+1:i_data+descriptor_mould_size) = transfer(spherical_harmonics(-l:l),descriptor_mould)
               endif

               if(my_do_grad_descriptor) then
                  do k = 1, 3
                     descriptor_out%x(i_desc)%grad_data(i_data+1:i_data+descriptor_mould_size,k,i_n) = &
                     transfer(grad_spherical_harmonics(k,-l:l),descriptor_mould)
                  enddo
               endif

               i_data = i_data + descriptor_mould_size

            enddo
         enddo
      enddo

      if(allocated(spherical_harmonics)) deallocate(spherical_harmonics)
      if(allocated(grad_spherical_harmonics)) deallocate(grad_spherical_harmonics)

      call system_timer('atom_real_space_calc')

   endsubroutine atom_real_space_calc

   subroutine power_so3_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(power_so3), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      type(cplx_1d), dimension(:), allocatable :: SphericalY_ij
      type(cplx_1d), dimension(:,:), allocatable :: fourier_so3

      type(cplx_2d), dimension(:), allocatable :: dSphericalY_ij
      type(cplx_2d), dimension(:,:,:), allocatable :: dfourier_so3

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, n, a, l, m, i_desc, i_pow, n_neighbours, n_i, n_descriptors, n_cross
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij
      real(dp), dimension(3) :: u_ij, d_ij
      real(dp), dimension(:), allocatable :: Rad_ij
      real(dp), dimension(:,:), allocatable :: dRad_ij
      integer, dimension(116) :: species_map

      INIT_ERROR(error)

      call system_timer('power_so3_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("power_so3_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call finalise(descriptor_out)

      d = power_so3_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            n_neighbours = atoms_n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      allocate(fourier_so3(0:this%l_max,this%n_max),SphericalY_ij(0:this%l_max),Rad_ij(this%n_max))
      do a = 1, this%n_max
         do l = 0, this%l_max
            allocate(fourier_so3(l,a)%m(-l:l))
            fourier_so3(l,a)%m(:) = CPLX_ZERO
         enddo
      enddo
      do l = 0, this%l_max
         allocate(SphericalY_ij(l)%m(-l:l))
      enddo

      if(my_do_grad_descriptor) then
         allocate( dRad_ij(3,this%n_max), dSphericalY_ij(0:this%l_max) )
         do l = 0, this%l_max
            allocate(dSphericalY_ij(l)%mm(3,-l:l))
         enddo
      endif

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         if(my_do_descriptor) descriptor_out%x(i_desc)%has_data = .true.
         do a = 1, this%n_max
            do l = 0, this%l_max
               fourier_so3(l,a)%m(:) = CPLX_ZERO
            enddo
         enddo

         if(my_do_grad_descriptor) then
            allocate( dfourier_so3(0:this%l_max,this%n_max,0:atoms_n_neighbours(at,i,max_dist=this%cutoff)) )
            do n = 0, atoms_n_neighbours(at,i,max_dist=this%cutoff)
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     allocate(dfourier_so3(l,a,n)%mm(3,-l:l))
                     dfourier_so3(l,a,n)%mm(:,:) = CPLX_ZERO
                  enddo
               enddo
            enddo
         endif

         n_i = 0
         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij > this%cutoff ) cycle

            n_i = n_i + 1
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif

            do a = 1, this%n_max
               Rad_ij(a) = RadialFunction(this%Radial, r_ij, a)
               if(my_do_grad_descriptor) dRad_ij(:,a) = GradRadialFunction(this%Radial, r_ij, a) * u_ij
            enddo

            do l = 0, this%l_max
               do m = -l, l
                  SphericalY_ij(l)%m(m) = SphericalYCartesian(l,m,d_ij)
                  if(my_do_grad_descriptor) dSphericalY_ij(l)%mm(:,m) = GradSphericalYCartesian(l,m,d_ij)
               enddo
            enddo

            do a = 1, this%n_max
               do l = 0, this%l_max
                  do m = -l, l
                     fourier_so3(l,a)%m(m) = fourier_so3(l,a)%m(m) + Rad_ij(a)*SphericalY_ij(l)%m(m)
                     if(my_do_grad_descriptor) then
                        dfourier_so3(l,a,n_i)%mm(:,m) = dfourier_so3(l,a,n_i)%mm(:,m) + &
                        dRad_ij(:,a) * SphericalY_ij(l)%m(m) + Rad_ij(a)*dSphericalY_ij(l)%mm(:,m)
                     endif
                  enddo
               enddo
            enddo

         enddo ! n

         if(my_do_descriptor) then
            i_pow = 0
            do a = 1, this%n_max
               do l = 0, this%l_max
                  i_pow = i_pow + 1

                  descriptor_out%x(i_desc)%data(i_pow) = dot_product(fourier_so3(l,a)%m,fourier_so3(l,a)%m)
               enddo
            enddo
         endif

         if(my_do_grad_descriptor) then
            do n = 1, atoms_n_neighbours(at,i,max_dist=this%cutoff)
               i_pow = 0
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     i_pow = i_pow + 1

                     descriptor_out%x(i_desc)%grad_data(i_pow,:,n) = 2.0_dp * matmul(conjg(dfourier_so3(l,a,n)%mm(:,:)),fourier_so3(l,a)%m(:))
                  enddo
               enddo
               descriptor_out%x(i_desc)%grad_data(:,:,0) = descriptor_out%x(i_desc)%grad_data(:,:,0) - descriptor_out%x(i_desc)%grad_data(:,:,n)
            enddo
         endif

         if(allocated(dfourier_so3)) then
            do n = lbound(dfourier_so3,3), ubound(dfourier_so3,3)
               do a = lbound(dfourier_so3,2), ubound(dfourier_so3,2)
                  do l = lbound(dfourier_so3,1), ubound(dfourier_so3,1)
                     deallocate(dfourier_so3(l,a,n)%mm)
                  enddo
               enddo
            enddo
            deallocate(dfourier_so3)
         endif

      enddo ! i

      if(allocated(Rad_ij)) deallocate(Rad_ij)
      if(allocated(dRad_ij)) deallocate(dRad_ij)

      if(allocated(fourier_so3)) then
         do a = lbound(fourier_so3,2), ubound(fourier_so3,2)
            do l = lbound(fourier_so3,1), ubound(fourier_so3,1)
               deallocate(fourier_so3(l,a)%m)
            enddo
         enddo
         deallocate(fourier_so3)
      endif

      if(allocated(SphericalY_ij)) then
         do l = lbound(SphericalY_ij,1), ubound(SphericalY_ij,1)
            deallocate(SphericalY_ij(l)%m)
         enddo
         deallocate(SphericalY_ij)
      endif

      if(allocated(dSphericalY_ij)) then
         do l = lbound(dSphericalY_ij,1), ubound(dSphericalY_ij,1)
            deallocate(dSphericalY_ij(l)%mm)
         enddo
         deallocate(dSphericalY_ij)
      endif

      call system_timer('power_so3_calc')

   endsubroutine power_so3_calc

   subroutine power_SO4_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(power_SO4), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      type(cplx_2d), dimension(:), allocatable :: U
      type(cplx_3d), dimension(:,:), allocatable :: dU

      real(dp), dimension(3) :: diff, u_ij
      real(dp) :: r
      integer :: i, n, n_i, ji, jn, k, j, i_desc, i_bisp, d, n_descriptors, n_cross, n_neighbours
      integer, dimension(3) :: shift
      integer, dimension(116) :: species_map
      logical :: my_do_descriptor, my_do_grad_descriptor

      INIT_ERROR(error)

      call system_timer('power_SO4_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("power_SO4_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      species_map = 0
      do i = 1, size(this%species_Z)
         if(this%species_Z(i) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i)) = i
         endif
      enddo

      call finalise(descriptor_out)

      d = power_SO4_dimensions(this,error)
      call descriptor_sizes(this,at,n_descriptors,n_cross,error)
      allocate(descriptor_out%x(n_descriptors))

      i_desc = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle

         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            descriptor_out%x(i_desc)%data = 0.0_dp
            descriptor_out%x(i_desc)%has_data = .true.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif

         if(my_do_grad_descriptor) then
            n_neighbours = atoms_n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:n_neighbours))
            descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif

      enddo

      i_desc = 0
      do i = 1, at%N

         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         i_desc = i_desc + 1

         if(my_do_descriptor) descriptor_out%x(i_desc)%has_data = .true.
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
         endif

         n_i = 0
         do n = 1, atoms_n_neighbours(at,i)
            ji = atoms_neighbour(at, i, n, jn=jn, distance=r, diff=diff, cosines=u_ij,shift=shift)
            if( r > this%cutoff ) cycle

            n_i = n_i + 1

            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = ji
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,ji) + matmul(at%lattice,shift)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif
         enddo

         if(my_do_grad_descriptor) then
            call fourier_SO4_calc(this%fourier_SO4,at,i,U,dU,error=error)
         else
            call fourier_SO4_calc(this%fourier_SO4,at,i,U,error=error)
         endif

         if(my_do_descriptor) then

            i_bisp = 0
            do j = 0, this%j_max
               i_bisp = i_bisp + 1
               descriptor_out%x(i_desc)%data(i_bisp) =  sum( conjg(U(j)%mm)*U(j)%mm )
            enddo
         endif

         if(my_do_grad_descriptor) then
            n_i = 0
            do n = 1, atoms_n_neighbours(at,i)
               ji = atoms_neighbour(at, i, n, distance=r)
               if( r > this%cutoff ) cycle
               n_i = n_i + 1
               i_bisp = 0
               do j = 0, this%j_max
                  i_bisp = i_bisp + 1
                  do k = 1, 3
                     descriptor_out%x(i_desc)%grad_data(i_bisp,k,n_i) = 2.0_dp * sum( conjg(U(j)%mm)*dU(j,n_i)%mm(k,:,:) )
                  enddo
               enddo 
            enddo
            descriptor_out%x(i_desc)%grad_data(:,:,0) = -sum(descriptor_out%x(i_desc)%grad_data(:,:,:), dim=3)
         endif

         call finalise(dU)
      enddo ! i

      ! clear U from the memory
      call finalise(U)

      call system_timer('power_SO4_calc')

   endsubroutine power_SO4_calc

   subroutine soap_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,error)
      type(soap), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      integer, optional, intent(out) :: error

      type(cplx_1d), dimension(:), allocatable :: SphericalY_ij
      type(cplx_2d), dimension(:), allocatable :: grad_SphericalY_ij

      !SPEED type(cplx_1d), dimension(:,:), allocatable :: fourier_so3
      !SPEED type(cplx_2d), dimension(:,:,:), allocatable :: grad_fourier_so3
      type(real_1d), dimension(:,:), allocatable :: fourier_so3_r, fourier_so3_i
      type(real_2d), dimension(:,:,:), allocatable :: grad_fourier_so3_r, grad_fourier_so3_i
      real(dp), allocatable :: t_g_r(:,:), t_g_i(:,:), t_f_r(:,:), t_f_i(:,:), t_g_f_rr(:,:), t_g_f_ii(:,:)
      integer :: alpha

      logical :: my_do_descriptor, my_do_grad_descriptor
      integer :: d, i, j, n, a, b, k, l, m, i_desc, i_pow, n_neighbours, n_i, n_descriptors, n_cross
      integer, dimension(3) :: shift_ij
      real(dp) :: r_ij, arg_bess, mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lm, mo_spher_bess_fi_ki_lmm, mo_spher_bess_fi_ki_lp, exp_p, exp_m, f_cut, df_cut, norm_descriptor_i
      real(dp), dimension(3) :: u_ij, d_ij
      real(dp), dimension(:,:), allocatable :: radial_fun, radial_coefficient, grad_radial_fun, grad_radial_coefficient, grad_descriptor_i
      real(dp), dimension(:), allocatable :: descriptor_i

      INIT_ERROR(error)

      call system_timer('soap_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("soap_calc: descriptor object not initialised", error)
      endif

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      call finalise(descriptor_out)

      d = soap_dimensions(this,error)
      allocate(descriptor_i(d))
      if(my_do_grad_descriptor) allocate(grad_descriptor_i(d,3))

      call descriptor_sizes(this,at,n_descriptors,n_cross,error)

      allocate(descriptor_out%x(n_descriptors))

      allocate(radial_fun(0:this%l_max, this%n_max), radial_coefficient(0:this%l_max, this%n_max)) 
      !SPEED allocate(fourier_so3(0:this%l_max,this%n_max), SphericalY_ij(0:this%l_max))
      allocate(fourier_so3_r(0:this%l_max,this%n_max), fourier_so3_i(0:this%l_max,this%n_max), SphericalY_ij(0:this%l_max))

      if(my_do_grad_descriptor) then
         allocate(grad_radial_fun(0:this%l_max, this%n_max), grad_radial_coefficient(0:this%l_max, this%n_max))
         allocate(grad_SphericalY_ij(0:this%l_max))
      endif

      do a = 1, this%n_max
         do l = 0, this%l_max
            !SPEED allocate(fourier_so3(l,a)%m(-l:l))
            !SPEED fourier_so3(l,a)%m(:) = CPLX_ZERO
            allocate(fourier_so3_r(l,a)%m(-l:l))
            allocate(fourier_so3_i(l,a)%m(-l:l))
            fourier_so3_r(l,a)%m(:) = 0.0_dp
            fourier_so3_i(l,a)%m(:) = 0.0_dp
         enddo
      enddo
      do l = 0, this%l_max
         allocate(SphericalY_ij(l)%m(-l:l))
         if(my_do_grad_descriptor) allocate(grad_SphericalY_ij(l)%mm(3,-l:l))
      enddo

      i_desc = 0
      do i = 1, at%N
         i_desc = i_desc + 1

         if(my_do_descriptor) then
            allocate(descriptor_out%x(i_desc)%data(d))
            !slow, no need
	    !descriptor_out%x(i_desc)%data = 0.0_dp
            descriptor_out%x(i_desc)%has_data = .false.
            descriptor_out%x(i_desc)%covariance_cutoff = 1.0_dp
         endif
         if(my_do_grad_descriptor) then
            n_neighbours = atoms_n_neighbours(at,i,max_dist=this%cutoff)

            allocate(descriptor_out%x(i_desc)%grad_data(d,3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%ii(0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%pos(3,0:n_neighbours))
            allocate(descriptor_out%x(i_desc)%has_grad_data(0:n_neighbours))
	    ! slow, no need
            ! descriptor_out%x(i_desc)%grad_data = 0.0_dp
            descriptor_out%x(i_desc)%grad_data(:,:,0) = 0.0_dp
            descriptor_out%x(i_desc)%ii = 0
            descriptor_out%x(i_desc)%pos = 0.0_dp
            descriptor_out%x(i_desc)%has_grad_data = .false.

            allocate(descriptor_out%x(i_desc)%grad_covariance_cutoff(3,0:n_neighbours))
            descriptor_out%x(i_desc)%grad_covariance_cutoff = 0.0_dp
         endif
      enddo

      i_desc = 0
      do i = 1, at%N

         i_desc = i_desc + 1

         if(my_do_descriptor) descriptor_out%x(i_desc)%has_data = .true.
         if(my_do_grad_descriptor) then
            descriptor_out%x(i_desc)%ii(0) = i
            descriptor_out%x(i_desc)%pos(:,0) = at%pos(:,i)
            descriptor_out%x(i_desc)%has_grad_data(0) = .true.
            !SPEED allocate( grad_fourier_so3(0:this%l_max,this%n_max,atoms_n_neighbours(at,i,max_dist=this%cutoff)) )
            allocate( grad_fourier_so3_r(0:this%l_max,this%n_max,atoms_n_neighbours(at,i,max_dist=this%cutoff)) )
            allocate( grad_fourier_so3_i(0:this%l_max,this%n_max,atoms_n_neighbours(at,i,max_dist=this%cutoff)) )
         endif

         !do a = 1, this%n_max
         !   radial_fun(0,a) = exp( -this%alpha * this%r_basis(a)**2 ) !* this%r_basis(a)
         !enddo
         !radial_coefficient(0,:) = matmul( radial_fun(0,:), this%transform_basis )
         radial_fun(0,:) = 0.0_dp
         radial_fun(0,1) = 1.0_dp
         radial_coefficient(0,:) = matmul( radial_fun(0,:), this%cholesky_overlap_basis)

         do a = 1, this%n_max
            !SPEED fourier_so3(0,a)%m(0) = radial_coefficient(0,a) * SphericalYCartesian(0,0,(/0.0_dp, 0.0_dp, 0.0_dp/))
            fourier_so3_r(0,a)%m(0) = real(radial_coefficient(0,a) * SphericalYCartesian(0,0,(/0.0_dp, 0.0_dp, 0.0_dp/)), dp)
            fourier_so3_i(0,a)%m(0) = aimag(radial_coefficient(0,a) * SphericalYCartesian(0,0,(/0.0_dp, 0.0_dp, 0.0_dp/)))
            do l = 1, this%l_max
               !SPEED fourier_so3(l,a)%m(:) = CPLX_ZERO
               fourier_so3_r(l,a)%m(:) = 0.0_dp
               fourier_so3_i(l,a)%m(:) = 0.0_dp
            enddo
         enddo

! soap_calc 20 takes 0.0052 s
! call system_timer("soap_calc 20")
         n_i = 0
         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, shift=shift_ij)
            if( r_ij > this%cutoff ) cycle

            n_i = n_i + 1
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc)%ii(n_i) = j
               descriptor_out%x(i_desc)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc)%has_grad_data(n_i) = .true.
            endif

            f_cut = coordination_function(r_ij,this%cutoff, this%cutoff_transition_width)
            if(my_do_grad_descriptor) then
               df_cut = dcoordination_function(r_ij,this%cutoff, this%cutoff_transition_width)
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     !SPEED allocate(grad_fourier_so3(l,a,n_i)%mm(3,-l:l))
                     !SPEED grad_fourier_so3(l,a,n_i)%mm(:,:) = CPLX_ZERO
                     allocate(grad_fourier_so3_r(l,a,n_i)%mm(3,-l:l))
                     allocate(grad_fourier_so3_i(l,a,n_i)%mm(3,-l:l))
                     grad_fourier_so3_r(l,a,n_i)%mm(:,:) = 0.0_dp
                     grad_fourier_so3_i(l,a,n_i)%mm(:,:) = 0.0_dp
                  enddo
               enddo
            endif

            do a = 1, this%n_max
               arg_bess = 2.0_dp * this%alpha * r_ij * this%r_basis(a)
               exp_p = exp( -this%alpha*( r_ij + this%r_basis(a) )**2 )
               exp_m = exp( -this%alpha*( r_ij - this%r_basis(a) )**2 )

               do l = 0, this%l_max
                  if( l == 0 ) then
                     if(arg_bess == 0.0_dp) then
                        !mo_spher_bess_fi_ki_l = 1.0_dp
                        mo_spher_bess_fi_ki_l = exp( -this%alpha * (this%r_basis(a)**2 + r_ij**2) )
                        if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = 0.0_dp
                     else
                        !mo_spher_bess_fi_ki_lm = cosh(arg_bess)/arg_bess
                        !mo_spher_bess_fi_ki_l = sinh(arg_bess)/arg_bess
                        mo_spher_bess_fi_ki_lm = 0.5_dp * (exp_m + exp_p) / arg_bess
                        mo_spher_bess_fi_ki_l  = 0.5_dp * (exp_m + exp_p) / arg_bess
                        if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                     endif
                  else
                     if(arg_bess == 0.0_dp) then
                        mo_spher_bess_fi_ki_l = 0.0_dp
                        if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = 0.0_dp
                     else
                        mo_spher_bess_fi_ki_lmm = mo_spher_bess_fi_ki_lm
                        mo_spher_bess_fi_ki_lm = mo_spher_bess_fi_ki_l
                        if(my_do_grad_descriptor) then
                           mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lp
                           mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                        else
                           mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm - (2*l-1)*mo_spher_bess_fi_ki_lm / arg_bess
                        endif
                     endif
                  endif

                  !radial_fun(l,a) = exp( -this%alpha * (this%r_basis(a)**2 + r_ij**2) ) * mo_spher_bess_fi_ki_l !* this%r_basis(a)
                  radial_fun(l,a) = mo_spher_bess_fi_ki_l !* this%r_basis(a)
                  if(my_do_grad_descriptor) grad_radial_fun(l,a) = -2.0_dp * this%alpha * r_ij * mo_spher_bess_fi_ki_l + &
                     l*mo_spher_bess_fi_ki_l / r_ij + mo_spher_bess_fi_ki_lp * 2.0_dp * this%alpha * this%r_basis(a)

               enddo
            enddo

            radial_coefficient = matmul( radial_fun, this%transform_basis )
            if(my_do_grad_descriptor) grad_radial_coefficient = matmul( grad_radial_fun, this%transform_basis ) * f_cut + radial_coefficient * df_cut
            radial_coefficient = radial_coefficient * f_cut

            do l = 0, this%l_max
               do m = -l, l
                  SphericalY_ij(l)%m(m) = SphericalYCartesian(l,m,d_ij)
                  if(my_do_grad_descriptor) grad_SphericalY_ij(l)%mm(:,m) = GradSphericalYCartesian(l,m,d_ij)
               enddo
            enddo

            do a = 1, this%n_max
               do l = 0, this%l_max
                  do m = -l, l
                     !SPEED fourier_so3(l,a)%m(m) = fourier_so3(l,a)%m(m) + radial_coefficient(l,a) * SphericalY_ij(l)%m(m)
                     !SPEED if(my_do_grad_descriptor) grad_fourier_so3(l,a,n_i)%mm(:,m) = grad_fourier_so3(l,a,n_i)%mm(:,m) + &
                        !SPEED grad_radial_coefficient(l,a) * SphericalY_ij(l)%m(m) * u_ij + radial_coefficient(l,a) * grad_SphericalY_ij(l)%mm(:,m)
                     fourier_so3_r(l,a)%m(m) = fourier_so3_r(l,a)%m(m) + real(radial_coefficient(l,a) * SphericalY_ij(l)%m(m), dp)
                     fourier_so3_i(l,a)%m(m) = fourier_so3_i(l,a)%m(m) + aimag(radial_coefficient(l,a) * SphericalY_ij(l)%m(m))
                     if(my_do_grad_descriptor) then
			grad_fourier_so3_r(l,a,n_i)%mm(:,m) = grad_fourier_so3_r(l,a,n_i)%mm(:,m) + &
			   real(grad_radial_coefficient(l,a) * SphericalY_ij(l)%m(m) * u_ij + radial_coefficient(l,a) * grad_SphericalY_ij(l)%mm(:,m), dp)
			grad_fourier_so3_i(l,a,n_i)%mm(:,m) = grad_fourier_so3_i(l,a,n_i)%mm(:,m) + &
			   aimag(grad_radial_coefficient(l,a) * SphericalY_ij(l)%m(m) * u_ij + radial_coefficient(l,a) * grad_SphericalY_ij(l)%mm(:,m))
		     endif
                  enddo
               enddo
            enddo

         enddo ! n
! call system_timer("soap_calc 20")

         i_pow = 0
         do a = 1, this%n_max
            do b = 1, this%n_max
               do l = 0, this%l_max
                  i_pow = i_pow + 1
                  !SPEED descriptor_i(i_pow) = real( dot_product(fourier_so3(l,a)%m, fourier_so3(l,b)%m) )
                  descriptor_i(i_pow) = dot_product(fourier_so3_r(l,a)%m, fourier_so3_r(l,b)%m) + dot_product(fourier_so3_i(l,a)%m, fourier_so3_i(l,b)%m)
               enddo !l
            enddo !b
         enddo !a
         norm_descriptor_i = sqrt(dot_product(descriptor_i,descriptor_i))

         if(my_do_descriptor) descriptor_out%x(i_desc)%data = descriptor_i / norm_descriptor_i

         if(my_do_grad_descriptor) then
! soap_calc 33 takes 0.047 s
! call system_timer("soap_calc 33")
	    allocate(t_g_r(this%n_max*3, 2*this%l_max+1), t_g_i(this%n_max*3, 2*this%l_max+1))
	    allocate(t_f_r(this%n_max, 2*this%l_max+1), t_f_i(this%n_max, 2*this%l_max+1))
	    allocate(t_g_f_rr(this%n_max*3, this%n_max), t_g_f_ii(this%n_max*3, this%n_max))
            do n_i = 1, atoms_n_neighbours(at,i,max_dist=this%cutoff)
               !SPEED i_pow = 0
               !SPEED do a = 1, this%n_max
                  !SPEED do b = 1, this%n_max
                     !SPEED do l = 0, this%l_max
                        !SPEED i_pow = i_pow + 1
                        !SPEED grad_descriptor_i(i_pow,:) = real( matmul(conjg(grad_fourier_so3(l,a,n_i)%mm),fourier_so3(l,b)%m) + matmul(grad_fourier_so3(l,b,n_i)%mm,conjg(fourier_so3(l,a)%m)) )
                     !SPEED enddo !l
                  !SPEED enddo !b
               !SPEED enddo !a
	       do l=0, this%l_max
		  do a = 1, this%n_max
		     do alpha=1, 3
			t_g_r(3*(a-1)+alpha, 1:2*l+1) = grad_fourier_so3_r(l,a,n_i)%mm(alpha,-l:l)
			t_g_i(3*(a-1)+alpha, 1:2*l+1) = grad_fourier_so3_i(l,a,n_i)%mm(alpha,-l:l)
		     end do
		  end do
		  do b = 1, this%n_max
		     t_f_r(b, 1:2*l+1) = fourier_so3_r(l,b)%m(-l:l)
		     t_f_i(b, 1:2*l+1) = fourier_so3_i(l,b)%m(-l:l)
		  end do
		  call dgemm('N','T',this%n_max*3, this%n_max, 2*l+1, 1.0_dp, &
		     t_g_r(1,1), size(t_g_r,1), t_f_r(1,1), size(t_f_r,1), 0.0_dp, t_g_f_rr(1,1), size(t_g_f_rr, 1))
		  call dgemm('N','T',this%n_max*3, this%n_max, 2*l+1, 1.0_dp, &
		     t_g_i(1,1), size(t_g_i,1), t_f_i(1,1), size(t_f_i,1), 0.0_dp, t_g_f_ii(1,1), size(t_g_f_ii, 1))

		  i_pow = l+1
		  do a = 1, this%n_max
		     do b = 1, this%n_max
			grad_descriptor_i(i_pow, 1:3) = t_g_f_rr(3*(a-1)+1:3*a,b) + t_g_f_ii(3*(a-1)+1:3*a,b) + &
			                                t_g_f_rr(3*(b-1)+1:3*b,a) + t_g_f_ii(3*(b-1)+1:3*b,a)
			i_pow = i_pow + this%l_max+1
		     end do
		  end do

	       end do !l
               descriptor_out%x(i_desc)%grad_data(:,:,n_i) = grad_descriptor_i / norm_descriptor_i
               do k = 1, 3
                  descriptor_out%x(i_desc)%grad_data(:,k,n_i) = descriptor_out%x(i_desc)%grad_data(:,k,n_i) - descriptor_i * dot_product(descriptor_i,grad_descriptor_i(:,k)) / norm_descriptor_i**3
               enddo

               descriptor_out%x(i_desc)%grad_data(:,:,0) = descriptor_out%x(i_desc)%grad_data(:,:,0) - descriptor_out%x(i_desc)%grad_data(:,:,n_i)
            enddo !ni
	    deallocate(t_f_r, t_f_i)
	    deallocate(t_g_r, t_g_i)
	    deallocate(t_g_f_rr, t_g_f_ii)
! call system_timer("soap_calc 33")

            do n_i = 1, atoms_n_neighbours(at,i,max_dist=this%cutoff)
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     !SPEED deallocate(grad_fourier_so3(l,a,n_i)%mm)
                     deallocate(grad_fourier_so3_r(l,a,n_i)%mm)
                     deallocate(grad_fourier_so3_i(l,a,n_i)%mm)
                  enddo
               enddo
            enddo
            !SPEED deallocate(grad_fourier_so3)
            deallocate(grad_fourier_so3_r)
            deallocate(grad_fourier_so3_i)
         endif

      enddo ! i

      !SPEED if(allocated(fourier_so3)) then
         !SPEED do a = lbound(fourier_so3,2), ubound(fourier_so3,2)
            !SPEED do l = lbound(fourier_so3,1), ubound(fourier_so3,1)
               !SPEED deallocate(fourier_so3(l,a)%m)
            !SPEED enddo
         !SPEED enddo
         !SPEED deallocate(fourier_so3)
      !SPEED endif
      if(allocated(fourier_so3_r)) then
         do a = lbound(fourier_so3_r,2), ubound(fourier_so3_r,2)
            do l = lbound(fourier_so3_r,1), ubound(fourier_so3_r,1)
               deallocate(fourier_so3_r(l,a)%m)
            enddo
         enddo
         deallocate(fourier_so3_r)
      endif
      if(allocated(fourier_so3_i)) then
         do a = lbound(fourier_so3_i,2), ubound(fourier_so3_i,2)
            do l = lbound(fourier_so3_i,1), ubound(fourier_so3_i,1)
               deallocate(fourier_so3_i(l,a)%m)
            enddo
         enddo
         deallocate(fourier_so3_i)
      endif

      if(allocated(SphericalY_ij)) then
         do l = lbound(SphericalY_ij,1), ubound(SphericalY_ij,1)
            deallocate(SphericalY_ij(l)%m)
         enddo
         deallocate(SphericalY_ij)
      endif

      if(allocated(grad_SphericalY_ij)) then
         do l = lbound(grad_SphericalY_ij,1), ubound(grad_SphericalY_ij,1)
            deallocate(grad_SphericalY_ij(l)%mm)
         enddo
         deallocate(grad_SphericalY_ij)
      endif

      if(allocated(radial_fun)) deallocate(radial_fun)
      if(allocated(radial_coefficient)) deallocate(radial_coefficient)
      if(allocated(grad_radial_fun)) deallocate(grad_radial_fun)
      if(allocated(grad_radial_coefficient)) deallocate(grad_radial_coefficient)
      if(allocated(descriptor_i)) deallocate(descriptor_i)
      if(allocated(grad_descriptor_i)) deallocate(grad_descriptor_i)

      call system_timer('soap_calc')

   endsubroutine soap_calc

   function descriptor_dimensions(this,error)
      type(descriptor), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: descriptor_dimensions

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            descriptor_dimensions = bispectrum_SO4_dimensions(this%descriptor_bispectrum_SO4,error)
         case(DT_BISPECTRUM_SO3)
            descriptor_dimensions = bispectrum_SO3_dimensions(this%descriptor_bispectrum_SO3,error)
         case(DT_BEHLER)
            descriptor_dimensions = behler_dimensions(this%descriptor_behler,error)
         case(DT_DISTANCE_2b)
            descriptor_dimensions = distance_2b_dimensions(this%descriptor_distance_2b,error)
         case(DT_COORDINATION)
            descriptor_dimensions = coordination_dimensions(this%descriptor_coordination,error)
         case(DT_ANGLE_3B)
            descriptor_dimensions = angle_3b_dimensions(this%descriptor_angle_3b,error)
         case(DT_CO_ANGLE_3B)
            descriptor_dimensions = co_angle_3b_dimensions(this%descriptor_co_angle_3b,error)
         case(DT_CO_DISTANCE_2b)
            descriptor_dimensions = co_distance_2b_dimensions(this%descriptor_co_distance_2b,error)
         case(DT_COSNX)
            descriptor_dimensions = cosnx_dimensions(this%descriptor_cosnx,error)
         case(DT_TRIHIS)
            descriptor_dimensions = trihis_dimensions(this%descriptor_trihis,error)
         case(DT_WATER_MONOMER)
            descriptor_dimensions = water_monomer_dimensions(this%descriptor_water_monomer,error)
         case(DT_WATER_DIMER)
            descriptor_dimensions = water_dimer_dimensions(this%descriptor_water_dimer,error)
         case(DT_A2_DIMER)
            descriptor_dimensions = A2_dimer_dimensions(this%descriptor_A2_dimer,error)
         case(DT_AB_DIMER)
            descriptor_dimensions = AB_dimer_dimensions(this%descriptor_AB_dimer,error)
         case(DT_BOND_REAL_SPACE)
            descriptor_dimensions = bond_real_space_dimensions(this%descriptor_bond_real_space,error)
         case(DT_ATOM_REAL_SPACE)
            descriptor_dimensions = atom_real_space_dimensions(this%descriptor_atom_real_space,error)
         case(DT_POWER_SO3)
            descriptor_dimensions = power_so3_dimensions(this%descriptor_power_so3,error)
         case(DT_POWER_SO4)
            descriptor_dimensions = power_so4_dimensions(this%descriptor_power_so4,error)
         case(DT_SOAP)
            descriptor_dimensions = soap_dimensions(this%descriptor_soap,error)
         case default
            RAISE_ERROR("descriptor_dimensions: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endfunction descriptor_dimensions

   function bispectrum_SO4_dimensions(this,error) result(i)
      type(bispectrum_SO4), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i
      integer :: j, j1, j2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO4_dimensions: descriptor object not initialised", error)
      endif

      i = 0
      do j1 = 0, this%j_max
         j2 = j1
         !do j2 = 0, this%j_max
            do j = abs(j1-j2), min(this%j_max,j1+j2)
               if( mod(j1+j2+j,2) == 1 ) cycle
               i = i + 1
            enddo
         !enddo
      enddo

   endfunction bispectrum_SO4_dimensions

   function bispectrum_SO3_dimensions(this,error) result(i)
      type(bispectrum_SO3), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i
      integer :: a, l1, l2, l

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO3_dimensions: descriptor object not initialised", error)
      endif

      i = 0
      do a = 1, this%n_max
         do l1 = 0, this%l_max
            l2 = l1
            !do l2 = 0, this%l_max
               do l = abs(l1-l2), min(this%l_max,l1+l2)
                  if( mod(l1,2)==1 .and. mod(l2,2)==1 .and. mod(l,2)==1 ) cycle
                  i = i + 1
               enddo
            !enddo
         enddo
      enddo

   endfunction bispectrum_SO3_dimensions

   function behler_dimensions(this,error) result(i)
      type(behler), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("behler_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_g2 + this%n_g3

   endfunction behler_dimensions

   function distance_2b_dimensions(this,error) result(i)
      type(distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_2b_dimensions: descriptor object not initialised", error)
      endif

      i = 1

   endfunction distance_2b_dimensions

   function coordination_dimensions(this,error) result(i)
      type(coordination), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("coordination_dimensions: descriptor object not initialised", error)
      endif

      i = 1

   endfunction coordination_dimensions

   function angle_3b_dimensions(this,error) result(i)
      type(angle_3b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("angle_3b_dimensions: descriptor object not initialised", error)
      endif

      i = 3

   endfunction angle_3b_dimensions

   function co_angle_3b_dimensions(this,error) result(i)
      type(co_angle_3b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_angle_3b_dimensions: descriptor object not initialised", error)
      endif

      i = 4

   endfunction co_angle_3b_dimensions

   function co_distance_2b_dimensions(this,error) result(i)
      type(co_distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_distance_2b_dimensions: descriptor object not initialised", error)
      endif

      i = 3

   endfunction co_distance_2b_dimensions

   function cosnx_dimensions(this,error) result(i)
      type(cosnx), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("cosnx_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_max*(this%l_max+1)

   endfunction cosnx_dimensions

   function trihis_dimensions(this,error) result(i)
      type(trihis), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("trihis_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_gauss

   endfunction trihis_dimensions

   function water_monomer_dimensions(this,error) result(i)
      type(water_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_monomer_dimensions: descriptor object not initialised", error)
      endif

      i = 3

   endfunction water_monomer_dimensions

   function water_dimer_dimensions(this,error) result(i)
      type(water_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_dimer_dimensions: descriptor object not initialised", error)
      endif

      i = 15

   endfunction water_dimer_dimensions

   function A2_dimer_dimensions(this,error) result(i)
      type(A2_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_dimensions: descriptor object not initialised", error)
      endif

      i = 6

   endfunction A2_dimer_dimensions

   function AB_dimer_dimensions(this,error) result(i)
      type(AB_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("AB_dimer_dimensions: descriptor object not initialised", error)
      endif

      i = 6

   endfunction AB_dimer_dimensions

   function bond_real_space_dimensions(this,error) result(i)
      type(bond_real_space), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bond_real_space_dimensions: descriptor object not initialised", error)
      endif

      i = 1

   endfunction bond_real_space_dimensions

   function atom_real_space_dimensions(this,error) result(i)
      type(atom_real_space), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("atom_real_space_dimensions: descriptor object not initialised", error)
      endif

      i = 2 * (this%l_max+1)**2 + 2

   endfunction atom_real_space_dimensions

   function power_so3_dimensions(this,error) result(i)
      type(power_so3), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_so3_dimensions: descriptor object not initialised", error)
      endif

      i = this%n_max*(this%l_max+1)

   endfunction power_so3_dimensions

   function power_SO4_dimensions(this,error) result(i)
      type(power_SO4), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_SO4_dimensions: descriptor object not initialised", error)
      endif

      i = this%j_max + 1

   endfunction power_SO4_dimensions

   function soap_dimensions(this,error) result(i)
      type(soap), intent(in) :: this
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_dimensions: descriptor object not initialised", error)
      endif

      i = (this%l_max+1) * this%n_max**2

   endfunction soap_dimensions

   function descriptor_cutoff(this,error)
      type(descriptor), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: descriptor_cutoff

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            descriptor_cutoff = cutoff(this%descriptor_bispectrum_SO4,error)
         case(DT_BISPECTRUM_SO3)
            descriptor_cutoff = cutoff(this%descriptor_bispectrum_SO3,error)
         case(DT_BEHLER)
            descriptor_cutoff = cutoff(this%descriptor_behler,error)
         case(DT_DISTANCE_2b)
            descriptor_cutoff = cutoff(this%descriptor_distance_2b,error)
         case(DT_COORDINATION)
            descriptor_cutoff = cutoff(this%descriptor_coordination,error)
         case(DT_ANGLE_3B)
            descriptor_cutoff = cutoff(this%descriptor_angle_3b,error)
         case(DT_CO_ANGLE_3B)
            descriptor_cutoff = cutoff(this%descriptor_co_angle_3b,error)
         case(DT_CO_DISTANCE_2b)
            descriptor_cutoff = cutoff(this%descriptor_co_distance_2b,error)
         case(DT_COSNX)
            descriptor_cutoff = cutoff(this%descriptor_cosnx,error)
         case(DT_TRIHIS)
            descriptor_cutoff = cutoff(this%descriptor_trihis,error)
         case(DT_WATER_MONOMER)
            descriptor_cutoff = cutoff(this%descriptor_water_monomer,error)
         case(DT_WATER_DIMER)
            descriptor_cutoff = cutoff(this%descriptor_water_dimer,error)
         case(DT_A2_DIMER)
            descriptor_cutoff = cutoff(this%descriptor_A2_dimer,error)
         case(DT_AB_DIMER)
            descriptor_cutoff = cutoff(this%descriptor_AB_dimer,error)
         case(DT_BOND_REAL_SPACE)
            descriptor_cutoff = cutoff(this%descriptor_bond_real_space,error)
         case(DT_ATOM_REAL_SPACE)
            descriptor_cutoff = cutoff(this%descriptor_atom_real_space,error)
         case(DT_POWER_SO3)
            descriptor_cutoff = cutoff(this%descriptor_power_so3,error)
         case(DT_POWER_SO4)
            descriptor_cutoff = cutoff(this%descriptor_power_so4,error)
         case(DT_SOAP)
            descriptor_cutoff = cutoff(this%descriptor_soap,error)
         case default
            RAISE_ERROR("descriptor_cutoff: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endfunction descriptor_cutoff

   function bispectrum_SO4_cutoff(this,error) 
      type(bispectrum_SO4), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: bispectrum_SO4_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO4_cutoff: descriptor object not initialised", error)
      endif

      bispectrum_SO4_cutoff = this%cutoff

   endfunction bispectrum_SO4_cutoff

   function bispectrum_SO3_cutoff(this,error) 
      type(bispectrum_SO3), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: bispectrum_SO3_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO3_cutoff: descriptor object not initialised", error)
      endif

      bispectrum_SO3_cutoff = this%cutoff

   endfunction bispectrum_SO3_cutoff

   function behler_cutoff(this,error) 
      type(behler), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: behler_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("behler_cutoff: descriptor object not initialised", error)
      endif

      behler_cutoff = this%cutoff

   endfunction behler_cutoff

   function distance_2b_cutoff(this,error) 
      type(distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: distance_2b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_2b_cutoff: descriptor object not initialised", error)
      endif

      distance_2b_cutoff = this%cutoff

   endfunction distance_2b_cutoff

   function co_distance_2b_cutoff(this,error) 
      type(co_distance_2b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: co_distance_2b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_distance_2b_cutoff: descriptor object not initialised", error)
      endif

      co_distance_2b_cutoff = this%cutoff

   endfunction co_distance_2b_cutoff

   function coordination_cutoff(this,error) 
      type(coordination), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: coordination_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("coordination_cutoff: descriptor object not initialised", error)
      endif

      coordination_cutoff = this%cutoff

   endfunction coordination_cutoff

   function angle_3b_cutoff(this,error) 
      type(angle_3b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: angle_3b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("angle_3b_cutoff: descriptor object not initialised", error)
      endif

      angle_3b_cutoff = this%cutoff

   endfunction angle_3b_cutoff

   function co_angle_3b_cutoff(this,error) 
      type(co_angle_3b), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: co_angle_3b_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_angle_3b_cutoff: descriptor object not initialised", error)
      endif

      co_angle_3b_cutoff = this%cutoff

   endfunction co_angle_3b_cutoff

   function cosnx_cutoff(this,error) 
      type(cosnx), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: cosnx_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("cosnx_cutoff: descriptor object not initialised", error)
      endif

      cosnx_cutoff = this%cutoff

   endfunction cosnx_cutoff

   function trihis_cutoff(this,error) 
      type(trihis), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: trihis_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("trihis_cutoff: descriptor object not initialised", error)
      endif

      trihis_cutoff = this%cutoff

   endfunction trihis_cutoff

   function water_monomer_cutoff(this,error) 
      type(water_monomer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: water_monomer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_monomer_cutoff: descriptor object not initialised", error)
      endif

      water_monomer_cutoff = this%cutoff

   endfunction water_monomer_cutoff

   function water_dimer_cutoff(this,error) 
      type(water_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: water_dimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_dimer_cutoff: descriptor object not initialised", error)
      endif

      water_dimer_cutoff = this%cutoff

   endfunction water_dimer_cutoff

   function A2_dimer_cutoff(this,error) 
      type(A2_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: A2_dimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_cutoff: descriptor object not initialised", error)
      endif

      A2_dimer_cutoff = this%cutoff

   endfunction A2_dimer_cutoff

   function AB_dimer_cutoff(this,error) 
      type(AB_dimer), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: AB_dimer_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("AB_dimer_cutoff: descriptor object not initialised", error)
      endif

      AB_dimer_cutoff = this%cutoff

   endfunction AB_dimer_cutoff

   function bond_real_space_cutoff(this,error)
      type(bond_real_space), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: bond_real_space_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bond_real_space_cutoff: descriptor object not initialised", error)
      endif

      bond_real_space_cutoff = max(this%cutoff, this%bond_cutoff)

   endfunction bond_real_space_cutoff

   function atom_real_space_cutoff(this,error)
      type(atom_real_space), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: atom_real_space_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("atom_real_space_cutoff: descriptor object not initialised", error)
      endif

      atom_real_space_cutoff = this%cutoff

   endfunction atom_real_space_cutoff

   function power_so3_cutoff(this,error) 
      type(power_so3), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: power_so3_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_so3_cutoff: descriptor object not initialised", error)
      endif

      power_so3_cutoff = this%cutoff

   endfunction power_so3_cutoff

   function power_so4_cutoff(this,error) 
      type(power_so4), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: power_so4_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_so4_cutoff: descriptor object not initialised", error)
      endif

      power_so4_cutoff = this%cutoff

   endfunction power_so4_cutoff

   function soap_cutoff(this,error) 
      type(soap), intent(in) :: this
      integer, optional, intent(out) :: error
      real(dp) :: soap_cutoff

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_cutoff: descriptor object not initialised", error)
      endif

      soap_cutoff = this%cutoff

   endfunction soap_cutoff

   subroutine descriptor_sizes(this,at,n_descriptors,n_cross,error)
      type(descriptor), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4)
            call bispectrum_SO4_sizes(this%descriptor_bispectrum_SO4,at,n_descriptors,n_cross,error=error)
         case(DT_BISPECTRUM_SO3)
            call bispectrum_SO3_sizes(this%descriptor_bispectrum_SO3,at,n_descriptors,n_cross,error=error)
         case(DT_BEHLER)
            call behler_sizes(this%descriptor_behler,at,n_descriptors,n_cross,error=error)
         case(DT_DISTANCE_2b)
            call distance_2b_sizes(this%descriptor_distance_2b,at,n_descriptors,n_cross,error=error)
         case(DT_COORDINATION)
            call coordination_sizes(this%descriptor_coordination,at,n_descriptors,n_cross,error=error)
         case(DT_ANGLE_3B)
            call angle_3b_sizes(this%descriptor_angle_3b,at,n_descriptors,n_cross,error=error)
         case(DT_CO_ANGLE_3B)
            call co_angle_3b_sizes(this%descriptor_co_angle_3b,at,n_descriptors,n_cross,error=error)
         case(DT_CO_DISTANCE_2b)
            call co_distance_2b_sizes(this%descriptor_co_distance_2b,at,n_descriptors,n_cross,error=error)
         case(DT_COSNX)
            call cosnx_sizes(this%descriptor_cosnx,at,n_descriptors,n_cross,error=error)
         case(DT_TRIHIS)
            call trihis_sizes(this%descriptor_trihis,at,n_descriptors,n_cross,error=error)
         case(DT_WATER_MONOMER)
            call water_monomer_sizes(this%descriptor_water_monomer,at,n_descriptors,n_cross,error=error)
         case(DT_WATER_DIMER)
            call water_dimer_sizes(this%descriptor_water_dimer,at,n_descriptors,n_cross,error=error)
         case(DT_A2_DIMER)
            call A2_dimer_sizes(this%descriptor_A2_dimer,at,n_descriptors,n_cross,error=error)
         case(DT_AB_DIMER)
            call AB_dimer_sizes(this%descriptor_AB_dimer,at,n_descriptors,n_cross,error=error)
         case(DT_BOND_REAL_SPACE)
            call bond_real_space_sizes(this%descriptor_bond_real_space,at,n_descriptors,n_cross,error=error)
         case(DT_ATOM_REAL_SPACE)
            call atom_real_space_sizes(this%descriptor_atom_real_space,at,n_descriptors,n_cross,error=error)
         case(DT_POWER_SO3)
            call power_so3_sizes(this%descriptor_power_so3,at,n_descriptors,n_cross,error=error)
         case(DT_POWER_SO4)
            call power_so4_sizes(this%descriptor_power_so4,at,n_descriptors,n_cross,error=error)
         case(DT_SOAP)
            call soap_sizes(this%descriptor_soap,at,n_descriptors,n_cross,error=error)
         case default
            RAISE_ERROR("descriptor_sizes: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endsubroutine descriptor_sizes

   subroutine bispectrum_SO4_sizes(this,at,n_descriptors,n_cross,error)
      type(bispectrum_SO4), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO4_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + atoms_n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine bispectrum_SO4_sizes

   subroutine bispectrum_SO3_sizes(this,at,n_descriptors,n_cross,error)
      type(bispectrum_SO3), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error
      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bispectrum_SO3_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + atoms_n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine bispectrum_SO3_sizes

   subroutine behler_sizes(this,at,n_descriptors,n_cross,error)
      type(behler), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("behler_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0
      do i = 1, at%N
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + atoms_n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo
   endsubroutine behler_sizes

   subroutine distance_2b_sizes(this,at,n_descriptors,n_cross,error)
      type(distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i, j, n
      logical :: Zi1, Zi2, Zj1, Zj2
      real(dp) :: r_ij

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("distance_2b_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance=r_ij)
            if(r_ij > this%cutoff) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            n_descriptors = n_descriptors + 1
         enddo
      enddo

      n_cross = n_descriptors*2

   endsubroutine distance_2b_sizes

   subroutine coordination_sizes(this,at,n_descriptors,n_cross,error)
      type(coordination), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("coordination_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0
      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + atoms_n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine coordination_sizes

   subroutine angle_3b_sizes(this,at,n_descriptors,n_cross,error)
      type(angle_3b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i, j, k, n, m
      real(dp) :: r_ij, r_ik
      logical :: Zk1, Zk2, Zj1, Zj2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("angle_3b_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle

         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij)
            if( r_ij > this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            do m = 1, atoms_n_neighbours(at,i)
               if( n == m ) cycle

               k = atoms_neighbour(at, i, m, distance = r_ik)
               if( r_ik > this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)
               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               n_descriptors = n_descriptors + 1
            enddo
         enddo
      enddo
      n_cross = n_descriptors * 3

   endsubroutine angle_3b_sizes

   subroutine co_angle_3b_sizes(this,at,n_descriptors,n_cross,error)
      type(co_angle_3b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i, j, k, n, m, n_neighbours_coordination
      real(dp) :: r_ij, r_ik
      logical :: Zk1, Zk2, Zj1, Zj2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_angle_3b_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( (this%Z /=0) .and. (at%Z(i) /= this%Z) ) cycle

         n_neighbours_coordination = atoms_n_neighbours(at,i,max_dist=this%coordination_cutoff)

         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at, i, n, distance = r_ij)
            if( r_ij > this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)

            do m = 1, atoms_n_neighbours(at,i)
               if( n == m ) cycle
               k = atoms_neighbour(at, i, m, distance = r_ik)
               if( r_ik > this%cutoff ) cycle

               Zk1 = (this%Z1 == 0) .or. (at%Z(k) == this%Z1)
               Zk2 = (this%Z2 == 0) .or. (at%Z(k) == this%Z2)
               if( .not. ( ( Zk1 .and. Zj2 ) .or. ( Zk2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

               n_descriptors = n_descriptors + 1
               n_cross = n_cross + 3 + n_neighbours_coordination
            enddo
         enddo
      enddo


   endsubroutine co_angle_3b_sizes

   subroutine co_distance_2b_sizes(this,at,n_descriptors,n_cross,error)
      type(co_distance_2b), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      real(dp) :: r_ij
      integer :: i, j, n
      logical :: Zi1, Zi2, Zj1, Zj2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("co_distance_2b_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         Zi1 = (this%Z1 == 0) .or. (at%Z(i) == this%Z1)
         Zi2 = (this%Z2 == 0) .or. (at%Z(i) == this%Z2)
         do n = 1, atoms_n_neighbours(at,i)
            j = atoms_neighbour(at,i,n,distance=r_ij)
            if( r_ij > this%cutoff ) cycle

            Zj1 = (this%Z1 == 0) .or. (at%Z(j) == this%Z1)
            Zj2 = (this%Z2 == 0) .or. (at%Z(j) == this%Z2)
            if( .not. ( ( Zi1 .and. Zj2 ) .or. ( Zi2 .and. Zj1 ) ) ) cycle ! this pair doesn't belong to the descriptor type

            n_descriptors = n_descriptors + 1
            n_cross = n_cross + 4 + atoms_n_neighbours(at,i,max_dist=this%coordination_cutoff) + atoms_n_neighbours(at,j,max_dist=this%coordination_cutoff)
         enddo
      enddo

   endsubroutine co_distance_2b_sizes

   subroutine cosnx_sizes(this,at,n_descriptors,n_cross,error)
      type(cosnx), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("cosnx_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + atoms_n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine cosnx_sizes

   subroutine trihis_sizes(this,at,n_descriptors,n_cross,error)
      type(trihis), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("trihis_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = at%N

      n_cross = 0

      do i = 1, at%N
         n_cross = n_cross + atoms_n_neighbours(at,i) + 1
      enddo

   endsubroutine trihis_sizes

   subroutine water_monomer_sizes(this,at,n_descriptors,n_cross,error)
      type(water_monomer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_monomer_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if(at%Z(i) == 8) then
            n_descriptors = n_descriptors + 1
            n_cross = n_cross + 3
         endif
      enddo

   endsubroutine water_monomer_sizes

   subroutine water_dimer_sizes(this,at,n_descriptors,n_cross,error)
      type(water_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i, j, n
      real(dp) :: r_ij

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("water_dimer_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if(at%Z(i) == 8) then
            do n = 1, atoms_n_neighbours(at,i)
               j = atoms_neighbour(at,i,n,distance=r_ij)
               if(at%Z(j) == 8 .and. r_ij < this%cutoff) then
                  n_descriptors = n_descriptors + 1
                  n_cross = n_cross + 6
               endif
            enddo
         endif
      enddo
   endsubroutine water_dimer_sizes

   subroutine A2_dimer_sizes(this,at,n_descriptors,n_cross,error)
      type(A2_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i, j, iA1, iA2, iB1, iB2
      integer, dimension(at%N) :: A2_monomer_index
      real(dp) :: r_A1_A2, r_B1_B2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_sizes: descriptor object not initialised", error)
      endif

      call find_A2_monomer(at,this%atomic_number, this%monomer_cutoff, A2_monomer_index)

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         iA1 = i
         iA2 = atoms_neighbour(at,i,A2_monomer_index(i),distance=r_A1_A2)
         if( iA1 > iA2 ) cycle

         do j = i + 1, at%N
            iB1 = j
            iB2 = atoms_neighbour(at,j,A2_monomer_index(j),distance=r_B1_B2)
            if( iB1 > iB2 ) cycle

            r_A1_B1 = distance_min_image(at,iA1,iB1)
            r_A1_B2 = distance_min_image(at,iA1,iB2)

            r_A2_B1 = distance_min_image(at,iA2,iB1)
            r_A2_B2 = distance_min_image(at,iA2,iB2)
            
            if( all( (/r_A1_A2,r_B1_B2,r_A1_B1,r_A1_B2,r_A2_B1,r_A2_B2/) < this%cutoff) ) then
               n_descriptors = n_descriptors + 1
               n_cross = n_cross + 4
            endif
         enddo
      enddo

   endsubroutine A2_dimer_sizes

   subroutine AB_dimer_sizes(this,at,n_descriptors,n_cross,error)
      type(AB_dimer), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i, j, n_monomers, iA1, iA2, iB1, iB2
      integer, dimension(:,:), allocatable :: AB_monomer_index
      real(dp) :: r_A1_A2, r_B1_B2, r_A1_B1, r_A1_B2, r_A2_B1, r_A2_B2

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("A2_dimer_sizes: descriptor object not initialised", error)
      endif

      if( count(at%Z == this%atomic_number1) == count(at%Z == this%atomic_number2) ) then
         n_monomers = count(at%Z == this%atomic_number1)
      else
         RAISE_ERROR("AB_dimer_sizes: number of monomer atoms 1 ("//count(at%Z == this%atomic_number1)//") not equal to number of monomer atoms 2 ("//count(at%Z == this%atomic_number1)//")",error)
      endif

      allocate(AB_monomer_index(2,n_monomers))
      call find_AB_monomer(at,(/this%atomic_number1,this%atomic_number2/), this%monomer_cutoff, AB_monomer_index)

      n_descriptors = 0
      n_cross = 0

      do i = 1, n_monomers
         iA1 = AB_monomer_index(1,i)
         iB1 = AB_monomer_index(2,i)
         do j = i + 1, n_monomers
            iA2 = AB_monomer_index(1,j)
            iB2 = AB_monomer_index(2,j)

            r_A1_B1 = distance_min_image(at,iA1,iB1)
            r_A2_B2 = distance_min_image(at,iA2,iB2)

            r_A1_A2 = distance_min_image(at,iA1,iA2)
            r_B1_B2 = distance_min_image(at,iB1,iB2)

            r_A1_B2 = distance_min_image(at,iA1,iB2)
            r_A2_B1 = distance_min_image(at,iA2,iB1)
            
            if( all( (/r_A1_A2,r_B1_B2,r_A1_B1,r_A1_B2,r_A2_B1,r_A2_B2/) < this%cutoff) ) then
               n_descriptors = n_descriptors + 1
               n_cross = n_cross + 4
            endif
         enddo
      enddo

      deallocate(AB_monomer_index)

   endsubroutine AB_dimer_sizes

   subroutine bond_real_space_sizes(this,at,n_descriptors,n_cross,error)
      type(bond_real_space), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      type(atoms) :: at_copy
      integer :: i, j, k, n, m, shift_j(3)

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("bond_real_space_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         n_descriptors = n_descriptors + atoms_n_neighbours(at, i, max_dist=this%bond_cutoff)

         do n = 1, atoms_n_neighbours(at, i)
            j = atoms_neighbour(at, i, n, shift=shift_j, max_dist=this%bond_cutoff)

            if(j == 0) cycle

            at_copy = at
            call add_atoms(at_copy, 0.5_dp * (at%pos(:,i) + at%pos(:,j) + matmul(at%lattice,shift_j)), 1)
            call calc_connect(at_copy)

            do m = 1, atoms_n_neighbours(at_copy, at%N + 1)
               k = atoms_neighbour(at_copy, at%N + 1, m, max_dist=this%cutoff)

               if(k == 0) cycle

               if(at_copy%pos(:,k) .feq. at_copy%pos(:,at%N + 1)) cycle

               n_cross = n_cross + 1
            enddo

            call finalise(at_copy)
         enddo
      enddo

   endsubroutine bond_real_space_sizes

   subroutine atom_real_space_sizes(this,at,n_descriptors,n_cross,error)
      type(atom_real_space), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("atom_real_space_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = at%N
      n_cross = 0

      do i = 1, at%N
         n_cross = n_cross + atoms_n_neighbours(at,i,max_dist=this%cutoff)*2 
      enddo

   endsubroutine atom_real_space_sizes

   subroutine power_so3_sizes(this,at,n_descriptors,n_cross,error)
      type(power_so3), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_so3_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + atoms_n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine power_so3_sizes

   subroutine power_SO4_sizes(this,at,n_descriptors,n_cross,error)
      type(power_SO4), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("power_SO4_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         if( at%Z(i) /= this%Z .and. this%Z /=0 ) cycle
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + atoms_n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine power_SO4_sizes

   subroutine soap_sizes(this,at,n_descriptors,n_cross,error)
      type(soap), intent(in) :: this
      type(atoms), intent(in) :: at
      integer, intent(out) :: n_descriptors, n_cross
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) then
         RAISE_ERROR("soap_sizes: descriptor object not initialised", error)
      endif

      n_descriptors = 0
      n_cross = 0

      do i = 1, at%N
         n_descriptors = n_descriptors + 1
         n_cross = n_cross + atoms_n_neighbours(at,i,max_dist=this%cutoff) + 1
      enddo

   endsubroutine soap_sizes

   function descriptor_n_permutations(this,error)
      type(descriptor), intent(in) :: this
      integer, optional, intent(out) :: error

      integer :: descriptor_n_permutations

      INIT_ERROR(error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4,DT_BISPECTRUM_SO3,DT_BEHLER,DT_DISTANCE_2b,DT_COORDINATION, &
            DT_ANGLE_3B,DT_CO_ANGLE_3B,DT_CO_DISTANCE_2b,DT_COSNX,DT_TRIHIS,DT_WATER_MONOMER,DT_BOND_REAL_SPACE,DT_ATOM_REAL_SPACE,DT_POWER_SO3,DT_POWER_SO4,DT_SOAP)

            descriptor_n_permutations = 1
            
         case(DT_WATER_DIMER)
            descriptor_n_permutations = NP_WATER_DIMER
         case(DT_A2_DIMER)
            descriptor_n_permutations = NP_A2_DIMER
         case(DT_AB_DIMER)
            descriptor_n_permutations = NP_AB_DIMER
         case default
            RAISE_ERROR("descriptor_permutations: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endfunction descriptor_n_permutations

   subroutine descriptor_permutations(this,permutations,error)
      type(descriptor), intent(in) :: this
      integer, dimension(:,:), intent(out) :: permutations
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      call check_size('permutations',permutations, &
           (/descriptor_dimensions(this,error),descriptor_n_permutations(this,error)/),'descriptor_permutations',error)

      selectcase(this%descriptor_type)
         case(DT_BISPECTRUM_SO4,DT_BISPECTRUM_SO3,DT_BEHLER,DT_DISTANCE_2b,DT_COORDINATION, &
            DT_ANGLE_3B,DT_CO_ANGLE_3B,DT_CO_DISTANCE_2b,DT_COSNX,DT_TRIHIS,DT_WATER_MONOMER,DT_BOND_REAL_SPACE,DT_ATOM_REAL_SPACE,DT_POWER_SO3,DT_POWER_SO4,DT_SOAP)
            
            permutations(:,1) = (/ (i, i = 1, size(permutations,1)) /)
         case(DT_WATER_DIMER)
            permutations(:,1) = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15/) ! original order
            permutations(:,2) = (/1, 3, 2, 4, 5, 7, 6, 8, 9, 10, 13, 14, 11, 12, 15/) ! swap Hs on monomer A
            permutations(:,3) = (/1, 2, 3, 5, 4, 6, 7, 9, 8, 10, 12, 11, 14, 13, 15/) ! swap Hs on monomer B
            permutations(:,4) = (/1, 3, 2, 5, 4, 7, 6, 9, 8, 10, 14, 13, 12, 11, 15/) ! swap Hs on both monomers
            permutations(:,5) = (/1, 8, 9, 6, 7, 4, 5, 2, 3, 15, 11, 13, 12, 14, 10/) ! swap monomers A and B
            permutations(:,6) = (/1, 9, 8, 6, 7, 5, 4, 2, 3, 15, 12, 14, 11, 13, 10/) ! swap monomers and Hs on monomer A
            permutations(:,7) = (/1, 8, 9, 7, 6, 4, 5, 3, 2, 15, 13, 11, 14, 12, 10/) ! swap monomers and Hs on monomer B
            permutations(:,8) = (/1, 9, 8, 7, 6, 5, 4, 3, 2, 15, 14, 12, 13, 11, 10/) ! swap monomers and Hs on both monomers

         case(DT_A2_DIMER)
            permutations(:,1) = (/1, 2, 3, 4, 5, 6/) ! original order
            permutations(:,2) = (/1, 2, 5, 6, 3, 4/) ! swap atoms on monomer A
            permutations(:,3) = (/1, 2, 4, 3, 6, 5/) ! swap atoms on monomer B
            permutations(:,4) = (/1, 2, 6, 5, 4, 3/) ! swap atoms on both monomers
            permutations(:,5) = (/2, 1, 3, 5, 4, 6/) ! swap monomers A and B
            permutations(:,6) = (/2, 1, 5, 3, 6, 4/) ! swap monomers and atoms on monomer A
            permutations(:,7) = (/2, 1, 4, 6, 3, 5/) ! swap monomers and atoms on monomer B
            permutations(:,8) = (/2, 1, 6, 4, 5, 3/) ! swap monomers and atoms on both monomers
            
         case(DT_AB_DIMER)
            permutations(:,1) = (/1, 2, 3, 4, 5, 6/) ! original order
            permutations(:,2) = (/2, 1, 3, 4, 6, 5/) ! swap monomers
            
         case default
            RAISE_ERROR("descriptor_permutations: unknown descriptor type "//this%descriptor_type,error)
      endselect

   endsubroutine descriptor_permutations

   subroutine real_space_fourier_coefficients(at,l_max,atom_coefficient)
      type(atoms), intent(in) :: at
      integer, intent(in) :: l_max
      type(neighbour_type), dimension(:), allocatable :: atom_coefficient

      integer :: i, j, n, l, m
      real(dp) :: r
      real(dp), dimension(3) :: d

      if(.not.allocated(atom_coefficient)) allocate(atom_coefficient(at%N))

      do i = 1, at%N
         if(.not. allocated(atom_coefficient(i)%neighbour)) allocate(atom_coefficient(i)%neighbour(atoms_n_neighbours(at,i)))
         do n = 1, atoms_n_neighbours(at,i)

            j = atoms_neighbour(at,i,n,distance = r, diff = d)
            atom_coefficient(i)%neighbour(n)%r = r
            atom_coefficient(i)%neighbour(n)%u = d / r

            if(.not. allocated(atom_coefficient(i)%neighbour(n)%spherical_harmonics)) allocate( atom_coefficient(i)%neighbour(n)%spherical_harmonics(0:l_max), &
            atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(0:l_max) )
            do l = 0, l_max
               if(.not. allocated(atom_coefficient(i)%neighbour(n)%spherical_harmonics(l)%m)) &
               allocate(atom_coefficient(i)%neighbour(n)%spherical_harmonics(l)%m(-l:l))
               if(.not. allocated(atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(l)%mm)) &
               allocate(atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(l)%mm(3,-l:l))

               atom_coefficient(i)%neighbour(n)%spherical_harmonics(l)%m = CPLX_ZERO
               atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(l)%mm = CPLX_ZERO

               do m = -l, l
                  atom_coefficient(i)%neighbour(n)%spherical_harmonics(l)%m(m) = SphericalYCartesian(l,m,d)
                  atom_coefficient(i)%neighbour(n)%grad_spherical_harmonics(l)%mm(:,m) = GradSphericalYCartesian(l,m,d)
               enddo
            enddo
         enddo
      enddo

   endsubroutine real_space_fourier_coefficients

   function real_space_covariance_coefficient(anc1,anc2,i1,i2,alpha,l_max,f1,f2)
      type(neighbour_type), dimension(:), intent(in) :: anc1, anc2
      real(dp), intent(in) :: alpha
      integer, intent(in) :: i1, i2, l_max
      real(dp), dimension(:,:), intent(out), optional :: f1, f2

      real(dp) :: real_space_covariance_coefficient

      complex(dp) :: real_space_covariance_in, I_lm1m2
      integer :: n1, n2, l, m1, m2, k
      real(dp) :: r1, r2, arg_bess, fac_exp, mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lm, mo_spher_bess_fi_ki_lmm, mo_spher_bess_fi_ki_lp, grad_mo_spher_bess_fi_ki_l
      real(dp), dimension(3) :: u1, u2, grad_arg_bess1, grad_fac_exp1, grad_arg_bess2, grad_fac_exp2
      type(cplx_2d), dimension(:), allocatable :: integral_r
      type(grad_spherical_harmonics_overlap_type), dimension(:), allocatable :: grad_integral_r1, grad_integral_r2

      logical :: do_derivative

      do_derivative = (present(f1) .or. present(f2)) 

      real_space_covariance_in = CPLX_ZERO

      allocate(integral_r(0:l_max))
      do l = 0, l_max
         allocate(integral_r(l)%mm(-l:l,-l:l))
         integral_r(l)%mm = CPLX_ZERO
      enddo

      if(present(f1)) then
         allocate(grad_integral_r1(0:size(anc1(i1)%neighbour)))
         do n1 = 0, size(anc1(i1)%neighbour)
            allocate(grad_integral_r1(n1)%grad_integral(0:l_max))
            do l = 0, l_max
               allocate(grad_integral_r1(n1)%grad_integral(l)%mm(3,-l:l,-l:l))
               grad_integral_r1(n1)%grad_integral(l)%mm = CPLX_ZERO
            enddo
         enddo
      endif

      if(present(f2)) then
         allocate(grad_integral_r2(0:size(anc2(i2)%neighbour)))
         do n2 = 0, size(anc2(i2)%neighbour)
            allocate(grad_integral_r2(n2)%grad_integral(0:l_max))
            do l = 0, l_max
               allocate(grad_integral_r2(n2)%grad_integral(l)%mm(3,-l:l,-l:l))
               grad_integral_r2(n2)%grad_integral(l)%mm = CPLX_ZERO
            enddo
         enddo
      endif
      do n1 = 1, size(anc1(i1)%neighbour)
         r1 = anc1(i1)%neighbour(n1)%r
         u1 = anc1(i1)%neighbour(n1)%u
         do n2 = 1, size(anc2(i2)%neighbour)
            r2 = anc2(i2)%neighbour(n2)%r

            u2 = anc2(i2)%neighbour(n2)%u

            arg_bess = alpha*r1*r2
            fac_exp = exp(-0.5_dp*alpha*(r1**2+r2**2))

            if(present(f1)) then
               grad_arg_bess1 = alpha*r2*u1
               grad_fac_exp1 = -fac_exp*alpha*r1*u1
            endif

            if(present(f2)) then
               grad_arg_bess2 = alpha*r1*u2
               grad_fac_exp2 = -fac_exp*alpha*r2*u2
            endif

            do l = 0, l_max
               if( l == 0 ) then
                  mo_spher_bess_fi_ki_lm = cosh(arg_bess)/arg_bess
                  mo_spher_bess_fi_ki_l = sinh(arg_bess)/arg_bess
                  if(do_derivative) mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
               else
                  mo_spher_bess_fi_ki_lmm = mo_spher_bess_fi_ki_lm
                  mo_spher_bess_fi_ki_lm = mo_spher_bess_fi_ki_l
                  if(do_derivative) then
                     mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lp
                     mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                  else
                     mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm - (2*l-1)*mo_spher_bess_fi_ki_lm / arg_bess
                  endif

               endif


               if(do_derivative) grad_mo_spher_bess_fi_ki_l = 0.5_dp * (mo_spher_bess_fi_ki_lp - mo_spher_bess_fi_ki_l / arg_bess + mo_spher_bess_fi_ki_lm)
                  
               do m1 = -l, l
                  do m2 = -l, l
                     I_lm1m2 = conjg(anc1(i1)%neighbour(n1)%spherical_harmonics(l)%m(m1)) * anc2(i2)%neighbour(n2)%spherical_harmonics(l)%m(m2) * mo_spher_bess_fi_ki_l*fac_exp
                     integral_r(l)%mm(m2,m1) = integral_r(l)%mm(m2,m1) + I_lm1m2
                     if(present(f1)) then
                        grad_integral_r1(n1)%grad_integral(l)%mm(:,m2,m1) = grad_integral_r1(n1)%grad_integral(l)%mm(:,m2,m1) + &
                        anc2(i2)%neighbour(n2)%spherical_harmonics(l)%m(m2) * &
                        ( conjg(anc1(i1)%neighbour(n1)%grad_spherical_harmonics(l)%mm(:,m1)) * mo_spher_bess_fi_ki_l*fac_exp + &
                        conjg(anc1(i1)%neighbour(n1)%spherical_harmonics(l)%m(m1)) * ( grad_mo_spher_bess_fi_ki_l * grad_arg_bess1 * fac_exp + mo_spher_bess_fi_ki_l * grad_fac_exp1 ) )
                     endif

                     if(present(f2)) then
                        grad_integral_r2(n2)%grad_integral(l)%mm(:,m2,m1) = grad_integral_r2(n2)%grad_integral(l)%mm(:,m2,m1) + &
                        conjg(anc1(i1)%neighbour(n1)%spherical_harmonics(l)%m(m1)) * &
                        ( anc2(i2)%neighbour(n2)%grad_spherical_harmonics(l)%mm(:,m2) * mo_spher_bess_fi_ki_l*fac_exp + &
                        anc2(i2)%neighbour(n2)%spherical_harmonics(l)%m(m2) * ( grad_mo_spher_bess_fi_ki_l * grad_arg_bess2 * fac_exp + mo_spher_bess_fi_ki_l * grad_fac_exp2 ) )
                     endif

                  enddo
               enddo
            enddo
         enddo
      enddo

      if(present(f1)) then
         f1 = 0.0_dp
         do n1 = 0, size(anc1(i1)%neighbour)
            do l = 0, l_max
               do k = 1, 3
                  f1(k,n1+1) = f1(k,n1+1) + real(sum(conjg(grad_integral_r1(n1)%grad_integral(l)%mm(k,:,:))*integral_r(l)%mm(:,:)))
               enddo
            enddo
         enddo
         f1 = 2.0_dp * f1
      endif

      if(present(f2)) then
         f2 = 0.0_dp
         do n2 = 0, size(anc2(i2)%neighbour)
            do l = 0, l_max
               do k = 1, 3
                  f2(k,n2+1) = f2(k,n2+1) + real(sum(conjg(grad_integral_r2(n2)%grad_integral(l)%mm(k,:,:))*integral_r(l)%mm(:,:)))
               enddo
            enddo
         enddo
         f2 = 2.0_dp * f2
      endif

      do l = 0, l_max
         real_space_covariance_in = real_space_covariance_in + sum(conjg(integral_r(l)%mm) * integral_r(l)%mm)
      enddo
      real_space_covariance_coefficient = real(real_space_covariance_in)

      do l = 0, l_max
         deallocate(integral_r(l)%mm)
      enddo
      deallocate(integral_r)

      if(present(f1)) then
         do n1 = 0, size(anc1(i1)%neighbour)
            do l = 0, l_max
               deallocate(grad_integral_r1(n1)%grad_integral(l)%mm)
            enddo
            deallocate(grad_integral_r1(n1)%grad_integral)
         enddo
         deallocate(grad_integral_r1)
      endif

      if(present(f2)) then
         do n2 = 0, size(anc2(i2)%neighbour)
            do l = 0, l_max
               deallocate(grad_integral_r2(n2)%grad_integral(l)%mm)
            enddo
            deallocate(grad_integral_r2(n2)%grad_integral)
         enddo
         deallocate(grad_integral_r2)
      endif

   endfunction real_space_covariance_coefficient

   function real_space_covariance(at1,at2,i1,i2,alpha,l_max,f1,f2)
      type(atoms), intent(in) :: at1, at2
      real(dp), intent(in) :: alpha
      integer, intent(in) :: i1, i2, l_max
      real(dp), dimension(:,:), intent(inout), optional :: f1, f2

      real(dp) :: real_space_covariance

      complex(dp) :: real_space_covariance_in, I_lm1m2
      integer :: j1, j2, n1, n2, l, m1, m2
      real(dp) :: r1, r2, arg_bess, fac_exp, mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lm, mo_spher_bess_fi_ki_lmm
      real(dp), dimension(3) :: d1, d2
      type(cplx_2d), dimension(:), allocatable :: integral_r

      logical :: do_derivative

      do_derivative = (present(f1) .or. present(f2)) 

      real_space_covariance_in = CPLX_ZERO

      allocate(integral_r(0:l_max))
      do l = 0, l_max
         allocate(integral_r(l)%mm(-l:l,-l:l))
         integral_r(l)%mm = CPLX_ZERO
      enddo

      do n1 = 1, atoms_n_neighbours(at1,i1)
         j1 = atoms_neighbour(at1,i1,n1,distance = r1, diff = d1)
         do n2 = 1, atoms_n_neighbours(at2,i2)
            j2 = atoms_neighbour(at2,i2,n2,distance = r2, diff = d2)

            arg_bess = alpha*r1*r2
            fac_exp = exp(-0.5_dp*alpha*(r1**2+r2**2))

            do l = 0, l_max
               if( l == 0 ) then
                  mo_spher_bess_fi_ki_lmm = sinh(arg_bess)/arg_bess
                  mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm
               elseif( l == 1 ) then
                  mo_spher_bess_fi_ki_lm = ( arg_bess*cosh(arg_bess) - sinh(arg_bess) ) / arg_bess**2
                  mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lm
               else
                  mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm - (2*l+1)*mo_spher_bess_fi_ki_lm / arg_bess
                  mo_spher_bess_fi_ki_lm = mo_spher_bess_fi_ki_l
                  mo_spher_bess_fi_ki_lmm = mo_spher_bess_fi_ki_lm
               endif
                  
               do m1 = -l, l
                  do m2 = -l, l
                     I_lm1m2 = conjg(SphericalYCartesian(l,m1,d1)) * SphericalYCartesian(l,m2,d2)*mo_spher_bess_fi_ki_l*fac_exp
                     integral_r(l)%mm(m2,m1) = integral_r(l)%mm(m2,m1) + I_lm1m2
                  enddo
               enddo
            enddo
         enddo
      enddo

      do l = 0, l_max
         real_space_covariance_in = real_space_covariance_in + sum(conjg(integral_r(l)%mm) * integral_r(l)%mm)
      enddo
      real_space_covariance = real(real_space_covariance_in)

      do l = 0, l_max
         deallocate(integral_r(l)%mm)
      enddo
      deallocate(integral_r)

   endfunction real_space_covariance

   function cutoff_function(r,cutoff_in)

      real(dp)             :: cutoff_function
      real(dp), intent(in) :: r, cutoff_in
      real(dp), parameter :: S = 0.25_dp

      if( r > cutoff_in ) then
          cutoff_function = 0.0_dp
      else
          cutoff_function = 0.5_dp * ( cos(PI*r/cutoff_in) + 1.0_dp )
      endif
      !if( r > cutoff_in ) then
      !    cutoff = 0.0_dp
      !elseif( r > (cutoff_in-S) ) then
      !    cutoff = 0.5_dp * ( cos(PI*(r-cutoff_in+S)/S) + 1.0_dp )
      !else
      !    cutoff = 1.0_dp
      !endif

   endfunction cutoff_function

   function dcutoff_function(r,cutoff_in)

      real(dp)             :: dcutoff_function
      real(dp), intent(in) :: r, cutoff_in
      real(dp), parameter :: S = 0.25_dp

      if( r > cutoff_in ) then
          dcutoff_function = 0.0_dp
      else
          dcutoff_function = - 0.5_dp * PI * sin(PI*r/cutoff_in) / cutoff_in
      endif
      !if( r > r_cut ) then
      !    dcutoff = 0.0_dp
      !elseif( r > (cutoff_in-S) ) then
      !    dcutoff = - 0.5_dp * PI * sin(PI*(r-cutoff_in+S)/S) / S
      !else
      !    dcutoff = 0.0_dp
      !endif

   endfunction dcutoff_function

   function coordination_function(r,cutoff_in,transition_width)

      real(dp)             :: coordination_function
      real(dp), intent(in) :: r, cutoff_in, transition_width

      if( r > cutoff_in ) then
          coordination_function = 0.0_dp
      elseif( r > (cutoff_in-transition_width) ) then
          coordination_function = 0.5_dp * ( cos(PI*(r-cutoff_in+transition_width)/transition_width) + 1.0_dp )
      else
          coordination_function = 1.0_dp
      endif

   endfunction coordination_function

   function dcoordination_function(r,cutoff_in,transition_width)

      real(dp)             :: dcoordination_function
      real(dp), intent(in) :: r, cutoff_in,transition_width

      if( r > cutoff_in ) then
          dcoordination_function = 0.0_dp
      elseif( r > (cutoff_in-transition_width) ) then
          dcoordination_function = - 0.5_dp * PI * sin(PI*(r-cutoff_in+transition_width)/transition_width) / transition_width
      else
          dcoordination_function = 0.0_dp
      endif

   endfunction dcoordination_function

   function RadialFunction(this,r,i)
      type(RadialFunction_type), intent(in) :: this
      real(dp), intent(in) :: r
      integer, intent(in) :: i

      real(dp) :: RadialFunction
   
      real(dp), dimension(this%n_max) :: h
      integer :: j
   
      if( r < this%cutoff ) then
         do j = 1, this%n_max
            h(j) = (this%cutoff-r)**(j+2) / this%NormFunction(j)
         enddo
         RadialFunction = dot_product(this%RadialTransform(:,i),h)
      else
         RadialFunction = 0.0_dp
      endif
   
   endfunction RadialFunction
   
   function GradRadialFunction(this,r,i)
      type(RadialFunction_type), intent(in) :: this
      real(dp), intent(in) :: r
      integer, intent(in) :: i

      real(dp) :: GradRadialFunction
   
      real(dp), dimension(this%n_max) :: h
      integer :: j
   
      if( r < this%cutoff ) then
         do j = 1, this%n_max
            h(j) = - (j+2) * (this%cutoff-r)**(j+1) / this%NormFunction(j)
         enddo
         GradRadialFunction = dot_product(this%RadialTransform(:,i),h)
      else
         GradRadialFunction = 0.0_dp
      endif
   
   endfunction GradRadialFunction

   !#################################################################################
   !#
   !% Solid Harmonic function using Cartesian coordinates
   !%
   !% $ R_{l m} = \sqrt{\frac{4 \pi}{2 l + 1}} r^l Y_{l m} $
   !#
   !#################################################################################

   function SolidRCartesian(l, m, x)

     complex(dp) :: SolidRCartesian
     integer, intent(in) :: l, m
     real(dp), intent(in) :: x(3)
     integer :: p, q, s

     SolidRCartesian = CPLX_ZERO

     do p = 0, l
        q = p - m
        s = l - p - q

        if ((q >= 0) .and. (s >= 0)) then
           SolidRCartesian = SolidRCartesian + ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**p) &
                                             * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**q) &
                                             * (x(3)**s) &
                                             / (factorial(p) * factorial(q) * factorial(s)))
        end if
     end do

     SolidRCartesian = SolidRCartesian * sqrt(factorial(l + m) * factorial(l - m))

   end function SolidRCartesian

   !#################################################################################
   !#
   !% Spherical Harmonic function using Cartesian coordinates
   !#
   !#################################################################################

   function SphericalYCartesian(l, m, x)

     complex(dp) :: SphericalYCartesian
     integer, intent(in) :: l, m
     real(dp), intent(in) :: x(3)

     SphericalYCartesian = SolidRCartesian(l, m, x) * sqrt(((2.0_dp * l) + 1) / (4.0_dp * PI)) &
                                                    * (normsq(x)**(-0.5_dp * l))

   end function SphericalYCartesian

    !#################################################################################
    !#
    !% Derivative of Spherical Harmonic function using Cartesian coordinates
    !#
    !#################################################################################

    function GradSphericalYCartesian(l, m, x)

      complex(dp) :: GradSphericalYCartesian(3)
      integer, intent(in) :: l, m
      real(dp), intent(in) :: x(3)
      integer :: p, q, s

      GradSphericalYCartesian = CPLX_ZERO

      do p = 0, l
         q = p - m
         s = l - p - q

         if ((p >= 1) .and. (q >= 0) .and. (s >= 0)) then
            GradSphericalYCartesian(1) = GradSphericalYCartesian(1) &
                                       - ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**(p - 1)) &
                                       * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**q) &
                                       * (x(3)**s) &
                                       * 0.5_dp &
                                       / (factorial(p - 1) * factorial(q) * factorial(s)))
            GradSphericalYCartesian(2) = GradSphericalYCartesian(2) &
                                       - ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**(p - 1)) &
                                       * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**q) &
                                       * (x(3)**s) &
                                       * 0.5_dp * cmplx(0.0_dp, 1.0_dp, dp) &
                                       / (factorial(p - 1) * factorial(q) * factorial(s)))
         end if

         if ((p >= 0) .and. (q >= 1) .and. (s >= 0)) then
            GradSphericalYCartesian(1) = GradSphericalYCartesian(1) &
                                       + ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**p) &
                                       * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**(q - 1)) &
                                       * (x(3)**s) &
                                       * 0.5_dp &
                                       / (factorial(p) * factorial(q - 1) * factorial(s)))
            GradSphericalYCartesian(2) = GradSphericalYCartesian(2) &
                                       - ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**p) &
                                       * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**(q - 1)) &
                                       * (x(3)**s) &
                                       * 0.5_dp * cmplx(0.0_dp, 1.0_dp, dp) &
                                       / (factorial(p) * factorial(q - 1) * factorial(s)))
         end if

         if ((p >= 0) .and. (q >= 0) .and. (s >= 1)) then
            GradSphericalYCartesian(3) = GradSphericalYCartesian(3) &
                                       + ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**p) &
                                       * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**q) &
                                       * (x(3)**(s - 1)) &
                                       / (factorial(p) * factorial(q) * factorial(s - 1)))
         end if
      end do

      GradSphericalYCartesian = GradSphericalYCartesian &
                              * sqrt(factorial(l + m) * factorial(l - m) * ((2.0_dp * l) + 1) / (4.0_dp * PI)) &
                              * (normsq(x)**(-0.5_dp * l))

      GradSphericalYCartesian = GradSphericalYCartesian &
                              - (l * x * SphericalYCartesian(l, m, x) / normsq(x))

    end function GradSphericalYCartesian
   subroutine cg_initialise(j,denom)

      integer, intent(in) :: j
      integer :: i_j1,i_m1,i_j2,i_m2,i_j,i_m
      integer, intent(in), optional :: denom

      integer :: my_denom

      if(cg_initialised) return

      my_denom = optional_default(1,denom)

      cg_j1_max = j
      cg_m1_max = j
      cg_j2_max = j
      cg_m2_max = j
      cg_j_max = j !(j1_max+j2_max)
      cg_m_max = j !(j1_max+j2_max)

      allocate( cg_array(0:cg_j1_max,-cg_m1_max:cg_m1_max,0:cg_j2_max,-cg_m2_max:cg_m2_max,&
      & 0:cg_j_max,-cg_j_max:cg_j_max) )
 
      cg_array = 0.0_dp

      do i_j1 = 0, cg_j1_max
      do i_m1 = -i_j1, i_j1, my_denom
      do i_j2 = 0, cg_j2_max
      do i_m2 = -i_j2, i_j2, my_denom
      do i_j = abs(i_j1-i_j2), min(cg_j_max,i_j1+i_j2)
      do i_m = -i_j, i_j, my_denom


         cg_array(i_j1,i_m1,i_j2,i_m2,i_j,i_m) = cg_calculate(i_j1,i_m1,i_j2,i_m2,i_j,i_m,denom)

      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      cg_initialised = .true.
    
   endsubroutine cg_initialise

   !#################################################################################
   !#
   !% Finalise global CG arrays
   !#
   !#################################################################################

   subroutine cg_finalise

      cg_j1_max = 0
      cg_m1_max = 0
      cg_j2_max = 0
      cg_m2_max = 0
      cg_j_max = 0
      cg_m_max = 0

      if(allocated(cg_array)) deallocate( cg_array )
      cg_initialised = .false.

   endsubroutine cg_finalise

   !#################################################################################
   !#
   !% Look up CG coefficient from CG array, previously calculated.
   !#
   !#################################################################################

   function cg_lookup(j1,m1,j2,m2,j,m,denom) result(cg)

     real(dp)            :: cg
     integer, intent(in) :: j1,m1,j2,m2,j,m
     integer, intent(in), optional :: denom

     cg=0.0_dp

     if ( .not. cg_check(j1,m1,j2,m2,j,m,denom) ) then
        return
     endif

     if( j1<=cg_j1_max .and. j2<=cg_j2_max .and. j<=cg_j_max .and. &
       abs(m1)<=cg_m1_max .and. abs(m2)<=cg_m2_max .and. abs(m) <= cg_m_max .and. cg_initialised ) then
         cg = cg_array(j1,m1,j2,m2,j,m)
     else
         cg = cg_calculate(j1,m1,j2,m2,j,m,denom)
     endif

   endfunction cg_lookup

   !#################################################################################
   !#
   !% Check if input variables for CG make sense.
   !% Source: http://mathworld.wolfram.com/Clebsch-GordanCoefficient.html \\
   !%                                                                    
   !% $ j_1 + j_2 \ge j $ \\
   !% $ j_1 - j_2 \ge -j $ \\
   !% $ j_1 - j_2 \le j $ \\
   !% $ j_1 \ge m_1 \ge j_1 $ \\
   !% $ j_2 \ge m_2 \ge j_2 $ \\
   !% $ j \ge m \ge j $ \\
   !% $ m_1 + m_2 = m $ \\
   !#
   !#################################################################################

   function cg_check(j1,m1,j2,m2,j,m,denom)

      logical :: cg_check
      integer, intent(in) :: j1,m1,j2,m2,j,m
      integer, intent(in), optional :: denom

      integer :: my_denom

      my_denom = optional_default(1,denom)
      cg_check = (j1>=0) .and. (j2>=0) .and. (j>=0) .and. &
               (abs(m1)<=j1) .and. (abs(m2)<=j2) .and. (abs(m)<=j) &
               .and. (m1+m2==m) .and. (j1+j2 >= j) .and. (abs(j1-j2) <= j) &
               .and. (mod(j1+j2+j,my_denom)==0)
      
   endfunction cg_check

   !#################################################################################
   !#
   !% Calculate a Clebsch-Gordan coefficient $\left< j_1 m_1 j_2 m_2 | j m \right>$
   !% Source: http://mathworld.wolfram.com/Clebsch-GordanCoefficient.html \\
   !% $ \left< j_1 m_1 j_2 m_2 | j m \right> = (-1)^{m+j_1-j_2) 
   !% \sqrt{2j+1} \left( \begin{array}{ccc}
   !% j_1 & j_2 & j \\
   !% m_1 & m_2 & -m \\
   !% \end{array} \right) $
   !% where the thing on the right-hand side is the Wigner 3J symbol.
   !#
   !#################################################################################

   function cg_calculate(j1,m1,j2,m2,j,m,denom) result(cg)

     real(dp)            :: cg
     integer, intent(in) :: j1,m1,j2,m2,j,m
     integer, intent(in), optional :: denom

     integer :: my_denom

     my_denom = optional_default(1,denom)

     cg=0.0_dp
     if ( .not. cg_check(j1,m1,j2,m2,j,m,denom) ) return

     cg = oscillate((m+j1-j2)/my_denom) * sqrt(2.0_dp*real(j,dp)/real(my_denom,dp)+1.0_dp) * &
     wigner3j(j1,m1,j2,m2,j,-m,denom)

   end function cg_calculate

   !#################################################################################
   !#
   !% Triangle coefficient
   !% Source: http://mathworld.wolfram.com/TriangleCoefficient.html
   !% 
   !% $ \Delta(a,b,c) = \frac{ (a+b-c)! (a-b+c)! (-a+b+c)! }{ (a+b+c+1)! } $
   !#
   !#################################################################################

   function tc(a,b,c,denom)

      real(dp) :: tc
      integer, intent(in) :: a, b, c
      integer, intent(in), optional :: denom
      integer :: my_denom

      my_denom = optional_default(1,denom)

      tc = factorial((a+b-c)/my_denom) * factorial((a-b+c)/my_denom) * &
      factorial((-a+b+c)/my_denom) / factorial((a+b+c)/my_denom+1)

   endfunction tc

   !#################################################################################
   !#
   !% Wigner 3J symbol
   !% Source: http://mathworld.wolfram.com/Wigner3j-Symbol.html
   !% 
   !% \[
   !% \left( \begin{array}{ccc}
   !% j_1 & j_2 & j \\
   !% m_1 & m_2 & m \\
   !% \end{array} \right) = (-1)^{j_1-j_2-m) \sqrt{ \Delta (j_1,j_2,j) }
   !% \sqrt{ (j_1+m_1)! (j_1-m_1)! (j_2+m_2)! (j_2-m_2)! (j+m)! (j-m)! }
   !% \sum_k \frac{ (-1)^k }{k! (j-j_2+k+m_1)! (j-j_1+k-m_2)! (j_1+j_2-j-k)!
   !% (j_1-k-m_1)! (j_2-k+m_2)! }
   !% \]
   !% the summation index k runs on all integers where none of the argument of
   !% factorials are negative
   !% $\Delta(a,b,c)$ is the triangle coefficient.
   !#
   !#################################################################################

   function wigner3j(j1,m1,j2,m2,j,m,denom)

       real(dp)            :: wigner3j
       integer, intent(in) :: j1,m1,j2,m2,j,m
       integer, intent(in), optional :: denom

       real(dp) :: pre_fac, triang_coeff, main_coeff, sum_coeff, sum_term
       integer  :: k, kmin, kmax, my_denom

       my_denom = optional_default(1,denom)

       pre_fac = oscillate((j1-j2-m)/my_denom)

       triang_coeff = sqrt( tc(j1,j2,j,denom) )

       main_coeff = sqrt( &
                  factorial((j1+m1)/my_denom) * factorial((j1-m1)/my_denom) * &
                  factorial((j2+m2)/my_denom) * factorial((j2-m2)/my_denom) * &
                  factorial((j+m)/my_denom) * factorial((j-m)/my_denom) )
                  
       sum_coeff = 0.0_dp

       kmin = max( j2-j-m1, j1+m2-j, 0 ) / my_denom
       kmax = min( j1+j2-j, j1-m1, j2+m2) / my_denom

       do k = kmin, kmax

          sum_term = 1.0_dp / ( factorial(k) * factorial((j-j2+m1)/my_denom+k) * &
                   factorial((j-j1-m2)/my_denom+k) * factorial((j1+j2-j)/my_denom-k) * &
                   factorial((j1-m1)/my_denom-k) * factorial((j2+m2)/my_denom-k) )

          sum_term = oscillate(k) * sum_term
          
          sum_coeff = sum_coeff + sum_term

       enddo

       wigner3j = pre_fac * triang_coeff * main_coeff * sum_coeff

   endfunction wigner3j

   !#################################################################################
   !#
   !% Factorial, real result
   !#
   !#################################################################################

    function factorial(n) result(res)

     ! factorial_real

     integer, intent(in) :: n
     real(dp)            :: res
     integer :: i

     if (n<0) then
        call system_abort('factorial: negative argument')
     elseif(n <= 16) then
        res = factorial_table(n)
     else
        res=1.0_dp
        do i=2,n
           res = res*i
        end do
     end if

   endfunction factorial

   !#################################################################################
   !#
   !% Factorial, integer result
   !#
   !#################################################################################

   function factorial_int(n) result(res)

     ! factorial_int

     integer, intent(in) :: n
     integer             :: res
     integer :: i

     if (n<0) then
        call system_abort('factorial_int: negative argument')
     else
        res=1
        do i=2,n
           res = res*i
        end do
     end if

   endfunction factorial_int

   !#################################################################################
   !#
   !% Double factorial, real result
   !#
   !#################################################################################

   recursive function factorial2(n) result(res)

     ! double factorial

     integer, intent(in) :: n
     real(dp)            :: res
     integer :: i

     if( n < -1 ) then
         call system_abort('factorial2: negative argument')
     else
         res = 1.0_dp
         do i = 2-mod(n,2), n, 2
            res = res*i
         enddo
     endif

   endfunction factorial2

   !#################################################################################
   !#
   !% Binomial coefficient, real result
   !#
   !#################################################################################

   recursive function binom(n,r) result(res)

     ! binomial coefficients

     integer, intent(in) :: n,r
     real(dp)            :: res
     integer             :: i

     if((n<0) .or. (r<0) .or. (n<r)) then
        res = 0.0_dp
     else
        res = 1.0_dp
        do i = 0, r-1
           res = real(n-i,dp)/real(r-i,dp) * res
        enddo
     endif

   endfunction binom

   !#################################################################################
   !#
   !% $ (-1)^n $ function.
   !#
   !#################################################################################

   function oscillate(m)

     integer, intent(in) :: m
     integer :: oscillate

     if( mod(m,2) == 0 ) then
         oscillate = 1
     else
         oscillate = -1
     endif

   endfunction oscillate

endmodule descriptors_module
