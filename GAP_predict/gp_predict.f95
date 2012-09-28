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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X Gaussian Process module
!X
!% Module for general GP function interpolations.
!% A gp object contains the training set (teaching points and function values),
!% important temporary matrices, vectors and parameters.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module gp_predict_module

   use iso_c_binding, only : C_NULL_CHAR
   use libatoms_module
   use descriptors_module
   use fox_wxml
   use FoX_sax, only: xml_t, dictionary_t, haskey, getvalue, parse, &
   open_xml_string, close_xml_t

   implicit none

   private

   integer, parameter :: besseli_max_n = 20

   real(dp), dimension(besseli_max_n), parameter :: besseli0_c = (/ &
   0.125_dp, &
   7.03125E-002_dp, &
   7.32421875E-002_dp, &
   0.112152099609375_dp, &    
   0.22710800170898438_dp, &    
   0.57250142097473145_dp, &    
   1.7277275025844574_dp, &    
   6.0740420012734830_dp, &    
   24.380529699556064_dp, &    
   110.01714026924674_dp, &    
   551.33589612202059_dp, &    
   3038.0905109223841_dp, &    
   18257.755474293175_dp, &    
   118838.42625678326_dp, &    
   832859.30401628942_dp, &    
   6252951.4934347980_dp, &    
   50069589.531988934_dp, &    
   425939216.50476694_dp, &    
   3836255180.2304339_dp, &    
   36468400807.065559_dp /)    

   real(dp), dimension(besseli_max_n), parameter :: besseli1_c = (/ &
   -0.375_dp, &    
   -0.1171875_dp, &    
   -0.1025390625_dp, &    
   -0.144195556640625_dp, &    
   -0.27757644653320313_dp, &    
   -0.67659258842468262_dp, &    
   -1.9935317337512970_dp, &    
   -6.8839142681099474_dp, &    
   -27.248827311268542_dp, &    
   -121.59789187653587_dp, &    
   -603.84407670507017_dp, &    
   -3302.2722944808525_dp, &    
   -19718.375912236628_dp, &    
   -127641.27264617461_dp, &    
   -890297.87670706783_dp, &    
   -6656367.7188176867_dp, &    
   -53104110.109685220_dp, &    
   -450278600.30503929_dp, &    
   -4043620325.1077542_dp, &    
   -38338575207.427895_dp /)

   real(dp), parameter :: besseli_max_x = 18.0_dp

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

   real(dp), parameter :: THETA_MIN = 1.0e-8_dp
   integer, parameter, public :: GP_SPARSE_RANDOM = 1
   integer, parameter, public :: GP_SPARSE_PIVOT = 2
   integer, parameter, public :: GP_SPARSE_CLUSTER = 3
   integer, parameter, public :: GP_SPARSE_UNIFORM = 4
   integer, parameter, public :: GP_SPARSE_KMEANS = 5
   integer, parameter, public :: GP_SPARSE_COVARIANCE = 6

   integer, parameter, public :: GP_COVARIANCE_FITC = 1
   integer, parameter, public :: GP_COVARIANCE_DTC = 2
   type gpCovariance_bond_real_space

      integer :: n
      real(dp) :: delta
      real(dp) :: atom_sigma

      integer, dimension(:), allocatable :: sparseX_N
      real(dp), dimension(:,:,:,:), allocatable :: sparseX_pos
      real(dp), dimension(:,:), allocatable :: sparseX_atom_covariance_cutoff
      real(dp), dimension(:), allocatable :: sparseX_bond_covariance_cutoff
      real(dp), dimension(:), allocatable :: sparseX_self_overlap
      logical :: sparseX_inisialised = .false.

      logical :: initialised = .false.

   endtype gpCovariance_bond_real_space

!   type gpCovariance_bond_real_space
!
!      real(dp) :: theta
!      real(dp) :: delta
!      integer :: integration_steps
!      real(dp) :: atom_sigma
!
!      real(dp), dimension(:,:,:,:), allocatable :: rot_matrix
!      integer, dimension(:), allocatable :: sparseX_N
!      real(dp), dimension(:,:,:,:,:), allocatable :: sparseX_pos_rot
!      real(dp), dimension(:,:), allocatable :: sparseX_atom_covariance_cutoff
!      real(dp), dimension(:), allocatable :: sparseX_bond_covariance_cutoff
!      real(dp), dimension(:), allocatable :: sparseX_invariant_buffer
!      logical :: sparseX_inisialised = .false.
!
!      logical :: initialised = .false.
!
!   endtype gpCovariance_bond_real_space

   type gpCovariance_atom_real_space

      integer :: l_max = 0
      real(dp) :: atom_sigma, delta, zeta
      logical :: sparseX_inisialised = .false.
      real(dp) :: cutoff, cutoff_transition_width

      logical :: initialised = .false.

   endtype gpCovariance_atom_real_space

   public :: gpCovariance_atom_real_space
   public :: gpcovariance_atom_real_space_calc

   type gpCoordinates

      integer :: d = 0, n_x, n_xPrime, n_sparseX, n_permutations
      ! dimension of descriptors, number of descriptors, number of derivatives of descriptors

      integer :: current_x, current_xPrime
      ! pointers to the last added values

      real(dp), dimension(:,:), allocatable :: x, xPrime
      ! descriptors (d,n_x), derivatives of descriptors (d, n_xPrime)
      ! for real space covariance descriptors (max(x_size),n_x), derivatives of descriptors (max(x_size),n_xPrime)
      real(dp), dimension(:), allocatable :: cutoff, cutoffPrime
      integer, dimension(:), allocatable :: x_size, xPrime_size
      real(dp), dimension(:), allocatable :: covarianceDiag_x_x, covarianceDiag_xPrime_xPrime

      real(dp), dimension(:,:), allocatable :: sparseX, covarianceDiag_x_xPrime
      real(dp), dimension(:), allocatable :: sparseCutoff
      ! sparse points stored as real array
      ! for real space covariance descriptors
      integer, dimension(:), allocatable :: sparseX_size
      real(dp), dimension(:), allocatable :: covarianceDiag_sparseX_sparseX

      real(dp), dimension(:), allocatable :: theta
      ! range parameters (d) for descriptors in each directions

      real(dp), dimension(:), allocatable :: alpha
      ! 

      real(dp) ::  delta, f0 = 0.0_dp
      ! range of GP (function value) and baseline of function

      integer, dimension(:), allocatable :: map_x_y, map_xPrime_yPrime, map_xPrime_x, config_type
      ! which descriptor is used for a given function value, which derivative descriptor is used for a given derivative function, which descriptor is differentiated

      integer, dimension(:), allocatable :: map_sparseX_globalSparseX
      ! sparse point in this descriptor type -> all sparse points in gpFull

      integer, dimension(:), allocatable :: sparseX_index
      ! sparse points stored as indices of the x array

      integer, dimension(:,:), allocatable :: permutations
      ! Lists the permutations symmetries of the coordinates

      type(gpCovariance_bond_real_space) :: bond_real_space_cov
      logical :: bond_real_space = .false.

      type(gpCovariance_atom_real_space) :: atom_real_space_cov
      logical :: atom_real_space = .false.

      logical :: soap = .false.

      type(extendable_str) :: descriptor_str

      type(LA_Matrix) :: LA_k_mm

      logical :: initialised = .false.
      logical :: sparsified = .false.
      logical :: sparseScore_initialised = .false.

   endtype gpCoordinates

   public :: gpCoordinates 

   type gpFull

      integer :: n_y, n_yPrime
      ! number of function values, number of derivative function values
      
      integer :: n_globalSparseX
      ! number of all sparse points in every descriptor type

      integer :: n_coordinate
      ! number of different descriptors

      integer :: current_y, current_yPrime

      real(dp) :: sparse_jitter = 1.0e-5_dp

      real(dp), dimension(:), allocatable :: y, yPrime
      ! function values, derivative function values

      real(dp), dimension(:), allocatable :: sigma_y, sigma_yPrime
      ! estimated error of function values, derivatives

      real(dp), dimension(:,:), allocatable :: covariance_subY_y, covariance_subY_subY, covariance_y_y, inverse_sparse_full
      ! covariance matrix

      real(dp), dimension(:), allocatable :: covarianceDiag_y_y, lambda, alpha
      ! covariance matrix

      integer, dimension(:), allocatable :: map_y_globalY, map_yPrime_globalY

      type(gpCoordinates), dimension(:), allocatable :: coordinate

      integer :: covariance_method = GP_COVARIANCE_DTC
      logical :: initialised = .false.

   endtype gpFull

   type gpSparse

      integer :: n_coordinate
      ! number of different descriptors

      type(gpCoordinates), dimension(:), allocatable :: coordinate

      logical :: initialised = .false.

   endtype gpSparse

   type cplx_1d_array
      complex(dp), dimension(:), allocatable :: value
   endtype cplx_1d_array

   type cplx_2d_array
      complex(dp), dimension(:,:), allocatable :: value
   endtype cplx_2d_array

   type neighbour_descriptor
      type(cplx_1d_array), dimension(:), allocatable :: spherical_harmonics
      real(dp) :: r
      integer :: n
   endtype neighbour_descriptor

   logical, save :: parse_matched_label, parse_in_gpCoordinates, parse_in_gpFull, parse_in_gpSparse, parse_in_sparseX, parse_sliced, parse_sparseX_separate_file
   integer, save :: parse_i_sparseX, parse_i_x, parse_i_xPrime, parse_i_permutation, parse_slice_start, parse_slice_end
   type(gpCoordinates), pointer :: parse_gpCoordinates
   type(gpFull), pointer :: parse_gpFull
   type(gpSparse), pointer :: parse_gpSparse
   type(extendable_str), save :: parse_cur_data
   integer, dimension(:,:), allocatable :: parse_in_permutations
   character(len=1024), save :: parse_gpCoordinates_label, parse_gpFull_label, parse_gpSparse_label

   public :: gpFull, gpSparse

   interface initialise
      module procedure gpSparse_initialise
   endinterface initialise
   public :: initialise

   interface finalise
      module procedure gpFull_Finalise, gpCoordinates_Finalise, gpSparse_finalise, gpNeighbourDescriptor_Finalise
   endinterface finalise
   public :: finalise

   interface gp_setTheta
      module procedure gpCoordinates_setTheta, gpFull_setTheta
   endinterface gp_setTheta
   public :: gp_setTheta

   interface gp_setThetaFactor
      module procedure gpFull_setTheta_thetaFactor !, gpFull_setTheta_thetaFactorArray, gpFull_setTheta_thetaFactorUniform
   endinterface gp_setThetaFactor
   public :: gp_setThetaFactor

   interface gp_setParameters
      module procedure gpFull_setParameters, gpFull_gpCoordinates_setParameters, gpCoordinates_setParameters, &
      gpCoordinates_setParameters_sparse, gpSparse_setParameters
   endinterface gp_setParameters
   public :: gp_setParameters

   interface gp_setPermutations
      module procedure gpCoordinates_setPermutations, gpFull_setPermutations, gpSparse_setPermutations
   endinterface gp_setPermutations
   public :: gp_setPermutations

   interface gp_addFunctionValue
      module procedure gpFull_addFunctionValue
   endinterface gp_addFunctionValue
   public :: gp_addFunctionValue

   interface gp_addFunctionDerivative
      module procedure gpFull_addFunctionDerivative
   endinterface gp_addFunctionDerivative
   public :: gp_addFunctionDerivative

   interface gp_addCoordinates
      module procedure gpFull_addCoordinates_1Darray, gpFull_addCoordinates_2Darray
   endinterface gp_addCoordinates
   public :: gp_addCoordinates

   interface gp_addCoordinateDerivatives
      module procedure gpFull_addCoordinateDerivatives_1Darray, gpFull_addCoordinateDerivatives_2Darray
   endinterface gp_addCoordinateDerivatives
   public :: gp_addCoordinateDerivatives

   interface gp_addDescriptor
      module procedure gpFull_addDescriptor
   endinterface gp_addDescriptor
   public :: gp_addDescriptor

   interface gp_printXML
      module procedure gpCoordinates_printXML, gpFull_printXML, gpSparse_printXML
   endinterface gp_PrintXML
   public :: gp_printXML

   interface gp_readXML
      module procedure gpCoordinates_readXML, gpFull_readXML, gpSparse_readXML, &
      gpCoordinates_readXML_string, gpFull_readXML_string, gpSparse_readXML_string
   endinterface gp_readXML
   public :: gp_readXML

   interface gp_covariance_sparse
      module procedure gpFull_covarianceMatrix_sparse
   endinterface gp_covariance_sparse
   public :: gp_covariance_sparse

   interface gp_covariance_full
      module procedure gpFull_covarianceMatrix
   endinterface gp_covariance_full
   public :: gp_covariance_full

   interface gp_Predict
      module procedure gpCoordinates_Predict
   endinterface gp_Predict
   public :: gp_Predict

   public :: gpCoordinates_Covariance, gpCoordinates_gpCovariance_bond_real_space_Initialise, gpCoordinates_gpCovariance_atom_real_space_Initialise
   public :: gpCoordinates_initialise_SparseScore

   contains

   subroutine gpFull_setParameters(this, n_coordinate, n_y, n_yPrime, sparse_jitter, error)

      type(gpFull), intent(inout) :: this
      integer, intent(in) :: n_coordinate, n_y, n_yPrime
      real(dp), intent(in) :: sparse_jitter
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(this%initialised) call finalise(this,error)

      this%n_coordinate = n_coordinate
      this%n_y = n_y
      this%n_yPrime = n_yPrime
      this%current_y = 0
      this%current_yPrime = 0
      this%sparse_jitter = sparse_jitter

      allocate( this%coordinate(n_coordinate) )
      allocate( this%y(n_y), this%yPrime(n_yPrime) )
      allocate( this%map_y_globalY(n_y), this%map_yPrime_globalY(n_yPrime) )
      allocate( this%sigma_y(n_y), this%sigma_yPrime(n_yPrime) )

      this%initialised = .true.

   endsubroutine gpFull_setParameters

   subroutine gpFull_gpCoordinates_setParameters(this, i, d, n_x, n_xPrime, delta, f0, bond_real_space, atom_real_space, soap, x_size_max, xPrime_size_max, error)

      type(gpFull), intent(inout) :: this
      integer, intent(in) :: i, d, n_x, n_xPrime
      real(dp), intent(in) :: delta, f0
      logical, optional, intent(in) :: bond_real_space, atom_real_space, soap
      integer, optional, intent(in) :: x_size_max, xPrime_size_max
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_set_gpCoordinates_parameters: object not initialised',error)
      endif

      if( i > this%n_coordinate ) then
         RAISE_ERROR( 'gpFull_set_gpCoordinates_parameters: access to descriptor '//i//' is not possible as number of descriptors is '//this%n_coordinate,error )
      endif

      call gpCoordinates_setParameters(this%coordinate(i), d, n_x, n_xPrime, delta, f0, bond_real_space=bond_real_space, atom_real_space=atom_real_space, soap = soap, x_size_max=x_size_max, xPrime_size_max=xPrime_size_max, error=error)

   endsubroutine gpFull_gpCoordinates_setParameters

   subroutine gpCoordinates_setParameters(this, d, n_x, n_xPrime, delta, f0, bond_real_space, atom_real_space, soap, x_size_max, xPrime_size_max, error)

      type(gpCoordinates), intent(inout) :: this
      integer, intent(in) :: d, n_x, n_xPrime
      real(dp), intent(in) :: delta, f0
      logical, optional, intent(in) :: bond_real_space, atom_real_space, soap
      integer, optional, intent(in) :: x_size_max, xPrime_size_max
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(this%initialised) call finalise(this,error)

      this%d = d
      this%n_x = n_x
      this%n_xPrime = n_xPrime
      this%delta = delta
      this%f0 = f0

      this%current_x = 0
      this%current_xPrime = 0
      this%n_sparseX = 0
      this%n_permutations = 1

      this%bond_real_space = optional_default(.false., bond_real_space)
      this%atom_real_space = optional_default(.false., atom_real_space)
      this%soap = optional_default(.false., soap)

      if(present(x_size_max)) then
         allocate( this%x(x_size_max,n_x) )
      else
         allocate( this%x(d,n_x) )
      endif
      this%x = 0.0_dp

      if(present(xPrime_size_max) .and. this%atom_real_space) then
         allocate( this%xPrime(xPrime_size_max,n_xPrime) )
      else
         allocate( this%xPrime(d,n_xPrime) )
      endif
      this%xPrime = 0.0_dp

      allocate(this%cutoff(n_x))
      this%cutoff = 1.0_dp
      allocate(this%cutoffPrime(n_xPrime))
      this%cutoffPrime = 0.0_dp

      allocate( this%config_type(n_x) )
      this%config_type = 0

      allocate( this%map_x_y(n_x), this%map_xPrime_yPrime(n_xPrime), this%map_xPrime_x(n_xPrime) )
      this%map_x_y = 0
      this%map_xPrime_yPrime = 0
      this%map_xPrime_x = 0

      allocate(this%covarianceDiag_x_x(n_x), this%covarianceDiag_xPrime_xPrime(n_xPrime))
      this%covarianceDiag_x_x = 1.0_dp
      this%covarianceDiag_xPrime_xPrime = 1.0_dp

      if(this%bond_real_space .or. this%atom_real_space) then
         allocate( this%x_size(n_x), this%xPrime_size(n_xPrime) )
         this%x_size = d
         this%xPrime_size = 0
         allocate( this%theta(1), this%permutations(1,1) )
         this%theta = 0.0_dp
         this%permutations = 1
      elseif( this%soap ) then
         allocate( this%theta(1), this%permutations(1,1) )
         this%theta = 1.0_dp
         this%permutations = 1
      else
         allocate( this%theta(d), this%permutations(d,1) )
         this%theta = 0.0_dp
         this%permutations(:,1) = (/ (i, i=1, d) /)
      endif

      this%sparsified = .false.
      this%initialised = .true.
   endsubroutine gpCoordinates_setParameters

   subroutine gpCoordinates_setParameters_sparse(this, d, n_sparseX, delta, f0, bond_real_space, atom_real_space, soap, sparseX_size_max, error)

      type(gpCoordinates), intent(inout) :: this
      integer, intent(in) :: d, n_sparseX
      real(dp), intent(in) :: delta, f0
      logical, optional, intent(in) :: bond_real_space, atom_real_space, soap
      integer, optional, intent(in) :: sparseX_size_max
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(this%initialised) call finalise(this,error)

      this%d = d
      this%n_x = 0
      this%n_xPrime = 0
      this%delta = delta
      this%f0 = f0

      this%current_x = 0
      this%current_xPrime = 0
      this%n_sparseX = n_sparseX
      this%n_permutations = 1

      this%bond_real_space = optional_default(.false., bond_real_space)
      this%atom_real_space = optional_default(.false., atom_real_space)
      this%soap = optional_default(.false., soap)

      if(present(sparseX_size_max)) then
         allocate( this%sparseX(sparseX_size_max,n_sparseX) )
      else
         allocate( this%sparseX(d,n_sparseX) )
      endif

      allocate( this%alpha(n_sparseX) )
      allocate( this%sparseCutoff(n_sparseX) )

      this%sparsified = .true.
      this%initialised = .true.

      allocate( this%covarianceDiag_sparseX_sparseX(n_sparseX) )
      this%covarianceDiag_sparseX_sparseX = 1.0_dp

      if(this%bond_real_space .or. this%atom_real_space) then
         allocate( this%sparseX_size(n_sparseX) )
         this%sparseX_size = 0
         allocate( this%theta(1), this%permutations(1,1) )
         this%theta = 0.0_dp
         this%permutations = 1
      elseif( this%soap ) then
         allocate( this%theta(1), this%permutations(1,1) )
         this%theta = 1.0_dp
         this%permutations = 1
      else
         allocate( this%theta(d), this%permutations(d,1) )
         this%theta = 0.0_dp
         this%permutations(:,1) = (/ (i, i=1, d) /)
      endif

   endsubroutine gpCoordinates_setParameters_sparse

   subroutine gpSparse_setParameters(this,n_coordinate,error)
      type(gpSparse), intent(inout) :: this
      integer, intent(in) :: n_coordinate
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(this%initialised) call finalise(this,error)
      this%n_coordinate = n_coordinate
      allocate( this%coordinate(this%n_coordinate) )

   endsubroutine gpSparse_setParameters

   subroutine gpCoordinates_setPermutations(this,permutations,error)
      type(gpCoordinates), intent(inout) :: this
      integer, dimension(:,:), intent(in) :: permutations
      integer, optional, intent(out) :: error

      real(dp), dimension(this%d) :: theta
      integer :: i

      INIT_ERROR(error)

      this%n_permutations = size(permutations,2)

      if(this%atom_real_space .or. this%bond_real_space .or. this%soap) then

      else
         call reallocate(this%permutations,this%d,this%n_permutations,zero=.true.)
         this%permutations = permutations
         ! Symmetrise theta wrt permutations
         theta = this%theta
         this%theta = 0.0_dp
         do i = 1, this%n_permutations
            this%theta = this%theta + theta(this%permutations(:,i))
         enddo
         this%theta = this%theta / real(this%n_permutations,kind=dp)
      endif


   endsubroutine gpCoordinates_setPermutations

   subroutine gpFull_setPermutations(this,i_coordinate,permutations,error)
      type(gpFull), intent(inout) :: this
      integer :: i_coordinate
      integer, dimension(:,:), intent(in) :: permutations
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( i_coordinate > this%n_coordinate ) then
         RAISE_ERROR( 'gpFull_setPermutations: access to descriptor '//i_coordinate//' is not possible as number of descriptors is set '//this%n_coordinate,error )
      endif

      call gpCoordinates_setPermutations(this%coordinate(i_coordinate),permutations,error)

   endsubroutine gpFull_setPermutations

   subroutine gpSparse_setPermutations(this,i_coordinate,permutations,error)
      type(gpSparse), intent(inout) :: this
      integer :: i_coordinate
      integer, dimension(:,:), intent(in) :: permutations
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( i_coordinate > this%n_coordinate ) then
         RAISE_ERROR( 'gpSparse_setPermutations: access to descriptor '//i_coordinate//' is not possible as number of descriptors is set '//this%n_coordinate,error )
      endif

      call gpCoordinates_setPermutations(this%coordinate(i_coordinate),permutations,error)

   endsubroutine gpSparse_setPermutations

   subroutine gpSparse_initialise(this, from, error)
      type(gpSparse), intent(inout) :: this
      type(gpFull), intent(in)  :: from
      integer, optional, intent(out) :: error

      integer :: i_coordinate, i_sparseX, i_global_sparseX, n_globalSparseX, n_globalY, i, j, i_y, i_yPrime, &
      i_globalY, i_global_yPrime
#ifdef HAVE_QR      
      real(qp), dimension(:,:), allocatable :: c_subYY_sqrtInverseLambda, factor_c_subYsubY, a
      real(qp), dimension(:), allocatable :: globalY, alpha
      type(LA_Matrix) :: LA_c_subYsubY, LA_q_subYsubY
#else
      real(qp), dimension(:,:), allocatable :: c_subYY_inverseLambda, c_subYY_inverseLambda_c_YsubY!, &
!      inverse_q_subYsubY, inverse_c_subYsubY
      real(qp), dimension(:), allocatable :: globalY, alpha
      type(LA_Matrix) :: LA_q_subYsubY
#endif

      INIT_ERROR(error)

      if( .not. from%initialised ) then
         RAISE_ERROR('gpSparse_initialise: gpFull object not initialised',error)
      endif

      if(this%initialised) call finalise(this,error)

      call gpSparse_setParameters(this, from%n_coordinate)

      do i_coordinate = 1, this%n_coordinate

         if(from%coordinate(i_coordinate)%bond_real_space .or. from%coordinate(i_coordinate)%atom_real_space) then
            call gpCoordinates_setParameters_sparse(this%coordinate(i_coordinate), &
                 from%coordinate(i_coordinate)%d, from%coordinate(i_coordinate)%n_sparseX, from%coordinate(i_coordinate)%delta, from%coordinate(i_coordinate)%f0, &
                 bond_real_space=from%coordinate(i_coordinate)%bond_real_space, atom_real_space=from%coordinate(i_coordinate)%atom_real_space, &
                 sparseX_size_max=maxval(from%coordinate(i_coordinate)%sparseX_size), &
                 error=error)
         else
            call gpCoordinates_setParameters_sparse(this%coordinate(i_coordinate), &
                 from%coordinate(i_coordinate)%d, from%coordinate(i_coordinate)%n_sparseX, from%coordinate(i_coordinate)%delta, from%coordinate(i_coordinate)%f0, soap = from%coordinate(i_coordinate)%soap, &
                 error=error)
         endif

         this%coordinate(i_coordinate)%sparseX = from%coordinate(i_coordinate)%sparseX
         this%coordinate(i_coordinate)%covarianceDiag_sparseX_sparseX = from%coordinate(i_coordinate)%covarianceDiag_sparseX_sparseX
         this%coordinate(i_coordinate)%atom_real_space_cov = from%coordinate(i_coordinate)%atom_real_space_cov

         if(from%coordinate(i_coordinate)%bond_real_space .or. from%coordinate(i_coordinate)%atom_real_space) then 
            this%coordinate(i_coordinate)%sparseX_size = from%coordinate(i_coordinate)%sparseX_size
         endif
         this%coordinate(i_coordinate)%theta = from%coordinate(i_coordinate)%theta
         this%coordinate(i_coordinate)%descriptor_str = from%coordinate(i_coordinate)%descriptor_str
         this%coordinate(i_coordinate)%sparseCutoff = from%coordinate(i_coordinate)%sparseCutoff

         call gpSparse_setPermutations(this,i_coordinate,from%coordinate(i_coordinate)%permutations,error)

      enddo
      
      n_globalSparseX = from%n_globalSparseX
      n_globalY = from%n_y + from%n_yPrime

#ifdef HAVE_QR
      allocate( c_subYY_sqrtInverseLambda(n_globalSparseX,n_globalY), factor_c_subYsubY(n_globalSparseX,n_globalSparseX), &
      a(n_globalY+n_globalSparseX,n_globalSparseX), globalY(n_globalY+n_globalSparseX), alpha(n_globalSparseX) )

      call matrix_product_vect_asdiagonal_sub(c_subYY_sqrtInverseLambda,from%covariance_subY_y,sqrt(1.0_qp/from%lambda)) ! O(NM)
      call initialise(LA_c_subYsubY,from%covariance_subY_subY)
      call LA_Matrix_Factorise(LA_c_subYsubY,factor_c_subYsubY,error=error)
      call finalise(LA_c_subYsubY)

      do i = 1, n_globalSparseX-1
         do j = i+1, n_globalSparseX
            factor_c_subYsubY(j,i) = 0.0_qp
         end do
      end do

      a(1:n_globalY,:) = transpose(c_subYY_sqrtInverseLambda)
      a(n_globalY+1:,:) = factor_c_subYsubY

      globalY = 0.0_qp
      do i_y = 1, from%n_y
         ! loop over all function values

         i_globalY = from%map_y_globalY(i_y)
         ! find unique function value/derivative identifier

         globalY(i_globalY) = from%y(i_y)*sqrt(1.0_qp/from%lambda(i_globalY))
      enddo

      do i_yPrime = 1, from%n_yPrime
         ! loop over all function values

         i_global_yPrime = from%map_yPrime_globalY(i_yPrime)
         ! find unique function value/derivative identifier

         globalY(i_global_yPrime) = from%yPrime(i_yPrime)*sqrt(1.0_qp/from%lambda(i_global_yPrime))
      enddo

      call initialise(LA_q_subYsubY,a)
      call LA_Matrix_QR_Solve_Vector(LA_q_subYsubY,globalY,alpha)
      call finalise(LA_q_subYsubY)

      do i_coordinate = 1, from%n_coordinate
         do i_sparseX = 1, from%coordinate(i_coordinate)%n_sparseX
            i_global_sparseX = from%coordinate(i_coordinate)%map_sparseX_globalSparseX(i_sparseX)
            this%coordinate(i_coordinate)%alpha(i_sparseX) = real(alpha(i_global_sparseX),kind=dp)
         enddo
      enddo

      if(allocated(c_subYY_sqrtInverseLambda)) deallocate(c_subYY_sqrtInverseLambda)
      if(allocated(factor_c_subYsubY)) deallocate(factor_c_subYsubY)
      if(allocated(a)) deallocate(a)
      if(allocated(globalY)) deallocate(globalY)
      if(allocated(alpha)) deallocate(alpha)
#else
      allocate( c_subYY_inverseLambda(n_globalSparseX,n_globalY), c_subYY_inverseLambda_c_YsubY(n_globalSparseX,n_globalSparseX), &
!      inverse_q_subYsubY(n_globalSparseX,n_globalSparseX), inverse_c_subYsubY(n_globalSparseX,n_globalSparseX), &
      alpha(n_globalSparseX), globalY(n_globalY))

      call matrix_product_vect_asdiagonal_sub(c_subYY_inverseLambda,from%covariance_subY_Y,1.0_qp/from%lambda) ! O(NM)

      c_subYY_inverseLambda_c_YsubY = matmul(c_subYY_inverseLambda,transpose(from%covariance_subY_Y))
      call initialise(LA_q_subYsubY,from%covariance_subY_subY + c_subYY_inverseLambda_c_YsubY)

      globalY = 0.0_qp
      do i_y = 1, from%n_y
         ! loop over all function values

         i_globalY = from%map_y_globalY(i_y)
         ! find unique function value/derivative identifier

         globalY(i_globalY) = from%y(i_y) !*sqrt(1.0_qp/from%lambda(i_globalY))
      enddo

      do i_yPrime = 1, from%n_yPrime
         ! loop over all function values

         i_global_yPrime = from%map_yPrime_globalY(i_yPrime)
         ! find unique function value/derivative identifier

         globalY(i_global_yPrime) = from%yPrime(i_yPrime) !*sqrt(1.0_qp/from%lambda(i_global_yPrime))
      enddo

      call Matrix_Solve(LA_q_subYsubY,matmul(c_subYY_inverseLambda, globalY),alpha)
      call finalise(LA_q_subYsubY)

      do i_coordinate = 1, from%n_coordinate
         do i_sparseX = 1, from%coordinate(i_coordinate)%n_sparseX
            i_global_sparseX = from%coordinate(i_coordinate)%map_sparseX_globalSparseX(i_sparseX)
            this%coordinate(i_coordinate)%alpha(i_sparseX) = real(alpha(i_global_sparseX),kind=dp)
         enddo
      enddo

      if(allocated(c_subYY_inverseLambda)) deallocate(c_subYY_inverseLambda)
      if(allocated(c_subYY_inverseLambda_c_YsubY)) deallocate(c_subYY_inverseLambda_c_YsubY)
!      if(allocated(inverse_q_subYsubY)) deallocate(inverse_q_subYsubY)
!      if(allocated(inverse_c_subYsubY)) deallocate(inverse_c_subYsubY)
      if(allocated(alpha)) deallocate(alpha)
      if(allocated(globalY)) deallocate(globalY)
#endif

      this%initialised = .true.

   endsubroutine gpSparse_initialise

   subroutine gpSparse_finalise(this,error)
      type(gpSparse), intent(inout) :: this
      integer, optional, intent(out) :: error

      integer :: i_coordinate

      INIT_ERROR(error)

      do i_coordinate = 1, this%n_coordinate
         call finalise(this%coordinate(i_coordinate), error)
      enddo
      deallocate(this%coordinate)

      this%n_coordinate = 0

   endsubroutine gpSparse_finalise

   subroutine gpFull_Finalise(this, error)
      type(gpFull), intent(inout) :: this
      integer, optional, intent(out) :: error

      integer :: i

      INIT_ERROR(error)

      if(.not. this%initialised) return

      if(allocated(this%coordinate)) then
         do i = 1, this%n_coordinate
            call finalise(this%coordinate(i))
         enddo
         deallocate( this%coordinate )
      endif

      if(allocated(this%y)) deallocate( this%y )
      if(allocated(this%yPrime)) deallocate( this%yPrime )
      if(allocated(this%sigma_y)) deallocate( this%sigma_y )
      if(allocated(this%sigma_yPrime)) deallocate( this%sigma_yPrime )
      if(allocated(this%map_y_globalY)) deallocate( this%map_y_globalY )
      if(allocated(this%map_yPrime_globalY)) deallocate( this%map_yPrime_globalY )
      if(allocated(this%covariance_subY_y)) deallocate( this%covariance_subY_y )
      if(allocated(this%covariance_subY_subY)) deallocate( this%covariance_subY_subY )
      if(allocated(this%covarianceDiag_y_y)) deallocate( this%covarianceDiag_y_y )
      if(allocated(this%lambda)) deallocate( this%lambda )
      if(allocated(this%alpha)) deallocate( this%alpha )
      if(allocated(this%inverse_sparse_full)) deallocate( this%inverse_sparse_full )


      this%n_coordinate = 0
      this%n_y = 0
      this%n_yPrime = 0
      this%current_y = 0
      this%current_yPrime = 0

      this%initialised = .false.

   endsubroutine gpFull_Finalise

   subroutine gpCoordinates_Finalise(this, error)
      type(gpCoordinates), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if(.not. this%initialised) return

      if(allocated(this%x)) deallocate( this%x )
      if(allocated(this%xPrime)) deallocate( this%xPrime )
      if(allocated(this%cutoff)) deallocate( this%cutoff )
      if(allocated(this%cutoffPrime)) deallocate( this%cutoffPrime )
      if(allocated(this%theta)) deallocate( this%theta )
      if(allocated(this%permutations)) deallocate(this%permutations)
      if(allocated(this%map_x_y)) deallocate( this%map_x_y )
      if(allocated(this%map_xPrime_yPrime)) deallocate( this%map_xPrime_yPrime )
      if(allocated(this%map_xPrime_x)) deallocate( this%map_xPrime_x )
      if(allocated(this%map_sparseX_globalSparseX)) deallocate( this%map_sparseX_globalSparseX )
      if(allocated(this%config_type)) deallocate( this%config_type )

      if(allocated(this%sparseX_index)) deallocate(this%sparseX_index)
      if(allocated(this%sparseX)) deallocate(this%sparseX)
      if(allocated(this%alpha)) deallocate(this%alpha)
      if(allocated(this%sparseCutoff)) deallocate(this%sparseCutoff)

      if(allocated(this%x_size)) deallocate( this%x_size )
      if(allocated(this%xPrime_size)) deallocate( this%xPrime_size )
      if(allocated(this%covarianceDiag_x_x)) deallocate( this%covarianceDiag_x_x )
      if(allocated(this%covarianceDiag_x_xPrime)) deallocate( this%covarianceDiag_x_xPrime )
      if(allocated(this%covarianceDiag_xPrime_xPrime)) deallocate( this%covarianceDiag_xPrime_xPrime )

      if(allocated(this%sparseX_size)) deallocate( this%sparseX_size )
      if(allocated(this%covarianceDiag_sparseX_sparseX)) deallocate( this%covarianceDiag_sparseX_sparseX )


      call finalise(this%descriptor_str)
      call gpCoordinates_finalise_SparseScore(this)

      this%d = 0
      this%n_x = 0
      this%n_xPrime = 0
      this%delta = 0.0_dp
      this%f0 = 0.0_dp

      this%current_x = 0
      this%current_xPrime = 0

      this%n_sparseX = 0
      this%n_permutations = 0

      this%sparsified = .false.
      this%initialised = .false.

      if(this%bond_real_space) call gpCovariance_bond_real_space_Finalise(this%bond_real_space_cov)
      this%bond_real_space = .false.

      if(this%atom_real_space) call gpCovariance_atom_real_space_Finalise(this%atom_real_space_cov)
      this%atom_real_space = .false.

      this%soap = .false.

   endsubroutine gpCoordinates_Finalise

   subroutine gpCovariance_bond_real_space_Finalise(this, error)
      type(gpCovariance_bond_real_space), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      this%n = 0
      this%delta = 0.0_dp
      this%atom_sigma = 0.0_dp

      if(allocated(this%sparseX_N)) deallocate(this%sparseX_N)
      if(allocated(this%sparseX_pos)) deallocate(this%sparseX_pos)
      if(allocated(this%sparseX_atom_covariance_cutoff)) deallocate(this%sparseX_atom_covariance_cutoff)
      if(allocated(this%sparseX_bond_covariance_cutoff)) deallocate(this%sparseX_bond_covariance_cutoff)
      if(allocated(this%sparseX_self_overlap)) deallocate(this%sparseX_self_overlap)
      this%sparseX_inisialised = .false.

      this%initialised = .false.

   endsubroutine gpCovariance_bond_real_space_Finalise

!   subroutine gpCovariance_bond_real_space_Finalise(this, error)
!      type(gpCovariance_bond_real_space), intent(inout) :: this
!      integer, optional, intent(out) :: error
!
!      INIT_ERROR(error)
!
!      this%theta = 0.0_dp
!      this%delta = 0.0_dp
!      this%integration_steps = 0
!      this%atom_sigma = 0.0_dp
!
!      if(allocated(this%rot_matrix)) deallocate(this%rot_matrix)
!      if(allocated(this%sparseX_N)) deallocate(this%sparseX_N)
!      if(allocated(this%sparseX_pos_rot)) deallocate(this%sparseX_pos_rot)
!      if(allocated(this%sparseX_atom_covariance_cutoff)) deallocate(this%sparseX_atom_covariance_cutoff)
!      if(allocated(this%sparseX_bond_covariance_cutoff)) deallocate(this%sparseX_bond_covariance_cutoff)
!      if(allocated(this%sparseX_invariant_buffer)) deallocate(this%sparseX_invariant_buffer)
!      this%sparseX_inisialised = .false.
!
!      this%initialised = .false.
!
!   endsubroutine gpCovariance_bond_real_space_Finalise

   subroutine gpCovariance_atom_real_space_Finalise(this, error)
      type(gpCovariance_atom_real_space), intent(inout) :: this
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      this%l_max = 0
      this%delta = 0.0_dp
      this%sparseX_inisialised = .false.

      this%initialised = .false.

   endsubroutine gpCovariance_atom_real_space_Finalise

   function gpFull_addFunctionValue(this,y,sigma_y, error)

      type(gpFull), intent(inout) :: this
      real(dp), intent(in) :: y, sigma_y           ! Function value
      integer :: gpFull_addFunctionValue  ! Which function value we added
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_addFunctionValue: object not initialised',error)
      endif

      if( this%current_y == this%n_y ) then
         RAISE_ERROR( 'gpFull_addFunctionValue: object full, no more function values can be added',error)
      endif

      this%current_y = this%current_y + 1
      this%y(this%current_y) = y
      this%sigma_y(this%current_y) = sigma_y

      gpFull_addFunctionValue = this%current_y

   endfunction gpFull_addFunctionValue

   function gpFull_addFunctionDerivative(this, yPrime, sigma_yPrime, error)
      type(gpFull), intent(inout) :: this
      real(dp), intent(in) :: yPrime, sigma_yPrime ! Function value
      integer :: gpFull_addFunctionDerivative  ! Which function value we added
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_addFunctionDerivative: object not initialised',error)
      endif

      if( this%current_yPrime == this%n_yPrime ) then
         RAISE_ERROR( 'gpFull_addFunctionDerivative: object full, no more function values can be added',error)
      endif

      this%current_yPrime = this%current_yPrime + 1
      this%yPrime(this%current_yPrime) = yPrime
      this%sigma_yPrime(this%current_yPrime) = sigma_yPrime

      gpFull_addFunctionDerivative = this%current_yPrime

   endfunction gpFull_addFunctionDerivative

   function gpFull_addCoordinates_2Darray(this,x,i_coordinate,cutoff_in, current_y, config_type, error) result(xLocation)
      type(gpFull), intent(inout) :: this
      real(dp), dimension(:,:), intent(in) :: x
      integer, intent(in) :: i_coordinate
      integer, optional, intent(in) :: current_y, config_type
      real(dp), dimension(:), intent(in), optional :: cutoff_in
      integer, optional, intent(out) :: error

      integer, dimension(:), pointer :: xLocation

      integer :: previous_x, i
      real(dp), dimension(:,:), allocatable :: new_x

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_addCoordinates: object not initialised',error)
      endif

      if( i_coordinate > this%n_coordinate ) then
         RAISE_ERROR( 'gpFull_addCoordinates: access to descriptor '//i_coordinate//' is not possible as number of descriptors is set '//this%n_coordinate ,error)
      endif

      if( .not. this%coordinate(i_coordinate)%initialised ) then
         RAISE_ERROR('gpFull_addCoordinates: '//i_coordinate//'th coordinate object is not initialised',error)
      endif

      if( this%coordinate(i_coordinate)%bond_real_space .or. this%coordinate(i_coordinate)%atom_real_space ) then
         if( size(x,1) > size(this%coordinate(i_coordinate)%x,1) ) then
            allocate( new_x(size(x,1),this%coordinate(i_coordinate)%n_x) )
            new_x = 0.0_dp
            new_x(1:size(this%coordinate(i_coordinate)%x,1),:) = this%coordinate(i_coordinate)%x
            deallocate( this%coordinate(i_coordinate)%x )
            allocate( this%coordinate(i_coordinate)%x(size(x,1),this%coordinate(i_coordinate)%n_x) )
            this%coordinate(i_coordinate)%x = new_x
            deallocate( new_x )
            this%coordinate(i_coordinate)%d = size(x,1)
         end if
      else
!         if( size(x,1) /= this%coordinate(i_coordinate)%d ) then
!            RAISE_ERROR('gpFull_addCoordinates: dimensionality of descriptors '//size(x,1)//' does not match what is given in the object '//this%coordinate(i_coordinate)%d,error)
!         endif
      endif

      previous_x = this%coordinate(i_coordinate)%current_x
      this%coordinate(i_coordinate)%current_x = previous_x + size(x,2)

      if( this%coordinate(i_coordinate)%current_x > this%coordinate(i_coordinate)%n_x ) then
         RAISE_ERROR('gpFull_addCoordinates: object full, no more descriptors can be added',error)
      endif

      if( this%coordinate(i_coordinate)%bond_real_space .or. this%coordinate(i_coordinate)%atom_real_space ) then
         this%coordinate(i_coordinate)%x(1:size(x,1),previous_x+1:this%coordinate(i_coordinate)%current_x) = x
         this%coordinate(i_coordinate)%x_size(previous_x+1:this%coordinate(i_coordinate)%current_x) = size(x,1)
      else
         this%coordinate(i_coordinate)%x(:,previous_x+1:this%coordinate(i_coordinate)%current_x) = x
      endif

      if(present(cutoff_in)) then
         this%coordinate(i_coordinate)%cutoff(previous_x+1:this%coordinate(i_coordinate)%current_x) = cutoff_in
      endif

      if(present(current_y)) &
         this%coordinate(i_coordinate)%map_x_y(previous_x+1:this%coordinate(i_coordinate)%current_x) = current_y

      if(present(config_type)) &
         this%coordinate(i_coordinate)%config_type(previous_x+1:this%coordinate(i_coordinate)%current_x) = config_type

      allocate(xLocation(size(x,2)))
      xLocation = (/ ( i, i = previous_x+1, this%coordinate(i_coordinate)%current_x ) /)

   endfunction gpFull_addCoordinates_2Darray

   function gpFull_addCoordinates_1Darray(this,x,i_coordinate,cutoff_in,current_y,config_type, error) result(xLocation)
      type(gpFull), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: x
      integer, intent(in) :: i_coordinate
      real(dp), optional, intent(in) :: cutoff_in
      integer, optional, intent(in) :: current_y, config_type
      integer, optional, intent(out) :: error

      integer :: xLocation

      integer, dimension(:), pointer :: xLocation_in

      INIT_ERROR(error)

      xLocation_in => gpFull_addCoordinates_2Darray(this,reshape(x,(/size(x),1/)),i_coordinate,(/cutoff_in/),current_y,config_type,error)

      xLocation =  xLocation_in(1)
      deallocate(xLocation_in)

   endfunction gpFull_addCoordinates_1Darray

   subroutine gpFull_addCoordinateDerivatives_2Darray(this,xPrime,i_coordinate,current_yPrime, xLocation, dcutoff_in, error)
      type(gpFull), intent(inout) :: this
      real(dp), dimension(:,:), intent(in) :: xPrime
      integer, intent(in) :: i_coordinate, current_yPrime
      integer, dimension(:), intent(in) :: xLocation
      real(dp), dimension(:), optional, intent(in) :: dcutoff_in
      integer, optional, intent(out) :: error

      integer :: previous_xPrime
      real(dp), dimension(:,:), allocatable :: new_xPrime

      integer :: i, first_nonzero, last_nonzero

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_addCoordinateDerivatives: object not initialised',error)
      endif

      if( i_coordinate > this%n_coordinate ) then
         RAISE_ERROR( 'gpFull_addCoordinateDerivatives: access to descriptor '//i_coordinate//' is not possible as number of descriptors is set '//this%n_coordinate,error )
      endif

      if( .not. this%coordinate(i_coordinate)%initialised ) then
         RAISE_ERROR('gpFull_addCoordinateDerivatives: '//i_coordinate//'th coordinate object is not initialised',error)
      endif

      if( this%coordinate(i_coordinate)%bond_real_space ) then
         if( size(xPrime,1) > size(this%coordinate(i_coordinate)%xPrime,1) ) then
            allocate( new_xPrime(size(xPrime,1),this%coordinate(i_coordinate)%n_xPrime) )
            new_xPrime = 0.0_dp
            new_xPrime(1:size(this%coordinate(i_coordinate)%xPrime,1),:) = this%coordinate(i_coordinate)%xPrime
            deallocate( this%coordinate(i_coordinate)%xPrime )
            allocate( this%coordinate(i_coordinate)%xPrime(size(xPrime,1),this%coordinate(i_coordinate)%n_xPrime) )
            this%coordinate(i_coordinate)%xPrime = new_xPrime
            deallocate( new_xPrime )
         end if
      elseif( this%coordinate(i_coordinate)%atom_real_space ) then

         ! Put non-zero part of data into the container. Check if everything is the right size. xPrime_size gets the first non-zero index so later the array can be reconstructed.
         previous_xPrime = this%coordinate(i_coordinate)%current_xPrime
         do i = 1, size(xPrime,2)
            call find_nonzero_chunk(xPrime(:,i), first = first_nonzero, last = last_nonzero)
            if( last_nonzero - first_nonzero + 1 > size(this%coordinate(i_coordinate)%xPrime,1) ) then ! something's wrong: more data, initialised with the wrong size...
               RAISE_ERROR('gpFull_addCoordinateDerivatives: wrong number of elements non-zero in case of atom real space descriptors',error)
            endif

            if( first_nonzero + size(this%coordinate(i_coordinate)%xPrime,1) - 1 > size(xPrime,1) ) & ! data is at the very end, and there are extra zeros in the front of the chunk, adjusting starting position to include them
               first_nonzero = size(xPrime,1) + 1 - size(this%coordinate(i_coordinate)%xPrime,1)
 
            if( last_nonzero - size(this%coordinate(i_coordinate)%xPrime,1) + 1 < 1 ) & ! data is at the very beginning, and there are extra zeros in the end of the chunk. adjusting finishing position
               last_nonzero = size(this%coordinate(i_coordinate)%xPrime,1)

            if( last_nonzero - first_nonzero + 1 < size(this%coordinate(i_coordinate)%xPrime,1) ) & ! data chunk is in the middle and there are zeros either at the front or at the end. adjusting finishing so array is padded with zeros at the end.
               last_nonzero = first_nonzero + size(this%coordinate(i_coordinate)%xPrime,1) - 1

            ! Now the array should be exactly the right size...

            this%coordinate(i_coordinate)%xPrime(:,previous_xPrime+i) = xPrime(first_nonzero:last_nonzero,i)
            this%coordinate(i_coordinate)%xPrime_size(previous_xPrime+i) = first_nonzero
         enddo
      else
         if( size(xPrime,1) /= this%coordinate(i_coordinate)%d ) then
            RAISE_ERROR('gpFull_addCoordinateDerivatives: dimensionality of descriptors '//size(xPrime,1)//' does not match what is given in the object '//this%coordinate(i_coordinate)%d,error)
         endif
      endif

      if( size(xPrime,2) /= size(xLocation) ) then
         RAISE_ERROR('gpFull_addCoordinateDerivatives: number of descriptors '//size(xPrime,2)//' has to match the dimensionality of the mapping array '//size(xLocation),error)
      endif

      previous_xPrime = this%coordinate(i_coordinate)%current_xPrime
      this%coordinate(i_coordinate)%current_xPrime = previous_xPrime + size(xPrime,2)

      if( this%coordinate(i_coordinate)%current_xPrime > this%coordinate(i_coordinate)%n_xPrime ) then
         RAISE_ERROR('gpFull_addCoordinateDerivatives: object full, no more descriptors can be added',error)
      endif

      if( this%coordinate(i_coordinate)%bond_real_space ) then
         this%coordinate(i_coordinate)%xPrime(1:size(xPrime,1),previous_xPrime+1:this%coordinate(i_coordinate)%current_xPrime) = xPrime
         this%coordinate(i_coordinate)%xPrime_size(previous_xPrime+1:this%coordinate(i_coordinate)%current_xPrime) = size(xPrime,1)
      elseif( this%coordinate(i_coordinate)%atom_real_space ) then
         ! everthing has been done before...
      else
         this%coordinate(i_coordinate)%xPrime(:,previous_xPrime+1:this%coordinate(i_coordinate)%current_xPrime) = xPrime
      endif

      if(present(dcutoff_in)) then
         this%coordinate(i_coordinate)%cutoffPrime(previous_xPrime+1:this%coordinate(i_coordinate)%current_xPrime) = dcutoff_in
      endif

      this%coordinate(i_coordinate)%map_xPrime_yPrime(previous_xPrime+1:this%coordinate(i_coordinate)%current_xPrime) = current_yPrime
      this%coordinate(i_coordinate)%map_xPrime_x(previous_xPrime+1:this%coordinate(i_coordinate)%current_xPrime) = xLocation

   endsubroutine gpFull_addCoordinateDerivatives_2Darray

   subroutine gpFull_addCoordinateDerivatives_1Darray(this,xPrime,i_coordinate,current_yPrime, xLocation, dcutoff_in, error)
      type(gpFull), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: xPrime
      integer, intent(in) :: i_coordinate, current_yPrime
      integer, intent(in) :: xLocation
      real(dp), optional, intent(in) :: dcutoff_in
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      call gpFull_addCoordinateDerivatives_2Darray(this, reshape(xPrime,(/size(xPrime),1/)),i_coordinate,current_yPrime,(/xLocation/),(/dcutoff_in/),error)

   endsubroutine gpFull_addCoordinateDerivatives_1Darray

   subroutine gpFull_addDescriptor(this,i_coordinate,descriptor_str,error)

      type(gpFull), intent(inout) :: this
      integer, intent(in) :: i_coordinate 
      character(len=*), intent(in) :: descriptor_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_addDescriptor: object not initialised',error)
      endif

      call gpCoordinates_addDescriptor(this%coordinate(i_coordinate),descriptor_str,error)

   endsubroutine gpFull_addDescriptor

   subroutine gpCoordinates_addDescriptor(this,descriptor_str,error)

      type(gpCoordinates), intent(inout) :: this
      character(len=*), intent(in) :: descriptor_str
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_addDescriptor: object not initialised',error)
      endif

      call initialise(this%descriptor_str)
      call zero(this%descriptor_str)

      call concat(this%descriptor_str,descriptor_str,keep_lf=.false.)

   endsubroutine gpCoordinates_addDescriptor

   subroutine gpCoordinates_setTheta(this, theta,error)
      type(gpCoordinates), intent(inout), target :: this
      real(dp), dimension(:), intent(in) :: theta
      integer, optional, intent(out) :: error

      integer :: i, x_i_size
      real(dp) :: delta
      real(dp), dimension(:), pointer :: xPrime_i

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_calculateTheta: object not initialised',error)
      endif

      call check_size('theta',theta,shape(this%theta),'gpCoordinates_setTheta',error)

      allocate(this%covarianceDiag_x_xPrime(size(this%x,1),this%n_x))
      this%covarianceDiag_x_xPrime = 0.0_dp

      if(this%bond_real_space) then
         if(.not. this%bond_real_space_cov%initialised) then
            call gpCoordinates_gpCovariance_bond_real_space_Initialise(this)
         endif

         delta = this%bond_real_space_cov%delta
         this%bond_real_space_cov%delta = 1.0_dp

         do i = 1, this%n_x
            this%covarianceDiag_x_x(i) = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, x_i = this%x(:,i), x_i_size = this%x_size(i), x_j = this%x(:,i), x_j_size = this%x_size(i))
         enddo

!         do i = 1, this%n_xPrime
!         enddo

         this%bond_real_space_cov%delta = delta

      elseif(this%atom_real_space) then
         if(.not. this%atom_real_space_cov%initialised) then
            call gpCoordinates_gpCovariance_atom_real_space_Initialise(this)
         endif

         do i = 1, this%n_x
            x_i_size = this%x_size(i)
            xPrime_i => this%covarianceDiag_x_xPrime(:,i)
            this%covarianceDiag_x_x(i) = gpCovariance_atom_real_space_Calc(this%atom_real_space_cov, x_i = this%x(:,i), x_i_size = x_i_size, x_j = this%x(:,i), x_j_size = x_i_size, xPrime_i=xPrime_i)
         enddo
      elseif(this%soap) then
         this%theta = theta
      else
         ! Symmetrise theta wrt permutations
         this%theta = 0.0_dp
         do i = 1, this%n_permutations
            this%theta = this%theta + theta(this%permutations(:,i))
         enddo
         this%theta = this%theta / real(this%n_permutations,kind=dp)
      endif

   endsubroutine gpCoordinates_setTheta

   subroutine gpCoordinates_setThetaFactor(this, thetaFactor,useSparseX,error)
      type(gpCoordinates), intent(inout) :: this
      real(dp), dimension(:), intent(in) :: thetaFactor
      logical, optional, intent(in) :: useSparseX
      integer, optional, intent(out) :: error

      integer :: i
      logical :: my_useSparseX
      real(dp), dimension(this%d) :: theta

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_calculateThetaFactor: object not initialised',error)
      endif

      my_useSparseX = .false.
      if( allocated(this%sparseX_index) ) then
         if( sum( this%sparseX_index ) > 0 ) my_useSparseX = optional_default(.true.,useSparseX)
      endif

      if( my_useSparseX ) then
         do i = 1, this%d
            theta(i) = ( maxval(this%x(i,this%sparseX_index)) - minval(this%x(i,this%sparseX_index)) ) * thetaFactor(i)
            if( theta(i) < THETA_MIN ) theta(i) = 1.0_dp
         enddo
      else
         do i = 1, this%d
            theta(i) = ( maxval(this%x(i,:)) - minval(this%x(i,:)) ) * thetaFactor(i)
            if( theta(i) < THETA_MIN ) theta(i) = 1.0_dp
         enddo
      endif

      call gpCoordinates_setTheta(this,theta,error)

   endsubroutine gpCoordinates_setThetaFactor

   subroutine gpFull_setTheta(this, i_coordinate, theta, error)
      type(gpFull), intent(inout) :: this
      integer, intent(in) :: i_coordinate
      real(dp), dimension(:), intent(in) :: theta
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_setTheta: object not initialised',error)
      endif

      call gpCoordinates_setTheta(this%coordinate(i_coordinate), theta, error)

   endsubroutine gpFull_setTheta

   subroutine gpFull_setTheta_thetaFactor(this, i_coordinate, thetaFactor, useSparseX, error)
      type(gpFull), intent(inout) :: this
      integer, intent(in) :: i_coordinate
      real(dp), dimension(:), intent(in) :: thetaFactor
      logical, optional, intent(in) :: useSparseX
      integer, optional, intent(out) :: error

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_setTheta_thetaFactor: object not initialised',error)
      endif

      call gpCoordinates_setThetaFactor(this%coordinate(i_coordinate), thetaFactor, useSparseX, error)

   endsubroutine gpFull_setTheta_thetaFactor

!   subroutine gpFull_setTheta_thetaFactorArray(this, thetaFactor, useSparseX, error)
!      type(gpFull), intent(inout) :: this
!      real(dp), dimension(:), intent(in) :: thetaFactor
!      logical, optional, intent(out) :: useSparseX
!      integer, optional, intent(out) :: error
!
!      integer :: i
!
!      INIT_ERROR(error)
!
!      if( .not. this%initialised ) then
!         RAISE_ERROR('gpFull_setTheta_thetaFactorArray: object not initialised',error)
!      endif
!
!      call check_size('thetaFactor',thetaFactor,(/this%n_coordinate/),'gpFull_setTheta_thetaFactorArray',error)
!
!      do i = 1, this%n_coordinate
!         call gpCoordinates_setThetaFactor(this%coordinate(i), thetaFactor(i), useSparseX, error)
!      enddo
!
!   endsubroutine gpFull_setTheta_thetaFactorArray
!
!   subroutine gpFull_setTheta_thetaFactorUniform(this, thetaFactor, useSparseX, error)
!      type(gpFull), intent(inout) :: this
!      real(dp), intent(in) :: thetaFactor
!      logical, optional, intent(out) :: useSparseX
!      integer, optional, intent(out) :: error
!
!      integer :: i
!
!      INIT_ERROR(error)
!
!      if( .not. this%initialised ) then
!         RAISE_ERROR('gpFull_setTheta_thetaFactorUniform: object not initialised',error)
!      endif
!
!      do i = 1, this%n_coordinate
!         call gpCoordinates_setThetaFactor(this%coordinate(i), thetaFactor, useSparseX, error)
!      enddo
!
!   endsubroutine gpFull_setTheta_thetaFactorUniform

   subroutine gpFull_covarianceMatrix_sparse(this,error)
      type(gpFull), intent(inout), target :: this
      integer, optional, intent(out) :: error

      integer :: i_coordinate, i_global_sparseX, j_global_sparseX, i_sparseX, j_sparseX, &
      n_globalY, i_globalY, i_global_yPrime, i_y, i_yPrime, i_x, j_x, n_x, i_xPrime, j_xPrime, n_xPrime, &
      n, i, first_nonzero, last_nonzero, first_nonzero_i, last_nonzero_i, first_nonzero_j, last_nonzero_j
      real(dp) :: covariance_xPrime_sparseX, covariance_x_x_single, covariance_xPrime_xPrime, fc_i, fc_j, dfc_i, dfc_j

      real(dp), dimension(:,:,:,:), allocatable :: grad2_Covariance
      real(dp), dimension(:,:,:), allocatable :: grad_Covariance_j
      real(dp), dimension(:,:), allocatable :: grad_Covariance_i, covariance_x_x
      real(dp), dimension(:), allocatable :: covariance_x_sparseX, covariance_subY_currentX_y, covariance_subY_currentX_suby

      real(dp), dimension(:), pointer :: xPrime_i, xPrime_j
      integer, dimension(:), allocatable :: xIndex, xPrime_Index, xPrime_x_Index
      type(LA_matrix) :: LA_covariance_subY_subY
      logical :: found_i_x

      INIT_ERROR(error)

      call system_timer('gpFull_covarianceMatrix_sparse')

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_covarianceMatrix: object not initialised',error)
      endif

      this%n_globalSparseX = 0
      i_global_sparseX = 0

      do i_coordinate = 1, this%n_coordinate
         if( .not. this%coordinate(i_coordinate)%initialised ) then
            RAISE_ERROR('gpFull_covarianceMatrix: '//i_coordinate//'th coordinate object not initialised',error)
         endif

         if( .not. this%coordinate(i_coordinate)%sparsified ) then
            RAISE_ERROR('gpFull_covarianceMatrix: '//i_coordinate//'th coordinate object not sparsified',error)
         endif

         if(this%coordinate(i_coordinate)%bond_real_space) then
            if((.not. this%coordinate(i_coordinate)%bond_real_space_cov%initialised) .or. (.not. this%coordinate(i_coordinate)%bond_real_space_cov%sparseX_inisialised)) then
               call gpCoordinates_gpCovariance_bond_real_space_Initialise(this%coordinate(i_coordinate))
            endif
         endif

         if(this%coordinate(i_coordinate)%atom_real_space) then
            if((.not. this%coordinate(i_coordinate)%atom_real_space_cov%initialised)) then
               call gpCoordinates_gpCovariance_atom_real_space_Initialise(this%coordinate(i_coordinate))
            endif
         endif

         this%n_globalSparseX = this%n_globalSparseX + this%coordinate(i_coordinate)%n_sparseX

         do i_sparseX = 1, this%coordinate(i_coordinate)%n_sparseX
            i_global_sparseX = i_global_sparseX + 1
            this%coordinate(i_coordinate)%map_sparseX_globalSparseX(i_sparseX) = i_global_sparseX
         enddo
      enddo

      n_globalY = this%n_y + this%n_yPrime

      i_globalY = 0
      do i_y = 1, this%n_y
         ! loop over all function values

         i_globalY = i_globalY + 1
         this%map_y_globalY(i_y) = i_globalY
      enddo
      do i_yPrime = 1, this%n_yPrime
         ! loop over all derivative values

         i_globalY = i_globalY + 1
         this%map_yPrime_globalY(i_yPrime) = i_globalY
      enddo

      call reallocate(this%covariance_subY_y, this%n_globalSparseX, n_globalY, zero = .true.)
      call reallocate(this%covariance_subY_subY, this%n_globalSparseX, this%n_globalSparseX, zero = .true.)
      call reallocate(this%covarianceDiag_y_y, n_globalY, zero = .true.)
      call reallocate(this%lambda, n_globalY, zero = .true.)
      call reallocate(this%inverse_sparse_full, this%n_globalSparseX, n_globalY, zero = .true.)

      allocate( covariance_subY_currentX_y(n_globalY),covariance_subY_currentX_suby(this%n_globalSparseX) )
      covariance_subY_currentX_y = 0.0_dp
      covariance_subY_currentX_suby = 0.0_dp

      do i_coordinate = 1, this%n_coordinate
         ! loop over different descriptor types
         call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate)

         call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate//'_sparse')
         do i_sparseX = 1, this%coordinate(i_coordinate)%n_sparseX
            ! loop over sparse points of each descriptor

            i_global_sparseX = this%coordinate(i_coordinate)%map_sparseX_globalSparseX(i_sparseX)
            ! find the unique number of the sparse point (to refer to it outside of the context of descriptor)

            allocate( grad_Covariance_i(this%coordinate(i_coordinate)%d,this%coordinate(i_coordinate)%n_x), &
                      covariance_x_sparseX(this%coordinate(i_coordinate)%n_x) )

            covariance_subY_currentX_y = 0.0_dp
            covariance_subY_currentX_suby = 0.0_dp
!$omp parallel do default(none) shared(this,i_coordinate,covariance_x_sparseX,grad_Covariance_i,i_sparseX) private(i_x,i_y,i_globalY) reduction(+:covariance_subY_currentX_y)
            do i_x = 1, this%coordinate(i_coordinate)%n_x
               ! loop over all data

               covariance_x_sparseX(i_x) = gpCoordinates_Covariance(this%coordinate(i_coordinate), &
		  i_x = i_x, j_sparseX = i_sparseX, grad_Covariance_i = grad_Covariance_i(:,i_x))

               i_y = this%coordinate(i_coordinate)%map_x_y(i_x)
               ! find which function value depends on the given descriptor

               if( i_y /= 0 ) then
                  i_globalY = this%map_y_globalY(i_y)
                  ! find unique function value/derivative identifier

                  !this%covariance_subY_y(i_global_sparseX, i_globalY) = this%covariance_subY_y(i_global_sparseX, i_globalY) + &
                  covariance_subY_currentX_y(i_globalY) = covariance_subY_currentX_y(i_globalY) + &
		     covariance_x_sparseX(i_x)*this%coordinate(i_coordinate)%cutoff(i_x)*this%coordinate(i_coordinate)%sparseCutoff(i_sparseX)
               endif
            enddo

!$omp parallel do default(none) shared(this,i_coordinate,i_sparseX,grad_Covariance_i,covariance_x_sparseX) private(i_xPrime,i_yPrime,i_x,i_global_yPrime,covariance_xPrime_sparseX,first_nonzero,last_nonzero) reduction(+:covariance_subY_currentX_y)
	    do i_xPrime = 1, this%coordinate(i_coordinate)%n_xPrime
	       ! loop over all derivative data

	       i_yPrime = this%coordinate(i_coordinate)%map_xPrime_yPrime(i_xPrime)
	       ! find which derivative depends on the given descriptor

	       i_x = this%coordinate(i_coordinate)%map_xPrime_x(i_xPrime)
	       if( i_yPrime /= 0 ) then
		  i_global_yPrime = this%map_yPrime_globalY(i_yPrime)
		  ! find unique function value/derivative identifier

		  if(this%coordinate(i_coordinate)%atom_real_space) then
		     first_nonzero = this%coordinate(i_coordinate)%xPrime_size(i_xPrime)
		  else
		     first_nonzero = 1
		  endif
		  last_nonzero = first_nonzero + size(this%coordinate(i_coordinate)%xPrime,1) - 1

		  ! on Xeon w/ ifort 12, sum is fastest .  ddot is close.  dot_product is terrible
		  ! on Opteron w/ ifort 12 acml 5.2, ddot is 14.95 s, dot_product is 22.5 s, and sum is 13.9 s
		  ! dgemv doesn't seem noticeably faster at Opterons (may be faster on Xeon for 'N' transpose setting)
		  ! covariance_xPrime_sparseX = ddot(size(this%coordinate(i_coordinate)%xPrime,1),grad_Covariance_i(first_nonzero,i_x),1,this%coordinate(i_coordinate)%xPrime(1,i_xPrime),1)*&
		  ! covariance_xPrime_sparseX = dot_product(grad_Covariance_i(first_nonzero:last_nonzero,i_x),this%coordinate(i_coordinate)%xPrime(:,i_xPrime))* &
		  covariance_xPrime_sparseX = sum(grad_Covariance_i(first_nonzero:last_nonzero,i_x)*this%coordinate(i_coordinate)%xPrime(:,i_xPrime))* &
		     this%coordinate(i_coordinate)%cutoff(i_x)*this%coordinate(i_coordinate)%sparseCutoff(i_sparseX) + &
		     covariance_x_sparseX(i_x)*this%coordinate(i_coordinate)%cutoffPrime(i_xPrime)*this%coordinate(i_coordinate)%sparseCutoff(i_sparseX)

		  !this%covariance_subY_y(i_global_sparseX, i_global_yPrime) = this%covariance_subY_y(i_global_sparseX, i_global_yPrime) + covariance_xPrime_sparseX
		  covariance_subY_currentX_y(i_global_yPrime) = covariance_subY_currentX_y(i_global_yPrime) + covariance_xPrime_sparseX
	       endif
	    enddo


            if(allocated(grad_Covariance_i)) deallocate(grad_Covariance_i)
            if(allocated(covariance_x_sparseX)) deallocate(covariance_x_sparseX)

!$omp parallel do default(none) shared(this,i_coordinate,covariance_x_sparseX,grad_Covariance_i,i_sparseX,i_global_sparseX) private(j_sparseX,j_global_sparseX) reduction(+:covariance_subY_currentX_suby)
            do j_sparseX = 1, this%coordinate(i_coordinate)%n_sparseX
               ! loop over sparse points of each descriptor

               j_global_sparseX = this%coordinate(i_coordinate)%map_sparseX_globalSparseX(j_sparseX)
               ! find the unique number of the sparse point (to refer to it outside of the context of descriptor)
               covariance_subY_currentX_suby(j_global_sparseX) = covariance_subY_currentX_suby(j_global_sparseX) + &
               gpCoordinates_Covariance(this%coordinate(i_coordinate), j_sparseX = j_sparseX, i_sparseX = i_sparseX) * this%coordinate(i_coordinate)%sparseCutoff(i_sparseX)*this%coordinate(i_coordinate)%sparseCutoff(j_sparseX)
            enddo

            this%covariance_subY_subY(i_global_sparseX,i_global_sparseX) = this%covariance_subY_subY(i_global_sparseX,i_global_sparseX) + this%sparse_jitter
            this%covariance_subY_y(i_global_sparseX,:) = this%covariance_subY_y(i_global_sparseX,:) + covariance_subY_currentX_y
            this%covariance_subY_subY(:,i_global_sparseX) = this%covariance_subY_subY(:,i_global_sparseX) + covariance_subY_currentX_suby

         enddo
         call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate//'_sparse')

         select case(this%covariance_method)
         case(GP_COVARIANCE_FITC)
            call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate//'_diag')
            do i_y = 1, this%n_y
               ! loop over all function values

               i_globalY = this%map_y_globalY(i_y)
               ! find unique function value/derivative identifier

               n_x = count( this%coordinate(i_coordinate)%map_x_y == i_y )
               allocate(xIndex(n_x))
               
               n = 1
               do i = 1, size(this%coordinate(i_coordinate)%map_x_y)
                  if(this%coordinate(i_coordinate)%map_x_y(i) == i_y) then
                     xIndex(n) = i
                     n = n+1
                  endif
               enddo

               ! construct index array of all contributions to a given function value

               do i_x = 1, n_x
                  do j_x = 1, n_x
                     covariance_x_x_single = gpCoordinates_Covariance( this%coordinate(i_coordinate), i_x = xIndex(i_x), j_x = xIndex(j_x) ) * &
                     this%coordinate(i_coordinate)%cutoff(xIndex(i_x))*this%coordinate(i_coordinate)%cutoff(xIndex(j_x))
                     this%covarianceDiag_y_y(i_globalY) = this%covarianceDiag_y_y(i_globalY) + covariance_x_x_single
                  enddo
               enddo
               deallocate(xIndex)
            enddo
            call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate//'_diag')

            call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate//'_grad_diag')

            do i_yPrime = 1, this%n_yPrime
               ! loop over all derivative values
               i_global_yPrime = this%map_yPrime_globalY(i_yPrime)
               ! find unique function value/derivative identifier

               !n_xPrime = count( this%coordinate(i_coordinate)%map_xPrime_yPrime == i_yPrime )
               allocate(xPrime_Index(size(this%coordinate(i_coordinate)%map_xPrime_yPrime)))
               
               ! construct index array of all contributions to a given derivative value
               n = 0
               do i = 1, size(this%coordinate(i_coordinate)%map_xPrime_yPrime)
                  if(this%coordinate(i_coordinate)%map_xPrime_yPrime(i) == i_yPrime) then
                     n = n + 1
                     xPrime_Index(n) = i
                  endif
               enddo
               n_xPrime = n
               call reallocate(xPrime_Index,n_xPrime,copy=.true.)

               ! xIndex: for each xPrime that contributes to the given i_yPrime, this contains the indices 
               ! of the x in the main database. Length is less or equal to the number of xPrimes.
               ! xPrime_x_Index: this is the inverse of xIndex, i.e. for each xPrime this contains the index
               ! of x in xIndex.
               allocate(xPrime_x_Index(n_xPrime),xIndex(n_xPrime))
               xIndex = -1
               n_x = 1
               xIndex(1) = this%coordinate(i_coordinate)%map_xPrime_x(xPrime_Index(1))
               xPrime_x_Index(1) = 1
               ! Search started off.
               do i_xPrime = 2, n_xPrime
                  i_x = this%coordinate(i_coordinate)%map_xPrime_x(xPrime_Index(i_xPrime))
                  found_i_x = .false.
                  do i = 1, n_x
                     if(xIndex(i) == i_x) then
                        xPrime_x_Index(i_xPrime) = i
                        found_i_x = .true.
                        exit
                     endif
                  enddo
                  if(.not. found_i_x) then
                     n_x = n_x + 1
                     xIndex(n_x) = i_x
                     xPrime_x_Index(i_xPrime) = n_x
                  endif

                  !if( any(xIndex(:n_x) == i_x) ) then
                  !   do i = 1, n_x
                  !      if(xIndex(i) == i_x) xPrime_x_Index(i_xPrime) = i
                  !   enddo
                  !   cycle
                  !else
                  !   n_x = n_x + 1
                  !   xIndex(n_x) = i_x
                  !   xPrime_x_Index(i_xPrime) = n_x
                  !endif
               enddo
               call reallocate(xIndex,n_x,copy=.true.)


               allocate( grad2_Covariance(this%coordinate(i_coordinate)%d,this%coordinate(i_coordinate)%d,n_x,n_x), &
              grad_Covariance_j(this%coordinate(i_coordinate)%d,n_x,n_x), covariance_x_x(n_x,n_x) ) !, &
   !            xPrime_i_grad2_Covariance(this%coordinate(i_coordinate)%d) )
               grad2_Covariance = 0.0_dp
               do i_x = 1, n_x
                  do j_x = 1, n_x
                     covariance_x_x(j_x,i_x) = gpCoordinates_Covariance( this%coordinate(i_coordinate), i_x = xIndex(i_x), j_x = xIndex(j_x), &
                     grad_Covariance_j = grad_Covariance_j(:,j_x,i_x), grad2_Covariance=grad2_Covariance(:,:,j_x,i_x) )
                  enddo
               enddo

               do i_xPrime = 1, n_xPrime
                  i_x = xPrime_x_Index(i_xPrime)
                  xPrime_i => this%coordinate(i_coordinate)%xPrime(:,xPrime_Index(i_xPrime))
                  
                  if(this%coordinate(i_coordinate)%atom_real_space) then
                     first_nonzero_i = this%coordinate(i_coordinate)%xPrime_size(xPrime_Index(i_xPrime))
                     last_nonzero_i = first_nonzero_i + size(this%coordinate(i_coordinate)%xPrime,1) - 1
                  else
                     first_nonzero_i = 1
                     last_nonzero_i = size(this%coordinate(i_coordinate)%xPrime_size,1)
                  endif

                  fc_i = this%coordinate(i_coordinate)%cutoff(xIndex(i_x))
                  dfc_i = this%coordinate(i_coordinate)%cutoffPrime(xPrime_Index(i_xPrime))
                  do j_xPrime = 1, n_xPrime
                     j_x = xPrime_x_Index(j_xPrime)
                     xPrime_j => this%coordinate(i_coordinate)%xPrime(:,xPrime_Index(j_xPrime))

                     if(this%coordinate(i_coordinate)%atom_real_space) then
                        first_nonzero_j = this%coordinate(i_coordinate)%xPrime_size(xPrime_Index(j_xPrime))
                        last_nonzero_j = first_nonzero_j + size(this%coordinate(i_coordinate)%xPrime,1) - 1
                     else
                        first_nonzero_j = 1
                        last_nonzero_j = size(this%coordinate(i_coordinate)%xPrime_size,1)
                     endif
                     
                     fc_j = this%coordinate(i_coordinate)%cutoff(xIndex(j_x))
                     dfc_j = this%coordinate(i_coordinate)%cutoffPrime(xPrime_Index(j_xPrime))

                     covariance_xPrime_xPrime = dot_product(matmul(xPrime_i,grad2_Covariance(first_nonzero_i:last_nonzero_i,first_nonzero_j:last_nonzero_j,j_x,i_x)),xPrime_j)*&
                     fc_i*fc_j + &
                     dot_product(grad_Covariance_j(first_nonzero_i:last_nonzero_i,j_x,i_x),xPrime_i)*fc_i*dfc_j + &
                     dot_product(grad_Covariance_j(first_nonzero_j:last_nonzero_j,j_x,i_x),xPrime_j)*fc_j*dfc_i + &
                     covariance_x_x(j_x,i_x)*dfc_i*dfc_j
                     !gpCoordinates_Covariance( this%coordinate(i_coordinate), i_xPrime = xPrime_Index(i_xPrime), j_xPrime = xPrime_Index(j_xPrime) )
                     this%covarianceDiag_y_y(i_global_yPrime) = this%covarianceDiag_y_y(i_global_yPrime) + covariance_xPrime_xPrime
                  enddo
               enddo

               if(allocated(xPrime_Index)) deallocate(xPrime_Index)
               if(allocated(xIndex)) deallocate(xIndex)
               if(allocated(xPrime_x_Index)) deallocate(xPrime_x_Index)
               if(allocated(grad2_Covariance)) deallocate(grad2_Covariance)
               if(allocated(grad_Covariance_j)) deallocate(grad_Covariance_j)
               if(allocated(covariance_x_x)) deallocate(covariance_x_x)
            enddo
            call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate//'_grad_diag')
         case(GP_COVARIANCE_DTC)
            this%covarianceDiag_y_y = 0.0_dp
         case default
            this%covarianceDiag_y_y = 0.0_dp
         endselect

         call system_timer('gpFull_covarianceMatrix_sparse_Coordinate'//i_coordinate)
      enddo

      call system_timer('gpFull_covarianceMatrix_sparse_LinearAlgebra')
      call initialise(LA_covariance_subY_subY,this%covariance_subY_subY)
      call Matrix_Solve(LA_covariance_subY_subY, this%covariance_subY_y, this%inverse_sparse_full,error=error)
      call finalise(LA_covariance_subY_subY)
      call system_timer('gpFull_covarianceMatrix_sparse_LinearAlgebra')

      call system_timer('gpFull_covarianceMatrix_sparse_FunctionValues')
      select case(this%covariance_method)
      case(GP_COVARIANCE_FITC)
         do i_globalY = 1, n_globalY
            ! loop over all function values

            this%lambda(i_globalY) = this%covarianceDiag_y_y(i_globalY) - &
               dot_product( this%inverse_sparse_full(:,i_globalY), this%covariance_subY_y(:,i_globalY) )
         enddo
      case(GP_COVARIANCE_DTC)
         this%lambda = 0.0_dp
      case default
         this%lambda = 0.0_dp
      endselect

      do i_y = 1, this%n_y
         ! loop over all function values

         i_globalY = this%map_y_globalY(i_y)
         ! find unique function value/derivative identifier

         this%lambda(i_globalY) = this%lambda(i_globalY) + &
            this%sigma_y(i_y)**2
      enddo

      do i_yPrime = 1, this%n_yPrime
         ! loop over all function values

         i_global_yPrime = this%map_yPrime_globalY(i_yPrime)
         ! find unique function value/derivative identifier

         this%lambda(i_global_yPrime) = this%lambda(i_global_yPrime) + &
            this%sigma_yPrime(i_yPrime)**2
      enddo
      call system_timer('gpFull_covarianceMatrix_sparse_FunctionValues')

      if(allocated(covariance_subY_currentX_y)) deallocate(covariance_subY_currentX_y)
      if(allocated(covariance_subY_currentX_suby)) deallocate(covariance_subY_currentX_suby)
      call system_timer('gpFull_covarianceMatrix_sparse')

   endsubroutine gpFull_covarianceMatrix_sparse

   subroutine gpFull_covarianceMatrix(this,error)
      type(gpFull), intent(inout) :: this
      integer, optional, intent(out) :: error

      integer :: i_coordinate, n_globalY, i_globalY, j_globalY, i_global_yPrime, j_global_yPrime, i_y, j_y, i_yPrime, j_yPrime, i_x, j_x, i_xPrime, j_xPrime, &
         first_nonzero, last_nonzero
      real(dp) :: covariance_x_x
      real(dp), dimension(:), allocatable :: globalY
      real(dp), dimension(:,:), allocatable :: grad_Covariance_i
      logical :: is_i_xPrime

      type(LA_matrix) :: LA_covariance_y_y

      INIT_ERROR(error)

      call system_timer('gpFull_covarianceMatrix')

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_covarianceMatrix: object not initialised',error)
      endif

      do i_coordinate = 1, this%n_coordinate
         if( .not. this%coordinate(i_coordinate)%initialised ) then
            RAISE_ERROR('gpFull_covarianceMatrix: '//i_coordinate//'th coordinate object not initialised',error)
         endif

         if(this%coordinate(i_coordinate)%bond_real_space) then
            if(.not. this%coordinate(i_coordinate)%bond_real_space_cov%initialised) then
               call gpCoordinates_gpCovariance_bond_real_space_Initialise(this%coordinate(i_coordinate))
            endif
         endif

         if(this%coordinate(i_coordinate)%atom_real_space) then
            if(.not. this%coordinate(i_coordinate)%atom_real_space_cov%initialised) then
               call gpCoordinates_gpCovariance_atom_real_space_Initialise(this%coordinate(i_coordinate))
            endif
         endif

      enddo

      n_globalY = this%n_y + this%n_yPrime

      i_globalY = 0
      do i_y = 1, this%n_y
         ! loop over all function values

         i_globalY = i_globalY + 1
         this%map_y_globalY(i_y) = i_globalY
      enddo
      do i_yPrime = 1, this%n_yPrime
         ! loop over all derivative values

         i_globalY = i_globalY + 1
         this%map_yPrime_globalY(i_yPrime) = i_globalY
      enddo

      call reallocate(this%covariance_y_y, n_globalY, n_globalY, zero = .true.)

      do i_coordinate = 1, this%n_coordinate
         ! loop over different descriptor types

!!$omp parallel default(none) private(j_x, j_y, j_globalY, i_x, i_y, i_globalY, covariance_x_x, i_xPrime, i_yPrime, i_global_yPrime, j_global_yPrime, j_xPrime, j_yPrime) shared(this,i_coordinate)
!!$omp do
         do j_x = 1, this%coordinate(i_coordinate)%n_x
            ! loop over all data

            j_y = this%coordinate(i_coordinate)%map_x_y(j_x)
            ! find which function value depends on the given descriptor
            if( j_y /= 0 ) then
               j_globalY = this%map_y_globalY(j_y)
               ! find unique function value/derivative identifier

               allocate( grad_Covariance_i(this%coordinate(i_coordinate)%d,this%coordinate(i_coordinate)%n_x) )
               do i_x = 1, this%coordinate(i_coordinate)%n_x
                  ! loop over all data

                  i_y = this%coordinate(i_coordinate)%map_x_y(i_x)
                  ! find which function value depends on the given descriptor

                  is_i_xPrime = any(this%coordinate(i_coordinate)%map_xPrime_x == i_x)

                  if( (i_y /= 0 .and. i_y <= j_y) .or. is_i_xPrime) then

                     if(is_i_xPrime) then
                        covariance_x_x = gpCoordinates_Covariance(this%coordinate(i_coordinate), i_x = i_x, j_x = j_x, grad_Covariance_i=grad_Covariance_i(:,i_x))
                     else
                        covariance_x_x = gpCoordinates_Covariance(this%coordinate(i_coordinate), i_x = i_x, j_x = j_x)
                     endif

                     if( i_y /= 0 ) then
                        i_globalY = this%map_y_globalY(i_y)
                        ! find unique function value/derivative identifier
                        this%covariance_y_y(i_globalY, j_globalY) = this%covariance_y_y(i_globalY, j_globalY) + covariance_x_x
                     endif
                  endif

               enddo

               do i_xPrime = 1, this%coordinate(i_coordinate)%n_xPrime
                  ! loop over all derivative data

                  i_yPrime = this%coordinate(i_coordinate)%map_xPrime_yPrime(i_xPrime)
                  ! find which derivative depends on the given descriptor

                  i_x = this%coordinate(i_coordinate)%map_xPrime_x(i_xPrime)

                  if( i_yPrime /= 0 ) then
                     i_global_yPrime = this%map_yPrime_globalY(i_yPrime)
                     ! find unique function value/derivative identifier

                     if(this%coordinate(i_coordinate)%atom_real_space) then
                        first_nonzero = this%coordinate(i_coordinate)%xPrime_size(i_xPrime)
                        last_nonzero = first_nonzero + size(this%coordinate(i_coordinate)%xPrime,1) - 1
                     else
                        first_nonzero = 1
                        last_nonzero = size(this%coordinate(i_coordinate)%xPrime,1)
                     endif

                     covariance_x_x = dot_product(grad_Covariance_i(first_nonzero:last_nonzero,i_x),this%coordinate(i_coordinate)%xPrime(:,i_xPrime))
                     !gpCoordinates_Covariance(this%coordinate(i_coordinate), i_xPrime = i_xPrime, j_x = j_x)
                     this%covariance_y_y(i_global_yPrime, j_globalY) = this%covariance_y_y(i_global_yPrime, j_globalY) + covariance_x_x
                  endif
               enddo
            endif
            if(allocated(grad_Covariance_i)) deallocate(grad_Covariance_i)
            
         enddo
!!$omp end do

!!$omp do
         do j_xPrime = 1, this%coordinate(i_coordinate)%n_xPrime
            ! loop over all derivative data

            j_yPrime = this%coordinate(i_coordinate)%map_xPrime_yPrime(j_xPrime)
            ! find which derivative depends on the given descriptor

            if( j_yPrime /= 0 ) then
               j_global_yPrime = this%map_yPrime_globalY(j_yPrime)
               ! find unique function value/derivative identifier

               do i_xPrime = 1, this%coordinate(i_coordinate)%n_xPrime
                  ! loop over all derivative data

                  i_yPrime = this%coordinate(i_coordinate)%map_xPrime_yPrime(i_xPrime)
                  ! find which derivative depends on the given descriptor

                  if( i_yPrime /= 0 .and. i_yPrime <= j_yPrime) then
                     i_global_yPrime = this%map_yPrime_globalY(i_yPrime)
                     ! find unique function value/derivative identifier

                     call system_abort('not implemented yet')
                     !covariance_x_x = gpCoordinates_Covariance(this%coordinate(i_coordinate), i_xPrime = i_xPrime, j_xPrime = j_xPrime)
                     this%covariance_y_y(i_global_yPrime, j_global_yPrime) = this%covariance_y_y(i_global_yPrime, j_global_yPrime) + covariance_x_x
                  endif
               enddo
            endif
         enddo
!!$omp end parallel
      enddo

      do j_y = 1, size(this%covariance_y_y,2)
         do i_y = j_y + 1, size(this%covariance_y_y,1)
            this%covariance_y_y(i_y,j_y) = this%covariance_y_y(j_y,i_y)
         enddo
      enddo

      allocate( globalY(n_globalY) )

      do i_y = 1, this%n_y
         ! loop over all function values

         i_globalY = this%map_y_globalY(i_y)
         ! find unique function value/derivative identifier

         this%covariance_y_y(i_globalY,i_globalY) = this%covariance_y_y(i_globalY,i_globalY) + &
            this%sigma_y(i_y)**2

         globalY(i_globalY) = this%y(i_y)
      enddo

      do i_yPrime = 1, this%n_yPrime
         ! loop over all function values

         i_global_yPrime = this%map_yPrime_globalY(i_yPrime)
         ! find unique function value/derivative identifier

         this%covariance_y_y(i_global_yPrime,i_global_yPrime) = this%covariance_y_y(i_global_yPrime,i_global_yPrime) + &
            this%sigma_yPrime(i_yPrime)**2

         globalY(i_global_yPrime) = this%y(i_yPrime)
      enddo

      call reallocate(this%alpha, n_globalY, zero = .true.)

      call initialise(LA_covariance_y_y,this%covariance_y_y)
      call Matrix_Solve(LA_covariance_y_y, globalY, this%alpha ,error=error)
      call finalise(LA_covariance_y_y)

      if(allocated(globalY)) deallocate(globalY)

      call system_timer('gpFull_covarianceMatrix')

   endsubroutine gpFull_covarianceMatrix

   function gpCoordinates_Covariance( this, i_x, j_x, i_sparseX, j_sparseX, grad_Covariance_i, grad_Covariance_j, grad2_Covariance, error )
      type(gpCoordinates), intent(inout), target :: this
      integer, intent(in), optional :: i_x, j_x, i_sparseX, j_sparseX
      real(dp), dimension(:), optional, intent(out), target :: grad_Covariance_i, grad_Covariance_j
      real(dp), dimension(:,:), optional, intent(out), target  :: grad2_Covariance
      integer, optional, intent(out) :: error

      real(dp) :: gpCoordinates_Covariance

      integer :: i_p, x_i_size, x_j_size, i
      real(dp) :: covarianceExp, covarianceDiag_x_x_i, covarianceDiag_x_x_j, normalisation, zeta, cov_delta
      real(dp), dimension(:), pointer :: x_i, x_j, xPrime_i, xPrime_j, grad_Covariance_Diag_i, grad_Covariance_Diag_j
      real(dp), dimension(:,:), pointer :: xPrime_ij
      real(dp), dimension(this%d) :: xI_xJ_theta2, xI_xJ

      INIT_ERROR(error)
      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_Covariance: object not initialised', error)
      endif

      if( count( (/ present(i_x), present(i_sparseX) /) ) /= 1 ) then
         RAISE_ERROR('gpCoordinates_Covariance: exactly one of i_x or i_sparseX can be present', error)
      endif

      if( count( (/ present(j_x), present(j_sparseX) /) ) /= 1 ) then
         RAISE_ERROR('gpCoordinates_Covariance: exactly one of j_x or j_sparseX can be present', error)
      endif

      x_i => null()
      x_j => null()
      grad_Covariance_Diag_i => null()
      grad_Covariance_Diag_j => null()
      xPrime_i => null()
      xPrime_j => null()
      xPrime_ij => null()

      x_i_size = 0
      x_j_size = 0

      if(present(i_x)) then
         x_i => this%x(:,i_x)
         covarianceDiag_x_x_i = this%covarianceDiag_x_x(i_x)
         grad_Covariance_Diag_i => this%covarianceDiag_x_xPrime(:,i_x)

         if(this%bond_real_space .or. this%atom_real_space) then
            x_i_size = this%x_size(i_x)
         endif
      endif

      if(present(j_x)) then
         x_j => this%x(:,j_x)
         covarianceDiag_x_x_j = this%covarianceDiag_x_x(j_x)
         grad_Covariance_Diag_j => this%covarianceDiag_x_xPrime(:,j_x)

         if(this%bond_real_space .or. this%atom_real_space) then
            x_j_size = this%x_size(j_x)
         endif
      endif

      if(present(i_sparseX)) then
         x_i => this%sparseX(:,i_sparseX)
         covarianceDiag_x_x_i = this%covarianceDiag_sparseX_sparseX(i_sparseX)

         if(this%bond_real_space .or. this%atom_real_space) then
            x_i_size = this%sparseX_size(i_sparseX)
         endif
      endif

      if(present(j_sparseX)) then
         x_j => this%sparseX(:,j_sparseX)
         covarianceDiag_x_x_j = this%covarianceDiag_sparseX_sparseX(j_sparseX)

         if(this%bond_real_space .or. this%atom_real_space) then
            x_j_size = this%sparseX_size(j_sparseX)
         endif
      endif

      if( .not. associated(x_i) .and. .not. associated(x_j) ) then
         RAISE_ERROR('gpCoordinates_Covariance: both i and j indices have to be present', error)
      endif

      gpCoordinates_Covariance = 0.0_dp
      if(present(grad_Covariance_i)) then
         grad_Covariance_i = 0.0_dp
         xPrime_i => grad_Covariance_i
      endif
      if(present(grad_Covariance_j)) then
         grad_Covariance_j = 0.0_dp
         xPrime_j => grad_Covariance_j
      endif
      if(present(grad2_Covariance)) then
         if(.not. associated(xPrime_i)) allocate(xPrime_i(size(x_i)))
         if(.not. associated(xPrime_j)) allocate(xPrime_j(size(x_j)))
         xPrime_i = 0.0_dp
         xPrime_j = 0.0_dp
         grad2_Covariance = 0.0_dp
         xPrime_ij => grad2_Covariance
      endif

      if(this%bond_real_space) then
         if(present(i_x)) then
            if(present(j_x)) then
               gpCoordinates_Covariance = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, x_i=this%x(:,i_x), x_i_size=this%x_size(i_x), &
               x_j=this%x(:,j_x), x_j_size=this%x_size(j_x)) / &
               sqrt(this%covarianceDiag_x_x(i_x) * this%covarianceDiag_x_x(j_x))
            elseif(present(j_sparseX)) then
               gpCoordinates_Covariance = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, x_i=this%x(:,i_x), x_i_size=this%x_size(i_x), &
               j_sparseX=j_sparseX) / sqrt(this%covarianceDiag_x_x(i_x) * this%covarianceDiag_sparseX_sparseX(j_sparseX))
!            elseif(present(j_xPrime)) then
            endif
         elseif(present(i_sparseX)) then
            if(present(j_x)) then
               gpCoordinates_Covariance = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, x_i=this%x(:,j_x), x_i_size=this%x_size(j_x), j_sparseX=i_sparseX) &
                                          / sqrt(this%covarianceDiag_sparseX_sparseX(i_sparseX) * this%covarianceDiag_x_x(j_x))
            elseif(present(j_sparseX)) then
               gpCoordinates_Covariance = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, i_sparseX=i_sparseX, j_sparseX=j_sparseX) &
                                          / sqrt(this%covarianceDiag_sparseX_sparseX(i_sparseX) * this%covarianceDiag_sparseX_sparseX(j_sparseX))
!            elseif(present(j_xPrime)) then
            endif
!         elseif(present(i_xPrime)) then
!            if(present(j_x)) then
!            elseif(present(j_sparseX)) then
!            elseif(present(j_xPrime)) then
!            endif
         endif
      elseif(this%atom_real_space) then
         zeta = this%atom_real_space_cov%zeta
         cov_delta = this%atom_real_space_cov%delta

         gpCoordinates_Covariance = gpCovariance_atom_real_space_Calc(this%atom_real_space_cov, x_i, x_i_size, x_j, x_j_size, &
         xPrime_i = xPrime_i, xPrime_j = xPrime_j, xPrime_ij = xPrime_ij)

         normalisation = sqrt(covarianceDiag_x_x_i * covarianceDiag_x_x_j)
         
         if(present(grad2_Covariance)) then
            !grad2_Covariance = xPrime_ij / normalisation - &
            !                   (xPrime_i .outer. grad_Covariance_Diag_j) / normalisation / covarianceDiag_x_x_j - &
            !                   (grad_Covariance_Diag_i .outer. xPrime_j) / normalisation / covarianceDiag_x_x_i + &
            !                   (grad_Covariance_Diag_i .outer. grad_Covariance_Diag_j) * gpCoordinates_Covariance / normalisation**3

            !grad2_Covariance = xPrime_ij / normalisation
            !do i = 1, size(grad2_Covariance,2)
            !   grad2_Covariance(:,i) = grad2_Covariance(:,i) - &
            !   xPrime_i*grad_Covariance_Diag_j(i) / normalisation / covarianceDiag_x_x_j - &
            !   grad_Covariance_Diag_i*xPrime_j(i) / normalisation / covarianceDiag_x_x_i + &
            !   grad_Covariance_Diag_i*grad_Covariance_Diag_j(i) * gpCoordinates_Covariance / normalisation**3
            !enddo
            grad2_Covariance = xPrime_ij * gpCoordinates_Covariance**(zeta-1.0_dp) / normalisation**zeta
            do i = 1, size(grad2_Covariance,2)
               grad2_Covariance(:,i) = grad2_Covariance(:,i) + &
               (zeta-1.0_dp) * gpCoordinates_Covariance**(zeta-2.0_dp) * xPrime_i*xPrime_j(i) / normalisation**zeta - &
               zeta * gpCoordinates_Covariance**(zeta-1.0_dp) * xPrime_i*grad_Covariance_Diag_j(i) / normalisation**zeta / covarianceDiag_x_x_j - &
               zeta * gpCoordinates_Covariance**(zeta-1.0_dp) * grad_Covariance_Diag_i*xPrime_j(i) / normalisation**zeta / covarianceDiag_x_x_i + &
               zeta * gpCoordinates_Covariance**zeta * grad_Covariance_Diag_i*grad_Covariance_Diag_j(i) / normalisation**(zeta+2.0_dp)
            enddo
            grad2_Covariance = zeta * grad2_Covariance * cov_delta**2
                              
         endif

         if(present(grad_Covariance_i)) then
            !grad_Covariance_i = grad_Covariance_i / normalisation - grad_Covariance_Diag_i * gpCoordinates_Covariance / normalisation / covarianceDiag_x_x_i
            grad_Covariance_i = gpCoordinates_Covariance**(zeta-1.0_dp) * grad_Covariance_i / normalisation**zeta - grad_Covariance_Diag_i * (gpCoordinates_Covariance / normalisation)**zeta / covarianceDiag_x_x_i
            grad_Covariance_i = zeta * grad_Covariance_i * cov_delta**2
         endif

         if(present(grad_Covariance_j)) then
            !grad_Covariance_j = grad_Covariance_j / normalisation - grad_Covariance_Diag_j * gpCoordinates_Covariance / normalisation / covarianceDiag_x_x_j
            grad_Covariance_j = gpCoordinates_Covariance**(zeta-1.0_dp) * grad_Covariance_j / normalisation**zeta - grad_Covariance_Diag_j * (gpCoordinates_Covariance / normalisation)**zeta / covarianceDiag_x_x_j
            grad_Covariance_j = zeta * grad_Covariance_j * cov_delta**2
         endif

         !gpCoordinates_Covariance = gpCoordinates_Covariance / normalisation
         gpCoordinates_Covariance = (gpCoordinates_Covariance / normalisation)**zeta * cov_delta**2

         if(.not. present(grad_Covariance_i) .and. associated(xPrime_i)) deallocate(xPrime_i)
         if(.not. present(grad_Covariance_j) .and. associated(xPrime_j)) deallocate(xPrime_j)

         xPrime_i => null()
         xPrime_j => null()
         xPrime_ij => null()
      elseif(this%soap) then
         gpCoordinates_Covariance = dot_product(x_i,x_j)

         if(present(grad_Covariance_i)) grad_Covariance_i = this%delta**2 * this%theta(1) * gpCoordinates_Covariance**(this%theta(1)-1.0_dp) * x_j
         if(present(grad_Covariance_j)) grad_Covariance_j = this%delta**2 * this%theta(1) * gpCoordinates_Covariance**(this%theta(1)-1.0_dp) * x_i

         gpCoordinates_Covariance = this%delta**2 * gpCoordinates_Covariance**this%theta(1)
      else
         do i_p = 1, this%n_permutations
            ! permute only i. theta should be symmetrised by now.

            xI_xJ = (x_i(this%permutations(:,i_p)) - x_j)
            xI_xJ_theta2 = xI_xJ / this%theta**2

            covarianceExp = this%delta**2 * exp( -0.5_dp * dot_product(xI_xJ_theta2,xI_xJ) )

            gpCoordinates_Covariance = gpCoordinates_Covariance + covarianceExp

            if(present(grad_Covariance_i)) then
               grad_Covariance_i(this%permutations(:,i_p)) = grad_Covariance_i(this%permutations(:,i_p)) - covarianceExp * xI_xJ_theta2
            endif

            if(present(grad_Covariance_j)) then
               grad_Covariance_j = grad_Covariance_j + covarianceExp * xI_xJ_theta2
            endif

            if(present(grad2_Covariance)) then
               do i = 1, this%d
                  grad2_Covariance(:,this%permutations(i,i_p)) = grad2_Covariance(:,this%permutations(i,i_p)) - covarianceExp * &
                  xI_xJ_theta2*xI_xJ_theta2(i)
                  grad2_Covariance(this%permutations(i,i_p),i) = grad2_Covariance(this%permutations(i,i_p),i) + covarianceExp/this%theta(i)**2
               enddo
            endif

            !if(present(i_xPrime) .and. .not. present(j_xPrime)) then

            !   gpCoordinates_Covariance = gpCoordinates_Covariance - covarianceExp * (dot_product(xI_xJ_theta,xPrime_i_theta(this%permutations(:,i_p)))*fc_i - dfc_i)*fc_j

            !elseif(.not. present(i_xPrime) .and. present(j_xPrime)) then

            !   gpCoordinates_Covariance = gpCoordinates_Covariance + covarianceExp * (dot_product(xI_xJ_theta,xPrime_j_theta)*fc_j + dfc_j)*fc_i

            !elseif(present(i_xPrime) .and. present(j_xPrime)) then

            !   gpCoordinates_Covariance = gpCoordinates_Covariance + covarianceExp * ( dot_product( xPrime_i_theta(this%permutations(:,i_p)), xPrime_j_theta )*fc_i*fc_j + &
            !   ( - dot_product( xI_xJ_theta, xPrime_i_theta(this%permutations(:,i_p)) )*fc_i + dfc_i ) * &
            !   ( dot_product( xI_xJ_theta, xPrime_j_theta )*fc_j + dfc_j ) )

            !else
            !   gpCoordinates_Covariance = gpCoordinates_Covariance + covarianceExp*fc_i*fc_j
            !endif
         enddo

      endif

      gpCoordinates_Covariance = gpCoordinates_Covariance + this%f0**2

   endfunction gpCoordinates_Covariance

   subroutine gpCoordinates_gpCovariance_bond_real_space_Initialise( this, error )
      type(gpCoordinates), intent(inout) :: this
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      integer :: i, j, k
      real(dp) :: r, r_x, x, y, z, buffer1, buffer2

      INIT_ERROR(error)

      if (.not. this%bond_real_space) then
         RAISE_ERROR('gpCoordinates_gpCovariance_bond_real_space_Initialise: covariance is not bond_real_space', error)
      endif

      call gpCovariance_bond_real_space_Finalise(this%bond_real_space_cov, error)

      call initialise(params)
      call param_register(params, 'n', '2', this%bond_real_space_cov%n, &
           help_string="Covariance degree for bond_real_space-type descriptors")
      call param_register(params, 'atom_sigma', '0.0', this%bond_real_space_cov%atom_sigma, &
           help_string="Atoms sigma for bond_real_space-type descriptors")

      if (.not. param_read_line(params, string(this%descriptor_str), ignore_unknown=.true., task='gpCoordinates_gpCovariance_bond_real_space_Initialise descriptor_str')) then
         RAISE_ERROR("gpCoordinates_gpCovariance_bond_real_space_Initialise failed to parse descriptor_str='"//trim(string(this%descriptor_str))//"'", error)
      endif
      call finalise(params)

      this%bond_real_space_cov%delta = this%delta

      if (this%sparsified) then
         allocate( this%bond_real_space_cov%sparseX_N(this%n_sparseX), &
                   this%bond_real_space_cov%sparseX_pos(0:1,3,maxval(this%sparseX_size) / 4,this%n_sparseX), &
                   this%bond_real_space_cov%sparseX_atom_covariance_cutoff(maxval(this%sparseX_size) / 4,this%n_sparseX), &
                   this%bond_real_space_cov%sparseX_bond_covariance_cutoff(this%n_sparseX), &
                   this%bond_real_space_cov%sparseX_self_overlap(this%n_sparseX) )

         this%bond_real_space_cov%sparseX_N = this%sparseX_size / 4

         do i = 1, this%n_sparseX
            do j = 1, this%bond_real_space_cov%sparseX_N(i)
               this%bond_real_space_cov%sparseX_pos(0,:,j,i) = this%sparseX((4 * (j - 1)) + 1:(4 * (j - 1)) + 3,i)
               this%bond_real_space_cov%sparseX_atom_covariance_cutoff(j,i) = this%sparseX(4 * j,i)
            enddo

            r = norm(this%bond_real_space_cov%sparseX_pos(0,:,1,i))
            r_x = norm(this%bond_real_space_cov%sparseX_pos(0,2:3,1,i))
            x = this%bond_real_space_cov%sparseX_pos(0,1,1,i)
            y = this%bond_real_space_cov%sparseX_pos(0,2,1,i)
            z = this%bond_real_space_cov%sparseX_pos(0,3,1,i)
            if(r_x .fne. 0.0_dp) then
               this%bond_real_space_cov%sparseX_pos(0,:,:,i) = matmul(reshape((/ r_x/r, 0.0_dp, x/r, &
                                                                                 -x*y/(r_x*r), z/r_x, y/r, &
                                                                                 -x*z/(r_x*r), -y/r_x, z/r /), (/ 3, 3 /)), &
                                                                      this%bond_real_space_cov%sparseX_pos(0,:,:,i))
            else
               this%bond_real_space_cov%sparseX_pos(0,:,:,i) = matmul(reshape((/ 0.0_dp, 0.0_dp, 1.0_dp, &
                                                                                 0.0_dp, 1.0_dp, 0.0_dp, &
                                                                                 -1.0_dp, 0.0_dp, 0.0_dp /), (/ 3, 3 /)), &
                                                                      this%bond_real_space_cov%sparseX_pos(0,:,:,i))
            endif

            this%bond_real_space_cov%sparseX_pos(1,:,:,i) = matmul(reshape((/ 1.0_dp, 0.0_dp, 0.0_dp, &
                                                                              0.0_dp, 1.0_dp, 0.0_dp, &
                                                                              0.0_dp, 0.0_dp, -1.0_dp /), (/ 3, 3 /)), &
                                                                   this%bond_real_space_cov%sparseX_pos(0,:,:,i))

            this%bond_real_space_cov%sparseX_bond_covariance_cutoff(i) = this%sparseX(0,i)

            buffer1 = 0.0_dp
            buffer2 = 0.0_dp

            do j = 1, this%bond_real_space_cov%sparseX_N(i)
               buffer1 = buffer1 + this%bond_real_space_cov%sparseX_atom_covariance_cutoff(j,i)**2
            enddo

            do j = 2, this%bond_real_space_cov%sparseX_N(i)
               do k = 1, (j - 1)
                  buffer2 = buffer2 + 2.0_dp * this%bond_real_space_cov%sparseX_atom_covariance_cutoff(j,i) * &
                                               this%bond_real_space_cov%sparseX_atom_covariance_cutoff(k,i) * &
                                               exp( -0.25_dp * normsq(this%bond_real_space_cov%sparseX_pos(0,:,j,i) - this%bond_real_space_cov%sparseX_pos(0,:,k,i)) / this%bond_real_space_cov%atom_sigma**2 )
               enddo
            enddo

            this%bond_real_space_cov%sparseX_self_overlap(i) = buffer1 + buffer2
         enddo

         this%bond_real_space_cov%sparseX_inisialised = .true.
      endif

      this%bond_real_space_cov%initialised = .true.

   endsubroutine gpCoordinates_gpCovariance_bond_real_space_Initialise

!   subroutine gpCoordinates_gpCovariance_bond_real_space_Initialise( this, error )
!      type(gpCoordinates), intent(inout) :: this
!      integer, optional, intent(out) :: error
!
!      type(Dictionary) :: params
!      integer :: i, j, k
!      real(dp) :: cos_z_angle, sin_z_angle, r, r_x, x, y, z, buffer1, buffer2
!
!      INIT_ERROR(error)
!
!      if (.not. this%bond_real_space) then
!         RAISE_ERROR('gpCoordinates_gpCovariance_bond_real_space_Initialise: covariance is not bond_real_space', error)
!      endif
!
!      call gpCovariance_bond_real_space_Finalise(this%bond_real_space_cov, error)
!
!      call initialise(params)
!      call param_register(params, 'integration_steps', '0', this%bond_real_space_cov%integration_steps, &
!           help_string="Number of integration steps when calculating covariance for bond_real_space-type descriptors")
!      call param_register(params, 'atom_sigma', '0.0', this%bond_real_space_cov%atom_sigma, &
!           help_string="Atoms sigma for bond_real_space-type descriptors")
!
!      if (.not. param_read_line(params, string(this%descriptor_str), ignore_unknown=.true., task='gpCoordinates_gpCovariance_bond_real_space_Initialise descriptor_str')) then
!         RAISE_ERROR("gpCoordinates_gpCovariance_bond_real_space_Initialise failed to parse descriptor_str='"//trim(string(this%descriptor_str))//"'", error)
!      endif
!      call finalise(params)
!
!      this%bond_real_space_cov%theta = this%theta(1)
!      this%bond_real_space_cov%delta = this%delta
!
!      allocate( this%bond_real_space_cov%rot_matrix(0:1,0:this%bond_real_space_cov%integration_steps - 1,3,3) )
!
!      do i = 0, (this%bond_real_space_cov%integration_steps - 1)
!         cos_z_angle = cos( 2.0_dp * PI * real(i, dp) / real(this%bond_real_space_cov%integration_steps, dp) )
!         sin_z_angle = sin( 2.0_dp * PI * real(i, dp) / real(this%bond_real_space_cov%integration_steps, dp) )
!
!         this%bond_real_space_cov%rot_matrix(0,i,:,:) = reshape((/ cos_z_angle, sin_z_angle, 0.0_dp, &
!                                                                   -sin_z_angle, cos_z_angle, 0.0_dp, &
!                                                                   0.0_dp, 0.0_dp, 1.0_dp /), (/ 3, 3 /))
!         this%bond_real_space_cov%rot_matrix(1,i,:,:) = reshape((/ cos_z_angle, sin_z_angle, 0.0_dp, &
!                                                                   -sin_z_angle, cos_z_angle, 0.0_dp, &
!                                                                   0.0_dp, 0.0_dp, -1.0_dp /), (/ 3, 3 /))
!      enddo
!
!      if (this%sparsified) then
!         allocate( this%bond_real_space_cov%sparseX_N(this%n_sparseX), &
!                   this%bond_real_space_cov%sparseX_pos_rot(0:1,0:this%bond_real_space_cov%integration_steps - 1,3,maxval(this%sparseX_size) / 4,this%n_sparseX), &
!                   this%bond_real_space_cov%sparseX_atom_covariance_cutoff(maxval(this%sparseX_size) / 4,this%n_sparseX), &
!                   this%bond_real_space_cov%sparseX_bond_covariance_cutoff(this%n_sparseX), &
!                   this%bond_real_space_cov%sparseX_invariant_buffer(this%n_sparseX) )
!
!         this%bond_real_space_cov%sparseX_N = this%sparseX_size / 4
!
!         do i = 1, this%n_sparseX
!            do j = 1, this%bond_real_space_cov%sparseX_N(i)
!               this%bond_real_space_cov%sparseX_pos_rot(0,0,:,j,i) = this%sparseX((4 * (j - 1)) + 1:(4 * (j - 1)) + 3,i)
!               this%bond_real_space_cov%sparseX_atom_covariance_cutoff(j,i) = this%sparseX(4 * j,i)
!            enddo
!
!            r = norm(this%bond_real_space_cov%sparseX_pos_rot(0,0,:,1,i))
!            r_x = norm(this%bond_real_space_cov%sparseX_pos_rot(0,0,2:3,1,i))
!            x = this%bond_real_space_cov%sparseX_pos_rot(0,0,1,1,i)
!            y = this%bond_real_space_cov%sparseX_pos_rot(0,0,2,1,i)
!            z = this%bond_real_space_cov%sparseX_pos_rot(0,0,3,1,i)
!            if(r_x .fne. 0.0_dp) then
!               this%bond_real_space_cov%sparseX_pos_rot(0,0,:,:,i) = matmul(reshape((/ r_x/r, 0.0_dp, x/r, &
!                                                                                       -x*y/(r_x*r), z/r_x, y/r, &
!                                                                                       -x*z/(r_x*r), -y/r_x, z/r /), (/ 3, 3 /)), &
!                                                                            this%bond_real_space_cov%sparseX_pos_rot(0,0,:,:,i))
!            else
!               this%bond_real_space_cov%sparseX_pos_rot(0,0,:,:,i) = matmul(reshape((/ 0.0_dp, 0.0_dp, 1.0_dp, &
!                                                                                       0.0_dp, 1.0_dp, 0.0_dp, &
!                                                                                       -1.0_dp, 0.0_dp, 0.0_dp /), (/ 3, 3 /)), &
!                                                                            this%bond_real_space_cov%sparseX_pos_rot(0,0,:,:,i))
!            endif
!
!            this%bond_real_space_cov%sparseX_pos_rot(1,0,:,:,i) = matmul(reshape((/ 1.0_dp, 0.0_dp, 0.0_dp, &
!                                                                                    0.0_dp, 1.0_dp, 0.0_dp, &
!                                                                                    0.0_dp, 0.0_dp, -1.0_dp /), (/ 3, 3 /)), &
!                                                                         this%bond_real_space_cov%sparseX_pos_rot(0,0,:,:,i))
!
!            do j = 1, (this%bond_real_space_cov%integration_steps - 1)
!               this%bond_real_space_cov%sparseX_pos_rot(0,j,:,:,i) = matmul(this%bond_real_space_cov%rot_matrix(0,j,:,:), this%bond_real_space_cov%sparseX_pos_rot(0,0,:,:,i))
!               this%bond_real_space_cov%sparseX_pos_rot(1,j,:,:,i) = matmul(this%bond_real_space_cov%rot_matrix(1,j,:,:), this%bond_real_space_cov%sparseX_pos_rot(0,0,:,:,i))
!            enddo
!
!            this%bond_real_space_cov%sparseX_bond_covariance_cutoff(i) = this%sparseX(0,i)
!
!            buffer1 = 0.0_dp
!            buffer2 = 0.0_dp
!
!            do j = 1, this%bond_real_space_cov%sparseX_N(i)
!               buffer1 = buffer1 + this%bond_real_space_cov%sparseX_atom_covariance_cutoff(j,i)**2
!            enddo
!
!            do j = 2, this%bond_real_space_cov%sparseX_N(i)
!               do k = 1, (j - 1)
!                  buffer2 = buffer2 + 2.0_dp * this%bond_real_space_cov%sparseX_atom_covariance_cutoff(j,i) * &
!                                               this%bond_real_space_cov%sparseX_atom_covariance_cutoff(k,i) * &
!                                               exp( -0.25_dp * normsq(this%bond_real_space_cov%sparseX_pos_rot(0,0,:,j,i) - this%bond_real_space_cov%sparseX_pos_rot(0,0,:,k,i)) / this%bond_real_space_cov%atom_sigma**2 )
!               enddo
!            enddo
!
!            this%bond_real_space_cov%sparseX_invariant_buffer(i) = buffer1 + buffer2
!         enddo
!
!         this%bond_real_space_cov%sparseX_inisialised = .true.
!      endif
!
!      this%bond_real_space_cov%initialised = .true.
!
!   endsubroutine gpCoordinates_gpCovariance_bond_real_space_Initialise

   subroutine gpCoordinates_gpCovariance_atom_real_space_Initialise( this, error )
      type(gpCoordinates), intent(inout) :: this
      integer, optional, intent(out) :: error

      type(Dictionary) :: params

!$omp barrier      
      INIT_ERROR(error)

      if (.not. this%atom_real_space) then
         RAISE_ERROR('gpCoordinates_gpCovariance_atom_real_space_Initialise: covariance is not atom_real_space', error)
      endif

      call gpCovariance_atom_real_space_Finalise(this%atom_real_space_cov, error)

      call initialise(params)
      call param_register(params, 'l_max', '0', this%atom_real_space_cov%l_max, &
           help_string="Cutoff of spherical harmonics expansion")
      call param_register(params, 'atom_sigma', '1.0', this%atom_real_space_cov%atom_sigma, &
           help_string="Atoms sigma for atom_real_space-type descriptors")
      call param_register(params, 'zeta', '1.0', this%atom_real_space_cov%zeta, &
      help_string="Exponent of covariance function")
      call param_register(params, 'delta', '1.0', this%atom_real_space_cov%delta, &
           help_string="delta for atom_real_space-type descriptors")
      call param_register(params, 'cutoff', '0.0', this%atom_real_space_cov%cutoff, &
           help_string="cutoff for atom_real_space-type descriptors")
      call param_register(params, 'cutoff_transition_width', '0.0', this%atom_real_space_cov%cutoff_transition_width, &
           help_string="cutoff_transition_width for atom_real_space-type descriptors")

      if (.not. param_read_line(params, string(this%descriptor_str), ignore_unknown=.true., task='gpCoordinates_gpCovariance_atom_real_space_Initialise descriptor_str')) then
         RAISE_ERROR("gpCoordinates_gpCovariance_bond_atom_space_Initialise failed to parse descriptor_str='"//trim(string(this%descriptor_str))//"'", error)
      endif
      call finalise(params)

      this%atom_real_space_cov%initialised = .true.

      this%atom_real_space_cov%atom_sigma = 0.5_dp / this%atom_real_space_cov%atom_sigma**2

   endsubroutine gpCoordinates_gpCovariance_atom_real_space_Initialise

   function gpCovariance_bond_real_space_Calc( this, x_i, x_i_size, x_j, x_j_size, i_sparseX, j_sparseX, error )
      type(gpCovariance_bond_real_space), intent(in) :: this
      real(dp), optional, intent(in) :: x_i(0:), x_j(0:)
      integer, optional, intent(in) :: x_i_size, x_j_size
      integer, optional, intent(in) :: i_sparseX, j_sparseX
      integer, optional, intent(out) :: error

      real(dp) :: gpCovariance_bond_real_space_Calc

      integer :: i, j, k, m, n
      integer, allocatable :: l(:,:)
      integer :: x_i_N, x_j_N
      real(dp), allocatable :: x_i_pos(:,:), x_j_pos(:,:,:), x_i_atom_covariance_cutoff(:), x_j_atom_covariance_cutoff(:)
      real(dp), allocatable :: ri2_rj2(:), rizrjz(:), rixrjx_riyrjy(:), riyrjx_rixrjy(:), RriRrj(:)
      real(dp), allocatable :: ri2_rj2_l(:), rizrjz_l(:), rixrjx_riyrjy_l(:), riyrjx_rixrjy_l(:), RriRrj_l(:)
      real(dp) :: x_i_bond_covariance_cutoff, x_j_bond_covariance_cutoff, x_i_self_overlap, x_j_self_overlap
      real(dp) :: r, r_x, x, y, z, buffer1, buffer2, C_iRj, exp_arg, exp_arg_mirror, bessel_arg

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCovariance_bond_real_space_Calc: object not initialised', error)
      endif

      if( count( (/ present(x_i), present(i_sparseX) /) ) /= 1 ) then
         RAISE_ERROR('gpCovariance_bond_real_space_Calc: exactly one of i_x, i_xPrime and i_sparseX and  has to be present', error)
      endif

      if( count( (/ present(x_j), present(j_sparseX) /) ) /= 1 ) then
         RAISE_ERROR('gpCovariance_bond_real_space_Calc: exactly one of j_x, j_xPrime and j_sparseX and  has to be present', error)
      endif

      if( present(x_i) ) then
         allocate( x_i_pos(3,x_i_size / 4), &
                   x_i_atom_covariance_cutoff(x_i_size / 4) )

         x_i_N = x_i_size / 4

         do i = 1, x_i_N
            x_i_pos(:,i) = x_i((4 * (i - 1)) + 1:(4 * (i - 1)) + 3)
            x_i_atom_covariance_cutoff(i) = x_i(4 * i)
         enddo

         r = norm(x_i_pos(:,1))
         r_x = norm(x_i_pos(2:3,1))
         x = x_i_pos(1,1)
         y = x_i_pos(2,1)
         z = x_i_pos(3,1)
         if(r_x .fne. 0.0_dp) then
            x_i_pos(:,:) = matmul(reshape((/ r_x/r, 0.0_dp, x/r, &
                                             -x*y/(r_x*r), z/r_x, y/r, &
                                             -x*z/(r_x*r), -y/r_x, z/r /), (/ 3, 3 /)), &
                                  x_i_pos(:,:))
         else
            x_i_pos(:,:) = matmul(reshape((/ 0.0_dp, 0.0_dp, 1.0_dp, &
                                             0.0_dp, 1.0_dp, 0.0_dp, &
                                             -1.0_dp, 0.0_dp, 0.0_dp /), (/ 3, 3 /)), &
                                  x_i_pos(:,:))
         endif

         x_i_bond_covariance_cutoff = x_i(0)

         buffer1 = 0.0_dp
         buffer2 = 0.0_dp

         do i = 1, x_i_N
            buffer1 = buffer1 + x_i_atom_covariance_cutoff(i)**2
         enddo

         do i = 2, x_i_N
            do j = 1, (i - 1)
               buffer2 = buffer2 + 2.0_dp * x_i_atom_covariance_cutoff(i) * &
                                            x_i_atom_covariance_cutoff(j) * &
                                            exp( -0.25_dp * normsq(x_i_pos(:,i) - x_i_pos(:,j)) / this%atom_sigma**2 )
            enddo
         enddo

         x_i_self_overlap = buffer1 + buffer2
      endif

      if( present(x_j) ) then
         allocate( x_j_pos(0:1,3,x_j_size / 4), &
                   x_j_atom_covariance_cutoff(x_j_size / 4) )

         x_j_N = x_j_size / 4

         do i = 1, x_j_N
            x_j_pos(0,:,i) = x_j((4 * (i - 1)) + 1:(4 * (i - 1)) + 3)
            x_j_atom_covariance_cutoff(i) = x_j(4 * i)
         enddo

         r = norm(x_j_pos(0,:,1))
         r_x = norm(x_j_pos(0,2:3,1))
         x = x_j_pos(0,1,1)
         y = x_j_pos(0,2,1)
         z = x_j_pos(0,3,1)
         if(r_x .fne. 0.0_dp) then
            x_j_pos(0,:,:) = matmul(reshape((/ r_x/r, 0.0_dp, x/r, &
                                               -x*y/(r_x*r), z/r_x, y/r, &
                                               -x*z/(r_x*r), -y/r_x, z/r /), (/ 3, 3 /)), &
                                    x_j_pos(0,:,:))
         else
            x_j_pos(0,:,:) = matmul(reshape((/ 0.0_dp, 0.0_dp, 1.0_dp, &
                                               0.0_dp, 1.0_dp, 0.0_dp, &
                                               -1.0_dp, 0.0_dp, 0.0_dp /), (/ 3, 3 /)), &
                                    x_j_pos(0,:,:))
         endif

         x_j_pos(1,:,:) = matmul(reshape((/ 1.0_dp, 0.0_dp, 0.0_dp, &
                                            0.0_dp, 1.0_dp, 0.0_dp, &
                                            0.0_dp, 0.0_dp, -1.0_dp /), (/ 3, 3 /)), &
                                 x_j_pos(0,:,:))

         x_j_bond_covariance_cutoff = x_j(0)

         buffer1 = 0.0_dp
         buffer2 = 0.0_dp

         do i = 1, x_j_N
            buffer1 = buffer1 + x_j_atom_covariance_cutoff(i)**2
         enddo

         do i = 2, x_j_N
            do j = 1, (i - 1)
               buffer2 = buffer2 + 2.0_dp * x_j_atom_covariance_cutoff(i) * &
                                            x_j_atom_covariance_cutoff(j) * &
                                            exp( -0.25_dp * normsq(x_j_pos(0,:,i) - x_j_pos(0,:,j)) / this%atom_sigma**2 )
            enddo
         enddo

         x_j_self_overlap = buffer1 + buffer2
      endif

      if( present(i_sparseX) ) then
         allocate( x_i_pos(3,this%sparseX_N(i_sparseX)), &
                   x_i_atom_covariance_cutoff(this%sparseX_N(i_sparseX)) )

         x_i_N = this%sparseX_N(i_sparseX)
         x_i_pos = this%sparseX_pos(0,:,:,i_sparseX)
         x_i_atom_covariance_cutoff = this%sparseX_atom_covariance_cutoff(:,i_sparseX)
         x_i_bond_covariance_cutoff = this%sparseX_bond_covariance_cutoff(i_sparseX)
         x_i_self_overlap = this%sparseX_self_overlap(i_sparseX)
      endif

      if( present(j_sparseX) ) then
         allocate( x_j_pos(0:1,3,this%sparseX_N(j_sparseX)), &
                   x_j_atom_covariance_cutoff(this%sparseX_N(j_sparseX)) )

         x_j_N = this%sparseX_N(j_sparseX)
         x_j_pos = this%sparseX_pos(:,:,:,j_sparseX)
         x_j_atom_covariance_cutoff = this%sparseX_atom_covariance_cutoff(:,j_sparseX)
         x_j_bond_covariance_cutoff = this%sparseX_bond_covariance_cutoff(j_sparseX)
         x_j_self_overlap = this%sparseX_self_overlap(j_sparseX)
      endif

      gpCovariance_bond_real_space_Calc = 0.0_dp

      if( this%n == 1 ) then
         n = x_i_N * x_j_N
      elseif( this%n == 2 ) then
         n = x_i_N * x_j_N * (x_i_N * x_j_N + 1) / 2
      else
         RAISE_ERROR('gpCovariance_bond_real_space_Calc: covariance degree n > 2 not implemented', error)
      endif

      allocate(l(x_i_N * x_j_N,n))
      allocate(ri2_rj2(x_i_N * x_j_N), rizrjz(x_i_N * x_j_N), rixrjx_riyrjy(x_i_N * x_j_N), riyrjx_rixrjy(x_i_N * x_j_N), RriRrj(x_i_N * x_j_N))
      allocate(ri2_rj2_l(x_i_N * x_j_N), rizrjz_l(x_i_N * x_j_N), rixrjx_riyrjy_l(x_i_N * x_j_N), riyrjx_rixrjy_l(x_i_N * x_j_N), RriRrj_l(x_i_N * x_j_N))

      l = 0

      if( this%n == 1 ) then
         do n = 1, x_i_N * x_j_N
            l(n,n) = 1
         enddo
      elseif( this%n == 2 ) then
         i = 1
         j = 2

         do n = 1, (x_i_N * x_j_N * (x_i_N * x_j_N + 1) / 2)
            if( n <= (x_i_N * x_j_N) ) then
               l(n,n) = 2
            else
               l(i,n) = 1
               l(j,n) = 1

               if( j < (x_i_N * x_j_N) ) then
                  j = j + 1
               else
                  i = i + 1
                  j = i + 1
               endif
            endif
         enddo
      else
         RAISE_ERROR('gpCovariance_bond_real_space_Calc: covariance degree n > 2 not implemented', error)
      endif

         ri2_rj2 = 0.0_dp
         rizrjz = 0.0_dp
         rixrjx_riyrjy = 0.0_dp
         riyrjx_rixrjy = 0.0_dp
         RriRrj = 1.0_dp

         do k = 1, (x_i_N * x_j_N)
            i = 1 + ((k - 1) / x_j_N)
            j = 1 + mod(k - 1,x_j_N)

            ri2_rj2(k) = normsq(x_i_pos(:,i)) + normsq(x_j_pos(m,:,j))
            rizrjz(k) = x_i_pos(3,i) * x_j_pos(m,3,j)
            rixrjx_riyrjy(k) = (x_i_pos(1,i) * x_j_pos(m,1,j)) + (x_i_pos(2,i) * x_j_pos(m,2,j))
            riyrjx_rixrjy(k) = (x_i_pos(2,i) * x_j_pos(m,1,j)) - (x_i_pos(1,i) * x_j_pos(m,2,j))
            RriRrj(k) = x_i_atom_covariance_cutoff(i) * x_j_atom_covariance_cutoff(j)
         enddo

         do n = 1, size(l, 2)
            ri2_rj2_l = 0.0_dp
            rizrjz_l = 0.0_dp
            rixrjx_riyrjy_l = 0.0_dp
            riyrjx_rixrjy_l = 0.0_dp
            RriRrj_l = 1.0_dp

            C_iRj = factorial(this%n)

            do k = 1, (x_i_N * x_j_N)
               if( l(k,n) /= 0 ) then
                  ri2_rj2_l(k) = real(l(k,n), dp) * ri2_rj2(k)
                  rizrjz_l(k) = real(l(k,n), dp) * rizrjz(k)
                  rixrjx_riyrjy_l(k) = real(l(k,n), dp) * rixrjx_riyrjy(k)
                  riyrjx_rixrjy_l(k) = real(l(k,n), dp) * riyrjx_rixrjy(k)
                  RriRrj_l(k) = RriRrj(k)**l(k,n)

                  C_iRj = C_iRj / factorial(l(k,n))
               endif
            enddo

            exp_arg = -0.25_dp * (sum(ri2_rj2_l) - (2.0_dp * sum(rizrjz_l))) / this%atom_sigma**2
            exp_arg_mirror = -0.25_dp * (sum(ri2_rj2_l) + (2.0_dp * sum(rizrjz_l))) / this%atom_sigma**2
            bessel_arg = 0.5_dp * sqrt(sum(rixrjx_riyrjy_l)**2 + sum(riyrjx_rixrjy_l)**2) / this%atom_sigma**2

            if( bessel_arg < besseli_max_x ) then
               C_iRj = C_iRj * product(RriRrj_l) * (exp(exp_arg) + exp(exp_arg_mirror)) * besseli0(bessel_arg)
            else
               buffer1 = 1.0_dp
               do i = 1, besseli_max_n
                  buffer1 = buffer1 + besseli0_c(i)/bessel_arg**i
               enddo

               C_iRj = C_iRj * product(RriRrj_l) * (exp(exp_arg + bessel_arg) + exp(exp_arg_mirror + bessel_arg)) * buffer1 / sqrt(2.0_dp * pi * bessel_arg)
            endif

            gpCovariance_bond_real_space_Calc = gpCovariance_bond_real_space_Calc + C_iRj
         enddo

      if( (x_i_self_overlap + x_j_self_overlap) .fne. 0.0_dp ) then
         gpCovariance_bond_real_space_Calc = gpCovariance_bond_real_space_Calc / (0.5_dp * (x_i_self_overlap + x_j_self_overlap))**this%n
         gpCovariance_bond_real_space_Calc = gpCovariance_bond_real_space_Calc * this%delta**2 * x_i_bond_covariance_cutoff * x_j_bond_covariance_cutoff / 2.0_dp
      else
         gpCovariance_bond_real_space_Calc = 0.0_dp
      endif

      deallocate(x_i_pos, x_j_pos, x_i_atom_covariance_cutoff, x_j_atom_covariance_cutoff, &
                 l, ri2_rj2, rizrjz, rixrjx_riyrjy, riyrjx_rixrjy, RriRrj, &
                 ri2_rj2_l, rizrjz_l, rixrjx_riyrjy_l, riyrjx_rixrjy_l, RriRrj_l)

   endfunction gpCovariance_bond_real_space_Calc

    function besseli0(x)

       real(dp), intent(in) :: x
       real(dp) :: besseli0

       real(dp) :: x2, r, k
       integer :: i

       x2 = x**2

       if(x == 0.0_dp) then
          besseli0 = 1.0_dp
       elseif( x < besseli_max_x ) then
          besseli0 = 1.0_dp
          r = 1.0_dp
          k = 1.0_dp
          do while ( abs(r/besseli0) > NUMERICAL_ZERO )
             r = 0.25_dp * r * x2 / k**2
             besseli0 = besseli0 + r
             k = k + 1.0_dp
          enddo
       else
          besseli0 = 1.0_dp
          do i = 1, besseli_max_n
             besseli0 = besseli0 + besseli0_c(i)/x**i
          enddo
          besseli0 = besseli0 * exp(x) / sqrt(2.0_dp*pi*x)
       endif

    endfunction besseli0

    function besseli1(x)

       real(dp), intent(in) :: x
       real(dp) :: besseli1

       real(dp) :: x2, r, k
       integer :: i

       x2 = x**2

       if(x == 0.0_dp) then
          besseli1 = 0.0_dp
       elseif( x < besseli_max_x ) then
          besseli1 = 1.0_dp
          r = 1.0_dp
          k = 1.0_dp
          do while ( abs(r/besseli1) > NUMERICAL_ZERO )
             r = 0.25_dp * r * x2 / (k*(k+1.0_dp))
             besseli1 = besseli1 + r
             k = k + 1.0_dp
          enddo
          besseli1 = besseli1 * 0.5_dp * x
       else
          besseli1 = 1.0_dp
          do i = 1, besseli_max_n
             besseli1 = besseli1 + besseli1_c(i)/x**i
          enddo
          besseli1 = besseli1 * exp(x) / sqrt(2.0_dp*pi*x)
       endif

    endfunction besseli1

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

!   function gpCovariance_bond_real_space_Calc( this, x_i, x_i_size, x_j, x_j_size, i_sparseX, j_sparseX, error )
!      type(gpCovariance_bond_real_space), intent(in) :: this
!      real(dp), optional, intent(in) :: x_i(0:), x_j(0:)
!      integer, optional, intent(in) :: x_i_size, x_j_size
!      integer, optional, intent(in) :: i_sparseX, j_sparseX
!      integer, optional, intent(out) :: error
!
!      real(dp) :: gpCovariance_bond_real_space_Calc
!
!      integer :: i, j, k, l
!      integer :: x_i_N, x_j_N
!      real(dp), allocatable :: x_i_pos_rot(:,:), x_j_pos_rot(:,:,:,:), x_i_atom_covariance_cutoff(:), x_j_atom_covariance_cutoff(:)
!      real(dp) :: x_i_bond_covariance_cutoff, x_j_bond_covariance_cutoff, x_i_invariant_buffer, x_j_invariant_buffer
!      real(dp) :: r, r_x, x, y, z, buffer1, buffer2, xI_xJ_2
!
!      INIT_ERROR(error)
!
!      if( .not. this%initialised ) then
!         RAISE_ERROR('gpCovariance_bond_real_space_Calc: object not initialised', error)
!      endif
!
!      if( count( (/ present(x_i), present(i_sparseX) /) ) /= 1 ) then
!         RAISE_ERROR('gpCovariance_bond_real_space_Calc: exactly one of i_x, i_xPrime and i_sparseX and  has to be present', error)
!      endif
!
!      if( count( (/ present(x_j), present(j_sparseX) /) ) /= 1 ) then
!         RAISE_ERROR('gpCovariance_bond_real_space_Calc: exactly one of j_x, j_xPrime and j_sparseX and  has to be present', error)
!      endif
!
!      if( present(x_i) ) then
!         allocate( x_i_pos_rot(3,x_i_size / 4), &
!                   x_i_atom_covariance_cutoff(x_i_size / 4) )
!
!         x_i_N = x_i_size / 4
!
!         do i = 1, x_i_N
!            x_i_pos_rot(:,i) = x_i((4 * (i - 1)) + 1:(4 * (i - 1)) + 3)
!            x_i_atom_covariance_cutoff(i) = x_i(4 * i)
!         enddo
!
!         r = norm(x_i_pos_rot(:,1))
!         r_x = norm(x_i_pos_rot(2:3,1))
!         x = x_i_pos_rot(1,1)
!         y = x_i_pos_rot(2,1)
!         z = x_i_pos_rot(3,1)
!         if(r_x .fne. 0.0_dp) then
!            x_i_pos_rot(:,:) = matmul(reshape((/ r_x/r, 0.0_dp, x/r, &
!                                                 -x*y/(r_x*r), z/r_x, y/r, &
!                                                 -x*z/(r_x*r), -y/r_x, z/r /), (/ 3, 3 /)), &
!                                      x_i_pos_rot(:,:))
!         else
!            x_i_pos_rot(:,:) = matmul(reshape((/ 0.0_dp, 0.0_dp, 1.0_dp, &
!                                                 0.0_dp, 1.0_dp, 0.0_dp, &
!                                                 -1.0_dp, 0.0_dp, 0.0_dp /), (/ 3, 3 /)), &
!                                      x_i_pos_rot(:,:))
!         endif
!
!         x_i_bond_covariance_cutoff = x_i(0)
!
!         buffer1 = 0.0_dp
!         buffer2 = 0.0_dp
!
!         do i = 1, x_i_N
!            buffer1 = buffer1 + x_i_atom_covariance_cutoff(i)**2
!         enddo
!
!         do i = 2, x_i_N
!            do j = 1, (i - 1)
!               buffer2 = buffer2 + 2.0_dp * x_i_atom_covariance_cutoff(i) * &
!                                            x_i_atom_covariance_cutoff(j) * &
!                                            exp( -0.25_dp * normsq(x_i_pos_rot(:,i) - x_i_pos_rot(:,j)) / this%atom_sigma**2 )
!            enddo
!         enddo
!
!         x_i_invariant_buffer = buffer1 + buffer2
!      endif
!
!      if( present(x_j) ) then
!         allocate( x_j_pos_rot(0:1,0:this%integration_steps - 1,3,x_j_size / 4), &
!                   x_j_atom_covariance_cutoff(x_j_size / 4) )
!
!         x_j_N = x_j_size / 4
!
!         do i = 1, x_j_N
!            x_j_pos_rot(0,0,:,i) = x_j((4 * (i - 1)) + 1:(4 * (i - 1)) + 3)
!            x_j_atom_covariance_cutoff(i) = x_j(4 * i)
!         enddo
!
!         r = norm(x_j_pos_rot(0,0,:,1))
!         r_x = norm(x_j_pos_rot(0,0,2:3,1))
!         x = x_j_pos_rot(0,0,1,1)
!         y = x_j_pos_rot(0,0,2,1)
!         z = x_j_pos_rot(0,0,3,1)
!         if(r_x .fne. 0.0_dp) then
!            x_j_pos_rot(0,0,:,:) = matmul(reshape((/ r_x/r, 0.0_dp, x/r, &
!                                                     -x*y/(r_x*r), z/r_x, y/r, &
!                                                     -x*z/(r_x*r), -y/r_x, z/r /), (/ 3, 3 /)), &
!                                          x_j_pos_rot(0,0,:,:))
!         else
!            x_j_pos_rot(0,0,:,:) = matmul(reshape((/ 0.0_dp, 0.0_dp, 1.0_dp, &
!                                                     0.0_dp, 1.0_dp, 0.0_dp, &
!                                                     -1.0_dp, 0.0_dp, 0.0_dp /), (/ 3, 3 /)), &
!                                          x_j_pos_rot(0,0,:,:))
!         endif
!
!         x_j_pos_rot(1,0,:,:) = matmul(reshape((/ 1.0_dp, 0.0_dp, 0.0_dp, &
!                                                  0.0_dp, 1.0_dp, 0.0_dp, &
!                                                  0.0_dp, 0.0_dp, -1.0_dp /), (/ 3, 3 /)), &
!                                       x_j_pos_rot(0,0,:,:))
!
!         do i = 1, (this%integration_steps - 1)
!            x_j_pos_rot(0,i,:,:) = matmul(this%rot_matrix(0,i,:,:), x_j_pos_rot(0,0,:,:))
!            x_j_pos_rot(1,i,:,:) = matmul(this%rot_matrix(1,i,:,:), x_j_pos_rot(0,0,:,:))
!         enddo
!
!         x_j_bond_covariance_cutoff = x_j(0)
!
!         buffer1 = 0.0_dp
!         buffer2 = 0.0_dp
!
!         do i = 1, x_j_N
!            buffer1 = buffer1 + x_j_atom_covariance_cutoff(i)**2
!         enddo
!
!         do i = 2, x_j_N
!            do j = 1, (i - 1)
!               buffer2 = buffer2 + 2.0_dp * x_j_atom_covariance_cutoff(i) * &
!                                            x_j_atom_covariance_cutoff(j) * &
!                                            exp( -0.25_dp * normsq(x_j_pos_rot(0,0,:,i) - x_j_pos_rot(0,0,:,j)) / this%atom_sigma**2 )
!            enddo
!         enddo
!
!         x_j_invariant_buffer = buffer1 + buffer2
!      endif
!
!      if( present(i_sparseX) ) then
!         allocate( x_i_pos_rot(3,this%sparseX_N(i_sparseX)), &
!                   x_i_atom_covariance_cutoff(this%sparseX_N(i_sparseX)) )
!
!         x_i_N = this%sparseX_N(i_sparseX)
!         x_i_pos_rot = this%sparseX_pos_rot(0,0,:,:,i_sparseX)
!         x_i_atom_covariance_cutoff = this%sparseX_atom_covariance_cutoff(:,i_sparseX)
!         x_i_bond_covariance_cutoff = this%sparseX_bond_covariance_cutoff(i_sparseX)
!         x_i_invariant_buffer = this%sparseX_invariant_buffer(i_sparseX)
!      endif
!
!      if( present(j_sparseX) ) then
!         allocate( x_j_pos_rot(0:1,0:this%integration_steps,3,this%sparseX_N(j_sparseX)), &
!                   x_j_atom_covariance_cutoff(this%sparseX_N(j_sparseX)) )
!
!         x_j_N = this%sparseX_N(j_sparseX)
!         x_j_pos_rot = this%sparseX_pos_rot(:,:,:,:,j_sparseX)
!         x_j_atom_covariance_cutoff = this%sparseX_atom_covariance_cutoff(:,j_sparseX)
!         x_j_bond_covariance_cutoff = this%sparseX_bond_covariance_cutoff(j_sparseX)
!         x_j_invariant_buffer = this%sparseX_invariant_buffer(j_sparseX)
!      endif
!
!      gpCovariance_bond_real_space_Calc = 0.0_dp
!
!      do i = 0, (this%integration_steps - 1)
!         do j = 0, 1
!            xI_xJ_2 = 0.0_dp
!
!            do k = 1, x_i_N
!               do l = 1, x_j_N
!                  xI_xJ_2 = xI_xJ_2 - 2.0_dp * x_i_atom_covariance_cutoff(k) * &
!                                               x_j_atom_covariance_cutoff(l) * &
!                                               exp( -0.25_dp * normsq(x_i_pos_rot(:,k) - x_j_pos_rot(j,i,:,l)) / this%atom_sigma**2 )
!               enddo
!            enddo
!
!            xI_xJ_2 = xI_xJ_2 + x_i_invariant_buffer + x_j_invariant_buffer
!            xI_xJ_2 = xI_xJ_2 * this%atom_sigma**3 * PI**(3.0_dp/2.0_dp)
!
!            gpCovariance_bond_real_space_Calc = gpCovariance_bond_real_space_Calc + exp( -0.5_dp * xI_xJ_2 / this%theta**2 )
!         enddo
!      enddo
!
!      gpCovariance_bond_real_space_Calc = gpCovariance_bond_real_space_Calc * this%delta**2 * &
!                                          x_i_bond_covariance_cutoff * x_j_bond_covariance_cutoff / &
!                                          (2.0_dp * real(this%integration_steps, dp))
!
!      deallocate(x_i_pos_rot, x_j_pos_rot, x_i_atom_covariance_cutoff, x_j_atom_covariance_cutoff)
!
!   endfunction gpCovariance_bond_real_space_Calc

   function gpCovariance_atom_real_space_Calc( this, x_i, x_i_size, x_j, x_j_size, xPrime_i, xPrime_j, xPrime_ij, error )
      type(gpCovariance_atom_real_space), intent(in) :: this
      real(dp), dimension(:), intent(in) :: x_i, x_j
      integer, intent(in) :: x_i_size, x_j_size
      real(dp), dimension(:), intent(out), optional, pointer :: xPrime_i, xPrime_j
      real(dp), dimension(:,:), intent(out), optional, pointer :: xPrime_ij
      integer, intent(out), optional :: error

      real(dp) :: gpCovariance_atom_real_space_Calc

      type(neighbour_descriptor), dimension(:), allocatable :: neighbour_i, neighbour_j, grad_spherical_i, grad_spherical_j
      type(neighbour_descriptor), dimension(:,:), allocatable :: grad_spherical_i_radial_j, grad_spherical_j_radial_i
      
      real(dp) :: r1, r2, arg_bess, fac_exp, mo_spher_bess_fi_ki_lmm, mo_spher_bess_fi_ki_lm, &
      mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lp, mo_spher_bess_fi_ki_lpp, &
      grad_mo_spher_bess_fi_ki_l, grad2_mo_spher_bess_fi_ki_l, &
      grad_arg_bess1, grad_arg_bess2, grad_fac_exp1, grad_fac_exp2, radial, grad_radial_i, grad_radial_j, grad2_radial_ij, &
      fcut1, fcut2, dfcut1, dfcut2, fac_r1r2, grad_fac_r1r2_1, grad_fac_r1r2_2, grad2_fac_exp, grad2_fac_r1r2

      real(dp), dimension(1) :: real_mould

      integer :: i, j, i_data, j_data, n1, n2, l, l1, l2, m1, m2, n_neighbour_i, n_neighbour_j, real_mould_size

      logical :: do_derivative, do_xPrime_i, do_xPrime_j, do_xPrime_ij

      complex(dp) :: I_lm1m2, tmp_complex
      type(cplx_2d_array), dimension(:), allocatable :: integral_r

      type grad_r_type
         type(cplx_2d_array), dimension(:), allocatable :: integral_r
      endtype grad_r_type

      type real_1d_array
         real(dp), dimension(:), allocatable :: value
      endtype real_1d_array

      type(grad_r_type), dimension(:), allocatable :: grad_ri, grad_rj
      type(grad_r_type), dimension(:,:), allocatable :: grad_rij

      type(real_1d_array), dimension(:,:), allocatable :: grad_spherical_ij

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCovariance_atom_real_space_Calc: object not initialised', error)
      endif

      do_xPrime_i = .false.
      do_xPrime_j = .false.
      do_xPrime_ij = .false.
      if(present(xPrime_i)) do_xPrime_i = associated(xPrime_i)
      if(present(xPrime_j)) do_xPrime_j = associated(xPrime_j)
      if(present(xPrime_ij)) do_xPrime_ij = associated(xPrime_ij)

      do_derivative = (do_xPrime_i .or. do_xPrime_j .or. do_xPrime_ij)

      call gpRealArray_NeighbourDescriptor(this,x_i,x_i_size,neighbour_i,n_neighbour_i)
      call gpRealArray_NeighbourDescriptor(this,x_j,x_j_size,neighbour_j,n_neighbour_j)

      if(do_xPrime_i .or. do_xPrime_ij) then
         allocate(grad_spherical_i(n_neighbour_i))
         allocate(grad_ri(n_neighbour_i))

         do i = 1, n_neighbour_i
            

            allocate(grad_ri(i)%integral_r(0:this%l_max))
            allocate(grad_spherical_i(i)%spherical_harmonics(0:this%l_max))

            do l = 0, this%l_max
               allocate(grad_spherical_i(i)%spherical_harmonics(l)%value(-l:l))
               grad_spherical_i(i)%spherical_harmonics(l)%value = CPLX_ZERO

               allocate(grad_ri(i)%integral_r(l)%value(-l:l,-l:l))
               grad_ri(i)%integral_r(l)%value = CPLX_ZERO
            enddo
         enddo
      endif

      if(do_xPrime_j .or. do_xPrime_ij) then
         allocate(grad_spherical_j(n_neighbour_j))
         allocate(grad_rj(n_neighbour_j))

         do i = 1, n_neighbour_j
            
            allocate(grad_rj(i)%integral_r(0:this%l_max))
            allocate(grad_spherical_j(i)%spherical_harmonics(0:this%l_max))

            do l = 0, this%l_max
               allocate(grad_spherical_j(i)%spherical_harmonics(l)%value(-l:l))
               grad_spherical_j(i)%spherical_harmonics(l)%value = CPLX_ZERO

               allocate(grad_rj(i)%integral_r(l)%value(-l:l,-l:l))
               grad_rj(i)%integral_r(l)%value = CPLX_ZERO
            enddo
         enddo
      endif

      if(do_xPrime_ij) then
         allocate(grad_rij(n_neighbour_j,n_neighbour_i))
         allocate(grad_spherical_ij(n_neighbour_j,n_neighbour_i))
         allocate(grad_spherical_i_radial_j(n_neighbour_j,n_neighbour_i))
         allocate(grad_spherical_j_radial_i(n_neighbour_j,n_neighbour_i))

         do i = 1, n_neighbour_i
            do j = 1, n_neighbour_j

               allocate(grad_spherical_ij(j,i)%value(0:this%l_max))
               grad_spherical_ij(j,i)%value = 0.0_dp
            
               allocate(grad_rij(j,i)%integral_r(0:this%l_max))

               allocate(grad_spherical_i_radial_j(j,i)%spherical_harmonics(0:this%l_max))
               allocate(grad_spherical_j_radial_i(j,i)%spherical_harmonics(0:this%l_max))

               do l = 0, this%l_max
                  allocate(grad_rij(j,i)%integral_r(l)%value(-l:l,-l:l))
                  grad_rij(j,i)%integral_r(l)%value = CPLX_ZERO

                  allocate(grad_spherical_i_radial_j(j,i)%spherical_harmonics(l)%value(-l:l))
                  grad_spherical_i_radial_j(j,i)%spherical_harmonics(l)%value = CPLX_ZERO

                  allocate(grad_spherical_j_radial_i(j,i)%spherical_harmonics(l)%value(-l:l))
                  grad_spherical_j_radial_i(j,i)%spherical_harmonics(l)%value = CPLX_ZERO
               enddo
            enddo
         enddo

      endif

      allocate(integral_r(0:this%l_max))

      do l = 0, this%l_max
         allocate(integral_r(l)%value(-l:l,-l:l))
         integral_r(l)%value = CPLX_ZERO
      enddo

      if(do_xPrime_i) xPrime_i = 0.0_dp
      if(do_xPrime_j) xPrime_j = 0.0_dp
      if(do_xPrime_ij) xPrime_ij = 0.0_dp

      ! Overlap of central atoms
      integral_r(0)%value(0,0) = 0.25_dp/PI !/ this%atom_sigma**1.5_dp

      ! Overlaps of central atom of first environment with the other atoms in the second.
      do n1 = 1, n_neighbour_i
         r1 = neighbour_i(n1)%r
         if(r1 > this%cutoff) cycle
         fcut1 = coordination_function(r1,this%cutoff,this%cutoff_transition_width)
         fac_exp = exp(-0.5_dp*this%atom_sigma*r1**2) * 0.25_dp/PI
         fac_r1r2 = fac_exp * fcut1

         integral_r(0)%value(0,0) = integral_r(0)%value(0,0) + fac_r1r2

         if(do_xPrime_i) then
            dfcut1 = dcoordination_function(r1,this%cutoff,this%cutoff_transition_width)
            grad_fac_exp1 = -fac_exp*this%atom_sigma*r1
            grad_fac_r1r2_1 = fac_exp * dfcut1 + grad_fac_exp1 * fcut1
            grad_ri(n1)%integral_r(0)%value(0,0) = grad_ri(n1)%integral_r(0)%value(0,0) + grad_fac_r1r2_1
         endif
      enddo

      ! Overlaps of central atom of second environment with the other atoms in the first.
      do n2 = 1, n_neighbour_j
         r2 = neighbour_j(n2)%r
         if(r2 > this%cutoff) cycle
         fcut2 = coordination_function(r2,this%cutoff,this%cutoff_transition_width)
         fac_exp = exp(-0.5_dp*this%atom_sigma*r2**2) * 0.25_dp/PI
         fac_r1r2 = fac_exp * fcut2

         integral_r(0)%value(0,0) = integral_r(0)%value(0,0) + fac_r1r2

         if(do_xPrime_j) then
            dfcut2 = dcoordination_function(r2,this%cutoff,this%cutoff_transition_width)
            grad_fac_exp2 = -fac_exp*this%atom_sigma*r2
            grad_fac_r1r2_2 = fac_exp * dfcut2 + grad_fac_exp2 * fcut2
            grad_rj(n2)%integral_r(0)%value(0,0) = grad_rj(n2)%integral_r(0)%value(0,0) + grad_fac_r1r2_2
         endif
      enddo

      ! Overlaps of non-central atoms.
      do n1 = 1, n_neighbour_i
         r1 = neighbour_i(n1)%r

         if(r1 > this%cutoff) cycle
         fcut1 = coordination_function(r1,this%cutoff,this%cutoff_transition_width)
         dfcut1 = dcoordination_function(r1,this%cutoff,this%cutoff_transition_width)
         do n2 = 1, n_neighbour_j
            r2 = neighbour_j(n2)%r

            if(r2 > this%cutoff) cycle
            fcut2 = coordination_function(r2,this%cutoff,this%cutoff_transition_width)
            dfcut2 = dcoordination_function(r2,this%cutoff,this%cutoff_transition_width)

            arg_bess = this%atom_sigma*r1*r2
            fac_exp = exp(-0.5_dp*this%atom_sigma*(r1**2+r2**2))
            fac_r1r2 = fac_exp * fcut1 * fcut2

            if(do_xPrime_i .or. do_xPrime_ij) then
               grad_arg_bess1 = this%atom_sigma*r2
               grad_fac_exp1 = -fac_exp*this%atom_sigma*r1
               grad_fac_r1r2_1 = (fac_exp * dfcut1 + grad_fac_exp1 * fcut1) * fcut2
            endif

            if(do_xPrime_j .or. do_xPrime_ij) then
               grad_arg_bess2 = this%atom_sigma*r1
               grad_fac_exp2 = -fac_exp*this%atom_sigma*r2
               grad_fac_r1r2_2 = (fac_exp * dfcut2 + grad_fac_exp2 * fcut2) * fcut1
            endif

            if(do_xPrime_ij) then
               grad2_fac_exp = fac_exp * this%atom_sigma**2 * r1*r2
               grad2_fac_r1r2 = grad2_fac_exp*fcut1*fcut2 + grad_fac_exp1*fcut1*dfcut2 + grad_fac_exp2*dfcut1*fcut2 + fac_exp*dfcut1*dfcut2
            endif

            do l = 0, this%l_max
               if( l == 0 ) then
                  mo_spher_bess_fi_ki_lm = cosh(arg_bess)/arg_bess
                  mo_spher_bess_fi_ki_l = sinh(arg_bess)/arg_bess
                  if(do_derivative) mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                  if(do_xPrime_ij) mo_spher_bess_fi_ki_lpp = mo_spher_bess_fi_ki_l - (2*l+3)*mo_spher_bess_fi_ki_lp / arg_bess
               else
                  mo_spher_bess_fi_ki_lmm = mo_spher_bess_fi_ki_lm
                  mo_spher_bess_fi_ki_lm = mo_spher_bess_fi_ki_l
                  if(do_derivative) then
                     mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lp
                     mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                     if(do_xPrime_ij) mo_spher_bess_fi_ki_lpp = mo_spher_bess_fi_ki_l - (2*l+3)*mo_spher_bess_fi_ki_lp / arg_bess
                  else
                     mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm - (2*l-1)*mo_spher_bess_fi_ki_lm / arg_bess
                  endif
               endif

               if(do_derivative) grad_mo_spher_bess_fi_ki_l = l * mo_spher_bess_fi_ki_l / arg_bess + mo_spher_bess_fi_ki_lp
               if(do_xPrime_ij) grad2_mo_spher_bess_fi_ki_l = ( l*(2*l+3)*(l-1) + (1.0_dp+2*l)*arg_bess**2 ) * &
               mo_spher_bess_fi_ki_lp / arg_bess**3 + &
               ( 1.0_dp + l*(l-1)/arg_bess**2 ) * mo_spher_bess_fi_ki_lpp

               !radial = mo_spher_bess_fi_ki_l*fac_exp
               radial = mo_spher_bess_fi_ki_l*fac_r1r2

               !if(do_xPrime_i .or. do_xPrime_ij) grad_radial_i = grad_mo_spher_bess_fi_ki_l * grad_arg_bess1 * fac_exp + mo_spher_bess_fi_ki_l * grad_fac_exp1
               if(do_xPrime_i .or. do_xPrime_ij) grad_radial_i = grad_mo_spher_bess_fi_ki_l * grad_arg_bess1 * fac_r1r2 + mo_spher_bess_fi_ki_l * grad_fac_r1r2_1

               !if(do_xPrime_j .or. do_xPrime_ij) grad_radial_j = grad_mo_spher_bess_fi_ki_l * grad_arg_bess2 * fac_exp + mo_spher_bess_fi_ki_l * grad_fac_exp2
               if(do_xPrime_j .or. do_xPrime_ij) grad_radial_j = grad_mo_spher_bess_fi_ki_l * grad_arg_bess2 * fac_r1r2 + mo_spher_bess_fi_ki_l * grad_fac_r1r2_2

               if(do_xPrime_ij) then
                  !grad2_radial_ij = fac_exp * this%atom_sigma**2 * r1 * r2 * mo_spher_bess_fi_ki_l + &
                  !grad_mo_spher_bess_fi_ki_l * grad_arg_bess1 * grad_fac_exp2 + &
                  !grad_mo_spher_bess_fi_ki_l * grad_arg_bess2 * grad_fac_exp1 + &
                  !fac_exp * ( this%atom_sigma * grad_mo_spher_bess_fi_ki_l + grad_arg_bess1*grad_arg_bess2*grad2_mo_spher_bess_fi_ki_l )
                  grad2_radial_ij = grad2_fac_r1r2 * mo_spher_bess_fi_ki_l + &
                  grad_mo_spher_bess_fi_ki_l * grad_arg_bess1 * grad_fac_r1r2_2 + &
                  grad_mo_spher_bess_fi_ki_l * grad_arg_bess2 * grad_fac_r1r2_1 + &
                  fac_r1r2 * ( this%atom_sigma * grad_mo_spher_bess_fi_ki_l + grad_arg_bess1*grad_arg_bess2*grad2_mo_spher_bess_fi_ki_l )

                  grad_spherical_ij(n2,n1)%value(l) = grad_spherical_ij(n2,n1)%value(l) + &
                  radial
               endif
                  
               do m1 = -l, l
                  if(do_xPrime_i .or. do_xPrime_ij) then
                     grad_spherical_i(n1)%spherical_harmonics(l)%value(m1) = &
                     grad_spherical_i(n1)%spherical_harmonics(l)%value(m1) + &
                     radial * neighbour_j(n2)%spherical_harmonics(l)%value(m1)
                  endif

                  if(do_xPrime_j .or. do_xPrime_ij) then
                     grad_spherical_j(n2)%spherical_harmonics(l)%value(m1) = &
                     grad_spherical_j(n2)%spherical_harmonics(l)%value(m1) + &
                     radial * conjg(neighbour_i(n1)%spherical_harmonics(l)%value(m1))
                  endif

                  if(do_xPrime_ij) then
                     grad_spherical_i_radial_j(n2,n1)%spherical_harmonics(l)%value(m1) = &
                     grad_spherical_i_radial_j(n2,n1)%spherical_harmonics(l)%value(m1) + &
                     grad_radial_j * neighbour_j(n2)%spherical_harmonics(l)%value(m1)

                     grad_spherical_j_radial_i(n2,n1)%spherical_harmonics(l)%value(m1) = &
                     grad_spherical_j_radial_i(n2,n1)%spherical_harmonics(l)%value(m1) + &
                     grad_radial_i * conjg(neighbour_i(n1)%spherical_harmonics(l)%value(m1))
                  endif

                  do m2 = -l, l
                     I_lm1m2 =  radial * &
                     conjg(neighbour_i(n1)%spherical_harmonics(l)%value(m1)) * &
                     neighbour_j(n2)%spherical_harmonics(l)%value(m2)

                     integral_r(l)%value(m2,m1) = integral_r(l)%value(m2,m1) + I_lm1m2

                     if(do_xPrime_i .or. do_xPrime_ij) then
                        grad_ri(n1)%integral_r(l)%value(m2,m1) = grad_ri(n1)%integral_r(l)%value(m2,m1) + &
                        grad_radial_i * &
                        conjg(neighbour_i(n1)%spherical_harmonics(l)%value(m1)) * &
                        neighbour_j(n2)%spherical_harmonics(l)%value(m2)
                     endif

                     if(do_xPrime_j .or. do_xPrime_ij) then
                        grad_rj(n2)%integral_r(l)%value(m2,m1) = grad_rj(n2)%integral_r(l)%value(m2,m1) + &
                        grad_radial_j * &
                        conjg(neighbour_i(n1)%spherical_harmonics(l)%value(m1)) * &
                        neighbour_j(n2)%spherical_harmonics(l)%value(m2)
                     endif

                     if(do_xPrime_ij) then
                        grad_rij(n2,n1)%integral_r(l)%value(m2,m1) = grad_rij(n2,n1)%integral_r(l)%value(m2,m1) + &
                        grad2_radial_ij * &
                        conjg(neighbour_i(n1)%spherical_harmonics(l)%value(m1)) * &
                        neighbour_j(n2)%spherical_harmonics(l)%value(m2)
                     endif

                  enddo
               enddo
            enddo
         enddo
      enddo

      gpCovariance_atom_real_space_Calc = 0.0_dp

      do l = 0, this%l_max
         gpCovariance_atom_real_space_Calc = gpCovariance_atom_real_space_Calc + sum(real(integral_r(l)%value)**2) + sum(aimag(integral_r(l)%value)**2)
      enddo

      if(do_xPrime_i) then
         i_data = 0
         do i = 1, n_neighbour_i
               
            i_data = i_data + 1
            do l = 0, this%l_max
               xPrime_i(i_data) = xPrime_i(i_data) + &
               sum(real(integral_r(l)%value)*real(grad_ri(i)%integral_r(l)%value)) + &
               sum(aimag(integral_r(l)%value)*aimag(grad_ri(i)%integral_r(l)%value)) 
            enddo
            i_data = i_data + 1
            xPrime_i(i_data) = 0.0_dp

            do l = 0, this%l_max
               real_mould_size = size(transfer(grad_spherical_i(i)%spherical_harmonics(l)%value(-l:l),real_mould))
               xPrime_i(i_data+1:i_data+real_mould_size) = transfer(matmul(grad_spherical_i(i)%spherical_harmonics(l)%value,conjg(integral_r(l)%value)),real_mould) 
               i_data = i_data + real_mould_size
            enddo
         enddo

          xPrime_i = xPrime_i * 2.0_dp
      endif

      if(do_xPrime_j) then
         i_data = 0
         do i = 1, n_neighbour_j
               
            i_data = i_data + 1
            do l = 0, this%l_max
               xPrime_j(i_data) = xPrime_j(i_data) + &
               sum(real(integral_r(l)%value)*real(grad_rj(i)%integral_r(l)%value)) + &
               sum(aimag(integral_r(l)%value)*aimag(grad_rj(i)%integral_r(l)%value)) 
            enddo
            i_data = i_data + 1
            xPrime_j(i_data) = 0.0_dp

            do l = 0, this%l_max
               real_mould_size = size(transfer(grad_spherical_j(i)%spherical_harmonics(l)%value(-l:l),real_mould))
               xPrime_j(i_data+1:i_data+real_mould_size) = transfer(matmul(integral_r(l)%value,conjg(grad_spherical_j(i)%spherical_harmonics(l)%value)),real_mould) 
               i_data = i_data + real_mould_size
            enddo
         enddo

         xPrime_j = xPrime_j * 2.0_dp
      endif

      if(do_xPrime_ij) then
         i_data = 0
         do i = 1, n_neighbour_i
            i_data = i_data + 1

            ! i-th neighbour, wrt r
            j_data = 0
            do j = 1, n_neighbour_j
               j_data = j_data + 1

               ! d r_i d r_j
               do l = 0, this%l_max
                  xPrime_ij(i_data,j_data) = xPrime_ij(i_data,j_data) + &
                  sum(real(grad_rj(j)%integral_r(l)%value)*real(grad_ri(i)%integral_r(l)%value)) + &
                  sum(aimag(grad_rj(j)%integral_r(l)%value)*aimag(grad_ri(i)%integral_r(l)%value)) + &
                  sum(real(integral_r(l)%value)*real(grad_rij(j,i)%integral_r(l)%value)) + &
                  sum(aimag(integral_r(l)%value)*aimag(grad_rij(j,i)%integral_r(l)%value)) 
               enddo
               j_data = j_data + 1
               xPrime_ij(i_data,j_data) = 0.0_dp

               ! d r_i d Y^{lm}_j
               do l = 0, this%l_max
                  real_mould_size = size(transfer(grad_spherical_j_radial_i(j,i)%spherical_harmonics(l)%value(-l:l),real_mould))
                  xPrime_ij(i_data,j_data+1:j_data+real_mould_size) = transfer( &
                  matmul(integral_r(l)%value, conjg(grad_spherical_j_radial_i(j,i)%spherical_harmonics(l)%value)) + &
                  matmul(grad_ri(i)%integral_r(l)%value, conjg(grad_spherical_j(j)%spherical_harmonics(l)%value)), real_mould)

                  j_data = j_data + real_mould_size
               enddo
            enddo

            i_data = i_data + 1
            xPrime_ij(i_data,:) = 0.0_dp

            do l1 = 0, this%l_max
               do m1 = -l1, l1
                  j_data = 0
                  i_data = i_data + 1

                  ! d Y^{lm}_i d r_j
                  do j = 1, n_neighbour_j
                     j_data = j_data + 1

                     tmp_complex = dot_product(grad_rj(j)%integral_r(l1)%value(-l1:l1,m1), grad_spherical_i(i)%spherical_harmonics(l1)%value(-l1:l1)) + &
                     dot_product(integral_r(l1)%value(-l1:l1,m1),grad_spherical_i_radial_j(j,i)%spherical_harmonics(l1)%value(-l1:l1))
                     xPrime_ij(i_data:i_data+1,j_data) = (/real(tmp_complex),aimag(tmp_complex)/)

                     j_data = j_data + 1
                     xPrime_ij(i_data:i_data+1,j_data) = 0.0_dp

                     do l2 = 0, this%l_max
                        real_mould_size = size(transfer(grad_spherical_j(j)%spherical_harmonics(l2)%value(-l2:l2),real_mould))
                        if(l1 == l2) then
                           xPrime_ij(i_data,j_data+1:j_data+real_mould_size) = transfer( &
                           grad_spherical_i(i)%spherical_harmonics(l1)%value*conjg(grad_spherical_j(j)%spherical_harmonics(l2)%value(m1)) + &
                           grad_spherical_ij(j,i)%value(l1)*integral_r(l1)%value(-l1:l1,m1), real_mould) 

                           xPrime_ij(i_data+1,j_data+1:j_data+real_mould_size) = transfer( &
                           -CPLX_IMAG*grad_spherical_i(i)%spherical_harmonics(l1)%value*conjg(grad_spherical_j(j)%spherical_harmonics(l2)%value(m1)) + &
                           CPLX_IMAG*grad_spherical_ij(j,i)%value(l1)*integral_r(l1)%value(-l1:l1,m1), real_mould) 
                        endif
                        j_data = j_data + real_mould_size
                     enddo

                  enddo
                  i_data = i_data + 1
               enddo
            enddo
         enddo
         xPrime_ij = xPrime_ij * 2.0_dp
      endif


      if(allocated(integral_r)) then
         do l = 0, this%l_max
            if(allocated(integral_r(l)%value)) deallocate(integral_r(l)%value)
         enddo
         deallocate(integral_r)
      endif
      if(allocated(grad_ri)) then
         do i = 1, size(grad_ri)
            if(allocated(grad_ri(i)%integral_r)) then
               do l = 0, this%l_max
                  if(allocated(grad_ri(i)%integral_r(l)%value)) deallocate(grad_ri(i)%integral_r(l)%value)
               enddo
               deallocate(grad_ri(i)%integral_r)
            endif
         enddo
         deallocate(grad_ri)
      endif

      if(allocated(grad_rj)) then
         do i = 1, size(grad_rj)
            if(allocated(grad_rj(i)%integral_r)) then
               do l = 0, this%l_max
                  if(allocated(grad_rj(i)%integral_r(l)%value)) deallocate(grad_rj(i)%integral_r(l)%value)
               enddo
               deallocate(grad_rj(i)%integral_r)
            endif
         enddo
         deallocate(grad_rj)
      endif

      if(do_xPrime_ij) then

         do i = 1, n_neighbour_i
            do j = 1, n_neighbour_j


               do l = 0, this%l_max
                  deallocate(grad_rij(j,i)%integral_r(l)%value)
                  deallocate(grad_spherical_i_radial_j(j,i)%spherical_harmonics(l)%value)
                  deallocate(grad_spherical_j_radial_i(j,i)%spherical_harmonics(l)%value)
               enddo

               deallocate(grad_spherical_ij(j,i)%value)
               deallocate(grad_rij(j,i)%integral_r)
               deallocate(grad_spherical_i_radial_j(j,i)%spherical_harmonics)
               deallocate(grad_spherical_j_radial_i(j,i)%spherical_harmonics)
            enddo
         enddo

         deallocate(grad_rij)
         deallocate(grad_spherical_ij)
         deallocate(grad_spherical_i_radial_j)
         deallocate(grad_spherical_j_radial_i)

      endif

      call finalise(neighbour_i)
      call finalise(neighbour_j)
      call finalise(grad_spherical_i)
      call finalise(grad_spherical_j)

   endfunction gpCovariance_atom_real_space_Calc

!   function gpCovariance_soap_Calc( this, x_i, x_j, xPrime_i, xPrime_j, xPrime_ij, error )
!      type(gpCovariance_soap), intent(in) :: this
!      real(dp), dimension(:), intent(in) :: x_i, x_j
!      real(dp), dimension(:), intent(out), optional, pointer :: xPrime_i, xPrime_j
!      real(dp), dimension(:,:), intent(out), optional, pointer :: xPrime_ij
!      integer, intent(out), optional :: error
!
!      real(dp) :: gpCovariance_soap_Calc
!
!      integer :: l, m, m1, m2, a, i
!      logical :: do_xPrime_i, do_xPrime_j, do_xPrime_ij, do_derivative
!
!      type(cplx_1d_array), dimension(:,:), allocatable :: fourier1_so3, fourier2_so3, dcov_dfourier1, dcov_dfourier2
!      type(cplx_2d_array), dimension(:), allocatable :: int_soap
!
!      INIT_ERROR(error)
!
!      if( .not. this%initialised ) then
!         RAISE_ERROR('gpCovariance_soap_Calc: object not initialised', error)
!      endif
!
!      do_xPrime_i = .false.
!      do_xPrime_j = .false.
!      do_xPrime_ij = .false.
!      if(present(xPrime_i)) do_xPrime_i = associated(xPrime_i)
!      if(present(xPrime_j)) do_xPrime_j = associated(xPrime_j)
!      if(present(xPrime_ij)) do_xPrime_ij = associated(xPrime_ij)
!
!      do_derivative = (do_xPrime_i .or. do_xPrime_j .or. do_xPrime_ij)
!
!      allocate( fourier1_so3(0:this%l_max,this%n_max), fourier2_so3(0:this%l_max,this%n_max), int_soap(0:this%l_max) )
!
!      if(do_xPrime_i) allocate( dcov_dfourier1(0:this%l_max,this%n_max) )
!      if(do_xPrime_j) allocate( dcov_dfourier2(0:this%l_max,this%n_max) )
!
!      do a = 1, this%n_max
!         do l = 0, this%l_max
!            allocate(fourier1_so3(l,a)%value(-l:l))
!            allocate(fourier2_so3(l,a)%value(-l:l))
!            if(do_xPrime_i) allocate(dcov_dfourier1(l,a)%value(-l:l))
!            if(do_xPrime_j) allocate(dcov_dfourier2(l,a)%value(-l:l))
!         enddo
!      enddo
!
!      do l = 0, this%l_max
!         allocate(int_soap(l)%value(-l:l,-l:l))
!         int_soap(l)%value = CPLX_ZERO
!      enddo
!
!
!      i = 0
!      do a = 1, this%n_max
!         do l = 0, this%l_max
!            do m = -l, l
!               fourier1_so3(l,a)%value(m) = cmplx(x_i(i+1), x_i(i+2))
!               fourier2_so3(l,a)%value(m) = cmplx(x_j(i+1), x_j(i+2))
!               i = i + 2
!            enddo
!         enddo
!      enddo
!
!      do a = 1, this%n_max
!         do l = 0, this%l_max
!            do m1 = -l, l
!               do m2 = -l, l
!                  int_soap(l)%value(m2,m1) = int_soap(l)%value(m2,m1) + &
!                     fourier1_so3(l,a)%value(m1) * conjg(fourier2_so3(l,a)%value(m2))
!               enddo
!            enddo
!         enddo
!      enddo
!
!      do a = 1, this%n_max
!         do l = 0, this%l_max
!            if(do_xPrime_i) dcov_dfourier1(l,a)%value = matmul(fourier2_so3(l,a)%value,int_soap(l)%value)
!            if(do_xPrime_j) dcov_dfourier2(l,a)%value = matmul(conjg(int_soap(l)%value),fourier1_so3(l,a)%value)
!         enddo
!      enddo
!
!      gpCovariance_soap_Calc = 0.0_dp
!      do l = 0, this%l_max
!         gpCovariance_soap_Calc = gpCovariance_soap_Calc + sum(real(int_soap(l)%value)**2+aimag(int_soap(l)%value)**2)
!      enddo
!
!      if(do_derivative) then
!         i = 0
!         do a = 1, this%n_max
!            do l = 0, this%l_max
!               do m = -l, l
!                  if(do_xPrime_i) then
!                     xPrime_i(i+1) = real(dcov_dfourier1(l,a)%value(m))
!                     xPrime_i(i+2) = aimag(dcov_dfourier1(l,a)%value(m))
!                  endif
!                  if(do_xPrime_j) then
!                     xPrime_j(i+1) = real(dcov_dfourier2(l,a)%value(m))
!                     xPrime_j(i+2) = aimag(dcov_dfourier2(l,a)%value(m))
!                  endif
!                  i = i + 2
!               enddo
!            enddo
!         enddo
!         if(do_xPrime_i) xPrime_i = xPrime_i*2.0_dp
!         if(do_xPrime_j) xPrime_j = xPrime_j*2.0_dp
!      endif
!
!      if(allocated(fourier1_so3)) then
!         do a = 1, this%n_max
!            do l = 0, this%l_max
!               deallocate(fourier1_so3(l,a)%value)
!            enddo
!         enddo
!         deallocate(fourier1_so3)
!      endif
!
!      if(allocated(fourier2_so3)) then
!         do a = 1, this%n_max
!            do l = 0, this%l_max
!               deallocate(fourier2_so3(l,a)%value)
!            enddo
!         enddo
!         deallocate(fourier2_so3)
!      endif
!
!      if(allocated(int_soap)) then
!         do l = 0, this%l_max
!            deallocate(int_soap(l)%value)
!         enddo
!         deallocate(int_soap)
!      endif
!
!      if(allocated(dcov_dfourier1)) then
!         do a = 1, this%n_max
!            do l = 0, this%l_max
!               deallocate(dcov_dfourier1(l,a)%value)
!            enddo
!         enddo
!         deallocate(dcov_dfourier1)
!      endif
!
!      if(allocated(dcov_dfourier2)) then
!         do a = 1, this%n_max
!            do l = 0, this%l_max
!               deallocate(dcov_dfourier2(l,a)%value)
!            enddo
!         enddo
!         deallocate(dcov_dfourier2)
!      endif
!
!   endfunction gpCovariance_soap_Calc

   subroutine gpRealArray_NeighbourDescriptor(this,x,x_size,neighbour,n_neighbour)
      type(gpCovariance_atom_real_space), intent(in) :: this
      real(dp), dimension(:), intent(in) :: x
      integer, intent(in) :: x_size
      type(neighbour_descriptor), dimension(:), allocatable, intent(out) :: neighbour
      integer, intent(out) :: n_neighbour

      integer :: l, i_data, i, real_mould_size
      real(dp), dimension(1) :: real_mould
      complex(dp), dimension(1) :: complex_mould

      n_neighbour = x_size / ( 2 * (this%l_max+1)**2 + 2 )
      
      call finalise(neighbour)

      allocate(neighbour(n_neighbour))

      i_data = 0
      do i = 1, n_neighbour
            
         i_data = i_data + 1
         neighbour(i)%r = x(i_data)
         i_data = i_data + 1
         neighbour(i)%n = abs(nint(x(i_data)))

         allocate(neighbour(i)%spherical_harmonics(0:this%l_max))
         do l = 0, this%l_max

            allocate(neighbour(i)%spherical_harmonics(l)%value(-l:l))

            real_mould_size = size(transfer(neighbour(i)%spherical_harmonics(l)%value(-l:l),real_mould))
            neighbour(i)%spherical_harmonics(l)%value = transfer(x(i_data+1:i_data+real_mould_size),complex_mould)
            i_data = i_data + real_mould_size
         enddo
      enddo

   endsubroutine gpRealArray_NeighbourDescriptor

   subroutine gpNeighbourDescriptor_Finalise(this)
      type(neighbour_descriptor), dimension(:), allocatable, intent(inout) :: this

      integer :: i, l

      if(allocated(this)) then
         do i = 1, size(this)
            do l = lbound(this(i)%spherical_harmonics,dim=1), ubound(this(i)%spherical_harmonics,dim=1)
               if(allocated(this(i)%spherical_harmonics(l)%value)) deallocate(this(i)%spherical_harmonics(l)%value)
            enddo
            if(allocated(this(i)%spherical_harmonics)) deallocate(this(i)%spherical_harmonics)
         enddo
         deallocate(this)
      endif

   endsubroutine gpNeighbourDescriptor_Finalise

   subroutine gp_atom_real_space_RealArray_XYZ(this,x_array,x_array_size,xyz_array,xyz_array_size)
      type(gpCovariance_atom_real_space), intent(in) :: this
      real(dp), dimension(:), intent(in), target :: x_array
      integer, intent(in) :: x_array_size
      real(dp), dimension(:), allocatable, intent(out) :: xyz_array
      integer, intent(out) :: xyz_array_size

      integer :: l, i_data, i, n_neighbour, n, n_data, xyz_start, xyz_end
      real(dp), dimension(:), pointer :: Y1_array
      real(dp), pointer :: Re_Y1m1, Im_Y1m1, Re_Y10
      real(dp) :: r, x, y, z

      real(dp), parameter :: xy_factor = 0.5_dp * sqrt(1.5_dp / PI)
      real(dp), parameter :: z_factor = 0.5_dp * sqrt(3.0_dp / PI)

      n_neighbour = x_array_size / ( 2 * (this%l_max+1)**2 + 2 )
      xyz_array_size = n_neighbour*3
      
      if(allocated(xyz_array)) deallocate(xyz_array)

      allocate(xyz_array(xyz_array_size))

      i_data = 0
      do i = 1, n_neighbour
            
         i_data = i_data + 1
         r = x_array(i_data)
         i_data = i_data + 1
         n = abs(nint(x_array(i_data)))

         do l = 0, this%l_max
            n_data = 2*(2*l + 1)

            if(l == 1) then
               Y1_array => x_array(i_data+1:i_data+n_data)
               Re_Y1m1 => Y1_array(1)
               Im_Y1m1 => Y1_array(2)
               Re_Y10 =>  Y1_array(3)
               !Im_Y10 =>  Y1_value(4)
               !Re_Y1p1 => Y1_value(5)
               !Im_Y1p1 => Y1_value(6)

               z = Re_Y10 * r / z_factor
               x = Re_Y1m1 * r / xy_factor
               y = -Im_Y1m1 * r / xy_factor
            endif

            i_data = i_data + n_data
         enddo

         xyz_start = (i-1)*3+1
         xyz_end = 3*i
         xyz_array(xyz_start:xyz_end) = (/x,y,z/)
      enddo

   endsubroutine gp_atom_real_space_RealArray_XYZ

   subroutine gp_atom_real_space_XYZ_RealArray(this,xyz_array,xyz_array_size,x_array,x_array_size)
      type(gpCovariance_atom_real_space), intent(in) :: this
      real(dp), dimension(:), intent(in), target :: xyz_array
      integer, intent(in) :: xyz_array_size
      real(dp), dimension(:), allocatable, intent(out) :: x_array
      integer, intent(out) :: x_array_size

      integer :: l, m, i_data, i, n_neighbour, xyz_start, xyz_end
      real(dp), dimension(:), pointer :: xyz
      complex(dp) :: Y_lm

      n_neighbour = xyz_array_size / 3
      x_array_size = n_neighbour * ( 2 * (this%l_max+1)**2 + 2 )
      
      if(allocated(x_array)) deallocate(x_array)

      allocate(x_array(x_array_size))

      i_data = 0
      do i = 1, n_neighbour
            
         xyz_start = (i-1)*3+1
         xyz_end = 3*i
         xyz => xyz_array(xyz_start:xyz_end)

         i_data = i_data + 1
         x_array(i_data) = norm(xyz)
         i_data = i_data + 1
         x_array(i_data) = real(i,dp)

         do l = 0, this%l_max
            do m = -l, l
               Y_lm = SphericalYCartesian(l,m,xyz)
               x_array(i_data+1:i_data+2) = (/real(Y_lm),aimag(Y_lm)/)
               i_data = i_data + 2
            enddo
         enddo

      enddo

   endsubroutine gp_atom_real_space_XYZ_RealArray

   function gpCoordinates_Predict( this, xStar, gradPredict, sparseScore, do_sparseScore, error )
      type(gpCoordinates), intent(inout), target :: this
      real(dp), dimension(:), intent(in) :: xStar
      real(dp), dimension(:), intent(out), optional :: gradPredict
      real(dp), intent(out), optional :: sparseScore
      logical, intent(in) :: do_sparseScore
      integer, optional, intent(out) :: error

      real(dp) :: gpCoordinates_Predict

      real(dp) :: covarianceExp, zeta
      real(dp), pointer :: fc_i
      real(dp), dimension(:), pointer :: x_i
      real(dp), dimension(this%d) :: xI_xJ_theta
      real(dp), dimension(this%n_sparseX) :: k
      integer :: i_sparseX, i_p
      real(dp) :: delta, covariance_x_x, covariance_x_xStar
      real(dp), allocatable :: covariance_x_xStars(:), alpha_scaled(:)
      real(dp), dimension(:), pointer :: xPrime_i
      real(dp), dimension(:), allocatable, target :: grad_kStar, k_mm_k
      real(dp), dimension(:,:), allocatable, target :: grad_k
      logical :: my_do_sparseScore


      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_Predict: object not initialised', error)
      endif

      if(this%bond_real_space) then
         if((.not. this%bond_real_space_cov%initialised) .or. (.not. this%bond_real_space_cov%sparseX_inisialised)) then
            call gpCoordinates_gpCovariance_bond_real_space_Initialise(this)
         endif
      endif

      if(present(sparseScore)) then
         my_do_sparseScore = optional_default(.true.,do_sparseScore)
      else
         my_do_sparseScore = .false.
      endif

      if(my_do_sparseScore) then
         allocate(k_mm_k(this%n_sparseX))
         if(.not.this%sparseScore_initialised) then
            RAISE_ERROR('gpCoordinates_Predict: sparseScore not initialised',error)
         endif
      endif

      if(this%atom_real_space) then
         if((.not. this%atom_real_space_cov%initialised)) then
            call gpCoordinates_gpCovariance_atom_real_space_Initialise(this)
         endif
      endif

      k = 0.0_dp
      xPrime_i => null()
      if(present(gradPredict)) then
         allocate(grad_k(size(xStar),this%n_sparseX))
      endif

      if(this%atom_real_space) then
         if(present(gradPredict)) then
            allocate(grad_kStar(size(xStar)))
            xPrime_i => grad_kStar
         endif
         if(this%atom_real_space) covariance_x_x = gpCovariance_atom_real_space_Calc(this%atom_real_space_cov, x_i=xStar, x_i_size=size(xStar), x_j=xStar, x_j_size=size(xStar), xPrime_i = xPrime_i)
      endif

      if (this%soap) then
	 allocate(covariance_x_xStars(this%n_sparseX))
	 call dgemv('T', size(this%sparseX,1), size(this%sparseX,2), 1.0_dp, this%sparseX(1,1), size(this%sparseX, 1), &
	            xStar(1), 1, 0.0_dp, covariance_x_xStars(1), 1)
	 k(:) = this%delta**2 * covariance_x_xStars(:)**this%theta(1)
	 if(present(gradPredict)) then
	    allocate(alpha_scaled(size(this%alpha)))
	    alpha_scaled(:) = this%alpha(:) * this%delta**2 * this%theta(1) * covariance_x_xStars(:)**(this%theta(1)-1.0_dp)
	 endif
	 deallocate(covariance_x_xStars)
      else
	 xPrime_i => null()
	 do i_sparseX = 1, this%n_sparseX
	    if(this%bond_real_space) then
	       delta = this%bond_real_space_cov%delta
	       this%bond_real_space_cov%delta = 1.0_dp
	       covariance_x_x = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, x_i=xStar, x_i_size=(size(xStar) - 1), x_j=xStar, x_j_size=(size(xStar) - 1))
	       this%bond_real_space_cov%delta = delta
	       k(i_sparseX) = gpCovariance_bond_real_space_Calc(this%bond_real_space_cov, x_i=xStar, x_i_size=(size(xStar) - 1), j_sparseX=i_sparseX) &
			      / sqrt(covariance_x_x * this%covarianceDiag_sparseX_sparseX(i_sparseX))
	    elseif(this%atom_real_space) then

	       if(present(gradPredict)) xPrime_i => grad_k(:,i_sparseX)

	       zeta = this%atom_real_space_cov%zeta
	       delta = this%atom_real_space_cov%delta

	       covariance_x_xStar = gpCovariance_atom_real_space_Calc(this%atom_real_space_cov, &
	       x_i=xStar, x_i_size=size(xStar), &
	       x_j=this%sparseX(:,i_sparseX), x_j_size=this%sparseX_size(i_sparseX), &
	       xPrime_i = xPrime_i)

	       if(present(gradPredict)) then
		  grad_k(:,i_sparseX) = covariance_x_xStar**(zeta-1.0_dp) * grad_k(:,i_sparseX) / sqrt(covariance_x_x * this%covarianceDiag_sparseX_sparseX(i_sparseX))**zeta - &
		  grad_kStar * ( covariance_x_xStar / sqrt(covariance_x_x * this%covarianceDiag_sparseX_sparseX(i_sparseX)) )**zeta / covariance_x_x
		  grad_k(:,i_sparseX) = grad_k(:,i_sparseX) * zeta * delta**2
	       endif

	       k(i_sparseX) = ( covariance_x_xStar / sqrt(covariance_x_x * this%covarianceDiag_sparseX_sparseX(i_sparseX)) )**zeta * delta**2

! now a single dgemv call outside the loop 
!	    elseif(this%soap) then
!	       covariance_x_xStar = dot_product(xStar,this%sparseX(:,i_sparseX))
!
!	       k(i_sparseX) = this%delta**2 * covariance_x_xStar**this%theta(1)
!
!	       if(present(gradPredict)) grad_k(:,i_sparseX) = this%delta**2 * this%theta(1) * covariance_x_xStar**(this%theta(1)-1.0_dp) * this%sparseX(:,i_sparseX)
	    else
	       x_i => this%sparseX(:,i_sparseX)
	       fc_i => this%sparseCutoff(i_sparseX)
	       do i_p = 1, this%n_permutations

		  xI_xJ_theta = (x_i(this%permutations(:,i_p)) - xStar(1:)) / this%theta

		  covarianceExp = this%delta**2 * exp( -0.5_dp * dot_product(xI_xJ_theta,xI_xJ_theta) )

		  if(present(gradPredict)) &
		  grad_k(:,i_sparseX) = grad_k(:,i_sparseX) + covarianceExp*xI_xJ_theta / this%theta
		  k(i_sparseX) = k(i_sparseX) + covarianceExp
	       enddo
	       k(i_sparseX) = k(i_sparseX) * fc_i
	       if(present(gradPredict)) grad_k(:,i_sparseX) = grad_k(:,i_sparseX) * fc_i
	    endif
	 enddo
      endif
      k = k + this%f0**2
      gpCoordinates_Predict = dot_product( k, this%alpha )

      if(my_do_sparseScore) then
         call Matrix_Solve(this%LA_k_mm,k,k_mm_k)
         sparseScore = (this%delta**2 + this%f0**2) - dot_product(k,k_mm_k)
         if(allocated(k_mm_k)) deallocate(k_mm_k)
      endif

      if (this%soap) then
	 if(present(gradPredict)) &
	    call dgemv('N', size(this%sparseX,1), size(this%sparseX,2), 1.0_dp, this%sparseX(1,1), size(this%sparseX,1), &
	       alpha_scaled(1), 1, 0.0_dp, gradPredict(1), 1)
      else
	 if(present(gradPredict)) &
	    call dgemv('N', size(grad_k,1), size(grad_k,2), 1.0_dp, grad_k(1,1), size(grad_k,1), &
	       this%alpha(1), 1, 0.0_dp, gradPredict(1), 1)
      endif

      if(allocated(alpha_scaled)) deallocate(alpha_scaled)
      if(allocated(grad_k)) deallocate(grad_k)
      if(allocated(grad_kStar)) deallocate(grad_kStar)

   endfunction gpCoordinates_Predict

   subroutine gpCoordinates_initialise_SparseScore(this,error)
      type(gpCoordinates), intent(inout), target :: this
      integer, intent(out), optional :: error

      real(dp), dimension(:,:), allocatable :: k_mm
      real(dp), dimension(:), pointer :: x_i, x_j
      real(dp), pointer :: fc_i, fc_j
      real(dp), dimension(this%d) :: xI_xJ_theta

      integer :: i, j, i_p

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_initialise_SparseScore: object not initialised', error)
      endif

      if(this%sparseScore_initialised) call gpCoordinates_finalise_SparseScore(this,error)

      allocate(k_mm(this%n_sparseX,this%n_sparseX))
      k_mm = 0.0_dp

      do i = 1, this%n_sparseX
         x_i => this%sparseX(:,i)
         fc_i => this%sparseCutoff(i)
         do j = i, this%n_sparseX
            if(this%bond_real_space) then
            elseif(this%atom_real_space) then
               k_mm(j,i) = ( gpCovariance_atom_real_space_Calc( this%atom_real_space_cov, x_i=this%sparseX(:,i), x_i_size=this%sparseX_size(i), &
               x_j=this%sparseX(:,j), x_j_size=this%sparseX_size(j) ) / &
               sqrt(this%covarianceDiag_sparseX_sparseX(i)*this%covarianceDiag_sparseX_sparseX(j)) )**this%atom_real_space_cov%zeta * this%atom_real_space_cov%delta**2
            elseif(this%soap) then
               k_mm(j,i) = this%delta**2 * dot_product( this%sparseX(:,i),this%sparseX(:,j) )**this%theta(1)
            else
               x_j => this%sparseX(:,j)
               fc_j => this%sparseCutoff(j)
               do i_p = 1, this%n_permutations
                  xI_xJ_theta = (x_i(this%permutations(:,i_p)) - x_j) / this%theta
                  k_mm(j,i) = k_mm(j,i) + this%delta**2 * exp( -0.5_dp * dot_product(xI_xJ_theta,xI_xJ_theta) ) * fc_i * fc_j
               enddo
            endif
            k_mm(i,j) = k_mm(j,i)
         enddo
      enddo

      k_mm = k_mm + this%f0**2

      call initialise(this%LA_k_mm, k_mm)
      if(allocated(k_mm)) deallocate(k_mm)

      this%sparseScore_initialised = .true.

   endsubroutine gpCoordinates_initialise_SparseScore

   subroutine gpCoordinates_finalise_SparseScore(this,error)
      type(gpCoordinates), intent(inout) :: this
      integer, intent(out), optional :: error
                     
      INIT_ERROR(error)

      if( .not. this%sparseScore_initialised) return

      call finalise(this%LA_k_mm)

      this%sparseScore_initialised = .false.

   endsubroutine gpCoordinates_finalise_SparseScore

   subroutine gpCoordinates_print_sparseX_file(this,label,sparseX_base_filename,error)
      type(gpCoordinates), intent(in) :: this
      character(len=*), intent(in), optional :: label
      character(len=*), intent(in), optional :: sparseX_base_filename
      integer, intent(out), optional :: error

      logical :: have_sparseX_base_filename

      INIT_ERROR(error)

      have_sparseX_base_filename = .false.
      if (present(sparseX_base_filename)) then
	 if (len_trim(sparseX_base_filename) > 0) have_sparseX_base_filename = .true.
      endif

      if (have_sparseX_base_filename) then
	 if (present(label)) then
	    call fwrite_array_d(size(this%sparseX), this%sparseX(1,1), trim(sparseX_base_filename)//"."//trim(label)//C_NULL_CHAR)
	 else
	    call fwrite_array_d(size(this%sparseX), this%sparseX(1,1), trim(sparseX_base_filename)//C_NULL_CHAR)
	 endif
      endif

   end subroutine gpCoordinates_print_sparseX_file

   subroutine gpCoordinates_printXML(this,xf,label,sparseX_base_filename,error)
      type(gpCoordinates), intent(in) :: this
      type(xmlf_t), intent(inout) :: xf
      character(len=*), intent(in), optional :: label
      character(len=*), intent(in), optional :: sparseX_base_filename
      integer, intent(out), optional :: error

      integer :: i, j, j_end
      real(dp), dimension(:), allocatable :: xyz_array
      integer :: xyz_array_size
      logical :: have_sparseX_base_filename

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_printXML: object not initialised', error)
      endif

      have_sparseX_base_filename = .false.
      if (present(sparseX_base_filename)) then
	 if (len_trim(sparseX_base_filename) > 0) have_sparseX_base_filename = .true.
      endif

      call xml_NewElement(xf,"gpCoordinates")

      if(present(label)) call xml_AddAttribute(xf,"label", trim(label))

      call xml_AddAttribute(xf,"dimensions", ""//this%d)
      call xml_AddAttribute(xf,"signal_variance", ""//this%delta)
      call xml_AddAttribute(xf,"signal_mean", ""//this%f0)
      call xml_AddAttribute(xf,"sparsified", ""//this%sparsified)
      call xml_AddAttribute(xf,"n_permutations", ""//this%n_permutations)
      if(this%bond_real_space) call xml_AddAttribute(xf,"bond_real_space", ""//this%bond_real_space)
      if(this%atom_real_space) call xml_AddAttribute(xf,"atom_real_space", ""//this%atom_real_space)
      if(this%soap) call xml_AddAttribute(xf,"soap", ""//this%soap)

      if(this%sparsified) then
         call xml_AddAttribute(xf,"n_sparseX",""//this%n_sparseX)
         if(this%bond_real_space .or. this%atom_real_space) call xml_AddAttribute(xf,"sparseX_size_max", ""//maxval(this%sparseX_size))
	 if (have_sparseX_base_filename) then
	    if (present(label)) then
	       call xml_AddAttribute(xf,"sparseX_filename",trim(sparseX_base_filename)//"."//trim(label))
	    else
	       call xml_AddAttribute(xf,"sparseX_filename",trim(sparseX_base_filename))
	    endif
	 endif
      else
         call xml_AddAttribute(xf,"n_x",""//this%n_x)
         call xml_AddAttribute(xf,"n_xPrime",""//this%n_xPrime)
         if(this%bond_real_space .or. this%atom_real_space) call xml_AddAttribute(xf,"x_size_max", ""//maxval(this%x_size))
         if(this%bond_real_space .or. this%atom_real_space) call xml_AddAttribute(xf,"xPrime_size_max", ""//maxval(this%xPrime_size))
      endif

      call xml_NewElement(xf,"theta")
      call xml_AddCharacters(xf, ""//this%theta//" ")
      call xml_EndElement(xf,"theta")

      call xml_NewElement(xf,"descriptor")
      call xml_AddCharacters(xf, string(this%descriptor_str))
      call xml_EndElement(xf,"descriptor")

      do i = 1, this%n_permutations
         call xml_NewElement(xf,"permutation")
         call xml_AddAttribute(xf,"i",""//i)
         call xml_AddCharacters(xf,""//this%permutations(:,i)//" ")
         call xml_EndElement(xf,"permutation")
      enddo

      if(this%sparsified) then
         do i = 1, this%n_sparseX
            call xml_NewElement(xf,"sparseX")
            call xml_AddAttribute(xf,"i", ""//i)
            call xml_AddAttribute(xf,"alpha", ""//this%alpha(i))
            call xml_AddAttribute(xf,"sparseCutoff", ""//this%sparseCutoff(i))
            if(this%bond_real_space .or. this%atom_real_space) then
               call xml_AddAttribute(xf,"covariance_sparseX_sparseX", ""//this%covarianceDiag_sparseX_sparseX(i))
            endif
            if(this%atom_real_space) then
               call gp_atom_real_space_RealArray_XYZ(this%atom_real_space_cov,this%sparseX(:,i),this%sparseX_size(i), &
               xyz_array,xyz_array_size)               
               call xml_AddAttribute(xf,"sparseX_size", ""//xyz_array_size)
               call xml_AddCharacters(xf, ""//xyz_array//"  ")
               if(allocated(xyz_array)) deallocate(xyz_array)
            elseif(this%bond_real_space) then
               call xml_AddAttribute(xf,"sparseX_size", ""//this%sparseX_size(i))
               call xml_AddCharacters(xf, ""//this%sparseX(:this%sparseX_size(i),i)//"  ")
            elseif (.not. have_sparseX_base_filename) then
	       if(this%d <= 50) then
		  call xml_AddCharacters(xf, ""//this%sparseX(:,i)//"  ")
	       else
		  call xml_AddAttribute(xf,"sliced", "T")
		  do j = 1, this%d, 50
		     j_end = min(j-1+50,this%d)
		     call xml_NewElement(xf,"sparseX_slice")
		     call xml_AddAttribute(xf,"start", ""//j)
		     call xml_AddAttribute(xf,"end", ""//j_end)
		     call xml_AddCharacters(xf, ""//this%sparseX(j:j_end,i)//"  ")
		     call xml_EndElement(xf,"sparseX_slice")
		  enddo
	       endif
            endif
            call xml_EndElement(xf,"sparseX")
         enddo
      else
         do i = 1, this%n_x
            call xml_NewElement(xf,"x")
            call xml_AddAttribute(xf,"i", ""//i)
            call xml_AddAttribute(xf,"map_x_y", ""//this%map_x_y(i))
            call xml_AddAttribute(xf,"cutoff", ""//this%cutoff(i))
            if(this%bond_real_space .or. this%atom_real_space) call xml_AddAttribute(xf,"x_size", ""//this%x_size(i))
            if(this%bond_real_space .or. this%atom_real_space) call xml_AddAttribute(xf,"covariance_x_x", ""//this%covarianceDiag_x_x(i))
            call xml_AddCharacters(xf, ""//this%x(:,i)//" ")
            call xml_EndElement(xf,"x")
         enddo
         do i = 1, this%n_xPrime
            call xml_NewElement(xf,"xPrime")
            call xml_AddAttribute(xf,"i", ""//i)
            call xml_AddAttribute(xf,"map_xPrime_yPrime", ""//this%map_xPrime_yPrime(i))
            call xml_AddAttribute(xf,"map_xPrime_x", ""//this%map_xPrime_x(i))
            call xml_AddAttribute(xf,"cutoffPrime", ""//this%cutoffPrime(i))
            if(this%bond_real_space .or. this%atom_real_space) call xml_AddAttribute(xf,"xPrime_size", ""//this%xPrime_size(i))
            if(this%bond_real_space .or. this%atom_real_space) call xml_AddAttribute(xf,"covariance_xPrime_xPrime", ""//this%covarianceDiag_xPrime_xPrime(i))
            call xml_AddCharacters(xf, ""//this%xPrime(:,i)//" ")
            call xml_EndElement(xf,"xPrime")
         enddo
      endif

      call xml_EndElement(xf,"gpCoordinates")

   endsubroutine gpCoordinates_printXML

   subroutine gpFull_printXML(this,xf,label,error)
      type(gpFull), intent(in) :: this
      type(xmlf_t), intent(inout) :: xf
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

      integer :: i

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_printXML: object not initialised', error)
      endif

      call xml_NewElement(xf,"gpFull")

      if(present(label)) call xml_AddAttribute(xf,"label", trim(label))

      call xml_AddAttribute(xf,"n_y", ""//this%n_y)
      call xml_AddAttribute(xf,"n_yPrime", ""//this%n_yPrime)
      call xml_AddAttribute(xf,"n_globalSparseX", ""//this%n_globalSparseX)
      call xml_AddAttribute(xf,"n_coordinate", ""//this%n_coordinate)
      call xml_AddAttribute(xf,"sparse_jitter", ""//this%sparse_jitter)

      do i = 1, this%n_y
         call xml_NewElement(xf,"y")
         call xml_AddAttribute(xf,"i", ""//i)
         call xml_AddAttribute(xf,"map_y_globalY", ""//this%map_y_globalY(i))
         call xml_AddAttribute(xf,"alpha", ""//this%alpha(this%map_y_globalY(i)) )
         call xml_EndElement(xf,"y")
      enddo

      do i = 1, this%n_yPrime
         call xml_NewElement(xf,"yPrime")
         call xml_AddAttribute(xf,"i", ""//i)
         call xml_AddAttribute(xf,"map_yPrime_globalY", ""//this%map_yPrime_globalY(i))
         call xml_AddAttribute(xf,"alpha", ""//this%alpha(this%map_yPrime_globalY(i)) )
         call xml_EndElement(xf,"yPrime")
      enddo

      do i = 1, this%n_coordinate
         call gpCoordinates_printXML(this%coordinate(i),xf,label=trim(optional_default("",label))//i,error=error)
      enddo

      call xml_EndElement(xf,"gpFull")

   endsubroutine gpFull_printXML

   subroutine gpSparse_printXML(this,xf,label,sparseX_base_filename,error)
      type(gpSparse), intent(in) :: this
      type(xmlf_t), intent(inout) :: xf
      character(len=*), intent(in), optional :: label
      character(len=*), intent(in), optional :: sparseX_base_filename
      integer, intent(out), optional :: error

      logical :: have_sparseX_base_filename
      integer :: i

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpSparse_printXML: object not initialised', error)
      endif

      have_sparseX_base_filename = .false.
      if (present(sparseX_base_filename)) then
	 if (len_trim(sparseX_base_filename) > 0) have_sparseX_base_filename = .true.
      endif

      call xml_NewElement(xf,"gpSparse")

      if(present(label)) call xml_AddAttribute(xf,"label", trim(label))

      call xml_AddAttribute(xf,"n_coordinate", ""//this%n_coordinate)

      do i = 1, this%n_coordinate
         call gpCoordinates_printXML(this%coordinate(i),xf,label=trim(optional_default("",label))//i,&
				     sparseX_base_filename=sparseX_base_filename, error=error)
	 if (have_sparseX_base_filename) &
	    call gpCoordinates_print_sparseX_file(this%coordinate(i),label=trim(optional_default("",label))//i,&
				                  sparseX_base_filename=sparseX_base_filename, error=error)
      enddo

      call xml_EndElement(xf,"gpSparse")

   endsubroutine gpSparse_printXML

   subroutine gpCoordinates_readXML(this,xp,label,error)
      type(gpCoordinates), intent(inout), target :: this
      type(xml_t), intent(inout) :: xp
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

      INIT_ERROR(error)

      if( this%initialised ) call finalise(this,error)

      parse_in_gpCoordinates = .false.
      parse_matched_label = .false.
      parse_gpCoordinates => this
      parse_gpCoordinates_label = optional_default("",label)

      call initialise(parse_cur_data)
      call parse(xp, &
         characters_handler = gpCoordinates_characters_handler, &
         startElement_handler = gpCoordinates_startElement_handler, &
         endElement_handler = gpCoordinates_endElement_handler)

      call finalise(parse_cur_data)

      this%initialised = .true.

   endsubroutine gpCoordinates_readXML

   subroutine gpFull_readXML(this,xp,label,error)
      type(gpFull), intent(inout), target :: this
      type(xml_t), intent(inout) :: xp
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

      integer :: i

      INIT_ERROR(error)

      if( this%initialised ) call finalise(this,error)

      parse_in_gpFull = .false.
      parse_matched_label = .false.
      parse_gpFull => this
      parse_gpFull_label = optional_default("",label)

      call initialise(parse_cur_data)

      call parse(xp, &
         characters_handler = gpFull_characters_handler, &
         startElement_handler = gpFull_startElement_handler, &
         endElement_handler = gpFull_endElement_handler)

      call finalise(parse_cur_data)

      do i = 1, this%n_coordinate
         call gpCoordinates_readXML(this%coordinate(i),xp,label=trim(parse_gpFull_label)//i,error=error)
      enddo

      this%initialised = .true.

   endsubroutine gpFull_readXML

   subroutine gpSparse_readXML(this,xp,label,error)
      type(gpSparse), intent(inout), target :: this
      type(xml_t), intent(inout) :: xp
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

!      integer :: i

      INIT_ERROR(error)

      if( this%initialised ) call finalise(this,error)

      parse_in_gpSparse = .false.
      parse_gpSparse => this
      parse_matched_label = .false.
      parse_gpSparse_label = optional_default("",label)

      call initialise(parse_cur_data)

      call parse(xp, &
         characters_handler = gpSparse_characters_handler, &
         startElement_handler = gpSparse_startElement_handler, &
         endElement_handler = gpSparse_endElement_handler)

      call finalise(parse_cur_data)

!      do i = 1, this%n_coordinate
!         call gpCoordinates_readXML(this%coordinate(i),xp,label=trim(parse_gpSparse_label)//i,error=error)
!      enddo

      this%initialised = .true.

   endsubroutine gpSparse_readXML

   subroutine gpFull_readXML_string(this,params_str,label,error)
      type(gpFull), intent(inout), target :: this
      character(len=*), intent(in) :: params_str
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

      type(xml_t) :: xp

      INIT_ERROR(error)

      call open_xml_string(xp, params_str)
      call gp_readXML(this,xp,label,error)
      call close_xml_t(xp)

   endsubroutine gpFull_readXML_string

   subroutine gpCoordinates_readXML_string(this,params_str,label,error)
      type(gpCoordinates), intent(inout), target :: this
      character(len=*), intent(in) :: params_str
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

      type(xml_t) :: xp

      INIT_ERROR(error)

      call open_xml_string(xp, params_str)
      call gp_readXML(this,xp,label,error)
      call close_xml_t(xp)

   endsubroutine gpCoordinates_readXML_string

   subroutine gpSparse_readXML_string(this,params_str,label,error)
      type(gpSparse), intent(inout), target :: this
      character(len=*), intent(in) :: params_str
      character(len=*), intent(in), optional :: label
      integer, intent(out), optional :: error

      type(xml_t) :: xp
      integer :: i

      INIT_ERROR(error)

      call open_xml_string(xp, params_str)
      call gp_readXML(this,xp,label,error)
      call close_xml_t(xp)

      do i = 1, this%n_coordinate
         call gp_readXML(this%coordinate(i),params_str,label=trim(parse_gpSparse_label)//i,error=error)
      enddo

   endsubroutine gpSparse_readXML_string

   subroutine gpCoordinates_startElement_handler(URI, localname, name, attributes)
      character(len=*), intent(in)   :: URI
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name
      type(dictionary_t), intent(in) :: attributes

      real(dp) :: delta, f0
      integer :: status, d, n_sparseX, n_x, n_xPrime, n_permutations, i, x_size_max, xPrime_size_max, sparseX_size_max
      logical :: sparsified, bond_real_space, atom_real_space, soap
      character(len=1024) :: value

      if(name == 'gpCoordinates') then ! new GP_data
         if(parse_in_gpCoordinates) then
            call system_abort("gpCoordinates_startElement_handler entered gpCoordinates with parse_in_gpCoordinates true. Probably a bug in FoX (4.0.1, e.g.)")
         endif

         if(parse_matched_label) return ! we already found an exact match for this label

         call GP_FoX_get_value(attributes, 'label', value, status)
         if (status /= 0) value = ''

         if(len(trim(parse_gpCoordinates_label)) > 0) then ! we were passed in a label
            if(trim(value) == trim(parse_gpCoordinates_label)) then
               parse_matched_label = .true.
               parse_in_gpCoordinates = .true.
            else ! no match
               parse_in_gpCoordinates = .false.
            endif
         else ! no label passed in
            parse_in_gpCoordinates = .true.
         endif

         if(parse_in_gpCoordinates) then
            if(parse_gpCoordinates%initialised) call finalise(parse_gpCoordinates)

            call GP_FoX_get_value(attributes, 'dimensions', value, status)
            if (status == 0) then
               read (value,*) d
            else
               call system_abort("gpCoordinates_startElement_handler did not find the dimensions attribute.")
            endif

            call GP_FoX_get_value(attributes, 'signal_variance', value, status)
            if (status == 0) then
               read (value,*) delta
            else
               call system_abort("gpCoordinates_startElement_handler did not find the signal_variance attribute.")
            endif

            call GP_FoX_get_value(attributes, 'signal_mean', value, status)
            if (status == 0) then
               read (value,*) f0
            else
               call system_abort("gpCoordinates_startElement_handler did not find the signal_variance attribute.")
            endif


            call GP_FoX_get_value(attributes, 'sparsified', value, status)
            if (status == 0) then
               read (value,*) sparsified
            else
               call system_abort("gpCoordinates_startElement_handler did not find the sparsified attribute.")
            endif

            call GP_FoX_get_value(attributes, 'n_permutations', value, status)
            if (status == 0) then
               read (value,*) n_permutations
            else
               call system_abort("gpCoordinates_startElement_handler did not find the n_permutations attribute.")
            endif

            call GP_FoX_get_value(attributes, 'bond_real_space', value, status)
            if (status == 0) then
               read (value,*) bond_real_space
            else
               bond_real_space = .false.
!               call system_abort("gpCoordinates_startElement_handler did not find the bond_real_space attribute.")
            endif

            call GP_FoX_get_value(attributes, 'atom_real_space', value, status)
            if (status == 0) then
               read (value,*) atom_real_space
            else
               atom_real_space = .false.
!               call system_abort("gpCoordinates_startElement_handler did not find the atom_real_space attribute.")
            endif

            call GP_FoX_get_value(attributes, 'soap', value, status)
            if (status == 0) then
               read (value,*) soap
            else
               soap = .false.
            endif

            call GP_FoX_get_value(attributes, 'x_size_max', value, status)
            if (status == 0) then
               read (value,*) x_size_max
            else
               if ((bond_real_space .or. atom_real_space) .and. (.not. sparsified)) call system_abort("gpCoordinates_startElement_handler did not find the x_size_max attribute.")
               x_size_max = 0
            endif

            call GP_FoX_get_value(attributes, 'xPrime_size_max', value, status)
            if (status == 0) then
               read (value,*) xPrime_size_max
            else
               if ((bond_real_space .or. atom_real_space) .and. (.not. sparsified)) call system_abort("gpCoordinates_startElement_handler did not find the xPrime_size_max attribute.")
               xPrime_size_max = 0
            endif

            call GP_FoX_get_value(attributes, 'sparseX_size_max', value, status)
            if (status == 0) then
               read (value,*) sparseX_size_max
            else
               if ((bond_real_space .or. atom_real_space) .and. sparsified) call system_abort("gpCoordinates_startElement_handler did not find the sparseX_size_max attribute.")
               sparseX_size_max = 0
            endif

            if(sparsified) then
               call GP_FoX_get_value(attributes, 'n_sparseX', value, status)
               if (status == 0) then
                  read (value,*) n_sparseX
               else
                  call system_abort("gpCoordinates_startElement_handler did not find the n_sparseX attribute.")
               endif

               if (bond_real_space .or. atom_real_space) then
                  call gpCoordinates_setParameters_sparse(parse_gpCoordinates,d,n_sparseX,delta,f0,bond_real_space=bond_real_space,atom_real_space=atom_real_space,sparseX_size_max=sparseX_size_max)
               else
                  call gpCoordinates_setParameters_sparse(parse_gpCoordinates,d,n_sparseX,delta,f0, soap = soap)
		  call GP_FoX_get_value(attributes, 'sparseX_filename', value, status)
		  if (status == 0) then
		     call fread_array_d(size(parse_gpCoordinates%sparseX), parse_gpCoordinates%sparseX(1,1), trim(value)//C_NULL_CHAR)
		     parse_sparseX_separate_file = .true.
		  else
		     parse_sparseX_separate_file = .false.
		  endif
               endif
            else
               call GP_FoX_get_value(attributes, 'n_x', value, status)
               if (status == 0) then
                  read (value,*) n_x
               else
                  call system_abort("gpCoordinates_startElement_handler did not find the n_x attribute.")
               endif

               call GP_FoX_get_value(attributes, 'n_xPrime', value, status)
               if (status == 0) then
                  read (value,*) n_xPrime
               else
                  call system_abort("gpCoordinates_startElement_handler did not find the n_xPrime attribute.")
               endif

               if (bond_real_space .or. atom_real_space) then
                  call gpCoordinates_setParameters(parse_gpCoordinates,d,n_x,n_xPrime,delta,f0,bond_real_space=bond_real_space,atom_real_space=atom_real_space,x_size_max=x_size_max,xPrime_size_max=xPrime_size_max)
               else
                  call gpCoordinates_setParameters(parse_gpCoordinates,d,n_x,n_xPrime,delta,f0,soap=soap)
               endif
            endif

            if (bond_real_space .or. atom_real_space .or. soap) then
               allocate(parse_in_permutations(1,n_permutations))
            else
               allocate(parse_in_permutations(d,n_permutations))
            endif

         endif

      elseif(parse_in_gpCoordinates .and. name == 'theta') then
         call zero(parse_cur_data)
      elseif(parse_in_gpCoordinates .and. name == 'descriptor') then
         call zero(parse_cur_data)
      elseif(parse_in_gpCoordinates .and. name == 'permutation') then

         call GP_FoX_get_value(attributes, 'i', value, status)
         if (status == 0) then
            read (value,*) i
         else
            call system_abort("gpCoordinates_startElement_handler did not find the i attribute.")
         endif

         parse_i_permutation = i

         call zero(parse_cur_data)

      elseif(parse_in_gpCoordinates .and. name == 'sparseX') then

         parse_in_sparseX = .true.

         if( .not. parse_gpCoordinates%sparsified ) then
            call system_abort("gpCoordinates_startElement_handler: not sparsified data and sparseX element found.")
         endif

         call GP_FoX_get_value(attributes, 'i', value, status)
         if (status == 0) then
            read (value,*) i
         else
            call system_abort("gpCoordinates_startElement_handler did not find the i attribute.")
         endif

         call GP_FoX_get_value(attributes, 'alpha', value, status)
         if (status == 0) then
            read (value,*) parse_gpCoordinates%alpha(i)
         else
            call system_abort("gpCoordinates_startElement_handler did not find the alpha attribute.")
         endif

         call GP_FoX_get_value(attributes, 'sparseCutoff', value, status)
         if (status == 0) then
            read (value,*) parse_gpCoordinates%sparseCutoff(i)
         else
            call system_abort("gpCoordinates_startElement_handler did not find the cutoff attribute.")
         endif

         call GP_FoX_get_value(attributes, 'sliced', value, status)
         if (status == 0) then
            read (value,*) parse_sliced
         else
            parse_sliced = .false.
         endif

         if( parse_gpCoordinates%bond_real_space .or. parse_gpCoordinates%atom_real_space ) then
            call GP_FoX_get_value(attributes, 'sparseX_size', value, status)
            if (status == 0) then
               read (value,*) parse_gpCoordinates%sparseX_size(i)
            else
               call system_abort("gpCoordinates_startElement_handler did not find the sparseX_size attribute.")
            endif
         endif

         if( parse_gpCoordinates%bond_real_space .or. parse_gpCoordinates%atom_real_space ) then
            call GP_FoX_get_value(attributes, 'covariance_sparseX_sparseX', value, status)
            if (status == 0) then
               read (value,*) parse_gpCoordinates%covarianceDiag_sparseX_sparseX(i)
            else
               call system_abort("gpCoordinates_startElement_handler did not find the covariance_sparseX_sparseX attribute.")
            endif
         endif

         parse_i_sparseX = i

         call zero(parse_cur_data)

      elseif(parse_in_gpCoordinates .and. parse_in_sparseX .and. name == 'sparseX_slice') then

         call GP_FoX_get_value(attributes, 'start', value, status)
         if (status == 0) then
            read (value,*) parse_slice_start
         else
            call system_abort("gpCoordinates_startElement_handler did not find the start attribute.")
         endif

         call GP_FoX_get_value(attributes, 'end', value, status)
         if (status == 0) then
            read (value,*) parse_slice_end
         else
            call system_abort("gpCoordinates_startElement_handler did not find the end attribute.")
         endif

         call zero(parse_cur_data)
      elseif(parse_in_gpCoordinates .and. name == 'x') then
         if( parse_gpCoordinates%sparsified ) then
            call system_abort("gpCoordinates_startElement_handler: sparsified=T but x element found.")
         endif

         call GP_FoX_get_value(attributes, 'i', value, status)
         if (status == 0) then
            read (value,*) i
         else
            call system_abort("gpCoordinates_startElement_handler did not find the i attribute.")
         endif

         call GP_FoX_get_value(attributes, 'map_x_y', value, status)
         if (status == 0) then
            read (value,*) parse_gpCoordinates%map_x_y(i)
         else
            call system_abort("gpCoordinates_startElement_handler did not find the map_x_y attribute.")
         endif

         if( parse_gpCoordinates%bond_real_space .or. parse_gpCoordinates%atom_real_space ) then
            call GP_FoX_get_value(attributes, 'x_size', value, status)
            if (status == 0) then
               read (value,*) parse_gpCoordinates%x_size(i)
            else
               call system_abort("gpCoordinates_startElement_handler did not find the x_size attribute.")
            endif
         endif

         if( parse_gpCoordinates%bond_real_space .or. parse_gpCoordinates%atom_real_space ) then
            call GP_FoX_get_value(attributes, 'covariance_x_x', value, status)
            if (status == 0) then
               read (value,*) parse_gpCoordinates%covarianceDiag_x_x(i)
            else
               call system_abort("gpCoordinates_startElement_handler did not find the covariance_x_x attribute.")
            endif
         endif

         parse_i_x = i

         call zero(parse_cur_data)

      elseif(parse_in_gpCoordinates .and. name == 'xPrime') then
         if( parse_gpCoordinates%sparsified ) then
            call system_abort("gpCoordinates_startElement_handler: sparsified=T but xPrime element found.")
         endif

         call GP_FoX_get_value(attributes, 'i', value, status)
         if (status == 0) then
            read (value,*) i
         else
            call system_abort("gpCoordinates_startElement_handler did not find the i attribute.")
         endif

         call GP_FoX_get_value(attributes, 'map_xPrime_yPrime', value, status)
         if (status == 0) then
            read (value,*) parse_gpCoordinates%map_xPrime_yPrime(i)
         else
            call system_abort("gpCoordinates_startElement_handler did not find the map_xPrime_yPrime attribute.")
         endif

         call GP_FoX_get_value(attributes, 'map_xPrime_x', value, status)
         if (status == 0) then
            read (value,*) parse_gpCoordinates%map_xPrime_x(i)
         else
            call system_abort("gpCoordinates_startElement_handler did not find the map_xPrime_x attribute.")
         endif

         if( parse_gpCoordinates%bond_real_space .or. parse_gpCoordinates%atom_real_space ) then
            call GP_FoX_get_value(attributes, 'xPrime_size', value, status)
            if (status == 0) then
               read (value,*) parse_gpCoordinates%xPrime_size(i)
            else
               call system_abort("gpCoordinates_startElement_handler did not find the xPrime_size attribute.")
            endif
         endif

         if( parse_gpCoordinates%bond_real_space .or. parse_gpCoordinates%atom_real_space ) then
            call GP_FoX_get_value(attributes, 'covariance_xPrime_xPrime', value, status)
            if (status == 0) then
               read (value,*) parse_gpCoordinates%covarianceDiag_xPrime_xPrime(i)
            else
               call system_abort("gpCoordinates_startElement_handler did not find the covariance_xPrime_xPrime attribute.")
            endif
         endif

         parse_i_xPrime = i

         call zero(parse_cur_data)

      endif

   endsubroutine gpCoordinates_startElement_handler

   subroutine gpCoordinates_endElement_handler(URI, localname, name)
      character(len=*), intent(in)   :: URI
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name

      real(dp), dimension(:), allocatable :: x_array, xyz_array
      integer :: x_array_size, xyz_array_size
      !character(len=100*parse_gpCoordinates%d) :: val
      !character(len=200*100) :: val

      if(parse_in_gpCoordinates) then
         if(name == 'gpCoordinates') then
            call gpCoordinates_setPermutations(parse_gpCoordinates,parse_in_permutations)
            deallocate(parse_in_permutations)
            parse_in_gpCoordinates = .false.
         elseif(name == 'theta') then
            !val = string(parse_cur_data)
            !read(val,*) parse_gpCoordinates%theta
            call string_to_numerical(string(parse_cur_data),parse_gpCoordinates%theta)
         elseif(name == 'descriptor') then
            parse_gpCoordinates%descriptor_str = parse_cur_data
            if(parse_gpCoordinates%atom_real_space) &
               call gpCoordinates_gpCovariance_atom_real_space_Initialise(parse_gpCoordinates)
         elseif(name == 'permutation') then
            
            if( parse_i_permutation > size(parse_in_permutations,2) ) then
               call system_abort("gpCoordinates_endElement_handler: parse_i_permutation ("//parse_i_permutation//") greater than n_permutations ("//size(parse_in_permutations,2)//")")
            endif

            !val = string(parse_cur_data)
            !read(val,*) parse_in_permutations(:,parse_i_permutation)
            call string_to_numerical(string(parse_cur_data),parse_in_permutations(:,parse_i_permutation))
         elseif(name == 'sparseX') then
            
            if( .not. allocated(parse_gpCoordinates%sparseX) ) then
               call system_abort("gpCoordinates_endElement_handler: sparseX not allocated")
            endif
            
            if( parse_i_sparseX > parse_gpCoordinates%n_sparseX ) then
               call system_abort("gpCoordinates_endElement_handler: parse_i_sparseX ("//parse_i_sparseX//") greater than n_sparseX ("//parse_gpCoordinates%n_sparseX//")")
            endif 

            !val = string(parse_cur_data)
            !read(val,*) parse_gpCoordinates%sparseX(:,parse_i_sparseX)
            if( parse_gpCoordinates%bond_real_space ) then
               parse_gpCoordinates%sparseX(:,parse_i_sparseX) = 0.0_dp
               call string_to_numerical(string(parse_cur_data),parse_gpCoordinates%sparseX(:parse_gpCoordinates%sparseX_size(parse_i_sparseX),parse_i_sparseX))
            elseif(parse_gpCoordinates%atom_real_space) then   
               xyz_array_size = parse_gpCoordinates%sparseX_size(parse_i_sparseX)
               allocate(xyz_array(xyz_array_size))

               call string_to_numerical(string(parse_cur_data),xyz_array)

               call gp_atom_real_space_XYZ_RealArray(parse_gpCoordinates%atom_real_space_cov, &
               xyz_array,xyz_array_size,x_array,x_array_size)

               parse_gpCoordinates%sparseX_size(parse_i_sparseX) = x_array_size
               parse_gpCoordinates%sparseX(:x_array_size,parse_i_sparseX) = x_array

               if(allocated(x_array)) deallocate(x_array)
               if(allocated(xyz_array)) deallocate(xyz_array)
            else
               if(.not. parse_sparseX_separate_file .and. .not. parse_sliced) call string_to_numerical(string(parse_cur_data),parse_gpCoordinates%sparseX(:,parse_i_sparseX))
            endif

            parse_in_sparseX = .false.
         elseif(name == 'sparseX_slice') then
            if(parse_slice_start < 1) then
               call system_abort("gpCoordinates_endElement_handler: slice start less than 1")
            endif

            if(parse_slice_end > parse_gpCoordinates%d) then
               call system_abort("gpCoordinates_endElement_handler: slice start greater than dimension")
            endif

            if(.not. parse_sparseX_separate_file .and. parse_sliced) call string_to_numerical(string(parse_cur_data),parse_gpCoordinates%sparseX(parse_slice_start:parse_slice_end,parse_i_sparseX))
         elseif(name == 'x') then

            if( .not. allocated(parse_gpCoordinates%x) ) then
               call system_abort("gpCoordinates_endElement_handler: x not allocated")
            endif
            
            if( parse_i_x > parse_gpCoordinates%n_x ) then
               call system_abort("gpCoordinates_endElement_handler: parse_i_x ("//parse_i_x//") greater than n_x ("//parse_gpCoordinates%n_x//")")
            endif

            !val = string(parse_cur_data)
            !read(val,*) parse_gpCoordinates%x(:,parse_i_x)
            call string_to_numerical(string(parse_cur_data),parse_gpCoordinates%x(:,parse_i_x))
         elseif(name == 'xPrime') then

            if( .not. allocated(parse_gpCoordinates%xPrime) ) then
               call system_abort("gpCoordinates_endElement_handler: xPrime not allocated")
            endif
            
            if( parse_i_xPrime > parse_gpCoordinates%n_xPrime ) then
               call system_abort("gpCoordinates_endElement_handler: parse_i_xPrime ("//parse_i_xPrime//") greater than n_xPrime ("//parse_gpCoordinates%n_xPrime//")")
            endif

            !val = string(parse_cur_data)
            !read(val,*) parse_gpCoordinates%xPrime(:,parse_i_xPrime)
            call string_to_numerical(string(parse_cur_data), parse_gpCoordinates%xPrime(:,parse_i_xPrime))
         endif
      endif

   endsubroutine gpCoordinates_endElement_handler

   subroutine gpCoordinates_characters_handler(in)
      character(len=*), intent(in) :: in

      if(parse_in_gpCoordinates) then
         call concat(parse_cur_data, in, keep_lf=.false.)
      endif
   endsubroutine gpCoordinates_characters_handler

   subroutine gpFull_startElement_handler(URI, localname, name, attributes)
      character(len=*), intent(in)   :: URI
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name
      type(dictionary_t), intent(in) :: attributes

      integer :: status, n_y, n_yPrime, n_coordinate, i
      real(dp) :: sparse_jitter
      character(len=1024) :: value

      if(name == 'gpFull') then ! new GP_data
         if(parse_in_gpFull) then
            call system_abort("gpFull_startElement_handler entered gpFull with parse_in_gpFull true. Probably a bug in FoX (4.0.1, e.g.)")
         endif

         if(parse_matched_label) return ! we already found an exact match for this label

         call GP_FoX_get_value(attributes, 'label', value, status)
         if (status /= 0) value = ''

         if(len(trim(parse_gpFull_label)) > 0) then ! we were passed in a label
            if(trim(value) == trim(parse_gpFull_label)) then
               parse_matched_label = .true.
               parse_in_gpFull = .true.
            else ! no match
               parse_in_gpFull = .false.
            endif
         else ! no label passed in
            parse_in_gpFull = .true.
         endif

         if(parse_in_gpFull) then
            if(parse_gpFull%initialised) call finalise(parse_gpFull)

            call GP_FoX_get_value(attributes, 'n_y', value, status)
            if (status == 0) then
               read (value,*) n_y
            else
               call system_abort("gpFull_startElement_handler did not find the n_y attribute.")
            endif

            call GP_FoX_get_value(attributes, 'n_yPrime', value, status)
            if (status == 0) then
               read (value,*) n_yPrime
            else
               call system_abort("gpFull_startElement_handler did not find the n_yPrime attribute.")
            endif

            call GP_FoX_get_value(attributes, 'n_coordinate', value, status)
            if (status == 0) then
               read (value,*) n_coordinate
            else
               call system_abort("gpFull_startElement_handler did not find the n_coordinate attribute.")
            endif

            call GP_FoX_get_value(attributes, 'sparse_jitter', value, status)
            if (status == 0) then
               read (value,*) sparse_jitter
            else
               call print_warning("gpFull_startElement_handler did not find the sparse_jitter attribute, using default value 1.0e-5.")
               sparse_jitter = 1.0e-5_dp
            endif
            call gpFull_setParameters(parse_gpFull,n_coordinate, n_y, n_yPrime, sparse_jitter)

         endif

      elseif(parse_in_gpFull .and. name == 'y') then

         call GP_FoX_get_value(attributes, 'i', value, status)
         if (status == 0) then
            read (value,*) i
         else
            call system_abort("gpFull_startElement_handler did not find the i attribute.")
         endif

         call GP_FoX_get_value(attributes, 'map_y_globalY', value, status)
         if (status == 0) then
            read (value,*) parse_gpFull%map_y_globalY(i)
         else
            call system_abort("gpFull_startElement_handler did not find the map_y_globalY attribute.")
         endif

         call GP_FoX_get_value(attributes, 'alpha', value, status)
         if (status == 0) then
            read (value,*) parse_gpFull%alpha(parse_gpFull%map_y_globalY(i))
         else
            call system_abort("gpFull_startElement_handler did not find the alpha attribute.")
         endif

      elseif(parse_in_gpFull .and. name == 'yPrime') then

         call GP_FoX_get_value(attributes, 'i', value, status)
         if (status == 0) then
            read (value,*) i
         else
            call system_abort("gpFull_startElement_handler did not find the i attribute.")
         endif

         call GP_FoX_get_value(attributes, 'map_yPrime_globalY', value, status)
         if (status == 0) then
            read (value,*) parse_gpFull%map_yPrime_globalY(i)
         else
            call system_abort("gpFull_startElement_handler did not find the map_yPrime_globalY attribute.")
         endif

         call GP_FoX_get_value(attributes, 'alpha', value, status)
         if (status == 0) then
            read (value,*) parse_gpFull%alpha(parse_gpFull%map_yPrime_globalY(i))
         else
            call system_abort("gpFull_startElement_handler did not find the alpha attribute.")
         endif

      endif

   endsubroutine gpFull_startElement_handler

   subroutine gpFull_endElement_handler(URI, localname, name)
      character(len=*), intent(in)   :: URI
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name

      if(parse_in_gpFull) then
         if(name == 'gpFull') then
            parse_in_gpFull = .false.
         endif
      elseif(name == 'y') then

      elseif(name == 'yPrime') then

      endif

   endsubroutine gpFull_endElement_handler

   subroutine gpFull_characters_handler(in)
      character(len=*), intent(in) :: in

      if(parse_in_gpFull) then
         call concat(parse_cur_data, in, keep_lf=.false.)
      endif
   endsubroutine gpFull_characters_handler

   subroutine gpSparse_startElement_handler(URI, localname, name, attributes)
      character(len=*), intent(in)   :: URI
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name
      type(dictionary_t), intent(in) :: attributes

      integer :: status, n_coordinate
      character(len=1024) :: value

      if(name == 'gpSparse') then ! new GP_data
         if(parse_in_gpSparse) then
            call system_abort("gpSparse_startElement_handler entered gpSparse with parse_in_gpSparse true. Probably a bug in FoX (4.0.1, e.g.)")
         endif

         if(parse_matched_label) return ! we already found an exact match for this label

         call GP_FoX_get_value(attributes, 'label', value, status)
         if (status /= 0) value = ''

         if(len(trim(parse_gpSparse_label)) > 0) then ! we were passed in a label
            if(trim(value) == trim(parse_gpSparse_label)) then
               parse_matched_label = .true.
               parse_in_gpSparse = .true.
            else ! no match
               parse_in_gpSparse = .false.
            endif
         else ! no label passed in
            parse_in_gpSparse = .true.
         endif

         if(parse_in_gpSparse) then
            if(parse_gpSparse%initialised) call finalise(parse_gpSparse)

            call GP_FoX_get_value(attributes, 'n_coordinate', value, status)
            if (status == 0) then
               read (value,*) n_coordinate
            else
               call system_abort("gpSparse_startElement_handler did not find the n_coordinate attribute.")
            endif

            call gpSparse_setParameters(parse_gpSparse,n_coordinate)

         endif

      endif

   endsubroutine gpSparse_startElement_handler

   subroutine gpSparse_endElement_handler(URI, localname, name)
      character(len=*), intent(in)   :: URI
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name

      if(parse_in_gpSparse) then
         if(name == 'gpSparse') then
            parse_in_gpSparse = .false.
         endif
      endif

   endsubroutine gpSparse_endElement_handler

   subroutine gpSparse_characters_handler(in)
      character(len=*), intent(in) :: in

      if(parse_in_gpSparse) then
         call concat(parse_cur_data, in, keep_lf=.false.)
      endif
   endsubroutine gpSparse_characters_handler

   subroutine gp_FoX_get_value(attributes, key, val, status)
     type(dictionary_t), intent(in) :: attributes
     character(len=*), intent(in) :: key
     character(len=*), intent(inout) :: val
     integer, intent(out), optional :: status
            
     if (HasKey(attributes,key)) then
       val = GetValue(attributes, trim(key))
       if (present(status)) status = 0
     else
       val = "" 
       if (present(status)) status = 1
     endif
   end subroutine gp_FoX_get_value

end module gp_predict_module
