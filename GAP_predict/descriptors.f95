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

module descriptors_module

   use libatoms_module

   implicit none

   real(dp), parameter :: QW_FP_ZERO = 1.0E-12_dp
   integer, parameter :: FOURIER_N_COMP = 10
   integer, parameter :: WATER_DIMER_D = 2*FOURIER_N_COMP+1

   type real_3
        real(dp), dimension(3) :: x
   endtype real_3

   type cplx_3
        complex(dp), dimension(3) :: x
   endtype cplx_3

   type cplx_2d_x
      complex(dp), dimension(:,:,:), allocatable :: mm
   endtype cplx_2d_x

   type cplx_2d
      complex(dp), dimension(:,:), allocatable :: mm
   endtype cplx_2d

   type cplx_1d_x
      complex(dp), dimension(:,:), allocatable :: m
   endtype cplx_1d_x

   type cplx_1d
      complex(dp), dimension(:), allocatable :: m
   endtype cplx_1d

   type fourier_so4
      integer :: j_max
      real(dp) :: z0, cutoff
      type(cplx_2d), dimension(:), allocatable :: U      
      logical                                  :: initialised = .false.
   endtype fourier_so4

   type grad_fourier_so4
      integer :: j_max
      real(dp) :: z0, cutoff
      type(cplx_2d_x), dimension(:), allocatable :: U      
      logical                                  :: initialised = .false.
   endtype grad_fourier_so4

   type fourier_so3
      integer :: l_max
      real(dp), dimension(:), allocatable :: cutoff, cutoff_r1
      integer, dimension(:), allocatable :: cutoff_f
      type(cplx_1d), dimension(:,:), allocatable :: Y_times_R
      logical :: initialised = .false.
   endtype fourier_so3

   type grad_fourier_so3
      integer :: l_max
      real(dp), dimension(:), allocatable :: cutoff, cutoff_r1
      integer, dimension(:), allocatable :: cutoff_f
      type(cplx_1d_x), dimension(:,:), allocatable :: Y_times_R
      logical :: initialised = .false.
   endtype grad_fourier_so3

   type magnetic
      complex(dp), dimension(:), allocatable :: mag
   endtype magnetic

   type bispectrum_so4
      integer :: j_max
      type(magnetic), dimension(:,:), allocatable :: coeff
      logical                                  :: initialised = .false.
   endtype bispectrum_so4

   type magnetic_x
      type(cplx_3), dimension(:), allocatable :: mag
   endtype magnetic_x

   type grad_bispectrum_so4
      integer :: j_max
      type(magnetic_x), dimension(:,:), allocatable :: coeff
      logical                                  :: initialised = .false.
   endtype grad_bispectrum_so4

   type qw_so3
      integer :: l_max, f_n
      logical :: do_q, do_w
      real(dp), dimension(:,:), allocatable :: q
      real(dp), dimension(:,:), allocatable :: w
      logical :: initialised = .false.
   endtype qw_so3

   type grad_qw_so3
      integer :: l_max, f_n
      logical :: do_q, do_w
      type(real_3), dimension(:,:), allocatable :: q
      type(real_3), dimension(:,:), allocatable :: w
      logical :: initialised = .false.
   endtype grad_qw_so3

   type angular_fourier
      type(magnetic), dimension(:), allocatable  :: ang
   endtype angular_fourier

   type radial_fourier
      type(angular_fourier), dimension(:), allocatable :: rad
   endtype radial_fourier

   type angular_bispectrum
      type(magnetic), dimension(:,:), allocatable :: ang
   endtype angular_bispectrum

   type radial_bispectrum
      type(angular_bispectrum), dimension(:,:), allocatable :: rad
   endtype radial_bispectrum

   type symm_radial
      real(dp) :: eta
      real(dp) :: rs
   endtype symm_radial

   type symm_angular
      real(dp) :: zeta
      real(dp) :: lambda
      real(dp) :: eta
      real(dp) :: rj
      real(dp) :: rk
   endtype symm_angular

   type symm
      logical                                       :: initialised = .false.
      real(dp)                                      :: cutoff
      type(symm_radial), dimension(:), allocatable  :: radial
      type(symm_angular), dimension(:), allocatable :: angular
      integer                                       :: n_radial = 0
      integer                                       :: n_angular = 0
   endtype symm

   type fourier_coefficient_atom
      logical                             :: initialised = .false.
      integer                             :: N_max = 0
      integer                             :: L_max = 0
      type(radial_fourier)                :: coeff
      real(dp), dimension(:), allocatable :: r0
      real(dp)                            :: sigma, cutoff_left, cutoff_right, cutoff_sigma
   endtype fourier_coefficient_atom

   type fourier_coefficient_pair
      logical                                    :: initialised = .false.
      integer                                    :: N_max = 0
      integer                                    :: L_max = 0
      complex(dp), dimension(:,:,:), allocatable :: coeff
      real(dp), dimension(:), allocatable        :: r0
      real(dp)                                   :: r_ij, r_cut, sigma, cutoff_left, cutoff_right, cutoff_sigma
   endtype fourier_coefficient_pair

   type bispectrum_coefficient_atom
      logical                             :: initialised = .false.
      integer                             :: N_max = 0
      integer                             :: L_max = 0
      type(radial_bispectrum)             :: coeff
   endtype bispectrum_coefficient_atom

   type bispectrum_coefficient_pair
      logical                                    :: initialised = .false.
      integer                                    :: N_max = 0
      integer                                    :: L_max = 0
      real(dp)                                   :: r_ij
      complex(dp), dimension(:,:,:), allocatable :: coeff
   endtype bispectrum_coefficient_pair

   type per
      complex(dp), dimension(:,:,:), allocatable :: f
      integer, dimension(3)                      :: k
      real(dp)                                   :: sig = 0.2_dp
      real(dp), dimension(:), allocatable        :: w
      logical                                    :: initialised = .false.
   endtype per

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

   interface initialise
      module procedure initialise_fourier_atom, initialise_bispectrum_atom, &
      & initialise_fourier_pair, initialise_bispectrum_pair, initialise_symm, &
      & initialise_fourier_so4, initialise_bispectrum_so4, &
      & initialise_grad_fourier_so4, initialise_grad_bispectrum_so4, &
      & initialise_fourier_periodic, &
      & initialise_fourier_so3, initialise_qw_so3, &
      & initialise_grad_fourier_so3, initialise_grad_qw_so3
   endinterface initialise

   interface finalise
      module procedure finalise_fourier_atom, finalise_bispectrum_atom, &
      & finalise_fourier_pair, finalise_bispectrum_pair, finalise_symm, &
      & finalise_fourier_so4, finalise_bispectrum_so4, finalise_grad_fourier_so4, &
      & finalise_grad_bispectrum_so4, finalise_fourier_periodic, &
      & finalise_fourier_so3, finalise_qw_so3, &
      & finalise_grad_fourier_so3, finalise_grad_qw_so3
   endinterface finalise

   interface print
      module procedure print_fourier_atom, print_bispectrum_atom, &
      & print_fourier_pair, print_bispectrum_pair, print_bispectrum_so4
   endinterface print

   interface fourier_transform
      module procedure fourier_transform_atom, fourier_transform_vec, fourier_transform_pair, &
      & fourier_transform_periodic, fourier_transform_so4, grad_fourier_transform_so4, &
      & fourier_transform_so3, grad_fourier_transform_so3
   endinterface fourier_transform

   interface calc_bispectrum
      module procedure calc_bispectrum1_atom, calc_bispectrum_pair, calc_bispectrum_so4, calc_grad_bispectrum_so4, &
      & bispectrum_periodic
   endinterface calc_bispectrum

   interface calc_qw
      module procedure calc_qw_so3, calc_grad_qw_so3
   endinterface calc_qw

   interface bispectrum2vec
      module procedure bispectrum2vec1_atom, bispectrum2vec2_atom, bispectrum2vec3_atom, &
      & bispectrum2vec1_pair, bispectrum2vec_so4, bispectrum2vec_grad_so4, bispectrum2vec_periodic
   endinterface bispectrum2vec

   interface qw2vec
      module procedure qw2vec_qw_so3, qw2vec_grad_qw_so3
   endinterface qw2vec

   interface qw2d
      module procedure qw2d_qw_so3, qw2d_grad_qw_so3
   end interface

   interface L_max2d
      module procedure L_max2d_bispectrum_atom, L_max2d_bispectrum_pair
   endinterface L_max2d

   interface cg
      module procedure cg_lookup !cg_calculate
   endinterface cg

   private

   public :: symm, atom2symm, symm_jacobian, per, get_weights

   public :: calc_qw_at
   public :: fourier_transform, initialise, finalise, calc_bispectrum, calc_qw, L_max2d, L_max2d_atom, L_max2d_pair, &
   & bispectrum2vec, qw2vec, qw2d, print

   public :: fourier_coefficient_atom, fourier_coefficient_pair, &
   & bispectrum_coefficient_atom, bispectrum_coefficient_pair
   !NB gfortran 4.3.4 doesn't like private type members of a public type
   public :: radial_bispectrum, magnetic, radial_fourier, cplx_1d, cplx_2d, magnetic_x, cplx_1d_x, &
     cplx_2d_x, real_3, cplx_3, angular_bispectrum, angular_fourier, symm_radial, symm_angular
   !NB

   public :: test_qw_gradient

   public :: order_lattice_vectors
   public :: transpose_bispectrum_pair, conjugate_bispectrum_pair
   public :: cg, cg_initialise, cg_finalise, wigner3j
   public :: gaussian, fermi_dirac
   public :: xyz2spherical, xyz2ell
   public :: SphericalY
   public :: LegendreP
   public :: SolidRCartesian
   public :: SphericalYCartesian
   public :: GradSphericalYCartesian
   public :: RadialFunction
   public :: GradRadialFunction
   public :: p_norm
   public :: factorial, factorial2, binom, oscillate
   public :: fourier_so4, bispectrum_so4, j_max2d, grad_fourier_so4, grad_bispectrum_so4
   public :: fourier_so3, grad_fourier_so3, qw_so3, grad_qw_so3
   public :: wigner_big_U, grad_wigner_big_U, fourier_transform_so4_old, reduce_lattice
   public :: kmax2d, suggest_kmax, suggest_resolution
   public :: water_monomer, water_dimer, WATER_DIMER_D

   contains

     subroutine initialise_fourier_so4(this,j_max,z0,cutoff)

        type(fourier_so4), intent(inout) :: this
        integer, intent(in) :: j_max
        real(dp), intent(in) :: z0, cutoff

        integer :: j

        if( this%initialised ) call finalise(this)

        this%z0 = z0
        this%cutoff = cutoff
        this%j_max = j_max

        allocate( this%U(0:j_max) )
        do j = 0, j_max
           allocate( this%U(j)%mm(-j:j,-j:j) )
           this%U(j)%mm = CPLX_ZERO
        enddo
        this%initialised = .true.
     
     endsubroutine initialise_fourier_so4

     subroutine initialise_grad_fourier_so4(this,j_max,z0,cutoff)

        type(grad_fourier_so4), intent(inout) :: this
        integer, intent(in) :: j_max
        real(dp), intent(in) :: z0, cutoff

        integer :: j

        if( this%initialised ) call finalise(this)

        this%z0 = z0
        this%cutoff = cutoff
        this%j_max = j_max

        allocate( this%U(0:j_max) )
        do j = 0, j_max
           allocate( this%U(j)%mm(3,-j:j,-j:j) )
           this%U(j)%mm = CPLX_ZERO
        enddo
        this%initialised = .true.
     
     endsubroutine initialise_grad_fourier_so4

     subroutine initialise_fourier_so3(this, l_max, cutoff, cutoff_f, cutoff_r1)

       type(fourier_so3), intent(inout) :: this
       integer, intent(in) :: l_max
       real(dp), dimension(:), intent(in) :: cutoff, cutoff_r1
       integer, dimension(:), intent(in) :: cutoff_f

       integer :: l, f

       if (mod(l_max, 2) /= 0) call system_abort('initialise_fourier_so3: l_max is odd')

       if (this%initialised) call finalise(this)

       this%l_max = l_max
       allocate(this%cutoff(size(cutoff)), this%cutoff_f(size(cutoff_f)), this%cutoff_r1(size(cutoff_r1)))
       this%cutoff = cutoff
       this%cutoff_f = cutoff_f
       this%cutoff_r1 = cutoff_r1

       allocate(this%Y_times_R(1:l_max / 2,1:size(cutoff)))
       do l = 2, l_max, 2
          do f = 1, size(cutoff)
             allocate(this%Y_times_R(l / 2,f)%m(-l:l))
             this%Y_times_R(l / 2,f)%m = CPLX_ZERO
          enddo
       enddo
       this%initialised = .true.
     
     endsubroutine initialise_fourier_so3

     subroutine initialise_grad_fourier_so3(this, l_max, cutoff, cutoff_f, cutoff_r1)

       type(grad_fourier_so3), intent(inout) :: this
       integer, intent(in) :: l_max
       real(dp), dimension(:), intent(in) :: cutoff, cutoff_r1
       integer, dimension(:), intent(in) :: cutoff_f

       integer :: l, f

       if (mod(l_max, 2) /= 0) call system_abort('initialise_fourier_so3: l_max is odd')

       if (this%initialised) call finalise(this)

       this%l_max = l_max
       allocate(this%cutoff(size(cutoff)), this%cutoff_f(size(cutoff_f)), this%cutoff_r1(size(cutoff_r1)))
       this%cutoff = cutoff
       this%cutoff_f = cutoff_f
       this%cutoff_r1 = cutoff_r1

       allocate(this%Y_times_R(1:l_max / 2,1:size(cutoff)))
       do l = 2, l_max, 2
          do f = 1, size(cutoff)
             allocate(this%Y_times_R(l / 2,f)%m(3,-l:l))
             this%Y_times_R(l / 2,f)%m = CPLX_ZERO
          enddo
       enddo
       this%initialised = .true.

     endsubroutine initialise_grad_fourier_so3

     subroutine initialise_bispectrum_so4(this,j_max)

        type(bispectrum_so4), intent(inout) :: this
        integer, intent(in) :: j_max

        integer :: j1, j2

        if( this%initialised ) call finalise(this)

        this%j_max = j_max

        allocate( this%coeff(0:this%j_max,0:this%j_max) )
        do j1 = 0, this%j_max
           do j2 = 0, this%j_max
               allocate( this%coeff(j2,j1)%mag(abs(j1-j2):min(j1+j2,this%j_max)) )
               this%coeff(j2,j1)%mag(:) = CPLX_ZERO
           enddo
        enddo

        this%initialised = .true.
     
     endsubroutine initialise_bispectrum_so4

     subroutine initialise_grad_bispectrum_so4(this,j_max)

        type(grad_bispectrum_so4), intent(inout) :: this
        integer, intent(in) :: j_max

        integer :: j1, j2, j

        if( this%initialised ) call finalise(this)

        this%j_max = j_max

        allocate( this%coeff(0:this%j_max,0:this%j_max) )
        do j1 = 0, this%j_max
           do j2 = 0, this%j_max
               allocate( this%coeff(j2,j1)%mag(abs(j1-j2):min(j1+j2,this%j_max)) )
               do j = abs(j1-j2), min(this%j_max,j1+j2)
                  this%coeff(j2,j1)%mag(j)%x(:) = CPLX_ZERO
               enddo
           enddo
        enddo

        this%initialised = .true.
     
     endsubroutine initialise_grad_bispectrum_so4

     subroutine initialise_qw_so3(this, l_max, f_n, do_q, do_w)

       type(qw_so3), intent(inout) :: this
       integer, intent(in) :: l_max, f_n
       logical, intent(in), optional :: do_q, do_w

       logical :: my_do_q, my_do_w

       if (mod(l_max, 2) /= 0) call system_abort('initialise_qw_so3: l_max is odd')
       if (f_n == 0) call system_abort('initialise_qw_so3: need at least 1 radial function')

       if (this%initialised) call finalise(this)

       my_do_q = optional_default(.true., do_q)
       my_do_w = optional_default(.true., do_w)

       this%l_max = l_max
       this%f_n = f_n
       this%do_q = my_do_q
       this%do_w = my_do_w

       if ((.not. my_do_q) .and. (.not. my_do_w)) call system_abort('initialise_qw_so3: nothing to initialise')

       if (my_do_q) allocate(this%q(1:l_max / 2,1:f_n))
       if (my_do_w) allocate(this%w(1:l_max / 2,1:f_n))

       if (my_do_q) this%q(:,:) = CPLX_ZERO
       if (my_do_w) this%w(:,:) = CPLX_ZERO

       this%initialised = .true.

     endsubroutine initialise_qw_so3

     subroutine initialise_grad_qw_so3(this, l_max, f_n, do_q, do_w)

       type(grad_qw_so3), intent(inout) :: this
       integer, intent(in) :: l_max, f_n
       logical, intent(in), optional :: do_q, do_w

       logical :: my_do_q, my_do_w
       integer :: l, f

       if (mod(l_max, 2) /= 0) call system_abort('initialise_grad_qw_so3: l_max is odd')
       if (f_n == 0) call system_abort('initialise_grad_qw_so3: need at least 1 radial function')

       if (this%initialised) call finalise(this)

       my_do_q = optional_default(.true., do_q)
       my_do_w = optional_default(.true., do_w)

       this%l_max = l_max
       this%f_n = f_n
       this%do_q = my_do_q
       this%do_w = my_do_w

       if ((.not. my_do_q) .and. (.not. my_do_w)) call system_abort('initialise_grad_qw_so3: nothing to initialise')

       if (my_do_q) allocate(this%q(1:l_max / 2,1:f_n))
       if (my_do_w) allocate(this%w(1:l_max / 2,1:f_n))

       do l = 2, (l_max / 2)
          do f = 1, f_n
             if (my_do_q) this%q(l,f)%x(:) = CPLX_ZERO
             if (my_do_w) this%w(l,f)%x(:) = CPLX_ZERO
          enddo
       enddo

       this%initialised = .true.

     endsubroutine initialise_grad_qw_so3

     subroutine finalise_fourier_so4(this)
        type(fourier_so4), intent(inout) :: this

        integer :: j

        if( .not. this%initialised ) return
        do j = 0, this%j_max
           deallocate( this%U(j)%mm )
        enddo
        deallocate( this%U )

        this%z0 = 0.0_dp
        this%cutoff = 0.0_dp
        this%j_max = 0
        this%initialised = .false.

     endsubroutine finalise_fourier_so4
     
     subroutine finalise_grad_fourier_so4(this)
        type(grad_fourier_so4), intent(inout) :: this

        integer :: j

        if( .not. this%initialised ) return
        do j = 0, this%j_max
           deallocate( this%U(j)%mm )
        enddo
        deallocate( this%U )

        this%z0 = 0.0_dp
        this%cutoff = 0.0_dp
        this%j_max = 0
        this%initialised = .false.

     endsubroutine finalise_grad_fourier_so4

     subroutine finalise_fourier_so3(this)
        type(fourier_so3), intent(inout) :: this

        integer :: l, f

        if (.not. this%initialised) return
        do l = 1, (this%l_max / 2)
           do f = 1, size(this%cutoff)
              deallocate(this%Y_times_R(l,f)%m)
           enddo
        enddo
        deallocate(this%Y_times_R)

        deallocate(this%cutoff, this%cutoff_f, this%cutoff_r1)
        this%l_max = 0
        this%initialised = .false.

     endsubroutine finalise_fourier_so3
     
     subroutine finalise_grad_fourier_so3(this)
        type(grad_fourier_so3), intent(inout) :: this

        integer :: l, f

        if (.not. this%initialised) return
        do l = 1, (this%l_max / 2)
           do f = 1, size(this%cutoff)
              deallocate(this%Y_times_R(l,f)%m)
           enddo
        enddo
        deallocate(this%Y_times_R)

        deallocate(this%cutoff, this%cutoff_f, this%cutoff_r1)
        this%l_max = 0
        this%initialised = .false.

     endsubroutine finalise_grad_fourier_so3
     
     subroutine finalise_bispectrum_so4(this)
        type(bispectrum_so4), intent(inout) :: this

        integer :: j1, j2

        if( .not. this%initialised ) return

        do j1 = 0, this%j_max
           do j2 = 0, this%j_max
               deallocate( this%coeff(j2,j1)%mag )
           enddo
        enddo
        deallocate( this%coeff )

        this%j_max = 0
        this%initialised = .false.

     endsubroutine finalise_bispectrum_so4
     
     subroutine finalise_grad_bispectrum_so4(this)
        type(grad_bispectrum_so4), intent(inout) :: this

        integer :: j1, j2

        if( .not. this%initialised ) return

        do j1 = 0, this%j_max
           do j2 = 0, this%j_max
               deallocate( this%coeff(j2,j1)%mag )
           enddo
        enddo
        deallocate( this%coeff )

        this%j_max = 0
        this%initialised = .false.

     endsubroutine finalise_grad_bispectrum_so4

     subroutine finalise_qw_so3(this)
        type(qw_so3), intent(inout) :: this

        if (.not. this%initialised) return

        if (allocated(this%q)) deallocate(this%q)
        if (allocated(this%w)) deallocate(this%w)

        this%l_max = 0
        this%f_n = 0
        this%do_q = .false.
        this%do_w = .false.
        this%initialised = .false.

     endsubroutine finalise_qw_so3
     
     subroutine finalise_grad_qw_so3(this)
        type(grad_qw_so3), intent(inout) :: this

        if (.not. this%initialised) return

        if (allocated(this%q)) deallocate(this%q)
        if (allocated(this%w)) deallocate(this%w)

        this%l_max = 0
        this%f_n = 0
        this%do_q = .false.
        this%do_w = .false.
        this%initialised = .false.

     endsubroutine finalise_grad_qw_so3
     
     subroutine fourier_transform_so4(this,at,i,w)
        type(fourier_so4), intent(inout) :: this
        type(atoms), intent(in) :: at
        integer, intent(in) :: i
        real(dp), dimension(:), intent(in), optional :: w

        integer :: j, n, jn, m1, m2
        real(dp) :: r0, r, cutoff, z0, theta0
        real(dp), dimension(3) :: diff
        complex(dp) :: z0_pls_Iz, z0_min_Iz, x_pls_Iy, x_min_Iy
        complex(dp), dimension(:,:), allocatable :: U, Up
        real(dp), dimension(:), allocatable :: my_w

        if( .not. this%initialised ) call system_abort('fourier_transform_so4: fourier_so4 not initialised')

        allocate(my_w(at%N))
        my_w = 1.0_dp
        if(present(w)) then
           if (size(w) /= at%N) call system_abort('fourier_transform_so4: size of w is not at%N')
           my_w = w
        endif

        do j = 0, this%j_max
           this%U(j)%mm(:,:) = CPLX_ZERO
           do m1 = -j, j, 2
              this%U(j)%mm(m1,m1) = CPLX_ONE * my_w(i)
           enddo
        enddo

        allocate( U(-this%j_max:this%j_max, -this%j_max:this%j_max), &
        & Up(-this%j_max:this%j_max, -this%j_max:this%j_max) )
        U = CPLX_ZERO
        Up = CPLX_ZERO

        do n = 1, atoms_n_neighbours(at,i)
           jn = atoms_neighbour(at, i, n, distance=r, diff=diff)
           if( r > this%cutoff ) cycle

!           r0 = sqrt( r**2 + this%z0**2 )
!
!           diff = diff * this%z0*sin( 2.0_dp * atan(r/this%z0) ) / r
!           z0 = this%z0*cos( 2.0_dp * atan(r/this%z0) )
!           r0 = sqrt( sum(diff**2) + z0**2 )
           theta0 = r / this%z0 ! + PI/3
           
           z0 = r / tan( theta0 )
           r0 = r / sin( theta0 )

           z0_pls_Iz = ( z0 + CPLX_IMAG*diff(3) ) / r0
           z0_min_Iz = ( z0 - CPLX_IMAG*diff(3) ) / r0
           x_pls_Iy = ( diff(1) + CPLX_IMAG*diff(2) ) / r0
           x_min_Iy = ( diff(1) - CPLX_IMAG*diff(2) ) / r0

!           z0_pls_Iz = ( this%z0 + CPLX_IMAG*diff(3) ) / r0
!           z0_min_Iz = ( this%z0 - CPLX_IMAG*diff(3) ) / r0
!           x_pls_Iy = ( diff(1) + CPLX_IMAG*diff(2) ) / r0
!           x_min_Iy = ( diff(1) - CPLX_IMAG*diff(2) ) / r0
           cutoff = coff(r,this%cutoff) * my_w(jn)
!           cutoff = 1.0_dp

           this%U(0)%mm(0,0) = this%U(0)%mm(0,0) + cutoff
           Up(0:0,0:0) = CPLX_ONE

           do j = 1, this%j_max
              U(-j:j,-j:j) = CPLX_ZERO

              do m1 = -j, j-2, 2
                 do m2 = -j, j, 2
                    if( (j-m2) /= 0 ) U(m2,m1) = U(m2,m1) + &
                    & sqrt( real(j-m2,dp)/real(j-m1,dp) ) * z0_pls_Iz * Up(m2+1,m1+1)

                    if( (j+m2) /= 0 ) U(m2,m1) = U(m2,m1) - &
                    & CPLX_IMAG * sqrt( real(j+m2,dp)/real(j-m1,dp) ) * x_min_Iy * &
                    & Up(m2-1,m1+1)
                 enddo
              enddo

              m1 = j
              do m2 = -j, j, 2
                 if( (j+m2) /= 0 ) U(m2,m1) = U(m2,m1) + &
                 & sqrt( real(j+m2,dp)/real(j+m1,dp) ) * z0_min_Iz * Up(m2-1,m1-1)
                 if( (j-m2) /= 0 ) U(m2,m1) = U(m2,m1) - &
                 & CPLX_IMAG * sqrt( real(j-m2,dp)/real(j+m1,dp) ) * x_pls_Iy * &
                 & Up(m2+1,m1-1)
              enddo

              this%U(j)%mm = this%U(j)%mm + U(-j:j,-j:j) * cutoff !* sqrt( 2.0_dp * j + 1.0_dp )
              Up(-j:j,-j:j) = U(-j:j,-j:j)
           enddo
        enddo

        deallocate( U, Up, my_w )

     endsubroutine fourier_transform_so4

     subroutine grad_fourier_transform_so4(this,at,i,ni,w)
        type(grad_fourier_so4), intent(inout) :: this
        type(atoms), intent(in) :: at
        integer, intent(in) :: i, ni
        real(dp), dimension(:), intent(in), optional :: w

        integer :: j, n, jn, m1, m2
        real(dp) :: r0, r, cutoff, theta0, z0
        real(dp), dimension(3) :: diff, u_ij, dcutoff, dz0, dr0
        complex(dp) :: z0_pls_Iz, z0_min_Iz, x_pls_Iy, x_min_Iy
        complex(dp), dimension(3) :: dz0_pls_Iz, dz0_min_Iz, dx_pls_Iy, dx_min_Iy
        complex(dp), dimension(:,:), allocatable :: U, Up
        complex(dp), dimension(:,:,:), allocatable :: dU, dUp
        real(dp), dimension(:), allocatable :: my_w

        if( .not. this%initialised ) call system_abort('fourier_transform_so4: fourier_so4 not initialised')

        allocate(my_w(at%N))
        my_w = 1.0_dp
        if(present(w)) then
           if (size(w) /= at%N) call system_abort('fourier_transform_so4: size of w is not at%N')
           my_w = w
        endif

        do j = 0, this%j_max
           this%U(j)%mm(:,:,:) = CPLX_ZERO
        enddo

        allocate( U(-this%j_max:this%j_max, -this%j_max:this%j_max), &
        & Up(-this%j_max:this%j_max, -this%j_max:this%j_max), &
        & dU(3,-this%j_max:this%j_max, -this%j_max:this%j_max), &
        & dUp(3,-this%j_max:this%j_max, -this%j_max:this%j_max) )

        U = CPLX_ZERO
        Up = CPLX_ZERO
        dU = CPLX_ZERO
        dUp = CPLX_ZERO

        if(ni==0) then
           do n = 1, atoms_n_neighbours(at,i)
              jn = atoms_neighbour(at, i, n, distance=r, diff=diff, cosines=u_ij )
              
              if( r > this%cutoff ) cycle
              !if( i==jn ) cycle

              !r0 = sqrt( r**2 + this%z0**2 )
              cutoff = coff(r,this%cutoff) * my_w(jn)
              dcutoff = -dcoff(r,this%cutoff)*u_ij * my_w(jn)

              theta0 = r / this%z0 !+ PI/3

              z0 = r / tan( theta0 )
              r0 = sin( theta0 ) / r

              dz0 = ( 1.0_dp / tan( theta0 ) - theta0 / sin(theta0)**2 ) * u_ij
              !dz0 = ( 1.0_dp / tan( theta0 ) - (theta0-pi/3) / sin(theta0)**2 ) * u_ij
              dr0 = ( cos( theta0 ) / (r*this%z0) - r0 / r ) * u_ij

              z0_pls_Iz = ( z0 + CPLX_IMAG*diff(3) ) * r0
              z0_min_Iz = ( z0 - CPLX_IMAG*diff(3) ) * r0
              x_pls_Iy = ( diff(1) + CPLX_IMAG*diff(2) ) * r0
              x_min_Iy = ( diff(1) - CPLX_IMAG*diff(2) ) * r0

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

              !dz0_pls_Iz = - z0_pls_Iz * diff / r0**2
              !dz0_pls_Iz(3) = dz0_pls_Iz(3) + CPLX_IMAG/r0
              !dz0_min_Iz = - z0_min_Iz * diff / r0**2
              !dz0_min_Iz(3) = dz0_min_Iz(3) - CPLX_IMAG/r0
              !dx_pls_Iy = - x_pls_Iy * diff / r0**2
              !dx_pls_Iy(1) = dx_pls_Iy(1) + 1.0_dp/r0
              !dx_pls_Iy(2) = dx_pls_Iy(2) + CPLX_IMAG/r0
              !dx_min_Iy = - x_min_Iy * diff / r0**2
              !dx_min_Iy(1) = dx_min_Iy(1) + 1.0_dp/r0
              !dx_min_Iy(2) = dx_min_Iy(2) - CPLX_IMAG/r0

              Up(0:0,0:0) = CPLX_ONE
              dUp(:,0:0,0:0) = CPLX_ZERO

              this%U(0)%mm(:,0,0) = this%U(0)%mm(:,0,0) + dcutoff
              do j = 1, this%j_max
                 U(-j:j,-j:j) = CPLX_ZERO
                 dU(:,-j:j,-j:j) = CPLX_ZERO

                 do m1 = -j, j-2, 2
                    do m2 = -j, j, 2
                       if( (j-m2) /= 0 ) then
                          U(m2,m1) = U(m2,m1) + &
                          & sqrt( real(j-m2,dp)/real(j-m1,dp) ) * z0_pls_Iz * Up(m2+1,m1+1)
                          dU(:,m2,m1) = dU(:,m2,m1) + &
                          & sqrt( real(j-m2,dp)/real(j-m1,dp) ) * &
                          & ( dz0_pls_Iz * Up(m2+1,m1+1) + z0_pls_Iz * dUp(:,m2+1,m1+1) )
                       endif

                       if( (j+m2) /= 0 ) then
                          U(m2,m1) = U(m2,m1) - &
                          & CPLX_IMAG * sqrt( real(j+m2,dp)/real(j-m1,dp) ) * x_min_Iy * &
                          & Up(m2-1,m1+1)
                          dU(:,m2,m1) = dU(:,m2,m1) - &
                          & CPLX_IMAG * sqrt( real(j+m2,dp)/real(j-m1,dp) ) * &
                          & ( dx_min_Iy * Up(m2-1,m1+1) + x_min_Iy * dUp(:,m2-1,m1+1) )
                       endif
                    enddo
                 enddo
                          
                 m1 = j
                 do m2 = -j, j, 2
                    if( (j+m2) /= 0 ) then
                       U(m2,m1) = U(m2,m1) + &
                       & sqrt( real(j+m2,dp)/real(j+m1,dp) ) * z0_min_Iz * Up(m2-1,m1-1)
                       dU(:,m2,m1) = dU(:,m2,m1) + &
                       & sqrt( real(j+m2,dp)/real(j+m1,dp) ) * &
                       & ( dz0_min_Iz * Up(m2-1,m1-1) + z0_min_Iz * dUp(:,m2-1,m1-1) )
                    endif
                    if( (j-m2) /= 0 ) then
                       U(m2,m1) = U(m2,m1) - &
                       & CPLX_IMAG * sqrt( real(j-m2,dp)/real(j+m1,dp) ) * x_pls_Iy * &
                       & Up(m2+1,m1-1)
                       dU(:,m2,m1) = dU(:,m2,m1) - &
                       & CPLX_IMAG * sqrt( real(j-m2,dp)/real(j+m1,dp) ) * &
                       & ( dx_pls_Iy * Up(m2+1,m1-1) + x_pls_Iy * dUp(:,m2+1,m1-1) )
                    endif
                 enddo

                 Up(-j:j,-j:j) = U(-j:j,-j:j)
                 dUp(:,-j:j,-j:j) = dU(:,-j:j,-j:j)
                 this%U(j)%mm = this%U(j)%mm - dU(:,-j:j,-j:j) * cutoff 
                 do m1 = -j, j, 2
                    do m2 = -j, j, 2
                       this%U(j)%mm(:,m2,m1) = this%U(j)%mm(:,m2,m1) &
                       & + U(m2,m1) * dcutoff
                    enddo
                 enddo
              enddo
           enddo
        else
           jn = atoms_neighbour(at,i,ni,diff=diff,distance=r,cosines=u_ij)
           if( r < this%cutoff ) then
              theta0 = r / this%z0 !+ PI/3

              z0 = r / tan( theta0 )
              r0 = sin( theta0 ) / r

              dz0 = ( 1.0_dp / tan( theta0 ) - theta0 / sin(theta0)**2 ) * u_ij
              !dz0 = ( 1.0_dp / tan( theta0 ) - (theta0-pi/3) / sin(theta0)**2 ) * u_ij
              dr0 = ( cos( theta0 ) / (r*this%z0) - r0 / r ) * u_ij

              z0_pls_Iz = ( z0 + CPLX_IMAG*diff(3) ) * r0
              z0_min_Iz = ( z0 - CPLX_IMAG*diff(3) ) * r0
              x_pls_Iy = ( diff(1) + CPLX_IMAG*diff(2) ) * r0
              x_min_Iy = ( diff(1) - CPLX_IMAG*diff(2) ) * r0

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

              cutoff = coff(r,this%cutoff) * my_w(jn)
              dcutoff = -dcoff(r,this%cutoff)*u_ij * my_w(jn)

              !dz0_pls_Iz = - z0_pls_Iz * diff / r0**2
              !dz0_pls_Iz(3) = dz0_pls_Iz(3) + CPLX_IMAG/r0
              !dz0_min_Iz = - z0_min_Iz * diff / r0**2
              !dz0_min_Iz(3) = dz0_min_Iz(3) - CPLX_IMAG/r0
              !dx_pls_Iy = - x_pls_Iy * diff / r0**2
              !dx_pls_Iy(1) = dx_pls_Iy(1) + 1.0_dp/r0
              !dx_pls_Iy(2) = dx_pls_Iy(2) + CPLX_IMAG/r0
              !dx_min_Iy = - x_min_Iy * diff / r0**2
              !dx_min_Iy(1) = dx_min_Iy(1) + 1.0_dp/r0
              !dx_min_Iy(2) = dx_min_Iy(2) - CPLX_IMAG/r0

              Up(0:0,0:0) = CPLX_ONE
              dUp(:,0:0,0:0) = CPLX_ZERO

              this%U(0)%mm(:,0,0) = this%U(0)%mm(:,0,0) - dcutoff
              do j = 1, this%j_max
                 U(-j:j,-j:j) = CPLX_ZERO
                 dU(:,-j:j,-j:j) = CPLX_ZERO

                 do m1 = -j, j-2, 2
                    do m2 = -j, j, 2
                       if( (j-m2) /= 0 ) then
                          U(m2,m1) = U(m2,m1) + &
                          & sqrt( real(j-m2,dp)/real(j-m1,dp) ) * z0_pls_Iz * Up(m2+1,m1+1)
                          dU(:,m2,m1) = dU(:,m2,m1) + &
                          & sqrt( real(j-m2,dp)/real(j-m1,dp) ) * &
                          & ( dz0_pls_Iz * Up(m2+1,m1+1) + z0_pls_Iz * dUp(:,m2+1,m1+1) )
                       endif

                       if( (j+m2) /= 0 ) then
                          U(m2,m1) = U(m2,m1) - &
                          & CPLX_IMAG * sqrt( real(j+m2,dp)/real(j-m1,dp) ) * x_min_Iy * &
                          & Up(m2-1,m1+1)
                          dU(:,m2,m1) = dU(:,m2,m1) - &
                          & CPLX_IMAG * sqrt( real(j+m2,dp)/real(j-m1,dp) ) * &
                          & ( dx_min_Iy * Up(m2-1,m1+1) + x_min_Iy * dUp(:,m2-1,m1+1) )
                       endif
                    enddo
                 enddo
                          
                 m1 = j
                 do m2 = -j, j, 2
                    if( (j+m2) /= 0 ) then
                       U(m2,m1) = U(m2,m1) + &
                       & sqrt( real(j+m2,dp)/real(j+m1,dp) ) * z0_min_Iz * Up(m2-1,m1-1)
                       dU(:,m2,m1) = dU(:,m2,m1) + &
                       & sqrt( real(j+m2,dp)/real(j+m1,dp) ) * &
                       & ( dz0_min_Iz * Up(m2-1,m1-1) + z0_min_Iz * dUp(:,m2-1,m1-1) )
                    endif
                    if( (j-m2) /= 0 ) then
                       U(m2,m1) = U(m2,m1) - &
                       & CPLX_IMAG * sqrt( real(j-m2,dp)/real(j+m1,dp) ) * x_pls_Iy * &
                       & Up(m2+1,m1-1)
                       dU(:,m2,m1) = dU(:,m2,m1) - &
                       & CPLX_IMAG * sqrt( real(j-m2,dp)/real(j+m1,dp) ) * &
                       & ( dx_pls_Iy * Up(m2+1,m1-1) + x_pls_Iy * dUp(:,m2+1,m1-1) )
                    endif
                 enddo

                 Up(-j:j,-j:j) = U(-j:j,-j:j)
                 dUp(:,-j:j,-j:j) = dU(:,-j:j,-j:j)
                 this%U(j)%mm = this%U(j)%mm + dU(:,-j:j,-j:j) * cutoff
                 do m1 = -j, j, 2
                    do m2 = -j, j, 2
                       this%U(j)%mm(:,m2,m1) = this%U(j)%mm(:,m2,m1) &
                       & - U(m2,m1) * dcutoff
                    enddo
                 enddo
              enddo
           endif
        endif
        deallocate( U, Up, dU, dUp, my_w )

     endsubroutine grad_fourier_transform_so4

     subroutine fourier_transform_so4_old(this,at,i)
        type(fourier_so4), intent(inout) :: this
        type(atoms), intent(in) :: at
        integer, intent(in) :: i

        integer :: j, n, jn, m1, m2
        real(dp) :: theta0, theta, phi, cutoff, r0, gs
        real(dp), dimension(3) :: polar, diff

        if( .not. this%initialised ) call system_abort('fourier_transform_so4: fourier_so4 not initialised')
        do j = 0, this%j_max
           this%U(j)%mm(:,:) = CPLX_ZERO
!           do m1 = -j, j, 2
!              this%U(j)%mm(m1,m1) = CPLX_ONE
!           enddo
        enddo

        do n = 1, atoms_n_neighbours(at,i)
           jn = atoms_neighbour(at,i,n,diff=diff)
           polar = xyz2spherical( diff )
           !theta0 = atan(polar(1)/this%z0)
           !theta0 = atan( ((erf((polar(1)-2.4_dp)/0.6_dp)+1.0_dp)*0.5_dp)/this%z0)
           theta0 = polar(1)/this%z0
           theta = polar(2)
           phi = polar(3)
           !cutoff = 1.0_dp 
           cutoff = 1.0_dp / (exp((polar(1)-2.15_dp)/0.25_dp) + 1.0_dp )
           !cutoff = 1.0_dp / (exp((polar(1)-3.05_dp)/0.30_dp) + 1.0_dp )
           do j = 0, this%j_max
!r0 = j*0.6_dp + 2.0_dp
!gs = exp( - (polar(1)-r0)**2/(1.2_dp**2) )
              do m1 = -j, j, 2
                 do m2 = -j, j, 2
                    this%U(j)%mm(m2,m1) = this%U(j)%mm(m2,m1) + &
                    & wigner_big_U(j,m2,m1,2.0_dp*theta0,theta,phi,2) * &
                    & cutoff !* sqrt( 2.0_dp * j + 1.0_dp ) !* gs
                 enddo
              enddo
           enddo
        enddo
     endsubroutine fourier_transform_so4_old

!     subroutine grad_fourier_transform_so4_old(this,at,i,ij)
!        type(grad_fourier_so4), intent(inout) :: this
!        type(atoms), intent(in) :: at
!        integer, intent(in) :: i, ij
!
!        integer :: j, n, jn, m1, m2
!        real(dp) :: theta0, theta, phi, cutoff, gs, r0, dgs
!        real(dp), dimension(3) :: polar, diff, u_ij
!        real(dp), dimension(3,3) :: jacobian
!        complex(dp), dimension(3) :: dw
!        logical :: is_neighbour
!
!        if( .not. this%initialised ) call system_abort('fourier_transform_so4: fourier_so4 not initialised')
!        do j = 0, this%j_max
!           do m1 = -j, j
!              do m2 = -j, j
!                 this%U(j)%mm(m2,m1)%x = CPLX_ZERO
!              enddo
!           enddo
!        enddo
!
!        if(i==ij) then
!           do n = 1, atoms_n_neighbours(at,i)
!              jn = atoms_neighbour(at,i,n,diff=diff,cosines=u_ij)
!              polar = xyz2spherical( diff )
!              theta0 = atan(polar(1)/this%z0)
!              theta = polar(2)
!              phi = polar(3)
!              !cutoff = 1.0_dp / (exp((polar(1)-2.15_dp)/0.25_dp) + 1.0_dp )
!              cutoff = 1.0_dp / (exp((polar(1)-3.05_dp)/0.30_dp) + 1.0_dp )
!
!              jacobian(:,1) = 1.0_dp/(1.0_dp + (polar(1)/this%z0)**2) / this%z0 * u_ij
!              jacobian(:,2) = diff(3)*u_ij/(polar(1)**2*sqrt(1.0_dp-diff(3)**2/polar(1)**2))
!              jacobian(3,2) = jacobian(3,2) - 1.0_dp/(polar(1)*sqrt(1.0_dp-diff(3)**2/polar(1)**2))
!              if( ( diff(1) == 0.0_dp ) .and. ( diff(2) == 0.0_dp ) ) then
!                 jacobian(:,3) = 0.0_dp
!              else
!                 jacobian(:,3) = (/-diff(2),diff(1),0.0_dp/) / ( diff(1)**2 + diff(2)**2 )
!              endif
!
!              do j = 0, this%j_max
!!r0 = j*0.6_dp + 2.0_dp
!!gs = exp( - (polar(1)-r0)**2/(1.2_dp**2) )
!!dgs = -2.0_dp*(polar(1)-r0)*exp( - (polar(1)-r0)**2/(1.2_dp**2) )/(1.2_dp**2)
!                 do m1 = -j, j, 2
!                    do m2 = -j, j, 2
!                       dw = grad_wigner_big_U(j,m2,m1,4.0_dp*theta0,theta,phi,2)
!                       dw(1) = 4.0_dp * dw(1)
!                       this%U(j)%mm(m2,m1)%x = this%U(j)%mm(m2,m1)%x - sqrt( 2.0_dp * j + 1.0_dp ) * &
!                       & ( matmul(jacobian,dw) * cutoff - &
!                       & wigner_big_U(j,m2,m1,4.0_dp*theta0,theta,phi,2) * &
!                       & cutoff**2 * exp((polar(1)-3.05_dp)/0.30_dp) / 0.30_dp * u_ij )
!                       !this%U(j)%mm(m2,m1)%x = this%U(j)%mm(m2,m1)%x - sqrt( 2.0_dp * j + 1.0_dp ) * &
!                       !& ( matmul(jacobian,dw) * cutoff * gs - &
!                       !& wigner_big_U(j,m2,m1,4.0_dp*theta0,theta,phi,2) * ( gs *  &
!                       !& cutoff**2 * exp((polar(1)-3.05_dp)/0.30_dp) / 0.30_dp - &
!                       !& cutoff * dgs ) * u_ij )
!                    enddo
!                 enddo
!              enddo
!           enddo
!        else
!           is_neighbour = .false.
!           do n = 1, atoms_n_neighbours(at,i)
!              if( ij==atoms_neighbour(at,i,n) ) is_neighbour = .true.
!           enddo
!           if( is_neighbour ) then
!              diff = diff_min_image(at,i,ij)
!              polar = xyz2spherical( diff )
!              theta0 = atan(polar(1)/this%z0)
!              theta = polar(2)
!              phi = polar(3)
!              !cutoff = 1.0_dp / (exp((polar(1)-2.15_dp)/0.25_dp) + 1.0_dp )
!              cutoff = 1.0_dp / (exp((polar(1)-3.05_dp)/0.30_dp) + 1.0_dp )
!              u_ij = diff / polar(1)
!
!              jacobian(:,1) = 1.0_dp/(1.0_dp + (polar(1)/this%z0)**2) / this%z0 * u_ij
!              jacobian(:,2) = diff(3)*u_ij/(polar(1)**2*sqrt(1.0_dp-diff(3)**2/polar(1)**2))
!              jacobian(3,2) = jacobian(3,2) - 1.0_dp/(polar(1)*sqrt(1.0_dp-diff(3)**2/polar(1)**2))
!              if( ( diff(1) == 0.0_dp ) .and. ( diff(2) == 0.0_dp ) ) then
!                 jacobian(:,3) = 0.0_dp
!              else
!                 jacobian(:,3) = (/-diff(2),diff(1),0.0_dp/) / ( diff(1)**2 + diff(2)**2 )
!              endif
!
!              do j = 0, this%j_max
!!r0 = j*0.6_dp + 2.0_dp
!!gs = exp( - (polar(1)-r0)**2/(1.2_dp**2) )
!!dgs = -2.0_dp*(polar(1)-r0)*exp( - (polar(1)-r0)**2/(1.2_dp**2) )/(1.2_dp**2)
!                 do m1 = -j, j, 2
!                    do m2 = -j, j, 2
!                       dw = grad_wigner_big_U(j,m2,m1,4.0_dp*theta0,theta,phi,2)
!                       dw(1) = 4.0_dp * dw(1)
!                       this%U(j)%mm(m2,m1)%x = sqrt( 2.0_dp * j + 1.0_dp ) * &
!                       & ( matmul(jacobian,dw) * cutoff - &
!                       & wigner_big_U(j,m2,m1,4.0_dp*theta0,theta,phi,2) * &
!                       & cutoff**2 * exp((polar(1)-3.05_dp)/0.30_dp) / 0.30_dp * u_ij )
!                       !this%U(j)%mm(m2,m1)%x = sqrt( 2.0_dp * j + 1.0_dp ) * &
!                       !& ( matmul(jacobian,dw) * cutoff * gs - &
!                       !& wigner_big_U(j,m2,m1,4.0_dp*theta0,theta,phi,2) * &
!                       !& ( cutoff**2 * exp((polar(1)-3.05_dp)/0.30_dp) / 0.30_dp * gs - &
!                       !& cutoff * dgs ) * u_ij )
!                    enddo
!                 enddo
!              enddo
!           endif
!        endif
!
!     endsubroutine grad_fourier_transform_so4_old

     subroutine fourier_transform_so3(this, at, i)
       type(fourier_so3), intent(inout) :: this
       type(atoms), intent(in) :: at
       integer, intent(in) :: i

       integer :: l, m, f, n, jn
       real(dp) :: x_ij(3), R

       if (.not. this%initialised) call system_abort('fourier_transform_so3: fourier_so3 not initialised')
       if (at%cutoff < maxval(this%cutoff)) call system_abort('fourier_transform_so3: radial function cutoff greater than atoms connectivity cutoff')

       do l = 2, this%l_max, 2
          do f = 1, size(this%cutoff)
             this%Y_times_R(l / 2,f)%m(:) = CPLX_ZERO
          enddo
       enddo

       do f = 1, size(this%cutoff)
          do n = 1, atoms_n_neighbours(at, i)
             jn = atoms_neighbour(at, i, n, diff = x_ij, max_dist = this%cutoff(f))

             if (jn /= 0) then
                R = RadialFunction(this%cutoff(f), x_ij, this%cutoff_f(f), this%cutoff_r1(f))

                do l = 2, this%l_max, 2
                   do m = -l, l
                      this%Y_times_R(l / 2,f)%m(m) = this%Y_times_R(l / 2,f)%m(m) &
                                                   + (SphericalYCartesian(l, m, x_ij) * R)
                   enddo
                enddo
             endif
          enddo
       enddo

     endsubroutine fourier_transform_so3

     subroutine grad_fourier_transform_so3(this, at, i, ni)
       type(grad_fourier_so3), intent(inout) :: this
       type(atoms), intent(in) :: at
       integer, intent(in) :: i, ni

       integer :: l, m, f, n, jn
       real(dp) :: x_ij(3), R, dR(3)

       if (.not. this%initialised) call system_abort('grad_fourier_transform_so3: grad_fourier_so3 not initialised')
       if (at%cutoff < maxval(this%cutoff)) call system_abort('grad_fourier_transform_so3: radial function cutoff greater than atoms connectivity cutoff')

       do l = 2, this%l_max, 2
          do f = 1, size(this%cutoff)
             this%Y_times_R(l / 2,f)%m(:,:) = CPLX_ZERO
          enddo
       enddo

       if (ni == 0) then
          do f = 1, size(this%cutoff)
             do n = 1, atoms_n_neighbours(at, i)
                jn = atoms_neighbour(at, i, n, diff = x_ij, max_dist = this%cutoff(f))

                if (jn /= 0) then
                   R = RadialFunction(this%cutoff(f), x_ij, this%cutoff_f(f), this%cutoff_r1(f))
                   dR = GradRadialFunction(this%cutoff(f), x_ij, this%cutoff_f(f), this%cutoff_r1(f))

                   do l = 2, this%l_max, 2
                      do m = -l, l
                         this%Y_times_R(l / 2,f)%m(:,m) = this%Y_times_R(l / 2,f)%m(:,m) &
                                                        - (GradSphericalYCartesian(l, m, x_ij) * R) &
                                                        - (SphericalYCartesian(l, m, x_ij) * dR)
                      enddo
                   enddo
                endif
             enddo
          enddo
       else
          do f = 1, size(this%cutoff)
             jn = atoms_neighbour(at, i, ni, diff = x_ij, max_dist = this%cutoff(f))

             if (jn /= 0) then
                R = RadialFunction(this%cutoff(f), x_ij, this%cutoff_f(f), this%cutoff_r1(f))
                dR = GradRadialFunction(this%cutoff(f), x_ij, this%cutoff_f(f), this%cutoff_r1(f))

                do l = 2, this%l_max, 2
                   do m = -l, l
                      this%Y_times_R(l / 2,f)%m(:,m) = this%Y_times_R(l / 2,f)%m(:,m) &
                                                     + (GradSphericalYCartesian(l, m, x_ij) * R) &
                                                     + (SphericalYCartesian(l, m, x_ij) * dR)
                   enddo
                enddo
             endif
          enddo
       endif

     endsubroutine grad_fourier_transform_so3

     subroutine calc_bispectrum_so4(this,f)
        type(bispectrum_so4), intent(inout) :: this
        type(fourier_so4), intent(in) :: f

        integer :: j_max, j, j1, j2, m1, m2, m11, m12, m21, m22
        complex(dp) :: sub

        if( .not. this%initialised ) call initialise(this,f%j_max)

        if( .not. cg_initialised ) then
           call cg_initialise(f%j_max,2)
        elseif( f%j_max > cg_j_max ) then
           call cg_finalise()
           call cg_initialise(f%j_max,2)
        endif

        j_max = min(this%j_max,f%j_max)

        do j1 = 0, j_max
            do j2 = 0, j_max
               this%coeff(j2,j1)%mag(:) = CPLX_ZERO
            enddo
        enddo

        do j1 = 0, j_max
           j2 = j1
           !do j2 = 0, j_max
              do j = abs(j1-j2), min(j_max,j1+j2)
                 do m1 = -j, j, 2
                    do m2 = -j, j, 2
                       sub = CPLX_ZERO
                       do m11 = max(-j1,m1-j2), min(j1,m1+j2), 2
                          do m12 = max(-j1,m2-j2), min(j1,m2+j2), 2
                             sub = sub + cg_array(j1,m11,j2,m1-m11,j,m1) &
                             & * cg_array(j1,m12,j2,m2-m12,j,m2) &
                             & * f%U(j1)%mm(m11,m12) * f%U(j2)%mm(m1-m11,m2-m12)
                          enddo
                       enddo
                       this%coeff(j1,j2)%mag(j) = this%coeff(j1,j2)%mag(j) + &
                       & sub*conjg(f%U(j)%mm(m1,m2))
                    enddo
                 enddo
              enddo
           !enddo
        enddo

     endsubroutine calc_bispectrum_so4

     subroutine calc_grad_bispectrum_so4(this,f,df)
        type(grad_bispectrum_so4), intent(inout) :: this
        type(fourier_so4), intent(in) :: f
        type(grad_fourier_so4), intent(in) :: df

        integer :: j_max, j, j1, j2, m1, m2, m11, m12, m21, m22
        complex(dp) :: sub
        complex(dp), dimension(3) :: dsub
        real(dp) :: tmp_cg

        if( .not. this%initialised ) call initialise(this,f%j_max)

        if( .not. cg_initialised ) then
           call cg_initialise(f%j_max,2)
        elseif( f%j_max > cg_j_max ) then
           call cg_finalise()
           call cg_initialise(f%j_max,2)
        endif

        j_max = min(this%j_max,f%j_max,df%j_max)

        do j1 = 0, j_max
            do j2 = 0, j_max
               do j = abs(j1-j2), min(j_max,j1+j2)
                  this%coeff(j2,j1)%mag(j)%x = CPLX_ZERO
               enddo
            enddo
        enddo

        do j1 = 0, j_max
           j2 = j1
           !do j2 = 0, j_max
              do j = abs(j1-j2), min(j_max,j1+j2)
                 do m1 = -j, j, 2
                    do m2 = -j, j, 2
                       sub = CPLX_ZERO
                       dsub = CPLX_ZERO
                       do m11 = max(-j1,m1-j2), min(j1,m1+j2), 2
                          do m12 = max(-j1,m2-j2), min(j1,m2+j2), 2

                             tmp_cg =  cg_array(j1,m11,j2,m1-m11,j,m1) &
                             & * cg_array(j1,m12,j2,m2-m12,j,m2) 
                             sub = sub + tmp_cg &
                             & * f%U(j1)%mm(m11,m12) * f%U(j2)%mm(m1-m11,m2-m12)
                             dsub = dsub + tmp_cg &
                             & * ( df%U(j1)%mm(:,m11,m12) * f%U(j2)%mm(m1-m11,m2-m12) + &
                             & f%U(j1)%mm(m11,m12) * df%U(j2)%mm(:,m1-m11,m2-m12) )
                          enddo
                       enddo
                       this%coeff(j1,j2)%mag(j)%x = this%coeff(j1,j2)%mag(j)%x + &
                       & dsub*conjg(f%U(j)%mm(m1,m2)) + sub*conjg(df%U(j)%mm(:,m1,m2))
                    enddo
                 enddo
              enddo
           !enddo
        enddo

     endsubroutine calc_grad_bispectrum_so4

     subroutine calc_qw_so3(this, f)
       type(qw_so3), intent(inout) :: this
       type(fourier_so3), intent(in) :: f

       real(dp) :: m_normsq_sum
       complex(dp) :: wc
       integer :: l, n, m1, m2, m3

       do n = 1, this%f_n
          do l = 2, this%l_max, 2
             m_normsq_sum = sum(real(f%Y_times_R(l / 2,n)%m * conjg(f%Y_times_R(l / 2,n)%m), dp))

             if (this%do_q) then
                if (abs(m_normsq_sum) < QW_FP_ZERO) then
                   this%q(l / 2,n) = 0.0_dp
                else
                   this%q(l / 2,n) = sqrt(4.0_dp * PI * m_normsq_sum / ((2.0_dp * l) + 1.0_dp))
                endif
             endif

             if (this%do_w) then
                wc = CPLX_ZERO

                do m1 = -l, l
                   do m2 = -l, l
                      do m3 = -l, l
                         if ((m1 + m2 + m3) /= 0 ) cycle

                         wc = wc + (wigner3j(l, m1, l, m2, l, m3) * &
                                    f%Y_times_R(l / 2,n)%m(m1) * &
                                    f%Y_times_R(l / 2,n)%m(m2) * &
                                    f%Y_times_R(l / 2,n)%m(m3))
                      enddo
                   enddo
                enddo

                if (abs(m_normsq_sum) < QW_FP_ZERO) then
                   wc = CPLX_ZERO
                else
                   wc = wc / (m_normsq_sum**1.5_dp)
                endif

                this%w(l / 2,n) = real(wc, dp)
             endif
          enddo
       enddo

     endsubroutine calc_qw_so3

     subroutine calc_grad_qw_so3(this, f, df)
       type(grad_qw_so3), intent(inout) :: this
       type(fourier_so3), intent(in) :: f
       type(grad_fourier_so3), intent(in) :: df

       real(dp) :: m_normsq_sum, dm_normsq_sum(3)
       complex(dp) :: wc, dwc(3)
       integer :: l, n, k, m1, m2, m3

       do n = 1, this%f_n
          do l = 2, this%l_max, 2
             m_normsq_sum = sum(real(f%Y_times_R(l / 2,n)%m * conjg(f%Y_times_R(l / 2,n)%m), dp))

             do k = 1, 3
                dm_normsq_sum(k) = sum(real((f%Y_times_R(l / 2,n)%m * conjg(df%Y_times_R(l / 2,n)%m(k,:))) + &
                                           (df%Y_times_R(l / 2,n)%m(k,:) * conjg(f%Y_times_R(l / 2,n)%m)), dp))
             end do

             if (this%do_q) then
                if (abs(m_normsq_sum) < QW_FP_ZERO) then
                   this%q(l / 2,n)%x = 0.0_dp
                else
                   this%q(l / 2,n)%x = 2.0_dp * PI * dm_normsq_sum / sqrt(4.0_dp * PI * m_normsq_sum * ((2.0_dp * l) + 1.0_dp))
                endif
             endif

             if (this%do_w) then
                wc = CPLX_ZERO
                dwc = CPLX_ZERO

                do m1 = -l, l
                   do m2 = -l, l
                      do m3 = -l, l
                         if ((m1 + m2 + m3) /= 0 ) cycle

                         wc = wc + (wigner3j(l, m1, l, m2, l, m3) * &
                                    f%Y_times_R(l / 2,n)%m(m1) * &
                                    f%Y_times_R(l / 2,n)%m(m2) * &
                                    f%Y_times_R(l / 2,n)%m(m3))

                         dwc = dwc + (wigner3j(l, m1, l, m2, l, m3) &
                                   * ((df%Y_times_R(l / 2,n)%m(:,m1) * f%Y_times_R(l / 2,n)%m(m2) * f%Y_times_R(l / 2,n)%m(m3)) &
                                    + (f%Y_times_R(l / 2,n)%m(m1) * df%Y_times_R(l / 2,n)%m(:,m2) * f%Y_times_R(l / 2,n)%m(m3)) &
                                    + (f%Y_times_R(l / 2,n)%m(m1) * f%Y_times_R(l / 2,n)%m(m2) * df%Y_times_R(l / 2,n)%m(:,m3))))
                      enddo
                   enddo
                enddo

                if (abs(m_normsq_sum) < QW_FP_ZERO) then
                   wc = CPLX_ZERO
                   dwc = CPLX_ZERO
                else
                   wc = wc / (m_normsq_sum**1.5_dp)
                   dwc = (dwc / (m_normsq_sum**1.5_dp)) - (1.5_dp * wc * dm_normsq_sum / m_normsq_sum)
                endif

                this%w(l / 2,n)%x = real(dwc, dp)
             endif
          enddo
       enddo

     endsubroutine calc_grad_qw_so3

     subroutine bispectrum2vec_so4(this, vec)
        type(bispectrum_so4), intent(in) :: this
        real(dp), dimension(:), intent(out) :: vec

        integer :: j1, j2, j, i

        if( size( vec ) < j_max2d(this%j_max) ) &
        & call system_abort('bispectrum2vec_so4: vec too small')

        vec = 0.0_dp
        i = 1
        do j1 = 0, this%j_max
           j2 = j1
           !do j2 = 0, this%j_max
              do j = abs(j1-j2), min(this%j_max,j1+j2)
                 if( mod(j1+j2+j,2) == 1 ) cycle
                 vec(i) = real(this%coeff(j1,j2)%mag(j))
                 i = i + 1
              enddo
           !enddo
        enddo

     endsubroutine bispectrum2vec_so4

     subroutine print_bispectrum_so4(this)
        type(bispectrum_so4), intent(in) :: this
        integer :: j1, j2, j, i

        do j1 = 0, this%j_max
           do j2 = 0, this%j_max
              do j = abs(j1-j2), min(this%j_max,j1+j2)
                 call print(j1//"    "//j2//"    "//j//"    "//this%coeff(j1,j2)%mag(j))
              enddo
           enddo
        enddo
     endsubroutine print_bispectrum_so4

     subroutine bispectrum2vec_grad_so4(this, vec)
        type(grad_bispectrum_so4), intent(in) :: this
        real(dp), dimension(:,:), intent(out) :: vec

        integer :: j1, j2, j, i

        if( size( vec ) < 3*j_max2d(this%j_max) ) &
        & call system_abort('bispectrum2vec_so4: vec too small')

        vec = 0.0_dp
        i = 1
        do j1 = 0, this%j_max
           j2 = j1
           !do j2 = 0, this%j_max
              do j = abs(j1-j2), min(this%j_max,j1+j2)
                 if( mod(j1+j2+j,2) == 1 ) cycle
                 vec(i,:) = real(this%coeff(j1,j2)%mag(j)%x(:))
                 i = i + 1
              enddo
           !enddo
        enddo

     endsubroutine bispectrum2vec_grad_so4

     subroutine qw2vec_qw_so3(this, vec)
       type(qw_so3), intent(in) :: this
       real(dp), dimension(:), intent(out) :: vec

       integer :: d, i, l, f

       d = (this%l_max / 2) * this%f_n
       if (this%do_q .and. this%do_w) d = d * 2

       if (size(vec) < d) call system_abort('qw2vec_so3: vec too small')

       vec = 0.0_dp
       i = 1
       do l = 2, this%l_max, 2
          do f = 1, this%f_n
             if (this%do_q .and. this%do_w) then
                vec(i) = this%q(l / 2,f)
                vec((d / 2) + i) = this%w(l / 2,f)
             elseif (this%do_q) then
                vec(i) = this%q(l / 2,f)
             elseif (this%do_w) then
                vec(i) = this%w(l / 2,f)
             endif
             i = i + 1
          enddo
       enddo
           
     endsubroutine qw2vec_qw_so3

     subroutine qw2vec_grad_qw_so3(this, vec)
       type(grad_qw_so3), intent(in) :: this
       real(dp), dimension(:,:), intent(out) :: vec

       integer :: d, i, l, f

       d = (this%l_max / 2) * this%f_n
       if (this%do_q .and. this%do_w) d = d * 2

       if (size(vec) < (3 * d)) call system_abort('qw2vec_grad_so3: vec too small')

       vec = 0.0_dp
       i = 1
       do l = 2, this%l_max, 2
          do f = 1, this%f_n
             if (this%do_q .and. this%do_w) then
                vec(i,:) = this%q(l / 2,f)%x(:)
                vec((d / 2) + i,:) = this%w(l / 2,f)%x(:)
             elseif (this%do_q) then
                vec(i,:) = this%q(l / 2,f)%x(:)
             elseif (this%do_w) then
                vec(i,:) = this%w(l / 2,f)%x(:)
             endif
             i = i + 1
          enddo
       enddo

     endsubroutine qw2vec_grad_qw_so3

     function j_max2d(j_max)
        integer :: j_max2d
        integer, intent(in) :: j_max
        integer :: j1, j2, j

        j_max2d = 0

        do j1 = 0, j_max
           j2 = j1
           !do j2 = 0, j_max
              do j = abs(j1-j2), min(j_max,j1+j2)
                 if( mod(j1+j2+j,2) == 1 ) cycle
                 j_max2d = j_max2d + 1
              enddo
           !enddo
        enddo
     
     endfunction j_max2d

     function qw2d_qw_so3(this)
       type(qw_so3), intent(in) :: this
       integer :: qw2d_qw_so3

       qw2d_qw_so3 = qw2d_so3(qw = this)

     endfunction qw2d_qw_so3

     function qw2d_grad_qw_so3(this)
       type(grad_qw_so3), intent(in) :: this
       integer :: qw2d_grad_qw_so3

       qw2d_grad_qw_so3 = qw2d_so3(dqw = this)

     endfunction qw2d_grad_qw_so3

     function qw2d_so3(qw, dqw)
       type(qw_so3), intent(in), optional :: qw
       type(grad_qw_so3), intent(in), optional :: dqw
       integer :: qw2d_so3

       if (present(qw)) then
          qw2d_so3 = (qw%l_max / 2) * qw%f_n
          if (qw%do_q .and. qw%do_w) qw2d_so3 = qw2d_so3 * 2
       elseif (present(dqw)) then
          qw2d_so3 = (dqw%l_max / 2) * dqw%f_n
          if (dqw%do_q .and. dqw%do_w) qw2d_so3 = qw2d_so3 * 2
       endif

     endfunction qw2d_so3

     function test_qw_gradient(f, df, qw, dqw, at, x0)
       logical :: test_qw_gradient
       type(fourier_so3), intent(inout) :: f
       type(grad_fourier_so3), intent(inout) :: df
       type(qw_so3), intent(inout) :: qw
       type(grad_qw_so3), intent(inout) :: dqw
       type(atoms), intent(inout) :: at
       real(dp), intent(in) :: x0(3)
       real(dp) :: x00(3)
       real(dp), allocatable :: vec(:), dvec(:,:)
       integer :: i, j, n

       allocate (vec(qw2d(qw)), dvec(qw2d(qw),3))

       test_qw_gradient = .true.

       call add_atoms(at, x0, 1)
       call calc_connect(at)

       do n = 0, atoms_n_neighbours(at, at%N)
          call print("neighbour: "//n)
          call print("")

          if (n == 0) then
             j = at%N
          else
             j = atoms_neighbour(at, at%N, n)
          end if

          x00 = at%pos(:,j)

          do i = 1, qw2d(qw)
             call print("qw dimension: "//i)
             call print("")

             !test_qw_gradient = test_qw_gradient .and. test_gradient(x00, test_qw, test_dqw)
             call print("")
          end do
       end do

       call remove_atoms(at, at%N)
       call calc_connect(at)

       if (test_qw_gradient) then
          call print("TEST_QW_GRADIENT: PASS")
       else
          call print("TEST_QW_GRADIENT: FAIL")
       end if

       contains

       function test_qw(x, data)
         real(dp) :: test_qw
         real(dp) :: x(:)
         character, optional :: data(:)
         real(dp) :: x_old(3)

         x_old = at%pos(:,j)
         at%pos(:,j) = x
         call calc_connect(at)
         call fourier_transform(f, at, at%N)
         call calc_qw(qw, f)
         call qw2vec(qw, vec)

         test_qw = vec(i)

         at%pos(:,j) = x_old
         call calc_connect(at)

       endfunction test_qw

       function test_dqw(x, data)
         real(dp) :: test_dqw(3)
         real(dp) :: x(:)
         character, optional :: data(:)
         real(dp) :: x_old(3)

         x_old = at%pos(:,j)
         at%pos(:,j) = x
         call calc_connect(at)
         call fourier_transform(f, at, at%N)
         call fourier_transform(df, at, at%N, n)
         call calc_qw(dqw, f, df)
         call qw2vec(dqw, dvec)

         test_dqw = dvec(i,:)

         at%pos(:,j) = x_old
         call calc_connect(at)

       endfunction test_dqw

     endfunction test_qw_gradient

     subroutine get_weights(this,w)

        type(atoms), intent(in) :: this
        real(dp), dimension(:), intent(out) :: w

        integer, dimension(:), allocatable :: Z
        real(dp), dimension(:), allocatable :: w_type
        integer :: i, Z_prev, n_type

        if( size(w) < maxval(this%Z) ) call system_abort('get_weights: size of w is smaller than &
        & the largest atomic number '// maxval(this%Z) //' in atoms object')

        allocate(Z(this%N)) 
        ! Make copy of atomic numbers
        Z = this%Z
        ! Sort them
        call sort_array(Z)

        ! Determine how many different atomtypes are present
        Z_prev = 0
        n_type = 0
        do i = 1, this%N
           if( Z_prev /= Z(i) ) then
              n_type = n_type + 1
              Z_prev = Z(i)
           endif
        enddo

        ! Each atomtype is assigned with a weight. These are determined such
        ! that their cube falls between 0.5 and 1.0
        allocate(w_type(n_type))
        !w_type = ( 0.5_dp * (/ (i, i=n_type,1,-1) /) / real(n_type,dp) + 0.5_dp )**(1.0_dp/3.0_dp)
        !w_type = ( 0.5_dp * (/ (i, i=n_type,1,-1) /) / real(n_type,dp) + 0.5_dp ) 
        if( n_type == 1 ) then
           w_type = 1.0_dp
        else
           w_type = ( 0.5_dp*( (/ (i, i=n_type,1,-1) /) -1 ) / real(n_type-1,dp) + 0.5_dp )
        endif
           
        ! Fill array w with weight for each atom. This is done via the index
        ! array from the sorting.
        Z_prev = 0
        n_type = 0
        w = 0.0_dp
        do i = 1, this%N
           if( Z_prev /= Z(i) ) then
              n_type = n_type + 1
              Z_prev = Z(i)
              w(Z(i)) = w_type(n_type)
           endif
        enddo
        deallocate(Z,w_type)

     endsubroutine get_weights
     
     subroutine initialise_fourier_periodic(this,k_max1,k_max2,k_max3,at)

        type(per), intent(inout) :: this
        integer, intent(in) :: k_max1
        integer, intent(in), optional :: k_max2, k_max3
        type(atoms), intent(in), optional :: at

        integer, dimension(3) :: k_max

        if( this%initialised ) call finalise(this)

        k_max = (/k_max1,optional_default(k_max1,k_max2),optional_default(k_max1,k_max3)/)

        allocate( this%f(-k_max(1):k_max(1),-k_max(2):k_max(2),-k_max(3):k_max(3)) )
        this%f = CPLX_ZERO

        this%k = k_max

        if( present(at) ) then
           allocate(this%w(maxval(at%Z)))
           call get_weights(at,this%w)
        else
           allocate(this%w(size(ElementName)))
           this%w = 1.0_dp
        endif

        this%initialised = .true.

     endsubroutine initialise_fourier_periodic

     subroutine finalise_fourier_periodic(this)

        type(per), intent(inout) :: this

        if( .not. this%initialised ) return

        deallocate( this%f, this%w )

        this%k = 0
        this%initialised = .false.

     endsubroutine finalise_fourier_periodic

     subroutine fourier_transform_periodic(at,f_hat)
        type(Atoms), intent(inout) :: at
        type(per), intent(inout) :: f_hat

        integer :: k1, k2, k3, i
        complex(dp) :: TWO_PI_I
        real(dp), dimension(:,:), allocatable :: frac

        real(dp), dimension(3) :: k

        TWO_PI_I = 2.0_dp * PI * CPLX_IMAG

        f_hat%f = CPLX_ZERO        

        if( size(f_hat%w) < maxval(at%Z) ) call system_abort('size of weights &
        & array is '// size(f_hat%w) //', which is smaller than the maximumal &
        & atomic number '// maxval(at%Z) //' in atoms object')

        allocate(frac(3,at%N))
        do i = 1, at%N
           frac(:,i) = matmul(at%g,at%pos(:,i))
        enddo
        
        do i = 1, at%N
           do k3 = -f_hat%k(3), f_hat%k(3)
              do k2 = -f_hat%k(2), f_hat%k(2)
                 do k1 = -f_hat%k(1), f_hat%k(1)
                    k = k1 * at%g(1,:) + k2 * at%g(2,:) + k3 * at%g(3,:)
                    f_hat%f(k1,k2,k3) = f_hat%f(k1,k2,k3) + exp(-0.5_dp * normsq(k) * f_hat%sig**2) * &
                    & exp( TWO_PI_I * dot_product( (/k1,k2,k3/), frac(:,i) ) ) * f_hat%w(at%Z(i))
                    !f_hat%f(k1,k2,k3) = f_hat%f(k1,k2,k3) + &
                    !& exp( TWO_PI_I * dot_product( (/k1,k2,k3/), frac(:,i) ) ) * f_hat%w(at%Z(i))
                 enddo
              enddo
           enddo
        enddo

        deallocate(frac)
        
     endsubroutine fourier_transform_periodic

     function suggest_resolution(at)
        type(Atoms), intent(in) :: at
        real(dp) :: suggest_resolution

        integer :: i, j
        real(dp) :: distance

        suggest_resolution = huge(1.0_dp)
        do i = 1, at%N - 1
           do j = i+1, at%N
              distance = distance_min_image(at,i,j)
              if( distance > 0.0_dp ) suggest_resolution=min(suggest_resolution,distance)
           enddo
        enddo
              
     endfunction suggest_resolution

     function suggest_kmax(at,resolution)

        type(Atoms), intent(in) :: at
        real(dp), intent(in), optional :: resolution
        integer, dimension(3) :: suggest_kmax
        real(dp) :: my_resolution
        integer :: i        

        if( present(resolution) ) then
           my_resolution = resolution
        else
           my_resolution = suggest_resolution(at)
        endif

        do i = 1, 3
           suggest_kmax(i) = ceiling(norm(at%lattice(:,i))/my_resolution)
        enddo

     endfunction suggest_kmax

     subroutine bispectrum_periodic(f_hat,bis)
       
        type(per), intent(inout) :: f_hat, bis
        integer :: k1, k2, k3

        if( .not. f_hat%initialised ) call system_abort('bispectrum_periodic: f_hat not initialised')

        if( .not. bis%initialised .or. any ( 2*bis%k > f_hat%k ) ) then
           call finalise(bis)
           call initialise(bis,f_hat%k(1)/2,f_hat%k(2)/2,f_hat%k(3)/2)
        endif

        do k3 = -bis%k(3), bis%k(3)
           do k2 = -bis%k(2), bis%k(2)
              do k1 = -bis%k(1), bis%k(1)
                 bis%f(k1,k2,k3) = f_hat%f(k1,k2,k3)**2 * conjg(f_hat%f(2*k1,2*k2,2*k3))
              enddo
           enddo
        enddo

     endsubroutine bispectrum_periodic

     subroutine bispectrum2vec_periodic(bis,vec,p)

        type(per), intent(in) :: bis
        complex(dp), dimension(:), intent(out) :: vec
        integer, dimension(3), intent(in), optional :: p
        
        integer :: k1, k2, k3, i
        integer, dimension(3) :: k
        integer, dimension(3), target :: kp
        integer, pointer :: kp1, kp2, kp3

        if( .not. bis%initialised ) call system_abort('bispectrum2vec_periodic: bis not initialised' )
        if( size(vec) < kmax2d(bis) ) call system_abort('bispectrum2vec_periodic: vec too small')
        if( present(p) .and. .not. ( (bis%k(1)==bis%k(2)) .and. (bis%k(1)==bis%k(3)) ) ) &
        & call system_abort('bispectrum2vec_periodic: permutation present, so k(1) = k(2) = k(3)')
        
        vec = CPLX_ZERO

        i = 0
        kp1 => kp(1)
        kp2 => kp(2)
        kp3 => kp(3)

        do k3 = -bis%k(3), bis%k(3)
           do k2 = -bis%k(2), bis%k(2)
              do k1 = -bis%k(1), bis%k(1)
                 i = i + 1
                 k = (/k1,k2,k3/)
                 if(present(p)) then
                    kp = k(p)
                 else
                    kp = k
                 endif
                 vec(i) = bis%f(kp1,kp2,kp3)
              enddo
           enddo
        enddo

     endsubroutine bispectrum2vec_periodic
     
     function kmax2d(this)
        type(per), intent(in) :: this
        integer :: kmax2d

        if( .not. this%initialised ) call system_abort('kmax2d: bispectrum not initialised' )

        kmax2d = product(2*this%k+1)
        
     endfunction kmax2d

     subroutine reduce_lattice(matrix)
        real(dp), dimension(3,3), intent(inout) :: matrix
        real(dp), dimension(3,3) :: matrix2

        real(dp) :: mn, mn1, mn2
        integer :: i, j, i_min, j_min
        logical :: o_min
        

        do
           mn = matrix_norm(matrix)
           mn1 = mn
           do i = 1, 3
              do j = 1, 3
                 if( i==j) cycle

                 matrix2 = matrix
                 matrix2(:,i) = matrix2(:,i) - matrix2(:,j)
                 mn2 = matrix_norm(matrix2)
                 if( mn2 < mn1 ) then
                    mn1 = mn2
                    i_min = i
                    j_min = j
                    o_min = .true.
                 endif
                 matrix2 = matrix
                 matrix2(:,i) = matrix2(:,i) + matrix2(:,j)
                 mn2 = matrix_norm(matrix2)
                 if( mn2 < mn1 ) then
                    mn1 = mn2
                    i_min = i
                    j_min = j
                    o_min = .false.
                 endif
              enddo
           enddo
           if( mn1 == mn ) exit
           if( o_min ) then
              matrix(:,i_min) = matrix(:,i_min) - matrix(:,j_min)
           else
              matrix(:,i_min) = matrix(:,i_min) + matrix(:,j_min)
           endif
        enddo


     endsubroutine reduce_lattice

     function matrix_norm(matrix)
        real(dp), dimension(:,:), intent(in) :: matrix
        real(dp) :: matrix_norm

        matrix_norm = sqrt( sum( matrix**2 ) )
     
     endfunction matrix_norm

     
     subroutine order_lattice_vectors(at)

        type(atoms), intent(inout) :: at

        real(dp), dimension(3)       :: bQb, bOb, bHb, com
        real(dp), dimension(3,3)     :: Q, new_lattice
        real(dp), dimension(3,3,3)   :: O
        real(dp), dimension(3,3,3,3) :: H

        integer :: i, a, b, c, d

        call add_property(at,'mass',ElementMass(at%Z))

        com = centre_of_mass(at)
        do i = 1, at%N
           at%pos(:,i) = at%pos(:,i) - com
        enddo

        Q = 0.0_dp
        O = 0.0_dp
        H = 0.0_dp

        do i = 1, at%N
           do a = 1, 3
              do b = 1, 3
                 Q(b,a) = Q(b,a) + at%pos(a,i) * at%pos(b,i) * at%mass(i)
                 do c = 1, 3
                    O(c,b,a) = O(c,b,a) + at%pos(a,i) * at%pos(b,i) * at%pos(c,i) * at%mass(i)
                    do d = 1, 3
                       H(d,c,b,a) = H(d,c,b,a) + at%pos(a,i) * at%pos(b,i) * at%pos(c,i) * at%pos(d,i) * at%mass(i)
                    enddo
                 enddo
              enddo
           enddo
        enddo

        bQb = 0.0_dp
        bOb = 0.0_dp
        bHb = 0.0_dp

        do i = 1, 3
           do a = 1, 3
              do b = 1, 3
                 bQb(i) = bQb(i) + at%g(i,a) * at%g(i,b) * Q(b,a)
                 do c = 1, 3
                    bOb(i) = bOb(i) + at%g(i,a) * at%g(i,b) * at%g(i,c) * O(c,b,a)
                    do d = 1, 3
                       bHb(i) = bHb(i) + at%g(i,a) * at%g(i,b) * at%g(i,c) * at%g(i,d) * H(d,c,b,a)
                    enddo
                 enddo
              enddo
           enddo
        enddo

        call print(bQb)
        call print(bOb)
        call print(bHb)
     
     endsubroutine order_lattice_vectors

     !#################################################################################
     !#
     !% Calculates Q and W (central symmetry parameters) for each atom in an atoms object.
     !#
     !#################################################################################

     subroutine calc_qw_at(at,l,q_global,w_global)
        type(Atoms), intent(inout) :: at     !% atoms object, q and w added as properties: q//l and w//l
        integer, intent(in), optional :: l   !% which L spherical harmonics to use
        real(dp), optional :: q_global, w_global

        real(dp), dimension(:), pointer :: at_q, at_w

        integer :: my_l
        real(dp), dimension(3) :: u_ij, polar
        real(dp), dimension(:), allocatable :: q

        integer :: i, j, n, m, m1, m2, m3, n_bond
        complex(dp) :: qm_bond
        complex(dp), dimension(:), allocatable :: w, qm_global
        complex(dp), dimension(:,:), allocatable :: qm

        my_l = optional_default(4,l)

        if( mod(my_l,2) /=0 ) call system_abort('Wisdom says l should be even number, let it be so.')

        call add_property(at,'q'//my_l,0.0_dp)
        call add_property(at,'w'//my_l,0.0_dp)
        if( .not. assign_pointer(at, 'q'//my_l, at_q) ) &
        & call system_abort('at_qw: could not assign pointer to atoms object')
        if( .not. assign_pointer(at, 'w'//my_l, at_w) ) &
        & call system_abort('at_qw: could not assign pointer to atoms object')

        allocate(qm(-my_l:my_l,at%N), q(at%N), w(at%N), qm_global(-my_l:my_l) )

        qm = CPLX_ZERO
        w = CPLX_ZERO

        n_bond = 0
        qm_global = CPLX_ZERO
        do i = 1, at%N
           do n = 1, atoms_n_neighbours(at,i)
              j = atoms_neighbour(at,i,n,cosines=u_ij)
              polar = xyz2spherical(u_ij)
              n_bond = n_bond + 1
              do m = -my_l, my_l
                 qm_bond = SphericalY(my_l,m,polar(2),polar(3))
                 qm(m,i) = qm(m,i) + qm_bond/atoms_n_neighbours(at,i)
                 qm_global(m) = qm_global(m) + qm_bond
              enddo
           enddo
        enddo
        n_bond = max(n_bond,1)
        qm_global = qm_global / real(n_bond, dp)

        q = sum( real( qm*conjg(qm) ), dim=1 )
        if ( present(q_global) ) q_global = sqrt( PI * 4.0_dp * sum( real( qm_global*conjg(qm_global) ) ) / (2.0_dp * my_l + 1.0_dp ) )

        do i = 1, at%N
           do m1 = -my_l, my_l
              do m2 = -my_l, my_l
                 do m3 = -my_l, my_l
                    if( m1+m2+m3 /= 0 ) cycle
                    w(i) = w(i) + wigner3j(my_l,m1,my_l,m2,my_l,m3) * qm(m1,i) * qm(m2,i) * qm(m3,i)
                 enddo
              enddo
           enddo
        enddo

        if ( present(w_global) ) then
           w_global = 0.0_dp
           do m1 = -my_l, my_l
              do m2 = -my_l, my_l
                 do m3 = -my_l, my_l
                    if( m1+m2+m3 /= 0 ) cycle
                    w_global = w_global + real( wigner3j(my_l,m1,my_l,m2,my_l,m3) * qm_global(m1) * qm_global(m2) * qm_global(m3) ) 
                 enddo
              enddo
           enddo
           w_global = w_global / ( sum( real( qm_global*conjg(qm_global) ) )**(3.0_dp/2.0_dp) )
        endif

        w = w / (q**(3.0_dp/2.0_dp))
        q = sqrt(4.0_dp*PI/(2.0_dp * my_l + 1.0_dp)*q)
        at_q = q
        at_w = real(w)

        deallocate(qm, q, w, qm_global)

        at_q=>null()
        at_w=>null()

     endsubroutine calc_qw_at
       
     subroutine initialise_symm(this) !,n_radial,n_angular)
        type(symm), intent(inout) :: this
        !integer, intent(in)              :: n_radial, n_angular

        if( this%initialised ) call finalise(this)

        !this%n_radial = n_radial
        !this%n_angular = n_angular

        allocate( this%radial(this%n_radial), this%angular(this%n_angular) )

        this%initialised = .true.

     endsubroutine initialise_symm

     subroutine finalise_symm(this)
        type(symm), intent(inout) :: this

        if( .not. this%initialised ) return
        this%n_radial = 0
        this%n_angular = 0

        deallocate( this%radial, this%angular)

        this%initialised = .false.

     endsubroutine finalise_symm

     subroutine atom2symm(at,coeff,vec,ii)

        type(atoms), intent(in)               :: at
        type(symm), intent(in)                :: coeff
        real(dp), dimension(:,:), intent(out) :: vec
        integer, intent(in), optional         :: ii

        integer :: a, n, m, i, j, k, start_i, end_i
        real(dp) :: r_ij, r_ik, r_jk, f_ij, f_ik, f_jk, cos_ijk
        real(dp), dimension(3) :: d_ij, d_ik, d_jk
        integer, dimension(3) :: shift_ij, shift_ik

        if( (coeff%n_radial + coeff%n_angular) > size(vec,1) ) &
        & call system_abort('atom2symm: array sizes do not conform')

        start_i = 1
        end_i   = at%N
        if(present(ii)) then
           start_i = ii
           end_i   = ii
        endif

        vec = 0.0_dp

        do i = start_i, end_i
           do n = 1, atoms_n_neighbours(at,i)
              j = atoms_neighbour(at,i,n,distance=r_ij,diff=d_ij,shift=shift_ij)
              if( r_ij > coeff%cutoff ) cycle

              f_ij = coff(r_ij,coeff%cutoff)
              do a = 1, coeff%n_radial
                 vec(a,i) = vec(a,i) + exp( - coeff%radial(a)%eta * (r_ij - coeff%radial(a)%rs)**2 ) * f_ij
              enddo

              do m = 1, atoms_n_neighbours(at,i)
                 k = atoms_neighbour(at,i,m,distance=r_ik,diff=d_ik,shift=shift_ik)
                 if( r_ik > coeff%cutoff ) cycle

                 d_jk = diff(at,j,k,shift_ik-shift_ij)
                 r_jk = norm(d_jk)
                 if( r_jk .feq. 0.0_dp ) cycle

                 cos_ijk = dot_product(d_ij,d_ik)/(r_ij*r_ik) !cosine(at,i,j,k)

                 f_ik = coff(r_ik,coeff%cutoff)
                 f_jk = coff(r_jk,coeff%cutoff)

                 do a = 1, coeff%n_angular
                    vec(a+coeff%n_radial,i) = vec(a+coeff%n_radial,i) + &
                    (1.0_dp + coeff%angular(a)%lambda*cos_ijk)**coeff%angular(a)%zeta * &
                    & exp(-coeff%angular(a)%eta*(r_ij**2+r_ik**2+r_jk**2))*f_ij*f_ik*f_jk
                 enddo
              enddo
           enddo
        enddo
              
     endsubroutine atom2symm

     subroutine symm_jacobian(at,coeff,jacobian)

        type(atoms), intent(in)  :: at
        type(symm), intent(in)   :: coeff
        real(dp), dimension(:,:,:) :: jacobian

        integer :: n_basis, a, i, j, k, n, m, alpha
        real(dp) :: r_ij, r_ik, r_jk, f_ij, f_ik, f_jk, cos_ijk, rad, ang, dang
        real(dp), dimension(3) :: u_ij, u_ik, u_jk, df_ij, df_ik, df_jk, d_ij, d_ik, d_jk, dcosijk_ij, dcosijk_ik, tmp, tmp_i, tmp_j, tmp_k
        integer, dimension(3) :: shift_ij, shift_ik

        jacobian = 0.0_dp
        n_basis = coeff%n_radial + coeff%n_angular

        do i = 1, at%N

           do n = 1, atoms_n_neighbours(at,i)
              j = atoms_neighbour(at,i,n,distance=r_ij,cosines=u_ij,diff=d_ij,shift=shift_ij)
              if( r_ij > coeff%cutoff ) cycle

              f_ij  = coff(r_ij,coeff%cutoff)
              df_ij = dcoff(r_ij,coeff%cutoff) * u_ij

              do a = 1, coeff%n_radial
                 tmp = exp( -coeff%radial(a)%eta * (r_ij-coeff%radial(a)%rs)**2) * & 
                        & (-2.0_dp * coeff%radial(a)%eta * (r_ij - coeff%radial(a)%rs) * f_ij * u_ij + df_ij )

                 alpha = a

                 jacobian(alpha,1:3,i) = jacobian(alpha,1:3,i) - tmp
                 jacobian(alpha,3*n+1:3*(n+1),i) = jacobian(alpha,3*n+1:3*(n+1),i) + tmp
              enddo

              do m = 1, atoms_n_neighbours(at,i)
                 k = atoms_neighbour(at,i,m,distance=r_ik,cosines=u_ik,diff=d_ik,shift=shift_ik)
                 if( r_ik > coeff%cutoff ) cycle
                 
                 d_jk = diff(at,j,k,shift_ik-shift_ij)
                 r_jk = norm(d_jk)
                 u_jk = d_jk / r_jk
                 if( r_jk .feq. 0.0_dp ) cycle

                 cos_ijk = dot_product(d_ij,d_ik)/(r_ij*r_ik)

                 f_ik = coff(r_ik,coeff%cutoff)
                 f_jk = coff(r_jk,coeff%cutoff)

                 df_ik = dcoff(r_ik,coeff%cutoff) * u_ik
                 df_jk = dcoff(r_jk,coeff%cutoff) * u_jk
                 dcosijk_ij = ( u_ik - cos_ijk * u_ij ) / r_ij
                 dcosijk_ik = ( u_ij - cos_ijk * u_ik ) / r_ik

                 do a = 1, coeff%n_angular

                   rad = exp(-coeff%angular(a)%eta*(r_ij**2+r_ik**2+r_jk**2)) 
                   ang = (1.0_dp + coeff%angular(a)%lambda*cos_ijk)**coeff%angular(a)%zeta
                   dang = coeff%angular(a)%zeta*coeff%angular(a)%lambda*(1.0_dp + coeff%angular(a)%lambda*cos_ijk)**(coeff%angular(a)%zeta-1.0_dp)
                   
                   tmp_i = rad * ( - dang * (dcosijk_ij+dcosijk_ik) * f_ij*f_ik*f_jk + &
                   & ang * 2.0_dp * coeff%angular(a)%eta * ( d_ij + d_ik ) * f_ij*f_ik*f_jk - &
                   & ang * ( df_ij*f_ik*f_jk + f_ij*df_ik*f_jk ) )

                   tmp_j = rad * ( dang*dcosijk_ij * f_ij*f_ik*f_jk + &
                   & ang * 2.0_dp * coeff%angular(a)%eta * ( - d_ij + d_jk ) * f_ij*f_ik*f_jk + &
                   & ang * ( df_ij*f_ik*f_jk - f_ij*f_ik*df_jk ) )

                   tmp_k = rad * ( dang*dcosijk_ik * f_ij*f_ik*f_jk + &
                   & ang * 2.0_dp * coeff%angular(a)%eta * ( - d_jk - d_ik ) * f_ij*f_ik*f_jk + &
                   & ang * ( f_ij*df_ik*f_jk + f_ij*f_ik*df_jk ) )

                   alpha = coeff%n_radial + a

                   jacobian(alpha,1:3,i) = jacobian(alpha,1:3,i) + tmp_i
                   jacobian(alpha,3*n+1:3*(n+1),i) = jacobian(alpha,3*n+1:3*(n+1),i) + tmp_j
                   jacobian(alpha,3*m+1:3*(m+1),i) = jacobian(alpha,3*m+1:3*(m+1),i) + tmp_k
                 enddo
              enddo
           enddo
        enddo

     endsubroutine symm_jacobian

     !#################################################################################
     !#
     !% Initialises an atom-centered fourier coefficient object
     !#
     !#################################################################################

     subroutine initialise_fourier_atom(f_hat,N_max,L_max,r_min,r_max,sigma, &
                & cutoff_left,cutoff_right,cutoff_sigma)

        type(fourier_coefficient_atom), intent(inout) :: f_hat   !% fourier coefficient object
        integer, intent(in)  :: L_max                       !% angular basis functions
        integer, intent(in)  :: N_max                       !% radial basis functions
        real(dp), intent(in) :: r_min                       !% centre of first radial basis function
        real(dp), intent(in) :: r_max                       !% centre of last radial basis function
        real(dp), intent(in) :: sigma                       !% width of radial basis function
        real(dp), intent(in), optional :: cutoff_left       !% centre of the 'left' cutoff fermi function
        real(dp), intent(in), optional :: cutoff_right      !% centre of the 'right' cutoff fermi function
        real(dp), intent(in), optional :: cutoff_sigma      !% width of the cutoff fermi functions
    
        integer  :: n, l
        real(dp) :: dr
    
        if( f_hat%initialised ) then
           call system_abort('initialise_fourier_array: f_hat already initialised')
        endif
    
        allocate( f_hat%coeff%rad(N_max), f_hat%r0(N_max) )
    
        do n = 1, N_max
           allocate( f_hat%coeff%rad(n)%ang(0:L_max) )
           do l = 0, L_max
              allocate( f_hat%coeff%rad(n)%ang(l)%mag( -L_max:L_max ) )
              f_hat%coeff%rad(n)%ang(l)%mag(:) = cplx_zero
           enddo
        enddo

        dr = (r_max-r_min)/(N_max-1)
        f_hat%r0(1) = r_min
        do n = 2, N_max
           f_hat%r0(n) = f_hat%r0(n-1)+dr
        enddo
        f_hat%sigma = sigma
        f_hat%cutoff_left  = optional_default(r_min,cutoff_left)
        f_hat%cutoff_right = optional_default(r_max,cutoff_right)
        f_hat%cutoff_sigma = optional_default(1.0_dp/sigma,cutoff_sigma)
    
        f_hat%N_max = N_max
        f_hat%L_max = L_max
        f_hat%initialised = .true.
    
     endsubroutine initialise_fourier_atom
    
     !#################################################################################
     !#
     !% Finalises an atom-centered fourier coefficient object
     !#
     !#################################################################################

     subroutine finalise_fourier_atom(f_hat)
         
        type(fourier_coefficient_atom), intent(inout) :: f_hat
        integer :: n, l
    
        if( .not. f_hat%initialised ) return
    
        do n = 1, f_hat%N_max
           do l = 0, f_hat%L_max
              deallocate( f_hat%coeff%rad(n)%ang(l)%mag )
           enddo
           deallocate( f_hat%coeff%rad(n)%ang )
        enddo
    
        deallocate( f_hat%coeff%rad, f_hat%r0 )
    
        f_hat%sigma = 0.0_dp
        f_hat%cutoff_left  = 0.0_dp
        f_hat%cutoff_right = 0.0_dp
        f_hat%cutoff_sigma = 0.0_dp

        f_hat%N_max = 0
        f_hat%L_max = 0
        f_hat%initialised = .false.
    
     endsubroutine finalise_fourier_atom
    
     !#################################################################################
     !#
     !% Initialises an atom-centered bispectrum object
     !#
     !#################################################################################

     subroutine initialise_bispectrum_atom(bispectrum,N_max,L_max) 
    
        type(bispectrum_coefficient_atom), intent(inout) :: bispectrum !% bispectrum
        integer, intent(in) :: N_max                                   !% number of radial basis functions
        integer, intent(in) :: L_max                                   !% number of angular basis functions
    
        integer :: l1, l2, n1, n2
    
        if( bispectrum%initialised ) call finalise(bispectrum)

!        if( mod(N_max,2) == 0 ) call system_abort('initialise_bispectrum: N_max must be odd number')
    
        !allocate( bispectrum%coeff%rad( (N_max-1)/2 ) )
        allocate( bispectrum%coeff%rad(N_max,N_max) )
    
        do n1 = 1, N_max
           do n2 = 1, N_max
               allocate( bispectrum%coeff%rad(n2,n1)%ang( 0:L_max, 0:L_max ) )
               do l1 = 0,L_max
                  do l2 = 0,L_max
                     allocate( bispectrum%coeff%rad(n2,n1)%ang(l2,l1)%mag( abs(l1-l2):min(l1+l2,L_max) ) )
                     bispectrum%coeff%rad(n2,n1)%ang(l2,l1)%mag(:) = cplx_zero
                  enddo
              enddo
           enddo
        enddo

        bispectrum%N_max = N_max
        bispectrum%L_max = L_max
        bispectrum%initialised = .true.
    
     endsubroutine initialise_bispectrum_atom
    
     !#################################################################################
     !#
     !% Finalises an atom-centered bispectrum object
     !#
     !#################################################################################

     subroutine finalise_bispectrum_atom(bispectrum)
         
        type(bispectrum_coefficient_atom), intent(inout) :: bispectrum !% bispectrum object
        integer :: l1, l2, n1, n2
    
        if( .not. bispectrum%initialised ) return
    
        do n1 = 1, bispectrum%N_max
           do n2 = 1, bispectrum%N_max
              do l1 = 0,bispectrum%L_max
                 do l2 = 0,bispectrum%L_max
                    deallocate( bispectrum%coeff%rad(n2,n1)%ang(l2,l1)%mag )
                 enddo
              enddo
              deallocate( bispectrum%coeff%rad(n2,n1)%ang )
           enddo
        enddo
    
        deallocate( bispectrum%coeff%rad )

        bispectrum%N_max = 0
        bispectrum%L_max = 0
        bispectrum%initialised = .false.
    
     endsubroutine finalise_bispectrum_atom

     !#################################################################################
     !#
     !% Calculates an atom-centered fourier spectrum
     !#
     !#################################################################################

     subroutine fourier_transform_vec(at,vec,f_hat,Z)

       type(atoms), intent(in)                       :: at       !% atoms object
       real(dp), dimension(3), intent(in)            :: vec      !% centre
       type(fourier_coefficient_atom), intent(inout) :: f_hat    !% fourier coefficents
       integer, intent(in), optional                 :: Z        !% which type of atoms to include

       real(dp), dimension(:,:), allocatable, target :: atoms_sph_in
       real(dp), dimension(3) :: diff
       real(dp), pointer :: r, theta, phi

       integer :: i, l, m, n, my_Z

       if( f_hat%initialised ) then
          do n = 1, f_hat%N_max
             do l = 0, f_hat%L_max
                f_hat%coeff%rad(n)%ang(l)%mag(:) = cplx_zero
             enddo
          enddo
       else
          call system_abort('fourier_transform: f_hat not initialised')
       endif

       my_Z = optional_default(at%Z(1),Z)

       if( .not. any( my_Z == at%Z(:) ) ) &
       & call print_warning('fourier_transform: no atoms '//Z//' in atoms object')

       allocate( atoms_sph_in( 3, at%N ) )

       atoms_sph_in = huge(1.0_dp)

       do i = 1, at%N
          if( my_Z /= at%Z(i) ) cycle
          diff = diff_min_image(at,i,vec)
          atoms_sph_in(:,i) = xyz2spherical( diff )
       enddo

       do n = 1, f_hat%N_max
          do l = 0, f_hat%L_max
             do m = -l, l
                do i = 1, at%N
                
                   if( my_Z /= at%Z(i) ) cycle

                   r => atoms_sph_in(1,i)
                   theta => atoms_sph_in(2,i)
                   phi => atoms_sph_in(3,i)
 
                   f_hat%coeff%rad(n)%ang(l)%mag(m) = f_hat%coeff%rad(n)%ang(l)%mag(m) + &
                   & gaussian(r,f_hat%r0(n),f_hat%sigma) * SphericalY( l, m, theta, phi ) 

                enddo
             enddo
          enddo
       enddo

       deallocate( atoms_sph_in )

     endsubroutine fourier_transform_vec

     !#################################################################################
     !#
     !% Calculates an atom-centered fourier spectrum
     !#
     !#################################################################################

     subroutine fourier_transform_atom(at,i,f_hat,Z)

       type(atoms), intent(in)                       :: at       !% atoms object
       integer, intent(in)                           :: i        !% central atom
       type(fourier_coefficient_atom), intent(inout) :: f_hat    !% fourier coefficents
       integer, intent(in), optional                 :: Z        !% which type of atoms to include

       real(dp), dimension(:,:), allocatable, target :: atoms_sph_in
       real(dp), dimension(:), allocatable           :: cutoff
       real(dp), dimension(3) :: diff
       real(dp) :: radial
       real(dp), pointer :: r, theta, phi

       integer :: j, l, n, nn, m, my_Z

       if( f_hat%initialised ) then
          do n = 1, f_hat%N_max
             do l = 0, f_hat%L_max
                f_hat%coeff%rad(n)%ang(l)%mag(:) = cplx_zero
             enddo
          enddo
       else
          call system_abort('fourier_transform: f_hat not initialised')
       endif

       my_Z = optional_default(at%Z(i),Z)

       if( .not. any( my_Z == at%Z(:) ) ) &
       & call print_warning('fourier_transform: no atoms '//Z//' in atoms object')

       allocate( atoms_sph_in( 3, atoms_n_neighbours(at,i) ), cutoff(atoms_n_neighbours(at,i)) )

       cutoff = 0.0_dp
       do nn = 1, atoms_n_neighbours(at,i)
          j = atoms_neighbour(at,i,nn,diff=diff)
          if( my_Z /= at%Z(j) ) cycle
          atoms_sph_in(:,nn) = xyz2spherical( diff )
          r => atoms_sph_in(1,nn)
          cutoff(nn) = cutoff_left(r,f_hat%cutoff_left,f_hat%cutoff_sigma) * &
                    & cutoff_right(r,f_hat%cutoff_right,f_hat%cutoff_sigma)
       enddo

       do n = 1, f_hat%N_max
          do l = 0, f_hat%L_max
             do m = -l, l
                do nn = 1, atoms_n_neighbours(at,i)
                   j = atoms_neighbour(at,i,nn)
                   if( my_Z /= at%Z(j) ) cycle
                   r => atoms_sph_in(1,nn)
                   theta => atoms_sph_in(2,nn)
                   phi => atoms_sph_in(3,nn)
 
                   radial = cutoff(nn) * gaussian(r,f_hat%r0(n),f_hat%sigma)

                   f_hat%coeff%rad(n)%ang(l)%mag(m) = f_hat%coeff%rad(n)%ang(l)%mag(m) + &
                   & radial * SphericalY( l, m, theta, phi ) 

                enddo
             enddo
          enddo
       enddo

       deallocate( atoms_sph_in, cutoff )

     endsubroutine fourier_transform_atom

     !#################################################################################
     !#
     !% Calculates an atom-centered bispectrum from a single fourier spectrum
     !#
     !#################################################################################

     subroutine calc_bispectrum1_atom(f_hat,bispectrum)

       type(fourier_coefficient_atom), intent(in) :: f_hat            !% fourier coefficients
       type(bispectrum_coefficient_atom), intent(inout) :: bispectrum !% bispectrum coefficients

       integer :: l1, l2, l, m, m1, N_max, L_max, n1, n2 
       complex(dp) :: sub_sum 

       if( .not. f_hat%initialised ) then
          call system_abort('calc_bispectrum1_atom: f_hat not initialised')
       endif

       if( .not. cg_initialised ) call cg_initialise(f_hat%L_max)

       if( bispectrum%initialised ) then
          do n1 = 1, bispectrum%N_max
             do n2 = 1, bispectrum%N_max
                do l1 = 0, bispectrum%L_max
                   do l2 = 0, bispectrum%L_max
                      bispectrum%coeff%rad(n2,n1)%ang(l2,l1)%mag(:) = cplx_zero
                   enddo
                enddo
             enddo
          enddo
       endif

       N_max = f_hat%N_max
       L_max = f_hat%L_max

       if( .not. bispectrum%initialised ) call initialise( bispectrum, N_max, L_max )

       if( bispectrum%L_max < L_max ) then
          call print_warning('calc_bispectrum1_atom: bispectrum%L_max < L_max, using bispectrum%L_max')
          L_max = bispectrum%L_max
       endif

       if( bispectrum%N_max < N_max ) then
          call print_warning('calc_bispectrum1_atom: bispectrum%N_max < N_max, using bispectrum%N_max')
          N_max = bispectrum%N_max
       endif

       do n1 = 1, N_max
          do n2 = 1, N_max
             do l1 = 0, L_max
                do l2 = 0, L_max
                   do l = abs(l1-l2), min(L_max,l1+l2)
                      do m = -l, l
                         sub_sum = cplx_zero
                         do m1 = max( -l1,m-l2), min(l1,m+l2)
       
                            sub_sum = sub_sum + &
                            & cg_array(l1,m1,l2,m-m1,l,m) * & 
                            & conjg(f_hat%coeff%rad(n1)%ang(l1)%mag(m1) * f_hat%coeff%rad(n2)%ang(l2)%mag(m-m1))
       
                         enddo
                      bispectrum%coeff%rad(n2,n1)%ang(l2,l1)%mag(l) = &
                      & bispectrum%coeff%rad(n2,n1)%ang(l2,l1)%mag(l) + f_hat%coeff%rad(n2)%ang(l)%mag(m) * sub_sum
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

     endsubroutine calc_bispectrum1_atom

     !#################################################################################
     !#
     !% Calculates an atom-centered bispectrum from two fourier spectrums
     !#
     !#################################################################################

!     subroutine calc_bispectrum2_atom(f_hat1,f_hat2,bispectrum)
!
!       type(fourier_coefficient_atom), intent(in) :: f_hat1            !% fourier spectrum #1
!       type(fourier_coefficient_atom), intent(in) :: f_hat2            !% fourier spectrum #2
!       type(bispectrum_coefficient_atom), intent(inout) :: bispectrum  !% bispectrum
!
!       integer :: l1, l2, l, m, m1, N_max, L_max, n 
!       complex(dp) :: sub_sum 
!
!       if( .not. (f_hat1%initialised.and.f_hat2%initialised) ) &
!       & call system_abort('calc_bispectrum2_atom: f_hat1 or f_hat2 not initialised')
!
!       if( (f_hat1%L_max/=f_hat2%L_max) .or. (f_hat1%N_max/=f_hat2%N_max) ) &
!       & call system_abort('calc_bispectrum2_atom: L_max or N_max does not match')
!
!       if( .not. cg_initialised ) call cg_initialise(f_hat1%L_max,f_hat1%L_max,f_hat1%L_max,f_hat1%L_max,f_hat1%L_max,f_hat1%L_max)
!
!       if( bispectrum%initialised ) then
!          do n = 1, (bispectrum%N_max-1)/2
!             do l1 = 0, bispectrum%L_max
!                do l2 = 0, bispectrum%L_max
!                   bispectrum%coeff%rad(n)%ang(l2,l1)%mag(:) = cplx_zero
!                enddo
!             enddo
!          enddo
!       endif
!
!       N_max = f_hat1%N_max
!       L_max = f_hat1%L_max
!
!       if( .not. bispectrum%initialised ) call initialise( bispectrum, N_max, L_max )
!
!       if( bispectrum%L_max < L_max ) then
!          call print_warning('calc_bispectrum2_atom: bispectrum%L_max < L_max, using bispectrum%L_max')
!          L_max = bispectrum%L_max
!       endif
!
!       if( bispectrum%N_max < N_max ) then
!          call print_warning('calc_bispectrum2_atom: bispectrum%N_max < N_max, using bispectrum%N_max')
!          N_max = bispectrum%N_max
!       endif
!
!       do n = 1, (N_max-1)/2
!          do l1 = 0, L_max
!             do l2 = 0, L_max
!                do l = abs(l1-l2), min(L_max,l1+l2)
!                   do m = -l, l
!                      sub_sum = cplx_zero
!                      do m1 = max( -l1,m-l2), min(l1,m+l2)
!       
!                         sub_sum = sub_sum + &
!                         & cg_array(l1,m1,l2,m-m1,l,m) * & 
!                         & conjg(f_hat1%coeff%rad(2*n-1)%ang(l1)%mag(m1) * &
!                         & f_hat2%coeff%rad(2*n)%ang(l2)%mag(m-m1))
!       
!                      enddo
!                      bispectrum%coeff%rad(n)%ang(l2,l1)%mag(l) = &
!                      & bispectrum%coeff%rad(n)%ang(l2,l1)%mag(l) + &
!                      & f_hat1%coeff%rad(2*n+1)%ang(l)%mag(m) * sub_sum
!                   enddo
!                enddo
!             enddo
!          enddo
!       enddo
!
!     endsubroutine calc_bispectrum2_atom

     !#################################################################################
     !#
     !% Converts a pair environment from xyz to cylindrical coordinates
     !#
     !#################################################################################

     subroutine xyz2cyl(at,i,j,cyl)

        type(atoms), intent(in)   :: at    !% atoms object
        integer, intent(in)       :: i, j  !% atom pair
        type(table), intent(out)  :: cyl   !% cylindrical coordinates

        integer :: k, n

        real(dp), dimension(3) :: diff_ij, diff_ik, diff_jk, diff_ref, cross_ref
        real(dp) :: r_ij, r_ik, r_jk, t, r_ik_dot_r_ij, r_jk_dot_r_ij, r, z, costheta, theta, r_ref, t_ref

        logical :: theta_flag, clockwise

        call allocate(cyl, Nint=2, Nreal=3, Nstr=0, Nlogical=0,max_length=0)

        diff_ij = diff_min_image(at, i, j)
        r_ij = distance_min_image(at, i, j)
        theta_flag = .false.

        do n = 1, atoms_n_neighbours(at,i)
           k = atoms_neighbour(at, i, n, distance=r_ik, diff=diff_ik)
           if( k == j ) cycle

           r_ik_dot_r_ij = dot_product(diff_ik,diff_ij)
           t = r_ik_dot_r_ij / r_ij**2
           r = sqrt( r_ik**2 + t**2 * r_ij**2 - 2.0_dp*t*r_ik_dot_r_ij )
           z = (0.5_dp - t) * r_ij

           if( r .fne. 0.0_dp ) then
               if( theta_flag ) then
                   cross_ref = (diff_ref .cross. diff_ik) + &
                   & ( (t_ref*diff_ik - t*diff_ref) .cross. diff_ij )
                   clockwise = (dot_product(cross_ref,diff_ij) > 0.0_dp )
                   costheta = ( dot_product(diff_ref,diff_ik) + t_ref*t*r_ij**2 - &
                   & t * dot_product(diff_ref,diff_ij) - t_ref * dot_product(diff_ik,diff_ij) ) / &
                   & (r_ref*r)
               else
                   theta_flag = .true.
                   clockwise = .true. 
                   costheta = 1.0_dp
                   diff_ref = diff_ik
                   t_ref = t
                   r_ref = r
               endif
           else
               clockwise = .true. 
               costheta = 1.0_dp
           endif
           theta = acos(costheta)
           if( .not. clockwise ) theta = 2*PI - theta
           call append(cyl,(/k,at%Z(k)/),(/r,z,theta/)) 
        enddo
           
        do n = 1, atoms_n_neighbours(at,j)
           k = atoms_neighbour(at, j, n, distance=r_jk, diff=diff_jk)

           if( k == i ) cycle
           if( any(k == cyl%int(1,:) ) ) cycle

           r_jk_dot_r_ij = dot_product(diff_jk,diff_ij)
           t = r_jk_dot_r_ij / r_ij**2
           r = sqrt( r_jk**2 + t**2 * r_ij**2 - 2.0_dp*t*r_jk_dot_r_ij )
           z = - (0.5_dp + t) * r_ij

           if( r .fne. 0.0_dp ) then
               if( theta_flag ) then
                   cross_ref = (diff_ref .cross. diff_jk) + &
                   & ( (t_ref*diff_jk - t*diff_ref) .cross. diff_ij )
                   clockwise = (dot_product(cross_ref,diff_ij) > 0.0_dp )
                   costheta = ( dot_product(diff_ref,diff_jk) + t_ref*t*r_ij**2 - &
                   & t * dot_product(diff_ref,diff_ij) - t_ref * dot_product(diff_jk,diff_ij) ) / &
                   & (r_ref*r)
               else
                   theta_flag = .true.
                   clockwise = .true.
                   costheta = 1.0_dp
                   diff_ref = diff_jk
                   t_ref = t
                   r_ref = r
               endif
           else
               clockwise = .true.
               costheta = 1.0_dp
           endif
           theta = acos(costheta)
           if( .not. clockwise ) theta = 2*PI - theta
           call append(cyl,(/k,at%Z(k)/),(/r,z,theta/)) 
        enddo
           
     endsubroutine xyz2cyl

     !#################################################################################
     !#
     !% Converts a pair environment from xyz to 'elliptical' or bispherical coordinates
     !#
     !#################################################################################

     subroutine xyz2ell(at,i,j,ell)

        type(atoms), intent(in)  :: at     !% atoms object
        integer, intent(in)      :: i, j   !% atom pair
        type(table), intent(out) :: ell    !% elliptical coordinates

        integer :: k, n

        real(dp), dimension(3) :: diff_ij, diff_ik, diff_jk, diff_pk, diff_ref, ref_cross_pk
        real(dp) :: r_ij, r_ik, r_jk, r_pk, t, costheta, theta, r_ref, t_ref

        logical :: theta_flag, clockwise

        call allocate(ell, Nint=2, Nreal=3, Nstr=0, Nlogical=0)

        diff_ij = diff_min_image(at, i, j)
        r_ij = distance_min_image(at, i, j)
        theta_flag = .false.

        do n = 1, atoms_n_neighbours(at,i)
           k = atoms_neighbour(at, i, n, distance=r_ik, diff=diff_ik)
           if( k == j ) cycle
           r_jk = distance_min_image(at,j,k)

           t = dot_product(diff_ik,diff_ij) / r_ij**2
           diff_pk = diff_ik - t * diff_ij
           r_pk = norm(diff_pk)

           if( r_pk .fne. 0.0_dp ) then
               if( theta_flag ) then
                   ref_cross_pk = (diff_ref .cross. diff_pk)
                   clockwise = (dot_product(ref_cross_pk,diff_ij) > 0.0_dp )

                   costheta = dot_product(diff_ref,diff_pk) / (r_ref*r_pk)
               else
                   theta_flag = .true.
                   clockwise = .true.
                   costheta = 1.0_dp
                   diff_ref = diff_pk
                   t_ref = t
                   r_ref = r_pk
               endif
           else
               clockwise = .true.
               costheta = 1.0_dp
           endif
           theta = acos(costheta)
           if( .not. clockwise ) theta = 2*PI - theta
           call append(ell,(/k,at%Z(k)/),(/r_ik,r_jk,theta/)) 
        enddo
           
        do n = 1, atoms_n_neighbours(at,j)
           k = atoms_neighbour(at, j, n, distance=r_jk, diff=diff_jk)

           if( k == i ) cycle
           if( any(k == ell%int(1,:) ) ) cycle

           r_ik = distance_min_image(at,i,k)

           t = dot_product(diff_jk,diff_ij) / r_ij**2
           diff_pk = diff_jk - t * diff_ij
           r_pk = norm(diff_pk)

           if( r_pk .fne. 0.0_dp ) then
               if( theta_flag ) then
                   ref_cross_pk = (diff_ref .cross. diff_pk)
                   clockwise = (dot_product(ref_cross_pk,diff_ij) > 0.0_dp )

                   costheta = dot_product(diff_ref,diff_pk) / (r_ref*r_pk)
               else
                   theta_flag = .true.
                   clockwise = .true.
                   costheta = 1.0_dp
                   diff_ref = diff_pk
                   t_ref = t
                   r_ref = r_pk
               endif
           else
               clockwise = .true.
               costheta = 1.0_dp
           endif
           theta = acos(costheta)
           if( .not. clockwise ) theta = 2*PI - theta
           call append(ell,(/k,at%Z(k)/),(/r_ik,r_jk,theta/)) 
        enddo
           
     endsubroutine xyz2ell

     !#################################################################################
     !#
     !% Initialise a fourier object for atom pairs.
     !#
     !#################################################################################

     subroutine initialise_fourier_pair(f_hat,N_max,L_max,r_cut,r_min,r_max,sigma, &
                & cutoff_left,cutoff_right,cutoff_sigma)

        type(fourier_coefficient_pair), intent(inout) :: f_hat  !% fourier coefficients
        integer, intent(in)  :: L_max                           !% angular basis functions
        integer, intent(in)  :: N_max                           !% radial basis functions
        real(dp), intent(in) :: r_cut                           !% elliptical cutoff value
        real(dp), intent(in) :: r_min                           !% centre of first radial basis function
        real(dp), intent(in) :: r_max                           !% centre of last radial basis function
        real(dp), intent(in) :: sigma                           !% width of radial basis function
        real(dp), intent(in) :: cutoff_left                     !% centre of the 'left' cutoff fermi function
        real(dp), intent(in) :: cutoff_right                    !% centre of the 'right' cutoff fermi function
        real(dp), intent(in) :: cutoff_sigma                    !% width of the cutoff fermi functions
    
        integer  :: n
        real(dp) :: dr

        if( mod(L_max,2) /= 0 ) &
        & call system_abort('initialise_fourier_pair: L_max must be even, here it is '//L_max//' clearly odd.')

        if( f_hat%initialised ) call finalise(f_hat)
        f_hat%r_cut = r_cut
        f_hat%r_ij = 0.0_dp
        allocate( f_hat%coeff(N_max,N_max,-L_max:L_max), f_hat%r0(N_max) )

        f_hat%coeff = CPLX_ZERO

        dr = (r_max-r_min)/(N_max-1)
        f_hat%r0(1) = r_min
        do n = 2, N_max
           f_hat%r0(n) = f_hat%r0(n-1)+dr
        enddo
        f_hat%sigma = sigma
        f_hat%cutoff_left  = cutoff_left
        f_hat%cutoff_right = cutoff_right
        f_hat%cutoff_sigma = cutoff_sigma
    
        f_hat%N_max = N_max
        f_hat%L_max = L_max
        f_hat%initialised = .true.
    
     endsubroutine initialise_fourier_pair
    
     !#################################################################################
     !#
     !% Finalise a fourier object for atom pairs.
     !#
     !#################################################################################

     subroutine finalise_fourier_pair(f_hat)

        type(fourier_coefficient_pair), intent(inout) :: f_hat !% fourier coefficients object

        if( .not. f_hat%initialised ) return
        f_hat%r_cut = 0.0_dp
        f_hat%r_ij = 0.0_dp

        deallocate( f_hat%coeff, f_hat%r0 )

        f_hat%sigma = 0.0_dp
        f_hat%cutoff_left  = 0.0_dp
        f_hat%cutoff_right = 0.0_dp
        f_hat%cutoff_sigma = 0.0_dp
    
        f_hat%N_max = 0
        f_hat%L_max = 0
        f_hat%initialised = .false.
    
     endsubroutine finalise_fourier_pair
    
     !#################################################################################
     !#
     !% Initialise a bispectrum object for an atom pair.
     !#
     !#################################################################################

     subroutine initialise_bispectrum_pair(this,N_max,L_max)

        type(bispectrum_coefficient_pair), intent(inout) :: this       !% bispectrum coefficients
        integer, intent(in)                              :: N_max      !% number of radial basis functions
        integer, intent(in)                              :: L_max      !% number of angular basis functions

        if( this%initialised ) call finalise(this)

        allocate(this%coeff(N_max,N_max,-L_max:L_max))
        this%L_max=L_max
        this%N_max=N_max
        this%initialised = .true.

     endsubroutine initialise_bispectrum_pair
     
     !#################################################################################
     !#
     !% Finalise a bispectrum object for an atom pair.
     !#
     !#################################################################################

     subroutine finalise_bispectrum_pair(this)

        type(bispectrum_coefficient_pair), intent(inout) :: this !% bispectrum coefficients

        if( .not. this%initialised ) return

        deallocate(this%coeff)
        this%L_max=0
        this%N_max=0
        this%initialised = .false.

     endsubroutine finalise_bispectrum_pair
     
     !#################################################################################
     !#
     !% Performs fourier transform on atom pairs.
     !#
     !#################################################################################

     subroutine fourier_transform_pair(at,i,j,f_hat)

        type(atoms), intent(in)                       :: at    !% atoms object
        integer, intent(in)                           :: i, j  !% atom pair
        type(fourier_coefficient_pair), intent(inout) :: f_hat !% fourier object

        type(table) :: ell
        real(dp), dimension(:), allocatable :: cutoff
        real(dp), dimension(:,:), allocatable :: g1, g2
        complex(dp), dimension(:,:), allocatable :: c
        real(dp) :: r_ij, r12, r1, r2, theta
        integer :: n, n1, n2, m, nn

        if( .not. f_hat%initialised ) call system_abort('fourier_transform_pair: f_hat not initialised')

        r_ij = distance_min_image(at,i,j)
        r12 = r_ij + 2.0_dp*f_hat%r_cut
        if( at%cutoff < r12/2.0_dp ) then
            !call print_warning('atoms cutoff to small, calling calc_connect')
            !call set_cutoff(at,r12/2.0_dp)
            !call calc_connect(at)
            call system_abort('fourier_transform_pair: atoms cutoff to small')
        endif
        call xyz2ell(at,i,j,ell)


        f_hat%r_ij = r_ij
        f_hat%coeff = CPLX_ZERO

        allocate( cutoff(ell%N),g1(ell%N,f_hat%N_max),g2(ell%N,f_hat%N_max),&
        & c(ell%N,-f_hat%L_max:f_hat%L_max) )

        do nn = 1, ell%N
           r1 = ell%real(1,nn)
           r2 = ell%real(2,nn)
           cutoff(nn) = cutoff_left(r1,f_hat%cutoff_left,f_hat%cutoff_sigma) * &
                      & cutoff_left(r2,f_hat%cutoff_left,f_hat%cutoff_sigma) * &
                      & cutoff_right(r1+r2,r12,f_hat%cutoff_sigma)
        enddo
        do n = 1, f_hat%N_max
           do nn = 1, ell%N
              r1 = ell%real(1,nn)
              r2 = ell%real(2,nn)
              g1(nn,n) = gaussian(r1,f_hat%r0(n),f_hat%sigma)
              g2(nn,n) = gaussian(r2,f_hat%r0(n),f_hat%sigma)
           enddo
        enddo
        do m = -f_hat%L_max, f_hat%L_max
           do nn = 1, ell%N
              theta = ell%real(3,nn)
              c(nn,m) = exp( CPLX_IMAG * m * theta )
           enddo
        enddo

        do m = -f_hat%L_max, f_hat%L_max
           do n2 = 1, f_hat%N_max
              do n1 = 1, f_hat%N_max
                 f_hat%coeff(n1,n2,m) = sum( cutoff(:) * g1(:,n1) * g2(:,n2) * c(:,m) )
              enddo
           enddo
        enddo
        deallocate(cutoff,g1,g2,c)
        call finalise(ell)

     endsubroutine fourier_transform_pair

     !#################################################################################
     !#
     !% Calculates bispectrum on a pair of atoms.
     !#
     !#################################################################################

     subroutine calc_bispectrum_pair(f_hat,bispectrum)

        type(fourier_coefficient_pair), intent(in)       :: f_hat      !% fourier coefficients
        type(bispectrum_coefficient_pair), intent(inout) :: bispectrum !% bispectrum coefficients

        integer :: m, n1, n2
        if( .not. f_hat%initialised ) then
          call system_abort('calc_bispectrum_pair: f_hat not initialised')
        endif

        if( bispectrum%initialised ) then
            if( (bispectrum%L_max == f_hat%L_max/2) .and. (bispectrum%N_max == f_hat%N_max) ) then
                bispectrum%coeff = CPLX_ZERO
	    else
                call finalise( bispectrum )
                call initialise( bispectrum, f_hat%N_max, f_hat%L_max/2 )
            endif
        else
            call initialise( bispectrum, f_hat%N_max, f_hat%L_max/2 )
        endif

        bispectrum%r_ij = f_hat%r_ij

        do m = -bispectrum%L_max, bispectrum%L_max
           do n2 = 1, bispectrum%N_max
              do n1 = 1, bispectrum%N_max
                 bispectrum%coeff(n1, n2, m) = conjg(f_hat%coeff(n1,n2,2*m))*&
                 & f_hat%coeff(n1,n2,m)*f_hat%coeff(n1,n2,m)
              enddo
           enddo
        enddo
     endsubroutine calc_bispectrum_pair

     !#################################################################################
     !#
     !% Transposes matrices in a bispectrum for each angular component.
     !#
     !#################################################################################

     subroutine transpose_bispectrum_pair(this,other)

        type(bispectrum_coefficient_pair), intent(inout) :: this, other !% bispectrum coefficients

        integer :: m

        if( .not. this%initialised ) call system_abort('transpose_bispectrum_pair: bispectrum not initialised')

        if( other%initialised ) then
            if( (this%N_max /= other%N_max) .or. (this%L_max /= other%L_max) ) then
                call finalise(other)
                call initialise(other, this%N_max, this%L_max)
            endif
        else
            call initialise(other, this%N_max, this%L_max)
        endif

        other%r_ij = this%r_ij
        do m = -this%L_max, this%L_max
           other%coeff(:,:,m) = transpose(this%coeff(:,:,m))
        enddo

     endsubroutine transpose_bispectrum_pair

     !#################################################################################
     !#
     !% Conjugates each bispectrum coefficient.
     !#
     !#################################################################################

     subroutine conjugate_bispectrum_pair(this,other)

        type(bispectrum_coefficient_pair), intent(inout) :: this, other !% bispectrum coefficients

        if( .not. this%initialised ) call system_abort('conjugate_bispectrum_pair: bispectrum not initialised')
        if( other%initialised ) then
            if( (this%N_max /= other%N_max) .or. (this%L_max /= other%L_max) ) then
                call finalise(other)
                call initialise(other, this%N_max, this%L_max)
            endif
        else
            call initialise(other, this%N_max, this%L_max)
        endif

        other%r_ij = this%r_ij
        other%coeff(:,:,:) = conjg(this%coeff(:,:,:))

     endsubroutine conjugate_bispectrum_pair

     !#################################################################################
     !#
     !% Prints fourier coefficients.
     !#
     !#################################################################################

     subroutine print_fourier_pair(this)

        type(fourier_coefficient_pair), intent(in) :: this !% fourier object

        integer  :: n1, n2, m

        if( .not. this%initialised ) call system_abort('print_fourier_pair: not initialised')

        print'("r_cut = ",f10.4)', this%r_cut
        print'("r_ij  = ",f10.4)', this%r_ij

        print'(3a5,2a10)','n1','n2','m','Re(x)','Im(x)'
        do m = -this%L_max, this%L_max
           do n2 = 1, this%N_max
              do n1 = 1, this%N_max
                 print'(3i5,2f10.4)',n1,n2,m,this%coeff(n1,n2,m)
              enddo
           enddo
        enddo

     endsubroutine print_fourier_pair
    
     !#################################################################################
     !#
     !% Prints fourier coefficients.
     !#
     !#################################################################################

     subroutine print_fourier_atom(this)

       type(fourier_coefficient_atom), intent(in) :: this
       integer :: l, m, n

       if( .not. this%initialised ) call system_abort('print_fourier_atom: not initialised')

       print'(2a5,2a10)','l','m','Re(x)','Im(x)'

       do n = 1, this%N_max
          do l = 0, this%L_max
             do m = -l, l
                print'(3i5,2f10.4)',n,l,m,this%coeff%rad(n)%ang(l)%mag(m)
             enddo
          enddo
       enddo

     endsubroutine print_fourier_atom

     !#################################################################################
     !#
     !% Prints bispectrum coefficients.
     !#
     !#################################################################################

     subroutine print_bispectrum_pair(this)

        type(bispectrum_coefficient_pair), intent(in) :: this ! bispectrum object

        integer  :: n1, n2, m

        if( .not. this%initialised ) call system_abort('print_bispectrum_pair: not initialised')

        print'(3a5,2a10)','n1','n2','m','Re(x)','Im(x)'
        do m = -this%L_max, this%L_max
           do n2 = 1, this%N_max
              do n1 = 1, this%N_max
                 print'(3i5,2f10.4)',n1,n2,m,this%coeff(n1,n2,m)
              enddo
           enddo
        enddo

     endsubroutine print_bispectrum_pair
    
     !#################################################################################
     !#
     !% Prints bispectrum coefficients.
     !#
     !#################################################################################

     subroutine print_bispectrum_atom(this)

       type(bispectrum_coefficient_atom), intent(in) :: this
       integer :: l1,l2,l, n1, n2

       if( .not. this%initialised ) call system_abort('print_bispectrum_atom: not initialised')

       print'(5a5,2a10)','n1','n2','l_1','l_2','l','Re(s)','Im(s)'

       do n1 = 1, this%N_max
          do n2 = 1, this%N_max
             do l1 = 0, this%L_max
                do l2 = 0, this%L_max
                   do l = abs(l1-l2), min(this%L_max,l1+l2)
                      print'(5i5,2f10.4)',n1,n2,l1,l2,l, this%coeff%rad(n2,n1)%ang(l2,l1)%mag(l)
                   enddo
                enddo
             enddo
          enddo
       enddo

     endsubroutine print_bispectrum_atom

     !#################################################################################
     !#
     !% Determines how many coefficients are in a bispectrum object (atom)
     !#
     !#################################################################################

     function L_max2d_atom(N_max,L_max)

       integer, intent(in) :: L_max !% number of angular basis functions
       integer, intent(in) :: N_max !% number of radial basis functions
       integer :: L_max2d_atom

       integer :: l1, l2, l, n1, n2

       L_max2d_atom = 0

       do n1 = 1, N_max
          do n2 = 1, N_max
             do l1 = 0, L_max
                do l2 = 0, L_max
                   do l = abs(l1-l2), min(L_max,l1+l2)
                      L_max2d_atom = L_max2d_atom + 2
                   enddo
                enddo
             enddo
          enddo
       enddo

     endfunction L_max2d_atom

     !#################################################################################
     !#
     !% Determines how many coefficients are in a bispectrum object (pair)
     !#
     !#################################################################################

     function L_max2d_pair(N_max,L_max)

       integer, intent(in) :: L_max !% number of angular basis functions
       integer, intent(in) :: N_max !% number of radial basis functions
       integer :: L_max2d_pair

       L_max2d_pair = 2 * (L_max*2+1) * N_max**2 + 1

     endfunction L_max2d_pair

     !#################################################################################
     !#
     !% Determines how many coefficients are in a bispectrum object (atom)
     !#
     !#################################################################################

     function L_max2d_bispectrum_atom(this1,this2,this3) 

       type(bispectrum_coefficient_atom), intent(in) :: this1           !% bispectrum object
       type(bispectrum_coefficient_atom), intent(in), optional :: this2 !% bispectrum object
       type(bispectrum_coefficient_atom), intent(in), optional :: this3 !% bispectrum object
       integer :: L_max2d_bispectrum_atom

       L_max2d_bispectrum_atom = L_max2d_atom(this1%N_max,this1%L_max)
       if( present(this2) ) L_max2d_bispectrum_atom = L_max2d_bispectrum_atom + L_max2d_atom(this2%N_max,this2%L_max)
       if( present(this3) ) L_max2d_bispectrum_atom = L_max2d_bispectrum_atom + L_max2d_atom(this3%N_max,this3%L_max)

     endfunction L_max2d_bispectrum_atom

     !#################################################################################
     !#
     !% Determines how many coefficients are in a bispectrum object (pair)
     !#
     !#################################################################################

     function L_max2d_bispectrum_pair(this1,this2,this3) 

       type(bispectrum_coefficient_pair), intent(in) :: this1           !% bispectrum object
       type(bispectrum_coefficient_pair), intent(in), optional :: this2 !% bispectrum object
       type(bispectrum_coefficient_pair), intent(in), optional :: this3 !% bispectrum object
       integer :: L_max2d_bispectrum_pair

       L_max2d_bispectrum_pair = L_max2d_pair(this1%N_max,this1%L_max)
       if( present(this2) ) L_max2d_bispectrum_pair = L_max2d_bispectrum_pair + L_max2d_pair(this2%N_max,this2%L_max)
       if( present(this3) ) L_max2d_bispectrum_pair = L_max2d_bispectrum_pair + L_max2d_pair(this3%N_max,this3%L_max)

     endfunction L_max2d_bispectrum_pair

     !#################################################################################
     !#
     !% Copies bispectrum to an ordinary vector
     !#
     !#################################################################################

     subroutine bispectrum2vec1_atom(this,vec)

       type(bispectrum_coefficient_atom), intent(in) :: this !% bispectrum
       real(dp), dimension(:), intent(out)      :: vec  !% array, dimension(L_max2d(this))

       integer :: l1, l2, l, n1, n2, i

       if( .not. this%initialised ) then
           call system_abort('bispectrum2vec1_atom: not initialised')
       endif

       if( L_max2d(this) > size(vec) ) then
          call system_abort('bispectrum2vec1_atom: vec too small')
       endif

       vec = 0.0_dp
       i = 1

       do n1 = 1, this%N_max
          do n2 = 1, this%N_max
             do l1 = 0, this%L_max
                do l2 = 0, this%L_max
                   do l = abs(l1-l2), min(this%L_max,l1+l2)
                      vec(i) = real(this%coeff%rad(n2,n1)%ang(l2,l1)%mag(l))
                      vec(i+1) = aimag(this%coeff%rad(n2,n1)%ang(l2,l1)%mag(l))
                      i=i+2
                   enddo
                enddo
             enddo
          enddo
       enddo

     endsubroutine bispectrum2vec1_atom

     !#################################################################################
     !#
     !% Copies two bispectrums to an ordinary vector
     !#
     !#################################################################################

     subroutine bispectrum2vec2_atom(this1,this2,vec)

       type(bispectrum_coefficient_atom), intent(in) :: this1, this2  !% bispectrums
       real(dp), dimension(:), intent(out)      :: vec           !% vector, dimension(L_max2d(this1,this2))

       if( .not. (this1%initialised.and.this2%initialised) ) then
           call system_abort('bispectrum2vec1_atom: not initialised')
       endif

       if( L_max2d(this1,this2) > size(vec) ) then
          call system_abort('bispectrum2vec1_atom: vec too small')
       endif

       vec = 0.0_dp

       call bispectrum2vec1_atom( this1,vec(1:L_max2d(this1)) )
       call bispectrum2vec1_atom( this2,vec(L_max2d(this1)+1:L_max2d(this1,this2)) )

     endsubroutine bispectrum2vec2_atom

     !#################################################################################
     !#
     !% Copies three bispectrums to an ordinary vector
     !#
     !#################################################################################

     subroutine bispectrum2vec3_atom(this1,this2,this3,vec)

       type(bispectrum_coefficient_atom), intent(in) :: this1, this2, this3 !% bispectrums
       real(dp), dimension(:), intent(out)      :: vec                 !% vector, dimension(L_max2d(this1,this2))

       if( .not. (this1%initialised.and.this2%initialised.and.this3%initialised) ) then
           call system_abort('bispectrum2vec3_atom: not initialised')
       endif

       if( L_max2d(this1,this2,this3) > size(vec) ) then
          call system_abort('bispectrum2vec3_atom: vec too small')
       endif

       vec = 0.0_dp
       call bispectrum2vec1_atom( this1,vec(1:L_max2d(this1)) )
       call bispectrum2vec1_atom( this2,vec(L_max2d(this1)+1:L_max2d(this1,this2)) )
       call bispectrum2vec1_atom( this3,vec(L_max2d(this1,this2)+1:L_max2d(this1,this2,this3)) )

     endsubroutine bispectrum2vec3_atom

     !#################################################################################
     !#
     !% Copies bispectrum to an ordinary vector (pair)
     !#
     !#################################################################################

     subroutine bispectrum2vec1_pair(this,vec)

       type(bispectrum_coefficient_pair), intent(in) :: this !% bispectrum
       real(dp), dimension(:), intent(out)      :: vec  !% array, dimension(L_max2d(this))

       integer :: n1, n2, m, i

       if( .not. this%initialised ) then
           call system_abort('bispectrum2vec1_pair: not initialised')
       endif

       if( L_max2d(this) > size(vec) ) then
          call system_abort('bispectrum2vec1_pair: vec too small')
       endif

       vec = 0.0_dp
       vec(1) = this%r_ij
       i = 2

       do m = -this%L_max, this%L_max
          do n2 = 1, this%N_max
             do n1 = 1, this%N_max
                vec(i) = real(this%coeff(n1,n2,m))
                vec(i+1) = aimag(this%coeff(n1,n2,m))
                i=i+2
             enddo
          enddo
       enddo

     endsubroutine bispectrum2vec1_pair
    
     !###############################################################
     !#
     !# special water descriptor functions
     !#
     !###############################################################

     function water_monomer(at,w) result(vec)
       type(atoms), intent(in) :: at
       integer, dimension(3), intent(in) :: w
       real(dp), dimension(3) :: vec, v1, v2
       real(dp) :: r1, r2
       integer :: iO, iH1, iH2

       iO = w(1)
       iH1 = w(2)
       iH2 = w(3)

       v1 = diff_min_image(at,iO,iH1)
       v2 = diff_min_image(at,iO,iH2)

       r1 = norm(v1)
       r2 = norm(v2)

       ! descriptors
       !vec(1) = normsq(v1+v2)
       !vec(2) = normsq(v1-v2)
       !vec(3) = ((v1+v2).dot.(v1-v2))**2
       vec(1) = r1 + r2
       vec(2) = (r1 - r2)**2
       vec(3) = dot_product(v1,v2)

     end function water_monomer

     subroutine water_dimer(at,w1,w2,cutoff, vec, dvec, rOHout, rHHout)
       type(atoms), intent(in) :: at
       integer, dimension(3), intent(in) :: w1, w2
       real(dp), dimension(WATER_DIMER_D), intent(out), optional     :: vec !, v
       real(dp), dimension(3,6,WATER_DIMER_D), intent(out), optional :: dvec !, v
       real(dp), intent(out), optional :: rOHout(8), rHHout(6)
       real(dp), intent(in) :: cutoff
       real(dp) :: rOH(8), rHH(6), drOH(3,6,8), drHH(3,6,6), &
       fOH(FOURIER_N_COMP), fHH(FOURIER_N_COMP), dfOH(3,6,FOURIER_N_COMP), dfHH(3,6,FOURIER_N_COMP), &
       arg, arg_r, fOH_ij, fHH_ij
       integer :: iAo, iAh1, iAh2, iBo, iBh1, iBh2, i, j, k
       logical, parameter :: DO_FOURIER = .true.

       real(dp), dimension(8), parameter :: r0_OH = (/0.92_dp, 0.95_dp, 0.98_dp, 1.01_dp, 1.90_dp, 2.50_dp, 3.25_dp, 4.00_dp/)
       real(dp), dimension(8), parameter :: sigma_OH = 1.0_dp / (/0.015_dp, 0.015_dp, 0.015_dp, 0.015_dp, 0.125_dp, 0.4_dp, 0.4_dp, 0.4_dp/)**2

       real(dp), dimension(6), parameter :: r0_HH = (/1.51_dp, 1.57_dp, 2.40_dp, 2.90_dp, 3.70_dp, 4.50_dp/)
       real(dp), dimension(6), parameter :: sigma_HH = 1.0_dp / (/0.05_dp, 0.05_dp, 0.25_dp, 0.40_dp, 0.60_dp, 0.60_dp/)**2

       !real(dp) :: rA1, rA2, rB1, rB2
       !real(dp), dimension(3) :: vA1, vA2, vB1, vB2, sA, sB, nA, nB, dA, dB, vAB

       ! Descriptors based on the interatomic distances. In theory, a
       ! configuration can be regenerated based on the interatomic distances. In
       ! this approach we take these and make them permutationally invariant.

       if( count( (/present(vec), present(dvec), present(rOHout), present(rHHout)/) ) == 0 ) return ! nothing to do

       ! Atomic indices
       iAo = w1(1)    ! 1
       iAh1 = w1(2)   ! 2
       iAh2 = w1(3)   ! 3
       iBo = w2(1)    ! 4
       iBh1 = w2(2)   ! 5
       iBh2 = w2(3)   ! 6

       ! All H-O distances
       rOH(1) = distance_min_image(at, iAo, iAh1)
       rOH(2) = distance_min_image(at, iAo, iAh2)
       rOH(3) = distance_min_image(at, iAo, iBh1)
       rOH(4) = distance_min_image(at, iAo, iBh2)
       rOH(5) = distance_min_image(at, iBo, iAh1)
       rOH(6) = distance_min_image(at, iBo, iAh2)
       rOH(7) = distance_min_image(at, iBo, iBh1)
       rOH(8) = distance_min_image(at, iBo, iBh2)
       

       ! All H-H distances
       rHH(1) = distance_min_image(at, iAh1, iAh2)
       rHH(2) = distance_min_image(at, iAh1, iBh1)
       rHH(3) = distance_min_image(at, iAh1, iBh2)
       rHH(4) = distance_min_image(at, iAh2, iBh1)
       rHH(5) = distance_min_image(at, iAh2, iBh2)
       rHH(6) = distance_min_image(at, iBh1, iBh2)

       if(present(rOHout)) then
          rOHout = rOH
       end if
       if(present(rHHout)) then
          rHHout = rHH
       end if
       
       if(present(dvec)) then
          drOH = 0.0_dp
          drHH = 0.0_dp

          drOH(:,1,1) = - diff_min_image(at, iAo, iAh1) / rOH(1) ! d r_{OH_1} / d r_1
          drOH(:,2,1) = - drOH(:,1,1)                            ! d r_{OH_1} / d r_2

          drOH(:,1,2) = - diff_min_image(at, iAo, iAh2) / rOH(2) ! d r_{OH_2} / d r_1
          drOH(:,3,2) = - drOH(:,1,2)                            ! d r_{OH_2} / d r_3

          drOH(:,1,3) = - diff_min_image(at, iAo, iBh1) / rOH(3) ! d r_{OH_3} / d r_1
          drOH(:,5,3) = - drOH(:,1,3)                            ! d r_{OH_3} / d r_5

          drOH(:,1,4) = - diff_min_image(at, iAo, iBh2) / rOH(4) ! d r_{OH_4} / d r_1
          drOH(:,6,4) = - drOH(:,1,4)                            ! d r_{OH_4} / d r_6

          drOH(:,4,5) = - diff_min_image(at, iBo, iAh1) / rOH(5) ! d r_{OH_5} / d r_4
          drOH(:,2,5) = - drOH(:,4,5)                            ! d r_{OH_5} / d r_2

          drOH(:,4,6) = - diff_min_image(at, iBo, iAh2) / rOH(6) ! d r_{OH_6} / d r_4
          drOH(:,3,6) = - drOH(:,4,6)                            ! d r_{OH_6} / d r_3

          drOH(:,4,7) = - diff_min_image(at, iBo, iBh1) / rOH(7) ! d r_{OH_7} / d r_4
          drOH(:,5,7) = - drOH(:,4,7)                            ! d r_{OH_7} / d r_5

          drOH(:,4,8) = - diff_min_image(at, iBo, iBh2) / rOH(8) ! d r_{OH_8} / d r_4
          drOH(:,6,8) = - drOH(:,4,8)                            ! d r_{OH_8} / d r_6


          drHH(:,2,1) = - diff_min_image(at, iAh1, iAh2) / rHH(1) ! d r_{HH_1} / d r_2
          drHH(:,3,1) = - drHH(:,2,1)                             ! d r_{HH_1} / d r_3

          drHH(:,2,2) = - diff_min_image(at, iAh1, iBh1) / rHH(2) ! d r_{HH_2} / d r_2
          drHH(:,5,2) = - drHH(:,2,2)                             ! d r_{HH_2} / d r_5

          drHH(:,2,3) = - diff_min_image(at, iAh1, iBh2) / rHH(3) ! d r_{HH_3} / d r_2
          drHH(:,6,3) = - drHH(:,2,3)                             ! d r_{HH_3} / d r_6

          drHH(:,3,4) = - diff_min_image(at, iAh2, iBh1) / rHH(4) ! d r_{HH_4} / d r_3
          drHH(:,5,4) = - drHH(:,3,4)                             ! d r_{HH_4} / d r_5

          drHH(:,3,5) = - diff_min_image(at, iAh2, iBh2) / rHH(5) ! d r_{HH_5} / d r_3
          drHH(:,6,5) = - drHH(:,3,5)                             ! d r_{HH_5} / d r_6

          drHH(:,5,6) = - diff_min_image(at, iBh1, iBh2) / rHH(6) ! d r_{HH_6} / d r_5
          drHH(:,6,6) = - drHH(:,5,6)                             ! d r_{HH_6} / d r_6
       endif



       if ( DO_FOURIER .eqv. .false. ) then

          ! Fit gaussians to the distance distribution functions

          call system_abort("Gaussian fit to water dimer distances is currently broken because the WATER_DIMER_D parameter is wrongly hardwired")

          if(present(vec) .or. present(dvec)) then
             fOH = 0.0_dp
             dfOH = 0.0_dp
             do i = 1, 8     ! basis function
                do j = 1, 8  ! atom pair
                   arg_r = rOH(j) - r0_OH(i)
                   fOH_ij = exp( - 0.5_dp * arg_r**2 * sigma_OH(i) )
                   if(present(vec)) fOH(i) = fOH(i) + fOH_ij
                   if(present(dvec)) then
                      do k = 1, 6
                         dfOH(:,k,i) = dfOH(:,k,i) - arg_r * sigma_OH(i) * fOH_ij * drOH(:,k,j)
                      enddo
                   endif
                enddo
             enddo
             
             fHH = 0.0_dp
             dfHH = 0.0_dp
             do i = 1, 6     ! basis function
                do j = 1, 6  ! atom pair
                   arg_r = rHH(j) - r0_HH(i)
                   fHH_ij = exp( - 0.5_dp * arg_r**2 * sigma_HH(i) )
                   if(present(vec)) fHH(i) = fHH(i) + fHH_ij
                   if(present(dvec)) then
                      do k = 1, 6
                         dfHH(:,k,i) = dfHH(:,k,i) - arg_r * sigma_HH(i) * fHH_ij * drHH(:,k,j)
                      enddo
                   endif
                enddo
             enddo
          endif
          
       else

          ! "Fourier Transform" distances. This should be completely reversible,
          ! i.e. we can reconstruct all the distances up to permutations.

          if(present(vec) .or. present(dvec)) then
             fOH = 0.0_dp
             dfOH = 0.0_dp
             do i = 1, FOURIER_N_COMP
                arg = PI*i/cutoff
                do j = 1, FOURIER_N_COMP
                   arg_r = arg * rOH(j)
                   if(present(vec)) fOH(i) = fOH(i) + cos( arg_r )
                   if(present(dvec)) then
                      do k = 1, 6
                         dfOH(:,k,i) = dfOH(:,k,i) - sin( arg_r ) * arg * drOH(:,k,j)
                      enddo
                   endif
                enddo
             enddo
             
             fHH = 0.0_dp
             dfHH = 0.0_dp
             do i = 1, FOURIER_N_COMP
                arg = PI*i/cutoff
                do j = 1, FOURIER_N_COMP
                   arg_r = arg * rHH(j)
                   if(present(vec)) fHH(i) = fHH(i) + cos( arg_r )
                   if(present(dvec)) then
                      do k = 1, 6
                         dfHH(:,k,i) = dfHH(:,k,i) - sin( arg_r ) * arg * drHH(:,k,j)
                      enddo
                   endif
                enddo
             enddo
          endif
          
       endif

       if (present(vec)) then
          vec(1) = distance_min_image(at, iAo, iBo)
          vec(2:FOURIER_N_COMP+1) = fOH
          vec(FOURIER_N_COMP+2:2*FOURIER_N_COMP+1) = fHH
       endif

       if (present(dvec)) then
          dvec = 0.0_dp

          dvec(:,1,1) = - diff_min_image(at, iAo, iBo) / distance_min_image(at, iAo, iBo)
          dvec(:,4,1) = - dvec(:,1,1)

          dvec(:,:,2:FOURIER_N_COMP+1)   = dfOH
          dvec(:,:,FOURIER_N_COMP+2:2*FOURIER_N_COMP+1) = dfHH
       endif

       !! O--H vectors
       !vA1 = diff_min_image(at,iAo,iAh1)
       !vA2 = diff_min_image(at,iAo,iAh2)
       !vB1 = diff_min_image(at,iBo,iBh1)
       !vB2 = diff_min_image(at,iBo,iBh2)
       !vAB = diff_min_image(at,iAo,iBo)

       !rA1 = norm(vA1)
       !rA2 = norm(vA2)
       !rB1 = norm(vB1)
       !rB2 = norm(vB2)
       !nA = vA1 .cross. vA2
       !nB = vB1 .cross. vB2
       !sA = (vA1+vA2)
       !dA = (vA1-vA2)
       !sB = (vB1+vB2)
       !dB = (vB1-vB2)
       !! descriptors
       !!vec(1) = normsq(sA)
       !!vec(2) = normsq(dA)
       !!vec(3) = (sA .dot. dA)**2
       !!vec(4) = normsq(sB)
       !!vec(5) = normsq(dB)
       !!vec(6) = (sB .dot. dB)**2
       !v(1) = rA1+rA2
       !v(2) = (rA1-rA2)**2
       !v(3) = vA1 .dot. vA2
       !v(4) = rB1+rB2
       !v(5) = (rB1-rB2)**2
       !v(6) = vB1 .dot. vB2

       !vec(1) = v(1)+v(4)
       !vec(2) = (v(1)-v(4))**2
       !vec(3) = v(2)+v(5)
       !vec(4) = (v(2)-v(5))**2
       !vec(5) = v(3)+v(6)
       !vec(6) = (v(3)-v(6))**2

       !vec(7) = sA .dot. sB

       !!v(8) = sA .dot. nB
       !!v(9) = sB .dot. nA
       !!v(10) = pA .dot. sB
       !!v(11) = pB .dot. sA

       !vec(8) = (dA .dot. dB)**2
       !vec(9) = (nA .dot. nB)**2
       !vec(10) = (vAB .dot. sB) - (vAB .dot. sA)
       !vec(11) = (vAB .dot. dB)**2 + (vAB .dot. dA)**2

       !vec(12) = distance_min_image(at, iAo, iBo)
    end subroutine water_dimer

    !#################################################################################
    !#
    !# Auxiliary functions start here.
    !#
    !#################################################################################

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

    !#################################################################################
    !#
    !% Gaussian: $ \exp( - \sigma(r-r_0)^2 )
    !#
    !#################################################################################

    function gaussian(r,r0,sigma)

      real(dp), intent(in) :: r             !% at r
      real(dp), intent(in) :: r0            !% centre of gaussian
      real(dp), intent(in) :: sigma         !% width of gaussian
      real(dp)             :: gaussian

      !gaussian = exp( -r**2 / (2.0_dp * sigma**2) ) / sqrt( 2.0_dp*PI ) / sigma
      gaussian = exp( -sigma * (r-r0)**2 )

    endfunction gaussian

    !#################################################################################
    !#
    !% Fermi-Dirac function: $ \frac{1}{ \exp( \frac{r-r_0}{\sigma} ) + 1 )
    !#
    !#################################################################################

    function fermi_dirac(r,r0,sigma)

      real(dp), intent(in) :: r             !% at r
      real(dp), intent(in) :: r0            !% centre of Fermi-Dirac
      real(dp), intent(in) :: sigma         !% width of Fermi-Dirac
      real(dp)             :: fermi_dirac

      fermi_dirac = 1.0_dp/( exp( (r-r0)/sigma ) + 1.0_dp )

    endfunction fermi_dirac

    !#################################################################################
    !#
    !% Fermi-Dirac function for cutoff left:
    !% $ \frac{\frac{r-r_0}{\sigma}}{ \exp( \frac{r-r_0}{\sigma} ) + 1 )
    !#
    !#################################################################################

    function cutoff_left(r,r0,w)

      real(dp), intent(in) :: r             !% at r
      real(dp), intent(in) :: r0            !% centre of Fermi-Dirac
      real(dp), intent(in) :: w         !% width of Fermi-Dirac
      real(dp)             :: cutoff_left
      real(dp)             :: arg

      arg = exp( (r-r0)/w )
      cutoff_left = arg / (arg+1.0_dp)

!      cutoff with a cos
!      r0w = r0 + w
!      if( r < r0 ) then
!          cutoff_left = 0.0_dp
!      elseif( r < r0w ) then
!          cutoff_left = 0.5_dp * cos( (r-r0w)*PI/w ) + 0.5_dp
!      else
!          cutoff_left = 1.0_dp
!      endif

    endfunction cutoff_left

    !#################################################################################
    !#
    !% Fermi-Dirac function for cutoff right:
    !% $ \frac{1}{ \exp( \frac{r-r_0}{\sigma} ) + 1 )
    !#
    !#################################################################################

    function cutoff_right(r,r0,w)

      real(dp), intent(in) :: r             !% at r
      real(dp), intent(in) :: r0            !% centre of Fermi-Dirac
      real(dp), intent(in) :: w         !% width of Fermi-Dirac
      real(dp)             :: cutoff_right

      cutoff_right = 1.0_dp / ( exp( (r-r0)/w ) + 1.0_dp )

!      cutoff with a cos
!      r0w = r0 - w
!      if( r < r0w ) then
!          cutoff_right = 1.0_dp
!      elseif( r < r0 ) then
!          cutoff_right = 0.5_dp * cos( (r-r0w)*PI/w ) + 0.5_dp
!      else
!          cutoff_right = 0.0_dp
!      endif

    endfunction cutoff_right


    function coff(r,r_cut)

       real(dp)             :: coff
       real(dp), intent(in) :: r, r_cut
       real(dp), parameter :: S = 0.25_dp

       if( r > r_cut ) then
           coff = 0.0_dp
       else
           coff = 0.5_dp * ( cos(PI*r/r_cut) + 1.0_dp )
       endif
       !if( r > r_cut ) then
       !    coff = 0.0_dp
       !elseif( r > (r_cut-S) ) then
       !    coff = 0.5_dp * ( cos(PI*(r-r_cut+S)/S) + 1.0_dp )
       !else
       !    coff = 1.0_dp
       !endif

    endfunction coff

    function dcoff(r,r_cut)

       real(dp)             :: dcoff
       real(dp), intent(in) :: r, r_cut
       real(dp), parameter :: S = 0.25_dp

       if( r > r_cut ) then
           dcoff = 0.0_dp
       else
           dcoff = - 0.5_dp * PI * sin(PI*r/r_cut) / r_cut
       endif
       !if( r > r_cut ) then
       !    dcoff = 0.0_dp
       !elseif( r > (r_cut-S) ) then
       !    dcoff = - 0.5_dp * PI * sin(PI*(r-r_cut+S)/S) / S
       !else
       !    dcoff = 0.0_dp
       !endif

    endfunction dcoff
    !#################################################################################
    !#
    !% Converts xyz coordinates to spherical
    !#
    !#################################################################################

    function xyz2spherical(xyz) result(spherical)

      real(dp), dimension(3) :: spherical            ! Spherical coordinates: $r$, $\theta$, $\phi$
      real(dp), dimension(3), intent(in) :: xyz      ! Cartesian coordinates

      spherical(1) = sqrt( sum( xyz**2 ) )

      if( spherical(1) == 0.0_dp ) then
          spherical(2) = 0.0_dp
          spherical(3) = 0.0_dp
          return
      endif

      spherical(2) = acos( xyz(3) / spherical(1) )

      if( xyz(1) /= 0.0_dp ) then
         spherical(3) = atan( xyz(2) / xyz(1) )
      elseif( xyz(2) > 0.0_dp ) then
         spherical(3) = PI / 2.0_dp
         return
      elseif( xyz(2) < 0.0_dp ) then
         spherical(3) = -PI / 2.0_dp
         return
      else
         spherical(3) = 0.0_dp
         return
      endif

      if( xyz(1) < 0.0_dp ) then
         spherical(3) = PI + spherical(3)
      endif

    endfunction xyz2spherical

    !#################################################################################
    !#
    !% Legendre polynomial
    !#
    !#################################################################################

    function LegendreP(l,m,x)

      real(dp), intent(in) :: x
      integer, intent(in)  :: l, m

      real(dp) :: LegendreP

      real(dp) :: pre_factor, m_over_two, sub_result_1, sub_result_2, sub_result_3
      integer  :: m_in, ll

      if( m>l .or. abs(x) > 1.0_dp ) then
          LegendreP = 0.0_dp
          return
      endif

      if( m>=0 ) then
          m_in = m
          pre_factor = 1.0_dp
      else
          m_in = -m
          pre_factor = oscillate(m) * factorial(l-m_in)/factorial(l+m_in)
      endif

      m_over_two = real(m_in,dp) / 2.0_dp

      sub_result_1 = oscillate(m_in) * factorial2(2*m_in-1) * (1.0_dp-x**2)**m_over_two
      sub_result_2 = x * (2.0_dp*m_in+1.0_dp) * sub_result_1

      if(l==m_in) then
         LegendreP = pre_factor * sub_result_1
         return
      elseif( l==m_in+1 ) then
         LegendreP = pre_factor * sub_result_2
         return
      endif

      do ll = m_in+2, l

         sub_result_3 = (x*(2.0_dp*ll-1)*sub_result_2 - (ll+m_in-1)*sub_result_1)/(ll-m_in)
         sub_result_1 = sub_result_2
         sub_result_2 = sub_result_3

      enddo

      LegendreP = pre_factor * sub_result_3

    endfunction LegendreP

    !#################################################################################
    !#
    !% Spherical Harmonic function
    !#
    !#################################################################################

    function SphericalY(l,m,theta,phi)

       integer,  intent(in) :: l, m
       real(dp), intent(in) :: theta, phi

       complex(dp) :: SphericalY
       complex(dp) :: arg

       arg = cplx_imag * m * phi

       SphericalY = sqrt( (2.0_dp*l+1) * factorial(l-m) / (4.0_dp*PI*factorial(l+m)) ) * &
                  & LegendreP(l,m,cos(theta)) * exp( arg )

    endfunction SphericalY

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

    !#################################################################################
    !#
    !% Radial functions used by fourier so3
    !#
    !#################################################################################

    function RadialFunction(cutoff, x, cutoff_f, cutoff_r1)

      real(dp) :: RadialFunction
      real(dp), intent(in) :: cutoff, x(3), cutoff_r1
      integer, intent(in) :: cutoff_f
      real(dp) :: r

      r = norm(x)

      if (cutoff_f == 1) then
         if ((r >= 0.0_dp) .and. (r < (cutoff - cutoff_r1))) then
            RadialFunction = 1.0_dp
         else if ((r >= (cutoff - cutoff_r1)) .and. (r <= cutoff)) then
            RadialFunction = (cos(0.5_dp * PI * (r - cutoff + cutoff_r1) / cutoff_r1))**2
         else
            call system_abort('RadialFunction: distance greater than cutoff')
         end if
      else if (cutoff_f == 2) then
         if ((r >= 0.0_dp) .and. (r < (cutoff - (2.0_dp * cutoff_r1)))) then
            RadialFunction = 0.0_dp
         else if ((r >= (cutoff - (2.0_dp * cutoff_r1))) .and. (r <= cutoff)) then
            RadialFunction = (cos(0.5_dp * PI * (r - cutoff + cutoff_r1) / cutoff_r1))**2
         else
            call system_abort('RadialFunction: distance greater than cutoff')
         end if
      else if (cutoff_f == 3) then
         if ((r >= 0.0_dp) .and. (r < (cutoff - cutoff_r1))) then
            RadialFunction = (cos(0.5_dp * PI * (r - cutoff + cutoff_r1) / (cutoff - cutoff_r1)))**2
         else if ((r >= (cutoff - cutoff_r1)) .and. (r <= cutoff)) then
            RadialFunction = (cos(0.5_dp * PI * (r - cutoff + cutoff_r1) / cutoff_r1))**2
         else
            call system_abort('RadialFunction: distance greater than cutoff')
         end if
      else
         call system_abort('RadialFunction: radial function type unknown')
      end if

    end function RadialFunction

    !#################################################################################
    !#
    !% Derivative of Radial functions used by fourier so3
    !#
    !#################################################################################

    function GradRadialFunction(cutoff, x, cutoff_f, cutoff_r1)
 
      real(dp) :: GradRadialFunction(3)
      real(dp), intent(in) :: cutoff, x(3), cutoff_r1
      integer, intent(in) :: cutoff_f
      real(dp) :: r

      r = norm(x)

      if (cutoff_f == 1) then
         if ((r >= 0.0_dp) .and. (r < (cutoff - cutoff_r1))) then
            GradRadialFunction = 0.0_dp
         else if ((r >= (cutoff - cutoff_r1)) .and. (r <= cutoff)) then
            GradRadialFunction = -PI * x * cos(0.5_dp * PI * (r - cutoff + cutoff_r1) / cutoff_r1) &
                                         * sin(0.5_dp * PI * (r - cutoff + cutoff_r1) / cutoff_r1) &
                                         / (cutoff_r1 * r)
         else
            call system_abort('GradRadialFunction: distance greater than cutoff')
         end if
      else if (cutoff_f == 2) then
         if ((r >= 0.0_dp) .and. (r < (cutoff - (2.0_dp * cutoff_r1)))) then
            GradRadialFunction = 0.0_dp
         else if ((r >= (cutoff - (2.0_dp * cutoff_r1))) .and. (r <= cutoff)) then
            GradRadialFunction = -PI * x * cos(0.5_dp * PI * (r - cutoff + cutoff_r1) / cutoff_r1) &
                                         * sin(0.5_dp * PI * (r - cutoff + cutoff_r1) / cutoff_r1) &
                                         / (cutoff_r1 * r)
         else
            call system_abort('GradRadialFunction: distance greater than cutoff')
         end if
      else if (cutoff_f == 3) then
         if ((r >= 0.0_dp) .and. (r < (cutoff - cutoff_r1))) then
            GradRadialFunction = -PI * x * cos(0.5_dp * PI * (r - cutoff + cutoff_r1) / (cutoff - cutoff_r1)) &
                                         * sin(0.5_dp * PI * (r - cutoff + cutoff_r1) / (cutoff - cutoff_r1)) &
                                         / ((cutoff - cutoff_r1) * r)
         else if ((r >= (cutoff - cutoff_r1)) .and. (r <= cutoff)) then
            GradRadialFunction = -PI * x * cos(0.5_dp * PI * (r - cutoff + cutoff_r1) / cutoff_r1) &
                                         * sin(0.5_dp * PI * (r - cutoff + cutoff_r1) / cutoff_r1) &
                                         / (cutoff_r1 * r)
         else
            call system_abort('GradRadialFunction: distance greater than cutoff')
         end if
      else
         call system_abort('GradRadialFunction: radial function type unknown')
      end if

    end function GradRadialFunction

    function wigner_big_U(j,m1,m2,omega,theta,phi,denom)
       complex(dp) :: wigner_big_U
       integer, intent(in) :: j,m1,m2
       real(dp), intent(in) :: omega,theta,phi
       integer, intent(in), optional :: denom

       complex(dp) :: v, u, tmp2
       real(dp) :: tmp1
       integer :: s, my_denom

       my_denom = optional_default(1,denom)

       v = sin(omega/2.0_dp)*sin(theta)*CPLX_ONE
       u = cos(omega/2.0_dp)*CPLX_ONE -sin(omega/2.0_dp)*cos(theta)*CPLX_IMAG
  
       tmp1 = 0.0_dp
  
       if( m1+m2 >= 0 ) then
          tmp2 = ( u / ( -CPLX_IMAG*v ) ) ** ((m1+m2)/my_denom)
          do s = max(0,-(m1+m2))/my_denom, min((j-m1),(j-m2))/my_denom
             tmp1 = tmp1 + (1.0_dp - 1.0_dp/v**2)**s / ( factorial(s) * &
             & factorial(s+(m1+m2)/my_denom) * factorial((j-m1)/my_denom-s) * factorial((j-m2)/my_denom-s) )
          enddo
       else
          tmp2 = ( conjg(u) / ( -CPLX_IMAG*v ) ) ** (-(m1+m2)/my_denom)
          do s = max(0,(m1+m2))/my_denom, min((j+m1),(j+m2))/my_denom
             tmp1 = tmp1 + (1.0_dp - 1.0_dp/v**2)**s / ( factorial(s) * &
             & factorial(s-(m1+m2)/my_denom) * factorial((j+m1)/my_denom-s) * factorial((j+m2)/my_denom-s) )
          enddo
       endif
  
       wigner_big_U = (-CPLX_IMAG*v)**j * tmp2 * exp( -CPLX_IMAG*((m1-m2)/my_denom)*phi ) * &
       & sqrt( factorial((j+m1)/my_denom) * factorial((j-m1)/my_denom) * factorial((j+m2)/my_denom) * &
       & factorial((j-m2)/my_denom) ) * tmp1

    endfunction wigner_big_U

    function grad_wigner_big_U(j,m1,m2,omega,theta,phi,denom)
       complex(dp), dimension(3) :: grad_wigner_big_U
       integer, intent(in) :: j,m1,m2
       real(dp), intent(in) :: omega,theta,phi
       integer, intent(in), optional :: denom

       complex(dp) :: u, tmp2, tmp5, tmp6, du_domega, du_dtheta, dU_dv, dU_du
       real(dp) :: tmp1, tmp3, tmp4, tmp7, v, dv_domega, dv_dtheta
       integer :: s, my_denom, my_2j

       my_denom = optional_default(1,denom)
       if( my_denom == 1 ) then
          my_2j = 2*j
       else
          my_2j = j
       endif

       v = sin(omega/2.0_dp)*sin(theta)
       u = cos(omega/2.0_dp)*CPLX_ONE - sin(omega/2.0_dp)*cos(theta)*CPLX_IMAG
  
       tmp1 = 0.0_dp
       tmp4 = 0.0_dp
  
       if( m1+m2 >= 0 ) then
          tmp2 = ( u / ( -CPLX_IMAG*v ) ) ** ((m1+m2)/my_denom)
          tmp5 = CPLX_IMAG*((m1+m2)/my_denom)*( u / ( -CPLX_IMAG*v ) ) ** ((m1+m2)/my_denom) / ( -CPLX_IMAG*v )
          tmp6 = ((m1+m2)/my_denom)*( u / ( -CPLX_IMAG*v ) ) ** ((m1+m2)/my_denom-1) / ( -CPLX_IMAG*v )
          
          du_domega = -sin(omega/2.0_dp)*CPLX_ONE/2.0_dp - cos(omega/2.0_dp)*cos(theta)*CPLX_IMAG / 2.0_dp
          du_dtheta = sin(omega/2.0_dp)*sin(theta)*CPLX_IMAG

          do s = max(0,-(m1+m2))/my_denom, min((j-m1),(j-m2))/my_denom
             tmp3 = 1.0_dp / ( factorial(s) * &
             & factorial(s+(m1+m2)/my_denom) * factorial((j-m1)/my_denom-s) * factorial((j-m2)/my_denom-s) )
             tmp1 = tmp1 + (1.0_dp - 1.0_dp/v**2)**s * tmp3
             tmp4 = tmp4 + 2.0_dp * s*(1.0_dp - 1.0_dp/v**2)**(s-1) / v**3 * tmp3
          enddo
       else
          tmp2 = ( conjg(u) / ( -CPLX_IMAG*v ) ) ** (-(m1+m2)/my_denom)
          tmp5 = -CPLX_IMAG*((m1+m2)/my_denom)*( conjg(u) / ( -CPLX_IMAG*v ) ) ** (-(m1+m2)/my_denom) / ( -CPLX_IMAG*v )
          tmp6 = (-(m1+m2)/my_denom)*( conjg(u) / ( -CPLX_IMAG*v ) ) ** (-(m1+m2)/my_denom-1) / ( -CPLX_IMAG*v )

          du_domega = -sin(omega/2.0_dp)*CPLX_ONE/2.0_dp + cos(omega/2.0_dp)*cos(theta)*CPLX_IMAG / 2.0_dp
          du_dtheta = -sin(omega/2.0_dp)*sin(theta)*CPLX_IMAG

          do s = max(0,(m1+m2))/my_denom, min((j+m1),(j+m2))/my_denom
             tmp3 = 1.0_dp / ( factorial(s) * &
             & factorial(s-(m1+m2)/my_denom) * factorial((j+m1)/my_denom-s) * factorial((j+m2)/my_denom-s) )
             tmp1 = tmp1 + (1.0_dp - 1.0_dp/v**2)**s * tmp3
             tmp4 = tmp4 + 2.0_dp * s*(1.0_dp - 1.0_dp/v**2)**(s-1) / v**3 * tmp3
          enddo
       endif
       tmp7 = sqrt( factorial((j+m1)/my_denom) * factorial((j-m1)/my_denom) * &
       & factorial((j+m2)/my_denom) * factorial((j-m2)/my_denom) )

       dU_dv = ( -my_2j*CPLX_IMAG*(-CPLX_IMAG*v)**(my_2j-1) * tmp2 * tmp1 + &
       & (-CPLX_IMAG*v)**my_2j * (tmp5 * tmp1 + tmp2 * tmp4) ) * exp( -CPLX_IMAG*((m1-m2)/my_denom)*phi ) * tmp7
       dU_du = (-CPLX_IMAG*v)**my_2j * tmp6 * tmp1 * exp( -CPLX_IMAG*((m1-m2)/my_denom)*phi ) * tmp7

       dv_domega = cos(omega/2.0_dp)*sin(theta) / 2.0_dp
       dv_dtheta = sin(omega/2.0_dp)*cos(theta)

       
       grad_wigner_big_U(1) = dU_dv*dv_domega + dU_du*du_domega
       grad_wigner_big_U(2) = dU_dv*dv_dtheta + dU_du*du_dtheta
       grad_wigner_big_U(3) = -CPLX_IMAG*((m1-m2)/my_denom) * (-CPLX_IMAG*v)**my_2j &
       & * tmp2 * exp( -CPLX_IMAG*((m1-m2)/my_denom)*phi ) * tmp7 * tmp1
  
     endfunction grad_wigner_big_U

    !#################################################################################
    !#
    !% Initialise Clebsch-Gordan coefficient for fast calculation
    !#
    !#################################################################################

    !subroutine cg_initialise(j1,m1,j2,m2,j,m)
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

!       do i_m = -m_max, m_max
!       do i_j = abs(j1_max-j2_max), j1_max+j2_max
!       do i_m2 = -m2_max, m2_max
!       do i_j2 = 0, j2_max
!       do i_m1 = -m1_max, m1_max
!       do i_j1 = 0, j1_max
!
!          cg_array(i_j1,i_m1,i_j2,i_m2,i_j,i_m) = cg_calculate(i_j1,i_m1,i_j2,i_m2,i_j,i_m,denom)
!
!       enddo
!       enddo
!       enddo
!       enddo
!       enddo
!       enddo
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
        & abs(m1)<=cg_m1_max .and. abs(m2)<=cg_m2_max .and. abs(m) <= cg_m_max .and. &
        & cg_initialised ) then
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
                & (abs(m1)<=j1) .and. (abs(m2)<=j2) .and. (abs(m)<=j) &
                & .and. (m1+m2==m) .and. (j1+j2 >= j) .and. (abs(j1-j2) <= j) &
                & .and. (mod(j1+j2+j,my_denom)==0)
       
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

      cg = oscillate((m+j1-j2)/my_denom) * sqrt(2.0_dp*real(j,dp)/real(my_denom,dp)+1.0_dp) &
      & * wigner3j(j1,m1,j2,m2,j,-m,denom)

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

       tc = factorial((a+b-c)/my_denom) * factorial((a-b+c)/my_denom) &
       & * factorial((-a+b+c)/my_denom) / factorial((a+b+c)/my_denom+1)

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
                   & factorial((j1+m1)/my_denom) * factorial((j1-m1)/my_denom) * &
                   & factorial((j2+m2)/my_denom) * factorial((j2-m2)/my_denom) * &
                   & factorial((j+m)/my_denom) * factorial((j-m)/my_denom) )
                   
        sum_coeff = 0.0_dp

        kmin = max( j2-j-m1, j1+m2-j, 0 ) / my_denom
        kmax = min( j1+j2-j, j1-m1, j2+m2) / my_denom

        do k = kmin, kmax

           sum_term = 1.0_dp / ( factorial(k) * factorial((j-j2+m1)/my_denom+k) * &
                    & factorial((j-j1-m2)/my_denom+k) * factorial((j1+j2-j)/my_denom-k)    * &
                    & factorial((j1-m1)/my_denom-k) * factorial((j2+m2)/my_denom-k) )

           sum_term = oscillate(k) * sum_term
           
           sum_coeff = sum_coeff + sum_term

        enddo

        wigner3j = pre_fac * triang_coeff * main_coeff * sum_coeff

    endfunction wigner3j

    !#################################################################################
    !#
    !% Variance of an array.
    !#
    !#################################################################################

    function p_norm(vec,p)

       real(dp) :: p_norm
       real(dp), dimension(:) :: vec
       integer, optional :: p
       integer :: my_p

       my_p = optional_default(2,p)

       p_norm = ( sum(vec**my_p) / size(vec) ) ** (1.0_dp/real(my_p,dp)) 

    endfunction p_norm

endmodule descriptors_module
