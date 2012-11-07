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

module gp_teach_module

   use iso_c_binding, only : C_NULL_CHAR
   use libatoms_module
   use gp_predict_module
   use clustering_module

   interface gp_sparsify
      module procedure gpFull_sparsify_array_config_type
   endinterface gp_sparsify
   public :: gp_sparsify

   contains

   subroutine gpCoordinates_sparsify_config_type(this, n_sparseX, default_all, sparseMethod, error)
      type(gpCoordinates), intent(inout) :: this
      integer, dimension(:), intent(in) :: n_sparseX
      logical, intent(in) :: default_all
      integer, optional, intent(in) :: sparseMethod
      integer, optional, intent(out) :: error

      integer :: my_sparseMethod, li, ui, i_config_type, n_config_type, d, n_x
      integer, dimension(:), allocatable :: config_type_index, sparseX_index
      real(dp), dimension(:,:), allocatable :: x
      integer, dimension(:), allocatable :: x_index
      logical, dimension(:), allocatable :: x_unique

      INIT_ERROR(error)

      my_sparseMethod = optional_default(GP_SPARSE_RANDOM,sparseMethod)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpCoordinates_sparsify: : object not initialised',error)
      endif

      d = size(this%x,1)
      n_x = size(this%x,2)

      if(my_sparseMethod == GP_SPARSE_UNIQ) then

         allocate(x(d,n_x))
         allocate(x_index(n_x))
         allocate(x_unique(n_x))

         x = this%x
         x_index = (/(i,i=1,n_x)/)

         call heap_sort(x,i_data=x_index)
         call uniq(x,unique=x_unique)
         this%n_sparseX = count(x_unique)

         call print('UNIQ type sparsification specified. The number of sparse point was changed to '//this%n_sparseX//' from '//n_sparseX//'.')
      else
         this%n_sparseX = sum(n_sparseX)
      endif

      call reallocate(this%sparseX, this%d,this%n_sparseX, zero = .true.)

      call reallocate(this%sparseX_index, this%n_sparseX, zero = .true.)
      call reallocate(this%map_sparseX_globalSparseX, this%n_sparseX, zero = .true.)
      call reallocate(this%alpha, this%n_sparseX, zero = .true.)
      call reallocate(this%sparseCutoff, this%n_sparseX, zero = .true.)
      this%sparseCutoff = 1.0_dp

      ui = 0
      do i_config_type = 1, size(n_sparseX)
         if(default_all) then
            allocate(config_type_index(n_x), sparseX_index(this%n_sparsex))
            config_type_index = (/(i,i=1,n_x)/)
            li = 1
            ui = this%n_sparseX
            n_config_type = this%n_sparsex
         else
            if( n_sparseX(i_config_type) == 0 ) cycle

            n_config_type = count(i_config_type == this%config_type)
            allocate(config_type_index(n_config_type),sparseX_index(n_sparseX(i_config_type)))
            config_type_index = find(i_config_type == this%config_type)

            li = ui + 1
            ui = ui + n_sparseX(i_config_type)
         endif
         
         select case(my_sparseMethod)
         case(GP_SPARSE_RANDOM)
            call fill_random_integer(sparseX_index, n_config_type)
         case(GP_SPARSE_PIVOT)
            call pivot(this%x(:,config_type_index), sparseX_index, theta = this%theta)
         case(GP_SPARSE_CLUSTER)
            call bisect_kmedoids(this%x(:,config_type_index), n_sparseX(i_config_type), med = sparseX_index, theta = this%theta)
         case(GP_SPARSE_UNIFORM)
            call select_uniform(this%x(:,config_type_index), sparseX_index)
         case(GP_SPARSE_KMEANS)
            call cluster_kmeans(this%x(:,config_type_index), sparseX_index, theta = this%theta)
         case(GP_SPARSE_COVARIANCE)
            call sparse_covariance(this,sparseX_index,config_type_index)
         case(GP_SPARSE_UNIQ)
            exit
         case default
            RAISE_ERROR('gpCoordinates_sparsify: '//my_sparseMethod//' method is unknown', error)
         endselect
         this%sparseX_index(li:ui) = config_type_index(sparseX_index)
         deallocate(config_type_index,sparseX_index)

         if(default_all) exit
      enddo

      if(my_sparseMethod == GP_SPARSE_UNIQ) this%sparseX_index = x_index(find_indices(x_unique))

      call sort_array(this%sparseX_index)
      if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
         if(allocated(this%sparseX)) deallocate(this%sparseX)
         allocate(this%sparseX(maxval(this%x_size(this%sparseX_index)),this%n_sparseX))
         if(allocated(this%sparseX_size)) deallocate(this%sparseX_size)
         allocate(this%sparseX_size(this%n_sparseX))
         this%sparseX(:,:) = this%x(1:maxval(this%x_size(this%sparseX_index)),this%sparseX_index)
         this%sparseX_size = this%x_size(this%sparseX_index)
      else
         this%sparseX(:,:) = this%x(:,this%sparseX_index)
      endif

      if(allocated(this%covarianceDiag_sparseX_sparseX)) deallocate(this%covarianceDiag_sparseX_sparseX)
      allocate(this%covarianceDiag_sparseX_sparseX(this%n_sparseX))
      this%covarianceDiag_sparseX_sparseX = this%covarianceDiag_x_x(this%sparseX_index)

      this%sparseCutoff = this%cutoff(this%sparseX_index)

      if(allocated(x)) deallocate(x)
      if(allocated(x_index)) deallocate(x_index)
      if(allocated(x_unique)) deallocate(x_unique)

      this%sparsified = .true.

   endsubroutine gpCoordinates_sparsify_config_type

   subroutine gpFull_sparsify_array_config_type(this, n_sparseX, default_all, sparseMethod, error)
      type(gpFull), intent(inout) :: this
      integer, dimension(:,:), intent(in) :: n_sparseX
      logical, dimension(:), intent(in) :: default_all
      integer, dimension(:), optional, intent(in) :: sparseMethod
      integer, optional, intent(out) :: error

      integer :: i
      integer, dimension(:), allocatable :: my_sparseMethod

      INIT_ERROR(error)

      if( .not. this%initialised ) then
         RAISE_ERROR('gpFull_sparsify_array: object not initialised',error)
      endif

      allocate(my_sparseMethod(this%n_coordinate))
      my_sparseMethod = optional_default((/ (GP_SPARSE_RANDOM, i=1,this%n_coordinate) /),sparseMethod)

      do i = 1, this%n_coordinate
         call gpCoordinates_sparsify_config_type(this%coordinate(i),n_sparseX(:,i), default_all(i), my_sparseMethod(i), error)
      enddo

      if(allocated(my_sparseMethod)) deallocate(my_sparseMethod)

   endsubroutine gpFull_sparsify_array_config_type

   subroutine sparse_covariance(this, index_out, config_type_index)
      type(gpCoordinates), intent(in) :: this
      integer, dimension(:), intent(out) :: index_out
      integer, dimension(:), intent(in), optional :: config_type_index
 
      real(dp), dimension(:), allocatable :: score, k_n, xI_xJ
      real(dp), dimension(:,:), allocatable :: k_mn, k_mm_k_m
      real(dp), dimension(1,1) :: k_mm
      integer :: m, n, i, ii, j, jj, i_p
      integer, dimension(1) :: j_loc
 
      type(LA_Matrix) :: LA_k_mm
 
      call system_timer('sparse_covariance')
      if(present(config_type_index)) then
         n = size(config_type_index)
      else
         n = size(this%x,2)
      endif
      m = size(index_out)
 
      allocate(k_n(n), k_mn(m,n), score(n), k_mm_k_m(m,n))
      k_mn = 0.0_dp
 
      allocate(xI_xJ(this%d))
 
      j = 1
      index_out(j) = 1 !ceiling(ran_uniform() * n)
 
      k_mm = 1.0_dp+1.0e-6_dp
      call initialise(LA_k_mm,k_mm)
 
      do j = 1, m-1
 
         if(present(config_type_index)) then
            jj = config_type_index(index_out(j))
         else
            jj = index_out(j)
         endif

!$omp parallel do default(none) shared(n,this,k_mn,jj,j,k_n,LA_k_mm,k_mm_k_m,score,config_type_index, index_out) private(i,i_p,ii)
         do i = 1, n

            if(present(config_type_index)) then
               ii = config_type_index(i)
            else
               ii = i
            endif
            if(this%covariance_type == COVARIANCE_BOND_REAL_SPACE) then
            elseif(this%covariance_type == COVARIANCE_DOT_PRODUCT) then
               k_mn(j,i) = dot_product( this%x(:,ii), this%x(:,jj) )**this%theta(1)
            elseif( this%covariance_type == COVARIANCE_ARD_SE ) then
               k_mn(j,i) = 0.0_dp
               do i_p = 1, this%n_permutations
                  !xI_xJ = (this%x(this%permutations(:,i_p),i) - this%x(:,j)) / 4.0_dp
                  k_mn(j,i) = k_mn(j,i) + exp( -0.5_dp * sum((this%x(this%permutations(:,i_p),ii) - this%x(:,jj))**2) / 16.0_dp )
               enddo
            endif
 
            call Matrix_Solve(LA_k_mm,k_mn(1:j,i),k_mm_k_m(1:j,i))
            score(i) = dot_product( k_mn(1:j,i), k_mm_k_m(1:j,i) )
         enddo
 
         j_loc = minloc(score)
         jj = j_loc(1)
         index_out(j+1) = jj

         if(j == 1) then
            call print('Initial score: '//score)
            call print('Min score: '//minval(score))
         endif
 
         !k_mm(1:j_i,j_i+1) = k_mn(1:j_i,j)
         !k_mm(j_i+1,1:j_i) = k_mn(1:j_i,j)
         !k_mm(j_i+1,j_i+1) = 1.0_dp
         call LA_Matrix_Expand_Symmetrically(LA_k_mm,(/k_mn(1:j,jj),1.0_dp+1.0e-6_dp/))
         !call initialise(LA_k_mm,k_mm(1:j_i+1,1:j_i+1))
 
      enddo
      call print('Final score: '//score)
      call print('Min score: '//minval(score))

      deallocate(k_n, k_mn, score, k_mm_k_m)
      if(allocated(xI_xJ)) deallocate(xI_xJ)
      call finalise(LA_k_mm)
      call system_timer('sparse_covariance')
 
   endsubroutine sparse_covariance

end module gp_teach_module
