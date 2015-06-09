! This Fortran95 Code/Module performs Machine Learning Force prediction based on XYZ structural Files;
! Sorting/Selecting Algorithms are used to construct sub training database during heavy machine learning task;
! Database can be updated according to an error_threshold, e.g. pred_err; 
! Internal Vector representation is used for describing the atomic structures; 
! Any issues or bugs regarding this code please contact Zhenwei Li (zhenwei.li@unibas.ch)
! Reference article : Molecular Dynamics with On-the-Fly Machine Learning of Quantum-Mechanical Forces, Zhenwei Li, James Kermode, Alessandro De Vita, PRL, 114, 096405 (2015) 
module force_machine_learning_module

  use libAtoms_module
  use potential_module
  implicit none

  type optimise_likelihood
      logical use_qr
      real(dp), dimension(:,:), pointer :: distance_matrix, covariance_matrix, force_covariance_matrix,  t
  endtype optimise_likelihood
 
public :: write_distance
interface  write_distance
   module procedure write_distance_file
   module procedure write_distance_file_real
end interface

  contains

  function likelihood(x,am_data)
      real(dp), dimension(:) :: x
      character(len=1), dimension(:), optional :: am_data
      real(dp) :: likelihood

      type(optimise_likelihood) :: am
      real(dp) :: sigma_error, sigma_covariance
      real(dp), dimension(:,:), allocatable :: ICF
      integer :: i
      type(LA_Matrix) :: LA_covariance

      am = transfer(am_data,am)
      sigma_error = x(1)
      sigma_covariance = x(2) 

      am%covariance_matrix = exp(-0.5_dp * am%distance_matrix/sigma_covariance**2)
      do i = 1, size(am%covariance_matrix,1)
         am%covariance_matrix(i,i) = am%covariance_matrix(i,i) + sigma_error**2
      enddo

      LA_covariance = am%covariance_matrix
      allocate(ICF(size(am%covariance_matrix,1),size(am%covariance_matrix,1)))

      if (am%use_qr) then
         call Matrix_QR_Solve(LA_covariance,am%force_covariance_matrix,ICF)
      else
         call Matrix_Solve(LA_covariance,am%force_covariance_matrix,ICF)
      end if

      likelihood = 1.5_dp * LA_Matrix_LogDet(LA_covariance) + 0.5_dp * trace(ICF)
!      write(*,*) "DEBUG Likelihood value: ", likelihood

      call finalise(LA_covariance)
      if(allocated(ICF)) deallocate(ICF)

   endfunction likelihood

   function dlikelihood(x,am_data)
      real(dp), dimension(:) :: x
      character(len=1), dimension(:), optional :: am_data
      real(dp), dimension(size(x)) :: dlikelihood

      type(optimise_likelihood) :: am
      real(dp) :: sigma_error, sigma_covariance
      real(dp), dimension(:,:), allocatable :: ICF, IC, ICDCDS, DCDS
      integer :: i, n
      type(LA_Matrix) :: LA_covariance

      am = transfer(am_data,am)
      sigma_error = x(1)
      sigma_covariance = x(2) 

      n = size(am%covariance_matrix,1)
      allocate(ICF(n,n),DCDS(n,n),IC(n,n),ICDCDS(n,n))

      am%covariance_matrix = exp(-0.5_dp * am%distance_matrix/sigma_covariance**2)
      DCDS = am%covariance_matrix * am%distance_matrix / sigma_covariance**3

      do i = 1, n
         am%covariance_matrix(i,i) = am%covariance_matrix(i,i) + sigma_error**2
      enddo
      
      LA_covariance = am%covariance_matrix

      if (am%use_qr) then
         call Matrix_QR_Solve(LA_covariance,am%force_covariance_matrix,ICF)
         call Matrix_QR_Solve(LA_covariance,DCDS,ICDCDS)
         call LA_Matrix_QR_Inverse(LA_covariance,IC)
      else
         call Matrix_Solve(LA_covariance,am%force_covariance_matrix,ICF)
         call Matrix_Solve(LA_covariance,DCDS,ICDCDS)
         call LA_Matrix_Inverse(LA_covariance,IC)
      end if

      dlikelihood(1) = sigma_error * (3*trace(IC) - sum(IC*ICF) )
      dlikelihood(2) = 1.5_dp * trace(ICDCDS) - 0.5_dp * sum(ICDCDS*transpose(ICF))
 !    write(*,*) "DEBUG dlikelihood", dlikelihood(1)**2, dlikelihood(2)**2
      
      call finalise(LA_covariance)
      if(allocated(ICF)) deallocate(ICF)
      if(allocated(DCDS)) deallocate(DCDS)
      if(allocated(IC)) deallocate(IC)
      if(allocated(ICDCDS)) deallocate(ICDCDS)

   endfunction dlikelihood

 function  weight_neighbour(r, r0, m, weight_index)

     real(dp), intent(in)                    :: r, r0, m
     integer, intent(in), optional           :: weight_index 
     integer                                 :: selector     
     real(dp)                                :: weight_neighbour, kernel

     if (present(weight_index)) then
           selector=weight_index
     else
         selector=0
     endif

     SELECT CASE (selector)
     CASE (0)
        weight_neighbour = exp(-((r/r0)**m)) /r
     CASE (1)
        weight_neighbour = exp(-((r/r0)**m)) 
     !CASE (2)
     !   kernel=sqrt(m*2_dp)*(r/r0)
     !   weight_neighbour= (1_dp/GAMMA(m) / (2_dp**(m-1)) ) * (kernel**m) * BESYN(int(m), kernel)  
     END SELECT

 end function weight_neighbour

 function ionic_charge(z)
     ! no partial charges 
     integer                           ::		  i, ionic_charge
     integer                           ::		  x(8)
     integer, allocatable              ::                 i_index(:)
     real(dp), allocatable             ::                 xx(:)
     integer, intent(in)               ::		  z
     !
     x=[0, 2, 10, 18, 36, 64, 86, 118]   
     allocate(xx(size(x)))
     allocate(i_index(size(x))) 
     do i=1, size(x)
		i_index(i) = i
		xx(i) = 1.0_dp*abs(x(i) - z)
     enddo
     call heap_sort(xx,i_data=i_index)
     ionic_charge = z- x(i_index(1))  
     !
     deallocate(xx, i_index)
 endfunction ionic_charge

function char_weight(z1, z2, do_atz_index)
  integer, intent(in)                    :: z1, z2
  integer, intent(in), optional          :: do_atz_index
  real(dp)                               :: char_weight
  integer                                :: selector

  if (present(do_atz_index) ) then
      selector=do_atz_index
  else
      selector=0
 endif
 
select case (selector)
case (0)
   char_weight=1.0_dp
case (1)
    char_weight=real(z1*z2, dp)
case (2) 
      char_weight=real(ionic_charge(z1)*ionic_charge(z2), dp)
end select
!
endfunction char_weight

!
function internal_vector(at, r0, m, species_number, n_center_atom, weight_index, do_atz_index)
      real(dp)                            :: internal_vector(3), delta_r(3), delta_r_len
      real(dp), intent(in)                :: r0, m
      type(Atoms), intent(in)             :: at
      integer, intent(in), optional       :: weight_index, do_atz_index
      integer, intent(in)                 :: species_number
      integer                             :: i, j, n_center_atom
      !
      internal_vector = 0.0_dp
      do i=1, n_neighbours(at, n_center_atom)
          j = neighbour(at, n_center_atom, i, distance=delta_r_len, diff=delta_r) 
          !
          if ((at%z(j) .eq. species_number) .or. (species_number.eq.0))  then
                 internal_vector = internal_vector + delta_r * weight_neighbour(delta_r_len, r0, m, weight_index) * char_weight(at%z(j), at%z(n_center_atom), do_atz_index=do_atz_index)
          endif
      enddo
      !
endfunction  internal_vector

 subroutine load_iv_params(data_dir, r_grid, m_grid, species, k, add_vector, n_data)

   ! reads r and m from iv_params file. The first line must contain the number
   ! of internal vector, and the following lines should be formatted in two columns 
   ! with commas as spacers (a common csv file)
   character(STRING_LENGTH), intent(in)                            :: data_dir                          
   real(dp), dimension(:), intent(out), allocatable                :: r_grid, m_grid
   integer, dimension(:), intent(out), allocatable                 :: species
   integer                                                         :: i, m
   integer, intent(out)                                            :: k, add_vector 
   integer, intent(out), optional                                  :: n_data 

   open (unit=22, file=trim(data_dir)//'grid.dat', status='old', action='read')
   read(22,*), k, add_vector, m 

   call print("Number of iv: "//k)
   if (present(n_data)) n_data=m
   if (.not. allocated(r_grid))  allocate(r_grid(k))
   if (.not. allocated(m_grid))  allocate(m_grid(k))
   if (.not. allocated(species)) allocate(species(k))

   call print(" Reading internal vectors parameters from file "//trim(data_dir)//"grid.dat")
   
   do i=1, k
      read(22,*) r_grid(i), m_grid(i), species(i)
      call print("Vector "//i//" : "//r_grid(i)//" "//m_grid(i)//" "//species(i))
   end do
   close(22)

 end subroutine load_iv_params

 subroutine write_iv_params(data_dir, r_grid, m_grid, species, k_in, add_vector, n_data)
  character(STRING_LENGTH), intent(in)                            :: data_dir 
  real(dp), dimension(:), intent(in)                              :: r_grid, m_grid
  integer,  dimension(:), intent(in)                              :: species
  integer                                                         :: i
  integer, intent(in)                                             :: k_in, n_data, add_vector

  open (unit=22, status='replace', file=trim(data_dir)//'grid.dat', form='FORMATTED')

   write(22, *) k_in,  add_vector, n_data
  do i=1, k_in
     write(22, *) r_grid(i), m_grid(i), species(i)
  enddo
  close(unit=22)
  
 end subroutine write_iv_params

 function cov_dij(sigma_cov, dist, func_type)  
    real(dp), intent(in)                  ::        sigma_cov, dist
    real(dp)                              ::        cov_dij, d_sq
    integer, intent(in), optional         ::        func_type
    integer                               ::        selector

    d_sq = dist

    if (present(func_type)) then
         selector=func_type
    else
         selector=0
    endif

    SELECT CASE (selector)
    CASE (0)
       cov_dij = exp(-0.5_dp*d_sq/sigma_cov**2)
    CASE (1)
       cov_dij = sqrt(d_sq)
    CASE (2)
       cov_dij = d_sq
    CASE (3)
       cov_dij = (sqrt(d_sq))**3
    CASE (4)
       cov_dij = d_sq**2
    END SELECT 

 end function cov_dij

 function cov(feature_matrix1, feature_matrix2, bent_space1, bent_space2, iv_weights, sigma_cov, distance, func_type)

    real(dp), intent(in)                  ::        feature_matrix1(:,:), feature_matrix2(:,:), bent_space1(:,:), bent_space2(:,:), iv_weights(:), sigma_cov
    real(dp), intent(out), optional       ::        distance
    real(dp)                              ::        cov, d_sq
    integer, intent(in), optional         ::        func_type
    integer                               ::        i, k, selector
    
    k = size(feature_matrix1(:,1))    
    d_sq = 0.0_dp    
    do i=1, k
       d_sq = d_sq + (distance_in_bent_space(feature_matrix1(i,:), feature_matrix2(i,:), bent_space1, bent_space2))**2 / (2.0_dp * (iv_weights(i)**2))
    enddo
 
    ! norm with the dimensionality k of the Internal Space
    d_sq = d_sq/real(k,dp)                            

    if (present(distance))  distance = d_sq  ! N.B. distance squared

    if (present(func_type)) then
         selector=func_type
    else
         selector=0
    endif

    SELECT CASE (selector)
    CASE (0)
       cov = exp(-0.5_dp*d_sq / sigma_cov**2)
    CASE (1)
       cov = sqrt(d_sq)
    CASE (2)
       cov = d_sq
    CASE (3)
       cov = (sqrt(d_sq))**3
    CASE (4)
       cov = d_sq**2
    CASE (5)  ! Laplace
       cov = exp(- sqrt(d_sq) / sigma_cov)
    END SELECT

 endfunction cov


 function  distance_in_bent_space(vect1, vect2, bent_space1, bent_space2)

    real(dp), intent(in)                    :: vect1(3), vect2(3), bent_space1(:,:), bent_space2(:,:)
    real(dp)                                :: distance_in_bent_space

    distance_in_bent_space = norm( (bent_space1 .mult. vect1) - (bent_space2 .mult. vect2))  
 
 endfunction distance_in_bent_space

! function eval_lhood(covariance_matrix,FCM)
!   real(dp), dimension(:,:), intent(in) :: covariance_matrix, FCM
!   real(dp) :: eval_lhood   
!   real(dp), dimension(:,:), allocatable :: ICF
!   type(LA_Matrix) :: LA_covariance
      
!   LA_covariance = covariance_matrix
!   allocate(ICF(size(FCM,1),size(FCM,1)))
   
!   call Matrix_QR_Solve(LA_covariance, FCM, ICF)
   !write(*,*) "DEBUG ICF: ", ICF
      
!   eval_lhood = 1.5_dp * LA_Matrix_LogDet(LA_covariance) + 0.5_dp * trace(ICF)
   
!   call finalise(LA_covariance)
!   if(allocated(ICF)) deallocate(ICF)
   
! endfunction eval_lhood


 subroutine matrix_statistics(in_matrix, max_value, mean_value, deviation_value, n)

   real(dp), intent(in)                    ::  in_matrix(:,:)
   integer, intent(in)                     ::  n
   integer                                 ::  i, j, t
   real(dp), intent(out)                   ::  max_value, mean_value, deviation_value
   real(dp), allocatable                   ::  data_array(:)

   allocate(data_array(n*(n-1)/2))

    t=0
    do i=1, n
       do j=1,n
          if (i>j) then
             t=t+1
             data_array(t) = in_matrix(i,j)  
          endif
       enddo
    enddo
    mean_value = sum(data_array)/real(size(data_array),dp)
    max_value  = maxval(data_array)
    
    deviation_value = 0.0_dp
    do i=1, size(data_array)
       deviation_value = deviation_value + ( (data_array(i) - mean_value)**2 )
    enddo
    deviation_value = sqrt(deviation_value/real(size(data_array),dp))

    deallocate(data_array)

 end subroutine  matrix_statistics


 subroutine sorting_configuration(matrix_predict, matrix_data, matrix_predict_norm, matrix_data_norm, sigma, distance_confs_unsorted, distance_index)
 
   ! returns a list of indexes corresponding to the configurations ordered according to their distance to target 
   real(dp), intent(in)                        ::  matrix_predict(:,:), matrix_data(:,:,:), matrix_predict_norm(:,:), matrix_data_norm(:,:,:), sigma(:)
   real(dp), intent(inout)                     ::  distance_confs_unsorted(:)
   real(dp), allocatable                       ::  distance_confs(:) ! only for sorting
   integer,  intent(inout)                     ::  distance_index(:)
   real(dp)                                    ::  cov_tmp  
   integer                                     ::  i
 
   call system_timer('Distance Calculation for sorting')

   !$omp parallel default(none) shared(matrix_data, matrix_data_norm, matrix_predict, matrix_predict_norm, sigma, distance_confs_unsorted, distance_index) private(i, cov_tmp)
   !$omp do  
   do i=1, size(matrix_data(1,1,:)) 
      cov_tmp = cov(matrix_predict, matrix_data(:,:,i), matrix_predict_norm, matrix_data_norm(:,:,i), sigma, 1.0_dp, distance=distance_confs_unsorted(i)) 
      distance_index(i) = i  
   enddo
   !$omp end parallel

   call system_timer('Distance Calculation for sorting')

   call system_timer('Sorting the DATABASE')

!  call insertion_sort(distance_confs_unsorted, idx=distance_index)
!  heap sorting below
   allocate(distance_confs(size(distance_confs_unsorted)))
   distance_confs = distance_confs_unsorted
   call heap_sort(distance_confs, i_data=distance_index)
   call system_timer('Sorting the DATABASE')
   deallocate(distance_confs)
 end subroutine sorting_configuration


 subroutine  real_expection_sampling_components(ivs_direction_matrix, internal_component_matrix, target_force)

    real(dp), intent(in)                        ::  ivs_direction_matrix(:,:), internal_component_matrix(:)
    real(dp), intent(out)                       ::  target_force(3)
    real(dp), allocatable                       ::  force_array(:,:)
    real(dp)                                    ::  const=0.1_dp, f_deviv(3), converting_matrix(3,3), converting_matrix_inv(3,3), force_value_matrix(3)
    integer                                     ::  i, j, t, f_counter, s

    f_counter=0

    s=size(internal_component_matrix)  ! the number of internal directions, equavalent to 'k'
    allocate(force_array(s*(s-1)*(s-2)/6,3)) ! number of combinations: C_{k}^{3}
    force_array=0.0_dp

    do i=1, s-2
       do j=i+1, s-1
          do t=j+1, s
               f_counter = f_counter + 1
               converting_matrix(:,1)=ivs_direction_matrix(i,:)
               converting_matrix(:,2)=ivs_direction_matrix(j,:)
               converting_matrix(:,3)=ivs_direction_matrix(t,:)
               force_value_matrix(1)=internal_component_matrix(i)
               force_value_matrix(2)=internal_component_matrix(j)
               force_value_matrix(3)=internal_component_matrix(t)

               if ((norm(converting_matrix(:,1) .cross. converting_matrix(:,2)) > const) .and. (norm(converting_matrix(:,2) &
                           .cross. converting_matrix(:,3)) > const) .and. (norm(converting_matrix(:,3) .cross. converting_matrix(:,1)) > const)) then
                 f_counter = f_counter + 1

                 write(*,*) "the converting matrix", converting_matrix
                 call inverse_svd_threshold(converting_matrix, 10.0_dp, result_inv=converting_matrix_inv)
                 force_array(f_counter,:)=transpose(converting_matrix_inv) .mult. force_value_matrix
                 write(*,*) "the inverse matrix : ", converting_matrix_inv
                 write(*,*) "the number of sequence of ivs:", i, j, t
                 call print("Target Force Distribution : "//force_array(f_counter,1)//" "//force_array(f_counter,2)//" "//force_array(f_counter,3))
               endif
           enddo
       enddo
     enddo

   do i=1, 3
      target_force(i) = sum(force_array(:,i)) / real(f_counter, dp)
   enddo
   f_deviv=0.0_dp
   do j=1, 3
      do i=1,f_counter
        f_deviv(j) = f_deviv(j) + (force_array(i,j)-target_force(j))**2 
     enddo
   enddo
   do j=1, 3
      f_deviv(j) = sqrt(f_deviv(j)/ real(f_counter,dp))
   enddo
   write(*,*) "deviation of the prediction: ", f_deviv
   deallocate(force_array)
 end subroutine real_expection_sampling_components

 subroutine internal_dimension_mark(in_matrix, TOL, mark)

    real(dp), intent(in)                ::  in_matrix(:,:), TOL
    integer, intent(out)                ::  mark(3) 
    integer                             ::  i

    mark(:)=1
    do i=1,3
       !if (sum(abs(in_matrix(:,i))/size(in_matrix(:,1)))< TOL)  then ! using Mean as Indicator
       if ( maxval(abs(in_matrix(:,i))) < TOL ) then ! using Max as indicator
              call print("dimension index: "//maxval(abs(in_matrix(:,i))) )
              mark(i)=0
              call print('the '//i//'-th dimension is zero')
       endif 
    enddo

 end subroutine  internal_dimension_mark 


 subroutine rank_lowered_inverse(in_matrix, out_inv)

   real(dp), intent(in)               ::    in_matrix(3,3)
   real(dp), intent(out)              ::    out_inv(3,3)
   real(dp)                           ::    a, b, c, d, m_factor

   a=in_matrix(1,1)
   b=in_matrix(1,2)
   c=in_matrix(2,1)
   d=in_matrix(2,2)
   m_factor = a*d-b*c
   out_inv=0.0_dp

   out_inv(1,1) =  d/m_factor
   out_inv(1,2) = -b/m_factor
   out_inv(2,1) = -c/m_factor
   out_inv(2,2) =  a/m_factor
 
 end subroutine rank_lowered_inverse

 subroutine do_optimise_likelihood(sigma_error, sigma_covariance)

   type(optimise_likelihood)           :: am_likelihood 
   logical                             :: use_qr
   integer                             :: am_data_size
   character(len=1), dimension(1)      :: am_mold
   character(len=1), allocatable       :: am_data(:)
   real(dp)                            :: x_sigma(2), covariance_tol
   real(dp), intent(inout)             :: sigma_error, sigma_covariance
 !  real(dp), dimension(:,:), target    :: distance_matrix, covariance_matrix, force_covariance_matrix
   

   covariance_tol = 1.0e-2
   use_qr=.true.

!   am_likelihood%use_qr = use_qr ! filling the am_likelihood container
!   am_likelihood%distance_matrix => distance_matrix
!   am_likelihood%covariance_matrix => covariance_matrix
!   am_likelihood%force_covariance_matrix => force_covariance_matrix

   am_data_size = size(transfer(am_likelihood,am_mold))
   allocate(am_data(am_data_size))
   am_data = transfer(am_likelihood,am_mold) ! pouring the content of am_likelihood into am_data
     
   ! Testing the implementation of likelihood and dlikelihood.
   if(.false.) then
      call print('starting test gradient')
      call verbosity_push(PRINT_VERBOSE)
      write (*,*) test_gradient((/sigma_error,sigma_covariance/),likelihood,dlikelihood,data=am_data)
      call verbosity_pop()
   endif
   call print_title('likelihood optimisation')

   if (use_qr) then
      call print('doing likelihood optimisation with QR decomposition')
   else
      call print('doing likelihood optimisation with Cholesky decomposition')
   end if
   x_sigma = (/sigma_error, sigma_covariance/)    ! starting points for the values
   ! n = minim( x_sigma, likelihood, dlikelihood, method='cg', convergence_tol=covariance_tol, max_steps = 100, data=am_data)
   ! x_sigma is the input vector, lh is the target function to minimise (using dlh), using conjugate gradient, with that tolerance, 
   ! with that max number of steps, and using the data contained in am_data
   write(*,*) "Fixed LIKELIHOOD", likelihood(x_sigma, am_data), dlikelihood(x_sigma, am_data)
   deallocate(am_data)

   sigma_error = x_sigma(1)  
   sigma_covariance = x_sigma(2)   ! return the values of sigma_error and sigma_covariance after the optimisation
   call print('Optimised hypers. ERROR = '//sigma_error//' COVARIANCE = '//sigma_covariance)
  
 end subroutine do_optimise_likelihood

subroutine sigma_factor_from_single_dimension(ivs_set, sigma, print_verbosity, TOL_REAL, data_dir)
	! calculation normalisation factor on IV_i
	integer                             :: t, i, j, n, k
	real(dp)                            :: dist_primitive, distance_ivs_stati, feature_len
	real(dp), intent(inout)             :: sigma(:)
real(dp), allocatable               :: distance_ivs(:,:), mean_value(:), max_value(:), deviation_value(:), ivs_set_norm(:,:,:)
	real(dp), intent(in)                :: ivs_set(:,:,:), TOL_REAL
	logical, intent(in)                 :: print_verbosity 
	character(STRING_LENGTH),intent(in) :: data_dir      

	n=size(ivs_set(1,1,:))
	k=size(ivs_set(:,1,1))
allocate(ivs_set_norm(k,3,n))

	do j=1, n
	do t=1, k
feature_len=norm(ivs_set(t,:,j))
	if (feature_len < TOL_REAL)  then 
	feature_len=1.0_dp
	call print("WARNING: Numerical Limit in getting the unit direction of IVs : Generating sigma.dat")
	endif
	ivs_set_norm(t,:,j) = ivs_set(t,:,j)/feature_len
	enddo
	enddo

allocate(distance_ivs(n,n), mean_value(k), max_value(k), deviation_value(k))

	call system_timer('Getting Hyper-parameters from the DATABASE')

	dist_primitive=0.0_dp
	do t = 1, k
	do i = 1, n
	do j=1, n
distance_ivs(i,j) = distance_in_bent_space(ivs_set(t,:,i), ivs_set(t,:,j), ivs_set_norm(:,:,i), ivs_set_norm(:,:,j))
	if ((i>j) .and. (print_verbosity)) call print("Dimensional Distance : "//distance_ivs(i, j)//"  dimension : "//t)
	enddo   ! i
	enddo     ! j

call matrix_statistics(distance_ivs, max_value(t), mean_value(t), deviation_value(t), n)
	call print("Statistics max value: "//max_value(t)//" mean value: "//mean_value(t)//" deviation_value: "//deviation_value(t))
	dist_primitive = dist_primitive + (mean_value(t)/deviation_value(t)/2.0_dp)**2 
	enddo

dist_primitive = sqrt(dist_primitive / real(k,dp)) 
	call print("primitive_distance : "//dist_primitive)

	if (print_verbosity) then
	! norm distance with hyperparameters derived from statistics analysis      
	do i=1,n
	do j=1,n
	distance_ivs_stati = 0.0_dp   
	do t=1, k                 
	distance_ivs_stati = distance_ivs_stati + &
	(distance_in_bent_space(ivs_set(t,:,i), ivs_set(t,:,j), ivs_set_norm(:,:,i), ivs_set_norm(:,:,j))/2.0_dp/deviation_value(t))**2  
	enddo
	distance_ivs_stati = sqrt(distance_ivs_stati/real(k,dp))          ! to take the square root and get a norm distance 
	if (i>j)   call print("Normalised Dimensional Distance : "//distance_ivs_stati//"  Normalised Value : "//distance_ivs_stati/dist_primitive)
	enddo
	enddo
	endif ! print_verbosity

    ! to write the derived sigma vector into file "sigma.dat" for later use
    ! sigma array is derived from the standard deviation on each single dimension of the IVs space
    sigma = deviation_value           
    open(unit=1, status='replace', file=trim(data_dir)//'sigma.dat')
    write(1,*) sigma
    close(unit=1)
     
    deallocate(distance_ivs, mean_value, max_value, deviation_value, ivs_set_norm)
    call system_timer('Getting hyper-parameters from the DATABASE')
end subroutine sigma_factor_from_single_dimension

subroutine write_distance_file_real(d, data_dir)

   character(STRING_LENGTH), intent(in)                            :: data_dir
   real(dp), intent(in)                                            :: d

!   OPEN(2, status='replace',file=trim(data_dir)//'Dist.dat',  form='UNFORMATTED')
!   !write(2)  d
!   close(2)
    call fwrite_array_d(1, 0.0_dp, trim(data_dir)//'Dist.dat')
end subroutine write_distance_file_real

!subroutine  refresh_info(data_dir, f_proj, f_norm, ivs_set)
!
!      character(STRING_LENGTH), intent(in)                            ::   data_dir
!      real(dp), intent(in)                                            ::   f_proj(:), ivs_set(:,:) 
!      real(dp), intent(in)                                            ::   f_norm  
!      integer                                                         ::   i
!      call write_distance(0.0_dp, data_dir) ! always refresh the Dist.dat whenever print_data
!
!      OPEN(2,status='replace',file=trim(data_dir)//'Force.dat', form='UNFORMATTED')
!      OPEN(3,status='replace',file=trim(data_dir)//'IV.dat', form='UNFORMATTED')
!
!      write(2)  f_proj(:), f_norm
!
!      do i=1,3
!             write(3)  ivs_set(:,i)
!      enddo
!      close(2)
!      close(3)
!end subroutine refresh_info

subroutine write_distance_file(distance_confs_unsorted, data_dir)

   character(STRING_LENGTH), intent(in)                            :: data_dir
   real(dp), intent(in), dimension(:)                              :: distance_confs_unsorted 
 
   call system_timer('Writing Distance File : Trivial Usually')
!   OPEN(2, status='old',file=trim(data_dir)//'Dist.dat', form='UNFORMATTED', position='APPEND')
!   write(2)  distance_confs_unsorted(:)
!   close(2)
   call fwrite_array_d(size(distance_confs_unsorted), distance_confs_unsorted(:), trim(data_dir)//'Dist.dat')
   call system_timer('Writing Distance File : Trivial Usually')

end subroutine write_distance_file

subroutine  read_distance_file(distance_matrix, data_dir)

       character(STRING_LENGTH), intent(in)                               :: data_dir
       real(dp), intent(inout), dimension(:, :)                           :: distance_matrix
       integer                                                            :: i, j

       call system_timer('Read Pair Distance and Bulid Covariance')
 !      OPEN(2,status='old',file=trim(data_dir)//'Dist.dat',form='UNFORMATTED')
       call fread_array_d(1, distance_matrix(1,1), trim(data_dir)//'Dist.dat') 

       do i=2, size(distance_matrix(:,1))
              call  fread_array_d(i-1, distance_matrix(i,:(i-1)), trim(data_dir)//'Dist.dat')
  !           read(2) distance_matrix(i,:(i-1))                 ! for row <= column
              distance_matrix(i,i)=0.0_dp
       enddo
  
       do i=1, size(distance_matrix(:,1))
            do j= i, size(distance_matrix(:,1)) 
                distance_matrix(i,j) = distance_matrix(j,i)  ! for row > column
            enddo
        enddo
  !      close(2)
       call system_timer('Read Pair Distance and Bulid Covariance')

end subroutine read_distance_file 
!
function fgp_calc(at_in) 

   type(Atoms)         :: at_in, fgp_calc
   type(Potential)     :: pot
   type(Dictionary)    :: params
   real(dp)            :: r_cut, feature_len, thresh, sigma_error, error_ls, error_gp, error_gp_sq, time, error_ls_max
   real(dp)            ::  sigma_covariance, sigma_cov, dist_shift_factor, dim_tol, error_frame, force_magntd_pred, dist
   real(dp), dimension(:), allocatable           :: r_grid, m_grid, sigma, covariance_pred,force_magntd_data,force_magntd_data_select 
   real(dp), dimension(:), allocatable           :: force_proj_test, force_proj_ivs_pred, distance_confs_unsorted
   real(dp), parameter                           :: TOL_REAL=1e-7_dp, SCALE_IVS=100.0_dp
   real(dp)                                      :: force(3), variance_xyz(3), feature_inner_matrix(3,3), feature_inv(3,3), kappa
   real(dp), dimension(:,:,:), allocatable       :: ivs_set_norm, ivs_set, ivs_set_pred, nearest_data
   real(dp), dimension(:,:), allocatable         :: force_proj_ivs, force_proj_ivs_select, inv_covariance, out_u, out_vt, geometry_matrix
   real(dp), dimension(:,:), allocatable, target :: covariance, distance_matrix, distance_matrix_dij, force_covariance_matrix, covariance_tiny
   real(dp), dimension(:,:), allocatable         :: ivs_set_norm_pred, ivs_set_norm_pred_t
   real(dp), dimension(:,:), pointer             :: force_ptr, force_ptr_mm, force_ptr_tff
   real(dp), dimension(:), pointer               :: pred_err
   integer                                       :: i,j, k, n,  t,  n_center_atom, error 
   integer                                       :: add_vector, n_teach, func_type, weight_index, temp_integer, ii, do_atz_index
   integer	                 		 :: local_ml_optim_size, dimension_mark(3)
   integer, dimension(:), allocatable            :: distance_index, mark_zero_ivs, species
   logical                         ::             do_gp, do_sigma, least_sq, do_svd, print_verbosity,  do_insert, do_dynamic_x, write_data_dij, read_data_dij 
   character(STRING_LENGTH)                      :: data_dir
   
   ! set by default exc. at_in
   call initialise(params)  
   call param_register(params, 'r_cut',  '8.0', r_cut, "the cutoff radius for the spherical atomic environment")
   call param_register(params, 'thresh', '10.0', thresh, "the threshold for doing the Singular Value Decompostion of the Covariance Matrix")
   call param_register(params, 'weight_index', '0', weight_index, "Type of weight function in deriving IVs")
   call param_register(params, 'sigma_error', '0.05', sigma_error, "the noise assumed on the teaching data")
   call param_register(params, 'func_type', '0', func_type, "which kernel function is used to build the covariance matrix")
   call param_register(params, 'do_gp',  'T', do_gp, "true for doing a gaussian processes")
   call param_register(params, 'n_teach', '1000', n_teach, "the number of relevant confs you would like to do machine learning with")
   call param_register(params, 'do_sigma', 'F', do_sigma, "set true to print out the distance on every single dimension of the IVs space")
   call param_register(params, 'do_dynamic_x', 'F', do_dynamic_x, "true for updating sigma.dat every test configuration based on nearest data points")
   !call param_register(params, 'top_teach', 'T', top_teach, "only the first configuration is considered when doing teaching")
   call param_register(params, 'least_sq',    'T', least_sq, "if true, the internal force components will be tranformed to real force using least squares")
   call param_register(params, 'do_svd', 'F', do_svd, "if true, doing inverting by SVD")
   call param_register(params, 'sigma_cov', '1.0', sigma_cov, "correlation length for pairs of data points") 
   call param_register(params, 'dist_shift_factor', '0.2', dist_shift_factor, "distance shifted with a given factor") 
   ! optional for updating sigma_cov
   call param_register(params, 'do_atz_index', '0', do_atz_index, "0: charge-weight disabled; 1: using at%Z, 2: using ionic charge instead")
   call param_register(params, 'dim_tol', '0.00005', dim_tol, "threshold for determing the dimensionality of the group of IVs")
   !call param_register(params, 'test_file', 'test.xyz', test_file, "file to read the testing configurations from")
   call param_register(params, 'data_dir', './info/', data_dir, "the directory where data files locate")
   call param_register(params, 'do_insert', 'F', do_insert, "INsert any addiitional vector as force_plus")
   call param_register(params, 'print_verbosity', 'F', print_verbosity, "if true, print out verbosity stuff")
   !call param_register(params, 'out_file', 'out.xyz', out_file, "output configuration file containing the predicted forces")
   call param_register(params, 'read_data_dij', 'F', read_data_dij, 'whether or not writing the distance file')
   call param_register(params, 'write_data_dij', 'F', write_data_dij, 'whether or not reading the pair distnace matrix') 
   call param_register(params, 'local_ml_optim_size', '0', local_ml_optim_size, "Optimise sigma_error and sigma_cov using local_ml_optim_size LS confs. If 0, no optimisation is performed")
   
   if (.not. param_read_args(params, task="fgp command line params")) then
      call print("Usage: fgp [options]")
      call system_abort("Confused by command line arguments")
   end if

   call print_title('params')
   call param_print(params)
   call finalise(params)

    call print_title('READING TEACHING INFORMATION')
    call system_timer('Collecting Information From Files')
    call load_iv_params(data_dir, r_grid, m_grid, species, k, add_vector, n_data=n)
    k = k + add_vector  ! adding the Classical Force Vectors

    allocate(force_proj_ivs(k,n))
    allocate(ivs_set(k,3,n))
    allocate(ivs_set_norm(k,3,n))
    allocate(force_magntd_data(n))
    allocate(sigma(k))

    OPEN(2,status='old',file=trim(data_dir)//'Force.dat',form='UNFORMATTED')
    OPEN(3,status='old',file=trim(data_dir)//'IV.dat',form='UNFORMATTED')

    do t=1, n
       read(2) force_proj_ivs(:, t), force_magntd_data(t)
       !write(*,*) "%", force_proj_ivs(:, t), force_magntd_data(t)
       do i=1, 3
          read(3) ivs_set(:,i,t)
       enddo
       !write(*,*) "$", ivs_set(:,:,t)
    enddo
    close(2)
    close(3)

  ! calculating the normalised vectors
  do t=1, n
    do j=1, k
     feature_len = norm(ivs_set(j,:, t))
     if (feature_len < TOL_REAL)  then
             feature_len=1.0_dp
             call print("WARNING: TEACHING, encountered the numerical limit in getting the unit direction of IVs")
     endif
     ivs_set_norm(j,:,t)=ivs_set(j,:, t)/feature_len
    enddo
 enddo
  
    ! reading sigma file
    open(unit=1, status='old', file=trim(data_dir)//'sigma.dat', form='formatted')
    read(1, *) sigma
    close(unit=1)
    call system_timer('Collecting Information From Files')
    call print('sigma is:    '//sigma)

!
if (n_teach > n) then
  n_teach = n
  call print("the Actual Number of Configs for Prediction: "//n_teach)
  call print("WARNNING: Trying to select more configurations than the database has")
endif
!
allocate(covariance_pred(n_teach))
allocate(force_proj_ivs_pred(k))
allocate(mark_zero_ivs(k))
allocate(geometry_matrix(3,k))
allocate(covariance(n_teach,n_teach)) 
allocate(distance_confs_unsorted(n), distance_index(n))
if (local_ml_optim_size > 0) then 
allocate(distance_matrix(local_ml_optim_size,local_ml_optim_size), force_covariance_matrix(local_ml_optim_size,local_ml_optim_size), covariance_tiny(local_ml_optim_size,local_ml_optim_size))
endif
allocate(inv_covariance(n_teach,n_teach))
allocate(force_proj_ivs_select(k, n_teach))
allocate(force_magntd_data_select(n_teach))
allocate(ivs_set_norm_pred(k,3))
if (read_data_dij) allocate(distance_matrix_dij(n, n))

if (read_data_dij) then
          call system_timer('Read Pair Distance and Bulid Covariance')
          call read_distance_file(distance_matrix_dij(:,:), data_dir)
          !write(*,*)  'Distance Read', distance_matrix_dij
          call system_timer('Read Pair Distance and Bulid Covariance')
endif

   if (.not. get_value(at_in%params, 'time', time))  then
     call print("TIME of FRAME : UNKOWN")
   else
     call print("TIME of FRAME :"//time)
   endif
   call add_property(at_in, 'force', 0.0_dp, n_cols=3, overwrite=.false.) ! false in case DFT force is reading in
   call assign_property_pointer(at_in, 'force', ptr=force_ptr)

   ! for writing out the predicted error
   call add_property(at_in, 'pred_err', 0.0_dp, n_cols=1, overwrite=.true.)
   call assign_property_pointer(at_in, 'pred_err', ptr=pred_err)

   if (add_vector > 0) then 
           if (do_insert) then
              call assign_property_pointer(at_in, 'force_plus', ptr=force_ptr_mm)
           else
             call Potential_Filename_Initialise(pot, args_str='IP SW', param_filename='SW.xml')
             call calc(pot, at_in, args_str='force=force_mm',  error=error)
             call assign_property_pointer(at_in, 'force_mm', ptr=force_ptr_mm)
             call finalise(pot)
           endif
   endif
   if (add_vector>1) call assign_property_pointer(at_in, 'force_i', ptr=force_ptr_tff)

   call set_cutoff(at_in, r_cut)
   call calc_connect(at_in)

   allocate(ivs_set_pred(k,3,at_in%N))

  error_frame=0.0_dp           ! for calculating the average error of the system.

  do n_center_atom=1, at_in%N     ! loop over atoms in the Frame, make prediction (incl. error)

     do j=1, k                  ! reset the mark for zero internal vectors
         mark_zero_ivs(j)=1
     enddo

     do j= 1, k-add_vector
        ivs_set_pred(j,:,n_center_atom) = internal_vector(at_in, r_grid(j), m_grid(j), species(j), n_center_atom, weight_index, do_atz_index) *SCALE_IVS
        call print("internal vectors ( "//r_grid(j)//" "//m_grid(j)//" "//species(j)//" ):  "//ivs_set_pred(j,1,n_center_atom)//"  "//ivs_set_pred(j,2,n_center_atom)//"  "//ivs_set_pred(j,3,n_center_atom))
       feature_len = norm(ivs_set_pred(j,:,n_center_atom))

        if (feature_len < TOL_REAL)  then
            feature_len=1.0_dp
            mark_zero_ivs(j)=0
            call print("WARNING: PREDICTION, encountered the numerical limit in getting the unit direction of IVs")
        endif

        ivs_set_norm_pred(j,:) = ivs_set_pred(j,:,n_center_atom)/feature_len
     enddo

     if (add_vector>0) then
         do j = k-add_vector+1, k
            temp_integer=k-j
            select case (temp_integer)
            case (0)
               ivs_set_pred(j,:,n_center_atom)=force_ptr_mm(:, n_center_atom)
            case (1)
               ivs_set_pred(j,:,n_center_atom)=force_ptr_tff(:, n_center_atom)
            end select
            feature_len = norm(ivs_set_pred(j,:,n_center_atom))

            if (feature_len < TOL_REAL) then
               feature_len=1.0_dp
               mark_zero_ivs(j)=0
               call print("WARNING: PREDICTION, encountered the numerical limit in getting the unit direction of IVs")
            endif
            
            ivs_set_norm_pred(j,:) = ivs_set_pred(j,:,n_center_atom)/feature_len
            if (print_verbosity) call print("added internal vectors : "//ivs_set_pred(j,1,n_center_atom)//"  "//ivs_set_pred(j,2,n_center_atom)//"  "//ivs_set_pred(j,3,n_center_atom))
         enddo
      endif


     if (print_verbosity) then
       do j=1, k 
          do ii=1, k   
           if(print_verbosity) then
              call print("ATOM : "//n_center_atom//" ivs_set: ("//j//" "//ii//") "//dot_product(ivs_set_pred(j,:,n_center_atom), ivs_set_norm_pred(ii,:)) )
           endif
          enddo   
       enddo
     endif

      
      ! initialise distance_index(:)
      do t=1, n
         distance_index(t) =t 
      enddo

        ! do the sorting and selection
      call sorting_configuration(ivs_set_pred(:,:,n_center_atom), ivs_set, ivs_set_norm_pred, ivs_set_norm, sigma, distance_confs_unsorted, distance_index)
      call print("Min and Max DISTANCE with Index after Sorting: "//distance_confs_unsorted(distance_index(1))//" and "// &
	                                    distance_confs_unsorted(distance_index(n_teach))//"  the INDEX: "//distance_index(1)//" and "//distance_index(n_teach))
       !!!!!!!!!!!!!
    if (do_dynamic_x) then 
        allocate(nearest_data(k,3,n_teach))
        do ii=1, n_teach
           nearest_data(:,:,ii)=ivs_set(:,:,distance_index(ii))
        enddo             
        call sigma_factor_from_single_dimension(nearest_data, sigma, print_verbosity, TOL_REAL, data_dir)
        ! update 'sigma.dat'
        deallocate(nearest_data) 
     endif     
   
     if ( distance_confs_unsorted( distance_index(1) ) > dist_shift_factor) then  
                ! distance_confs_unsorted: the pair distance betwen test and database.
                ! sigma_cov: the original input for covariance length 
                sigma_covariance = sigma_cov*distance_confs_unsorted(distance_index(1))/dist_shift_factor    
                ! this intrinsically changes the value of SIGMA_COVARIANCE and enters next loop
                call print("sigma_covariance are shifted by a factor of :"//distance_confs_unsorted(distance_index(1))/dist_shift_factor)
                call print("Modified sigma_covariance :"//sigma_covariance)
     else
               sigma_covariance = sigma_cov
               ! the sigma_covariance used in GP
     endif

     if (.false.) then
         ! massive output, only for testing use
         call print("Detail on the teaching information")
         do ii = 1, n_teach
            call print("Index: "//distance_index(ii)//" Distance: "//distance_confs_unsorted(ii))
            call print("Force in IVs space: "//force_proj_ivs(:,distance_index(ii)))
         enddo
     endif
   
     ! max marginal likelihood 
      if (local_ml_optim_size > 0) call print("Optimise: "//local_ml_optim_size)
      !!! ML optimisation here: first calculates the necessary input matrices, then optimises
      if (local_ml_optim_size > 0) then
         call print("Local Maximum Likelihood selected")
         do t = 1, local_ml_optim_size  ! loop to generate the force_cov_mat and the distance_matrix of the restricted subset of the LS
            do j = t, local_ml_optim_size
               force_covariance_matrix(t,j) = dot_product(force_proj_ivs(:,distance_index(t)), force_proj_ivs(:,distance_index(j))) 
               ! check if possible, may be a wreck because IVs make absolute representation of the force impossible 
               force_covariance_matrix(j,t) = force_covariance_matrix(t,j)
               covariance_tiny(t,j) = cov(ivs_set(:,:,distance_index(t)), ivs_set(:,:,distance_index(j)), &
                    ivs_set_norm(:,:,distance_index(t)), ivs_set_norm(:,:,distance_index(j)), sigma, sigma_covariance, func_type=func_type, &
                    distance=distance_matrix(t,j))
               covariance_tiny(j,t) = covariance_tiny(t,j)
               distance_matrix(j,t) = distance_matrix(t,j)
            enddo
         enddo
         call do_optimise_likelihood(sigma_error, sigma_covariance) 
      endif

      if (read_data_dij) then
          do t = 1, n_teach  ! loop to generate the covariance matrix of the learning set
                 dist=distance_matrix_dij(distance_index(t), distance_index(t))
                 covariance(t,t) = cov_dij(sigma_covariance, dist, func_type)
           !      write(*,*) 'diagoanl value', covariance(t,t)
                 if (do_gp) covariance(t,t) = covariance(t,t) + sigma_error**2
                 do j = t+1, n_teach ! above diagonal
                    dist=distance_matrix_dij(distance_index(t),distance_index(j))
                    covariance(t,j) = cov_dij(sigma_covariance, dist, func_type)
                    covariance(j,t) = covariance(t,j)
                 enddo
          enddo
      else
	      call system_timer('Evaluating the Cov')
	      !Generate covariance matrix  
	      !$omp parallel default(none) shared(covariance, n_teach, ivs_set, ivs_set_norm, distance_index, sigma, sigma_error, sigma_covariance, do_gp, func_type), private(j)
	      !$omp do 
	      do t = 1, n_teach  ! loop to generate the covariance matrix of the learning set
		 covariance(t,t) = cov(ivs_set(:,:,distance_index(t)), ivs_set(:,:,distance_index(t)), &
			ivs_set_norm(:,:,distance_index(t)), ivs_set_norm(:,:,distance_index(t)), sigma, sigma_covariance, func_type=func_type)
		 if (do_gp) covariance(t,t) = covariance(t,t) + sigma_error**2
		 do j = t+1, n_teach ! above diagonal
		    covariance(t,j) = cov(ivs_set(:,:,distance_index(t)), ivs_set(:,:,distance_index(j)), &
			 ivs_set_norm(:,:,distance_index(t)), ivs_set_norm(:,:,distance_index(j)), sigma, sigma_covariance, func_type=func_type)
		    covariance(j,t) = covariance(t,j)  
		enddo
	       enddo
	      !$omp end parallel
	       call system_timer('Evaluating the Cov')
      endif ! not read_data_dij

     call system_timer('Inverting the Covariance Matrix')
      if (do_gp) then
         call inverse(covariance, inv_covariance)
      else
         ! To Do Sigular Value Decomposition (SVD): A = U*SIGMA*VT
         call inverse_svd_threshold(covariance, thresh, result_inv=inv_covariance)
      endif
      call system_timer('Inverting the Covariance Matrix')
      

      do t=1, n_teach
         force_proj_ivs_select(:,t)=force_proj_ivs(:,distance_index(t))
         force_magntd_data_select(t) = force_magntd_data(distance_index(t))  
      enddo
      
      ! the test configuration covariance vector
      do t= 1, n_teach
        covariance_pred(t) = cov(ivs_set_pred(:,:,n_center_atom), ivs_set(:,:,distance_index(t)), ivs_set_norm_pred, &
                                         ivs_set_norm(:,:,distance_index(t)), sigma, sigma_covariance, func_type=func_type)
      enddo
      
      ! the predicted force components on each of the internal directions
      force_proj_ivs_pred(:) = matmul(covariance_pred, matmul(inv_covariance, transpose(force_proj_ivs_select(:,:)) )) 
        
      ! Predict the Force Magnitude
      force_magntd_pred =  dot_product(covariance_pred, matmul(inv_covariance, force_magntd_data_select))

      ! we need to mannually set the force components to be zero, projected on the zero internal-vector directions
      do j=1, k
         force_proj_ivs_pred(j)=force_proj_ivs_pred(j)*real(mark_zero_ivs(j), dp) 
      enddo

     if (print_verbosity) then 
         do j=1, k
               call print("ivs_set_norm_pred"//ivs_set_norm_pred(j,:))
         enddo
     endif

     allocate(out_u(3,3), out_vt(3,3))      ! to obtain the transformation matrix with the principal axis
     call inverse_svd_threshold(transpose(ivs_set_norm_pred) .mult. ivs_set_norm_pred, thresh, u_out=out_u, vt_out=out_vt)

     allocate(ivs_set_norm_pred_t(k,3))

     call internal_dimension_mark(ivs_set_pred(:,:,n_center_atom) .mult. out_u, dim_tol, dimension_mark)

     ivs_set_norm_pred_t = ivs_set_norm_pred .mult. out_u
 
  ! ivs_set_pred_t contains the Internal-directions after PCA transformation.
  if (print_verbosity) then
       do j=1, k
           write(*,*) "ivs_set_norm_pred_t : ", ivs_set_norm_pred_t(j, :)
       enddo
  endif

  select case (sum(dimension_mark(:)))
    case(0)
         force= 0.0_dp
    case(1)
         force=(transpose(ivs_set_norm_pred) .mult. force_proj_ivs_pred )/real(sum(mark_zero_ivs),dp)
    case(2) 
         feature_inner_matrix=transpose(ivs_set_norm_pred_t) .mult. ivs_set_norm_pred_t
         call rank_lowered_inverse(feature_inner_matrix, feature_inv) 
         force = out_u .mult. feature_inv .mult. transpose(ivs_set_norm_pred_t) .mult. force_proj_ivs_pred 
    case(3)

         feature_inner_matrix=transpose(ivs_set_norm_pred) .mult. ivs_set_norm_pred
         if (do_svd) then   
             if (print_verbosity) write(*,*)  "feature inner matrix :", feature_inner_matrix
             call inverse_svd_threshold(feature_inner_matrix, thresh, result_inv=feature_inv)
             if (print_verbosity) write(*,*)  "inverse feature inner matrix :", feature_inv
         else
             if (print_verbosity) write(*,*)  "feature_inner_matrix :", feature_inner_matrix
             call inverse(feature_inner_matrix, feature_inv)
             if (print_verbosity) write(*,*)  "inverse feature inner matrix :", feature_inv
         endif

         if (least_sq) then
            force = feature_inv .mult. transpose(ivs_set_norm_pred) .mult. force_proj_ivs_pred
         else
            call real_expection_sampling_components(ivs_set_norm_pred, force_proj_ivs_pred, force)
         endif
   end select

  deallocate(out_u)
  deallocate(out_vt)
  deallocate(ivs_set_norm_pred_t)

   call print("Predicted Force Magnitude: "//force_magntd_pred)
   call print("force in external space : "//force//" Magnitude: "//norm(force))
   call print("the original force:"//force_ptr(:, n_center_atom)//" Magnitude: "//norm(force_ptr(:, n_center_atom)))
   call print("max error :    "//maxval(abs(force_ptr(:,n_center_atom)-force))//" norm  error :  "//norm(force_ptr(:,n_center_atom)-force))

   ! measure the consistence between the least-square inference and the components obtained from GP processes
   if (least_sq) then
       error_ls=0.0_dp
       do j=1, k
            error_ls = error_ls + (dot_product(ivs_set_norm_pred(j,:), force(:)) - force_proj_ivs_pred(j))**2
       enddo
       error_ls = sqrt(error_ls / real(k,dp))
       call print("deviation of least-square inference : "//error_ls)
   error_ls_max=maxval(abs( matmul(ivs_set_norm_pred, force) - force_proj_ivs_pred(:) ))
   call print("Max deviation Component of LS : "//error_ls_max)
   endif 

   kappa = cov(ivs_set_pred(:,:,n_center_atom), ivs_set_pred(:,:,n_center_atom), ivs_set_norm_pred, ivs_set_norm_pred, sigma, sigma_covariance, func_type=func_type) + sigma_error**2
   ! the passing of  uncertainty from  Gaussian Processes to the Force Vector
   error_gp_sq = kappa - covariance_pred .dot. matmul(inv_covariance, covariance_pred)
   error_gp = sqrt(abs(error_gp_sq))
   geometry_matrix = feature_inv .mult. transpose(ivs_set_norm_pred) 

   do j=1, 3  ! Deviation in External Cartesian Coordinate 
      variance_xyz(j)= sqrt(dot_product(geometry_matrix(j,:), geometry_matrix(j,:))) * (error_gp) 
      !call print("Predicted Error Factor :"//sum(abs(geometry_matrix(j,:)))//" Dimension :"//j)
   enddo

   call print("GP Uncertainty in Cartesian Coordinate: "//variance_xyz//" GP variance: "//error_gp//" variance Squre :"//error_gp_sq)   
   pred_err(n_center_atom)=norm(variance_xyz) 
   force_ptr(:,n_center_atom)=force
!
enddo  ! loop over positions

error_frame = maxval(pred_err) ! Cautious. could generate crazy numbers
call print("Error Bar of the Frame :"//error_frame)

deallocate(ivs_set_pred)
deallocate(ivs_set_norm_pred)  
deallocate(geometry_matrix)
deallocate(r_grid, m_grid, species, sigma)
deallocate(force_proj_ivs_pred)
deallocate(ivs_set_norm)
deallocate(force_magntd_data, force_magntd_data_select)
deallocate(ivs_set)
deallocate(mark_zero_ivs)
deallocate(force_proj_ivs)
deallocate(force_proj_ivs_select)
deallocate(covariance, inv_covariance)
deallocate(covariance_pred)
deallocate(distance_index, distance_confs_unsorted)
if(allocated(force_covariance_matrix)) deallocate(force_covariance_matrix)
if(allocated(covariance_tiny)) deallocate(covariance_tiny)
if(allocated(distance_matrix)) deallocate(distance_matrix)
if (allocated(distance_matrix_dij)) deallocate(distance_matrix_dij)
!
fgp_calc=at_in  ! output the updated atomic object
!
end function fgp_calc ! return updated at_in
!
!subroutine update_database(at) ! update database based on all info in at
!!
!          type(Atoms)                              ::   at
!          real(dp), allocatable, dimension(:)      ::   r_grid, m_grid, force_proj_test
!          real(dp), allocatable, dimension(:,:,:)  ::   ivs_set_pred
!          character(STRING_LEN)                    ::		  data_dir
!          integer                                  ::		  i, j, k, add_vector, n, n_center_atom, n0
!          logical                                  ::		  write_data_dij 
!          real(dp), dimension(:,:), pointer             ::  force_ptr
!          integer, dimension(:), allocatable            ::  species 
!
!          call system_timer('Updating the DATABSE : grid.dat force.dat ivs.dat, exclud. dist.dat')
!          ! probably first load the grid file
!          call load_iv_params(data_dir, r_grid, m_grid, species, k, add_vector, n_data=n_0)
!          ! the current k is number of defined internal vectors 
!          if (add_vector > 0)   k=k+add_vector
!          ! becomes the number of internal dimensions
!          allocate(force_proj_test(k))
!      
!          do i=1, at%N
!    
!          do j = 1, k
!               force_proj_test(j) = dot_product(force_ptr(:,n_center_atom), ivs_set_norm_pred(j,:) )
!          enddo 
!           
!          call  write_iv_params(data_dir, r_grid, m_grid, species, k-add_vector, add_vector, n+1)
!
!          !if (write_data_dij) call write_distance(distance_confs_unsorted, data_dir) ! append distance file
!          !write(*,*) 'Distance', distance_confs_unsorted(:)
!
!          OPEN(2,status='old',file=trim(data_dir)//'Force.dat', form='UNFORMATTED', position='APPEND')
!          OPEN(3,status='old',file=trim(data_dir)//'IV.dat', form='UNFORMATTED',  position='APPEND')
!          write(2)  force_proj_test(:), norm(force_ptr(:,n_center_atom))
!          do j=1, 3
!                write(3)  ivs_set_pred(:,j,n_center_atom)
!          enddo
!       
!          enddo ! loop on positions
!          close(2)
!          close(3)
!          deallocate(force_proj_test)
!          call system_timer('Updating the DATABSE : dist.dat force.dat ivs.dat')
!!
!end subroutine update_database
!
subroutine optimisation_ivs_subset(data_dir, x_index) 
   ! giving x_index in the order of decreasing weight
   !
   integer                                         ::  k, t, add_vector, n, i, j
   real(dp), dimension(:,:,:), allocatable         ::  ivs_set
   real(dp), dimension(:,:), allocatable           ::  force_proj_ivs
   real(dp), dimension(:), allocatable             ::  force_magntd_data,  uni_ave, norm_ave, r_grid, m_grid 
   integer, dimension(:), allocatable              ::  species
   integer, dimension(:), allocatable, intent(out)             ::  x_index
   real(dp)                                        ::  max_corr
   type(Dictionary)                                ::  params
   character(STRING_LENGTH), intent(in)            ::  data_dir
  
   call   load_iv_params(data_dir, r_grid, m_grid, species, k, add_vector, n_data=n)
 
   k=k+add_vector  ! adding the Classical Force Vectors
   call print_title('Starting optimise the IVs directions')

   allocate(force_proj_ivs(k,n), ivs_set(k,3,n), force_magntd_data(n))
   allocate(uni_ave(k), norm_ave(k), x_index(k))

   OPEN(2,status='old',file=trim(data_dir)//'Force.dat',form='UNFORMATTED')
   OPEN(3,status='old',file=trim(data_dir)//'IV.dat',form='UNFORMATTED')

   do t=1, n
       read(2) force_proj_ivs(:, t), force_magntd_data(t)
       do i=1, 3
          read(3) ivs_set(:,i,t)
       enddo
   enddo
   close(2)
   close(3)
   
   ! output force_proj_ivs, force_magntd_data, norm of ivs 
   do i=1, k
      uni_ave(i)  = -1.0_dp* abs(sum(force_proj_ivs(i,:)) / real(n, dp)) 
      ! universal average. to be confirmed
      norm_ave(i) = sum(abs(force_proj_ivs(i,:))) / real(n, dp)
      x_index(i) =i
   enddo
   !
   call heap_sort(uni_ave, i_data=x_index) 
   ! what if reversed order is required.
   !
   deallocate(force_proj_ivs, ivs_set, force_magntd_data, uni_ave, norm_ave) 
   call print_title('Starting optimise the IVs directions')
   !
end subroutine optimisation_ivs_subset 
!
end module force_machine_learning_module

