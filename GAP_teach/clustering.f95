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

module clustering_module

  ! use libatoms_module
  use error_module
  use system_module ! , only : dp, optional_default, ran_uniform, reallocate
  use linearalgebra_module

  implicit none
  private

  public :: pivot, bisect_kmedoids, cluster_kmeans, select_uniform

  integer, parameter  :: n_trial = 10
  integer, parameter  :: n_trial_k_med = 100
  real(dp), parameter :: cluster_jitter = 1.0e-7_dp

  type lst
     integer, dimension(:), allocatable :: object
     integer :: medoid
     real(dp) :: sse
     integer :: N
  endtype lst
  
  type clstr
     type(lst), dimension(:), allocatable :: cluster
     real(dp), dimension(:,:), allocatable :: dm
     integer :: N
  endtype clstr
  
  contains

  subroutine distance_matrix(x,dm,theta_fac,theta)
     real(dp), dimension(:,:), intent(in) :: x
     real(dp), dimension(:,:), intent(out) :: dm
     real(dp), intent(in), optional :: theta_fac
     real(dp), dimension(:), intent(in), target, optional :: theta

     real(dp), dimension(:), pointer :: my_theta => null()
     real(dp) :: my_theta_fac
     integer :: i, j, d, n

     my_theta_fac = optional_default(1.0_dp, theta_fac)
     d = size(x,1)
     n = size(x,2)

     if( present(theta) ) then
        my_theta => theta
     else
        allocate(my_theta(d))

        do i = 1, d
           my_theta(i) = ( maxval(x(i,:)) - minval(x(i,:)) ) 
   !        theta(i) = sqrt( & !take square root
   !                         & sum( x(i,:)**2 ) / size(x(i,:)) - &
   !                         & (sum( x(i,:) ) / size(x(i,:)))**2 )
           if( my_theta(i) .feq. 0.0_dp ) my_theta(i) = 1.0_dp
        enddo
        my_theta = my_theta * my_theta_fac
     endif

     do i = 1, n
        do j = i, n
           dm(j,i) = normsq( (x(:,j) - x(:,i)) / my_theta ) + cluster_jitter*ran_uniform()
           dm(i,j) = dm(j,i)
        enddo
        dm(i,i) = 0.0_dp
     enddo

     if( present(theta) ) then
        my_theta => null()
     else
        deallocate(my_theta)
     endif

  endsubroutine distance_matrix

  subroutine pca(x,x_mean,v)

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(:), intent(out) :: x_mean
    real(dp), dimension(:,:), intent(out) :: v

    real(dp), dimension(:), allocatable :: diag_c
    real(dp), dimension(:,:), allocatable :: cov
    integer :: i, j, d, n

    d = size(x,1)
    n = size(x,2)
    allocate(cov(d,d),diag_c(d))

    x_mean = sum(x,dim=2) / n ! empirical mean

    do i = 1, d
       do j = 1, d
          cov(j,i) = dot_product(x(i,:),x(j,:)) / n - x_mean(i)*x_mean(j)
       enddo
    enddo

    call diagonalise(cov,diag_c, evects=v)

    deallocate(cov, diag_c)

  endsubroutine pca

  subroutine pivot(x,pivout,theta_fac,theta)
     real(dp), dimension(:,:), intent(in) :: x
     integer, dimension(:), intent(out) :: pivout
     real(dp), intent(in), optional :: theta_fac
     real(dp), dimension(:), intent(in), optional :: theta
     
     real(dp), dimension(:,:), allocatable :: knn
     real(dp), dimension(:), allocatable :: ktmp
     integer, dimension(:), allocatable :: pivin
     
     integer :: stat, i, j, k, d, m, n, jtmp, jmax
     real(dp) :: dmax
     
     d = size(x,1)
     n = size(x,2)
     
     m = size(pivout)

     if( m > n ) call system_abort('pivot: required number of changes ('//m//') greater than possible number of changes ('//n//')')
     
     allocate(knn(n,n),stat=stat)
     if(stat /=0 ) call system_abort('pivot: could not allocate knn matrix.')
     
     allocate(pivin(n),ktmp(n))
     
     call distance_matrix(x,knn,theta_fac=theta_fac,theta=theta)
     do i = 1, n
        do j = 1, n
           knn(j,i) = exp(-0.5_dp*knn(j,i))
        enddo
     enddo
     
     pivin = (/ (i, i=1,n) /)
     
     do k = 1, m
        dmax = 0.0_dp
        do j = k, n
           if( dmax < knn(j,j) ) then
              jmax = j
              dmax = knn(j,j)
           endif
        enddo
        if( jmax /= k ) then
            jtmp = pivin(jmax)
            pivin(jmax) = pivin(k)
            pivin(k) = jtmp

            ktmp = knn(k,:)
            knn(k,:) = knn(jmax,:)
            knn(jmax,:) = ktmp

            ktmp = knn(:,k)
            knn(:,k) = knn(:,jmax)
            knn(:,jmax) = ktmp
         endif

         knn(k,k) = sqrt(knn(k,k))

         knn(k+1:n,k) = knn(k+1:n,k)/knn(k,k)
         do j = k+1, n
            knn(j:n,j) = knn(j:n,j) - knn(j:n,k)*knn(j,k)
         enddo

         do j = 1, n
            do i = j+1,n
               knn(j,i) = knn(i,j)
            enddo
         enddo
      enddo
        
      pivout = pivin(1:m)
     
      deallocate(knn,pivin,ktmp)
  
  endsubroutine pivot

  subroutine bisect_kmedoids(x,n_clusters_in,c,med,theta_fac,theta)
     real(dp), dimension(:,:), intent(in) :: x
     integer, intent(in) :: n_clusters_in
     integer, dimension(:), intent(out),optional :: c, med
     real(dp), intent(in), optional :: theta_fac
     real(dp), dimension(:), intent(in), optional :: theta

     type(clstr) :: my_cluster, tmp

     real(dp), dimension(:), allocatable :: dv
     real(dp) :: max_sse, min_sse, sse
  
     integer, dimension(:), allocatable :: sub_cluster1, sub_cluster2, sub_cluster1_min, sub_cluster2_min
     integer, dimension(1) :: ml
     integer  :: stat, i, j, k, km, m, n, nc, &
     & lo_med, hi_med, lo_med_new, hi_med_new, lo_med_min, hi_med_min, n1, n2, n1_min, n2_min

     n = size(x,2)

     if( n_clusters_in > n ) call system_abort('bisect_kmedoids: required number of cluster greater than total number of data points')

     if(present(c) ) c = 0
     allocate(my_cluster%dm(n,n),stat=stat)
     if(stat /=0 ) call system_abort('bisect_kmedoids: could not allocate dm matrix.')

     call distance_matrix(x,my_cluster%dm,theta_fac=theta_fac,theta=theta)

     ! start clustering
     my_cluster%N = 1                               ! start with one big cluster
     allocate( my_cluster%cluster(1) )              
     my_cluster%cluster(1)%N = n                    ! put every object in the initial cluster
     allocate( my_cluster%cluster(1)%object(n) )
     my_cluster%cluster(1)%object = (/(i,i=1,n)/)

     allocate(dv(n)) ! distance vector, the sum of square of distances of points from central object
     dv = sum(my_cluster%dm,dim=1)
     my_cluster%cluster(1)%sse = minval( dv )  ! determine initial medoid, the object that is the
     ml = minloc( dv )                         ! closest to any other object in cluster
     my_cluster%cluster(1)%medoid = ml(1)
     deallocate(dv)

     ! main loop starts here, bisects initial clusters until desired number of
     ! clusters are found
                                                       
     do
        if( my_cluster%N == n_clusters_in )  exit
        max_sse = -1.0_dp                                 ! select cluster with greatest sse
        do j = 1, my_cluster%N
           if( max_sse < my_cluster%cluster(j)%sse ) then
              i = j
              max_sse = my_cluster%cluster(j)%sse
           endif
        enddo
        nc = my_cluster%cluster(i)%N
        if( nc==1 ) cycle
        allocate( sub_cluster1(nc), sub_cluster2(nc), sub_cluster1_min(nc),sub_cluster2_min(nc) )

        min_sse = huge(1.0_dp)
        do j = 1, n_trial
           m = ceiling( ran_uniform()*(nc-1) ) ! choose a bisecting point randomly
           ml = minloc( sum( my_cluster%dm( my_cluster%cluster(i)%object(:m), my_cluster%cluster(i)%object(:m) ), dim=1) )
           lo_med_new = my_cluster%cluster(i)%object(ml(1))

           ml = minloc( sum( my_cluster%dm( my_cluster%cluster(i)%object(m+1:), my_cluster%cluster(i)%object(m+1:) ), dim=1) )

           hi_med_new = my_cluster%cluster(i)%object(ml(1) + m)

           ! the median of the 2 subclusters determined
           lo_med = 0
           hi_med = 0

           ! perform k-medoid clustering on the two subclusters
           do km = 1, n_trial_k_med
              if( (lo_med_new == lo_med) .and. (hi_med_new == hi_med) ) exit
              lo_med = lo_med_new
              hi_med = hi_med_new
              n1 = 0
              n2 = 0
              !n1 = 1
              !n2 = 1
              !sub_cluster1(n1) = lo_med
              !sub_cluster1(n2) = hi_med

              do k = 1, my_cluster%cluster(i)%N
                 if( my_cluster%dm(lo_med,my_cluster%cluster(i)%object(k)) < &
                 & my_cluster%dm(hi_med,my_cluster%cluster(i)%object(k)) ) then
                    n1 = n1 + 1
                    sub_cluster1(n1) = my_cluster%cluster(i)%object(k)
                 else
                    n2 = n2 + 1
                    sub_cluster2(n2) = my_cluster%cluster(i)%object(k)
                 endif
              enddo

              ml = minloc( sum( my_cluster%dm( sub_cluster1(:n1), sub_cluster1(:n1) ), dim=1) )
              lo_med_new = sub_cluster1(ml(1))
              ml = minloc( sum( my_cluster%dm( sub_cluster2(:n2), sub_cluster2(:n2) ), dim=1) )
              hi_med_new = sub_cluster2(ml(1))
           enddo
           sse = sum( my_cluster%dm(lo_med_new,sub_cluster1(:n1)) ) + sum( my_cluster%dm(hi_med_new,sub_cluster2(:n2)) )

           ! choose the clustering that resulted the smallest sse
           if( sse < min_sse ) then
              min_sse = sse
              sub_cluster1_min = sub_cluster1
              sub_cluster2_min = sub_cluster2
              n1_min = n1
              n2_min = n2
              lo_med_min = lo_med_new
              hi_med_min = hi_med_new
           endif
        enddo

        ! now update the the clusters with the two new subclusters
        tmp = my_cluster

        do j = 1, my_cluster%N
           deallocate( my_cluster%cluster(j)%object )
        enddo
        deallocate( my_cluster%cluster )
        my_cluster%N = my_cluster%N + 1
        allocate( my_cluster%cluster( my_cluster%N ) )

        do j = 1, my_cluster%N - 1
           if( i == j ) then
              allocate( my_cluster%cluster(j)%object(n1_min) )
              my_cluster%cluster(j)%N = n1_min
              my_cluster%cluster(j)%object = sub_cluster1_min(:n1_min)
              my_cluster%cluster(j)%sse = sum( my_cluster%dm(lo_med_min,sub_cluster1_min(:n1_min)) )
              my_cluster%cluster(j)%medoid = lo_med_min
           else
              my_cluster%cluster(j) = tmp%cluster(j)
           endif
        enddo
        allocate( my_cluster%cluster(my_cluster%N)%object(n2_min) )
        my_cluster%cluster(my_cluster%N)%N = n2_min
        my_cluster%cluster(my_cluster%N)%object = sub_cluster2_min(:n2_min)
        my_cluster%cluster(my_cluster%N)%sse = sum( my_cluster%dm(hi_med_min,sub_cluster2_min(:n2_min)) )
        my_cluster%cluster(my_cluster%N)%medoid = hi_med_min
     
        do j = 1, tmp%N
           deallocate( tmp%cluster(j)%object )
        enddo
        deallocate( tmp%cluster, sub_cluster1, sub_cluster2, sub_cluster1_min, sub_cluster2_min )

        call kmedoid(my_cluster)
     enddo

     if( present(c) ) then
        do j = 1, my_cluster%N
           do k = 1, my_cluster%cluster(j)%N
              i = my_cluster%cluster(j)%object(k)
              c(i) = j
           enddo
        enddo
     endif

     if( present(med) ) then
        do j = 1, my_cluster%N
           med(j) =  my_cluster%cluster(j)%medoid
        enddo
     endif

     do j = 1, my_cluster%N
        deallocate( my_cluster%cluster(j)%object )
     enddo
     deallocate( my_cluster%cluster, my_cluster%dm )
           
  endsubroutine bisect_kmedoids
     
  subroutine kmedoid(this)
     type(clstr), intent(inout) :: this

     type(clstr) :: tmp
     integer, dimension(:), allocatable :: medoids
     integer, dimension(1) :: ml
     integer :: n, j, k
     logical :: refined

     ! k-medoid-refinement
     n = size(this%dm,1)
     ! n: total number of objects

     tmp%N = this%N
     allocate( tmp%cluster(tmp%N), medoids(tmp%N) )
     do j = 1, tmp%N
        allocate( tmp%cluster(j)%object(n) )
        medoids(j) = this%cluster(j)%medoid
     enddo

     ! main loop starts here, perfom k-medoid clustering until medoids don't
     ! change anymore
     do 
        do j = 1, tmp%N
           tmp%cluster(j)%N = 0
        enddo
        do j = 1, n
           ml = minloc( this%dm(j,medoids) ) ! determine to which medoid each object belongs
           k = ml(1)
           tmp%cluster(k)%N = tmp%cluster(k)%N + 1
           tmp%cluster(k)%object(tmp%cluster(k)%N) = j
        enddo

        ! re-determine the medoid in each cluster
        do j = 1, tmp%N
           ml = minloc( sum( this%dm( tmp%cluster(j)%object(:tmp%cluster(j)%N), &
           & tmp%cluster(j)%object(:tmp%cluster(j)%N) ), dim=1) )
           tmp%cluster(j)%medoid = tmp%cluster(j)%object(ml(1))
        enddo
        
        refined = .true.

        ! check whether medoids have changed
        do j = 1, tmp%N
           refined = refined .and. (tmp%cluster(j)%medoid == medoids(j))
           medoids(j) = tmp%cluster(j)%medoid
        enddo
        if(refined) exit
     enddo

     ! write results
     do j = 1, tmp%N
        deallocate( this%cluster(j)%object )
        allocate( this%cluster(j)%object( tmp%cluster(j)%N ) )
        this%cluster(j)%object = tmp%cluster(j)%object(:tmp%cluster(j)%N)
        this%cluster(j)%N = tmp%cluster(j)%N
        this%cluster(j)%medoid = tmp%cluster(j)%medoid
        this%cluster(j)%sse = sum( this%dm(this%cluster(j)%medoid,&
        & this%cluster(j)%object ) )

        deallocate( tmp%cluster(j)%object )
     enddo
     deallocate( tmp%cluster, medoids )
     
  endsubroutine kmedoid

  subroutine cluster_kmeans(x,cluster_index,theta_fac,theta)
     real(dp), dimension(:,:), intent(in) :: x
     integer, dimension(:), intent(out) :: cluster_index
     real(dp), intent(in), optional :: theta_fac
     real(dp), dimension(:), intent(in), target, optional :: theta

     real(dp), dimension(:), pointer :: my_theta => null()
     real(dp) :: my_theta_fac, d_min, d_ij

     real(dp), dimension(:,:), allocatable :: cluster_centre
     integer, dimension(:), allocatable :: cluster_info
     integer :: d, n, m, i, j, k, cluster_info_old, iter
     logical :: cluster_same

     d = size(x,1)
     n = size(x,2)
     m = size(cluster_index)
     if( m > n ) call system_abort('cluster_kmeans: required number of clusters ('//m//') greater than total number of points ('//n//')')

     my_theta_fac = optional_default(1.0_dp, theta_fac)
     if( present(theta) ) then
        if( size(theta) == d) then
           my_theta => theta
        else
           allocate(my_theta(d))
           my_theta = theta(1)
        endif
     else
        allocate(my_theta(d))
        do i = 1, d
           my_theta(i) = ( maxval(x(i,:)) - minval(x(i,:)) )
           if( my_theta(i) .feq. 0.0_dp ) my_theta(i) = 1.0_dp
        enddo
        my_theta = my_theta * my_theta_fac
     endif

     allocate(cluster_centre(d,m),cluster_info(n))

     call fill_random_integer(cluster_index, n) !choose random points as cluster centres.

     cluster_centre = x(:,cluster_index)
     cluster_info = 0

     iter = 0
     do 
        iter = iter + 1
        call print("iteration: "//iter,verbosity=PRINT_NERD)
        cluster_same = .true.

!$omp parallel do default(none) shared(n,m,x,cluster_info,cluster_centre,my_theta) &
!$omp reduction(.and.:cluster_same) &
!$omp private(i,j,d_min,d_ij,cluster_info_old)
        do i = 1, n
           d_min = huge(0.0_dp)
           cluster_info_old = cluster_info(i)
           do j = 1, m
              d_ij = sum(( (cluster_centre(:,j) - x(:,i))/my_theta )**2)
              if( d_ij < d_min ) then
                 d_min = d_ij
                 cluster_info(i) = j
              endif
           enddo
           if( cluster_info_old /= cluster_info(i) ) cluster_same = cluster_same .and. .false.
        enddo
!$omp end parallel do        

!$omp parallel do default(none) shared(x,cluster_centre,cluster_info,m,d) private(j,k)
        do j = 1, m
           do k = 1, d
              cluster_centre(k,j) = sum(x(k,:),mask=(cluster_info==j)) / count(cluster_info==j)
           enddo
        enddo
!$omp end parallel do
        if( cluster_same ) exit
     enddo

     do j = 1, m
        d_min = huge(0.0_dp)
        do i = 1, n
           d_ij = sum(( (cluster_centre(:,j) - x(:,i))/my_theta )**2)
           if( d_ij < d_min ) then
              d_min = d_ij
              cluster_index(j) = i
           endif
        enddo
     enddo

     deallocate(cluster_centre, cluster_info)

     if(present(theta)) then
        my_theta => null()
     else
        deallocate(my_theta)
     endif

  endsubroutine cluster_kmeans

  subroutine select_uniform(x,index_out)
     real(dp), dimension(:,:), intent(in) :: x
     integer, dimension(:), intent(out) :: index_out

     integer :: i, j, d, n, m, n_grid, i_global, i_index_out
     integer, dimension(:), allocatable :: p_grid, i_hist, histogram, x_histogram, index_out_histogram
     real(dp), dimension(:), allocatable :: lower_bound, upper_bound, x_range

     d = size(x,1)
     n = size(x,2)
     m = size(index_out)

     if( n < m ) call system_abort('select_uniform: n = '//n//' < m = '//m)

     allocate(lower_bound(d), upper_bound(d), x_range(d), p_grid(d), i_hist(d))

     lower_bound = minval(x,dim=2)
     upper_bound = maxval(x,dim=2)
     x_range = upper_bound - lower_bound

     n_grid = ceiling( real(m,kind=dp)**(1.0_dp/real(d,kind=dp)) )
     p_grid = (/ ( n_grid**(i-1), i = 1, d ) /)

     allocate(histogram(n_grid**d))
     allocate(x_histogram(n),index_out_histogram(m))

     histogram = 0

     do i = 1, n
        do j = 1, d
           i_hist(j) = nint( ( x(j,i) - lower_bound(j) ) / x_range(j) * (n_grid-1) ) + 1
        enddo
        i_global = sum((i_hist-1)*p_grid)+1

        histogram(i_global) = histogram(i_global) + 1
        x_histogram(i) = i_global
     enddo

     index_out = 0
     i_index_out = 0
     index_out_histogram = 0

     do i = 1, n
        if( all(index_out /= 0) ) exit

        if( all(x_histogram(i) /= index_out_histogram) ) then
           i_index_out = i_index_out + 1
           index_out(i_index_out) = i
           index_out_histogram(i_index_out) = x_histogram(i)
        elseif( count( histogram(index_out_histogram(1:i_index_out)) > 0 ) >= count(histogram > 0) ) then
           i_index_out = i_index_out + 1
           index_out(i_index_out) = i
           index_out_histogram(i_index_out) = x_histogram(i)
        endif
     enddo
     deallocate(lower_bound, upper_bound, x_range, p_grid, i_hist, histogram, x_histogram,index_out_histogram)

  endsubroutine select_uniform

endmodule clustering_module
