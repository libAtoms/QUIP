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

module k_means_clustering_module
implicit none
private

   integer, parameter :: dp = 8

   public :: k_means_clustering_pick
   interface k_means_clustering_pick
      module procedure k_means_clustering_pick_1d, k_means_clustering_pick_nd
   end interface

contains

   subroutine k_means_clustering_pick_1d(x, cluster_indices)
      real(dp), intent(in) :: x(:)
      integer, intent(out) :: cluster_indices(:)

      integer :: k
      real(dp), allocatable :: k_means(:)
      integer :: n, rv, i, j, closest_i, closest_j
      real(dp) :: rrv
      real(dp) :: closest_r, r, t_sum
      integer, allocatable :: a(:), prev_a(:)

      k = size(cluster_indices)
      n = size(x)
      allocate(k_means(k))

      ! random initial guess
      cluster_indices = 0
      do i=1, k
	 call random_number(rrv); rv = 1+floor(rrv*n)
	 do while (any(cluster_indices == rv))
	    call random_number(rrv); rv = 1+floor(rrv*n)
	 end do
	 cluster_indices(i) = rv
      end do
      k_means = x(cluster_indices)

      allocate(a(n), prev_a(n))

      ! evaluate initial assignments
      do j=1, n
	 closest_i = 0
	 closest_r = HUGE(1.0_dp)
	 do i=1, k
	    r = sqrt((x(j)-k_means(i))**2)
	    if (r < closest_r) then
	       closest_r = r
	       closest_i = i
	    endif
	 end do
	 a(j) = closest_i
      end do

      prev_a = 0
      do while (any(a /= prev_a))
	 prev_a = a

	 ! update positions
	 do i=1, k
	    t_sum = 0.0_dp
	    do j=1, n
	       if (a(j) == i) t_sum = t_sum + x(j)
	    end do
	    k_means(i) = t_sum / real(count(a == i),dp)
	 end do

	 ! update assignments
	 do j=1, n
	    closest_i = 0
	    closest_r = HUGE(1.0_dp)
	    do i=1, k
	       r = sqrt((x(j)-k_means(i))**2)
	       if (r < closest_r) then
		  closest_r = r
		  closest_i = i
	       endif
	    end do
	    a(j) = closest_i
	 end do
      end do

      do i=1, k
	 closest_i = 0
	 closest_r = HUGE(1.0_dp)
	 do j=1, n
	    r = sqrt((x(j)-k_means(i))**2)
	    if (r < closest_r) then
	       closest_r = r
	       closest_j = j
	    endif
	 end do
	 cluster_indices(i) = closest_j
      end do

      deallocate(a, prev_a, k_means)

   end subroutine k_means_clustering_pick_1d

   subroutine k_means_clustering_pick_nd(x, cluster_indices)
      real(dp), intent(in) :: x(:,:)
      integer, intent(out) :: cluster_indices(:)

      integer :: k
      real(dp), allocatable :: k_means(:,:)
      integer :: n, rv, i, j, closest_i, closest_j, nd
      real(dp) :: rrv
      real(dp) :: closest_r, r, t_sum(size(x,1))
      integer, allocatable :: a(:), prev_a(:)

      k = size(cluster_indices)
      n = size(x,2)
      nd = size(x,1)
      allocate(k_means(nd,k))

      ! random initial guess
      cluster_indices = 0
      do i=1, k
	 call random_number(rrv); rv = 1+floor(rrv*n)
	 do while (any(cluster_indices == rv))
	    call random_number(rrv); rv = 1+floor(rrv*n)
	 end do
	 cluster_indices(i) = rv
      end do
      k_means(:,:) = x(:,cluster_indices(:))

      allocate(a(n), prev_a(n))

      ! evaluate initial assignments
      do j=1, n
	 closest_i = 0
	 closest_r = HUGE(1.0_dp)
	 do i=1, k
	    r = sqrt(sum((x(:,j)-k_means(:,i))**2))
	    if (r < closest_r) then
	       closest_r = r
	       closest_i = i
	    endif
	 end do
	 a(j) = closest_i
      end do

      prev_a = 0
      do while (any(a /= prev_a))
	 prev_a = a

	 ! update positions
	 do i=1, k
	    t_sum = 0.0_dp
	    do j=1, n
	       if (a(j) == i) t_sum(:) = t_sum(:) + x(:,j)
	    end do
	    k_means(:,i) = t_sum(:) / real(count(a == i),dp)
	 end do

	 ! update assignments
	 do j=1, n
	    closest_i = 0
	    closest_r = HUGE(1.0_dp)
	    do i=1, k
	       r = sqrt(sum((x(:,j)-k_means(:,i))**2))
	       if (r < closest_r) then
		  closest_r = r
		  closest_i = i
	       endif
	    end do
	    a(j) = closest_i
	 end do
      end do

      do i=1, k
	 closest_i = 0
	 closest_r = HUGE(1.0_dp)
	 do j=1, n
	    r = sqrt(sum((x(:,j)-k_means(:,i))**2))
	    if (r < closest_r) then
	       closest_r = r
	       closest_j = j
	    endif
	 end do
	 cluster_indices(i) = closest_j
      end do

      deallocate(a, prev_a, k_means)

   end subroutine k_means_clustering_pick_nd

end module k_means_clustering_module
