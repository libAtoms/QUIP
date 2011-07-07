! based on Varshney et al IEEE Comp. Graphics and Appl. v. 14 p 19 (1994)
! Not sure why they think that intersection of power-cell half spaces and enclosing
!   tetrahedron half spaces is null for internal atoms - seems to be just something
!   like Voronoi cell of atom (if all (covalent) radii are equal, anyway)
! instead, look for intersection of power cell with half spaces that don't include
!   system (iterate over 6 - could be done with 4 with a tetrahedron).  If any 
!   intersections are not null, atom has an unbounded power cell and is a surface atom.
module find_surface_atoms_module
use periodictable_module
use atoms_module
implicit none

private
public :: find_surface_atoms_atoms

contains

subroutine find_surface_atoms_atoms(at, probe_r)
   type(Atoms), intent(inout) :: at
   real(dp), intent(in) :: probe_r

   integer, pointer :: surf(:)
   integer :: i, ji, j
   real(dp) :: dist_ij, cos_ij(3), p(3), v(3), di, dj
   real(dp), allocatable :: plane_constraints(:,:)
   real(dp) :: bounds(3,2)
   real(dp) :: bounds_constraints(6,6)
   real(dp) :: nn_cutoff
   integer :: nn

   nn_cutoff = 2.0_dp*probe_r + 2.0_dp*maxval(ElementCovRad(at%Z))
   call set_cutoff(at, nn_cutoff)
   call calc_connect(at)

   call add_property(at, "surf", 0, n_cols=1, ptr=surf)

   bounds(1,1) = minval(at%pos(1,:))-(maxval(at%pos(1,:))-minval(at%pos(1,:)))
   bounds(1,2) = maxval(at%pos(1,:))+(maxval(at%pos(1,:))-minval(at%pos(1,:)))
   bounds(2,1) = minval(at%pos(2,:))-(maxval(at%pos(2,:))-minval(at%pos(2,:)))
   bounds(2,2) = maxval(at%pos(2,:))+(maxval(at%pos(2,:))-minval(at%pos(2,:)))
   bounds(3,1) = minval(at%pos(3,:))-(maxval(at%pos(3,:))-minval(at%pos(3,:)))
   bounds(3,2) = maxval(at%pos(3,:))+(maxval(at%pos(3,:))-minval(at%pos(3,:)))
   bounds_constraints(1:6,1) = (/ bounds(1,1), 0.0_dp, 0.0_dp,  -1.0_dp,  0.0_dp,  0.0_dp /)
   bounds_constraints(1:6,2) = (/ bounds(1,2), 0.0_dp, 0.0_dp,   1.0_dp,  0.0_dp,  0.0_dp /)
   bounds_constraints(1:6,3) = (/ 0.0_dp, bounds(2,1), 0.0_dp,   0.0_dp, -1.0_dp,  0.0_dp /)
   bounds_constraints(1:6,4) = (/ 0.0_dp, bounds(2,2), 0.0_dp,   0.0_dp,  1.0_dp,  0.0_dp /)
   bounds_constraints(1:6,5) = (/ 0.0_dp, 0.0_dp, bounds(3,1),   0.0_dp,  0.0_dp, -1.0_dp /)
   bounds_constraints(1:6,6) = (/ 0.0_dp, 0.0_dp, bounds(3,2),   0.0_dp,  0.0_dp,  1.0_dp /)

   bounds(1,1) = minval(at%pos(1,:))-10.0_dp*(maxval(at%pos(1,:))-minval(at%pos(1,:)))
   bounds(1,2) = maxval(at%pos(1,:))+10.0_dp*(maxval(at%pos(1,:))-minval(at%pos(1,:)))
   bounds(2,1) = minval(at%pos(2,:))-10.0_dp*(maxval(at%pos(2,:))-minval(at%pos(2,:)))
   bounds(2,2) = maxval(at%pos(2,:))+10.0_dp*(maxval(at%pos(2,:))-minval(at%pos(2,:)))
   bounds(3,1) = minval(at%pos(3,:))-10.0_dp*(maxval(at%pos(3,:))-minval(at%pos(3,:)))
   bounds(3,2) = maxval(at%pos(3,:))+10.0_dp*(maxval(at%pos(3,:))-minval(at%pos(3,:)))

   do i=1, at%N

      ! count close enough neighbours
      nn = 0
      do ji=1, atoms_n_neighbours(at, i)
	 j = atoms_neighbour(at, i, ji, distance=dist_ij)
	 if (dist_ij <= ElementCovRad(at%Z(i)) + ElementCovRad(at%Z(j)) + 2.0_dp*probe_R) then
	    nn = nn + 1
	 endif
      end do

      if (allocated(plane_constraints)) then
	 if (size(plane_constraints,2) /= nn+1) deallocate(plane_constraints)
      endif
      if (.not. allocated(plane_constraints)) allocate(plane_constraints(6,nn+1))

      ! do constraints for close enough neighbours
      nn = 0
      do ji=1, atoms_n_neighbours(at, i)
	 j = atoms_neighbour(at, i, ji, distance=dist_ij, cosines=cos_ij)
	 if (dist_ij <= ElementCovRad(at%Z(i)) + ElementCovRad(at%Z(j)) + 2.0_dp*probe_R) then
	    nn = nn + 1
	    ! pi = di^2 - (ElementCovRad(i)+probe_r)^2
	    ! pj = dj^2 - (ElementCovRad(j)+probe_r)^2
	    ! along r_ij, di+dj=dij, dj=dij-di

	    ! pi == pj 
	    ! di^2 - (Ci+pr)^2 = (dij-di)^2 - (Cj+pr)^2
	    ! di^2 - (Ci+pr)^2 = dij^2 - 2 di dij + di^2 - (Cj+pr)^2
	    ! -(Ci+pr)^2 = dij^2 - 2 di dij - (Cj+pr)^2
	    ! -2 di dij = (Cj+pr)^2-(Ci+pr)^2-dij^2
	    ! di = ((Cj+pr)^2-(Ci+pr)^2-dij^2)/(-2 dij)
	    di = ((ElementCovRad(at%Z(j))+probe_r)**2 - (ElementCovRad(at%Z(i))+probe_r)**2 - dist_ij**2)/(-2.0_dp*dist_ij)
	    dj = dist_ij - di
	    if (abs(di) < dist_ij .and. abs(dj) < dist_ij) then ! in between
	       p(:) = at%pos(:,i) + abs(di)*cos_ij(:)
	    else if (abs(di) > abs(dj)) then ! not in between, and closer to j than to i
	       p(:) = at%pos(:,i) + abs(di)*cos_ij(:)
	    else ! not in between and closer to i than to j
	       p(:) = at%pos(:,i) - abs(di)*cos_ij(:)
	    endif
	    v(:) = -cos_ij(:)
	    plane_constraints(1:3,nn) = p
	    plane_constraints(4:6,nn) = v
	 endif
      end do
      do j=1,6
	 plane_constraints(1:6,nn+1) = bounds_constraints(1:6,j)
	 if (feasible_cell_has_solution(plane_constraints, bounds)) then
	    surf(i) = 1
	    exit
	 endif
      end do
   end do

end subroutine find_surface_atoms_atoms

function feasible_cell_has_solution(plane_constraints, bounds)
   real(dp), intent(in) :: plane_constraints(:,:), bounds(3,2)
   real(dp) :: x(3)
   logical :: feasible_cell_has_solution
   logical :: constraint_check

   integer i

   call print("has_solution "//shape(plane_constraints))
   do i=1, size(plane_constraints,2)
      call print("p "//plane_constraints(1:3,i)//" v "//plane_constraints(4:6,i))
   end do

   x = seidel_solve(a=-plane_constraints(4:6,:),b=-sum(plane_constraints(1:3,:)*plane_constraints(4:6,:),1),&
      bounds=bounds, has_solution=feasible_cell_has_solution)
   call print(" has_solution " // feasible_cell_has_solution)
   if (feasible_cell_has_solution) then
      call print("solution x " //x)
      constraint_check=.true.
      do i=1, size(plane_constraints,2)
	 call print(i//" a.x "//sum(-plane_constraints(4:6,i)*x)//" <=? "//sum(-plane_constraints(1:3,i)*plane_constraints(4:6,i)) // &
	    " " // (sum(-plane_constraints(4:6,i)*x) <= sum(-plane_constraints(1:3,i)*plane_constraints(4:6,i))))
	 constraint_check = constraint_check .and. (sum(-plane_constraints(4:6,i)*x) <= sum(-plane_constraints(1:3,i)*plane_constraints(4:6,i)))
      end do
      if (.not. constraint_check) then 
	 call print("WARNING: has_solution=T but some constraint violated")
	 do i=1, size(plane_constraints, 2)
	    if (.not. (sum(-plane_constraints(4:6,i)*x) <= sum(-plane_constraints(1:3,i)*plane_constraints(4:6,i)))) then
	       call print(i//" a.x "//sum(-plane_constraints(4:6,i)*x)//" <=? "//sum(-plane_constraints(1:3,i)*plane_constraints(4:6,i)) // &
		  " " // (sum(-plane_constraints(4:6,i)*x) <= sum(-plane_constraints(1:3,i)*plane_constraints(4:6,i))))
	    endif
	 end do
      endif
   endif
end function feasible_cell_has_solution

! solve a.x <= b
recursive function seidel_solve(a, b, bounds, has_solution) result(x)
   real(dp) :: a(:,:), b(:), bounds(:,:)
   logical :: has_solution
   real(dp) :: x(size(a,1))

   real(dp), allocatable :: a_skip_m(:,:), b_skip_m(:), x_skip_m(:)
   real(dp), allocatable :: a_skip_dm(:,:), b_skip_dm(:), x_skip_dm(:), bounds_skip_dm(:,:)
   integer :: d, m, i, j, ii, jj, i_skip, j_skip
   real(dp) :: x_min, x_max
   integer, allocatable :: skip_m_index(:), skip_d_index(:)
   character(len=1024) :: prefix

   d = size(a,1)
   m = size(a,2)
! call print("seidel_solve d "//d//" m "//m)

   if (size(b) /= m) call system_abort("Bad size of b "//size(b))
   if (size(bounds,1) /= d) call system_abort("Bad size of bounds dim 1 "//size(bounds,1))
   if (size(bounds,2) /= 2) call system_abort("Bad size of bounds dim 2 "//size(bounds,2))

   has_solution = .true.
   if (d == 1) then ! 1-D, just find range of minima, maxima
      x_min = bounds(1,1)
      x_max = bounds(1,2)
      do i=1, m
	 if (a(1,i) .feq. 0.0_dp) then 
	    if (0.0_dp <= b(i)) then
	       has_solution = .true.
	    else
	       has_solution = .false.
	    endif
	 else if (a(1,i) > 0) then ! x <= b/a
	    x_max = min(x_max, b(i)/a(1,i))
	 else ! x >= b/a
	    x_min = max(x_min, b(i)/a(1,i))
	 endif
      end do ! i=1,m
      if (x_max >= x_min) then
	 has_solution = has_solution .and. .true.
      else
	 has_solution = has_solution .and. .false.
      endif
      x = 0.5_dp*(x_max+x_min)
! call print("1-D returning has_solution "//has_solution//" "//x_min//" "//x_max)
      return
   endif ! d == 1

   if (m == 0) then ! no constraints
      call print("WARNING: should never get here", PRINT_ALWAYS)
      has_solution = .true.
      x = 0.5_dp*(bounds(:,1)+bounds(:,2))
      return
   endif

   ! omit a random constraint
   i_skip = mod(ran(),m)+1
   allocate(x_skip_m(d))
   allocate(skip_m_index(m-1))
! call print("skip constraint "//i_skip)
   if (m == 1) then ! next would be a call with no constraints
! call print("faking m=0 call")
      x_skip_m = 0.5_dp*(bounds(:,1)+bounds(:,2))
   else ! some constraints left to solve
      allocate(a_skip_m(d,m-1), b_skip_m(m-1))
      ii=1
      do i=1, m
	 if (i /= i_skip) then
	    skip_m_index(ii) = i
	    ii = ii + 1
	 endif
      end do
      do j=1, d
	 a_skip_m(j,:) = a(j,skip_m_index(:))
      end do
      b_skip_m(:) = b(skip_m_index(:))
prefix=trim(mainlog%prefix)
mainlog%prefix=trim(mainlog%prefix)//"S"
      x_skip_m = seidel_solve(a_skip_m, b_skip_m, bounds, has_solution)
mainlog%prefix=trim(prefix)
! call print("got has_solution "//has_solution//" "//x_skip_m)
      deallocate(a_skip_m, b_skip_m)
      if (.not. has_solution) then
	 deallocate(x_skip_m)
	 deallocate(skip_m_index)
	 return
      endif
   endif
   if (sum(x_skip_m(:)*a(:,i_skip)) <= b(i_skip)) then
! call print("skipped constraint is satisfied (a.x="//sum(x_skip_m(:)*a(:,i_skip))//" <= b="//b(i_skip)//"), returning")
      has_solution=.true.
      x = x_skip_m
      deallocate(x_skip_m)
      deallocate(skip_m_index)
      return
   endif
! call print("skipped constraint wasn't satisfied, continuing")
   deallocate(x_skip_m)
   ! still have skip_m_index allocated

   j_skip=0
   do j=1,d
      if (a(j,i_skip) .fne. 0.0_dp) then
	 j_skip=j
	 exit
      endif
   end do
! call print("eliminating dim "//j_skip)
   if (j_skip == 0) then
! call print("can't eliminate, assuming no solution")
      has_solution = .false.
      deallocate(skip_m_index)
      return
   endif
   allocate(x_skip_dm(d-1))
   allocate(skip_d_index(d-1))
   jj=1
   do j=1, d
      if (j /= j_skip) then
	 skip_d_index(jj) = j
	 jj = jj + 1
      endif
   end do
   if (m == 1) then ! fake m=0 call
! call print("faking m=0 call")
      has_solution = .true.
      x_skip_dm(1:d-1) = 0.5_dp*(bounds(skip_d_index(:),1)+bounds(skip_d_index(:),2))
   else
      allocate(a_skip_dm(d-1,m-1), b_skip_dm(m-1))
      do i=1, m-1
	 do j=1, d-1
	    a_skip_dm(j,i) = a(skip_d_index(j),skip_m_index(i)) - a(j_skip,skip_m_index(i))*a(skip_d_index(j),i_skip)/a(j_skip,i_skip)
	 end do
	 b_skip_dm(i) = b(skip_m_index(i)) - a(j_skip,skip_m_index(i))*b(i_skip)/a(j_skip,i_skip)
      end do
      allocate(bounds_skip_dm(d-1,2))
      bounds_skip_dm(:,1) = bounds(skip_d_index(:),1)
      bounds_skip_dm(:,2) = bounds(skip_d_index(:),2)
prefix=trim(mainlog%prefix)
mainlog%prefix=trim(mainlog%prefix)//"S"
      x_skip_dm = seidel_solve(a_skip_dm, b_skip_dm, bounds_skip_dm, has_solution)
mainlog%prefix=trim(prefix)
! call print("got has_solution "//has_solution//" "//x_skip_dm)
! if (has_solution) then
! call print("check raw solution")
! do i=1, m-1
!    call print("  a.x="//sum(a_skip_dm(:,i)*x_skip_dm(:))//" <=? "//b_skip_dm(i)//" "// (sum(a_skip_dm(:,i)*x_skip_dm(:)) <=b_skip_dm(i)))
! end do
! endif
      deallocate(bounds_skip_dm)
      deallocate(a_skip_dm, b_skip_dm)
   endif
   deallocate(skip_m_index)
   if (.not. has_solution) then
! call print("no solution, returning")
      deallocate(x_skip_dm)
      deallocate(skip_d_index)
      return
   endif
! call print("has solution, back substituting")
   x(skip_d_index(:)) = x_skip_dm(:)
   x(j_skip) = (b(i_skip) - sum(a(skip_d_index(:),i_skip)*x_skip_dm(:)))/a(j_skip,i_skip)
! call print("check backsub solution")
! do i=1, m
!    call print("  a.x="//sum(a(:,i)*x(:))//" <=? "//b(i)//" "// (sum(a(:,i)*x(:)) <=b(i)))
! end do
   has_solution = .true.
   deallocate(x_skip_dm)
   deallocate(skip_d_index)

end function seidel_solve

end module find_surface_atoms_module
