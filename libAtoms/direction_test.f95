! Directionality tester:
!
! Sets of some random atoms, selects atom 1 as the centre, 2 and 3 as current splines, then tries to
! Add other atoms to make the splines more evenly distributed
!
! Look at the structure growing in file direction.xyz
!
program dtest

  include "lotf.h"

  implicit none

  integer, parameter :: N_ATOMS = 5000
  integer, parameter :: N_EXTRA_ATOMS = 50
  type(Atoms)    :: at, movie_at
  type(table)    :: indices
  type(inoutput) :: movie
  real(dp)    :: evals(3), evects(3,3), largest, smallest, small_evect(3), ratio, best_angle, this_angle
  integer     :: i, j, n, best_atom, origin

  call System_Initialise(seed = 172846274)

  call Atoms_Initialise(at,N_ATOMS,Make_Lattice(20.0_dp))
  call Set_Atoms(at,'O')

  call Initialise(movie,'direction.xyz')

  call Table_Allocate(indices,1,0,N_EXTRA_ATOMS)

  !Randomise atom positions
  at%pos = 0.0_dp
  call randomise(at%pos,5.0_dp)
  
  !Use atoms 1 as the origin atom
  origin = 1
  at%Z(origin) = 7
  at%pos(:,origin) = 0.0_dp
  !Start with just atoms 2 and 3 in the indices list
  call Append(indices,(/2/))
  call Append(indices,(/3/))

  call CalcConnectFast(at)

  !Create an atoms object of these three atoms
  call Initialise(movie_at,indices%N+1,at%lattice)
  call Set_Atoms(movie_at,'O')
  do i = 1, indices%N
     movie_at%pos(:,i) = at%pos(:,indices%int(1,i))
  end do
  movie_at%pos(:,indices%N+1) = at%pos(:,origin)
  movie_at%Z(indices%N+1) = 7 ! Highlight origin atom
  call Print_xyz(movie_at,movie)
  call Finalise(movie_at)

  !See how things change as atoms are added
  do n = 1, N_EXTRA_ATOMS

     call Print('Calling Directionality:')

     call Directionality(at,origin,indices,evals,evects,1) !1 = Directionality Ellipsoid, 2 = SVD

     !Use the ratio of smallest to largest eigenvalues as measure of 'goodness' (0 -> 1)
     largest = -huge(0.0_dp)
     smallest = huge(0.0_dp)

     call Print('e-values : e-vectors')
     do i = 1,3
        write(line,'(f8.5,a,3(f8.5,1x))') evals(i),' : ',evects(:,i)
        call Print(line)
        if (evals(i) > largest) largest = evals(i)
        if (evals(i) < smallest) then 
           smallest = evals(i)
           small_evect = evects(:,i)
        end if
     end do

     ratio = smallest / largest
     
     call Print('')
     write(line,'(a,f0.5)') 'Ratio: ',ratio
     call Print(line)
     call Print('')

     !Now find which of the remaining atoms lies directionally closest to the weakest eigenvector
     best_angle = 0.0_dp
     best_atom = 0
     do i = 1, at%N
        if (i == origin) cycle                         ! Don't re-add the origin atom
        if (Is_In_Array(Int_Part(indices,1),i)) cycle  ! Don't re-add atoms which are already selected

        this_angle = CosAngle_To_Line(at,origin,small_evect,i)
        
        if (this_angle > best_angle) then
           best_angle = this_angle
           best_atom = i
        end if

     end do

     !Add the best atom to the indices list
     if (best_atom==0) call System_Abort('Error in adding atom to indices list: No atom selected')
     write(line,'(a,i0)') 'Adding atom ',best_atom
     call Print(line)
     call Print('')
     call Append(indices,(/best_atom/))

     !Create an atoms object of all the atoms in the group so far and print it to the movie file
     call Initialise(movie_at,indices%N+1,at%lattice)
     call Set_Atoms(movie_at,'O')
     do i = 1, indices%N
        movie_at%pos(:,i) = at%pos(:,indices%int(1,i))
     end do
     movie_at%pos(:,indices%N+1) = at%pos(:,origin)
     movie_at%Z(indices%N+1) = 7 ! Highlight origin atom
     !Draw a line in the direction of the smallest evector
     do i = -5, 5
        if (i==0) cycle
        call AddAtoms(movie_at,at%pos(:,origin)+i*small_evect,1)
     end do
     call Print_xyz(movie_at,movie)
     call Finalise(movie_at)

  end do

  !Finally, print the whole atoms object to the file to show all the atoms we could have picked from
  call Print_xyz(at,movie)

  call Free(movie)
  call Free(indices)
  call Finalise(at)

  call System_Finalise

end program dtest
