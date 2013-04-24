!!$
!!$------------Permutation Generator---Alan Nichol---------------------------------

!!$ Generate permutations of the interatomic distance vector of
!!$ a number of atoms. Information about the symmetries present in the cluster
!!$ is specified as an array called 'equivalents' - this is generated automatically
!!$ when a permutation_data_type is initialised. At the moment the 'signature' you
!!$ give as input has a maximum of 8 atoms, but this can trivially be increased.
!!$--------------------------------------------------------------------------------
!!$--------------------------------------------------------------------------------


#include "error.inc"

module permutation_maker_module
  use error_module
  use system_module, only : dp, print, optional_default, system_timer, operator(//)

  implicit none

type permutation_data_type
   integer :: perm_number, n_perms
   integer, dimension(:), allocatable :: signature_one, signature_two, signature_three, counter, rank, dist_vec !dist_vec is internal name for descriptor
   integer, dimension(:,:), allocatable :: dist_vec_permutations
   integer, dimension(:,:,:), allocatable :: perm_array
   logical :: internal_swaps_only, initialised
endtype permutation_data_type

!% Overloaded assigment operators for permutation data objects. 
private :: permutation_data_assignment
interface assignment(=)
   module procedure permutation_data_assignment
end interface assignment(=)

contains

subroutine permutation_data_assignment(to,from)
!this is the slow way to copy this object, because it goes through the motions of initialisation again
! call permutation_data_copy for a much faster alternative
  type(permutation_data_type), intent(inout) :: to
  type(permutation_data_type), intent(in)    :: from

  ! We do not fail if *from* is unitialised, since overloaded operator
  ! routines are outside scope of error handling mechanism.
  if(.not. from%initialised) then
     call permutation_data_finalise(to)
     return
  end if 
  
  call permutation_data_initialise(to,from%signature_one,signature_two=from%signature_two,signature_three=from%signature_three,internal_swaps_only=from%internal_swaps_only)

end subroutine permutation_data_assignment

subroutine permutation_data_copy(to,from)
  type(permutation_data_type), intent(inout) :: to
  type(permutation_data_type), intent(in)    :: from

  if(.not. from%initialised) then
     call permutation_data_finalise(to)
     return
  end if 


  allocate(to%counter(size(from%counter)))
  allocate(to%rank(size(from%rank)))
  allocate(to%dist_vec(size(from%dist_vec)))
  allocate(to%dist_vec_permutations(size(from%dist_vec_permutations,1),size(from%dist_vec_permutations,2)))
  allocate(to%perm_array(size(from%perm_array,1),size(from%perm_array,2),size(from%perm_array,3)))

  allocate(to%signature_one(size(from%signature_one)))
  to%signature_one = from%signature_one
  if (allocated(from%signature_two)) then
    allocate(to%signature_two(size(from%signature_two)))
    to%signature_two = from%signature_two    
    if (allocated(from%signature_three)) then 
      allocate(to%signature_three(size(from%signature_three)))
      to%signature_three = from%signature_three
    end if
  end if
  to%counter = from%counter
  to%rank = from%rank
  to%dist_vec = from%dist_vec
  to%dist_vec_permutations = from%dist_vec_permutations
  to%perm_array = from%perm_array

  to%n_perms = from%n_perms
  to%internal_swaps_only = from%internal_swaps_only

  to%initialised = .true.
  to%perm_number = 1

end subroutine permutation_data_copy

subroutine equivalents_row_atoms(equivalents_row,signature,offset,N)
  implicit none
  integer, dimension(:), allocatable, intent(inout) :: equivalents_row
  integer, dimension(:), intent(in) :: signature
  integer, dimension(:), allocatable :: scratch_row, equivalents_temp
  integer, intent(in) :: offset, N
  integer :: z_index, i, j, repeats

  allocate(scratch_row(N))
  allocate(equivalents_temp(1))

  do z_index=1,maxval(signature)
     repeats=0
     scratch_row=0
     do i=1,size(signature)
        if (signature(i) == z_index) then
           do j=repeats,0,-1
             scratch_row(i+offset)=scratch_row(i+offset)+10**(j)
           end do
           repeats = repeats+1
        end if
     end do
     if (repeats .le. 1) cycle

     if (.not. allocated(equivalents_row)) then
         allocate(equivalents_row(N))
         equivalents_row=scratch_row
     else
         deallocate(equivalents_temp)
         allocate(equivalents_temp(size(equivalents_row)))
         equivalents_temp=equivalents_row
         deallocate(equivalents_row)
         allocate(equivalents_row(size(equivalents_temp)+N))
         equivalents_row=(/equivalents_temp,scratch_row/)
     end if           
  end do

end subroutine equivalents_row_atoms

subroutine equivalents_row_monomers(equivalents_row,N,signature,pos_a,pos_b,pos_c)

  integer,dimension(:), allocatable, intent(inout) :: equivalents_row
  integer, intent(in) :: N, pos_a, pos_b
  integer, intent(in), optional :: pos_c
  integer :: i
  integer, dimension(:), intent(in) :: signature
  integer, dimension(:), allocatable :: scratch_row, equivalents_temp

  allocate(scratch_row(N))
  allocate(equivalents_temp(1))
  scratch_row=0

  do i=1,size(signature)
    scratch_row(i+pos_a) = i
    scratch_row(i+pos_b) = i*11
    if (present(pos_c)) then
      scratch_row(i+pos_c) = i*111
    end if
  end do

  if (.not. allocated(equivalents_row)) then
      allocate(equivalents_row(N))
      equivalents_row=scratch_row
  else
      deallocate(equivalents_temp)
      allocate(equivalents_temp(size(equivalents_row)))
      equivalents_temp=equivalents_row
      deallocate(equivalents_row)
      allocate(equivalents_row(size(equivalents_temp)+N))
      equivalents_row=(/equivalents_temp,scratch_row/)
  end if     
      
end subroutine equivalents_row_monomers

subroutine permutation_data_initialise(this,signature_one,signature_two,signature_three,internal_swaps_only,error)
   implicit none
   type(permutation_data_type) :: this
   integer, dimension(:), allocatable :: counter, rank, dist_vec, equivalents_row, scratch_row, &
                                          equivalents_temp, atoms, group
   integer, dimension(:), allocatable :: signature
   integer, dimension(:), optional :: signature_two, signature_three
   integer, dimension(:) :: signature_one
   integer, dimension(:,:), allocatable :: group_array, dist_vec_permutations, equivalents
   integer, dimension(:,:,:), allocatable :: perm_array
   integer, optional, intent(out) :: error
   integer :: repeats, num_groups, N, dist_vec_n_perms, i,j,z_index,max_rank,num_distances, &
              num_perms, offset_one, offset_two, offset_three
   logical, optional :: internal_swaps_only
   logical :: two_monomers_given, three_monomers_given, my_internal_swaps_only
   real(dp) :: cutoff

   INIT_ERROR(error)
   call permutation_data_finalise(this)

   two_monomers_given=.false.
   three_monomers_given=.false.
   if(present(signature_two) .and. size(signature_two) .ge. 1) then
     two_monomers_given= .true.
     if(present(signature_three) .and. size(signature_three) .ge. 1) three_monomers_given= .true.
   end if
   my_internal_swaps_only = optional_default(.true., internal_swaps_only)

   if (three_monomers_given) then
     N = size(signature_one)+size(signature_two)+size(signature_three)
   else if (two_monomers_given) then
     N = size(signature_one)+size(signature_two)
   else
     N = size(signature_one)
   end if

   allocate(scratch_row(N))
   allocate(equivalents_temp(1))

    if (two_monomers_given .and. .not. three_monomers_given) then
      if (minval(signature_one-signature_two) .eq. 0 .and. maxval(signature_one-signature_two) .eq. 0) then
        offset_one=0
        offset_two=size(signature_one)
        call equivalents_row_monomers(equivalents_row,N,signature_one,offset_one,offset_two)
      end if
    else if (three_monomers_given) then
      if (minval(signature_one-signature_two) .eq. 0 .and. maxval(signature_one-signature_two) .eq. 0) then
       ! one and two are equivalent
        if (minval(signature_one-signature_three) .eq. 0 .and. maxval(signature_one-signature_three) .eq. 0) then
         ! all three are equivalent
          offset_one=0
          offset_two=size(signature_one)
          offset_three =size(signature_one)+size(signature_two)
          call equivalents_row_monomers(equivalents_row,N,signature_one,offset_one,offset_two,offset_three)
        else
         ! only one and two are equivalent
          offset_one=0
          offset_two=size(signature_one)
          call equivalents_row_monomers(equivalents_row,N,signature_one,offset_one,offset_two)
        end if
      else if (minval(signature_one-signature_three) .eq. 0 .and. maxval(signature_one-signature_three) .eq. 0) then
       ! only one and three are equivalent
          offset_one=0
          offset_two=size(signature_one)+size(signature_two)
          call equivalents_row_monomers(equivalents_row,N,signature_one,offset_one,offset_two)
      else if (minval(signature_two-signature_three) .eq. 0 .and. maxval(signature_two-signature_three) .eq. 0) then
       ! only two and three are equivalent
          offset_one=size(signature_one)
          offset_two=size(signature_one)+size(signature_two)
          call equivalents_row_monomers(equivalents_row,N,signature_two,offset_one,offset_two)
      end if 
    end if 

 ! if more than one signature is given, and swapping atoms between monomers is allowed
 ! then the signatures get concatenated to 'signature'. Otherwise 'signature' is just set
 ! to refer to signature_one
   if(two_monomers_given) then
     if (my_internal_swaps_only) then
        allocate(signature(size(signature_one)))
        signature = signature_one
     else if (three_monomers_given) then
        allocate(signature(size(signature_one)+size(signature_two)+size(signature_three)))
        signature = (/signature_one,signature_two,signature_three/)
     else
        allocate(signature(size(signature_one)+size(signature_two)))
        signature = (/signature_one,signature_two/)
     end if
   else
     allocate(signature(size(signature_one)))
     signature = signature_one     
   end if

   offset_one=0
   call equivalents_row_atoms(equivalents_row,signature,offset_one,N) 

   if (two_monomers_given .and. my_internal_swaps_only) then
     offset_one=size(signature)
     call equivalents_row_atoms(equivalents_row,signature_two,offset_one,N) 
     if (three_monomers_given .and. my_internal_swaps_only) then
       offset_one =size(signature)+size(signature_two)
       call equivalents_row_atoms(equivalents_row,signature_three,offset_one,N) 
     end if
   end if

   num_groups=size(equivalents_row)/N
   allocate(equivalents(num_groups,N))
 ! make the equivalents array
   equivalents =transpose(reshape(equivalents_row,(/ size(equivalents, 2), size(equivalents, 1) /)))

!!$  do i=1,size(equivalents,1)
!!$    write(*,*) equivalents(i,:)
!!$  end do


!--------- Further Array allocations and Initialisation --------------------------!
  allocate(atoms(N))
  allocate(group(N))
  do i=1,size(atoms)
    atoms(i)=i
  end do

  num_distances = N*(N-1)/2
  allocate(dist_vec(num_distances))


  num_groups = size(equivalents,1)
  allocate(counter(num_groups))
  allocate(rank(num_groups))

  !make rank vector
  do i=1,num_groups
    group(:) = equivalents(i,:)
    num_perms = num_group_perms(group)
    rank(i) = num_perms
  end do

  max_rank = maxval(rank)

  !get total number of permutations
  dist_vec_n_perms = 1
  do i=1,size(rank)
    dist_vec_n_perms = dist_vec_n_perms*rank(i)
  end do

  !initialise counter
  counter=1

  allocate(dist_vec_permutations(size(dist_vec),dist_vec_n_perms))
  allocate(group_array(max_rank,N))
  allocate(perm_array(num_groups,N,max_rank))

  !initialise arrays permutations to zero
  perm_array = 0
  dist_vec_permutations=0


!-------------------------------------------------------------------------!
!make 2D array of permutations of each group and add to 3D array perm_array
!-------------------------------------------------------------------------!
  do i = 1,num_groups 
    group(:) = equivalents(i,:)
    group_array = permute_atoms(atoms,group,N,max_rank)!this padded with zeroes in case group is of less than max_rank
    do j=1,size(group_array, 1)
      perm_array(i,:,j) = group_array(j,:)   
    end do
  end do


!-------------------------------------------------------------------------!
!Now assign relevant stuff to the permutation_data_type
!-------------------------------------------------------------------------!


  allocate(this%signature_one(size(signature_one)))
  this%signature_one=signature_one
  if (two_monomers_given) then
    allocate(this%signature_two(size(signature_two)))
    this%signature_two=signature_two
  end if


  allocate(this%counter(size(counter)))
  allocate(this%rank(size(rank)))
  allocate(this%perm_array(size(perm_array,1),size(perm_array,2),size(perm_array,3)))
  allocate(this%dist_vec(num_distances))
  allocate(this%dist_vec_permutations(size(dist_vec_permutations,1),size(dist_vec_permutations,2)))


  this%counter=counter
  this%rank=rank
  this%perm_array=perm_array
  this%dist_vec_permutations=dist_vec_permutations
  this%perm_number=1
  this%n_perms=dist_vec_n_perms
  this%initialised=.true.

end subroutine permutation_data_initialise

subroutine permutation_data_finalise(this)
  implicit none
  type(permutation_data_type) :: this

  if (.not. this%initialised) return

  if(allocated(this%counter)) deallocate(this%counter)
  if(allocated(this%rank)) deallocate(this%rank)
  if(allocated(this%perm_array)) deallocate(this%perm_array)
  if(allocated(this%dist_vec)) deallocate(this%dist_vec)
  this%initialised = .false.

end subroutine permutation_data_finalise
 
  subroutine add_combined_permutation (counter, perm_array, dist_vec_permutations,perm_number)
    implicit none
  ! this gets called by the subroutine next, it should receive a vector 'counter' from which it 
  ! figures out which permutations to combine. It then asks combine_perms to do so and 
  ! gets the dist_vec vector from do_swaps
    integer ::  i, num_distances, N, perm_number
    integer, dimension(:), intent(inout) :: counter
    integer, dimension(:), allocatable :: combo, next_perm, dist_vec
    integer, dimension(:,:,:), intent(in) :: perm_array
    integer, dimension(:,:) :: dist_vec_permutations

    N=size(perm_array,2)
    num_distances = N*(N-1)/2
    allocate(dist_vec(num_distances))
    allocate(combo(N))
    allocate(next_perm(N))

    combo = perm_array(1,:,counter(1))
    do i=1, size(counter)-1
      next_perm = perm_array(i+1,:,counter(i+1))
      combo = combine_perms(combo,next_perm)
    end do
    call do_swaps(combo, dist_vec)

    dist_vec_permutations(:,perm_number)=dist_vec
 
    deallocate(dist_vec)
    deallocate(combo)
    deallocate(next_perm)

  end subroutine add_combined_permutation
 
   recursive subroutine next(this, m)
    implicit none
    type(permutation_data_type), intent(inout) :: this
    integer :: m, num_groups

    num_groups = size(this%counter)

    if (m .gt. num_groups) then

      call add_combined_permutation(this%counter, this%perm_array, this%dist_vec_permutations, this%perm_number)
      this%perm_number=this%perm_number+1

    else
      do while (this%counter(m) .lt. this%rank(m))
        call next(this, m+1)
        this%counter(m+1:) = 1
        this%counter(m) = this%counter(m) + 1
      end do
      this%counter(m+1:) = 1
      call next(this,m+1)

    end if
  end subroutine next

  function num_group_perms(group)
    implicit none
    integer :: num_group_perms, n_members
    integer, dimension(:) :: group
    integer, dimension(8) :: factorial

    factorial = (/ 1,2,6,24,120,720,5040,40320 /)
    n_members = ceiling(log10(real(maxval(group))))
    num_group_perms = factorial(n_members)
    return
  end function num_group_perms

  ! This modified from Rosetta Code
  recursive subroutine update_matrix(std_perms,n_members,position,i_perm,perm_vec)
    implicit none
    integer :: n_members, value, position, i_perm
    integer, dimension(:,:) :: std_perms
    integer, dimension(:) :: perm_vec

    if (position > n_members) then
      std_perms(i_perm,:) = perm_vec
      i_perm=i_perm+1
    else
      do value = 1, n_members
        if (.not. any (perm_vec(:position - 1) == value)) then
          perm_vec(position)= value
          call update_matrix(std_perms,n_members,position+1,i_perm,perm_vec)
        end if
      end do
    end if
  end subroutine update_matrix

  function permute_atoms(atoms,group,N,max_rank)
    implicit none
    integer :: i, j, k, i_perm, n_members, num_perms, atoms_per_monomer
    integer, intent(IN) :: N, max_rank
    integer, dimension(N) :: atoms, group
    integer, dimension(:,:), allocatable :: permute_atoms, std_perms
    integer, dimension(:), allocatable :: group_vec, perm_vec, indices, offsets
    integer, dimension(1) ::  p, q, temp

    n_members = ceiling(log10(real(maxval(group))))
    num_perms = num_group_perms(group)
    allocate(group_vec(n_members))
    allocate(perm_vec(n_members))
    !allocate(indices(n_members))
    allocate(std_perms(num_perms,n_members))
    allocate(permute_atoms(max_rank,N))
    permute_atoms = 0

    if (num_perms .eq. 2) then
       !just a pair of equivalent atoms or monomers
       permute_atoms(1,:) = atoms
       permute_atoms(2,:) = atoms
       do i=1,count(group .lt. 10 .and. group .gt. 0)
          p = minloc(group, mask=group .ge. i)
          q = minloc(group, mask=group .gt. 10*i)
          !write(*,*) "equivalent pair"
          !write(*,'(2I3)') p,q 
          temp = permute_atoms(2,p(1))
          permute_atoms(2,p(1)) = permute_atoms(2,q(1))
          permute_atoms(2,q(1)) = temp(1)
       end do
    else
       !Permutations of groups of >2 atoms, no support for >2 monomers yet
       i_perm=1
       call update_matrix(std_perms,n_members,1,i_perm,perm_vec)


       if (.not. any(group .eq. 2)) then ! permutations of identical atoms
         allocate(indices(n_members))
         ! get indices of equivalent atoms
         do i=1,n_members
           temp =minloc(group, mask = group .ge. 10**(i-1))
           indices(i) = temp(1)
         end do

         do i=1,size(std_perms,1)
            perm_vec = std_perms(i,:)
            group_vec = indices(perm_vec)
            do j=1,n_members
               permute_atoms(i,indices(j)) = group_vec(j)
            end do
            do j=1,N
              if (permute_atoms(i,j) ==0) permute_atoms(i,j) = j
            end do
         end do

       else ! permutations of identical monomers
         allocate(offsets(n_members))
         atoms_per_monomer = maxval(group, mask = group .lt. 10)
!         allocate(indices(n_members*atoms_per_monomer))
        
         do i=1,n_members
           ! find the 1, 11, 111, etc. with which the monomer starts
           temp =minloc(group, mask = group .ge. 10**(i-1))
           offsets(i) = temp(1)-1
!!$           do j=1,atoms_per_monomer
!!$             indices(j+(i-1)*n_members) = j + offsets(i)
!!$           end do
         end do

         do i=1,size(std_perms,1)
            perm_vec = std_perms(i,:)
            group_vec = offsets(perm_vec)
            do j=1,n_members
              do k=1,atoms_per_monomer
                permute_atoms(i,k+offsets(j)) = k+group_vec(j)
              end do
            end do
            do j=1,N
              if (permute_atoms(i,j) ==0) permute_atoms(i,j) = j
            end do
         end do

       end if
    end if

    return

  end function permute_atoms

  function combine_perms(vec1,vec2)
  implicit none
  integer, dimension(:), intent(in) :: vec1, vec2
  integer, dimension(:), allocatable :: combine_perms
  integer :: j
  allocate(combine_perms(size(vec1)))

  if (size(vec1) /= size(vec2)) then
    write(*,*) "combine_perms received vectors of mismatched lengths"
    call exit(1)
  end if

  do j=1,size(vec1)
    combine_perms(j) =  vec1(vec2(j))
  end do
  return
  end function combine_perms

  subroutine do_swaps(atom_vec, dist_vec)
  implicit none
  integer :: N, start, finish, length, temp, j, i
  integer, dimension(:), intent(in) :: atom_vec
  integer, dimension(:) :: dist_vec 
  integer, dimension(1) :: i_vec
  integer, dimension(:), allocatable :: temp_vec, scratch_vec!, do_swaps
  integer, dimension(:,:), allocatable :: dist_mat, dist_mat_upper

  !initialise vector and matrix
  N = size(atom_vec)
  allocate(scratch_vec(N))
  allocate(temp_vec(N))
  do i=1,N
    temp_vec(i)=i
  end do

  do i=1,size(dist_vec)
    dist_vec(i)=i
  end do


  allocate(dist_mat(N,N))
  allocate(dist_mat_upper(N,N))
  dist_mat=0
  dist_mat_upper = 0

  start = 1
  do i=1,N
    finish=start + N-i
    dist_mat_upper(i,i+1:N) = dist_vec(start:finish-1)
    start = finish
  end do

  dist_mat = dist_mat_upper + transpose(dist_mat_upper)

  do while (any(temp_vec .ne. atom_vec))
    i_vec = minloc(temp_vec, temp_vec .ne. atom_vec)
    i=i_vec(1)
 
!  do i=1,N
!    if (temp_vec(i) .ne. atom_vec(i)) then
      ! keep track of swaps
      temp = temp_vec(i)
      temp_vec(i) = temp_vec(atom_vec(i))
      temp_vec(atom_vec(i)) = temp
      ! now swap in array - rows then columns
      scratch_vec = dist_mat(i,:)
      dist_mat(i,:) = dist_mat(atom_vec(i),:)
      dist_mat(atom_vec(i),:) = scratch_vec

      scratch_vec = dist_mat(:,i)
      dist_mat(:,i) = dist_mat(:,atom_vec(i))
      dist_mat(:,atom_vec(i)) = scratch_vec
!    end if
  end do
  !convert back into vector
  start = 1
  do i=1,N
    finish=start + N-i
    dist_vec(start:finish) = dist_mat(i,i+1:N)  
    start = finish
  end do
  
  deallocate(temp_vec)
  deallocate(scratch_vec)
  deallocate(dist_mat)
  deallocate(dist_mat_upper)
  return
  end subroutine do_swaps

end module permutation_maker_module
