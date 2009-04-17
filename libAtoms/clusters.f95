module clusters_module
  use Atoms_module
  use Table_module
  use ParamReader_module
  use system_module
  use atoms_module
  use dynamicalsystem_module
  use table_module


implicit none
private

integer, parameter, public :: &
   HYBRID_ACTIVE_MARK = 1, &
   HYBRID_BUFFER_MARK = 2, &
   HYBRID_TRANS_MARK = 3, &
   HYBRID_TERM_MARK = 4, &
   HYBRID_BUFFER_OUTER_LAYER_MARK = 5, &
   HYBRID_FIT_MARK = 6, &
   HYBRID_NO_MARK = 0


character(len=TABLE_STRING_LENGTH), parameter :: hybrid_mark_name(0:6) = &
  (/ "h_none    ", &
     "h_active  ", &
     "h_buffer  ", &
     "h_trans   ", &
     "h_term    ", &
     "h_outer_l ", &
     "h_fit     " /)

public :: create_cluster, create_cluster_hybrid_mark, create_hybrid_weights, bfs_grow, bfs_step, &
     multiple_images, discard_non_min_images, make_convex, create_embed_and_fit_lists, add_cut_hydrogens, & 
     construct_buffer, select_hysteretic_quantum_region

interface create_hybrid_weights
   module procedure create_hybrid_weights_args
   module procedure create_hybrid_weights_args_str
end interface

!% Grow a selection list by bond hopping.
interface bfs_grow
   module procedure bfs_grow_single
   module procedure bfs_grow_list
end interface

contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X   Cluster carving routines
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% On exit, 'list' will contain 'atom' (with shift '000') 
  !% plus the atoms within 'n' bonds hops of it.
  subroutine bfs_grow_single(this, list, atom, n, nneighb_only, min_images_only)
    type(Atoms), intent(in)  :: this
    type(Table), intent(out) :: list
    integer, intent(in) :: atom, n
    logical, optional, intent(in)::nneighb_only, min_images_only

    call append(list, (/atom, 0,0,0/)) ! Add atom with shift 000
    call bfs_grow_list(this, list, n, nneighb_only, min_images_only)

  end subroutine bfs_grow_single


  !% On exit, 'list' will have been grown by 'n' bond hops.
  subroutine bfs_grow_list(this, list, n, nneighb_only, min_images_only)
    type(Atoms), intent(in)   ::this
    type(Table), intent(inout)::list
    integer, intent(in) :: n
    logical, optional, intent(in)::nneighb_only, min_images_only
    type(Table)::tmplist
    integer::i

    do i=1,n
       call bfs_step(this, list, tmplist, nneighb_only, min_images_only)
       call append(list, tmplist)
    end do

    if (n >= 1) call finalise(tmplist)
  end subroutine bfs_grow_list


  !% Execute one Breadth-First-Search move on the atomic connectivity graph.
  subroutine bfs_step(this,input,output,nneighb_only, min_images_only)
    type(Atoms),        intent(in)      :: this  !% The atoms structure to perform the step on.
    type(Table),        intent(in)      :: input !% Table with intsize 4. First integer column is indices of atoms
                                                 !% already in the region, next 3 are shifts.
    type(Table),        intent(out)     :: output !% Table with intsize 4, containing the new atomic 
                                                  !% indices and shifts.

    logical, optional,  intent(in)      :: nneighb_only
    !% If present and true, sets whether only neighbours
    !% within the sum of the two respective covalent radii (multiplied by the atom's nneightol) are included,
    !% irrespective of the cutoff in the atoms structure
    !% (default is true).

    logical, optional, intent(in)       :: min_images_only 
    !% If true, there will be no repeated atomic indices in final list - only the
    !% minimum shift image of those found will be included. Default is false.

    !local
    logical                             :: do_nneighb_only, do_min_images_only
    integer                             :: i, j, n, m, jshift(3), ishift(3)

    integer :: n_i, keep_row(4), in_i, min_image
    integer, allocatable, dimension(:) :: repeats
    real(dp), allocatable, dimension(:) :: norm2shift

    if (.not.this%connect%initialised) &
         call system_abort('BFS_Step: Atomic structure has no connectivity data')

    do_nneighb_only = .true.
    if(present(nneighb_only)) do_nneighb_only = nneighb_only

    do_min_images_only = .false.
    if (present(min_images_only)) do_min_images_only = min_images_only

    if(input%intsize /= 4 .or. input%realsize /= 0) &
         call system_abort("bfs_step: input table must have intsize=4.")

    call table_allocate(output, 4, 0, 0, 0)

    ! Now go though the atomic indices
    do m = 1, input%N
       i = input%int(1,m)
       ishift = input%int(2:4,m)

       ! Loop over i's neighbours
       do n = 1, atoms_n_neighbours(this,i)
          j = atoms_neighbour(this,i,n,shift=jshift)

          ! Look at next neighbour if j with correct shift is already in the cluster
          ! Must check input AND output tables
          if (find(input,(/j,ishift+jshift/)) > 0) cycle
          if (find(output,(/j,ishift+jshift/)) > 0) cycle

          if (do_nneighb_only .and. .not. is_nearest_neighbour(this, i, n)) cycle

          ! Everything checks out ok, so add j to the output table
          ! with correct shift
          call append(output,(/j,ishift+jshift/)) 
       end do
    end do

    if (do_min_images_only) then

       ! If there are repeats of any atomic index, 
       ! we want to keep the one with smallest norm2(shift)

       ! must check input and output

       n = 1
       do while (n <= output%N)

          ! How many occurances in the output list?
          n_i = count(int_part(output,1) == output%int(1,n))

          ! Is there one in the input list as well?
          in_i = find_in_array(int_part(input,1),output%int(1,n))

          ! If there's only one, and we've not got one already then move on
          if (n_i == 1 .and. in_i == 0) then 
             n = n + 1
             cycle
          end if

          ! otherwise, things are more complicated...
          ! we want to keep the one with the smallest shift, bearing
          ! in mind that it could well be the one in the input list

          allocate(repeats(n_i), norm2shift(n_i))

          ! Get indices of repeats of this atomic index
          repeats = pack((/ (j, j=1,output%N) /), &
               int_part(output,1) == output%int(1,n))

          if (in_i /= 0) then
             ! atom is in input list, remove all new occurances
             call delete_multiple(output, repeats)
          else
             ! Find row with minimum norm2(shift)
             norm2shift = norm2(real(output%int(2:4,repeats),dp),1)
             min_image = repeats(minloc(norm2shift,dim=1))
             keep_row = output%int(:,min_image)

             ! keep the minimum image
             call delete_multiple(output, &
                  pack(repeats, repeats /= min_image))
          end if

          ! don't increment n, since delete_multiple copies items from
          ! end of list over deleted items, so we need to retest

          deallocate(repeats, norm2shift)

       end do

    end if

  end subroutine bfs_step


  !% Check if the list of indices and shifts 'list' contains
  !% any repeated atomic indices.
  function multiple_images(list)
    type(Table), intent(in) :: list
    logical :: multiple_images


    integer :: i
    multiple_images = .false.

    if (list%N == 0) return

    do i = 1, list%N
       if (count(int_part(list,1) == list%int(1,i)) /= 1) then
          multiple_images = .true.
          return
       end if
    end do

  end function multiple_images

  !% Given an input list with 4 integer columns for atomic indices
  !% and shifts, keep only the minimum images of each atom
  subroutine discard_non_min_images(list)
    type(Table), intent(inout) :: list

    integer :: n, n_i, j, min_image, keep_row(4)
    integer, allocatable, dimension(:) :: repeats
    real(dp), allocatable, dimension(:) :: norm2shift

    n = 1
    do while (n <= list%N)

       ! How many occurances in the list list?
       n_i = count(int_part(list,1) == list%int(1,n))

       ! If there's only one, and we've not got one already then move on
       if (n_i == 1) then 
          n = n + 1
          cycle
       end if

       ! otherwise, things are more complicated...
       ! we want to keep the one with the smallest shift

       allocate(repeats(n_i), norm2shift(n_i))

       ! Get indices of repeats of this atomic index
       repeats = pack((/ (j, j=1,list%N) /), &
            int_part(list,1) == list%int(1,n))

       ! Find row with minimum norm2(shift)
       norm2shift = norm2(real(list%int(2:4,repeats),dp),1)
       min_image = repeats(minloc(norm2shift,dim=1))
       keep_row = list%int(:,min_image)

       ! keep the minimum image
       call delete_multiple(list, &
            pack(repeats, repeats /= min_image))

       ! don't increment n, since delete_multiple copies items from
       ! end of list over deleted items, so we need to retest

       deallocate(repeats, norm2shift)

    end do

  end subroutine discard_non_min_images


  !% Add atoms to 'list' to make the selection region convex, i.e. if $i$ and
  !% $j$ are nearest neighbours, with $i$ in the list and not $j$ then $j$ will be added
  !% if more than half its nearest neighbours are in the list.
  subroutine make_convex(this, list)
    type(Atoms), intent(in) :: this
    type(Table), intent(inout) :: list
    type(Table)::tmplist
    do while(make_convex_step(this, list, tmplist) /= 0)
       call append(list, tmplist)
    end do
    call finalise(tmplist)
  end subroutine make_convex

  !% OMIT
  ! Given an input list, what needs to be added to make the selection region convex?
  function make_convex_step(this, input, output) result(newatoms)
    type(Atoms), intent(in) :: this
    type(Table), intent(in) :: input
    type(Table), intent(out) :: output
    integer::newatoms

    integer :: n, i, j, k, p, m, n_in, nn, ishift(3), jshift(3), kshift(3)
    real(dp) :: r_ij, r_kj

    ! Check table size
    if(input%intsize /= 4 .or. input%realsize /= 0) &
         call system_abort("make_convex_step: input table must have intsize=4")

    if(input%intsize /= 4 .or. input%realsize /= 0) &
         call system_abort("bfs_step: input table must have intsize=4.")

    call table_allocate(output, 4, 0, 0, 0)

    do n=1,input%N
       i = input%int(1,n)
       ishift = input%int(2:4,n)

       !Loop over neighbours
       do m = 1, atoms_n_neighbours(this,i)
          j = atoms_neighbour(this,i,m, r_ij,shift=jshift)

          ! Look for nearest neighbours of i not in input list
          if (find(input,(/j,ishift+jshift/)) == 0 .and. is_nearest_neighbour(this, i, m)) then

             n_in = 0
             nn = 0
             ! Count number of nearest neighbours of j, and how many
             ! of them are in input list
             do p = 1, atoms_n_neighbours(this,j)
                k = atoms_neighbour(this,j,p, r_kj, shift=kshift)
                if (is_nearest_neighbour(this, j, p)) then
                   nn = nn + 1
                   if (find(input,(/k,ishift+jshift+kshift/)) /= 0) n_in = n_in + 1
                end if
             end do

             !If more than half of  j's nearest neighbours are in then add it to output
             if (real(n_in,dp)/real(nn,dp) > 0.5_dp) &
                  call append(output, (/j,ishift+jshift/))

          end if

       end do
    end do

    newatoms = output%N

  end function  make_convex_step


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Create Cluster:
  !% Returns an Atoms object (cluster) which contains the atoms whose
  !% indices are given in atomlist, possibly with some extras for consistency,
  !% and optionally terminated with Hydrogens.
  !% The output cluster contains all properties of the initial atoms object, and
  !% some additional columns, which are:
  !% "index" : index of the cluster atoms into the initial atoms object.
  !% "termindex": nonzero for termination atoms, and is an index into the cluster atoms specifiying which atom
  !% is being terminated, it is used in collecting the forces.
  !% "rescale" : a real number which for nontermination atoms is 1.0, for termination atoms records
  !% the scaling applied to termination bond lengths
  !% "shift" : the shift of each atom
  !% 
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function create_cluster(this, atomlist, terminate, periodic, same_lattice, even_hydrogens, vacuum_size, &
       cut_bonds, allow_cluster_modification) result(cluster)
    type(Atoms),               intent(in)    :: this           !% Input Atoms object
    type(Table),               intent(in)    :: atomlist       !% List of atoms to include in cluster. This should be
                                                               !% either 1 column with indices, or 4 columns with indices
                                                               !% and shifts relative to first atom in list.
    logical,     optional,     intent(in)    :: terminate      !% Should Hydrogens be added to cut bonds (default .true.)
    logical,     optional,     intent(in)    :: periodic(3)    !% Should cluster be periodic in each direction.
                                                               !% Default is '(/false,false,false/)'.
                                                               !% Number of true entries must be zero, one or three.
    logical,	 optional,     intent(in)    :: same_lattice   !% Should lattice be left the same, overrides periodic variable
                                                               !% Default false
    logical,     optional,     intent(in)    :: even_hydrogens !% If true, then cluster will always be terminated with
                                                               !% an even number of hydrogens to prevent an imbalance
                                                               !% of spin up and spin down electrons. If a hydrogen has to
                                                               !% be removed it will be taken from an atom with as many 
                                                               !% termination hydrogens as possible.
    real(dp),    optional                    :: vacuum_size    !% Amount of vacuum around clusters, in \AA{}.
                                                               !% Default is 10 \AA{}.
    type(Table), optional,     intent(out)   :: cut_bonds      !% Return a list of the bonds cut when making
                                                               !% the cluster. Table with 8 'int' columns,
                                                               !% for $i$, $j$, 'shift_i' and 'shift_j'.
                                                               !% for the atom indices at each end of the cut bonds.
    logical,     optional,     intent(in)    :: allow_cluster_modification  !% if false, don't try to fix cluster surface
    type(Atoms)                              :: cluster      ! this is the output

    type(Table)                              :: cluster_temp,  n_term
    integer                                  :: i, j, k, m, n, p, lookup(3)

    real(dp),    dimension(3)                :: diff_ik
    real(dp),    dimension(3)                :: dhat_ij, dhat_jk, H1, H2, maxlen, sep
    real(dp),    dimension(3)                :: lat_maxlen, lat_sep
    real(dp)                                 :: r_ij, r_jk, cluster_vacuum, rescale
    logical                                  :: all_in
    logical                                  :: do_terminate, do_periodic(3), do_even_hydrogens, do_same_lattice
    integer                                  :: ishift(3), jshift(3), kshift(3), oldN, most_hydrogens
    logical                                  :: atom_mask(6)
    logical allow_cluster_mod
    integer, allocatable, dimension(:,:)     :: periodic_shift

    integer, pointer :: hybrid_mark(:)
    ! optional defaults

    allow_cluster_mod = optional_default(.true., allow_cluster_modification)
    cluster_vacuum    = optional_default(10.0_dp, vacuum_size)
    do_terminate      = optional_default(.true., terminate)
    do_same_lattice   = optional_default(.false., same_lattice)
    do_even_hydrogens = optional_default(.false., even_hydrogens)
    
    do_periodic = (/.false.,.false.,.false./)
    if (present(periodic)) do_periodic = periodic

    if (.not. assign_pointer(this, "hybrid_mark", hybrid_mark)) &
      call system_abort("create_cluster impossible failure to assing hybrid_mark pointer")

    ! 
    ! Validate arguments
    !

    ! check for consistency in optional arguments

    if (.not. (count(do_periodic) == 0 .or. count(do_periodic) == 1 .or. count(do_periodic) == 3)) &
         call system_abort('count(periodic) must be zero, one or three.')

    if (do_same_lattice) do_periodic = .true.

    if (any(do_periodic) .and. multiple_images(atomlist)) &
         call system_abort("create_cluster: can't make a periodic cluster since atomlist contains repeats")

    ! check for empty list

    if(atomlist%N == 0) then
       call print('create_cluster: empty atomlist', NORMAL)
       return
    end if

    call print('create_cluster: Entering create_cluster', NERD)

    if(.not.(atomlist%intsize == 1 .or. atomlist%intsize == 4) .or. atomlist%realsize /= 0) &
         call system_abort("create_cluster: atomlist table must have intsize=1 or 4 and realsize=0.")


    ! Cluster_temp is a temporary, extensible storage for the cluster
    ! It stores atomic indices and shifts (4 ints)
    ! atomic number (1 int)
    ! termination index (1 int): for termination atoms, which atom is being terminated?
    ! and atomic positions (3 reals)
    ! It's length will be at least atomlist%N

    call print('create_cluster: Creating temporary cluster table', NERD)
    call table_allocate(cluster_temp,6,4,1,0,atomlist%N)


    ! First, put all the marked atoms into cluster_temp, storing their positions and shifts
    call print('create_cluster: Adding specified atoms to the cluster', NERD)
    do i = 1, atomlist%N
       if(atomlist%intsize == 4) then
          ! we have shifts
          ishift = atomlist%int(2:4,i)
       else
          ! no incoming shifts
          ishift = (/0,0,0/)
       end if
       call append(cluster_temp, (/atomlist%int(1,i),ishift,this%Z(atomlist%int(1,i)),0/),&
	    (/this%pos(:,atomlist%int(1,i)),1.0_dp/), (/ hybrid_mark_name(hybrid_mark(atomlist%int(1,i))) /) )
    end do

    call print("create_cluster: cluster list:", NERD)
    call print(cluster_temp, NERD)

    ! Next, check for various gotchas

    ! this mask is used to match with atoms already in the cluster_temp table
    ! if we are periodic in a direction, we don't care about the shifts in that direction when matching
    atom_mask = (/.true.,.not.do_periodic, .true., .true./)

    if(allow_cluster_mod) then
       ! Gotcha 1: Hollow sections
       ! OUT and IN refers to the list in cluster_temp
       ! Look at the OUT nearest neighbours of IN atoms. If all the nearest neighbours of the OUT
       ! atom are IN, then make the OUT atom IN.
       n = 1
       ! Loop over cluster atoms (including ones that may get added in this loop)
       call print('create_cluster: Checking for hollow sections', NERD)
       do while(n <= cluster_temp%N)
          i = cluster_temp%int(1,n)
          ishift = cluster_temp%int(2:4,n)
          call print('create_cluster: i = '//i//' ['//ishift//'] Looping over '//atoms_n_neighbours(this,i)//' neighbours...',ANAL)

          !Loop over neighbours
          do m = 1, atoms_n_neighbours(this,i)
             j = atoms_neighbour(this,i,m, r_ij,shift=jshift)

             if (find(cluster_temp,(/j,ishift+jshift,this%Z(j),0/), atom_mask) == 0 .and. &
                  is_nearest_neighbour(this, i, m)) then

                call print('create_cluster:   checking j = '//j//" ["//jshift//"]",ANAL)

                ! We have an OUT nearest neighbour, loop over it's nearest neighbours to see if they
                ! are all IN

                all_in = .true.
                do p = 1, atoms_n_neighbours(this,j)
                   k = atoms_neighbour(this,j,p, shift=kshift)
                   if (find(cluster_temp,(/k,ishift+jshift+kshift,this%Z(k),0/), atom_mask) == 0 .and. &
                        is_nearest_neighbour(this, j, p)) then
                      all_in = .false.
                      exit
                   end if
                end do

                !If all j's nearest neighbours are IN then add it
                if (all_in) then
                   call append(cluster_temp, (/j,ishift+jshift,this%Z(j),0/), (/this%pos(:,j), 1.0_dp/), (/ "hollow    "/) )
                   if(current_verbosity() .ge. NERD) then
		      call print('create_cluster:  Adding atom ' //j//' ['//(ishift+jshift)//'] to cluster. Atoms = ' // cluster_temp%N, NERD)
                   end if
                end if

             end if

          end do
          n = n + 1
       end do

       call print('create_cluster: Finished checking for hollow sections',NERD)
       call print("create_cluster: cluster list:", NERD)
       call print(cluster_temp, NERD)

       ! Gotcha 2: In-out-in structures
       ! Find cases where two IN atoms have a common
       ! OUT nearest neighbour, and see if termination would cause the hydrogen
       ! atoms to be too close together. If so, include the OUT nearest neighbour
       ! in the cluster

       call print('create_cluster: Checking for termination clashes', NERD)

       !Loop over atoms in the cluster
       n = 1
       if(do_terminate) then
         do while (n <= cluster_temp%N)
            i = cluster_temp%int(1,n)     ! index of atom in the cluster
            ishift = cluster_temp%int(2:4,n)
            call print('create_cluster: i = '//i//'. Looping over '//atoms_n_neighbours(this,i)//' neighbours...',ANAL)
            !Loop over atom i's neighbours
            do m = 1, atoms_n_neighbours(this,i)
               j = atoms_neighbour(this,i,m, shift=jshift, diff=dhat_ij, distance=r_ij)
               dhat_ij = dhat_ij/r_ij
        
               !If j is IN the cluster, or not a nearest neighbour then try the next neighbour
        
               if(find(cluster_temp,(/j,ishift+jshift,this%Z(j),0/), atom_mask) /= 0) then
                  call print('create_cluster:   j = '//j//" ["//jshift//"] is in cluster",ANAL)
                  cycle
               end if
               if(.not. is_nearest_neighbour(this,i, m)) then
                  call print('create_cluster:   j = '//j//" ["//jshift//"] not nearest neighbour",ANAL)
                  cycle
               end if
        
               ! So j is an OUT nearest neighbour of i.
               call print('create_cluster:   checking j = '//j//" ["//jshift//"]",ANAL)
        
               !Determine the position of the would-be hydrogen along i--j
               call print('create_cluster:  Finding i--j hydrogen position',ANAL)
        
               H1 = this%pos(:,i) + (this%lattice .mult. (ishift)) + &
                    termination_bond_rescale(this%Z(i), this%Z(j)) * r_ij * dhat_ij
        
               !Do a loop over j's nearest neighbours
               call print('create_cluster:  Looping over '//atoms_n_neighbours(this,j)//' neighbours of j',&
                    ANAL)
        
               do p = 1, atoms_n_neighbours(this,j)
        
                  k = atoms_neighbour(this,j,p, shift=kshift, diff=dhat_jk, distance=r_jk)
                  dhat_jk = dhat_jk/r_jk
        
                  !If k is OUT of the cluster or k == i or it is not a nearest neighbour of j
                  !then try the next neighbour
        
                  if(find(cluster_temp,(/k,ishift+jshift+kshift,this%Z(k),0/), atom_mask) == 0) then
                     call print('create_cluster:   k = '//k//" ["//kshift//"] not in cluster",ANAL)
                     cycle
                  end if
                  if(k == i .and. all( jshift+kshift == 0 )) cycle
                  if(.not. is_nearest_neighbour(this,j, p)) then
                     call print('create_cluster:   k = '//k//" ["//kshift//"] not nearest neighbour",ANAL)
                     cycle
                  end if
        
                  call print('create_cluster: testing k = '//k//" ["//kshift//"]", ANAL)
                  !Determine the position of the would-be hydrogen along k--j
                  call print('create_cluster:   Finding k--j hydrogen position',ANAL)
        
                  diff_ik = r_ij * dhat_ij + r_jk * dhat_jk
        
                  H2 = this%pos(:,i) + (this%lattice .mult. (ishift)) + diff_ik - &
                       termination_bond_rescale(this%Z(k), this%Z(j)) * r_jk * dhat_jk
                  call print('create_cluster:   Checking i--k distance and hydrogen distance',ANAL)
        
                  call print("create_cluster: proposed hydrogen positions:", ANAL)
                  call print(H1, ANAL)
                  call print(H2, ANAL)
                  call print("create_cluster: hydrogen distance would be "//norm(H1-H2), ANAL)
                  ! If i and k are nearest neighbours, or the terminating hydrogens would be very close, then
                  ! include j in the cluster. The H--H checking is conservative, hence the extra factor of 1.2
        
        
                  if ((norm(diff_ik) < bond_length(this%Z(i),this%Z(k))*this%nneightol) .or. &
                       (norm(H1-H2) < bond_length(1,1)*this%nneightol*1.2_dp) ) then
        
                     call append(cluster_temp,(/j,ishift+jshift,this%Z(j),0/),(/this%pos(:,j),1.0_dp/), (/ "clash     "/) )
                     call print('create_cluster:  Atom '//j//' added to cluster. Atoms = '//cluster_temp%N, &
                          NERD)
                     ! j is now included in the cluster, so we can exit this do loop (over p)
                     exit
                  end if
        
               end do
        
            end do
            n = n + 1
         end do
       endif

       call print('create_cluster: Finished checking for termination clashes',NERD)
       call print("create_cluster: cluster list:", NERD)
       call print(cluster_temp, NERD)
    end if ! allow_cluster_mod

    !So now cluster_temp contains all the atoms that are going to be in the cluster.
    !If do_terminate is set, we need to add terminating hydrogens along nearest neighbour bonds
    if (do_terminate) then
       call print('create_cluster: Terminating cluster with hydrogens',NERD)

       call table_allocate(n_term, 5, 0, 0, 0)
       oldN = cluster_temp%N

       !Loop over atoms in the cluster
       do n = 1, oldN

          i = cluster_temp%int(1,n)
          ishift = cluster_temp%int(2:4,n)
          !Loop over atom i's neighbours
          do m = 1, atoms_n_neighbours(this,i)

             j = atoms_neighbour(this,i,m, r_ij, diff=dhat_ij, shift=jshift)
             dhat_ij = dhat_ij / r_ij

             if (find(cluster_temp,(/j,ishift+jshift,this%Z(j),0/), atom_mask) == 0 .and. &
                  is_nearest_neighbour(this, i, m)) then

                ! If j is an OUT atom, and it is close enough, put a terminating hydrogen
                ! at the scaled distance between i and j
               
                rescale = termination_bond_rescale(this%Z(i), this%Z(j))
                H1 = this%pos(:,i) + rescale * r_ij * dhat_ij

                ! Label term atom with indices into original atoms structure.
                ! j is atom it's generated from and n is index into cluster table of atom it's attached to
                call append(cluster_temp,(/j,ishift,1,n/),(/H1, rescale/), (/ "term      " /)) 

                ! Keep track of how many termination atoms each cluster atom has
                p = find_in_array(int_part(n_term,(/1,2,3,4/)),(/n,ishift/))
                if (p == 0) then
                   call append(n_term, (/n,ishift,1/))
                else
                   n_term%int(5,p) = n_term%int(5,p) + 1
                end if

                ! optionally keep a record of the bonds that we have cut
                if (present(cut_bonds)) &
                     call append(cut_bonds, (/i,j,ishift,jshift/))

                if(current_verbosity() .ge. NERD) then
                   write(line,'(a,i0,a,i0,a)')'create_cluster: Replacing bond ',i,'--',j,' with hydrogen'
                   call print(line, NERD)
                end if
             end if
          end do

       end do

       ! Do we need to remove a hydrogen atom to ensure equal n_up and n_down electrons?
       if (do_even_hydrogens .and. mod(count(int_part(cluster_temp,5) == 1),2) == 1) then

          ! Find first atom with a maximal number of terminating hydrogens
          most_hydrogens = maxloc(int_part(n_term,5),dim=1)
          n = n_term%int(1,most_hydrogens)
          ishift = n_term%int(2:4,most_hydrogens)

          ! Loop over termination atoms
          do j=oldN,cluster_temp%N
             ! Remove first H atom attached to atom i
             if (all(cluster_temp%int(2:6,j) == (/ishift,1,n/))) then
                call delete(cluster_temp, j)
                call print('create_cluster: removed one of atom '//i//"'s "//maxval(int_part(n_term,5))// &
                     ' terminating hydrogens to zero total spin', VERBOSE)
                exit 
             end if
          end do
       end if

       call finalise(n_term)

       call print('create_cluster: Finished terminating cluster',NERD)

       if (allocated(periodic_shift)) deallocate(periodic_shift)
    end if

    !Now turn the cluster_temp table into an atoms structure
    call print('create_cluster: Copying atomic data to output object',NERD)
    call print('List of atoms in cluster:', NERD)
    call print(int_part(cluster_temp,1), NERD)


    ! first pick up an atoms structure with the right number of atoms and copy the properties
    call select(cluster, this, list=int_part(cluster_temp,1))
    ! then reset the positions species and Z (latter two needed because termination atoms have Z=1)
    ! unfold the positions to real positions using the stored shifts, this is neede because
    ! next we might make the unit cell much smaller
    do i=1,cluster_temp%N
       cluster%pos(:,i) = cluster_temp%real(1:3, i)+(this%lattice .mult. cluster_temp%int(2:4, i))
       cluster%Z(i) = cluster_temp%int(5,i)
       cluster%species(i) = ElementName(cluster_temp%int(5,i))
    end do
    ! add properties to cluster
    call add_property(cluster, 'index', int_part(cluster_temp,1))
    call add_property(cluster, 'shift', 0, n_cols=3, lookup=lookup)
    cluster%data%int(lookup(2):lookup(3),1:cluster%N) = cluster_temp%int(2:4,1:cluster_temp%N)
    call add_property(cluster, 'termindex', int_part(cluster_temp,6))
    call add_property(cluster, 'rescale', real_part(cluster_temp,4))
    call add_property(cluster,"cluster_ident", cluster_temp%str(1,:))

    ! Find smallest bounding box for cluster
    ! Find boxes aligned with xyz (maxlen) and with a1 a2 a3 (lat_maxlen)
    maxlen = 0.0_dp
    lat_maxlen = 0.0_dp
    do i=1,cluster%N
       do j=1,cluster%N
          sep = cluster%pos(:,i)-cluster%pos(:,j)
	  lat_sep = cluster%g .mult. sep
          do k=1,3
             if (abs(sep(k)) > maxlen(k)) maxlen(k) = abs(sep(k))
	     if (abs(lat_sep(k)) > lat_maxlen(k)) lat_maxlen(k) = abs(lat_sep(k))
          end do
       end do
    end do

    ! renormalize lat_maxlen to real dist units
    do k=1,3
      lat_maxlen(k) = lat_maxlen(k) * norm(cluster%lattice(:,k))
    end do

    ! Round up maxlen to be divisible by 3 A, so that it does not fluctuate too much
    forall (k=1:3) maxlen(k) = 3.0_dp*ceiling(maxlen(k)/3.0_dp)
    forall (k=1:3) lat_maxlen(k) = 3.0_dp*ceiling(lat_maxlen(k)/3.0_dp)

    ! vacuum pad cluster (if didn't set same_lattice)
    ! if not periodic at all, just do vacuum padding
    ! if periodic along some dir, keep supercell vector directions, and set
    !    extent in each direction to lesser of cluster extent + vacuum or original extent
    if (do_same_lattice) then
      cluster%lattice = this%lattice
    else
      if (any(do_periodic)) then
	do k=1,3
	  if (do_periodic(k)) then
	    if (lat_maxlen(k)+cluster_vacuum >= norm(this%lattice(:,k))) then
	      cluster%lattice(:,k) = this%lattice(:,k)
	    else
	      cluster%lattice(:,k) = (lat_maxlen(k)+cluster_vacuum)*this%lattice(:,k)/norm(this%lattice(:,k))
	    endif
	  else
	    cluster%lattice(:,k) = (lat_maxlen(k)+cluster_vacuum)*this%lattice(:,k)/norm(this%lattice(:,k))
	  endif
	end do
      else
	cluster%lattice = 0.0_dp
	do k=1,3
	  cluster%lattice(k,k) = maxlen(k) + cluster_vacuum
	end do
      endif
    endif

    call atoms_set_lattice(cluster, cluster%lattice)

    ! Remap positions so any image atoms end up inside the cell
    call map_into_cell(cluster)

    write (line, '(a,f10.3,f10.3,f10.3)') &
         'create_cluster: Cluster dimensions are ', cluster%lattice(1,1), &
         cluster%lattice(2,2), cluster%lattice(3,3)
    call print(line, VERBOSE)

    !Deallocate the table to tidy up
    call print('create_cluster: Freeing temporary storage',NERD)
    call finalise(cluster_temp)

    call print ('create_cluster: carved cluster with '//cluster%N//' atoms', VERBOSE)

    call print ('Exiting create_cluster', NERD)

  end function create_cluster

  !% Create a cluster using the 'hybrid_mark' property and options in 'args_str'.
  !% All atoms that are marked with anything other than 'HYBRID_NO_MARK' will
  !% be included in the cluster; this includes active, transition and buffer
  !% atoms. 
  function create_cluster_hybrid_mark(at, args_str) result(cluster)
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in) :: args_str
    type(Atoms) :: cluster

    type(Dictionary) :: params
    logical :: terminate, periodic_x, periodic_y, periodic_z, &
       even_hydrogens, do_calc_connect, do_periodic(3), cluster_nneighb_only, cluster_allow_modification, &
       do_rescale_r, print_clusters, randomise_buffer, in_outer_layer
    real(dp) :: cluster_vacuum, r_scale, r, r_min, centre(3)
    type(Table) :: cluster_list, outer_layer, currentlist, nextlist, activelist, bufferlist
    integer :: i, j, jj, n, first_active, old_n, n_cluster, shift(3)
    integer, pointer :: hybrid_mark(:), cluster_index(:), cluster_hybrid_mark(:)
    integer, allocatable, dimension(:) :: uniqed, tmp_index

    type(Inoutput)                    :: clusterfile
    character(len=255)                :: clusterfilename

#ifdef _MPI
    integer::mpi_size, mpi_rank, error
    include "mpif.h"
    integer :: mpi_force_size
    real(dp), allocatable, dimension(:)  :: mpi_force

    call get_mpi_size_rank(MPI_COMM_WORLD, mpi_size, mpi_rank)
#endif _MPI

    call print('create_cluster_hybrid_mark got args_str "'//trim(args_str)//'"', VERBOSE)

    call initialise(params)
    call param_register(params, 'terminate', 'T', terminate)
    call param_register(params, 'cluster_periodic_x', 'F', periodic_x)
    call param_register(params, 'cluster_periodic_y', 'F', periodic_y)
    call param_register(params, 'cluster_periodic_z', 'F', periodic_z)
    call param_register(params, 'even_hydrogens', 'F', even_hydrogens)
    call param_register(params, 'cluster_calc_connect', 'F', do_calc_connect)
    call param_register(params, 'cluster_nneighb_only', 'T', cluster_nneighb_only)
    call param_register(params, 'cluster_vacuum', '10.0', cluster_vacuum)
    call param_register(params, 'cluster_allow_modificiation', 'T', cluster_allow_modification)
    call param_register(params, 'do_rescale_r', 'F', do_rescale_r)
    call param_register(params, 'r_scale', '1.0', r_scale)
    call param_register(params, 'randomise_buffer', 'T', randomise_buffer)
    call param_register(params, 'print_clusters', 'F', print_clusters)
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.) ) &
      call system_abort("create_cluster_hybrid_mark failed to parse args_str='"//trim(args_str)//"'")
    call finalise(params)

    do_periodic = (/periodic_x,periodic_y,periodic_z/)

    if (.not. has_property(at, 'hybrid_mark')) &
         call system_abort('create_cluster_hybrid_mark: atoms structure has no "hybrid_mark" property')
    
    if (.not. assign_pointer(at, 'hybrid_mark', hybrid_mark)) &
         call system_abort('create_cluster_hybrid_mark passed atoms structure with no hybrid_mark property')


    ! Calculate centre of cluster
    call allocate(cluster_list, 1,0,0,0)
    call append(cluster_list, find(hybrid_mark /= HYBRID_NO_MARK))
    do i=1,cluster_list%N
       centre = centre + at%pos(:,cluster_list%int(1,i))
    end do
    centre = centre / cluster_list%N
!!$    call print('centre = '//centre)

    ! Find atom closest to centre of cluster, using min image convention  
    r_min = huge(1.0_dp)
    do i=1,cluster_list%N
       r = distance_min_image(at, centre, cluster_list%int(1,i))
!!$       call print(cluster_list%int(1,i)//' '//r)
       if (r < r_min) then
          first_active = cluster_list%int(1,i)
          r_min = r
       end if
    end do

    n_cluster = cluster_list%N
    call wipe(cluster_list)
    call allocate(cluster_list, 4,0,0,0)

    ! Add first marked atom to cluster_list. shifts will be relative to this atom
    call print('Growing cluster starting from atom '//first_active//', n_cluster='//n_cluster, VERBOSE)
    call append(cluster_list, (/first_active,0,0,0/))
    call append(currentlist, cluster_list)

    ! Add other active atoms using bond hopping from the central cluster atom
    ! to find the other cluster atoms and hence to determine the correct 
    ! periodic shifts. 
    !
    ! This will fail if marked atoms do not form a single connected cluster
    old_n = cluster_list%N
    do 
       call BFS_step(at, currentlist, nextlist, nneighb_only = cluster_nneighb_only, min_images_only = any(do_periodic))
       do j=1,nextlist%N
          jj = nextlist%int(1,j)
          shift = nextlist%int(2:4,j)
          if (hybrid_mark(jj) /= HYBRID_NO_MARK) &
               call append(cluster_list, nextlist%int(:,j))
       end do
       call append(currentlist, nextlist)

       ! check exit condition
       allocate(tmp_index(cluster_list%N))
       tmp_index = int_part(cluster_list,1)
       call sort_array(tmp_index)
       call uniq(tmp_index, uniqed)
       call print('cluster hopping: got '//cluster_list%N//' atoms, of which '//size(uniqed)//' are unique.', VERBOSE)
       if (size(uniqed) == n_cluster) exit !got them all
       deallocate(uniqed, tmp_index)
       
       ! check that cluster is still growing
       if (cluster_list%N == old_n) &
            call system_abort('create_cluster_hybrid_mark: cluster stopped growing before all marked atoms found - check for split QM region')
       old_n = cluster_list%N
    end do
    deallocate(tmp_index, uniqed)

    ! partition cluster_list so that active atoms come first
    do i=1,cluster_list%N
       if (hybrid_mark(cluster_list%int(1,i)) == HYBRID_ACTIVE_MARK) then
          call append(activelist, cluster_list%int(:,i))
       else
          call append(bufferlist, cluster_list%int(:,i))
       end if
    end do

    call wipe(cluster_list)
    call append(cluster_list, activelist)
    call append(cluster_list, bufferlist)
    call finalise(activelist)
    call finalise(bufferlist)
    
    cluster = create_cluster(at, cluster_list, terminate=terminate, &
         periodic=do_periodic, even_hydrogens=even_hydrogens, &
         vacuum_size=cluster_vacuum, allow_cluster_modification=cluster_allow_modification)

    ! reassign pointers
    if (.not. assign_pointer(at, 'hybrid_mark', hybrid_mark)) &
         call system_abort('cannot reassign hybrid_mark property')
    
    ! rescale cluster positions and lattice 
    if (do_rescale_r) then
       cluster%pos = r_scale * cluster%pos
       call set_lattice(cluster, r_scale * cluster%lattice)
    end if

    if (randomise_buffer .and. .not. any(hybrid_mark == HYBRID_BUFFER_OUTER_LAYER_MARK)) &
         do_calc_connect = .true.

    if (do_calc_connect) then
       call print('create_cluster_hybrid_mark: doing calc_connect', VERBOSE)
       ! Does QM force model need connectivity information?
       if (at%use_uniform_cutoff) then
          call atoms_set_cutoff(cluster, at%cutoff)
       else
          call atoms_set_cutoff_factor(cluster, at%cutoff)
       end if
       call calc_connect(cluster)
    end if

    if (randomise_buffer) then
       ! If outer buffer layer is not marked do so now. This will be
       ! the case if we're using hysteretic_buffer selection option.
       ! In this case we consider any atoms connected to terminating
       ! hydrogens to be in the outer layer - obviously this breaks down 
       ! if we're not terminating so we abort if that's the case.
       if (.not. assign_pointer(cluster, 'index', cluster_index)) &
            call system_abort('create_cluster_hybrid_mark: cluster is missing index property')

       if (.not. any(hybrid_mark == HYBRID_BUFFER_OUTER_LAYER_MARK)) then
          if (.not. terminate) call system_abort('cannot determine which buffer atoms to randomise if terminate=F and hysteretic_buffer=T')

          if (.not. assign_pointer(cluster, 'hybrid_mark', cluster_hybrid_mark)) &
               call system_abort('hybrid_mark property not found in cluster')

          do i=1,cluster_list%N
             if (hybrid_mark(cluster_list%int(1,i)) /= HYBRID_BUFFER_MARK) cycle
             in_outer_layer = .false.
             do n=1,atoms_n_neighbours(cluster,i)
                if (.not. is_nearest_neighbour(cluster,i,n)) cycle
                j = atoms_neighbour(cluster,i,n)
                if (j > cluster_list%N .and. cluster%Z(j) == 1) then
                   in_outer_layer = .true.
                   exit
                end if
             end do
             if (in_outer_layer) then
                hybrid_mark(cluster_list%int(1,i)) = HYBRID_BUFFER_OUTER_LAYER_MARK
                cluster_hybrid_mark(i) = HYBRID_BUFFER_OUTER_LAYER_MARK
             end if
          end do
       end if

      if (any(hybrid_mark == HYBRID_BUFFER_OUTER_LAYER_MARK)) then
	 ! Slightly randomise the positions of the outermost layer
	 ! of the buffer region in order to avoid systematic errors
	 ! in QM forces.

	 call initialise(outer_layer, 1,0,0,0)
	 call append(outer_layer, find(hybrid_mark == HYBRID_BUFFER_OUTER_LAYER_MARK))

	 do i=1,outer_layer%N
	    ! Shift atom by randomly distributed vector of magnitude 0.05 A
	    j = find_in_array(cluster_index, outer_layer%int(1,i))
	    cluster%pos(:,j) = cluster%pos(:,j) + 0.05_dp*random_unit_vector()
	 end do

	 call finalise(outer_layer)
      end if

    end if

    if (value(mainlog%verbosity_stack) >= VERBOSE .or. print_clusters) then
#ifdef _MPI
       write (clusterfilename, '(a,i3.3,a)') 'clusters.',mpi_rank,'.xyz'
#else
       clusterfilename = 'clusters.xyz'
#endif _MPI
       call initialise(clusterfile, clusterfilename, append=.true., action=OUTPUT)
       call inoutput_mpi_all_inoutput(clusterfile, .true.)
       call print_xyz(cluster, clusterfile, all_properties=.true.)
       call finalise(clusterfile)
    end if

    call finalise(currentlist)
    call finalise(nextlist)
    call finalise(cluster_list)
    call finalise(outer_layer)

  end function create_cluster_hybrid_mark

  
  subroutine create_hybrid_weights_args_str(at, args_str)
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in) :: args_str

    type(Dictionary) :: params
    logical :: min_images_only, mark_buffer_outer_layer, nneighb_only, hysteretic_buffer
    real(dp) :: hysteretic_buffer_inner_radius, hysteretic_buffer_outer_radius
    integer :: buffer_hops, transition_hops
    character(FIELD_LENGTH) :: weight_interpolation
    
    call initialise(params)
    call param_register(params, 'min_images_only', 'F', min_images_only)
    call param_register(params, 'mark_buffer_outer_layer', 'T', mark_buffer_outer_layer)
    call param_register(params, 'transition_hops', '0', transition_hops)
    call param_register(params, 'buffer_hops', '3', buffer_hops)
    call param_register(params, 'nneighb_only', 'T', nneighb_only)
    call param_register(params, 'weight_interpolation', 'hop_ramp', weight_interpolation)
    call param_register(params, 'hysteretic_buffer', 'F', hysteretic_buffer)
    call param_register(params, 'hysteretic_buffer_inner_radius', '5.0', hysteretic_buffer_inner_radius)
    call param_register(params, 'hysteretic_buffer_outer_radius', '7.0', hysteretic_buffer_outer_radius)
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.) ) &
      call system_abort("create_hybrid_weights_args_str failed to parse args_str='"//trim(args_str)//"'")
    call finalise(params)

    call create_hybrid_weights(at, trans_width=transition_hops, buffer_width=buffer_hops, &
         weight_interpolation=weight_interpolation, nneighb_only=nneighb_only, &
         min_images_only = min_images_only, mark_buffer_outer_layer=mark_buffer_outer_layer, &
         hysteretic_buffer=hysteretic_buffer, hysteretic_buffer_inner_radius=hysteretic_buffer_inner_radius, &
         hysteretic_buffer_outer_radius=hysteretic_buffer_outer_radius)

    
  end subroutine create_hybrid_weights_args_str


  !% Given an atoms structure with a 'hybrid_mark' property, this routine
  !% creates a 'weight_region1' property, whose values are between 0 and
  !% 1. Atoms marked with HYBRID_ACTIVE_MARK in 'hybrid_mark' get weight
  !% 1, the neighbourhopping is done trans_width times, during which
  !% the weight linearly decreases to zero (either with hop count if 'weight_interpolation=hop_ramp'
  !% or with distance from the centre of mass of the embed region if 
  !% 'weight_interpolation=distance_ramp') decreases and atoms are marked with
  !% HYBRID_TRANS_MARK. Further hopping is done 'buffer_width' times
  !% and atoms are marked with HYBRID_BUFFER_MARK and given weight
  !% zero.
  !%
  !% Optionally we return a Table 'clusterlist' which includes the indices
  !% and periodic shifts relative to the first marked atom of all atoms in
  !% the cluster. This is suitable for passing to create_cluster.
  subroutine create_hybrid_weights_args(at, trans_width, buffer_width, weight_interpolation, &
       nneighb_only, min_images_only, mark_buffer_outer_layer, hysteretic_buffer, &
       hysteretic_buffer_inner_radius, hysteretic_buffer_outer_radius)
    type(Atoms), intent(inout) :: at
    integer, intent(in) :: trans_width, buffer_width
    character(len=*), optional, intent(in) :: weight_interpolation
    logical, optional, intent(in) :: nneighb_only, min_images_only, mark_buffer_outer_layer, hysteretic_buffer
    real(dp), optional, intent(in) :: hysteretic_buffer_inner_radius, hysteretic_buffer_outer_radius
    !logical :: terminate

    integer, pointer :: hybrid_mark(:)
    real(dp), pointer :: weight_region1(:)
    integer :: n_region1, n_trans, n_region2 !, n_term
    integer :: i, j, jj, first_active, shift(3)
    logical :: dummy
    logical :: do_nneighb_only, do_min_images_only, do_mark_buffer_outer_layer, distance_ramp, hop_ramp, do_hysteretic_buffer
    character(FIELD_LENGTH) :: do_weight_interpolation
    type(Table) :: activelist, currentlist, nextlist, distances, oldbuffer, bufferlist, buffer_inner, buffer_outer, embedlist
    real(dp) :: core_CoM(3), core_mass, mass, do_hysteretic_buffer_inner_radius, do_hysteretic_buffer_outer_radius

    do_nneighb_only = optional_default(.false., nneighb_only)
    do_min_images_only = optional_default(.true., min_images_only)
    do_mark_buffer_outer_layer = optional_default(.true., mark_buffer_outer_layer)
    do_weight_interpolation = optional_default('hop_ramp', weight_interpolation)
    do_hysteretic_buffer = optional_default(.false., hysteretic_buffer)
    do_hysteretic_buffer_inner_radius = optional_default(5.0_dp, hysteretic_buffer_inner_radius)
    do_hysteretic_buffer_outer_radius = optional_default(7.0_dp, hysteretic_buffer_outer_radius)

    hop_ramp = .false.
    distance_ramp = .false.
    if (trim(do_weight_interpolation) == 'hop_ramp') then
       hop_ramp = .true.
    else if (trim(do_weight_interpolation) == 'distance_ramp') then
       distance_ramp = .true.
    else
       call system_abort('create_hybrid_weights_args: unknown weight_interpolation value: '//trim(do_weight_interpolation))
    end if

    call print('create_hybrid_weights: trans_width='//trans_width//' buffer_width='//buffer_width//' weight_interpolation='//do_weight_interpolation, VERBOSE)
    call print('  nneighb_only='//do_nneighb_only//' min_images_only='//do_min_images_only//' mark_buffer_outer_layer='//do_mark_buffer_outer_layer, VERBOSE)
    call print('  hysteretic_buffer='//do_hysteretic_buffer//' hysteretic_buffer_inner_radius='//do_hysteretic_buffer_inner_radius, VERBOSE)
    call print('  hysteretic_buffer_outer_radius='//do_hysteretic_buffer_outer_radius, VERBOSE)

    ! check to see if atoms has a 'weight_region1' property already, if so, check that it is compatible, if not present, add it
    if(assign_pointer(at, 'weight_region1', weight_region1)) then
       weight_region1 = 0.0_dp
    else
       call add_property(at, 'weight_region1', 0.0_dp)
       dummy = assign_pointer(at, 'weight_region1', weight_region1)
    end if

    ! check for a compatible hybrid_mark property. it must be present
    if(.not.assign_pointer(at, 'hybrid_mark', hybrid_mark)) &
         call system_abort('create_hybrid_weights: atoms structure has no "hybrid_mark" property')

    ! Add first marked atom to activelist. shifts will be relative to this atom
    first_active = find_in_array(hybrid_mark, HYBRID_ACTIVE_MARK)
    call append(activelist, (/first_active,0,0,0/))
    call append(currentlist, activelist)
    weight_region1(first_active) = 1.0_dp

    if (distance_ramp) then
       if (has_property(at, 'mass')) then
          core_mass = at%mass(first_active)
       else
          core_mass = ElementMass(at%Z(first_active))
       end if
       core_CoM =core_mass*at%pos(:,first_active) ! we're taking this as reference atom so no shift needed here
    end if

    n_region1 = count(hybrid_mark == HYBRID_ACTIVE_MARK)

    ! Add other active atoms using bond hopping from the first atom
    ! in the cluster to find the other active atoms and hence to determine the correct 
    ! periodic shifts
    do while (activelist%N < n_region1)
       call BFS_step(at, currentlist, nextlist, nneighb_only = do_nneighb_only, min_images_only = do_min_images_only)
       do j=1,nextlist%N
          jj = nextlist%int(1,j)
          shift = nextlist%int(2:4,j)
          if (hybrid_mark(jj) == HYBRID_ACTIVE_MARK) then
             call append(activelist, nextlist%int(:,j))
             weight_region1(jj) = 1.0_dp

             if (distance_ramp) then
                if (has_property(at, 'mass')) then
                   mass = at%mass(jj)
                else
                   mass = ElementMass(at%Z(jj))
                end if
                core_CoM = core_CoM + mass*(at%pos(:,jj) + (at%lattice .mult. (shift)))
                core_mass = core_mass + mass
             end if

          end if
       end do
       call append(currentlist, nextlist)
    end do

    if (distance_ramp) core_CoM = core_CoM/core_mass

    call wipe(currentlist)
    call append(currentlist, activelist)

    ! create transition region
    if (distance_ramp) call initialise(distances, 1, 1, 0, 0)
    n_trans = 0
    do i = 0,trans_width-2
       call BFS_step(at, currentlist, nextlist, nneighb_only = do_nneighb_only, min_images_only = do_min_images_only)
       call wipe(currentlist)
       do j = 1,nextlist%N
          jj = nextlist%int(1,j)
          if(hybrid_mark(jj) == HYBRID_NO_MARK) then
             call append(currentlist, nextlist%int(:,j))
             if (hop_ramp) weight_region1(jj) = 1.0_dp - real(i+1)/real(trans_width) ! linear transition
             if (distance_ramp) &        ! Save distance, weight will be calculated later
                  call append(distances, jj, distance_min_image(at, core_CoM, jj))
             hybrid_mark(jj) = HYBRID_TRANS_MARK
             n_trans = n_trans+1
          end if
       end do
    end do


    if (distance_ramp) then
       ! Normalise distances from core_CoM so that maximum distance is equal to 1
       distances%real(1,1:distances%N) = distances%real(1,1:distances%N)/maxval(real_part(distances,1))
       ! Fill in the weights
       do i=1,distances%N
          weight_region1(jj) = 1.0_dp - distances%real(1,i)
       end do
       call finalise(distances)
    end if

    call append(embedlist, currentlist)

    ! create region2 (buffer region)
    if (.not. do_hysteretic_buffer) then
       n_region2 = 0
       do i = 0,buffer_width-1
          call BFS_step(at, currentlist, nextlist, nneighb_only = do_nneighb_only, min_images_only = do_min_images_only)
          call wipe(currentlist)
          do j = 1,nextlist%N
             jj = nextlist%int(1,j)
             if(hybrid_mark(jj) == HYBRID_NO_MARK) then
                call append(currentlist, nextlist%int(:,j))
                weight_region1(jj) = 0.0_dp
                if (i==buffer_width-1 .and. do_mark_buffer_outer_layer) then
                   hybrid_mark(jj) = HYBRID_BUFFER_OUTER_LAYER_MARK
                else
                   hybrid_mark(jj) = HYBRID_BUFFER_MARK
                end if
                n_region2 = n_region2+1
             end if
          end do
       end do
    else

       call wipe(oldbuffer)
       call append(oldbuffer, find(hybrid_mark /= HYBRID_NO_MARK))

       call wipe(currentlist)

       ! Add first marked atom to embedlist. shifts will be relative to this atom
       first_active = find_in_array(hybrid_mark, HYBRID_ACTIVE_MARK)
       call append(bufferlist, (/first_active,0,0,0/))
       call append(currentlist, bufferlist)

       n_region2 = count(hybrid_mark /= HYBRID_NO_MARK)

       ! Find old embed + buffer atoms using bond hopping from the first atom
       ! in the cluster to find the other buffer atoms and hence to determine the correct 
       ! periodic shifts
       do while (bufferlist%N < n_region2)
          call BFS_step(at, currentlist, nextlist, nneighb_only = do_nneighb_only, min_images_only = do_min_images_only)
          do j=1,nextlist%N
             jj = nextlist%int(1,j)
             shift = nextlist%int(2:4,j)
             if (hybrid_mark(jj) == HYBRID_ACTIVE_MARK) call append(bufferlist, nextlist%int(:,j))
          end do
          call append(currentlist, nextlist)
       end do
    
       ! Remove marks on all buffer atoms
       do i=1,oldbuffer%N
          if (hybrid_mark(oldbuffer%int(1,i)) == HYBRID_BUFFER_MARK .or. &
              hybrid_mark(oldbuffer%int(1,i)) == HYBRID_BUFFER_OUTER_LAYER_MARK) &
              hybrid_mark(oldbuffer%int(1,i)) = HYBRID_NO_MARK
       end do

       ! Construct inner and outer buffer lists
       call construct_buffer(at, embedlist, do_hysteretic_buffer_inner_radius, buffer_inner, has_property(at, 'avgpos'))
       call construct_buffer(at, embedlist, do_hysteretic_buffer_outer_radius, buffer_outer, has_property(at, 'avgpos'))

       call print('bufferlist=')
       call print(bufferlist)

       call print('buffer_inner=')
       call print(buffer_inner)

       call print('buffer_outer=')
       call print(buffer_outer)

       call select_hysteretic_quantum_region(at, buffer_inner, buffer_outer, bufferlist)

       ! Mark new buffer region, leaving core QM region alone
       do i=1,bufferlist%N
          if (hybrid_mark(bufferlist%int(1,i)) == HYBRID_NO_MARK) & 
               hybrid_mark(bufferlist%int(1,i)) = HYBRID_BUFFER_MARK

          ! Marking the outer layer with  HYBRID_BUFFER_OUTER_LAYER_MARK is
          ! dealt with in create_cluster_hybrid_mark.

       end do

       call finalise(buffer_inner)
       call finalise(buffer_outer)
       call finalise(bufferlist)
       call finalise(oldbuffer)
       call finalise(embedlist)
       
    end if

       ! this is commented out for now, terminations are controlled in create_cluster only
       ! create terminations
       !n_term = 0
       !if(terminate) then
       !   call BFS_step(at, currentlist, nextlist, nneighb_only = .true., min_images_only = .true.)
       !   do j = 1,nextlist%N
       !      jj = nextlist%int(1,j)
       !      if(hybrid_mark(jj) == HYBRID_NO_MARK) then
       !         weight_region1(jj) = 0.0_dp
       !         hybrid_mark(jj) = HYBRID_TERM_MARK
       !         n_term = n_term+1
       !      end if
       !   end do
       !end do

    call print('create_hybrid_weights: '//n_region1//' region 1, '//n_trans//' transition, '//n_region2//&
         ' region 2, '//count(hybrid_mark /= HYBRID_NO_MARK)//' in total', VERBOSE)

    call finalise(activelist)
    call finalise(currentlist)
    call finalise(nextlist)
    call finalise(distances)
    
  end subroutine create_hybrid_weights_args


  !% Given an Atoms structure with an active region marked in the 'hybrid_mark'
  !% property using 'HYBRID_ACTIVE_MARK', grow the embed region by 'fit_hops'
  !% bond hops to form a fit region. Returns the embedlist and fitlist with correct
  !% periodic shifts.
  subroutine create_embed_and_fit_lists(at, fit_hops, embedlist, fitlist, nneighb_only, min_images_only)

    type(Atoms), intent(inout) :: at
    integer :: fit_hops
    type(Table), intent(out) :: embedlist, fitlist
    logical, intent(in), optional :: nneighb_only, min_images_only

    integer, pointer :: hybrid_mark(:)
    integer :: n_region1, n_region2 !, n_term
    integer :: i, j, jj, first_active, shift(3)
    logical :: do_nneighb_only, do_min_images_only
    type(Table) :: currentlist, nextlist, tmpfitlist
    
    do_nneighb_only = optional_default(.false., nneighb_only)
    do_min_images_only = optional_default(.true., min_images_only)

    ! check for a compatible hybrid_mark property. it must be present
    if(.not.assign_pointer(at, 'hybrid_mark', hybrid_mark)) &
       call system_abort('create_fit_region: atoms structure has no "hybrid_mark" property')

    call wipe(embedlist)
    call wipe(fitlist)

    ! Add first marked atom to embedlist. shifts will be relative to this atom
    first_active = find_in_array(hybrid_mark, HYBRID_ACTIVE_MARK)
    call append(embedlist, (/first_active,0,0,0/))
    call append(currentlist, embedlist)

    n_region1 = count(hybrid_mark == HYBRID_ACTIVE_MARK)

    ! Add other active atoms using bond hopping from the first atom
    ! in the cluster to find the other active atoms and hence to determine the correct 
    ! periodic shifts
    do while (embedlist%N < n_region1)
       call BFS_step(at, currentlist, nextlist, nneighb_only = do_nneighb_only, min_images_only = do_min_images_only)
       do j=1,nextlist%N
          jj = nextlist%int(1,j)
          shift = nextlist%int(2:4,j)
          if (hybrid_mark(jj) == HYBRID_ACTIVE_MARK) call append(embedlist, nextlist%int(:,j))
       end do
       call append(currentlist, nextlist)
    end do
    
    call wipe(currentlist)
    call append(currentlist, embedlist)

!    call append(fitlist, embedlist)

    ! create region2 (fit region)
    n_region2 = 0
    do i = 0,fit_hops-1
       call BFS_step(at, currentlist, nextlist, nneighb_only = do_nneighb_only, min_images_only = do_min_images_only)
       do j = 1,nextlist%N
          jj = nextlist%int(1,j)
          if(hybrid_mark(jj) /= HYBRID_ACTIVE_MARK) then
             call append(tmpfitlist, nextlist%int(:,j))
             n_region2 = n_region2+1
          end if
       end do
       call append(currentlist, nextlist)
    end do

    call print('create_embed_and_fit_lists: '//n_region1//' embed, '//n_region2//' fit', VERBOSE)

    ! Sort order to we are stable to changes in neighbour ordering introduced
    ! by calc_connect. 
    call sort(embedlist)
    call sort(tmpfitlist)

    ! fitlist consists of sorted embedlist followed by sorted list of remainder of fit atoms
    call append(fitlist, embedlist)
    call append(fitlist, tmpfitlist)

    call finalise(currentlist)
    call finalise(nextlist)
    call finalise(tmpfitlist)

  end subroutine create_embed_and_fit_lists

  !% Return the atoms in a hysteretic quantum region:
  !% To become part of the quantum region, atoms must drift within
  !% 'inner' of any core atom. To leave the quantum region, an atom
  !% must drift further than 'outer; from any core atom.
  !% Optionally use time averaged positions
  !

  subroutine select_hysteretic_quantum_region(at,inner,outer,list,verbosity)

    type(Atoms),       intent(in)    :: at
    type(Table),       intent(inout) :: inner, outer
    type(Table),       intent(inout) :: list
    integer, optional, intent(in)    :: verbosity

    integer                          :: n, my_verbosity, atom(list%intsize)

    my_verbosity = optional_default(NORMAL,verbosity)

    if (my_verbosity == ANAL) then
       call print('In Select_Hysteretic_Quantum_Region:')
       call print('List currently contains '//list%N//' atoms')
    end if

    if ((list%intsize /= inner%intsize) .or. (list%intsize /= outer%intsize)) &
       call system_abort('select_hysteretic_quantum_region: inner, outer and list must have the same intsize')

    ! Speed up searching
    call sort(outer)
    call sort(list)

    !Check for atoms in 'list' and not in 'outer'
    do n = list%N, 1, -1
       atom = list%int(:,n) !the nth atom in the list
       if (search(outer,atom)==0) then
          call delete(list,n)
          if (my_verbosity > NORMAL) call print('Removed atom ('//atom//') from quantum list')
       end if
    end do

    call sort(list)

    !Check for new atoms in 'inner' cluster and add them to list
    do n = 1, inner%N
       atom = inner%int(:,n)
       if (search(list,atom)==0) then
          call append(list,atom)
          call sort(list)
          if (my_verbosity > NORMAL) call print('Added atom ('//atom//') to quantum list')
       end if
    end do

    if (my_verbosity >= NORMAL) call print(list%N//' atoms selected for quantum treatment')

  end subroutine select_hysteretic_quantum_region

  !
  !% Given an atoms object, and a list of core atoms in the first integer of
  !% the 'core' table, fill the 'buffer' table with all atoms within 'radius'
  !% of any core atom (which can be reached by connectivity hopping). 
  !% Optionally use the time averaged positions. 
  !
  subroutine construct_buffer(at,core,radius,buffer,use_avgpos,verbosity)

    type(Atoms),       intent(in)  :: at
    type(Table),       intent(in)  :: core
    real(dp),          intent(in)  :: radius
    type(Table),       intent(out) :: buffer
    logical, optional, intent(in)  :: use_avgpos
    integer, optional, intent(in)  :: verbosity

    logical                             :: do_use_avgpos, more_atoms, add_atom
    integer                             :: buffer_atoms_guess, i(4), j(4), n, nn, my_verbosity, lookup(3)
    type(Table)                         :: cluster, extra_atoms
    integer,  allocatable, dimension(:) :: append_row_int
    real(dp), allocatable, dimension(:) :: append_row_real
    
    !Check the input arguments
    if (.not.at%connect%initialised) call system_abort('Construct_Buffer: Atomic connectivity data required')
    if (core%intsize < 1) call system_abort('Construct_Buffer: Core table must have at least one integer')
    if (radius < 0.0_dp) call system_abort('Construct_Buffer: Buffer radius must be positive')
    do_use_avgpos = .false.
    if (present(use_avgpos)) do_use_avgpos = use_avgpos
    my_verbosity = NORMAL
    if (present(verbosity)) my_verbosity = verbosity
    if (do_use_avgpos .and. .not. get_value(at%properties, 'avgpos', lookup)) &
         call system_abort('Construct_Buffer: Average positions are not present in atoms structure')

    if (my_verbosity == ANAL) call print('In Construct_Buffer:')

    !Guess the number of buffer atoms we'll add from the average number density and radius
    buffer_atoms_guess = int(10 * at%N * radius*radius*radius / Cell_Volume(at))

    !Set up the output table and append rows
    call allocate(buffer,core%intsize,core%realsize,0,0,buffer_atoms_guess)
    allocate(append_row_int(core%intsize),append_row_real(core%realsize))
    append_row_int = 0
    append_row_real = 0.0_dp

    !Add the core atoms
    do n = 1, core%N
       append_row_int = core%int(:,n)
       if (core%realsize > 0) then
          call append(buffer,append_row_int,append_row_real)
       else
          call append(buffer,append_row_int)
       end if       
    end do

    !Loop through the neighbours and add them to the cluster if they are close enough.
    !If we add any atoms this way, trigger an additional neighbour search       

    more_atoms = .true.
    
    do while(more_atoms)
       
       more_atoms = .false.
       
       call wipe(cluster)
       call append(cluster,core)
       call append(cluster,buffer)

       !Find the neighbours of all the current cluster atoms 
       if (my_verbosity == ANAL) call print('Searching for neighbours...')
       call bfs_step(at,cluster,extra_atoms,nneighb_only=.false.,min_images_only=.true.)

       if (my_verbosity == ANAL) then
          write(line,'(a,i0,a)')'Found ',extra_atoms%N,' neighbours'
          call print(line)
       end if
          
       !Now go through the found atoms and see if they are within 'radius' of ANY core atom
       do n = 1, extra_atoms%N
          
          add_atom = .false.

          i = extra_atoms%int(:,n)
          
          do nn = 1, core%N
             
             j = core%int(:,nn)
             
             if (do_use_avgpos) then
                if (norm(diff_min_image(at,at%avgpos(:,i(1)),at%avgpos(:,j(1))) + &
                          (at%lattice .mult. (j(2:4)-i(2:4)))) < radius) then
                   append_row_int = i
                   add_atom = .true.
                end if
             else
                if (distance(at,i(1),j(1),j(2:4)-i(2:4)) < radius) then
                   append_row_int = i
                   add_atom = .true.
                end if
             end if

             if (add_atom) then
                if (core%realsize > 0) then
                   call append(buffer,append_row_int,append_row_real)
                else
                   call append(buffer,append_row_int)
                end if
                more_atoms = .true. !Trigger another neighbour search
                if (my_verbosity == ANAL) then                 
                   call print('Adding atom '//i)
                end if
                exit !Don't check other core atoms
             end if

          end do
          
       end do

    end do

    call finalise(cluster)
    call finalise(extra_atoms)

    if (my_verbosity == ANAL) call print('Leaving Construct_Buffer')

  end subroutine construct_buffer

  
  subroutine update_active(this, nneightol, avgpos, reset)
    type(Atoms) :: this
    real(dp), optional :: nneightol
    logical, optional :: avgpos, reset

    type(Atoms), save :: nn_atoms
    integer, pointer, dimension(:) :: nn, old_nn, active
    integer  :: i
    real(dp) :: use_nn_tol
    logical  :: use_avgpos, do_reset

    use_nn_tol = optional_default(this%nneightol, nneightol)
    use_avgpos = optional_default(.true., avgpos)
    do_reset   = optional_default(.false., reset)
    
    ! First time copy entire atoms structure into nn_atoms and
    ! add properties for nn, old_nn and active
    if (reset .or. .not. nn_atoms%initialised) then
       nn_atoms = this
       call add_property(this, 'nn', 0)
       call add_property(this, 'old_nn', 0)
       call add_property(this, 'active', 0)
       
       call atoms_set_cutoff_factor(nn_atoms, use_nn_tol)
    end if

    if (this%N /= nn_atoms%N) &
         call system_abort('update_actives: Number mismatch between this%N ('//this%N// &
         ') and nn_atoms%N ('//nn_atoms%N//')')

    if (.not. assign_pointer(this, 'nn', nn)) &
         call system_abort('update_actives: Atoms is missing "nn" property')

    if (.not. assign_pointer(this, 'old_nn', old_nn)) &
         call system_abort('update_actives: Atoms is missing "old_nn" property')

    if (.not. assign_pointer(this, 'active', active)) &
         call system_abort('update_actives: Atoms is missing "active" property')

    call print('update_actives: recalculating nearest neighbour table', VERBOSE)
    if (use_avgpos .and. associated(this%avgpos)) then
       call print('update_actives: using time averaged atomic positions', VERBOSE)
       nn_atoms%pos = this%avgpos
    else
       call print('update_actives: using instantaneous atomic positions', VERBOSE)
       nn_atoms%pos = this%pos
    end if
    call calc_connect(nn_atoms)
    
    nn = 0
    do i = 1,nn_atoms%N
       nn(i) = atoms_n_neighbours(nn_atoms, i)
    end do

    if (all(old_nn == 0)) old_nn = nn ! Special case for first time

    ! Decrement the active counts
    where (active /= 0) active = active -1

    ! Find newly active atoms
    where (nn /= old_nn) active = 1

    if (count(active == 1) /= 0) then
       call print('update_actives: '//count(active == 1)//' atoms have become active.')
    end if

    old_nn = nn

  end subroutine update_active


  !
  !% Given an atoms structure and a list of quantum atoms, find X-H
  !% bonds which have been cut and include the other atom of 
  !% the pair in the quantum list.
  !
  subroutine add_cut_hydrogens(this,qmlist,verbosity)

    type(Atoms),       intent(in)    :: this
    type(Table),       intent(inout) :: qmlist
    integer, optional, intent(in)    :: verbosity
    
    type(Table)                :: neighbours, bonds, centre
    logical                    :: more_atoms
    integer                    :: i, j, n, nn, added
    
    more_atoms = .true.
    added = 0
    call allocate(centre,4,0,0,0,1)

    !Repeat while the search trigger is true
    do while(more_atoms)

       more_atoms = .false.

       !Find nearest neighbours of the cluster
       call bfs_step(this,qmlist,neighbours,nneighb_only=.true.,min_images_only=.true.)

       !Loop over neighbours
       do n = 1, neighbours%N

          i = neighbours%int(1,n)
          
          call wipe(centre)
          call append(centre,(/i,0,0,0/))

          ! Find atoms bonded to this neighbour
          call bfs_step(this,centre,bonds,nneighb_only=.true.,min_images_only=.true.)
          
          !Loop over these bonds
          do nn = 1, bonds%N
             
             j = bonds%int(1,nn)

             !Try to find j in the qmlist
             if (find_in_array(int_part(qmlist,1),j)/=0) then
                !We have a cut bond. If i or j are hydrogen then add i to the 
                !quantum list and trigger another search after this one
                if (this%Z(i)==1 .or. this%Z(j)==1) then
                   call append(qmlist,(/i,0,0,0/))
                   if (present(verbosity)) then
                      if (verbosity >= NORMAL) call print('Add_Cut_Hydrogens: Added atom '//i//', neighbour of atom '//j)
                   end if
                   more_atoms = .true.
                   added = added + 1
                   exit !Don't add the same atom more than once
                end if
             end if
             
          end do

       end do

    end do

    !Report findings
    if (present(verbosity)) then
       if (verbosity >= NORMAL) then
          write(line,'(a,i0,a)')'Add_Cut_Hydrogens: Added ',added,' atoms to quantum region'
          call print(line)
       end if
    end if

    call finalise(centre)
    call finalise(bonds)
    call finalise(neighbours)

  end subroutine add_cut_hydrogens

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X GROUP CONVERSION
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !
  !% Convert a constrained atom's group to a normal group prior to quantum treatment
  !% The constraints are not deleted, just ignored. The number of degrees of freedom
  !% of the system is also updated.
  !
  subroutine constrained_to_quantum(this,i)

    type(DynamicalSystem), intent(inout) :: this
    integer,               intent(in)    :: i
    integer                              :: g!, n, j, a1, a2

    g = this%group_lookup(i)

    if (this%group(g)%type == TYPE_CONSTRAINED) then
       write(line,'(2(a,i0),a)')' INFO: Converting group ',g,' from Constrained to Quantum'
       call print(line)
       
       !Change the group type
       call set_type(this%group(g),TYPE_ATOM)
       
       !Update the number of degrees of freedom (add 1 per removed constraint)
       this%Ndof = this%Ndof + Group_N_Objects(this%group(g))

       !Thermalize the released constraints
!       do n = 1, group_n_objects(this%group(g))
!
!          j = group_nth_object(this%group(g),n)
!
!          !If not a two body constraint then go to the next one.
!          if (this%constraint(j)%N /= 2) cycle
!
!          !Get the indices of the atoms in this constraint
!          a1 = this%constraint(j)%atom(1)
!          a2 = this%constraint(j)%atom(2)
!
!          call thermalize_bond(this,a1,a2,this%sim_temp)
!
!       end do

!       call thermalize_group(this,g,this%sim_temp)

    end if

  end subroutine constrained_to_quantum

  !
  !% When a constraint is released as a molecule drifts into the quantum region the
  !% new, unconstrained degree of freedom has zero temperature. This routine adds
  !% velocities in the directions between atoms which make up two body constraints
  !% at a specified temperature
  !
  subroutine thermalize_bond(this,i,j,T)

    type(DynamicalSystem), intent(inout) :: this         !% The dynamical system
    integer,               intent(in)    :: i,j          !% The two atoms
    real(dp),              intent(in)    :: T            !% Temperature
    real(dp)                             :: meff         !  effective mass
    real(dp), dimension(3)               :: u_hat        !  unit vector between the atoms
    real(dp)                             :: w, wi, wj, & ! magnitude of velocity component in bond direction
                                            mi, mj       ! masses of the atoms
    real(dp)                             :: E0,E1        ! kinetic energy in bond before and after
    real(dp), save                       :: Etot = 0.0_dp ! cumulative energy added to the system
    real(dp), save                       :: Eres = 0.0_dp ! cumulative energy already in bonds (QtoC error)
    real(dp)                             :: r, v         ! relative distance and velocity

    call print('Thermalize_Bond: Thermalizing bond '//i//'--'//j)

    ! Find unit vector along bond and current relative velocity
    call distance_relative_velocity(this%atoms,i,j,r,v)

    ! Find unit vector between atoms
    u_hat = diff_min_image(this%atoms,i,j)
    u_hat = u_hat / norm(u_hat)

    ! Calculate effective mass
    mi = ElementMass(this%atoms%Z(i))
    mj = ElementMass(this%atoms%Z(j))
    meff = mi * mj / (mi + mj)

    ! Calculate energy currently stored in this bond (should be zero)
    E0 = 0.5_dp * meff * v * v
    Eres = Eres + E0

    call print('Thermalize_Bond: Bond currently contains '//round(E0,8)//' eV')
    call print('Thermalize_Bond: Residual energy encountered so far: '//round(Eres,8)//' eV')

    ! Draw a velocity for this bond
    w = gaussian_velocity_component(meff,T)
    w = w - v ! Adjust for any velocity which the bond already had

    ! Distribute this between the two atoms, conserving momentum
    wi = (meff / mi) * w
    wj = -(meff / mj) * w
    
    this%atoms%velo(:,i) = this%atoms%velo(:,i) + wi * u_hat
    this%atoms%velo(:,j) = this%atoms%velo(:,j) + wj * u_hat

    call distance_relative_velocity(this%atoms,i,j,r,v)
    E1 = 0.5_dp * meff * v * v
        
    Etot = Etot + E1 - E0

    call print('Thermalize_Bond: Added '//round(E1-E0,8)//' eV to bond')
    call print('Thermalize_Bond: Added '//round(Etot,8)//' eV to the system so far')

  end subroutine thermalize_bond

  !
  !% Convert a quantum atom's group to a constrained group. The group is searched for
  !% two body constraints, which are assumed to be bond length or time dependent bond length
  !% constraints. The required final length is extracted, and a new time dependent bond length
  !% constraint is imposed, taking the current length and relative velocity along the bond
  !% to the required values, in the time given by smoothing_time.
  !%
  !% Any other constraints are left untouched.
  !
  subroutine quantum_to_constrained(this,i,smoothing_time)

    type(DynamicalSystem), intent(inout) :: this
    integer,               intent(in)    :: i
    real(dp),              intent(in)    :: smoothing_time
    integer, save                        :: CUBIC_FUNC
    logical, save                        :: first_call = .true.
    integer                              :: g, j, n, a1, a2, datalength
    real(dp)                             :: x0,y0,y0p,x1,y1,y1p, coeffs(4), t, dE, m1, m2, meff
    real(dp), save                       :: Etot = 0.0_dp

    !Register the constraint function if this is the first call
    if (first_call) then
       CUBIC_FUNC = Register_Constraint(CUBIC_BOND)
       first_call = .false.
    end if

    g = this%group_lookup(i)

    if (this%group(g)%type == TYPE_ATOM) then
       write(line,'(2(a,i0),a)')' INFO: Converting group ',g,' from Quantum to Constrained'
       call print(line)
       
       !Change the group type
       call set_type(this%group(g),TYPE_CONSTRAINED)
       
       !Update the number of degrees of freedom (subtract 1 per added constraint)
       this%Ndof = this%Ndof - Group_N_Objects(this%group(g))

       !Loop over constraints
       do n = 1, Group_N_Objects(this%group(g))

          j = Group_Nth_Object(this%group(g),n)

          !If not a two body constraint then go to the next one.
          if (this%constraint(j)%N /= 2) cycle

          !Get the indices of the atoms in this constraint
          a1 = this%constraint(j)%atom(1)
          a2 = this%constraint(j)%atom(2)

          !Detect the type of constraint (simple bondlength or time dependent bond length) 
          !from the size of the data array (1 or 6)
          datalength = size(this%constraint(j)%data)
          select case(datalength)
          case(1) 
             !Simple bond length: the required final length is the only entry
             y1 = this%constraint(j)%data(1)
          case(6)
             !Time dependent bond length: the final length is the polynomial evaluated at the end point
             coeffs = this%constraint(j)%data(1:4)
             t = this%constraint(j)%data(6)
             y1 = ((coeffs(1)*t + coeffs(2))*t + coeffs(3))*t + coeffs(4)
          case default
             !This isn't a known bond length constraint. Set y1 to a negative value.
             y1 = -1.0_dp
          end select
          
          !If an existing bond length constraint has been found then set up a time dependent constraint
          if (y1 > 0.0_dp) then

             call print('Quantum_to_Constrained: Reconstraining bond '//a1//'--'//a2)
             
             x0 = this%t
             call distance_relative_velocity(this%atoms,a1,a2,y0,y0p)
             x1 = this%t + smoothing_time
             y1p = 0.0_dp

             !Find the coefficients of the cubic which fits these boundary conditions
             call fit_cubic(x0,y0,y0p,x1,y1,y1p,coeffs)

             !Change the constraint
             call ds_amend_constraint(this,j,CUBIC_FUNC,(/coeffs,x0,x1/))

             !Work out the amount of energy which will be removed
             m1 = ElementMass(this%atoms%Z(a1))
             m2 = ElementMass(this%atoms%Z(a2))
             meff = m1 * m2 / (m1 + m2)

             dE = 0.5_dp * meff * y0p * y0p
             Etot = Etot + dE

             call print('Quantum_to_Constrained: Gradually removing '//round(dE,8)//' eV from bond')
             call print('Quantum_to_Constrained: Will have removed approximately '//round(Etot,8)//' eV in total')

          end if

       end do

    end if

  end subroutine quantum_to_constrained
 

  subroutine thermalize_group(this,g,T)

    type(DynamicalSystem), intent(inout) :: this
    integer,               intent(in)    :: g
    real(dp),              intent(in)    :: T

    integer  :: i, j, k
    real(dp) :: mi, mj, mk, mij, mik, mjk, vij, vik, vjk, rij, rik, rjk, gij, gik, gjk
    real(dp) :: ri(3), rj(3), rk(3), Rcom(3)
    real(dp) :: A(9,9), b(9), A_inv(9,9), x(9)

    real(dp) :: dEij, dEik, dEjk
    real(dp), save :: Etot = 0.0_dp

    ! Get atoms
    i = this%group(g)%atom(1)
    j = this%group(g)%atom(2)
    k = this%group(g)%atom(3)

    ! Get shifted positions
    ri = 0.0_dp
    rj = diff_min_image(this%atoms,i,j)
    rk = diff_min_image(this%atoms,i,k)  

    ! Get separations and current bond velocities
    call distance_relative_velocity(this%atoms,i,j,rij,vij)
    call distance_relative_velocity(this%atoms,i,k,rik,vik)
    call distance_relative_velocity(this%atoms,j,k,rjk,vjk)

    ! Calculate masses
    mi = ElementMass(this%atoms%Z(i))
    mj = ElementMass(this%atoms%Z(j))
    mk = ElementMass(this%atoms%Z(k))   
    mij = mi * mj / (mi + mj)
    mik = mi * mk / (mi + mk)
    mjk = mj * mk / (mj + mk)

    ! Calculate centre of mass
    Rcom = (mi*ri + mj*rj + mk*rk) / (mi+mj+mk)

    ! Select a gaussian distributed velocity adjustment for each bond
    gij = gaussian_velocity_component(mij,T)
    gik = gaussian_velocity_component(mik,T)
    gjk = gaussian_velocity_component(mjk,T)
    vij = gij - vij
    vik = gik - vik
    vjk = gjk - vjk

    ! Build matrix...
    A = 0.0_dp
    b = 0.0_dp

    ! Momentum conservation
    A(1,1) = mi; A(2,2) = mi; A(3,3) = mi
    A(1,4) = mj; A(2,5) = mj; A(3,6) = mj
    A(1,7) = mk; A(2,8) = mk; A(3,9) = mk

    ! Angular momentum conservation
    A(4,2) = -mi * (ri(3) - Rcom(3)); A(4,3) = mi * (ri(2) - Rcom(2))
    A(4,5) = -mj * (rj(3) - Rcom(3)); A(4,6) = mj * (rj(2) - Rcom(2))
    A(4,8) = -mk * (rk(3) - Rcom(3)); A(4,9) = mk * (rk(2) - Rcom(2))

    A(5,1) = mi * (ri(3) - Rcom(3)); A(5,3) = -mi * (ri(1) - Rcom(1))
    A(5,4) = mj * (rj(3) - Rcom(3)); A(5,6) = -mj * (rj(1) - Rcom(1))
    A(5,7) = mk * (rk(3) - Rcom(3)); A(5,9) = -mk * (rk(1) - Rcom(1))

    A(6,1) = -mi * (ri(2) - Rcom(2)); A(6,2) = mi * (ri(1) - Rcom(1))
    A(6,4) = -mj * (rj(2) - Rcom(2)); A(6,5) = mj * (rj(1) - Rcom(1))
    A(6,7) = -mk * (rk(2) - Rcom(2)); A(6,8) = mk * (rk(1) - Rcom(1))

    ! Bond length velocities
    A(7,1) = ri(1) - rj(1); A(7,2) = ri(2) - rj(2); A(7,3) = ri(3) - rj(3)
    A(7,4) = rj(1) - ri(1); A(7,5) = rj(2) - ri(2); A(7,6) = rj(3) - ri(3)
    b(7) = rij * vij

    A(8,1) = ri(1) - rk(1); A(8,2) = ri(2) - rk(2); A(8,3) = ri(3) - rk(3)
    A(8,7) = rk(1) - ri(1); A(8,8) = rk(2) - ri(2); A(8,9) = rk(3) - ri(3)
    b(8) = rik * vik

    A(9,4) = rj(1) - rk(1); A(9,5) = rj(2) - rk(2); A(9,6) = rj(3) - rk(3)
    A(9,7) = rk(1) - rj(1); A(9,8) = rk(2) - rj(2); A(9,9) = rk(3) - rj(3)
    b(9) = rjk * vjk

    ! Invert
    call matrix_inverse(A,A_inv)
    x = A_inv .mult. b

    call print('Momentum before update = '//momentum(this,(/i,j,k/)))

    call print('Angular momentum before update = '//&
         (mi*((ri-Rcom) .cross. this%atoms%velo(:,i)) + &
          mj*((rj-Rcom) .cross. this%atoms%velo(:,j)) + &
          mk*((rk-Rcom) .cross. this%atoms%velo(:,k))))

    dEij = -bond_energy(this,i,j)
    dEik = -bond_energy(this,i,k)
    dEjk = -bond_energy(this,j,k)

    this%atoms%velo(:,i) = this%atoms%velo(:,i) + x(1:3)
    this%atoms%velo(:,j) = this%atoms%velo(:,j) + x(4:6)
    this%atoms%velo(:,k) = this%atoms%velo(:,k) + x(7:9)

    dEij = dEij + bond_energy(this,i,j)
    dEik = dEik + bond_energy(this,i,k)
    dEjk = dEjk + bond_energy(this,j,k)

    call print('Momentum after update = '//momentum(this,(/i,j,k/)))

    call print('Angular momentum after update = '//&
         (mi*((ri-Rcom) .cross. this%atoms%velo(:,i)) + &
          mj*((rj-Rcom) .cross. this%atoms%velo(:,j)) + &
          mk*((rk-Rcom) .cross. this%atoms%velo(:,k))))

    ! Now check the result

    call print('Required velocity for bond '//i//'--'//j//': '//gij)
    call print('Required velocity for bond '//i//'--'//k//': '//gik)
    call print('Required velocity for bond '//j//'--'//k//': '//gjk)

    call distance_relative_velocity(this%atoms,i,j,rij,vij)
    call distance_relative_velocity(this%atoms,i,k,rik,vik)
    call distance_relative_velocity(this%atoms,j,k,rjk,vjk)
    
    call print('')

    call print('Actual velocity for bond '//i//'--'//j//': '//vij)
    call print('Actual velocity for bond '//i//'--'//k//': '//vik)
    call print('Actual velocity for bond '//j//'--'//k//': '//vjk)

    call print('Energy added to bond '//i//'--'//j//': '//dEij)
    call print('Energy added to bond '//i//'--'//k//': '//dEik)
    call print('Energy added to bond '//j//'--'//k//': '//dEjk)

    Etot = Etot + dEij + dEik + dEjk

    call print('Total energy added so far = '//Etot)

  end subroutine thermalize_group
  
  function bond_energy(this,i,j)

    type(DynamicalSystem), intent(in) :: this
    integer,               intent(in) :: i, j
    real(dp)                          :: bond_energy

    real(dp)                          :: mi, mj, vij, rij

    mi = ElementMass(this%atoms%Z(i))
    mj = ElementMass(this%atoms%Z(j))

    call distance_relative_velocity(this%atoms,i,j,rij,vij)

    bond_energy = 0.5_dp * mi * mj * vij * vij / (mi + mj)

  end function bond_energy

end module clusters_module
