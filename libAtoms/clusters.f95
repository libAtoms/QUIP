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

#include "error.inc"
module clusters_module
  use system_module
  use linearalgebra_module
  use Table_module
  use periodictable_module
  use Atoms_types_module
  use Atoms_module
  use dictionary_module
  use ParamReader_module
  use group_module
  use constraints_module
  use dynamicalsystem_module
  use cinoutput_module

implicit none
private

integer, parameter :: &
   HYBRID_ACTIVE_MARK = 1, &
   HYBRID_BUFFER_MARK = 2, &
   HYBRID_TRANS_MARK = 3, &
   HYBRID_TERM_MARK = 4, &
   HYBRID_BUFFER_OUTER_LAYER_MARK = 5, &
   HYBRID_ELECTROSTATIC_MARK = 6, &
   HYBRID_NO_MARK = 0
!   HYBRID_FIT_MARK = 6, &

public :: HYBRID_ACTIVE_MARK, HYBRID_BUFFER_MARK, HYBRID_TRANS_MARK, HYBRID_TERM_MARK, &
   HYBRID_BUFFER_OUTER_LAYER_MARK, HYBRID_ELECTROSTATIC_MARK, HYBRID_NO_MARK
! HYBRID_FIT_MARK

integer, parameter :: MAX_CUT_BONDS = 6
public :: MAX_CUT_BONDS

character(len=TABLE_STRING_LENGTH), parameter :: hybrid_mark_name(0:6) = &
  (/ "h_none    ", &
     "h_active  ", &
     "h_buffer  ", &
     "h_trans   ", &
     "h_term    ", &
     "h_outer_l ", &
     "h_fit     " /)

public :: create_cluster_info_from_mark, carve_cluster, create_cluster_simple, create_hybrid_weights, &
    bfs_grow, bfs_step, multiple_images, discard_non_min_images, make_convex, create_embed_and_fit_lists, &
    create_embed_and_fit_lists_from_cluster_mark, &
    add_cut_hydrogens, construct_hysteretic_region, &
    create_pos_or_list_centred_hybrid_region, get_hybrid_list, &
    estimate_origin_extent
    !, construct_region, select_hysteretic_quantum_region

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
  subroutine bfs_grow_single(this, list, atom, n, nneighb_only, min_images_only, alt_connect)
    type(Atoms), intent(in)  :: this
    type(Table), intent(out) :: list
    integer, intent(in) :: atom, n
    logical, optional, intent(in)::nneighb_only, min_images_only
    type(Connection), intent(in), optional :: alt_connect

    call append(list, (/atom, 0,0,0/)) ! Add atom with shift 000
    call bfs_grow_list(this, list, n, nneighb_only, min_images_only, alt_connect)

  end subroutine bfs_grow_single


  !% On exit, 'list' will have been grown by 'n' bond hops.
  subroutine bfs_grow_list(this, list, n, nneighb_only, min_images_only, alt_connect)
    type(Atoms), intent(in)   ::this
    type(Table), intent(inout)::list
    integer, intent(in) :: n
    logical, optional, intent(in)::nneighb_only, min_images_only
    type(Connection), intent(in), optional :: alt_connect

    type(Table)::tmplist
    integer::i

    do i=1,n
       call bfs_step(this, list, tmplist, nneighb_only, min_images_only, alt_connect=alt_connect)
       call append(list, tmplist)
    end do

    if (n >= 1) call finalise(tmplist)
  end subroutine bfs_grow_list


  !% Execute one Breadth-First-Search move on the atomic connectivity graph.
  subroutine bfs_step(this,input,output,nneighb_only, min_images_only, max_r, alt_connect, property, debugfile, error)
    type(Atoms),        intent(in), target      :: this  !% The atoms structure to perform the step on.
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

    real(dp), optional, intent(in)      :: max_r
    !% if present, only neighbors within this range will be included

    type(Connection), intent(in), optional, target :: alt_connect
    integer, intent(in), optional:: property(:)

    type(inoutput), optional :: debugfile

    integer, optional, intent(out) :: error

    !local
    logical                             :: do_nneighb_only, do_min_images_only
    integer                             :: i, j, n, m, jshift(3), ishift(3)

    integer :: n_i, keep_row(4), in_i, min_image
    integer, allocatable, dimension(:) :: repeats
    real(dp), allocatable, dimension(:) :: normsqshift
    type(Connection), pointer :: use_connect

    INIT_ERROR(error)

    if (present(debugfile)) call print("bfs_step", file=debugfile)
    if (present(alt_connect)) then
       use_connect => alt_connect
    else
      use_connect => this%connect
    endif

    if (present(debugfile)) call print("bfs_step cutoff " // this%cutoff // " " // this%use_uniform_cutoff, file=debugfile)

    if (.not.use_connect%initialised) then
         RAISE_ERROR('BFS_Step: Atomic structure has no connectivity data', error)
    endif

    do_nneighb_only = optional_default(.true., nneighb_only)

    do_min_images_only = optional_default(.false., min_images_only)

    if (present(debugfile)) call print('bfs_step: do_nneighb_only = ' // do_nneighb_only // ' do_min_images_only = '//do_min_images_only, file=debugfile)
    call print('bfs_step: do_nneighb_only = ' // do_nneighb_only // ' do_min_images_only = '//do_min_images_only, PRINT_NERD)

    if(input%intsize /= 4 .or. input%realsize /= 0) then
         RAISE_ERROR("bfs_step: input table must have intsize=4.", error)
    endif

    call allocate(output, 4, 0, 0, 0)

    ! Now go though the atomic indices
    do m = 1, input%N
       i = input%int(1,m)
       ishift = input%int(2:4,m)

       if (present(debugfile)) call print("bfs_step check atom " // m // " " // i // " with n_neigh " // n_neighbours(this, i, alt_connect=use_connect), file=debugfile)
       ! Loop over i's neighbours
       do n = 1, n_neighbours(this,i, alt_connect=use_connect)
	  j = neighbour(this,i,n,shift=jshift,max_dist=max_r, alt_connect=use_connect)
	  if (present(debugfile)) call print("bfs_step   check neighbour " // n // " " // j, file=debugfile)
	  if (j == 0) cycle

          ! Look at next neighbour if j with correct shift is already in the cluster
          ! Must check input AND output tables
          if (find(input,(/j,ishift+jshift/)) > 0) cycle
	  if (present(debugfile)) call print("bfs_step   not in input", file=debugfile)
          if (find(output,(/j,ishift+jshift/)) > 0) cycle
	  if (present(debugfile)) call print("bfs_step   not in output", file=debugfile)

          if (do_nneighb_only .and. .not. is_nearest_neighbour(this, i, n, alt_connect=use_connect)) cycle
	  if (present(debugfile)) call print("bfs_step   acceptably near neighbor", file=debugfile)

          if (present(property)) then
             if (property(j) == 0) cycle
          endif
	  if (present(debugfile)) call print("bfs_step   property matches OK", file=debugfile)

          ! Everything checks out ok, so add j to the output table
          ! with correct shift
	  if (present(debugfile)) call print("bfs_step   appending", file=debugfile)
          call append(output,(/j,ishift+jshift/)) 
       end do
    end do

    if (present(debugfile)) call print("bfs_step   raw output", file=debugfile)
    if (present(debugfile)) call print(output, file=debugfile)

    if (do_min_images_only) then

       ! If there are repeats of any atomic index, 
       ! we want to keep the one with smallest normsq(shift)

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

          allocate(repeats(n_i), normsqshift(n_i))

          ! Get indices of repeats of this atomic index
          repeats = pack((/ (j, j=1,output%N) /), &
               int_part(output,1) == output%int(1,n))

          if (in_i /= 0) then
             ! atom is in input list, remove all new occurances
             call delete_multiple(output, repeats)
          else
             ! Find row with minimum normsq(shift)
             normsqshift = normsq(real(output%int(2:4,repeats),dp),1)
             min_image = repeats(minloc(normsqshift,dim=1))
             keep_row = output%int(:,min_image)

             ! keep the minimum image
             call delete_multiple(output, &
                  pack(repeats, repeats /= min_image))
          end if

          ! don't increment n, since delete_multiple copies items from
          ! end of list over deleted items, so we need to retest

          deallocate(repeats, normsqshift)

       end do

    end if

    if (present(debugfile)) call print("bfs_step   final output", file=debugfile)
    if (present(debugfile)) call print(output, file=debugfile)

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
    real(dp), allocatable, dimension(:) :: normsqshift

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

       allocate(repeats(n_i), normsqshift(n_i))

       ! Get indices of repeats of this atomic index
       repeats = pack((/ (j, j=1,list%N) /), &
            int_part(list,1) == list%int(1,n))

       ! Find row with minimum normsq(shift)
       normsqshift = normsq(real(list%int(2:4,repeats),dp),1)
       min_image = repeats(minloc(normsqshift,dim=1))
       keep_row = list%int(:,min_image)

       ! keep the minimum image
       call delete_multiple(list, &
            pack(repeats, repeats /= min_image))

       ! don't increment n, since delete_multiple copies items from
       ! end of list over deleted items, so we need to retest

       deallocate(repeats, normsqshift)

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
  function make_convex_step(this, input, output, error) result(newatoms)
    type(Atoms), intent(in) :: this
    type(Table), intent(in) :: input
    type(Table), intent(out) :: output
    integer, optional, intent(out) :: error
    integer::newatoms

    integer :: n, i, j, k, p, m, n_in, nn, ishift(3), jshift(3), kshift(3)
    real(dp) :: r_ij, r_kj

    INIT_ERROR(error)

    ! Check table size
    if(input%intsize /= 4 .or. input%realsize /= 0) then
         RAISE_ERROR("make_convex_step: input table must have intsize=4", error)
    endif

    if(input%intsize /= 4 .or. input%realsize /= 0) then
         RAISE_ERROR("bfs_step: input table must have intsize=4.", error)
    endif

    call allocate(output, 4, 0, 0, 0)

    do n=1,input%N
       i = input%int(1,n)
       ishift = input%int(2:4,n)

       !Loop over neighbours
       do m = 1, n_neighbours(this,i)
          j = neighbour(this,i,m, r_ij,shift=jshift)

          ! Look for nearest neighbours of i not in input list
          if (find(input,(/j,ishift+jshift/)) == 0 .and. is_nearest_neighbour(this, i, m)) then

             n_in = 0
             nn = 0
             ! Count number of nearest neighbours of j, and how many
             ! of them are in input list
             do p = 1, n_neighbours(this,j)
                k = neighbour(this,j,p, r_kj, shift=kshift)
                if (is_nearest_neighbour(this, j, p)) then
                   nn = nn + 1
                   if (find(input,(/k,ishift+jshift+kshift/)) /= 0) n_in = n_in + 1
                end if
             end do

             !If more than half of  j's nearest neighbours are in then add it to output
             if (find(output, (/j,ishift+jshift/)) == 0 .and. real(n_in,dp)/real(nn,dp) > 0.5_dp) &
                  call append(output, (/j,ishift+jshift/))

          end if

       end do
    end do

    newatoms = output%N

  end function  make_convex_step


! Gotcha 1: Hollow sections 
! NB equivalent to reduce_n_cut_bonds when new number of bonds is 0 
! JRK: reduce_n_cut_bonds() is not equivalent to cluster_in_out_in()
! for silica, where adding IN-OUT-IN atoms doesn't change number of
! cut bonds OUT and IN refers to the list in cluster_info Look at the
! OUT nearest neighbours of IN atoms. If mode='all', then if ALL the
! nearest neighbours of the OUT atom are IN, then make the OUT atom
! IN. If mode='any', then add the OUT atom if ANY of it's nearest
! neighours (apart from the original IN atom)!  are IN. This stronger
! condition is required for cp2k.

function cluster_in_out_in(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, mode) result(cluster_changed)
  type(Atoms), intent(in) :: this
  type(Table), intent(inout) :: cluster_info
  logical, intent(in) :: connectivity_just_from_connect
  type(Connection), intent(in) :: use_connect
  logical, intent(in) :: atom_mask(6)
  character(len=*), intent(in) :: mode
  logical :: cluster_changed

  integer :: n, i, ishift(3), m, j, jshift(3), p, k, kshift(3), n_nearest, n_in
  logical :: add_atom

  cluster_changed = .false.
  if (trim(mode) /= 'all' .and. trim(mode) /= 'any') &
       call system_abort('cluster_in_out_in: bad mode '//mode//' - should be "all" or "any"')

  n = 1
  ! Loop over cluster atoms (including ones that may get added in this loop)
  call print('cluster_in_out_in: Checking for hollow sections', PRINT_NERD)
  do while (n <= cluster_info%N)
    i = cluster_info%int(1,n)
    ishift = cluster_info%int(2:4,n)
    call print('cluster_in_out_in: i = '//i//' ['//ishift//'] Looping over '//n_neighbours(this,i,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)

    !Loop over neighbours
    do m = 1, n_neighbours(this,i,alt_connect=use_connect)
    j = neighbour(this,i,m, shift=jshift,alt_connect=use_connect)

    if (find(cluster_info,(/j,ishift+jshift,this%Z(j),0/), atom_mask) == 0 .and. &
        (connectivity_just_from_connect .or. is_nearest_neighbour(this, i, m, alt_connect=use_connect)) ) then
      ! j is out and is nearest neighbour

      call print('cluster_in_out_in:   checking j = '//j//" ["//jshift//"]",PRINT_ANAL)

      ! We have an OUT nearest neighbour, loop over its nearest neighbours to see if they are IN

      n_nearest = 0
      n_in = 0
      do p = 1, n_neighbours(this,j,alt_connect=use_connect)
         k = neighbour(this,j,p, shift=kshift,alt_connect=use_connect)
         if (k == i) cycle
         if (connectivity_just_from_connect .or. is_nearest_neighbour(this, j, p,alt_connect=use_connect)) then
            n_nearest = n_nearest + 1
            if (find(cluster_info,(/k,ishift+jshift+kshift,this%Z(k),0/), atom_mask) /= 0) then
               n_in = n_in + 1
            end if
         end if
      end do

      if (trim(mode) == "all") then
         !If all of j's nearest neighbours are IN then add it
         add_atom = n_nearest == n_in
      else
         !If any of j's nearest neighbours (apart from i) are in then add it
         add_atom = n_in /= 0
      end if

      if (add_atom) then
        call append(cluster_info, (/j,ishift+jshift,this%Z(j),0/), (/this%pos(:,j), 1.0_dp/), (/ "inoutin   "/) )
        cluster_changed = .true.
        call print('cluster_in_out_in:  Added atom ' //j//' ['//(ishift+jshift)//'] to cluster. Atoms = ' // cluster_info%N, PRINT_NERD)
      end if

    end if

    end do ! m
    n = n + 1
  end do ! while (n <= cluster_info%N)

  call print('cluster_in_out_in: Finished checking',PRINT_NERD)
  call print("cluster_in_out_in: cluster list:", PRINT_NERD)
  call print(cluster_info, PRINT_NERD)
end function cluster_in_out_in

  !% Find cases where two IN atoms have a common
  !% OUT nearest neighbour, and see if termination would cause the hydrogen
  !% atoms to be too close together. If so, include the OUT nearest neighbour
  !% in the cluster
  !% returns true if cluster was changed
  function cluster_fix_termination_clash(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, termination_clash_factor) result(cluster_changed)
    type(Atoms), intent(in) :: this !% atoms structure 
    type(Table), intent(inout) :: cluster_info !% table of cluster info, modified if necessary on output
    logical, intent(in) :: connectivity_just_from_connect !% if true, we're doing hysterestic connect and should rely on the connection object completely
    type(Connection), intent(in) :: use_connect !% connection object to use for connectivity info
    logical, intent(in) :: atom_mask(6) !% which fields in int part of table to compare when checking for identical atoms
    real(dp), intent(in) :: termination_clash_factor
    logical :: cluster_changed

    integer :: n, i, ishift(3), m, j, jshift(3), p, k, kshift(3)
    real(dp) :: dhat_ij(3), dhat_jk(3), r_ij, r_jk, H1(3), H2(3), diff_ik(3)
    ! real(dp) :: t_norm

    cluster_changed = .false.

    call print('doing cluster_fix_termination_clash', PRINT_NERD)

    !Loop over atoms in the cluster
    n = 0
    do while (n < cluster_info%N)
      n = n + 1
      i = cluster_info%int(1,n)     ! index of atom in the cluster
      ishift = cluster_info%int(2:4,n)
      call print('cluster_fix_termination_clash: i = '//i//'. Looping over '//n_neighbours(this,i,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)
      !Loop over atom i's neighbours
      do m = 1, n_neighbours(this,i,alt_connect=use_connect)
	j = neighbour(this,i,m, shift=jshift, diff=dhat_ij, distance=r_ij,alt_connect=use_connect)
	dhat_ij = dhat_ij/r_ij

	!If j is IN the cluster, or not a nearest neighbour then try the next neighbour
	if(find(cluster_info,(/j,ishift+jshift,this%Z(j),0/), atom_mask) /= 0) then
	  call print('cluster_fix_termination_clash:   j = '//j//" ["//jshift//"] is in cluster",PRINT_ANAL)
	  cycle
	end if
	if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, m, alt_connect=use_connect))) then
	  call print('cluster_fix_termination_clash:   j = '//j//" ["//jshift//"] not nearest neighbour",PRINT_ANAL)
	  cycle
	end if

	! So j is an OUT nearest neighbour of i.
	call print('cluster_fix_termination_clash:   checking j = '//j//" ["//jshift//"]",PRINT_ANAL)

	!Determine the position of the would-be hydrogen along i--j
	call print('cluster_fix_termination_clash:  Finding i--j hydrogen position',PRINT_ANAL)

	H1 = this%pos(:,i) + (this%lattice .mult. (ishift)) + &
	     termination_bond_rescale(this%Z(i), this%Z(j)) * r_ij * dhat_ij

	!Do a loop over j's nearest neighbours
	call print('cluster_fix_termination_clash:  Looping over '//n_neighbours(this,j, alt_connect=use_connect)//' neighbours of j', PRINT_ANAL)

	do p = 1, n_neighbours(this,j, alt_connect=use_connect)
	  k = neighbour(this,j,p, shift=kshift, diff=dhat_jk, distance=r_jk, alt_connect=use_connect)
	  dhat_jk = dhat_jk/r_jk

	  !If k is OUT of the cluster or k == i or it is not a nearest neighbour of j
	  !then try the next neighbour

	  if(find(cluster_info,(/k,ishift+jshift+kshift,this%Z(k),0/), atom_mask) == 0) then
	    call print('cluster_fix_termination_clash:   k = '//k//" ["//kshift//"] not in cluster",PRINT_ANAL)
	    cycle
	  end if
	  if(k == i .and. all( jshift+kshift == 0 )) cycle
	  if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,j, p, alt_connect=use_connect))) then
	    call print('cluster_fix_termination_clash:   k = '//k//" ["//kshift//"] not nearest neighbour",PRINT_ANAL)
	    cycle
	  end if

	  call print('cluster_fix_termination_clash: testing k = '//k//" ["//kshift//"]", PRINT_ANAL)
	  !Determine the position of the would-be hydrogen along k--j
	  call print('cluster_fix_termination_clash:   Finding k--j hydrogen position',PRINT_ANAL)

	  diff_ik = r_ij * dhat_ij + r_jk * dhat_jk

	  H2 = this%pos(:,i) + (this%lattice .mult. (ishift)) + diff_ik - &
	    termination_bond_rescale(this%Z(k), this%Z(j)) * r_jk * dhat_jk
	  call print('cluster_fix_termination_clash:   Checking i--k distance and hydrogen distance',PRINT_ANAL)

	  call print("cluster_fix_termination_clash: proposed hydrogen positions:", PRINT_ANAL)
	  call print(H1, PRINT_ANAL)
	  call print(H2, PRINT_ANAL)
	  !NB workaround for pgf90 bug (as of 9.0-1)
	  ! t_norm = norm(H1-H2); call print("cluster_fix_termination_clash: hydrogen distance would be "//t_norm, PRINT_ANAL)
	  !NB end of workaround for pgf90 bug (as of 9.0-1)
	  ! If i and k are nearest neighbours, or the terminating hydrogens would be very close, then
	  ! include j in the cluster. The H--H checking is conservative, hence the extra factor of 1.2
	  if ((norm(diff_ik) < bond_length(this%Z(i),this%Z(k))*this%nneightol) .or. &
	      (norm(H1-H2) < bond_length(1,1)*this%nneightol*termination_clash_factor) ) then

	    call append(cluster_info,(/j,ishift+jshift,this%Z(j),0/),(/this%pos(:,j),1.0_dp/), (/ "clash     "/) )
	    cluster_changed = .true.
	    call print('cluster_fix_termination_clash:  Atom '//j//'with shift '//(ishift+jshift)//' added to cluster. Atoms = '//cluster_info%N, PRINT_NERD)
	    ! j is now included in the cluster, so we can exit this do loop (over p)
	    exit
	  end if 

	end do ! p

      end do ! m

    end do ! while (n <= cluster_info%N)

    call print('cluster_fix_termination_clash: Finished checking',PRINT_NERD)
    call print("cluster_fix_termination_clash: cluster list:", PRINT_NERD)
    call print(cluster_info, PRINT_NERD)
  end function cluster_fix_termination_clash

  !% keep each whole residue (using atom_res_number(:) property) if any bit of it is already included
  function cluster_keep_whole_residues(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, present_keep_whole_residues) result(cluster_changed)
    type(Atoms), intent(in) :: this !% atoms structure 
    type(Table), intent(inout) :: cluster_info !% table of cluster info, modified if necessary on output
    logical, intent(in) :: connectivity_just_from_connect !% if true, we're doing hysterestic connect and should rely on the connection object completely
    type(Connection), intent(in) :: use_connect !% connection object to use for connectivity info
    logical, intent(in) :: atom_mask(6) !% which fields in int part of table to compare when checking for identical atoms
    logical, intent(in) :: present_keep_whole_residues !% if true, keep_whole_residues was specified by user, not just a default value
    logical :: cluster_changed
    
    integer :: n, i, ishift(3), m, j, jshift(3)
    integer, pointer :: atom_res_number(:)

    call print('doing cluster_keep_whole_residues', PRINT_NERD)

    cluster_changed = .false.
    if (.not. assign_pointer(this, 'atom_res_number', atom_res_number)) then
      if (present_keep_whole_residues) then
	call print("WARNING: cluster_keep_whole_residues got keep_whole_residues requested explicitly, but no proper atom_res_number property available", PRINT_ALWAYS)
      endif
      return
    endif

    n = 0
    do while (n < cluster_info%N)
      n = n + 1
      i = cluster_info%int(1,n)
      ishift = cluster_info%int(2:4,n)
      if (atom_res_number(i) < 0) cycle
      call print('cluster_keep_whole_residues: i = '//i//' residue # = ' // atom_res_number(i) //'. Looping over '//n_neighbours(this,i,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)
      do m=1, n_neighbours(this, i, alt_connect=use_connect)
	j = neighbour(this, i, m, shift=jshift, alt_connect=use_connect)

	if (atom_res_number(i) /= atom_res_number(j)) then
	  call print("cluster_keep_whole_residues:   j = "//j//" ["//jshift//"] has different res number " // atom_res_number(j), PRINT_ANAL)
	  cycle
	endif
	if(find(cluster_info,(/j,ishift+jshift,this%Z(j),0/), atom_mask) /= 0) then
	  call print("cluster_keep_whole_residues:   j = "//j//" ["//jshift//"] is in cluster",PRINT_ANAL)
	  cycle
	end if
	if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, m, alt_connect=use_connect))) then
	  call print("cluster_keep_whole_residues:   j = "//j//" ["//jshift//"] not nearest neighbour",PRINT_ANAL)
	  cycle
	end if

	call append(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), (/ this%pos(:,j), 1.0_dp /), (/"res_num   "/) )
	cluster_changed = .true.
	call print('cluster_keep_whole_residues:  Added atom ' //j//' ['//(ishift+jshift)//'] to cluster. Atoms = ' // cluster_info%N, PRINT_NERD)

      end do ! m
    end do ! while (n <= cluster_info%N)

    call print('cluster_keep_whole_residues: Finished checking',PRINT_NERD)
    call print("cluster_keep_whole_residues: cluster list:", PRINT_NERD)
    call print(cluster_info, PRINT_NERD)

  end function cluster_keep_whole_residues

  !% keep each whole proline (using atom_subgroup_number(:) and atom_res_number(:) property) if any bit of it is already included
  function cluster_keep_whole_prolines(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, present_keep_whole_prolines) result(cluster_changed)
    type(Atoms), intent(in) :: this !% atoms structure 
    type(Table), intent(inout) :: cluster_info !% table of cluster info, modified if necessary on output
    logical, intent(in) :: connectivity_just_from_connect !% if true, we're doing hysterestic connect and should rely on the connection object completely
    type(Connection), intent(in) :: use_connect !% connection object to use for connectivity info
    logical, intent(in) :: atom_mask(6) !% which fields in int part of table to compare when checking for identical atoms
    logical, intent(in) :: present_keep_whole_prolines !% if true, keep_whole_prolines was specified by user, not just a default value
    logical :: cluster_changed
    
    integer :: n, i, ishift(3), m, j, jshift(3)
    integer, pointer :: atom_res_number(:), atom_subgroup_number(:)
    logical :: is_proline

    call print('doing cluster_keep_whole_prolines', PRINT_NERD)

    cluster_changed = .false.
    if (.not. assign_pointer(this, 'atom_res_number', atom_res_number) .or. .not. assign_pointer(this, 'atom_subgroup_number', atom_subgroup_number)) then
      if (present_keep_whole_prolines) then
	call print("WARNING: cluster_keep_whole_prolines got keep_whole_prolines requested explicitly, but no proper atom_res_number or atom_subgroup_number properties available", PRINT_ALWAYS)
      endif
      return
    endif

    n = 0
    do while (n < cluster_info%N)
      n = n + 1
      i = cluster_info%int(1,n)
      ishift = cluster_info%int(2:4,n)
      if (atom_res_number(i) < 0) cycle
      call print('cluster_keep_whole_prolines: i = '//i//' residue # = ' // atom_res_number(i) //' subgroup # = '// atom_subgroup_number(i) // &
                 '. Looping over '//n_neighbours(this,i,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)
      do m=1, n_neighbours(this, i, alt_connect=use_connect)
	j = neighbour(this, i, m, shift=jshift, alt_connect=use_connect)
	is_proline = .false.

	if (atom_res_number(i) /= atom_res_number(j)) then
	  call print("cluster_keep_whole_prolines:   j = "//j//" ["//jshift//"] has different res number " // atom_res_number(j), PRINT_ANAL)
	  cycle
	endif
	if (atom_subgroup_number(i) <= -10 .or. atom_subgroup_number(j) <= -10) then
	  call print("cluster_keep_whole_prolines:   j = "//j//" ["//jshift//"] has proline subgroup number " // atom_subgroup_number(j), PRINT_ANAL)
	  is_proline = .true.
	endif
	if(find(cluster_info,(/j,ishift+jshift,this%Z(j),0/), atom_mask) /= 0) then
	  call print("cluster_keep_whole_prolines:   j = "//j//" ["//jshift//"] is in cluster",PRINT_ANAL)
	  cycle
	end if
	if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, m, alt_connect=use_connect))) then
	  call print("cluster_keep_whole_prolines:   j = "//j//" ["//jshift//"] not nearest neighbour",PRINT_ANAL)
	  cycle
	end if

	if (is_proline) then
	  call append(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), (/ this%pos(:,j), 1.0_dp /), (/"subgp_num "/) )
	  cluster_changed = .true.
	  call print('cluster_keep_whole_prolines:  Added atom ' //j//' ['//(ishift+jshift)//'] to cluster. Atoms = ' // cluster_info%N, PRINT_NERD)
	endif

      end do ! m

    end do ! while (n <= cluster_info%N)

    call print('cluster_keep_whole_prolines: Finished checking',PRINT_NERD)
    call print("cluster_keep_whole_prolines: cluster list:", PRINT_NERD)
    call print(cluster_info, PRINT_NERD)

  end function cluster_keep_whole_prolines

  !% keep each proline sidechain and directly connected bit of backbone (using atom_subgroup_number and atom_res_number(:) property) if any bit of it is already included
  function cluster_keep_whole_proline_sidechains(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, present_keep_whole_proline_sidechains) result(cluster_changed)
    type(Atoms), intent(in) :: this !% atoms structure 
    type(Table), intent(inout) :: cluster_info !% table of cluster info, modified if necessary on output
    logical, intent(in) :: connectivity_just_from_connect !% if true, we're doing hysterestic connect and should rely on the connection object completely
    type(Connection), intent(in) :: use_connect !% connection object to use for connectivity info
    logical, intent(in) :: atom_mask(6) !% which fields in int part of table to compare when checking for identical atoms
    logical, intent(in) :: present_keep_whole_proline_sidechains !% if true, keep_whole_proline_sidechains was specified by user, not just a default value
    logical :: cluster_changed
    
    integer :: n, i, ishift(3), m, j, jshift(3)
    integer, pointer :: atom_res_number(:), atom_subgroup_number(:)
    logical :: is_proline

    call print('doing cluster_keep_whole_proline_sidechains', PRINT_NERD)

    cluster_changed = .false.
    if (.not. assign_pointer(this, 'atom_res_number', atom_res_number) .or. .not. assign_pointer(this, 'atom_subgroup_number', atom_subgroup_number)) then
      if (present_keep_whole_proline_sidechains) then
	call print("WARNING: cluster_keep_whole_proline_sidechains got keep_whole_proline_sidechains requested explicitly, but no proper atom_res_number or atom_subgroup_number properties available", PRINT_ALWAYS)
      endif
      return
    endif

    n = 0
    do while (n < cluster_info%N)
      n = n + 1
      i = cluster_info%int(1,n)
      ishift = cluster_info%int(2:4,n)
      if (atom_res_number(i) < 0) cycle
      call print('cluster_keep_whole_proline_sidechains: i = '//i//' residue # = ' // atom_res_number(i) //' subgroup # = '// atom_subgroup_number(i) // &
                 '. Looping over '//n_neighbours(this,i,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)
      do m=1, n_neighbours(this, i, alt_connect=use_connect)
	j = neighbour(this, i, m, shift=jshift, alt_connect=use_connect)
	is_proline = .false.

	if (atom_res_number(i) /= atom_res_number(j)) then
	  call print("cluster_keep_whole_proline_sidechains:   j = "//j//" ["//jshift//"] has different res number " // atom_res_number(j), PRINT_ANAL)
	  cycle
	endif
	if (atom_subgroup_number(i) == -1000 .or. atom_subgroup_number(j) == -1000) then
	  call print("cluster_keep_proline_sidechains:   i or j = "//j//" ["//jshift//"] has proline sidechain subgroup number " // atom_subgroup_number(j), PRINT_ANAL)
	  is_proline = .true.
	endif
	if(find(cluster_info,(/j,ishift+jshift,this%Z(j),0/), atom_mask) /= 0) then
	  call print("cluster_keep_proline_sidechains:   j = "//j//" ["//jshift//"] is in cluster",PRINT_ANAL)
	  cycle
	end if
	if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, m, alt_connect=use_connect))) then
	  call print("cluster_keep_whole_proline_sidechains:   j = "//j//" ["//jshift//"] not nearest neighbour",PRINT_ANAL)
	  cycle
	end if

	if (is_proline) then
	  call append(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), (/ this%pos(:,j), 1.0_dp /), (/"subgp_num "/) )
	  cluster_changed = .true.
	  call print('cluster_keep_whole_proline_sidechains:  Added atom ' //j//' ['//(ishift+jshift)//'] to cluster. Atoms = ' // cluster_info%N, PRINT_NERD)
	endif

      end do ! m

    end do ! while (n <= cluster_info%N)

    call print('cluster_keep_whole_proline_sidechains: Finished checking',PRINT_NERD)
    call print("cluster_keep_whole_proline_sidechains: cluster list:", PRINT_NERD)
    call print(cluster_info, PRINT_NERD)

  end function cluster_keep_whole_proline_sidechains

  !% keep each whole subgroup (using atom_subgroup_number and atom_res_number(:) property) if any bit of it is already included
  function cluster_keep_whole_subgroups(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, present_keep_whole_subgroups) result(cluster_changed)
    type(Atoms), intent(in) :: this !% atoms structure 
    type(Table), intent(inout) :: cluster_info !% table of cluster info, modified if necessary on output
    logical, intent(in) :: connectivity_just_from_connect !% if true, we're doing hysterestic connect and should rely on the connection object completely
    type(Connection), intent(in) :: use_connect !% connection object to use for connectivity info
    logical, intent(in) :: atom_mask(6) !% which fields in int part of table to compare when checking for identical atoms
    logical, intent(in) :: present_keep_whole_subgroups !% if true, keep_whole_subgroups was specified by user, not just a default value
    logical :: cluster_changed
    
    integer :: n, i, ishift(3), m, j, jshift(3)
    integer, pointer :: atom_res_number(:), atom_subgroup_number(:)

    call print('doing cluster_keep_whole_subgroups', PRINT_NERD)

    cluster_changed = .false.
    if (.not. assign_pointer(this, 'atom_res_number', atom_res_number) .or. .not. assign_pointer(this, 'atom_subgroup_number', atom_subgroup_number)) then
      if (present_keep_whole_subgroups) then
	call print("WARNING: cluster_keep_whole_subgroups got keep_whole_subgroups requested explicitly, but no proper atom_res_number or atom_subgroup_number properties available", PRINT_ALWAYS)
      endif
      return
    endif

    n = 0
    do while (n < cluster_info%N)
      n = n + 1
      i = cluster_info%int(1,n)
      ishift = cluster_info%int(2:4,n)
      if (atom_res_number(i) < 0) cycle
      if (atom_subgroup_number(i) == 0) cycle
      call print('cluster_keep_whole_subgroups: i = '//i//' residue # = ' // atom_res_number(i) //' subgroup # = '// atom_subgroup_number(i) // '. Looping over '//n_neighbours(this,i,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)
      do m=1, n_neighbours(this, i, alt_connect=use_connect)
	j = neighbour(this, i, m, shift=jshift, alt_connect=use_connect)
	call print("cluster_keep_whole_subgroups:    m = " // m //" j = " // j, PRINT_ANAL)
	if (atom_res_number(j) < 0) cycle
	if (atom_subgroup_number(j) == 0) cycle

	if (atom_res_number(i) /= atom_res_number(j)) then
	  call print("cluster_keep_whole_subgroups:   j = "//j//" ["//jshift//"] has different res number " // atom_res_number(j), PRINT_ANAL)
	  cycle
	endif
	if (atom_subgroup_number(i) /= atom_subgroup_number(j)) then
	  call print("cluster_keep_whole_subgroups:   j = "//j//" ["//jshift//"] has different subgroup number " // atom_subgroup_number(j), PRINT_ANAL)
	  cycle
	endif
	if(find(cluster_info,(/j,ishift+jshift,this%Z(j),0/), atom_mask) /= 0) then
	  call print("cluster_keep_whole_subgroups:   j = "//j//" ["//jshift//"] is in cluster",PRINT_ANAL)
	  cycle
	end if
	if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, m, alt_connect=use_connect))) then
	  call print("cluster_keep_whole_subgroups:   j = "//j//" ["//jshift//"] not nearest neighbour",PRINT_ANAL)
	  cycle
	end if

	call append(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), (/ this%pos(:,j), 1.0_dp /), (/"subgp_num "/) )
	cluster_changed = .true.
	call print('cluster_keep_whole_subgroups:  Added atom ' //j//' ['//(ishift+jshift)//'] to cluster. Atoms = ' // cluster_info%N, PRINT_NERD)

      end do ! m

    end do ! while (n <= cluster_info%N)

    call print('cluster_keep_whole_subgroups: Finished checking',PRINT_NERD)
    call print("cluster_keep_whole_subgroups: cluster list:", PRINT_NERD)
    call print(cluster_info, PRINT_NERD)

  end function cluster_keep_whole_subgroups

  !% keep whole silica tetrahedra -- that is, for each silicon atom, keep all it's oxygen nearest neighbours
  function cluster_keep_whole_silica_tetrahedra(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask) result(cluster_changed)
    type(Atoms), intent(in) :: this !% atoms structure 
    type(Table), intent(inout) :: cluster_info !% table of cluster info, modified if necessary on output
    logical, intent(in) :: connectivity_just_from_connect !% if true, we're doing hysterestic connect and should rely on the connection object completely
    type(Connection), intent(in) :: use_connect !% connection object to use for connectivity info
    logical, intent(in) :: atom_mask(6) !% which fields in int part of table to compare when checking for identical atoms
    logical :: cluster_changed
    
    integer :: n, i, ishift(3), m, j, jshift(3)

    call print('doing cluster_keep_whole_silica_tetrahedra', PRINT_NERD)

    cluster_changed = .false.
    n = 0
    do while (n < cluster_info%N)
      n = n + 1
      i = cluster_info%int(1,n)
      ishift = cluster_info%int(2:4,n)
      if (this%z(i) /= 14) then
         ! Consider only silicon atoms, which form centres of tetrahedra
         cycle  
      end if

      call print('cluster_keep_whole_silica_tetrahedra: i = '//i//'. Looping over '//n_neighbours(this,i,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)
      do m=1, n_neighbours(this, i, alt_connect=use_connect)
	j = neighbour(this, i, m, shift=jshift, alt_connect=use_connect)

	if (this%z(j) /= 8) then
	  call print("cluster_keep_whole_silica_tetrahedra:   j = "//j//" ["//jshift//"] is not oxygen ", PRINT_ANAL)
	  cycle
	endif
	if(find(cluster_info,(/j,ishift+jshift,this%Z(j),0/), atom_mask) /= 0) then
	  call print("cluster_keep_whole_silica_tetrahedra:   j = "//j//" ["//jshift//"] is in cluster",PRINT_ANAL)
	  cycle
	end if
	if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, m, alt_connect=use_connect))) then
	  call print("cluster_keep_whole_silica_tetrahedra:   j = "//j//" ["//jshift//"] not nearest neighbour",PRINT_ANAL)
	  cycle
	end if

	call append(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), (/ this%pos(:,j), 1.0_dp /), (/"tetra     "/) )
	cluster_changed = .true.
	call print('cluster_keep_whole_silica_tetrahedra:  Added atom ' //j//' ['//(ishift+jshift)//'] to cluster. Atoms = ' // cluster_info%N, PRINT_NERD)

      end do ! m

    end do ! while (n <= cluster_info%N)

    call print('cluster_keep_whole_silica_tetrahedra: Finished checking',PRINT_NERD)
    call print("cluster_keep_whole_silica_tetrahedra: cluster list:", PRINT_NERD)
    call print(cluster_info, PRINT_NERD)

  end function cluster_keep_whole_silica_tetrahedra

  !% keep whole octahedra for titania that is, for each Ti atom, keep all it's oxygen nearest neighbours but the O atoms with only 1 neighbour.
  function cluster_keep_whole_titania_octahedra(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask) result(cluster_changed)
    type(Atoms), intent(in) :: this !% atoms structure 
    type(Table), intent(inout) :: cluster_info !% table of cluster info, modified if necessary on output
    logical, intent(in) :: connectivity_just_from_connect !% if true, we're doing hysterestic connect and should rely on the connection object completely
    type(Connection), intent(in) :: use_connect !% connection object to use for connectivity info
    logical, intent(in) :: atom_mask(6) !% which fields in int part of table to compare when checking for identical atoms
    logical :: cluster_changed
    
    integer :: n, i, ishift(3), m, j, jshift(3)

    call print('doing cluster_keep_whole_titania_octahedra', PRINT_NERD)

    cluster_changed = .false.
    n = 0
    do while (n < cluster_info%N)
      n = n + 1
      i = cluster_info%int(1,n)
      ishift = cluster_info%int(2:4,n)
      if (this%z(i) /= 22) then
         ! Consider only Ti atoms, which form centres of tetrahedra
         cycle  
      end if

      call print('cluster_keep_whole_titania_octahedra: i = '//i//'. Looping over '//n_neighbours(this,i,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)
      do m=1, n_neighbours(this, i, alt_connect=use_connect)
        j = neighbour(this, i, m, shift=jshift, alt_connect=use_connect)

        if (this%z(j) /= 8) then
          call print("cluster_keep_whole_titania_octahedra:   j = "//j//" ["//jshift//"] is not oxygen ", PRINT_ANAL)
          cycle
        endif
        if(find(cluster_info,(/j,ishift+jshift,this%Z(j),0/), atom_mask) /= 0) then
          call print("cluster_keep_whole_titania_octahedra:   j = "//j//" ["//jshift//"] is in cluster",PRINT_ANAL)
          cycle
        end if
        if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, m, alt_connect=use_connect))) then
          call print("cluster_keep_whole_titania_octahedra:   j = "//j//" ["//jshift//"] not nearest neighbour",PRINT_ANAL)
          cycle
        end if

        !check if the oxygen has more than 1 neighbour in the cluster.
!       neigh_O = 0
!       do l=1,n_neighbours(this, j, alt_connect=use_connect)
!           k = neighbour(this, j, l, distance=rdist, shift=kshift, alt_connect=use_connect)
!           if(find(cluster_info,(/k,ishift+jshift+kshift,this%Z(k),0/), atom_mask) == 0) then
!             call print(" j : "// j // " Z(j) :" // this%Z(j) // ' shift ' // jshift)
!             call print(" k : "// k // " Z(k) :" // this%Z(k) // ' shift ' // kshift)
!             call print(" i : "// i // " Z(i) :" // this%Z(i) // ' shift ' // ishift)
!             call print('distance between k and j '// rdist)
!             call print("cluster_keep_whole_titania_octahedra:   k = "//k//" ["//kshift//"] is not in cluster",PRINT_ANAL)
!             cycle
!           end if
!           if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,j, l, alt_connect=use_connect))) then
!             call print(" j : "// j // " Z(j) :" // this%Z(j) // ' shift ' // jshift)
!             call print(" k : "// k // " Z(k) :" // this%Z(k) // ' shift ' // kshift)
!             call print(" i : "// i // " Z(i) :" // this%Z(i) // ' shift ' // ishift)
!             call print('distance between k and j '// rdist)
!             call print("cluster_keep_whole_titania_octahedra:   k = "//k//" ["//kshift//"] not nearest neighbour",PRINT_ANAL)
!             cycle
!           end if
!           neigh_O = neigh_O + 1 
!       end do ! l

        call append(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), (/ this%pos(:,j), 1.0_dp /), (/"octa      "/) )
        cluster_changed = .true.
        call print('cluster_keep_whole_titania_octahedra  Added atom ' //j//' ['//(ishift+jshift)//'] to cluster. Atoms = ' // cluster_info%N, PRINT_NERD)
        
      end do ! m


    end do ! while (n <= cluster_info%N)

    call print('cluster_keep_whole_titania_octahedra: Finished checking',PRINT_NERD)
    call print("cluster_keep_whole_titania_octahedra: cluster list:", PRINT_NERD)
    call print(cluster_info, PRINT_NERD)

  end function cluster_keep_whole_titania_octahedra

  !% adding each neighbour outside atom if doing so immediately reduces the number of cut bonds
  function cluster_reduce_n_cut_bonds(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask) result(cluster_changed)
    type(Atoms), intent(in) :: this !% atoms structure 
    type(Table), intent(inout) :: cluster_info !% table of cluster info, modified if necessary on output
    logical, intent(in) :: connectivity_just_from_connect !% if true, we're doing hysterestic connect and should rely on the connection object completely
    type(Connection), intent(in) :: use_connect !% connection object to use for connectivity info
    logical, intent(in) :: atom_mask(6) !% which fields in int part of table to compare when checking for identical atoms
    logical :: cluster_changed

    integer :: n, i, ishift(3), m, j, jshift(3), p, k, kshift(3)
    integer :: n_bonds_in, n_bonds_out

    call print('doing cluster_reduce_n_cut_bonds', PRINT_NERD)

    cluster_changed = .false.

    ! loop over each atom in the cluster already
    n = 0
    do while (n < cluster_info%N)
      n = n + 1
      i = cluster_info%int(1,n)
      ishift = cluster_info%int(2:4,n)
      call print('cluster_reduce_n_cut_bonds: i = '//i//'. Looping over '//n_neighbours(this,i,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)

      ! loop over every neighbour, looking for outside neighbours
      do m=1, n_neighbours(this, i, alt_connect=use_connect)
	j = neighbour(this, i, m, shift=jshift, alt_connect=use_connect)
	if (find(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), atom_mask ) /= 0) cycle
	if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, m, alt_connect=use_connect))) then
	  call print("cluster_reduce_n_cut_bonds:   j = "//j//" ["//jshift//"] not nearest neighbour",PRINT_ANAL)
	  cycle
	end if
	! if we're here, j must be an outside neighbour of i
	call print('cluster_reduce_n_cut_bonds: j = '//j//'. Looping over '//n_neighbours(this,j,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)

	n_bonds_in = 0
	n_bonds_out = 0
	do p=1, n_neighbours(this, j, alt_connect=use_connect)
	  k = neighbour(this, j, p, shift=kshift, alt_connect=use_connect)
	  if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,j, p, alt_connect=use_connect))) then
	    call print("cluster_reduce_n_cut_bonds:   k = "//k//" ["//kshift//"] not nearest neighbour of j = " // j,PRINT_ANAL)
	    cycle
	  end if
	  ! count how many bonds point in vs. out
	  if (find(cluster_info, (/ k, ishift+jshift+kshift, this%Z(k), 0 /), atom_mask ) /= 0) then
	    n_bonds_in = n_bonds_in + 1
	  else
	    n_bonds_out = n_bonds_out + 1
	  endif
	end do ! p

	if (n_bonds_out < n_bonds_in) then ! adding this one would reduce number of cut bonds
	  call append(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), (/ this%pos(:,j), 1.0_dp /), (/"n_cut_bond" /) )
	  cluster_changed = .true.
	  call print('cluster_reduce_n_cut_bonds:  Added atom ' //j//' ['//(ishift+jshift)//'] n_bonds_in ' // n_bonds_in // &
		     ' out ' // n_bonds_out // '  to cluster. Atoms = ' // cluster_info%N, PRINT_NERD)
	endif
      end do ! m

    end do ! while (n <= cluster_info%N)

    call print('cluster_reduce_n_cut_bonds: Finished checking',PRINT_NERD)
    call print("cluster_reduce_n_cut_bonds: cluster list:", PRINT_NERD)
    call print(cluster_info, PRINT_NERD)
  end function cluster_reduce_n_cut_bonds

  !% Go through the cluster atoms and find cut bonds. If the bond is to
  !% hydrogen, then include the hydrogen in the cluster list (since AMBER
  !% doesn't cap a bond to hydrogen with a link atom)
  function cluster_protect_X_H_bonds(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask) result(cluster_changed)
    type(Atoms), intent(in) :: this !% atoms structure 
    type(Table), intent(inout) :: cluster_info !% table of cluster info, modified if necessary on output
    logical, intent(in) :: connectivity_just_from_connect !% if true, we're doing hysterestic connect and should rely on the connection object completely
    type(Connection), intent(in) :: use_connect !% connection object to use for connectivity info
    logical, intent(in) :: atom_mask(6) !% which fields in int part of table to compare when checking for identical atoms
    logical :: cluster_changed

    integer :: n, i, ishift(3), m, j, jshift(3)

    call print('doing cluster_protect_X_H_bonds', PRINT_NERD)
    cluster_changed = .false.

    ! loop over each atom in the cluster already
    n = 0
    do while (n < cluster_info%N)
      n = n + 1
      i = cluster_info%int(1,n)
      ishift = cluster_info%int(2:4,n)
      call print('cluster_protect_X_H_bonds: i = '//i//'. Looping over '//n_neighbours(this,i,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)

      ! loop over every neighbour, looking for outside neighbours
      do m=1, n_neighbours(this, i, alt_connect=use_connect)
	j = neighbour(this, i, m, shift=jshift, alt_connect=use_connect)

	! if j is in, or not a nearest neighbour, go on to next
	if (find(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), atom_mask ) /= 0) cycle
	if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, m, alt_connect=use_connect))) then
	  call print("cluster_protect_X_H_bonds:   j = "//j//" ["//jshift//"] not nearest neighbour",PRINT_ANAL)
	  cycle
	end if
	! if we're here, j must be an outside neighbour of i

	if (this%Z(i) == 1 .or. this%Z(j) == 1) then
	  call append(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), (/ this%pos(:,j), 1.0_dp /), (/"X_H_bond  " /) )
	  cluster_changed = .true.
	  call print('cluster_protect_X_H_bonds:  Added atom ' //j//' ['//(ishift+jshift)//'] to cluster. Atoms = ' // cluster_info%N, PRINT_NERD)
	endif
      end do ! m

    end do ! while (n <= cluster_info%N)

    call print('cluster_protect_X_H_bonds: Finished checking',PRINT_NERD)
    call print("cluster_protect_X_H_bonds: cluster list:", PRINT_NERD)
    call print(cluster_info, PRINT_NERD)
  end function cluster_protect_X_H_bonds

  !% Go through the cluster atoms and find cut bonds. 
  !% Include the bonded-to atom in the cluster list.
  !% cheap alternative to keep_residues_whole
  function cluster_keep_whole_molecules(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask) result(cluster_changed)
    type(Atoms), intent(in) :: this !% atoms structure 
    type(Table), intent(inout) :: cluster_info !% table of cluster info, modified if necessary on output
    logical, intent(in) :: connectivity_just_from_connect !% if true, we're doing hysterestic connect and should rely on the connection object completely
    type(Connection), intent(in) :: use_connect !% connection object to use for connectivity info
    logical, intent(in) :: atom_mask(6) !% which fields in int part of table to compare when checking for identical atoms
    logical :: cluster_changed

    integer :: n, i, ishift(3), m, j, jshift(3)

    call print('doing cluster_keep_whole_molecules', PRINT_NERD)
    cluster_changed = .false.

    ! loop over each atom in the cluster already
    n = 0
    do while (n < cluster_info%N)
      n = n + 1
      i = cluster_info%int(1,n)
      ishift = cluster_info%int(2:4,n)
      call print('cluster_keep_whole_molecules: i = '//i//'. Looping over '//n_neighbours(this,i,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)

      ! loop over every neighbour, looking for outside neighbours
      do m=1, n_neighbours(this, i, alt_connect=use_connect)
	j = neighbour(this, i, m, shift=jshift, alt_connect=use_connect)

	! if j is in, or not a nearest neighbour, go on to next
	if (find(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), atom_mask ) /= 0) cycle
	if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, m, alt_connect=use_connect))) then
	  call print("cluster_keep_whole_molecules:   j = "//j//" ["//jshift//"] not nearest neighbour",PRINT_ANAL)
	  cycle
	end if
	! if we're here, j must be an outside neighbour of i

	call append(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), (/ this%pos(:,j), 1.0_dp /), (/"molec     " /) )
	cluster_changed = .true.
	call print('cluster_keep_whole_molecules:  Added atom ' //j//' ['//(ishift+jshift)//'] to cluster. Atoms = ' // cluster_info%N, PRINT_NERD)
      end do ! m

    end do ! while (n <= cluster_info%N)

    call print('cluster_keep_whole_molecules: Finished checking',PRINT_NERD)
    call print("cluster_keep_whole_molecules: cluster list:", PRINT_NERD)
    call print(cluster_info, PRINT_NERD)
  end function cluster_keep_whole_molecules

  !if by accident a single atom is included, that is part of a larger entity (not ion)
  !then delete it from the list as it will cause error in SCF (e.g. C=O oxygen)
  !call allocate(to_delete,4,0,0,0,1)
  !do n = 1, cluster%N
  !        i = cluster%int(1,n)
  !        is_alone = .true.
  !        do j = 1, n_neighbours(at,i)
  !                k = neighbour(at,i,j)
  !                if(find(cluster_info,(/k,shift/)) /= 0) then
  !                   if(is_nearest_neighbour(at,i,j)) GOTO 120
  !                else
  !                   if(is_nearest_neighbour(at,i,j)) is_alone = .false.
  !                end if
  !        end do
  !        if(.not. is_alone) call append(to_delete,(/i,shift/))
  !   120  cycle
  !end do

  !do i = 1, to_delete%N
  ! call delete(cluster_info,to_delete%int(:,i))
  !end do

!OUTDATED   function cluster_biochem_in_out_in(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask) result(cluster_changed)
!OUTDATED     type(Atoms), intent(in) :: this
!OUTDATED     type(Table), intent(inout) :: cluster_info
!OUTDATED     logical, intent(in) :: connectivity_just_from_connect
!OUTDATED     type(Connection), intent(in) :: use_connect
!OUTDATED     logical, intent(in) :: atom_mask(6)
!OUTDATED     logical :: cluster_changed
!OUTDATED 
!OUTDATED     integer :: n, i, ishift(3), m, j, jshift(3), p, k, kshift(3)
!OUTDATED 
!OUTDATED     call print('doing cluster_biochem_in_out_in', PRINT_NERD)
!OUTDATED     cluster_changed = .false.
!OUTDATED 
!OUTDATED     ! loop over each atom in the cluster already
!OUTDATED     n = 1
!OUTDATED     do while (n <= cluster_info%N)
!OUTDATED       i = cluster_info%int(1,n)
!OUTDATED       ishift = cluster_info%int(2:4,n)
!OUTDATED       call print('cluster_reduce_n_cut_bonds: i = '//i//'. Looping over '//n_neighbours(this,i,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)
!OUTDATED 
!OUTDATED       ! loop over every neighbour, looking for outside neighbours
!OUTDATED       do m=1, n_neighbours(this, i, alt_connect=use_connect)
!OUTDATED 	j = neighbour(this, i, m, shift=jshift, alt_connect=use_connect)
!OUTDATED 
!OUTDATED 	! if j is in, or not a neareste neighbour, go on to next
!OUTDATED 	if (find(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), atom_mask ) /= 0) cycle
!OUTDATED 	if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, m, alt_connect=use_connect))) then
!OUTDATED 	  call print("cluster_reduce_n_cut_bonds:   j = "//j//" ["//jshift//"] not nearest neighbour",PRINT_ANAL)
!OUTDATED 	  cycle
!OUTDATED 	end if
!OUTDATED 	! if we're here, j must be an outside neighbour of i
!OUTDATED 	call print('cluster_reduce_n_cut_bonds: j = '//j//'. Looping over '//n_neighbours(this,j,alt_connect=use_connect)//' neighbours...',PRINT_ANAL)
!OUTDATED 
!OUTDATED 	do p=1, n_neighbours(this, j, alt_connect=use_connect)
!OUTDATED 	  k = neighbour(this, j, p, shift=kshift, alt_connect=use_connect)
!OUTDATED 	  if (k == i .and. all(jshift + kshift == 0)) cycle
!OUTDATED 	  if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,j, p, alt_connect=use_connect))) then
!OUTDATED 	    call print("cluster_reduce_n_cut_bonds:   k = "//k//" ["//kshift//"] not nearest neighbour of j = " // j,PRINT_ANAL)
!OUTDATED 	    cycle
!OUTDATED 	  end if
!OUTDATED 	  ! if k is in, then j has 2 in neighbours, and so we add it
!OUTDATED 	  if (find(cluster_info, (/ k, ishift+jshift+kshift, this%Z(k), 0 /), atom_mask ) /= 0) then
!OUTDATED 	    call append(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), (/ this%pos(:,j), 1.0_dp /), (/"bio_IOI   " /) )
!OUTDATED 	    cluster_changed = .true.
!OUTDATED 	    call print('cluster_biochem_in_out_in:  Added atom ' //j//' ['//(ishift+jshift)//'] to cluster. Atoms = ' // cluster_info%N, PRINT_NERD)
!OUTDATED 	    exit
!OUTDATED 	  end if 
!OUTDATED 	end do ! p
!OUTDATED       end do ! m
!OUTDATED 
!OUTDATED       n = n + 1
!OUTDATED     end do ! while (n <= cluster_info%N)
!OUTDATED 
!OUTDATED     call print('cluster_biochem_in_out_in: Finished checking',PRINT_NERD)
!OUTDATED     call print("cluster_biochem_in_out_in: cluster list:", PRINT_NERD)
!OUTDATED     call print(cluster_info, PRINT_NERD)
!OUTDATED   end function cluster_biochem_in_out_in

  !% if an in non-H atom is in a residue and not at right coordination, add all of its neighbours
  function cluster_protect_double_bonds(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, present_protect_double_bonds) result(cluster_changed)
    type(Atoms), intent(in) :: this !% atoms structure 
    type(Table), intent(inout) :: cluster_info !% table of cluster info, modified if necessary on output
    logical, intent(in) :: connectivity_just_from_connect !% if true, we're doing hysterestic connect and should rely on the connection object completely
    type(Connection), intent(in) :: use_connect !% connection object to use for connectivity info
    logical, intent(in) :: atom_mask(6) !% which fields in int part of table to compare when checking for identical atoms
    logical, intent(in) :: present_protect_double_bonds !% if true, protect_double_bonds was specified by user, not just a default value
    logical :: cluster_changed

    integer :: n, i, ishift(3), m, j, jshift(3)
    integer ::  n_nearest_neighbours
    integer, pointer :: atom_res_number(:)

    call print('doing cluster_protect_double_bonds', PRINT_NERD)
    cluster_changed = .false.

    if (.not. assign_pointer(this, 'atom_res_number', atom_res_number)) then
      if (present_protect_double_bonds) then
	call print("WARNING: cluster_protect_double_bonds got protect_double_bonds requested explicitly, but no proper atom_res_number property available", PRINT_ALWAYS)
      endif
      return
    endif

    n = 0
    do while (n < cluster_info%N)
      n = n + 1
      i = cluster_info%int(1,n)
      ishift = cluster_info%int(2:4,n)

      if (this%Z(i) /= 0 .and. atom_res_number(i) >= 0) then

	! count nearest neighbours of i
	n_nearest_neighbours = 0
	do m=1, n_neighbours(this, i, alt_connect=use_connect)
	  j = neighbour(this, i, m, alt_connect=use_connect)

	  ! if j is not nearest neighbour, go on to next one
	  if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, m, alt_connect=use_connect))) then
	    call print("cluster_protect_double_bonds:   j = "//j//" not nearest neighbour",PRINT_ANAL)
	    cycle
	  end if
	  ! if we're here, j must be an outside neighbour of i
	  n_nearest_neighbours = n_nearest_neighbours + 1
	end do ! m

	if (ElementValence(this%Z(i)) /= -1 .and. (ElementValence(this%Z(i)) /= n_nearest_neighbours)) then ! i has a valence, and it doesn't match nn#
	  do m=1, n_neighbours(this, i)
	    j = neighbour(this, i, m, shift=jshift, alt_connect=use_connect)

	    ! if j is in, or isn't nearest neighbour, go on to next
	    if (find(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), atom_mask ) /= 0) cycle
	    if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, m, alt_connect=use_connect))) then
	      call print("cluster_protect_double_bonds:   j = "//j//" ["//jshift//"] not nearest neighbour",PRINT_ANAL)
	      cycle
	    end if

	    ! if we're here, j must be an outside neighbour of i
	    call append(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), (/ this%pos(:,j), 1.0_dp /), (/"dbl_bond  " /) )
	    cluster_changed = .true.
	    call print('cluster_protect_double_bonds:  Added atom ' //j//' ['//(ishift+jshift)//'] to cluster. Atoms = ' // cluster_info%N, PRINT_NERD)
	  end do ! m
	endif ! valence defined
      end if !  atom_res_number(i) >= 0

    end do ! while (n <= cluster_info%N)

    call print('cluster_protect_double_bonds: Finished checking',PRINT_NERD)
    call print("cluster_protect_double_bonds: cluster list:", PRINT_NERD)
    call print(cluster_info, PRINT_NERD)
  end function cluster_protect_double_bonds

  !% protect bond between R-C ... -N-R
  !%                        |      |
  !%                        O      R (usually H or C (proline))
  function cluster_protect_peptide_bonds(this, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, present_protect_peptide_bonds) result(cluster_changed)
    type(Atoms), intent(in) :: this !% atoms structure 
    type(Table), intent(inout) :: cluster_info !% table of cluster info, modified if necessary on output
    logical, intent(in) :: connectivity_just_from_connect !% if true, we're doing hysterestic connect and should rely on the connection object completely
    type(Connection), intent(in) :: use_connect !% connection object to use for connectivity info
    logical, intent(in) :: atom_mask(6) !% which fields in int part of table to compare when checking for identical atoms
    logical, intent(in) :: present_protect_peptide_bonds !% if true, protect_peptide_bonds was specified by user, not just a default value
    logical :: cluster_changed

    integer :: n, i, ishift(3), ji, j, ki, k, jshift(3)
    integer, pointer :: atom_res_number(:)
    integer :: found_other_j, found_other_jshift(3)
    integer :: this_neighbour_Z, other_neighbour_Z, other_Z
    logical :: found_this_neighbour, found_other_neighbour

    call print('doing cluster_protect_peptide_bonds', PRINT_NERD)
    cluster_changed = .false.

    if (.not. assign_pointer(this, 'atom_res_number', atom_res_number)) then
      if (present_protect_peptide_bonds) then
	call print("WARNING: cluster_protect_peptide_bonds got protect_peptide_bonds requested explicitly, but no proper atom_res_number property available", PRINT_ALWAYS)
      endif
      return
    endif

    n = 0
    do while (n < cluster_info%N)
      n = n + 1
      i = cluster_info%int(1,n)
      ishift = cluster_info%int(2:4,n)

      if (this%Z(i) /= 0 .and. atom_res_number(i) >= 0) then

	call print("cluster_protect_peptide_bonds: i = "//i//" Z(i) " // this%Z(i), PRINT_ANAL)

	! peptide bond is C-N
	!                 |
	!                 O
	if (this%Z(i) == 6) then
	  other_Z = 7
	  this_neighbour_Z = 8
	  other_neighbour_Z = 0
	else if (this%Z(i) == 7) then
	  other_Z = 6
	  this_neighbour_Z = 0
	  other_neighbour_Z = 8
	else
	  call print("cluster_protect_peptide_bonds:   i = "//i//" not C or N",PRINT_ANAL)
	  cycle
	endif

	found_other_j = 0
	found_this_neighbour = (this_neighbour_Z == 0)
	found_other_neighbour = (other_neighbour_Z == 0)

	call print("cluster_protect_peptide_bonds: other_Z " // other_Z // " this_neighbour_Z " // this_neighbour_Z // " other_neighbour_z " // other_neighbour_Z, PRINT_ANAL)
	call print("cluster_protect_peptide_bonds: found_this_neighbour " // found_this_neighbour // " found_other_neighbour " // found_other_neighbour, PRINT_ANAL)

	! must have 3 neighbours
	if (n_neighbours(this, i, alt_connect=use_connect) /= 3) then
	  call print("cluster_protect_peptide_bonds:   i = "//i//" doesn't have 3 neighbours",PRINT_ANAL)
	  cycle
	endif
	do ji=1, n_neighbours(this, i, alt_connect=use_connect) ! look for correct neighbour
	  j=neighbour(this, i, ji, shift=jshift, alt_connect=use_connect)
	  call print("cluster_protect_peptide_bonds:   j = " //j//" Z(j) " // this%Z(j), PRINT_ANAL)
	  ! only nearest neighbours count, usually
	  if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,i, ji, alt_connect=use_connect))) then
	    call print("cluster_protect_peptide_bonds:   neighbour j = "//j//" not nearest neighbour",PRINT_ANAL)
	    cycle
	  end if
	  if (this%Z(j) == this_neighbour_Z .and. n_neighbours(this, j) == 1) then ! neighbour has right Z and coordination
	    found_this_neighbour = .true.
	    call print("cluster_protect_peptide_bonds:   neighbour j = "//j//" is neighbor (O) we're looking for",PRINT_ANAL)
	  endif
	  if (this%Z(j) == other_Z) then ! found other
	    ! if j is already in, this isn't a peptide bond in need of protection from being cut
	    if (find(cluster_info, (/ j, ishift+jshift, this%Z(j), 0 /), atom_mask) /= 0) then
	      call print("cluster_protect_peptide_bonds:   other_Z neighbour j = "//j//" already in",PRINT_ANAL)
	      cycle
	    endif
	    ! save info on other
	    call print("cluster_protect_peptide_bonds: found_other_j " // j // " Z(found_other_j) " // this%Z(j), PRINT_ANAL)
	    found_other_j = j
	    found_other_jshift = jshift

	    ! if (.not. found_this_neighbour .or. found_other_j == 0) then ! some part of motif missing
	      ! call print("cluster_protect_peptide_bonds:   failed to find this_neighbour_Z or other_j",PRINT_ANAL)
	      ! cycle
	    ! endif

	    ! other must have 3 neighbours
	    if (n_neighbours(this, found_other_j, alt_connect=use_connect) /= 3) then
	      call print("cluster_protect_peptide_bonds:   found_other_j = "//found_other_j//" doesn't have 3 neighbours", PRINT_ANAL)
	      cycle
	    endif
	    ! perhaps we should require that other be in a different residue, but not doing that now
	    ! if (atoms_res_number(i) == atom_res_number(found_other_j)) then
	      ! call print("cluster_protect_peptide_bonds:   i and found_other_j are in the same residue", PRINT_ANAL)
	      ! cycle
	    ! endif

	    ! check neighbours of other
	    do ki=1, n_neighbours(this, found_other_j, alt_connect=use_connect)
	      k=neighbour(this, found_other_j, ki, alt_connect=use_connect)
	      ! only nearest neighbours count
	      if(.not. (connectivity_just_from_connect .or. is_nearest_neighbour(this,found_other_j, ki, alt_connect=use_connect))) then
		call print("cluster_protect_peptide_bonds:   found_other_j's neighbour k = "//k//" not nearest neighbour",PRINT_ANAL)
		cycle
	      end if

	      if (this%Z(k) == other_neighbour_Z .and. n_neighbours(this, k) == 1) then ! neighbour has right Z and coordination
		found_other_neighbour = .true.
		call print("cluster_protect_peptide_bonds:   neighbour k = "//k//" is neighbor (O) we're looking for",PRINT_ANAL)
	      endif
	    end do ! ki

	    ! some part of motif is missing
	    if (.not. found_other_neighbour) then
	      call print("cluster_protect_peptide_bonds:   couldn't find other_neighbour",PRINT_ANAL)
	      cycle
	    endif

	    call append(cluster_info, (/ found_other_j, ishift+found_other_jshift, this%Z(found_other_j), 0 /), &
			(/ this%pos(:,found_other_j), 1.0_dp /), (/"pep_bond  " /) )
	    cluster_changed = .true.
	    call print('cluster_protect_peptide_bonds:  Added atom ' //found_other_j//' ['//(ishift+found_other_jshift)//'] to cluster. Atoms = ' // cluster_info%N, PRINT_NERD)
	  endif ! Z(j) == other_Z
	end do ! ji

      end if ! Z(i) /= 0 .and. atom_res_number(i) /= 0

    end do ! while (n <= cluster_info%N)

    call print('cluster_protect_peptide_bonds: Finished checking',PRINT_NERD)
    call print("cluster_protect_peptide_bonds: cluster list:", PRINT_NERD)
    call print(cluster_info, PRINT_NERD)
  end function cluster_protect_peptide_bonds

  subroutine create_cluster_simple(at, args_str, cluster, mark_name, error)
     type(Atoms), intent(inout) :: at
     character(len=*), intent(in) :: args_str
     type(Atoms), intent(out) :: cluster
     character(len=*), intent(in), optional :: mark_name
     integer, intent(out), optional :: ERROR

     type(Table) :: cluster_info

     INIT_ERROR(error)

     call calc_connect(at)
     cluster_info = create_cluster_info_from_mark(at, args_str, mark_name=mark_name, error=error)
     PASS_ERROR(error)
     call carve_cluster(at, args_str, cluster_info, cluster=cluster, mark_name=mark_name, error=error)
     PASS_ERROR(error)
     call finalise(cluster_info)
     
  end subroutine create_cluster_simple

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Create cluster from atoms and cluster information table (which
  !% which should be generated by a call to :func:`create_cluster_info_from_hybrid_mark`.
  !%
  !% The output cluster contains all properties of the initial atoms object, and
  !% some additional columns, which are:
  !%
  !%  ``index``
  !%      index of the cluster atoms into the initial atoms object.
  !%
  !%  ``termindex`` 
  !%     nonzero for termination atoms, and is an index into
  !%     the cluster atoms specifiying which atom is being terminated,
  !%     it is used in collecting the forces.
  !%
  !%  ``rescale`` 
  !%      a real number which for nontermination atoms is 1.0,
  !%      for termination atoms records the scaling applied to
  !%      termination bond lengths
  !%  
  !%  ``shift``
  !%      the shift of each atom
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  subroutine carve_cluster(at, args_str, cluster_info, cluster, mark_name, error)
    type(Atoms), intent(in), target :: at
    character(len=*), intent(in) :: args_str
    type(Table), intent(in) :: cluster_info
    type(Atoms), intent(out) :: cluster
    character(len=*), optional, intent(in) :: mark_name
    integer, optional, intent(out) :: error

    type(Dictionary) :: params
    logical :: do_rescale_r
    real(dp) :: r_scale, cluster_vacuum, cluster_fix_lattice(3), cluster_round_lattice
    logical :: terminate, randomise_buffer, print_clusters, do_calc_connect, do_same_lattice, do_fix_lattice
    logical :: hysteretic_connect
    integer :: i, j, k, m
    real(dp) :: maxlen(3), sep(3), lat_maxlen(3), lat_sep(3)
    logical :: do_periodic(3)
    integer, pointer :: hybrid_mark(:), cluster_index(:), cluster_hybrid_mark(:), cluster_shift(:,:)
    type(CInoutput)                    :: clusterfile
    character(len=255)                :: clusterfilename
    character(len=255)                :: my_mark_name
    character(STRING_LENGTH)          :: run_suffix
    type(Table) :: outer_layer
    logical :: in_outer_layer, do_map_into_cell

#ifdef _MPI
    integer::mpi_size, mpi_rank, PRINT_ALWAYS
    include "mpif.h"
    integer :: mpi_force_size
    real(dp), allocatable, dimension(:)  :: mpi_force

    call get_mpi_size_rank(MPI_COMM_WORLD, mpi_size, mpi_rank)
#endif /* _MPI */

    INIT_ERROR(error)

    call print('carve_cluster got args_str "'//trim(args_str)//'"', PRINT_VERBOSE)

    call initialise(params)
    call param_register(params, 'terminate', 'T', terminate,&
      help_string="do hydrogen termination")
    call param_register(params, 'do_rescale_r', 'F', do_rescale_r,&
      help_string="do rescale cluster lattice positions and lattice (for matching lattice constants of different models")
    call param_register(params, 'r_scale', '1.0', r_scale,&
      help_string="cluster position/lattice rescaling factor")
    call param_register(params, 'randomise_buffer', 'T', randomise_buffer,&
      help_string="randomly perturb positions of buffer atoms")
    call param_register(params, 'print_clusters', 'F', print_clusters,&
      help_string="Print out clusters into clusters.xyz file (clusters.MPI_RANK.xyz if in MPI)")
    call param_register(params, 'cluster_calc_connect', 'F', do_calc_connect,&
      help_string="run calc_connect before doing stuff (passivation, etc)")
    call param_register(params, 'cluster_same_lattice', 'F', do_same_lattice,&
      help_string="Don't modify cluster cell vectors from incoming cell vectors")
    call param_register(params, 'cluster_fix_lattice', '0. 0. 0.', cluster_fix_lattice, &
      has_value_target=do_fix_lattice, &
      help_string="fix the cluster lattice to be cubic with given 3 lengths. useful with DFT codes where one doesn't pay for vacuum.")
    call param_register(params, 'cluster_round_lattice', '3.0', cluster_round_lattice, &
      help_string="distance to which to round cluster lattice (default is 3.0 A)")
    call param_register(params, 'cluster_periodic_x', 'F', do_periodic(1),&
      help_string="do not add vaccum along x, so as to make cluster's cell vectors periodic")
    call param_register(params, 'cluster_periodic_y', 'F', do_periodic(2),&
      help_string="do not add vaccum along x, so as to make cluster's cell vectors periodic")
    call param_register(params, 'cluster_periodic_z', 'F', do_periodic(3),&
      help_string="do not add vaccum along x, so as to make cluster's cell vectors periodic")
    call param_register(params, 'cluster_vacuum', '10.0', cluster_vacuum,&
      help_string="Amount of vacuum to add around cluster (minimum - lattice vectors are rounded up to integers)")
    call param_register(params, 'hysteretic_connect', 'F', hysteretic_connect,&
      help_string="Use hysteretic connection object for identifying bonds")
    call param_register(params, 'map_into_cell', 'T', do_map_into_cell,&
      help_string="Call map_into_cell (lattice coordinates=[-0.5,0.5)) on cluster")
    call param_register(params, 'run_suffix', '', run_suffix, &
         help_string="string to append to names of extra properties added to cluster. NOTE: does not override mark_name argument")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='carve_cluster arg_str') ) then
      RAISE_ERROR("carve_cluster failed to parse args_str='"//trim(args_str)//"'", error)
    endif
    call finalise(params)

    ! first pick up an atoms structure with the right number of atoms and copy the properties
    !Now turn the cluster_temp table into an atoms structure
    call print('carve_cluster: Copying atomic data to output object',PRINT_NERD)
    call print('List of atoms in cluster:', PRINT_NERD)
    call print(int_part(cluster_info,1), PRINT_NERD)

    call select(cluster, at, list=int_part(cluster_info,1), error=error)
    PASS_ERROR(error)
    ! then reset the positions species and Z (latter two needed because termination atoms have Z=1)
    ! unfold the positions to real positions using the stored shifts, at is neede because
    ! next we might make the unit cell much smaller
    do i=1,cluster_info%N
       cluster%pos(:,i) = cluster_info%real(1:3, i)+(at%lattice .mult. cluster_info%int(2:4, i))
       cluster%Z(i) = cluster_info%int(5,i)
       cluster%species(:,i) = ' '
       cluster%species(1:3,i) = s2a(ElementName(cluster_info%int(5,i)))
    end do
    ! add properties to cluster
    call add_property(cluster, 'index'//trim(run_suffix), int_part(cluster_info,1))
    call add_property(cluster, 'shift'//trim(run_suffix), 0, n_cols=3, ptr2=cluster_shift)
    cluster_shift(:,1:cluster%N) = cluster_info%int(2:4,1:cluster_info%N)
    call add_property(cluster, 'termindex'//trim(run_suffix), int_part(cluster_info,6))
    call add_property(cluster, 'rescale'//trim(run_suffix), real_part(cluster_info,4))
    call add_property(cluster, 'cluster_ident'//trim(run_suffix), cluster_info%str(1,1:cluster_info%N))

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

    ! Round up maxlen to be divisible by cluster_round_lattice (default 3 A), so that it does not fluctuate too much
    forall (k=1:3) maxlen(k) = cluster_round_lattice*ceiling(maxlen(k)/cluster_round_lattice)
    forall (k=1:3) lat_maxlen(k) = cluster_round_lattice*ceiling(lat_maxlen(k)/cluster_round_lattice)

    ! vacuum pad cluster (if didn't set do_same_lattice or cluster_fix_lattice)
    ! if not periodic at all, just do vacuum padding
    ! if cluster_fix_lattice was given, just set the diagonal lattice elements to those values
    ! if periodic along some dir, keep supercell vector directions, and set
    !    extent in each direction to lesser of cluster extent + vacuum or original extent
    if (do_same_lattice) then
      cluster%lattice = at%lattice
    else if (do_fix_lattice) then
      cluster%lattice = 0.0_dp
      do k=1,3
         cluster%lattice(k,k) = cluster_fix_lattice(k)
      end do
    else
      if (any(do_periodic)) then
	do k=1,3
	  if (do_periodic(k)) then
	    if (lat_maxlen(k)+cluster_vacuum >= norm(at%lattice(:,k))) then
	      cluster%lattice(:,k) = at%lattice(:,k)
	    else
	      cluster%lattice(:,k) = (lat_maxlen(k)+cluster_vacuum)*at%lattice(:,k)/norm(at%lattice(:,k))
	    endif
	  else
	    cluster%lattice(:,k) = (lat_maxlen(k)+cluster_vacuum)*at%lattice(:,k)/norm(at%lattice(:,k))
	  endif
	end do
      else
	cluster%lattice = 0.0_dp
	do k=1,3
	  cluster%lattice(k,k) = maxlen(k) + cluster_vacuum
	end do
      endif
    endif

    call set_lattice(cluster, cluster%lattice, scale_positions=.false.)

    ! Remap positions so any image atoms end up inside the cell
    if (do_map_into_cell) then
       call map_into_cell(cluster)
    end if

    call print ('carve_cluster: carved cluster with '//cluster%N//' atoms', PRINT_VERBOSE)

    write (line, '(a,f10.3,f10.3,f10.3)') &
         'carve_cluster: Cluster dimensions are ', cluster%lattice(1,1), &
         cluster%lattice(2,2), cluster%lattice(3,3)
    call print(line, PRINT_VERBOSE)

    my_mark_name = optional_default('hybrid_mark', mark_name)
    ! reassign pointers
    if (.not. assign_pointer(at, trim(my_mark_name), hybrid_mark)) then
         RAISE_ERROR('cannot reassign "'//trim(my_mark_name)//'" property', error)
    endif

    ! rescale cluster positions and lattice 
    if (do_rescale_r) then
       call print('carve_cluster: rescaling cluster positions (and lattice if not fixed) by factor '//r_scale, PRINT_VERBOSE)
       cluster%pos = r_scale * cluster%pos
       if (.not. do_fix_lattice) then
          call set_lattice(cluster, r_scale * cluster%lattice, scale_positions=.false.)
       end if
    end if

    if (randomise_buffer .and. .not. any(hybrid_mark == HYBRID_BUFFER_OUTER_LAYER_MARK)) &
         do_calc_connect = .true.

    if (do_calc_connect) then
       call print('carve_cluster: doing calc_connect', PRINT_VERBOSE)
       ! Does QM force model need connectivity information?
       if (at%use_uniform_cutoff) then
          call set_cutoff(cluster, at%cutoff)
       else
          call set_cutoff_factor(cluster, at%cutoff)
       end if
       call calc_connect(cluster)
    end if

    if (randomise_buffer) then
       ! If outer buffer layer is not marked do so now. This will be
       ! the case if we're using hysteretic_buffer selection option.
       ! In this case we consider any atoms connected to terminating
       ! hydrogens to be in the outer layer - obviously this breaks down 
       ! if we're not terminating so we abort if that's the case.
       if (.not. assign_pointer(cluster, 'index', cluster_index)) then
            RAISE_ERROR('carve_cluster: cluster is missing index property', error)
       endif

       if (.not. any(hybrid_mark == HYBRID_BUFFER_OUTER_LAYER_MARK)) then
          if (.not. terminate) then
	    RAISE_ERROR('cannot determine which buffer atoms to randomise if terminate=F and hysteretic_buffer=T', error)
	  endif

          if (.not. assign_pointer(cluster, trim(my_mark_name), cluster_hybrid_mark)) then
               RAISE_ERROR(trim(my_mark_name)//' property not found in cluster', error)
	  endif

          do i=1,cluster%N
             if (hybrid_mark(cluster_info%int(1,i)) /= HYBRID_BUFFER_MARK) cycle
             in_outer_layer = .false.
             do m = 1, cluster_info%N
                if (cluster_info%int(6,m).eq. i) then
                   in_outer_layer = .true.
                   exit
                end if
             enddo
             if (in_outer_layer) then
                hybrid_mark(cluster_info%int(1,i)) = HYBRID_BUFFER_OUTER_LAYER_MARK
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

    if (value(mainlog%verbosity_stack) >= PRINT_VERBOSE .or. print_clusters) then
#ifdef _MPI
       write (clusterfilename, '(a,i3.3,a)') 'clusters.',mpi_rank,'.xyz'
#else
       clusterfilename = 'clusters.xyz'
#endif /* _MPI */
       call initialise(clusterfile, clusterfilename, append=.true., action=OUTPUT)
       call write(clusterfile, cluster)
       call finalise(clusterfile)
    end if

 end subroutine carve_cluster

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Create a cluster using the mark_name (optional arg, default 'hybrid_mark') 
  !% property and options in 'args_str'.
  !% All atoms that are marked with anything other than 'HYBRID_NO_MARK' will
  !% be included in the cluster; this includes active, transition and buffer
  !% atoms.   Atoms added (to satisfy heuristics) will be marked in new property
  !% modified_//trim(mark_name) as 'HYBRID_BUFFER_MARK'.
  !% Returns a Table object (cluster_info) which contains info on atoms whose
  !% indices are given in atomlist, possibly with some extras for consistency,
  !% and optionally terminated with Hydrogens, that can be used by carve_cluster().
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function create_cluster_info_from_mark(at, args_str, cut_bonds, mark_name, error) result(cluster_info)
    type(Atoms), intent(inout), target :: at
    character(len=*), intent(in) :: args_str
    type(Table), optional, intent(out)   :: cut_bonds !% Return a list of the bonds cut when making
    !% the cluster.  See create_cluster() documentation.
    character(len=*), intent(in), optional :: mark_name
    integer, optional, intent(out) :: error
    type(Table) :: cluster_info

    type(Dictionary) :: params
    logical :: terminate, periodic_x, periodic_y, periodic_z, &
         even_electrons, do_periodic(3), cluster_hopping_nneighb_only, cluster_heuristics_nneighb_only, &
         cluster_allow_modification, hysteretic_connect, same_lattice, &
         fix_termination_clash, keep_whole_residues, keep_whole_subgroups, keep_whole_prolines, keep_whole_proline_sidechains, &
	 keep_whole_silica_tetrahedra, keep_whole_titania_octahedra, terminate_octahedra,reduce_n_cut_bonds, in_out_in, &
         protect_X_H_bonds, protect_double_bonds, protect_peptide_bonds, keep_whole_molecules, has_termination_rescale, &
	 combined_protein_heuristics, force_no_fix_termination_clash
    character(STRING_LENGTH) :: in_out_in_mode
    logical :: keep_whole_residues_has_value, keep_whole_subgroups_has_value, keep_whole_prolines_has_value, keep_whole_proline_sidechains_has_value, &
               protect_double_bonds_has_value, protect_peptide_bonds_has_value, keep_whole_molecules_has_value
    real(dp) :: r, r_min, centre(3), termination_rescale, termination_clash_factor
    type(Table) :: cluster_list, currentlist, nextlist, activelist, bufferlist
    integer :: i, j, jj, first_active, old_n, n_cluster
    integer, pointer :: hybrid_mark(:), modified_hybrid_mark(:)
    integer :: prev_cluster_info_n
    integer, allocatable, dimension(:) :: uniqed, tmp_index

    type(Table)                              :: n_term, in_n_term, sorted_n_term
    integer                                  :: m, n, p, p_in, in_n

    real(dp)                                 :: H1(3)
    real(dp)                                 :: dhat_ij(3), d(3)
    real(dp)                                 :: r_ij, rescale
    real(dp)                                 :: orig_length
    real(dp), dimension(3)                   :: bond1, bond2, new_bond 
    integer                                  :: ishift(3), jshift(3), oldN, most_hydrogens
    logical                                  :: atom_mask(6)
    integer, allocatable, dimension(:) :: idx

    type(Connection), pointer :: use_connect
    logical :: connectivity_just_from_connect
    logical :: cluster_changed, cluster_hopping
    character(len=STRING_LENGTH) :: my_mark_name

    INIT_ERROR(error)

    my_mark_name = optional_default('hybrid_mark', mark_name)

    call print('create_cluster_info_from_mark got args_str "'//trim(args_str)//'"', PRINT_VERBOSE)

    call initialise(params)
    call param_register(params, 'terminate', 'T', terminate, &
      help_string="Do hydrogen termination for dangling bonds")
    call param_register(params, 'cluster_periodic_x', 'F', periodic_x,&
      help_string="do not add vaccum along x, so as to make cluster's cell vectors periodic")
    call param_register(params, 'cluster_periodic_y', 'F', periodic_y,&
      help_string="do not add vaccum along y, so as to make cluster's cell vectors periodic")
    call param_register(params, 'cluster_periodic_z', 'F', periodic_z,&
      help_string="do not add vaccum along z, so as to make cluster's cell vectors periodic")
    call param_register(params, 'even_electrons', 'F', even_electrons,&
      help_string="Remove passivation hydrogen to keep total number of electrons even")
    call param_register(params, 'cluster_hopping_nneighb_only', 'T', cluster_hopping_nneighb_only,&
      help_string="Apply is_nearest_neighbor test to bonds when doing bond hops to find cluster.  If false, uses connectivity object as is")
    call param_register(params, 'cluster_heuristics_nneighb_only', 'T', cluster_heuristics_nneighb_only,&
      help_string="Apply is_nearest_neighbor test to bonds when deciding what's bonded for heuristics.  If false, uses connectivity object as is")
    call param_register(params, 'cluster_allow_modification', 'T', cluster_allow_modification,&
      help_string="Allow cluster to be modified using various heuristics")
    call param_register(params, 'hysteretic_connect', 'F', hysteretic_connect,&
      help_string="Use hysteretic connect object to decide what's bonded.  Usually goes with cluster_*_nneighb_only=F.")
    call param_register(params, 'cluster_same_lattice', 'F', same_lattice,&
      help_string="Don't modify cluster cell vectors from incoming cell vectors")
    call param_register(params, 'fix_termination_clash','T', fix_termination_clash,&
      help_string="Apply termination clash (terimnation H too close together) heureistic")
    call param_register(params, 'force_no_fix_termination_clash','F', force_no_fix_termination_clash,&
      help_string="Overrides values of (terminate .or. termination_clash) to explicitly disable termination clash heuristic")
    call param_register(params, 'combined_protein_heuristics','F', combined_protein_heuristics, &
      help_string="Apply heuristics for proteins: overrides keep_whole_molecules=F, keep_whole_residues=F, keep_whole_subgroups=T, keep_whole_prolines=T, keep_whole_proline_sidechains=F, protect_double_bonds=F, protect_peptide_bonds=T")
    call param_register(params, 'keep_whole_residues','T', keep_whole_residues, has_value_target=keep_whole_residues_has_value,&
      help_string="Apply heuristic keeping identified residues whole")
    call param_register(params, 'keep_whole_prolines','F', keep_whole_prolines, has_value_target=keep_whole_prolines_has_value,&
      help_string="Apply heuristic keeping identified prolines whole")
    call param_register(params, 'keep_whole_proline_sidechains','F', keep_whole_proline_sidechains, has_value_target=keep_whole_proline_sidechains_has_value,&
      help_string="Apply heuristic keeping identified prolines sidechain+directly connected backbone part whole")
    call param_register(params, 'keep_whole_subgroups','F', keep_whole_subgroups, has_value_target=keep_whole_subgroups_has_value,&
      help_string="Apply heuristic keeping identified subgroups within a residue whole")
    call param_register(params, 'keep_whole_silica_tetrahedra','F', keep_whole_silica_tetrahedra,&
      help_string="Apply heureistic keeping SiO2 tetrahedra whole")
    call param_register(params, 'keep_whole_titania_octahedra','F', keep_whole_titania_octahedra,&
      help_string="Apply heureistic keeping TiO2 tetragonal structure whole")
    call param_register(params, 'terminate_octahedra','F', terminate_octahedra,&
      help_string="Apply heureistic terminating TiO2 octahedra")
    call param_register(params, 'reduce_n_cut_bonds','T', reduce_n_cut_bonds,&
      help_string="Apply heuristic that adds atoms if doing so reduces number of broken bonds")
    call param_register(params, 'in_out_in', 'F', in_out_in, &
      help_string="Apply heuristic than adds atoms to avoid IN-OUT-IN configurations. Usually equivalent to reduce_n_cut_bonds.")
    call param_register(params, 'in_out_in_mode', 'all', in_out_in_mode, &
      help_string="Mode to use for in_out_in heuristic - can be either 'all' or 'any'. Default is 'all'.")
    call param_register(params, 'protect_X_H_bonds','T', protect_X_H_bonds,&
      help_string="Apply heuristic protecting X-H bonds - no point H passivating bonds with an H")
    call param_register(params, 'protect_double_bonds','T', protect_double_bonds, has_value_target=protect_double_bonds_has_value,&
      help_string="Apply heuristic protecting double bonds from being cut (based one some heuristic idea of what's a double bond")
    call param_register(params, 'protect_peptide_bonds','F', protect_peptide_bonds, has_value_target=protect_peptide_bonds_has_value,&
      help_string="Apply heuristic protecting peptide bonds from being cut (RO-C...N-RR)")
    call param_register(params, 'keep_whole_molecules','F', keep_whole_molecules, has_value_target=keep_whole_molecules_has_value, &
      help_string="Apply heuristic keeping identified molecules (i.e. contiguous bonded groups of atoms) whole")
    call param_register(params, 'termination_rescale', '0.0', termination_rescale, has_value_target=has_termination_rescale,&
      help_string="rescale factor for X-H passivation bonds")
    call param_register(params, 'termination_clash_factor', '1.2', termination_clash_factor, &
      help_string="Factor to multiply H-H bond-length by when deciding if two Hs are too close. Default 1.2")
    call param_register(params, 'cluster_hopping', 'T', cluster_hopping,&
      help_string="Identify cluster region by hopping from a seed atom.  Needed so that shifts will be correct.")

    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='create_cluster_info_from_mark args_str') ) then
         RAISE_ERROR("create_cluster_info_from_mark failed to parse args_str='"//trim(args_str)//"'", error)
    endif
    call finalise(params)

    if (combined_protein_heuristics) then
      if (keep_whole_molecules_has_value) call print("WARNING: create_cluster_info_from_mark got combined_protein_heuristics, overriding keep_whole_molecules")
      if (keep_whole_residues_has_value) call print("WARNING: create_cluster_info_from_mark got combined_protein_heuristics, overriding keep_whole_residues")
      if (keep_whole_subgroups_has_value) call print("WARNING: create_cluster_info_from_mark got combined_protein_heuristics, overriding keep_whole_subgroups")
      if (keep_whole_prolines_has_value) call print("WARNING: create_cluster_info_from_mark got combined_protein_heuristics, overriding keep_whole_prolines")
      if (keep_whole_proline_sidechains_has_value) call print("WARNING: create_cluster_info_from_mark got combined_protein_heuristics, overriding keep_whole_proline_sidechains")
      if (protect_double_bonds_has_value) call print("WARNING: create_cluster_info_from_mark got combined_protein_heuristics, overriding protect_double_bonds_has_value")
      if (protect_peptide_bonds_has_value) call print("WARNING: create_cluster_info_from_mark got combined_protein_heuristics, overriding protect_peptide_bonds_has_value")
      keep_whole_molecules = .false.
      keep_whole_molecules_has_value = .true.
      keep_whole_residues = .false.
      keep_whole_residues_has_value = .true.
      keep_whole_subgroups = .true.
      keep_whole_subgroups_has_value = .true.
      keep_whole_prolines = .true.
      keep_whole_prolines_has_value = .true.
      keep_whole_proline_sidechains = .true.
      keep_whole_proline_sidechains_has_value = .true.
      protect_double_bonds = .false.
      protect_peptide_bonds = .true.
      protect_peptide_bonds_has_value = .true.
    endif

    do_periodic = (/periodic_x,periodic_y,periodic_z/)

    if (.not. has_property(at, trim(my_mark_name))) then
       RAISE_ERROR('create_cluster_info_from_mark: atoms structure has no "'//trim(my_mark_name)//'" property', error)
    endif

    call assign_property_pointer(at, trim(my_mark_name), hybrid_mark, error=error)
    PASS_ERROR_WITH_INFO('create_cluster_info_from_mark passed atoms structure with no hybrid_mark property', error)

    ! Calculate centre of cluster
    call allocate(cluster_list, 1,0,0,0)
    call append(cluster_list, find(hybrid_mark /= HYBRID_NO_MARK))
    ! unreliable when cluster extends across periodic boundaries
    ! centre = 0.0_dp
    ! do i=1,cluster_list%N
    ! centre = centre + at%pos(:,cluster_list%int(1,i)) + (at%lattice .mult. cluster_list%int(2:4,i))
    ! end do
    ! centre = centre / cluster_list%N
    centre = at%pos(:,cluster_list%int(1,1))
!!$    call print('centre = '//centre)

    ! Find ACTIVE atom closest to centre of cluster, using min image convention  
    r_min = huge(1.0_dp)
    do i=1,cluster_list%N
       if (hybrid_mark(cluster_list%int(1,i)) /= HYBRID_ACTIVE_MARK) cycle
       r = distance_min_image(at, centre, cluster_list%int(1,i))
!!$       call print(cluster_list%int(1,i)//' '//r)
       if (r < r_min) then
          first_active = cluster_list%int(1,i)
          r_min = r
       end if
    end do

    if (cluster_hopping) then
       call system_timer('cluster_hopping')

       n_cluster = cluster_list%N
       call wipe(cluster_list)
       call allocate(cluster_list, 4,0,0,0)

       ! Add first marked atom to cluster_list. shifts will be relative to this atom
       call print('Growing cluster starting from atom '//first_active//', n_cluster='//n_cluster, PRINT_VERBOSE)
       call append(cluster_list, (/first_active,0,0,0/))
       call append(currentlist, cluster_list)

       ! Add other active atoms using bond hopping from the central cluster atom
       ! to find the other cluster atoms and hence to determine the correct 
       ! periodic shifts. 
       !
       ! This will fail if marked atoms do not form a single connected cluster
       old_n = cluster_list%N
       do 
          call wipe(nextlist)
          if (hysteretic_connect) then
             call BFS_step(at, currentlist, nextlist, nneighb_only=cluster_hopping_nneighb_only, &
                  min_images_only = any(do_periodic) .or. same_lattice , alt_connect=at%hysteretic_connect)
          else
             call BFS_step(at, currentlist, nextlist, nneighb_only=cluster_hopping_nneighb_only, &
                  min_images_only = any(do_periodic) .or. same_lattice)
          endif
          do j=1,nextlist%N
             jj = nextlist%int(1,j)
             !          shift = nextlist%int(2:4,j)
             if (hybrid_mark(jj) /= HYBRID_NO_MARK .and. find(cluster_list, nextlist%int(:,j)) == 0) &
                  call append(cluster_list, nextlist%int(:,j))
          end do
          call append(currentlist, nextlist)

          ! check exit condition
          allocate(tmp_index(cluster_list%N))
          tmp_index = int_part(cluster_list,1)
          call sort_array(tmp_index)
          call uniq(tmp_index, uniqed)
          call print('cluster hopping: got '//cluster_list%N//' atoms, of which '//size(uniqed)//' are unique.', PRINT_VERBOSE)
          if (size(uniqed) == n_cluster) exit !got them all
          deallocate(uniqed, tmp_index)

          ! check that cluster is still growing
          if (cluster_list%N == old_n) then
             call write(at, 'create_cluster_abort.xyz')
             call print(cluster_list)
             RAISE_ERROR('create_cluster_info_from_mark: cluster stopped growing before all marked atoms found - check for split QM region', error)
          end if
          old_n = cluster_list%N
       end do
       deallocate(tmp_index, uniqed)
       call finalise(nextlist)
       call finalise(currentlist)
       call system_timer('cluster_hopping')
    else
       call system_timer('cluster_diff_min_image')
       call allocate(cluster_list, 4,0,0,0)

       call append(cluster_list, (/first_active,0,0,0/))
       do i=1, at%n
          if (hybrid_mark(i) == HYBRID_NO_MARK .or. i == first_active) cycle
          ! shifts relative to first_active
          d = distance_min_image(at, first_active, i, shift=ishift)
          call append(cluster_list, (/i,ishift/))
       end do
       call system_timer('cluster_diff_min_image')
    end if
    ! partition cluster_list so that active atoms come first
    call wipe(activelist)
    call wipe(bufferlist)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! 
    ! Validate arguments
    !

    ! check for consistency in optional arguments

    if (.not. (count(do_periodic) == 0 .or. count(do_periodic) == 1 .or. count(do_periodic) == 3)) then
         RAISE_ERROR('count(periodic) must be zero, one or three.', error)
    endif

    if (same_lattice) do_periodic = .true.

    if (any(do_periodic) .and. multiple_images(cluster_list)) then
         RAISE_ERROR("create_cluster: can't make a periodic cluster since cluster_list contains repeats", error)
    endif

    ! check for empty list

    if(cluster_list%N == 0) then
       call print('create_cluster: empty cluster_list', PRINT_NORMAL)
       return
    end if

    call print('create_cluster: Entering create_cluster', PRINT_NERD)

    if(.not.(cluster_list%intsize == 1 .or. cluster_list%intsize == 4) .or. cluster_list%realsize /= 0) then
         RAISE_ERROR("create_cluster: cluster_list table must have intsize=1 or 4 and realsize=0.", error)
    endif


    ! Cluster_info is extensible storage for the cluster
    ! It stores atomic indices and shifts (4 ints)
    ! atomic number (1 int)
    ! termination index (1 int): for termination atoms, which atom is being terminated?
    ! and atomic positions (3 reals)
    ! It's length will be at least cluster_list%N

    call print('create_cluster: Creating temporary cluster table', PRINT_NERD)
    call allocate(cluster_info,6,4,1,0,cluster_list%N)

    ! First, put all the marked atoms into cluster_info, storing their positions and shifts
    call print('create_cluster: Adding specified atoms to the cluster', PRINT_NERD)
    do i = 1, cluster_list%N
       if(cluster_list%intsize == 4) then
          ! we have shifts
          ishift = cluster_list%int(2:4,i)
       else
          ! no incoming shifts
          ishift = (/0,0,0/)
       end if
       call append(cluster_info, (/cluster_list%int(1,i),ishift,at%Z(cluster_list%int(1,i)),0/),&
	    (/at%pos(:,cluster_list%int(1,i)),1.0_dp/), (/ hybrid_mark_name(hybrid_mark(cluster_list%int(1,i))) /) )
    end do

    call print("create_cluster: cluster list:", PRINT_NERD)
    call print(cluster_info, PRINT_NERD)

    call system_timer('cluster_consistency')
    ! Next, check for various gotchas

    ! at mask is used to match with atoms already in the cluster_info table
    ! if we are periodic in a direction, we don't care about the shifts in that direction when matching
    atom_mask = (/.true.,.not.do_periodic, .true., .true./)

    connectivity_just_from_connect = .not. cluster_heuristics_nneighb_only
    if (hysteretic_connect) then
       use_connect => at%hysteretic_connect
       ! will also pass true for connectivity_just_from_connect to cluster_...() routines, so that
       ! is_nearest_neighbour() won't be used
       connectivity_just_from_connect = .true.
    else
       use_connect => at%connect
    endif

    if (cluster_allow_modification) then
       call add_property(at, 'modified_'//trim(my_mark_name), 0, ptr=modified_hybrid_mark)
       cluster_changed = .true.
       modified_hybrid_mark = hybrid_mark
       prev_cluster_info_n = cluster_info%N

       do while (cluster_changed) 
          cluster_changed = .false.
          call print("fixing up cluster according to heuristics keep_whole_residues " // keep_whole_residues // &
               ' keep_whole_prolines ' // keep_whole_prolines // &
               ' keep_whole_proline_sidechains ' // keep_whole_proline_sidechains // &
               ' keep_whole_subgroups ' // keep_whole_subgroups // &
               ' keep_whole_silica_tetrahedra ' // keep_whole_silica_tetrahedra // &
               ' keep_whole_titania_octahedra ' // keep_whole_titania_octahedra // &
               ' terminate_octahedra ' // terminate_octahedra // &
               ' reduce_n_cut_bonds ' // reduce_n_cut_bonds // &
               ' in_out_in ' // in_out_in // &
               ' in_out_in_mode ' // trim(in_out_in_mode) // &
               ' protect_X_H_bonds ' // protect_X_H_bonds // &
               ' protect_double_bonds ' // protect_double_bonds // &
               ' protect_peptide_bonds ' // protect_peptide_bonds // &
               ' terminate .or. fix_termination_clash ' // (terminate .or. fix_termination_clash) // &
               ' force_no_fix_termination_clash '// force_no_fix_termination_clash, verbosity=PRINT_NERD)
          if (keep_whole_proline_sidechains) then
             cluster_changed = cluster_changed .or. cluster_keep_whole_proline_sidechains(at, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, keep_whole_proline_sidechains_has_value)
          endif
          if (keep_whole_prolines) then
             cluster_changed = cluster_changed .or. cluster_keep_whole_prolines(at, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, keep_whole_prolines_has_value)
          endif
          if (keep_whole_subgroups) then
             cluster_changed = cluster_changed .or. cluster_keep_whole_subgroups(at, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, keep_whole_subgroups_has_value)
          endif
          if (keep_whole_residues) then
             cluster_changed = cluster_changed .or. cluster_keep_whole_residues(at, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, keep_whole_residues_has_value)
          endif
          if (keep_whole_silica_tetrahedra) then
             cluster_changed = cluster_changed .or. cluster_keep_whole_silica_tetrahedra(at, cluster_info, connectivity_just_from_connect, use_connect, atom_mask)
          end if
          if (reduce_n_cut_bonds) then
             cluster_changed = cluster_changed .or. cluster_reduce_n_cut_bonds(at, cluster_info, connectivity_just_from_connect, use_connect, atom_mask)
          endif
          if (keep_whole_titania_octahedra) then
             cluster_changed = cluster_changed .or. cluster_keep_whole_titania_octahedra(at, cluster_info, connectivity_just_from_connect, use_connect, atom_mask)
          end if
          if (in_out_in) then
             cluster_changed = cluster_changed .or. cluster_in_out_in(at, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, in_out_in_mode)
          endif
          if (protect_X_H_bonds) then
             cluster_changed = cluster_changed .or. cluster_protect_X_H_bonds(at, cluster_info, connectivity_just_from_connect, use_connect, atom_mask)
          endif
          if (protect_double_bonds) then
             cluster_changed = cluster_changed .or. cluster_protect_double_bonds(at, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, protect_double_bonds_has_value)
          endif
          if (protect_peptide_bonds) then
             cluster_changed = cluster_changed .or. cluster_protect_peptide_bonds(at, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, protect_peptide_bonds_has_value)
          endif
          if ((terminate .or. fix_termination_clash) .and. .not. force_no_fix_termination_clash) then
             cluster_changed = cluster_changed .or. cluster_fix_termination_clash(at, cluster_info, connectivity_just_from_connect, use_connect, atom_mask, termination_clash_factor)
          endif
          if (keep_whole_molecules) then
             cluster_changed = cluster_changed .or. cluster_keep_whole_molecules(at, cluster_info, connectivity_just_from_connect, use_connect, atom_mask)
          endif
       end do ! while cluster_changed

       if (prev_cluster_info_n < cluster_info%N) modified_hybrid_mark(cluster_info%int(1,prev_cluster_info_n+1:cluster_info%N)) = HYBRID_BUFFER_MARK

       call print('create_cluster: Finished fixing cluster for various heuristic pathologies',PRINT_NERD)
       call print("create_cluster: cluster list:", PRINT_NERD)
       call print(cluster_info, PRINT_NERD)
    end if ! allow_cluster_mod

    call system_timer('cluster_consistency')

    call system_timer('cluster_termination')

    !So now cluster_info contains all the atoms that are going to be in the cluster.
    !If terminate is set, we need to add terminating hydrogens along nearest neighbour bonds
    if (terminate.or.terminate_octahedra) then
       call print('create_cluster: Terminating cluster with hydrogens',PRINT_NERD)

       call allocate(n_term, 5, 0, 0, 0)
       call allocate(in_n_term, 5, 0, 0, 0)
       oldN = cluster_info%N

       !Loop over atoms in the cluster
       do n = 1, oldN

          i = cluster_info%int(1,n)
          ishift = cluster_info%int(2:4,n)
          if(terminate_octahedra .and. at%z(i).ne.8 ) cycle 
          !Loop over atom i's neighbours
          do m = 1, n_neighbours(at,i, alt_connect=use_connect)

             j = neighbour(at,i,m, r_ij, diff=dhat_ij, shift=jshift, alt_connect=use_connect)
             dhat_ij = dhat_ij / r_ij

             if (find(cluster_info,(/j,ishift+jshift,at%Z(j),0/), atom_mask) == 0 .and. &
                  (hysteretic_connect .or. is_nearest_neighbour(at, i, m, alt_connect=use_connect))) then

                ! If j is an OUT atom, and it is close enough, put a terminating hydrogen
                ! at the scaled distance between i and j

                if (.not. has_termination_rescale) then
                   rescale = termination_bond_rescale(at%Z(i), at%Z(j))
                else
                   rescale = termination_rescale
                end if
                H1 = at%pos(:,i) + rescale * r_ij * dhat_ij

                ! Label term atom with indices into original atoms structure.
                ! j is atom it's generated from and n is index into cluster table of atom it's attached to
                call append(cluster_info,(/j,ishift,1,n/),(/H1, rescale/), (/ "term      " /)) 
                 
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

                if(current_verbosity() .ge. PRINT_NERD) then
                   write(line,'(a,i0,a,i0,a)')'create_cluster: Replacing bond ',i,'--',j,' with hydrogen'
                   call print(line, PRINT_NERD)
                end if
             end if

             ! Keep track of how many other neighbours has
             if (terminate_octahedra.and.find(cluster_info,(/j,ishift+jshift,at%Z(j),0/), atom_mask) /= 0 .and. &
                  (hysteretic_connect .or. is_nearest_neighbour(at, i, m, alt_connect=use_connect))) then
                p_in = find_in_array(int_part(in_n_term,(/1,2,3,4/)),(/n,ishift/))
                if (p_in == 0) then
                   call append(in_n_term, (/n,ishift,1/))
                else
                   in_n_term%int(5,p_in) = in_n_term%int(5,p_in) + 1
                end if
             endif
          end do

       end do

       !If terminate octahedra is set we have to adjust how many H atoms are added to the oxygens.
       !Since with Ti, the oxygens are threefold coordinated, we have to remove some of the added O
       !to restore the correct twofold coordination.
       if(terminate_octahedra) then

          j = cluster_info%N 
          do i=n_term%n,1,-1
             n = n_term%int(1,i) !Reference atom
             in_n = find_in_array(int_part(in_n_term,(/1/)),(/n/))
             ishift = n_term%int(2:4,i)

             if (n_term%int(5,i) .eq. 2) then
               if (all(cluster_info%int(2:6,j) == (/ishift,1,n/))) then
                 bond1 = cluster_info%real(1:3,j)   - cluster_info%real(1:3,n)
                 bond2 = cluster_info%real(1:3,j-1) - cluster_info%real(1:3,n)
                 orig_length = sqrt(sum(bond1 * bond1))
                 
                 new_bond = bond1 + bond2 
                 new_bond = new_bond / sqrt(sum(new_bond*new_bond))
                 new_bond = new_bond * orig_length 
                 call delete(cluster_info, j) !Delete second H atom
                 call print('create_cluster: removed one of atom '//cluster_info%int(1,n)//" "//maxval(int_part(n_term,5))// &
                     ' terminating hydrogens ', PRINT_VERBOSE)
                 cluster_info%real(1:3,j-1) = cluster_info%real(1:3,n) + new_bond 
               endif
               j = j - 2
             elseif (n_term%int(5,i) .eq. 1 .and. in_n_term%int(5,in_n) .le. 1) then  !Do not remove the H atom
               j = j - 1
             elseif (n_term%int(5,i) .eq. 1 .and. in_n_term%int(5,in_n) .gt. 1) then  !Remove the H atom
               if (all(cluster_info%int(2:6,j) == (/ishift,1,n/))) then
                 call delete(cluster_info, j)
                 call print('create_cluster: removed one of atom '//cluster_info%int(1,n)//" "//maxval(int_part(n_term,5))// &
                     ' terminating hydrogens ', PRINT_VERBOSE)
               endif
               j = j - 1
             endif
          end do

       endif

       ! Do we need to remove a hydrogen atom to ensure equal n_up and n_down electrons?
       if (even_electrons .and. mod(sum(int_part(cluster_info,5)),2) == 1) then

          ! Find first atom with a maximal number of terminating hydrogens

          do i=1,n_term%n
             call append(sorted_n_term, (/n_term%int(5,i), n_term%int(1,i)/))
          end do
          allocate(idx(n_term%n))
          call sort(sorted_n_term, idx)

          n = sorted_n_term%int(2, sorted_n_term%n)
          most_hydrogens = idx(sorted_n_term%n)
          ishift = n_term%int(2:4,most_hydrogens)

          ! Loop over termination atoms
          do j=oldN,cluster_info%N
             ! Remove first H atom attached to atom i
             if (all(cluster_info%int(2:6,j) == (/ishift,1,n/))) then
                call delete(cluster_info, j)
                call print('create_cluster: removed one of atom '//cluster_info%int(1,n)//" "//maxval(int_part(n_term,5))// &
                     ' terminating hydrogens to zero total spin', PRINT_VERBOSE)
                exit 
             end if
          end do

          deallocate(idx)
       end if

       call finalise(n_term)

       call print('create_cluster: Finished terminating cluster',PRINT_NERD)
    end if
    call system_timer('cluster_termination')

    call print ('Exiting create_cluster_info', PRINT_NERD)

    call finalise(cluster_list)

  end function create_cluster_info_from_mark


  !% Given an atoms structure with a 'hybrid_mark' property, set to 
  !% 'HYBRID_ACTIVE_MARK' for some atoms, this routine adds transition 
  !% and buffer regions, and creates a 'weight_region1' property, 
  !% whose values are between 0 and 1. 
  !%
  !% Atoms marked with 'HYBRID_ACTIVE_MARK' in 'hybrid_mark' get weight
  !% 1.  Neighbour hopping is done up to transition_hops times, over
  !% which the weight linearly decreases to zero with hop count if
  !% 'weight_interpolation=hop_ramp'.  If
  !% 'weight_interpolation=distance_ramp', weight is 1 up to
  !% 'distance_ramp_inner_radius', and goes linearly to 0 by
  !% 'distance_ramp_outer_radius', where distance is calculated from
  !% 'distance_ramp_center(1:3)'.  'distance_ramp' makes the most sense if
  !% the 'HYBRID_ACTIVE_MARK' is set on a sphere with radius
  !% distance_ramp_inner_radius from distance_ramp_center, no
  !% hysteresis, and with 'transition_hops' < 0 so hops continue until
  !% no more atoms are added.  Transition atoms are marked with
  !% 'HYBRID_TRANS_MARK'. If 'hysteretic_buffer=F' (default), buffer
  !% atoms are selected by hopping 'buffer_hops' times, otherwise
  !% hysteretically with radii 'hysteretic_buffer_inner_radius' and
  !% 'hysteretic_buffer_outer_radius'.  Buffer atoms are marked with
  !% 'HYBRID_BUFFER_MARK' and given weight 0.  If 'hysteretic_connect' is
  !% set, all hops are done on bonds in the hysteretic_connect
  !% structure.

  subroutine create_hybrid_weights(at, args_str, error)
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in) :: args_str
    integer, optional, intent(out) :: error

    type(Dictionary) :: params
    logical :: has_distance_ramp_inner_radius, has_distance_ramp_outer_radius, has_distance_ramp_center
    real(dp) :: distance_ramp_inner_radius, distance_ramp_outer_radius, distance_ramp_center(3)
    logical :: min_images_only, mark_buffer_outer_layer, cluster_hopping_nneighb_only, hysteretic_buffer, hysteretic_connect, hysteretic_buffer_nneighb_only
    real(dp) :: hysteretic_buffer_inner_radius, hysteretic_buffer_outer_radius
    real(dp) :: hysteretic_connect_cluster_radius, hysteretic_connect_inner_factor, hysteretic_connect_outer_factor
    integer :: buffer_hops, transition_hops
    character(STRING_LENGTH) :: weight_interpolation, run_suffix
    logical :: construct_buffer_use_only_heavy_atoms, cluster_hopping_skip_unreachable

    integer, pointer :: hybrid_mark(:)
    real(dp), pointer :: weight_region1(:)
    integer :: n_region1, n_trans, n_region2 !, n_term
    integer :: i, j, jj, first_active, shift(3)
    logical :: dummy
    type(Table) :: activelist, currentlist, nextlist, distances, oldbuffer, bufferlist
    real(dp) :: core_CoM(3), core_mass, mass
    integer :: list_1, hybrid_number 
    type(Table) :: total_embedlist 

    real(dp) :: origin(3), extent(3,3)
    real(dp) :: save_cutoff, save_cutoff_break
    logical :: save_use_uniform_cutoff
    integer :: old_n
    integer, allocatable :: uniq_Z(:)

    logical :: distance_ramp, hop_ramp
    logical :: add_this_atom, more_hops
    integer :: cur_trans_hop
    real(dp) :: bond_len, d
    logical :: ignore_silica_residue
    integer :: res_num_silica

    INIT_ERROR(error)

! only one set of defaults now, not one in args_str and one in arg list 
    call initialise(params)
    call param_register(params, 'run_suffix', '', run_suffix, help_string="string to append to hyrbid_mark for proper mark property name")
    call param_register(params, 'transition_hops', '0', transition_hops, help_string="Number of hops to do transition region over")
    call param_register(params, 'buffer_hops', '3', buffer_hops, help_string="Number of hops for buffer region")
    call param_register(params, 'weight_interpolation', 'hop_ramp', weight_interpolation, help_string="How to ramp down weights for transition region: hop_ramp or distance_ramp")
    call param_register(params, 'distance_ramp_inner_radius', '0', distance_ramp_inner_radius, has_value_target=has_distance_ramp_inner_radius, help_string="If distance_ramp, inner radius to start reducing weight from 1")
    call param_register(params, 'distance_ramp_outer_radius', '0', distance_ramp_outer_radius, has_value_target=has_distance_ramp_outer_radius, help_string="If distance_ramp, outer radius by which to reach weight of 0")
    call param_register(params, 'distance_ramp_center', '0 0 0', distance_ramp_center, has_value_target=has_distance_ramp_center, help_string="If present, origin for distances of distance ramp (otherwise cluster center of mass)")
    call param_register(params, 'cluster_hopping_nneighb_only', 'T', cluster_hopping_nneighb_only, help_string="If true, only hop to atom pairs that are is_nearest_neighbor().")
    call param_register(params, 'cluster_hopping_skip_unreachable', 'F', cluster_hopping_skip_unreachable, help_string='If true, remove buffer atoms no longer reachable by bond hopping')
    call param_register(params, 'min_images_only', 'F', min_images_only, help_string="If true, consider only minimum images, not multiple periodic images")
    call param_register(params, 'mark_buffer_outer_layer', 'T', mark_buffer_outer_layer, help_string="If true, mark outermost buffer layer")
    call param_register(params, 'hysteretic_buffer', 'F', hysteretic_buffer, help_string="If true, do hysteretic buffer")
    call param_register(params, 'hysteretic_buffer_inner_radius', '5.0', hysteretic_buffer_inner_radius, help_string="Inner radius for hysteretic buffer")
    call param_register(params, 'hysteretic_buffer_outer_radius', '7.0', hysteretic_buffer_outer_radius, help_string="Outer radius for hysteretic buffer")
    call param_register(params, 'hysteretic_buffer_nneighb_only', 'F', hysteretic_buffer_nneighb_only, help_string="If true, only hop to nearest neighbour pairs when constructing hysteretic buffer")
    call param_register(params, 'hysteretic_connect', 'F', hysteretic_connect, help_string="If true, use hysteretic connect object")
    call param_register(params, 'hysteretic_connect_cluster_radius', '1.2', hysteretic_connect_cluster_radius, help_string="If using hysteretic connect, use as cluster radius for origin and extent estimates")
    call param_register(params, 'hysteretic_connect_inner_factor', '1.2', hysteretic_connect_inner_factor, help_string="Inner factor (multiplier for covalent bond distance) for hysteretic connect")
    call param_register(params, 'hysteretic_connect_outer_factor', '1.5', hysteretic_connect_outer_factor, help_string="Outer factor (multiplier for covalent bond distance) for hysteretic connect")
    call param_register(params, 'construct_buffer_use_only_heavy_atoms', 'F', construct_buffer_use_only_heavy_atoms, help_string="If true, use only non-H atoms for constructing buffer")
    call param_register(params, 'ignore_silica_residue', 'F', ignore_silica_residue, help_string="If true, do special things for silica") !lam81
    call param_register(params, 'res_num_silica', '1', res_num_silica, help_string="Residue number for silica") !lam81
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='create_hybrid_weights args_str') ) then
      RAISE_ERROR("create_hybrid_weights failed to parse args_str='"//trim(args_str)//"'", error)
    endif
    call finalise(params)

    hop_ramp = .false.
    distance_ramp = .false.
    if (trim(weight_interpolation) == 'hop_ramp') then
       hop_ramp = .true.
    else if (trim(weight_interpolation) == 'distance_ramp') then
       distance_ramp = .true.
    else
       RAISE_ERROR('create_hybrid_weights: unknown weight_interpolation value: '//trim(weight_interpolation), error)
    end if

    call print('create_hybrid_weights: transition_hops='//transition_hops//' buffer_hops='//buffer_hops//' weight_interpolation='//weight_interpolation, PRINT_VERBOSE)
    call print('  cluster_hopping_nneighb_only='//cluster_hopping_nneighb_only//'  min_images_only='//min_images_only//' mark_buffer_outer_layer='//mark_buffer_outer_layer, PRINT_VERBOSE)
    call print('  hysteretic_buffer='//hysteretic_buffer//' hysteretic_buffer_inner_radius='//hysteretic_buffer_inner_radius, PRINT_VERBOSE)
    call print('  hysteretic_buffer_outer_radius='//hysteretic_buffer_outer_radius//' hysteretic_buffer_nneighb_only='//hysteretic_buffer_nneighb_only, PRINT_VERBOSE)
    call print('  hysteretic_connect='//hysteretic_connect//' hysteretic_connect_cluster_radius='//hysteretic_connect_cluster_radius, PRINT_VERBOSE)
    call print('  hysteretic_connect_inner_factor='//hysteretic_connect_inner_factor //' hysteretic_connect_outer_factor='//hysteretic_connect_outer_factor, PRINT_VERBOSE)


    ! check to see if atoms has a 'weight_region1' property already, if so, check that it is compatible, if not present, add it
    if(assign_pointer(at, 'weight_region1'//trim(run_suffix), weight_region1)) then
       weight_region1 = 0.0_dp
    else
       call add_property(at, 'weight_region1'//trim(run_suffix), 0.0_dp)
       dummy = assign_pointer(at, 'weight_region1'//trim(run_suffix), weight_region1)
    end if

    ! check for a compatible hybrid_mark property. it must be present
    if(.not.assign_pointer(at, 'hybrid_mark'//trim(run_suffix), hybrid_mark)) then
       RAISE_ERROR('create_hybrid_weights: atoms structure has no "hybrid_mark'//trim(run_suffix)//'" property', error)
    endif

    ! Fast implementation of trivial case where buffer_hops=0 and transition_hops=0
    if (buffer_hops == 0 .and. transition_hops == 0) then

       where (hybrid_mark /= HYBRID_ACTIVE_MARK)
          hybrid_mark = HYBRID_NO_MARK
       end where

       weight_region1 = 0.0_dp
       where (hybrid_mark == HYBRID_ACTIVE_MARK)
          weight_region1 = 1.0_dp
       end where

       return
    end if

    ! Add first marked atom to activelist. shifts will be relative to this atom
    first_active = find_in_array(hybrid_mark, HYBRID_ACTIVE_MARK)

    n_region1 = count(hybrid_mark == HYBRID_ACTIVE_MARK)
    list_1 = 0  
    call allocate(activelist, 4,0,0,0)
    call allocate(total_embedlist, 4,0,0,0) 

    do while (total_embedlist%N < n_region1)

       call wipe(currentlist)
       call wipe(activelist)
       call wipe(nextlist)    

       call append(activelist, (/first_active,0,0,0/))
       call append(currentlist, activelist)
       weight_region1(first_active) = 1.0_dp

       if (distance_ramp) then ! we might need ACTIVE_MARK center of mass
          if (has_property(at, 'mass')) then
             core_mass = at%mass(first_active)
          else
             core_mass = ElementMass(at%Z(first_active))
          end if
          core_CoM = core_mass*at%pos(:,first_active) ! we're taking this as reference atom so no shift needed here
       end if

       if (hysteretic_connect) then
         call system_timer('hysteretic_connect')
         call estimate_origin_extent(at, hybrid_mark == HYBRID_ACTIVE_MARK, hysteretic_connect_cluster_radius, origin, extent)
         save_use_uniform_cutoff = at%use_uniform_cutoff
         save_cutoff = at%cutoff
         save_cutoff_break = at%cutoff_break
         call set_cutoff_factor(at, hysteretic_connect_inner_factor, hysteretic_connect_outer_factor)
         call calc_connect_hysteretic(at, at%hysteretic_connect, origin, extent)
         if (save_use_uniform_cutoff) then
           call set_cutoff(at, save_cutoff, save_cutoff_break)
         else
           call set_cutoff_factor(at, save_cutoff, save_cutoff_break)
         endif
         call system_timer('hysteretic_connect')
       endif

       ! Add other active atoms using bond hopping from the first atom
       ! in the cluster to find the other active atoms and hence to determine the correct 
       ! periodic shifts
       
       hybrid_number = 1
       do while (hybrid_number .ne. 0)
          if (hysteretic_connect) then
            call BFS_step(at, currentlist, nextlist, nneighb_only=cluster_hopping_nneighb_only, min_images_only = min_images_only, alt_connect=at%hysteretic_connect)
          else
            call BFS_step(at, currentlist, nextlist, nneighb_only=cluster_hopping_nneighb_only, min_images_only = min_images_only, property=hybrid_mark)
          endif
          hybrid_number = 0 
          do j=1,nextlist%N
             jj = nextlist%int(1,j) 
             shift = nextlist%int(2:4,j)
             if (hybrid_mark(jj) == HYBRID_ACTIVE_MARK .and. find_in_array(activelist%int(1,1:activelist%N), jj) == 0) then
                hybrid_number = hybrid_number+1 
                call append(activelist, nextlist%int(:,j))
                weight_region1(jj) = 1.0_dp

                if (distance_ramp) then ! we might need ACTIVE_MARK center of mass
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
       enddo
       list_1 = list_1 + activelist%N
       call append(total_embedlist, activelist)    

       if (distance_ramp) then ! calculate actual distance ramp parameters 
	 ! distance_ramp center is as specified, otherwise ACTIVE_MARK center of mass
	 core_CoM = core_CoM/core_mass
	 if (.not. has_distance_ramp_center) distance_ramp_center = core_CoM

	 ! distance_ramp_inner_radius is as specified, otherwise distance of atom furthest from distance_ramp_center
	 if (.not. has_distance_ramp_inner_radius) then
	   call initialise(distances, 1, 1, 0, 0)
	   do i=1, activelist%N
	     jj = activelist%int(1,i)
	     call append(distances, jj, distance_min_image(at, distance_ramp_center, jj))
	   end do
	   distance_ramp_inner_radius = maxval(distances%real(1,1:distances%N))
	   call finalise(distances)
	 endif

	 if (.not. has_distance_ramp_outer_radius) then
	   distance_ramp_outer_radius = 0.0_dp
	   call uniq(at%Z, uniq_Z)
	   do i=1, size(uniq_Z)
	   do j=1, size(uniq_Z)
	     bond_len = bond_length(uniq_Z(i), uniq_Z(j))
	     if (bond_len > distance_ramp_outer_radius) distance_ramp_outer_radius = bond_len 
	   end do
	   end do
	   if (transition_hops <= 0) &
	     call print("WARNING: using transition_hops for distance_ramp outer radius, but transition_hops="//transition_hops,PRINT_ALWAYS)
	   distance_ramp_outer_radius = distance_ramp_inner_radius + distance_ramp_outer_radius*transition_hops
	 endif
       endif ! distance_ramp

       call wipe(currentlist)
       call append(currentlist, activelist)

       ! create transition region
       call print('create_hybrid_mark: creating transition region',PRINT_VERBOSE)
       n_trans = 0

       cur_trans_hop = 1
       if (distance_ramp) transition_hops = -1 ! distance ramp always does as many hops as needed to get every atom within outer radius
       more_hops = (transition_hops < 0 .or. transition_hops >= 2)
       do while (more_hops)
	 more_hops = .false.
         if (hysteretic_connect) then
           call BFS_step(at, currentlist, nextlist, nneighb_only = cluster_hopping_nneighb_only .and. (transition_hops > 0), min_images_only = min_images_only, alt_connect=at%hysteretic_connect)
         else
           call BFS_step(at, currentlist, nextlist, nneighb_only = cluster_hopping_nneighb_only .and. (transition_hops > 0), min_images_only = min_images_only)
         endif

         call wipe(currentlist)
         do j = 1,nextlist%N
            jj = nextlist%int(1,j)
            if(hybrid_mark(jj) == HYBRID_NO_MARK) then
	       add_this_atom = .true.
	       if (distance_ramp) then
                  d = distance_min_image(at, distance_ramp_center, jj)
		  if (d >= distance_ramp_outer_radius) add_this_atom = .false.
	       endif

	       if (add_this_atom) then
		 if (hop_ramp) weight_region1(jj) = 1.0_dp - real(cur_trans_hop,dp)/real(transition_hops,dp) ! linear transition
		 if (distance_ramp) then        ! Save distance, weight will be calculated later
		   if (d <= distance_ramp_inner_radius) then
		     weight_region1(jj) = 1.0_dp
		   else
		     weight_region1(jj) = 1.0_dp - (d-distance_ramp_inner_radius)/(distance_ramp_outer_radius-distance_ramp_inner_radius)
		   endif
		 endif
		 call append(currentlist, nextlist%int(:,j))
		 hybrid_mark(jj) = HYBRID_TRANS_MARK
		 n_trans = n_trans+1
		 more_hops = .true.
	       end if
            end if
         end do ! do j=1,nextlist%N

	 cur_trans_hop = cur_trans_hop + 1
	 if (transition_hops >= 0 .and. cur_trans_hop >= transition_hops) more_hops = .false.
       end do ! more_hops

       if (list_1 < n_region1) then
          call print('searching for a new quantum zone as found '//list_1//' atoms, need to get to '//n_region1, PRINT_VERBOSE)
          do i =1, at%N
             if (hybrid_mark(i) == HYBRID_ACTIVE_MARK .and. .not. Is_in_Array(total_embedlist%int(1,1:total_embedlist%N), i)) then
                first_active = i
                exit
             endif
          enddo
       endif

    enddo ! while (total_embedlist%N < n_region1)

    ! create region2 (buffer region) 
    if (.not. hysteretic_buffer) then
       ! since no hysteresis, safe to reset non ACTIVE/TRANS atoms to HYBRID_NO_MARK
       where (hybrid_mark /= HYBRID_ACTIVE_MARK .and. hybrid_mark /= HYBRID_TRANS_MARK) hybrid_mark = HYBRID_NO_MARK

       n_region2 = 0
       do i = 0,buffer_hops-1
	  if (hysteretic_connect) then
	    call BFS_step(at, currentlist, nextlist, nneighb_only = cluster_hopping_nneighb_only, min_images_only = min_images_only, alt_connect=at%hysteretic_connect)
	  else
	    call BFS_step(at, currentlist, nextlist, nneighb_only = cluster_hopping_nneighb_only, min_images_only = min_images_only)
	  endif
          call wipe(currentlist)
          do j = 1,nextlist%N
             jj = nextlist%int(1,j)
             if(hybrid_mark(jj) == HYBRID_NO_MARK) then
                call append(currentlist, nextlist%int(:,j))
                weight_region1(jj) = 0.0_dp
                if (i==buffer_hops-1 .and. mark_buffer_outer_layer) then
                   hybrid_mark(jj) = HYBRID_BUFFER_OUTER_LAYER_MARK
                else
                   hybrid_mark(jj) = HYBRID_BUFFER_MARK
                end if
                n_region2 = n_region2+1
             end if
          end do
       end do
    else 
       ! hysteretic buffer here

       call initialise(oldbuffer, 1,0,0,0)
       !call wipe(oldbuffer)
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
       !
       ! This will fail if marked atoms do not form connected clusters around the active atoms
       old_n = bufferlist%N
       do while (bufferlist%N < n_region2)
          if (hysteretic_connect) then
             call BFS_step(at, currentlist, nextlist, nneighb_only = cluster_hopping_nneighb_only, min_images_only = min_images_only, alt_connect=at%hysteretic_connect)
          else
             call BFS_step(at, currentlist, nextlist, nneighb_only = cluster_hopping_nneighb_only, min_images_only = min_images_only, property =hybrid_mark)
          endif
          do j=1,nextlist%N
             jj = nextlist%int(1,j)
             shift = nextlist%int(2:4,j)
             if (hybrid_mark(jj) /= HYBRID_NO_MARK) call append(bufferlist, nextlist%int(:,j))
          end do
          call append(currentlist, nextlist)

          ! check that cluster is still growing
          if (bufferlist%N == old_n) then
             if (cluster_hopping_skip_unreachable) then
                call print('WARNING: skipping buffer atoms no longer connected')
                exit
             else
                RAISE_ERROR('create_hybrid_weights: buffer cluster stopped growing before all marked atoms found - check for split QM or buffer region', error)
             end if
          end if
          old_n = bufferlist%N
       end do
       n_region2 = bufferlist%N ! possibly excluding some no-longer-reachable atoms

       ! Remove marks on all buffer atoms
       do i=1,oldbuffer%N
          if (hybrid_mark(oldbuffer%int(1,i)) == HYBRID_BUFFER_MARK .or. &
              hybrid_mark(oldbuffer%int(1,i)) == HYBRID_BUFFER_OUTER_LAYER_MARK) &
              hybrid_mark(oldbuffer%int(1,i)) = HYBRID_NO_MARK
       end do

       !construct the hysteretic buffer region:
       if (hysteretic_connect) then
	 call print("create_hybrid_weights calling construct_hysteretic_region", verbosity=PRINT_NERD)
	 call construct_hysteretic_region(region=bufferlist,at=at,core=total_embedlist,loop_atoms_no_connectivity=.false., &
	   inner_radius=hysteretic_buffer_inner_radius,outer_radius=hysteretic_buffer_outer_radius,use_avgpos=.false., &
	   add_only_heavy_atoms=construct_buffer_use_only_heavy_atoms, cluster_hopping_nneighb_only=cluster_hopping_nneighb_only, &
           force_hopping_nneighb_only=hysteretic_buffer_nneighb_only, min_images_only=min_images_only, &
	   alt_connect=at%hysteretic_connect, &
           ignore_silica_residue=ignore_silica_residue, res_num_silica=res_num_silica, error=error) !NB, debugfile=mainlog) lam81
       else
	 call print("create_hybrid_weights calling construct_hysteretic_region", verbosity=PRINT_NERD)
	 call construct_hysteretic_region(region=bufferlist,at=at,core=total_embedlist,loop_atoms_no_connectivity=.false., &
	   inner_radius=hysteretic_buffer_inner_radius,outer_radius=hysteretic_buffer_outer_radius,use_avgpos=.false., &
	   add_only_heavy_atoms=construct_buffer_use_only_heavy_atoms, cluster_hopping_nneighb_only=cluster_hopping_nneighb_only, &
           force_hopping_nneighb_only=hysteretic_buffer_nneighb_only, min_images_only=min_images_only, &
           ignore_silica_residue=ignore_silica_residue, res_num_silica=res_num_silica, error=error) !NB, debugfile=mainlog) lam81
       endif

       call print('bufferlist=',PRINT_VERBOSE)
       call print(bufferlist,PRINT_VERBOSE)

       ! Mark new buffer region, leaving core QM region alone
       ! at the moment only ACTIVE and NO marks are present, because hybrid_mark was set to =hybrid
       do i=1,bufferlist%N
          if (hybrid_mark(bufferlist%int(1,i)) == HYBRID_NO_MARK) & 
               hybrid_mark(bufferlist%int(1,i)) = HYBRID_BUFFER_MARK

          ! Marking the outer layer with  HYBRID_BUFFER_OUTER_LAYER_MARK is
          ! dealt with in create_cluster_from_mark.

       end do

       call finalise(bufferlist)
       call finalise(oldbuffer)
!       call finalise(embedlist)

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

    call print('create_hybrid_weights: '//list_1//' region 1, '//n_trans//' transition, '//n_region2//&
         ' region 2, '//count(hybrid_mark /= HYBRID_NO_MARK)//' in total', PRINT_VERBOSE)
    !call print('create_hybrid_weights: '//n_region1//' region 1, '//n_trans//' transition, '//n_region2//&
    !     ' region 2, '//count(hybrid_mark /= HYBRID_NO_MARK)//' in total', PRINT_VERBOSE)
    !    call sort(total_embedlist)

    call finalise(activelist)
    call finalise(currentlist)
    call finalise(nextlist)
    call finalise(distances)
    call finalise(total_embedlist)

  end subroutine create_hybrid_weights

  !% Return estimated `(origin, extent)` for hysteretic connectivity 
  !% calculator to include all atoms within a distance `cluster_radius`
  !% of an atom in the `active` array.
  subroutine estimate_origin_extent(at, active, cluster_radius, origin, extent)
    type(Atoms), intent(in) :: at
    logical, intent(in) :: active(:)
    real(dp), intent(in) :: cluster_radius
    real(dp), intent(out) :: origin(3), extent(3,3)

    real(dp) :: center(3), low_corner(3), high_corner(3), dr(3)
    integer :: i, n_active, first_active
    logical :: found_first_active

    found_first_active = .false.
    n_active = 0
    do i=1, at%N
      if (.not. active(i)) cycle
      n_active = n_active + 1
      if (found_first_active) then
	center = center + diff_min_image(at, first_active, i)
      else
	center = 0.0_dp
	first_active = i
	found_first_active = .true.
      endif
    end do
    center = center / real(n_active, dp)
    center = center + at%pos(:,first_active)

    call print("estimate_origin_extent: got center" // center, verbosity=PRINT_VERBOSE)

    low_corner = 1.0e38_dp
    high_corner = -1.0e38_dp
    do i=1, at%N
      if (.not. active(i)) cycle
      dr = diff_min_image(at, at%pos(:,i), center)
      low_corner = min(low_corner, dr)
      high_corner = max(high_corner, dr)
    end do
    call print("estimate_origin_extent: got relative low_corner" // low_corner, verbosity=PRINT_NERD)
    call print("estimate_origin_extent: got relative high_corner" // high_corner, verbosity=PRINT_NERD)
    low_corner = low_corner + center
    high_corner = high_corner + center

    call print("estimate_origin_extent: got low_corner" // low_corner, verbosity=PRINT_NERD)
    call print("estimate_origin_extent: got high_corner" // high_corner, verbosity=PRINT_NERD)

    origin = low_corner - cluster_radius
    extent = 0.0_dp
    extent(1,1) = (high_corner(1)-low_corner(1))+2.0_dp*cluster_radius
    extent(2,2) = (high_corner(2)-low_corner(2))+2.0_dp*cluster_radius
    extent(3,3) = (high_corner(3)-low_corner(3))+2.0_dp*cluster_radius

    call print("estimate_origin_extent: got origin" // origin, verbosity=PRINT_VERBOSE)
    call print("estimate_origin_extent: got extent(1,:) " // extent(1,:), verbosity=PRINT_VERBOSE)
    call print("estimate_origin_extent: got extent(2,:) " // extent(2,:), verbosity=PRINT_VERBOSE)
    call print("estimate_origin_extent: got extent(3,:) " // extent(3,:), verbosity=PRINT_VERBOSE)

  end subroutine estimate_origin_extent

  !% Given an Atoms structure with an active region marked in the mark_name (default 'hybrid_mark')
  !% property using 'HYBRID_ACTIVE_MARK', grow the embed region by 'fit_hops'
  !% bond hops to form a fit region. Returns the embedlist and fitlist with correct
  !% periodic shifts.
  subroutine create_embed_and_fit_lists(at, fit_hops, embedlist, fitlist, cluster_hopping_nneighb_only, min_images_only, mark_name, error)

    type(Atoms), intent(inout) :: at
    integer :: fit_hops
    type(Table), intent(out) :: embedlist, fitlist
    logical, intent(in), optional :: cluster_hopping_nneighb_only, min_images_only
    character(len=*), intent(in), optional :: mark_name
    integer, optional, intent(out) :: error

    character(len=255)                :: my_mark_name
    integer, pointer :: hybrid_mark(:)
    integer :: n_region1, n_region2 !, n_term
    integer :: i, j, jj, first_active, shift(3)
    logical :: do_hopping_nneighb_only, do_min_images_only
    type(Table) :: currentlist, nextlist, tmpfitlist
    integer :: n, hybrid_number
    type(Table) :: totallist
    integer :: list_1

    integer :: old_n

    INIT_ERROR(error)

    my_mark_name = optional_default('hybrid_mark', mark_name)

    call print('Entered create_embed_and_fit_lists.',PRINT_VERBOSE)
    do_hopping_nneighb_only = optional_default(.false., cluster_hopping_nneighb_only)
    do_min_images_only = optional_default(.true., min_images_only)

    ! check for a compatible hybrid_mark property. it must be present
    if(.not.assign_pointer(at, trim(my_mark_name), hybrid_mark)) then
       RAISE_ERROR('create_fit_region: atoms structure has no "'//trim(my_mark_name)//'"', error)
    endif

    call wipe(embedlist)
    call wipe(fitlist)

    ! Add first marked atom to embedlist. shifts will be relative to this atom
    first_active = find_in_array(hybrid_mark, HYBRID_ACTIVE_MARK)

    n_region1 = count(hybrid_mark == HYBRID_ACTIVE_MARK)
    list_1 = 0
    call allocate(totallist,4,0,0,0)
    call allocate(currentlist,4,0,0,0)

    ! Add other active atoms using bond hopping from the first atom
    ! in the cluster to find the other active atoms and hence to determine the correct 
    ! periodic shifts
    !
    old_n = embedlist%N
    do while (embedlist%N < n_region1)

      call wipe(currentlist)
      call wipe(nextlist)
      call wipe(totallist)

      call append(totallist, (/first_active,0,0,0/))
      call print('create_embed_and_fit_lists: expanding quantum region starting from '//first_active//' atom', PRINT_VERBOSE)
      call append(currentlist, totallist)

      hybrid_number = 1 
      do while (hybrid_number .ne. 0) 
!!!CHECK THIS
       call BFS_step(at, currentlist, nextlist, nneighb_only = .false., min_images_only = do_min_images_only, property =hybrid_mark)

       hybrid_number = 0 
       do j=1,nextlist%N
          jj = nextlist%int(1,j)
          shift = nextlist%int(2:4,j)
          if (hybrid_mark(jj) == HYBRID_ACTIVE_MARK) then
                hybrid_number = hybrid_number + 1 
                call append(totallist, nextlist%int(:,j))
          endif
       end do
       call append(currentlist, nextlist)
      enddo

      list_1 = list_1 + totallist%N
      call append(embedlist, totallist)
      call print('create_embed_and_fit_lists: number of atoms in the embedlist is '//embedlist%N, PRINT_VERBOSE)

      if (list_1 .lt. n_region1) then
         n = 0
         do i =1, at%N
            if (hybrid_mark(i) == HYBRID_ACTIVE_MARK) then
               n = n+1
               if (.not. Is_in_Array(int_part(embedlist,1), i)) then
                  first_active = i
                  exit
               endif
            endif
         enddo
      endif

       ! check that cluster is still growing
       !if (embedlist%N == old_n) then
       !     RAISE_ERROR('create_embed_and_fit_lists: (embedlist) cluster stopped growing before all marked atoms found - check for split QM region', error)
       !endif
       !old_n = embedlist%N
    end do


    call wipe(currentlist)
    call append(currentlist, embedlist)


    ! create region2 (fit region)
    call initialise(tmpfitlist,4,0,0,0,0)
    n_region2 = 0
    do i = 0,fit_hops-1
       call BFS_step(at, currentlist, nextlist, nneighb_only = do_hopping_nneighb_only, min_images_only = do_min_images_only)
       do j = 1,nextlist%N
          jj = nextlist%int(1,j)
          if(hybrid_mark(jj) /= HYBRID_ACTIVE_MARK) then
            call append(tmpfitlist, nextlist%int(:,j))
            n_region2 = n_region2+1
          end if
       end do
       call append(currentlist, nextlist)
    end do

    call print('create_embed_and_fit_lists: '//list_1//' embed, '//n_region2//' fit', PRINT_VERBOSE)

    ! Sort order to we are stable to changes in neighbour ordering introduced
    ! by calc_connect. 
    call sort(embedlist)
    call sort(tmpfitlist)

    ! fitlist consists of sorted embedlist followed by sorted list of remainder of fit atoms
    call append(fitlist, embedlist)
    call append(fitlist, tmpfitlist)

    if (do_min_images_only) then
       if (multiple_images(embedlist)) then
          call discard_non_min_images(embedlist)
          call print('create_embed_and_fits_lists: multiple images discarded from embedlist', PRINT_VERBOSE)
       endif

       if (multiple_images(fitlist)) then
          call discard_non_min_images(fitlist)
          call print('create_embed_and_fits_lists: multiple images discarded from fitlist', PRINT_VERBOSE)
       endif
    endif

    call finalise(currentlist)
    call finalise(nextlist)
    call finalise(tmpfitlist)

  end subroutine create_embed_and_fit_lists

  !% Given an Atoms structure with an active region marked in the mark_name (default 'cluster_mark')
  !% property using 'HYBRID_ACTIVE_MARK', and buffer region marked using 'HYBRID_BUFFER_MARK',
  !% 'HYBRID_TRANS_MARK' or 'HYBRID_BUFFER_OUTER_LAYER_MARK', simply returns the embedlist and fitlist
  !% according to 'cluster_mark'. It does not take into account periodic shifts.
  subroutine create_embed_and_fit_lists_from_cluster_mark(at,embedlist,fitlist, mark_name, error)

    type(Atoms), intent(in)  :: at
    type(Table), intent(out) :: embedlist, fitlist
    character(len=*), intent(in), optional :: mark_name
    integer, optional, intent(out) :: error

    character(len=255) :: my_mark_name
    type(Table)              :: tmpfitlist

    INIT_ERROR(error)

    my_mark_name = optional_default('cluster_mark', mark_name)

    call print('Entered create_embed_and_fit_lists_from_cluster_mark.',PRINT_VERBOSE)
    call wipe(embedlist)
    call wipe(fitlist)
    call wipe(tmpfitlist)

    !build embed list from ACTIVE atoms
    call list_matching_prop(at,embedlist,trim(my_mark_name),HYBRID_ACTIVE_MARK)

    !build fitlist from BUFFER and TRANS atoms
    call list_matching_prop(at,tmpfitlist,trim(my_mark_name),HYBRID_BUFFER_MARK)
    call append(fitlist,tmpfitlist)
    call list_matching_prop(at,tmpfitlist,trim(my_mark_name),HYBRID_TRANS_MARK)
    call append(fitlist,tmpfitlist)
    call list_matching_prop(at,tmpfitlist,trim(my_mark_name),HYBRID_BUFFER_OUTER_LAYER_MARK)
    call append(fitlist,tmpfitlist)

    call wipe(tmpfitlist)
    call append(tmpfitlist,fitlist)

    ! Sort in order to we are stable to changes in neighbour ordering introduced
    ! by calc_connect. 
    call sort(embedlist)
    call sort(tmpfitlist)

    ! fitlist consists of sorted embedlist followed by sorted list of remainder of fit atoms
    call wipe(fitlist)
    call append(fitlist, embedlist)
    call append(fitlist, tmpfitlist)

    call print('Embedlist:',PRINT_ANAL)
    call print(int_part(embedlist,1),PRINT_ANAL)
    call print('Fitlist:',PRINT_ANAL)
    call print(int_part(fitlist,1),PRINT_ANAL)

    call finalise(tmpfitlist)
    call print('Leaving create_embed_and_fit_lists_from_cluster_mark.',PRINT_VERBOSE)

  end subroutine create_embed_and_fit_lists_from_cluster_mark

  !% Return the atoms in a hysteretic region:
  !% To become part of the 'list' region, atoms must drift within the
  !% 'inner' list. To leave the 'list' region, an atom
  !% must drift further than 'outer'.
  !% Optionally use time averaged positions
  !
  subroutine update_hysteretic_region(at,inner,outer,list,min_images_only,error)

    type(Atoms),       intent(in)    :: at
    type(Table),       intent(inout) :: inner, outer
    type(Table),       intent(inout) :: list
    logical, optional, intent(in)    :: min_images_only
    integer, optional, intent(out)   :: error

    integer n, my_verbosity
    logical do_min_images_only
    integer, allocatable :: cols(:)

    INIT_ERROR(error)

    my_verbosity = current_verbosity()
    do_min_images_only = optional_default(.false., min_images_only)
    call print('update_hystertic_region: do_min_images_only='//do_min_images_only, PRINT_VERBOSE)

    if (do_min_images_only) then
       ! only compare atom indices - first column in table
       allocate(cols(1))
       cols(1) = 1
    else
       ! compare atom indices and shifts - all four columns in table
       allocate(cols(4))
       cols(:) = (/1, 2, 3, 4/)
    end if

    if (my_verbosity == PRINT_ANAL) then
       call print('In Select_Hysteretic_Quantum_Region:')
       call print('List currently contains '//list%N//' atoms')
    end if

    if ((list%intsize /= inner%intsize) .or. (list%intsize /= outer%intsize)) then
       RAISE_ERROR('update_hysteretic_region: inner, outer and list must have the same intsize', error)
    endif

    ! Speed up searching
    call sort(outer)
    call sort(list)

    !Check for atoms in 'list' and not in 'outer'
    do n = list%N, 1, -1
       if (search(outer,list%int(cols,n))==0) then
          call delete(list,n)
          if (my_verbosity > PRINT_NORMAL) call print('Removed atom ('//list%int(cols,n)//') from quantum list')
       end if
    end do

    call sort(list)

    !Check for new atoms in 'inner' cluster and add them to list
    do n = 1, inner%N
       if (search(list,inner%int(cols,n))==0) then
          call append(list,inner%int(:,n)) ! copy shifts columns as well, even if we're not comparing them
          call sort(list)
          if (my_verbosity > PRINT_NORMAL) call print('Added atom ('//inner%int(cols,n)//') to quantum list')
       end if
    end do

    if (my_verbosity >= PRINT_NORMAL) call print(list%N//' atoms selected for quantum treatment')
    deallocate(cols)

  end subroutine update_hysteretic_region


  !
  !% Given an atoms object, and a point 'centre' or a list of atoms 'core',
  !% fill the 'region' table with all atoms hysteretically within 'inner_radius -- outer_radius'
  !% of the 'centre' or any 'core' atom, which can be reached by connectivity hopping.
  !% In the case of building the 'region' around 'centre', simply_loop_over_atoms loops instead of
  !% bond hopping.
  !% Optionally use the time averaged positions.
  !% Optionally use only heavy atom selection.
  !
  subroutine construct_hysteretic_region(region,at,core,centre,loop_atoms_no_connectivity,inner_radius,outer_radius,use_avgpos,add_only_heavy_atoms,cluster_hopping_nneighb_only, &
                                         force_hopping_nneighb_only, min_images_only,alt_connect,debugfile, ignore_silica_residue, res_num_silica, error) !lam81

    type(Table),           intent(inout) :: region
    type(Atoms),           intent(in)    :: at
    type(Table), optional, intent(in)    :: core
    real(dp),    optional, intent(in)    :: centre(3)
    logical,     optional, intent(in)    :: loop_atoms_no_connectivity
    real(dp),              intent(in)    :: inner_radius
    real(dp),              intent(in)    :: outer_radius
    logical,     optional, intent(in)    :: use_avgpos
    logical,     optional, intent(in)    :: add_only_heavy_atoms
    logical,     optional, intent(in)  :: cluster_hopping_nneighb_only
    logical,     optional, intent(in)  :: force_hopping_nneighb_only
    logical,     optional, intent(in)  :: min_images_only
    type(Connection), optional,intent(in) :: alt_connect
    type(inoutput), optional :: debugfile
    integer, optional, intent(out) :: error

    type(Table)                          :: inner_region
    type(Table)                          :: outer_region
    logical                              :: no_hysteresis

    logical,     optional, intent(in)    :: ignore_silica_residue ! lam81
    integer,     optional, intent(in)    :: res_num_silica        ! lam81

    INIT_ERROR(error)

    if (present(debugfile)) call print("construct_hysteretic_region radii " // inner_radius // " " // outer_radius, file=debugfile)
    !check the input arguments only in construct_region.
    if (inner_radius.lt.0._dp) then
      RAISE_ERROR('inner_radius must be > 0 and it is '//inner_radius, error)
    endif
    if (outer_radius.lt.inner_radius) then
      RAISE_ERROR('outer radius ('//outer_radius//') must not be smaller than inner radius ('//inner_radius//').', error)
    endif

    no_hysteresis = .false.
    if ((outer_radius-inner_radius).lt.epsilon(0._dp)) no_hysteresis = .true.
    if (inner_radius .feq. 0.0_dp) call print('WARNING: hysteretic region with inner_radius .feq. 0.0', PRINT_ALWAYS)
    if (no_hysteresis) call print('WARNING! construct_hysteretic_region: inner_buffer=outer_buffer. no hysteresis applied: outer_region = inner_region used.',PRINT_ALWAYS)
    if (present(debugfile)) call print("   no_hysteresis " // no_hysteresis, file=debugfile)

    if (present(debugfile)) call print("   constructing inner region", file=debugfile)
    call construct_region(region=inner_region,at=at,core=core,centre=centre,loop_atoms_no_connectivity=loop_atoms_no_connectivity, &
      radius=inner_radius,use_avgpos=use_avgpos,add_only_heavy_atoms=add_only_heavy_atoms,cluster_hopping_nneighb_only=cluster_hopping_nneighb_only, &
      force_hopping_nneighb_only=force_hopping_nneighb_only, &
      min_images_only=min_images_only,alt_connect=alt_connect,debugfile=debugfile, &
      ignore_silica_residue=ignore_silica_residue, res_num_silica=res_num_silica) !lam81
    if (no_hysteresis) then
       call initialise(outer_region,4,0,0,0,0)
       call append(outer_region,inner_region)
    else
    if (present(debugfile)) call print("   constructing outer region", file=debugfile)
       call construct_region(region=outer_region,at=at,core=core,centre=centre,loop_atoms_no_connectivity=loop_atoms_no_connectivity, &
	radius=outer_radius, use_avgpos=use_avgpos,add_only_heavy_atoms=add_only_heavy_atoms,cluster_hopping_nneighb_only=cluster_hopping_nneighb_only, &
        force_hopping_nneighb_only=force_hopping_nneighb_only, &
	min_images_only=min_images_only,alt_connect=alt_connect,debugfile=debugfile, &
        ignore_silica_residue=ignore_silica_residue, res_num_silica=res_num_silica) ! lam81
    endif

    if (present(debugfile)) call print("   orig inner_region list", file=debugfile)
    if (present(debugfile)) call print(inner_region, file=debugfile)
    if (present(debugfile)) call print("   orig outer_region list", file=debugfile)
    if (present(debugfile)) call print(outer_region, file=debugfile)
    if (present(debugfile)) call print("   old region list", file=debugfile)
    if (present(debugfile)) call print(region, file=debugfile)

     call print('construct_hysteretic_region: old region list:',PRINT_VERBOSE)
     call print(region,PRINT_VERBOSE)

    call update_hysteretic_region(at,inner_region,outer_region,region, min_images_only=min_images_only)

    call print('construct_hysteretic_region: inner_region list:',PRINT_VERBOSE)
    call print(inner_region,PRINT_VERBOSE)
    call print('construct_hysteretic_region: outer_region list:',PRINT_VERBOSE)
    call print(outer_region,PRINT_VERBOSE)

    if (present(debugfile)) call print("   new inner region list", file=debugfile)
    if (present(debugfile)) call print(inner_region, file=debugfile)
    if (present(debugfile)) call print("   new outer region list", file=debugfile)
    if (present(debugfile)) call print(outer_region, file=debugfile)

    call print('construct_hysteretic_region: new region list:',PRINT_VERBOSE)
    call print(region,PRINT_VERBOSE)

    if (present(debugfile)) call print("   new region list", file=debugfile)
    if (present(debugfile)) call print(region, file=debugfile)

    call finalise(inner_region)
    call finalise(outer_region)

  end subroutine construct_hysteretic_region

  !
  !% Given an atoms object, and a point or a list of atoms in the first integer of
  !% the 'core' table, fill the 'buffer' table with all atoms within 'radius'
  !% of any core atom (which can be reached by connectivity hopping) or the point
  !%(with bond hopping or simply looping once over the atoms).
  !% Optionally use the time averaged positions (only for the radius).
  !% Optionally use a heavy atom based selection (applying to both the core and the region atoms).
  !% Alternatively use the hysteretic connection, only nearest neighbours and/or min_images (only for the n_connectivity_hops).
  !
  subroutine construct_region(region,at,core,centre,loop_atoms_no_connectivity,radius,n_connectivity_hops,use_avgpos,add_only_heavy_atoms,cluster_hopping_nneighb_only,&
                              force_hopping_nneighb_only, min_images_only,alt_connect,debugfile, ignore_silica_residue, res_num_silica, error) !lam81

    type(Table),           intent(out) :: region
    type(Atoms),           intent(in)  :: at
    type(Table), optional, intent(in)  :: core
    real(dp),    optional, intent(in)  :: centre(3)
    logical,     optional, intent(in)  :: loop_atoms_no_connectivity
    real(dp),    optional, intent(in)  :: radius
    integer,     optional, intent(in)  :: n_connectivity_hops
    logical,     optional, intent(in)  :: use_avgpos
    logical,     optional, intent(in)  :: add_only_heavy_atoms
    logical,     optional, intent(in)  :: cluster_hopping_nneighb_only
    logical,     optional, intent(in)  :: force_hopping_nneighb_only
    logical,     optional, intent(in)  :: min_images_only
    type(Connection), optional,intent(in) :: alt_connect
type(inoutput), optional :: debugfile
    integer, optional, intent(out) :: error

    logical                             :: do_use_avgpos
    logical                             :: do_loop_atoms_no_connectivity
    logical                             :: do_add_only_heavy_atoms
    logical                             :: do_hopping_nneighb_only, do_force_hopping_nneighb_only
    logical                             :: do_min_images_only
    integer                             :: i, j, ii, ij
    type(Table)                         :: nextlist
    logical                             :: more_hops, add_i
    integer                             :: cur_hop
    real(dp), pointer                   :: use_pos(:,:)
    integer                             ::  shift_i(3)
    logical :: do_ignore_silica_residue, got_atom_res_number

    integer, pointer :: atom_res_number(:) !lam81
    logical,     optional, intent(in)  :: ignore_silica_residue !lam81
    integer,     optional, intent(in)  :: res_num_silica

    INIT_ERROR(error)

    do_loop_atoms_no_connectivity = optional_default(.false.,loop_atoms_no_connectivity)
    do_ignore_silica_residue = optional_default(.false., ignore_silica_residue)

    do_add_only_heavy_atoms = optional_default(.false.,add_only_heavy_atoms)
    if (do_add_only_heavy_atoms .and. .not. has_property(at,'Z')) then
      RAISE_ERROR("construct_region: atoms has no Z property", error)
    endif

    got_atom_res_number = assign_pointer(at, 'atom_res_number', atom_res_number)

    do_use_avgpos = optional_default(.false.,  use_avgpos)

    if (count((/ present(centre), present(core) /)) /= 1) then
      RAISE_ERROR("Need either centre or core, but not both present(centre) " // present(centre) // " present(core) "// present(core), error)
    endif

    ! if we have radius, we'll have to decide what positions to use
    if (present(radius)) then
      if (do_use_avgpos) then
        if (.not. assign_pointer(at, "avgpos", use_pos)) then
	  RAISE_ERROR("do_use_avgpos is true, but no avgpos property", error)
	endif
      else
        if (.not. assign_pointer(at, "pos", use_pos)) then
	  RAISE_ERROR("do_use_avgpos is false, but no pos property", error)
	endif
      endif
    endif

    call initialise(region,4,0,0,0,0)

    if (do_loop_atoms_no_connectivity) then
      if (present(debugfile)) call print("   loop over atoms", file=debugfile)
      if (.not. present(radius)) then
	RAISE_ERROR("do_loop_atoms_no_connectivity=T requires radius", error)
      endif
      if (present(n_connectivity_hops) .or. present(min_images_only) .or. present(cluster_hopping_nneighb_only)) &
	call print("WARNING: do_loop_atoms_no_connectivity, but specified unused arg n_connectivity_hops " // present(n_connectivity_hops) // &
	  " min_images_only " // present(min_images_only) // " cluster_hopping_nneighb_only " // present(cluster_hopping_nneighb_only), PRINT_ALWAYS)
       ! jrk33: disabled this warning - should be fine with small cells since distance_min_image is used.
       ! call print('WARNING: check if your cell is greater than the radius, looping only works in that case.',PRINT_ALWAYS)
       !if (any((/at%lattice(1,1),at%lattice(2,2),at%lattice(3,3)/) < radius)) then
       !RAISE_ERROR('too small cell', error)
       !endifx
       do i = 1,at%N
	 if (present(debugfile)) call print("check atom i " // i // " Z " // at%Z(i) // " pos " // at%pos(:,i), file=debugfile)
	 if (do_add_only_heavy_atoms .and. at%Z(i) == 1) cycle
	 if (present(centre)) then ! distance from centre is all that matters
	   if (present(debugfile)) call print(" distance_min_image use_pos i " // use_pos(:,i) // " centre " // centre // " dist " // distance_min_image(at, use_pos(:,i), centre(1:3), shift_i(1:3)), file=debugfile)
	   if (distance_min_image(at,use_pos(:,i),centre(1:3),shift_i(1:3)) < radius) then
	     if (present(debugfile)) call print("   adding i " // i  // " at " // use_pos(:,i) // " center " // centre // " dist " // distance_min_image(at, use_pos(:,i), centre(1:3), shift_i(1:3)), file=debugfile)
	     call append(region,(/i,shift_i(1:3)/))
	    endif
	  else ! no centre, check distance from each core list atom
	    do ij=1, core%N
	      j = core%int(1,ij)
	      if (present(debugfile)) call print(" core atom ij " // ij // " j " // j // " distance_min_image use_pos i " // use_pos(:,i) // " use_pos j " // use_pos(:,j) // " dist " // distance_min_image(at, use_pos(:,i), use_pos(:,j), shift_i(1:3)), file=debugfile)
	      if (distance_min_image(at,use_pos(:,i),use_pos(:,j),shift_i(1:3)) < radius) then
		if (present(debugfile)) call print("   adding i " // i // " at " // use_pos(:,i) // " near j " // j // " at " // use_pos(:,j) // " dist " // distance_min_image(at, use_pos(:,i), use_pos(:,j), shift_i(1:3)), file=debugfile)
	        call append(region,(/i,shift_i(1:3)/))
		exit
	      endif
	    end do ! j
	  endif ! no centre, use core
       enddo ! i

    else ! do_loop_atoms_no_connectivity

      if (.not. present(core)) then
	RAISE_ERROR("do_loop_atoms_no_connectivity is false, trying connectivity hops, but no core list is specified", error)
      endif

      if (present(debugfile)) call print("   connectivity hopping", file=debugfile)
      if (present(debugfile) .and. present(radius)) call print("    have radius " // radius, file=debugfile)
      if (present(debugfile)) call print("   present cluster_hopping_nneighb_only " // present(cluster_hopping_nneighb_only), file=debugfile)
      if (present(debugfile) .and. present(cluster_hopping_nneighb_only)) call print("   cluster_hopping_nneighb_only " // cluster_hopping_nneighb_only, file=debugfile)
      do_hopping_nneighb_only = optional_default(.true., cluster_hopping_nneighb_only)
      do_force_hopping_nneighb_only = optional_default(.false., force_hopping_nneighb_only)
      do_min_images_only = optional_default(.true., min_images_only)
      if (do_use_avgpos) then
	RAISE_ERROR("can't use avgpos with connectivity hops - make sure your connectivity is based on pos instead", error)
      endif

      if (do_ignore_silica_residue) then
         call print_warning('Overriding min_images_only since ignore_silica_residue=T')
         do_min_images_only = .true. !lam81 
      end if
      ! start with core
      call append(region,core)

      more_hops = .true.
      cur_hop = 1
      do while (more_hops)
	if (present(debugfile)) call print('   construct_region do_hopping_nneighb_only = ' // do_hopping_nneighb_only // ' do_min_images_only = '//do_min_images_only, file=debugfile)
	if (present(debugfile)) call print('   doing hop ' // cur_hop, file=debugfile)
	if (present(debugfile)) call print('   cutoffs ' // at%cutoff // ' ' // at%use_uniform_cutoff, file=debugfile)
	more_hops = .false.
	call bfs_step(at, region, nextlist, nneighb_only=(do_hopping_nneighb_only .and. .not. present(radius)) .or. do_force_hopping_nneighb_only, &
                      min_images_only=do_min_images_only, max_r=radius, alt_connect=alt_connect, debugfile=debugfile)
	if (present(debugfile)) call print("   bfs_step returned nextlist%N " // nextlist%N, file=debugfile)
	if (present(debugfile)) call print(nextlist, file=debugfile)
	if (nextlist%N /= 0) then ! go over things in next hop
	  do ii=1, nextlist%N
	    i = nextlist%int(1,ii)
	    add_i = .true.
	    if (do_add_only_heavy_atoms) then
	      if (at%Z(i) == 1) cycle
	    endif
            if (do_ignore_silica_residue .and. got_atom_res_number) then     ! lam81
               if (atom_res_number(i) == res_num_silica) cycle ! lam81
            endif
	    if (present(debugfile)) call print("  i " // i // " is heavy or all atoms requested", file=debugfile)
	    if (present(radius)) then ! check to make sure we're close enough to core list
	      if (present(debugfile)) call print("  radius is present, check distance from core list", file=debugfile)
	      add_i = .false.
	      do ij=1, core%N
		j = core%int(1,ij)
		if (present(debugfile)) call print("  distance of ii " // ii // " i " // i // " at  " // use_pos(:,i) // " from ij " // ij // " j " // j // " at " // use_pos(:,j) // " is " // distance_min_image(at,use_pos(:,i),use_pos(:,j),shift_i(1:3)), file=debugfile)
		if (distance_min_image(at,use_pos(:,i),use_pos(:,j),shift_i(1:3)) <= radius) then
		  if (present(debugfile)) call print("  decide to add i", file=debugfile)
		  add_i = .true.
		  exit ! from ij loop
		endif
	      end do ! ij
	    endif ! present(radius)
	    if (add_i) then
	      if (present(debugfile)) call print("  actually adding i", file=debugfile)
	      more_hops = .true.
	      call append(region, nextlist%int(1:4,ii))
	    endif
	  end do
	else ! nextlist%N == 0
	  if (present(n_connectivity_hops)) then
	    if (n_connectivity_hops > 0) &
	      call print("WARNING: n_connectivity_hops = " // n_connectivity_hops // " cur_hop " // cur_hop // " but no atoms added")
	  endif
	endif ! nextlist%N

	cur_hop = cur_hop + 1
	if (present(n_connectivity_hops)) then
	  if (cur_hop > n_connectivity_hops) more_hops = .false.
	endif
      end do ! while more_hops

    endif ! do_loop_atoms_no_connectivity
    if (present(debugfile)) call print("leaving construct_region()", file=debugfile)

  end subroutine construct_region


  subroutine update_active(this, nneightol, avgpos, reset, error)
    type(Atoms) :: this
    real(dp), optional :: nneightol
    logical, optional :: avgpos, reset
    integer, optional, intent(out) :: error

    type(Atoms), save :: nn_atoms
    integer, pointer, dimension(:) :: nn, old_nn, active
    integer  :: i
    real(dp) :: use_nn_tol
    logical  :: use_avgpos, do_reset

    INIT_ERROR(error)

    use_nn_tol = optional_default(this%nneightol, nneightol)
    use_avgpos = optional_default(.true., avgpos)
    do_reset   = optional_default(.false., reset)

    ! First time copy entire atoms structure into nn_atoms and
    ! add properties for nn, old_nn and active
    if (reset .or. .not. is_initialised(nn_atoms)) then
       nn_atoms = this
       call add_property(this, 'nn', 0)
       call add_property(this, 'old_nn', 0)
       call add_property(this, 'active', 0)

       call set_cutoff_factor(nn_atoms, use_nn_tol)
    end if

    if (this%N /= nn_atoms%N) then
         RAISE_ERROR('update_actives: Number mismatch between this%N ('//this%N//') and nn_atoms%N ('//nn_atoms%N//')', error)
    endif

    if (.not. assign_pointer(this, 'nn', nn)) then
         RAISE_ERROR('update_actives: Atoms is missing "nn" property', error)
    endif

    if (.not. assign_pointer(this, 'old_nn', old_nn)) then
         RAISE_ERROR('update_actives: Atoms is missing "old_nn" property', error)
    endif

    if (.not. assign_pointer(this, 'active', active)) then
         RAISE_ERROR('update_actives: Atoms is missing "active" property', error)
    endif

    call print('update_actives: recalculating nearest neighbour table', PRINT_VERBOSE)
    if (use_avgpos .and. associated(this%avgpos)) then
       call print('update_actives: using time averaged atomic positions', PRINT_VERBOSE)
       nn_atoms%pos = this%avgpos
    else
       call print('update_actives: using instantaneous atomic positions', PRINT_VERBOSE)
       nn_atoms%pos = this%pos
    end if
    call calc_connect(nn_atoms)

    nn = 0
    do i = 1,nn_atoms%N
       nn(i) = n_neighbours(nn_atoms, i)
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
  subroutine add_cut_hydrogens(this,qmlist,heuristics_nneighb_only,verbosity,alt_connect)

    type(Atoms),       intent(in),          target :: this
    type(Table),       intent(inout)               :: qmlist
    integer, optional, intent(in)                  :: verbosity
    logical, optional, intent(in)                  :: heuristics_nneighb_only
    type(Connection), intent(in), optional, target :: alt_connect

    type(Table)                :: neighbours, bonds, centre
    logical                    :: more_atoms
    integer                    :: i, j, n, nn, added
    type(Connection), pointer :: use_connect
    logical :: my_heuristics_nneighb_only

    ! Check for atomic connectivity
    if (present(alt_connect)) then
      use_connect => alt_connect
    else
      use_connect => this%connect
    endif

    my_heuristics_nneighb_only = optional_default(.true., heuristics_nneighb_only)

    more_atoms = .true.
    added = 0
    call allocate(centre,4,0,0,0,1)

    !Repeat while the search trigger is true
    do while(more_atoms)

       more_atoms = .false.

       !Find nearest neighbours of the cluster
       call bfs_step(this,qmlist,neighbours,nneighb_only=my_heuristics_nneighb_only,min_images_only=.true.,alt_connect=use_connect)

       !Loop over neighbours
       do n = 1, neighbours%N

          i = neighbours%int(1,n)

          call wipe(centre)
          call append(centre,(/i,0,0,0/))

          ! Find atoms bonded to this neighbour
          call bfs_step(this,centre,bonds,nneighb_only=my_heuristics_nneighb_only,min_images_only=.true.,alt_connect=use_connect)

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
                      if (verbosity >= PRINT_NORMAL) call print('Add_Cut_Hydrogens: Added atom '//i//', neighbour of atom '//j)
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
       if (verbosity >= PRINT_NORMAL) then
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
    integer, save                        :: CUBIC_BONDLENGTH_SQ_FUNC
    logical, save                        :: first_call = .true.
    integer                              :: g, j, n, a1, a2, datalength
    real(dp)                             :: x0,y0,y0p,x1,y1,y1p, coeffs(4), t, dE, m1, m2, meff
    real(dp), save                       :: Etot = 0.0_dp

    !Register the constraint function if this is the first call
    if (first_call) then
       CUBIC_BONDLENGTH_SQ_FUNC = Register_Constraint(CUBIC_BONDLENGTH_SQ)
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
             call ds_amend_constraint(this,j,CUBIC_BONDLENGTH_SQ_FUNC,(/coeffs,x0,x1/))

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
    call inverse(A,A_inv)
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

  !% Updates the core QM flags saved in 'hybrid' and 'hybrid_mark'
  !% properties.  Do this hysteretically, from 'R_inner' to 'R_outer'
  !% around 'origin' or 'atomlist', that is the centre of the QM
  !% region (a position in space or a list of atoms).  Also saves the
  !% old 'hybrid_mark' property in 'old_hybrid_mark' property.
  !% optionally correct selected region with heuristics (as coded in
  !% :func:`create_cluster_info_from_mark`())
  !
  subroutine create_pos_or_list_centred_hybrid_region(my_atoms,R_inner,R_outer,origin, atomlist,use_avgpos,add_only_heavy_atoms, &
	     cluster_hopping_nneighb_only,cluster_heuristics_nneighb_only,min_images_only,use_create_cluster_info, create_cluster_info_args, list_changed, mark_postfix, &
             ignore_silica_residue, res_num_silica, error) ! lam81

    type(Atoms),        intent(inout) :: my_atoms
    real(dp),           intent(in)    :: R_inner
    real(dp),           intent(in)    :: R_outer
    real(dp), optional, intent(in)    :: origin(3)
    type(Table), optional, intent(in)    :: atomlist !the seed of the QM region
    logical,  optional, intent(in)   :: use_avgpos, add_only_heavy_atoms, cluster_hopping_nneighb_only, cluster_heuristics_nneighb_only, min_images_only, use_create_cluster_info
    character(len=*), optional, intent(in) :: create_cluster_info_args
    logical,  optional, intent(out)   :: list_changed
    character(len=*),  optional, intent(in)   :: mark_postfix
    integer, optional, intent(out) :: error

    logical ::  do_hopping_nneighb_only, do_heuristics_nneighb_only
    logical :: my_use_create_cluster_info
    type(Atoms) :: atoms_for_add_cut_hydrogens
    type(Table) :: core, old_core, old_all_but_term
    integer, pointer :: hybrid_p(:), hybrid_mark_p(:), old_hybrid_mark_p(:)
    character(len=STRING_LENGTH) :: my_create_cluster_info_args, my_mark_postfix
    integer, pointer :: hybrid_region_core_tmp_p(:)
    type(Table) :: padded_cluster_info

    logical,  optional, intent(in)   :: ignore_silica_residue ! lam81
    integer,  optional, intent(in)   :: res_num_silica        ! lam81

    INIT_ERROR(error)

    my_use_create_cluster_info = optional_default(.false., use_create_cluster_info)

    do_hopping_nneighb_only = optional_default(.false., cluster_hopping_nneighb_only)
    do_heuristics_nneighb_only = optional_default(.true., cluster_heuristics_nneighb_only)

    if (count((/present(origin),present(atomlist)/))/=1) then
      RAISE_ERROR('create_pos_or_list_centred_hybrid_mark: Exactly 1 of origin and atomlist must be present.', error)
    endif

    my_mark_postfix = optional_default("", trim(mark_postfix))

    call map_into_cell(my_atoms)

    ! save various subsets of atoms based on old hybrid marks
    call allocate(core,4,0,0,0,0)
    call allocate(old_core,4,0,0,0,0)
    call allocate(old_all_but_term,4,0,0,0,0)
    call get_hybrid_list(my_atoms,old_core,active_trans_only=.true., int_property='hybrid_mark'//trim(my_mark_postfix), error=error)
    PASS_ERROR_WITH_INFO("create_pos_or_list_centred_hybrid_region getting old core", error)
    call get_hybrid_list(my_atoms,core,active_trans_only=.true., int_property='hybrid_mark'//trim(my_mark_postfix), error=error)
    PASS_ERROR_WITH_INFO("create_pos_or_list_centred_hybrid_region getting core", error)
    call get_hybrid_list(my_atoms,old_all_but_term, all_but_term=.true., int_property='hybrid_mark'//trim(my_mark_postfix), error=error)
    PASS_ERROR_WITH_INFO("create_pos_or_list_centred_hybrid_region getting all but term", error)

   ! Build the new hysteretic QM core based on distances from point/list
   if (present(atomlist)) then
     call print("create_pos_or_list_centred_hybrid_region calling construct_hysteretic_region", verbosity=PRINT_NERD)
     call construct_hysteretic_region(region=core,at=my_atoms,core=atomlist,loop_atoms_no_connectivity=.false., &
       inner_radius=R_inner,outer_radius=R_outer, use_avgpos=use_avgpos, add_only_heavy_atoms=add_only_heavy_atoms, &
       cluster_hopping_nneighb_only=cluster_hopping_nneighb_only, min_images_only=min_images_only, &
       ignore_silica_residue=ignore_silica_residue, res_num_silica=res_num_silica, error=error) !NB , debugfile=mainlog) lam81
   else !present origin
     call print("create_pos_or_list_centred_hybrid_region calling construct_hysteretic_region", verbosity=PRINT_NERD)
     call construct_hysteretic_region(region=core,at=my_atoms,centre=origin,loop_atoms_no_connectivity=.true., &
       inner_radius=R_inner,outer_radius=R_outer, use_avgpos=use_avgpos, add_only_heavy_atoms=add_only_heavy_atoms, &
       cluster_hopping_nneighb_only=cluster_hopping_nneighb_only, min_images_only=min_images_only, &
       ignore_silica_residue=ignore_silica_residue, res_num_silica=res_num_silica, error=error) !NB , debugfile=mainlog) lam81
   endif
   PASS_ERROR_WITH_INFO("create_pos_or_list_centred_hybrid_region constructing hysteretic region", error)

   ! call create_cluster_info_from_mark, if requested, for various chemical intuition fixes to core
   if (my_use_create_cluster_info) then
      call add_property(my_atoms, "hybrid_region_core_tmp", HYBRID_NO_MARK, ptr=hybrid_region_core_tmp_p)
      hybrid_region_core_tmp_p(int_part(core,1)) = HYBRID_ACTIVE_MARK
      if (present(create_cluster_info_args) .and. (present(cluster_hopping_nneighb_only) .or. present(cluster_heuristics_nneighb_only))) then
	 RAISE_ERROR("Got both create_cluster_info_args and (do_hopping_nneighb_only or do_heuristics_nneighb_only), but these conflict", error)
      endif
      my_create_cluster_info_args = optional_default("terminate=F cluster_hopping_nneighb_only="//do_hopping_nneighb_only// &
	 " cluster_heuristics_nneighb_only="//do_heuristics_nneighb_only//" cluster_allow_modification", create_cluster_info_args)
      padded_cluster_info = create_cluster_info_from_mark(my_atoms, my_create_cluster_info_args, mark_name="hybrid_region_core_tmp")
      call finalise(core)
      call initialise(core, Nint=4, Nreal=0, Nstr=0, Nlogical=0, max_length=padded_cluster_info%N)
      call append(core, blank_rows=padded_cluster_info%N)
      core%int(1:4,1:padded_cluster_info%N) = padded_cluster_info%int(1:4,1:padded_cluster_info%N)
      call remove_property(my_atoms,"hybrid_region_core_tmp")
   endif

!TO BE OPTIMIZED : add avgpos to add_cut_hydrogen
   ! add cut hydrogens, according to avgpos
    atoms_for_add_cut_hydrogens = my_atoms
    atoms_for_add_cut_hydrogens%oldpos = my_atoms%avgpos
    atoms_for_add_cut_hydrogens%avgpos = my_atoms%avgpos
    atoms_for_add_cut_hydrogens%pos = my_atoms%avgpos

    call set_cutoff_factor(atoms_for_add_cut_hydrogens,DEFAULT_NNEIGHTOL)
    call calc_connect(atoms_for_add_cut_hydrogens)

    call add_cut_hydrogens(atoms_for_add_cut_hydrogens,core)
    call finalise(atoms_for_add_cut_hydrogens)

    ! check changes in QM list and set the new QM list
    if (present(list_changed)) then
       list_changed = check_list_change(old_list=old_core,new_list=core)
       if (list_changed)  call print('QM list around the origin  has changed')
    endif

    ! sadd old hybrid_mark property
    call add_property(my_atoms,'old_hybrid_mark'//trim(my_mark_postfix),HYBRID_NO_MARK)
    if (.not. assign_pointer(my_atoms,'old_hybrid_mark'//trim(my_mark_postfix),old_hybrid_mark_p)) then
      RAISE_ERROR("create_pos_or_list_centred_hybrid_region couldn't get old_hybrid_mark"//trim(my_mark_postfix)//" property", error)
    endif
    ! update QM_flag of my_atoms
    if (.not. assign_pointer(my_atoms,'hybrid_mark'//trim(my_mark_postfix),hybrid_mark_p)) then
      RAISE_ERROR("create_pos_or_list_centred_hybrid_region couldn't get hybrid_mark"//trim(my_mark_postfix)//" property", error)
    endif
    ! save old_hybrid_mark
    old_hybrid_mark_p(1:my_atoms%N) = hybrid_mark_p(1:my_atoms%N)
    ! default to NO_MARK
    hybrid_mark_p(1:my_atoms%N) = HYBRID_NO_MARK
    ! anything which was marked before (except termination) is now BUFFER (so hysteretic buffer, later, has correct previous values)
    hybrid_mark_p(int_part(old_all_but_term,1)) = HYBRID_BUFFER_MARK
    ! anything returned in the core list is now ACTIVE
    hybrid_mark_p(int_part(core,1)) = HYBRID_ACTIVE_MARK

   ! update hybrid property of my_atoms
    if (.not. assign_pointer(my_atoms,'hybrid'//trim(my_mark_postfix),hybrid_p)) then
      RAISE_ERROR("create_pos_or_list_centred_hybrid_region couldn't get hybrid"//trim(my_mark_postfix)//" property", error)
    endif
    hybrid_p(1:my_atoms%N) = 0
    hybrid_p(int_part(core,1)) = 1

    call finalise(core)
    call finalise(old_core)
    call finalise(old_all_but_term)

  end subroutine create_pos_or_list_centred_hybrid_region

!  !% Returns a $hybridlist$ table with the atom indices whose $cluster_mark$
!  !% (or optionally any $int_property$) property takes no greater than $hybridflag$ positive value.
!  !
!  subroutine get_hybrid_list(my_atoms,hybridflag,hybridlist,int_property,get_up_to_mark_value)
!
!    type(Atoms), intent(in)  :: my_atoms
!    integer,     intent(in)  :: hybridflag
!    type(Table), intent(out) :: hybridlist
!    character(len=*), optional, intent(in) :: int_property
!    logical, intent(in), optional :: get_up_to_mark_value
!
!    integer, pointer :: mark_p(:)
!    integer              :: i
!    logical :: do_get_up_to_mark_value
!
!    do_get_up_to_mark_value = optional_default(.false., get_up_to_mark_value)
!
!    if (present(int_property)) then
!       if (.not. assign_pointer(my_atoms, trim(int_property), mark_p)) then
!	 RAISE_ERROR("get_hybrid_list_int couldn't get int_property='"//trim(int_property)//"'", error)
!       endif
!    else
!       if (.not. assign_pointer(my_atoms, 'cluster_mark', mark_p)) then
!	 RAISE_ERROR("get_hybrid_list_int couldn't get default int_property='cluster_mark'", error)
!       endif
!    endif
!
!    call initialise(hybridlist,4,0,0,0,0)      !1 int, 0 reals, 0 str, 0 log, num_hybrid_atoms entries
!    if (do_get_up_to_mark_value) then
!       do i=1,my_atoms%N
!          ! if (my_atoms%data%int(hybrid_flag_index,i).gt.0.and. &
!             ! my_atoms%data%int(hybrid_flag_index,i).le.hybridflag) &
!          if (mark_p(i) > 0 .and.  mark_p(i) <= hybridflag) &
!               call append(hybridlist,(/i,0,0,0/))
!       enddo
!    else
!       do i=1,my_atoms%N
!          ! if (my_atoms%data%int(hybrid_flag_index,i).eq.hybridflag) &
!          if (mark_p(i) == hybridflag) call append(hybridlist,(/i,0,0,0/))
!       enddo
!    endif
!
!    if (hybridlist%N.eq.0) call print('Empty QM list with cluster_mark '//hybridflag,verbosity=PRINT_SILENT)
!
!  end subroutine get_hybrid_list

  !% Checks and reports the changes between two tables, $old_qmlist$ and $new_qmlist$.
  !
  function check_list_change(old_list,new_list) result(list_changed)

    type(Table), intent(in) :: old_list, new_list
    integer :: i
    logical :: list_changed

    list_changed = .false.
    if (old_list%N.ne.new_list%N) then
       call print ('list has changed: new number of atoms is: '//new_list%N//', was: '//old_list%N)
       list_changed = .true.
    else
       if (any(old_list%int(1,1:old_list%N).ne.new_list%int(1,1:new_list%N))) then
           do i=1,old_list%N
              if (.not.find_in_array(int_part(old_list,1),(new_list%int(1,i))).gt.0) then
                 call print('list has changed: atom '//new_list%int(1,i)//' has entered the region')
                 list_changed = .true.
              endif
              if (.not.find_in_array(int_part(new_list,1),(old_list%int(1,i))).gt.0) then
                 call print('list has changed: atom '//old_list%int(1,i)//' has left the region')
                 list_changed = .true.
              endif
           enddo
       endif
    endif

  end function check_list_change

  !% return list of atoms that have various subsets of hybrid marks set
  subroutine get_hybrid_list(at,hybrid_list,all_but_term,active_trans_only,int_property, error)
    type(Atoms), intent(in)  :: at !% object to scan for marked atoms
    type(Table), intent(out) :: hybrid_list !% on return, list of marked atoms
    logical, intent(in), optional :: all_but_term !% if present and true, select all marked atoms that aren't TERM
    logical, intent(in), optional :: active_trans_only !% if present and true, select all atoms marked ACTIVE or TRANS
    !% exactly one of all_but_term and active_trans_only must be present and true
    character(len=*), optional, intent(in) :: int_property !% if present, property to check, default cluster_mark
    integer, optional, intent(out) :: error

    integer :: i
    integer, pointer :: hybrid_mark(:)
    logical              :: my_all_but_term, my_active_trans_only
    character(STRING_LENGTH) :: my_int_property

    INIT_ERROR(error)

    if (.not. present(all_but_term) .and. .not. present(active_trans_only)) then
      RAISE_ERROR("get_hybrid_list called with neither all_but_term nor active_trans_only present", error)
    endif

    my_all_but_term = optional_default(.false., all_but_term)
    my_active_trans_only = optional_default(.false., active_trans_only)

    if ((my_all_but_term .and. my_active_trans_only) .or. (.not. my_all_but_term .and. .not. my_active_trans_only)) then
      RAISE_ERROR("get_hybrid_list needs exactly one of all_but_term=" // all_but_term // " and active_trans_only="//my_active_trans_only, error)
    endif

    my_int_property = ''
    if (present(int_property)) then
       my_int_property = trim(int_property)
    else
       my_int_property = "cluster_mark"
    endif
    if (.not.(assign_pointer(at, trim(my_int_property), hybrid_mark))) then
      RAISE_ERROR("get_hybrid_list couldn't find "//trim(my_int_property)//" field", error)
    endif

    call initialise(hybrid_list,4,0,0,0,0)      !1 int, 0 reals, 0 str, 0 log, num_qm_atoms entries
    do i=1, at%N
      if (my_all_but_term) then
	if (hybrid_mark(i) /= HYBRID_NO_MARK .and. hybrid_mark(i) /= HYBRID_TERM_MARK) call append(hybrid_list,(/i,0,0,0/))
      else if (my_active_trans_only) then
	if (hybrid_mark(i) == HYBRID_ACTIVE_MARK .or. hybrid_mark(i) == HYBRID_TRANS_MARK) call append(hybrid_list,(/i,0,0,0/))
      else
	RAISE_ERROR("impossible! get_hybrid_list has no selection mode set", error)
      endif
    end do

    if (hybrid_list%N.eq.0) call print('get_hybrid_list returns empty hybrid list with field '//trim(my_int_property)// &
                                   ' all_but_term ' // my_all_but_term // ' active_trans_only ' // my_active_trans_only ,PRINT_ALWAYS)
  end subroutine get_hybrid_list


end module clusters_module
