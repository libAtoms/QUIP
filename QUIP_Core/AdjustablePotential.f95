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

!X
!X     Learn-on-the-fly (LOTF) hybrid molecular dynamics code
!X
!X    
!X     Authors: Gabor Csanyi, Alessio Commisso, Steven Winfield
!X     James Kermode, Gianpietro Moras, Michael Payne, Alessandro De Vita
!X     
!X     Copyright 2005, All Rights Reserved 
!X
!X     This source code is confidential, all distribution is 
!X     prohibited. Making unauthorized copies is also prohibited
!X
!X     When using this software, the following should be referenced:
!X
!X     Gabor Csanyi, Tristan Albaret, Mike C. Payne and Alessandro De Vita
!X     "Learn on the fly": a hybrid classical and quantum-mechanical
!X         molecular dynamics simulation
!X     Physical Review Letters 93 p. 175503 (2004) >>PDF [626 KB]
!X
!X     Gabor Csanyi, T. Albaret, G. Moras, M. C. Payne, A. De Vita
!X     Multiscale hybrid simulation methods for material systems
!X     J. Phys. Cond. Mat. 17 R691-R703 Topical Review (2005)
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X  AdjustablePotential_linear module
!X  
!X  a set of linear potentials that can be adapted a reproduce 
!X  given forces on atoms
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


module AdjustablePotential_module
  use libAtoms_module
  implicit none
  SAVE
  

  ! ratio of smallest to largest eigenvalue of spring space spanning matrix
  ! for an atom to be considered not well described
  real(dp), parameter :: AP_bad_atom_threshold = 0.1_dp

  ! Threshold for spring constants above which a warning will be issued
  real(dp), parameter :: adjustable_potential_max_spring_constant = 1.0_dp

  integer, parameter :: adjustable_potential_nparams = 1

  ! are we initialised or not
  logical::adjustable_potential_initialised = .false.

  ! have we got some parameters saved
  logical::adjustable_potential_parameters_saved = .false.

  ! this table contains the pairs of atoms for which we have a two-body
  ! adjustable spring potential, together with their distance
  ! and the potential parameters. the parameters are just the y values of the spring
  ! !!!the integer indices refer to the order of atoms in the fitlist
  ! argument of the init() routine!!!
  !
  !                    int               real
  ! 
  !fitlist indexes ->  fi,fj,shift       rij, fij <- fij: variable parameter
  type(table)::twobody

  ! this table is used as a store for old data which can be used to initialise
  ! parameters to values from the last interpolation step. It also stores
  ! parameters from one end point during interpolation
  type(table)::twobody_old 

  real(dp), dimension(:,:), allocatable:: fitforce_old

  !The exclusion list table is used to prevent springs being created between
  !certain pairs of atoms, e.g. if there is a bond length constraint between
  !them.
  type(table) :: exclusion_list
  logical     :: exclusions = .false.

  ! Operation                                            twobody   twobody_old
  !
  ! Extrapolate from t1 -> t2 using parameters P1           P1  --,     --
  ! Compute parameters P2' at t2                            P2'   `-->  P1
  ! Interpolate from t1 -> t2 using P1 and P2'              P2' --,     P1
  ! Call CalcConnect                                        --    `-->  P2'
  ! Create parameters P2 using P2' as initial values        P2          P2'
  ! Extrapolate from t2 -> t3 using parameters P2           P2          --
  ! ...

  ! local copy of target force to fit to: F1x,F1yF1z,F2x,F2y,F2z...
  real(dp), allocatable, dimension(:)   ::target_force

  ! local copy of fitlist indices
  ! atomlist%int(1,fi)   = i,    fi: index in fitlist, i: index in Atoms structure
  type(Table) :: atomlist, atomlist_old
!!$  type(Table) :: fitlist_old

  ! forcematrix multiplied onto the vector of
  ! parameters (fij...) 
  ! gives the force vector on the atoms. the number of columns
  ! is the number of parameters, the number of rows is 3*the number of atoms
  !
  ! force = forcematrix * params 
  !
  ! ...0..y1_x y2_x...0......       
  ! ...0..y1_y y2_y...0......
  ! ...0..y1_z y2_z...0......      
  type(sparse)                          ::forcematrix ! linear part

contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X Initialize adjustable potential
  !X
  !X inputs are an atoms structure and a table with a single list integers 
  !X of atoms for which target forces are to be fitted.
  !X
  !X The adjustable potential will be fitted to the difference between
  !X a quantum mechanical force and a classical force.
  !X
  !X An adjustable potential will be created between every pair of neighbouring
  !X atoms that are given in the fitlist.
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! 

  ! this is called with a table
  subroutine adjustable_potential_init(atoms_in, fitlist, directionN, spring_hops, &
       map, method, nnonly)
    type(Atoms),       intent(in) :: atoms_in
    type(Table),       intent(in) :: fitlist
    integer,           intent(in) :: directionN
    integer, optional, intent(in) :: spring_hops
    logical, optional, intent(in) :: map, nnonly
    character(len=*),optional,intent(in):: method


    real(dp):: ratio, evals(3), evecs(3,3), small_evec(3), r_ij, u_ij(3), cos_theta, oldratio
    real(dp), allocatable, dimension(:)::default_row
    real(dp), allocatable, dimension(:,:) :: df
    integer::Nfit_components, i, j, fi, fj, nj, n, hops, atomcount, oldn, &
         shift(3), k, best_j, best_fj, best_shift(3), sloc(1), sx, sy, sz
    type(Table) :: springs, newsprings, tmpsprings,cand_springs
    logical :: do_map, do_nnonly
    logical :: do_refit

    do_map = .false.
    if (present(map)) do_map = map

    do_nnonly = .true.
    if (present(nnonly)) do_nnonly = nnonly

    hops = 2
    if(present(spring_hops)) then
       hops = spring_hops
    end if

    ! If QM region hasn't changed since we last made out spring list, then
    ! there's nothing to be done now.

    if(do_map) then
       ! this save makes the twobody_old equal to twobody
       call adjustable_potential_save_parameters
       call print("Saving old parameters ....", PRINT_VERBOSE)
    end if

!!$    if (first_time) then
!!$       first_time = .false.
!!$    else if (fitlist%N == fitlist_old%N) then
!!$
!!$       allocate(s1(fitlist%N),s2(fitlist%N))
!!$       s1 = fitlist%int(1,1:fitlist%N)
!!$       s2 = atomlist%int(1,1:atomlist%N)
!!$       call sort_array(s1)
!!$       call sort_array(s2)
!!$
!!$       if (all(int_part(fitlist) == int_part(fitlist_old))) then
!!$          call print('fitlist unchanged, nothing to be done in adjustable_potential_init', PRINT_VERBOSE)
!!$          deallocate(s1)
!!$          deallocate(s2)
          
          ! Atom set matches, but it's possible some of the shifts have changed
          ! if any fit atoms have crossed a periodic boundary since last call to 
          ! adjustable_potential_init. Let's overwrite with new shifts. 
!!$          do n=1,twobody%N
!!$             i = atomlist%int(1,twobody%int(1,n))
!!$             j = atomlist%int(1,twobody%int(2,n))
!!$             ! overwrite old shift in twobody%int(3:5,n) with correct shift
!!$             r = distance_min_image(atoms_in, i, j, shift=shift)
!!$
!!$             if (abs(r-twobody%real(1,i)) < 0.1_dp) then
!!$                call print('shift update: '//shift//' '//twobody%int(3:5,i)//' '//r//' '//twobody%real(1,i))
!!$                twobody%int(3:5,i) = shift
!!$             end if
!!$
!!$             !! for the thin crack system, not all springs are min
!!$             !! image springs; some of them span periodic boundary.
!!$             !! how can we get round this?
!!$
!!$          end do
!!$
!!$          return
!!$       end if
!!$       
!!$       deallocate(s1)
!!$       deallocate(s2)
!!$    end if
!!$    call wipe(fitlist_old)
!!$    call append(fitlist_old, fitlist)

    call print("Setting up adjustable potential....", PRINT_VERBOSE)

    ! test table size
    if(fitlist%intsize /= 4 .and. fitlist%intsize /= 1)  &
         call system_abort("adjustable_potential_init: fitlist%intsize must be 1 or 4")
    if(fitlist%realsize /= 0) &
         call system_abort("adjustable_potential_init: fitlist%realsize /= 0")

    if (multiple_images(fitlist)) &
         call system_abort('adjustable_potential_init: fitlist contains repeated atom indices')

    call print("Fitting on atoms:", PRINT_VERBOSE)
    call print(fitlist, PRINT_VERBOSE)
    call print("", PRINT_VERBOSE)

    ! allocate globals
    call print("Allocating target_force", PRINT_NERD)
    Nfit_components = fitlist%N*3

    call reallocate(target_force,Nfit_components)
    target_force = 0.0_dp

    ! copy atom indices from argument to global
    
    call wipe(atomlist)
    call table_allocate(atomlist,1,0,0,0,fitlist%N)
    call append(atomlist, fitlist%int(1,1:fitlist%N))

    ! set up twobody table that holds pairs of atom indices and shiftts
    ! and the associated potential parameters (initialized to 0)

    ! wipe previous table
    call print("Wiping twobody table", PRINT_NERD)
    call wipe(twobody)
    call Table_Allocate(twobody, 5, 1+2*adjustable_potential_nparams,0,0)

    call table_allocate(springs, 4, 0, 0, 0)
    call table_allocate(newsprings, 4, 0, 0, 0)
    call table_allocate(tmpsprings, 4, 0, 0, 0)


    ! default parameters
    allocate(default_row(1+2*adjustable_potential_nparams))
    default_row = 0.0_dp


    ! First we try to add springs between pairs of atoms in fitlist within 'hops'
    ! bond hops of one another
    call print('Adding default springs')
    do fi=1,atomlist%N ! loop through fit atoms
       i = atomlist%int(1,fi)

       call wipe(newsprings)
       call append(newsprings, (/i,0,0,0/))
       call bfs_grow(atoms_in, newsprings, hops, nneighb_only = do_nnonly)

       do nj = 1, newsprings%N

          j = newsprings%int(1,nj)
          shift = newsprings%int(2:4,nj)

          fj = find(atomlist, j)
          if (fj == 0) cycle ! not in fitlist
          if (i == j) cycle  ! don't join to self or images of self

          ! Have we already got this spring?
          ! always store in order fi, fj with fi < fj  
          if (find(twobody, (/min(fi,fj), max(fi,fj), sign(1,fj-fi)*shift/)) /= 0) cycle

          ! Is this pair excluded?
          if (exclusions) then
             if (is_excluded(i,j,shift)) cycle
          end if

          default_row(1) = distance(atoms_in, i, j, shift)
          call append(twobody, (/min(fi,fj), max(fi,fj), sign(1,fj-fi)*shift/), default_row)
       end do

    end do
    write(line, '(a,i0,a)') 'Got ', twobody%N, ' springs' ; call print(line)
    write(line, '(a,f0.2)') 'Searching for difficult atoms, anisotropy threshold: ', AP_bad_atom_threshold
    call print(line, PRINT_VERBOSE)

    oldn = twobody%N		
    atomcount = 0
    ! For those fit atoms that we require directionality for, check how well the 
    ! springs connected to it span 3D space. If necessary add more springs for those atoms.
    do fi=1,directionN
       i = atomlist%int(1, fi)

       ! Find who atom i is connected to with springs, and get associated shifts
       call wipe(springs)
       do j = 1, twobody%N
          if (twobody%int(1,j) == fi) &
               call append(springs, (/atomlist%int(1,twobody%int(2,j)), twobody%int(3:5,j)/))
          if (twobody%int(2,j) == fi) &
               call append(springs, (/atomlist%int(1,twobody%int(1,j)),-twobody%int(3:5,j)/))
       end do

       ! How well do these springs span 3D space?
       if (springs%N > 2) then
          call directionality(atoms_in, i, springs, evals, evecs, method=2)
          ratio = minval(evals)/maxval(evals)
       else
          ratio = 0.0_dp ! less than three springs - really bad!
       end if

       sloc = minloc(evals)
       small_evec = evecs(:,sloc(1))

       if (ratio < AP_bad_atom_threshold) then

          write (line,'(a,i6,a,i6,a,f10.3)'), 'Atom ', i, ' with ', springs%N, &
               ' springs has ratio ', ratio
          call print(line, PRINT_VERBOSE)
          call print('Smallest evect ='//small_evec, PRINT_VERBOSE)
          
          ! Not well enough
          atomcount = atomcount + 1

          ! Tabulate candidate springs
          call allocate(cand_springs, 5, 2, 0, 0)
          do fj = 1,atomlist%N
             j = atomlist%int(1,fj)

             if (i == j) cycle ! Don't join to self

             do sx=-1,1
                do sy=-1,1
                   do sz=-1,1

                      shift = (/sx,sy,sz/)
                      
                      if (find(springs, (/j,shift/)) /= 0) cycle ! We've tried this spring already

                      ! Don't add if this pair is excluded
                      if (exclusions) then
                         if (is_excluded(i,j,shift)) cycle
                      end if

                      r_ij = distance(atoms_in,i,j,shift)
                      u_ij = diff(atoms_in,i,j,shift)
                      cos_theta = abs(small_evec .dot. u_ij) / (norm(small_evec)*norm(u_ij))

                      call append(cand_springs, (/j,fj,shift/), (/r_ij,cos_theta/))

                   end do
                end do
             end do
          end do

          ! Try to add springs until ratio exceeds thresholds
          n = 0
          do while (ratio < AP_bad_atom_threshold .and. n < 8)

             ! Find candidate spring with largest cos_theta
             sloc = maxloc(cand_springs%real(2,1:cand_springs%N))
             k = sloc(1)

             best_j     = cand_springs%int(1,k)
             best_fj    = cand_springs%int(2,k)
             best_shift = cand_springs%int(3:5,k)
             r_ij       = cand_springs%real(1,k)
             
             call print('Trying spring to atom '//best_j//' (index '//best_fj//') shift '//best_shift//' length '//r_ij)

             call delete(cand_springs, k)

             call append(springs, (/best_j,best_shift/))

             oldratio = ratio

             call directionality(atoms_in, i, springs, evals, evecs, 2)
             ratio = minval(evals)/maxval(evals)

             ! Add this spring to twobody only if it improves ratio
             if (ratio > oldratio) then
                default_row(1) = r_ij
                call append(twobody,(/min(fi,best_fj),max(fi,best_fj),sign(1,best_fj-fi)*best_shift/), default_row)
                
                write (line, '(a,i0,a,f0.3)') '  Now got ', springs%N,' springs. Ratio improved to ', ratio
                call print(line, PRINT_VERBOSE)
             end if

             n = n + 1

!!$     newsprings = springs
!!$
!!$          call bfs_step(atoms_in, newsprings, tmpsprings, nneighb_only = .true.)
!!$          call append(newsprings, tmpsprings)
!!$
!!$          if (tmpsprings%N == 0) then
!!$             ! Didn't manage to add any atoms, so try nneigb_only = false, since
!!$             ! these are the 'bad' atoms and their neighbours are 
!!$             ! probably some distance away
!!$
!!$             call bfs_step(atoms_in, newsprings, tmpsprings, nneighb_only = .false.)
!!$             call append(newsprings, tmpsprings)
!!$          end if
!!$
!!$          do nj = 1, newsprings%N
!!$             j = newsprings%int(1,nj)
!!$             shift = newsprings%int(2:4,nj)
!!$             fj = find(atomlist, j)
!!$             if (fj == 0) cycle ! not in fitlist
!!$             if (i == j) cycle ! don't join to self or images of self
!!$
!!$             ! already got this spring
!!$             if (find(twobody, (/min(fi,fj), max(fi,fj), sign(1,fj-fi)*shift/)) /= 0) cycle
!!$
!!$             ! Don't add if this pair is excluded
!!$             if (exclusions) then
!!$                if (is_excluded(i,j,shift)) cycle
!!$             end if
!!$
!!$             default_row(1) = distance(atoms_in, i, j, shift)
!!$             call append(twobody,(/min(fi,fj),max(fi,fj),sign(1,fj-fi)*shift/), default_row)
!!$          end do
!!$
!!$          ! Work out the springs again
!!$          call wipe(springs)
!!$          do j = 1, twobody%N
!!$             if (twobody%int(1,j) == fi) &
!!$                  call append(springs, (/atomlist%int(1,twobody%int(2,j)), twobody%int(3:5,j)/))
!!$             if (twobody%int(2,j) == fi) &
!!$                  call append(springs, (/atomlist%int(1,twobody%int(1,j)),-twobody%int(3:5,j)/))
!!$          end do
!!$
!!$          ! Did we fail to add any springs at all? If so system has probably exploded
!!$          if (springs%N == 0) then
!!$             write(line, '(a,i0)') 'adjustable_potential_init: can''t add any springs for atom ', i
!!$             call system_abort(line)
!!$          end if
!!$
!!$          call directionality(atoms_in, i, springs, evals, evecs, 2)
!!$          ratio = minval(evals)/maxval(evals)
!!$
!!$          write (line, '(a,i0,a,f0.3)') '  Now got ', springs%N,' springs. Ratio improved to ', ratio
!!$          call print(line, PRINT_VERBOSE)
!!$
!!$          n = n + 1
          end do
          if (n == 8 .and. ratio < AP_bad_atom_threshold) then
             call system_abort('adjustable_potential_init: failed to sufficiently improve bad spring space spanning')
          end if
          
       end if

    end do

    if(atomcount > 0) then
       write(line, '(a,i0,a,i0,a)') 'Added ', twobody%N-oldn, ' springs to fix ', atomcount, ' atoms' 
       call print(line)
    end if

    call finalise(springs)
    call finalise(newsprings)
    call finalise(tmpsprings)
    deallocate(default_row)

    call print("Allocating twobody table to proper size", PRINT_NERD)
    call table_allocate(twobody) ! reduce table to proper length

    write(line, '(a,i0)') "Number of force components: ", size(target_force) ; call print(line)
    write(line, '(a,i0)') "Number of parameters:       ", twobody%N*adjustable_potential_nparams ; call print(line)

    adjustable_potential_initialised = .true.

    if(do_map) then
       ! this maps from the twobody_old with the old structure into twobody with the new structure
       call adjustable_potential_map_parameters

       ! Only need to refit if spline set has changed since last time
       do_refit = .false.
       ! if sizes don't match definitely need to refit
       if (twobody%N /= twobody_old%N) then 
          do_refit = .true.
       else
          do i=1,twobody%N
             ! As soon as we can't find a new spring in twobody_old, know we have to refit
             if (find(twobody_old,twobody%int(:,i)) == 0) then
                do_refit = .true.
                exit
             end if
          end do
       end if

       ! reoptimise with new springs to reproduce old forces at new positions
       ! if fit region has changed, just copy the target force components that
       ! are in both the old and new regions
       if (do_refit) then
          call print('Springs changed, refitting old forces with new spring set')
          allocate(df(3,atomlist%N))
          df = 0.0_dp 
          do i=1,atomlist_old%N
             fi = find(atomlist, atomlist_old%int(:,i)) ! find old fit atom in new fitlist
             if (fi /= 0) df(:,fi) = fitforce_old(:,i)         ! if found, copy target force
          end do
          call adjustable_potential_optimise(atoms_in, df, method=method)
          deallocate(df)
       end if

    endif

    ! this save makes the twobody_old equal to twobody (both new structure)
    call adjustable_potential_save_parameters

  end subroutine adjustable_potential_init

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X adjustable_potential_finalise()
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  subroutine adjustable_potential_finalise()
    if(allocated(target_force)) deallocate(target_force)
    call finalise(atomlist)
    call finalise(atomlist_old)
    call finalise(twobody)
    call finalise(twobody_old)
    call adjustable_potential_delete_saved_parameters
    call sparse_finalise(forcematrix)
    if (exclusions) call finalise(exclusion_list)
  end subroutine adjustable_potential_finalise

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X adjustable_potential_delete_saved_parameters
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine adjustable_potential_delete_saved_parameters
    call finalise(twobody_old)
    call finalise(atomlist_old)
    adjustable_potential_parameters_saved = .false.
  end subroutine adjustable_potential_delete_saved_parameters

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X adjustable_potential_save_parameters
  !X
  !X Makes a copy of the twobody and atomlist data, avoiding a
  !X reallocation if possible
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine adjustable_potential_save_parameters
    if (.not.adjustable_potential_initialised) &
         call system_abort('Adjustable_Potential_Save_Parameters: No parameters exist to save')
    call wipe(twobody_old)
    call append(twobody_old,twobody)

    call wipe(atomlist_old)
    call append(atomlist_old, atomlist)

    adjustable_potential_parameters_saved = .true.

  end subroutine adjustable_potential_save_parameters

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X adjustable_potential_map_parameters
  !X
  !X Given an old set of data (twobody_old, atomlist_old), fill in the
  !X y1 and y2 values which were used previously
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine adjustable_potential_map_parameters
    integer :: i, j, fi, fj, n, m, k, shift(3)

    if(.not.adjustable_potential_initialised) &
         call system_abort('Adjustable_Potential_Map_Parameters: No parameters to map!')
    if (twobody_old%N == 0) &
         call system_abort('Adjustable_Potential_Map_Parameters: No parameters are saved')
    
    if (twobody%N >= twobody_old%N) then  ! loop over the smallest table

       !Loop over the entries in the old table
       do n = 1, twobody_old%N
          
          !look up the atom indices
          i     = atomlist_old%int(1,twobody_old%int(1,n))
          j     = atomlist_old%int(1,twobody_old%int(2,n))
          shift = twobody_old%int(3:5,n)
          
          !try to find the corresponding atoms in the new list
          fi = find(atomlist,i)
          if (fi == 0) cycle !atom i is not present in the new list
          fj = find(atomlist,j)
          if (fj == 0) cycle !atom j is not present in the new list
          
          !find the entry in the new twobody table
          m = find(twobody,(/min(fi,fj),max(fi,fj),sign(1,fj-fi)*shift/))
          if (m == 0) cycle ! i and j are no longer connected
          
          ! copy the parameter values across
          do k=1,2*adjustable_potential_nparams
             twobody%real(1+k,m) = twobody_old%real(1+k,n)
          end do
          
       end do

    else

       !Loop over the entries in the new table
       do n = 1, twobody%N
          
          !look up the atom indices
          i     = atomlist%int(1,twobody%int(1,n))
          j     = atomlist%int(1,twobody%int(2,n))
          shift = twobody%int(3:5,n)
          
          !try to find the corresponding atoms in the old list
          fi = find(atomlist_old, i)
          if (fi == 0) cycle !atom i is not present in the old list
          fj = find(atomlist_old, j)
          if (fj == 0) cycle !atom j is not present in the old list
          
          !find the entry in the old twobody table
          m = find(twobody_old,(/min(fi,fj),max(fi,fj),sign(1,fj-fi)*shift/))
          if (m == 0) cycle ! i and j were not neighbours
          
          !copy the param values across
          do k=1,2*adjustable_potential_nparams
             twobody%real(1+k,n) = twobody_old%real(1+k,m)
          end do

       end do

    end if

  end subroutine adjustable_potential_map_parameters


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X adjustable_potential_force(atoms, force)
  !X
  !X force and energy from the adjustable potential, using the interatomic
  !X distance information from atoms. If interp is given, then interpolate
  !X the parameters between twobody_old and twobody
  !X
  !X the force(:) array is the same size as the stored atomlist.
  !X 
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine adjustable_potential_force(atoms_in, force, interp, interp_space, interp_order, energy, power)
    type(Atoms),                            intent(in)  :: atoms_in
    real(dp), dimension(3,atomlist%N), intent(out) :: force
    real(dp), optional,                     intent(in)  :: interp
    logical, optional,                      intent(in)  :: interp_space
    character(len=*),optional,              intent(in)  :: interp_order
    real(dp), optional, intent(out) :: energy, power

    ! local
    integer:: i, a1, a2, fi, fj, shift(3), my_interp_order
    real(dp), dimension(3)::f
    real(dp)::fij, r, interp_point, dalpha
    logical::interpolate, do_interp_space   
    
    do_interp_space = optional_default(.false., interp_space)

    interpolate = present(interp) .or. do_interp_space
    if(interpolate .and. .not.adjustable_potential_parameters_saved) &
         call System_Abort('Adjustable_Potential_Force: No parameters saved to interpolate with')

    my_interp_order = 1
    if(present(interp_order)) then
       if(trim(interp_order) .eq. 'linear') then
          my_interp_order = 1
       else if(trim(interp_order) .eq. 'quadratic') then
          my_interp_order = 2
       else if(trim(interp_order) .eq. 'cubic') then
          my_interp_order = 3
       else
          call System_Abort("Invalid interpolation order `"//interp_order//"'in Adjustable_Potential_Force")
       end if
    end if

    force = 0.0_dp

    if(present(energy)) energy = 0.0_dp
    if(present(power))   power = 0.0_dp

    do i=1,twobody%N       

       fi = twobody%int(1,i)
       fj = twobody%int(2,i)
       shift = twobody%int(3:5,i)

       a1 = atomlist%int(1,fi)
       a2 = atomlist%int(1,fj)
       r = distance(atoms_in, a1, a2, shift)

       if(interpolate) then
          if(do_interp_space) then
             interp_point = (r - twobody_old%real(1,i))/(twobody%real(1,i)-twobody_old%real(1,i))
             fij = linear_interpolate(0.0_dp,twobody_old%real(2,i),&   ! x0,y0
               1.0_dp,twobody%real(2,i),interp_point) !x1, y1, x
          else
             select case(my_interp_order)
                case(1) ! linear interpolation
                   fij = linear_interpolate(0.0_dp,twobody_old%real(2,i),&   ! x0,y0
                        1.0_dp,twobody%real(2,i),interp) !x1, y1, x

                   dalpha = (twobody%real(2,i) - twobody_old%real(2,i))
                case(2) ! quadratic using the derivative info at start point
                   fij = twobody_old%real(2,i)+interp*twobody_old%real(3,i)&
                        +interp*interp*(twobody%real(2,i)-twobody_old%real(2,i)-twobody_old%real(3,i))

                   dalpha = 0.0_dp
                case(3) ! cubic with zero derivatives at endpoints
                   fij = cubic_interpolate(0.0_dp,twobody_old%real(2,i),&   ! x0,y0
                        1.0_dp,twobody%real(2,i),interp) ! x1,y1, x
                   
                   dalpha = 0.0_dp                
                case default
                   call System_Abort("Invalid interpolation order in Adjustable_Potential_Force")
             end select
          end if
       else
          ! Extrapolate
          fij = twobody%real(2,i)
          dalpha = 0.0_dp ! Parameters don't change during extrap
       end if

       ! optionally calculate energy and power
       if(present(energy)) energy = energy + r*fij
       if(present(power))  power  = power  + r*dalpha ! (since dV_dalpha = r, and dt = 1 here)

       if (value(mainlog%verbosity_stack) >= PRINT_NERD) &
            call print(i//' '//a1//' '//a2//' '//r//' '//fij//' <-- SPRING', PRINT_NERD)

       ! get force
       f = direction_cosines(atoms_in,a1,a2,shift)*fij
       force(:,fi) = force(:,fi) + f
       force(:,fj) = force(:,fj) - f
    end do

  end subroutine adjustable_potential_force



  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X adjustable_potential_force_sparse(params)
  !X
  !X force from the adjustable potential, using the precomputed
  !X sparse matrix
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! 
  function adjustable_potential_force_sparse(params) result(force) 
    real(dp)::params(:)
    real(dp)::force(size(target_force))

    if(.not.adjustable_potential_initialised) & 
         call system_abort("adjustable_potential_force_sparse: not initialised!")

    force = (forcematrix .mult. params) 
  end function adjustable_potential_force_sparse

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X adjustable_potential_force_error(params)
  !X
  !X normsq() of difference between target and adjustable force 
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  function adjustable_potential_force_error(params,data) result(error)
    real(dp)::params(:)
    character,optional :: data(:)
    real(dp)::error

    if(.not.adjustable_potential_initialised) & 
         call system_abort("adjustable_potential_force_error: not initialised!")

    if(value(mainlog%verbosity_stack) .ge. PRINT_NERD) then
       write(line, *) "adjustable_potential_force_error: params(", size(params), ")"
       call print(line, PRINT_NERD)
       call print(params, PRINT_NERD)
    end if

    error = normsq(target_force - (forcematrix .mult. params))

  end function adjustable_potential_force_error

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X adjustable_potential_force_error_deriv(params)
  !X
  !X derivative of the force_error wrt variable parameters
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! 
  function adjustable_potential_force_error_deriv(params,data) result(deriv)
    real(dp)::params(:)
    character,optional :: data(:)
    real(dp)::force(size(target_force))
    real(dp)::deriv(size(params))

    if(.not.adjustable_potential_initialised) & 
         call system_abort("adjustable_potential_force_error_deriv: not initialised!")

    if(value(mainlog%verbosity_stack) .ge. PRINT_NERD) then
       write(line, *) "adjustable_potential_force_error_deriv: params(", size(params), ")"
       call print(line, PRINT_NERD)
       call print(params, PRINT_NERD)
    end if

    force = (forcematrix .mult. params)
    deriv = (-2.0_dp * (target_force - force)) .mult. forcematrix
  end function adjustable_potential_force_error_deriv

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X adjustable_potential_create_forcematrix()
  !X
  !X creates the sparse matrix forcematrix that we
  !X use to compute the adjustable forces quickly for different
  !X parameters, but unchanging atomic positions
  !X
  !X forcematrix has as many rows as the 3*number of atoms to be fitted
  !X on, so same as size(target_force), and as many columns as the
  !X number of adjustable parameters, so twobody%N*adjustable_potential_nparams
  !X 
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !   

  subroutine adjustable_potential_create_forcematrix(atoms_in, precondition)
    type(Atoms), intent(IN)::atoms_in
    real(dp), intent(in) :: precondition(:)

    ! temporary hack to go via a full matrix. to be redone
    real(dp), allocatable, dimension(:,:)::tmpmatrix
    integer, parameter::nparams = adjustable_potential_nparams
    integer::i, j, k, a1, a2, fi, fj, shift(3)
    real(dp)::cosi(3)
    integer, allocatable::colcounter(:), tmpcolumns(:,:)


    write(line, *) "Setting up adjustable potential forcematrix.." ; call print(line, PRINT_VERBOSE)
    write(line, *) ; call print(line, PRINT_VERBOSE)

    allocate(tmpmatrix(size(target_force), twobody%N*nparams))
    allocate(colcounter(size(target_force)))
    allocate(tmpcolumns(size(target_force), twobody%N*nparams))
    tmpmatrix = 0.0_dp
    colcounter = 0
    tmpcolumns = 0

    ! loop through the springs
    do i=1,twobody%N       
       ! we take and cache the distance from atoms_in, as it might have changed since time
       ! we initialised the twobody table
       fi = twobody%int(1,i)
       fj = twobody%int(2,i)
       shift = twobody%int(3:5,i)
       a1 = atomlist%int(1,fi)
       a2 = atomlist%int(1,fj)

       twobody%real(1,i) = distance(atoms_in, a1, a2, shift)
       cosi = direction_cosines(atoms_in, a1, a2, shift)

       fi = (twobody%int(1,i)-1)  ! subtract 1 so we can use it as index
       do k=1,3
          colcounter(fi*3+k) = colcounter(fi*3+k) + 1
          ! Multiply by preconditioning factor for this spring
          tmpmatrix (fi*3+k, colcounter(fi*3+k)) = cosi(k)*precondition(i)
          tmpcolumns(fi*3+k, colcounter(fi*3+k)) = (i-1)*nparams+1
       end do
       fj = (twobody%int(2,i)-1) ! subtract 1 so we can use it as index
       do k=1,3
          colcounter(fj*3+k) = colcounter(fj*3+k) + 1
          ! Multiply by preconditioning factor for this spring
          tmpmatrix (fj*3+k, colcounter(fj*3+k)) = -cosi(k)*precondition(i)
          tmpcolumns(fj*3+k, colcounter(fj*3+k)) = (i-1)*nparams+1
       end do
    end do ! i=1,twobody%N


    ! fill in sparse matrix from values in tmpmatrix and tmpcolumns
    call sparse_init(forcematrix, size(target_force), twobody%N*nparams)
    do i=1,size(target_force)
       forcematrix%rows(i) = forcematrix%table%N+1
       do j=1,colcounter(i)
          call append(forcematrix%table, tmpcolumns(i,j), tmpmatrix(i,j))
       end do
    end do
    ! put in last row value
    forcematrix%rows(i) = forcematrix%table%N+1

    deallocate(tmpmatrix)
    deallocate(tmpcolumns)
    deallocate(colcounter)

    ! print it
    if(value(mainlog%verbosity_stack) .ge. PRINT_NERD) then
       call print('Forcematrix:',PRINT_NERD)
       call print(forcematrix, PRINT_NERD)
       call print('',PRINT_NERD)
    end if
  end subroutine adjustable_potential_create_forcematrix

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X adjustable_potential_optimise()
  !X
  !X the whole point of this module: vary the adjustable parameters
  !X to minimize the force_error() 
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine adjustable_potential_optimise(atoms_in, fitforce, fromscratch, method, prec_func, parallel_svd)
    type(Atoms),              intent(in):: atoms_in
    real(dp), dimension(:,:), intent(in):: fitforce
    logical, optional                   :: fromscratch

    logical                             :: do_fromscratch
    character(len=*),optional,intent(in):: method
    logical, optional, intent(in)       :: parallel_svd
    logical::do_svd, do_parallel_svd
    integer, parameter::nparam = adjustable_potential_nparams
    integer::nparamtot, bad_spring
    real(dp), allocatable::params(:)
    integer::steps, i, shift(3), k
    real(dp):: err, ferr(3, atomlist%N)
    real(dp), parameter::eigenvalue_threshold=1e-10_dp
    real(dp), allocatable::full_forcematrix(:,:), X(:,:), S(:), WORK(:), precondition(:)
    integer::LWORK, INFO, maxmn, minmn, M,N, RANK, a1, a2, error_code
    optional :: prec_func

#ifdef _MPI
    include 'mpif.h'
#endif

    interface 
       function prec_func(r)
         use System_module
         real(dp), intent(in) :: r
         real(dp) :: prec_func
       end function prec_func
    end interface

    if(twobody%N == 0) then
       call print("There are no parameters to optimise")
       return
    end if
    
    do_fromscratch = .true.
    if(present(fromscratch)) do_fromscratch = fromscratch

    do_svd = .true.
    if(present(method)) then
       if(trim(method) .eq. "SVD") then
          do_svd = .true.
       else if(trim(method) .eq. "minim") then
          do_svd = .false.
       else
          call system_abort("Adjustable_potential_optimise: invalid method '" // trim(method) // "'")
       end if
    end if
    
    !If the chosen method is SVD then the default is to run serially to
    !avoid memory bottlenecks when all processors do SVD.
    do_parallel_svd = optional_default(.false.,parallel_svd)

    if(size(target_force) /= size(fitforce)) &
         call system_abort("adjustable_potential_optimise: target force is wrong size")
    
    ! Calculate preconditioning factors for each spring
    allocate(precondition(twobody%N))
    if (present(prec_func)) then
       do i=1,twobody%N
          a1 = atomlist%int(1,twobody%int(1,i))
          a2 = atomlist%int(1,twobody%int(2,i))
          shift = twobody%int(3:5,i)
          precondition(i) = prec_func(distance(atoms_in, a1, a2, shift))
       end do
    else
       precondition = 1.0_dp
    end if

    ! create matrices for fast evaluation
    call adjustable_potential_create_forcematrix(atoms_in, precondition)
    

    target_force = reshape(fitforce, (/size(target_force)/))

    nparamtot = twobody%N*nparam
    allocate(params(nparamtot))

    if(.not.adjustable_potential_initialised) & 
         call system_abort("adjustable_potential_optimise: not initialised!")

    call print("Target force:", PRINT_NERD)
    call print(target_force, PRINT_NERD)

    call print('Optimising '//nparamtot//' adjustable parameters')
    call print("RMS force component error before optimisation : "// &
         sqrt(normsq(target_force)/(real(atomlist%N,dp)*3.0_dp)))
    call print("Max force component error before optimisation : "//&
         maxval(abs(target_force)))

    if(do_fromscratch) then
       params = 0.0_dp
    else
       params = reshape(twobody%real(2:1+nparam,1:twobody%N), (/nparamtot/))
    endif

    if(do_svd) then

#ifdef _MPI

       if (do_parallel_svd .or. mpi_id()==0) then

#endif

       ! Solution by SVD
       call print('Using SVD for least squares fit, eigenvalue threshold = '//eigenvalue_threshold)
       allocate(full_forcematrix(size(target_force), size(params)))
       full_forcematrix = forcematrix
       M = size(full_forcematrix,1)
       N = size(full_forcematrix,2)
       minmn = minval((/M,N/))
       maxmn = maxval((/M,N/))
       allocate(X(maxmn,1))
       allocate(S(minmn))

       allocate(WORK(1))

       X(1:M,1) = target_force

       ! Do a workspace query first
       call dgelss(M, N, 1, full_forcematrix, M, X,maxmn,S, eigenvalue_threshold, RANK, WORK, -1, INFO )

       ! Allocate optimal size workspace
       LWORK = WORK(1)
       deallocate(WORK)
       allocate(WORK(LWORK))

       ! Do the SVD
       call dgelss(M, N, 1, full_forcematrix, M, X,maxmn,S, eigenvalue_threshold, RANK, WORK, LWORK, INFO )

       if (INFO /= 0) then
          call system_abort('adjustable_potential_optimise: DGELSS failed with error code '//INFO)
       end if

       params = X(1:N,1)
       deallocate(WORK)
       deallocate(S)
       deallocate(X)
       deallocate(full_forcematrix)
       
#ifdef _MPI
       
       endif

       !Broadcast the params to each process if SVD only done on proc 0
       if (.not.do_parallel_svd) then
          call MPI_BCAST(params,size(params),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,error_code)
          !check for errors
          call abort_on_mpi_error(error_code, "Adjustable_Potentia_Optimise: Params MPI_Bcast()")
       end if

#endif

    else
       !! Solution by minimisation
       steps = minim(params,&
            adjustable_potential_force_error, &
            adjustable_potential_force_error_deriv,&
            "cg",1e-10_dp,1000,"FAST_LINMIN")
    end if
       

    ! compute RMS after optimisation
    err = sqrt(adjustable_potential_force_error(params)/(real(atomlist%N,dp)*3.0_dp))

    ! Print all magnitudes of force errors if we're a NERD
    call print('Force errors', PRINT_NERD)
    ferr = reshape(target_force - (forcematrix .mult. params), &
         (/3, atomlist%N/))

!!$    call print('sparse force=')
!!$    call Print(forcematrix .mult. params)
!!$    call print('')
    
    do i=1,atomlist%N
       write (line, '(i4,f12.6)') i, sqrt(normsq(ferr(:,i)))
       call print(line, PRINT_NERD)
    end do
    call Print('', PRINT_NERD)

    write(line, '(a,e10.2)') "RMS force component error after  optimisation : ", err
    call print(line)
    write(line, '(a,e10.2)') "Max force component error after  optimisation : ", &
         maxval(abs(target_force - (forcematrix .mult. params)))
    call print(line)

    if (value(mainlog%verbosity_stack) >= PRINT_NERD) then
       call Print('Final Springs:', PRINT_NERD)
       do i=1,twobody%N
          write (line, '(3i8,4f16.8)') i, atomlist%int(1,twobody%int(1,i)), atomlist%int(1,twobody%int(2,i)), twobody%real(1,i), precondition(i), params(i), params(i)*precondition(i)
          call Print(line, PRINT_NERD)
       end do
       call Print('', PRINT_NERD)
    end if

    ! Get physical spring constants by multiplying by precondition factor
    do i=1,twobody%N
       params(i) = params(i)*precondition(i)
    end do

    twobody%real(2:1+nparam,1:twobody%N) = reshape(params, (/nparam,twobody%N/))

    write(line, '(a,e10.2)') "Max abs spring constant   after  optimisation : ", maxval(abs(params))
    call print(line)

    if (maxval(abs(params)) > adjustable_potential_max_spring_constant) then
       write(line, '(a,f10.2)') "WARNING: Max spring constant value after optimisation is large : ", &
            maxval(abs(params))
       call print(line)
       bad_spring = maxloc(abs(params),1)
       write(line,'(2(a,i0))')'         Max spring constant belongs to spring ', &
                               atomlist%int(1,twobody%int(1,bad_spring)),'--',atomlist%int(1,twobody%int(2,bad_spring))
       call print(line)
    end if

    call print('')

    deallocate(params)
    deallocate(precondition)

    ! calculate new parameter derivatives
    do i=1,min(twobody%N,twobody_old%N)
       do k=1,nparam
          twobody%real(1+nparam+k,i) = twobody_old%real(1+nparam+k, i)+2.0_dp*(twobody%real(1+k,i)-twobody_old%real(1+k,i)-twobody_old%real(1+nparam+k,i))
       end do
    end do

    ! Save the old fit force
    if (allocated(fitforce_old)) deallocate(fitforce_old)
    allocate(fitforce_old(size(fitforce,1),size(fitforce,2)))
    fitforce_old = fitforce

  end subroutine adjustable_potential_optimise

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X adjustable_potential_print_params()
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine adjustable_potential_print_params()
    real(dp)::f(size(target_force))
    integer::np

    if(.not.adjustable_potential_initialised) & 
         call system_abort("adjustable_potential_print_params: not initialised!")


    write(line,*) "Atoms              Shift                   Parameters" ; call print(line)
    write(line,*) "-----------------------------------------------------" 
    call print(line)
    call print(twobody)
    write(line, *) ; call print(line)


    write(line, *) "Force matrix: " ; call print(line)
    write(line, *) ; call print(line)
    call print_full(forcematrix)
    write(line, *) ; call print(line)

    write(line,*) "Forces from force matrix:" ; call print( line)
    np = twobody%N*adjustable_potential_nparams
    f = adjustable_potential_force_sparse(reshape(twobody%real(2:1+adjustable_potential_nparams,1:twobody%N), (/np/)))
    call print(transpose(reshape(f, (/3,atomlist%N/))))
    write(line, *) ; call print(line)

  end subroutine adjustable_potential_print_params

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X Set up, alter and query the exclusion list
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine adjustable_potential_add_exclusion(atom1,atom2,shift)

    integer, intent(in) :: atom1,atom2,shift(3)

    if (atom1 == atom2) then
       call print_warning('Adjustable_Potential_Add_Exclusion: Atom1 = Atom2')
       return
    end if

    if (.not.exclusions) then !Set up the exclusion table
       call allocate(exclusion_list,5,0,0,0)
       exclusions = .true.
    end if

    !Append (/a, b, shift/) where a < b
    if (atom1 < atom2) then
       call append(exclusion_list,(/atom1,atom2,shift/))
    else
       call append(exclusion_list,(/atom2,atom1,-shift/))
    end if

  end subroutine adjustable_potential_add_exclusion

  subroutine adjustable_potential_delete_exclusion(atom1,atom2,shift)

    integer, intent(in) :: atom1, atom2, shift(3)
    
    if (atom1 < atom2) then
       call delete(exclusion_list,(/atom1,atom2,shift/))
    else
       call delete(exclusion_list,(/atom2,atom1,-shift/))
    end if

    if (exclusion_list%N == 0) then
       call finalise(exclusion_list)
       exclusions = .false.
    end if

  end subroutine adjustable_potential_delete_exclusion

  function is_excluded(atom1,atom2, shift)

    integer, intent(in) :: atom1, atom2, shift(3)
    logical             :: Is_Excluded

    if (atom1 < atom2) then
       if (find(exclusion_list,(/atom1,atom2,shift/))==0) then
          Is_Excluded = .false.
       else
          Is_Excluded = .true.
       end if
    else
       if (find(exclusion_list,(/atom2,atom1,-shift/))==0) then
          Is_Excluded = .false.
       else
          Is_Excluded = .true.
       end if
       
    end if

  end function is_excluded

end module AdjustablePotential_module

