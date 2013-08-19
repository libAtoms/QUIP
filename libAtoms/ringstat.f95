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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X Ring statistics module
!X
!% Functions to compute ring statistics according to
!% D.S. Franzblau, Phys. Rev. B 44, 4925 (1991)
!%
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


#include "error.inc"

module ringstat_module
  use error_module
  use system_module, only : dp, operator(//)
  use linearalgebra_module
  use atoms_module

  implicit none
  private

  public :: count_sp_rings, distance_map

contains

  !% Look for the shortest distance starting at atom *root*,
  !% look only for atoms where *mask* is true
  subroutine find_shortest_distances(at, root, dist, mask, error)
    implicit none

    type(Atoms), intent(in)           :: at
    integer, intent(in)               :: root
    integer, intent(out)              :: dist(at%N)
    logical, intent(in), optional     :: mask(at%N)
    integer, intent(out), optional  :: error

    ! ---

    integer, parameter  :: MAX_WALKER = 8192

    ! ---

    integer   :: i, ni, j, k
    
    integer   :: n_walker
    integer   :: walker(MAX_WALKER)
    
    integer   :: n_new_walker
    integer   :: new_walker(MAX_WALKER)

    ! ---

    INIT_ERROR(error)

    n_walker     = 1
    walker(1)    = root

    do while (n_walker > 0)

       n_new_walker  = 0

       do k = 1, n_walker
          i  = walker(k)

!          do ni = nl%seed(i), nl%last(i)
          do ni = 1, n_neighbours(at, i)
!             j  = nl%neighbors(ni)
             j = neighbour(at, i, ni)
             if ( ( .not. present(mask) .or. mask(j) ) .and. dist(j) == 0 ) then

                n_new_walker              = n_new_walker+1
                if (n_new_walker > MAX_WALKER) then
                   RAISE_ERROR("MAX_WALKER ("//MAX_WALKER//") exceeded.", error)
                endif

                new_walker(n_new_walker)  = j

                dist(j)                   = dist(i)+1

             endif

          enddo

       enddo

       n_walker            = n_new_walker
       walker(1:n_walker)  = new_walker(1:n_walker)

    enddo

    dist(root)  = 0

  endsubroutine find_shortest_distances


  !% Look for the shortest distance starting at atom *root*,
  !% look only for atoms where *mask* is true
  subroutine distance_map(at, dist, mask, diameter, error)
    implicit none

    type(Atoms), intent(in)           :: at
    integer, intent(out)              :: dist(at%N, at%N)
    logical, intent(in), optional     :: mask(at%N)
    integer, intent(out), optional    :: diameter
    integer, intent(out), optional  :: error

    ! ---

    integer   :: i
    
    ! ---

    INIT_ERROR(error)

    if (present(diameter)) then
       diameter  = 0
    endif

    dist  = 0

    do i = 1, at%N
       if ( .not. present(mask) .or. mask(i) ) then

          call find_shortest_distances(at, i, dist(:, i), mask, error)
          PASS_ERROR(error)

          if (present(diameter)) then
             diameter = max(diameter, maxval(dist(:, i)))
          endif

       endif
    enddo

    ! Check that matrix is symmetric
    if (any(dist /= transpose(dist))) then
       RAISE_ERROR("*dist* matrix not symmetric.", error)
    endif

  endsubroutine distance_map


  !% Look for the "shortest path" ring starting at atom *at*,
  !% look only for atoms where *mask* is true
  subroutine count_sp_rings(at, cutoff, dist, max_ring_len, stat, mask, rings_out, dr_out, error)
    implicit none

    type(Atoms), intent(in)                     :: at
    real(DP), intent(in)                        :: cutoff
    integer, intent(in)                         :: dist(at%N, at%N)
    integer, intent(in)                         :: max_ring_len
    integer, intent(inout)                      :: stat(max_ring_len)
    logical, intent(in), optional               :: mask(at%N)
    integer, intent(out), optional              :: rings_out(:, :)
    real(dp), intent(out), optional             :: dr_out(:, :, :)
    integer, intent(out), optional            :: error

    ! ---

    integer, parameter   :: MAX_WALKER  = 16384

    real(DP), parameter  :: EPS  = 0.0001_DP

    ! ---

    integer   :: i, ni, j, k, m, n, a, na, b, dn

    real(DP)  :: d(3), abs_dr_sq, cutoff_sq

    logical, allocatable  :: done(:)
    logical   :: is_sp

    integer   :: n_walker
    integer   :: walker(MAX_WALKER)
    integer   :: last(MAX_WALKER)

    integer   :: ring_len(MAX_WALKER)
    integer   :: rings(max_ring_len, MAX_WALKER)
    real(DP)  :: dr(3, max_ring_len, MAX_WALKER)

    integer   :: n_new_walker
    integer   :: new_walker(MAX_WALKER)
    integer   :: new_last(MAX_WALKER)

    integer   :: new_ring_len(MAX_WALKER)
    integer   :: new_rings(max_ring_len, MAX_WALKER)
    real(DP)  :: new_dr(3, max_ring_len, MAX_WALKER)

    integer   :: no(at%N), n_rings_out

    ! ---

!    allocate(done(at%connect%neighbour1%t%N+at%connect%neighbour2%t%N))

    INIT_ERROR(error)

    if (present(rings_out)) then
       if (size(rings_out,1) < max_ring_len+1) then
          RAISE_ERROR('size(rings_out,1) < max_ring_len+1', error)
       end if
       if (.not. present(dr_out)) then
          RAISE_ERROR('rings_out is present but dr_out is not', error)
       end if
       rings_out(:,:) = 0
       dr_out(:,:,:) = 0.0_dp
       n_rings_out = 0
    end if

    cutoff_sq  = cutoff**2

    no(1) = 0
    do i = 2, at%N
       no(i) = no(i-1) + n_neighbours(at, i-1)
    enddo

    allocate(done(no(at%N) + n_neighbours(at, at%N)))

    done  = .false.

    bonds_loop1: do a = 1, at%N
       if ( .not. present(mask) .or. mask(a) ) then
!          bonds_loop2: do na = nl%seed(a), nl%last(a)
          bonds_loop2: do na = 1, n_neighbours(at, a)
!             b  = nl%neighbors(na)
             b = neighbour(at, a, na)

             if ( a < b .and. ( .not. present(mask) .or. mask(b) ) ) then

                done(no(a)+na)  = .true.
!                do ni = nl%seed(b), nl%last(b)
!                   if (nl%neighbors(ni) == a) then
!                      done(ni)  = .true.
!                   endif
!                enddo
                do ni = 1, n_neighbours(at, b)
                   if (neighbour(at, b, ni) == a) then
                      done(no(b)+ni)  = .true.
                   endif
                enddo

                n_walker     = 1
                walker(1)    = b
                last(1)      = a

                ring_len     = 1
                rings(1, 1)  = b
!                dr(:, 1, 1)  = GET_DRJ(p, nl, a, b, na)
                b = neighbour(at, a, na, diff=dr(:, 1, 1))

                while_walkers_exist: do while (n_walker > 0)

                   n_new_walker  = 0

                   do k = 1, n_walker
                      i  = walker(k)

                      if ( i > 0 ) then

!                         ni_away_from_root_loop: do ni = nl%seed(i), nl%last(i)
                         ni_away_from_root_loop: do ni = 1, n_neighbours(at, i)
                            ni_away_from_root_if: if ( .not. done(no(i)+ni) ) then

!                               DISTJ_SQ(p, nl, i, ni, j, d, abs_dr_sq)
                               j = neighbour(at, i, ni, diff=d)
                               abs_dr_sq = d .dot. d

                               if ( abs_dr_sq < cutoff_sq .and. ( .not. present(mask) .or. mask(j) ) .and. j /= last(k) ) then

                                  if ( dist(j, a) == dist(i, a)+1 ) then

                                     if ( ring_len(k) < (max_ring_len-1)/2 ) then

                                        n_new_walker              = n_new_walker+1
                                        if (n_new_walker > MAX_WALKER) then
                                           RAISE_ERROR("MAX_WALKER ("//MAX_WALKER//") exceeded.", error)
                                        endif

                                        new_walker(n_new_walker)  = j
                                        new_last(n_new_walker)    = i

                                        new_rings(1:ring_len(k), n_new_walker)               = rings(1:ring_len(k), k)
                                        new_ring_len(n_new_walker)                           = ring_len(k)+1
                                        new_rings(new_ring_len(n_new_walker), n_new_walker)  = j

                                        new_dr(:, 1:ring_len(k), n_new_walker)               = dr(:, 1:ring_len(k), n_new_walker)
                                        new_dr(:, new_ring_len(n_new_walker), n_new_walker)  = dr(:, ring_len(k), k) + d

                                     endif

                                  else

                                     if ( dist(j, a) == dist(i, a) .or. dist(j, a) == dist(i, a)-1) then

                                        ! Reverse search direction
                                        n_new_walker              = n_new_walker+1
                                        if (n_new_walker > MAX_WALKER) then
                                           RAISE_ERROR("MAX_WALKER ("//MAX_WALKER//") exceeded.", error)
                                        endif

                                        new_walker(n_new_walker)  = -j
                                        new_last(n_new_walker)    = i

                                        new_rings(1:ring_len(k), n_new_walker)               = rings(1:ring_len(k), k)
                                        new_ring_len(n_new_walker)                           = ring_len(k)+1
                                        new_rings(new_ring_len(n_new_walker), n_new_walker)  = j

                                        new_dr(:, 1:ring_len(k), n_new_walker)               = dr(:, 1:ring_len(k), k)
                                        new_dr(:, new_ring_len(n_new_walker), n_new_walker)  = dr(:, ring_len(k), k) + d

                                     else

                                        RAISE_ERROR("Something is wrong with the distance map.", error)

                                     endif

                                  endif

                               endif

                            endif ni_away_from_root_if
                         enddo ni_away_from_root_loop

                      else

                         i = -i

!                         ni_towards_root_loop: do ni = nl%seed(i), nl%last(i)
                         ni_towards_root_loop: do ni = 1, n_neighbours(at, i)
                            ni_towards_root_if: if ( .not. done(no(i)+ni) ) then

!                               DISTJ_SQ(p, nl, i, ni, j, d, abs_dr_sq)
                               j = neighbour(at, i, ni, diff=d)
                               abs_dr_sq = d .dot. d

                               if ( abs_dr_sq < cutoff_sq .and. ( .not. present(mask) .or. mask(j) ) .and. j /= last(k) ) then

                                  if ( j == a ) then

                                     d = d + dr(:, ring_len(k), k)

                                     if ( dot_product(d, d) < EPS ) then

                                        ! Now we need to check whether this ring is SP

                                        is_sp  = .true.

                                        ring_len(k)            = ring_len(k)+1
                                        rings(ring_len(k), k)  = a

                                        do m = 1, ring_len(k)
                                           do n = m+2, ring_len(k)

                                              dn  = n-m
                                              if (dn > ring_len(k)/2)  dn = ring_len(k)-dn

                                              if (dist(rings(n, k), rings(m, k)) /= dn) then 
                                                 is_sp  = .false.
                                              endif
                                           enddo
                                        enddo

                                        if (is_sp) then
                                           if (present(rings_out)) then
                                              n_rings_out = n_rings_out + 1
                                              if (n_rings_out > size(rings_out,2)) then
                                                 RAISE_ERROR('Too many rings to store in rings_out', error)
                                              end if
                                              ! save atoms around the ring
                                              rings_out(:, n_rings_out) = (/ring_len(k), rings(1:ring_len(k), k) /) 
                                              ! save displacement vectors around the ring
                                              dr_out(:, 1:ring_len(k), n_rings_out) = dr(:, 1:ring_len(k), k)
                                           end if
                                           stat(ring_len(k))  = stat(ring_len(k))+1
                                        endif

                                     endif

                                  else if (dist(j, a) == dist(i, a)-1) then

                                     if (all(rings(1:ring_len(k), k) /= j)) then

                                        n_new_walker              = n_new_walker+1
                                        if (n_new_walker > MAX_WALKER) then
                                           RAISE_ERROR("MAX_WALKER ("//MAX_WALKER//") exceeded.", error)
                                        endif

                                        new_walker(n_new_walker)  = -j
                                        new_last(n_new_walker)    = i

                                        new_rings(1:ring_len(k), n_new_walker)               = rings(1:ring_len(k), k)
                                        new_ring_len(n_new_walker)                           = ring_len(k)+1
                                        new_rings(new_ring_len(n_new_walker), n_new_walker)  = j

                                        new_dr(:, 1:ring_len(k), n_new_walker)               = dr(:, 1:ring_len(k), k)
                                        new_dr(:, new_ring_len(n_new_walker), n_new_walker)  = dr(:, ring_len(k), k) + d

                                     endif

                                  endif

                               endif

                            endif ni_towards_root_if
                         enddo ni_towards_root_loop

                      endif
                   enddo

                   n_walker              = n_new_walker
                   walker(1:n_walker)    = new_walker(1:n_walker)
                   last(1:n_walker)      = new_last(1:n_walker)

                   ring_len(1:n_walker)  = new_ring_len(1:n_walker)
                   rings(:, 1:n_walker)  = new_rings(:, 1:n_walker)

                   dr(:, :, 1:n_walker)  = new_dr(:, :, 1:n_walker)

                enddo while_walkers_exist

             endif

          enddo bonds_loop2
       endif
    enddo bonds_loop1

    deallocate(done)

  endsubroutine count_sp_rings

endmodule ringstat_module
