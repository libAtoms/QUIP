! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   This module written by Marco Caccin
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

module partition_module

  use iso_c_binding
  use error_module
  use system_module, only: dp, operator(//)
  use atoms_types_module
  use atoms_module
  use connection_module
  use clusters_module, only: HYBRID_ACTIVE_MARK
  use ringstat_module, only: distance_map

  implicit none
  private

  public :: kway_partition_atoms, partition_qm_list

  interface
     function METIS_SetDefaultOptions(options) bind(c)
       use iso_c_binding
       integer(C_INT),  dimension(0:40) :: options
     end function METIS_SetDefaultOptions

     function METIS_PartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, vsize, &
          adjwgt, nparts, tpwgts, ubvec, options, objval, part) bind(c)
       use iso_c_binding
       integer(C_INT) :: nvtxs, ncon, nparts, objval
       integer(C_INT), dimension(*) :: xadj, adjncy, part
       type(C_PTR) :: vwgt, vsize, adjwgt, tpwgts, ubvec
       integer(C_INT),  dimension(0:40) :: options
       integer(C_INT) :: METIS_PartGraphKway
     end function METIS_PartGraphKway
  end interface

  ! There is no METIS Fortran API, so we have to declare constants. The following are valid for versions >= 5. 
  ! Tested with 5.1.0_3, Nov 2014
  integer(C_INT), parameter :: METIS_OPTION_PTYPE=0, METIS_OPTION_OBJTYPE=1, METIS_OPTION_CTYPE=2, METIS_OPTION_IPTYPE=3, &
         METIS_OPTION_RTYPE=4, METIS_OPTION_DBGLVL=5, METIS_OPTION_NITER=6, METIS_OPTION_NCUTS=7, METIS_OPTION_SEED=8, &
         METIS_OPTION_NO2HOP=9, METIS_OPTION_MINCONN=10, METIS_OPTION_CONTIG=11, METIS_OPTION_COMPRESS=12, &
         METIS_OPTION_CCORDER=13, METIS_OPTION_PFACTOR=14, METIS_OPTION_NSEPS=15, METIS_OPTION_UFACTOR=16, &
         METIS_OPTION_NUMBERING=17, METIS_OPTION_HELP=18, METIS_OPTION_TPWGTS=19, METIS_OPTION_NCOMMON=20, &
         METIS_OPTION_NOOUTPUT=21, METIS_OPTION_BALANCE=22, METIS_OPTION_GTYPE=23, METIS_OPTION_UBVEC=24

contains

!   subroutine atoms_to_metis_graph(at, mask, error, filename)
!     implicit none
!     type(Atoms), intent(in) :: at
!     integer :: distance_matrix(at%N, at%N)
!     integer :: i, j
!     integer :: nnodes, nedges
!     logical, intent(in), optional :: mask(at%N)
!     integer, intent(out), optional  :: error
!     character(len=*), intent(in) :: filename
!     character(len=STRING_LENGTH), dimension(nnodes) :: lines
!     character(len=STRING_LENGTH) :: tempstr

!     call distance_map(at, distance_matrix, mask, error)
!     ! convert a bondhop distance matrix into lines of METIS graph file (line 1+1 := list of edges for node i)
!     do i = 1, nnodes
!        lines(i) = " "
!        do j = 1, nnodes
!           if (dist(i,j) == 1) then
!              if (j>i) then
!                 nedges = nedges + 1
!              end if
!              write( tempstr, '(i10)' ) j ! convert j to string
!              lines(i) = trim(lines(i))//" "//trim(tempstr) ! append it to list of neighbours of i on the same line
!           end if
!        end do
!     end do
!     ! write the METIS graph file corresponding to at
!     open(12, file=filename, status="replace", action="write")
!     write(12,*) nnodes, nedges
!     do i = 1, nnodes
!        write(12,*) trim(lines(i))
!     end do
!     close(12)
!   end subroutine atoms_to_metis_graph


  subroutine preprocess_atoms_to_csr_adjacency_matrix(at, nneightol, nedges)
    implicit none
    type(Atoms), intent(in) :: at
    integer :: i
    real(dp), intent(in) :: nneightol
    integer, intent(out) :: nedges
    type(Atoms) :: at_copy

    call atoms_copy_without_connect(at_copy, at) 
    at_copy%nneightol = nneightol
    !% Use hysteretic connect so we can work with relative cutoffs
    call calc_connect_hysteretic(at_copy, at_copy%nneightol, at_copy%nneightol)

    nedges = 0
    do i = 1, at%N
       nedges = nedges + n_neighbours(at_copy, i)
    end do

    call finalise(at_copy)

  end subroutine preprocess_atoms_to_csr_adjacency_matrix


  subroutine atoms_to_csr_adjacency_matrix(at, nedges, xadj, adjncy)
    implicit none
    type(Atoms), intent(in) :: at
    integer, intent(in) :: nedges
    integer :: i, j, edge
    integer :: nneighs, neigh_idx
    integer(C_INT), intent(out) :: xadj(at%N+1), adjncy(2*nedges)

    edge = 0
    do i = 1, at%N
       xadj(i) = edge
       nneighs = n_neighbours(at, i)
       do j = 1, nneighs
          edge = edge + 1
          adjncy(edge) = neighbour(at, i, j)
       end do
    end do
    xadj(at%N+1) = nedges

  end subroutine atoms_to_csr_adjacency_matrix


  subroutine kway_partition_atoms(at, k, nneightol, part, objval, status)

    ! Generate a balanced and connected partition of an atomic configuration using METIS' PartGraphKway method. 
    ! INPUT Atoms object at : configuration to be partitioned
    ! INPUT integer k : number of partitions
    ! OUTPUT integer array part, size at%N : partition array, assigning each atom of at to part i \in [1,k]

    implicit none
    type(Atoms), intent(in) :: at
    integer, intent(in) :: k
    real(dp), intent(in) :: nneightol
    integer, intent(out), dimension(:) :: part
    integer, intent(out) :: objval, status

    integer :: nedges
    integer(C_INT) :: nvert, ncon, nparts, c_status, c_objval
    type(C_PTR) :: vwgt, vsize, adjwgt, tpwgts, ubvec
    integer(C_INT), allocatable :: options(:), adjncy(:), xadj(:)

#ifdef HAVE_METIS
    integer(C_INT), allocatable :: cpart(:)

    allocate(options(0:40))
    c_status = METIS_SetDefaultOptions(options)
    options(:) = -1
    options(METIS_OPTION_NUMBERING) = 1 ! METIS_OPTION_NUMBERING: fortran indexing
    options(METIS_OPTION_NCUTS) = 5 ! number of different partitionings to compute (then choose the best edgecut) 
    options(METIS_OPTION_CONTIG) = 1 ! forces contiguous partitions

    nvert = at%n
    ncon = 1
    ! from at create csr adjacency matrix
    call preprocess_atoms_to_csr_adjacency_matrix(at, nneightol, nedges)
    allocate(adjncy(2*nedges), xadj(at%n + 1), cpart(at%n))
    call atoms_to_csr_adjacency_matrix(at, nedges, xadj, adjncy)

    ! perform partitioning using METIS
    vwgt = C_NULL_PTR
    vsize = C_NULL_PTR
    adjwgt = C_NULL_PTR
    tpwgts = C_NULL_PTR
    ubvec = C_NULL_PTR

    c_status = METIS_PartGraphKway(nvert, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, &
         tpwgts, ubvec, options, c_objval, part)

    ! Copy output
    objval = c_objval
    status = c_status
    if (size(cpart) > size(part)) then
       part(:) = cpart
    else
       call system_abort('size of part array ('//(size(part))//') insufficient to store results')
    end if

    deallocate(options, adjncy, cpart)
#else
    call system_abort('kway_partition_atoms requires HAVE_METIS=1 at compile time')
#endif

  end subroutine kway_partition_atoms

  function submatrix(matrix, mask)
    ! select submatrix choosing only rows and columns in mask. only works for square matrices
    implicit none
    integer, intent(in) :: matrix(:,:)
    logical, intent(in) :: mask(:)
    integer, allocatable :: submatrix(:,:), my_mask(:)
    integer :: i

    if (size(matrix(1,:)) /= size(matrix(:,1)) .or. size(matrix(1,:)) /= size(mask)) then
       write(*,*) "Inconsistent input, submatrix function now exiting..."
    end if
    do i = 1, size(mask)
       if (mask(i)) then
          my_mask(i) = i
       else
          my_mask(i) = 0
       end if
    end do
    submatrix = matrix(pack(my_mask, my_mask > 0), pack(my_mask, my_mask > 0))
  end function submatrix


  function costG(dm, part, idx, kappa)
    ! Calculate 'free energy' of region idx using distance matrix dm and assignment assign
    ! kappa measures the relative importance of surface and 1/diameter components of G (arbitrarily)
    integer, intent(in), dimension(:,:) :: dm
    integer, intent(in), dimension(:)   :: part
    integer, intent(in) :: idx
    real(dp), intent(in), optional :: kappa
    logical :: mask(size(part))
    integer, allocatable :: sub_dm(:,:)
    real(dp) :: costG, surface_energy, r, k
    integer :: i, j

    k = 1.0
    if(present(kappa)) k = kappa

    mask = (part == idx)
    sub_dm = submatrix(dm, mask)
    surface_energy = sum(sub_dm**2)

    if (size(shape(sub_dm)) == 1) then
       r = 1.0E-6_dp
    else
       r = maxval(sub_dm) / real(count(mask), dp)**2
    end if
    costG = surface_energy + k/r

  end function costG


  subroutine connect2subconnects(connect, nodelist, part, sub1, sub2, nodelist1, nodelist2, cut_edges)
    implicit none
    integer, intent(in) :: connect(:,:), nodelist(:)
    logical, intent(in) :: part(:)
    integer, intent(out), allocatable :: sub1(:,:), sub2(:,:), nodelist1(:), nodelist2(:), cut_edges(:,:)
    integer, allocatable :: cut_edges_temp(:,:)
    integer, dimension(size(nodelist)) :: part_mask1, part_mask2
    integer :: i, j, k, ncuts

    ! FIXME compiler doesn't like logical comparison
!     if ((size(connect(1,:) /= size(nodelist)) .or. size(nodelist) /= size(part) &
!          .or. size(connect(1,:)) /= size(part)) then
!        call system_abort('inconsistent sizes for connectivity matrix, nodelist and part')
!     end if

    do i = 1, size(part)
       if (part(i)) then
          part_mask1(i) = i
          part_mask2(i) = 0
       else
          part_mask1(i) = 0
          part_mask2(i) = i
       end if
    end do

    ! slice the connectivity matrix
    sub1 = connect(pack(part_mask1, part_mask1 > 0), pack(part_mask1, part_mask1 > 0))
    sub2 = connect(pack(part_mask2, part_mask2 > 0), pack(part_mask2, part_mask2 > 0))
    ! also return list of nodes for each new subgraph
    nodelist1 = nodelist(pack(part_mask1, part_mask1 > 0))
    nodelist2 = nodelist(pack(part_mask2, part_mask2 > 0))     

    ! number of edges is the sum of upper triangle "1"s
    ncuts = (sum(connect) - sum(sub1) - sum(sub2)) / 2
    write(*,*) "number of cut edges:", ncuts
    allocate(cut_edges_temp(ncuts, 2))     
    ! find edges that are missing from the two new subgraphs
    k = 1
    do i = 1, size(connect(1,:))
       do j = i+1, size(connect(1,:))
          if (connect(i,j) == 1) then
             if (any(nodelist1 .eq. nodelist(i)) .and. any(nodelist2 .eq. nodelist(j))) then
                cut_edges(k,:) = (/ nodelist(i), nodelist(j) /)
                k = k +1
             end if
          end if
       end do
    end do
    cut_edges = cut_edges_temp

    deallocate(cut_edges_temp)     
  end subroutine connect2subconnects

  function region_neighbouring_atoms(connect, nodelist, part, which)
    implicit none
    integer, intent(in) :: connect(:,:), part(:), nodelist(:), which
    integer, allocatable :: sub1(:,:), sub2(:,:), dummy1(:), dummy2(:), cut_edges(:,:)
    integer, dimension(2) :: edge, regions
    integer, allocatable :: neigh_temp(:,:), region_neighbouring_atoms(:,:)
    integer i
    ! generate a partition of the whole graph using part == which as the carving mask.
    ! returns a list of cut edges
    call connect2subconnects(connect, nodelist, part==which, &
         sub1, sub2, dummy1, dummy2, cut_edges)

    allocate(neigh_temp(size(cut_edges(:,1)),3))
    do i = 1, size(cut_edges(:,1))
       edge = cut_edges(i,:)
       regions = (/part(minloc(abs(nodelist - edge(1)), 1)), &
            part(minloc(abs(nodelist - edge(2)), 1))/)
       ! which region other node belongs to
       if (regions(1) == which) then
          neigh_temp(i,1) = regions(2)
       else
          neigh_temp(i,1) = regions(1)
       end if
       ! which edge
       neigh_temp(i,2:3) = edge
    end do
    region_neighbouring_atoms = neigh_temp

    deallocate(neigh_temp)
  end function region_neighbouring_atoms


  subroutine partition_qm_list(at, k, nneightol, ripening, error)

    implicit none
    type(Atoms), intent(inout) :: at
    integer, intent(in) :: k
    real(dp), intent(in) :: nneightol
    logical, intent(in) :: ripening
    integer, optional, intent(out) :: error

    type(Atoms) :: qm_at
    integer, pointer :: hybrid_mark(:), orig_index(:), qm_cluster_ptr(:)
    logical, allocatable :: my_hybrid_mask(:)
    integer :: i, objval, status
    integer, dimension(:), allocatable :: part

    INIT_ERROR(error)

    call assign_property_pointer(at, 'hybrid_mark', hybrid_mark, error)
    PASS_ERROR(error)

    allocate(my_hybrid_mask(at%n))
    my_hybrid_mask = (hybrid_mark == HYBRID_ACTIVE_MARK)

    ! select qm region here. NB: qm_list is now qm_at%orig_index
    call select(qm_at, at, mask=my_hybrid_mask, orig_index=.true.)
    allocate(part(qm_at%N))

    ! then generate the part vector (was called assign, unfortunate choice for a fortran code)
    call kway_partition_atoms(qm_at, k, nneightol, part, objval, status)

    ! now try writing down the ripening stuff in fortran

    ! add property qm_cluster
    call add_property(at, 'qm_cluster', 0, ptr=qm_cluster_ptr, overwrite=.true., error=error)
    PASS_ERROR(error)

    ! zero padding qm_cluster property
    do i = 1, at%N
       qm_cluster_ptr(i) = 0
    end do
    ! update values in qm region

    call assign_property_pointer(at, 'orig_index', orig_index, error)
    PASS_ERROR(error)

    do i = 1, qm_at%N
       qm_cluster_ptr(orig_index(i)) = part(i)
    end do

    deallocate(part)
    deallocate(my_hybrid_mask)
    ! WARNING: this routine does not return medoids, which are anyway useless at the current stage of partitioning method

  end subroutine partition_qm_list

end module partition_module
