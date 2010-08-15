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

program domdecomp_example
  use libAtoms_module
  use QUIP_module

  implicit none

  real(DP), parameter :: a0 = 3.0_DP
  real(DP), parameter :: dt = 1.0_DP

  character(len=10240) :: lj_str
  type(Potential) :: pot
  type(Atoms) :: at, cell
  type(MPI_context) :: mpi
  type(DomainDecomposition) :: dd
  type(DynamicalSystem) :: ds
  type(InOutput) :: ener

  real(dp) :: e, ekin, dr, dr2, v(3)
  real(dp), allocatable :: f(:,:)
  real(DP), pointer :: local_e(:)
  integer :: it, i, j

  integer :: error = ERROR_NONE

  call system_initialise()

  lj_str = '<LJ_params n_types="1">' // &
       '<per_type_data type="1" atomic_num="29" />' // &
       '<per_pair_data type1="1" type2="1" sigma="2.0" eps6="2.0" eps12="2.0" cutoff="2.5" energy_shift="T" linear_force_shift="T" />' // &
       '</LJ_params>'

  call print("dt = " // dt)

  call Initialise(pot, "IP LJ", param_str=lj_str)

  call print(pot)

  call fcc(cell, a0, 29)
!  call sc(cell, a0, 29)
!  cell%pos(1, :) = cell%pos(1, :) + a0/4
!  cell%pos(2, :) = cell%pos(2, :) + a0/4
!  cell%pos(3, :) = cell%pos(3, :) + a0/4
  cell%pos(1, :) = cell%pos(1, :) + 0.001_DP
  cell%pos(2, :) = cell%pos(2, :) + 0.001_DP
  cell%pos(3, :) = cell%pos(3, :) + 0.001_DP
!  cell%pos(1, :) = cell%pos(1, :) + a0/2
!  cell%pos(2, :) = cell%pos(2, :) + a0/2
!  cell%pos(3, :) = cell%pos(3, :) + a0/2
  call supercell(at, cell, 10, 10, 10)
!  call supercell(at, cell, 4, 1, 1)

  call print("cell = " // at%lattice(1, :))
  call print("       " // at%lattice(2, :))
  call print("       " // at%lattice(3, :))

  do i = 1, at%N
     at%pos(:, i) = at%pos(:, i) + 0.1_DP*sin(2*PI*(at%g .mult. at%pos(:, i)))
  enddo
!  at%pos(1, 1) = at%pos(1, 1) - a0/4
!  at%pos(1, 2) = at%pos(1, 2) + a0/4
!  at%pos(1, 3) = at%pos(1, 3) - a0/4
!  at%pos(1, 4) = at%pos(1, 4) + a0/4

  call print("pos(1) = " // at%pos(:, 1))

  call initialise(ds, at)

  ds%atoms%velo(1, 1) = -5e-4_DP

  call add_property(ds%atoms, "local_e", 0.0_DP, ptr=local_e)

!  call verbosity_push(PRINT_VERBOSE)

  call initialise(mpi)
#ifdef _MPI
  call initialise(dd, ds%atoms, mpi, (/ 2, 2, 2 /), error=error)
#else
  call initialise(dd, ds%atoms, mpi, (/ 1, 1, 1 /), error=error)
#endif
  HANDLE_ERROR(error)
  call allocate(dd, ds%atoms, error=error)
  HANDLE_ERROR(error)

  call set_border(dd, ds%atoms, cutoff(pot))

  call enable(dd, ds%atoms)

  call print("After enabling domain decomposition.")
  call print("N = " // ds%atoms%N // ", Ndomain  = " // ds%atoms%Ndomain // ", Ntotal   = " // dd%Ntotal)

  call set_cutoff(ds%atoms, cutoff(pot))


  allocate(f(3,ds%atoms%Nbuffer))
  !    call rescale_velo(ds, 300.0_dp)

  call calc_connect(ds%atoms, error=error)
  HANDLE_ERROR(error)

  call calc(pot, ds%atoms, local_e = local_e, f = f, error=error)
  HANDLE_ERROR(error)

  e = sum(local_e(1:ds%atoms%Ndomain))

  call initialise(ener, "ener.out", action=OUTPUT)
  call print("0 0 " // e // " " // e, file=ener)

  call ds_print_status(ds, 'D', e, &
       instantaneous  = .true., &
       mpi_obj        = dd%mpi, &
       error          = error)
  HANDLE_ERROR(error)

  do it = 1, 1000
     call advance_verlet1(ds, dt, f, error=error)
     HANDLE_ERROR(error)

     call communicate_domain(dd, ds%atoms, error=error)
     HANDLE_ERROR(error)
     call communicate_ghosts(dd, ds%atoms, .true., error=error)
     HANDLE_ERROR(error)

     call calc_connect(ds%atoms, error=error)
     HANDLE_ERROR(error)

!     do i = 1, atoms_n_neighbours(ds%atoms, 1)
!        j = atoms_neighbour(ds%atoms, 1, i, distance=dr)
!        v = ds%atoms%pos(:, i) - ds%atoms%pos(:, j)
!        v = ds%atoms%lattice .mult. ((ds%atoms%g .mult. v) - nint(ds%atoms%g .mult. v))
!        dr2 = norm(v)
!        call print("neighbour = " // j // ", dr = " // dr // ", dr2 = " // dr2)
!        call print("posi = " // ds%atoms%pos(:, i))
!        call print("posj = " // ds%atoms%pos(:, j))
!     enddo

     call calc(pot, ds%atoms, local_e = local_e, f = f, error=error)
     HANDLE_ERROR(error)

     call print("pos(1) = " // ds%atoms%pos(:, 1))
     call print("vel(1) = " // ds%atoms%velo(:, 1))
     call print("for(1) = " // f(:, 1))

!     if (it == 2) then
!     write (2000+dd%mpi%my_proc, '(A)')  "---"
!     do i = 1, ds%atoms%N
!        write (2000+dd%mpi%my_proc, '(I10,9ES20.10)') &
!             dd%local_to_global(i), &
!             ds%atoms%pos(:, i), ds%atoms%velo(:, i), f(:, i)
!     enddo
!     endif

     e = sum(local_e(1:ds%atoms%Ndomain))

     call advance_verlet2(ds, dt, f, error=error)
     HANDLE_ERROR(error)

!     write (1000+dd%mpi%my_proc, '(A,3F20.10)') "r(:, 1)    = ", ds%atoms%pos(:, 1)
!     write (1000+dd%mpi%my_proc, '(A,3F20.10)') "v(:, 1)    = ", ds%atoms%velo(:, 1)
!     write (1000+dd%mpi%my_proc, '(A,3F20.10)') "f(:, 1)    = ", f(:, 1)
!     write (1000+dd%mpi%my_proc, '(A,3F20.10)') "r(:, 2)    = ", ds%atoms%pos(:, 2)
!     write (1000+dd%mpi%my_proc, '(A,3F20.10)') "v(:, 2)    = ", ds%atoms%velo(:, 2)
!     write (1000+dd%mpi%my_proc, '(A,3F20.10)') "f(:, 2)    = ", f(:, 2)
!     write (1000+dd%mpi%my_proc, '(A,3F20.10)') "r(:, 2001) = ", ds%atoms%pos(:, 2001)
!     write (1000+dd%mpi%my_proc, '(A,3F20.10)') "v(:, 2001) = ", ds%atoms%velo(:, 2001)
!     write (1000+dd%mpi%my_proc, '(A,3F20.10)') "f(:, 2001) = ", f(:, 2001)
!     write (1000+dd%mpi%my_proc, '(A,3F20.10)') "r(:, 2002) = ", ds%atoms%pos(:, 2002)
!     write (1000+dd%mpi%my_proc, '(A,3F20.10)') "v(:, 2002) = ", ds%atoms%velo(:, 2002)
!     write (1000+dd%mpi%my_proc, '(A,3F20.10)') "f(:, 2002) = ", f(:, 2002)


!     call write(ds%atoms, "proc_" // dd%mpi%my_proc // ".xyz", append=.true.)

     call ds_print_status(ds, 'D', e, &
          instantaneous  = .true., &
          mpi_obj        = dd%mpi, &
          error          = error)
     HANDLE_ERROR(error)

     call sum_in_place(dd%mpi, e, error=error)
     HANDLE_ERROR(error)
     ekin = kinetic_energy(ds, dd%mpi, error=error)
     HANDLE_ERROR(error)
     call print(it // " " // ekin // " " // e // " " // (ekin+e), file=ener)

     call print("N = " // ds%atoms%N // ", Ndomain  = " // ds%atoms%Ndomain // ", Ntotal   = " // dd%Ntotal)
  end do
  call finalise(ener)

  call finalise(dd)
  call finalise(mpi)

endprogram domdecomp_example
