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

  call system_initialise(verbosity=PRINT_ANAL)
!  call system_initialise

  lj_str = '<LJ_params n_types="1">' // &
       '<per_type_data type="1" atomic_num="29" />' // &
       '<per_pair_data type1="1" type2="1" sigma="2.0" eps6="2.0" eps12="2.0" cutoff="2.5" energy_shift="T" linear_force_shift="T" />' // &
       '</LJ_params>'

  call print("dt = " // dt)

  call Initialise(pot, "IP LJ", param_str=lj_str)

  call print(pot)

#if 0
  call fcc(cell, a0, 29)
  cell%pos(1, :) = cell%pos(1, :) + 0.001_DP
  cell%pos(2, :) = cell%pos(2, :) + 0.001_DP
  cell%pos(3, :) = cell%pos(3, :) + 0.001_DP
  call supercell(at, cell, 10, 10, 10)

  call print("cell = " // at%lattice(1, :))
  call print("       " // at%lattice(2, :))
  call print("       " // at%lattice(3, :))
  do i = 1, at%N
     at%pos(:, i) = at%pos(:, i) + 0.1_DP*sin(2*PI*(at%g .mult. at%pos(:, i)))
  enddo
#endif

!  call verbosity_push(PRINT_VERBOSE)

  call initialise(mpi)
#ifdef _MPI
  call initialise(dd, mpi, (/ 2, 2, 2 /), error=error)
#else
  call initialise(dd, mpi, (/ 1, 1, 1 /), error=error)
#endif
  HANDLE_ERROR(error)

#if 0
  if (mpi%my_proc == 0) then
     call write(ds%atoms, "01.xyz", error=error)
     HANDLE_ERROR(error)
  endif
#endif

!#define SERIAL_READ

#ifdef SERIAL_READ
  ! Read first, enable later
  call read(at, "initial.xyz", error=error)
  HANDLE_ERROR(error)
  call initialise(ds, at, error=error)
  HANDLE_ERROR(error)
  call set_border(dd, ds%atoms, cutoff(pot))
  call add_property(ds%atoms, "local_e", 0.0_DP, ptr=local_e, error=error)
  HANDLE_ERROR(error)
  call enable(dd, ds%atoms, error=error)
  HANDLE_ERROR(error)
#else
  ! Enable first, parallel read later
  call print("...enable...")
  call enable(dd, error=error)
  HANDLE_ERROR(error)
  call print("...read...")
  call read(at, "initial.xyz", domain=dd, error=error)
  HANDLE_ERROR(error)

  call print("...initialise ds...")
  call initialise(ds, at, error=error)
  HANDLE_ERROR(error)
  call print("...set_border...")
  call set_border(dd, ds%atoms, cutoff(pot))
  call print("...add_property...")
  call add_property(ds%atoms, "local_e", 0.0_DP, ptr=local_e, error=error)
  HANDLE_ERROR(error)
#endif

  call print("...disable...")
  call disable(dd, ds%atoms, error=error)
  HANDLE_ERROR(error)

  call print("...write...")
  if (mpi%my_proc == 0) then
     call write(ds%atoms, "02.xyz", error=error)
     HANDLE_ERROR(error)
  endif

  call print("...enable...")
  call enable(dd, ds%atoms)

  call print("After enabling domain decomposition.")
  call print("N = " // ds%atoms%N // ", Ndomain  = " // ds%atoms%Ndomain // ", Ntotal   = " // dd%Ntotal)

  call set_cutoff(ds%atoms, cutoff(pot))


  allocate(f(3,ds%atoms%Nbuffer))
  !    call rescale_velo(ds, 300.0_dp)

  call calc_connect(ds%atoms, error=error)
  HANDLE_ERROR(error)

  call calc(pot, ds%atoms, local_energy = local_e, force = f, error=error)
  HANDLE_ERROR(error)

  e = sum(local_e(1:ds%atoms%Ndomain))

  call initialise(ener, "ener.out", action=OUTPUT)
  call print("0 0 " // e // " " // e, file=ener)

  call ds_print_status(ds, 'D', e, &
       instantaneous  = .true., &
       mpi_obj        = dd%mpi, &
       error          = error)
  HANDLE_ERROR(error)

  do it = 1, 10
     call advance_verlet1(ds, dt, f, error=error)
     HANDLE_ERROR(error)

     call communicate_domain(dd, ds%atoms, error=error)
     HANDLE_ERROR(error)
     call communicate_ghosts(dd, ds%atoms, .true., error=error)
     HANDLE_ERROR(error)

     call calc_connect(ds%atoms, error=error)
     HANDLE_ERROR(error)

     call calc(pot, ds%atoms, local_energy = local_e, force = f, error=error)
     HANDLE_ERROR(error)

     e = sum(local_e(1:ds%atoms%Ndomain))

     call advance_verlet2(ds, dt, f, error=error)
     HANDLE_ERROR(error)

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
  end do
  call finalise(ener)

  call print("...disable...")
  call disable(dd, ds%atoms, error=error)
  HANDLE_ERROR(error)

  call print("...")

  call print("Ndomain = " // at%Ndomain)

  call print("...finalise(dd)...")
  call finalise(dd)
  call print("...finalise(mpi)...")
  call finalise(mpi)
  call print("...done...")

endprogram domdecomp_example
