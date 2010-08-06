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

    character(len=10240) :: lj_str
    type(Potential) :: pot
    type(Atoms) :: at, cell
    type(MPI_context) :: mpi
    type(DomainDecomposition) :: dd
    type(DynamicalSystem) :: ds

    real(dp) :: e
    real(dp), allocatable :: f(:,:)
    real(DP), pointer :: local_e(:)
    integer :: it

    integer :: error = ERROR_NONE

    call system_initialise()

    lj_str = '<LJ_params n_types="1">' // &
         '<per_type_data type="1" atomic_num="29" />' // &
         '<per_pair_data type1="1" type2="1" sigma="2.0" eps6="2.0" eps12="2.0" cutoff="4.0" energy_shift="T" linear_force_shift="T" />' // &
         '</LJ_params>'

    call Initialise(pot, "IP LJ", param_str=lj_str)

    call print(pot)

    call fcc(cell, a0, 29)
    cell%pos(1, :) = cell%pos(1, :) + a0/4
    cell%pos(2, :) = cell%pos(2, :) + a0/4
    cell%pos(3, :) = cell%pos(3, :) + a0/4
    call supercell(at, cell, 10, 10, 10)

    call initialise(ds, at)

    call add_property(ds%atoms, "local_e", 0.0_DP, ptr=local_e)

    !    call verbosity_push(PRINT_VERBOSE)

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
    call print("N        = " // ds%atoms%N)
    call print("Ndomain  = " // ds%atoms%Ndomain)
    call print("Ntotal   = " // dd%Ntotal)

    call set_cutoff(ds%atoms, cutoff(pot))


    allocate(f(3,ds%atoms%Nbuffer))
    !    call rescale_velo(ds, 300.0_dp)

    call calc_connect(ds%atoms, error=error)
    HANDLE_ERROR(error)

    call calc(pot, ds%atoms, local_e = local_e, f = f, error=error)
    HANDLE_ERROR(error)

    e = sum(local_e(1:ds%atoms%Ndomain))

    call ds_print_status(ds, 'D', e, &
         instantaneous  = .true., &
         mpi_obj        = dd%mpi, &
         error          = error)
    HANDLE_ERROR(error)

    do it=1, 3
       call advance_verlet1(ds, 1.0_dp, error=error)
       HANDLE_ERROR(error)

       call calc_connect(ds%atoms, error=error)
       HANDLE_ERROR(error)

       call calc(pot, ds%atoms, local_e = local_e, f = f, error=error)
       HANDLE_ERROR(error)

       e = sum(local_e(1:ds%atoms%Ndomain))

       call advance_verlet2(ds, 1.0_dp, f, error=error)
       HANDLE_ERROR(error)

       call ds_print_status(ds, 'D', e, &
            instantaneous  = .true., &
            mpi_obj        = dd%mpi, &
            error          = error)
       HANDLE_ERROR(error)

       call communicate_domain(dd, ds%atoms, error=error)
       HANDLE_ERROR(error)
       call communicate_ghosts(dd, ds%atoms, .true., error=error)
       HANDLE_ERROR(error)
    end do

    call finalise(dd)
    call finalise(mpi)

  end program domdecomp_example
