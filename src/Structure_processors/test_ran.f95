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

program test_ran
  use libatoms_module

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  integer, parameter :: nrand = 100000

  type(Histogram1D), save :: h
  !$omp threadprivate(h)
  integer :: i, n
  real(DP) :: d

  ! ---

  call system_initialise(PRINT_NORMAL)

  call print("Start generating random numbers on proc 0...")

  if (mpi_id() == 0) then
     !$omp parallel default(none) private(i, n, d)
     do i = 1, nrand
        d = ran_normal()
     enddo
     !$omp end parallel
  endif

  call print("...done")

  call print("Synching rngs.")
  call system_resync_rng

  call print("Start generating random numbers on all procs...")

  !$omp parallel default(none) private(i, n, d)
  call initialise(h, 1000, -10.0_DP, 10.0_DP)

  do i = 1, nrand
     call add(h, ran_normal())
  enddo

#ifdef _OPENMP
  n = omp_get_thread_num()
#else
  n = 0
#endif

  call write(h, "normal_"//mpi_id()//"_"//n//".out")
  call finalise(h)
  !$omp end parallel

  call print("...done")

  call print("Please check that all normal_*_1.out are identical.")

  call system_finalise()

end program test_ran
