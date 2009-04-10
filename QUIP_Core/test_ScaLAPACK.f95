program test_Matrix

use System_module
use Matrix_module

implicit none

  type(MatrixD) md, evecs, t
  real(dp), allocatable :: evals(:)
  type(MPI_context) mpi_glob
  type(ScaLAPACK) scalapack_glob

  integer i, j

  call system_initialise()

  call Initialise(mpi_glob)
  call Initialise(scalapack_glob, mpi_glob)

  call Print("Initialise MatrixD 8x8 with scalapack")
  call Initialise(md, 8, 8, 2, 2, scalapack_obj = scalapack_glob)

  call Print("initial md desc " // md%ScaLAPACK_Info_obj%desc)
  do i=1, md%l_N
  do j=1, md%l_M
    md%data(i,j) = (real(i+j,dp)+sqrt(real(i+j,dp)))/20.0_dp
  end do
  end do
  call Print("Post initialise")
  call Print(md)

  call system_finalise()
  stop

  call Print("Call add_identity")
  call add_identity(md)
  call Print("Post add_identity")
  call Print(md)

  call Initialise(evecs, 8, 8, 2, 2, scalapack_obj = scalapack_glob)
  call Print("initial evecs desc " // evecs%ScaLAPACK_Info_obj%desc)
  allocate(evals(8))

  call Print("Call diagonalise")
  call diagonalise(md, evals, evecs)

  call Print("evals " // evals)
  call Print(evecs,name="evecs")

  call Initialise(t, 8, 8, 2, 2, scalapack_obj = scalapack_glob)

  call Print(md, name="md")

  call mult(t, md, evecs)

  call Print(t, name="md evecs")

  t%data = t%data / evecs%data

  call Print(t, name="md evecs / evecs")


end program
