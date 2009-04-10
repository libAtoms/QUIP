program test_Matrix

use System_module
use Matrix_module

implicit none

  type(MatrixD) md
  type(matrixZ) mz

  call system_initialise()

  call Print(logger, "Initialise MatrixD 5")
  call Initialise(md, 5)
  md%data = 1.0_dp
  call Print(logger, "Post initialise")
  call Print(md)

  call Print(logger, "Wipe MatrixD")
  call Wipe(md)
  call Print(logger, "Post wipe")
  call Print(md)

  call Print(logger, "Initialise MatrixD 3 5")
  call Initialise(md, 3, 5)
  md%data = 2.0_dp
  call Print(logger, "Post initialise")
  call Print(md)

  call Print(logger, "Initialise MatrixZ 5")
  call Initialise(mz, 5)
  mz%data = 1.0_dp
  call Print(logger, "Post initialise")
  call Print(mz)

  call Print(logger, "Wipe MatrixZ")
  call Wipe(mz)
  call Print(logger, "Post wipe")
  call Print(mz)

  call Print(logger, "Initialise MatrixZ 3 5")
  call Initialise(mz, 3, 5)
  mz%data = 2.0_dp
  call Print(logger, "Post initialise")
  call Print(mz)

  call system_finalise()

end program
