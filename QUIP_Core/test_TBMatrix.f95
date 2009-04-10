program test_TBMatrix

use System_module
use Atoms_module
use TBMatrix_module

implicit none
  
  type (TBMatrix) tbm
  type (TBVector) tbv

  call system_initialise()

  call Print (logger, 'Initialise TBMatrix N=8')
  call Initialise(tbm, 8)
  call Print (logger, "Post initialise")
  call Print(tbm)

  call Print (logger, "Wipe")
  call Wipe(tbm)
  call Print (logger, "Post wipe")
  call Print(tbm)

  call Print (logger, 'Initialise TBMatrix N=8 n_matrices=2, is_complex')
  call Initialise(tbm, 8, 2, .true.)
  call Print (logger, "Post initialise")
  call Print(tbm)

  call Print (logger, "Wipe")
  call Wipe(tbm)
  call Print (logger, "Post wipe")
  call Print(tbm)

  call Print (logger, "Initialise TBVector")
  call Initialise(tbv, 8)
  call Print (logger, "Post initialise")
  call Print(tbv)

  call Print (logger, "Wipe")
  call Wipe(tbv)
  call Print (logger, "Post wipe")
  call Print(tbv)

  call system_finalise()

end program
