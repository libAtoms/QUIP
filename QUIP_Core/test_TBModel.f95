program test_TBModel

use System_module
use TBModel_module

implicit none

  type(TBModel) tbm

  call system_initialise()

  call Print ("Initialize TBModel NRL-TB")
  call Initialise(tbm, "NRL-TB")
  call Print ("Post initialize")
  call Print(tbm)

  call Print ("Finalise")
  call Finalise(tbm)

  call system_finalise()

end program
