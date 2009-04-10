program test_KP

use System_module
use KPoints_module

implicit none

  type(KPoints) kp

  call system_initialise()

  call Print ("Initialise KP")
  call Initialise(kp)
  call Print ("Post initialise")
  call Print(kp)

  call Print ("")
  call Print ("Finalise KP")
  call Finalise(kp)
  call Print ("Post finalise")
  call Print(kp)
  call Print("")

  call Print ("Initialise KP from 4x4x1 mesh, no shift")
  call Initialise(kp, (/4,4,1/), .false.)
  call Print(kp)
  call Finalise(kp)
  call Print("")

  call Print ("Initialise KP from 4x4x1 mesh, shift")
  call Initialise(kp, (/4,4,1/), .true.)
  call Print(kp)
  call Finalise(kp)
  call Print("")

  call Print ("Initialise KP from 5x5x1 mesh, no shift")
  call Initialise(kp, (/5,5,1/), .false.)
  call Print(kp)
  call Finalise(kp)
  call Print("")

  call Print ("Initialise KP from 5x5x1 mesh, shift")
  call Initialise(kp, (/5,5,1/), .true.)
  call Print(kp)
  call Finalise(kp)
  call Print("")

  call Print ("Initialise KP from 4x5x1 mesh, no shift")
  call Initialise(kp, (/4,5,1/), .false.)
  call Print(kp)
  call Finalise(kp)
  call Print("")

  call Print ("Initialise KP from 4x5x1 mesh, shift")
  call Initialise(kp, (/4,5,1/), .true.)
  call Print(kp)
  call Finalise(kp)
  call Print("")

  call Print ("Initialise KP from 0x0x0 mesh, shift")
  call Initialise(kp, (/0,0,0/), .true.)
  call Print(kp)
  call Finalise(kp)
  call Print("")


  call system_finalise()
end program

