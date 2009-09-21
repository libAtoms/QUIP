program cxyzreader
use libatoms_module
implicit none
  type(Atoms) :: prim, at
  type(CInoutput) :: io
  integer i

  call system_initialise()

  call system_command("rm -f test.xyz test.xyz.idx")
  call print("writing test.xyz", ERROR)
  call initialise(io,"test.xyz", OUTPUT)
  call diamond(prim, 5.43_dp, (/ 14 /) )
  call supercell(at, prim, 3, 3, 3)
  do i=1, 200
    call write(at, io)
    at%pos(1,1) = at%pos(1,1) + 0.001
  end do
  call finalise(io)
  call finalise(at)

  call print("reading test.xyz", ERROR)
  call initialise(io,"test.xyz", INPUT)
  call read(io, at)
  call print("first read at%pos(1,1) " // at%pos(1,1))
  call read(io, at)
  call print("second read at%pos(1,1) " // at%pos(1,1))
  call finalise(io)
  call finalise(at)

  call initialise(io,"test.xyz", INPUT)
  call read(io, at, frame=50)
  call print("frame=50 read at%pos(1,1) " // at%pos(1,1))
  call read(io, at, frame=25)
  call print("frame=25 read at%pos(1,1) " // at%pos(1,1))
  call finalise(io)
  call finalise(at)

  call system_finalise()
end program cxyzreader
