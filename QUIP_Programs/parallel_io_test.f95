program parallel_io_test

  use libAtoms_module
  implicit none

  type(InOutput) :: init
  type(CInOutput) :: in, out
  type(Atoms)     :: at
  type(Dictionary) :: dict
  type(MPI_Context) :: mpi
  type(Table) :: tab
  logical :: texist
  integer :: iunit

  call system_initialise(mpi_all_inoutput=.true.)

  call initialise(mpi)
  
  call print_title('Test dict_bcast')
  call initialise(dict)
  if (mpi%my_proc == 0) then
     call set_value(dict, 'a', 0)
     call set_value(dict, 'b', 1.0_dp)
     call set_value(dict, 'c', (/.false., .false., .false./))
     call set_value(dict, 'd', (/3, 1, 4, 1, 5, 9/))
     call set_value(dict, 'e', (/3.141459_dp/))
     call set_value(dict, 'f', (1.0_dp, 2.0_dp))
  end if
  call bcast(mpi, dict)
  call print(dict)

  call print_title('Test table_bcast')
  call initialise(tab, 1, 0, 0, 0)
  if (mpi%my_proc == 0) then
     call append(tab, (/1, 2, 3, 4/))
  end if
  call bcast(mpi, tab)
  call print(tab)

  if (mpi%my_proc == 0) then
     ! Write input file
     call initialise(init, 'input.xyz', OUTPUT)
     call print('8', file=init)
     call print('Lattice="5.440000 0.000000 0.000000 0.000000 5.440000 0.000000 0.000000 0.000000 5.440000" Properties=species:S:1:pos:R:3:Z:I:1:log:L:1', file=init)
     call print('Si              0.00000000      0.00000000      0.00000000      14    F', file=init)
     call print('Si              1.36000000      1.36000000      1.36000000      14    F', file=init)
     call print('Si              2.72000000      2.72000000      0.00000000      14    F', file=init)
     call print('Si              4.08000000      4.08000000      1.36000000      14    F', file=init)
     call print('Si              2.72000000      0.00000000      2.72000000      14    F', file=init)
     call print('Si              4.08000000      1.36000000      4.08000000      14    F', file=init)
     call print('Si              0.00000000      2.72000000      2.72000000      14    F', file=init)
     call print('Si              1.36000000      4.08000000      4.08000000      14    F', file=init)
     call finalise(init)

     ! Remove index if present
     inquire (file='input.xyz.idx',exist=texist)
     if (texist) then
        iunit = pick_up_unit()
        open (iunit,file='input.xyz.idx',status='old')
        close (iunit,status='delete')
     endif
  end if
  
  call print_title('Test atoms_bcast')
  call initialise(in, 'input.xyz', INPUT, mpi=mpi)
  call initialise(out, 'output.xyz', OUTPUT, mpi=mpi)

  call read(in, at)
  call print(at)
  call write(out, at)

  call system_finalise

end program parallel_io_test
