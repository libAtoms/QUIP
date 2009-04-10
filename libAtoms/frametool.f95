program frametool
  use System_module
  use linearalgebra_module
  use Atoms_module

  implicit none
  
  type(Atoms)::at
  type(inoutput)::stdin, stdout
  type(table)::aux
  type(table)::aux_out
  character(len=1024):: comment
  integer::status

  ! initialise program
  call system_initialise(SILENT)
  call initialise(stdin, "stdin")
  call initialise(stdout, "stdout", verbosity=NORMAL)

  call allocate(aux_out, 0, 12, 0, 0)
  status=0
  do while(status == 0)
     call read_ext_xyz(at, comment, aux, stdin, status)
     if(status /= 0) cycle
     call table_allocate(aux_out,max_length=aux%N)
     aux_out%N = aux%N
     aux_out%real(1:10,:) = aux%real(1:10,:)
     aux_out%real(11,:) = norm(aux%real(5:7,:)-aux%real(8:10,:), 1)
     aux_out%real(12,:) = norm(aux%real(2:4,:)-aux%real(5:7,:), 1)
     call print_ext_xyz(at, extras=aux_out, xyzfile=stdout)
  end do

  call Atoms_Finalise(at)
  call system_finalise()
end program frametool

