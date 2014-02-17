program socktest2

  use libAtoms_module
  use Potential_module
  use SocketTools_module

  implicit none

  type(Atoms) :: at
  type(Potential) :: pot
  type(MPI_context) :: mpi

  character(STRING_LENGTH) :: ip, tmp
  integer :: port, client_id, buffsize, n_atoms, i, n, label, old_n
  real(dp) :: e

  call system_initialise

  call initialise(mpi)

  if (mpi%my_proc == 0) then
  
     call get_cmd_arg(1, ip)
     call get_cmd_arg(2, tmp)
     port = string_to_int(tmp)
     call get_cmd_arg(3, tmp)
     client_id = string_to_int(tmp)
     buffsize = 1000000
     call get_cmd_arg(4, tmp)
     label = string_to_int(tmp)
     
     mainlog%prefix = 'CLIENT '//client_id

     ! Read initial structure from .xyz
     call read(at, 'atoms.'//client_id//'.xyz')
     old_n = at%n

     call potential_filename_initialise(pot, 'IP SW', 'params.xml')

     call print('Connecting to QUIP server on host '//trim(ip)//':'//port//' as client '//client_id)

     n = 0
     do while (.true.)

        call set_cutoff(at, cutoff(pot))
        call calc_connect(at)
        call calc(pot, at, args_str="energy force virial")
        if (.not. get_value(at%params, 'energy', e)) &
             call system_abort('energy value missing from atoms after calc() call')

        call set_value(at%params, 'label', label)
        call print('completed calculation n='//n//' label= '//label//' on '//at%n//' atoms. energy='//e)
        call socket_send_xyz(ip, port, client_id, at)

        n = n + 1

        old_n = at%n
        call socket_recv_xyz(ip, port, client_id, buffsize, at)
        if (.not. get_value(at%params, 'label', label)) call system_abort('missing label param in received atoms')

        if (at%n == 0 .or. old_n /= at%n) then
           call print('old_n='//old_n//', at%n='//at%n//' - shutting down QUIP server')
           exit
        end if


     end do

     call finalise(at)
     call finalise(pot)

  end if
  call system_finalise

end program socktest2
