program socktest

  use libAtoms_module
  use Potential_module
  use SocketTools_module

  implicit none

  type(Atoms) :: at
  type(Potential) :: pot

  character(STRING_LENGTH) :: ip, tmp
  integer :: port, client_id, buffsize, z(100), n_atoms, i, n
  real(dp) :: e, v(3,3), frac_pos(3,100), lattice(3,3)
  real(dp), allocatable :: f(:,:)

  call system_initialise
  
  call get_cmd_arg(1, ip)
  call get_cmd_arg(2, tmp)
  port = string_to_int(tmp)
  call get_cmd_arg(3, tmp)
  client_id = string_to_int(tmp)
  buffsize = 10000

  mainlog%prefix = 'CLIENT '//client_id

  ! Read initial structure from .xyz
  call read(at, 'atoms.'//client_id//'.xyz')
  call potential_filename_initialise(pot, 'IP SW', 'params.xml')

  call print('Connecting to QUIP server on host '//trim(ip)//':'//port//' as client '//client_id)

  n = 0
  do while (.true.)

     allocate(f(3,at%n))

     call set_cutoff(at, cutoff(pot))
     call calc_connect(at)
     call calc(pot, at, energy=e, force=f, virial=v)

     call print('completed calculation '//n//' on '//at%n//' atoms')
     n = n + 1

     call print('Sending data')
     call socket_send_data(ip, port, client_id, at%n, e, f, v)
     call print('Finished sending data')
  
     deallocate(f)

     call print('Waiting to receive data...')
     call socket_recv_data(ip, port, client_id, buffsize, n_atoms, z, lattice, frac_pos)
     call print('Received data')

     if (n_atoms == 0 .or. n_atoms /= at%n) then
        call print('n_atoms='//n_atoms//', at%n='//at%n//' - shutting down QUIP server')
        exit
     end if

     do i=1, at%n
        at%pos(:,i) = at%lattice .mult. frac_pos(:, i)
     end do
     call set_atoms(at, z(1:at%n))

     call print('Setup atoms')

  end do

  call finalise(at)
  call finalise(pot)

  call system_finalise

end program socktest
