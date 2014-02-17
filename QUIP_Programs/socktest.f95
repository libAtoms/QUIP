program socktest

  use libAtoms_module
  use Potential_module
  use SocketTools_module

  implicit none

  type(Atoms) :: at
  type(Potential) :: pot
  type(MPI_context) :: mpi

  character(STRING_LENGTH) :: ip, tmp
  integer :: port, client_id, buffsize, n_atoms, i, n, label
  real(dp) :: e, v(3,3), lattice(3,3)
  real(dp), allocatable :: frac_pos(:,:), f(:,:)

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

     allocate(frac_pos(3, at%n))

     call potential_filename_initialise(pot, 'IP SW', 'params.xml')

     call print('Connecting to QUIP server on host '//trim(ip)//':'//port//' as client '//client_id)

     n = 0
     do while (.true.)

        allocate(f(3,at%n))

        call set_cutoff(at, cutoff(pot))
        call calc_connect(at)
        do i=1,50
           call calc(pot, at, energy=e, force=f, virial=v)
        end do

        call print('completed calculation n='//n//' label= '//label//' on '//at%n//' atoms. energy='//e)
        call socket_send_reftraj(ip, port, client_id, label, at%n, e, f, v)

        n = n + 1

        deallocate(f)

        call socket_recv_reftraj(ip, port, client_id, buffsize, label, n_atoms, lattice, frac_pos)

        if (n_atoms == 0 .or. n_atoms /= at%n) then
           call print('n_atoms='//n_atoms//', at%n='//at%n//' - shutting down QUIP server')
           exit
        end if

        do i=1, at%n
           at%pos(:,i) = at%lattice .mult. frac_pos(:, i)
        end do

     end do

     call finalise(at)
     call finalise(pot)

     deallocate(frac_pos)

  end if
  call system_finalise

end program socktest
