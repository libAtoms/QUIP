program socktest2

  use libAtoms_module
  use Potential_module
  use SocketTools_module

  implicit none

  type(Atoms) :: at, cluster, qm_cluster
  type(Potential) :: pot
  type(MPI_context) :: mpi

  character(STRING_LENGTH) :: ip, tmp, cluster_args
  integer :: port, client_id, buffsize, n_atoms, i, n, label, old_n
  real(dp) :: e

  cluster_args = 'cluster_calc_connect=F min_images_only hysteretic_connect_cluster_radius=10.0 hysteretic_buffer cluster_periodic_z hysteretic_buffer_inner_radius=6.0 single_cluster=F hysteretic_connect=F cluster_hopping=F randomise_buffer=F hysteretic_connect_outer_factor=1.5 hysteretic_buffer_outer_radius=8.0 cluster_hopping_nneighb_only=F hysteretic_connect_inner_factor=1.3'

  call system_initialise(enable_timing=.true.)

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
        call system_timer('step')

        call system_timer('make_cluster')
        call set_cutoff(at, cutoff(pot))
        call calc_connect(at)
        call create_hybrid_weights(at, args_str=cluster_args)
        call create_cluster_simple(at, args_str=cluster_args, cluster=cluster)
        call system_timer('make_cluster')

        call system_timer('socket_send')
        call set_value(cluster%params, 'label', label)
        call socket_send_xyz(ip, port, client_id, cluster)
        call system_timer('socket_send')

        ! do MM calc while we're waiting for QM
        call system_timer('MM calc')
        call calc(pot, at, args_str="energy force virial")
        if (.not. get_value(at%params, 'energy', e)) &
             call system_abort('energy value missing from atoms after calc() call')
        call print('completed MM calculation n='//n//' label= '//label//' on '//at%n//' atoms. energy='//e)
        call system_timer('MM calc')

        call system_timer('socket_recv')
        n = n + 1
        old_n = at%n
        call socket_recv_xyz(ip, port, client_id, buffsize, qm_cluster)
        call system_timer('socket_recv')

        if (.not. get_value(qm_cluster%params, 'label', label)) call system_abort('missing label param in received atoms')
        call print('Received qm_cluster label='//label)

        ! fill in QM forces from qm_cluster
        ! advance the dynamics
        
        call system_timer('step')
     end do

     call finalise(at)
     call finalise(cluster)
     call finalise(qm_cluster)
     call finalise(pot)

  end if
  call system_finalise

end program socktest2
