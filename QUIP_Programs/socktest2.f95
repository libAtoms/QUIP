program socktest2

  use libAtoms_module
  use Potential_module
  use SocketTools_module

  implicit none

  type(Atoms) :: at, cluster, qm_cluster
  type(DynamicalSystem) :: ds
  type(Potential) :: mmpot, sockpot, fmpot
  type(MPI_context) :: mpi

  character(STRING_LENGTH) :: ip, tmp, cluster_args
  integer :: port, client_id, n_atoms, i, n, label
  real(dp) :: dt = 1.0_dp
  real(dp), pointer :: force(:,:)

  cluster_args = 'cluster_calc_connect=F min_images_only hysteretic_connect_cluster_radius=10.0 hysteretic_buffer cluster_periodic_z hysteretic_buffer_inner_radius=6.0 single_cluster=T hysteretic_connect=F cluster_hopping=F randomise_buffer=F hysteretic_connect_outer_factor=1.5 hysteretic_buffer_outer_radius=8.0 cluster_hopping_nneighb_only=F hysteretic_connect_inner_factor=1.3'

  call system_initialise(PRINT_NORMAL, enable_timing=.true.)

  call initialise(mpi)

  if (mpi%my_proc == 0) then
     ! usage: socktest2 <ip> <port> <client_id> <label>
     call get_cmd_arg(1, ip)
     call get_cmd_arg(2, tmp)
     port = string_to_int(tmp)
     call get_cmd_arg(3, tmp)
     client_id = string_to_int(tmp)
     call get_cmd_arg(4, tmp)
     label = string_to_int(tmp)
     
     ! Read initial structure from .xyz
     call read(at, 'atoms.'//client_id//'.xyz')
     call set_value(at%params, 'label', label)

     call initialise(ds, at)
     call rescale_velo(ds, 300.0_dp)
     call finalise(at)

     call potential_filename_initialise(mmpot, 'IP SW', 'params.xml')
     call print(mmpot)

     call initialise(sockpot, 'SocketPot server_ip="'//trim(ip)//'" server_port='//port//' client_id='//client_id//' property_list=species:pos:hybrid_mark:hybrid_mark_1:hybrid_mark_2:hybrid_mark_3:qm_medoids:qm_cluster')
     call print(sockpot)

     call initialise(fmpot, 'ForceMixing '//trim(cluster_args)//' qm_args_str={'//trim(cluster_args)//'}', mmpot, sockpot)
     call print(fmpot)

     call set_cutoff(ds%atoms, cutoff(mmpot), cutoff_skin=1.0_dp)
     call calc_connect(ds%atoms)

     n = 0
     do while (.true.)
        call system_timer('step')
        call calc_connect(ds%atoms) ! only if atom has moved > skin will anything be done
      
        ! FIXME would be nice to do MM calc while we're waiting for QM, but breaks encapsualtion of SocketPot
        call system_timer('FM_calc')
        call calc(fmpot, ds%atoms, args_str='force')
        call system_timer('FM_calc')

        if (mod(n, 10) == 0) then
           call write(ds%atoms, 'traj.xyz', append=.true.)
        end if

        ! advance the dynamics using QM/MM forces
        call assign_property_pointer(ds%atoms, 'force', force)
        call advance_verlet(ds, dt, force)
        call ds_print_status(ds, 'D')
        
        n = n + 1
        call system_timer('step')
     end do

     call finalise(ds)
     call finalise(mmpot)
     call finalise(sockpot)
     call finalise(fmpot)

  end if
  call system_finalise

end program socktest2
