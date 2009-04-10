program test
use libatoms_module
use structures_module


implicit none
  type(Atoms) at, cluster
  type(Table) embedlist
  real(dp) :: axes(3,3)
  integer :: Z(2), i
  real(dp) :: width, height, a
  Type(inoutput) ::xyzfile

  call system_initialise(verbosity=ANAL)
  call initialise(xyzfile, "cluster.xyz")


  ! axes are  (111)[11b0], so
  axes(:,2) = (/ 1, 1, 1 /)
  axes(:,3) = (/1, -1, 0 /)
  axes(:,1) = axes(:,2) .cross. axes(:,3)


  Z(1) = 14
  Z(2) = 6
  width =20.0_dp
  height=20.0_dp
  a=4.385

  call slab(at, axes,  a, width, height, 1, Z)
  call randomise(at%pos, 0.3_dp)
  call set_cutoff(at, 4.0_dp)
  call calc_connect(at)

  call print_xyz(at, xyzfile)

  do i = 1,at%N
     call append(embedlist, (/i,0,0,0/))
     call BFS_grow(at, embedlist, 3, nneighb_only=.true., min_images_only=.true.)
     
     cluster = create_cluster(at, embedlist, terminate = .true., &
          periodic = (/.false., .false., .true./), even_hydrogens = .false., &
          vacuum_size = 10.0_dp)
     
     call print_xyz(cluster, xyzfile)

     call wipe(embedlist)
  end do

  call system_finalise()
end program
