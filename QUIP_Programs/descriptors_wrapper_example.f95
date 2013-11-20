program descriptors_wrapper_example

  implicit none

  integer, parameter :: n = 8
  real(8), dimension(3,3) :: lattice
  character(len=3), dimension(n) :: symbol
  real(8), dimension(3,n) :: coord
  character(len=10240) :: descriptor_str

  real(8), dimension(n,n) :: distances

  lattice(1,1) = 5.34224017d0; lattice(1,2) =        0.0d0; lattice(1,3) =        0.0d0;
  lattice(2,1) =        0.0d0; lattice(2,2) = 5.34224017d0; lattice(2,3) =        0.0d0;
  lattice(3,1) =        0.0d0; lattice(3,2) =        0.0d0; lattice(3,3) = 5.34224017d0;

  symbol = "Si "

  coord(1,1) =  -0.52654131d0; coord(2,1) =    -2.00406449d0; coord(3,1) =     1.19474806d0;
  coord(1,2) =   1.99352351d0; coord(2,2) =     0.42020741d0; coord(3,2) =     1.15114713d0;
  coord(1,3) =   1.21646828d0; coord(2,3) =     2.61257833d0; coord(3,3) =    -1.34327040d0;
  coord(1,4) =  -1.42817712d0; coord(2,4) =    -1.89404678d0; coord(3,4) =    -1.94681516d0;
  coord(1,5) =   0.10104488d0; coord(2,5) =    -0.15357825d0; coord(3,5) =    -0.33804307d0;
  coord(1,6) =   2.08539850d0; coord(2,6) =    -0.02403592d0; coord(3,6) =    -1.84183528d0;
  coord(1,7) =   1.79283079d0; coord(2,7) =    -2.48879868d0; coord(3,7) =     1.55589928d0;
  coord(1,8) =  -1.54238274d0; coord(2,8) =     1.54656737d0; coord(3,8) =    -0.12887274d0;

  descriptor_str ="soap n_max=10 l_max=10 atom_sigma=0.5 cutoff=2.7 cutoff_transition_width=0.5"

  call descriptors_wrapper(n,lattice,symbol,coord,descriptor_str,4,distances)

  print*,'Distances'
  print*, distances

  coord(1,4) =  -1.42817712d0 - 0.0221501781132952d0; coord(2,4) =    -1.89404678d0 + 0.00468815192049838d0; coord(3,4) =    -1.94681516d0 + 0.0457506835434298d0;
  call descriptors_wrapper(n,lattice,symbol,coord,descriptor_str,4,distances)

  print*,'Distances'
  print*, distances

  call descriptors_wrapper(n,lattice,symbol,coord,descriptor_str,0,distances)

  print*,'Distances'
  print*, distances

endprogram descriptors_wrapper_example
