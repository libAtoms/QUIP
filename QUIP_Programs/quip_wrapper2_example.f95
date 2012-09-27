program quip_wrapper2_example

  implicit none

  integer, parameter :: dp = selected_real_kind(10)
  integer, parameter :: n = 6
  real(dp), dimension(3,3) :: lattice
  integer, dimension(n) :: Z
  real(dp), dimension(3,n) :: coord
  character(len=10240) :: args_str

  real(dp) :: energy
  real(dp), dimension(3,n) :: force
  real(dp), dimension(3,3) :: stress

  lattice(1,1) = 6.0d0; lattice(1,2) =  0.0d0; lattice(1,3) =  0.0d0;
  lattice(2,1) =  0.0d0; lattice(2,2) = 6.0d0; lattice(2,3) =  0.0d0;
  lattice(3,1) =  0.0d0; lattice(3,2) =  0.0d0; lattice(3,3) = 6.0d0;

  Z = (/8,1,1,8,1,1/)

  coord(1,1) = -7.110371d0; coord(2,1) = -3.533572d0; coord(3,1) =  2.147261d0;
  coord(1,2) = -7.933029d0; coord(2,2) = -3.234956d0; coord(3,2) =  2.573383d0;
  coord(1,3) = -6.492180d0; coord(2,3) = -3.628907d0; coord(3,3) =  2.852075d0;
  coord(1,4) = -8.691981d0; coord(2,4) = -6.070463d0; coord(3,4) = -0.236430d0;
  coord(1,5) = -8.188371d0; coord(2,5) = -5.684074d0; coord(3,5) = -0.942021d0;
  coord(1,6) = -8.927810d0; coord(2,6) = -6.879853d0; coord(3,6) = -0.621913d0;
  
  coord = coord /6.0d0 ! fractional coords

  args_str="xml_label=water"

  call quip_wrapper2(n,Z, coord, lattice,args_str,0,energy,force,stress)

  print*,'Energy = ', energy
  print*,'Force  = ', force
  print*,'Stress = ', stress

endprogram quip_wrapper2_example
