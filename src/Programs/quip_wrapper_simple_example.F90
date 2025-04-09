program quip_wrapper_simple_example

  implicit none

  integer, parameter :: n = 2
  real(8), dimension(3,3) :: lattice
  integer, dimension(n) :: Z
  real(8), dimension(3,n) :: coord


  real(8) :: energy
  real(8), dimension(3,n) :: force
  real(8), dimension(3,3) :: virial



  lattice(1,1) = 20.0d0; lattice(1,2) =  0.0d0; lattice(1,3) =  0.0d0;
  lattice(2,1) =  0.0d0; lattice(2,2) = 20.0d0; lattice(2,3) =  0.0d0;
  lattice(3,1) =  0.0d0; lattice(3,2) =  0.0d0; lattice(3,3) = 20.0d0;

  Z = (/29,29/)

  coord(1,1) = -7.110371d0; coord(2,1) = -3.533572d0; coord(3,1) =  2.147261d0;
  coord(1,2) = -7.933029d0; coord(2,2) = -3.234956d0; coord(3,2) =  2.573383d0;
 
   
  call quip_wrapper_simple(n,lattice,Z,coord,energy,force,virial)
  
  
  print*,'Energy = ', energy
  print*, "forces:"
  print*, force
  print*, "virial:"
  print*, virial
endprogram quip_wrapper_simple_example
