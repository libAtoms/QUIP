program quip_wrapper_example

  implicit none

  integer, parameter :: n = 6
  real(8), dimension(3,3) :: lattice
  character(len=3), dimension(n) :: symbol
  real(8), dimension(3,n) :: coord
  character(len=1023) :: args_str
  integer :: args_str_length

  real(8) :: energy
  real(8), dimension(3,n) :: force
  real(8), dimension(3,3) :: virial

  lattice(1,1) = 20.0d0; lattice(1,2) =  0.0d0; lattice(1,3) =  0.0d0;
  lattice(2,1) =  0.0d0; lattice(2,2) = 20.0d0; lattice(2,3) =  0.0d0;
  lattice(3,1) =  0.0d0; lattice(3,2) =  0.0d0; lattice(3,3) = 20.0d0;

  symbol = (/"O  ","H  ","H  ","O  ","H  ","H  "/)

  coord(1,1) = -7.110371d0; coord(2,1) = -3.533572d0; coord(3,1) =  2.147261d0;
  coord(1,2) = -7.933029d0; coord(2,2) = -3.234956d0; coord(3,2) =  2.573383d0;
  coord(1,3) = -6.492180d0; coord(2,3) = -3.628907d0; coord(3,3) =  2.852075d0;
  coord(1,4) = -8.691981d0; coord(2,4) = -6.070463d0; coord(3,4) = -0.236430d0;
  coord(1,5) = -8.188371d0; coord(2,5) = -5.684074d0; coord(3,5) = -0.942021d0;
  coord(1,6) = -8.927810d0; coord(2,6) = -6.879853d0; coord(3,6) = -0.621913d0;

  args_str="IP PartridgeSchwenke force_using_fd"
  args_str_length = len_trim(args_str)

  
  call quip_wrapper(n,lattice,symbol,coord,args_str,args_str_length,energy,force,virial,.true.,.true.,.false.)
  
  
  print*,'Energy = ', energy
  print*, "forces:"
  print*, force
  
endprogram quip_wrapper_example
