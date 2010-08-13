subroutine lattice_abc_to_xyz(cell_lengths, cell_angles, lattice)
use system_module
use units_module
use linearalgebra_module
implicit none
   real(dp), intent(in) :: cell_lengths(3), cell_angles(3)
   real(dp), intent(out) :: lattice(3,3)

   real(dp) a, b, c, alpha, beta, gamma, cos_alpha, cos2_alpha, cos_beta, cos2_beta, &
      cos_gamma, cos2_gamma, sin_gamma, sin2_gamma

   a = cell_lengths(1);  b = cell_lengths(2); c = cell_lengths(3)
   alpha = cell_angles(1)*PI/180.0_dp; beta = cell_angles(2)*PI/180.0_dp; gamma = cell_angles(3)*PI/180.0_dp;    
   cos_alpha = cos(alpha); cos2_alpha = cos_alpha*cos_alpha
   cos_beta  = cos(beta);  cos2_beta  = cos_beta *cos_beta
   cos_gamma = cos(gamma); cos2_gamma = cos_gamma*cos_gamma
   sin_gamma = sin(gamma); sin2_gamma = sin_gamma*sin_gamma

   lattice(1,1) = a
   lattice(2,1) = 0.0_dp
   lattice(3,1) = 0.0_dp

   lattice(1,2) = b * cos_gamma
   lattice(2,2) = b * sin_gamma
   lattice(3,2) = 0.0_dp

   lattice(1,3) = c * cos_beta
   lattice(2,3) = c * (cos_alpha - cos_beta*cos_gamma) / sin_gamma
   lattice(3,3) = c * sqrt(1.0_dp - (cos2_alpha + cos2_beta - 2.0*cos_alpha*cos_beta*cos_gamma)/ sin2_gamma)

end subroutine lattice_abc_to_xyz

subroutine lattice_xyz_to_abc(lattice, cell_lengths, cell_angles)
use system_module
use units_module
use linearalgebra_module
implicit none
   real(dp), intent(in) :: lattice(3,3)
   real(dp), intent(out) :: cell_lengths(3), cell_angles(3)

   cell_lengths = sqrt(sum(lattice**2,1))
   cell_angles(1) = 180.0_dp/PI * acos( (lattice(:,2) .dot. lattice(:,3)) / (cell_lengths(2)*cell_lengths(3)) )
   cell_angles(2) = 180.0_dp/PI * acos( (lattice(:,1) .dot. lattice(:,3)) / (cell_lengths(1)*cell_lengths(3)) )
   cell_angles(3) = 180.0_dp/PI * acos( (lattice(:,2) .dot. lattice(:,1)) / (cell_lengths(2)*cell_lengths(1)) )

end subroutine lattice_xyz_to_abc
