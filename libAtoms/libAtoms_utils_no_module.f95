subroutine lattice_abc_to_xyz(cell_lengths, cell_angles, lattice)
use system_module
use atoms_module
implicit none
   real(dp), intent(in) :: cell_lengths(3), cell_angles(3)
   real(dp), intent(out) :: lattice(3,3)

   lattice = make_lattice(cell_lengths(1), cell_lengths(2), cell_lengths(3), &
			  cell_angles(1), cell_angles(2), cell_angles(3))

end subroutine lattice_abc_to_xyz

subroutine lattice_xyz_to_abc(lattice, cell_lengths, cell_angles)
use system_module
use atoms_module
implicit none
   real(dp), intent(in) :: lattice(3,3)
   real(dp), intent(out) :: cell_lengths(3), cell_angles(3)

   call get_lattice_params(lattice, cell_lengths(1), cell_lengths(2), cell_lengths(3), &
				    cell_angles(1), cell_angles(2), cell_angles(3))

end subroutine lattice_xyz_to_abc
