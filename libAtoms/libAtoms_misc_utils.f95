module libatoms_misc_utils_module
use libatoms_module
implicit none
private

public :: dipole_moment

contains

function dipole_moment(pos, net_charge)
  real(dp), intent(in) :: pos(:,:), net_charge(:)
  real(dp) :: dipole_moment(size(pos,1))

  integer :: i

  if (size(pos,2) /= size(net_charge)) &
    call system_abort("dipole_moment called with incompatible size(pos) " // size(pos) //&
                      ", size(net_charge) " // size(net_charge))

  dipole_moment = 0.0_dp
  do i=1, size(pos,1)
    dipole_moment(:) = dipole_moment(:) - pos(:,i)*net_charge(i)
  end do
end function dipole_moment

end module libatoms_misc_utils_module
