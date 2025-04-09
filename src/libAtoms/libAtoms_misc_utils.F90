! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   Portions of this code were written by
! H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! H0 X
! H0 X   Copyright 2006-2010.
! H0 X
! H0 X   These portions of the source code are released under the GNU General
! H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
! H0 X
! H0 X   If you would like to license the source code under different terms,
! H0 X   please contact Gabor Csanyi, gabor@csanyi.net
! H0 X
! H0 X   Portions of this code were written by Noam Bernstein as part of
! H0 X   his employment for the U.S. Government, and are not subject
! H0 X   to copyright in the USA.
! H0 X
! H0 X
! H0 X   When using this software, please cite the following reference:
! H0 X
! H0 X   http://www.libatoms.org
! H0 X
! H0 X  Additional contributions by
! H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! H0 X
! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

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
