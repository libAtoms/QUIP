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

module frametools_module
use system_module
use linearalgebra_module
use atoms_module
use quaternions_module
implicit none

!interface rotate
!  module procedure ft_rotate
!end interface rotate

contains

subroutine atoms_mark(this, mark_f, f_i_data, f_r_data)
  type(atoms), intent(inout) :: this
  interface
    function mark_f(at, i, i_data, r_data)
    use system_module
    use atoms_module
      type(atoms) :: at
      integer :: i
      integer, optional :: i_data(:)
      real(dp), optional :: r_data(:)
      logical mark_f
    end function mark_f
  end interface
  integer, intent(in), optional :: f_i_data(:)
  real(dp), intent(in), optional :: f_r_data(:)

  integer i
  integer, pointer :: mark(:)
  logical dummy

  call print ("called atoms_mark, this%N " // this%N)

  nullify(mark)
  if (.not. assign_pointer(this, "mark", mark)) then
    call add_property(this, "mark", 0)
    dummy = atoms_assign_pointer(this, "mark", mark)
  endif
  call print("associated(mark) " // associated(mark))
  call print("size(mark) " // size(mark))

  do i=1, this%N
    if (mark_f(this, i, f_i_data, f_r_data)) then
      mark(i) = 1
    endif
  end do
end subroutine atoms_mark


function is_in_cylinder(this, i, i_data, r_data)
  type(atoms), intent(in) :: this
  integer, intent(in) :: i
  integer, intent(in), optional :: i_data(:)
  real(dp), intent(in), optional :: r_data(:)
  logical is_in_cylinder

  real(dp) :: delta_p(3), v(3), r

  if (present(i_data)) call system_abort("is_in_cylinder doesn't take i_data")
  if (.not.present(r_data)) call system_abort("is_in_cylinder requires r_data")

  if (size(r_data) /= 7) call system_abort("is_in_cylinder requires [px py pz vx vy vz r] in cylinder")

  delta_p = this%pos(:,i) - r_data(1:3)
  v = r_data(4:6)
  v = v/norm(v)
  r = r_data(7)

  delta_p = delta_p - v*(v .dot. delta_p)
  is_in_cylinder = (norm(delta_p) < r)
end function is_in_cylinder

subroutine ft_rotate(field, axis, angle, origin)
  real(dp), intent(inout) :: field(:,:)
  real(dp), intent(in) :: axis(3)
  real(dp), intent(in) :: angle
  real(dp), intent(in), optional :: origin(3)

  real(dp) :: dr(3)
  type(Quaternion) :: Q
  integer i, N

  Q = rotation(axis, angle)
  
  N = size(field,2)

  do i=1, N
    if (present(origin)) then
      dr = field(:,i)-origin
    else
      dr = field(:,i)
    endif
    field(:,i) = rotate(dr, Q)
  end do
end subroutine ft_rotate


end module frametools_module
