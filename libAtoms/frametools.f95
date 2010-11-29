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

#include "error.inc"
module frametools_module
use system_module
use linearalgebra_module
use atoms_module
use quaternions_module
implicit none


public :: atoms_mark

public :: mark_cylinder, mark_sphere

private :: ft_rotate
public :: rotate
interface rotate
   module procedure ft_rotate
end interface rotate

contains

subroutine atoms_mark(this, mark_f, f_i_data, f_r_data, periodic, mark_name, mark_value, intersection, error)
  type(atoms), intent(inout) :: this
  interface
    function mark_f(at, i, periodic, i_data, r_data, error)
    use system_module
    use atoms_module
      type(atoms), intent(inout) :: at
      integer, intent(in) :: i
      logical, intent(in) :: periodic
      integer, intent(in), optional :: i_data(:)
      real(dp), intent(in), optional :: r_data(:)
      integer, intent(out), optional :: error
      logical mark_f
    end function mark_f
  end interface
  integer, intent(in), optional :: f_i_data(:)
  real(dp), intent(in), optional :: f_r_data(:)
  logical, intent(in), optional :: periodic
  character(len=*), intent(in), optional :: mark_name
  integer, intent(in), optional :: mark_value
  logical, intent(in), optional :: intersection
  integer, intent(out), optional :: error

  integer :: i
  character(len=128) :: my_mark_name
  integer :: my_mark_value
  logical :: my_intersection, my_periodic
  integer, pointer :: mark(:)

  INIT_ERROR(error)

  my_mark_name = optional_default("mark", mark_name)
  my_mark_value = optional_default(1, mark_value)
  my_intersection = optional_default(.false., intersection)
  my_periodic = optional_default(.true., periodic)

  if (.not. assign_pointer(this, trim(my_mark_name), mark)) then
      call add_property(this, trim(my_mark_name), 0, ptr=mark)
  else
      call assign_property_pointer(this, trim(my_mark_name), mark)
  endif

  do i=1, this%N
    if (my_intersection .and. mark(i) == 0) cycle ! if we're doing intersection and this one isn't already marked, skip it
    if (my_intersection) mark(i) = 0 ! if we're doing intersection, set mark to .false., and it'll will be set to true if it's marked in the current shape
    if (mark_f(this, i, my_periodic, f_i_data, f_r_data, error=error)) then
      mark(i) = my_mark_value
    endif
    PASS_ERROR(error)
  end do
end subroutine atoms_mark

function is_in_cylinder_f(this, i, periodic, i_data, r_data, error)
  type(atoms), intent(inout) :: this
  integer, intent(in) :: i
  logical, intent(in) :: periodic
  integer, intent(in), optional :: i_data(:)
  real(dp), intent(in), optional :: r_data(:)
  integer, intent(out), optional :: error
  logical is_in_cylinder_f

  logical :: my_periodic
  real(dp) :: delta_p(3), v(3), r

  INIT_ERROR(error)

  if (present(i_data)) then
     RAISE_ERROR("is_in_cylinder_f doesn't take i_data", error)
  endif
  if (.not.present(r_data)) then
     RAISE_ERROR("is_in_cylinder_f requires r_data", error)
  endif
  if (size(r_data) /= 7) then
     RAISE_ERROR("is_in_cylinder_f requires [px py pz vx vy vz r] in cylinder", error)
  endif

  if (periodic) then
     delta_p = diff_min_image(this, i, r_data(1:3))
  else
     delta_p = this%pos(:,i) - r_data(1:3)
  endif
  v = r_data(4:6)
  v = v/norm(v)
  r = r_data(7)

  delta_p = delta_p - v*(v .dot. delta_p)
  if (r > 0.0_dp) then
     is_in_cylinder_f = (norm(delta_p) < r)
  else
     is_in_cylinder_f = (norm(delta_p) >= -r)
  endif
end function is_in_cylinder_f

function is_in_sphere_f(this, i, periodic, i_data, r_data, error)
  type(atoms), intent(inout) :: this
  integer, intent(in) :: i
  logical, intent(in) :: periodic
  integer, intent(in), optional :: i_data(:)
  real(dp), intent(in), optional :: r_data(:)
  integer, intent(out), optional :: error
  logical is_in_sphere_f

  real(dp) :: dr, r

  INIT_ERROR(error)

  if (present(i_data)) then
     RAISE_ERROR("is_in_sphere_f doesn't take i_data", error)
  endif
  if (.not.present(r_data)) then
     RAISE_ERROR("is_in_sphere_f requires r_data", error)
  endif
  if (size(r_data) /= 4) then
     RAISE_ERROR("is_in_sphere_f requires [px py pz r] in sphere", error)
  endif

  if (periodic) then
     dr = distance_min_image(this, i, r_data(1:3))
  else
     dr = norm(this%pos(:,i)-r_data(1:3))
  endif
  r = r_data(4)

  if (r > 0.0_dp) then
     is_in_sphere_f = (dr < r)
  else
     is_in_sphere_f = (dr >= -r)
  endif
end function is_in_sphere_f

subroutine mark_cylinder(this, p, v, r, periodic, mark_name, mark_value, intersection)
  type(atoms), intent(inout) :: this
  real(dp), intent(in) :: p(3), v(3), r
  character(len=*), intent(in), optional :: mark_name
  integer, intent(in), optional :: mark_value
  logical, intent(in), optional :: periodic, intersection

  call atoms_mark(this, is_in_cylinder_f, f_r_data=(/ p, v, r /), periodic=periodic, mark_name=mark_name, mark_value=mark_value, intersection=intersection)
end subroutine mark_cylinder

subroutine mark_sphere(this, p, r, periodic, mark_name, mark_value, intersection)
  type(atoms), intent(inout) :: this
  real(dp), intent(in) :: p(3), r
  character(len=*), intent(in), optional :: mark_name
  integer, intent(in), optional :: mark_value
  logical, intent(in), optional :: periodic, intersection

  call atoms_mark(this, is_in_sphere_f, f_r_data=(/ p, r /), periodic=periodic, mark_name=mark_name, mark_value=mark_value, intersection=intersection)
end subroutine mark_sphere

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
