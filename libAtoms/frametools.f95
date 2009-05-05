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
