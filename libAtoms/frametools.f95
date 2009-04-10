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



subroutine elastic_fields(at, a, C11, C12, C44)
  type(Atoms), intent(inout) :: at
  real(dp), intent(in) :: a, C11, C12, C44

  real(dp) :: C(6,6), strain(6), stress(6)
  real(dp), dimension(3) :: n1,n2,n3, d
  real(dp), dimension(3,3) :: rotXYZ, E123, EEt, V, S, R, SS, Sig, RSigRt, RtE, SigEvecs
  integer :: i, j, m, nn, ngood

  real(dp), pointer, dimension(:) :: S_xx_sub1, S_yy_sub1, S_zz_sub1, S_yz, S_xz, S_xy, &
       Sig_xx, Sig_yy, Sig_zz, Sig_yz, Sig_xz, Sig_xy, SigEval1, SigEval2, SigEval3, &
       von_mises_stress, von_mises_strain
  
  real(dp), pointer, dimension(:,:) :: SigEvec1, SigEvec2, SigEvec3
  logical :: dum
  real(dp), dimension(4) :: dist
  real(dp), dimension(3,4) :: ndiff


  rotXYZ = 0.0_dp
  rotXYZ(1,2) = 1.0_dp
  rotXYZ(2,3) = 1.0_dp
  rotXYZ(3,1) = 1.0_dp

  ! Elastic constants matrix C_{ij}
  C = 0.0_dp
  C(1,1) = C11; C(2,2) = C11; C(3,3) = C11;
  C(4,4) = C44; C(5,5) = C44; C(6,6) = C44;
  C(1,2) = C12; C(1,3) = C12; C(2,3) = C12;
  C(2,1) = C12; C(3,1) = C12; C(3,2) = C12;


  ! Create properties and assign pointers

  call add_property(at, 'S_xx_sub1', 0.0_dp)
  call add_property(at, 'S_yy_sub1', 0.0_dp)
  call add_property(at, 'S_zz_sub1', 0.0_dp)
  call add_property(at, 'S_yz', 0.0_dp)
  call add_property(at, 'S_xz', 0.0_dp)
  call add_property(at, 'S_xy', 0.0_dp)

  call add_property(at, 'Sig_xx', 0.0_dp)
  call add_property(at, 'Sig_yy', 0.0_dp)
  call add_property(at, 'Sig_zz', 0.0_dp)
  call add_property(at, 'Sig_yz', 0.0_dp)
  call add_property(at, 'Sig_xz', 0.0_dp)
  call add_property(at, 'Sig_xy', 0.0_dp)

  call add_property(at, 'SigEval1', 0.0_dp)
  call add_property(at, 'SigEval2', 0.0_dp)
  call add_property(at, 'SigEval3', 0.0_dp)
  call add_property(at, 'SigEvec1', 0.0_dp, n_cols=3)
  call add_property(at, 'SigEvec2', 0.0_dp, n_cols=3)
  call add_property(at, 'SigEvec3', 0.0_dp, n_cols=3)

  call add_property(at, 'von_mises_stress', 0.0_dp)
  call add_property(at, 'von_mises_strain', 0.0_dp)

  dum = assign_pointer(at, 'S_xx_sub1', S_xx_sub1)
  dum = assign_pointer(at, 'S_yy_sub1', S_yy_sub1)
  dum = assign_pointer(at, 'S_zz_sub1', S_zz_sub1)
  dum = assign_pointer(at, 'S_yz', S_yz)
  dum = assign_pointer(at, 'S_xz', S_xz)
  dum = assign_pointer(at, 'S_xy', S_xy)

  dum = assign_pointer(at, 'Sig_xx', Sig_xx)
  dum = assign_pointer(at, 'Sig_yy', Sig_yy)
  dum = assign_pointer(at, 'Sig_zz', Sig_zz)
  dum = assign_pointer(at, 'Sig_yz', Sig_yz)
  dum = assign_pointer(at, 'Sig_xz', Sig_xz)
  dum = assign_pointer(at, 'Sig_xy', Sig_xy)

  dum = assign_pointer(at, 'SigEval1', SigEval1)
  dum = assign_pointer(at, 'SigEval2', SigEval2)
  dum = assign_pointer(at, 'SigEval3', SigEval3)
  dum = assign_pointer(at, 'SigEvec1', SigEvec1)
  dum = assign_pointer(at, 'SigEvec2', SigEvec2)
  dum = assign_pointer(at, 'SigEvec3', SigEvec3)

  dum = assign_pointer(at, 'von_mises_stress', von_mises_stress)
  dum = assign_pointer(at, 'von_mises_strain', von_mises_strain)

  call calc_connect(at)

  ! Loop over all atoms
  ngood = 0
  do i=1,at%N

     call print('Atom '//i)

     ! Count number of nearest neighbours
     nn = 0
     do m=1,atoms_n_neighbours(at, i) 
        if (is_nearest_neighbour(at,i,m)) then
           nn = nn + 1
           j = atoms_neighbour(at, i, m, distance=dist(nn), diff=ndiff(:,nn))

           if (at%Z(j) == 1) then ! Skip hydrogen neighbours
              nn = nn - 1
              cycle
           end if

           call print(nn//': '//dist(nn)//ndiff(:,nn), VERBOSE)
           if (nn > 5) exit
        end if
     end do

     if (nn == 4) then

        ngood = ngood + 1

        ! Find cubic axes from neighbours
        n1 = ndiff(:,2)-ndiff(:,1)
        n2 = ndiff(:,3)-ndiff(:,1)
        n3 = ndiff(:,4)-ndiff(:,1)

        call print('n1 '//norm(n1)//n1, VERBOSE)
        call print('n2 '//norm(n2)//n2, VERBOSE)
        call print('n3 '//norm(n3)//n3, VERBOSE)

        e123(:,1) = (n1 + n2 - n3)/a
        e123(:,2) = (n2 + n3 - n1)/a
        e123(:,3) = (n3 + n1 - n2)/a

        ! Kill near zero elements
        where (abs(e123) < 1.0e-6_dp) e123 = 0.0_dp

        if (all(e123 < 0.0_dp)) e123 = -e123

        call print('det(e) = '//matrix3x3_det(e123), VERBOSE)
        call print(e123, VERBOSE)

        if (matrix3x3_det(e123) < 0) then
           e123(:,3) = -e123(:,3)
           call print('reflected axis 3', VERBOSE)
           call print(e123, VERBOSE)
        end if

        ! Find polar decomposition: e123 = S*R where S is symmetric, 
        ! and R is a rotation
        !
        !  EEt = E*E', EEt = VDV' D diagonal, S = V D^1/2 V', R = S^-1*E

        EEt = e123 .mult. transpose(e123) ! Normal
        call diagonalise(EEt, D, V)

        ! Check positive definite
        if (any(D < 0)) then
           call print(e123, VERBOSE)
           call system_abort("EE' is not positive definite")
        end if

        S = V .mult. diag(sqrt(D)) .mult. transpose(V)
        R = V .mult. diag(D ** (-0.5_dp)) .mult. transpose(V) .mult. e123

        call print('S:', VERBOSE); call print(S, VERBOSE)
        call print('R:', VERBOSE); call print(R, VERBOSE)

        RtE = transpose(R) .mult. e123

        ! Check for permutations - which way does x point?
        if (RtE(2,1) > RtE(1,1) .and. RtE(2,1) > RtE(3,1)) then
           ! y direction
           R = rotXYZ .mult. R
        else if (RtE(3,1) > RtE(1,1) .and. RtE(3,1) > RtE(2,1)) then
           ! z direction
           R = transpose(rotXYZ) .mult. R
        end if

        SS = transpose(R) .mult. S .mult. R

        call print('R:', VERBOSE);    call print(R, VERBOSE)
        call print('RtE:', VERBOSE);  call print(RtE, VERBOSE)
        call print('RtSR:', VERBOSE); call print(SS, VERBOSE)

        ! Strain(1:6) = (/eps11,eps22,eps33,eps12,eps13,eps23/)
        strain(1) = SS(1,1) - 1.0_dp
        strain(2) = SS(2,2) - 1.0_dp
        strain(3) = SS(3,3) - 1.0_dp
        strain(4) = SS(1,2)
        strain(5) = SS(1,3)
        strain(6) = SS(2,3)
     else
        strain = 0.0_dp
        S = 0.0_dp
        S(1,1) = 1.0_dp; S(2,2) = 1.0_dp; S(3,3) = 1.0_dp
     end if

     stress = C .mult. strain

     ! Now stress(1:6) = (/sig11,sig22,sig33,sig12,sig13,sig23/)

     sig = 0.0_dp
     sig(1,1) = stress(1)
     sig(2,2) = stress(2)
     sig(3,3) = stress(3)
     sig(1,2) = stress(4)
     sig(1,3) = stress(5)
     sig(2,3) = stress(6)
     sig(2,1) = stress(4)
     sig(3,1) = stress(5)
     sig(3,2) = stress(6)

     von_mises_stress(i) = sqrt(0.5*((stress(1) - stress(2))**2.0_dp + &
          (stress(2) - stress(3))**2.0_dp + &
          (stress(3) - stress(1))**2.0_dp))

     von_mises_strain(i) = sqrt(strain(6)**2.0_dp + strain(5)**2.0_dp + &
          strain(4)**2.0_dp + &
          1.0_dp/6.0_dp*((strain(2) - strain(3))**2.0_dp + &
          (strain(1) - strain(3))**2.0_dp + &
          (strain(1) - strain(2))**2.0_dp))

     RSigRt = R .mult. sig .mult. transpose(R)

     call symmetrise(RSigRt)
     call diagonalise(RSigRt,D,SigEvecs)

     S_xx_sub1(i) = S(1,1) - 1.0_dp
     S_yy_sub1(i) = S(2,2) - 1.0_dp
     S_zz_sub1(i) = S(3,3) - 1.0_dp
     S_yz(i)      = S(2,3)
     S_xz(i)      = S(1,3)
     S_xy(i)      = S(1,2)

     Sig_xx(i)    = RSigRt(1,1)
     Sig_yy(i)    = RSigRt(2,2)
     Sig_zz(i)    = RSigRt(3,3)
     Sig_yz(i)    = RSigRt(2,3)
     Sig_xz(i)    = RSigRt(1,3)
     Sig_xy(i)    = RSigRt(1,2)

     SigEval1(i)  = D(1)
     SigEval2(i)  = D(2)
     SigEval3(i)  = D(3)

     SigEvec1(:,i) = SigEvecs(:,1)
     SigEvec2(:,i) = SigEvecs(:,2)
     SigEvec3(:,i) = SigEvecs(:,3)
  end do

  call print('Processed '//ngood//' of '//at%N//' atoms.')

end subroutine elastic_fields

end module frametools_module
