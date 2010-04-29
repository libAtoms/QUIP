module qw_continous_module

  use libatoms_module

  implicit none

  real(dp), parameter :: factorial_table(0:16) = (/ &
       1_dp, &
       1_dp, &
       2_dp, &
       6_dp, &
       24_dp, &
       120_dp, &
       720_dp, &
       5040_dp, &
       40320_dp, &
       362880_dp, &
       3628800_dp, &
       39916800_dp, &
       479001600_dp, &
       6227020800_dp, &
       87178291200_dp, &
       1307674368000_dp, &
       20922789888000_dp /)

  interface calc_qw
     module procedure calc_qw_atom, calc_qw_position
  end interface

  private

  public :: calc_qw_at
  public :: calc_qw_array
  public :: test_qw_grad
  public :: calc_qw
  public :: SolidHarmonicsCartesian, SphericalHarmonicsCartesian, GradSphericalHarmonicsCartesian
  public :: WeightFunction, GradWeightFunction
  public :: factorial, oscillate, tc, wigner3j

  contains

  function calc_qw_array(array, l, cutoff, wf_char_length, wf_type, do_weight)

    real(dp), intent(in) :: array(:,:)                ! 3,N array of atomic positions
    integer, intent(in) :: l                          ! degree of spherical harmonics used
    real(dp), intent(in) :: cutoff                    ! value of the radial cutoff
    real(dp), intent(in), optional :: wf_char_length  ! value of the characteristic length of the weight function
    integer, intent(in), optional :: wf_type          ! type of the weight function
    logical, intent(in), optional :: do_weight        ! include normalisation
    real(dp) :: calc_qw_array(2,size(array, 2))       ! 2,N array of qw parameters

    type(Atoms) :: at
    real(dp), pointer :: at_q(:), at_w(:)

    call initialise(at, size(array, 2), reshape((/ 100.0_dp * (maxval(array(1,:)) - minval(array(1,:))),0.0_dp,0.0_dp, &
                                                   0.0_dp,100.0_dp * (maxval(array(2,:)) - minval(array(2,:))),0.0_dp, &
                                                   0.0_dp,0.0_dp,100.0_dp * (maxval(array(3,:)) - minval(array(3,:))) /), (/3,3/)))
    at%pos = array
    call set_cutoff(at, cutoff)
    call calc_connect(at)

    call calc_qw_at(at, l = l, cutoff = cutoff, wf_char_length = wf_char_length, wf_type = wf_type, do_weight = do_weight)

    if (.not. assign_pointer(at, 'q' // l, at_q)) call system_abort('calc_qw_at: could not assign pointer to atoms object')
    if (.not. assign_pointer(at, 'w' // l, at_w)) call system_abort('calc_qw_at: could not assign pointer to atoms object')

    calc_qw_array(1,:) = at_q
    calc_qw_array(2,:) = at_w

    call finalise(at)

  end function calc_qw_array

  subroutine calc_qw_at(at, l, cutoff, wf_char_length, wf_type, do_weight)

    type(Atoms), intent(inout) :: at
    integer, intent(in), optional :: l
    real(dp), intent(in), optional :: cutoff
    real(dp), intent(in), optional :: wf_char_length
    integer, intent(in), optional :: wf_type
    logical, intent(in), optional :: do_weight

    integer :: my_l
    real(dp), pointer :: at_q(:), at_w(:)
    integer :: i

    my_l = optional_default(4, l)

    call add_property(at, 'q' // my_l, 0.0_dp)
    call add_property(at, 'w' // my_l, 0.0_dp)

    if (.not. assign_pointer(at, 'q' // my_l, at_q)) call system_abort('calc_qw_at: could not assign pointer to atoms object')
    if (.not. assign_pointer(at, 'w' // my_l, at_w)) call system_abort('calc_qw_at: could not assign pointer to atoms object')

    do i = 1, at%N
       call calc_qw(at, i, my_l, q = at_q(i), w = at_w(i), cutoff = cutoff, wf_char_length = wf_char_length, wf_type = wf_type, do_weight = do_weight)
    end do

  end subroutine calc_qw_at

  function test_qw_grad(at, x0, l, cutoff, wf_char_length, wf_type, do_weight)

    logical :: test_qw_grad
    type(Atoms), intent(inout) :: at
    real(dp), intent(in), optional :: x0(3)
    integer, intent(in), optional :: l
    real(dp), intent(in), optional :: cutoff
    real(dp), intent(in), optional :: wf_char_length
    integer, intent(in), optional :: wf_type
    logical, intent(in), optional :: do_weight

    real(dp) :: my_x0(3)
    integer :: my_l

    my_x0 = optional_default((/0.0_dp, 0.0_dp, 0.0_dp/), x0)
    my_l = optional_default(4, l)

    test_qw_grad = test_gradient(my_x0, q, grad_q) .and. test_gradient(my_x0, w, grad_w)

  contains

    function q(x, data)

      real(dp) :: q
      real(dp) :: x(:)
      character, optional :: data(:)

      call calc_qw(at, x, my_l, q = q, cutoff = cutoff, wf_char_length = wf_char_length, wf_type = wf_type, do_weight = do_weight)

    endfunction q

    function w(x, data)

      real(dp) :: w
      real(dp) :: x(:)
      character, optional :: data(:)

      call calc_qw(at, x, my_l, w = w, cutoff = cutoff, wf_char_length = wf_char_length, wf_type = wf_type, do_weight = do_weight)

    endfunction w

    function grad_q(x, data)

      real(dp) :: grad_q(3)
      real(dp) :: x(:)
      character, optional :: data(:)

      call calc_qw(at, x, my_l, grad_q = grad_q, cutoff = cutoff, wf_char_length = wf_char_length, wf_type = wf_type, do_weight = do_weight)

    endfunction grad_q

    function grad_w(x, data)

      real(dp) :: grad_w(3)
      real(dp) :: x(:)
      character, optional :: data(:)

      call calc_qw(at, x, my_l, grad_w = grad_w, cutoff = cutoff, wf_char_length = wf_char_length, wf_type = wf_type, do_weight = do_weight)

    endfunction grad_w

  endfunction test_qw_grad

  subroutine calc_qw_atom(at, i, l, q, w, grad_q, grad_w, d, cutoff, wf_char_length, wf_type, do_weight)

    type(Atoms), intent(in) :: at
    integer, intent(in) :: i, l
    real(dp), intent(out), optional :: q, w
    real(dp), intent(out), optional :: grad_q(3), grad_w(3)
    integer, intent(in), optional :: d
    real(dp), intent(in), optional :: cutoff
    real(dp), intent(in), optional :: wf_char_length
    integer, intent(in), optional :: wf_type
    logical, intent(in), optional :: do_weight

    integer :: my_d
    real(dp) :: my_cutoff
    real(dp) :: my_wf_char_length
    integer :: my_wf_type
    logical :: my_do_weight
    integer :: n, j, m, m1, m2, m3
    real(dp) :: x_ij(3)
    complex(dp), allocatable :: qm(:), grad_qm(:,:)
    real(dp) :: weights, grad_weights(3)
    real(dp) :: qm_norm2_sum, grad_qm_norm2_sum(3)
    complex(dp) :: wc, grad_wc(3)
    real(dp) :: temporary_weight, temporary_gradweight(3)

    if (.not. at%connect%initialised) &
       call system_abort('calc_qw_atom: atoms structure has no connectivity data - call calc_connect first')

    my_wf_char_length = optional_default(1.0_dp, wf_char_length)
    my_wf_type = optional_default(1, wf_type)

    if ((my_wf_char_length .feq. 0.0_dp) .and. (present(grad_q) .or. present(grad_w))) &
       call system_abort('calc_qw_atom: qw not differentiable - cannot calculate gradient')

    if (my_wf_char_length .feq. 0.0_dp) then
       my_do_weight = optional_default(.true., do_weight)
    else
       my_do_weight = optional_default(.false., do_weight)
    end if

    if (my_do_weight .and. (present(grad_q) .or. present(grad_w))) &
       call system_abort('calc_qw_atom: qw not differentiable - cannot calculate gradient')

    my_cutoff = optional_default(at%cutoff, cutoff)

    if (my_cutoff > at%cutoff) &
       call system_abort('calc_qw_atom: weight function cutoff greater than connectivity cutoff distance - call calc_connect first')

    allocate(qm(-l:l))

    qm = CPLX_ZERO
    weights = 0.0_dp

    do n = 1, atoms_n_neighbours(at, i)
       j = atoms_neighbour(at, i, n, diff = x_ij, max_dist = my_cutoff)

       if (j /= 0) then
          temporary_weight = WeightFunction(my_cutoff, x_ij, my_wf_char_length, my_wf_type)

          do m = -l, l
             qm(m) = qm(m) + (SphericalHarmonicsCartesian(l, m, x_ij) * temporary_weight)
          end do

          if (my_do_weight) weights = weights + temporary_weight
       end if
    end do

    if (my_do_weight) then
       if (weights .feq. 0.0_dp) then
          qm = CPLX_ZERO
       else
          qm = qm / weights
       end if
    end if

    qm_norm2_sum = sum(real(qm * conjg(qm), dp))

    if (present(q)) q = sqrt(4.0_dp * PI * qm_norm2_sum / ((2.0_dp * l) + 1.0_dp))

    if (present(w) .or. present(grad_w)) then
       wc = CPLX_ZERO

       do m1 = -l, l
          do m2 = -l, l
             do m3 = -l, l
                if ((m1 + m2 + m3) /= 0 ) cycle

                wc = wc + (wigner3j(l, m1, l, m2, l, m3) * qm(m1) * qm(m2) * qm(m3))
             end do
          end do
       end do

       if (qm_norm2_sum .feq. 0.0_dp) then
          wc = CPLX_ZERO
       else
          wc = wc / (qm_norm2_sum**1.5_dp)
       end if

       if (present(w)) w = real(wc, dp)
    end if

    if (present(grad_q) .or. present(grad_w)) then
       my_d = optional_default(i, d)

       allocate(grad_qm(3,-l:l))

       grad_qm = CPLX_ZERO
       grad_weights = 0.0_dp

       if (my_d == i) then
          do n = 1, atoms_n_neighbours(at, i)
             j = atoms_neighbour(at, i, n, diff = x_ij, max_dist = my_cutoff)

             if (j /= 0) then
                temporary_weight = WeightFunction(my_cutoff, x_ij, my_wf_char_length, my_wf_type)
                temporary_gradweight = GradWeightFunction(my_cutoff, x_ij, my_wf_char_length, my_wf_type)

                do m = -l, l
                   grad_qm(:,m) = grad_qm(:,m) - (GradSphericalHarmonicsCartesian(l, m, x_ij) * temporary_weight) &
                                               - (SphericalHarmonicsCartesian(l, m, x_ij) * temporary_gradweight)
                end do

                if (my_do_weight) grad_weights = grad_weights - temporary_gradweight
             end if
          end do
       else
          x_ij = diff_min_image(at, i, my_d)

          if (norm(x_ij) <= my_cutoff) then
             do m = -l, l
                temporary_weight = WeightFunction(my_cutoff, x_ij, my_wf_char_length, my_wf_type)
                temporary_gradweight = GradWeightFunction(my_cutoff, x_ij, my_wf_char_length, my_wf_type)

                grad_qm(:,m) = grad_qm(:,m) + (GradSphericalHarmonicsCartesian(l, m, x_ij) * temporary_weight) &
                                            + (SphericalHarmonicsCartesian(l, m, x_ij) * temporary_gradweight)
             end do

             if (my_do_weight) grad_weights = grad_weights + temporary_gradweight
          end if
       end if

       if (my_do_weight) then
          if (weights .feq. 0.0_dp) then
             grad_qm = CPLX_ZERO
          else
             do j = 1, 3
                grad_qm(j,:) = (grad_qm(j,:) - (grad_weights(j) * qm)) / weights
             end do
          end if
       end if

       do j = 1, 3
          grad_qm_norm2_sum(j) = sum(real((qm * conjg(grad_qm(j,:))) + (grad_qm(j,:) * conjg(qm)), dp))
       end do

       if (present(grad_q)) then
          if (qm_norm2_sum .feq. 0.0_dp) then
             grad_q = 0.0_dp
          else
             grad_q = 2.0_dp * PI * grad_qm_norm2_sum / sqrt(4.0_dp * PI * qm_norm2_sum * ((2.0_dp * l) + 1.0_dp))
          end if
       end if

       if (present(grad_w)) then
          grad_wc = CPLX_ZERO

          do m1 = -l, l
             do m2 = -l, l
                do m3 = -l, l
                   if ((m1 + m2 + m3) /= 0 ) cycle

                   grad_wc = grad_wc + (wigner3j(l, m1, l, m2, l, m3) * ((grad_qm(:,m1) * qm(m2) * qm(m3)) &
                                                                       + (qm(m1) * grad_qm(:,m2) * qm(m3)) &
                                                                       + (qm(m1) * qm(m2) * grad_qm(:,m3))))
                end do
             end do
          end do

          if (qm_norm2_sum .feq. 0.0_dp) then
             grad_wc = CPLX_ZERO
          else
             grad_wc = (grad_wc / (qm_norm2_sum**1.5_dp)) - (1.5_dp * wc * grad_qm_norm2_sum / qm_norm2_sum)
          end if

          grad_w = real(grad_wc, dp)
       end if

       deallocate(grad_qm)
    end if

    deallocate(qm)

  end subroutine calc_qw_atom

  subroutine calc_qw_position(at, x, l, q, w, grad_q, grad_w, d, cutoff, wf_char_length, wf_type, do_weight)

    type(Atoms), intent(inout) :: at
    real(dp), intent(in) :: x(3)
    integer, intent(in) :: l
    real(dp), intent(out), optional :: q, w
    real(dp), intent(out), optional :: grad_q(3), grad_w(3)
    integer, intent(in), optional :: d
    real(dp), intent(in), optional :: cutoff
    real(dp), intent(in), optional :: wf_char_length
    integer, intent(in), optional :: wf_type
    logical, intent(in), optional :: do_weight

    call add_atoms(at, x, 1)
    call calc_connect(at)

    call calc_qw_atom(at, at%N, l, q = q, w = w, grad_q = grad_q, grad_w = grad_w, d = d, cutoff = cutoff, wf_char_length = wf_char_length, wf_type = wf_type, do_weight = do_weight)

    call remove_atoms(at, at%N)
    call calc_connect(at)

  end subroutine calc_qw_position

  function SolidHarmonicsCartesian(l, m, x)

    complex(dp) :: SolidHarmonicsCartesian
    integer, intent(in) :: l, m
    real(dp), intent(in) :: x(3)
    integer :: p, q, s

    SolidHarmonicsCartesian = CPLX_ZERO

    do p = 0, l
       q = p - m
       s = l - p - q

       if ((q >= 0) .and. (s >= 0)) then
          SolidHarmonicsCartesian = SolidHarmonicsCartesian + ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**p) &
                                                             * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**q) &
                                                             * (x(3)**s) &
                                                             / (factorial(p) * factorial(q) * factorial(s)))
       end if
    end do

    SolidHarmonicsCartesian = SolidHarmonicsCartesian * sqrt(factorial(l + m) * factorial(l - m) * ((2.0_dp * l) + 1) / (4.0_dp * PI))

  end function SolidHarmonicsCartesian

  function SphericalHarmonicsCartesian(l, m, x)

    complex(dp) :: SphericalHarmonicsCartesian
    integer, intent(in) :: l, m
    real(dp), intent(in) :: x(3)

    SphericalHarmonicsCartesian = SolidHarmonicsCartesian(l, m, x) * (norm2(x)**(-0.5_dp * l))

  end function SphericalHarmonicsCartesian

  function GradSphericalHarmonicsCartesian(l, m, x)

    complex(dp) :: GradSphericalHarmonicsCartesian(3)
    integer, intent(in) :: l, m
    real(dp), intent(in) :: x(3)
    integer :: p, q, s

    GradSphericalHarmonicsCartesian = CPLX_ZERO

    do p = 0, l
       q = p - m
       s = l - p - q

       if ((p >= 1) .and. (q >= 0) .and. (s >= 0)) then
          GradSphericalHarmonicsCartesian(1) = GradSphericalHarmonicsCartesian(1) - ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**(p - 1)) &
                                                                                   * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**q) &
                                                                                   * (x(3)**s) &
                                                                                   * 0.5_dp &
                                                                                   / (factorial(p - 1) * factorial(q) * factorial(s)))
          GradSphericalHarmonicsCartesian(2) = GradSphericalHarmonicsCartesian(2) - ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**(p - 1)) &
                                                                                   * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**q) &
                                                                                   * (x(3)**s) &
                                                                                   * 0.5_dp * cmplx(0.0_dp, 1.0_dp, dp) &
                                                                                   / (factorial(p - 1) * factorial(q) * factorial(s)))
       end if

       if ((p >= 0) .and. (q >= 1) .and. (s >= 0)) then
          GradSphericalHarmonicsCartesian(1) = GradSphericalHarmonicsCartesian(1) + ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**p) &
                                                                                   * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**(q - 1)) &
                                                                                   * (x(3)**s) &
                                                                                   * 0.5_dp &
                                                                                   / (factorial(p) * factorial(q - 1) * factorial(s)))
          GradSphericalHarmonicsCartesian(2) = GradSphericalHarmonicsCartesian(2) - ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**p) &
                                                                                   * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**(q - 1)) &
                                                                                   * (x(3)**s) &
                                                                                   * 0.5_dp * cmplx(0.0_dp, 1.0_dp, dp) &
                                                                                   / (factorial(p) * factorial(q - 1) * factorial(s)))
       end if

       if ((p >= 0) .and. (q >= 0) .and. (s >= 1)) then
          GradSphericalHarmonicsCartesian(3) = GradSphericalHarmonicsCartesian(3) + ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**p) &
                                                                                   * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**q) &
                                                                                   * (x(3)**(s - 1)) &
                                                                                   / (factorial(p) * factorial(q) * factorial(s - 1)))
       end if
    end do

    GradSphericalHarmonicsCartesian = GradSphericalHarmonicsCartesian * sqrt(factorial(l + m) * factorial(l - m) * ((2.0_dp * l) + 1) / (4.0_dp * PI)) &
                                                                      * (norm2(x)**(-0.5_dp * l))

    GradSphericalHarmonicsCartesian = GradSphericalHarmonicsCartesian - (l * x * SphericalHarmonicsCartesian(l, m, x) / norm2(x))

  end function GradSphericalHarmonicsCartesian

  function WeightFunction(cutoff, x, char_length, type)

    real(dp) :: WeightFunction
    real(dp), intent(in) :: cutoff, x(3), char_length
    integer, intent(in) :: type
    real(dp) :: r

    r = norm(x)

    if (type == 1) then
       if ((r >= 0.0_dp) .and. (r < (cutoff - char_length))) then
          WeightFunction = 1.0_dp
       else if ((r >= (cutoff - char_length)) .and. (r <= cutoff)) then
          WeightFunction = (cos(0.5_dp * PI * (r - cutoff + char_length) / char_length))**2
       else
          call system_abort('WeightFunction: distance greater than cutoff')
       end if
    else if (type == 2) then
       if ((r >= 0.0_dp) .and. (r < (cutoff - (2.0_dp * char_length)))) then
          WeightFunction = 0.0_dp
       else if ((r >= (cutoff - (2.0_dp * char_length))) .and. (r <= cutoff)) then
          WeightFunction = (cos(0.5_dp * PI * (r - cutoff + char_length) / char_length))**2
       else
          call system_abort('WeightFunction: distance greater than cutoff')
       end if
    else if (type == 3) then
       if ((r >= 0.0_dp) .and. (r < (cutoff - char_length))) then
          WeightFunction = (cos(0.5_dp * PI * (r - cutoff + char_length) / (cutoff - char_length)))**2
       else if ((r >= (cutoff - char_length)) .and. (r <= cutoff)) then
          WeightFunction = (cos(0.5_dp * PI * (r - cutoff + char_length) / char_length))**2
       else
          call system_abort('WeightFunction: distance greater than cutoff')
       end if
    else
       call system_abort('WeightFunction: weightfunction type unknown')
    end if

  endfunction WeightFunction

  function GradWeightFunction(cutoff, x, char_length, type)

    real(dp) :: GradWeightFunction(3)
    real(dp), intent(in) :: cutoff, x(3), char_length
    integer, intent(in) :: type
    real(dp) :: r

    r = norm(x)

    if (type == 1) then
       if ((r >= 0.0_dp) .and. (r < (cutoff - char_length))) then
          GradWeightFunction = 0.0_dp
       else if ((r >= (cutoff - char_length)) .and. (r <= cutoff)) then
          GradWeightFunction = -PI * x * cos(0.5_dp * PI * (r - cutoff + char_length) / char_length) &
                                       * sin(0.5_dp * PI * (r - cutoff + char_length) / char_length) &
                                       / (char_length * r)
       else
          call system_abort('GradWeightFunction: distance greater than cutoff')
       end if
    else if (type == 2) then
       if ((r >= 0.0_dp) .and. (r < (cutoff - (2.0_dp * char_length)))) then
          GradWeightFunction = 0.0_dp
       else if ((r >= (cutoff - (2.0_dp * char_length))) .and. (r <= cutoff)) then
          GradWeightFunction = -PI * x * cos(0.5_dp * PI * (r - cutoff + char_length) / char_length) &
                                       * sin(0.5_dp * PI * (r - cutoff + char_length) / char_length) &
                                       / (char_length * r)
       else
          call system_abort('GradWeightFunction: distance greater than cutoff')
       end if
    else if (type == 3) then
       if ((r >= 0.0_dp) .and. (r < (cutoff - char_length))) then
          GradWeightFunction = -PI * x * cos(0.5_dp * PI * (r - cutoff + char_length) / (cutoff - char_length)) &
                                       * sin(0.5_dp * PI * (r - cutoff + char_length) / (cutoff - char_length)) &
                                       / ((cutoff - char_length) * r)
       else if ((r >= (cutoff - char_length)) .and. (r <= cutoff)) then
          GradWeightFunction = -PI * x * cos(0.5_dp * PI * (r - cutoff + char_length) / char_length) &
                                       * sin(0.5_dp * PI * (r - cutoff + char_length) / char_length) &
                                       / (char_length * r)
       else
          call system_abort('GradWeightFunction: distance greater than cutoff')
       end if
    else
       call system_abort('WeightFunction: weightfunction type unknown')
    end if

  end function GradWeightFunction

    function factorial(n) result(res)

      ! factorial_real

      integer, intent(in) :: n
      real(dp)            :: res
      integer :: i

      if (n<0) then
         call system_abort('factorial: negative argument')
      elseif(n <= 16) then
         res = factorial_table(n)
      else
         res=1.0_dp
         do i=2,n
            res = res*i
         end do
      end if

    endfunction factorial

    function oscillate(m)

      integer, intent(in) :: m
      integer :: oscillate

      if( mod(m,2) == 0 ) then
          oscillate = 1
      else
          oscillate = -1
      endif

    endfunction oscillate

    function tc(a,b,c)

       ! Triangle coefficient
       ! http://mathworld.wolfram.com/TriangleCoefficient.html
       !
       ! \[
       ! \Delta(a,b,c) = \frac{ (a+b-c)! (a-b+c)! (-a+b+c)! }{ (a+b+c+1)! }
       ! \]

       real(dp) :: tc
       integer, intent(in) :: a, b, c

       tc = factorial(a+b-c) * factorial(a-b+c) * factorial(-a+b+c) / factorial(a+b+c+1)

    endfunction tc

    function wigner3j(j1,m1,j2,m2,j,m)

        ! Wigner 3J symbol
        ! Source: http://mathworld.wolfram.com/Wigner3j-Symbol.html
        ! 
        ! \[
        ! \left( \begin{array}{ccc}
        ! j_1 & j_2 & j \\
        ! m_1 & m_2 & m \\
        ! \end{array} \right) = (-1)^{j_1-j_2-m) \sqrt{ \Delta (j_1,j_2,j) }
        ! \sqrt{ (j_1+m_1)! (j_1-m_1)! (j_2+m_2)! (j_2-m_2)! (j+m)! (j-m)! }
        ! \sum_k \frac{ (-1)^k }{k! (j-j_2+k+m_1)! (j-j_1+k-m_2)! (j_1+j_2-j-k)!
        ! (j_1-k-m_1)! (j_2-k+m_2)! }
        ! \]
        ! the summation index k runs on all integers where none of the argument of
        ! factorials are negative
        ! $\Delta(a,b,c)$ is the triangle coefficient.

        real(dp)            :: wigner3j
        integer, intent(in) :: j1, m1, j2, m2, j, m

        real(dp) :: pre_fac, triang_coeff, main_coeff, sum_coeff, sum_term
        integer  :: k, kmin, kmax

        pre_fac = oscillate(j1-j2-m)

        triang_coeff = sqrt( tc(j1,j2,j) )

        main_coeff = sqrt( &
                   & factorial(j1+m1) * factorial(j1-m1) * &
                   & factorial(j2+m2) * factorial(j2-m2) * &
                   & factorial(j+m) * factorial(j-m) )
                   
        sum_coeff = 0.0_dp

        kmin = max( j2-j-m1, j1+m2-j, 0 )
        kmax = min( j1+j2-j, j1-m1, j2+m2)

        do k = kmin, kmax

           sum_term = 1.0_dp / ( factorial(k) * factorial(j-j2+k+m1) * &
                    & factorial(j-j1+k-m2) * factorial(j1+j2-j-k)    * &
                    & factorial(j1-k-m1) * factorial(j2-k+m2) )

           sum_term = oscillate(k) * sum_term
           !if (mod(k,2)==1) sum_term = -sum_term
           
           sum_coeff = sum_coeff + sum_term

        enddo

        wigner3j = pre_fac * triang_coeff * main_coeff * sum_coeff

    endfunction wigner3j

end module qw_continous_module
