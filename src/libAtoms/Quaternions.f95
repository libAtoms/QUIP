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

!% Implementation of quaternion algebra. Primarily for use with 
!% rigid body dynamics, although the 'Orientation', 'Rotation' 
!% and 'Rotate' routines can be useful in many other situations.
!%
!% \begin{displaymath}
!% \begin{array}{c}
!% i^2 = j^2 = k^2 = ijk = -1 \\ ij = k \quad jk = i \quad ki = j
!% \end{array}
!% \end{displaymath}

module quaternions_module

  use system_module        ! for dp
  use units_module         ! for PI
  use linearalgebra_module ! for .feq./.fne

  implicit none
  private

  public :: quaternion, initialise, finalise, print
  public :: rotate, rotation, rotation_matrix, orientation, assignment(=)
  public :: rotation_parameters
  public :: operator(.dot.), operator(*), operator(+), operator(.conj.)

  type Quaternion

     real(dp) :: a = 0.0_dp   !
     real(dp) :: b = 0.0_dp   !
     real(dp) :: c = 0.0_dp   ! = a*1 + b*i + c*j + d*k
     real(dp) :: d = 0.0_dp   !

  end type Quaternion

  private :: quaternion_initialise
  interface initialise
     module procedure quaternion_initialise
  end interface initialise

  private :: quaternion_finalise
  interface finalise
     module procedure quaternion_finalise
  end interface finalise

  private :: quaternion_norm
  interface norm
     module procedure quaternion_norm
  end interface norm

  private :: quaternion_normsq
  interface normsq
     module procedure quaternion_normsq
  end interface normsq

  private :: quat_plus_quat, quat_plus_vect
  interface operator(+)
     module procedure quat_plus_quat, quat_plus_vect
  end interface

  private :: quat_minus_quat, quat_minus_vect
  interface operator(-)
     module procedure quat_minus_quat, quat_minus_vect
  end interface

  private :: quat_mult_real, quat_mult_quat, real_mult_quat
  interface operator(*)
     module procedure quat_mult_real, quat_mult_quat, real_mult_quat
  end interface

  private :: quat_divide_real, quat_divide_quat
  interface operator(/)
     module procedure quat_divide_real, quat_divide_quat
  end interface

  private :: quat_assign_vect, vect_assign_quat, quat_assign_real
  interface assignment(=)
     module procedure quat_assign_vect, vect_assign_quat, quat_assign_real
  end interface

  private :: quaternion_conjugate
  interface operator(.conj.)
     module procedure quaternion_conjugate
  end interface

  private :: quat_eq_quat, quat_eq_vect
  interface operator(.feq.)
     module procedure quat_eq_quat, quat_eq_vect
  end interface

  private :: quat_ne_quat, quat_ne_vect
  interface operator(.fne.)
     module procedure quat_ne_quat, quat_ne_vect
  end interface

  private :: quat_dot_quat
  interface operator(.dot.)
     module procedure quat_dot_quat
  end interface

  private :: rotate_vect, rotate_quat
  interface rotate
     module procedure rotate_vect, rotate_quat
  end interface rotate

  private :: quaternion_print
  interface print
     module procedure quaternion_print
  end interface print

contains

  subroutine quaternion_initialise(this, a, b, c, d)
    type(Quaternion), intent(out) :: this
    real(dp), intent(in), optional :: a, b, c, d
    
    this%a = optional_default(0.0_dp, a)
    this%b = optional_default(0.0_dp, b)
    this%c = optional_default(0.0_dp, c)
    this%d = optional_default(0.0_dp, d)

  end subroutine quaternion_initialise

  subroutine quaternion_finalise(this)
    type(Quaternion), intent(inout) :: this

    this%a = 0.0_dp
    this%b = 0.0_dp
    this%c = 0.0_dp    
    this%d = 0.0_dp
    
  end subroutine quaternion_finalise

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X Assignment
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine quat_assign_vect(this,vect)

    type(Quaternion),       intent(out) :: this
    real(dp), dimension(:), intent(in)  :: vect
    
    select case(size(vect))

       case(3)
          this%a = 0.0_dp
          this%b = vect(1)
          this%c = vect(2)
          this%d = vect(3)

       case(4)
          this%a = vect(1)
          this%b = vect(2)
          this%c = vect(3)
          this%d = vect(4)

       case default
          call system_abort('Quat_Assign_Vect: Vector must have 3 or 4 components')

    end select

  end subroutine quat_assign_vect

  subroutine vect_assign_quat(vect,q)

    real(dp), dimension(:), intent(out) :: vect
    type(Quaternion),       intent(in)  :: q

    select case(size(vect))

       case(3)
          vect = (/q%b,q%c,q%d/)

       case(4)
          vect = (/q%a,q%b,q%c,q%d/)

       case default
          call system_abort('Vect_Assign_Quat: Vector must have 3 or 4 components')
          
    end select

  end subroutine vect_assign_quat

  subroutine quat_assign_real(q,r)
    
    type(Quaternion), intent(out) :: q
    real(dp),         intent(in)  :: r

    q%a = r
    q%b = 0.0_dp
    q%c = 0.0_dp
    q%d = 0.0_dp

  end subroutine quat_assign_real

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X Simple arithmetic
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function quat_plus_quat(q1,q2) result(res)

    type(Quaternion), intent(in) :: q1,q2
    type(Quaternion)             :: res

    res%a = q1%a + q2%a
    res%b = q1%b + q2%b
    res%c = q1%c + q2%c
    res%d = q1%d + q2%d

  end function quat_plus_quat

  function quat_plus_vect(q,v) result(res)

    type(Quaternion),       intent(in) :: q
    real(dp), dimension(4), intent(in) :: v
    type(Quaternion)                   :: res

    res%a = q%a + v(1)
    res%b = q%b + v(2)
    res%c = q%c + v(3)
    res%d = q%d + v(4)

  end function quat_plus_vect

  function quat_minus_quat(q1,q2) result(res)

    type(Quaternion), intent(in) :: q1,q2
    type(Quaternion)             :: res

    res%a = q1%a - q2%a
    res%b = q1%b - q2%b
    res%c = q1%c - q2%c
    res%d = q1%d - q2%d

  end function quat_minus_quat

  function quat_minus_vect(q,v) result(res)

    type(Quaternion),       intent(in) :: q
    real(dp), dimension(4), intent(in) :: v
    type(Quaternion)                   :: res

    res%a = q%a - v(1)
    res%b = q%b - v(2)
    res%c = q%c - v(3)
    res%d = q%d - v(4)

  end function quat_minus_vect

  function quat_mult_quat(q1,q2) result(res)

    type(Quaternion), intent(in) :: q1,q2
    type(Quaternion)             :: res

    res%a = q1%a * q2%a  -  q1%b * q2%b  -  q1%c * q2%c - q1%d * q2%d
    res%b = q1%a * q2%b  +  q1%b * q2%a  +  q1%c * q2%d - q1%d * q2%c
    res%c = q1%a * q2%c  -  q1%b * q2%d  +  q1%c * q2%a + q1%d * q2%b
    res%d = q1%a * q2%d  +  q1%b * q2%c  -  q1%c * q2%b + q1%d * q2%a

  end function quat_mult_quat

  function quat_mult_real(q,r) result(res)

    type(Quaternion), intent(in) :: q
    real(dp),         intent(in) :: r
    type(Quaternion)             :: res

    res%a = q%a * r
    res%b = q%b * r
    res%c = q%c * r
    res%d = q%d * r

  end function quat_mult_real

  function real_mult_quat(r,q) result(res)

    type(Quaternion), intent(in) :: q
    real(dp),         intent(in) :: r
    type(Quaternion)             :: res

    res%a = q%a * r
    res%b = q%b * r
    res%c = q%c * r
    res%d = q%d * r

  end function real_mult_quat

  function quat_divide_real(q,r) result(res)

    type(Quaternion), intent(in) :: q
    real(dp),         intent(in) :: r
    type(Quaternion)             :: res

    res%a = q%a / r
    res%b = q%b / r
    res%c = q%c / r
    res%d = q%d / r

  end function quat_divide_real

  function quat_divide_quat(q1,q2) result(res)

    type(Quaternion), intent(in) :: q1,q2
    type(Quaternion)             :: res, q2c
    real(dp)                     :: q2n

    q2c = .conj. q2
    q2n = normsq(q2)

    res = (q1 * q2c) / q2n

  end function quat_divide_quat

  !Dot product as if they were 4 component vectors
  function quat_dot_quat(q1,q2) result(res)

    type(Quaternion), intent(in) :: q1,q2
    real(dp)                     :: res

    res = q1%a * q2%a  +  q1%b * q2%b  +  q1%c * q2%c + q1%d * q2%d

  end function quat_dot_quat

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X Comparison
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function quat_eq_quat(q1,q2) result(res)

    type(Quaternion), intent(in) :: q1,q2
    logical                      :: res

    res = .false.

    if (q1%a .fne. q2%a) return
    if (q1%b .fne. q2%b) return
    if (q1%c .fne. q2%c) return
    if (q1%d .fne. q2%d) return

    res = .true.

  end function quat_eq_quat
  
  function quat_ne_quat(q1,q2) result(res)

    type(Quaternion), intent(in) :: q1,q2
    logical                      :: res

    res = .not. (q1 .feq. q2)

  end function quat_ne_quat

  function quat_eq_vect(q,v) result(res)

    type(Quaternion), intent(in)       :: q
    real(dp), dimension(:), intent(in) :: v
    logical                            :: res
    type(Quaternion)                   :: qv

    qv = v

    res = .false.

    if (q%a .fne. qv%a) return
    if (q%b .fne. qv%b) return
    if (q%c .fne. qv%c) return
    if (q%d .fne. qv%d) return

    res = .true.

  end function quat_eq_vect

  function quat_ne_vect(q,v) result(res)

    type(Quaternion),       intent(in) :: q
    real(dp), dimension(:), intent(in) :: v
    logical                            :: res

    res = .not. (q .feq. v)

  end function quat_ne_vect

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X Norms
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function quaternion_normsq(this) result(normsq)

    type(Quaternion), intent(in) :: this
    real(dp)                     :: normsq

    normsq =   this%a*this%a &
            + this%b*this%b &
            + this%c*this%c &
            + this%d*this%d

  end function quaternion_normsq

  function quaternion_norm(this) result(norm)

    type(Quaternion), intent(in) :: this
    real(dp)                     :: norm

    norm = sqrt(quaternion_normsq(this))

  end function quaternion_norm

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X Quaternion conjugate
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  
  function quaternion_conjugate(this) result(res)

    type(Quaternion), intent(in) :: this
    type(Quaternion)             :: res

    res%a =  this%a
    res%b = -this%b
    res%c = -this%c
    res%d = -this%d

  end function quaternion_conjugate

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X Printing
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine quaternion_print(this,file)

    type(Quaternion), intent(in) :: this
    type(Inoutput),  optional, intent(in) :: file
    !local variables
    character                    :: s2,s3,s4

    s2 = '+'; s3 = '+'; s4 = '+'
    if (this%b < 0.0_dp) s2 = '-'
    if (this%c < 0.0_dp) s3 = '-'
    if (this%d < 0.0_dp) s4 = '-'
    
    write(line,'(4(f0.5,a))')     this%a,' '//s2//' ',   &
                              abs(this%b),'i '//s3//' ', &
                              abs(this%c),'j '//s4//' ', &
                              abs(this%d),'k'
    call Print(line,file=file)

  end subroutine quaternion_print


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X Rotation routines
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% Many rotations can be collapsed into one, e.g.
  !% to rotate using $q_1$ then $q_2$ :
  !% rotate using $q_3 = q_2 \times q_1$
  !%
  !%
  !% Construct a rotation quaternion from axis and angle (in a righthanded sense)
  !
  function rotation(axis,angle) result(res)

    real(dp), dimension(3), intent(in) :: axis
    real(dp),               intent(in) :: angle  ! in radians
    type(Quaternion)                   :: res
    !local variables
    real(dp)                           :: n
    
    n = norm(axis)
    
    res = (/ cos(angle/2.0_dp), sin(angle/2.0_dp)*axis/n /)

  end function rotation

  !
  !% Rotate a vector using the given quaternion
  !
  subroutine rotate_vect(v,qr)

    real(dp), dimension(3), intent(inout) :: v
    type(Quaternion),       intent(in)    :: qr
    !local variables
    type(Quaternion)                      :: qv

    qv = v
    call rotate_quat(qv,qr)
    v = qv

  end subroutine rotate_vect

  !
  !% Rotate a vector represented as a quaternion ('qv'), using a given rotation quaternion ('qr')
  !% This assumes that the rotation quaternion is properly normalised
  !%
  !% 'qv = qr * qv * (.conj.qr)'
  subroutine rotate_quat(qv,qr)

    type(Quaternion), intent(inout)  :: qv
    type(Quaternion), intent(in)     :: qr

    qv = qr * qv * (.conj. qr)

  end subroutine rotate_quat

  !
  !% Given two vectors ('a1' and 'b1'), calculate the rotation quaternion 
  !% that must be used to rotate them into two other vectors ('a2' and 'b2').
  !% The angle between 'a1' and 'b1' must be the same as 'a2' and 'b2',
  !% unless 'correct_angle' is present and true
  !
  function orientation(a1, b1, a2, b2, correct_angle)

    real(dp), dimension(3), intent(in) :: a1,b1,a2,b2
    logical, optional,      intent(in) :: correct_angle
    type(Quaternion)                   :: Orientation
    !local variables
    type(Quaternion)                   :: r1, r2 !The two parts of the orientation
    real(dp), dimension(3)             :: n, p, q, a1_hat, b1_hat, a2_hat, b2_hat, b3_hat, zero3
    real(dp)                           :: theta
    logical                            :: do_correct_angle

    do_correct_angle = .false.
    if (present(correct_angle)) do_correct_angle = correct_angle
    zero3 = (/0.0_dp,0.0_dp,0.0_dp/)

    !Normalise the vectors
    a1_hat = a1 / norm(a1)
    b1_hat = b1 / norm(b1)
    a2_hat = a2 / norm(a2)
    b2_hat = b2 / norm(b2)

    !Check the angles a1-b1 and a2-b2 match
    if (.not. do_correct_angle .and. ((a1_hat .dot. b1_hat) .fne. (a2_hat .dot. b2_hat))) then
       call system_abort('Orientation: Angle a1--b1 is not equal to angle a2--b2')
    end if

    !Correct the a1--b1 angle to be the same as the a2--b2 angle by moving b1
    if (do_correct_angle) then
       n = b1_hat - (a1_hat .dot. b1_hat)*a1_hat !Perpendicular to a1 in a1--b1 plane
       n = n / norm(n)
       theta = angle(a2_hat, b2_hat)
       b1_hat = a1_hat * cos(theta) + n * sin(theta)
    end if

    !Now construct the rotation quaternion needed to take a1 to a2 around an axis perpendicular to both
    n = a1_hat .cross. a2_hat

    !Special case: The vectors could be almost (anti)parallel, in which case it we can choose any vector
    !perpendicular to them as the first rotation axis
    if (n .feq. zero3) then
       do while(n .feq. zero3)       !This is for the TINY possiblility that the random vector we choose
                                     !is STILL (anti)parallel to a1
       !Create a random unit vector
          p = random_unit_vector()
          n = a1_hat .cross. p
       end do
    end if

    n = n / norm(n)
    
    theta = angle(a1_hat,a2_hat)
    r1 = Rotation(n,theta)

    !Acting on b1 with r1 will take it to a new position, b3
    b3_hat = b1_hat
    call Rotate(b3_hat,r1)

    !Now figure out the angle which we need to rotate around a2 (or -a2) to take b3 to b2
    ! -> Construct vectors perpendicular to a2, in the a2-b3 and a2-b2 planes, then find angle between them   
    p = b3_hat - (a2_hat .dot. b3_hat)*a2_hat
    p = p / norm(p)
    q = b2_hat - (a2_hat .dot. b2_hat)*a2_hat
    q = q / norm(q)
    n = p .cross. q   ! This will be either in the a2 or the -a2 direction

    !Another Special case: Miraculously, p and q could be (anti) parallel, in which case we can just
    !use a2_hat, then is doesn't matter if we rotate  by 0 (or PI) clockwise or anticlockwise
    if (n .feq. zero3) n = a2_hat

    n = n / norm(n)
    theta = angle(p,q)
    r2 = rotation(n,theta)

    !Combine the two rotations into one using quaternion multiplication
    orientation = r2 * r1 

  end function orientation

  !
  !% Returns the angle and rotation axis of a (normalised) quaternion
  !
  subroutine rotation_parameters(q,theta,axis)
    
    type(Quaternion),       intent(in)  :: q
    real(dp),               intent(out) :: theta
    real(dp), dimension(3), intent(out) :: axis

    theta = 2.0_dp * acos(q%a)

    axis = (/q%b,q%c,q%d/) / (sqrt(1 - q%a*q%a))

  end subroutine rotation_parameters

  !
  !% Return the equivalent rotation matrix of a (unit) quaternion
  !
  function rotation_matrix(q_in) result(A)
    
    type(Quaternion), intent(in) :: q_in
    real(dp), dimension(3,3)     :: A
    type(Quaternion)             :: q

    !Normalise q_in
    q = q_in / norm(q_in)

    !Now build A
    A(1,1) = q%a*q%a+q%b*q%b-q%c*q%c-q%d*q%d; A(1,2) = 2.0_dp*(q%b*q%c-q%a*q%d);        A(1,3) = 2.0_dp*(q%b*q%d+q%a*q%c)
    A(2,1) = 2.0_dp*(q%b*q%c+q%a*q%d);        A(2,2) = q%a*q%a-q%b*q%b+q%c*q%c-q%d*q%d; A(2,3) = 2.0_dp*(q%c*q%d-q%a*q%b)
    A(3,1) = 2.0_dp*(q%b*q%d-q%a*q%c);        A(3,2) = 2.0_dp*(q%c*q%d+q%a*q%b);        A(3,3) = q%a*q%a-q%b*q%b-q%c*q%c+q%d*q%d

  end function rotation_matrix

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X Testing
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !%OMIT
  subroutine quaternion_test

    type(Quaternion) :: q1, q2, q3
    real(dp)         :: v3(3), v4(4), a(3), b(3), c(3), d(3), M(3,3)

    !Assignment
    call Print('Testing assignment...')
    call Print('q1 <- 4vector:')
    q1 = (/1.0_dp,1.0_dp,1.0_dp,1.0_dp/)
    call Print(q1)
    call Print('q2 <- 3vector:')
    q2 = (/2.0_dp,2.0_dp,2.0_dp/)
    call Print(q2)
    call Print('q1 -> 3vector:')
    v3 = q1
    call Print(v3)
    call Print('q2 -> 4vector:')
    v4 = q2
    call Print(v4)
    call Print('')

    !Algebra
    call Print('Testing alegbra...')
    call Print('q1 + q2:')
    q3 = q1 + q2
    call Print(q3)
    call Print('q1 - q2:')
    q3 = q1 - q2
    call Print(q3)
    call Print('q1 * q2:')
    q3 = q1 * q2
    call Print(q3)
    call Print('q1 / q2:')
    q3 = q1 / q2
    call Print(q3)
    call Print('')

    !Norms
    call Print('Testing norms...')
    call Print('normsq(q1):')
    call Print(normsq(q1))
    call Print('normsq(q2):')
    call Print(normsq(q2))
    call Print('norm(q1):')
    call Print(norm(q1))
    call Print('norm(q2):')
    call Print(norm(q2))
    call Print('')

    !Conjugate
    call Print('Testing conjugates...')
    call Print('q1*:')
    call Print(.conj.q1)
    call Print('q2*:')
    call Print(.conj.q2)
    call Print('')

    !Rotation
    call Print('Testing Rotations...')
    call Print('Setting q1 to x axis...')
    q1 = (/1.0_dp,0.0_dp,0.0_dp/)
    call Print('q1:')
    call Print(q1)
    call Print('Setting up q2 for rotation about z axis by 90deg...')
    q2 = Rotation(unit_vector(0.0_dp,0.0_dp),PI/2.0_dp)
    call Print('q2:')
    call Print(q2)
    call Print('Rotating...')
    call Rotate(q1,q2)
    q3 = q1
    call Print(q3)
    v3 = q3
    call Print('Rotated vector:')
    call Print(v3)
    call Print('')
    call Print('Rotating this by a further 90deg around z axis')
    call Rotate(v3,q2)
    call Print('Rotated vector:')
    call Print(v3)
    call Print('')
    call Print('Creating random rotation quaternion...')
    q1 = Rotation(random_unit_vector(),ran_uniform()*PI)
    call Print(q1)
    call Print('Calculating equivalent rotation matrix...')
    M = Rotation_matrix(q1)
    call Print(M)
    call Print('Creating random vector...')
    v3 = random_unit_vector()
    call Print(v3)
    call Print('Rotating with quaternion...')
    a = v3
    call Rotate(a,q1)
    call Print(a)
    call Print('Rotating with matrix...')
    b = M .mult. v3
    call Print(b)
    call Print('')    

    !Orientation
    !NOTE: remember that q and -q represent the same rotation!
    call Print('Testing orientation...')
    call Print('Creating two random unit vectors...')
    a = random_unit_vector()
    b = random_unit_vector()
    call Print('a:')
    call Print(a)
    call Print('b:')
    call Print(b)
    call Print('Creating random rotation quaternion...')
    q1 = Rotation(random_unit_vector(),ran_uniform()*PI)
    call Print(q1)
    call Print('Rotating random vectors with random quaternion...')
    c = a
    call Rotate(c,q1)
    d = b
    call Rotate(d,q1)
    call Print('a -> c:')
    call Print(c)
    call Print('b -> d:')
    call Print(b)
    call Print('Calculating rotation quaternion based on a,b,c,d:')
    q2 = Orientation(a,b,c,d,.true.)
    call Print(q2)
    call Print('')

  end subroutine quaternion_test

end module quaternions_module
