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

!% Types and routines for rigid body dynamics.
!%
!% \textbf{The implementation of the NO_SQUISH algorithm in AdvanceVerlet still requires
!% some work, so rigid bodies should be constructed from constraints for the time
!% being.}
!
module rigidbody_module

  use system_module
  use periodictable_module
  use atoms_module
  use linearalgebra_module
  use quaternions_module

  implicit none

  type RigidBodyModel
  !% Object used as a reference by the RigidBody type.
  !% Many RigidBodys can point to the same model
     
     type(Atoms)            :: at                   !% The actual atoms
     integer                :: d1a, d1b, d2a, d2b   !% Atoms used as reference directions, 'd1a$>$d1b' and 'd2a$>$d2b'
     real(dp), dimension(3) :: I                    !% Moments of inertia along x,y,z axes (after reorientation)
     real(dp)               :: Mtot                 !% The total mass of the body
     logical                :: initialised = .false.

  end type RigidBodyModel


  type RigidBody
  !% Holds one rigid body (e.g. one water molecule)

     integer                               :: N = 0   !% Number of atoms in the rigid body
     integer,  allocatable, dimension(:)   :: indices !% Atoms that make up the rigid body
     type(RigidBodyModel), pointer         :: model   !% A model of the rigid body in its correct conformation
     type(Quaternion)                      :: q       !% The orientation of the body
     type(Quaternion)                      :: p       !% The momentum of the orientation
     real(dp), dimension(3)                :: RCoM    !% The centre of mass position
     real(dp), dimension(3)                :: VCoM    !% The centre of mass velocity
     logical                               :: initialised = .false.

  end type RigidBody

  integer, parameter :: M_ROT = 10 !Increasing this number integrates the rigid rotations more accurately

  interface initialise
     module procedure rigidbody_initialise
  end interface initialise
  
  interface finalise
     module procedure rigidbody_finalise, rigidbodies_finalise
  end interface finalise

  interface assignment(=)
     module procedure rigidbody_assignment, rigidbodies_assignment
  end interface

  interface print
     module procedure rigidbodymodel_print, rigidbody_print, rigidbodies_print
  end interface print

  private P0, P1, P2, P3

contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X INITIALISATION / FINALISATION
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !
  !% Initialise a RigidBodyModel: 
  !% \begin{itemize}
  !% \item Calculates total mass.
  !% \item Shifts the centre of mass to the origin.
  !% \item Calculates moments of inertia and rotates
  !%       the atoms so that the principal axes.
  !%       lie along cartesian directions.
  !% \item Finds two reference directions which are
  !%       not (anti)parallel.
  !% \end{itemize}
  !
  subroutine rigidbodymodel_initialise(this,at)

    type(RigidBodyModel), intent(inout) :: this
    type(Atoms),          intent(in)    :: at
    !local variables
    real(dp), dimension(3,3)            :: I_matrix, I_vects
    real(dp), dimension(3)              :: XCoM, d1, d2
    real(dp)                            :: m
    integer                             :: i
    type(Quaternion)                    :: q

    !A rigid body must have at least 2 atoms
    if (at%N < 3) call system_abort('RigidBodyModel_Initialise: A rigid body must have at least 3 atoms')

    if (this%initialised) call rigidbodymodel_finalise(this)

    !Copy the atomic info
    call initialise(this%at,at%N,at%lattice)
    this%at%Z = at%Z
    this%at%pos = at%pos

    !Shift the model so that the centre of mass is at the origin (note: no periodic wrapping is assumed)
    XCoM = 0.0_dp
    this%Mtot = 0.0_dp

    do i = 1, this%at%N
       m = this%at%mass(i)
       XCoM = XCoM + m * this%at%pos(:,i)
       this%Mtot = this%Mtot + m
    end do

    XCoM = XCoM / this%Mtot

    do i = 1, this%at%N
       this%at%pos(:,i) = this%at%pos(:,i) - XCoM
    end do

    !Calculate the model's inertia tensor
    I_matrix = Inertia_Tensor(this%at)
    call diagonalise(I_matrix,this%I,I_vects) !Stores the moments of inertia in this%I

    !Rotate the model so that the principle axes point along the x,y,z directions.
    !The positions are then stored in body-fixed coordinates
    !Eigenvalues are already normalised and orthogonal.
    q = Orientation(I_vects(:,1),            I_vects(:,2), &                 !Go from where the eigenvalues point...
                    (/1.0_dp,0.0_dp,0.0_dp/),(/0.0_dp,1.0_dp,0.0_dp/),.true.)!... to the x and y axes

    do i = 1, this%at%N
       this%at%pos(:,i) = Rotate(this%at%pos(:,i),q)
    end do

    !Find two non-parallel reference directions within the structure
    this%d1a = 1; this%d1b = 2
    this%d2a = 1; this%d2b = 3
    do
       d1 = this%at%pos(:,this%d1b) - this%at%pos(:,this%d1a)
       d2 = this%at%pos(:,this%d2b) - this%at%pos(:,this%d2a)
       d1 = d1 / norm(d1)
       d2 = d2 / norm(d2)
       if (abs(d1 .dot. d2) .fne. 1.0_dp) exit
       !Directions are parallel... try a different set.
       this%d2b = this%d2b + 1
       if (this%d2b > this%at%N) then
          this%d1b = this%d1b + 1
          this%d2b = this%d1b + 1
          if (this%d2b == this%at%N) call System_Abort('RigidBodyModel_Initialise: Structure is Linear!')
       end if
    end do

    !Now the object is initialised
    this%initialised = .true.

  end subroutine rigidbodymodel_initialise

  subroutine rigidbodymodel_finalise(this)

    type(RigidBodyModel), intent(inout) :: this

    call finalise(this%at)
    this%d1a = 0; this%d1b = 0
    this%d2a = 0; this%d2b = 0
    this%I = 0.0_dp
    this%Mtot = 0.0_dp

    this%initialised = .false.
    
  end subroutine rigidbodymodel_finalise

  !
  !% Initialise: The atoms in 'indices' must be in the same order as in the RigidModel
  !
  subroutine rigidbody_initialise(this,at,indices,RigidModel)

    type(RigidBody),              intent(inout) :: this       !% The rigid body object we are initialising
    type(Atoms),                  intent(inout) :: at         !% The atoms object to which the indices refer
    integer, dimension(:),        intent(in)    :: indices    !% The indices of the atoms which make up the rigid body
    type(RigidBodyModel), target, intent(in)    :: RigidModel !% A model of the rigid body in its proper conformation
    !local variables
    integer                              :: i

    !Check the model has been initialised
    if (.not.RigidModel%initialised) call system_abort('RigidBody_Initialise: Model is not initialised')

    !Check for the correct number of atoms
    if (size(indices) /= RigidModel%at%N) then
       write(line,'(a,i0,a,i0,a)')'RigidBody_Initialise: The model contains ',RigidModel%at%N,' atoms but ', &
                                  size(indices),' indices have been provided'
       call system_abort(line)
    end if

    !Check the atomic indices and atomic numbers
    do i = 1, size(indices)
       if (indices(i) > at%N) then
          write(line,'(a,i0,a,i0,a)')'RigidBody_Initialise: Atom ',indices(i),' is out of range (',at%N,')'
          call system_abort(line)
       end if
       if (RigidModel%at%Z(i) /= at%Z(indices(i))) then
          write(line,'(4(a,i0),a)')'RigidBody_Initialise: Atom ',i, &
                                     ' of model ('//trim(ElementName(RigidModel%at%Z(i)))//', Z=', &
                                     RigidModel%at%Z(i),') does not match atom ',indices(i), &
                                     ' of atoms object ('//trim(ElementName(at%Z(indices(i))))//', Z=', &
                                     at%Z(indices(i)),')'
          call system_abort(line)
       end if
    end do

    call reallocate(this%indices,size(indices))
    this%indices = indices
    this%N = size(indices)
    
    this%model => RigidModel

    this%q = 0.0_dp
    this%p = 0.0_dp
    this%RCoM = 0.0_dp
    this%VCoM = 0.0_dp

    !Now we are initialised!
    this%initialised = .true.

  end subroutine rigidbody_initialise

  subroutine rigidbody_finalise(this)

    type(RigidBody), intent(inout) :: this

    this%N = 0
    if (allocated(this%indices)) deallocate(this%indices)
    this%model => null()

    this%q = 0.0_dp
    this%p = 0.0_dp
    this%RCoM = 0.0_dp
    this%VCoM = 0.0_dp

    this%initialised = .false.

  end subroutine rigidbody_finalise

  !
  !% Finalise an allocatable array of RigidBodys
  !
  subroutine rigidbodies_finalise(this)

    type(RigidBody), dimension(:), allocatable, intent(inout) :: this
    integer                                                   :: i

    if (allocated(this)) then
       do i = 1, size(this)
          call finalise(this(i))
       end do
       deallocate(this)
    end if
    
  end subroutine rigidbodies_finalise

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X ASSIGNMENT
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine rigidbody_assignment(to,from)

    type(RigidBody), intent(inout) :: to
    type(RigidBody), intent(in)    :: from

    call finalise(to)

    if (from%initialised) then

       allocate(to%indices(size(from%indices)))

       to%N           =  from%N
       to%indices     =  from%indices
       to%model       => from%model
       to%q           =  from%q
       to%p           =  from%p
       to%RCoM        =  from%RCoM
       to%VCoM        =  from%VCoM
       to%initialised =  from%initialised
       
    end if

  end subroutine rigidbody_assignment

  subroutine rigidbodies_assignment(to,from)

    type(RigidBody), dimension(:), intent(inout) :: to
    type(RigidBody), dimension(:), intent(in)    :: from
    integer                                      :: i
    
    if (size(to) /= size(from)) call system_abort('RigidBodies_Assignment: Target array has different length to source array')

    do i = 1, size(from)
       to(i) = from(i)
    end do

  end subroutine rigidbodies_assignment

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X I/O
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine rigidbodymodel_print(this, file)

    type(RigidBodyModel), intent(inout) :: this
    type(Inoutput), intent(inout), optional :: file

    call print('================', file=file)
    call print(' RigidBodyModel', file=file)
    call print('================', file=file)
    call print('', file=file)
    if (this%initialised) then
       call print('Atoms:', file=file)
       call print('------------', file=file)
       call print(this%at, file=file)
       call print('------------', file=file)
       call print('', file=file)
       write(line,'(2(a,i0))') 'Reference direction 1: atom ',this%d1a,' --> atom ',this%d1b
       call print(line, file=file)
       write(line,'(2(a,i0))') 'Reference direction 2: atom ',this%d2a,' --> atom ',this%d2b
       call print(line, file=file)
       call print('', file=file)
       call print('Inertia Tensor:', file=file)
       call print(inertia_tensor(this%at), file=file)
       call print('Stored moments of inertia:')
       call print(this%I, file=file)
       call print('', file=file)
       call print('Total Mass:', file=file)
       call print(this%Mtot, file=file)
    else
       call print('(uninitialised)', file=file)
    end if
    call print('', file=file)
    call print('================', file=file)

  end subroutine rigidbodymodel_print

  subroutine rigidbody_print(this, file)

    type(RigidBody), intent(inout) :: this
    type(Inoutput), intent(inout), optional :: file

    call print('XXXXXXXXXXXXXXXX', file=file)
    call print('   Rigid Body', file=file)
    call print('XXXXXXXXXXXXXXXX', file=file)
    call print('', file=file)
    if (this%initialised) then
       write(line,'(a,i0)')'Number of atoms = ',this%N
       call print(line, file=file)
       call print('', file=file)
       call print('Atomic indices:', file=file)
       call print(this%indices, file=file)
       call print('', file=file)
       call print('Reference Model =>', file=file)
       call print(this%model, file=file)
       call print('', file=file)
       call print('Orientation quaternion:', file=file)
       call print(this%q, file=file)
       call print('Conjugate momentum quaternion:', file=file)
       call print(this%p, file=file)
       call print('', file=file)
       call print('Centre of mass position:', file=file)
       call print(this%RCoM, file=file)
       call print('Centre of mass velocity:', file=file)
       call print(this%VCoM, file=file)
    else
       call print('(uninitialised)', file=file)
    end if
    call print('', file=file)
    call print('XXXXXXXXXXXXXXXX', file=file)

  end subroutine rigidbody_print

  subroutine rigidbodies_print(this,file)

    type(Rigidbody), dimension(:), intent(inout) :: this
    type(inoutput), optional,      intent(inout) :: file
    integer :: i

    do i = 1, size(this)
       call print(this(i),file)
    end do

  end subroutine rigidbodies_print

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X FUNCTIONS FOR WORKING WITH RIGID BODIES
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function inertia_tensor(this) result(I)

    type(Atoms), intent(in)  :: this
    real(dp), dimension(3,3) :: I
    !local variables
    integer                  :: j
    real(dp)                 :: r(3), m, x2, xy, xz, y2, yz, z2

    I = 0.0_dp
    do j = 1, this%N
       r = this%pos(:,j)
       m = this%mass(j)
       x2 = r(1)*r(1); y2 = r(2)*r(2); z2 = r(3)*r(3)
       xy = r(1)*r(2); xz = r(1)*r(3); yz = r(2)*r(3)
       I(1,1) = I(1,1) + m * (y2 + z2)
       I(2,2) = I(2,2) + m * (x2 + z2)
       I(3,3) = I(3,3) + m * (x2 + y2)
       I(1,2) = I(1,2) - m * xy
       I(1,3) = I(1,3) - m * xz
       I(2,3) = I(2,3) - m * yz
    end do
    !Fill in symmetric bits
    I(2,1) = I(1,2)
    I(3,1) = I(1,3)
    I(3,2) = I(2,3)    
    
  end function inertia_tensor

  !
  ! Functions and routines required by NO_SQUISH
  !
  ! From Miller et al.:
  ! P_0 q = {q0,q1,q2,q3}   P_1 q = {-q1,q0,q3,-q2}   P_2 q = {-q2,-q3,q0,q1}   P_3 q = {-q3,q2,-q1,q0}
  !
  !               ( (  |  ) (  |  ) (  |  ) (  |  ) )
  ! Matrix S(q) = ( (P_0 q) (P_1 q) (P_2 q) (P_3 q) )
  !               ( (  |  ) (  |  ) (  |  ) (  |  ) )
  !

  function P0(q)  !does nothing, but here for completeness
    type(Quaternion), intent(in) :: q
    type(Quaternion)             :: P0
    P0 = q
  end function P0

  function P1(q)
    type(Quaternion), intent(in) :: q
    type(Quaternion)             :: P1
    P1 = (/-q%b,q%a,q%d,-q%c/)
  end function P1

  function P2(q)
    type(Quaternion), intent(in) :: q
    type(Quaternion)             :: P2
    P2 = (/-q%c,-q%d,q%a,q%b/)
  end function P2

  function P3(q)
    type(Quaternion), intent(in) :: q
    type(Quaternion)             :: P3
    P3 = (/-q%d,q%c,-q%b,q%a/)
  end function P3

  function no_squish_S(q) result(S)
    type(Quaternion), intent(in) :: q
    real(dp), dimension(4,4)     :: S
    S(1,1) = +q%a; S(1,2) = -q%b; S(1,3) = -q%c; S(1,4) = -q%d
    S(2,1) = +q%b; S(2,2) = +q%a; S(2,3) = -q%d; S(2,4) = +q%c
    S(3,1) = +q%c; S(3,2) = +q%d; S(3,3) = +q%a; S(3,4) = -q%b
    S(4,1) = +q%d; S(4,2) = -q%c; S(4,3) = +q%b; S(4,4) = +q%a
  end function no_squish_S

  function no_squish_A_dot_transpose(q,q_dot) result(A_t)

    type(Quaternion), intent(in) :: q,q_dot
    real(dp), dimension(3,3)     :: A_t
   
    A_t(1,1) = q%a*q_dot%a + q%b*q_dot%b - q%c*q_dot%c - q%d*q_dot%d
    A_t(2,1) = q%b*q_dot%c + q%c*q_dot%b + q%a*q_dot%d + q%d*q_dot%a
    A_t(3,1) = q%b*q_dot%d + q%d*q_dot%b - q%a*q_dot%c - q%c*q_dot%a
    A_t(1,2) = q%b*q_dot%c + q%c*q_dot%b - q%a*q_dot%d - q%d*q_dot%a
    A_t(2,2) = q%a*q_dot%a - q%b*q_dot%b + q%c*q_dot%c - q%d*q_dot%d
    A_t(3,2) = q%c*q_dot%d + q%d*q_dot%c + q%a*q_dot%b + q%b*q_dot%a
    A_t(1,3) = q%b*q_dot%d + q%d*q_dot%b + q%a*q_dot%c + q%c*q_dot%a
    A_t(2,3) = q%c*q_dot%d + q%d*q_dot%c - q%a*q_dot%b - q%b*q_dot%a
    A_t(3,3) = q%a*q_dot%a - q%b*q_dot%b - q%c*q_dot%c + q%d*q_dot%d

    A_t = 2.0_dp * A_t
    
  end function no_squish_A_dot_transpose

  !
  ! Free_Rotor: part of the NO_SQUISH algorithm
  !
  subroutine no_squish_free_rotor(q,p,I,dt)

    type(Quaternion),       intent(inout) :: q,p
    real(dp), dimension(3), intent(in)    :: I
    real(dp)                              :: dt
    !local variables
    real(dp)                              :: zetadt_2, zetadt, delta_t
    integer                               :: j

    !Free rotation in M_ROT steps
    delta_t = dt / real(M_ROT,dp)
    do j = 1, M_ROT

       zetadt_2 = delta_t * (p .dot. P3(q)) / (8.0_dp * I(3))
       q = cos(zetadt_2)*q + sin(zetadt_2)*P3(q)
       p = cos(zetadt_2)*p + sin(zetadt_2)*P3(p)

       zetadt_2 = delta_t * (p .dot. P2(q)) / (8.0_dp * I(2))
       q = cos(zetadt_2)*q + sin(zetadt_2)*P2(q)
       p = cos(zetadt_2)*p + sin(zetadt_2)*P2(p)

       zetadt   = delta_t * (p .dot. P1(q)) / (4.0_dp * I(1))
       q = cos(zetadt_2)*q + sin(zetadt_2)*P1(q)
       p = cos(zetadt_2)*p + sin(zetadt_2)*P1(p)

       zetadt_2 = delta_t * (p .dot. P2(q)) / (8.0_dp * I(2))
       q = cos(zetadt_2)*q + sin(zetadt_2)*P2(q)
       p = cos(zetadt_2)*p + sin(zetadt_2)*P2(p)

       zetadt_2 = delta_t * (p .dot. P3(q)) / (8.0_dp * I(3))
       q = cos(zetadt_2)*q + sin(zetadt_2)*P3(q)
       p = cos(zetadt_2)*p + sin(zetadt_2)*P3(p)

    end do   

  end subroutine no_squish_free_rotor

end module RigidBody_module

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !X
   !X RIGID BODY ROUTINES
   !X
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
   !This is based on the NO_SQUISH algorithm
   !Miller et al., JChemPhys 116(20) 8649 (2002)
   
   !
   !% From Positions and Velocities of individual atoms, calculate the
   !% centre of mass, centre of mass velocity,  the orientation, and the
   !% orientation conjugate momentum.
   !
!    subroutine rigidbody_calculate_variables(this,i)
!      
!      type(DynamicalSystem), intent(inout) :: this
!      integer,               intent(in)    :: i
!      integer                              :: j,k, d1s_a, d1s_b, d2s_a, d2s_b
!      real(dp), dimension(3)               :: VCoM,L,r,v, d1s,d2s,d1b,d2b
!      real(dp), dimension(4,4)             :: S
!      real(dp), dimension(4)               :: L4
!      real(dp)                             :: m
!      type(Quaternion)                     :: q
!      
!      if (i < 1 .or. i > this%Nrigid) then
!         write(line,'(a,i0,a)')'RigidBody_CalculateVariables: Rigid Body ',i,' is out of range'
!         call system_abort(line)
!      end if
!      
!      !Work out which atoms the reference directions go bewteen in the full atomic system
!      d1s_a = this%rigidbody(i)%indices(this%rigidbody(i)%model%d1a)
!      d1s_b = this%rigidbody(i)%indices(this%rigidbody(i)%model%d1b)
!      d2s_a = this%rigidbody(i)%indices(this%rigidbody(i)%model%d2a)
!      d2s_b = this%rigidbody(i)%indices(this%rigidbody(i)%model%d2b)
!      !Find the directions in body and space coordinates
!      d1s = diff_min_image(this%atoms,d1s_a,d1s_b)
!      d2s = diff_min_image(this%atoms,d2s_a,d2s_b)
!      d1b = this%rigidbody(i)%model%at%pos(:,this%rigidbody(i)%model%d1b) - &
!           this%rigidbody(i)%model%at%pos(:,this%rigidbody(i)%model%d1a)
!      d2b = this%rigidbody(i)%model%at%pos(:,this%rigidbody(i)%model%d2b) - &
!           this%rigidbody(i)%model%at%pos(:,this%rigidbody(i)%model%d2a)
!      !Calculate the rotation quaternion to go from body to space coords
!      q = Orientation(d1b,d2b,d1s,d2s,correct_angle = .true.)
!      this%rigidbody(i)%q = q
!      !Calculate its angular momentum in the body frame of reference
!      VCoM = centre_of_mass_velo(this,this%rigidbody(i)%indices)
!      L = 0.0_dp
!      q = this%rigidbody(i)%q
!      do j = 1, this%rigidbody(i)%N
!         
!         k = this%rigidbody(i)%indices(j) !index of the jth atom in rigid body i
!         m = this%atoms%Z(k) !the mass of the jth atom
!         r = this%rigidbody(i)%model%at%pos(:,j) !the body fixed coordinates of jth atom 
!         v = Rotate((this%atoms%velo(:,k) - VCoM),.conj.q) !body fixed velocity of jth atom
!         L = L + m * (r .cross. v)
!         
!      end do
!      
!      !Calculate the momentum conjugate to the orientation quaternion
!      L4 = (/0.0_dp,L/) !this is D^-1 omega(4) in equation 2.15 of Miller et al.
!      S = no_squish_S(q)
!      this%rigidbody(i)%p = 2.0_dp * (S .mult. L4)     
!      
!      !Store the centre of mass position and velocity
!      this%rigidbody(i)%RCoM = centre_of_mass(this%atoms,this%rigidbody(i)%indices)
!      this%rigidbody(i)%VCoM = VCoM
!      
!    end subroutine rigidbody_calculate_variables
!    
!    subroutine rigidbody_write_positions(this,i)
!      
!      type(DynamicalSystem), intent(inout) :: this
!      integer,               intent(in)    :: i
!      !local variables
!      integer                              :: j, k
!      real(dp), dimension(3)               :: r
!      type(Quaternion)                     :: q
!      
!      if (i < 1 .or. i > this%Nrigid) then
!         write(line,'(a,i0,a)')'RigidBody_WritePositions: Rigid Body ',i,' is out of range'
!         call system_abort(line)
!      end if
!      
!      q = this%rigidbody(i)%q
!      
!      do j = 1, this%rigidbody(i)%model%at%N
!         k = this%rigidbody(i)%indices(j)        !jth atom of rigid body
!         r = this%rigidbody(i)%model%at%pos(:,j) !body coords of jth atom
!         this%atoms%pos(:,k) = this%rigidbody(i)%RCoM + Rotate(r,this%rigidbody(i)%q)
!      end do
!      
!    end subroutine rigidbody_write_positions
!    
!    subroutine rigidbody_write_velocities(this,i)
!      
!      type(DynamicalSystem), intent(inout) :: this
!      integer,               intent(in)    :: i
!      !local variables
!      integer                              :: j, k
!      real(dp), dimension(3)               :: r, omega
!      type(Quaternion)                     :: q
!      
!      if (i < 1 .or. i > this%Nrigid) then
!         write(line,'(a,i0,a)')'RigidBody_WriteVelocities: Rigid Body ',i,' is out of range'
!         call system_abort(line)
!      end if
!      
!      q = this%rigidbody(i)%q
!      omega = RigidBody_angular_velocity(this,i)
!      
!      do j = 1, this%rigidbody(i)%model%at%N
!         k = this%rigidbody(i)%indices(j)        !jth atom of rigid body
!         r = this%rigidbody(i)%model%at%pos(:,j) !body coords of jth atom
!         this%atoms%velo(:,k) = this%rigidbody(i)%VCoM + Rotate((omega .cross. r),q)
!      end do
!      
!    end subroutine rigidbody_write_velocities
!
!    !
!    !% Reconstitute a rigid body that has been split over periodic boundaries
!    !
!    subroutine rigidbody_reconstitute(this,i)
!      
!      type(DynamicalSystem), intent(inout) :: this
!      integer,               intent(in)    :: i
!      integer                              :: j, k, n, t(3)
!      real(dp), dimension(3)               :: r
!      
!      if (i < 1 .or. i > this%Nrigid) then
!         write(line,'(a,i0,a)')'RigidBody_Reconstitute: Rigid Body ',i,' is out of range'
!         call system_abort(line)
!      end if
!      
!      !Map everything back into the cell of the first atom in the body
!      j = this%rigidbody(i)%indices(1)
!      r = this%atoms%pos(:,j)
!      t = this%atoms%travel(:,j)
!      do n = 2, this%rigidbody(i)%N
!         k = this%rigidbody(i)%indices(n)
!         !Alter the position
!         this%atoms%pos(:,k) = r + diff_min_image(this%atoms,j,k)
!         !Make the travel values agree... RIGID BODIES CANNOT SPAN MORE THAN ONE CELL (why would they?)
!         this%atoms%travel(:,k) = t
!      end do
!      
!    end subroutine rigidbody_reconstitute
!    
!    !
!    !% Calculate the current angular velocity of rigid body 'i', in body frame
!    !
!    function rigidbody_angular_velocity(this,i) result(omega)
!      
!      type(DynamicalSystem), intent(in) :: this
!      integer,               intent(in) :: i
!      real(dp), dimension(3)            :: omega
!      real(dp), dimension(4)            :: p4, sTp
!      
!      if (i < 1 .or. i > this%Nrigid) then
!         write(line,'(a,i0,a)')'RigidBody_AngularVelocity: Rigid Body ',i,' is out of range'
!         call System_Abort(line)
!      end if
!      
!      p4 = this%rigidbody(i)%p
!      !sTp(1) should always be zero
!      sTp = transpose(no_squish_S(this%rigidbody(i)%q)) .mult. p4
!      omega = 0.5_dp * sTp(2:4) / this%rigidbody(i)%model%I
!      
!    end function rigidbody_angular_velocity
!    
!    !
!    !% Calculate the current angular velocity of rigid body 'i', in body frame
!    !
!    function rigidbody_angular_momentum(this,i) result(L)
!      
!      type(DynamicalSystem), intent(in) :: this
!      integer,               intent(in) :: i
!      real(dp), dimension(3)            :: L
!      
!      if (i < 1 .or. i > this%Nrigid) then
!         write(line,'(a,i0,a)')'RigidBody_AngularMomentum: Rigid Body ',i,' is out of range'
!         call system_abort(line)
!      end if
!      
!      L = this%rigidbody(i)%model%I * RigidBody_Angular_Velocity(this,i) !elementwise multiplication since 
!                                                                         !I is a vector, not a matrix
!      
!    end function rigidbody_angular_momentum
!    
!    !
!    !% Calculate the current torque on the rigid body 'i', in body frame
!    !
!    function rigidbody_torque(this,i) result(tau)
!      
!      type(DynamicalSystem), intent(in) :: this
!      integer,               intent(in) :: i
!      real(dp), dimension(3)            :: tau, ACoM, r, a
!      real                              :: m
!      integer                           :: j, k
!      type(Quaternion)                  :: q
!      
!      if (i < 1 .or. i > this%Nrigid) then
!         write(line,'(a,i0,a)')'RigidBody_Torque: Rigid Body ',i,' is out of range'
!         call system_abort(line)
!      end if
!      
!      ACoM = centre_of_mass_acc(this,this%rigidbody(i)%indices)
!      q = this%rigidbody(i)%q
!      tau = 0.0_dp
!      
!      do j = 1, this%rigidbody(i)%model%at%N
!         k = this%rigidbody(i)%indices(j)        !jth atom of rigid body i
!         m = this%rigidbody(i)%model%at%mass(j)
!         r = this%rigidbody(i)%model%at%pos(:,j) !body coords of jth atom
!         a = Rotate((this%atoms%acc(:,k) - ACoM), .conj.q) !Rotate the non-CoM acceleration into the body frame      
!         tau = tau + m * (r .cross. a)
!      end do
!      
!    end function rigidbody_torque
!    
!    function rigidbody_rotational_energy(this,i) result(Erot)
!      
!      type(DynamicalSystem), intent(in) :: this
!      integer,               intent(in) :: i
!      real(dp)                          :: Erot
!      real(dp), dimension(3)            :: Im, omega
!      
!      
!      if (i < 1 .or. i > this%Nrigid) then
!         write(line,'(a,i0,a)')'RigidBody_RotationalEnergy: Rigid Body ',i,' is out of range'
!         call system_abort(line)
!      end if
!      
!      omega = RigidBody_Angular_Velocity(this,i)
!      Im = this%rigidbody(i)%model%I    
!      Erot = 0.5_dp * (Im(1)*omega(1)*omega(1) + Im(2)*omega(2)*omega(2) + Im(3)*omega(3)*omega(3))
!      
!    end function rigidbody_rotational_energy
!    
!    function rigidbody_translational_energy(this,i) result(Etrans)
!      
!      type(DynamicalSystem), intent(in) :: this
!      integer,               intent(in) :: i
!      real(dp)                          :: Etrans
!      
!      if (i < 1 .or. i > this%Nrigid) then
!         write(line,'(a,i0,a)')'RigidBody_TranslationalEnergy: Rigid Body ',i,' is out of range'
!         call system_abort(line)
!      end if
!      
!      Etrans = 0.5_dp*this%rigidbody(i)%model%Mtot * normsq(centre_of_mass_velo(this,this%rigidbody(i)%indices))
!      
!    end function rigidbody_translational_energy
!
!   !
!   !% Add rigid bodies to the DynamicalSystem
!   !
!   subroutine ds_add_rigidbody(this,indices,model)
!     
!     type(DynamicalSystem),        intent(inout) :: this
!     integer, dimension(:),        intent(in)    :: indices
!     type(RigidBodyModel), target, intent(in)    :: model
!     !local variables
!     integer                                     :: i, g1, g2
!
!     !Check if rigid bodies have been set up
!     if (.not.allocated(this%rigidbody)) call system_abort('DS_AddRigidBody: Rigid Bodies have not been allocated')
!
!     !Have we got space for another rigid body?
!     if (this%Nrigid == size(this%rigidbody)) then
!        write(line,'(a,i0,a)')'DS_AddRigidBody: Maximum number of rigid bodies reached (',this%Nrigid,')'
!        call system_abort(line)
!     end if
!
!     !Initialise the rigid body and update the group variables
!     this%Nrigid = this%Nrigid + 1
!     g1 = this%group_lookup(indices(1))
!     do i = 1, size(indices)
!        if (Atom_Type(this,indices(i)) /= TYPE_ATOM) then
!           write(line,'(a,i0,a)')'DS_AddRigidBody: Atom ',indices(i),&
!                ' is already being treated by another integration method'
!           call system_abort(line)
!        end if
!        g2 = this%group_lookup(indices(i))
!        call merge_groups(this%group(g1),this%group(g2))
!     end do
!
!     !Make atoms belong to this rigid body
!     call group_add_object(this%group(g1),this%Nrigid)
!     call initialise(this%rigidbody(this%Nrigid),this%atoms,indices,model)
!
!     !Calculate the dynamical variables
!     call rigidbody_calculate_variables(this,this%Nrigid)
!
!     !Reduce the number of degrees of freedom:
!     !Subtract 3 translational DoF for each atom, then add on the 3 trans + 3 rot DoF for the body
!     this%Ndof = this%Ndof - 3*size(indices) + 6
!     
!   end subroutine ds_add_rigidbody
!
!   !
!   !% Correct the positions of all the atoms in rigid bodies
!   !
!   subroutine ds_correct_rigid_positions(this)
!
!     type(DynamicalSystem), intent(inout) :: this
!     integer                              :: i,j,n
!
!     do i = 1, this%Nrigid
!        
!        call rigidbody_calculate_variables(this,i)
!
!        do j = 1, this%rigidbody(i)%N
!           n = this%rigidbody(i)%indices(j)
!           this%atoms%pos(:,n) = Rotate(this%rigidbody(i)%model%at%pos(:,j), this%rigidbody(i)%q) &
!                                 + this%rigidbody(i)%RCoM
!        end do
!        
!     end do
!
!   end subroutine ds_correct_rigid_positions
