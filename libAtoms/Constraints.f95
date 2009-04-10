!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     libAtoms: atomistic simulation library
!X     
!X     Copyright 2006-2007.
!X
!X     Authors: Gabor Csanyi, Steven Winfield, James Kermode
!X     Contributors: Noam Bernstein, Alessio Comisso
!X
!X     The source code is released under the GNU General Public License,
!X     version 2, http://www.gnu.org/copyleft/gpl.html
!X
!X     If you would like to license the source code under different terms,
!X     please contact Gabor Csanyi, gabor@csanyi.net
!X
!X     When using this software, please cite the following reference:
!X
!X     http://www.libatoms.org
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X  Constraints module
!X  
!%  Contains the Constraint type and constraint subroutines
!%  
!%  Constraints are applied using the RATTLE algorithm:
!%  H. C. Andersen, JCompPhys 52 24-34 (1983)
!%  This itself is based on the SHAKE algorithm:
!%  J. P. Ryckaert et al., JCompPhys 23 327-341 (1977)
!%
!%  To integrate the system with constraints Newtons equations of motion have to be modified to
!%  \begin{displaymath}
!%  m\mathbf{a} = \mathbf{F} - \lambda\nabla\sigma
!%  \end{displaymath}
!%  where $\sigma$ is a constraint function (see below) and
!%  $\lambda$ is a Lagrange multiplier which is dependent upon atomic positions
!%  and velocities.
!%  
!%  When the equations of motion are discretised two separate approximation to $\lambda$ must be
!%  made so that the positions and velocities both obey the constraint exactly (or to within a
!%  specified precision) at each time step.
!%  
!%  During the position updates, a Lagrange multiplier $\lambda_j$ attached to a constraint
!%  $\sigma_j$, which involves particles $i=1 \ldots N$ is solved for iteratively by the following:
!%  \begin{displaymath}
!%  \lambda_j^S = \frac{2\sigma_j}{(\Delta t)^2 \sum_i \frac{1}{m_i} \nabla_i \sigma_j(t+\Delta t) \cdot 
!%                                                                       \nabla_i \sigma_j(t)}
!%  \end{displaymath}
!%  
!%  The velocity updates also include the $\lambda_j^S\nabla\sigma_j$ term, plus an additional
!%  $\lambda_j^R\nabla\sigma_j$ term where $\lambda_j^R$ is solved for via:
!%  \begin{displaymath}
!%  \lambda_j^R = \frac{2\dot{\sigma}_j}{\Delta t \sum_i \frac{1}{m_i} |\nabla_i \sigma_j(t+\Delta t)|^2}
!%  \end{displaymath}
!%
!%  The full, non-specialised form of the iterative solver is present in the code for added flexibility and
!%  to allow users to write and use their own constraint functions, which can be time dependent.
!%
!%
!% \textbf{Writing extra constraint subroutines}
!% 
!% A constraint subroutine evaluates a constraint function ('C'), its derivatives with respect to all
!% coordinates ('dC_dr') and its full time derivative ('dC_dt'). A constraint function is a function which 
!% evaluates to zero when the constraint is satisfied.
!% 
!% Probably the easiest way of constructing a new subroutine is to copy the following and fill in the gaps:
!% 
!%>   subroutine CONSTRAINT(pos, velo, t, data, C, dC_dr, dC_dt)
!%> 
!%>     real(dp), dimension(:),         intent(in)  :: pos, velo, data
!%>     real(dp),                       intent(in)  :: t
!%>     real(dp),                       intent(out) :: C
!%>     real(dp), dimension(size(pos)), intent(out) :: dC_dr
!%>     real(dp),                       intent(out) :: dC_dt
!%>     if(size(pos) /= ?) &
!%>         call System_Abort('CONSTRAINT: Exactly ? atoms must be specified')
!%>     if(size(data) /= ?) &
!%>         call System_Abort('CONSTRAINT: "data" must contain exactly ? value(s)')
!%>     C = ?
!%>     dC_dr = ?
!%>     dC_dt = ?
!%> 
!%>   end subroutine CONSTRAINT
!% 
!% A few notes:
!% \begin{itemize}
!%  \item The subroutine must make its own checks on the sizes of the input arguments (as shown)
!% 
!%  \item dC_dr should not be zero when the constraint is satisfied, since then no force can be applied to 
!%    enforce the constraint. For example, the plane constraint could have been defined as the squared 
!%    distance from the plane, but that would cause it to suffer from this problem, so the signed 
!%    distance was used instead.
!% 
!%  \item Currently, a constraint function can only depend explicitly on positions and time, so the kinetic
!%    energy of a group of particles cannot be constrained for example.
!% 
!%  \item dC_dt will usually be dC_dr .dot. velo, except in the case where the constraint depends explicitly 
!%    on time; then there will be an extra partial d/dt term.
!%
!%  \item If a constraint contains one particle only, then the absolute position of that particle is passed
!%    in 'pos'. If more than one particle is in the constraint then the relative positions are passed
!%    (with the first particle at the origin), so periodic boundary conditions should not be a worry.
!%
!%  \item When a constraint is added to a DynamicalSystem its number of degrees of freedom is decreased by 1.
!%    A constraint should only attempt to remove one degree of freedom MAXIMUM. Constraints which try to,
!%    for examples, keep a particle on a line (removing 2 degrees of freedom) are doomed to failure since
!%    the particle will not always be able to return to the line by the application of forces in the direction
!%    from the particle to the nearest point on the line. Instead, implement the line as the intersection of
!%    two planes, which makes the removal of two degrees of freedom more explicit, and allows the algorithm
!%    to use two directions to construct the constraint force.
!%
!%  \item Some constraints will be extraneous (e.g. constraining a bond length which, by the action of many
!%    other constraints, is already fixed implicitly), and the decrease of ds%Ndof will be wrong, in which 
!%    case it should be correctly set by hand. Also, the algorithm will fail if the 'explicit' constraint 
!%    is not the same as the 'implicit' one.
!%  \end{itemize}
!% 
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! $Id: Constraints.f95,v 1.13 2008-02-04 13:31:19 saw44 Exp $

! $Log: Constraints.f95,v $
! Revision 1.13  2008/02/04 13:31:19  saw44
! Fixed comment typo
!
! Revision 1.12  2007/10/08 15:13:09  saw44
! Split lambda into lambdaR and lambdaV which are not zeroed between calls
!
! Revision 1.11  2007/04/17 15:42:39  jrk33
! Standardised subroutine and function references and printing argument order.
!
! Revision 1.10  2007/04/17 09:57:19  gc121
! put copyright statement in each file
!
! Revision 1.9  2007/04/11 15:44:17  saw44
! Updated/Added comments for the documentation generator
!
! Revision 1.8  2007/03/12 16:52:52  jrk33
! Reformatted documentation
!
! Revision 1.7  2007/03/01 13:51:45  jrk33
! Documentation comments reformatted and edited throughout. Anything starting "!(no space)%" 
! is picked up by the documentation generation script
!
! Revision 1.6  2007/02/28 15:41:06  saw44
! Converted SMOOTH_BOND to CUBIC_BOND
!
! Revision 1.5  2007/01/24 11:20:49  saw44
! Removed unused variables
!
! Revision 1.4  2007/01/19 15:18:44  saw44
! Reduced WARNING_TOLERANCE to allow for numerical errors in lapack matrix inversion
!
! Revision 1.3  2007/01/09 16:56:49  saw44
! Added CONSTRAINT_WARNING_TOLERANCE parameter (used in DynamicalSystem)
!
! Revision 1.2  2006/12/13 14:06:28  saw44
! Made data and func optional in Constraint_Amend. Increased MAX_ITERATIONS and reduced convergence tolerance
!
! Revision 1.1.1.1  2006/12/04 11:11:30  gc121
! Imported sources
!
! Revision 1.5  2006/11/27 14:30:56  saw44
! Added SMOOTH_BOND and RELAX_BOND constraint subroutines. Added Constraint_Amend routine. Made the number of registered constraint subroutines a global to allow checking in Constraint_Amend
!
! Revision 1.4  2006/11/24 12:19:26  saw44
! Updated some comments
!
! Revision 1.3  2006/11/21 13:37:56  saw44
! Added time dependence to constraints. Printing a constraint now gives more info
!
! Revision 1.2  2006/11/17 13:05:33  saw44
! Added assignment subroutines. Fixed ReadB interface clash.
!
! Revision 1.1  2006/11/13 15:23:49  saw44
! Constraints added to the repository
!
module constraints_module

  use system_module
  use linearalgebra_module 
  use atoms_module
  use group_module

  implicit none

  integer, parameter  :: MAX_CONSTRAINT_SUBS = 20         !% This must be consistent with Constraint_Pointers.c
  integer             :: REGISTERED_CONSTRAINTS = 0       !% Stores the number of registered constraint functions

  integer, parameter  :: RATTLE_MAX_ITERATIONS = 3000      !% Max iterations before failure (same as AMBER)
  real(dp), parameter :: CONSTRAINT_TOLERANCE = 1.0E-10_dp !% Default tolerance of constraints
  real(dp), parameter :: CONSTRAINT_WARNING_TOLERANCE = 1.0E-3_dp !% Warn if a constraint is added which is not obeyed
                                                                  !% to within this tolerance
  type Constraint

     integer                             :: N = 0     !% Number of atoms in the constraint
     integer,  allocatable, dimension(:) :: atom      !% Indices of the atoms in the constraint
     real(dp), allocatable, dimension(:) :: data      !% Data to be passed to the constraint function
     integer                             :: func      !% Pointer to the constraint function
     real(dp)                            :: lambdaR   !% Lagrange multipler for position constraint forces
     real(dp)                            :: dlambdaR  !% Last computed correction to lambdaR
     real(dp)                            :: lambdaV   !% Lagrange multipler for velocity constraint forces
     real(dp)                            :: dlambdaV  !% Last computed correction to lambdaV
     real(dp)                            :: C         !% Value of the constraint
     real(dp), allocatable, dimension(:) :: dC_dr     !% Gradient of the constraint at the start of the timestep
     real(dp)                            :: dC_dt     !% Time derivative of the constraint
     real(dp), allocatable, dimension(:) :: old_dC_dr !% Value of 'dC_dr' at the beginning of the algorithm
     real(dp)                            :: tol       !% Convergence tolerance of the constraint
     logical                             :: initialised  = .false.

  end type Constraint

  interface initialise
     module procedure constraint_initialise
  end interface initialise

  interface finalise
     module procedure constraint_finalise, constraints_finalise
  end interface finalise

  interface assignment(=)
     module procedure constraint_assignment, constraints_assignment
  end interface assignment(=)
     
  interface print
     module procedure constraint_print, constraints_print
  end interface print

  interface write_binary
     module procedure write_binary_constraint, write_binary_constraints
  end interface write_binary

  interface read_binary
     module procedure read_binary_constraint, read_binary_constraints
  end interface read_binary

  !
  ! interface blocks for the C functions
  !
  !% OMIT
  interface 
     subroutine register_constraint_sub(sub)
       interface 
          subroutine sub(pos,velo,t,data,C,dC_dr,dC_dt)
            use system_module  !for dp definition
            real(dp), dimension(:),         intent(in)  :: pos,velo,data
            real(dp),                       intent(in)  :: t
            real(dp),                       intent(out) :: C, dC_dt
            real(dp), dimension(size(pos)), intent(out) :: dC_dr
          end subroutine sub
       end interface
     end subroutine register_constraint_sub
  end interface 

  !% OMIT
  interface 
     subroutine call_constraint_sub(i,pos,velo,t,data,C,dC_dr,dC_dt)
       use system_module
       integer                                     :: i
       real(dp), dimension(:),         intent(in)  :: pos, velo, data
       real(dp),                       intent(in)  :: t
       real(dp),                       intent(out) :: C, dC_dt
       real(dp), dimension(size(pos)), intent(out) :: dC_dr
     end subroutine call_constraint_sub
  end interface 

contains
  
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X INITIALISE
  !% 
  !% Initialise this Constraint.
  !%
  !% NOTE: the atomic indices are not sorted into ascending order, because
  !%       a three body constraint (for example) may depend on the order
  !%       (and the wrong angle could be constrained).
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  
  subroutine constraint_initialise(this,indices,func,data)

    type(Constraint),       intent(inout) :: this
    integer,  dimension(:), intent(in)    :: indices
    integer,                intent(in)    :: func
    real(dp), dimension(:), intent(in)    :: data

    if (this%initialised) call constraint_finalise(this)
    
    if (size(indices) == 0) &
         call system_abort('Constraint_Initialise: There must be at least one particle in the constraint')

    this%N = size(indices)
    allocate(this%atom(this%N))
    this%atom = indices

    allocate(this%data(size(data)))
    this%data = data

    if (func < 0 .or. func >= REGISTERED_CONSTRAINTS) then
       write(line,'(a,i0,a)')'Constraint_Initialised: Invalid constraint subroutine (',func,')'
       call System_Abort(line)
    end if

    this%func = func
    this%lambdaR = 0.0_dp
    this%lambdaV = 0.0_dp
    this%C = 0.0_dp
    this%dC_dt = 0.0_dp
    this%tol = CONSTRAINT_TOLERANCE

    allocate(this%dC_dr(3*size(indices)))
    this%dC_dr = 0.0_dp
    allocate(this%old_dC_dr(3*size(indices)))
    this%old_dC_dr = 0.0_dp

    !initialised
    this%initialised = .true.

  end subroutine constraint_initialise

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X FINALISE
  !X
  !% Finalise this Constraint object.
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine constraint_finalise(this)

    type(Constraint), intent(inout) :: this
    
    this%N = 0
    this%func = 0
    this%lambdaR = 0.0_dp
    this%lambdaV = 0.0_dp
    this%C = 0.0_dp
    this%dC_dt = 0.0_dp
    this%tol = 0.0_dp
    if (allocated(this%atom))  deallocate(this%atom)
    if (allocated(this%data))  deallocate(this%data)
    if (allocated(this%dC_dr)) deallocate(this%dC_dr)
    if (allocated(this%old_dC_dr)) deallocate(this%old_dC_dr)
    this%initialised = .false.

  end subroutine constraint_finalise

  subroutine constraints_finalise(this)

    type(Constraint), dimension(:), allocatable, intent(inout) :: this
    integer                                                    :: i

    if (allocated(this)) then
       do i = 1, size(this)
          call finalise(this(i))
       end do
       deallocate(this)
    end if

  end subroutine constraints_finalise

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X AMEND A CONSTRAINT: 
  !% Change constraint data and/or constraint function
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine constraint_amend(this,func,data)
    
    type(Constraint),                 intent(inout) :: this
    integer,                optional, intent(in)    :: func
    real(dp), dimension(:), optional, intent(in)    :: data

    if (present(data)) then
       call reallocate(this%data,size(data))
       this%data = data
    end if

    if (present(func)) then
       if (func < 0 .or. func >= REGISTERED_CONSTRAINTS) then
          write(line,'(a,i0,a)')'Constraint_Amend: Invalid constraint subroutine (',func,')'
          call System_Abort(line)
       end if
       
       this%func = func
    end if

  end subroutine constraint_amend

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X ASSIGNMENT
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine constraint_assignment(to,from)

    type(Constraint), intent(inout) :: to
    type(Constraint), intent(in)    :: from

    !First wipe the target of the assignment
    call Finalise(to)

    if (from%initialised) then
       
       !Allocate target arrays
       allocate(to%atom(size(from%atom)),   &
                to%data(size(from%data)),   &
                to%dC_dr(size(from%dC_dr)), &
                to%old_dC_dr(size(from%old_dC_dr)))
    
       !Then copy over the members one by one
       to%N           = from%N
       to%atom        = from%atom
       to%data        = from%data
       to%func        = from%func
       to%lambdaR     = from%lambdaR
       to%lambdaV     = from%lambdaV
       to%C           = from%C
       to%dC_dr       = from%dC_dr
       to%dC_dt       = from%dC_dt
       to%old_dC_dr   = from%old_dC_dr
       to%tol         = from%tol
       to%initialised = from%initialised

    end if

  end subroutine constraint_assignment

  subroutine constraints_assignment(to,from)

    type(Constraint), dimension(:), intent(inout) :: to
    type(Constraint), dimension(:), intent(in)    :: from
    integer                                       :: i

    do i = 1, size(from)
       to(i) = from(i)
    end do

  end subroutine constraints_assignment

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X CALCULATE CONSTRAINT VALUES
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine constraint_calculate_values(this,pos,velo,t)

    type(Constraint),       intent(inout) :: this
    real(dp), dimension(:), intent(in)    :: pos, velo
    real(dp),               intent(in)    :: t

    call call_constraint_sub(this%func,pos,velo,t,this%data,this%C,this%dC_dr,this%dC_dt)

  end subroutine constraint_calculate_values

  subroutine constraint_calculate_values_at(this,at,t)

    type(constraint), intent(inout) :: this
    type(atoms),      intent(in)    :: at
    real(dp),         intent(in)    :: t

    real(dp) :: pos(3*this%N), velo(3*this%N)
    integer  :: n, i, o

    o = this%atom(1)
    velo(1:3) = at%velo(:,o)
    if (this%N > 1) then
       pos(1:3) = 0.0_dp
       do n = 2, this%N
          i = this%atom(n)
          pos(3*n-2:3*n) = diff_min_image(at,o,i)
          velo(3*n-2:3*n) = at%velo(:,i)
       end do
    else
       pos = at%pos(:,o)
    end if

    call call_constraint_sub(this%func,pos,velo,t,this%data,this%C,this%dC_dr,this%dC_dt)

  end subroutine constraint_calculate_values_at
  

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X STORE THE GRADIENT AT THE START OF THE ALGORITHM
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine constraint_store_gradient(this)

    type(Constraint), intent(inout) :: this

    this%old_dC_dr = this%dC_dr

  end subroutine constraint_store_gradient

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X REGISTERING A CONSTRAINT
  !X
  !% Wrapper for constraint registering. Allows error catching
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function register_constraint(sub) result(p)

    integer       :: p
    interface
       subroutine sub(pos, velo, t, data, C, dC_dr, dC_dt)
         use system_module
         real(dp), dimension(:),         intent(in)  :: pos, velo, data
         real(dp),                       intent(in)  :: t
         real(dp),                       intent(out) :: C
         real(dp), dimension(size(pos)), intent(out) :: dC_dr
         real(dp),                       intent(out) :: dC_dt           
       end subroutine sub
    end interface

    if ((REGISTERED_CONSTRAINTS+1) == MAX_CONSTRAINT_SUBS) &
         call system_abort('Register_Constraint: Pointer table is full')
    call register_constraint_sub(sub)
    p = REGISTERED_CONSTRAINTS
    REGISTERED_CONSTRAINTS = REGISTERED_CONSTRAINTS + 1

  end function register_constraint

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X I/O
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine constraint_print(this,verbosity,file,index)

    type(Constraint),  intent(in) :: this
    integer, optional, intent(in) :: verbosity
    type(Inoutput), optional, intent(inout) :: file
    integer, optional, intent(in) :: index
    integer                       :: i

    call print('================', verbosity, file)
    if (present(index)) then
       write(line,'(a,i0)')' Constraint ',index
       call print(line, verbosity, file)
    else
       call print('Constraint', verbosity, file)
    end if
    call print('================', verbosity, file)  
    call print('', verbosity, file)
    if (this%initialised) then       
       write(line,'(a,i0)')'Number of atoms = ',this%N
       call print(line, verbosity, file)      
       call print('Atoms:', verbosity, file)
       call print(this%atom, verbosity, file)
       call print('Data:', verbosity, file)
       call print(this%data, verbosity, file)
       write(line,'(a,i0)')'Constraint function = ',this%func
       call print(line, verbosity, file)      
       write(line,'(a,f0.8)')'Constraint value    = ',this%C
       call print(line, verbosity, file)
       call print('Constraint gradients :')
       do i = 1, this%N
          write(line,'(3(f0.8,1x))'),this%dC_dr(3*i-2:3*i)
          call print(line, verbosity, file)
       end do
       write(line,'(a,f0.8)')'Time Derivative     = ', this%dC_dt
       call print(line, verbosity, file)
       write(line,'(a,f0.8)')'Lagrange multiplier(R) = ', this%lambdaR
       call print(line, verbosity, file)
       write(line,'(a,f0.8)')'Lagrange multiplier(V) = ', this%lambdaV
       call print(line, verbosity, file)
    else
       call print('(uninitialised)', verbosity, file)
    end if
    call print('', verbosity, file)     
    call print('================', verbosity, file)

  end subroutine constraint_print

  !
  !% Print an array of constraints, with optional first index
  !
  subroutine constraints_print(this,verbosity,file,first)

    type(Constraint), dimension(:), intent(in) :: this
    integer,              optional, intent(in) :: verbosity
    type(Inoutput),    optional, intent(inout) :: file
    integer,          optional,     intent(in) :: first
    integer                                    :: my_first, i

    my_first = 1
    if (present(first)) my_first = first

    do i = 1, size(this)
       call print(this(i),verbosity,file,i-1+my_first)
       call print('',verbosity,file)
    end do

  end subroutine constraints_print

  subroutine write_binary_constraint(this,file)

    type(Constraint), intent(in)    :: this
    type(Inoutput),   intent(inout) :: file

    !Header
    call write_binary('Constraint',file)

    !Initialisation status
    call write_binary(this%initialised,file)
    
    if(this%initialised) then

       !Other members
       call write_binary(this%N,file)          !This is size(atom), and size(dC_dr)/3
       call write_binary(this%atom,file)
       call write_binary(size(this%data),file)
       call write_binary(this%data,file)
       call write_binary(this%func,file)
       call write_binary(this%lambdaR,file)
       call write_binary(this%lambdaV,file)
       call write_binary(this%C,file)
       call write_binary(this%dC_dr,file)
       call write_binary(this%dC_dt,file)
       call write_binary(this%old_dC_dr,file)
       call write_binary(this%tol,file)

    end if

  end subroutine write_binary_constraint

  subroutine read_binary_constraint(this,file)

    type(Constraint), intent(inout) :: this
    type(Inoutput),   intent(inout) :: file
    character(10)                   :: id
    integer                         :: datasize

    !Header
    call read_binary(id,file)
    if (id /= 'Constraint') call System_Abort('read_binary_Constraint: This is not a constraint')

    call Finalise(this)

    !Initialisation status
    call read_binary(this%initialised,file)
    
    if(this%initialised) then

       !Other members
       call read_binary(this%N,file)
       allocate(this%atom(this%N),&
                this%dC_dr(3*this%N),&
                this%old_dC_dr(3*this%N))
       call read_binary(this%atom,file)
       call read_binary(datasize,file)
       allocate(this%data(datasize))
       call read_binary(this%data,file)
       call read_binary(this%func,file)
       call read_binary(this%lambdaR,file)
       call read_binary(this%lambdaV,file)
       call read_binary(this%C,file)
       call read_binary(this%dC_dr,file)
       call read_binary(this%dC_dt,file)
       call read_binary(this%old_dC_dr,file)
       call read_binary(this%tol,file)

    end if

  end subroutine read_binary_constraint

  subroutine write_binary_constraints(this,file)

    type(Constraint), dimension(:), intent(in)    :: this
    type(Inoutput),                 intent(inout) :: file
    integer                                       :: i
    
    !Header
    call write_binary('ConstArray',file)

    !Number of constraints
    call write_binary(size(this),file)

    !Write the constraints
    do i = 1, size(this)
       call write_binary(this(i),file)
    end do
    
  end subroutine write_binary_constraints

  !
  ! Allow reading into an allocatable array, as for Groups
  !
  subroutine read_binary_constraints(alloc,file)

    type(Constraint), dimension(:), allocatable, intent(inout) :: alloc
    type(Inoutput),                              intent(inout) :: file
    integer                                                    :: i, consize
    character(10)                                              :: id
    
    !Header
    call read_binary(id,file)
    if (id /= 'ConstArray') call System_Abort('read_binary_constraints: This is not an array of constraints')

    !Number of constraints
    call read_binary(consize,file)

    !Allocate the array to the correct size
    allocate(alloc(consize))
    
    do i = 1, consize
       call read_binary(alloc(i),file)
    end do
    
  end subroutine read_binary_constraints

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X SHAKE/RATTLE ROUTINES
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% This is the SHAKE part of the SHAKE/RATTLE algorithm. It
  !% adds forces to each atom which are in the directions of the 
  !% gradients of the constraint functions such that the constraints
  !% are obeyed to within the preset tolerance at the new positions.
  subroutine shake(at,g,constraints,t,dt,store_constraint_force)

    type(atoms),       intent(inout) :: at !% the atoms object
    type(group),       intent(in)    :: g  !% the constrained group
    type(constraint),  intent(inout) :: constraints(:) !% an array of all the dynamicalsystem's constraints
    real(dp),          intent(in)    :: t, dt !% the current time, and the time step
    logical, optional, intent(in)    :: store_constraint_force

    integer                             :: i,j,n, nn, Nobj, Nat, iterations, o, pos_indices(3), cf_index
    real(dp)                            :: da(3), m, df(3)
    logical                             :: converged, do_store

    Nobj = group_n_objects(g)
    Nat = group_n_atoms(g)
    do_store = optional_default(.false.,store_constraint_force)

    !Check for "constraint_force" property
    if (do_store) then
       if (get_value(at%properties,'constraint_force',pos_indices)) then
          cf_index = pos_indices(2)
          !zero the constraint forces
          do n = 1, Nat
             i = group_nth_atom(g,n)
             at%data%real(cf_index:cf_index+2,i) = 0.0_dp
          end do
       else
          call system_abort('shake: cannot find "constraint_force" property')
       end if
    end if

    !Calculate the constraint functions and store their values - the derivatives at this point will be
    !used as the directions of the forces

    do n = 1, Nobj
       i = group_nth_object(g,n)
       call constraint_calculate_values_at(constraints(i),at,t)
       constraints(i)%old_dC_dr = constraints(i)%dC_dr
       constraints(i)%dlambdaR = 0.0_dp
       constraints(i)%lambdaR = 0.0_dp
    end do

    !Make the usual velocity Verlet step
    do n = 1, Nat
       i = group_nth_atom(g,n)
       at%pos(:,i) = at%pos(:,i) + at%velo(:,i)*dt
    end do

    iterations = 0
    do  
       ! update the  accelerations, velocities, and positions using dlambdaR
       do n = 1, Nobj
          i = group_nth_object(g,n)
          do nn = 1, constraints(i)%N
             j = constraints(i)%atom(nn)
             da = -constraints(i)%dlambdaR*constraints(i)%old_dC_dr(3*nn-2:3*nn)/at%mass(j)
             at%acc(:,j) = at%acc(:,j) + da 
             at%velo(:,j) = at%velo(:,j) + 0.5_dp*dt*da
             at%pos(:,j) = at%pos(:,j) + 0.5_dp*dt*dt*da
          end do
       end do

       ! recalculate constraints at the new positions (and new time) and test for convergence
       converged = .true.
       do n = 1, Nobj
          i = group_nth_object(g,n)
          call constraint_calculate_values_at(constraints(i),at,t+dt)
          if (abs(constraints(i)%C) > constraints(i)%tol) converged = .false.
       end do

       if (converged .or. iterations > RATTLE_MAX_ITERATIONS) exit

       !Not converged, so calculate an update to the Lagrange multiplers (dlambdaR)
       do n = 1, Nobj
          i = group_nth_object(g,n)
          m = 0.0_dp
          do nn = 1, constraints(i)%N       !Loop over each particle in the constraint
             j = constraints(i)%atom(nn)
             m = m + (constraints(i)%dC_dr(3*nn-2:3*nn) .dot. constraints(i)%old_dC_dr(3*nn-2:3*nn)) / at%mass(j)
          end do
          constraints(i)%dlambdaR =  2.0_dp * constraints(i)%C / (m * dt * dt)
          constraints(i)%lambdaR = constraints(i)%lambdaR + constraints(i)%dlambdaR
       end do

       iterations = iterations + 1

    end do

    ! Store first part of the constraint force
    if (do_store) then
       do n = 1, Nobj
          i = group_nth_object(g,n)
          do nn = 1, constraints(i)%N
             j = constraints(i)%atom(nn)
             df = -constraints(i)%lambdaR*constraints(i)%old_dC_dr(3*nn-2:3*nn)
             at%data%real(cf_index:cf_index+2,j) = at%data%real(cf_index:cf_index+2,j) + df
          end do
       end do
    end if

    if (.not.converged) then
       call print('shake: could not converge in '//RATTLE_MAX_ITERATIONS//&
            ' iterations for this group:')
       call print(g)
       call system_abort('shake convergence problem')
    end if

  end subroutine shake

  !% This is the RATTLE part of the SHAKE/RATTLE algorithm. It
  !% adds forces to each atom which are in the directions of the 
  !% gradients of the constraint functions (evaluated at the NEW positions)
  !% such that the constraints are obeyed to within the preset tolerance
  !% by the new velocities.
  subroutine rattle(at,g,constraints,t,dt,store_constraint_force)

    type(atoms),       intent(inout) :: at
    type(group),       intent(in)    :: g
    type(constraint),  intent(inout) :: constraints(:)
    real(dp),          intent(in)    :: t
    real(dp),          intent(in)    :: dt
    logical, optional, intent(in)    :: store_constraint_force 

    integer                          :: i,j,n, nn, Nobj, Nat, iterations, o, cf_index, pos_indices(3)
    real(dp)                         :: da(3), m, df(3)
    logical                          :: converged, do_store

    Nobj = group_n_objects(g)
    Nat = group_n_atoms(g)
    do_store = optional_default(.false.,store_constraint_force)

    !Check for "constraint_force" property
    if (do_store) then
       if (get_value(at%properties,'constraint_force',pos_indices)) then
          cf_index = pos_indices(2)
       else
          call system_abort('rattle: cannot find "constraint_force" property')
       end if
    end if

    !Make the usual velocity Verlet step
    do n = 1, Nat
       i = group_nth_atom(g,n)
       at%velo(:,i) = at%velo(:,i) + 0.5_dp*at%acc(:,i)*dt
    end do

    !Calculate the constraint functions and store their values - the derivatives at this point will be
    !used as the directions of the forces

    do n = 1, Nobj
       i = group_nth_object(g,n)
       call constraint_calculate_values_at(constraints(i),at,t+dt)
       constraints(i)%old_dC_dr = constraints(i)%dC_dr
       constraints(i)%dlambdaV = 0.0_dp
       constraints(i)%lambdaV = 0.0_dp
    end do

    iterations = 0
    do  
       ! update the  accelerations, and velocities using dlambdaV
       do n = 1, Nobj
          i = group_nth_object(g,n)
          do nn = 1, constraints(i)%N
             j = constraints(i)%atom(nn)
             da = -constraints(i)%dlambdaV*constraints(i)%old_dC_dr(3*nn-2:3*nn)/at%mass(j)
             at%acc(:,j) = at%acc(:,j) + da 
             at%velo(:,j) = at%velo(:,j) + 0.5_dp*dt*da
          end do
       end do

       ! recalculate constraints using new velocities and test for convergence
       converged = .true.
       do n = 1, Nobj
          i = group_nth_object(g,n)
          call constraint_calculate_values_at(constraints(i),at,t+dt)
          if (abs(constraints(i)%dC_dt) > constraints(i)%tol) converged = .false.
       end do

       if (converged .or. iterations > RATTLE_MAX_ITERATIONS) exit

       !Not converged, so calculate an update to the Lagrange multiplers (dlambdaV)
       do n = 1, Nobj
          i = group_nth_object(g,n)
          m = 0.0_dp
          do nn = 1, constraints(i)%N       !Loop over each particle in the constraint
             j = constraints(i)%atom(nn)
             m = m + norm2(constraints(i)%dC_dr(3*nn-2:3*nn)) / at%mass(j)
          end do
          constraints(i)%dlambdaV =  2.0_dp * constraints(i)%dC_dt / (m * dt)
          constraints(i)%lambdaV = constraints(i)%lambdaV + constraints(i)%dlambdaV
       end do

       iterations = iterations + 1

    end do

    ! Store second poart of the constraint force
    if (do_store) then
       do n = 1, Nobj
          i = group_nth_object(g,n)
          do nn = 1, constraints(i)%N
             j = constraints(i)%atom(nn)
             df = -constraints(i)%lambdaV*constraints(i)%old_dC_dr(3*nn-2:3*nn)
             at%data%real(cf_index:cf_index+2,j) = at%data%real(cf_index:cf_index+2,j) + df
          end do
       end do
    end if

    if (.not.converged) then
       call print('rattle: could not converge in '//RATTLE_MAX_ITERATIONS//&
            ' iterations for this group:')
       call print(g)
       call system_abort('rattle convergence problem')
    end if

  end subroutine rattle

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X START OF CONSTRAINT SUBROUTINES
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X BONDLENGTH:
  !X
  !% Constrain a bond length.
  !% 'data' should contain the required bond length
  !% The minimum image convention is used.
  !%
  !% The function used is $C = |\mathbf{r}_1 - \mathbf{r}_2|^2 - d^2$
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine BONDLENGTH(pos, velo, t, data, C, dC_dr, dC_dt)

    real(dp), dimension(:),         intent(in)  :: pos, velo, data
    real(dp),                       intent(in)  :: t
    real(dp),                       intent(out) :: C
    real(dp), dimension(size(pos)), intent(out) :: dC_dr
    real(dp),                       intent(out) :: dC_dt
    !local variables                             
    real(dp), dimension(3)                       :: d

    if(size(pos) /= 6) call system_abort('BONDLENGTH: Exactly 2 atoms must be specified')
    if(size(data) /= 1) call system_abort('BONDLENGTH: "data" must contain exactly one value')

    d = pos(1:3)-pos(4:6)

    C = norm2(d) - data(1)*data(1)

    dC_dr(1:3) = 2.0_dp * d
    dC_dr(4:6) = -2.0_dp * d

    dC_dt = dC_dr .dot. velo

  end subroutine BONDLENGTH

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X RELAX_BOND:
  !X
  !% Exponentially decay a bond length towards a final value.
  !%
  !% \begin{itemize}
  !% \item data(1) = initial bond length
  !% \item data(2) = final bond length
  !% \item data(3) = initial time
  !% \item data(4) = relaxation time
  !% \end{itemize}
  !%
  !% Constraint function is $C = |\mathbf{r}_1 - \mathbf{r}_2|^2 - d^2$, where
  !% \begin{displaymath}
  !% d = d_{final} + (d_{init} - d_{final})\exp\left(-\frac{t-t_{init}}{t_{relax}}\right)
  !% \end{displaymath}
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine RELAX_BOND(pos, velo, t, data, C, dC_dr, dC_dt)

    real(dp), dimension(:),         intent(in)  :: pos, velo, data
    real(dp),                       intent(in)  :: t
    real(dp),                       intent(out) :: C
    real(dp), dimension(size(pos)), intent(out) :: dC_dr
    real(dp),                       intent(out) :: dC_dt
    !local variables                             
    real(dp)                                    :: r(3), d, diff, efact

    if(size(pos) /= 6) call system_abort('RELAX_BOND: Exactly 2 atoms must be specified')
    if(size(data) /= 4) call system_abort('RELAX_BOND: "data" must contain exactly four value')

    r = pos(1:3)-pos(4:6)
    diff = data(1) - data(2)
    efact = exp(-(t-data(3))/data(4))
    d = data(2) + diff * efact

    C = norm2(r) - d*d
    dC_dr(1:3) = 2.0_dp * r
    dC_dr(4:6) = -2.0_dp * r
    dC_dt = dC_dr .dot. velo + 2.0_dp * d * diff * efact / data(4)

  end subroutine RELAX_BOND

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X CUBIC_BOND:
  !X
  !% Constrain a bond length to the time dependent quantity
  !% \begin{displaymath}
  !%            at^3 + bt^2 + ct + d
  !% \end{displaymath}
  !% where $t$ is clamped to the range $[t_{init},t_{final}]$
  !% 
  !% \begin{itemize}
  !% \item 'data(1:4) = (/ a,b,c,d /)'
  !% \item 'data(5) =' initial time $t_{init}$
  !% \item 'data(6) =' final time $t_{final}$
  !% \end{itemize}
  !%
  !% The constraint function is:
  !% \begin{displaymath}
  !% \begin{array}{l}
  !% C = |\mathbf{r}_1 - \mathbf{r}_2|^2 - l^2 \\
  !% l = at^3 + bt^2 + ct + d
  !% \end{array}
  !% \end{displaymath}
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine CUBIC_BOND(pos, velo, t, data, C, dC_dr, dC_dt)

    real(dp), dimension(:),         intent(in)  :: pos, velo, data
    real(dp),                       intent(in)  :: t
    real(dp),                       intent(out) :: C
    real(dp), dimension(size(pos)), intent(out) :: dC_dr
    real(dp),                       intent(out) :: dC_dt
    !local variables                             
    real(dp)                                    :: r(3), t_clamped, l, dl_dt, x, x2, x3

    if(size(pos) /= 6) call system_abort('CUBIC_BOND: Exactly 2 atoms must be specified')
    if(size(data) /= 6) call system_abort('CUBIC_BOND: "data" must contain exactly six values')

    r = pos(1:3)-pos(4:6)
    t_clamped = t
    if (t_clamped > data(6)) then
       t_clamped = data(6)
    else if (t_clamped < data(5)) then
       t_clamped = data(5)
    end if

    x = t_clamped; x2 = x*x; x3 = x2 * x
    l = data(1)*x3 + data(2)*x2 + data(3)*x + data(4)
    dl_dt = 3.0_dp*data(1)*x2 + 2.0_dp*data(2)*x + data(3)

    C = norm2(r) - l*l
    dC_dr(1:3) = 2.0_dp * r
    dC_dr(4:6) = -2.0_dp * r
    dC_dt = dC_dr .dot. velo - 2.0_dp * l * dl_dt

  end subroutine CUBIC_BOND

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X PLANE
  !X
  !% Constrain a particle to the plane $ax+by+cz=d$.
  !% 'data' should contain '(/a, b, c, d/)'. Periodic boundary conditions 
  !% are not applied (the plane only cuts the 
  !% cell once, if at all).
  !%
  !% The function used is $C = \frac{\mathbf{r} \cdot \mathbf{n} - d}{|\mathbf{n}|}$
  !% where $\mathbf{n} = (a,b,c)$
  !% 
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine PLANE(pos, velo, t, data, C, dC_dr, dC_dt)

    real(dp), dimension(:),         intent(in)  :: pos, velo, data
    real(dp),                       intent(in)  :: t
    real(dp),                       intent(out) :: C
    real(dp), dimension(size(pos)), intent(out) :: dC_dr
    real(dp),                       intent(out) :: dC_dt
    !local variables
    real(dp)                                    :: d, n(3), n_hat(3)

    if(size(data) /= 4 ) call system_abort('PLANE: "data" must be of length 4')
    if(size(pos) /= 3) call system_abort('PLANE: Exactly 1 atom must be specified')

    n = data(1:3)              !Normal to the plane
    n_hat = n / norm(n)        !Unit normal
    d = data(4) / norm(n)      !Distance from plane to origin
    C = (pos .dot. n_hat) - d  !Distance of atom from plane
    dC_dr = n_hat
    dC_dt = dC_dr .dot. velo
    
  end subroutine PLANE

end module constraints_module
