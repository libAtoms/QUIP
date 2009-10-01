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
!X  Dynamical System module
!X  
!% A DynamicalSystem object contains an Atoms object, which holds infomation about
!% velocities and accelerations of each atom, scalar quantities such as
!% thermostat settings, and logical masks so that thermostatting can be applied
!% to selected atoms etc.
!% 
!% Initialise a DynamicalSystem object like this:
!%> 	call initialise(MyDS, MyAtoms)
!% which copies MyAtoms into the internal atoms structure (and so MyAtoms is not
!% required by MyDS after this call). 
!%
!% DynamicalSystem has an integrator,
!%> 	call advance_verlet(MyDS,dt,forces)
!% which takes a set of forces and integrates the equations of motion forward
!% for a time 'dt'.
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! $Id: DynamicalSystem.f95,v 1.70 2008-07-10 12:59:25 nb326 Exp $

! $Log: not supported by cvs2svn $
! Revision 1.69  2008/07/02 13:37:58  jrk33
! Fixed bug in ds_add_atom_single - missing optional arguments cannot be reshaped
!
! Revision 1.68  2008/06/16 13:03:32  jrk33
! Updated for changes to add_atoms
!
! Revision 1.67  2008/05/29 16:48:08  jrk33
! Typo for MPI_ defined
!
! Revision 1.66  2008/05/29 15:25:28  jrk33
! Bug fix to constraint force calculation - thanks Steve
!
! Revision 1.65  2008/05/15 17:27:34  jrk33
! Add travel in ds_initialise
!
! Revision 1.64  2008/05/13 15:48:36  jrk33
! Add mass property when initialising dynamical system
!
! Revision 1.63  2008/05/13 10:39:57  saw44
! Added constraint_force subroutine to return constraint forces used in the last advance_verlet call
!
! Revision 1.62  2008/05/08 14:21:06  saw44
! Fixed temperature calc in adaptive_thermostat
!
! Revision 1.61  2008/05/07 15:51:07  nb326
! Clean up ifort warnings, mostly unused variables
!
! Revision 1.60  2008/05/02 14:13:47  saw44
! Added more functionality to adaptive thermostat: thermal_taus can now vary from region to region, regions can be named, regions can be made non-adaptive. Removed some unused variables. Added some docs
!
! Revision 1.59  2008/05/01 18:48:48  saw44
! Improved adaptive thermostat - now reacts to changes in number of atoms in each region. Fixed bug in optimum_alpha. Fixed ds_assignment
!
! Revision 1.58  2008/04/29 11:30:50  saw44
! Adaptive thermostat no prints the number of atoms in each region
!
! Revision 1.57  2008/04/28 15:45:13  saw44
! Adaptive thermostat is now usable
!
! Revision 1.56  2008/04/25 14:25:50  saw44
! fixed bug in adaptive_thermostat_add_to_region
!
! Revision 1.55  2008/04/25 13:57:33  saw44
! Updated some of the adaptive thermostat routines
!
! Revision 1.54  2008/04/22 19:25:36  saw44
! Fixed bug in adaptive thermstat (atoms%mass added in wrong place)
!
! Revision 1.53  2008/04/14 18:29:25  saw44
! Added add_group_members subroutine, and now advance_verlet stores the thermostat forces if a property called "thermo_force" is present
!
! Revision 1.52  2008/03/07 23:49:29  ab686
! reordered arguments in call to add_atoms(atoms,...)
!
! Revision 1.51  2008/03/07 23:23:51  nb326
! Fix ambiguous interface in add_atoms
!
! Revision 1.50  2008/02/12 18:18:52  jrk33
! intent(in) -> intent(inout) for things that ultimately call atoms_print_xyz
!
! Revision 1.49  2008/02/04 14:07:34  jrk33
! Decoupled atomic number and mass to allow mixed isotope simulations. New mass property must be set
!
! Revision 1.48  2008/02/04 13:33:46  saw44
! Added adaptive thermostat
!
! Revision 1.47  2007/11/30 21:40:19  nb326
! make various things operate on atoms instead of dynamical system.  Should this be done more widely?
!
! Revision 1.46  2007/11/29 16:56:12  nb326
! No need for true_pos, true_velo
!
! Revision 1.45  2007/11/28 22:52:17  nb326
! Save in-sync positions and velocities into true_pos and true_velo (both corresponding to time before advance_verlet() was called if such fields exist in atoms struct
!
! Revision 1.44  2007/11/28 18:40:56  nb326
! Remove roqgue print from torque()
!
! Revision 1.43  2007/11/28 18:33:54  nb326
! Fix operator precendence in torque()
!
! Revision 1.42  2007/11/28 15:35:21  nb326
! Add optional origin argument to torque()
!
! Revision 1.41  2007/11/27 18:20:46  nb326
! Add angular_momentum() and torque().  kinetic_energy() and momentum() can operate on DS, Atoms, or raw arrays
!
! Revision 1.40  2007/11/26 17:20:04  nb326
! Add momentum() that operates on Z and velo
!
! Revision 1.39  2007/11/21 15:00:41  jrk33
! Fixed mistake in per atom average kinetic energy
!
! Revision 1.38  2007/11/21 14:55:13  jrk33
! Added per atom average kinetic energy
!
! Revision 1.37  2007/11/19 16:01:22  nb326
! Add angular_momentum() and zero_angular_momentum().  Add mass_weighted and zero_L as optional arguments to rescale_velo
!
! Revision 1.36  2007/11/14 21:42:29  nb326
! Move guts of kinetic_energy into raw_kinetic energy, which takes arrays of Z and velo
!
! Revision 1.35  2007/11/07 10:09:05  gc121
! printing dW
!
! Revision 1.34  2007/10/25 17:25:41  saw44
! decode_mpi_error -> abort_on_mpi_error
!
! Revision 1.33  2007/10/25 16:54:41  saw44
! Fixed MPI bug when using a mixture of constrained and free atoms. Removed printing of dW every time advance_verlet is called. Changed Sum(p) to Norm(p) in ds_print_status. Deallocate mpi variables in advance_verlet only when parallel=.true.
!
! Revision 1.32  2007/10/11 18:34:00  nb326
! Replace decode_mpi_error with abort_on_mpi_error
!
! Revision 1.31  2007/10/08 15:14:15  saw44
! Parallelised over groups in Advance_Verlet. Changed SHAKE/RATTLE to not zero Lagrange multipliers between calls (slight speed up)
!
! Revision 1.30  2007/10/05 14:05:15  saw44
! Allowed passing through of data table in ds_add_atoms_multiple
!
! Revision 1.29  2007/08/28 08:49:52  jrk33
! Use Trapezium rule rather than Simpsons rule to integrate work done since it is piecewise linear
!
! Revision 1.28  2007/08/22 17:21:32  jrk33
! Resync the RNG before doing Langevein thermostatting
!
! Revision 1.27  2007/08/21 09:06:53  jrk33
! Added calculation and logging of work done by thermostat and extended energy to DS
!
! Revision 1.26  2007/07/18 13:18:36  nb326
! Use new verbosity system (intro.tex v 1.13)
!
! Revision 1.25  2007/04/30 14:34:40  jrk33
! Removed call to symmetrise_postions
!
! Revision 1.24  2007/04/27 14:41:53  jrk33
! Symmetry removed
!
! Revision 1.23  2007/04/27 11:58:56  gc121
! Massive change to the way we deal with neighbour connections. We no longer follow minimum image conventions, but store all images of neighbouring atoms within a (perhaps huge) cutoff. Much code is expected to be now broken.
!
! Revision 1.22  2007/04/20 09:49:37  saw44
! Removed PURE attribute from gaussian_velocity(_component) functions
!
! Revision 1.21  2007/04/19 12:33:41  saw44
! Added gaussian_velocity and gaussian_velocity_component functions to draw a random velocity from the correct gaussian distrbution for a degree of freedom with a specified effective mass and temperature
!
! Revision 1.20  2007/04/18 12:41:09  jrk33
! Fixed a couple of doc comments
!
! Revision 1.19  2007/04/18 01:31:03  gc121
! updated to reflect changes in printing and other naming conventions
!
! Revision 1.18  2007/04/17 16:43:21  jrk33
! Standardised subroutine and function references and printing argument order.
!
! Revision 1.17  2007/04/17 09:57:19  gc121
! put copyright statement in each file
!
! Revision 1.16  2007/03/30 16:46:11  jrk33
! Modified print argument order to conform with changes to System
!
! Revision 1.15  2007/03/28 17:44:14  saw44
! Cleaned up unused variables. Added Distance_Relative_Velocity subroutine
!
! Revision 1.14  2007/03/27 14:28:40  jrk33
! Changed to use Atoms properties to store all dynamical variables. Masks changed from logical type to integer - your code will break!
!
! Revision 1.13  2007/03/12 16:55:09  jrk33
! Reformatted documentation. Corrected BAERENDS to BERENDSEN thermostat
!
! Revision 1.12  2007/03/01 13:51:45  jrk33
! Documentation comments reformatted and edited throughout. Anything starting "!(no space)%"
!  is picked up by the documentation generation script
!
! Revision 1.11  2007/02/28 15:44:30  saw44
! Added Constrain_Bond subroutine for easy use. Allowed getMomentum and ZeroMomentum to work on a subset of the atoms
!
! Revision 1.10  2007/01/24 11:21:49  saw44
! Removed unused variables
!
! Revision 1.9  2007/01/17 14:23:44  jrk33
! Added optional parallel argument to CalcDists. Defaults to serial operation. Set parallel to true when called from DynamicalSystem.advanceVerlet
!
! Revision 1.8  2007/01/09 16:57:48  saw44
! Adding constraints now checks CONSTRAINT_WARNING_TOLERANCE
!
! Revision 1.7  2006/12/14 16:54:25  saw44
! Fixed bug in handling of time-dependent constraints
!
! Revision 1.6  2006/12/13 12:07:05  saw44
! Copied over RATTLE fix from old repository
!
! Revision 1.5  2006/12/12 14:53:29  gc121
! DS assignment did not work properly, now it does
!
! Revision 1.4  2006/12/12 00:14:44  gc121
! fixed errors, now compiles
!
! Revision 1.3  2006/12/11 23:27:51  gc121
! Moved (:,N) variables from DynamicalSystem to Atoms
!
! Revision 1.2  2006/12/04 18:06:26  nb326
! Doesn't seem to need Sparse_module to compile, so remove it
!
! Revision 1.1.1.1  2006/12/04 11:11:30  gc121
! Imported sources
!
! Revision 1.46  2006/11/27 14:29:16  saw44
! Added a ZeroMomentum call to RescaleVelo for the case where the temperature is initially zero. Added DS_Amend_Constraint routine
!
! Revision 1.45  2006/11/24 12:20:29  saw44
! If a constraint is not initially obeyed then DS_Add_Constraint now prints a more helpful warning
!
! Revision 1.44  2006/11/22 19:40:05  saw44
! Added comments about velocity Verlet and for what time the velocity corresponds to throughout the routine
!
! Revision 1.43  2006/11/21 13:42:02  saw44
! Added time dependence to constraints. Moved the position of RATTLE algorithm in AdvanceVerlet, hopefully for improved stability
!
! Revision 1.42  2006/11/17 13:32:30  saw44
! ADDED CONSTRAINTS. Made lots of changes to ReadB, WriteB, the DynamicalSystem type, partitioned AdvanceVerlet into a few smaller routines and integration now occurs by group rather than by atom
!
! Revision 1.41  2006/08/07 11:00:46  saw44
! Fixed bug in AdvanceVerlet: all components of random force in Langevin thermostat were equal
!
! Revision 1.40  2006/06/29 11:08:34  jrk33
! Removed unused variables, changed real(1e-4,dp) to 1e-4_dp in initialisation statement
!
! Revision 1.39  2006/06/28 17:26:50  saw44
! Added constrained dynamics: New object DSConstraints. See comments for use and constraint.f95 for sample program. Added centre_of_mass_velocity function.
!
! Revision 1.38  2006/06/22 18:54:53  jrk33
! Bug fixed in advanceVerlet for symmetric systems: atoms on the edge of cells with non zero travels were not being correctly remapped
!
! Revision 1.37  2006/06/20 17:23:18  gc121
! added new copyright notice to include James, Gian, Mike and Alessandro
!
! Revision 1.36  2006/06/14 10:30:41  saw44
! Renamed Constraint_Tol to Within_Constraint_Tol (name clash with new variables)
!
! Revision 1.35  2006/06/08 13:58:22  saw44
! Added Initialise and Finalise interfaces (as with other modules). Changed calculation of forcesum to sum(force,dim = 2) for possible compiler optimisation
!
! Revision 1.34  2006/06/08 13:28:55  jrk33
! Bug fix to advanceVerlet for symmetric systems: atoms that had been remapped by May_Into_Cell were not having their displacements correctly taken into account. Now we use realpos instead of pos
!
! Revision 1.33  2006/05/30 11:10:13  jrk33
! Removed declarations for unused variables
!
! Revision 1.32  2006/05/25 11:02:17  jrk33
! Shortened long lines to allow compilation with Sun f95
!
! Revision 1.31  2006/05/12 13:58:12  jrk33
! Moved force symmetrisation out of DS to hybrid_model
!
! Revision 1.30  2006/05/10 10:42:00  jrk33
! Added symmetrisation to advanceVerlet if atoms%symmetrise /= 0
!
! Revision 1.29  2006/05/05 09:02:06  saw44
! Fixed bug found by Gian in DS_File_Read
!
! Revision 1.28  2006/03/29 10:51:10  saw44
! Included dissipative term in fix of langevin thermostat
!
! Revision 1.27  2006/03/28 09:54:04  saw44
! Fixed langevin thermostat non-zero random momenta problem. needs testing
!
! Revision 1.26  2006/03/10 00:41:44  gc121
! fixed bugs in addatoms, Remove_Atoms: forgot to update this%N, flipped one loop exit condition, made do loop go from 1,N, it was not terminated before.
!
! Revision 1.25  2006/03/03 16:35:52  saw44
! Removed SAVE - not needed
!
! Revision 1.24  2006/02/28 17:00:32  saw44
! Typo correction
!
! Revision 1.23  2006/02/15 14:16:54  saw44
! Removed commas before first entry in write statements
!
! Revision 1.22  2006/02/06 16:48:22  saw44
! General Code clean-up: routine names changed to match the changes in System and linearalgebra
!
! Revision 1.21  2006/01/31 15:41:34  gc121
! added a ds_print() to initialise
!
! Revision 1.20  2006/01/31 14:28:19  saw44
! Updated ReadB and WriteB argument order
!
! Revision 1.19  2006/01/31 12:37:17  saw44
! Pointer target now properly deallocated
!
! Revision 1.18  2006/01/30 13:04:36  gc121
! correct format statement
!
! Revision 1.17  2006/01/30 12:05:43  gc121
! removeds *s from some writes
!
! Revision 1.16  2006/01/30 11:41:13  gc121
! swapped arguments of is_in_array
!
! Revision 1.15  2006/01/30 10:38:52  gc121
! removed logger from print lines, its the default
!
! Revision 1.14  2006/01/27 16:12:30  gc121
! added status line. rescaleVelo works with zero velocities. added zeroMomentum. fixed default timescales to be in our units
!
! Revision 1.13  2006/01/26 16:10:44  gc121
! added verbosity to printing, fixed function names
!
! Revision 1.12  2006/01/26 01:53:28  gc121
! corrected typo
!
! Revision 1.11  2006/01/25 17:35:27  gc121
! changed vector_abs and vector_abs2 to norm() and norm2()
!
! Revision 1.10  2006/01/25 16:09:44  gc121
! changes some argument names
!
! Revision 1.9  2006/01/24 12:07:31  saw44
! DS now compiles with some of the discussed changes
!
! Revision 1.8  2006/01/18 16:07:59  gc121
! cleanup started by gc121 and saw44
!

module dynamicalsystem_module
 
   use system_module
   use table_module
   use atoms_module
   use periodictable_module
   use units_module
   use linearalgebra_module
   use rigidbody_module
   use group_module
   use constraints_module
   use thermostat_module
   
   implicit none

   !Different integration types. Stored in group%type. Used by advance_verlet to integrate
   !the equations of motion in different ways. 
   integer, parameter :: TYPE_IGNORE      = 0
   integer, parameter :: TYPE_ATOM        = 1
   integer, parameter :: TYPE_CONSTRAINED = 2
   integer, parameter :: TYPE_RIGID       = 3
   real(dp), parameter :: min_ext_p = 1.0e-4_dp / GPA

   type DynamicalSystem 

      ! Scalar members
      integer                               :: N = 0                    !% Number of atoms
      integer                               :: nSteps = 0               !% Number of integration steps
      integer                               :: Nrigid = 0               !% Number of rigid bodies
      integer                               :: Nconstraints = 0         !% Number of constraints
      integer                               :: Ndof = 0                 !% Number of degrees of freedom
      real(dp)                              :: t = 0.0_dp               !% Time
      real(dp)                              :: avg_temp = 0.0_dp        !% Time-averaged temperature
      real(dp)                              :: cur_temp = 0.0_dp        !% Current temperature
      real(dp)                              :: avg_time = 100.0_dp      !% Averaging time, in fs
      real(dp)                              :: dW = 0.0_dp              !% Increment of work done this time step
      real(dp)                              :: work = 0.0_dp            !% Total work done
      real(dp)                              :: Epot = 0.0_dp            !% Total potential energy
      real(dp)                              :: ext_energy = 0.0_dp      !% Extended energy
      real(dp)                              :: thermostat_dW = 0.0_dp   !% Increment of work done by thermostat
      real(dp)                              :: thermostat_work = 0.0_dp !% Total work done by thermostat
      logical         			    :: initialised = .false.

      integer                               :: random_seed              !% RNG seed, used by 'ds_save_state' and 'ds_restore_state' only.

      ! Array members
      integer,  allocatable, dimension(:)   :: group_lookup !% Stores which group atom $i$ is in

      ! Derived type members
      type(Atoms)                           :: atoms

      ! Derived Type array members
      type(Constraint), allocatable, dimension(:) :: constraint
      type(RigidBody),  allocatable, dimension(:) :: rigidbody
      type(Group),      allocatable, dimension(:) :: group
      type(thermostat), allocatable, dimension(:) :: thermostat

   end type DynamicalSystem

   interface assignment(=)
      module procedure DS_Assignment
   end interface

   !% Add one or more atoms to this DynamicalSystem
   interface add_atoms
      module procedure ds_add_atom_single, ds_add_atom_multiple
   end interface add_atoms

   !% Remove one or more atoms from this DynamicalSystem
   interface remove_atoms
      module procedure ds_remove_atom_single, ds_remove_atom_multiple
   end interface remove_atoms

   !% Write this DynamicalSystem to a binary file
   interface write_binary
      module procedure ds_file_write
   end interface

   !% Read this DynamicalSystem from a binary file
   interface read_binary
      module procedure ds_file_read
   end interface

   interface print
      module procedure ds_print
   end interface

   !% Initialise this DynamicalSystem from an Atoms object
   interface initialise
      module procedure ds_initialise
   end interface initialise

   !% Free up the memory used by this DynamicalSystem
   interface finalise
      module procedure ds_finalise
   end interface finalise

   interface kinetic_energy
     module procedure raw_kinetic_energy, atoms_kinetic_energy, ds_kinetic_energy
   end interface kinetic_energy

   interface angular_momentum
     module procedure raw_angular_momentum, atoms_angular_momentum, ds_angular_momentum
   end interface angular_momentum

   interface momentum
     module procedure raw_momentum, atoms_momentum, ds_momentum
   end interface momentum

   interface add_thermostat
      module procedure ds_add_thermostat
   end interface add_thermostat

contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X INITIALISE AND FINALISE
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine ds_initialise(this,atoms_in,velocity,acceleration,constraints,rigidbodies)
   
      type(DynamicalSystem),              intent(inout) :: this
      type(Atoms),                        intent(in)    :: atoms_in
      real(dp), dimension(:,:), optional, intent(in)    :: velocity
      real(dp), dimension(:,:), optional, intent(in)    :: acceleration
      integer,                  optional, intent(in)    :: constraints
      integer,                  optional, intent(in)    :: rigidbodies
      integer                                           :: i,N
      
      ! Check to see if the object has already been initialised

      if (this%initialised) call ds_finalise(this)  

      N = atoms_in%N
      this%N = N

      ! Check the atomic numbers
      do i = 1, atoms_in%N
         if (atoms_in%Z(i) < 1) then
            write(line,'(a,i0,a)')'DS_Initialise: Atom ',i,' has not had its atomic number set'
            call system_abort(line)
         end if
      end do

      ! allocate local arrays
      allocate(this%group_lookup(N),   &
               this%group(N)) 

      ! Now copy Atoms   
      this%atoms = atoms_in

      ! Add properties for the dynamical variables to this%atoms if they don't
      ! already exist. This will update pointers in this%atoms so we can use 
      ! this%atoms%velo etc. afterwards.

      call add_property(this%atoms, 'mass', ElementMass(this%atoms%Z))
      call add_property(this%atoms, 'travel', 0, n_cols=3)

      call add_property(this%atoms, 'move_mask', 1)
      call add_property(this%atoms, 'damp_mask', 1)
      call add_property(this%atoms, 'thermostat_region', 1)

      call add_property(this%atoms, 'avg_ke', 0.0_dp)
      call add_property(this%atoms, 'velo', 0.0_dp, n_cols=3)
      call add_property(this%atoms, 'acc', 0.0_dp, n_cols=3)
      call add_property(this%atoms, 'avgpos', 0.0_dp, n_cols=3)
      call add_property(this%atoms, 'oldpos', 0.0_dp, n_cols=3)

      ! The input arrays must have 3N components if they are present
      ! if not, local arrays are set to zero
      
      if(present(velocity)) then 
         call check_size('Velocity',velocity,(/3,N/),'DS_Initialise')
         this%atoms%velo = velocity
      end if

      if(present(acceleration)) then
         call check_size('Acceleration',acceleration,(/3,N/),'DS_Initialise')
         this%atoms%acc = acceleration
      end if

      ! Initialize oldpos and avgpos
      do i = 1, this%N
         this%atoms%avgpos(:,i) = realpos(this%atoms,i)
      end do
      this%atoms%oldpos = this%atoms%avgpos
      do i=1, this%atoms%N
	this%atoms%avg_ke(i) = 0.5_dp*this%atoms%mass(i)*norm2(this%atoms%velo(:,i))
      end do

      ! Check for constraints
      if (present(constraints)) then
         if (constraints > 0) then
            allocate(this%constraint(constraints))
            this%Nconstraints = 0
         end if
      end if

      ! Check for rigid bodies
      if (present(rigidbodies)) then
         if (rigidbodies > 0) then
            allocate(this%rigidbody(rigidbodies))
            this%Nrigid = 0
         end if
      end if

      ! Initialise all the groups, default to TYPE_ATOM
      do i = 1, this%N
         call initialise(this%group(i),TYPE_ATOM, atoms = (/i/))
         this%group_lookup(i) = i
      end do

      this%Ndof = 3 * this%N

      allocate(this%thermostat(0:0))
      call initialise(this%thermostat(0),NONE,0.0_dp) !a dummy thermostat, which is turned into
                                                      !a langevin thermostat if damping is needed
      this%cur_temp = temperature(this, include_all=.true., instantaneous=.true.)
      this%avg_temp = this%cur_temp

      this%initialised = .true.
      
      call verbosity_push_decrement(ANAL)
      call print(this)
      call verbosity_pop()
         
   end subroutine ds_initialise


   subroutine ds_finalise(this)

      type(DynamicalSystem), intent(inout) :: this

      if (this%initialised) then      
         call finalise(this%atoms)
         call finalise(this%group)
         deallocate(this%group_lookup)
         
         !Finalise constraints
         if (allocated(this%constraint)) call finalise(this%constraint)

         !Finalise rigid bodies
         if (allocated(this%rigidbody)) call finalise(this%rigidbody)

         call finalise(this%thermostat)

         this%N = 0
         this%nSteps = 0
         this%Nrigid = 0
         this%Nconstraints = 0
         this%Ndof = 0
         this%t = 0.0_dp
         this%avg_temp = 0.0_dp
         this%dW = 0.0_dp
         this%work = 0.0_dp
         this%thermostat_dW = 0.0_dp
         this%thermostat_work = 0.0_dp
         this%ext_energy = 0.0_dp
         this%Epot = 0.0_dp

         this%initialised = .false.

      end if

   end subroutine ds_finalise

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !X
   !X FREE EXCESS MEMORY USED BY GROUPS
   !X
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   !% Free up excess memory used by Groups
   subroutine ds_free_groups(this)

     type(DynamicalSystem), intent(inout) :: this
     
     call tidy_groups(this%group)
     call groups_create_lookup(this%group,this%group_lookup)

   end subroutine ds_free_groups

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !X
   !X ASSIGNMENT
   !X
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   !% Overloaded assignment operator for DynamicalSystems
   subroutine ds_assignment(to,from)

      type(DynamicalSystem), intent(inout) :: to
      type(DynamicalSystem), intent(in)    :: from

      integer :: from_constraint_size, from_rigidbody_size

      ! Re-initialise the destination object to be the correct size.
      ! Note: size(from%constraint) in general is not equal to from%Nconstraints etc.
      ! Also, the constraints aren't added using add_constraint etc. since this
      ! will automatically evaluate the constraints at the current positions, which may not
      ! correspond to the values currently stored. Same for rigid bodies.

      from_constraint_size = 0
      if (allocated(from%constraint)) from_constraint_size = size(from%constraint)
      from_rigidbody_size = 0
      if (allocated(from%rigidbody)) from_rigidbody_size = size(from%rigidbody)

      call initialise(to, from%atoms, constraints=from_constraint_size, rigidbodies=from_rigidbody_size)

      ! Copy over scalar members
      ! to%N is set in the initialisation         1
      to%nSteps          = from%nSteps           !2
      to%Nrigid          = from%Nrigid           !4
      to%Nconstraints    = from%Nconstraints     !5
      to%Ndof            = from%Ndof             !6

      to%t               = from%t                !8
      to%avg_temp        = from%avg_temp          !11
      to%avg_time        = from%avg_time         !12
      to%dW              = from%dW               !18
      to%work            = from%work             !19
      to%Epot            = from%Epot             !20
      to%ext_energy      = from%ext_energy       !21
      to%thermostat_dW   = from%thermostat_dW    !22
      to%thermostat_work = from%thermostat_work  !23
      ! to%initialised is set in initialisation   

      ! Copy over array members
      to%group_lookup    = from%group_lookup

      ! Derived type members
      ! to%atoms already set
      if (from%Nconstraints /= 0) to%constraint      = from%constraint
      if (from%Nrigid /= 0)       to%rigidbody       = from%rigidbody
      if (size(to%group)/=size(from%group)) then
         deallocate(to%group)
         allocate(to%group(size(from%group)))
      end if
      to%group           = from%group

      ! Copy thermostats - will have been wiped by ds_initialise
      ! (uses overloaded routine to copy array of thermostats)
      ! to%thermostat = from%thermostat
      call thermostat_array_assignment(to%thermostat, from%thermostat)

   end subroutine ds_assignment

   !% Save the state of a DynamicalSystem. The output object
   !% cannot be used as an initialised DynamicalSystem since
   !% connectivity and group information is not copied to save
   !% memory. Only scalar members and the 'ds%atoms' object
   !% (minus 'ds%atoms%connect') are copied. The current
   !% state of the random number generator is also saved.
   subroutine ds_save_state(to, from)

      type(DynamicalSystem), intent(inout) :: to
      type(DynamicalSystem), intent(in)    :: from

      ! Copy over scalar members
      to%N               = from%N                !1
      to%nSteps          = from%nSteps           !2
      to%Nrigid          = from%Nrigid           !4
      to%Nconstraints    = from%Nconstraints     !5
      to%Ndof            = from%Ndof             !6

      to%t               = from%t                !8
      to%avg_temp        = from%avg_temp         !11
      to%avg_time        = from%avg_time         !12
      to%dW              = from%dW               !18
      to%work            = from%work             !19
      to%Epot            = from%Epot             !20
      to%ext_energy      = from%ext_energy       !21
      to%thermostat_dW   = from%thermostat_dW    !22
      to%thermostat_work = from%thermostat_work  !23

      to%random_seed     = system_get_random_seed()

      call atoms_copy_without_connect(to%atoms, from%atoms)

    end subroutine ds_save_state

    !% Restore a DynamicalSystem to a previously saved state.
    !% Only scalar members and 'ds%atoms' (minus 'ds%atoms%connect')
    !% are copied back; 'to' should be a properly initialised
    !% DynamicalSystem object. The saved state of the random
    !% number generator is also restored. 'calc_dists()' is
    !% called on the restored atoms object.
    subroutine ds_restore_state(to, from)

      type(DynamicalSystem), intent(inout) :: to
      type(DynamicalSystem), intent(in)    :: from

      ! Copy over scalar members
      to%N               = from%N                !1
      to%nSteps          = from%nSteps           !2
      to%Nrigid          = from%Nrigid           !4
      to%Nconstraints    = from%Nconstraints     !5
      to%Ndof            = from%Ndof             !6

      to%t               = from%t                !8
      to%avg_temp        = from%avg_temp          !11
      to%avg_time        = from%avg_time         !12
      to%dW              = from%dW               !18
      to%work            = from%work             !19
      to%Epot            = from%Epot             !20
      to%ext_energy      = from%ext_energy       !21
      to%thermostat_dW   = from%thermostat_dW    !22
      to%thermostat_work = from%thermostat_work  !23

      call system_reseed_rng(from%random_seed)

      call atoms_copy_without_connect(to%atoms, from%atoms)      
      call calc_dists(to%atoms)
      
    end subroutine ds_restore_state


   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !X
   !X ADDING / REMOVING ATOMS
   !X
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine ds_add_atom_single(this,z,mass,p,v,a,t)
     
     type(DynamicalSystem),            intent(inout) :: this
     integer,                          intent(in)    :: z
     real(dp), optional,               intent(in)    :: mass
     real(dp), optional, dimension(3), intent(in)    :: p,v,a
     integer,  optional, dimension(3), intent(in)    :: t

     real(dp) :: my_mass, my_p(3), my_v(3), my_a(3)
     integer :: my_t(3)
          
     my_mass = optional_default(ElementMass(Z),mass) 
     my_p = optional_default((/0.0_dp,0.0_dp,0.0_dp/),p)
     my_v = optional_default((/0.0_dp,0.0_dp,0.0_dp/),v)
     my_a = optional_default((/0.0_dp,0.0_dp,0.0_dp/),a)
     my_t = optional_default((/0,0,0/),t)

     call ds_add_atom_multiple(this, (/z/), (/my_mass/), &
                               reshape(my_p, (/3,1/)), &
                               reshape(my_v, (/3,1/)), &
                               reshape(my_a, (/3,1/)), &
                               reshape(my_t, (/3,1/)))

   end subroutine ds_add_atom_single
 
   !
   ! Updated to work with groups. NEEDS TESTING
   !
   subroutine ds_add_atom_multiple(this,z,mass,p,v,a,t,data)
     
     type(DynamicalSystem),              intent(inout) :: this
     integer,  dimension(:),             intent(in)    :: z
     real(dp), dimension(:),   optional, intent(in)    :: mass
     real(dp), dimension(:,:), optional, intent(in)    :: p,v,a
     integer, dimension(:,:),  optional, intent(in)    :: t
     type(table),              optional, intent(in)    :: data
     
     integer                                 :: oldN, newN, n, f, i
     integer,     dimension(this%N)          :: tmp_group_lookup
     type(Group), dimension(:), allocatable  :: tmp_group
     
     oldN = this%N

     ! Check the sizes are ok
     if (present(data)) then

        n = data%N
        newN = oldN + n

        if (this%atoms%data%intsize /= data%intsize) &
             call system_abort('Add_Atoms: this%data%intsize /= data%intsize')
        if (this%atoms%data%realsize /= data%realsize) &
             call system_abort('Add_Atoms: this%data%realsize /= data%realsize')
     else        

        n = size(z)
        call check_size('Position',p,(/3,n/),'DS_Add_Atoms') 
        call check_size('Velocity',v,(/3,n/),'DS_Add_Atoms') 
        call check_size('Acceleration',a,(/3,n/),'DS_Add_Atoms') 
        
        newN = oldN + n

     end if       
     
     !Copy all non-scalar data into the temps
     tmp_group_lookup = this%group_lookup
     
     !Adjust the size of the dynamical system's arrays
     deallocate(this%group_lookup)
     allocate(this%group_lookup(newN))
     
     ! Implement the changes in Atoms
     call add_atoms(this%atoms,pos=p,Z=z,mass=mass,velo=v,acc=a,travel=t,data=data)
     
     !update the scalars (if needed)
     this%N = newN
     this%Ndof = this%Ndof + 3*n
     
     !See if there are enough free groups to accommodate the new atoms
     f = Num_Free_Groups(this%group)
     if (f < n) then
        !We need more groups: 
        !Make a copy of the current groups
        allocate(tmp_group(size(this%group)))
        tmp_group = this%group
        !Resize the group array
        call finalise(this%group)
        allocate(this%group(size(tmp_group) + n - f))
        !Copy the groups back
        this%group(1:size(tmp_group)) = tmp_group
        call finalise(tmp_group)
     end if
     
     !Add the new groups
     do i = oldN+1, newN
        f = Free_Group(this%group)
        call group_add_atom(this%group(f),i)
        call set_type(this%group(f),TYPE_ATOM) !default to TYPE_ATOM
     end do
     
     !Rebuild the lookup table
     call groups_create_lookup(this%group,this%group_lookup)
     
   end subroutine ds_add_atom_multiple


   subroutine ds_remove_atom_single(this,i)

      type(DynamicalSystem), intent(inout) :: this
      integer,               intent(in)    :: i

      call ds_remove_atom_multiple(this,(/i/))

   end subroutine ds_remove_atom_single



   subroutine ds_remove_atom_multiple(this,atomlist_in)

      type(DynamicalSystem),  intent(inout)    :: this
      integer,  dimension(:), intent(in)       :: atomlist_in

      integer,  dimension(size(atomlist_in))   :: atomlist
      integer                                  :: oldN, newN, g, i, copysrc

      !Make our own copy of the indices so that we can sort them
      atomlist = atomlist_in
      call insertion_sort(atomlist)

      !Check for repeated indices, and non-TYPE_ATOM atoms
      do i = 1, size(atomlist)
         if (Atom_Type(this,atomlist(i)) /= TYPE_ATOM) then
            write(line,'(a,i0,a)')'Remove_Atoms: Atom ',atomlist(i),' is not a normal atom'
            call system_abort(line)
         end if
         if (i > 1) then
            if (atomlist(i) == atomlist(i-1)) then
               write(line,'(a,i0,a)')'Remove_Atoms: Tried to remove the same atom twice (',atomlist(i),')'
               call system_abort(line)
            end if
         end if
      end do

      oldN = this%N
      newN = this%N - size(atomlist)

      ! Implement the data changes in the atoms structure
      ! Since we're mangling the atom indices and some atoms are being removed, the
      ! connection data will no longer be valid after this
      call remove_atoms(this%atoms,atomlist)

      !Make sure the group lookup table is up-to-date
      call groups_create_lookup(this%group,this%group_lookup)      
 
      !Delete all atoms in atomlist from their respective groups
      do i = 1, size(atomlist)
         g = this%group_lookup(atomlist(i))
         call group_delete_atom(this%group(g),atomlist(i))
      end do
 
      ! Algorithm: Find first atom to be removed and last atom to not be removed
      !            and swap them. Repeat until all atoms to be removed are at the end.
      ! note: this loop must be logically identical to the corresponding one in 
      ! Atoms
      copysrc = oldN
      do i=1,size(atomlist)

         do while(Is_in_array(atomlist,copysrc))
            copysrc = copysrc - 1
         end do

         if (atomlist(i) > copysrc) exit
                  
         !Relabel copysrc to atomlist(i) in copysrc's group
         
         g = this%group_lookup(copysrc)
         call group_delete_atom(this%group(g),copysrc)
         call group_add_atom(this%group(g),atomlist(i))
         
         copysrc = copysrc - 1

      end do

      !Resize dynamical system's arrays
      deallocate(this%group_lookup)
      allocate(this%group_lookup(newN))

      !Update scalars
      this%N = newN
      this%Ndof = this%Ndof - 3*size(atomlist)

      !Rebuild the group lookup table
      call groups_create_lookup(this%group,this%group_lookup)

   end subroutine ds_remove_atom_multiple

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !X 
   !X Atom - Group lookup
   !X
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   !
   !% Look up the type of the group that atom i belongs to
   !
   pure function atom_type(this,i)

     type(DynamicalSystem), intent(in) :: this
     integer,               intent(in) :: i
     integer                           :: atom_type

     atom_type = this%group(this%group_lookup(i))%type

   end function atom_type

   !
   !% If an atom is in the list, add the rest of its group
   !
   subroutine add_group_members(this,list)
     
     type(dynamicalsystem), intent(in)    :: this
     type(table),           intent(inout) :: list
     
     integer :: i, j, n, g, nn, jn
     
     do n = 1, list%N
        
        i = list%int(1,n)
        g = this%group_lookup(i)
        do nn = 1, Group_N_Atoms(this%group(g))
           j = Group_Nth_Atom(this%group(g),nn)
           jn = find_in_array(int_part(list,1),j)
           if (jn == 0) call append(list,(/j,list%int(2:,n)/))
        end do
        
     end do

   end subroutine add_group_members

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !X 
   !X Adding a thermostat
   !X
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine ds_add_thermostat(this,type,T,gamma,Q,tau,p)

     type(dynamicalsystem), intent(inout) :: this
     integer,               intent(in)    :: type
     real(dp),              intent(in)    :: T
     real(dp), optional,    intent(in)    :: gamma
     real(dp), optional,    intent(in)    :: Q
     real(dp), optional,    intent(in)    :: tau
     real(dp), optional,    intent(in)    :: p

     real(dp) :: w_p, cell_gamma, mass1,mass2
     real(dp) :: gamma_eff

     if (count( (/present(gamma), present(tau) /) ) /= 1 ) call system_abort('ds_add_thermostat: exactly one of gamma, tau must be present')

     if (present(gamma)) then
       gamma_eff = gamma
     else
       gamma_eff = 1.0_dp/tau
     endif

     if(present(p)) then
        cell_gamma = gamma_eff * 0.1_dp

        mass1 = 9.0_dp*abs(p)*cell_volume(this%atoms)/((cell_gamma*2*PI)**2)
        mass2 = (this%Ndof+3.0_dp)*T/((cell_gamma*2*PI)**2)

        w_p = max(mass1,mass2)
     endif

     call add_thermostat(this%thermostat,type,T,gamma_eff,Q,p,0.1_dp*gamma_eff,w_p)
     
   end subroutine ds_add_thermostat

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !X 
   !X enable/disable damping
   !X
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine enable_damping(this,damp_time)

     type(dynamicalsystem), intent(inout) :: this
     real(dp),              intent(in)    :: damp_time
     
     if (damp_time <= 0.0_dp) call system_abort('enable_damping: damp_time must be > 0')
     call initialise(this%thermostat(0),LANGEVIN,0.0_dp,gamma=1.0_dp/damp_time)
     
   end subroutine enable_damping

   subroutine disable_damping(this)

     type(dynamicalsystem), intent(inout) :: this    
     call initialise(this%thermostat(0),NONE)
     
   end subroutine disable_damping

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !X
   !X Procedures for returning/changing the state of the system
   !X
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   !% Return the total momentum $\mathbf{p} = \sum_i \mathbf{m_i} \mathbf{v_i}$.
   !% Optionally only include the contribution of a subset of atoms.
   pure function ds_momentum(this,indices) result(p) ! sum(mv) 
     type(DynamicalSystem),           intent(in) :: this
     integer, optional, dimension(:), intent(in) :: indices
     real(dp)                                    :: p(3)

     p = momentum(this%atoms, indices)
   end function ds_momentum

   !% Return the total momentum $\mathbf{p} = \sum_i \mathbf{m_i} \mathbf{v_i}$.
   !% Optionally only include the contribution of a subset of atoms.
   pure function atoms_momentum(this,indices) result(p) ! sum(mv) 
     type(Atoms),                     intent(in) :: this
     integer, optional, dimension(:), intent(in) :: indices
     real(dp)                                    :: p(3)

     p = momentum(this%mass, this%velo, indices)
   end function atoms_momentum

   !% Return the total momentum $\mathbf{p} = \sum_i \mathbf{m_i} \mathbf{v_i}$.
   !% Optionally only include the contribution of a subset of atoms.
   pure function raw_momentum(mass, velo, indices) result(p) ! sum(mv) 
     real(dp), intent(in)                         :: mass(:)
     real(dp), intent(in)                        :: velo(:,:)
     integer, optional, dimension(:), intent(in) :: indices

     real(dp)                                    :: p(3)
     integer                                     :: i, N

     N = size(mass)
     p = 0.0_dp
     if (present(indices)) then
        do i = 1,size(indices)
           p = p + velo(:,indices(i)) * mass(indices(i))
        end do
     else
        do i = 1,N
           p = p + velo(:,i) * mass(i)
        end do
     end if
   end function raw_momentum

   !% Return the angular momentum of all the atoms in this DynamicalSystem, defined by
   !% $\mathbf{L} = \sum_{i} \mathbf{r_i} \times \mathbf{v_i}$.
   pure function ds_angular_momentum(this, origin, indices) result(L)
     type(DynamicalSystem), intent(in) :: this
     real(dp), intent(in), optional :: origin(3)
     integer, intent(in), optional :: indices(:)
     real(dp) :: L(3)

     L = angular_momentum(this%atoms, origin, indices)
   end function ds_angular_momentum

   !% Return the angular momentum of all the atoms in this DynamicalSystem, defined by
   !% $\mathbf{L} = \sum_{i} \mathbf{r_i} \times \mathbf{v_i}$.
   pure function atoms_angular_momentum(this, origin, indices) result(L)
     type(Atoms), intent(in) :: this
     real(dp), intent(in), optional :: origin(3)
     integer, intent(in), optional :: indices(:)
     real(dp) :: L(3)

     L = angular_momentum(this%mass, this%pos, this%velo, origin, indices)
   end function atoms_angular_momentum

   !% Return the angular momentum of all the atoms in this DynamicalSystem, defined by
   !% $\mathbf{L} = \sum_{i} \mathbf{r_i} \times \mathbf{v_i}$.
   pure function raw_angular_momentum(mass, pos, velo, origin, indices) result(L)
     real(dp), intent(in) :: mass(:)
     real(dp), intent(in) :: pos(:,:), velo(:,:)
     real(dp), intent(in), optional    :: origin(3)
     integer, intent(in), optional :: indices(:)
     real(dp) :: L(3)

     integer                           :: ii, i, N
     
     if (present(indices)) then
       N = size(indices)
     else
       N = size(mass)
     endif

     L = 0.0_dp
     do ii = 1, N
	if (present(indices)) then
	  i = indices(ii)
	else
	  i = ii
	end if
	if (present(origin)) then
	  L = L + mass(i) * ((pos(:,i)-origin) .cross. velo(:,i))
	else
	  L = L + mass(i) * (pos(:,i) .cross. velo(:,i))
	endif
     end do

   end function raw_angular_momentum

   !% Return the moment of inertia of all the atoms in this DynamicalSystem about a particular axis
   !% placed at (optional) origin
   !% $I = \sum_{i} m_i |\mathbf{r_i}-\mathbf{o}|^2 $.
   pure function moment_of_inertia(this,axis,origin) result (MoI)
     type(Atoms), intent(in) :: this
     real(dp), intent(in)              :: axis(3)
     real(dp), intent(in), optional    :: origin(3)
     real(dp) :: MoI

     real(dp) :: my_origin(3)
     integer i
     real(dp) :: dr_proj(3), dr(3), dr_normal(3)
     real(dp) :: axis_hat(3)

     if (present(origin)) then
       my_origin = origin
     else
       my_origin = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
     end if

     axis_hat = axis/norm(axis)

     MoI = 0.0_dp
     do i=1, this%N
       dr = (this%pos(:,i)-my_origin)
       dr_proj = (dr.dot.axis_hat)*axis_hat
       dr_normal = dr - dr_proj
       MoI = MoI + this%mass(i)*norm2(dr_normal)
     end do
   end function moment_of_inertia

   pure function moment_of_inertia_tensor(this,origin) result (MoI)
     type(Atoms), intent(in) :: this
     real(dp), intent(in), optional    :: origin(3)
     real(dp) :: MoI(3,3)

     real(dp) :: dr(3)
     real(dp) :: m
     integer i, ii, jj

     MoI = 0.0_dp
     do i=1, this%N
       m = this%mass(i)
       if (present(origin)) then
	 dr = this%pos(:,i)-origin
       else
	 dr = this%pos(:,i)
       endif
       do ii=1,3
       do jj=1,3
	 if (ii == jj) MoI(ii,jj) = MoI(ii,jj) + m*norm2(dr)
	 MoI(ii,jj) = MoI(ii,jj) - m*dr(ii)*dr(jj)
       end do
       end do
     end do
   end function moment_of_inertia_tensor

   !% Return the total kinetic energy $E_k = \sum_{i} \frac{1}{2} m v^2$
   function DS_kinetic_energy(this) result(ke)       ! sum(0.5mv^2)
      type(DynamicalSystem), intent(in) :: this
      real(dp)                          :: ke

      ke = kinetic_energy(this%atoms)
    end function DS_kinetic_energy

   !% Return the total kinetic energy $E_k = \sum_{i} \frac{1}{2} m v^2$
   function atoms_kinetic_energy(this) result(ke)       ! sum(0.5mv^2)
      type(atoms), intent(in) :: this
      real(dp)                :: ke

      if (.not. associated(this%mass)) call system_abort("atoms_kinetic_energy called on atoms without mass property")
      if (.not. associated(this%velo)) call system_abort("atoms_kinetic_energy called on atoms without velo property")
      ke = kinetic_energy(this%mass, this%velo)
    end function atoms_kinetic_energy

   !% Return the total kinetic energy given atomic numbers and velocities
   pure function raw_kinetic_energy(mass, velo) result(ke)
     real(dp), intent(in) :: mass(:)
     real(dp), intent(in) :: velo(:,:)
     real(dp) :: ke

     integer i, N

     N = size(mass)

     ke = 0.0_dp
     do i = 1,N
        ke = ke + mass(i) * norm2(velo(:,i))
     end do
     ke = 0.5_dp * ke      
   end function raw_kinetic_energy

   pure function torque(pos, force, origin) result(tau)
     real(dp), intent(in) :: pos(:,:), force(:,:)
     real(dp), intent(in), optional :: origin(3)
     real(dp) :: tau(3)
 
     integer i, N
 
     N = size(pos,2)
     tau = 0.0_dp
     do i=1, N
       if (present(origin)) then
	 tau = tau + ((pos(:,i)-origin).cross.force(:,i))
       else
	 tau = tau + (pos(:,i).cross.force(:,i))
       endif
     end do
   end function torque


   !% Return the temperature, assuming each degree of freedom contributes
   !% $\frac{1}{2}kT$. By default only moving and thermostatted atoms are
   !% included --- this can be overriden by setting 'include_all' to true.
   function temperature(this, region, include_all, instantaneous)

      type(DynamicalSystem), intent(in) :: this
      integer, intent(in), optional  :: region
      logical, intent(in), optional  :: include_all
      logical, intent(in), optional  :: instantaneous
      real(dp)                          :: temperature

      logical ::  my_include_all, my_instantaneous
      integer                           :: i, N
      real(dp)                          :: Ndof

      my_instantaneous = optional_default(.false., instantaneous)
      my_include_all = optional_default(.false., include_all)

      if (my_instantaneous) then
	temperature = 0.0_dp
	N = 0
	Ndof = 0.0_dp
	do i = 1,this%N
	   if (present(region)) then
	      if (this%atoms%thermostat_region(i) == region .and. this%atoms%move_mask(i) == 1) then
		 temperature = temperature + this%atoms%mass(i) * norm2(this%atoms%velo(:,i))
		 N = N + 1
		 Ndof = Ndof + degrees_of_freedom(this,i)
	      end if
	   else
	      if (my_include_all .or. (this%atoms%move_mask(i) == 1)) then
		 temperature = temperature + this%atoms%mass(i) * norm2(this%atoms%velo(:,i))
		 N = N + 1
		 Ndof = Ndof + degrees_of_freedom(this,i)
	      end if
	   end if
	end do

	if (N /= 0) temperature = temperature / ( Ndof * BOLTZMANN_K )         
      else
	temperature = this%cur_temp
      endif

    end function temperature

   !% Add or remove heat from the system by scaling the atomic velocities
   subroutine add_heat(this, heat, Ekin)

      type(DynamicalSystem), intent(inout) :: this
      real(dp),              intent(in)    :: heat, Ekin
      real(dp)                             :: r

      if ( Ekin + heat < 0.0_dp ) call System_Abort('AddHeat: Tried to remove too much heat')

      r = sqrt((Ekin+heat)/Ekin)

      this%atoms%velo = this%atoms%velo * r

   end subroutine add_heat

   !% Rescale the atomic velocities to temperature 'temp'. If the current
   !% temperature is zero, we first randomise the velocites.
   subroutine rescale_velo(this, temp, mass_weighted, zero_L)
      type(DynamicalSystem), intent(inout) :: this
      real(dp),              intent(in)    :: temp
      logical, intent(in), optional        :: mass_weighted, zero_L

      logical                              ::  my_mass_weighted, my_zero_L
      real(dp)                             :: r, currTemp

      my_mass_weighted = optional_default(.false., mass_weighted)
      my_zero_L = optional_default(.false., zero_L)
      currTemp = temperature(this, instantaneous=.true.)
      write(line, '(a,f8.1,a,f8.1,a)')"Rescaling velocities from ",currTemp," K to ",temp," K"; call print(line)

      if ( currTemp .feq. 0.0_dp ) then
         call print('Randomizing velocities')
	 if (my_mass_weighted) then
	   call randomise(this%atoms%velo, 1.0_dp/sqrt(this%atoms%mass))
	 else
	   call randomise(this%atoms%velo, 1.0_dp)
	 endif
         call print('Zeroing total momentum')
         call zero_momentum(this)
	 if (my_zero_L) then
	   call print('Zeroing angular momentum')
	   call zero_angular_momentum(this%atoms)
	 endif
         currTemp = temperature(this, instantaneous=.true.)
      end if

      r = sqrt(temp/currTemp)
      this%atoms%velo = this%atoms%velo * r
   end subroutine rescale_velo

   !% Reinitialise the atomic velocities to temperature 'temp'.
   !% We first randomise the velocites.
   subroutine reinitialise_velo_normal(this, temp, mass_weighted, zero_L)
      type(DynamicalSystem), intent(inout) :: this
      real(dp),              intent(in)    :: temp
      logical, intent(in), optional        :: mass_weighted, zero_L

      logical                              ::  my_mass_weighted, my_zero_L
      real(dp)                             :: r, currTemp
      integer                              :: i

      my_mass_weighted = optional_default(.false., mass_weighted)
      my_zero_L = optional_default(.false., zero_L)
      currTemp = temperature(this, instantaneous=.true.)
      ! write(line, '(a,f8.1,a,f8.1,a)')"Reinitialising velocities from ",currTemp," K to ",temp," K"; call print(line)

      !call print('Randomizing velocities')
      if (my_mass_weighted) then
        do i=1,this%atoms%N
           this%atoms%velo(1:3,i) = ran_normal3()/sqrt(this%atoms%mass(i))
        enddo
      else
        do i=1,this%atoms%N
           this%atoms%velo(1:3,i) = ran_normal3()
        enddo
      endif
      !call print('Zeroing total momentum')
      call zero_momentum(this)
      if (my_zero_L) then
        call print('Zeroing angular momentum')
        call zero_angular_momentum(this%atoms)
      endif
      currTemp = temperature(this, instantaneous=.true.)

      r = sqrt(temp/currTemp)
      this%atoms%velo = this%atoms%velo * r
   end subroutine reinitialise_velo_normal

   !% Draw a velocity component from the correct Gaussian distribution for
   !% a degree of freedom with (effective) mass 'm' at temperature 'T'
   function gaussian_velocity_component(m,T) result(v)

     real(dp), intent(in) :: m, T
     real(dp)             :: v

     v = sqrt(BOLTZMANN_K * T / m) * ran_normal()

   end function gaussian_velocity_component

   !% Construct a velocity vector for a particle of mass 'm'
   !% at temperature 'T' from Gaussian distributed components
   function gaussian_velocity(m,T) result(v)

     real(dp), intent(in)   :: m,T
     real(dp), dimension(3) :: v

     v(1) = gaussian_velocity_component(m,T)
     v(2) = gaussian_velocity_component(m,T)
     v(3) = gaussian_velocity_component(m,T)

   end function gaussian_velocity

   !% Change velocities to those that the system would have in the zero momentum frame.
   !% Optionalally zero the total momentum of a subset of atoms, specified by 'indices'.
   subroutine zero_momentum(this,indices)

     type(DynamicalSystem),           intent(inout) :: this
     integer, optional, dimension(:), intent(in)    :: indices
     real(dp)                                       :: p(3)
     integer                                        :: i

     if (present(indices)) then        
        p = momentum(this,indices)/size(indices)
        forall(i=1:size(indices)) this%atoms%velo(:,indices(i)) = &
                                  this%atoms%velo(:,indices(i)) - p/this%atoms%mass(indices(i))
     else
        p = momentum(this)/real(this%N,dp)
        forall(i=1:this%N) this%atoms%velo(:,i) = this%atoms%velo(:,i)-p/this%atoms%mass(i)
     end if

   end subroutine zero_momentum

   !% give the system a rigid body rotation so as to zero the angular momentum about the centre of mass
   subroutine zero_angular_momentum(this)
     type(Atoms) :: this
 
     real(dp) :: CoM(3), L(3), MoI(3,3), MoI_inv(3,3), angular_vel(3)
     real(dp) :: axis(3), angular_vel_mag
     real(dp) :: dr(3), dr_parallel(3), dr_normal(3), v(3)
 
     integer i
 
     CoM = centre_of_mass(this)
     L = angular_momentum(this,CoM)
     MoI = moment_of_inertia_tensor(this,CoM)
     call inverse(MoI, MoI_inv)
     angular_vel = matmul(MoI_inv,L)
 
     angular_vel_mag = norm(angular_vel)
     axis = angular_vel/angular_vel_mag
 
     if (angular_vel_mag .fne. 0.0_dp) then
       do i=1, this%N
	 dr = this%pos(:,i)-CoM
	 dr_parallel = axis*(dr.dot.axis)
	 dr_normal = dr - dr_parallel
	 v = angular_vel_mag*(axis .cross. dr_normal)
	 this%velo(:,i) = this%velo(:,i) - v
       end do
     end if
   end subroutine zero_angular_momentum

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !X
   !X Centre of Mass calculations
   !X
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   !% Calculates the velocity of the centre of mass c.f. 'centre_of_mass' in Atoms module. No origin atom required.

   pure function centre_of_mass_velo(this,index_list) result(V_CoM)
     
     type(DynamicalSystem),            intent(in) :: this
     integer,  dimension(:), optional, intent(in) :: index_list
     real(dp), dimension(3)                       :: V_CoM
     !local variables
     integer                                      :: i
     real(dp)                                     :: M_tot
     
     V_CoM = 0.0_dp
     M_tot = 0.0_dp
     if (present(index_list)) then
        do i = 1, size(index_list)
           V_CoM = V_CoM + this%atoms%velo(:,index_list(i)) * this%atoms%mass(index_list(i))
           M_tot = M_tot + this%atoms%mass(index_list(i))
        end do
     else
        do i = 1, this%N
           V_CoM = V_CoM + this%atoms%velo(:,i) * this%atoms%mass(i)
           M_tot = M_tot + this%atoms%mass(i)
        end do
     end if

     V_CoM = V_CoM / M_tot

   end function centre_of_mass_velo

   !% Calculates the acceleration of the centre of mass c.f. 'centre_of_mass' in Atoms module. No origin atom required

   pure function centre_of_mass_acc(this,index_list) result(A_CoM)
     
     type(DynamicalSystem),            intent(in) :: this
     integer,  dimension(:), optional, intent(in) :: index_list
     real(dp), dimension(3)                       :: A_CoM
     !local variables
     integer                                      :: i
     real(dp)                                     :: M_tot
     
     A_CoM = 0.0_dp
     M_tot = 0.0_dp

     if (present(index_list)) then
        do i = 1, size(index_list)
           A_CoM = A_CoM + this%atoms%acc(:,index_list(i)) * this%atoms%mass(index_list(i))
           M_tot = M_tot + this%atoms%mass(index_list(i))
        end do
     else
        do i = 1, this%N
           A_CoM = A_CoM + this%atoms%acc(:,i) * this%atoms%mass(i)
           M_tot = M_tot + this%atoms%mass(i)
        end do
     end if

     A_CoM = A_CoM / M_tot

   end function centre_of_mass_acc

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !XXX
   !XXX  ADVANCE VERLET 1
   !XXX
   !XXX  This is the first step of the integration algorithm.
   !XXX  
   !XXX  On entry we have r(t), v(t) and a(t). 
   !XXX  
   !XXX  For normal atoms this routine performs the following steps of velocity
   !XXX  Verlet:
   !XXX
   !XXX  p(t+dt/2) = p(t) + F(t)dt/2          ->  v(t+dt/2) = v(t) + a(t)dt/2
   !XXX
   !XXX  r(t+dt)   = r(t) + (1/m)p(t+dt/2)dt  ->  r(t+dt)   = r(t) + v(t+dt/2)dt
   !XXX
   !XXX  Integration algorithms for other types of atoms fit around this. After this
   !XXX  routine, the user code must calculate F(t+dt) and call advance_verlet2 to
   !XXX  complete the integration step.
   !XXX  
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
   subroutine advance_verlet1(this,dt,f,virial,parallel,store_constraint_force)

     type(dynamicalsystem), intent(inout)  :: this
     real(dp),              intent(in)     :: dt
     real(dp),              intent(in)     :: f(:,:)
     real(dp),dimension(3,3), intent(in), optional :: virial
     logical, optional,     intent(in)     :: parallel
     logical, optional,     intent(in)     :: store_constraint_force

     logical                               :: do_parallel, do_store
     integer                               :: i, j, g, n, ntherm
     real(dp), allocatable                 :: therm_ndof(:)

#ifdef _MPI
     include 'mpif.h'
     real(dp), dimension(:,:), allocatable :: mpi_pos, mpi_velo, mpi_acc, mpi_constraint_force
     integer                               :: error_code, pos_indices(3), cf_indices(2)
#endif
     
     do_parallel = optional_default(.false.,parallel)
     do_store = optional_default(.false.,store_constraint_force)
     call check_size('Force',f,(/3,this%N/),'advance_verlet1')

     this%dW = 0.0_dp
     ntherm = size(this%thermostat)-1
     allocate(therm_ndof(ntherm))
    
#ifdef _MPI
     if (do_parallel) then
        allocate(mpi_pos(3,this%N), mpi_velo(3,this%N), mpi_acc(3,this%N))
        mpi_pos = 0.0_dp
        mpi_velo = 0.0_dp
        mpi_acc = 0.0_dp
        if (do_store) then
           allocate(mpi_constraint_force(3,this%N))
           mpi_constraint_force = 0.0_dp
           if (get_value(this%atoms%properties,'constraint_force',pos_indices)) then
              cf_indices = pos_indices(2:3)
           else
              call system_abort('advance_verlet1: no constraint_force property found')
           end if
        end if
     end if
#endif

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X CALCULATE AND SET DEGREES OF FREEDOM FOR EACH THERMOSTAT
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     therm_ndof = 0.0_dp
     do i = 1, this%N
        j = this%atoms%thermostat_region(i)
        if (j>0 .and. j<=ntherm) therm_ndof(j) = therm_ndof(j) + degrees_of_freedom(this,i)
     end do
     do i = 1, ntherm
        call set_degrees_of_freedom(this%thermostat(i),therm_ndof(i))
     end do
     deallocate(therm_ndof)

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X DAMPING
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     if (this%thermostat(0)%type==LANGEVIN) then
        call thermostat1(this%thermostat(0),this%atoms,f,dt,'damp_mask',1)
     end if

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X THERMOSTATTING
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     do i = 1, ntherm
        call thermostat1(this%thermostat(i),this%atoms,f,dt,'thermostat_region',i,virial)
     end do

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X VELOCITY UPDATE
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     do g = 1, size(this%group)

#ifdef _MPI
        if (do_parallel .and. mod(g,mpi_n_procs()) /= mpi_id()) cycle
#endif
        
        select case(this%group(g)%type)

        case(TYPE_ATOM, TYPE_CONSTRAINED) ! Constrained atoms undergo the usual verlet step here
           
           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
           !X
           !X v(t+dt/2) = v(t) + a(t)dt/2
           !X
           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
           do n = 1, group_n_atoms(this%group(g))
              i = group_nth_atom(this%group(g),n)
#ifdef _MPI
              if (do_parallel) then
                 if (this%atoms%move_mask(i)==0) then
                    mpi_velo(:,i) = 0.0_dp
                 else
                    mpi_velo(:,i) = this%atoms%velo(:,i) + 0.5_dp*this%atoms%acc(:,i)*dt
                 end if
              else
                 if (this%atoms%move_mask(i)==0) then
                    this%atoms%velo(:,i) = 0.0_dp
                 else
                    this%atoms%velo(:,i) = this%atoms%velo(:,i) + 0.5_dp*this%atoms%acc(:,i)*dt
                 end if
              end if
#else
              if (this%atoms%move_mask(i)==0) then
                 this%atoms%velo(:,i) = 0.0_dp
              else
                 this%atoms%velo(:,i) = this%atoms%velo(:,i) + 0.5_dp*this%atoms%acc(:,i)*dt
              end if
#endif              
           end do
           
        end select

     end do

#ifdef _MPI
     ! Broadcast the new velocities
     if (do_parallel) then
        call MPI_ALLREDUCE(mpi_velo,this%atoms%velo,size(mpi_velo),MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,error_code)
        call abort_on_mpi_error(error_code,'advance_verlet1 - velocity update')       
     end if
#endif

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X DAMPING
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     if (this%thermostat(0)%type==LANGEVIN) then
        call thermostat2(this%thermostat(0),this%atoms,f,dt,'damp_mask',1)
     end if

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X THERMOSTATTING
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     do i=1, ntherm
        call thermostat2(this%thermostat(i),this%atoms,f,dt,'thermostat_region',i)
     end do
 
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X POSITION UPDATE
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     do g = 1, size(this%group)

#ifdef _MPI
        if (do_parallel .and. mod(g,mpi_n_procs()) /= mpi_id()) cycle
#endif

        select case(this%group(g)%type)

        case(TYPE_ATOM)
           
           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
           !X
           !X r(t+dt) = r(t) + v(t+dt/2)dt
           !X
           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
           do n = 1, group_n_atoms(this%group(g))
              i = group_nth_atom(this%group(g),n)
#ifdef _MPI
              if (do_parallel) then
                 if (this%atoms%move_mask(i)==0) then
                    mpi_pos(:,i) =  this%atoms%pos(:,i)
                    mpi_velo(:,i) = 0.0_dp
                    mpi_acc(:,i) = 0.0_dp
                 else
                    mpi_pos(:,i) =  this%atoms%pos(:,i) + this%atoms%velo(:,i)*dt
                    mpi_velo(:,i) = this%atoms%velo(:,i)
                    mpi_acc(:,i) = this%atoms%acc(:,i)
                 end if
              else
                 if (this%atoms%move_mask(i)==0) then
                    this%atoms%velo(:,i) = 0.0_dp
                    this%atoms%acc(:,i) = 0.0_dp
                 else
                    this%atoms%pos(:,i) = this%atoms%pos(:,i) + this%atoms%velo(:,i)*dt
                 end if
              end if
#else
              if (this%atoms%move_mask(i)==0) then
                 this%atoms%velo(:,i) = 0.0_dp
                 this%atoms%acc(:,i) = 0.0_dp
              else
                 this%atoms%pos(:,i) = this%atoms%pos(:,i) + this%atoms%velo(:,i)*dt
              end if
#endif                 
           end do
           
        case(TYPE_CONSTRAINED)

           !We don't test move_mask here: constrained atoms are fixed by extra constraints.
           call shake(this%atoms,this%group(g),this%constraint,this%t,dt,store_constraint_force)
#ifdef _MPI
           if (do_parallel) then
              do n = 1, group_n_atoms(this%group(g))
                 i = group_nth_atom(this%group(g),n)
                 mpi_pos(:,i) =  this%atoms%pos(:,i)
                 mpi_velo(:,i) = this%atoms%velo(:,i)
                 mpi_acc(:,i) = this%atoms%acc(:,i)
                 if (do_store) mpi_constraint_force(:,i) = this%atoms%data%real(cf_indices,i)
              end do
           end if
#endif                 

        end select

     end do

#ifdef _MPI
     ! Broadcast the new positions, velocities, accelerations and possibly constraint forces
     if (do_parallel) then
        call MPI_ALLREDUCE(mpi_pos,this%atoms%pos,size(mpi_pos),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error_code)
        call abort_on_mpi_error(error_code,'advance_verlet1 - position update')       
        call MPI_ALLREDUCE(mpi_velo,this%atoms%velo,size(mpi_velo),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error_code)
        call abort_on_mpi_error(error_code,'advance_verlet1 - velocity update 2')       
        call MPI_ALLREDUCE(mpi_acc,this%atoms%acc,size(mpi_acc),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error_code)
        call abort_on_mpi_error(error_code,'advance_verlet1 - acceleration update')       
        if (do_store) then
           call MPI_ALLREDUCE(mpi_constraint_force,this%atoms%data%real(cf_indices,1:this%N),size(mpi_constraint_force),&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error_code)
           call abort_on_mpi_error(error_code,'advance_verlet1 - constraint force update')
        end if
     end if
#endif


     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X BOOKKEEPING
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     
#ifdef _MPI
     if (do_parallel) then
        deallocate(mpi_pos, mpi_velo, mpi_acc)
        if (do_store) deallocate(mpi_constraint_force)
     end if
#endif

     call calc_dists(this%atoms,parallel=do_parallel)

   end subroutine advance_verlet1

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !XXX
   !XXX  ADVANCE VERLET 2
   !XXX
   !XXX  This is the second step of the integration algorithm.
   !XXX  
   !XXX  On entry we have r(t+dt), v(t+dt/2) and a(t). F(t+dt) is calculated by the
   !XXX  user code between calls to advance_verlet1 and advance_verlet2, and passed
   !XXX  in as the argument 'f'
   !XXX  
   !XXX  For normal atoms this routine performs the last step of velocity Verlet:
   !XXX
   !XXX  p(t+dt) = p(t+dt/2) + F(t+dt)dt/2  ->  v(t+dt) = v(t+dt/2) + a(t+dt)dt/2
   !XXX
   !XXX  where a(t+dt) = (1/m) F(t+dt)
   !XXX
   !XXX  Integration algorithms for other types of atoms fit around this.
   !XXX  
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine advance_verlet2(this,dt,f,virial,parallel,store_constraint_force)

     type(dynamicalsystem), intent(inout)  :: this
     real(dp),              intent(in)     :: dt
     real(dp),              intent(in)     :: f(:,:)
     real(dp),dimension(3,3), intent(in), optional :: virial
     logical, optional,     intent(in)     :: parallel
     logical, optional,     intent(in)     :: store_constraint_force

     logical                               :: do_parallel, do_store
     integer                               :: i, g, n, ntherm
     real(dp)                              :: decay
     
#ifdef _MPI
     include 'mpif.h'
     real(dp), dimension(:,:), allocatable :: mpi_velo, mpi_acc, mpi_constraint_force
     integer                               :: error_code, pos_indices(3), cf_indices(2)
#endif

     do_parallel = optional_default(.false.,parallel)
     do_store = optional_default(.false.,store_constraint_force)
     ntherm = size(this%thermostat)-1
     call check_size('Force',f,(/3,this%N/),'advance_verlet2')

#ifdef _MPI
     if (do_parallel) then
        allocate(mpi_velo(3,this%N), mpi_acc(3,this%N))
        mpi_velo = 0.0_dp
        mpi_acc = 0.0_dp
        if (do_store) then
           allocate(mpi_constraint_force(3,this%N))
           mpi_constraint_force = 0.0_dp
           if (get_value(this%atoms%properties,'constraint_force',pos_indices)) then
              cf_indices = pos_indices(2:3)
           else
              call system_abort('advance_verlet1: no constraint_force property found')
           end if
        end if
     end if
#endif
    
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X CONVERT FORCES TO ACCELERATIONS
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     forall (i=1:this%N) this%atoms%acc(:,i) = f(:,i) / this%atoms%mass(i)

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X DAMPING
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     if (this%thermostat(0)%type==LANGEVIN) then
        call thermostat3(this%thermostat(0),this%atoms,f,dt,'damp_mask',1)
     end if

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X THERMOSTATTING
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     do i = 1, ntherm
        call thermostat3(this%thermostat(i),this%atoms,f,dt,'thermostat_region',i)
     end do

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X VELOCITY UPDATE
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     do g = 1, size(this%group)

#ifdef _MPI
        if (do_parallel .and. mod(g,mpi_n_procs()) /= mpi_id()) cycle
#endif

        select case(this%group(g)%type)

        case(TYPE_ATOM)
           
           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
           !X
           !X v(t+dt) = v(t+dt/2) + a(t+dt)dt/2
           !X
           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
           do n = 1, group_n_atoms(this%group(g))
              i = group_nth_atom(this%group(g),n)
#ifdef _MPI
              if (do_parallel) then
                 if (this%atoms%move_mask(i)==0) then
                    mpi_velo(:,i) = 0.0_dp
                    mpi_acc(:,i) = 0.0_dp
                 else
                    mpi_velo(:,i) = this%atoms%velo(:,i) + 0.5_dp*this%atoms%acc(:,i)*dt
                    mpi_acc(:,i) = this%atoms%acc(:,i)
                 end if
              else
                 if (this%atoms%move_mask(i)==0) then
                    this%atoms%velo(:,i) = 0.0_dp
                    this%atoms%acc(:,i) = 0.0_dp
                 else
                    this%atoms%velo(:,i) = this%atoms%velo(:,i) + 0.5_dp*this%atoms%acc(:,i)*dt
                 end if
              end if
#else
              if (this%atoms%move_mask(i)==0) then
                 this%atoms%velo(:,i) = 0.0_dp
                 this%atoms%acc(:,i) = 0.0_dp
              else
                 this%atoms%velo(:,i) = this%atoms%velo(:,i) + 0.5_dp*this%atoms%acc(:,i)*dt
              end if
#endif              
           end do
           
        case(TYPE_CONSTRAINED)

           !As with shake, we don't test move_mask here
           call rattle(this%atoms,this%group(g),this%constraint,this%t,dt,store_constraint_force)
#ifdef _MPI
           if (do_parallel) then
              do n = 1, group_n_atoms(this%group(g))
                 i = group_nth_atom(this%group(g),n)
                 mpi_velo(:,i) = this%atoms%velo(:,i)
                 mpi_acc(:,i) = this%atoms%acc(:,i)
                 if (do_store) mpi_constraint_force(:,i) = this%atoms%data%real(cf_indices,i)
              end do
           end if
#endif

        end select

     end do

#ifdef _MPI
     ! Broadcast the new velocities, accelerations and possibly constraint forces
     if (do_parallel) then
        call MPI_ALLREDUCE(mpi_velo,this%atoms%velo,size(mpi_velo),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error_code)
        call abort_on_mpi_error(error_code,'advance_verlet1 - velocity update 2')       
        call MPI_ALLREDUCE(mpi_acc,this%atoms%acc,size(mpi_acc),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error_code)
        call abort_on_mpi_error(error_code,'advance_verlet1 - acceleration update')       
        if (do_store) then
           call MPI_ALLREDUCE(mpi_constraint_force,this%atoms%data%real(cf_indices,1:this%N),size(mpi_constraint_force),&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error_code)
           call abort_on_mpi_error(error_code,'advance_verlet1 - constraint force update')
        end if
     end if
#endif

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X DAMPING
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     if (this%thermostat(0)%type==LANGEVIN) then
        call thermostat4(this%thermostat(0),this%atoms,f,dt,'damp_mask',1)
     end if

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X THERMOSTATTING
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     do i = 1, ntherm
        call thermostat4(this%thermostat(i),this%atoms,f,dt,'thermostat_region',i,virial)
     end do

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X BOOKKEEPING
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     this%t = this%t + dt
     this%nSteps = this%nSteps + 1
     this%work = this%work + this%dW ! Not currently calculated

     this%cur_temp = temperature(this, instantaneous=.true.)

     if (this%avg_time > 0.0_dp) then
        decay = dt/this%avg_time
        call update_exponential_average(this%avg_temp,decay,this%cur_temp)
        do i = 1, this%N
           call update_exponential_average(this%atoms%avgpos(:,i),decay,realpos(this%atoms,i))
           call update_exponential_average(this%atoms%avg_ke(i),decay,0.5_dp*this%atoms%mass(i)*norm2(this%atoms%velo(:,i)))
        end do
     end if

#ifdef _MPI
     if (do_parallel) then
        deallocate(mpi_velo, mpi_acc)
        if (do_store) deallocate(mpi_constraint_force)
     end if
#endif

   end subroutine advance_verlet2

   !% Calls advance_verlet2 followed by advance_verlet1. Outside this routine the
   !% velocities will be half-stepped.
   subroutine advance_verlet(ds,dt,f,virial,parallel,store_constraint_force)

     type(dynamicalsystem), intent(inout) :: ds
     real(dp),              intent(in)    :: dt
     real(dp),              intent(in)    :: f(:,:)
     real(dp),dimension(3,3), intent(in), optional :: virial
     logical, optional,     intent(in)    :: parallel
     logical, optional,     intent(in)    :: store_constraint_force
     logical, save                        :: first_call = .true.

     if (first_call) then
        call print_title('SINGLE STEP VERLET IN USE')
        call print('Consider changing to the two-step integrator')
        call print_title('=')
        first_call = .false.
     end if
     call advance_verlet2(ds,dt,f,virial,parallel,store_constraint_force)
     call advance_verlet1(ds,dt,f,virial,parallel,store_constraint_force)

   end subroutine advance_verlet


!   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   !XXX
!   !XXX  ADVANCE VERLET
!   !XXX
!   !XXX  On entry we have r(t), v(t-dt) and a(t-dt)
!   !XXX  'force' contains f(t) = f[r(t)] = m a(t)
!   !XXX
!   !XXX  The algorithm is:
!   !XXX
!   !XXX  v(t-dt/2) = v(t-dt) + 1/2 a(t-dt) dt
!   !XXX
!   !XXX  v(t)      = v(t-dt/2) + 1/2 a(t) dt
!   !XXX
!   !XXX  r(t+dt)   = r(t) + v(t) dt + 1/2 a(t) dt^2
!   !XXX
!   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   
!   
!   ! Advance the positions, accelerations and velocites by one timestep 'dt'.
!   ! On entry we have $\mathbf{r}(t)$, $\mathbf{v}(t-\mathrm{d}t)$ and $\mathbf{a}(t-\mathrm{d}t)$
!   ! 'force' contains $\mathbf{f}(t) = \mathbf{f}(\mathbf{r}(t)) = m \mathbf{a}(t)$. The
!   ! algorithm used is:
!   ! \begin{displaymath}
!   ! \begin{array}{l}
!   ! \mathbf{v}\left(t-\frac{\mathrm{d}t}{2}\right) = 
!   !    \mathbf{v}\left(t-\mathrm{d}t\right) + \frac{1}{2} \mathbf{a}\left(t-\mathrm{d}t\right) \mathrm{d}t \\
!   ! \mathbf{v}\left(t\right)      = 
!   !   \mathbf{v}\left(t-\frac{\mathrm{d}t}{2}\right) + \frac{1}{2} \mathbf{a}\left(t\right) \mathrm{d}t \\
!   ! \mathbf{r}\left(t+\mathrm{d}t\right)   = 
!   !    \mathbf{r}\left(t\right) + \mathbf{v}\left(t\right) \mathrm{d}t + 
!   !    \frac{1}{2} \mathbf{a}\left(t\right) \mathrm{d}t^2
!   ! \end{array}
!   ! \end{displaymath}
!
!   subroutine advance_verlet(this,dt,force,dV_dt,parallel,constraint_force)
!
!     !ARGUMENTS  
!     type(DynamicalSystem), target, intent(inout)         :: this
!     real(dp),                      intent(in)            :: dt
!     real(dp), dimension(:,:),      intent(in)            :: force
!     real(dp), optional, intent(in)                       :: dV_dt
!     logical,  optional, intent(in)                       :: parallel
!     real(dp), dimension(:,:), optional, intent(out)      :: constraint_force
!     !GENERAL
!     integer                                              :: i, j, g, n, nn, o, type, pos_indices(3)
!     logical                                              :: do_parallel
!     !THERMOSTATTING, AVERAGING
!     real(dp)                                             :: currTemp, av_fact1, av_fact2
!     real(dp), dimension(3)                               :: tmp_realpos
!     real(dp), dimension(:,:), allocatable                :: thermostat_force
!     !CONSTRAINTS
!     real(dp), dimension(:), allocatable                  :: pos, velo, oldacc, oldpos, oldvelo
!     real(dp)                                             :: dlambda, m
!     logical                                              :: converged, store_constraint_forces
!     integer                                              :: iterations
!     !RIGID BODIES
!!     real(dp), dimension(3)                               :: ACoM, tau      
!!     real(dp), dimension(4)                               :: tau4, F4
!
!     ! Extended energy
!     real(dp) :: my_dV_dt
!     real(dp), save :: dV_dt_old = 0.0_dp
!
!#ifdef _MPI
!     include 'mpif.h'
!
!     real(dp), dimension(:,:), allocatable :: mpi_pos, mpi_velo, mpi_acc, mpi_constraint_force
!     integer                               :: error_code
!#endif
!
!     do_parallel = optional_default(.false.,parallel)
!     store_constraint_forces = present(constraint_force)
!     
!     if (store_constraint_forces) constraint_force = 0.0_dp
!
!#ifdef _MPI
!
!     if (do_parallel) then
!        allocate(mpi_pos(3,this%N), mpi_velo(3,this%N), mpi_acc(3,this%N))
!        mpi_pos = 0.0_dp
!        mpi_velo = 0.0_dp
!        mpi_acc = 0.0_dp
!        if (store_constraint_forces) then
!           allocate(mpi_constraint_force(3,this%N))
!           mpi_constraint_force = 0.0_dp
!        end if
!     end if
!
!#endif
!
!     !Complain if the size of force is wrong
!     call check_size('Force',force,(/3,this%N/),'AdvanceVerlet')
!
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     !X Check the temperature
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!     currTemp = temperature(this, instantaneous=.true.)
!     if ((this%old_thermostat /= NO_THERM) .and. (abs(currTemp - this%sim_temp) > 0.5*this%sim_temp)) then
!        write(line,'(3(a,f0.3))') 'advance_verlet: Temperature of tempered zone (',currTemp,  &
!             ') is very different from target (',this%sim_temp,') at time t = ',this%t
!        call print_warning(line)
!     end if
!     
!     !store the sum of the forces
!     this%forcesum = sum(force,dim = 2)
!     
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     !X
!     !X Advance the velocities by 0.5*dt: v(t-dt/2) = v(t-dt) + 1/2 a(t-dt) dt
!     !X
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!     do g = 1, size(this%group)
!
!#ifdef _MPI
!        if (do_parallel .and. mod(g,mpi_n_procs()) /= mpi_id()) cycle
!#endif
!        
!        type = this%group(g)%type
!        
!        select case(type)
!           
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!           !X NORMAL ATOMS
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX            
!        case(TYPE_ATOM)
!           do n = 1, Group_N_Atoms(this%group(g))  !Loop over all atoms in this group
!              i = Group_Nth_Atom(this%group(g),n)
!#ifdef _MPI
!              if (do_parallel) then
!                 mpi_velo(:,i) = this%atoms%velo(:,i) + 0.5_dp * dt * this%atoms%acc(:,i)
!              else
!                 this%atoms%velo(:,i) = this%atoms%velo(:,i) + 0.5_dp * dt * this%atoms%acc(:,i)
!              end if
!#else
!              this%atoms%velo(:,i) = this%atoms%velo(:,i) + 0.5_dp * dt * this%atoms%acc(:,i)
!#endif
!           end do
!
!
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!           !X CONSTRAINED ATOMS
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!        case(TYPE_CONSTRAINED)
!           !Integrate each atom as usual.
!           do n =1, Group_N_Atoms(this%group(g))
!              i = Group_Nth_Atom(this%group(g),n)
!#ifdef _MPI
!              if (do_parallel) then
!                 mpi_velo(:,i) = this%atoms%velo(:,i) + 0.5_dp * dt * this%atoms%acc(:,i)
!              else
!                 this%atoms%velo(:,i) = this%atoms%velo(:,i) + 0.5_dp * dt * this%atoms%acc(:,i)
!              end if
!#else
!              this%atoms%velo(:,i) = this%atoms%velo(:,i) + 0.5_dp * dt * this%atoms%acc(:,i)
!#endif
!           end do
!
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!           !X RIGID ATOMS
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX           
!        case(TYPE_RIGID)
!
!           cycle
!
!!           j = Group_Nth_Object(this%group(g),1)    !This is the rigid body index
!!           call RigidBody_CalculateVariables(this,j)
!!
!!           !Integrate the centre of mass variables
!!           ACoM = Centre_Of_Mass_Acceleration(this,this%rigidbody(j)%indices)
!!           this%rigidbody(j)%VCoM = this%rigidbody(j)%VCoM + 0.5_dp * ACoM * dt
!!           
!!           !Step the quaternion momentum by 0.5dt:
!!           tau = Rigidbody_Torque(this,j)
!!           tau4 = (/0.0_dp,tau/)
!!           F4 = 2.0_dp * (no_squish_S(this%rigidbody(j)%q) .mult. tau4)
!!           this%rigidbody(j)%p = this%rigidbody(j)%p + 0.5_dp * F4 * dt
!!           !Do a free rotation
!!           call no_squish_Free_Rotor(this%rigidbody(j)%q, &
!!                this%rigidbody(j)%p, &
!!                this%rigidbody(j)%model%I, &
!!                0.5_dp * dt)
!!           !Write back the positions and velocities
!!           call RigidBody_WritePositionsVelocities(this,j)
!!
!        end select
!
!     end do
!
!#ifdef _MPI
!
!     ! Update the velocities on all processes
!     if (do_parallel) then
!
!        call MPI_ALLREDUCE(mpi_velo,this%atoms%velo,size(mpi_velo),MPI_DOUBLE_PRECISION,&
!                           MPI_SUM,MPI_COMM_WORLD,error_code)
!        call abort_on_mpi_error(error_code,'Advance_verlet - pre thermostat1')
!        
!     end if
!
!#endif
!
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     !X Get new accelerations 
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!     forall(i = 1:this%N) this%atoms%acc(:,i) = force(:,i) / this%atoms%mass(i)
!
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     !X Do first part of thermostatting
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!        
!     allocate(thermostat_force(3,this%N))
!     thermostat_force = 0.0_dp
!
!#ifdef _MPI
!     ! Update the velocities on all processes
!     if (do_parallel) then
!
!        call MPI_ALLREDUCE(mpi_velo,this%atoms%velo,size(mpi_velo),MPI_DOUBLE_PRECISION,&
!                           MPI_SUM,MPI_COMM_WORLD,error_code)
!        call abort_on_mpi_error(error_code,'Advance_verlet - pre thermostat1')
!        
!     end if
!#endif
!
!     call old_thermostat1(this, dt, thermostat_force)
!
!     !Store thermostat force
!     if (get_value(this%atoms%properties,'thermo_force',pos_indices)) then
!        this%atoms%data%real(pos_indices(2):pos_indices(3),:) = thermostat_force
!     end if
!
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     !X
!     !X Advance the velocities again by 0.5*dt: v(t) = v(t-dt/2) + 1/2 a(t) dt
!     !X
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!     do g = 1, size(this%group)
!
!#ifdef _MPI
!        if (do_parallel .and. mod(g,mpi_n_procs()) /= mpi_id()) cycle
!#endif
!
!        type = this%group(g)%type
!
!        select case(type)
!
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!           !X NORMAL ATOMS
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX           
!        case(TYPE_ATOM)
!           do n = 1, Group_N_Atoms(this%group(g))
!              i = Group_Nth_Atom(this%group(g),n)
!#ifdef _MPI
!              if (do_parallel) then
!                 mpi_velo(:,i) = this%atoms%velo(:,i) + 0.5_dp * dt * this%atoms%acc(:,i)
!                 mpi_acc(:,i) = this%atoms%acc(:,i)
!              else
!                 this%atoms%velo(:,i) = this%atoms%velo(:,i) + 0.5_dp * dt * this%atoms%acc(:,i)
!              end if
!#else
!              this%atoms%velo(:,i) = this%atoms%velo(:,i) + 0.5_dp * dt * this%atoms%acc(:,i)
!#endif
!           end do
!
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!           !X CONSTRAINED ATOMS
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX           
!        case(TYPE_CONSTRAINED)
!           !RATTLE
!           !Store the current accelerations and velocities
!           allocate(oldacc(3*size(this%group(g)%atom)),&
!                    oldvelo(3*size(this%group(g)%atom)))
!           do n = 1, Group_N_Atoms(this%group(g))
!              i = Group_Nth_Atom(this%group(g),n)
!              oldacc(3*n-2:3*n) = this%atoms%acc(:,i)
!              oldvelo(3*n-2:3*n) = this%atoms%velo(:,i)
!           end do
!
!           !Calculate and store the gradients of the constraints at the current positions and
!           !initialise lambdas to zero
!           do n = 1, Group_N_Objects(this%group(g)) !Loop over all constraints in this group
!
!              i = Group_Nth_Object(this%group(g),n)
!              call reallocate(pos,3*this%constraint(i)%N)
!              call reallocate(velo,3*this%constraint(i)%N)
!
!              !Use minimum image convention for constraints with more than one atom
!              o = this%constraint(i)%atom(1) !'origin' atom for this constraint
!              if (this%constraint(i)%N /= 1) then
!                 do nn = 1, this%constraint(i)%N
!                    j = this%constraint(i)%atom(nn)
!                    pos(3*nn-2:3*nn) = diff_min_image(this%atoms,o,j)
!                    velo(3*nn-2:3*nn) = this%atoms%velo(:,j)
!                 end do
!              else
!                 !Otherwise, use the absolute position of that single atom
!                 pos = this%atoms%pos(:,o)
!                 velo = this%atoms%velo(:,o)
!              end if
!
!              call constraint_calculate_values(this%constraint(i),pos,velo,this%t)
!
!              call constraint_store_gradient(this%constraint(i))
!
!!              this%constraint(i)%lambdaV = 0.0_dp
!
!           end do           
!
!           iterations = 0
!           do
!              !Return the velocities and accelerations to their initial values
!              do n = 1, Group_N_Atoms(this%group(g))
!                 i = Group_Nth_Atom(this%group(g),n)
!                 this%atoms%velo(:,i) = oldvelo(3*n-2:3*n)
!                 this%atoms%acc(:,i) = oldacc(3*n-2:3*n)
!              end do
!
!              !Update accelerations via a -> a - lambda * grad(C)_initial / mass
!              do n = 1, Group_N_Objects(this%group(g))
!                 i = Group_Nth_Object(this%group(g),n)
!                 do nn = 1, this%constraint(i)%N
!                    j = this%constraint(i)%atom(nn)
!                    this%atoms%acc(:,j) = this%atoms%acc(:,j) - this%constraint(i)%lambdaV &
!                         * this%constraint(i)%old_dC_dr(3*nn-2:3*nn) &
!                         / this%atoms%mass(j)                   
!                 end do
!              end do
!              
!              !Calculate the new velocities with these accelerations
!              do n = 1, Group_N_Atoms(this%group(g))
!                 i = Group_Nth_Atom(this%group(g),n)
!                 this%atoms%velo(:,i) = this%atoms%velo(:,i) + 0.5_dp * dt * this%atoms%acc(:,i)
!              end do
!
!              !Calculate constraints and test for convergence
!              converged = .true.
!              do n = 1, Group_N_Objects(this%group(g))
!
!                 i = Group_Nth_Object(this%group(g),n)
!                 call reallocate(pos,3*this%constraint(i)%N)
!                 call reallocate(velo,3*this%constraint(i)%N)
!                 o = this%constraint(i)%atom(1)
!                 if (this%constraint(i)%N /= 1) then
!
!                    do nn = 1, this%constraint(i)%N
!                       j = this%constraint(i)%atom(nn)
!                       pos(3*nn-2:3*nn) = diff_min_image(this%atoms,o,j)
!                       velo(3*nn-2:3*nn) = this%atoms%velo(:,j)
!                    end do
!                 else
!                    pos = this%atoms%pos(:,o)
!                    velo = this%atoms%velo(:,o)
!                 end if
!
!                 call constraint_calculate_values(this%constraint(i),pos,velo,this%t) !this%t is the end point
!                                                                                      !of this integration step
!                                                                                      !for the velocities
!
!                 if (abs(this%constraint(i)%dC_dt) > this%constraint(i)%tol) converged = .false.
!
!              end do
!
!              iterations = iterations + 1
!
!              if (converged) exit
!              if (iterations > RATTLE_MAX_ITERATIONS) then
!                 call print_warning('advance_verlet: RATTLE did not converge for this group')
!                 call print(this%group(g),g)
!                 call print('')
!                 call print('Constraints:')
!                 call print('')
!                 do n = 1, Group_N_Objects(this%group(g))
!                    i = Group_Nth_Object(this%group(g),n)
!                    call print(this%constraint(i),i)
!                 end do
!                 call system_abort('AdvanceVerlet: RATTLE reached maximum iterations')
!              end if
!             
!              !Update lambdas
!              do n = 1, group_n_objects(this%group(g)) !Loop over constraints in the group
!                 i = group_nth_object(this%group(g),n)
!                 m = 0.0_dp
!                 do nn = 1, this%constraint(i)%N       !Loop over each particle in the constraint
!                    j = this%constraint(i)%atom(nn)
!                    m = m + norm2(this%constraint(i)%dC_dr(3*nn-2:3*nn)) / this%atoms%mass(j)
!                 end do
!                 dlambda = 2.0_dp * this%constraint(i)%dC_dt / (m * dt)
!                 this%constraint(i)%lambdaV = this%constraint(i)%lambdaV + dlambda
!              end do
!
!           end do
!
!           !Update stored constraint forces
!           if (store_constraint_forces) then
!              do n = 1, Group_N_Objects(this%group(g))
!                 i = Group_Nth_Object(this%group(g),n)
!                 do nn = 1, this%constraint(i)%N
!                    j = this%constraint(i)%atom(nn)
!                    constraint_force(:,j) = constraint_force(:,j) - this%constraint(i)%lambdaV &
!                         * this%constraint(i)%old_dC_dr(3*nn-2:3*nn)
!                 end do
!              end do
!           end if
!
!           deallocate(pos,velo,oldacc,oldvelo)
!
!#ifdef _MPI
!           if (do_parallel) then
!              
!              ! Gather the updated velocities and accelerations into the send buffers
!              do n = 1, Group_N_Atoms(this%group(g))
!                 i = Group_Nth_Atom(this%group(g),n)
!                 mpi_velo(:,i) = this%atoms%velo(:,i)
!                 mpi_acc(:,i) = this%atoms%acc(:,i)
!                 if (store_constraint_forces) mpi_constraint_force(:,i) = constraint_force(:,i)
!              end do
!
!           end if
!#endif           
!           
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!           !X RIGID ATOMS
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX           
!        case(TYPE_RIGID)
!
!!           j = group_nth_object(this%group(g),1)
!!
!!           !VCoM(t) = VCoM(t-dt/2) + 0.5 ACoM(t) dt
!!           ACoM = centre_of_mass_acc(this,this%rigidbody(j)%indices)
!!           this%rigidbody(j)%VCoM = this%rigidbody(j)%VCoM + 0.5_dp * ACoM * dt
!!
!!           !Step the quaternion momentum by 0.5dt:
!!           !P(t) = P(t-dt/2) + 0.5 F(t) dt
!!           tau = rigidbody_torque(this,j)
!!           tau4 = (/0.0_dp,tau/)
!!           F4 = 2.0_dp * (no_squish_S(this%rigidbody(j)%q) .mult. tau4)
!!           this%rigidbody(j)%p = this%rigidbody(j)%p + 0.5_dp * F4 * dt
!!
!!           !NOTE: calling WriteVelocities at this point HERE would give v(t)
!!
!!           !VCoM(t+dt/2) = VCoM(t) + 0.5 ACoM(t) dt
!!           this%rigidbody(j)%VCoM = this%rigidbody(j)%VCoM + 0.5_dp * ACoM * dt
!!           !P'(t+dt/2) = P(t) + 0.5 F(t) dt
!!           this%rigidbody(j)%p = this%rigidbody(j)%p + 0.5_dp * F4 * dt
!           
!        end select
!
!     end do
!
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     !X Do second part of thermostatting
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!#ifdef _MPI
!     if (do_parallel) then
!
!        call MPI_ALLREDUCE(mpi_velo,this%atoms%velo,size(mpi_velo),MPI_DOUBLE_PRECISION,&
!                           MPI_SUM,MPI_COMM_WORLD,error_code)
!        call abort_on_mpi_error(error_code,'Advance_verlet - pre thermostat2 velocities')
!
!        call MPI_ALLREDUCE(mpi_acc,this%atoms%acc,size(mpi_acc),MPI_DOUBLE_PRECISION,&
!                           MPI_SUM,MPI_COMM_WORLD,error_code)
!        call abort_on_mpi_error(error_code,'Advance_verlet - pre thermostat2 accelerations')
!               
!     end if
!#endif
!
!     call old_thermostat2(this, dt)
!
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     !X Zero the velocities and accelerations of the fixed atoms
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!     forall( i = 1:this%N, this%atoms%move_mask(i) == 0)
!        this%atoms%velo(:,i) = 0.0_dp
!        this%atoms%acc(:,i) = 0.0_dp
!     end forall
!
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     !X
!     !X Advance the positions by dt: r(t+dt) = r(t) + v(t) dt + 1/2 a(t) dt^2
!     !X
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!     do g = 1, size(this%group)
!
!#ifdef _MPI
!        if (do_parallel .and. mod(g,mpi_n_procs()) /= mpi_id()) cycle
!#endif
!
!        type = this%group(g)%type
!
!        select case(type)
!
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!           !X NORMAL ATOMS
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX           
!        case(TYPE_ATOM)
!           do n = 1, group_n_atoms(this%group(g))
!              i = group_nth_atom(this%group(g),n)
!#ifdef _MPI
!              if (do_parallel) then
!                 mpi_pos(:,i) = this%atoms%pos(:,i) + this%atoms%velo(:,i) * dt + 0.5_dp * dt * dt * this%atoms%acc(:,i)
!                 mpi_acc(:,i) = this%atoms%acc(:,i)
!              else
!                 this%atoms%pos(:,i) = this%atoms%pos(:,i) + this%atoms%velo(:,i) * dt + 0.5_dp * dt * dt * this%atoms%acc(:,i)
!              end if
!#else
!              this%atoms%pos(:,i) = this%atoms%pos(:,i) + this%atoms%velo(:,i) * dt + 0.5_dp * dt * dt * this%atoms%acc(:,i)
!#endif
!
!           end do
!
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!           !X CONSTRAINED ATOMS
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX           
!        case(TYPE_CONSTRAINED)
!           ! This is the SHAKE-like part of the RATTLE routine
!           ! It modifies the accelerations, so that in the next AdvanceVerlet call
!           ! this step: v(t-dt/2) = v(t-dt) + 1/2 a(t-dt) dt
!           ! automatically applies some constraint forces to the velocities.
!           ! however, a small fix is needed to make the velocities obey the constraints
!           ! exactly in the v(t) = ... step
!
!           !Store the current accelerations and positions
!           allocate(oldacc(3*size(this%group(g)%atom)),&
!                    oldpos(3*size(this%group(g)%atom)))
!           do n = 1, group_n_atoms(this%group(g))
!              i = group_nth_atom(this%group(g),n)
!              oldacc(3*n-2:3*n) = this%atoms%acc(:,i)
!              oldpos(3*n-2:3*n) = this%atoms%pos(:,i)
!           end do
!
!           !Constraints cannot currently be velocity dependent, so the gradients stored
!           !before the last velocity update step should still be valid. If velocity dependent
!           !constraints are added in future the constraint gradients will require re-calculating
!           !here.
!
!           !Zero the Lagrange multipliers
!!           do n = 1, group_n_objects(this%group(g))
!!              i = group_nth_object(this%group(g),n)
!!              this%constraint(i)%lambdaR = 0.0_dp
!!           end do
!
!           !ITERATE TO CONVERGENCE
!           iterations = 0
!           do !loop is escaped by convergence or max iterations below
!              !Reset the positions and accelerations to their initial values
!              do n = 1, group_n_atoms(this%group(g))
!                 i = group_nth_atom(this%group(g),n)
!                 this%atoms%acc(:,i) = oldacc(3*n-2:3*n)
!                 this%atoms%pos(:,i) = oldpos(3*n-2:3*n)
!              end do
!              !Update accelerations via a -> a - lambda * grad(C)_initial / mass
!              do n = 1, group_n_objects(this%group(g))
!                 i = group_nth_object(this%group(g),n)
!                 do nn = 1, this%constraint(i)%N
!                    j = this%constraint(i)%atom(nn)
!                    this%atoms%acc(:,j) = this%atoms%acc(:,j) - this%constraint(i)%lambdaR &
!                         * this%constraint(i)%old_dC_dr(3*nn-2:3*nn) &
!                         / this%atoms%mass(j)
!                 end do
!              end do
!
!              !Update positions according to new accelerations
!              do n = 1, group_n_atoms(this%group(g))
!                 i = group_nth_atom(this%group(g),n)
!                 this%atoms%pos(:,i) = this%atoms%pos(:,i) + this%atoms%velo(:,i)*dt + 0.5_dp*dt*dt*this%atoms%acc(:,i)
!              end do
!
!              !Calculate new values of constraints and gradients, and test for convergence
!              converged = .true.
!              do n = 1, group_n_objects(this%group(g))
!
!                 i = group_nth_object(this%group(g),n)
!                 call reallocate(pos,3*this%constraint(i)%N)
!                 call reallocate(velo,3*this%constraint(i)%N)
!
!                 o = this%constraint(i)%atom(1)
!                 if (this%constraint(i)%N /= 1) then
!                    do nn = 1, this%constraint(i)%N
!                       j = this%constraint(i)%atom(nn)
!                       pos(3*nn-2:3*nn) = diff_min_image(this%atoms,o,j)
!                       velo(3*nn-2:3*nn) = this%atoms%velo(:,j)
!                    end do
!                 else
!                    pos = this%atoms%pos(:,o)
!                    velo = this%atoms%velo(:,o)
!                 end if
!
!                 call constraint_calculate_values(this%constraint(i),pos,velo,this%t+dt)
!
!                 if (abs(this%constraint(i)%C) > this%constraint(i)%tol) converged = .false.
!
!              end do
!
!              iterations = iterations + 1
!
!              if (converged) exit
!              if (iterations > RATTLE_MAX_ITERATIONS) then
!                 call print_warning('AdvanceVerlet: RATTLE (SHAKE) did not converge for this group')
!                 call print(this%group(g),g)
!                 call print('')
!                 call print('Forces on atoms in this group before constraints:')
!                 do n = 1, group_n_atoms(this%group(g))
!                    i = group_nth_atom(this%group(g),n)
!                    write(line,'(i7,a,3(f0.5,1x))') i,' : ',this%atoms%mass(i)*oldacc(3*n-2:3*n)
!                    call print(line)                                        
!                 end do
!                 call print('')
!                 call print('Constraints:')
!                 call print('')
!                 do n = 1, group_n_objects(this%group(g))
!                    i = group_nth_object(this%group(g),n)
!                    call print(this%constraint(i),i)
!                 end do
!                 call system_abort('AdvanceVerlet: RATTLE (SHAKE) reached maximum iterations')
!              end if
!
!              !Update lambdas
!              do n = 1, group_n_objects(this%group(g)) !Loop over constraints in the group
!                 i = group_nth_object(this%group(g),n)
!                 m = 0.0_dp
!                 do nn = 1, this%constraint(i)%N       !Loop over each particle in the constraint
!                    j = this%constraint(i)%atom(nn)
!                    m = m+(this%constraint(i)%dC_dr(3*nn-2:3*nn).dot.this%constraint(i)%old_dC_dr(3*nn-2:3*nn)) &
!                         / this%atoms%mass(j)
!                 end do
!                 dlambda =  2.0_dp * this%constraint(i)%C / (m * dt * dt)
!                 this%constraint(i)%lambdaR = this%constraint(i)%lambdaR + dlambda
!              end do
!
!           end do
!
!           !Update stored constraint_forces
!           if (store_constraint_forces) then
!              do n = 1, group_n_objects(this%group(g))
!                 i = group_nth_object(this%group(g),n)
!                 do nn = 1, this%constraint(i)%N
!                    j = this%constraint(i)%atom(nn)
!                    constraint_force(:,j) = constraint_force(:,j) - this%constraint(i)%lambdaR &
!                         * this%constraint(i)%old_dC_dr(3*nn-2:3*nn)
!                 end do
!              end do
!           end if
!
!           deallocate(pos,velo,oldacc,oldpos)
!
!#ifdef _MPI
!           do n = 1, group_n_atoms(this%group(g))
!              i = group_nth_atom(this%group(g),n)
!              mpi_pos(:,i) = this%atoms%pos(:,i)
!              mpi_acc(:,i) = this%atoms%acc(:,i)
!              mpi_constraint_force(:,i) = constraint_force(:,i)
!           end do
!#endif          
!
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!           !X RIGID ATOMS
!           !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX           
!        case(TYPE_RIGID)
!
!!           j = group_nth_object(this%group(g),1)
!!
!!           !RCoM(t+dt) = RCoM(t) + VCoM(t+dt/2) dt
!!           this%rigidbody(j)%RCoM = this%rigidbody(j)%RCoM + this%rigidbody(j)%VCoM * dt 
!!
!!           !Do a free rotation: P'(t+dt/2), Q(t) -> P(t+dt/2), Q(t+dt)
!!           call no_squish_Free_Rotor(this%rigidbody(j)%q, this%rigidbody(j)%p, this%rigidbody(j)%model%I, dt)
!!
!!           !Write back all the positions (t+dt) and velocities (t+dt/2)
!!           call rigidbody_write_positions(this,j)
!!           call rigidbody_write_velocities(this,j)
!
!        end select
!
!     end do
!
!#ifdef _MPI
!     if (do_parallel) then
!
!        call MPI_ALLREDUCE(mpi_pos,this%atoms%pos,size(mpi_pos),MPI_DOUBLE_PRECISION,&
!                           MPI_SUM,MPI_COMM_WORLD,error_code)
!        call abort_on_mpi_error(error_code,'Advance_verlet - post pos update positions')
!
!        call MPI_ALLREDUCE(mpi_acc,this%atoms%acc,size(mpi_acc),MPI_DOUBLE_PRECISION,&
!                           MPI_SUM,MPI_COMM_WORLD,error_code)
!        call abort_on_mpi_error(error_code,'Advance_verlet - post pos update accelerations')
!        
!        if (store_constraint_forces) then
!           call MPI_ALLREDUCE(mpi_constraint_force,constraint_force,size(mpi_constraint_force),MPI_DOUBLE_PRECISION,&
!                MPI_SUM,MPI_COMM_WORLD,error_code)
!           call abort_on_mpi_error(error_code,'Advance_verlet - post pos update constraint forces')
!        end if
!        
!     end if
!#endif
!
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     !X Calculate work done and update the old positions
!     !X Also calculate work done by the thermostat
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!     ! Accumulate work done using trapezium rule
!
!     if (this%nSteps > 0) then
!        ! top up old finishing term to full term
!        this%work = this%work + 1.0_dp/2.0_dp*this%dW*dt
!        this%thermostat_work = this%thermostat_work + 1.0_dp/2.0_dp*this%thermostat_dW*dt
!     end if
!        
!     ! Calculate increment of work done this step using
!     ! dW = F(t) .dot. v(t) * dt
!     ! Note JRK 20/8/06: using dW = F(t) .dot. [r(t+dt)-r(t)] gives a dW
!     ! that is consistently too large resulting in a drift in total work
!     ! I don't really understand why this should be...
!     this%dW = 0.0_dp
!     this%thermostat_dW = 0.0_dp
!     do i = 1, this%N
!        tmp_realpos = realpos(this%atoms,i)
!        this%dW = this%dW + &
!             (this%atoms%acc(:,i) .dot. this%atoms%velo(:,i)) * dt * this%atoms%mass(i)
!        this%thermostat_dW = this%thermostat_dW + &
!             (thermostat_force(:,i) .dot. this%atoms%velo(:,i)) * dt
!
!        this%atoms%oldpos(:,i) = tmp_realpos
!     end do
!
!     call print('dW = '//this%dW, NERD)
!
!     ! Add new finishing terms to this%work and this%thermostat_work
!     this%work = this%work + 1.0_dp/2.0_dp*this%dW*dt
!     this%thermostat_work = this%thermostat_work + 1.0_dp/2.0_dp*this%thermostat_dW*dt
!
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     !X Update extended energy 
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!     my_dV_dt = 0.0_dp ! If not specified then assume potential is not time-dependant
!     if (present(dV_dt)) my_dV_dt = dV_dt
!     this%ext_energy = this%ext_energy - dt/2.0_dp*(dV_dt_old - my_dV_dt)
!     dV_dt_old = my_dV_dt
!     
!
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     !X Update the averages
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!     if (this%avg_time /= 0.0_dp) then
!        av_fact1 = exp(-dt/this%avg_time)
!        av_fact2 = 1.0_dp / ( this%avg_time/dt + 0.5_dp)
!        do i = 1, this%N
!           this%atoms%avgpos(:,i) = this%atoms%avgpos(:,i) * av_fact1
!           this%atoms%avgpos(:,i) = this%atoms%avgpos(:,i) + av_fact2 * realpos(this%atoms,i)
!
!           this%atoms%avg_ke(i) = this%atoms%avg_ke(i) * av_fact1
!           this%atoms%avg_ke(i) = this%atoms%avg_ke(i) + av_fact2 * &
!                0.5_dp*this%atoms%mass(i)*norm2(this%atoms%velo(:,i))
!        end do
!        this%avg_temp = this%avg_temp * av_fact1
!        this%avg_temp = this%avg_temp + av_fact2 * temperature(this, instantaneous=.true.)
!     end if
!
!
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     !X Recalculate the distance tables
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!     ! Do this calc_dists in parallel if we're using MPI
!     call calc_dists(this%atoms, parallel=.true.)
!
!     
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     !X Advance the timer
!     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     
!     this%t = this%t + dt
!     this%nSteps = this%nSteps + 1
!     
!     deallocate(thermostat_force)
!
!#ifdef _MPI
!     if (do_parallel) then
!        deallocate(mpi_pos, mpi_velo, mpi_acc)
!        if (store_constraint_forces) deallocat(mpi_constraint_force)
!     end if
!#endif
!
!     if (this%adaptive_thermostat) call update_thermostat(this,dt)
!
!   end subroutine advance_verlet
     
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !X
    !X Thermostat routines and other things called by AdvanceVerlet
    !X
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    !
    ! This goes after we have the new accelerations in acc, but before the velocities
    ! are updated
    !
    !%OMIT
!    subroutine thermostat1(this, dt, thermostat_force)
!      
!      type(DynamicalSystem), intent(inout)  :: this
!      real(dp),              intent(in)     :: dt
!      real(dp), dimension(:,:), intent(out) :: thermostat_force
!
!      !local variables
!      integer                               :: i, j, n_move_therm
!      real(dp)                              :: langevin_factor, langevin_fluct_factor, langevin_correct_factor
!      real(dp), allocatable, dimension(:,:) :: deltaF
!
!      integer :: pos_indices(3), thermo_region_index, region
!      real(dp) :: scale
!
!#ifdef _MPI
!      ! Make sure random numbers are the same in all processes
!      call system_resync_rng()
!#endif
!
!      ! Possibly do some thermostat stuff
!
!      if (.not.this%adaptive_thermostat) then
!
!         select case(this%thermostat)
!            
!         case(NOSE_THERM)
!            forall(i=1:this%N,this%atoms%move_mask(i) == 1 .and.this%atoms%thermostat_mask(i) == 1)  
!               this%atoms%acc(:,i) = this%atoms%acc(:,i) - this%nose*this%atoms%velo(:,i)
!               thermostat_force(:,i) = -this%nose*this%atoms%velo(:,i)*this%atoms%mass(i)
!            end forall
!            
!         case(LANGEVIN_THERM)
!            !Formula from D Quigley & MIJ Probert JChemPhys 120 p11432
!            langevin_factor         = 1.0_dp / this%thermal_tau
!            langevin_fluct_factor   = sqrt(2.0_dp * BOLTZMANN_K * this%sim_temp * langevin_factor / dt)
!            langevin_correct_factor = 1.0_dp / (1.0_dp + dt / (2.0_dp * this%thermal_tau))        
!            j = 0
!            ! Count the moving and thermostatted atoms
!            n_move_therm = count(this%atoms%move_mask == 1 .and. this%atoms%thermostat_mask == 1)      
!            ! allocate space for the fluctuating part (random forces)
!            allocate( deltaF(3,n_move_therm) )
!            
!            do i=1,this%N           
!               if (this%atoms%move_mask(i) == 1 .and. this%atoms%thermostat_mask(i) == 1) then
!                  j = j + 1
!                  deltaF(:,j) = - langevin_factor * this%atoms%mass(i) * this%atoms%velo(:,i) & !Dissipative term
!                       + sqrt(this%atoms%mass(i)) * langevin_fluct_factor * ran_normal3() !Fluctuating term
!               end if
!            end do
!            
!            ! Now zero the sum of the added forces
!            call zero_sum(deltaF)
!            
!            ! Add these forces (divided by mass) to the accelerations
!            j = 0
!            do i = 1, this%N           
!               if(this%atoms%move_mask(i) == 1.and.this%atoms%thermostat_mask(i) == 1) then
!                  j = j + 1
!                  
!                  thermostat_force(:,i) = this%atoms%mass(i)*this%atoms%acc(:,i)*&
!                       (langevin_correct_factor-1.0_dp) + langevin_correct_factor*deltaF(:,j)
!                  
!                  this%atoms%acc(:,i) = this%atoms%acc(:,i) + deltaF(:,j) / this%atoms%mass(i)
!                  ! Correction factor, from Castep
!                  this%atoms%acc(:,i) = this%atoms%acc(:,i) * langevin_correct_factor
!                  
!               end if
!            end do
!            deallocate(deltaF)
!            
!         end select
!
!      else
!         
!         if (get_value(this%atoms%properties,'thermo_region',pos_indices)) then
!            thermo_region_index = pos_indices(2)
!         end if
!
!         !Formula from D Quigley & MIJ Probert JChemPhys 120 p11432
!         j = 0
!         ! Count the moving and thermostatted atoms
!         n_move_therm = count(this%atoms%move_mask == 1 .and. this%atoms%thermostat_mask == 1)      
!         ! allocate space for the fluctuating part (random forces)
!         allocate( deltaF(3,n_move_therm) )
!         
!         do i=1,this%N           
!
!            region = this%atoms%data%int(thermo_region_index,i)
!            if (region == 0) then
!               scale = 1.0_dp
!               langevin_factor         = 1.0_dp / this%thermal_tau
!               langevin_fluct_factor   = sqrt(2.0_dp * BOLTZMANN_K * this%sim_temp * langevin_factor / dt)         
!               langevin_correct_factor = 1.0_dp / (1.0_dp + dt / (2.0_dp * this%thermal_tau))        
!            else
!               scale = this%adaptive_scale_factor(region)
!               langevin_factor         = 1.0_dp / this%adaptive_region_tau(region)
!               langevin_fluct_factor   = sqrt(2.0_dp * BOLTZMANN_K * this%sim_temp * langevin_factor / dt)         
!               langevin_correct_factor = 1.0_dp / (1.0_dp + dt / (2.0_dp * this%adaptive_region_tau(region)))   
!            end if
!
!            if(this%atoms%move_mask(i) == 1 .and. this%atoms%thermostat_mask(i) == 1) then
!               j = j + 1
!               deltaF(:,j) = - langevin_factor * this%atoms%mass(i) * this%atoms%velo(:,i) & !Dissipative term
!                    + sqrt(this%atoms%mass(i)) * langevin_fluct_factor * ran_normal3() * scale !Fluctuating term
!            end if            
!         end do
!         
!         ! Now zero the sum of the added forces
!         call zero_sum(deltaF)
!         
!         ! Add these forces (divided by mass) to the accelerations
!         j = 0
!         do i = 1, this%N           
!            if(this%atoms%move_mask(i) == 1.and.this%atoms%thermostat_mask(i) == 1) then
!               j = j + 1
!               
!               thermostat_force(:,i) = this%atoms%mass(i)*this%atoms%acc(:,i)*&
!                    (langevin_correct_factor-1.0_dp) + langevin_correct_factor*deltaF(:,j)
!
!               this%atoms%acc(:,i) = this%atoms%acc(:,i) + deltaF(:,j) / this%atoms%mass(i)
!
!               ! Correction factor, from Castep
!               this%atoms%acc(:,i) = this%atoms%acc(:,i) * langevin_correct_factor
!               
!            end if
!         end do
!         deallocate(deltaF)
!         
!      end if
!      
!    end subroutine thermostat1
!
!    !
!    ! This goes after we have updated the velocities but before we update the positions
!    !
!    !%OMIT
!    subroutine thermostat2(this, dt)
!      
!      type(DynamicalSystem), intent(inout) :: this
!      real(dp),              intent(in)    :: dt
!      !local variables
!      integer                              :: i
!      real(dp)                             :: damp_factor, currTemp
!      
!      ! Damping
!      if (this%damping) then
!         damp_factor = (1.0_dp - this%damp) / (1.0_dp + this%damp)
!         forall (i=1:this%N, this%atoms%move_mask(i) == 1.and. this%atoms%damp_mask(i) == 1) 
!            this%atoms%velo(:,i) = this%atoms%velo(:,i) * damp_factor
!         end forall
!      end if
!
!      select case(this%thermostat)
!
!      case(NOSE_THERM)
!         this%v_nose = (temperature(this, instantaneous=.true.) - this%sim_temp) / this%Q_nose
!         this%nose = this%nose + dt * this%v_nose
!         
!      case(BERENDSEN_THERM)
!         currTemp = temperature(this, instantaneous=.true.)
!         if (currTemp < NUMERICAL_ZERO) &
!          call system_abort('Cannot rescale velocities. Temperature of thermostatted atoms is zero')
!         if (this%thermal_tau /= 0.0_dp) then
!            this%rescalefactor = sqrt(1.0_dp + (dt/this%thermal_tau)*(this%sim_temp/currTemp - 1.0_dp))
!         else
!            this%rescalefactor = sqrt(this%sim_temp/currTemp)
!         end if
!         forall(i=1:this%N,this%atoms%move_mask(i) == 1.and.this%atoms%thermostat_mask(i) == 1)
!            this%atoms%velo(:,i) = this%atoms%velo(:,i) * this%rescalefactor
!         end forall
!
!      end select
!
!    end subroutine thermostat2

      
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !X
    !X I/O
    !X
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    !% Print a status line showing the current time, temperature, the mean temperature
    !% the total energy and the total momentum for this
    !% DynamicalSystem. If present, the optional
    !% 'label' parameter should be a one character label for the log lines and is printed
    !% in the first column of the output.
   subroutine ds_print_status(this, label, epot, instantaneous)
     type(DynamicalSystem),  intent(in):: this
     character(*), optional, intent(in)::label
     character(2)                      ::string
     real(dp),     optional, intent(in)::epot
     logical, optional :: instantaneous
     logical, save :: firstcall = .true.
     real(dp) :: temp

     string = " "
     if(present(label)) string = label

     if(firstcall) then
        if(present(epot)) then
             write(line, '(a14,a12,a12,a12,5a20)') "Time", "Temp", "Mean temp", &
                  "Norm(p)","Total work","Thermo work", &
                  "Ekin", "Etot", "Eext"; call print(line)
          else
             write(line, '(a14,a12,a12,a12,a12,2a20)') "Time", "Temp", "Mean temp", &
                  "Norm(p)",  "Total work","Thermo work"; call print(line)
          end if
        firstcall = .false.
     end if

     temp = temperature(this, instantaneous=instantaneous)

     if(present(epot)) then
        write(line, '(a,f12.2,2f12.4,e12.2,5e20.8)') string, this%t, &
             temp, this%avg_temp, norm(momentum(this)),&
             this%work, this%thermostat_work, kinetic_energy(this), &
             epot+kinetic_energy(this),epot+kinetic_energy(this)+this%ext_energy
     else
        write(line, '(a,f12.2,2f12.4,e12.2,2e20.8)') string, this%t, &
             temp, this%avg_temp,norm(momentum(this)), &
             this%work, this%thermostat_work

     end if
     call print(line)
   end subroutine ds_print_status

   !% Print lots of information about this DynamicalSystem in text format.
   subroutine ds_print(this,file)

      type(DynamicalSystem),  intent(inout) :: this
      type(Inoutput),optional,intent(inout) :: file
      integer                               :: i

      write(line,'(a)') 'DynamicalSystem:'
      call Print(line,file=file)

      if (.not.this%initialised) then
         write(line,'(a)')' (not initialised)'
         call Print(line,file=file)
         return
      end if

      write(line,'(a,i0)')        'Number of atoms = ',this%N ; call Print(line,file=file)
      write(line,'(a,f12.4)')     'Simulation Time = ',this%t ; call Print(line,file=file)
      write(line,'(a,f8.4)')      'Averaging Time  = ',this%avg_time ; call print(line,file=file)
      call print("")

      call print('-------------------------------------------------------------------------------', VERBOSE, file)
      call print('| Index  |   Average Position   |      Velocity        |     Acceleration     |', VERBOSE, file)
      call print('-------------------------------------------------------------------------------', VERBOSE, file)
      do i = 1, this%N
	 call print('| '//i//' |'//this%atoms%avgpos(:,i)//' |'//this%atoms%velo(:,i)//' |'//this%atoms%acc(:,i)//' |', VERBOSE, file)
      end do
      call print('-------------------------------------------------------------------------------', VERBOSE,file)

      call verbosity_push_decrement()
      call print(this%atoms,file)
      call verbosity_pop()

   end subroutine ds_print
   
   subroutine ds_file_write(this,file)

      type(Inoutput),        intent(inout) :: file
      type(DynamicalSystem), intent(inout) :: this
      
      !Header: used to make sure reading occurs from the correct place
      call write_binary('DS',file)

      !Scalar logicals
      call write_binary(this%initialised,file)
      if (.not.this%initialised) return


      !Scalar members
      call write_binary(this%N,file)                 !1
      call write_binary(this%nSteps,file)            !2
      call write_binary(this%Nrigid,file)            !3
      call write_binary(this%Nconstraints,file)      !4
      call write_binary(this%Ndof,file)              !5

      call write_binary(this%t,file)                 !6
      call write_binary(this%avg_temp,file)          !7
      call write_binary(this%avg_time,file)          !8
      call write_binary(this%dW,file)                !9
      call write_binary(this%work,file)              !10
      call write_binary(this%Epot,file)              !11
      call write_binary(this%ext_energy,file)        !12     
      call write_binary(this%thermostat_dW,file)     !13
      call write_binary(this%thermostat_work,file)   !14
      !this%initialised was written at start

      !Array members
      call write_binary(this%group_lookup,file)

      !Derived Type members
      call write_binary(this%atoms,file)

      !Derived Type Array members
      call write_binary(this%constraint,file)
      !call write_binary(this%rigidbody,file) not yet implemented... issue with pointer
      call write_binary(this%group,file)
      call write_binary(this%thermostat,file)

   end subroutine ds_file_write

   subroutine ds_file_read(this,file)

      type(Inoutput),        intent(inout) :: file
      type(DynamicalSystem), intent(inout) :: this
      logical                              :: initialised
      character(2)                         :: id

      !Make sure we are reading a dynamical system structure
      call read_binary(id,file)
      if (id /= 'DS') call system_abort('DS_File_Read: Bad DynamicalSystem structure in file '//trim(file%filename))

      !Make sure the dynamical system is uninitialised
      call finalise(this)    

      !Header
      call read_binary(initialised,file)            
      if (.not.initialised) return

      !Scalar members
      call read_binary(this%N,file)                 !1
      call read_binary(this%nSteps,file)            !2
      call read_binary(this%Nrigid,file)            !3
      call read_binary(this%Nconstraints,file)      !4
      call read_binary(this%Ndof,file)              !5

      call read_binary(this%t,file)                 !6
      call read_binary(this%avg_temp,file)          !7
      call read_binary(this%avg_time,file)          !8
      call read_binary(this%dW,file)                !9
      call read_binary(this%work,file)              !10
      call read_binary(this%Epot,file)              !11
      call read_binary(this%ext_energy,file)        !12
      call read_binary(this%thermostat_dW,file)     !13
      call read_binary(this%thermostat_work,file)   !14
      this%initialised = .true.               

      !Array members
      allocate(this%group_lookup(this%N))
      call read_binary(this%group_lookup,file)
    
      !Derived Type members
      call read_binary(this%atoms,file)

      !Derived Type Array members
      call read_binary(this%constraint,file)
      !call read_binary(this%rigidbody,file) not yet implemented... issue with pointer
      call read_binary(this%group,file)
      call read_binary(this%thermostat,file)

   end subroutine ds_file_read

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !X
   !X CONSTRAINED DYNAMICS ROUTINES
   !X
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   
   !
   ! Routines to make adding constraints easier
   !
   !% Constrain the bond between atoms i and j
   subroutine constrain_bond(this,i,j)

     type(DynamicalSystem), intent(inout) :: this
     integer,               intent(in)    :: i,j
     logical, save                        :: first_call = .true.
     integer, save                        :: BOND_FUNC
     real(dp)                             :: d

     !Do nothing for i==j
     if (i==j) then
        write(line,'(2(a,i0))')'Constrain_Bond: Tried to constrain bond ',i,'--',j
        call print_warning(line)
        return
     end if
     
     !Report bad atom indices
     if ( (i>this%N) .or. (i<1) .or. (j>this%N) .or. (j<1) ) then
        write(line,'(a,3(i0,a))')'Constrain_Bond: Cannot constrain bond ',i,'--',j,&
                                 ': Atom out of range (N=',this%N,')'
        call system_abort(line)
     end if
     
     !Register the constraint function if this is the first call
     if (first_call) then
        BOND_FUNC = register_constraint(BONDLENGTH)
        first_call = .false.
     end if
     
     !Add the constraint
     d = distance_min_image(this%atoms,i,j)
     call ds_add_constraint(this,(/i,j/),BOND_FUNC,(/d/))

   end subroutine constrain_bond

   !
   ! Routines to make adding constraints easier
   !
   !% Constrain the difference of bond length between atoms i--j and j--k
   subroutine constrain_bond_diff(this,i,j,k)

     type(DynamicalSystem), intent(inout) :: this
     integer,               intent(in)    :: i,j,k
     logical, save                        :: first_call = .true.
     integer, save                        :: BOND_DIFF_FUNC
     real(dp)                             :: d

     !Do nothing for i==j or i==k or j==k
     if (i==j.or.i==k.or.j==k) then
        write(line,'(3(a,i0))')'Constrain_Bond_Diff: Tried to constrain bond ',i,'--',j,'--',k
        call print_warning(line)
        return
     end if
     
     !Report bad atom indices
     if ( (i>this%N) .or. (i<1) .or. (j>this%N) .or. (j<1) .or. (k>this%N) .or. (k<0) ) then
        write(line,'(a,4(i0,a))')'Constrain_Bond_Diff: Cannot constrain bond ',i,'--',j,'--',k,&
                                 ': Atom out of range (N=',this%N,')'
        call system_abort(line)
     end if
     
     !Register the constraint function if this is the first call
     if (first_call) then
        BOND_DIFF_FUNC = register_constraint(BONDLENGTH_DIFF)
        first_call = .false.
     end if
     
     !Add the constraint
     d = abs(distance_min_image(this%atoms,i,j) - distance_min_image(this%atoms,j,k))
     call ds_add_constraint(this,(/i,j,k/),BOND_DIFF_FUNC,(/d/))

   end subroutine constrain_bond_diff

   !% Add a constraint to the DynamicalSystem and reduce the number of degrees of freedom,
   !% unless 'update_Ndof' is present and false.
   subroutine ds_add_constraint(this,atoms,func,data,update_Ndof)

     type(DynamicalSystem),  intent(inout) :: this
     integer,  dimension(:), intent(in)    :: atoms
     integer,                intent(in)    :: func
     real(dp), dimension(:), intent(in)    :: data
     logical,  optional,     intent(in)    :: update_Ndof
     integer                               :: i, type, g1, g2, n, new_constraint
     logical                               :: do_update_Ndof
    
     do_update_Ndof = .true. !default
     if (present(update_Ndof)) do_update_Ndof = update_Ndof

     !Have constraints been declared when the dynamicalsystem was initialised?
     if (.not.allocated(this%constraint)) &
          call system_abort('ds_add_constraint: Constraints have not been initialised')
     
     !First make sure the atoms are of TYPE_ATOM or TYPE_CONSTRAINED. If they are of 
     !TYPE_ATOM they must be in a group on their own, otherwise all the other atoms will
     !be added to the CONSTRAINED group, but never be integrated.
     do i = 1, size(atoms)
        type = Atom_Type(this,atoms(i))
        if ((type /= TYPE_ATOM) .and. (type /= TYPE_CONSTRAINED)) then
           write(line,'(2(a,i0),a)')'ds_add_constraint: Atom ',i, &
                                    ' is not a normal or constrained atom (type = ',type,')'
           call system_abort(line)
        end if
        if (type == TYPE_ATOM) then
           g1 = this%group_lookup(atoms(i))
           if (Group_N_Atoms(this%group(g1)) > 1) call system_abort(&
                'ds_add_constraint: A normal atom must be in its own group when added to a constraint')
        end if
     end do

     !Merge all the groups containing the atoms in this constraint
     g1 = this%group_lookup(atoms(1))
     call set_type(this%group(g1),TYPE_CONSTRAINED)
     do i = 2, size(atoms)
        g2 = this%group_lookup(atoms(i))
        
        call set_type(this%group(g2),TYPE_CONSTRAINED)
        call merge_groups(this%group(g1),this%group(g2))
     end do

     !Increase the constraint counter and initialise the new constraint
     this%Nconstraints = this%Nconstraints + 1
     new_constraint = this%Nconstraints
     if (this%Nconstraints > size(this%constraint)) call system_abort('ds_add_constraint: Constraint array full')
     call initialise(this%constraint(new_constraint),atoms,func,data)

     !Update the group_lookups
     do n = 1, Group_N_atoms(this%group(g1))
        i = Group_Nth_Atom(this%group(g1),n)
        this%group_lookup(i) = g1
     end do

     !Add the constraint as an object for this group
     call group_add_object(this%group(g1),new_constraint)

     !Reduce the number of degrees of freedom
     if (do_update_Ndof) this%Ndof = this%Ndof - 1

     !Call the constraint subroutine to fill in variables:
     call constraint_calculate_values_at(this%constraint(new_constraint),this%atoms,this%t)

     !Give warning if they constraint is not being satisfied
     if (abs(this%constraint(new_constraint)%C) > CONSTRAINT_WARNING_TOLERANCE) &
          call print_warning('ds_add_constraint: This constraint ('//new_constraint//') is not currently obeyed: C = '// &
          round(this%constraint(new_constraint)%C,5))

     call constraint_store_gradient(this%constraint(new_constraint))

   end subroutine ds_add_constraint

   !
   !% Replace a constraint involving some atoms with a different constraint involving
   !% the same atoms
   !
   subroutine ds_amend_constraint(this,constraint,func,data)

     type(DynamicalSystem),  intent(inout) :: this
     integer,                intent(in)    :: constraint, func
     real(dp), dimension(:), intent(in)    :: data

     call constraint_amend(this%constraint(constraint),func,data)

   end subroutine ds_amend_constraint

   !
   !% Return the distance between two atoms and the relative velocity
   !% between them projected along the bond direction.
   !% This is useful for time dependent constraints.
   !
   subroutine distance_relative_velocity(at,i,j,dist,rel_velo)
     
     type(Atoms), intent(in)  :: at
     integer,     intent(in)  :: i,j
     real(dp),    intent(out) :: dist,rel_velo
     real(dp), dimension(3)   :: r,v
     
     r = diff_min_image(at,i,j) !rj - ri
     dist = norm(r)
     
     v = at%velo(:,j) - at%velo(:,i)
     
    rel_velo = (v .dot. r) / dist
    
  end subroutine distance_relative_velocity
  
!   subroutine accumulate_kinetic_energy(this,i,ek,Ndof)
!
!     type(dynamicalsystem), intent(in)    :: this
!     integer,               intent(in)    :: i
!     real(dp),              intent(inout) :: ek, Ndof
!
!     integer :: g
!
!     g = this%group_lookup(i)
!     
!     select case(this%group(g)%type)
!        
!     case(TYPE_ATOM)
!        ek = eK + 0.5_dp*this%atoms%mass(i)*norm2(this%atoms%velo(:,i))
!        Ndof = Ndof + 3.0_dp
!        
!     case(TYPE_CONSTRAINED)
!        ek = ek + 0.5_dp*this%atoms%mass(i)*norm2(this%atoms%velo(:,i))
!        Ndof = Ndof + 3.0_dp - (real(group_n_objects(this%group(g)),dp) / real(group_n_atoms(this%group(g)),dp))
!        
!     case default
!        call system_abort('accumulate_kinetic_energy: cannot handle atoms with type = '//this%group(g)%type)
!        
!     end select
!
!   end subroutine accumulate_kinetic_energy

   function degrees_of_freedom(ds,i) result(Ndof)

     type(dynamicalsystem), intent(in) :: ds
     integer,               intent(in) :: i
     real(dp)                          :: Ndof

     integer :: g

     if (ds%atoms%move_mask(i) == 0) then
        Ndof = 0.0_dp
        return
     end if

     g = ds%group_lookup(i)

     select case(ds%group(g)%type)
        
     case(TYPE_ATOM)
        Ndof = 3.0_dp
        
     case(TYPE_CONSTRAINED)
        Ndof = 3.0_dp - (real(group_n_objects(ds%group(g)),dp) / real(group_n_atoms(ds%group(g)),dp))
        
     case default
        call system_abort('degrees_of_freedom: cannot handle atoms with type = '//ds%group(g)%type)
        
     end select

   end function degrees_of_freedom
  
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !X
   !X FIXING ATOMS
   !X
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   !% This subroutine fixes a list of atoms. If the atoms are unconstrained then it
   !% simply sets move_mask to 0 and advance_verlet does nothing with them.
   !% If the atom is constrained then the situation is more complicated: if only one
   !% atom in a constrained group is fixed then this must be implemented as 3
   !% constraints (the atom is constrained to the intersection of 3 planes). This is
   !% required so that the shake and rattle algorithms work correctly with any further
   !% constraints of which the atom is a part. Move_mask is also zeroed in this case.

   subroutine fix_atoms(this,index)

     type(dynamicalsystem), intent(inout) :: this
     integer, dimension(:), intent(in)    :: index

     integer  :: i, g, n
     logical, save :: plane_registered = .false.
     integer, save :: PLANE_FUNC

     do n = 1, size(index)

        !Zero the move_mask
        i = index(n)

        if (this%atoms%move_mask(i) == 0) cycle

        this%atoms%move_mask(i) = 0

        !Update the number of DoF
        this%Ndof = this%Ndof - 3

        !See if constraints need adding
        g = this%group_lookup(i)
        if (this%group(g)%type == TYPE_CONSTRAINED) then
           
           ! Make sure the PLANE constraint function is registered
           if (.not.plane_registered) then
              PLANE_FUNC = register_constraint(PLANE)
              plane_registered = .true.
           end if

           !Add the constraints for three intersecting planes
           call ds_add_constraint(this,(/i/),PLANE_FUNC,(/1.0_dp,0.0_dp,0.0_dp,this%atoms%pos(1,i)/))
           call ds_add_constraint(this,(/i/),PLANE_FUNC,(/0.0_dp,1.0_dp,0.0_dp,this%atoms%pos(2,i)/))
           call ds_add_constraint(this,(/i/),PLANE_FUNC,(/0.0_dp,0.0_dp,1.0_dp,this%atoms%pos(3,i)/))

        end if

     end do

   end subroutine fix_atoms

end module DynamicalSystem_module
