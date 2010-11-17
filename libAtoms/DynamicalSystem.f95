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
!% which (shallowly) copies MyAtoms into the internal atoms structure (and so 
!% MyAtoms is not required by MyDS after this call and can be finalised).
!%
!% DynamicalSystem has an integrator,
!%> 	call advance_verlet(MyDS,dt,forces)
!% which takes a set of forces and integrates the equations of motion forward
!% for a time 'dt'.
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module dynamicalsystem_module
 
   use system_module
   use units_module
   use error_module
   use linearalgebra_module
   use table_module
   use periodictable_module
   use mpi_context_module
   use atoms_module
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
      integer                               :: Nrestraints = 0         !% Number of restraints
      integer                               :: Ndof = 0                 !% Number of degrees of freedom
      real(dp)                              :: t = 0.0_dp               !% Time
      real(dp)                              :: dt = 1.0_dp              !% Last time step
      real(dp)                              :: avg_temp = 0.0_dp        !% Time-averaged temperature
      real(dp)                              :: cur_temp = 0.0_dp        !% Current temperature
      real(dp)                              :: avg_time = 100.0_dp      !% Averaging time, in fs
      real(dp)                              :: dW = 0.0_dp              !% Increment of work done this time step
      real(dp)                              :: work = 0.0_dp            !% Total work done
      real(dp)                              :: Epot = 0.0_dp            !% Total potential energy
      real(dp)                              :: Ekin = 0.0_dp            !% Current kinetic energy
      real(dp)                              :: Wkin(3, 3) = 0.0_dp      !% Current kinetic contribution to the virial
      real(dp)                              :: ext_energy = 0.0_dp      !% Extended energy
      real(dp)                              :: thermostat_dW = 0.0_dp   !% Increment of work done by thermostat
      real(dp)                              :: thermostat_work = 0.0_dp !% Total work done by thermostat
      logical         			    :: initialised = .false.

      integer                               :: random_seed              !% RNG seed, used by 'ds_save_state' and 'ds_restore_state' only.

      ! Array members
      integer,  pointer, dimension(:)       :: group_lookup !% Stores which group atom $i$ is in

      ! Derived type members
      type(Atoms), pointer                  :: atoms => NULL()

      ! Derived Type array members
      type(Constraint), allocatable, dimension(:) :: constraint, restraint
      type(RigidBody),  allocatable, dimension(:) :: rigidbody
      type(Group),      allocatable, dimension(:) :: group
      type(thermostat), allocatable, dimension(:) :: thermostat
      logical :: print_thermostat_temps = .true.

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
     module procedure single_kinetic_energy, arrays_kinetic_energy, atoms_kinetic_energy, ds_kinetic_energy
   end interface kinetic_energy

   interface kinetic_virial
      module procedure single_kinetic_virial, arrays_kinetic_virial
     module procedure atoms_kinetic_virial, ds_kinetic_virial
   endinterface kinetic_virial

   interface angular_momentum
     module procedure arrays_angular_momentum, atoms_angular_momentum, ds_angular_momentum
   end interface angular_momentum

   interface momentum
     module procedure arrays_momentum, atoms_momentum, ds_momentum
   end interface momentum

   interface add_thermostat
      module procedure ds_add_thermostat
   end interface add_thermostat

   interface update_thermostat
      module procedure ds_update_thermostat
   end interface update_thermostat

contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X INITIALISE AND FINALISE
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine ds_initialise(this,atoms_in,velocity,acceleration,constraints,restraints,rigidbodies,error)
   
      type(DynamicalSystem),              intent(inout) :: this
      type(Atoms),                        target        :: atoms_in
      real(dp), dimension(:,:), optional, intent(in)    :: velocity
      real(dp), dimension(:,:), optional, intent(in)    :: acceleration
      integer,                  optional, intent(in)    :: constraints, restraints
      integer,                  optional, intent(in)    :: rigidbodies
      integer                                           :: i,N
      integer,                  optional, intent(out)   :: error

      logical :: added_avgpos, added_avgke

      ! ---

      INIT_ERROR(error)

      ! Check to see if the object has already been initialised
      if (this%initialised) call ds_finalise(this)  

      N = atoms_in%Nbuffer
      this%N = N

      ! Check the atomic numbers
      do i = 1, atoms_in%N
         if (atoms_in%Z(i) < 1) then
            RAISE_ERROR('DS_Initialise: Atom '//i//' has not had its atomic number set', error)
         end if
      end do

      ! 'group_lookup' needs to be in the Atoms object in order to be
      ! distributed properly with domain decomposition enabled.
      ! XXX FIXME this crashes in quippy
      !call add_property(this%atoms, 'group_lookup', 0, &
      !     ptr=this%group_lookup, error=error)
      !PASS_ERROR(error)
      ! XXX FIXME this is the temporary fix
      allocate(this%group_lookup(N))
      ! allocate local arrays
      allocate(this%group(N))

      ! Now point to the atoms 
      call shallowcopy(this%atoms, atoms_in)

      ! Add single valued properties to this%atoms
      ! XXX FIXME also not quippy compatible, check
      !call set_value(this%atoms%properties, 'time', 0.0_DP)
      !if (.not. assign_pointer(this%atoms%properties, 'time', this%t)) then
      !   RAISE_ERROR("Could not assign property 'time'.", error)
      !endif

      ! Add properties for the dynamical variables to this%atoms if they don't
      ! already exist
      call add_property(this%atoms, 'mass', 0.0_DP, error=error)
      PASS_ERROR(error)
      call add_property(this%atoms, 'travel', 0, n_cols=3, error=error)
      PASS_ERROR(error)

      call add_property(this%atoms, 'move_mask', 1, error=error)
      PASS_ERROR(error)
      call add_property(this%atoms, 'damp_mask', 1, error=error)
      PASS_ERROR(error)
      call add_property(this%atoms, 'thermostat_region', 1, error=error)
      PASS_ERROR(error)

      added_avgke = .not. has_property(this%atoms, 'avgke')
      call add_property(this%atoms, 'avg_ke', 0.0_dp, error=error)
      PASS_ERROR(error)
      call add_property(this%atoms, 'velo', 0.0_dp, n_cols=3, error=error)
      PASS_ERROR(error)
      call add_property(this%atoms, 'acc', 0.0_dp, n_cols=3, error=error)
      PASS_ERROR(error)
      added_avgpos = .not. has_property(this%atoms, 'avgpos')
      call add_property(this%atoms, 'avgpos', 0.0_dp, n_cols=3, error=error)
      PASS_ERROR(error)
      call add_property(this%atoms, 'oldpos', 0.0_dp, n_cols=3, error=error)
      PASS_ERROR(error)

      ! Update pointers in this%atoms so we can use this%atoms%velo etc.
      call atoms_repoint(this%atoms)

      ! Set mass
      this%atoms%mass = ElementMass(this%atoms%Z)

      ! The input arrays must have 3N components if they are present
      ! if not, local arrays are set to zero
      
      if(present(velocity)) then 
         call check_size('Velocity',velocity,(/3,N/),'DS_Initialise',error)
         PASS_ERROR(error)
         this%atoms%velo = velocity
      end if

      if(present(acceleration)) then
         call check_size('Acceleration',acceleration,(/3,N/),'DS_Initialise',error)
         PASS_ERROR(error)
         this%atoms%acc = acceleration
      end if

      ! Initialize oldpos and avgpos
      if (added_avgpos) then
	do i = 1, this%N
	   this%atoms%avgpos(:,i) = realpos(this%atoms,i)
	end do
      endif
      this%atoms%oldpos = this%atoms%avgpos
      if (added_avgke) then
	do i=1, this%atoms%N
	  this%atoms%avg_ke(i) = 0.5_dp*this%atoms%mass(i)*norm2(this%atoms%velo(:,i))
	end do
      endif

      ! Check for constraints
      if (present(constraints)) then
         if (constraints > 0) then
            allocate(this%constraint(constraints))
            this%Nconstraints = 0
         end if
      end if

      ! Check for restraints
      if (present(restraints)) then
         if (restraints > 0) then
            allocate(this%restraint(restraints))
            this%Nrestraints = 0
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
      
      call verbosity_push_decrement(PRINT_ANAL)
      call print(this)
      call verbosity_pop()
         
   end subroutine ds_initialise


   subroutine ds_finalise(this, error)

      type(DynamicalSystem),           intent(inout)  :: this
      integer,               optional, intent(out)    :: error

      ! ---

      INIT_ERROR(error)

      if (this%initialised) then      
         call finalise_ptr(this%atoms)
         call finalise(this%group)
         deallocate(this%group_lookup)
         
         !Finalise constraints
         if (allocated(this%constraint)) call finalise(this%constraint)

         !Finalise restraints
         if (allocated(this%restraint)) call finalise(this%restraint)

         !Finalise rigid bodies
         if (allocated(this%rigidbody)) call finalise(this%rigidbody)

         call finalise(this%thermostat)

         this%N = 0
         this%nSteps = 0
         this%Nrigid = 0
         this%Nconstraints = 0
         this%Nrestraints = 0
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

      integer :: from_constraint_size, from_restraint_size, from_rigidbody_size

      ! Re-initialise the destination object to be the correct size.
      ! Note: size(from%constraint) in general is not equal to from%Nconstraints etc.
      ! Also, the constraints aren't added using add_constraint etc. since this
      ! will automatically evaluate the constraints at the current positions, which may not
      ! correspond to the values currently stored. Same for rigid bodies.

      from_constraint_size = 0
      if (allocated(from%constraint)) from_constraint_size = size(from%constraint)
      from_restraint_size = 0
      if (allocated(from%restraint)) from_restraint_size = size(from%restraint)
      from_rigidbody_size = 0
      if (allocated(from%rigidbody)) from_rigidbody_size = size(from%rigidbody)

      call initialise(to, from%atoms, constraints=from_constraint_size, restraints=from_restraint_size, rigidbodies=from_rigidbody_size)

      ! Copy over scalar members
      ! to%N is set in the initialisation         1
      to%nSteps          = from%nSteps           !2
      to%Nrigid          = from%Nrigid           !4
      to%Nconstraints    = from%Nconstraints     !5
      to%Nrestraints     = from%Nrestraints     !5
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
      if (from%Nrestraints /= 0)  to%restraint       = from%restraint
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
   subroutine ds_save_state(to, from, error)

      type(DynamicalSystem),           intent(inout)  :: to
      type(DynamicalSystem),           intent(in)     :: from
      integer,               optional, intent(out)    :: error

      INIT_ERROR(error)

      if (associated(to%atoms, from%atoms)) then
         ! Both DynamicalSystem object point to the same Atoms object, so we
         ! need to allocate our own buffer.
         allocate(to%atoms)
         to%atoms%own_this = .true.
      endif

      ! Copy over scalar members
      to%N               = from%N                !1
      to%nSteps          = from%nSteps           !2
      to%Nrigid          = from%Nrigid           !4
      to%Nconstraints    = from%Nconstraints     !5
      to%Nrestraints     = from%Nrestraints      !6
      to%Ndof            = from%Ndof             !7

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

      call atoms_copy_without_connect(to%atoms, from%atoms, error=error)
      PASS_ERROR(error)

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
      to%Nrestraints     = from%Nrestraints      !6
      to%Ndof            = from%Ndof             !7

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

   subroutine ds_add_atom_single(this,z,mass,p,v,a,t, error)
     
     type(DynamicalSystem),            intent(inout) :: this
     integer,                          intent(in)    :: z
     real(dp), optional,               intent(in)    :: mass
     real(dp), optional, dimension(3), intent(in)    :: p,v,a
     integer,  optional, dimension(3), intent(in)    :: t
     integer, optional, intent(out) :: error

     real(dp) :: my_mass, my_p(3), my_v(3), my_a(3)
     integer :: my_t(3)

     INIT_ERROR(error)
          
     my_mass = optional_default(ElementMass(Z),mass) 
     my_p = optional_default((/0.0_dp,0.0_dp,0.0_dp/),p)
     my_v = optional_default((/0.0_dp,0.0_dp,0.0_dp/),v)
     my_a = optional_default((/0.0_dp,0.0_dp,0.0_dp/),a)
     my_t = optional_default((/0,0,0/),t)

     call ds_add_atom_multiple(this, (/z/), (/my_mass/), &
                               reshape(my_p, (/3,1/)), &
                               reshape(my_v, (/3,1/)), &
                               reshape(my_a, (/3,1/)), &
                               reshape(my_t, (/3,1/)), error)
     PASS_ERROR(error)

   end subroutine ds_add_atom_single
 
   !
   ! Updated to work with groups. NEEDS TESTING
   !
   subroutine ds_add_atom_multiple(this,z,mass,p,v,a,t,error)
     
     type(DynamicalSystem),              intent(inout) :: this
     integer,  dimension(:),             intent(in)    :: z
     real(dp), dimension(:),   optional, intent(in)    :: mass
     real(dp), dimension(:,:), optional, intent(in)    :: p,v,a
     integer, dimension(:,:),  optional, intent(in)    :: t
     integer, optional, intent(out) :: error
     
     integer                                 :: oldN, newN, n, f, i
     integer,     dimension(this%N)          :: tmp_group_lookup
     type(Group), dimension(:), allocatable  :: tmp_group
     
     INIT_ERROR(error)

     oldN = this%N

     ! Check the sizes are ok
     n = size(z)
     call check_size('Position',p,(/3,n/),'DS_Add_Atoms',error) 
     PASS_ERROR(error)

     call check_size('Velocity',v,(/3,n/),'DS_Add_Atoms', error) 
     PASS_ERROR(error)

     call check_size('Acceleration',a,(/3,n/),'DS_Add_Atoms', error) 
     PASS_ERROR(error)

     newN = oldN + n
     
     !Copy all non-scalar data into the temps
     tmp_group_lookup = this%group_lookup
     
     !Adjust the size of the dynamical system's arrays
     deallocate(this%group_lookup)
     allocate(this%group_lookup(newN))
     
     ! Implement the changes in Atoms
     call add_atoms(this%atoms,pos=p,Z=z,mass=mass,velo=v,acc=a,travel=t)
     
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


   subroutine ds_remove_atom_single(this,i,error)

      type(DynamicalSystem), intent(inout) :: this
      integer,               intent(in)    :: i
      integer, intent(out), optional       :: error

      INIT_ERROR(error)
      call ds_remove_atom_multiple(this,(/i/),error)
      PASS_ERROR(error)

   end subroutine ds_remove_atom_single



   subroutine ds_remove_atom_multiple(this,atomlist_in,error)

      type(DynamicalSystem),  intent(inout)    :: this
      integer,  dimension(:), intent(in)       :: atomlist_in
      integer, intent(out), optional           :: error

      integer,  dimension(size(atomlist_in))   :: atomlist
      integer                                  :: oldN, newN, g, i, copysrc

      INIT_ERROR(error)

      !Make our own copy of the indices so that we can sort them
      atomlist = atomlist_in
      call insertion_sort(atomlist)

      !Check for repeated indices, and non-TYPE_ATOM atoms
      do i = 1, size(atomlist)
         if (Atom_Type(this,atomlist(i)) /= TYPE_ATOM) then
            RAISE_ERROR('Remove_Atoms: Atom '//atomlist(i)//' is not a normal atom', error)
         end if
         if (i > 1) then
            if (atomlist(i) == atomlist(i-1)) then
               RAISE_ERROR('Remove_Atoms: Tried to remove the same atom twice ('//atomlist(i)//')', error)
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

   subroutine ds_add_thermostat(this,type,T,gamma,Q,tau,tau_cell,p)

     type(dynamicalsystem), intent(inout) :: this
     integer,               intent(in)    :: type
     real(dp),              intent(in)    :: T
     real(dp), optional,    intent(in)    :: gamma
     real(dp), optional,    intent(in)    :: Q
     real(dp), optional,    intent(in)    :: tau
     real(dp), optional,    intent(in)    :: tau_cell
     real(dp), optional,    intent(in)    :: p

     real(dp) :: w_p, gamma_cell, mass1, mass2, volume_0
     real(dp) :: gamma_eff

     if (count( (/present(gamma), present(tau) /) ) /= 1 ) call system_abort('ds_add_thermostat: exactly one of gamma, tau must be present')

     if (present(gamma)) then
       gamma_eff = gamma
     else
       gamma_eff = 1.0_dp/tau
     endif

     if(present(p)) then
        if(present(tau_cell)) then
           gamma_cell = 1.0_dp / tau_cell
        else
           gamma_cell = gamma_eff * 0.1_dp
        endif
        select case(type)
        case(LANGEVIN_NPT,NPH_ANDERSEN)
           mass1 = 9.0_dp*abs(p)*cell_volume(this%atoms)/((gamma_cell*2*PI)**2)
           mass2 = (this%Ndof+3.0_dp)*BOLTZMANN_K*max(T,MIN_TEMP)/((gamma_cell*2*PI)**2)
           w_p = max(mass1,mass2)
        case(LANGEVIN_PR,NPH_PR)
           w_p = (this%Ndof+3.0_dp)*BOLTZMANN_K*max(T,MIN_TEMP)/((gamma_cell*2*PI)**2)/3.0_dp
        endselect
        volume_0 = cell_volume(this%atoms)
     endif
     call add_thermostat(this%thermostat,type,T,gamma_eff,Q,p,gamma_cell,w_p,volume_0)
     
   end subroutine ds_add_thermostat

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !X 
   !X Updating a thermostat
   !X
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine ds_update_thermostat(this,T,p,i)

     type(dynamicalsystem), intent(inout) :: this
     real(dp), optional,    intent(in)    :: T
     real(dp), optional,    intent(in)    :: p
     integer,  optional,    intent(in)    :: i

     real(dp) :: w_p, mass1, mass2, my_T
     integer :: my_i

     my_i = optional_default(1,i)

     if(present(p)) then
        
        my_T = optional_default(this%thermostat(my_i)%T,T)
        
        select case(this%thermostat(my_i)%type)
        case(LANGEVIN_NPT,NPH_ANDERSEN)
           mass1 = 9.0_dp*abs(p)*cell_volume(this%atoms)/((this%thermostat(my_i)%gamma_p*2*PI)**2)
           mass2 = (this%Ndof+3.0_dp)*BOLTZMANN_K*max(my_T,MIN_TEMP)/((this%thermostat(my_i)%gamma_p*2*PI)**2)
           w_p = max(mass1,mass2)
        case(LANGEVIN_PR,NPH_PR)
           w_p = (this%Ndof+3.0_dp)*BOLTZMANN_K*max(my_T,MIN_TEMP)/((this%thermostat(my_i)%gamma_p*2*PI)**2)/3.0_dp
        case default
           call print_warning('Pressure passed but thermostat does not have barostat.')
        endselect
     endif

     call update_thermostat(this%thermostat(my_i),T=T,p=p,w_p=w_p)
     
   end subroutine ds_update_thermostat

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
   pure function arrays_momentum(mass, velo, indices) result(p) ! sum(mv) 
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
   end function arrays_momentum

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
   pure function arrays_angular_momentum(mass, pos, velo, origin, indices) result(L)
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

   end function arrays_angular_momentum

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
   function DS_kinetic_energy(this, mpi_obj, error) result(ke)    ! sum(0.5mv^2)
      type(DynamicalSystem),       intent(in)   :: this
      type(MPI_context), optional, intent(in)   :: mpi_obj
      integer,           optional, intent(out)  :: error
      real(dp)                                  :: ke

      INIT_ERROR(error)
      ke = kinetic_energy(this%atoms, mpi_obj, error=error)
      PASS_ERROR(error)
    end function DS_kinetic_energy

   !% Return the total kinetic energy $E_k = \sum_{i} \frac{1}{2} m v^2$
   function atoms_kinetic_energy(this, mpi_obj, error) result(ke) ! sum(0.5mv^2)
      type(atoms),                 intent(in)   :: this
      type(MPI_context), optional, intent(in)   :: mpi_obj
      integer,           optional, intent(out)  :: error
      real(dp)                                  :: ke

      INIT_ERROR(error)
      if (.not. associated(this%mass)) call system_abort("atoms_kinetic_energy called on atoms without mass property")
      if (.not. associated(this%velo)) call system_abort("atoms_kinetic_energy called on atoms without velo property")
      ke = kinetic_energy(this%mass(1:this%Ndomain), this%velo(1:3, 1:this%Ndomain))

      if (present(mpi_obj)) then
         call sum_in_place(mpi_obj, ke, error=error)
         PASS_ERROR(error)
      endif
    end function atoms_kinetic_energy

   !% Return the kinetic energy given a mass and a velocity
   pure function single_kinetic_energy(mass, velo) result(ke)
     real(dp), intent(in) :: mass, velo(3)
     real(dp) :: ke

     ke = 0.5_dp*mass * norm2(velo)
   end function single_kinetic_energy

   !% Return the total kinetic energy given atomic masses and velocities
   pure function arrays_kinetic_energy(mass, velo) result(ke)
     real(dp), intent(in) :: mass(:)
     real(dp), intent(in) :: velo(:,:)
     real(dp) :: ke

     ke = 0.5_dp * sum(sum(velo**2,dim=1)*mass)
   end function arrays_kinetic_energy

   !% Return the total kinetic virial $w_ij = \sum_{k} \frac{1}{2} m v_i v_j$
   function DS_kinetic_virial(this, mpi_obj, error) result(kv)    ! sum(0.5mv^2)
      type(DynamicalSystem),       intent(in)   :: this
      type(MPI_context), optional, intent(in)   :: mpi_obj
      integer,           optional, intent(out)  :: error
      real(dp),          dimension(3,3)         :: kv

      INIT_ERROR(error)
      kv = kinetic_virial(this%atoms, mpi_obj, error=error)
      PASS_ERROR(error)
    end function DS_kinetic_virial

   !% Return the total kinetic virial $w_ij = \sum_{k} \frac{1}{2} m v_i v_j$
   function atoms_kinetic_virial(this, mpi_obj, error) result(kv) ! sum(0.5mv^2)
      type(atoms),                 intent(in)   :: this
      type(MPI_context), optional, intent(in)   :: mpi_obj
      integer,           optional, intent(out)  :: error
      real(dp),          dimension(3,3)         :: kv

      INIT_ERROR(error)
      if (.not. associated(this%mass)) call system_abort("atoms_kinetic_virial called on atoms without mass property")
      if (.not. associated(this%velo)) call system_abort("atoms_kinetic_virial called on atoms without velo property")
      kv = kinetic_virial(this%mass(1:this%Ndomain), this%velo(1:3, 1:this%Ndomain))

      if (present(mpi_obj)) then
         call sum_in_place(mpi_obj, kv, error=error)
         PASS_ERROR(error)
      endif
    end function atoms_kinetic_virial

   !% Return the kinetic virial given a mass and a velocity
   pure function single_kinetic_virial(mass, velo) result(kv)
     real(dp), intent(in)      :: mass, velo(3)
     real(dp), dimension(3,3)  :: kv

     kv = 0.5_dp*mass * (velo .outer. velo)
   end function single_kinetic_virial

   !% Return the total kinetic virial given atomic masses and velocities
   pure function arrays_kinetic_virial(mass, velo) result(kv)
     real(dp), intent(in)      :: mass(:)
     real(dp), intent(in)      :: velo(:,:)
     real(dp), dimension(3,3)  :: kv

     kv = 0.5_dp * sum(sum(velo**2,dim=1)*mass)
   end function arrays_kinetic_virial

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

   subroutine thermostat_temperatures(this, temps)
      type(DynamicalSystem), intent(in) :: this
      real(dp), intent(out) :: temps(:)

      integer :: i

      if (size(temps) /= size(this%thermostat)) &
	call system_abort("thermostat_temperatures needs a temps array to match size of this%thermostat() " //size(this%thermostat))

      temps = -1.0_dp
      if (this%thermostat(0)%type == LANGEVIN) then
	temps(1) = temperature(this, property='damp_mask', value=1, instantaneous=.true.)
      endif
      do i=1, size(this%thermostat)-1
	if (this%thermostat(i)%type /= NONE) temps(i+1) = temperature(this, property='thermostat_region', value=i, instantaneous=.true.)
      end do
   end subroutine thermostat_temperatures

   !% Return the temperature, assuming each degree of freedom contributes
   !% $\frac{1}{2}kT$. By default only moving and thermostatted atoms are
   !% included --- this can be overriden by setting 'include_all' to true.
   function temperature(this, property, value, include_all, instantaneous, &
        mpi_obj, error)
      type(DynamicalSystem), intent(in) :: this
      character(len=*), intent(in), optional  :: property
      integer, intent(in), optional  :: value
      logical, intent(in), optional  :: include_all
      logical, intent(in), optional  :: instantaneous
      type(MPI_context), intent(in), optional  :: mpi_obj
      integer, intent(out), optional  :: error
      real(dp)                          :: temperature

      logical ::  my_include_all, my_instantaneous
      integer                           :: i, N
      real(dp)                          :: Ndof
      integer, pointer :: property_p(:)

      my_instantaneous = optional_default(.false., instantaneous)
      my_include_all = optional_default(.false., include_all)

      if (my_instantaneous) then
         nullify(property_p)
         if (present(property)) then
            if (.not. present(value)) then
               RAISE_ERROR("temperature called with property but no value to match", error)
            endif
            if (.not. assign_pointer(this%atoms, trim(property), property_p)) then
               RAISE_ERROR("temperature failed to assign integer pointer for property '"//trim(property)//"'", error)
            endif
         endif

	temperature = 0.0_dp
	N = 0
	Ndof = 0.0_dp
	do i = 1,this%atoms%Ndomain
	   if (associated(property_p)) then
	      if (property_p(i) == value .and. this%atoms%move_mask(i) == 1) then
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

        if (present(mpi_obj)) then
           call sum_in_place(mpi_obj, temperature, error=error)
           PASS_ERROR(error)
           call sum_in_place(mpi_obj, N, error=error)
           PASS_ERROR(error)
           call sum_in_place(mpi_obj, Ndof, error=error)
           PASS_ERROR(error)
        endif

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
      this%Wkin = kinetic_virial(this)
      this%Ekin = trace(this%Wkin)/2

   end subroutine add_heat

   !% Rescale the atomic velocities to temperature 'temp'. If the current
   !% temperature is zero, we first randomise the velocites.
   subroutine rescale_velo(this, temp, mass_weighted, zero_L)
      type(DynamicalSystem), intent(inout) :: this
      real(dp),              intent(in)    :: temp
      logical, intent(in), optional        :: mass_weighted, zero_L

      logical                              ::  my_mass_weighted, my_zero_L
      real(dp)                             :: r, currTemp

      my_mass_weighted = optional_default(.true., mass_weighted)
      my_zero_L = optional_default(.false., zero_L)
      currTemp = temperature(this, instantaneous=.true.)
      call print("Rescaling velocities from "//currTemp//" K to "//temp//" K")

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

      this%Wkin = kinetic_virial(this)
      this%Ekin = trace(this%Wkin)/2
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

      my_mass_weighted = optional_default(.true., mass_weighted)
      my_zero_L = optional_default(.false., zero_L)
      currTemp = temperature(this, instantaneous=.true.)

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
    
   subroutine advance_verlet1(this,dt,virial,parallel,store_constraint_force,do_calc_dists,mpi_obj,error)

     type(dynamicalsystem), intent(inout)  :: this
     real(dp),              intent(in)     :: dt
!NB     real(dp),              intent(in)     :: f(:,:)
     real(dp),dimension(3,3), intent(in), optional :: virial
     logical, optional,     intent(in)     :: parallel
     logical, optional,     intent(in)     :: store_constraint_force, do_calc_dists
     type(MPI_context), optional, intent(in)   :: mpi_obj
     integer,           optional, intent(out)  :: error

     logical                               :: do_parallel, do_store, my_do_calc_dists
     integer                               :: i, j, g, n, ntherm
     real(dp), allocatable                 :: therm_ndof(:)

#ifdef _MPI
     include 'mpif.h'
     real(dp), dimension(:,:), allocatable :: mpi_pos, mpi_velo, mpi_acc, mpi_constraint_force
     real(dp), dimension(:,:), pointer :: constraint_force
     integer                               :: error_code
#endif
     
     INIT_ERROR(error)

     this%dt = dt

     do_parallel = optional_default(.false.,parallel)
     do_store = optional_default(.false.,store_constraint_force)
     my_do_calc_dists = optional_default(.true., do_calc_dists)
!NB     call check_size('Force',f,(/3,this%N/),'advance_verlet1')

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
           if (.not. assign_pointer(this%atoms,'constraint_force', constraint_force)) then
              RAISE_ERROR('advance_verlet1: no constraint_force property found', error)
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
     do i = 1, this%atoms%Ndomain
        j = this%atoms%thermostat_region(i)
        if (j>0 .and. j<=ntherm) therm_ndof(j) = therm_ndof(j) + degrees_of_freedom(this,i)
     end do
     if (present(mpi_obj)) then
        call sum_in_place(mpi_obj, therm_ndof, error=error)
        PASS_ERROR(error)
     endif
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
        call thermostat1(this%thermostat(0),this%atoms,dt,'damp_mask',1)
     end if

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X THERMOSTATTING
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     do i = 1, ntherm
        call thermostat1(this%thermostat(i),this%atoms,dt,'thermostat_region',i,virial)
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

              if (i <= this%atoms%Ndomain) then
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
              endif
           end do
           
        end select

     end do
!     call print("advance_verlet1 <delta vel> "//(0.5_dp*dt*sum(norm(this%atoms%acc,1))/real(this%atoms%N,dp)))

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
        call thermostat2(this%thermostat(0),this%atoms,dt,'damp_mask',1)
     end if

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X THERMOSTATTING
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     do i=1, ntherm
        call thermostat2(this%thermostat(i),this%atoms,dt,'thermostat_region',i)
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

              if (i <= this%atoms%Ndomain) then
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
              endif
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
                 if (do_store) mpi_constraint_force(:,i) = constraint_force(:,i)
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
           call MPI_ALLREDUCE(mpi_constraint_force, constraint_force, size(mpi_constraint_force),&
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

     if (my_do_calc_dists) then
        call calc_dists(this%atoms,parallel=do_parallel,error=error)
        PASS_ERROR(error)
     endif

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

   subroutine advance_verlet2(this,dt,f,virial,E,parallel,store_constraint_force,error)

     type(dynamicalsystem), intent(inout)  :: this
     real(dp),              intent(in)     :: dt
     real(dp),              intent(inout)     :: f(:,:)
     real(dp),dimension(3,3), intent(in), optional :: virial
     real(dp), optional, intent(inout)     :: E
     logical, optional,     intent(in)     :: parallel
     logical, optional,     intent(in)     :: store_constraint_force
     integer, optional, intent(out) :: error

     logical                               :: do_parallel, do_store
     integer                               :: i, g, n, ntherm
     real(dp)                              :: decay
     
#ifdef _MPI
     include 'mpif.h'
     real(dp), dimension(:,:), allocatable :: mpi_velo, mpi_acc, mpi_constraint_force
     real(dp), dimension(:,:), pointer     :: constraint_force
     integer                               :: error_code
#endif

     INIT_ERROR(error)

     this%dt = dt

     do_parallel = optional_default(.false.,parallel)
     do_store = optional_default(.false.,store_constraint_force)
     ntherm = size(this%thermostat)-1
     call check_size('Force',f,(/3,this%N/),'advance_verlet2', error=error)
     PASS_ERROR(error)

#ifdef _MPI
     if (do_parallel) then
        allocate(mpi_velo(3,this%N), mpi_acc(3,this%N))
        mpi_velo = 0.0_dp
        mpi_acc = 0.0_dp
        if (do_store) then
           allocate(mpi_constraint_force(3,this%N))
           mpi_constraint_force = 0.0_dp
           if (.not. assign_pointer(this%atoms,'constraint_force', constraint_force)) then
              RAISE_ERROR('advance_verlet2: no constraint_force property found', error)
           end if
        end if
     end if
#endif

     call add_restraint_forces(this%atoms, this%Nrestraints, this%restraint, this%t, f, E, store_restraint_force=store_constraint_force)

     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !X
     !X CONVERT FORCES TO ACCELERATIONS
     !X
     !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     forall (i=1:this%atoms%Ndomain) this%atoms%acc(:,i) = f(:,i) / this%atoms%mass(i)

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

              if (i <= this%atoms%Ndomain) then
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
              endif
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
                 if (do_store) mpi_constraint_force(:,i) = constraint_force(:,i)
              end do
           end if
#endif

        end select

     end do
!     call print("advance_verlet2 <delta vel> "//(0.5_dp*dt*sum(norm(this%atoms%acc,1))/real(this%atoms%N,dp)))

#ifdef _MPI
     ! Broadcast the new velocities, accelerations and possibly constraint forces
     if (do_parallel) then
        call MPI_ALLREDUCE(mpi_velo,this%atoms%velo,size(mpi_velo),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error_code)
        call abort_on_mpi_error(error_code,'advance_verlet1 - velocity update 2')       
        call MPI_ALLREDUCE(mpi_acc,this%atoms%acc,size(mpi_acc),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error_code)
        call abort_on_mpi_error(error_code,'advance_verlet1 - acceleration update')       
        if (do_store) then
           call MPI_ALLREDUCE(mpi_constraint_force,constraint_force,size(mpi_constraint_force),&
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
        do i = 1, this%atoms%Ndomain
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

     this%Wkin = kinetic_virial(this, error=error)
     PASS_ERROR(error)
     this%Ekin = trace(this%Wkin)/2

   end subroutine advance_verlet2

   !% Calls advance_verlet2 followed by advance_verlet1. Outside this routine the
   !% velocities will be half-stepped.
   subroutine advance_verlet(ds,dt,f,virial,E,parallel,store_constraint_force,do_calc_dists,error)

     type(dynamicalsystem), intent(inout) :: ds
     real(dp),              intent(in)    :: dt
     real(dp),              intent(inout)    :: f(:,:)
     real(dp),dimension(3,3), intent(in), optional :: virial
     real(dp),              intent(inout), optional    :: E
     logical, optional,     intent(in)    :: parallel
     logical, optional,     intent(in)    :: store_constraint_force, do_calc_dists
     logical, save                        :: first_call = .true.
     integer, optional, intent(out) :: error

     INIT_ERROR(error)
     if (first_call) then
        call print_title('SINGLE STEP VERLET IN USE')
        call print('Consider changing to the two-step integrator')
        call print_title('=')
        first_call = .false.
     end if
     call advance_verlet2(ds,dt,f,virial,E,parallel,store_constraint_force,error=error)
     PASS_ERROR(error)
     call advance_verlet1(ds,dt,virial,parallel,store_constraint_force,do_calc_dists,error=error)
     PASS_ERROR(error)

   end subroutine advance_verlet

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
   subroutine ds_print_status(this, label, epot, instantaneous, mpi_obj, error)
     type(DynamicalSystem),       intent(in)   :: this
     character(*),      optional, intent(in)   :: label
     real(dp),          optional, intent(in)   :: epot
     logical,           optional, intent(in)   :: instantaneous
     type(MPI_context), optional, intent(in)   :: mpi_obj
     integer,           optional, intent(out)  :: error

     character(len=1023) :: string
     logical, save :: firstcall = .true.
     real(dp) :: temp, my_epot, my_ekin
     real(dp) :: region_temps(size(this%thermostat))
     integer :: i

     INIT_ERROR(error)

     string = ""
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

     if (present(epot)) then
        if (present(mpi_obj)) then
           my_epot = sum(mpi_obj, epot, error=error)
           PASS_ERROR(error)
        else
           my_epot = epot
        endif
     endif

     temp = temperature(this, instantaneous=instantaneous, &
          mpi_obj=mpi_obj, error=error)
     PASS_ERROR(error)
     call thermostat_temperatures(this, region_temps)

     my_ekin = kinetic_energy(this, mpi_obj, error=error)
     PASS_ERROR(error)

     if(present(epot)) then
        write(line, '(a,f12.2,2f12.4,e12.2,5e20.8)') trim(string), this%t, &
             temp, this%avg_temp, norm(momentum(this)),&
             this%work, this%thermostat_work, my_ekin, &
             my_epot+my_ekin,my_epot+my_ekin+this%ext_energy
     else
        write(line, '(a,f12.2,2f12.4,e12.2,2e20.8)') trim(string), this%t, &
             temp, this%avg_temp, norm(momentum(this)), &
             this%work, this%thermostat_work

     end if
     call print(line)

     if (this%print_thermostat_temps) then
       if (any(region_temps >= 0.0_dp)) then
	 call print("T", nocr=.true.)
	 do i=0, size(region_temps)-1
	   if (this%thermostat(i)%type /= NONE) then
	     call print(" "// i // " " // round(this%thermostat(i)%T,2) // " " // round(region_temps(i+1),2), nocr=.true.)
	   else
	     call print(" "// i // " type=NONE", nocr=.true.)
	   endif
	 end do
	 call print("")
       endif
     endif

   end subroutine ds_print_status

   !% Print lots of information about this DynamicalSystem in text format.
   subroutine ds_print(this,file)

      type(DynamicalSystem),  intent(inout) :: this
      type(Inoutput),optional,intent(inout) :: file
      integer                               :: i

      call Print('DynamicalSystem:',file=file)

      if (.not.this%initialised) then
         
         call Print(' (not initialised)',file=file)
         return
      end if

      call print('Number of atoms = '//this%N)
      call print('Simulation Time = '//this%t)
      call print('Averaging Time  = '//this%avg_time)
      call print("")

      call print('Thermostat')
      call print(this%thermostat)

      call print('-------------------------------------------------------------------------------', PRINT_VERBOSE, file)
      call print('| Index  |   Average Position   |      Velocity        |     Acceleration     |', PRINT_VERBOSE, file)
      call print('-------------------------------------------------------------------------------', PRINT_VERBOSE, file)
      do i = 1, this%N
	 call print('| '//i//' |'//this%atoms%avgpos(:,i)//' |'//this%atoms%velo(:,i)//' |'//this%atoms%acc(:,i)//' |', PRINT_VERBOSE, file)
      end do
      call print('-------------------------------------------------------------------------------', PRINT_VERBOSE,file)

      call verbosity_push_decrement()
      call print(this%atoms,file)
      call verbosity_pop()

   end subroutine ds_print
   
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !X
   !X CONSTRAINED DYNAMICS ROUTINES
   !X
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   
   !
   ! Routines to make adding constraints easier
   !
   !% Constrain the bond angle cosine between atoms i, j, and k
   subroutine constrain_bondanglecos(this,i,j,k,c,restraint_k,tol)

     type(DynamicalSystem), intent(inout) :: this
     integer,               intent(in)    :: i,j,k
     real(dp), intent(in), optional       :: c, restraint_k, tol

     logical, save                        :: first_call = .true.
     integer, save                        :: BONDANGLECOS_FUNC
     real(dp)                             :: use_c
     real(dp) :: r_ji(3), r_jk(3)

     !Do nothing for i==j
     if (i==j .or. i == k .or. j == k) then
	call print_warning('Constrain_bondanglecos: Tried to constrain bond angle cosine '//i//'--'//j//'--'//k)
        return
     end if

     !Report bad atom indices
     if ( i>this%N .or. i<1 .or. j>this%N .or. j<1  .or. k > this%N .or. k<1) then
        call system_abort('Constrain_bondanglecos: Cannot constrain bond angle cosine '//i//'--'//j//'--'//k// &
                                 ' : Atom out of range (N='//this%N//')')
     end if

     !Register the constraint function if this is the first call
     if (first_call) then
        BONDANGLECOS_FUNC = register_constraint(BONDANGLECOS)
        first_call = .false.
     end if

     !Add the constraint
     if (present(c)) then
	use_c = c
     else
        r_ji = diff_min_image(this%atoms,j,i)
        r_jk = diff_min_image(this%atoms,j,k)
	use_c = (r_ji .dot. r_jk)/(norm(r_ji)*norm(r_jk))
     endif
     call ds_add_constraint(this,(/i,j,k/),BONDANGLECOS_FUNC,(/use_c/), restraint_k=restraint_k, tol=tol)

   end subroutine constrain_bondanglecos

   !% Constrain the bond between atoms i and j
   subroutine constrain_bondlength(this,i,j,d,restraint_k,tol)

     type(DynamicalSystem), intent(inout) :: this
     integer,               intent(in)    :: i,j
     real(dp), intent(in), optional       :: d, restraint_k, tol

     logical, save                        :: first_call = .true.
     integer, save                        :: BOND_FUNC
     real(dp)                             :: use_d

     !Do nothing for i==j
     if (i==j) then
        call print_warning('Constrain_bondlength: Tried to constrain bond '//i//'--'//j)
        return
     end if

     !Report bad atom indices
     if ( (i>this%N) .or. (i<1) .or. (j>this%N) .or. (j<1) ) then
        call system_abort('Constrain_bondlength: Cannot constrain bond '//i//'--'//j//&
                                 ': Atom out of range (N='//this%N//')')
     end if

     !Register the constraint function if this is the first call
     if (first_call) then
        BOND_FUNC = register_constraint(BONDLENGTH)
        first_call = .false.
     end if

     !Add the constraint
     use_d = optional_default(distance_min_image(this%atoms,i,j), d)
     call ds_add_constraint(this,(/i,j/),BOND_FUNC,(/use_d/), restraint_k=restraint_k, tol=tol)

   end subroutine constrain_bondlength

   !% Constrain the bond between atoms i and j
   subroutine constrain_relax_bondlength(this,i,j,t0,tau,df,di,restraint_k,tol)

     type(DynamicalSystem), intent(inout) :: this
     integer,               intent(in)    :: i,j
     real(dp), intent(in)                 :: t0, tau, df
     real(dp), intent(in), optional       :: di, restraint_k, tol

     logical, save                        :: first_call = .true.
     integer, save                        :: BOND_FUNC
     real(dp)                             :: use_di

     !Do nothing for i==j
     if (i==j) then
        call print_warning('Constrain_relax_bondlength: Tried to constrain bond '//i//'--'//j)
        return
     end if

     !Report bad atom indices
     if ( (i>this%N) .or. (i<1) .or. (j>this%N) .or. (j<1) ) then
        call system_abort('Constrain_relax_bondlength: Cannot constrain bond '//i//'--'//j//&
                                 ': Atom out of range (N='//this%N//')')
     end if

     !Register the constraint function if this is the first call
     if (first_call) then
        BOND_FUNC = register_constraint(RELAX_BONDLENGTH)
        first_call = .false.
     end if

     !Add the constraint
     use_di = optional_default(distance_min_image(this%atoms,i,j), di)
     call ds_add_constraint(this,(/i,j/),BOND_FUNC,(/use_di,df,t0,tau/), restraint_k=restraint_k, tol=tol)

   end subroutine constrain_relax_bondlength

   !% Constrain the bond between atoms i and j
   subroutine constrain_bondlength_sq(this,i,j,d,restraint_k,tol)

     type(DynamicalSystem), intent(inout) :: this
     integer,               intent(in)    :: i,j
     real(dp), intent(in), optional       :: d, restraint_k, tol

     logical, save                        :: first_call = .true.
     integer, save                        :: BOND_FUNC
     real(dp)                             :: use_d

     !Do nothing for i==j
     if (i==j) then
        call print_warning('Constrain_bondlength_sq: Tried to constrain bond '//i//'--'//j)
        return
     end if

     !Report bad atom indices
     if ( (i>this%N) .or. (i<1) .or. (j>this%N) .or. (j<1) ) then
        call system_abort('Constrain_bondlength_sq: Cannot constrain bond '//i//'--'//j//&
                                 ': Atom out of range (N='//this%N//')')
     end if

     !Register the constraint function if this is the first call
     if (first_call) then
        BOND_FUNC = register_constraint(BONDLENGTH_SQ)
        first_call = .false.
     end if

     !Add the constraint
     use_d = optional_default(distance_min_image(this%atoms,i,j), d)
     call ds_add_constraint(this,(/i,j/),BOND_FUNC,(/use_d/), restraint_k=restraint_k, tol=tol)

   end subroutine constrain_bondlength_sq

   !
   ! Routines to make adding constraints easier
   !
   !% Constrain the difference of bond length between atoms i--j and j--k
   subroutine constrain_bondlength_diff(this,i,j,k,d,restraint_k,tol)

     type(DynamicalSystem), intent(inout) :: this
     integer,               intent(in)    :: i,j,k
     real(dp), intent(in), optional       :: d
     real(dp), intent(in), optional       :: restraint_k, tol

     logical, save                        :: first_call = .true.
     integer, save                        :: BOND_DIFF_FUNC
     real(dp)                             :: use_d

     !Do nothing for i==j or i==k or j==k
     if (i==j.or.i==k.or.j==k) then
        
        call print_warning('Constrain_bondlength_Diff: Tried to constrain bond '//i//'--'//j//'--'//k)
        return
     end if
     
     !Report bad atom indices
     if ( (i>this%N) .or. (i<1) .or. (j>this%N) .or. (j<1) .or. (k>this%N) .or. (k<0) ) then
        call system_abort('Constrain_bondlength_Diff: Cannot constrain bond '//i//'--'//j//'--'//k//&
                                 ': Atom out of range (N='//this%N//')')
     end if
     
     !Register the constraint function if this is the first call
     if (first_call) then
        BOND_DIFF_FUNC = register_constraint(BONDLENGTH_DIFF)
        first_call = .false.
     end if
     
     !Add the constraint
     use_d = optional_default(abs(distance_min_image(this%atoms,i,j) - distance_min_image(this%atoms,j,k)),d)
     call ds_add_constraint(this,(/i,j,k/),BOND_DIFF_FUNC,(/use_d/), restraint_k=restraint_k, tol=tol)

   end subroutine constrain_bondlength_diff

   !
   ! Routines to make adding constraints easier
   !
   !% Constrain the energy gap of two resonance structures
   subroutine constrain_gap_energy(this,d,restraint_k,tol)

     type(DynamicalSystem), intent(inout) :: this
     real(dp), intent(in)                 :: d
     real(dp), intent(in), optional       :: restraint_k, tol

     logical, save                        :: first_call = .true.
     integer, save                        :: GAP_ENERGY_FUNC
     real(dp)                             :: use_d
     integer, allocatable                 :: all_atoms(:)
     integer                              :: i

     !Register the constraint function if this is the first call
     if (first_call) then
        GAP_ENERGY_FUNC = register_constraint(GAP_ENERGY)
        first_call = .false.
     end if
     
     !Add the constraint
     allocate(all_atoms(this%atoms%N))
     do i=1,this%atoms%N
        all_atoms(i) = i
     enddo
     call ds_add_constraint(this,all_atoms,GAP_ENERGY_FUNC,(/d/), restraint_k=restraint_k, tol=tol)
     deallocate(all_atoms)

   end subroutine constrain_gap_energy

   !% Constrain an atom to lie in a particluar plane
   subroutine constrain_atom_plane(this,i,plane_n,d,restraint_k,tol)

     type(DynamicalSystem), intent(inout) :: this
     integer,               intent(in)    :: i
     real(dp), intent(in)                 :: plane_n(3)
     real(dp), intent(in), optional       :: d, restraint_k, tol

     logical, save                        :: first_call = .true.
     integer, save                        :: PLANE_FUNC
     real(dp)                             :: use_d, plane_n_hat(3)

     !Report bad atom indices
     if ( (i>this%N) .or. (i<1)) then
        call system_abort('Constrain_atom_plane: Cannot constrain atom '//i// &
                                 ': Atom out of range (N='//this%N//')')
     end if

     !Register the constraint function if this is the first call
     if (first_call) then
        PLANE_FUNC = register_constraint(PLANE)
        first_call = .false.
     end if

     plane_n_hat = plane_n/norm(plane_n)
     use_d = optional_default((this%atoms%pos(:,i) .dot. plane_n_hat),d)

     !Add the constraint
     call ds_add_constraint(this,(/i/),PLANE_FUNC,(/plane_n,use_d/), restraint_k=restraint_k, tol=tol)

   end subroutine constrain_atom_plane

   !% Add a constraint to the DynamicalSystem and reduce the number of degrees of freedom,
   !% unless 'update_Ndof' is present and false.
   subroutine ds_add_constraint(this,atoms,func,data,update_Ndof, restraint_k, tol)

     type(DynamicalSystem),  intent(inout) :: this
     integer,  dimension(:), intent(in)    :: atoms
     integer,                intent(in)    :: func
     real(dp), dimension(:), intent(in)    :: data
     logical,  optional,     intent(in)    :: update_Ndof
     real(dp), optional,      intent(in)    :: restraint_k, tol

     integer                               :: i, type, g1, g2, n, new_constraint
     logical                               :: do_update_Ndof
    
     do_update_Ndof = .true. !default
     if (present(update_Ndof)) do_update_Ndof = update_Ndof

     !Have constraints been declared when the dynamicalsystem was initialised?

     if (present(restraint_k)) then ! restraint

	if (.not.allocated(this%restraint)) &
	     call system_abort('ds_add_constraint: Restraints array has not been allocated')

	this%Nrestraints = this%Nrestraints + 1
	new_constraint = this%Nrestraints
	if (this%Nrestraints > size(this%restraint)) call system_abort('ds_add_constraint: Constraint array full')
        if (present(tol)) call print("ds_add_constraint: Restraints will ignore the explicitly given constraint tolerance.",PRINT_ALWAYS)
	call initialise(this%restraint(new_constraint),atoms,func,data,restraint_k)
	call constraint_calculate_values_at(this%restraint(new_constraint),this%atoms,this%t)

     else ! constraint

	if (.not.allocated(this%constraint)) &
	     call system_abort('ds_add_constraint: Constraints array has not been allocated')
	
	!First make sure the atoms are of TYPE_ATOM or TYPE_CONSTRAINED. If they are of 
	!TYPE_ATOM they must be in a group on their own, otherwise all the other atoms will
	!be added to the CONSTRAINED group, but never be integrated.
	do i = 1, size(atoms)
	   type = Atom_Type(this,atoms(i))
	   if ((type /= TYPE_ATOM) .and. (type /= TYPE_CONSTRAINED)) then
	      call system_abort('ds_add_constraint: Atom '//i// &
				       ' is not a normal or constrained atom (type = '//type//')')
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
	call initialise(this%constraint(new_constraint),atoms,func,data,tol=tol)

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
     endif

   end subroutine ds_add_constraint

   !
   !% Replace a constraint involving some atoms with a different constraint involving
   !% the same atoms
   !
   subroutine ds_amend_constraint(this,constraint,func,data,k)

     type(DynamicalSystem),  intent(inout) :: this
     integer,                intent(in)    :: constraint, func
     real(dp), dimension(:), intent(in)    :: data
     real(dp), optional, intent(in) :: k

     call constraint_amend(this%constraint(constraint),func,data,k)

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

   subroutine fix_atoms(this,index,restraint_k)

     type(dynamicalsystem), intent(inout) :: this
     integer, dimension(:), intent(in)    :: index
     real(dp), optional, intent(in)       :: restraint_k

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
           call ds_add_constraint(this,(/i/),PLANE_FUNC,(/1.0_dp,0.0_dp,0.0_dp,this%atoms%pos(1,i)/),restraint_k=restraint_k)
           call ds_add_constraint(this,(/i/),PLANE_FUNC,(/0.0_dp,1.0_dp,0.0_dp,this%atoms%pos(2,i)/),restraint_k=restraint_k)
           call ds_add_constraint(this,(/i/),PLANE_FUNC,(/0.0_dp,0.0_dp,1.0_dp,this%atoms%pos(3,i)/),restraint_k=restraint_k)

        end if

     end do

   end subroutine fix_atoms

end module DynamicalSystem_module
