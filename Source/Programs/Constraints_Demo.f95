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

program Constraints_Demo

! This program demonstrates constrained dynamics using LOTF. It shows how to 
! apply constraints using function which are already available and also how
! to create your own.

  use libAtoms_module
  use LOTF_module

  implicit none

  type(Atoms)                           :: at
  type(DynamicalSystem)                 :: ds
  type(Inoutput)                        :: movie
  real(dp), dimension(:,:), allocatable :: forces
  integer                               :: i, BOND_FUNC, SPHERE_FUNC

  external SPHERE ! The user defined constraint subroutine

  call System_Initialise(seed = 5671232)

  !Set up some atoms
  call Initialise(at,500,Make_Lattice(50.0_dp))
  !The first 250 will be oxygen atoms constrained to the surface of a sphere
  !centred on the origin, radius 10.0
  do i = 1, 250
     at%Z(i) = Atomic_Number('O')
     at%pos(:,i) = random_unit_vector() * 10.0_dp
  end do
  !The last 250 will be randomly placed hydrogen molecules (bond length = 0.7413 A)
  do i = 251,500,2
     at%Z(i) = Atomic_Number('H')
     at%Z(i+1) = Atomic_Number('H')
     at%pos(:,i) = 0.0_dp
     call randomise(at%pos(:,i),50.0_dp)
     at%pos(:,i+1) = at%pos(:,i) + random_unit_vector() * 0.7413_dp
  end do

  !Initialise the dynamical system with the constraints:
  !We have 250 single particle constraints (atoms on the sphere)
  !      + 125 two particle constraints    (bonds)
  !      = 375 in total

  call Initialise(ds,at,constraints = 375)
  
  !Now we create the constraints: Start by registering the constraint functions
  BOND_FUNC = Register_Constraint(BONDLENGTH) !This is pre-coded in Constraints.f95
  SPHERE_FUNC = Register_Constraint(SPHERE)   !This is user defined below, and declared
                                              !as external above. 
                                              !Up to 20 different types of constraint can be 
                                              !registered without editing the source code.
  
  !Add the constraints one by one
  !The arguments are: The DynamicalSystem, the atoms involved, the subroutine pointer, constraint data.
  !See the comments in Constraints.f95 for other function, possible atom orders, and required data.
  do i = 1, 250
     call DS_Add_Constraint(ds,(/i/),SPHERE_FUNC,(/0.0_dp,0.0_dp,0.0_dp,10.0_dp/)) !Centre = 0,0,0 radius = 10
  end do
  do i = 251, 500, 2
     call DS_Add_Constraint(ds,(/i,i+1/),BOND_FUNC,(/0.7413_dp/)) !Bond length = 0.7413 A
  end do

  !Give the atoms some initial velocities
  call rescale_velo(ds,300.0_dp)

  !Free up memory used by empty groups
  call DS_Free_Groups(ds)

  !Set up a movie file so we can see what happens
  call Initialise(movie,'Constraints_Demo.xyz')  

  !Now we can do dynamics as usual
  allocate(forces(3,ds%N))
  forces = 0.0_dp     

  do i = 0, 10000

     if (mod(i,50)==0) then
        call Calc_Connect(ds%atoms)     
        call Print_xyz(ds%atoms,movie)
        call ds_print_status(ds)
     end if

     call advance_verlet(ds,1.0_dp,forces)

  end do

  !Finalise all the objects
  deallocate(forces)
  call Finalise(movie)
  call Finalise(ds)
  call Finalise(at)

  call System_Finalise

end program Constraints_Demo

!
!User defined constraint to keep an atom on a sphere.
!data(1:3) contains the coordinates of the centre of the sphere
!data(4) contains the radius of the sphere
!
! The constraint function is C = |r - r0|^2 - a^2
!
! where r0 is the centre, a is the radius
!
! To find out more about writing your own constraint, look at the comments near
! the end of Constraints.f95
!
subroutine SPHERE(pos, velo, t, data, C, dC_dr, dC_dt)

  use System_module !for dp definition
  use linearalgebra_module !for .dot. definition

  real(dp), dimension(:),         intent(in)  :: pos, velo, data
  real(dp),                       intent(in)  :: t !not used, since this is not time dependent
  real(dp),                       intent(out) :: C
  real(dp), dimension(size(pos)), intent(out) :: dC_dr
  real(dp),                       intent(out) :: dC_dt
  !local variables
  real(dp)                                    :: d(3), radius2

  if(size(pos) /= 3) call System_Abort('SPHERE: Exactly 1 atom must be specified')
  if(size(data) /= 4) call System_Abort('SPHERE: "data" must contain exactly 4 values')

  d = pos - data(1:3)
  radius2 = data(4)*data(4)
  
  C = normsq(d) - radius2
  dC_dr = 2.0_dp * d
  dC_dt = dC_dr .dot. velo

end subroutine SPHERE
