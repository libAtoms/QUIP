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
!X  Structures module
!X
!%  A collection of utility functions that generate useful Atoms structures
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! $Id: Structures.f95,v 1.17 2008-06-17 11:03:47 jrk33 Exp $
! $Log: not supported by cvs2svn $
! Revision 1.16  2008/05/14 14:19:22  jrk33
! Trivial typo fix
!
! Revision 1.15  2008/05/13 15:50:29  jrk33
! Only set masses if mass property exists
!
! Revision 1.14  2008/05/05 13:39:59  jrk33
! Changed table_allocate calls to reflect changes to Table.f95
!
! Revision 1.13  2008/04/29 11:04:22  jrk33
! Larger supercell required in unit_slab() for some orientations
!
! Revision 1.12  2008/04/09 10:32:02  jrk33
! Added tube_radius function for nanotube average radius
!
! Revision 1.11  2008/02/04 14:57:21  jrk33
! Explicitly set atoms%mass for various structures
!
! Revision 1.10  2008/02/04 14:07:34  jrk33
! Decoupled atomic number and mass to allow mixed isotope simulations. New mass property must be set
!
! Revision 1.9  2007/10/22 14:40:01  gc121
! added the transform() function to apply a linear transformation to the atoms and lattice of a structure
!
! Revision 1.8  2007/09/21 14:31:00  nb326
! Move structures_module into libatoms_module
!
! Revision 1.7  2007/09/06 14:25:17  gc121
! put in some missing finalise statements
!
! Revision 1.6  2007/08/30 14:29:07  nb326
! Moved diamond, fcc, supercell, make_structure from Atoms.f95, support for multicomponent in diamond, unit_slab, slab_...
!
! Revision 1.5  2007/08/20 15:19:04  saw44
! Added water() function to create TIP3P water molecules
!
! Revision 1.4  2007/07/19 15:56:31  gc121
! disloc_malc still doesnt work sensibly for an edge disloc
!
! Revision 1.3  2007/07/19 12:45:12  jrk33
! Added slab and graphene structure generation routines
!
! Revision 1.2  2007/07/17 08:44:23  gc121
! added screw component to edge_disloc
!
! Revision 1.1  2007/07/16 21:35:00  gc121
! some more complex structures. disloc_malc does not yet work for edge dislocs
!


module  structures_module
  use linearalgebra_module
  use atoms_module
  use clusters_module
  implicit none
  public

  interface slab
     module procedure slab_width_height_nz, slab_nx_ny_nz
  end interface

contains

  !X 
  !X Edge dislocation generator
  !X
  subroutine fcc_11b2_edge_disloc(at, a0, n1, n2, n3)
    type(Atoms), intent(out) :: at
    real(dp),    intent(in)  :: a0 ! cubic lattice constant
    integer,     intent(in)  :: n1, n2, n3 ! supercell dimensions

    type(atoms) :: at1
    type(table) :: removelist
    integer :: i

    call fcc_z111_ortho(at1, a0)
    call supercell(at, at1, n1, n2, n3)

    call map_into_cell(at)

    do i= 1,at%N
       if( at%pos(3,i) > 0.0_dp) then
          if(abs(at%pos(1,i)) < 0.1_dp) then
             call append(removelist, i)
          else 
             at%pos(1,i) = at%pos(1,i)-sign(at%pos(1,i))*min(a0/sqrt(2.0_dp)/4.0_dp, at%pos(3,i)/10.0_dp)
             at%pos(2,i) = at%pos(2,i)-sign(at%pos(1,i))*min(a0*sqrt(3.0_dp/2.0_dp)/4.0_dp, at%pos(3,i)/10.0_dp)
          end if
       end if
    end do

    if(removelist%N > 0) then
       call allocate(removelist) ! make the table size exact
       call print('Removing '//removelist%N//' atoms')
       call remove_atoms(at, int_part(removelist,1))
    end if


    call set_lattice(at, at%lattice+reshape((/ &
         20.0_dp,0.0_dp,0.0_dp, &
         0.0_dp, 0.0_dp,0.0_dp, &
         0.0_dp, 0.0_dp,20.0_dp/), (/3, 3/)))
         
    call finalise(at1)
    call finalise(removelist)

  end subroutine fcc_11b2_edge_disloc


  !X 
  !X Dislocation generator, formulas from Malcolm Heggie
  !X
  subroutine fcc_disloc_malc(at, a0, nu, n1, n2, n3, d, type)
    type(Atoms), intent(out) :: at
    real(dp),    intent(in)  :: a0, nu ! cubic lattice constant, poisson ratio
    integer,     intent(in)  :: n1, n2, n3, d ! supercell dimensions, distance between partials
    character(*) :: type

    type(atoms) :: at1
    type(table) :: removelist
    real(dp) :: x, y, z, ux, uy, uz, r2, b_screw, b_edge, b_edge_p, theta, ox(2), r_ij
    integer :: p, i, n, j

    call fcc_z111_ortho(at1, a0)
    call supercell(at, at1, n1, n2, n3)

    ox(1) = -real(d)*sqrt(3.0_dp/2.0_dp)*a0
    ox(2) =  real(d)*sqrt(3.0_dp/2.0_dp)*a0

    if(type == 'edge') then
       b_screw = a0*sqrt(3.0_dp/2.0_dp)/2.0_dp
       b_edge = a0/sqrt(2.0_dp)/2.0_dp
    else if (type == 'screw') then
       b_screw = a0/2.0_dp/sqrt(2.0_dp)
       b_edge = 0.0_dp
    else
       call system_abort("screw_disloc(): unknown type `"//type//"'")
    end if

    do p=1,2
       b_edge_p = b_edge*0.5_dp + a0/2.0_dp/sqrt(6.0_dp)*sign(p-1.5_dp)
       do i=1,at%N
          x = at%pos(1,i)-ox(p)-0.1
          y = at%pos(2,i)
          z = at%pos(3,i)-0.1
          theta = atan2(z,x)
          r2 = x*x+z*z
          ux = b_edge_p/(2.0_dp*PI)*(theta+0.5_dp*x*z/r2/(1.0_dp-nu))
          uy = b_screw/(2.0_dp*PI)*theta
          uz = -0.25_dp*b_edge_p/(2.0_dp*PI)*((1.0_dp-2.0_dp*nu)*log(r2)+(x*x-z*z)/r2)/(1.0_dp-nu)

          at%pos(:,i) = at%pos(:,i)+(/ux, uy, uz/)
       end do
    end do

    call set_cutoff(at, 2.1_dp)
    call calc_connect(at)

    do i= 1,at%N
       do n = 1,atoms_n_neighbours(at, i)
          j = atoms_neighbour(at, i, n, distance=r_ij)
          if(r_ij < 2.0_dp .and. i > j) then
             call print('i= '//i//' j= '//j//' r_ij = '//r_ij)
             call append(removelist, j)
          end if
       end do
    end do

    !call print(at)

    if(removelist%N > 0) then
       call allocate(removelist) ! make the table size exact
       call print('Removing '//removelist%N//' atoms')
       call remove_atoms(at, int_part(removelist,1))
       call calc_connect(at)
    end if

    call set_lattice(at, at%lattice+reshape((/ &
         20.0_dp,0.0_dp,0.0_dp, &
         0.0_dp, 0.0_dp,0.0_dp, &
         0.0_dp, 0.0_dp,20.0_dp/), (/3, 3/)))

    call finalise(at1)
    call finalise(removelist)
  end subroutine fcc_disloc_malc

  subroutine disloc_noam(at, p, l, b, close_threshold)
    type(Atoms), intent(inout) :: at
    real(dp) :: p(3), l(3), b(3)
    real(dp), optional :: close_threshold

    real(dp) :: l_hat(3), b_hat(3), r1_hat(3), r2_hat(3)
    real(dp) :: theta
    real(dp) :: delta_p(3), delta_p_rot(3)
    real(dp) :: r

    integer :: i_at, j_at, j
    logical, allocatable :: atoms_remove(:)
    integer, allocatable :: remove_list(:)

    l_hat = l/norm(l)
    b_hat = b/norm(b)

    if ((l_hat .dot. b_hat) .feq. 1.0_dp) then ! screw
      r1_hat = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
      r1_hat = r1_hat - l_hat * (r1_hat .dot. l_hat)
      if (norm(r1_hat) .feq. 0.0_dp) then
	r1_hat = (/ 0.0_dp, 1.0_dp, 0.0_dp /)
	r1_hat = r1_hat - l_hat * (r1_hat .dot. l_hat)
      endif
      r1_hat = r1_hat / norm(r1_hat)
      r2_hat = l_hat .cross. r1_hat
    else
      r2_hat = l_hat .cross. b_hat
      r2_hat = r2_hat/norm(r2_hat)
      r1_hat = r2_hat .cross. l_hat
    endif

    call print("disloc_noam: r1_hat " // r1_hat, VERBOSE)
    call print("disloc_noam: r2_hat " // r2_hat, VERBOSE)
    call print("disloc_noam: l_hat " // l_hat, VERBOSE)

    do i_at = 1, at%N
      delta_p = at%pos(:,i_at) - p
      if (norm(delta_p) .feq. 0.0_dp) then
	theta = 0.0_dp
      else
	delta_p_rot(1) = delta_p .dot. r1_hat
	delta_p_rot(2) = delta_p .dot. r2_hat
	delta_p_rot(3) = delta_p .dot. l_hat
	theta = atan2(delta_p_rot(2), delta_p_rot(1))
      endif
      call print("atom " // i_at // " pos " // at%pos(:,i_at) // " delta_p " // delta_p // " theta " // theta, ANAL)
      delta_p = b*theta/(2.0_dp*PI)
      at%pos(:,i_at) = at%pos(:,i_at) + delta_p
    end do

    allocate(atoms_remove(at%N))
    atoms_remove = .false.

    call calc_connect(at)
    do i_at = 1, at%N
      if (atoms_remove(i_at)) cycle
      do j=1, atoms_n_neighbours(at, i_at)
	j_at = atoms_neighbour(at, i_at, j, distance=r)
	if (atoms_remove(j_at)) cycle
	if (present(close_threshold)) then
	  if (r < close_threshold) then
	    if (j_at == i_at) then
	      call print("WARNING: disloc_noam found atom too close to itself", ERROR)
	    else
	      atoms_remove(j_at) = .true.
	    endif
	  endif
	else
	  if (r < 0.5_dp*bond_length(at%Z(i_at), at%Z(j_at))) then
	    if (i_at == j_at) then
	      call print("WARNING: disloc_noam found atom too close to itself", ERROR)
	    else
	      atoms_remove(j_at) = .true.
	    endif
	  endif
	endif ! present(close_threshold)
      end do ! j
    end do ! i_at

    if (count(atoms_remove) > 0) then
      allocate(remove_list(count(atoms_remove)))
      j_at = 0
      do i_at=1, at%N
	if (atoms_remove(i_at)) then
	  j_at = j_at + 1
	  remove_list(j_at) = i_at
	endif
      end do
      call remove_atoms(at, remove_list)

      deallocate(remove_list)
    end if

    deallocate(atoms_remove)

  end subroutine disloc_noam

  subroutine fcc_z111_ortho(at, a0)
    type(Atoms), intent(out) :: at
    real(dp),    intent(in)  :: a0 ! cubic lattice constant
    real(dp) :: lat(3,3), b

    b = a0*sqrt(2.0_dp)/2.0_dp
    lat(:,1) = (/b, 0.0_dp, 0.0_dp/)
    lat(:,2) = (/0.0_dp, b*sqrt(3.0_dp), 0.0_dp/)
    lat(:,3) = (/0.0_dp, 0.0_dp, a0*sqrt(3.0_dp)/)

    call initialise(at, 6, lat)
    at%pos(:,1) = lat .mult. (/0.0_dp,0.0_dp,0.0_dp/)
    at%pos(:,2) = lat .mult. (/0.5_dp,0.5_dp,0.0_dp/)
    at%pos(:,3) = lat .mult. (/0.0_dp, 2.0_dp/3.0_dp, 1.0_dp/3.0_dp/)
    at%pos(:,4) = lat .mult. (/0.5_dp, 1.0_dp/6.0_dp, 1.0_dp/3.0_dp/)
    at%pos(:,5) = lat .mult. (/0.0_dp, 1.0_dp/3.0_dp, 2.0_dp/3.0_dp/)
    at%pos(:,6) = lat .mult. (/0.5_dp, 5.0_dp/6.0_dp, 2.0_dp/3.0_dp/)
  end subroutine fcc_z111_ortho

  subroutine fcc_z111(at, a0)
    type(Atoms), intent(out) :: at
    real(dp),    intent(in)  :: a0 ! cubic lattice constant
    real(dp) :: lat(3,3), b

    b = a0*sqrt(2.0_dp)/2.0_dp
    lat(:,1) = (/b, 0.0_dp, 0.0_dp/)
    lat(:,2) = (/b*0.5_dp, b*sqrt(3.0_dp)/2.0_dp, 0.0_dp/)
    lat(:,3) = (/0.0_dp, 0.0_dp, a0*sqrt(3.0_dp)/)

    call initialise(at, 3, lat)
    at%pos(:,1) = lat .mult. (/0.0_dp,0.0_dp,0.0_dp/)
    at%pos(:,2) = lat .mult. (/1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp/)
    at%pos(:,3) = lat .mult. (/2.0_dp/3.0_dp, 2.0_dp/3.0_dp, 2.0_dp/3.0_dp/)
  end subroutine fcc_z111


 !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 ! 
 ! Make a slab of diamond structure, facing in any direction. Axes is 3x3
 ! matrix with column vectors being axes of slab. For a crack (abc)[def], 
 ! axes(:,2)=[abc], i.e. plane that is opened by crack and axes(:,3)=[def]
 ! is crack tip line, perpendicular to surface and to plane in which 
 ! crack propagates. Lattice constant is 'a'. Slab axes returned as
 ! columns in the matrix 'axes'.
 !
 !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

 subroutine unit_slab(myatoms, axes, a, atnum, lat_type)
   type(Atoms), intent(OUT) :: myatoms
   real(dp), intent(IN), dimension(3,3) :: axes
   real(dp), intent(IN) :: a ! Lattice vector
   integer, intent(in), optional :: atnum(:) ! atomic numbers
   character(len=*), intent(inout), optional :: lat_type  ! lattice type (diamond, bcc, fcc)

   integer :: i, Nrep(3) ! Number of repeats in x,y and z
   real(dp), dimension(3) :: x, y, z, t, d
   real(dp), dimension(3,3) :: rot
   type(Atoms) :: at
   character(20) :: my_lat_type
   
   i = 0
   Nrep = (/ 1,1,1 /)

   my_lat_type = optional_default('diamond', lat_type)
   call print("unit_slab : Lattice type " // trim(my_lat_type))    

   ! Need special cases for some surfaces to ensure we only get one unit cell
   if (all(axes(:,2) == (/ 1,1,1 /)) .and. all(axes(:,3) == (/ 1,-1,0 /))) & 
        Nrep = (/ 2,1,2 /)
   if (all(axes(:,2) == (/ 1,1,0 /)) .and. all(axes(:,3) == (/ 0,0,-1 /))) & 
        Nrep = (/ 2,2,1 /)
   if (all(axes(:,2) == (/ 1,1,0 /)) .and. all(axes(:,3) == (/ 1,-1,0 /))) &
        Nrep = (/ 2,1,2 /)

   if(trim(my_lat_type).eq.'diamond') then
      call diamond(at, a, atnum)
   elseif(trim(my_lat_type).eq.'bcc') then
      call bcc(at, a, atnum(1))
   elseif(trim(my_lat_type).eq.'fcc') then
      call fcc(at, a, atnum(1))
   else
      call system_abort('unit_slab: unknown lattice type '//my_lat_type)
   endif
   x = axes(:,1);   y = axes(:,2);    z = axes(:,3)
   rot(1,:) = x/norm(x);  rot(2,:) = y/norm(y);  rot(3,:) = z/norm(z)

   ! Rotate atom positions and lattice
   at%pos = rot .mult. at%pos
   at%lattice = rot .mult. at%lattice

   call supercell(myatoms, at, 5, 5, 5)
   call zero_sum(myatoms%pos)


   ! Project rotated lattice onto <xyz> axes
   do i=1,3
      myatoms%lattice(:,i) = (axes(1,i)*at%lattice(:,1) + &
           axes(2,i)*at%lattice(:,2) + &
           axes(3,i)*at%lattice(:,3))/Nrep(i)
   end do
   call set_lattice(myatoms, myatoms%lattice)

   ! Form primitive cell by discarding atoms with 
   ! lattice coordinates outside range [-0.5,0.5]
   d = (/ 0.01_dp,0.02_dp,0.03_dp /) ! Small shift to avoid conincidental alignments
   i = 1
   do
      t = myatoms%g .mult. (myatoms%pos(:,i) + d)
      if (any(t <= -0.5_dp) .or. any(t >= 0.5_dp)) then
         call remove_atoms(myatoms, i)
         i = i - 1 ! Retest since we've removed an atom
      end if
      if (i == myatoms%N) exit
      i = i + 1
   end do

   call map_into_cell(myatoms)

   call finalise(at)

 end subroutine unit_slab


 !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 ! 
 ! Slab with dimensions width x height in xy plane and nz layers deep.
 ! Origin is at centre of slab.
 !
 !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

 subroutine slab_width_height_nz(myslab, axes, a, width, height, nz, atnum, lat_type)
   type(Atoms), intent(out) :: myslab
   real(dp), intent(in), dimension(3,3) :: axes
   real(dp), intent(in) :: a, width, height
   integer, intent(in) :: nz ! Number of layers
   integer, intent(in), optional :: atnum(:) ! atomic numbers to use
   character(len=*),   optional  ::  lat_type 

   type(Atoms) :: unit, layer
   integer nx, ny

   call unit_slab(unit, axes, a, atnum, lat_type)

   nx = int(floor(width/unit%lattice(1,1)))
   ny = int(floor(height/unit%lattice(2,2)))

   if (.not. (nx > 0 .and. ny > 0 .and. nz > 0)) then
      write(line,'(a, i0, i0, i0)') 'Error in slab: nx,ny,nz = ', nx, ny, nz
      call system_abort(line)
   end if

   call supercell(layer, unit, nx, ny, 1)
   ! z layers last for symmetrise
   call supercell(myslab, layer, 1, 1, nz)
   call zero_sum(myslab%pos)
   call finalise(unit, layer)

 end subroutine slab_width_height_nz


 !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 ! 
 ! Slab of size nx x ny x nx primitive cells.
 ! Origin is at centre of slab.
 !
 !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

 subroutine slab_nx_ny_nz(myslab, axes, a, nx, ny, nz, atnum, lat_type)
   type(Atoms), intent(out) :: myslab
   real(dp), intent(in), dimension(3,3) :: axes
   real(dp), intent(in) :: a
   integer, intent(in) :: nx, ny, nz
   integer, intent(in), optional :: atnum(:)
   character(len=*), intent(inout), optional ::   lat_type

   type(Atoms) :: unit, layer

   call unit_slab(unit, axes, a, atnum, lat_type)
   call supercell(layer, unit, nx, ny, 1)
   ! z layers last for symmetrise
   call supercell(myslab, layer, 1, 1, nz)
   call zero_sum(myslab%pos)
   call finalise(unit, layer)

 end subroutine slab_nx_ny_nz


 !% Cubic graphene unit cell with lattice parameter 'a'.
  function Graphene_Cubic(a) result(cube)
    real(dp), intent(in) :: a
    type(Atoms) :: cube

    call Atoms_Initialise(cube, 4, &
         reshape((/ 3.0_dp*a, 0.0_dp,   0.0_dp, &
                    0.0_dp,   sqrt(3.0_dp)*a, 0.0_dp, &
                    0.0_dp,   0.0_dp,   10.0_dp/), (/3, 3/)))

    cube%Z = 6
    cube%pos(:,1) = a*(/0.0_dp, sqrt(3.0_dp)/2.0, 0.0_dp/) 
    cube%pos(:,2) = a*(/0.5_dp, 0.0_dp, 0.0_dp/) 
    cube%pos(:,3) = a*(/1.5_dp, 0.0_dp, 0.0_dp/) 
    cube%pos(:,4) = a*(/2.0_dp, sqrt(3.0_dp)/2.0, 0.0_dp/)
    
  end function Graphene_Cubic

  !% Construct a slab of graphene of a given with and height, at a given angle.
  !% 'a' is lattice parameter.
  subroutine Graphene_Slab(slab, a, theta, width, height)
    real(dp), intent(in) :: a, theta, width, height
    type(Atoms), intent(out) :: slab

    type(Atoms) :: unit, tmp_slab
    type(Table) :: keep_list
    integer :: i, nx, ny, nz
    real(dp) :: rot(3,3), lattice(3,3), d(3), t(3)

    unit = Graphene_Cubic(a)

    rot = reshape((/ cos(theta), -sin(theta), 0.0_dp, &
                     sin(theta),  cos(theta), 0.0_dp, &
                     0.0_dp,      0.0_dp,     1.0_dp /), (/3,3/))

    unit%pos = rot .mult. unit%pos
    lattice = rot .mult. unit%lattice
    call Atoms_Set_Lattice(unit, lattice)

    ! Choose nx and ny so we can cut a width by height square out of the
    ! rotated slab
    call Fit_Box_in_Cell(width, height, norm(unit%lattice(:,3)), unit%lattice, nx, ny, nz)
    call Supercell(tmp_slab, unit, nx, ny, 1)
    call Zero_Sum(tmp_slab%pos)

    ! Final lattice for unitslab is width by height
    lattice = 0.0_dp
    lattice(1,1) = width
    lattice(2,2) = height
    lattice(3,3) = tmp_slab%lattice(3,3)
    call Atoms_Set_Lattice(tmp_slab, lattice)

    ! Form primitive cell by discarding atoms with 
    ! lattice coordinates outside range [-0.5,0.5]
    d = (/ 0.01, 0.02, 0.03 /)
    call Table_Allocate(keep_list, 1, 0, 0, 0, tmp_slab%N)

    do i=1,tmp_slab%N
       t = tmp_slab%g .mult. (tmp_slab%pos(:,i) + d)
       if (.not. (any(t < -0.5_dp) .or. any(t >= 0.5_dp))) &
            call Append(keep_list, i)
    end do

    call Atoms_Initialise(slab, keep_list%N, tmp_slab%lattice)
    slab%Z = tmp_slab%Z(keep_list%int(1,1:keep_list%N))

    if (has_property(slab, 'mass')) &
         forall (i=1:slab%N) slab%mass(i) = ElementMass(slab%Z(i))

    slab%pos = tmp_slab%pos(:,keep_list%int(1,1:keep_list%N))

    call Finalise(unit)
    call Finalise(tmp_slab)
    call finalise(keep_list)
  end subroutine Graphene_Slab

  !% Construct a graphene sheet of index $(n,m)$ with lattice constant 'a' with 'rep_x' 
  !% repeat units in the $x$ direction and 'rep_y' in the $y$ direction.
  subroutine Graphene_Sheet(sheet, a, n, m, rep_x, rep_y)
    type(Atoms), intent(out) :: sheet
    real(dp), intent(in) :: a
    integer, intent(in) :: n, m, rep_x, rep_y

    type(Atoms) :: unit, unitsheet
    real(dp) :: a1(2), a2(2), x(2), y(2), lattice(3,3), theta, rot(3,3), d(3), t(3), normx, normy
    integer :: nx, ny, nz, i

    unit = Graphene_Cubic(a)

    a1 = a*(/3.0_dp/2.0_dp,  sqrt(3.0_dp)/2.0/)
    a2 = a*(/3.0_dp/2.0_dp, -sqrt(3.0_dp)/2.0/)

    ! x gets mapped to circumference of tube
    x = n*a1 + m*a2
    normx = norm(x)

    ! y is smallest vector perpendicular to x that can be written
    ! in form a*a1 + b*a2 where a and b are integers, i.e. it is
    ! the periodicity of nanotube along tube length
    y = (2.0_dp*n + m)*a1 - (n + 2.0_dp*m)*a2
    y = y / real(gcd(2*n + m, n + 2*m), dp)
    normy = norm(y)

    ! Rotate unit cell to align x axis with vector x
    theta = -atan(x(2)/x(1))
    rot = reshape((/ cos(theta), -sin(theta), 0.0_dp, &
                     sin(theta),  cos(theta), 0.0_dp, &
                     0.0_dp,      0.0_dp,     1.0_dp /), (/3,3/))

    unit%pos = rot .mult. unit%pos
    lattice = rot .mult. unit%lattice
    call Atoms_Set_Lattice(unit, lattice)

    ! Choose nx and ny so we can cut an x by y square out of the
    ! rotated unitsheet
    call Fit_Box_in_Cell(normx, normy, norm(unit%lattice(:,3)), unit%lattice, nx, ny, nz)
    call Supercell(unitsheet, unit, nx, ny, 1)
    call Zero_Sum(unitsheet%pos)

    ! Final lattice for unitsheet is norm(x) by norm(y)
    lattice = 0.0_dp
    lattice(1,1) = normx
    lattice(2,2) = normy
    lattice(3,3) = unitsheet%lattice(3,3)

    call Atoms_Set_Lattice(unitsheet, lattice)

    ! Form primitive cell by discarding atoms with 
    ! lattice coordinates outside range [-0.5,0.5]
    d = (/ 0.01, 0.02, 0.03 /)
    i = 1
    do
       t = unitsheet%g .mult. (unitsheet%pos(:,i) + d)
       if (any(t < -0.5_dp) .or. any(t >= 0.5_dp)) then
          call Remove_Atoms(unitsheet, i)
          i = i - 1 ! Retest since we've removed an atom
       end if
       if (i == unitsheet%N) exit
       i = i + 1
    end do

    call Supercell(sheet, unitsheet, rep_x, rep_y, 1)

    call Finalise(unit)
    call Finalise(unitsheet)
  end subroutine Graphene_Sheet


  !% Calcualte average radius of a nanotube
  function Tube_Radius(tube) result(r)
    type(Atoms), intent(in) :: tube
    real(dp) :: r

    real(dp), dimension(2,tube%N) :: pos
    real(dp) :: r_sum
    integer :: i

    ! Collapse in z direction
    pos(:,:) = tube%pos(1:2,:)

    ! Centre on (0,0)
    call Zero_Sum(pos)

    r_sum = 0.0_dp
    do i=1,tube%N
       r_sum = r_sum + sqrt(pos(1,i)**2 + pos(2,i)**2)
    end do
    r = r_sum/tube%N ! Average radius
   
  end function Tube_Radius



  !% Euclidean algorithm for greatest common divisor of 'a' and 'b'
  !% needed by 'graphene_tube'
  function gcd(ai, bi)
    integer, intent(in) :: ai, bi
    integer :: gcd

    integer :: a, b, c

    a = ai
    b = bi
    do while (b /= 0)
       c = b
       b = mod(a,b)
       a = c
    end do
    gcd = a

  end function gcd

  !% Construct a $(n,m)$ nanotube with lattice parameter 'a' and 'nz'
  !% unit cells along the tube length. Returns the radius of the tube.
  function Graphene_Tube(tube, a, n, m, nz) result(r)
    type(Atoms), intent(out) :: tube
    real(dp), intent(in) :: a
    integer, intent(in) :: n, m, nz
    real(dp) :: r

    type(Atoms) :: unit, sheet, unittube
    real(dp) :: a1(2), a2(2), x(2), lattice(3,3), normx
    integer :: i
    
    real, parameter :: NANOTUBE_VACCUUM = 10.0_dp

    a1 = a*(/3.0_dp/2.0_dp,  sqrt(3.0_dp)/2.0/)
    a2 = a*(/3.0_dp/2.0_dp, -sqrt(3.0_dp)/2.0/)

    call Graphene_Sheet(sheet, a, n, m, 1, 1)

    ! x gets mapped to circumference of tube
    x = n*a1 + m*a2
    normx = norm(x)

    ! Circumference of tube is norm(x)
    r = normx/(2.0_dp*pi)

    ! Tube lattice 
    lattice = 0.0_dp
    lattice(1,1) = 2*r + NANOTUBE_VACCUUM
    lattice(2,2) = lattice(1,1)
    lattice(3,3) = sheet%lattice(2,2)
    
    call Initialise(unittube, sheet%N, lattice)
    unittube%Z = sheet%Z

    if (has_property(unittube,'mass')) &
         forall(i=1:unittube%N) unittube%mass(i) = ElementMass(unittube%Z(i))

    ! Roll up sheet to form tube
    do i = 1,unittube%N
       unittube%pos(1,i) = r*cos(2.0_dp*pi*sheet%pos(1,i)/normx)
       unittube%pos(2,i) = r*sin(2.0_dp*pi*sheet%pos(1,i)/normx)
       unittube%pos(3,i) = sheet%pos(2,i)
    end do

    call Supercell(tube, unittube, 1, 1, nz)
    call Finalise(unit)
    call Finalise(sheet)
    call Finalise(unittube)
    
  end function Graphene_Tube

  !
  !% Return an atoms object containing one TIP3P water molecule in a box 
  !% giving the correct density at 300K
  !
  function water()
    
    type(Atoms) :: water

    call initialise(water,3,make_lattice(3.129058333_dp))

    water%Z(1) = 8
    water%Z(2) = 1
    water%Z(3) = 1

    water%pos(:,1) = (/0.0_dp,0.0_dp,0.0_dp/)
    water%pos(:,2) = (/0.9572_dp,0.0_dp,0.0_dp/)
    water%pos(:,3) = (/-0.2399872084_dp,0.9266272065_dp,0.0_dp/)
    
  end function water

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! supercell(lotsofatoms, atoms, n1, n2, n3)
  !
  !% Replicates the unit cell 'n1*n2*n3' times along the lattice vectors.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine supercell(aa, a, n1, n2, n3)
    type(Atoms), intent(out)::aa  !% Output (big) cell
    type(Atoms), intent(in)::a    !% Input cell
    integer, intent(in)::n1, n2, n3
    real(dp)::lattice(3,3), p(3)
    integer::i,j,k,n,nn
    type(Table) :: big_data

    call allocate(big_data,a%data%intsize,a%data%realsize,&
         a%data%strsize, a%data%logicalsize, a%N*n1*n2*n3)

    ! Replicate atomic data n1*n2*n3 times
    do i=1,n1*n2*n3
       call append(big_data,a%data)
    end do

    lattice(:,1) = a%lattice(:,1)*n1
    lattice(:,2) = a%lattice(:,2)*n2
    lattice(:,3) = a%lattice(:,3)*n3
    call atoms_initialise(aa, a%N*n1*n2*n3, lattice, data=big_data, properties=a%properties)
    aa%use_uniform_cutoff = a%use_uniform_cutoff
    aa%cutoff = a%cutoff

    do i = 0,n1-1
       do j = 0,n2-1
          do k = 0,n3-1
             p = a%lattice .mult. (/i,j,k/)
             do n = 1,a%N
                nn = ((i*n2+j)*n3+k)*a%n+n 
                ! overwrite position with shifted pos
                aa%pos(:,nn) = a%pos(:,n)+p
             end do
          end do
       end do
    end do

    call finalise(big_data)

  end subroutine supercell


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! diamond(myatoms, a)
  !
  !% Creates an 8-atom diamond-structure with cubic lattice constant of 'a'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine diamond(myatoms, a, Z)
    type(Atoms), intent(out)      :: myatoms
    real(dp), intent(in)          :: a
    integer, intent(in), optional :: Z(:)

    integer :: i

    call atoms_initialise(myatoms, 8, &
         reshape((/a,0.0_dp,0.0_dp,0.0_dp,a,0.0_dp,0.0_dp,0.0_dp,a/), (/3,3/)))
    
    myatoms%pos(:,1) = a*(/0.00_dp, 0.00_dp, 0.00_dp/)
    myatoms%pos(:,2) = a*(/0.25_dp, 0.25_dp, 0.25_dp/)
    myatoms%pos(:,3) = a*(/0.50_dp, 0.50_dp, 0.00_dp/)
    myatoms%pos(:,4) = a*(/0.75_dp, 0.75_dp, 0.25_dp/)
    myatoms%pos(:,5) = a*(/0.50_dp, 0.00_dp, 0.50_dp/)
    myatoms%pos(:,6) = a*(/0.75_dp, 0.25_dp, 0.75_dp/)
    myatoms%pos(:,7) = a*(/0.00_dp, 0.50_dp, 0.50_dp/)
    myatoms%pos(:,8) = a*(/0.25_dp, 0.75_dp, 0.75_dp/)

    if (present(Z)) then
      if (size(Z) == 1) then
	myatoms%Z = Z(1)
        if (has_property(myatoms, 'mass')) &
             myatoms%mass = ElementMass(Z(1))
      else if (size(Z) == 2) then
	myatoms%Z(1:7:2) = Z(1)
        if (has_property(myatoms, 'mass')) &
             myatoms%mass(1:7:2) = ElementMass(Z(1))

	myatoms%Z(2:8:2) = Z(2)
        if (has_property(myatoms, 'mass')) &
             myatoms%mass(2:8:2) = ElementMass(Z(2))
      else
	call system_abort("diamond got passed Z with size="//size(Z)//" which isn't 1 or 2")
      endif

      do i=1,myatoms%N
         myatoms%species(i) = ElementName(myatoms%Z(i))
      end do

    endif
  end subroutine diamond

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! fcc(myatoms, a, Z)
  !
  !% Creates a 4-atom fcc-structure with cubic lattice constant of 'a'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine fcc(myatoms, a, Z)
    type(Atoms), intent(out)      :: myatoms
    real(dp)                      :: a
    integer, intent(in), optional :: Z

    call atoms_initialise(myatoms, 4, &
         reshape((/a,0.0_dp,0.0_dp,0.0_dp,a,0.0_dp,0.0_dp,0.0_dp,a/), (/3,3/)))
    
    myatoms%pos(:,1) = a*(/0.00_dp, 0.00_dp, 0.00_dp/)
    myatoms%pos(:,2) = a*(/0.50_dp, 0.50_dp, 0.00_dp/)
    myatoms%pos(:,3) = a*(/0.50_dp, 0.00_dp, 0.50_dp/)
    myatoms%pos(:,4) = a*(/0.00_dp, 0.50_dp, 0.50_dp/)

    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine fcc

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! bcc(myatoms, a, Z)
  !
  !% Creates a 2-atom bcc-structure with cubic lattice constant of 'a'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 
  subroutine bcc(myatoms, a, Z)
    type(Atoms), intent(out) :: myatoms
    real(dp)                 :: a
    integer, intent(in), optional :: Z

    call atoms_initialise(myatoms, 2, &
         reshape((/a,0.0_dp,0.0_dp,0.0_dp,a,0.0_dp,0.0_dp,0.0_dp,a/), (/3,3/)))
    
    myatoms%pos(:,1) = a*(/0.00_dp, 0.00_dp, 0.00_dp/)
    myatoms%pos(:,2) = a*(/0.50_dp, 0.50_dp, 0.50_dp/)
 
    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine bcc

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! a15(myatoms, a, Z)
  !
  !% Creates an 8-atom a15-structure with cubic lattice constant of 'a'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 
  subroutine a15(myatoms, a, Z)
    type(Atoms), intent(out)      :: myatoms
    real(dp)                      :: a
    integer, intent(in), optional :: Z

    call atoms_initialise(myatoms, 8, &
         reshape((/a,0.0_dp,0.0_dp,0.0_dp,a,0.0_dp,0.0_dp,0.0_dp,a/), (/3,3/)))

    myatoms%pos(:,1) = a*(/0.00_dp, 0.00_dp, 0.00_dp/)
    myatoms%pos(:,2) = a*(/0.50_dp, 0.50_dp, 0.50_dp/)
    myatoms%pos(:,3) = a*(/0.25_dp, 0.50_dp, 0.00_dp/)
    myatoms%pos(:,4) = a*(/0.75_dp, 0.50_dp, 0.00_dp/)
    myatoms%pos(:,5) = a*(/0.00_dp, 0.25_dp, 0.50_dp/)
    myatoms%pos(:,6) = a*(/0.00_dp, 0.75_dp, 0.50_dp/)
    myatoms%pos(:,7) = a*(/0.50_dp, 0.00_dp, 0.25_dp/)
    myatoms%pos(:,8) = a*(/0.50_dp, 0.00_dp, 0.75_dp/)

    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine a15

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! graphite(myatoms, a, c, Z)
  !
  !% Creates a 4-atom graphite lattice with lattice constants of 'a' and 'c'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine graphite(myatoms, a, c, Z)
    type(Atoms), intent(out)      :: myatoms
    real(dp), intent(in)          :: a, c
    integer, intent(in), optional :: Z

    call atoms_initialise(myatoms, 4, &
         reshape( (/0.5_dp*a,-0.5_dp*sqrt(3.0_dp)*a, 0.0_dp, &
                  & 0.5_dp*a, 0.5_dp*sqrt(3.0_dp)*a, 0.0_dp, &
                  & 0.0_dp,   0.0_dp,                c/),(/3,3/)))
    
    myatoms%pos(:,1) = matmul(myatoms%lattice, (/ 0.0_dp,        0.0_dp,        0.25_dp /))
    myatoms%pos(:,2) = matmul(myatoms%lattice, (/ 0.0_dp,        0.0_dp,        0.75_dp /))
    myatoms%pos(:,3) = matmul(myatoms%lattice, (/ 1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.25_dp /))
    myatoms%pos(:,4) = matmul(myatoms%lattice, (/ 2.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.75_dp /))

    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine graphite

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Given a lattice type ('P','I','F','A','B','C') and a motif,
  !% creates an atoms object which consists of the motif
  !% applied to each lattice point in one unit cell.
  !% This object can then be supercelled.
  !%
  !% e.g. The following code creates the YBCO superconductor
  !%      structure, and allows the oxygen positions to be changed
  !%      by altering delta_O:
  !%
  !%>   type(atoms)              :: ybco
  !%>   type(atoms)              :: big_ybco
  !%>   type(table)              :: motif
  !%>   real(dp), dimension(3,3) :: lattice
  !%>   real(dp)                 :: delta_O = 0.03_dp
  !%>
  !%>   lattice = make_lattice(3.8227_dp,3.8827_dp,11.6802_dp)
  !%>
  !%>   call allocate(motif,1,3,0,0,13)
  !%>
  !%>   call append(motif,Atomic_Number('Y'), &
  !%>               (/real(0.5,dp),real(0.5,dp),real(0.5,dp)/))
  !%>   call append(motif,Atomic_Number('Ba'),
  !%>               (/real(0.5,dp),real(0.5,dp),real((1.0/6.0),dp)/))
  !%>   call append(motif,Atomic_Number('Ba'),
  !%>               (/real(0.5,dp),real(0.5,dp),real((5.0/6.0),dp)/))
  !%>   call append(motif,Atomic_Number('Cu'),
  !%>               (/real(0.0,dp),real(0.0,dp),real(0.0,dp)/))
  !%>   call append(motif,Atomic_Number('Cu'),
  !%>               (/real(0.0,dp),real(0.0,dp),real((1.0/3.0),dp)/))
  !%>   call append(motif,Atomic_Number('Cu'),
  !%>               (/real(0.0,dp),real(0.0,dp),real((2.0/3.0),dp)/))
  !%>   call append(motif,Atomic_Number('O'),
  !%>               (/real(0.0,dp),real(0.0,dp),real((1.0/6.0),dp)/))
  !%>   call append(motif,Atomic_Number('O'),
  !%>               (/real(0.0,dp),real(0.5,dp),real(0.0,dp)/))
  !%>   call append(motif,Atomic_Number('O'), &
  !%>               (/real(0.5,dp),real(0.0,dp),real((1.0/3.0),dp)+delta_O/))
  !%>   call append(motif,Atomic_Number('O'), &
  !%>               (/real(0.0,dp),real(0.5,dp),real((1.0/3.0),dp)+delta_O/))
  !%>   call append(motif,Atomic_Number('O'), &
  !%>               (/real(0.5,dp),real(0.0,dp),real((2.0/3.0),dp)-delta_O/))
  !%>   call append(motif,Atomic_Number('O'), &
  !%>               (/real(0.0,dp),real(0.5,dp),real((2.0/3.0),dp)-delta_O/))
  !%>   call append(motif,Atomic_Number('O'), &
  !%>               (/real(0.0,dp),real(0.0,dp),real((5.0/6.0),dp)/))
  !%>
  !%>   ybco = make_structure(lattice,'P',motif) 
  !%>
  !%>   call supercell(big_ybco,ybco,2,2,2)
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function make_structure(lattice, type, motif) result(structure)

    real(dp), dimension(3,3), intent(in) :: lattice
    type(table),              intent(in) :: motif
    character(1),             intent(in) :: type     ! P,I,F,A,B or C
    type(Atoms)                          :: structure
    !local variables
    integer                              :: i, j, n
    type(table)                          :: points

    !Check the motif
    if (motif%N < 1) call system_abort('Make_Structure: Must have at least one atom in motif')
    if (motif%intsize /= 1 .or. motif%realsize /= 3) &
         call system_abort('Make_Structure: Motif must have 1 int and 3 reals')

    !Set up the lattice points
    call allocate(points,0,3,0,0,4)
    call append(points, realpart=(/0.0_dp,0.0_dp,0.0_dp/))

    select case(type)
       case('P')
          !Primitive. No extra points
       case('I')
          !Body centred
          call append(points, realpart=(/0.5_dp,0.5_dp,0.5_dp/))
       case('F')
          !Face centred
          call append(points, realpart=(/0.5_dp,0.5_dp,0.0_dp/))
          call append(points, realpart=(/0.5_dp,0.0_dp,0.5_dp/))
          call append(points, realpart=(/0.0_dp,0.5_dp,0.5_dp/))
       case('A')
          call append(points, realpart=(/0.0_dp,0.5_dp,0.5_dp/))
       case('B')
          call append(points, realpart=(/0.5_dp,0.0_dp,0.5_dp/))
       case('C')
          call append(points, realpart=(/0.5_dp,0.5_dp,0.0_dp/))
       case default
          call system_abort('Make_Structure: Unknown lattice type: '//type)
    end select

    call atoms_initialise(structure,(motif%N * points%N),lattice)

    do j = 1, points%N
       do i = 1, motif%N
          
          if (motif%int(1,i) < 1) call system_abort('Make_Structure: Atom '//i//' not recognised')

          n = (j-1)*motif%N + i

          structure%Z(n) = motif%int(1,i)
          if (has_property(structure, 'mass')) &
               structure%mass(n) = ElementMass(structure%Z(n))
          structure%pos(:,n) = lattice .mult. (motif%real(:,i) + points%real(:,j))
          
       end do
    end do
    
    call finalise(points)

  end function make_structure

  subroutine transform(at_out, at_in, t)
    type(Atoms), intent(out)::at_out  !% Output 
    type(Atoms), intent(in) ::at_in   !% Input 
    real(dp), intent(in)::t(3,3)
    integer::i

    at_out = at_in
    call set_lattice(at_out,(at_in%lattice .mult. transpose(t)))

    do i=1,at_out%N
       at_out%pos(:,i) = t .mult. at_in%pos(:,i)
    end do
  end subroutine transform


  !%
  !% Subgraph isomorphism identifier based on J.R. Ullmann, JACM 23(1) 31-42 (1976)
  !%
  !% Slight modifications are that if we include two unconnected vertices in the subgraph
  !% then the corresponding vertices in the isomorphism MUST also be disconnected.
  !% i.e. if we don\'t include a bond between two atoms it is because it really isn\'t there.
  !%
  !% Pattern matching problems are combinatorial in time required, so the routine itself
  !% has seven different escape routes; the first six try to quickly fail and prevent the 
  !% main part of the algorithm from executing at all. 
  !%
  !% 'at' is the atoms structure with connectivity data precalculated to at least first nearest neighbours
  !%
  !% 'motif' is an integer matrix describing the connectivity of the region you wish to match.
  !% it has dimension (number of atoms, max number of neighbours + 1)
  !% motif(i,1) is the atomic number of an atom
  !% motif(i,2), motif(i,3) etc. are the indices of atoms to which this atom connects in this motif
  !% or zero if there are no more neighbours.
  !%
  !% E.g. to match a water molecule we could use:
  !%
  !%> water_motif  = reshape( (/  8, 1, 1, &
  !%>                             2, 1, 1, &
  !%>                             3, 0, 0  /), (/3,3/) )
  !%
  !% or, alternatively:
  !%> water_motif2 = reshape( (/ 1, 8, 1, &
  !%>                            2, 1, 2, &
  !%>                            0, 3, 0/), (/3,3/) )
  !%
  !% and for an alpha carbon 
  !%>          O   
  !%>          |   
  !%>  N - C - C   
  !%>      | 
  !%>      H  
  !%>
  !%> c_alpha = reshape( (/ 6,6,7,8,1, &
  !%>                       2,1,1,2,1, &
  !%>                       3,4,0,0,0, &
  !%>                       5,0,0,0,0/), (/5,4/) )
  !%
  !% The routine will identify an optimum atom in the motif which it will try to find in the
  !% atoms structure before doing any further matching. The optimum atom is the one which
  !% minimises the total number of bond hops required to include all atoms in the motif
  !%
  !% 'matches' is a table containing one line for each match found in the atoms structure
  !% or for optimum atoms with indices between 'start' and 'end'. The integers in each line
  !% give the indices of the atoms, in the same order as in the motif, which consitute a single
  !% match.
  !%
  !% 'mask' allows individual atoms to be selected for searching, e.g. for preventing a water
  !% molecule from being re-identified as an OH, and then later as two hydrogens.
  !%
  subroutine find_motif(at,motif,matches,start,end,mask)

    type(atoms),       intent(in)  :: at           !% The atoms structure to search
    integer,           intent(in)  :: motif(:,:)   !% The motif to search for
    type(table),       intent(out) :: matches      !% All matches
    integer, optional, intent(in)  :: start, end   !% Start and End atomic indices for search
    logical, optional, intent(in)  :: mask(:)      !% If present only masked atoms are searched

    character(*), parameter  :: me = 'find_motif: '

    integer,     allocatable :: A(:,:),B(:,:),C(:,:),depth(:), neighbour_Z(:), Z(:), M0(:,:), M(:,:), &
                                match_indices(:), depth_real(:)
    integer                  :: i,j,k,p,q,N,max_depth,max_depth_real, my_start,my_end, opt_atom
    logical                  :: match
    type(table)              :: neighbours, core, num_species_at, num_species_motif

    integer                  :: discards(7)

    discards = 0

    !XXXXXXXXXXXXXXXX
    !X INITIALISATION
    !XXXXXXXXXXXXXXXX

    ! Check atomic numbers
    if (any(motif(:,1)<1)) call system_abort(me//'Bad atomic numbers ('//motif(:,1)//')')
    ! Check for atomic connectivity
    if (.not.at%connect%initialised) call system_abort(me//'No connectivity data present in atoms structure')    
    if (size(motif,1).eq.0) call system_abort(me//'zero motif size')

    ! Build adjacency matrix from motif
    N = size(motif,1)
    allocate(A(N,N),C(N,N),depth(N))
    A = 0
    do i = 1, N
       do p = 2, size(motif,2)
          j = motif(i,p)
          if (j>N) call system_abort(me//'Undefined atom referenced in motif ('//j//')')
          if (j>0) A(i,j) = 1
       end do
    end do

    ! Make array of atomic numbers in motif
    allocate(Z(N))
    Z = motif(:,1)

    ! Check symmetry 
    if (.not.is_symmetric(A)) then
       call print_title('Non-symmetic adjacency matrix')
       call print(A)
       call print('bad rows/columns:')
       do i = 1, size(A,1)
          if (.not.all(A(:,i)==A(i,:))) call print(i)
       end do
       call system_abort(me//'Adjacency matrix for given motif is not symmetric!')
    end if

    ! Allocate matches table: the ints are atomic indices which match, in order, the connectivity info
    call allocate(matches,N,0,0,0,1)

    ! Find the best atom (opt_atom) to identify and grow clusters around in the real structure
    call find_optimum_central_atom(A,opt_atom,depth)
    max_depth = maxval(depth)

    !XXXXXXXXXXXXXXX
    !X THE ALGORITHM
    !XXXXXXXXXXXXXXX

    !Loop over all atoms looking for candidates for the optimum atom 'opt_atom'
    call allocate(core,4,0,0,0,1)
    allocate(neighbour_Z(count(A(opt_atom,:)==1)))
    neighbour_Z = pack(Z,(A(opt_atom,:)==1))

    my_start = 1
    if (present(start)) then
       if (start < 1 .or. start > at%N) call system_abort(me//'Starting atom is out of range ('//start//')')
       my_start = start
    end if
    my_end = at%N
    if (present(end)) then
       if (end < 1 .or. end > at%N) call system_abort(me//'Ending atom is out of range ('//start//')')
       my_end = end
    end if

    ! if end < start then no matches will be returned, a la do loops

    do i = my_start, my_end

       if (present(mask)) then
          if (.not.mask(i)) cycle
       end if

       !XXXXXXXXXXXXXXXXXXXX
       !X ESCAPE ROUTE 1
       !XXXXXXXXXXXXXXXXXXXX

       ! Discard atoms which don't match the atomic number of the optimum atom

       if (at%Z(i) /= Z(opt_atom)) then
          discards(1) = discards(1) + 1
          cycle
       end if

       !---------------------



       !XXXXXXXXXXXXXXXXXXXX
       !X ESCAPE ROUTE 2
       !XXXXXXXXXXXXXXXXXXXX

       ! Check if the real atom's neighbours can be matched with those in the motif

       if (.not.atoms_compatible(at,i,A,Z,opt_atom)) then
          discards(2) = discards(2) + 1
          cycle
       end if

       !---------------------

       ! Grow a cluster around the real atom i which is  max_depth hops deep
       call wipe(core)
       call append(core,(/i,0,0,0/))
       call bfs_grow(at,core,max_depth,nneighb_only = .true., min_images_only = .true.)

       !XXXXXXXXXXXXXXXXXXXX
       !X ESCAPE ROUTE 3
       !XXXXXXXXXXXXXXXXXXXX

       ! If the number of atoms in the cluster is less than those in the motif we can't have a match
       if (core%N < N) then
          discards(3) = discards(3) + 1
          cycle
       end if

       !---------------------

       ! Count the number of each atomic species in the motif and in the cluster
       call allocate(num_species_motif,2,0,0,0,6)
       call allocate(num_species_at,2,0,0,0,6)
       do k = 1, N
          p = find_in_array(int_part(num_species_motif,1),Z(k))
          if (p == 0) then
             call append(num_species_motif,(/Z(k),1/))
          else
             num_species_motif%int(2,p) = num_species_motif%int(2,p) + 1
          end if
       end do
       do k = 1, core%N
          p = find_in_array(int_part(num_species_at,1),at%Z(core%int(1,k)))
          if (p == 0) then
             call append(num_species_at,(/at%Z(core%int(1,k)),1/))
          else
             num_species_at%int(2,p) = num_species_at%int(2,p) + 1
          end if
       end do

       !XXXXXXXXXXXXXXXXXXXX
       !X ESCAPE ROUTE 4
       !XXXXXXXXXXXXXXXXXXXX

       !Compare the tables: if there are more atoms of one type in motif than core then a fit is impossible
       match = .true.
       do k = 1, num_species_motif%N
          p = find_in_array(int_part(num_species_at,1),num_species_motif%int(1,k))
          if (p == 0) then
             match = .false.
             exit
          else if (num_species_at%int(2,p) < num_species_motif%int(2,k)) then
             match = .false.
             exit
          end if
       end do
       call finalise(num_species_motif)
       call finalise(num_species_at)
       if (.not.match) then
          discards(4) = discards(4) + 1
          cycle
       end if

       !---------------------

       ! Create adjacency matrix
       allocate(B(core%N,core%N))
       B = 0
       do p = 1, core%N
          k = core%int(1,p)
          do q = p+1, core%N
             j = core%int(1,q)
             if (distance_min_image(at,k,j) < bond_length(at%Z(k),at%Z(j))*at%nneightol) then
                B(p,q) = 1
                B(q,p) = 1
             end if
          end do
       end do

       ! Find depth of connectivity for real atoms
       allocate(depth_real(core%N))
       call find_connectivity_depth(B,1,depth_real)
       max_depth_real = maxval(depth_real)

       !XXXXXXXXXXXXXXXXXX
       !X ESCAPE ROUTE 5
       !XXXXXXXXXXXXXXXXXX

       ! If there are atoms which aren't deep enough in the real structure then a match is impossible
       if (max_depth > max_depth_real) then
          deallocate(B,depth_real)
          discards(5) = discards(5) + 1
          cycle
       end if

       !---------------------

    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !X THE MAIN PART OF THE ALGORITHM
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

       ! Create possibility matrix M0, with the optimum atom already matched

       ! ***** THINGS CAN BE SPED UP CONSIDERABLY BY LIMITING THE NUMBER OF 1'S IN M0 *****

       allocate(M0(N,core%N))
       M0 = 0
       M0(opt_atom,1) = 1
       do q = 2, core%N
          do p = 1, N
             if (p==opt_atom) cycle
             if (depth_real(q) == depth(p)) then
                if (atoms_compatible(at,core%int(1,q),A,Z,p)) M0(p,q) = 1
             end if
          end do
       end do

       deallocate(depth_real)

       !XXXXXXXXXXXXXXXXXX
       !X ESCAPE ROUTE 6
       !XXXXXXXXXXXXXXXXXX

       ! Check to see if any atoms have no possibilities
       if (any(sum(M0,dim=2)==0)) then
          deallocate(M0,B)
          discards(6) = discards(6) + 1
          cycle
       end if

       !---------------------

       ! Loop over all permuations
       allocate(M(N,core%N))
       call first_trial_matrix(M0,M)
       match = .false.
       do 

          ! For each trial matrix create the permuted adjacency matrix
          C = matmul(M,transpose(matmul(M,B))) ! ***** THIS LINE IS ANOTHER CANDIDATE FOR OPTIMISATION *****
                                               ! Use sparse matrices maybe?

          if (all(C==A)) then
             match = .true.
             exit
          end if

          call next_trial_matrix(M0,M)

          !XXXXXXXXXXXXXXXXXX
          !X ESCAPE ROUTE 7
          !XXXXXXXXXXXXXXXXXX

          ! If next_trial_matrix deletes the only atom we know to be fitted then all
          ! permutations have been exhausted
          if (M(opt_atom,1)==0) then
             discards(7) = discards(7) + 1
             exit
          end if
          !---------------------

       end do

       if (match) then
          !decode trial matrix into atomic indices
          allocate(match_indices(N))
          do j = 1, N
             do k = 1, core%N
                if (M(j,k)==1) then
                   match_indices(j) = core%int(1,k)
                   exit
                end if
             end do
          end do
          call append(matches,match_indices)
         deallocate(match_indices)
       end if

       deallocate(B,M0,M)

    end do

    deallocate(A,C,depth,Z,neighbour_Z)
    call finalise(core)
    call finalise(neighbours)

    !When debugging/improving this algorithm uncomment this line
    !as it gives an indication of how successful the predetection
    !of bad matches is

    !call print('find_motif: Discards by level: '//discards)

  end subroutine find_motif

  !% OMIT
  ! Used by find_motif to generate the first permutation matrix
  subroutine first_trial_matrix(M0,M)

    integer, dimension(:,:), intent(in)    :: M0
    integer, dimension(:,:), intent(inout) :: M
    integer, dimension(size(M0,2))         :: taken
    integer :: A,B,i,j

    A = size(M0,1)
    B = size(M0,2)
    M = 0
    taken = 0

    do i = 1, A

       do j = 1, B

          if (M0(i,j) == 1 .and. taken(j) == 0) then
             M(i,j) = 1
             taken(j) = 1
             exit
          end if

       end do

    end do

  end subroutine first_trial_matrix

  !% OMIT
  ! Used by find_motif to generate the next permutation matrix given the current one
  subroutine next_trial_matrix(M0,M)

    integer, dimension(:,:), intent(in)    :: M0
    integer, dimension(:,:), intent(inout) :: M
    integer, dimension(size(M0,2))         :: taken

    integer :: A,B,i,j
    logical :: found, placed, l1, l2

    A = size(M,1)
    B = size(M,2)

    i = A

    taken = sum(M,dim=1)
    do while (i<=A .and. i>0)

       !try to move the 1 in the ith row along, or place it if it isn't there

       !find the 1
       found = .false.
       j = 1
       do while(j<=B)
          if (M(i,j)==1) then
             found = .true.
             exit
          else
             j = j + 1
          end if
       end do
       if (found) then
          M(i,j) = 0
          j = j + 1
       else
          j = 1
       end if

       !place it
       placed = .false.
       taken = sum(M,dim=1)
       do while (j<=B)
          l1 = M0(i,j) == 1
          l2 = taken(j) == 0
          if (l1 .and. l2) then
             M(i,j) = 1
             taken(j) = 1
             placed = .true.
             exit
          end if
          j = j + 1
       end do
       if (.not.placed) then
          i = i - 1
       else
          i = i + 1
       end if

    end do

  end subroutine next_trial_matrix

  !% OMIT
  ! Used by find_motif
  subroutine find_optimum_central_atom(A,i,depth)

    integer, dimension(:,:),       intent(in)  :: A ! The adjacency matrix
    integer,                       intent(out) :: i ! The optimum central atom
    integer, dimension(size(A,1)), intent(out) :: depth ! The connectivity depth using this atom

    integer                       :: n, max_depth, new_max_depth
    integer, dimension(size(A,1)) :: new_depth

    ! Starting with each atom in turn, calculate the depth of the connectivity
    ! The best is the one with the smallest maximum depth

    max_depth = huge(1)

    do n = 1, size(A,1)

       call find_connectivity_depth(A,n,new_depth)
       new_max_depth = maxval(new_depth)

       if (new_max_depth < max_depth) then
          max_depth = new_max_depth
          depth = new_depth
          i = n
       end if

    end do

  end subroutine find_optimum_central_atom

  !% OMIT
  ! Used by find_motif
  subroutine find_connectivity_depth(A,i,depth)

    integer, dimension(:,:),       intent(in)  :: A ! The adjacency matrix
    integer,                       intent(in)  :: i ! The central atom
    integer, dimension(size(A,1)), intent(out) :: depth ! The connectivity depth using this atom

    integer                       :: p, q, j
    logical                       :: atom_marked

    depth = -1
    depth(i) = 0
    atom_marked = .true.
    p = -1

    do while(atom_marked)

       atom_marked = .false.
       p = p + 1

       do q = 1, size(A,1)

          if (depth(q) == p) then !...mark its neighbours with p+1 if unmarked

             do j = 1, size(A,1)

                if (A(q,j) == 1 .and. depth(j) == -1) then
                   depth(j) = p + 1
                   atom_marked = .true.
                end if

             end do

          end if

       end do

    end do

    if (any(depth == -1)) call system_abort('find_conectivity_depth: Adjacency matrix contains dijoint regions')

  end subroutine find_connectivity_depth

  !% OMIT
  ! Used by find_motif
  function atoms_compatible(at,i,A,Z,j)

    type(atoms),             intent(in) :: at !the atoms structure
    integer,                 intent(in) :: i  !the atom in the atoms structure we are testing
    integer, dimension(:,:), intent(in) :: A  !the adjacency matrix of the motif
    integer, dimension(:),   intent(in) :: Z  !the atomic number of the motif
    integer,                 intent(in) :: j  !the atom in the motif we are testing
    logical                             :: atoms_compatible

    type(table)                        :: core, neighbours
    logical, allocatable, dimension(:) :: neighbour_used
    integer, allocatable, dimension(:) :: neighbour_Z
    integer                            :: p,q
    logical                            :: match, found

    ! First check the atomic numbers of the two atoms to test
    if (at%Z(i) /= Z(j)) then
       atoms_compatible = .false.
       return
    end if

    ! find nearest neighbours of real atom
    call allocate(core,4,0,0,0,1)
    call append(core,(/i,0,0,0/))
    call bfs_step(at,core,neighbours,nneighb_only=.true.,min_images_only=.true.)

    allocate(neighbour_used(neighbours%N))
    neighbour_used = .false.

    ! find the atomic numbers of the neighbours of atom j in the motif
    allocate(neighbour_Z(count(A(j,:)==1)))
    neighbour_Z = pack(Z,(A(j,:)==1))

    ! see if a tentative match between the neighbours of atom i in at and atom j in
    ! the adjacency matrix is possible, just using atomic numbers

    match = .true.

    do p = 1, size(neighbour_Z)

       found = .false.

       do q = 1, neighbours%N

          if (neighbour_used(q)) cycle

          if (at%Z(neighbours%int(1,q)) == neighbour_Z(p)) then
             neighbour_used(q) = .true.
             found = .true.
             exit
          end if
       end do

       if (.not.found) then
          match = .false.
          exit
       end if

    end do

    atoms_compatible = match

    deallocate(neighbour_used,neighbour_Z)
    call finalise(core)
    call finalise(neighbours)

  end function atoms_compatible


end module structures_module
