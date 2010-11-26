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
!X  Structures module
!X
!%  A collection of utility functions that generate useful Atoms structures
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module  structures_module
  use system_module
  use linearalgebra_module
  use atoms_module
  use cinoutput_module
  use paramreader_module
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
         0.0_dp, 0.0_dp,20.0_dp/), (/3, 3/)), scale_positions=.false.)
         
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
         0.0_dp, 0.0_dp,20.0_dp/), (/3, 3/)), scale_positions=.false.)

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

    call print("disloc_noam: r1_hat " // r1_hat, PRINT_VERBOSE)
    call print("disloc_noam: r2_hat " // r2_hat, PRINT_VERBOSE)
    call print("disloc_noam: l_hat " // l_hat, PRINT_VERBOSE)

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
      call print("atom " // i_at // " pos " // at%pos(:,i_at) // " delta_p " // delta_p // " theta " // theta, PRINT_ANAL)
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
	      call print("WARNING: disloc_noam found atom too close to itself", PRINT_ALWAYS)
	    else
	      atoms_remove(j_at) = .true.
	    endif
	  endif
	else
	  if (r < 0.5_dp*bond_length(at%Z(i_at), at%Z(j_at))) then
	    if (i_at == j_at) then
	      call print("WARNING: disloc_noam found atom too close to itself", PRINT_ALWAYS)
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

 subroutine unit_slab(myatoms, axes, a, atnum, lat_type, c, u, x, y, z)
   type(Atoms), intent(out) :: myatoms
   real(dp), intent(in), dimension(3,3) :: axes
   real(dp), intent(in) :: a ! Lattice vector
   real(dp), intent(in), optional :: c, u, x, y, z ! Lattice vector
   integer, intent(in), optional :: atnum(:) ! atomic numbers
   character(len=*), intent(in), optional :: lat_type  ! lattice type (diamond, bcc, fcc)

   integer :: i, Nrep(3) ! Number of repeats in x,y and z
   real(dp), dimension(3) :: a1, a2, a3, t, d
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
   elseif(trim(my_lat_type) .eq. 'alpha_quartz') then
      if (.not. present(c) .or. .not. present(u) .or. .not. present(x) &
           .or. .not. present(y) .or. .not. present(z)) call system_abort('unit_slab: alpha_quartz missing c, u, x, y or z')
      call alpha_quartz_cubic(at, a, c, u, x, y, z)
   elseif(trim(my_lat_type) .eq. 'anatase') then
      if (.not. present(c) .or. .not. present(u) ) call system_abort('unit_slab: anatase missing c, u')
      call anatase(at, a, c, u)
   elseif(trim(my_lat_type) .eq. 'rutile') then
      if (.not. present(c) .or. .not. present(u) ) call system_abort('unit_slab: rutile missing c, u')
      call rutile(at, a, c, u)
   else
      call system_abort('unit_slab: unknown lattice type '//my_lat_type)
   endif
   a1 = axes(:,1);   a2 = axes(:,2);    a3 = axes(:,3)
   rot(1,:) = a1/norm(a1);  rot(2,:) = a2/norm(a2);  rot(3,:) = a3/norm(a3)

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
   call set_lattice(myatoms, myatoms%lattice, scale_positions=.false.)

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

 subroutine slab_width_height_nz(myslab, axes, a, width, height, nz, atnum, lat_type, c, u, x, y, z, even_nx, even_ny)
   type(Atoms), intent(out) :: myslab
   real(dp), intent(in), dimension(3,3) :: axes
   real(dp), intent(in) :: a, width, height
   real(dp), intent(IN), optional :: c, u, x, y, z ! Lattice vector
   integer, intent(in) :: nz ! Number of layers
   integer, intent(in), optional :: atnum(:) ! atomic numbers to use
   character(len=*),   optional  ::  lat_type 
   logical, optional, intent(in) :: even_nx, even_ny ! round nx or ny to even numbers

   type(Atoms) :: unit, layer
   integer nx, ny
   logical :: do_even_nx, do_even_ny

   ! even_nx and even_ny default to true for lat_type="alpha_quartz", otherwise they default to false
   if (present(lat_type)) then
      if (trim(lat_type) == 'alpha_quartz') then
         do_even_nx = optional_default(.true., even_nx)
         do_even_ny = optional_default(.true., even_ny)
      else
         do_even_nx = optional_default(.false., even_nx)
         do_even_ny = optional_default(.false., even_ny)
      end if
   else
      do_even_nx = optional_default(.false., even_nx)
      do_even_ny = optional_default(.false., even_ny)
   end if

   call unit_slab(unit, axes, a, atnum, lat_type, c, u, x, y, z)

   nx = int(floor(width/unit%lattice(1,1)))
   ny = int(floor(height/unit%lattice(2,2)))

   if (do_even_nx .and. mod(nx, 2) == 1) nx = nx - 1
   if (do_even_ny .and. mod(ny, 2) == 1) ny = ny - 1

   if (.not. (nx > 0 .and. ny > 0 .and. nz > 0)) then
      write(line,'(a, i0, i0, i0)') 'Error in slab: nx,ny,nz = ', nx, ny, nz
      call system_abort(line)
   end if

   call print('slab_width_height_nz: nx='//nx//' ny='//ny//' nz='//nz, PRINT_VERBOSE)

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

 subroutine slab_nx_ny_nz(myslab, axes, a, nx, ny, nz, atnum, lat_type, c, u, x, y, z)
   type(Atoms), intent(out) :: myslab
   real(dp), intent(in), dimension(3,3) :: axes
   real(dp), intent(in) :: a
   real(dp), intent(IN), optional :: c, u, x, y, z
   integer, intent(in) :: nx, ny, nz
   integer, intent(in), optional :: atnum(:)
   character(len=*), intent(inout), optional ::   lat_type

   type(Atoms) :: unit, layer

   call unit_slab(unit, axes, a, atnum, lat_type, c, u, x, y, z)
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

    call Initialise(cube, 4, &
         reshape((/ 3.0_dp*a, 0.0_dp,   0.0_dp, &
                    0.0_dp,   sqrt(3.0_dp)*a, 0.0_dp, &
                    0.0_dp,   0.0_dp,   10.0_dp/), (/3, 3/)))

    call set_atoms(cube, 6)
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
    call Set_Lattice(unit, lattice, scale_positions=.false.)

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
    call Set_Lattice(tmp_slab, lattice, scale_positions=.false.)

    ! Form primitive cell by discarding atoms with 
    ! lattice coordinates outside range [-0.5,0.5]
    d = (/ 0.01, 0.02, 0.03 /)
    call Table_Allocate(keep_list, 1, 0, 0, 0, tmp_slab%N)

    do i=1,tmp_slab%N
       t = tmp_slab%g .mult. (tmp_slab%pos(:,i) + d)
       if (.not. (any(t < -0.5_dp) .or. any(t >= 0.5_dp))) &
            call Append(keep_list, i)
    end do

    call Initialise(slab, keep_list%N, tmp_slab%lattice)
    slab%Z = tmp_slab%Z(keep_list%int(1,1:keep_list%N))
    slab%species = tmp_slab%species(:,keep_list%int(1,1:keep_list%N))

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
    call Set_Lattice(unit, lattice, scale_positions=.false.)

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

    call Set_Lattice(unitsheet, lattice, scale_positions=.false.)

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
    unittube%species = sheet%species

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

    water%pos(:,1) = (/0.0_dp,0.0_dp,0.0_dp/)
    water%pos(:,2) = (/0.9572_dp,0.0_dp,0.0_dp/)
    water%pos(:,3) = (/-0.2399872084_dp,0.9266272065_dp,0.0_dp/)

    call set_atoms(water, (/8, 1, 1/))
    
  end function water

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! supercell(lotsofatoms, atoms, n1, n2, n3)
  !
  !% Replicates the unit cell 'n1*n2*n3' times along the lattice vectors.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine supercell(aa, a, n1, n2, n3, error)
    type(Atoms), intent(out)::aa  !% Output (big) cell
    type(Atoms), intent(in)::a    !% Input cell
    integer, intent(in)::n1, n2, n3
    integer, intent(out), optional :: error

    real(dp)::lattice(3,3), p(3)
    integer::i,j,k,n,nn

    integer, pointer :: a_int_ptr(:), a_int2_ptr(:,:), aa_int_ptr(:), aa_int2_ptr(:,:)
    real(dp), pointer :: a_real_ptr(:), a_real2_ptr(:,:), aa_real_ptr(:), aa_real2_ptr(:,:)
    logical, pointer :: a_logical_ptr(:), aa_logical_ptr(:)
    character, pointer :: a_char_ptr(:,:), aa_char_ptr(:,:)

    INIT_ERROR(error)
    lattice(:,1) = a%lattice(:,1)*n1
    lattice(:,2) = a%lattice(:,2)*n2
    lattice(:,3) = a%lattice(:,3)*n3
    call initialise(aa, a%N*n1*n2*n3, lattice)
    if (a%use_uniform_cutoff) then
       call set_cutoff(aa, a%cutoff)
    else
       call set_cutoff_factor(aa, a%cutoff)
    end if

    do n=1,a%properties%n
       select case(a%properties%entries(n)%type)

       case(T_INTEGER_A)
          if (.not. assign_pointer(a, string(a%properties%keys(n)), a_int_ptr)) then
             RAISE_ERROR('supercell: cannot assign pointer to property '//a%properties%keys(n), error)
          end if
          call add_property(aa, string(a%properties%keys(n)), 0, ptr=aa_int_ptr, error=error, overwrite=.true.)
          PASS_ERROR(error)
          do i=0,n1-1
             do j=0,n2-1
                do k=0,n3-1
                   aa_int_ptr(((i*n2+j)*n3+k)*a%n+1:((i*n2+j)*n3+k)*a%n+a%n) = a_int_ptr(:)
                end do
             end do
          end do

       case(T_REAL_A)
          if (.not. assign_pointer(a, string(a%properties%keys(n)), a_real_ptr)) then
             RAISE_ERROR('supercell: cannot assign pointer to property '//a%properties%keys(n), error)
          end if
          call add_property(aa, string(a%properties%keys(n)), 0.0_dp, ptr=aa_real_ptr, error=error, overwrite=.true.)
          PASS_ERROR(error)
          do i=0,n1-1
             do j=0,n2-1
                do k=0,n3-1
                   aa_real_ptr(((i*n2+j)*n3+k)*a%n+1:((i*n2+j)*n3+k)*a%n+a%n) = a_real_ptr(:)
                end do
             end do
          end do

       case(T_LOGICAL_A)
          if (.not. assign_pointer(a, string(a%properties%keys(n)), a_logical_ptr)) then
             RAISE_ERROR('supercell: cannot assign pointer to property '//a%properties%keys(n), error)
          end if
          call add_property(aa, string(a%properties%keys(n)), .false., ptr=aa_logical_ptr, error=error, overwrite=.true.)
          PASS_ERROR(error)
          do i=0,n1-1
             do j=0,n2-1
                do k=0,n3-1
                   aa_logical_ptr(((i*n2+j)*n3+k)*a%n+1:((i*n2+j)*n3+k)*a%n+a%n) = a_logical_ptr(:)
                end do
             end do
          end do

       case(T_INTEGER_A2)
          if (.not. assign_pointer(a, string(a%properties%keys(n)), a_int2_ptr)) then
             RAISE_ERROR('supercell: cannot assign pointer to property '//a%properties%keys(n), error)
          end if
          call add_property(aa, string(a%properties%keys(n)), 0, n_cols=size(a_int2_ptr,1), &
               ptr2=aa_int2_ptr, error=error, overwrite=.true.)
          PASS_ERROR(error)
          do i=0,n1-1
             do j=0,n2-1
                do k=0,n3-1
                   aa_int2_ptr(:,((i*n2+j)*n3+k)*a%n+1:((i*n2+j)*n3+k)*a%n+a%n) = a_int2_ptr(:,:)
                end do
             end do
          end do

       case(T_REAL_A2)
          if (.not. assign_pointer(a, string(a%properties%keys(n)), a_real2_ptr)) then
             RAISE_ERROR('supercell: cannot assign pointer to property '//a%properties%keys(n), error)
          end if
          call add_property(aa, string(a%properties%keys(n)), 0.0_dp, n_cols=size(a_real2_ptr,1), &
               ptr2=aa_real2_ptr, error=error, overwrite=.true.)
          PASS_ERROR(error)
          do i=0,n1-1
             do j=0,n2-1
                do k=0,n3-1
                   aa_real2_ptr(:,((i*n2+j)*n3+k)*a%n+1:((i*n2+j)*n3+k)*a%n+a%n) = a_real2_ptr(:,:)
                end do
             end do
          end do

       case(T_CHAR_A)
          if (.not. assign_pointer(a, string(a%properties%keys(n)), a_char_ptr)) then
             RAISE_ERROR('supercell: cannot assign pointer to property '//a%properties%keys(n), error)
          end if
          call add_property(aa, string(a%properties%keys(n)), repeat(' ', TABLE_STRING_LENGTH), ptr=aa_char_ptr, error=error, overwrite=.true.)
          PASS_ERROR(error)
          do i=0,n1-1
             do j=0,n2-1
                do k=0,n3-1
                   aa_char_ptr(:,((i*n2+j)*n3+k)*a%n+1:((i*n2+j)*n3+k)*a%n+a%n) = a_char_ptr(:,:)
                end do
             end do
          end do

       case default
          RAISE_ERROR('supercell: bad property type '//a%properties%entries(i)%type, error)
       end select
    end do

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

    call initialise(myatoms, 8, &
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

      call set_atoms(myatoms,myatoms%Z)
    endif
  end subroutine diamond

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! diamond2(myatoms, a, Z)
  !
  !% Creates a 2-atom diamond-structure with cubic lattice constant of 'a'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine diamond2(myatoms, a, Z1, Z2)
    type(Atoms), intent(out)      :: myatoms
    real(dp)                      :: a
    integer, intent(in), optional :: Z1, Z2
    integer :: my_Z1, my_Z2

    my_Z1 = optional_default(6, Z1)
    my_Z2 = optional_default(6, Z2)

    call initialise(myatoms, 2, &
         reshape((/a/2,a/2,0.0_dp,0.0_dp,a/2,a/2,a/2,0.0_dp,a/2/), (/3,3/)))
    
    myatoms%pos(:,1) = a*(/0.00_dp, 0.00_dp, 0.00_dp/)
    myatoms%pos(:,2) = a*(/0.25_dp, 0.25_dp, 0.25_dp/)

    call set_atoms(myatoms,(/my_Z1,my_Z2/))

  end subroutine diamond2

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! beta_tin(myatoms, a, c, Z)
  !
  !% Creates a 2-atom beta-tin structure with lattice constants of a and c
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine beta_tin(myatoms, a, c, Z)
    type(Atoms), intent(out)      :: myatoms
    real(dp), intent(in)          :: a, c
    integer, intent(in), optional :: Z

    call initialise(myatoms, 2, &
         reshape((/        a,   0.0_dp,   0.0_dp, &
                  &   0.0_dp,        a,   0.0_dp, &
                  & 0.5_dp*a, 0.5_dp*a, 0.5_dp*c /), (/3,3/) ))
    
    myatoms%pos(:,1) = matmul(myatoms%lattice, (/ -0.125_dp, -0.375_dp,  0.25_dp /))
    myatoms%pos(:,2) = matmul(myatoms%lattice, (/  0.125_dp,  0.375_dp, -0.25_dp /))
    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine beta_tin

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! beta_tin4(myatoms, a, c, Z)
  !
  !% Creates a 4-atom beta-tin structure with lattice constants of a and c
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine beta_tin4(myatoms, a, c, Z)
    type(Atoms), intent(out)      :: myatoms
    real(dp), intent(in)          :: a, c
    integer, intent(in), optional :: Z

    call initialise(myatoms, 4, &
         reshape((/a, 0.0_dp, 0.0_dp, 0.0_dp, a, 0.0_dp, 0.0_dp, 0.0_dp, c/), (/3,3/)))
    
    myatoms%pos(:,1) = matmul(myatoms%lattice, (/ 0.00_dp,        0.00_dp,        0.00_dp /))
    myatoms%pos(:,2) = matmul(myatoms%lattice, (/ 0.50_dp,        0.50_dp,        0.50_dp /))
    myatoms%pos(:,3) = matmul(myatoms%lattice, (/ 0.00_dp,        0.50_dp,        0.25_dp /))
    myatoms%pos(:,4) = matmul(myatoms%lattice, (/ 0.50_dp,        0.00_dp,        0.75_dp /))
    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine beta_tin4

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

    call initialise(myatoms, 4, &
         reshape((/a,0.0_dp,0.0_dp,0.0_dp,a,0.0_dp,0.0_dp,0.0_dp,a/), (/3,3/)))
    
    myatoms%pos(:,1) = a*(/0.00_dp, 0.00_dp, 0.00_dp/)
    myatoms%pos(:,2) = a*(/0.50_dp, 0.50_dp, 0.00_dp/)
    myatoms%pos(:,3) = a*(/0.50_dp, 0.00_dp, 0.50_dp/)
    myatoms%pos(:,4) = a*(/0.00_dp, 0.50_dp, 0.50_dp/)

    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine fcc

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! hcp(myatoms, a, Z)
  !
  !% Creates a 2-atom hcp lattice with lattice constants of 'a'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine hcp(myatoms, a, Z)
    type(Atoms), intent(out)      :: myatoms
    real(dp), intent(in)          :: a
    integer, intent(in), optional :: Z

    call initialise(myatoms, 2, &
         reshape( (/0.5_dp*a,-0.5_dp*sqrt(3.0_dp)*a, 0.0_dp, &
                  & 0.5_dp*a, 0.5_dp*sqrt(3.0_dp)*a, 0.0_dp, &
                  & 0.0_dp,   0.0_dp,                a*sqrt(8.0_dp/3.0_dp) /),(/3,3/)))
    
    myatoms%pos(:,1) = matmul(myatoms%lattice, (/        0.0_dp,        0.0_dp, 0.0_dp /))
    myatoms%pos(:,2) = matmul(myatoms%lattice, (/ 1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.5_dp /))

    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine hcp

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! sh(myatoms, a, Z)
  !
  !% Creates a 1-atom simple hexagonal lattice with lattice constants of 'a' and 'c'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine sh(myatoms, a, c, Z)
    type(Atoms), intent(out)      :: myatoms
    real(dp), intent(in)          :: a, c
    integer, intent(in), optional :: Z

    call initialise(myatoms, 1, &
         reshape( (/0.5_dp*a,-0.5_dp*sqrt(3.0_dp)*a, 0.0_dp, &
                  & 0.5_dp*a, 0.5_dp*sqrt(3.0_dp)*a, 0.0_dp, &
                  & 0.0_dp,   0.0_dp,                c /),(/3,3/)))
    
    myatoms%pos(:,1) = matmul(myatoms%lattice, (/ 0.0_dp, 0.0_dp, 0.0_dp /))

    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine sh

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! sh2(myatoms, a, Z)
  !
  !% Creates a 2-atom simple hexagonal lattice with lattice constants of 'a' and 'c'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine sh2(myatoms, a, c, Z)
    type(Atoms), intent(out)      :: myatoms
    real(dp), intent(in)          :: a, c
    integer, intent(in), optional :: Z

    call initialise(myatoms, 2, &
         reshape( (/     a,         0.0_dp, 0.0_dp, &
                  & 0.0_dp, sqrt(3.0_dp)*a, 0.0_dp, &
                  & 0.0_dp,         0.0_dp, c /),(/3,3/)))
    
    myatoms%pos(:,1) = matmul(myatoms%lattice, (/ 0.0_dp, 0.0_dp, 0.0_dp /))
    myatoms%pos(:,2) = matmul(myatoms%lattice, (/ 0.5_dp, 0.5_dp, 0.0_dp /))

    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine sh2

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! imma(myatoms, a, b, c, u, Z)
  !
  !% Creates a 2-atom Imma cell
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine imma(myatoms, a, b, c, u, Z)
    type(Atoms), intent(out)      :: myatoms
    real(dp), intent(in)          :: a, b, c, u
    integer, intent(in), optional :: Z

    call initialise(myatoms, 2, &
         0.5_dp * reshape( (/ a,  b, -c, &
                             -a,  b,  c, &
                              a, -b,  c /),(/3,3/)))
    
    myatoms%pos(:,1) =  matmul(myatoms%lattice, (/ 0.25_dp, 0.25_dp+u, u /))
    myatoms%pos(:,2) = -matmul(myatoms%lattice, (/ 0.25_dp, 0.25_dp+u, u /))

    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine imma

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! sc(myatoms, a, Z)
  !
  !% Creates a 1-atom simple cubic lattice with lattice constants of 'a'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine sc(myatoms, a, Z)
    type(Atoms), intent(out)      :: myatoms
    real(dp), intent(in)          :: a
    integer, intent(in), optional :: Z

    call initialise(myatoms, 1, &
         reshape( (/     a, 0.0_dp, 0.0_dp, &
                  & 0.0_dp,      a, 0.0_dp, &
                  & 0.0_dp, 0.0_dp,      a /),(/3,3/)))
    
    myatoms%pos(:,1) = matmul(myatoms%lattice, (/ 0.0_dp, 0.0_dp, 0.0_dp /))

    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine sc

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

    call initialise(myatoms, 2, &
         reshape((/a,0.0_dp,0.0_dp,0.0_dp,a,0.0_dp,0.0_dp,0.0_dp,a/), (/3,3/)))
    
    myatoms%pos(:,1) = a*(/0.00_dp, 0.00_dp, 0.00_dp/)
    myatoms%pos(:,2) = a*(/0.50_dp, 0.50_dp, 0.50_dp/)
 
    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine bcc

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! bcc1(myatoms, a, Z)
  !
  !% Creates a 1-atom primitive bcc-structure with cubic lattice constant of 'a'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 
  subroutine bcc1(myatoms, a, Z)
    type(Atoms), intent(out) :: myatoms
    real(dp)                 :: a
    integer, intent(in), optional :: Z

    call initialise(myatoms, 1, &
    & 0.5_dp*a*reshape( (/1.0_dp,-1.0_dp, 1.0_dp, &
                        & 1.0_dp, 1.0_dp,-1.0_dp, &
                        &-1.0_dp, 1.0_dp, 1.0_dp/),(/3,3/)))
    
    myatoms%pos(:,1) = (/0.00_dp, 0.00_dp, 0.00_dp/)
 
    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine bcc1

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! zincblende2(myatoms, a, Z1, Z2)
  !
  !% Creates a 2-atom zincblende lattice with lattice constants of 'a'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine zincblende(myatoms, a, Z1, Z2)
    type(Atoms), intent(out)      :: myatoms
    real(dp), intent(in)          :: a
    integer, intent(in), optional :: Z1, Z2

    call initialise(myatoms, 2, &
         reshape( (/0.0_dp,   0.5_dp*a, 0.5_dp*a, &
                  & 0.5_dp*a, 0.0_dp,   0.5_dp*a, &
                  & 0.5_dp*a, 0.5_dp*a, 0.0_dp /),(/3,3/)))
    
    myatoms%pos(:,1) = matmul(myatoms%lattice, (/ 0.00_dp, 0.00_dp, 0.00_dp /))
    myatoms%pos(:,2) = matmul(myatoms%lattice, (/ 0.25_dp, 0.25_dp, 0.25_dp /))

    if (present(Z1).and.present(Z2)) call set_atoms(myatoms,(/Z1,Z2/))
    if (present(Z1).and. .not.present(Z2)) call set_atoms(myatoms,Z1)
    if (.not. present(Z1).and. present(Z2)) call set_atoms(myatoms,Z2)

  end subroutine zincblende

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! wurtzite(myatoms, a, Z1, Z2)
  !
  !% Creates a 4-atom wurtzite lattice with lattice constants of 'a' and 'c'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine wurtzite(myatoms, a, c, Z1, Z2, u)
    type(Atoms), intent(out)       :: myatoms
    real(dp), intent(in)           :: a
    real(dp), intent(in), optional :: c, u
    integer, intent(in), optional  :: Z1, Z2

    real(dp) :: my_c, my_u

    my_c = optional_default(sqrt(8.0_dp/3.0_dp)*a,c)
    my_u = optional_default(3.0_dp/8.0_dp,u)

    call initialise(myatoms, 4, &
         reshape( (/0.5_dp*a,-0.5_dp*sqrt(3.0_dp)*a, 0.0_dp, &
                  & 0.5_dp*a, 0.5_dp*sqrt(3.0_dp)*a, 0.0_dp, &
                  & 0.0_dp,   0.0_dp,                my_c/),(/3,3/)))
    
    myatoms%pos(:,1) = matmul(myatoms%lattice, (/ 1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.00_dp /))
    myatoms%pos(:,2) = matmul(myatoms%lattice, (/ 2.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.50_dp /))
    myatoms%pos(:,3) = matmul(myatoms%lattice, (/ 1.0_dp/3.0_dp, 2.0_dp/3.0_dp, my_u    /))
    myatoms%pos(:,4) = matmul(myatoms%lattice, (/ 2.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.50_dp+my_u /))

    if (present(Z1).and.present(Z2)) call set_atoms(myatoms,(/Z1,Z1,Z2,Z2/))
    if (present(Z1).and. .not.present(Z2)) call set_atoms(myatoms,Z1)
    if (.not. present(Z1).and. present(Z2)) call set_atoms(myatoms,Z2)

  end subroutine wurtzite

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

    call initialise(myatoms, 8, &
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

    call initialise(myatoms, 4, &
         reshape( (/0.5_dp*a,-0.5_dp*sqrt(3.0_dp)*a, 0.0_dp, &
                  & 0.5_dp*a, 0.5_dp*sqrt(3.0_dp)*a, 0.0_dp, &
                  & 0.0_dp,   0.0_dp,                c/),(/3,3/)))
    
    myatoms%pos(:,1) = matmul(myatoms%lattice, (/ 0.0_dp,        0.0_dp,        0.25_dp /))
    myatoms%pos(:,2) = matmul(myatoms%lattice, (/ 0.0_dp,        0.0_dp,        0.75_dp /))
    myatoms%pos(:,3) = matmul(myatoms%lattice, (/ 1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.25_dp /))
    myatoms%pos(:,4) = matmul(myatoms%lattice, (/ 2.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.75_dp /))

    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine graphite

  subroutine graphite_rhombohedral(myatoms, a, c, Z)
    type(Atoms), intent(out)      :: myatoms
    real(dp), intent(in)          :: a, c
    integer, intent(in), optional :: Z

    call initialise(myatoms, 6, &
         reshape( (/0.5_dp*a,-0.5_dp*sqrt(3.0_dp)*a, 0.0_dp, &
                  & 0.5_dp*a, 0.5_dp*sqrt(3.0_dp)*a, 0.0_dp, &
                  & 0.0_dp,   0.0_dp,                c/),(/3,3/)))
    
    myatoms%pos(:,1) = matmul(myatoms%lattice, (/ 0.0_dp,        0.0_dp,        0.0_dp /))
    myatoms%pos(:,2) = matmul(myatoms%lattice, (/ 1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.0_dp /))
    myatoms%pos(:,3) = matmul(myatoms%lattice, (/ 0.0_dp,        0.0_dp,        1.0_dp/3.0_dp /))
    myatoms%pos(:,4) = matmul(myatoms%lattice, (/ 2.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp /))
    myatoms%pos(:,5) = matmul(myatoms%lattice, (/ 1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 2.0_dp/3.0_dp /))
    myatoms%pos(:,6) = matmul(myatoms%lattice, (/ 2.0_dp/3.0_dp, 1.0_dp/3.0_dp, 2.0_dp/3.0_dp /))

    if (present(Z)) call set_atoms(myatoms,Z)

  end subroutine graphite_rhombohedral

  !%  Primitive 9-atom trigonal alpha quartz cell
  subroutine alpha_quartz(at, a, c, u, x, y, z)
    type(Atoms), intent(out) :: at
    real(dp), intent(in) :: a, c, u, x, y, z

    real(dp) :: lattice(3,3), a1(3), a2(3), a3(3)

    lattice = 0.0_dp
    a1 = (/ 0.5_dp*a, -0.5_dp*sqrt(3.0_dp)*a, 0.0_dp /)
    a2 = (/ 0.5_dp*a,  0.5_dp*sqrt(3.0_dp)*a, 0.0_dp /)
    a3 = (/ 0.0_dp,    0.0_dp,                c /)
    lattice(:,1) = a1 
    lattice(:,2) = a2 
    lattice(:,3) = a3
   
    call initialise(at, n=9, lattice=lattice)
    call set_atoms(at, (/14,14,14,8,8,8,8,8,8 /))

    at%pos(:,1) =  u*a1 + 2.0_dp/3.0_dp*a3
    at%pos(:,2) =  u*a2 + 1.0_dp/3.0_dp*a3
    at%pos(:,3) = -u*a1 - u*a2
    at%pos(:,4) =  x*a1 + y*a2 + (z + 2.0_dp/3.0_dp)*a3
    at%pos(:,5) = -y*a1 + (x-y)*a2  + (4.0_dp/3.0_dp + z)*a3
    at%pos(:,6) = (y-x)*a1 - x*a2   + (1.0_dp+ z)*a3
    at%pos(:,7) = y*a1 + x*a2 - (2.0_dp/3.0_dp + z)*a3
    at%pos(:,8) = -x*a1 + (y-x)*a2  - z*a3
    at%pos(:,9) = (x - y)*a1 - y*a2 - (1.0_dp/3.0_dp + z)*a3

    call add_property(at, 'primitive_index', (/ 1,2,3,4,5,6,7,8,9 /))
    
 end subroutine alpha_quartz

 !%  Non-primitive 18-atom cubic quartz cell
 subroutine alpha_quartz_cubic(at, a, c, u, x, y, z)
   type(Atoms), intent(out) :: at
   real(dp), intent(in) :: a, c, u, x, y, z

   type(Atoms) :: a0, a1
   real(dp) :: lattice(3,3), g(3,3), t(3)
   logical, allocatable, dimension(:) :: unit_cell
   integer :: i

   call alpha_quartz(a0, a, c, u, x, y, z)
   call supercell(a1, a0, 4, 4, 1)
   call map_into_cell(a1)

   lattice = 0.0_dp
   lattice(1,1) = a0%lattice(1,1)*2.0_dp
   lattice(2,2) = a0%lattice(2,2)*2.0_dp
   lattice(3,3) = a0%lattice(3,3)
   call matrix3x3_inverse(lattice,g)

   allocate(unit_cell(a1%n))
   unit_cell = .false.
   do i=1,a1%n
      t = g .mult. a1%pos(:,i)
      if (all(t >= -0.5) .and. all(t < 0.5)) unit_cell(i) = .true.
   end do

   call select(at, a1, mask=unit_cell)
   call set_lattice(at, lattice, scale_positions=.false.)

   deallocate(unit_cell)

 end subroutine alpha_quartz_cubic

 subroutine rutile(at, a, c, u)
!For atomic coordinates see Ref. PRB 65, 224112 (2002).
!Experimental values: a=4.587, c=2.954, u=0.305 (Ref. JACS 109, 3639 (1987)).
    type(Atoms), intent(out) :: at
    real(dp), intent(in) :: a, c, u

    real(dp) :: lattice(3,3), a1(3), a2(3), a3(3)

    lattice = 0.0_dp
    a1 = (/ a, 0.0_dp, 0.0_dp /)
    a2 = (/ 0.0_dp, a, 0.0_dp /)
    a3 = (/ 0.0_dp, 0.0_dp, c /)
    lattice(:,1) = a1
    lattice(:,2) = a2
    lattice(:,3) = a3

    call initialise(at, n=6, lattice=lattice)
    call set_atoms(at, (/22,22,8,8,8,8 /))

    at%pos(:,1) =  0.0_dp 
    at%pos(:,2) =  0.5_dp * a1 + 0.5_dp*a2 + 0.5_dp*a3
    at%pos(:,3) =  u*a1 + u*a2
    at%pos(:,4) = -u*a1 - u*a2 
    at%pos(:,5) = (0.5_dp + u)*a1 + (0.5_dp - u)*a2 + 0.5_dp*a3 
    at%pos(:,6) = (0.5_dp - u)*a1 + (0.5_dp + u)*a2 + 0.5_dp*a3

    call add_property(at, 'primitive_index', (/ 1,2,3,4,5,6/))

 end subroutine rutile 

 subroutine anatase(at, a, c, u)
!Conventional 12 atoms unit cell defined as a function of a,c,u.
!Experimental values: a=3.782, c=9.502, u=0.208 (Ref. JACS 109, 3639 (1987)).
    type(Atoms), intent(out) :: at
    real(dp), intent(in) :: a, c, u

    real(dp) :: lattice(3,3), a1(3), a2(3), a3(3)

    lattice = 0.0_dp
    a1 = (/ a, 0.0_dp, 0.0_dp /)
    a2 = (/ 0.0_dp, a, 0.0_dp /)
    a3 = (/ 0.0_dp, 0.0_dp, c  /)
    lattice(:,1) = a1
    lattice(:,2) = a2
    lattice(:,3) = a3

    call initialise(at, n=12, lattice=lattice)
    call set_atoms(at, (/22,22,22,22,8,8,8,8,8,8,8,8 /))

    at%pos(:,1) = (/0.0_dp*a, 0.5_dp*a,  0.0_dp   /)
    at%pos(:,2) = (/0.0_dp,   0.0_dp  ,  0.25_dp*c/) 
    at%pos(:,3) = (/0.5_dp*a, 0.5_dp*a, -0.25_dp*c/) 
    at%pos(:,4) = (/0.5_dp*a, 0.0_dp  ,  0.50_dp*c/)  
    at%pos(:,5) = (/0.0_dp  , 0.5_dp*a,  u * c   /) 
    at%pos(:,6) = (/0.0_dp  , 0.5_dp*a, -u * c   /) 
    at%pos(:,7) = (/0.0_dp  , 0.0_dp  ,  (0.25_dp+u)*c /) 
    at%pos(:,8) = (/0.0_dp  , 0.0_dp  ,  (0.25_dp-u)*c /) 
    at%pos(:,9) = (/0.5_dp*a, 0.5_dp*a,  (0.75_dp-u)*c /) 
    at%pos(:,10)= (/0.5_dp*a, 0.0_dp  ,  (0.50_dp-u)*c /) 
    at%pos(:,11)= (/0.5_dp*a, 0.0_dp  , -(0.50_dp-u)*c /) 
    at%pos(:,12)= (/0.5_dp*a, 0.5_dp*a, -(0.25_dp-u)*c /) 

    call add_property(at, 'primitive_index', (/ 1,2,3,4,5,6,7,8,9,10,11,12/))

 end subroutine anatase 

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

    call initialise(structure,(motif%N * points%N),lattice)

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

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% create a supercell of a structure, read from a file, with chosen
  !% volume per atom or volume per unit cell, with desired supercell
  !% repeats, and specified Z values.
  !% file may contain default Z values as a property Z_values='Z1 Z2 ...'
  !% structures that begin with . or / are searched for as paths,
  !% and everything else is searched for in 'QUIP_ARCH/structures/struct.xyz'
  !% or in 'HOME/share/quip_structures/struct.xyz'
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function structure_from_file(struct, vol_per_atom, vol_per_unit_cell, repeat, Z_values_str, error) result(dup_cell)
    character(len=*), intent(in) :: struct
    real(dp), intent(in), optional :: vol_per_atom, vol_per_unit_cell
    integer, intent(in), optional :: repeat(3)
    character(len=*), intent(in), optional :: Z_values_str
    integer, intent(out), optional :: error
    type(Atoms) :: dup_cell

    character(len=FIELD_LENGTH) :: struct_file
    character(len=FIELD_LENGTH) :: quip_structs_dir
    type(Atoms) :: cell
    real(dp) :: vol, scale
    integer :: stat
    logical :: cell_has_Z_values, args_have_Z_values
    integer :: Z_values_i
    integer, allocatable :: Z_values(:), new_Z(:)
    integer :: i, n_types
    real(dp) :: u_vol_per_atom, u_vol_per_unit_cell
    type(Dictionary) :: Z_values_params
    integer :: l_error

    INIT_ERROR(error)

    if (minval(repeat) <= 0) &
      call system_abort("all elements of repeat='"//repeat//"' must be >= 1")

    u_vol_per_atom = optional_default(-1.0_dp, vol_per_atom)
    u_vol_per_unit_cell = optional_default(-1.0_dp, vol_per_unit_cell)

    if (u_vol_per_atom > 0.0_dp .and. u_vol_per_unit_cell > 0.0_dp) &
      call system_abort("Only one of u_vol_per_atom="//u_vol_per_atom//" and u_vol_per_unit="//u_vol_per_unit_cell//" cell can be specified")
    if (u_vol_per_atom <= 0.0_dp .and. u_vol_per_unit_cell <= 0.0_dp) &
      call system_abort("One of u_vol_per_atom ="//u_vol_per_atom//"and u_vol_per_unit ="//u_vol_per_unit_cell//"cell must be specified")

    if (struct(1:1) == '.' .or. struct(1:1) == '/') then
      struct_file = trim(struct)
    else
      call get_env_var("QUIP_STRUCTS_DIR", quip_structs_dir, stat)
      if (stat /= 0) then
	call get_env_var("QUIP_ROOT", quip_structs_dir, stat)
	if (stat /= 0) then
	  call get_env_var("HOME", quip_structs_dir, stat)
	  if (stat /= 0) &
	    call system_abort("Could not get QUIP_STRUCTS_DIR or QUIP_ROOT or HOME env variables")
	  quip_structs_dir = trim(quip_structs_dir) // "/share/quip_structures"
	else
	  quip_structs_dir = trim(quip_structs_dir) // "/structures"
	endif
      endif
      struct_file = trim(quip_structs_dir)//"/"//trim(struct)//".xyz"
    endif
    call read(cell, struct_file, error=l_error)
    if (l_error /= 0) then
      call print("structure_from_file looks for structure files in this order:", PRINT_ALWAYS)
      call print("  $QUIP_STRUCTS_DIR", PRINT_ALWAYS)
      call print("  $QUIP_ROOT/structures", PRINT_ALWAYS)
      call print("  $HOME/share/quip_structures", PRINT_ALWAYS)
      RAISE_ERROR("file not found", error)
    endif

    n_types=maxval(cell%Z)-minval(cell%Z) + 1
    allocate(Z_values(n_types))
    allocate(new_Z(cell%N))

    if (n_types == 1) then
      cell_has_Z_values = get_value(cell%params, "z_values", Z_values_i)
      Z_values(1) = Z_values_i
    else
      cell_has_Z_values = get_value(cell%params, "z_values", Z_values)
    endif
    args_have_Z_values = .false.
    if (present(Z_values_str)) then
      if (len_trim(Z_values_str) > 0) then
	call initialise(Z_values_params)
	if (n_types == 1) then
	  call param_register(Z_values_params, "Z_values", PARAM_MANDATORY, Z_values_i, help_string="No help yet.  This source file was $LastChangedBy$")
	else
	  call param_register(Z_values_params, "Z_values", PARAM_MANDATORY, Z_values, help_string="No help yet.  This source file was $LastChangedBy$")
	endif
	if (.not. param_read_line(Z_values_params, 'Z_values="'//trim(Z_values_str)//'"', ignore_unknown=.true.)) then
	  call system_abort("Z_values '"//trim(Z_values_str)//"' specified on command line, but can't be parsed to find "//n_types//" Z values")
	else ! succeeded in parsing
	  if (n_types == 1) Z_values(1) = Z_values_i
	endif
	call finalise(Z_values_params)
      endif
    endif
    if (.not. cell_has_Z_values .and. .not. args_have_Z_values) &
      call system_abort("structure_from_file could not find appropriate Z values in cell file or Z_values_str argument")

    do i=1, n_types
      where (cell%Z == i)
	new_Z = Z_values(i)
      end where
    end do
    call set_atoms(cell, new_Z)

    vol = cell_volume(cell)
    if (u_vol_per_atom > 0.0_dp) then
      scale = 1.0_dp/vol**(1.0_dp/3.0_dp) * (cell%N*u_vol_per_atom)**(1.0_dp/3.0_dp)
    else
      scale = 1.0_dp/vol**(1.0_dp/3.0_dp) * u_vol_per_atom**(1.0_dp/3.0_dp)
    endif
    call set_lattice(cell, cell%lattice*scale, scale_positions=.true.)
    call supercell(dup_cell, cell, repeat(1), repeat(2), repeat(3))

    deallocate(Z_values)
    deallocate(new_Z)

  end function structure_from_file

  subroutine transform(at_out, at_in, t)
    type(Atoms), intent(out)::at_out  !% Output 
    type(Atoms), intent(in) ::at_in   !% Input 
    real(dp), intent(in)::t(3,3)

    at_out = at_in
    call set_lattice(at_out,(t .mult. at_in%lattice), scale_positions=.true.)
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
  !% if find_all_possible_matches is true, all possible matches, not just non-overlapping ones,
  !% are returned.  Useful for situations where things are ambiguous and need to be resolved
  !% with more information outside this routine
  !%
  !% The routine could optionally find hysteretically defined connections between neighbours,
  !% if the alt_connect's cutoff were the same as at%cutoff(_break)
  !%
  subroutine find_motif(at,motif,matches,start,end,mask,find_all_possible_matches,nneighb_only,alt_connect) !,hysteretic_neighbours)

    type(atoms),       intent(in), target  :: at           !% The atoms structure to search
    integer,           intent(in)  :: motif(:,:)   !% The motif to search for
    type(table),       intent(out) :: matches      !% All matches
    integer, optional, intent(in)  :: start, end   !% Start and End atomic indices for search
    logical, optional, intent(in)  :: mask(:)      !% If present only masked atoms are searched
    logical, optional, intent(in)  :: find_all_possible_matches !% if true, don't exclude matches that overlap
    logical, intent(in), optional :: nneighb_only
    type(Connection), intent(in), optional, target :: alt_connect
!    logical, optional, intent(in) :: hysteretic_neighbours

    character(*), parameter  :: me = 'find_motif: '

    integer,     allocatable :: A(:,:),B(:,:),C(:,:),depth(:), neighbour_Z(:), Z(:), M0(:,:), M(:,:), &
                                match_indices(:), depth_real(:)
    integer                  :: i,j,k,p,q,N,max_depth,max_depth_real, my_start,my_end, opt_atom, ji
    logical                  :: match
    type(table)              :: neighbours, core, num_species_at, num_species_motif
    logical, allocatable :: assigned_to_motif(:)
    logical :: do_append
    logical :: do_find_all_possible_matches
    logical :: my_nneighb_only

    integer                  :: discards(7)
!    logical :: use_hysteretic_neighbours
!    real(dp) :: cutoff

    call print("find_motif", verbosity=PRINT_ANAL)
    call print(motif, verbosity=PRINT_ANAL)

    my_nneighb_only = optional_default(.true., nneighb_only)

    discards = 0
    do_find_all_possible_matches = optional_default(.false., find_all_possible_matches)

    !XXXXXXXXXXXXXXXX
    !X INITIALISATION
    !XXXXXXXXXXXXXXXX

    ! Check atomic numbers
    if (any(motif(:,1)<1)) call system_abort(me//'Bad atomic numbers ('//motif(:,1)//')')

    ! Check for atomic connectivity
    if (present(alt_connect)) then
      if (.not.alt_connect%initialised) call system_abort(me//'got alt_connect uninitialised')
    else
      if (.not.at%connect%initialised) call system_abort(me//'got at%connect uninitialised')
    endif
!    use_hysteretic_neighbours = optional_default(.false.,hysteretic_neighbours)
!
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

    allocate(assigned_to_motif(at%N))
    assigned_to_motif = .false.

    !Loop over all atoms looking for candidates for the optimum atom 'opt_atom'
    call allocate(core,4,0,0,0,1)
    allocate(neighbour_Z(count(A(opt_atom,:)==1)))
    neighbour_Z = pack(Z,(A(opt_atom,:)==1))

    call print("opt_atom " // opt_atom // " Z " // Z(opt_atom), verbosity=PRINT_ANAL)

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

       call print("check atom " // i // " Z " // at%Z(i), verbosity=PRINT_ANAL)

       if (present(mask)) then
          if (.not.mask(i)) then
	    call print("   i not in mask, skipping", verbosity=PRINT_ANAL)
	    cycle
	   endif
       end if

       !XXXXXXXXXXXXXXXXXXXX
       !X ESCAPE ROUTE 1
       !XXXXXXXXXXXXXXXXXXXX

       ! Discard atoms which don't match the atomic number of the optimum atom

       if (at%Z(i) /= Z(opt_atom)) then
	  call print("   i has wrong Z to match opt_atom, skipping", verbosity=PRINT_ANAL)
          discards(1) = discards(1) + 1
          cycle
       end if

       !---------------------

       !XXXXXXXXXXXXXXXXXXXX
       !X ESCAPE ROUTE 2
       !XXXXXXXXXXXXXXXXXXXX

       ! Check if the real atom's neighbours can be matched with those in the motif

       if (.not.atoms_compatible(at,i,A,Z,opt_atom,nneighb_only=nneighb_only,alt_connect=alt_connect)) then
          discards(2) = discards(2) + 1
	  call print("   i has incompatible neighbour list, skipping", verbosity=PRINT_ANAL)
          cycle
       end if

       !---------------------

       ! Grow a cluster around the real atom i which is  max_depth hops deep
       call wipe(core)
       call append(core,(/i,0,0,0/))
       call bfs_grow(at,core,max_depth,nneighb_only = nneighb_only, min_images_only = .true.,alt_connect=alt_connect)

       !XXXXXXXXXXXXXXXXXXXX
       !X ESCAPE ROUTE 3
       !XXXXXXXXXXXXXXXXXXXX

       ! If the number of atoms in the cluster is less than those in the motif we can't have a match
       if (core%N < N) then
	  call print("   i's cluster is too small, skipping", verbosity=PRINT_ANAL)
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
	  call print("   i's cluster has mismatching atom numbers, skipping", verbosity=PRINT_ANAL)
          discards(4) = discards(4) + 1
          cycle
       end if

       !---------------------

       ! Create adjacency matrix
       allocate(B(core%N,core%N))
       B = 0
       do p = 1, core%N
          k = core%int(1,p)
	  do ji=1, atoms_n_neighbours(at, k, alt_connect=alt_connect)
	    j = atoms_neighbour(at, k, ji, alt_connect=alt_connect)
	    if (my_nneighb_only) then
	       if (.not. is_nearest_neighbour(at, k, ji, alt_connect=alt_connect)) cycle
	    endif
	    q = find_in_array(core%int(1,:),j)
	    if (q > 0) then
	       B(p,q) = 1
	       B(q,p) = 1
	    endif
	 end do ! ji
       end do !p

       call print("core", verbosity=PRINT_ANAL)
       call print(core, verbosity=PRINT_ANAL)
       call print("adjacency mat", verbosity=PRINT_ANAL)
       call print(B, verbosity=PRINT_ANAL)

       ! Find depth of connectivity for real atoms
       allocate(depth_real(core%N))
       call find_connectivity_depth(B,1,depth_real)
       max_depth_real = maxval(depth_real)

       !XXXXXXXXXXXXXXXXXX
       !X ESCAPE ROUTE 5
       !XXXXXXXXXXXXXXXXXX

       ! If there are atoms which aren't deep enough in the real structure then a match is impossible
       if (max_depth > max_depth_real) then
	  call print("   i's cluster isn't actually keep enough, skipping", verbosity=PRINT_ANAL)
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
                if (atoms_compatible(at,core%int(1,q),A,Z,p,nneighb_only=nneighb_only,alt_connect=alt_connect)) M0(p,q) = 1
             end if
          end do
       end do

       deallocate(depth_real)

       !XXXXXXXXXXXXXXXXXX
       !X ESCAPE ROUTE 6
       !XXXXXXXXXXXXXXXXXX

       ! Check to see if any atoms have no possibilities
       if (any(sum(M0,dim=2)==0)) then
	  call print("   i has some neighbor that doesn't match, skipping", verbosity=PRINT_ANAL)
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

	  call print("  check permutation loop start", verbosity=PRINT_ANAL)
          ! For each trial matrix create the permuted adjacency matrix
          C = matmul(M,transpose(matmul(M,B))) ! ***** THIS LINE IS ANOTHER CANDIDATE FOR OPTIMISATION *****
                                               ! Use sparse matrices maybe?

          if (all(C==A)) then ! match
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

	     do_append = .true.
	     ! if mask is present, make sure all atoms are in mask, otherwise skip
	     if (present(mask)) then
		 if (.not. all(mask(match_indices))) do_append = .false.
	     end if
	     if (do_find_all_possible_matches) then
	       ! if we're finding all possible matches, skip only if _all_ atoms are already in a motif
	       if(all(assigned_to_motif(match_indices))) do_append = .false.
	     else
	       ! if we're not, skip if _any_ atom is already in a motif
	       if (any(assigned_to_motif(match_indices))) do_append = .false.
	     endif

	     if (do_append) then
		call append(matches,match_indices)
		assigned_to_motif(match_indices) = .true.
		call print("  found match, indices " // match_indices, verbosity=PRINT_ANAL)
	     endif
	     deallocate(match_indices)

	     if (.not. do_find_all_possible_matches) then
	       call print("  not looking for _all_ matches, finished", verbosity=PRINT_ANAL)
	       exit
	     endif
          end if ! match

          call next_trial_matrix(M0,M)

          !XXXXXXXXXXXXXXXXXX
          !X ESCAPE ROUTE 7
          !XXXXXXXXXXXXXXXXXX

          ! If next_trial_matrix deletes the only atom we know to be fitted then all
          ! permutations have been exhausted
          if (M(opt_atom,1)==0) then
	     call print("   i has no more premutations, leaving loop", verbosity=PRINT_ANAL)
             discards(7) = discards(7) + 1
             exit
          end if
          !---------------------

       end do


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
  ! Added alt_connect.
  function atoms_compatible(at,i,A,Z,j,nneighb_only,alt_connect)

    type(atoms),             intent(in), target :: at !the atoms structure
    integer,                 intent(in) :: i  !the atom in the atoms structure we are testing
    integer, dimension(:,:), intent(in) :: A  !the adjacency matrix of the motif
    integer, dimension(:),   intent(in) :: Z  !the atomic number of the motif
    integer,                 intent(in) :: j  !the atom in the motif we are testing
    logical, intent(in), optional :: nneighb_only
    type(Connection), intent(in), optional, target :: alt_connect
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
    call bfs_step(at,core,neighbours,nneighb_only=nneighb_only,min_images_only=.true.,alt_connect=alt_connect)

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
