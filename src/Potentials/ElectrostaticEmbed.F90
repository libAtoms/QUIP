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

!% ElectrostaticEmbed - tools for electrostatic embedding of potentials

#include "error.inc"

module ElectrostaticEmbed_module

  use error_module
  use system_module, only : dp, inoutput, optional_default, initialise, OUTPUT, operator(//), system_timer
  use units_module
  use periodictable_module
  use linearalgebra_module
  use dictionary_module
  use connection_module
  use atoms_types_module
  use atoms_module
  use clusters_module

  use QUIP_Common_module
  use Functions_module
  use Potential_module
#ifdef HAVE_THIRDPARTY
  use cube_tools
#endif

  implicit none
  private

  integer, private, parameter :: nspins = 2

  public :: extent_to_ngrid, assign_grid_coordinates, calc_electrostatic_potential, &
#ifdef HAVE_THIRDPARTY
       write_electrostatic_potential_cube, &
#endif
       make_periodic_potential

contains

  subroutine extent_to_ngrid(extent, fine_gmax, ngrid)
    real(dp), intent(in) :: extent(3), fine_gmax
    integer, intent(out) :: ngrid(3)

    integer i

    do i=1,3
       ngrid(i) = ceiling(fine_gmax*extent(i)/BOHR/PI)
       call round_prime_factors(ngrid(i))
    end do
    
  end subroutine extent_to_ngrid

  subroutine assign_grid_coordinates(ngrid, origin, extent, grid_size, r_real_grid)
    integer, intent(in) :: ngrid(3)
    real(dp), intent(in) :: origin(3), extent(3)
    real(dp), intent(out) :: grid_size(3)
    real(dp), dimension(:,:), intent(out) :: r_real_grid

    real(dp) :: lattice(3,3)
    integer :: nx, ny, nz, point

    lattice(:,:) = 0.0_dp
    lattice(1,1) = extent(1)
    lattice(2,2) = extent(2)
    lattice(3,3) = extent(3)

    point = 0
    do nz=1,ngrid(3)
       do ny=1,ngrid(2)
          do nx=1,ngrid(1)

             point=point+1

             r_real_grid(1,point) = real(nx-1,dp)*lattice(1,1)/real(ngrid(1),dp) +&
                  & real(ny-1,dp)*lattice(2,1)/real(ngrid(2),dp) +&
                  & real(nz-1,dp)*lattice(3,1)/real(ngrid(3),dp) + origin(1)
             r_real_grid(2,point) = real(nx-1,dp)*lattice(1,2)/real(ngrid(1),dp) +&
                  & real(ny-1,dp)*lattice(2,2)/real(ngrid(2),dp) +&
                  & real(nz-1,dp)*lattice(3,2)/real(ngrid(3),dp) + origin(2)
             r_real_grid(3,point) = real(nx-1,dp)*lattice(1,3)/real(ngrid(1),dp) +&
                  & real(ny-1,dp)*lattice(2,3)/real(ngrid(2),dp) +&
                  & real(nz-1,dp)*lattice(3,3)/real(ngrid(3),dp) + origin(3)

          end do
       end do
    end do

    grid_size(1) = extent(1)/real(ngrid(1), dp)
    grid_size(2) = extent(2)/real(ngrid(2), dp)
    grid_size(3) = extent(3)/real(ngrid(3), dp)

  end subroutine assign_grid_coordinates

#ifdef HAVE_THIRDPARTY
  subroutine write_electrostatic_potential_cube(at, filename, ngrid, origin, extent, pot, &
       write_efield, flip_sign, convert_to_atomic_units, error)
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ngrid(3)
    real(dp), intent(in) :: origin(3), extent(3)
    real(dp), dimension(:,:,:), intent(in) :: pot
    logical, optional, intent(in) :: write_efield, flip_sign, convert_to_atomic_units
    integer, optional, intent(out) :: error

    type(InOutput) :: out
    type(cube_type) :: cube
    type(Dictionary) :: metadata
    integer i, j, n, z!, na, nb, nc
    real(dp), pointer, dimension(:) :: charge, es_pot
    real(dp), pointer, dimension(:,:) :: es_efield
    integer, dimension(size(ElementMass)) :: nsp
    logical do_efield, do_atomic_units, do_flip_sign

    INIT_ERROR(error)

    do_efield = optional_default(.true., write_efield)
    do_flip_sign = optional_default(.false., flip_sign)
    do_atomic_units = optional_default(.true., convert_to_atomic_units)

    if (do_efield) then
       call assign_property_pointer(at, 'es_pot', es_pot, error)
       PASS_ERROR(error)

       call assign_property_pointer(at, 'es_efield', es_efield, error)
       PASS_ERROR(error)
    else
       call assign_property_pointer(at, 'charge', charge, error)
       PASS_ERROR(error)
    end if

    call cube_clear(cube)
    cube%comment1 = trim(filename)
    cube%N = at%n
    cube%r0 = real(origin)
    cube%na = ngrid(1)+1
    cube%nb = ngrid(2)+1
    cube%nc = ngrid(3)+1
    
    ! Size of voxels
    cube%da = 0.0_dp; cube%da(1) = real(extent(1)/ngrid(1))
    cube%db = 0.0_dp; cube%db(2) = real(extent(2)/ngrid(2))
    cube%dc = 0.0_dp; cube%dc(3) = real(extent(3)/ngrid(3))

    ! We sort atoms by atomic number for compatibility with CASTEP
    
    ! Count number of each species
    nsp(:) = 0
    do i=1,at%n
       nsp(at%z(i)) = nsp(at%z(i)) + 1
    end do

    ! Fill in atomic information
    allocate(cube%atoms(cube%n))
    n = 0
    do z = 1,size(nsp)
       if (nsp(z) == 0) cycle

       do j=1,at%n
          if (at%z(j) /= z) cycle

          n = n + 1
          cube%atoms(n)%number = at%z(j)
          if (do_efield) then
             ! Potential on ions and electric field
             cube%atoms(n)%unknown = real(es_pot(j))
             if (do_flip_sign) cube%atoms(n)%unknown = -cube%atoms(n)%unknown
             cube%atoms(n)%r(:) = real(es_efield(:,j))
             if (do_flip_sign) cube%atoms(n)%r = -cube%atoms(n)%r
          else
             ! Charge and position of ions
             cube%atoms(n)%unknown = real(charge(j))
             cube%atoms(n)%r(:) = real(at%pos(:,j))
          end if
       end do
    end do

    ! Fill in volumetric information - note conversion to single precision
    cube%NMOs=1
    allocate(cube%voxels(cube%na,cube%nb,cube%nc,cube%NMOs))
    cube%voxels(1:cube%na-1,1:cube%nb-1,1:cube%nc-1,1) = real(pot)

    ! periodic boundary conditions: in cube file format, values of
    ! n{a,b,c}=1 entries are repeated at n{a,b,c} = cube%n{a,b,c}
    cube%voxels(cube%na,:,:,1) = real(cube%voxels(1,:,:,1))
    cube%voxels(:,cube%nb,:,1) = real(cube%voxels(:,1,:,1))
    cube%voxels(:,:,cube%nc,1) = real(cube%voxels(:,:,1,1))

    if (do_flip_sign) cube%voxels = -cube%voxels

    if (do_atomic_units) then
       cube%da = cube%da/BOHR
       cube%db = cube%db/BOHR
       cube%dc = cube%dc/BOHR
       cube%r0 = cube%r0/BOHR
       cube%voxels = cube%voxels/HARTREE ! convert potential to Hartree
       if (do_efield) then
          ! convert potential and electric field to atomic units
          do n=1,cube%n
             cube%atoms(n)%unknown = cube%atoms(n)%unknown/HARTREE
             cube%atoms(n)%r(:) = cube%atoms(n)%r(:)/(HARTREE/BOHR)
          end do
       else
          ! convert positions to bohr
          do n=1,cube%n
             cube%atoms(n)%r = cube%atoms(n)%r/BOHR
          end do
       end if
    end if

    ! Write  meta-data into second line of .cube file
    call initialise(metadata)
    call set_value(metadata, 'volumetric_data', 'electrostatic_potential')
    if (do_atomic_units) then
       call set_value(metadata, 'units', 'atomic')
    else
       call set_value(metadata, 'units', 'eV_A_fs')
    end if
    if (do_efield) then
       call set_value(metadata, 'atom_data', 'z:pot:efield')
    else
       call set_value(metadata, 'atom_data', 'z:charge:pos')
    end if
    cube%comment2 = write_string(metadata)
    call finalise(metadata)

    call initialise(out, filename, OUTPUT)
    call cube_write(cube, out%unit, error)
    PASS_ERROR(error)
    call finalise(out)
    deallocate(cube%atoms)
    deallocate(cube%voxels)

  end subroutine write_electrostatic_potential_cube
#endif

  subroutine make_periodic_potential(at, real_grid, ngrid, is_periodic, cutoff_radius, cutoff_width, mark_name, pot, error)
    type(Atoms), intent(inout) :: at
    real(dp), dimension(:,:), intent(in) :: real_grid
    integer, intent(in) :: ngrid(3)
    logical, intent(in) :: is_periodic(3)
    real(dp), intent(in) :: cutoff_radius, cutoff_width
    character(len=*), intent(in) :: mark_name
    real(dp), intent(inout) :: pot(:,:,:)
    integer, optional, intent(out) :: error

    integer i,j,k,a,igrid, nsurface, cell_image_Na, cell_image_Nb, cell_image_Nc
    logical save_is_periodic(3), dir_mask(3)
    real(dp) :: cutoff, surface_avg, d(3), fc, dfc, masked_d
    logical, dimension(:), allocatable :: atom_mask
    integer, pointer, dimension(:) :: mark

    INIT_ERROR(error)
    
    call print('is_periodic = '//is_periodic)
    save_is_periodic = at%is_periodic
    at%is_periodic = is_periodic
    dir_mask = .not. is_periodic
    call calc_connect(at)
    call print('at%is_periodic = '//at%is_periodic//' dir_mask='//dir_mask)
    
    ! Calculate average value of potential on open surfaces
    nsurface = 0
    surface_avg = 0.0_dp
    if (dir_mask(1)) then
       surface_avg = surface_avg + sum(pot(1,:,:)) + sum(pot(ngrid(1),:,:))
       nsurface = nsurface + 2*ngrid(2)*ngrid(3)
    end if
    if (dir_mask(2)) then
       surface_avg = surface_avg + sum(pot(:,1,:)) + sum(pot(:,ngrid(2),:))
       nsurface = nsurface + 2*ngrid(1)*ngrid(3)
    end if
    if (dir_mask(3)) then
       surface_avg = surface_avg + sum(pot(:,:,1)) + sum(pot(:,:,ngrid(3)))
       nsurface = nsurface + 2*ngrid(1)*ngrid(2)
    end if
    surface_avg = surface_avg/nsurface

    if (.not. assign_pointer(at, mark_name, mark)) then
       RAISE_ERROR('make_periodic_potential: failed to assign pointer to "'//trim(mark_name)//'" property', error)
    end if
    allocate(atom_mask(at%n))
    atom_mask = mark /= HYBRID_ELECTROSTATIC_MARK

    call fit_box_in_cell(cutoff_radius+cutoff_width, cutoff_radius+cutoff_width, cutoff_radius+cutoff_width, &
         at%lattice, cell_image_Na, cell_image_Nb, cell_image_Nc)
    cell_image_Na = max(1,(cell_image_Na+1)/2)
    cell_image_Nb = max(1,(cell_image_Nb+1)/2)
    cell_image_Nc = max(1,(cell_image_Nc+1)/2)

    ! Smoothly interpolate potential to surface_avg away from atoms
    igrid = 0
    do k=1,ngrid(3)
       do j=1,ngrid(2)
          do i=1,ngrid(1)
             igrid = igrid + 1
             a = closest_atom(at, real_grid(:,igrid), cell_image_Na, cell_image_Nb, cell_image_Nc, mask=atom_mask, diff=d)

             ! project difference vector to mask out periodic directions
             masked_d = norm(pack(d, dir_mask))

             ! if this grid point is on an open surface, check we're at least cutoff_radius + cutoff_width away from closest atom
             if (dir_mask(1) .and. i==1 .or. dir_mask(2) .and. j==1 .or. dir_mask(3) .and. k==1) then
                if (masked_d < cutoff_radius+cutoff_width) then
                   RAISE_ERROR('atom '//a//' pos='//at%pos(:,a)//' too close to surface cell=['//i//' '//j//' '//k//'] pos=['//real_grid(:,igrid)//'] d='//masked_d, error)
                end if
             end if

             call smooth_cutoff(masked_d, cutoff_radius, cutoff_width, fc, dfc)
             pot(i,j,k) = pot(i,j,k)*fc + surface_avg*(1.0_dp-fc)
          end do
       end do
    end do

    ! restore periodic directions
    at%is_periodic = save_is_periodic
    call calc_connect(at)

    deallocate(atom_mask)

  end subroutine make_periodic_potential

  !% Evaluate electrostatic potential on a grid by adding test atoms at grid points
  subroutine calc_electrostatic_potential(this, at, cluster, mark_name, ngrid, origin, extent, real_grid, pot, args_str, error)
#ifdef _OPENMP
    use omp_lib
#endif

    type(Potential), intent(inout) :: this
    type(Atoms), intent(in) :: at
    type(Atoms), intent(inout) :: cluster
    character(len=*), intent(in) :: args_str, mark_name
    integer, intent(in) :: ngrid(3)
    real(dp), intent(in) :: origin(3), extent(3)
    real(dp), intent(inout) :: real_grid(:,:), pot(:,:,:)
    integer, optional, intent(out) :: error

    type(Atoms) :: at_copy
    type(Atoms), save :: private_at_copy
    real(dp) :: grid_size(3), tmp_cutoff
    real(dp), pointer :: at_charge(:), cluster_charge(:), cluster_es_pot(:), cluster_es_efield(:,:), es_pot(:), es_efield(:,:)
    real(dp), allocatable :: local_e(:)
    integer, pointer, dimension(:) :: mark
    integer, allocatable, dimension(:) :: z
    integer, allocatable, dimension(:,:) :: lookup
    logical, pointer, dimension(:) :: atom_mask, source_mask
    logical, allocatable, dimension(:) :: select_mask
#ifdef _OPENMP
    logical old_nested
#endif
    integer i, j, k, n, igrid, offset, stride, npoint, select_n

    !$omp threadprivate(private_at_copy)

    INIT_ERROR(error)
    call system_timer('calc_electrostatic_potential')

    call assign_grid_coordinates(ngrid, origin, extent, grid_size, real_grid)


    call assign_property_pointer(at, mark_name, mark, error)
    PASS_ERROR(error)

    ! Make copy of atoms within cutoff distance of a marked atom
    ! so we can add test atoms at grid points

!     allocate(select_mask(at%n))
!     tmp_cutoff = cutoff(this)
!     select_mask = .false.
!     do i=1,at%n
!        dr = diff_min_image(at, at%pos(:,i), origin+extent/2.0_dp)
!        if (all(abs(dr(:)) < extent/2.0_dp+tmp_cutoff)) select_mask(i) = .true.
!     end do
!     call select(at_copy, at, mask=select_mask)
!     select_n = at_copy%n
!     deallocate(select_mask)

    allocate(select_mask(at%n))
    tmp_cutoff = cutoff(this)
    select_mask = .false.
    do i=1,at%n
       if (mark(i) == HYBRID_NO_MARK) cycle
       select_mask(i) = .true.
       do j=1,at%n
          if (mark(j) /= HYBRID_NO_MARK) cycle
          if (distance_min_image(at, i, j) < tmp_cutoff) select_mask(j) = .true.
       end do
    end do
    call select(at_copy, at, mask=select_mask)
    select_n = at_copy%n

    call assign_property_pointer(at_copy, mark_name, mark, error)
    PASS_ERROR(error)

    call assign_property_pointer(at_copy, 'charge', at_charge, error)
    PASS_ERROR(error)

    ! choose balance between number of evaluations and number of simulataneous grid points
    npoint = select_n
    stride = max(1,size(real_grid, 2)/npoint)
    npoint = size(real_grid(:,1:size(real_grid,2):stride), 2)
    call print('npoint = '//npoint//' stride = '//stride)

    allocate(z(npoint))
    z(:) = 0
    call add_atoms(at_copy, real_grid(:,1:size(real_grid,2):stride), z)

    ! added atoms so must reassign pointers
    if (.not. assign_pointer(at_copy, mark_name, mark)) then
       RAISE_ERROR('calc_electrostatic_potential_calc failed to assign pointer to "'//trim(mark_name)//'" property', error)
    end if
    if (.not. assign_pointer(at_copy, 'charge', at_charge)) then
       RAISE_ERROR('calc_electrostatic_potential_calc failed to assign pointer to "charge" property', error)
    endif
    at_charge(select_n+1:at_copy%n) = 1.0_dp

    ! atom_mask is true for atoms for which we want the local potential energy, i.e grid point atoms
    call add_property(at_copy, 'atom_mask', .false., overwrite=.true., ptr=atom_mask, error=error)
    PASS_ERROR(error)
    atom_mask = .false. 
    atom_mask(select_n+1:at_copy%n) = .true.

    ! source_mask is true for atoms which should act as electrostatic sources
    call add_property(at_copy, 'source_mask', .false., overwrite=.true., ptr=source_mask, error=error)
    PASS_ERROR(error)
    source_mask(1:select_n) = mark(1:select_n) == HYBRID_ELECTROSTATIC_MARK
    source_mask(select_n+1:at_copy%n) = .false.

    call print('calc_electrostatic_potential: evaluating electrostatic potential due to '//count(source_mask)//' sources')

    ! Construct mapping from igrid to (i,j,k)
    igrid = 0
    allocate(lookup(3,size(real_grid,2)))
    do k=1,ngrid(3)
       do j=1,ngrid(2)
          do i=1,ngrid(1)
             igrid = igrid + 1
             lookup(:,igrid) = (/i,j,k/)
          end do
       end do
    end do

#ifdef _OPENMP
    old_nested = omp_get_nested()
    call omp_set_nested(.false.)
#endif
    !$omp parallel default(none) shared(at_copy,stride,args_str,grid_size,this,pot,real_grid,select_n,at_charge,lookup) private(local_e,npoint,igrid)

    call atoms_copy_without_connect(private_at_copy, at_copy)
    call set_cutoff(private_at_copy, cutoff(this))
    
    allocate(local_e(private_at_copy%n))

    !$omp do
    do offset=1,stride
       npoint = size(real_grid(:,offset:size(real_grid,2):stride),2)
       private_at_copy%n = select_n + npoint
       private_at_copy%nbuffer = select_n + npoint
       private_at_copy%pos(:,select_n+1:private_at_copy%n) = real_grid(:,offset:size(real_grid,2):stride)
       call calc_connect(private_at_copy, skip_zero_zero_bonds=.true.)

       local_e = 0.0_dp
       call calc(this, private_at_copy, local_energy=local_e(1:private_at_copy%n), &
            args_str=trim(args_str)//' atom_mask_name=atom_mask source_mask_name=source_mask grid_size='//minval(grid_size))
       do n=1,npoint
          igrid = (n-1)*stride + offset
          pot(lookup(1,igrid),lookup(2,igrid),lookup(3,igrid)) = 2.0_dp*local_e(select_n+n)/at_charge(select_n+n)
       end do
    end do

    call finalise(private_at_copy)
    deallocate(local_e)

    !$omp end parallel
#ifdef _OPENMP
    call omp_set_nested(old_nested)       
#endif

    ! Now do a single calculation of potential and field at ion positions in cluster
    if (cluster%n > (size(at_copy%z)-select_n)) then
       RAISE_ERROR("cluster contains too many atoms!", error)
    end if
    at_copy%n = select_n + cluster%n
    at_copy%nbuffer = select_n + cluster%n
    at_copy%pos(:,select_n+1:at_copy%n) = cluster%pos(:,:)
    if (.not. assign_pointer(cluster, 'charge', cluster_charge)) then
       RAISE_ERROR('calc_electrostatic_potential failed to assign pointer to cluster "charge" property', error)
    end if
    at_charge(select_n+1:at_copy%n) = 1.0_dp
    at_copy%z(select_n+1:at_copy%n) = 0 ! treat as non-interacting dummy atoms
    call calc_connect(at_copy, skip_zero_zero_bonds=.true.)

    call add_property(at_copy, 'es_efield', 0.0_dp, n_cols=3, overwrite=.true., ptr2=es_efield,  error=error)
    PASS_ERROR(error)
    call add_property(at_copy, 'es_pot', 0.0_dp, overwrite=.true., ptr=es_pot, error=error)
    PASS_ERROR(error)
    call add_property(cluster, 'es_efield', 0.0_dp, n_cols=3, overwrite=.true., ptr2=cluster_es_efield, error=error)
    PASS_ERROR(error)
    call add_property(cluster, 'es_pot', 0.0_dp, overwrite=.true., ptr=cluster_es_pot, error=error)
    PASS_ERROR(error)

    call calc(this, at_copy, args_str=trim(args_str)//' efield=es_efield local_energy=es_pot atom_mask_name=atom_mask source_mask_name=source_mask', error=error)
    PASS_ERROR(error)

    ! Electrostatic potential at ion positions is 2*local_energy/charge
    cluster_es_pot = 2.0_dp*es_pot(select_n+1:at_copy%n)/cluster_charge
    cluster_es_efield = es_efield(:,select_n+1:at_copy%n)

    deallocate(select_mask)
    deallocate(lookup,z)
    call finalise(at_copy)
    call system_timer('calc_electrostatic_potential')

  end subroutine calc_electrostatic_potential


end module ElectrostaticEmbed_module
