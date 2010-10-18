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

  use libatoms_module
  use QUIP_Common_module
  use Functions_module
  use Potential_module

  implicit none
  private

  integer, private, parameter :: nspins = 2

  public :: calc_electrostatic_potential, write_electrostatic_potential, &
       calc_electric_field, write_electric_field, make_periodic_potential

contains

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
                  & real(nz-1,dp)*lattice(3,1)/real(ngrid(3),dp) + origin(1) - extent(1)/2.0_dp
             r_real_grid(2,point) = real(nx-1,dp)*lattice(1,2)/real(ngrid(1),dp) +&
                  & real(ny-1,dp)*lattice(2,2)/real(ngrid(2),dp) +&
                  & real(nz-1,dp)*lattice(3,2)/real(ngrid(3),dp) + origin(2) - extent(2)/2.0_dp
             r_real_grid(3,point) = real(nx-1,dp)*lattice(1,3)/real(ngrid(1),dp) +&
                  & real(ny-1,dp)*lattice(2,3)/real(ngrid(2),dp) +&
                  & real(nz-1,dp)*lattice(3,3)/real(ngrid(3),dp) + origin(3) - extent(3)/2.0_dp

          end do
       end do
    end do

    grid_size(1) = extent(1)/real(ngrid(1), dp)
    grid_size(2) = extent(2)/real(ngrid(2), dp)
    grid_size(3) = extent(3)/real(ngrid(3), dp)

  end subroutine assign_grid_coordinates

  function gaussian_charge_pot_func(r, s, grid_size) result(V)
    real(dp) :: r, s, grid_size(3)
    real(dp) :: V

    real(dp), parameter :: sqrt_2 = 1.41421356237309504880168872421_dp

    ! Truncate r so smallest value is one grid point
    r = max(r, maxval(grid_size))

    ! point charge part
    V = 1.0_dp/r
    if (s > 0.0_dp .and. r < 9.0_dp*s) then ! screening if width > 0, correction for r >= 9s is < 1e-16
       V = V * erf(r/(sqrt_2*s))
    endif
  end function gaussian_charge_pot_func

  function gaussian_charge_pot_force_func(r, rsq, s, grid_size) result(Vf)
    real(dp) :: r(3), rsq, s, grid_size(3)
    real(dp) :: Vf(3)

    real(dp), parameter :: sqrt_2 = 1.41421356237309504880168872421_dp
    real(dp), parameter :: sqrt_pi = 1.77245385090551602729816748334_dp
    real(dp) :: r_mag, V_pt, erf_val, erf_deriv(3)
    integer :: i

    ! V = 1/r erf(r/(sqrt_2 s))
    ! derf(x)/dx = (2/sqrt(pi)) exp(-x^2)

    do i=1,3
       r(i) = max(r(i), grid_size(i))
    end do
    rsq = max(r(1)*r(1)+r(2)*r(2)+r(3)*r(3), rsq)

    ! point charge part
    Vf = 1.0_dp/(rsq*sqrt(rsq)) * r
    if (s > 0.0_dp .and. rsq < 81.0_dp*s*s) then ! screening if width > 0, correction for r >= 9s is <= 1e-16
       r_mag = sqrt(rsq)
       V_pt = 1.0_dp/r_mag
       erf_val = erf(r_mag/(sqrt_2 * s))
       erf_deriv = -(sqrt_2/(sqrt_pi*s)) * exp(-rsq/(2.0*s*s)) * r/r_mag
       Vf = Vf * erf_val + V_pt * erf_deriv
    endif
  end function gaussian_charge_pot_force_func

  subroutine write_electrostatic_header(out, lattice, ngrid)
    type(InOutput), intent(inout) :: out
    real(dp), intent(in) :: lattice(3,3)
    integer, intent(in) :: ngrid(3)
    real(dp) :: a, b, c, alpha, beta, gamma

    call get_lattice_params(lattice, a, b, c, alpha, beta, gamma)

    call print('BEGIN header', file=out)
    call print('', file=out)
    call print('           Real Lattice(A)               Lattice parameters(A)    Cell Angles', file=out)
    write (out%unit, '(3f12.7,5x,"a =",f12.6,2x,"alpha =",f12.6)') lattice(:,1), a, alpha
    write (out%unit, '(3f12.7,5x,"b =",f12.6,2x,"beta  =",f12.6)') lattice(:,2), b, beta
    write (out%unit, '(3f12.7,5x,"c =",f12.6,2x,"gamma =",f12.6)') lattice(:,3), c, gamma
    call print('', file=out)
    write(out%unit,'(i4,T30,a)') nspins,'   ! nspins'
    write(out%unit,'(3(i4,2x),T30,a)') ngrid, '   ! fine FFT grid along <a,b,c>'

  end subroutine write_electrostatic_header

  subroutine write_electric_field(at, charge, filename, ngrid, extent, efield, atom_mask_name)
    type(Atoms), intent(inout) :: at
    real(dp), intent(in) :: charge(:)
    character(len=*), intent(in) :: filename, atom_mask_name
    integer, intent(in) :: ngrid(3)
    real(dp), intent(in) :: extent(3)
    real(dp), dimension(:,:), intent(in) :: efield

    logical, dimension(:), pointer :: atom_mask
    type(InOutput) :: out
    integer i
    real(dp) :: lattice(3,3)

    lattice(:,:) = 0.0_dp
    lattice(1,1) = extent(1)
    lattice(2,2) = extent(2)
    lattice(3,3) = extent(3)

    call assign_property_pointer(at, atom_mask_name, atom_mask)

    ! Write electric field to file
    call initialise(out, filename, OUTPUT)
    call write_electrostatic_header(out, lattice, ngrid)
    call print('END header: data is "sp nsp charge efield_x efield_y efield_z"', file=out)
    do i=1,at%N
       if (.not. atom_mask(i)) cycle
       write (out%unit,'(a6,i6,4f12.6)'), a2s(at%species(:,i)), index_to_z_index(at, i), charge(i), efield(:,i)
    end do
    call finalise(out)

  end subroutine write_electric_field

  subroutine write_electrostatic_potential(at, filename, ngrid, extent, pot)
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ngrid(3)
    real(dp), intent(in) :: extent(3)
    real(dp), dimension(:,:,:), intent(in) :: pot

    type(InOutput) :: out
    integer igrid, spin, i, j, k
    real(dp) :: lattice(3,3)

    lattice(:,:) = 0.0_dp
    lattice(1,1) = extent(1)
    lattice(2,2) = extent(2)
    lattice(3,3) = extent(3)

    call initialise(out, filename, OUTPUT)
    call write_electrostatic_header(out, lattice, ngrid)
    call print('END header: data is "i j k pot"', file=out)  
    do spin=1,nspins
       igrid = 0
       do k=1,ngrid(3)
          do j=1,ngrid(2)
             do i=1,ngrid(1)
                igrid = igrid + 1
                write (out%unit,'(3i6,2f20.6)') i, j, k, pot(i,j,k)
             end do
          end do
       end do
    end do
    call finalise(out)

  end subroutine write_electrostatic_potential

  subroutine make_periodic_potential(at, real_grid, ngrid, is_periodic, cutoff_radius, cutoff_width, mark_name, pot, error)
    type(Atoms), intent(inout) :: at
    real(dp), dimension(:,:), intent(in) :: real_grid
    integer, intent(in) :: ngrid(3)
    logical, intent(in) :: is_periodic(3)
    real(dp), intent(in) :: cutoff_radius, cutoff_width
    character(len=*), intent(in) :: mark_name
    real(dp), intent(inout) :: pot(:,:,:)
    integer, optional, intent(out) :: error

    integer i,j,k,a,igrid, nsurface, cell_image_na, cell_image_nb, cell_image_nc
    logical save_is_periodic(3), dir_mask(3)
    real(dp) :: cutoff, surface_avg, f, d(3), fc, dfc, masked_d
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

  !% Evaluate electric field on real atoms, i.e. those not marked with HYBRID_ELECTROSTATIC_MARK
  subroutine calc_electric_field(this, at, mark_name, args_str, error)
    type(Potential), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in) :: args_str, mark_name
    integer, optional, intent(out) :: error

    integer, dimension(:), pointer :: mark
    logical, dimension(:), pointer :: atom_mask, source_mask

    INIT_ERROR(error)

    call assign_property_pointer(at, mark_name, mark, error)
    PASS_ERROR(error)

    call add_property(at, 'atom_mask', .false., overwrite=.true., ptr=atom_mask, error=error)
    PASS_ERROR(error)
    atom_mask = mark /= HYBRID_ELECTROSTATIC_MARK

    call add_property(at, 'source_mask', .false., overwrite=.true., ptr=source_mask, error=error)
    PASS_ERROR(error)
    source_mask = mark == HYBRID_ELECTROSTATIC_MARK

    call add_property(at, 'es_efield', 0.0_dp, overwrite=.true., error=error)
    PASS_ERROR(error)

    call calc(this, at, args_str=trim(args_str)//' efield=es_efield atom_mask_name=atom_mask source_mask_name=source_mask', error=error)
    PASS_ERROR(error)

  end subroutine calc_electric_field

  !% Evaluate electrostatic potential on a grid by adding test atoms at grid points
  subroutine calc_electrostatic_potential(this, at, mark_name, ngrid, origin, extent, real_grid, pot, args_str, error)
    type(Potential), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in) :: args_str, mark_name
    integer, intent(in) :: ngrid(3)
    real(dp), intent(in) :: origin(3), extent(3)
    real(dp), intent(inout) :: real_grid(:,:), pot(:,:,:)
    integer, optional, intent(out) :: error

    type(Atoms) :: at_copy
    type(Table) :: remove_list
    real(dp) :: r, charge, sigma, grid_size(3)
    real(dp), pointer :: at_charge(:)
    real(dp), allocatable :: local_e(:)
    integer, pointer, dimension(:) :: mark
    integer, allocatable, dimension(:) :: z
    integer, allocatable, dimension(:,:) :: lookup
    logical, pointer, dimension(:) :: atom_mask, source_mask
    integer i, j, k, n, igrid, offset, stride, npoint

    INIT_ERROR(error)
    call system_timer('calc_electrostatic_potential')

    call assign_grid_coordinates(ngrid, origin, extent, grid_size, real_grid)

    ! Make copy of atoms, so we can add test atoms at grid points
    call atoms_copy_without_connect(at_copy, at)
    call set_cutoff(at_copy, cutoff(this))

    call assign_property_pointer(at_copy, mark_name, mark, error)
    PASS_ERROR(error)

    call assign_property_pointer(at_copy, 'charge', at_charge, error)
    PASS_ERROR(error)

    ! choose balance between number of evaluations and number of simulataneous grid points
    npoint = 150
    stride = size(real_grid, 2)/npoint
    npoint = size(real_grid(:,1:size(real_grid,2):stride), 2)
    call print('npoint = '//npoint//' stride = '//stride)

    allocate(z(npoint))
    z(:) = 0
    call add_atoms(at_copy, real_grid(:,1:size(real_grid,2):stride), z)

    allocate(local_e(at_copy%n))

    ! added atoms so must reassign pointers
    if (.not. assign_pointer(at_copy, mark_name, mark)) then
       RAISE_ERROR('IPModel_ASAP2_calc failed to assign pointer to "'//trim(mark_name)//'" property', error)
    end if
    if (.not. assign_pointer(at_copy, 'charge', at_charge)) then
       RAISE_ERROR('IPModel_ASAP2_calc failed to assign pointer to "charge" property', error)
    endif
    at_charge(at%n+1:at_copy%n) = 1.0_dp

    ! atom_mask is true for atoms for which we want the local potential energy, i.e dummy atoms
    call add_property(at_copy, 'atom_mask', .false., overwrite=.true., ptr=atom_mask, error=error)
    PASS_ERROR(error)
    atom_mask = .false. 
    atom_mask(at%n+1:at_copy%n) = .true.

    call print('atom_mask='//atom_mask)

    ! source_mask is true for atoms which should act as electrostatic sources
    call add_property(at_copy, 'source_mask', .false., overwrite=.true., ptr=source_mask, error=error)
    PASS_ERROR(error)
    source_mask(1:at%n) = mark(1:at%n) == HYBRID_ELECTROSTATIC_MARK
    source_mask(at%n+1:at_copy%n) = .false.

    call print('source_mask='//source_mask)

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

    do offset=1,stride
       npoint = size(real_grid(:,offset:size(real_grid,2):stride),2)
       at_copy%n = at%n + npoint
       at_copy%nbuffer = at%nbuffer + npoint
       at_copy%pos(:,at%n+1:at_copy%n) = real_grid(:,offset:size(real_grid,2):stride)

       !call write(at_copy, 'embed.'//offset//'.xyz')

       call calc_connect(at_copy, skip_zero_zero_bonds=.true.)

       local_e = 0.0_dp
       call calc(this, at_copy, local_energy=local_e(1:at_copy%n), &
            args_str=trim(args_str)//' atom_mask_name=atom_mask source_mask_name=source_mask grid_size='//minval(grid_size))
       
       do n=1,npoint
          igrid = (n-1)*stride + offset
          pot(lookup(1,igrid),lookup(2,igrid),lookup(3,igrid)) = 2.0_dp*local_e(at%n+n)/at_charge(at%n+n)/HARTREE
       end do
    end do

    deallocate(lookup)
    call finalise(at_copy)
    deallocate(local_e, z)
    call system_timer('calc_electrostatic_potential')

  end subroutine calc_electrostatic_potential

end module ElectrostaticEmbed_module
