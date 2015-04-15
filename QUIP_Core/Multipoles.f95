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


#include "error.inc"

module Multipoles_module

use error_module
!use system_module, only : dp, print, inoutput, optional_default, system_timer, operator(//), print_warning
use system_module
use units_module
use dictionary_module
use paramreader_module
use periodictable_module
use linearalgebra_module
use atoms_types_module
use atoms_module
use topology_module
use extendable_str_module  
use table_module           
use connection_module           
use clusters_module        
use structures_module      
use multipole_interactions_module
use partridge_schwenke_dipole_module


#ifdef HAVE_GAP
use descriptors_module
#endif

use mpi_context_module
use QUIP_Common_module


implicit none
private
public :: multipole_sites_setup, atomic_forces_from_sites, electrostatics_calc,finalise, build_polarisation_matrix, calc_induced_dipoles,ewald_setup

integer, parameter :: Multipole_Position_Atomic = 1
integer, parameter :: Multipole_Position_Centre_of_Mass = 2
integer, parameter :: Multipole_Position_M_Site = 3

integer, parameter :: Charge_Method_None = 0
integer, parameter :: Charge_Method_Fixed = 1
integer, parameter :: Charge_Method_GAP = 2
integer, parameter :: Charge_Method_Partridge_Schwenke = 3

integer, parameter :: Dipole_Method_None = 0
integer, parameter :: Dipole_Method_Partridge_Schwenke = 1
integer, parameter :: Dipole_Method_GAP = 2

integer, parameter :: Polarisation_Method_None = 0
integer, parameter :: Polarisation_Method_FPI = 1
integer, parameter :: Polarisation_Method_GMRES = 2
integer, parameter :: Polarisation_Method_QR = 3

integer, parameter :: Damping_None = 0
integer, parameter :: Damping_Exp = 1
integer, parameter :: Damping_Erf = 2
integer, parameter :: Damping_Erf_Uniform = 3

integer, parameter :: Screening_None = 0
integer, parameter :: Screening_Yukawa = 1
integer, parameter :: Screening_Erfc_Uniform = 2

real(dp), parameter :: reciprocal_time_by_real_time = 1.0_dp / 3.0_dp

type Monomers
  type(Multipole_Interactions_Site), dimension(:), allocatable :: site_types ! this just a dummy list of site types with no actual positions
  integer, dimension(:), allocatable :: signature
  integer, dimension(:,:), allocatable :: excluded_pairs,monomer_indices
  real(dp), dimension(:), allocatable :: masses
  real(dp) :: monomer_cutoff = 0.0_dp, step ! step size (Angstroms) used for finite difference gradient calculations
  real(dp) :: gammaM = 0.426706882_dp ! param for determining M-site position in water models
end type Monomers

public :: Multipole_Moments
type Multipole_Moments
  real(dp) :: cutoff = 0.0_dp, dipole_tolerance=1e-7_dp, polarisation_cutoff=0.0_dp
  type(Monomers), dimension(:), allocatable :: monomer_types 
  type(Multipole_Calc_Opts) :: calc_opts
  type(Multipole_Interactions_Site), dimension(:), allocatable :: sites ! 1 to 1 correspondence w dummy atom sites
  logical ::  strict=.true., initialised=.false.,intermolecular_only = .true.

  integer, dimension(:), allocatable :: polarisable_map
  logical, dimension(:), allocatable :: site_polarisable
  integer, dimension(:,:), allocatable :: exclude_list
end type Multipole_Moments

public :: Ewald_arrays
type Ewald_arrays
  real(dp), dimension(:,:), allocatable :: Q
  logical :: do_ewald = .false.
end type Ewald_arrays

interface Finalise
  module procedure Multipole_Moments_Finalise
end interface Finalise

contains


subroutine Multipole_Moments_Finalise(this)

  type(Multipole_Moments),intent(inout) :: this

  if (allocated(this%monomer_types)) deallocate(this%monomer_types)
  if (allocated(this%polarisable_map)) deallocate(this%polarisable_map)
  if (allocated(this%site_polarisable)) deallocate(this%site_polarisable)
  this%initialised = .False.

end subroutine Multipole_Moments_Finalise 

subroutine Monomers_Finalise(this)
  type(Monomers),intent(inout) :: this
  this%monomer_cutoff = 0.0_dp

  ! deallocate arrays
end subroutine Monomers_Finalise

subroutine multipole_sites_setup(at,multipoles,dummy_atoms,do_grads,strict)
  type(Atoms)             :: at,dummy_atoms
  type(Multipole_Moments) :: multipoles
  logical, optional       :: do_grads,strict

  type(Monomers), dimension(:), allocatable :: monomer_lists
  logical,dimension(:), allocatable :: associated_to_monomer
  integer, dimension(3) :: shift
  integer :: n_sites, n_mono_types, position,i,j
  logical :: my_do_grads,my_strict

  

  my_do_grads=optional_default(.false.,do_grads)
  my_strict=optional_default(.true.,strict)
  allocate(associated_to_monomer(at%N))
  associated_to_monomer = .false.

  n_sites=0
  do i=1,n_mono_types     ! can introcduce a case switch here if want to support more monomer finding routines
    if(allocated(monomer_lists(i)%monomer_indices))deallocate(monomer_lists(i)%monomer_indices)
    call find_general_monomer(at,monomer_lists(i)%monomer_indices,monomer_lists(i)%signature,associated_to_monomer,monomer_lists(i)%monomer_cutoff)
    n_sites = n_sites + size(monomer_lists(i)%site_types) * size(monomer_lists(i)%monomer_indices,2)
  end do
!  if (.not. all(associated_to_monomer)) then
!    RAISE_ERROR('Multipoles: Not all atoms are assigned to a monomer', error)     
!  end if 
  allocate(multipoles%sites(n_sites))

  position=0
  do i=1,n_mono_types    ! fill in list of sites and exclude list
    do j=1,size(multipoles%monomer_types(i)%monomer_indices,2)
      call add_sites_for_monomer(at,multipoles%sites,position,multipoles%exclude_list,multipoles%monomer_types(i),multipoles%monomer_types(i)%monomer_indices(:,j),my_do_grads)
    end do
  end do

  call initialise(dummy_atoms,0,at%lattice)   ! set up dummy atoms obj
  do i=1,n_sites
    call add_atoms(dummy_atoms,multipoles%sites(i)%position,1) ! all dummy atoms are Hydrogens
  end do
  call set_cutoff(dummy_atoms,multipoles%cutoff)
  call calc_connect(dummy_atoms)

  ! if want to implement this way for efficiency, will probably have to create a second Connection object, which doesn't contain excluded interactions
  ! this is because for calculated induced dipoles we typically don't ignore intramolecular interactions, so we shuold ignore the exclude list.
!!$
!!$  do i=1,size(exclude_list,2)   ! delete bonds in exclude list
!!$    call distance_min_image(dummy_atoms,exclude_list(1,i),exclude_list(2,i),shift=shift)
!!$    call remove_bond(dummy_atoms%connection,exclude_list(1,i),exclude_list(2,i),shift=shift)
!!$  end do

  deallocate(associated_to_monomer)

end subroutine multipole_sites_setup

subroutine ewald_setup(dummy_atoms,multipoles,ewald,my_ewald_error)
  type(Atoms) :: dummy_atoms
  type(Multipole_Moments) :: multipoles
  type(Ewald_arrays) :: ewald

  real(dp) :: r_ij, erfc_ar, arg, my_ewald_error, alpha, kmax, kmax2, prefac, infac, two_alpha_over_sqrt_pi, v, &
   & ewald_precision, ewald_cutoff, my_cutoff, my_smooth_coulomb_cutoff, smooth_arg, smooth_f, dsmooth_f

  ewald_precision = -log(my_ewald_error)
  ewald_cutoff = sqrt(ewald_precision/PI) * reciprocal_time_by_real_time**(1.0_dp/6.0_dp) * &
  & minval(sqrt( sum(dummy_atoms%lattice(:,:)**2,dim=1) )) / dummy_atoms%N**(1.0_dp/6.0_dp)
  call print('Ewald cutoff = '//ewald_cutoff,PRINT_ANAL)
  multipoles%cutoff = ewald_cutoff

  write(*,*) "more to do to fill these arrays"

end subroutine ewald_setup

subroutine add_sites_for_monomer(at,sites,offset,exclude_list,monomer_type,atom_indices,do_grads)
  type(Atoms) :: at
  type(Multipole_Interactions_Site), dimension(:) :: sites
  integer :: offset
  integer, dimension(:,:), allocatable :: exclude_list
  type(Monomers) :: monomer_type
  integer, dimension(:) :: atom_indices
  logical :: do_grads
  
  real(dp),dimension(:,:),allocatable :: atomic_positions
  integer :: monomer_size,sites_per_mono,exclude_offset,n_exclude
  integer :: i_atom_site,i_site_glob,i_atom,i_site,i
  !integer, dimension(3) :: water_signature = (/1,1,8/)
  real(dp), dimension(3,3) :: identity3x3
  real(dp), dimension(3)   :: com_pos

  real(dp) :: monomer_mass
  integer,dimension(:), allocatable :: signature_copy,site_atom_map


  identity3x3 = 0.0_dp
  call add_identity(identity3x3)
  
  monomer_size = size(monomer_type%signature)
  sites_per_mono = size(monomer_type%site_types)

  allocate(atomic_positions(3,monomer_size))
  do i=1,monomer_size
    monomer_type%masses(i) =ElementMass(monomer_type%signature(i))
  end do
  monomer_mass = sum(monomer_type%masses)

  allocate(site_atom_map(sites_per_mono))
  allocate(signature_copy(monomer_size))
  signature_copy = monomer_type%signature
  site_atom_map=0

  do i_site=1,sites_per_mono
    if (monomer_type%site_types(i_site)%pos_type .eq. Multipole_Position_Atomic) then 
      i_atom_site=find_in_array(signature_copy,monomer_type%site_types(i_site)%atomic_number)
      signature_copy(i_atom_site)=0 ! make sure we don't put two sites on the same atom
      site_atom_map(i_site) =i_atom_site 
    end if
  end do

  atomic_positions=0.0_dp
  com_pos=0.0_dp

  ! find positions of atoms and centre of mass 
  do i=1,monomer_size
    i_atom = atom_indices(i)
    atomic_positions(:,i) = at%pos(:,i_atom)
    com_pos = com_pos + monomer_type%masses(i) * at%pos(:,i_atom)
  end do

  com_pos = com_pos * (1.0_dp/monomer_mass)

  do i_site=1,sites_per_mono 
    i_site_glob=offset+i_site

    sites(i_site_glob) = monomer_type%site_types(i_site) ! copy info from template site

    sites(i_site_glob)%charge = 0.0_dp
    sites(i_site_glob)%dipole = 0.0_dp

    allocate(sites(i_site_glob)%atom_indices(monomer_size))
    sites(i_site_glob)%atom_indices = atom_indices 

    selectcase(monomer_type%site_types(i_site)%pos_type) ! assign positions
      case(Multipole_Position_Centre_of_Mass)
        sites(i_site_glob)%position = com_pos
      case(Multipole_Position_Atomic)
        sites(i_site_glob)%position = atomic_positions(:,site_atom_map(i_site))
      case(Multipole_Position_M_Site)
        call m_site_position(atomic_positions,monomer_type%gammaM,sites(i_site_glob)%position)
      case default
        call print("site position param : "//sites(i_site_glob)%pos_type)
        call system_abort("Multipole_Moments_Assign doesn't know where to put this multipole moment")
    end select


    selectcase(sites(i_site_glob)%charge_method) ! assign charges
      case(Charge_Method_Fixed)
        sites(i_site_glob)%charge = monomer_type%site_types(i_site)%charge           
      case(Charge_Method_GAP)
          call system_abort("GAP moments not yet implemented")  ! call charge_GAP()
      case(Charge_Method_Partridge_Schwenke)
!#ifndef HAVE_FX
!  RAISE_ERROR('Multipoles: PS water as charges requested but FX model was not compiled in. Check the HAVE_FX flag in the Makefiles.', error)
!#endif
        call charges_PS(atomic_positions,sites(i_site_glob)%charge,i_site,monomer_type%gammaM)
        call test_charge_grads_PS(atomic_positions,monomer_type%gammaM)
    end select
    selectcase(sites(i_site_glob)%dipole_method) ! assign dipoles
      case(Dipole_Method_GAP)
        continue     ! call dipole_moment_GAP()
      case(Dipole_Method_Partridge_Schwenke)
        call dipole_moment_PS(atomic_positions,sites(i_site_glob)%dipole) 
    end select

    if (do_grads) then

      allocate(sites(i_site_glob)%pos_grad_positions(3,3,monomer_size)) !  gradients of position of site wrt atomic positions
      allocate(sites(i_site_glob)%charge_grad_positions(1,3,monomer_size)) !  gradients of charge wrt atomic positions
      allocate(sites(i_site_glob)%dipole_grad_positions(3,3,monomer_size)) !  gradients of dipole components wrt atomic positions

      sites(i_site_glob)%e_grad_pos = 0.0_dp 
      sites(i_site_glob)%e_grad_charge = 0.0_dp
      sites(i_site_glob)%e_grad_dipole = 0.0_dp

      sites(i_site_glob)%charge_grad_positions = 0.0_dp 
      sites(i_site_glob)%dipole_grad_positions = 0.0_dp 
      sites(i_site_glob)%pos_grad_positions=0.0_dp

      ! calc gradient of this site's multipole components with respect to atomic positions
      selectcase(sites(i_site_glob)%charge_method) ! assign charge gradients
        case(Charge_Method_Fixed)
          continue ! gradients are zero
        case(Charge_Method_GAP)
          call system_abort("GAP moments not yet implemented")
          ! something like call charge_gradients_GAP()
        case(Charge_Method_Partridge_Schwenke)
          call charge_gradients_PS(atomic_positions,sites(i_site_glob)%charge_grad_positions,monomer_type%gammaM,i_site)
      end select
      selectcase(sites(i_site_glob)%dipole_method) 
        case(Dipole_Method_Partridge_Schwenke)
          call dipole_moment_gradients_PS(atomic_positions,sites(i_site_glob)%dipole_grad_positions,step=monomer_type%step) 
        case(Dipole_Method_GAP)
          call system_abort("GAP moments not yet implemented")
          ! something like call dipole_gradients_GAP()
      end select

      ! calc gradients of this site's position with respect to atomic positions
      selectcase(sites(i_site_glob)%pos_type) 
        case(Multipole_Position_Centre_of_Mass)
          do i=1,monomer_size
            sites(i_site_glob)%pos_grad_positions(:,:,i) =  ( monomer_type%masses(i) / monomer_mass) * identity3x3
          end do
        case(Multipole_Position_Atomic)
          sites(i_site_glob)%pos_grad_positions(:,:,site_atom_map(i_site)) = identity3x3
        case(Multipole_Position_M_Site)
          call m_site_position_grads(monomer_type%gammaM,sites(i_site_glob)%pos_grad_positions)
        case default
          call system_abort("Multipole_Moments_Assign doesn't know how multipole position varies with atomic positions")
      end select
    end if

  end do

  n_exclude=size(monomer_type%excluded_pairs,2)
  exclude_offset=size(exclude_list,2)
  call reallocate(exclude_list,2,exclude_offset+n_exclude,copy=.true.)
  do i=1,n_exclude
    exclude_list(1,exclude_offset+i)=offset+monomer_type%excluded_pairs(1,i)
    exclude_list(2,exclude_offset+i)=offset+monomer_type%excluded_pairs(2,i)    
  end do

end subroutine add_sites_for_monomer

subroutine electrostatics_calc(at,multipoles,ewald,do_pot,do_field,do_force,e)
  type(Atoms) :: at
  type(Multipole_Moments) :: multipoles
  type(Ewald_arrays) :: ewald
  logical,optional :: do_pot,do_field,do_force
  real(dp),optional :: e

  real(dp) :: my_energy,r_ij,site_site_energy
  integer :: i_site,n,j_site
  real(dp), dimension(3) :: diff
  type(Multipole_Interactions_Site) :: site1,site2


  multipoles%calc_opts%do_pot = optional_default(.false.,do_pot)
  multipoles%calc_opts%do_field = optional_default(.false.,do_field)
  multipoles%calc_opts%do_force = optional_default(.false.,do_force)


  ! real space part
  do i_site=1,at%N
    site1 = multipoles%sites(i_site)
    do n = 1, n_neighbours(at,i_site)
      j_site = neighbour(at,i_site,n,distance=r_ij,diff=diff)

      if( r_ij > multipoles%cutoff )  cycle
      if ( find_in_array(multipoles%exclude_list,(/i_site,j_site/)) + find_in_array(multipoles%exclude_list,(/j_site,i_site/)) .gt. 0 ) cycle 

      site2 = multipoles%sites(j_site)
      site2%position = site1%position + diff         ! make a copy of neighbour site and move it to the correct position

      call Multipole_Moments_Site_Site_Interaction(site_site_energy,site1,site2, multipoles%calc_opts, cutoff=multipoles%cutoff)
      my_energy = my_energy + site_site_energy
    end do
  end do

  ! if ewald
  ! recip part
  ! self part

  if (present(e))e=e+my_energy
   
end subroutine electrostatics_calc


subroutine atomic_forces_from_sites(sites,f)
  type(Multipole_Interactions_Site), dimension(:) :: sites
  real(dp), intent(out), dimension(:,:),optional :: f

  integer :: i, j, k, a, n_atoms,i_site,i_atom


  do i_site=1,size(sites)  ! add up force contributions in most transparent way possible

    n_atoms = size(sites(i_site)%atom_indices)

    do a=1,n_atoms ! a is atom number, i component of atomic position, j component of dipole, k component of site position

      i_atom=sites(i_site)%atom_indices(a)
      do i=1,3
        f(i,i_atom) = f(i,i_atom) - sites(i_site)%e_grad_charge * sites(i_site)%charge_grad_positions(1,i,a) ! NB force is minus grad

        do j=1,3
          f(i,i_atom) = f(i,i_atom) - sites(i_site)%e_grad_dipole(j) * sites(i_site)%dipole_grad_positions(j,i,a) 
        end do

        do k=1,3
          f(i,i_atom) = f(i,i_atom) - sites(i_site)%e_grad_pos(k) * sites(i_site)%pos_grad_positions(k,i,a)
        end do
      end do

    end do

  end do

end subroutine atomic_forces_from_sites

subroutine build_polarisation_matrix(at,multipoles,pol_matrix,ewald)
  type(Atoms) :: at
  type(Multipole_Moments) :: multipoles
  real(dp), dimension(:,:) :: pol_matrix 
  type(Ewald_arrays) :: ewald

  write(*,*) "not implemented yet"
end subroutine build_polarisation_matrix

subroutine  calc_induced_dipoles(pol_matrix,multipoles,polarisation)
  real(dp), dimension(:,:) :: pol_matrix
  type(Multipole_Moments) :: multipoles
  integer :: polarisation
  write(*,*) "not implemented yet" 
end subroutine  calc_induced_dipoles

end module Multipoles_module
