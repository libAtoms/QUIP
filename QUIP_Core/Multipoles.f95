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
use system_module, only : dp, print, inoutput, optional_default, system_timer, operator(//), print_warning
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
public :: Multipole_Moments_Pre_Calc, Multipole_Moments_Calc,Multipole_Moments_Make_Dummy_Atoms,Multipole_Moments_Forces_From_Dummy_Atoms,Finalise

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


type Monomers
  type(Multipole_Interactions_Site), dimension(:), allocatable :: site_types ! this just a dummy list of site types with no actual positions
  integer, dimension(:), allocatable :: signature
  integer, dimension(:,:) allocatable :: excluded_pairs
  real(dp), dimension(:), allocatable :: masses
  real(dp) :: monomer_cutoff = 0.0_dp, step ! step size (Angstroms) used for finite difference gradient calculations
  real(dp) :: gammaM = 0.426706882_dp ! param for determining M-site position in water models
end type Monomer_List

public :: Multipole_Moments
type Multipole_Moments
  real(dp) :: cutoff = 0.0_dp, dipole_tolerance=1e-7_dp, polarisation_cutoff=0.0_dp
  type(Monomers), dimension(:), allocatable :: monomer_types 
  type(Multipole_Calc_Opts) :: calc_opts
  type(Multipole_Interactions_Site), dimension(:), allocatable :: sites ! 1 to 1 correspondence w dummy atom sites
  logical ::  strict=.true., initialised=.false.
  character(len=STRING_LENGTH) label
  integer :: n_monomer_types,polarisation
  integer, dimension(:), allocatable :: polarisable_map
  logical, dimension(:), allocatable :: site_polarisable
end type Multipole_Moments


interface Finalise
  module procedure Multipole_Moments_Finalise
end interface Finalise

contains


subroutine Multipole_Moments_Finalise(this)

  type(Multipole_Moments),intent(inout) :: this
  this%label=''
  this%n_monomer_types = 0
  if (allocated(this%monomer_types)) deallocate(this%monomer_types)
  if (allocated(this%polarisable_map)) deallocate(this%polarisable_map)
  if (allocated(this%site_polarisable)) deallocate(this%site_polarisable)
  this%initialised = .False.

end subroutine Multipole_Moments_Finalise 

subroutine Monomer_List_Finalise(this)
  type(Monomer_list),intent(inout) :: this
  this%monomer_cutoff = 0.0_dp
  this%n_monomers = 0

  if (allocated(this%sites)) deallocate(this%sites)
  if (allocated(this%signature)) deallocate(this%signature)
  if (allocated(this%masses)) deallocate(this%signature)
  if (allocated(this%monomer_indices)) deallocate(this%monomer_indices)

end subroutine Monomer_List_Finalise

subroutine multipole_sites_setup(at,multipoles,dummy_atoms)
  type(Atoms)             :: at,dummy_atoms
  type(Multipole_Moments) :: multipoles
  !
  type(Monomers), dimension(:), allocatable :: monomer_lists
  integer,dimension(:,:) allocatable :: idx
  integer,dimension(:) allocatable :: associated_to_monomer
  ! want to end up with:
  !   a list of sites
  !   dummy atoms object with correct connectivity
  !   
  ! Topology first
  ! loop over monomer types
  ! figure out number of total sites
  ! allocate sites array
  ! find all the monomers in the Atoms object and fill out the monomer lists in mono_lists
  n_sites=0
  do i=1,n_mono_types
    if(allocated(monomer_lists(i)%monomer_indices))deallocate(monomer_lists(i)%monomer_indices)
    call find_general_monomer(at,monomer_lists(i)%monomer_indices,monomer_lists(i)%signature,associated_to_monomer,monomer_lists(i)%monomer_cutoff,atom_ordercheck,use_smooth_cutoff,error)
    n_sites = n_site + size(monomer_lists(i)%site_types) * size(monomer_lists(i)%monomer_indices,2)
  end do

  allocate(multipoles%sites(n_sites))
  offset=0
  do i=1,n_mono_types  
    do j=1,size(monomer_lists(i)%monomer_indices,2)
      call add_sites_for_monomer(at,multipoles%sites,exclude_list,offset,this%monomer_types(i),monomer_lists(i)%monomer_indices(:,j))
    end do
  end do
  ! loop over monomer types
  ! fill in sites w correct info
  ! append to exclude list

  ! loop over sites
  ! create atoms in atoms object
  
  ! call calc_connect
  ! loop over exclude list
  ! delete bonds
end subroutine multipole_sites_setup

subroutine Multipole_Moments_Pre_Calc(at, this, intramolecular_factor, e, mpi, cutoff, error)
  type(Atoms), intent(inout) :: at
  type(Multipole_Moments), intent(inout) :: this
  integer ::  intramolecular_factor
  real(dp), optional :: e
  type(MPI_Context), intent(in), optional :: mpi
  real(dp), intent(in), optional :: cutoff 
  integer, intent(out), optional :: error

  type(Monomer_List),dimension(:), allocatable :: mono_lists  ! these hold the actual lists of monomers
  logical, dimension(:), allocatable :: associated_to_monomer
  integer, dimension(:,:), allocatable :: monomer_pairs
  integer :: n_mono_types, i, j
  type(Dictionary) :: params
  logical :: use_smooth_cutoff=.False., atom_ordercheck=.True.

  n_mono_types = size(this%monomer_types)
  allocate(associated_to_monomer(at%N))
  associated_to_monomer=.false.

  ! find all the monomers in the Atoms object and fill out the monomer lists in mono_lists
  do i=1,n_mono_types
    if (allocated(this%monomer_types(i)%monomer_indices) )deallocate(this%monomer_types(i)%monomer_indices)
    call find_general_monomer(at,this%monomer_types(i)%monomer_indices,this%monomer_types(i)%signature,associated_to_monomer,this%monomer_types(i)%monomer_cutoff,atom_ordercheck,use_smooth_cutoff,error)
    this%monomer_types(i)%n_monomers = size(this%monomer_types(i)%monomer_indices,2)
  end do
  ! if topology_error_check try assigning in different order

  if(.not. all(associated_to_monomer)) then
     RAISE_ERROR("Multipole_Moments_Calc: not all atoms assigned to a monomer", error)
  endif
  if (present(e)) e=0.0_dp
  !% specify the atoms in each monomer and a set of multipoles ( according to chosen method )
  !% this also calculates derivatives of the moments' magnitude and position with respect to atomic positions
  !% and pre-calculates interactions between sites on the same molecule

  do i=1,n_mono_types
    call Multipole_Moments_Assign(at, this%monomer_types(i),intramolecular_factor,this%calc_opts, e=e, cutoff=cutoff, error=error)
  end do

end subroutine Multipole_Moments_Pre_Calc

subroutine Multipole_Moments_Assign(at, mono_list, intramolecular_factor, calc_opts,e, cutoff,error)
  type(Atoms), intent(in) :: at
  type(Monomer_List), intent(inout) :: mono_list
  integer :: intramolecular_factor 
  type(Multipole_Calc_Opts) :: calc_opts
  real(dp),optional :: e, cutoff
  integer, optional :: error

  real(dp), dimension(:,:), allocatable :: atomic_positions
  integer :: monomer_size, sites_per_mono, n_monomers, i, j, i_site, j_site, i_atom_site, i_atom, d
  integer, dimension(1) :: unit_array
  integer, dimension(3) :: water_signature = (/1,1,8/)
  integer, dimension(:), allocatable :: signature_copy, site_atom_map
  real(dp), dimension(3) :: ps_charges, com_pos ! position of monomer centre of mass
  real(dp), dimension(3,3,3) :: ps_charge_grads
  real(dp) :: monomer_mass, site_site_energy, my_cutoff
  type(Multipole_Interactions_Site) :: dummy_i, dummy_j
  logical :: do_e

  real(dp), dimension(3,3) :: identity3x3
  real(dp), dimension(3) :: distance_vector

  identity3x3 = 0.0_dp
  call add_identity(identity3x3)
  calc_opts%do_energy = present(e)
  my_cutoff = optional_default(at%cutoff,cutoff)
  
  monomer_size = size(mono_list%signature)
  sites_per_mono = size(mono_list%site_types)

  do i=1,sites_per_mono
    if (mono_list%site_types(i)%dipole_method == Dipole_Method_Partridge_Schwenke .or. mono_list%site_types(i)%charge_method == Charge_Method_Partridge_Schwenke) then
      if (allocated(mono_list%signature)) then
        if(size(mono_list%signature) /= 3 ) then
          call system_abort("invalid signature given for water monomer, cannot use Partridge Schwenke model, signature must be "//water_signature)
        else if (any(mono_list%signature /= water_signature)) then
          call system_abort("signature for water monomer incompatible with Partridge Schwenke, has to be "//water_signature)
        end if
      else
        allocate(mono_list%signature(3))
        mono_list%signature = water_signature
      end if
    end if
  end do

  allocate(atomic_positions(3,monomer_size))
  do i=1,monomer_size
    mono_list%masses(i) =ElementMass(mono_list%signature(i))
  end do
  monomer_mass = sum(mono_list%masses)


  allocate(site_atom_map(sites_per_mono))
  allocate(signature_copy(monomer_size))
  signature_copy = mono_list%signature
  site_atom_map=0

  do i_site=1,sites_per_mono
    if (mono_list%site_types(i_site)%atomic_number > 0) then 
      unit_array = maxloc(signature_copy, signature_copy .eq. mono_list%site_types(i_site)%atomic_number)
      i_atom_site = unit_array(1)
      signature_copy(i_atom_site)=0 ! make sure we don't put two sites on the same atom
      site_atom_map(i_site) =i_atom_site 
    end if
  end do

  if (allocated(mono_list%sites)) deallocate(mono_list%sites)
  allocate(mono_list%sites(sites_per_mono,mono_list%n_monomers))

  do j=1,mono_list%n_monomers
    atomic_positions=0.0_dp
    com_pos=0.0_dp

    ! find positions of atoms and centre of mass !! TODO do this correctly with min img vecs instead of abs positions
    do i=1,monomer_size
      i_atom = mono_list%monomer_indices(i,j)
      atomic_positions(:,i) = at%pos(:,i_atom)
      com_pos = com_pos + mono_list%masses(i) * at%pos(:,i_atom)
    end do

    com_pos = com_pos * (1.0_dp/monomer_mass)

    do i_site=1,sites_per_mono 

      mono_list%sites(i_site,j) = mono_list%site_types(i_site) ! copy info from template site

      selectcase(mono_list%sites(i_site,j)%pos_type) ! assign positions
        case(Multipole_Position_Centre_of_Mass)
          mono_list%sites(i_site,j)%position = com_pos
        case(Multipole_Position_Atomic)
          mono_list%sites(i_site,j)%position = atomic_positions(:,site_atom_map(i_site))
        case(Multipole_Position_M_Site)
          call m_site_position(atomic_positions,mono_list%gammaM,mono_list%sites(i_site,j)%position)
        case default
          call print("site position param : "//mono_list%sites(i_site,j)%pos_type)
          call system_abort("Multipole_Moments_Assign doesn't know where to put this multipole moment")
      end select

      mono_list%sites(i_site,j)%charge = 0.0_dp
      mono_list%sites(i_site,j)%dipole = 0.0_dp
      mono_list%sites(i_site,j)%quadrupole = 0.0_dp

      selectcase(mono_list%sites(i_site,j)%charge_method) ! assign charges
        case(Charge_Method_Fixed)
          mono_list%sites(i_site,j)%charge = mono_list%site_types(i_site)%charge           
        case(Charge_Method_GAP)
            call system_abort("GAP moments not yet implemented")  ! call charge_GAP()
        case(Charge_Method_Partridge_Schwenke)
#ifndef HAVE_FX
  RAISE_ERROR('Multipoles: PS water as charges requested but FX model was not compiled in. Check the HAVE_FX flag in the Makefiles.', error)
#endif
          call charges_PS(atomic_positions,mono_list%sites(i_site,j)%charge,i_site,mono_list%gammaM)
          call test_charge_grads_PS(atomic_positions,mono_list%gammaM)
      end select
      selectcase(mono_list%sites(i_site,j)%dipole_method) ! assign dipoles
        case(Dipole_Method_GAP)
          continue     ! call dipole_moment_GAP()
        case(Dipole_Method_Partridge_Schwenke)
          call dipole_moment_PS(atomic_positions,mono_list%sites(i_site,j)%dipole) 
      end select

      if (calc_opts%do_force) then

        allocate(mono_list%sites(i_site,j)%pos_grad_positions(3,3,monomer_size)) !  gradients of position of site wrt atomic positions
        allocate(mono_list%sites(i_site,j)%charge_grad_positions(1,3,monomer_size)) !  gradients of charge wrt atomic positions
        allocate(mono_list%sites(i_site,j)%dipole_grad_positions(3,3,monomer_size)) !  gradients of dipole components wrt atomic positions

        mono_list%sites(i_site,j)%e_grad_pos = 0.0_dp ! just initialising these to zero
        mono_list%sites(i_site,j)%e_grad_charge = 0.0_dp
        mono_list%sites(i_site,j)%e_grad_dipole = 0.0_dp

        mono_list%sites(i_site,j)%charge_grad_positions = 0.0_dp 
        mono_list%sites(i_site,j)%dipole_grad_positions = 0.0_dp 
        mono_list%sites(i_site,j)%pos_grad_positions=0.0_dp

      ! calc gradient of this site's multipole components with respect to atomic positions
        selectcase(mono_list%sites(i_site,j)%charge_method) ! assign charge gradients
          case(Charge_Method_Fixed)
            continue ! gradients are zero
          case(Charge_Method_GAP)
            call system_abort("GAP moments not yet implemented")
            ! something like call charge_gradients_GAP()
          case(Charge_Method_Partridge_Schwenke)
            call charge_gradients_PS(atomic_positions,mono_list%sites(i_site,j)%charge_grad_positions,mono_list%gammaM,i_site)
        end select
        selectcase(mono_list%sites(i_site,j)%dipole_method) 
          case(Dipole_Method_Partridge_Schwenke)
            call dipole_moment_gradients_PS(atomic_positions,mono_list%sites(i_site,j)%dipole_grad_positions,step=mono_list%step) 
          case(Dipole_Method_GAP)
            call system_abort("GAP moments not yet implemented")
            ! something like call dipole_gradients_GAP()
        end select

     ! calc gradients of this site's position with respect to atomic positions
        selectcase(mono_list%sites(i_site,j)%pos_type) 
          case(Multipole_Position_Centre_of_Mass)
            do i=1,monomer_size
              mono_list%sites(i_site,j)%pos_grad_positions(:,:,i) =  ( mono_list%masses(i) / monomer_mass) * identity3x3
            end do
          case(Multipole_Position_Atomic)
            mono_list%sites(i_site,j)%pos_grad_positions(:,:,site_atom_map(i_site)) = identity3x3
          case(Multipole_Position_M_Site)
            call m_site_position_grads(mono_list%gammaM,mono_list%sites(i_site,j)%pos_grad_positions)
          case default
            call system_abort("Multipole_Moments_Assign doesn't know how multipole position varies with atomic positions")
        end select
      end if

    end do

    if (sites_per_mono > 1 .and. intramolecular_factor /= 0) then
      do i_site=1,sites_per_mono ! no double counting here
        do j_site=i_site+1,sites_per_mono

            ! make a couple of dummy sites
            dummy_i = mono_list%sites(i_site,j)
            dummy_j = mono_list%sites(j_site,j)
            
            distance_vector = diff_min_image(at, dummy_i%position,dummy_j%position) ! within same molecule so always minimum image position
            dummy_j%position = dummy_i%position + distance_vector
            if (norm(distance_vector) > my_cutoff) cycle
            site_site_energy=0.0_dp
            if (calc_opts%do_force) then
              dummy_i%e_grad_pos = 0.0_dp
              dummy_i%e_grad_charge = 0.0_dp
              dummy_i%e_grad_dipole = 0.0_dp

              dummy_j%e_grad_pos = 0.0_dp
              dummy_j%e_grad_charge = 0.0_dp
              dummy_j%e_grad_dipole = 0.0_dp
            endif

            call Multipole_Moments_Site_Site_Interaction(site_site_energy,dummy_i,dummy_j,calc_opts, cutoff=cutoff, error=error,test=.false.)

            if (calc_opts%do_energy) e = e +intramolecular_factor * site_site_energy
            if (calc_opts%do_force) then
              mono_list%sites(i_site,j)%e_grad_pos    = mono_list%sites(i_site,j)%e_grad_pos       + intramolecular_factor * dummy_i%e_grad_pos
              mono_list%sites(i_site,j)%e_grad_charge = mono_list%sites(i_site,j)%e_grad_charge    + intramolecular_factor * dummy_i%e_grad_charge
              mono_list%sites(i_site,j)%e_grad_dipole = mono_list%sites(i_site,j)%e_grad_dipole    + intramolecular_factor * dummy_i%e_grad_dipole

              mono_list%sites(j_site,j)%e_grad_pos    = mono_list%sites(j_site,j)%e_grad_pos       + intramolecular_factor * dummy_j%e_grad_pos
              mono_list%sites(j_site,j)%e_grad_charge = mono_list%sites(j_site,j)%e_grad_charge    + intramolecular_factor * dummy_j%e_grad_charge
              mono_list%sites(j_site,j)%e_grad_dipole = mono_list%sites(j_site,j)%e_grad_dipole    + intramolecular_factor * dummy_j%e_grad_dipole
            end if
        end do
      end do
    end if
  end do


end subroutine Multipole_Moments_Assign

subroutine Multipole_Moments_Make_Dummy_Atoms(at_out,at,multipoles,error)
  type(Atoms), intent(out) :: at_out
  type(Atoms), intent(in) :: at
  type(Multipole_Moments), intent(inout) :: multipoles
  integer,optional,intent(inout) :: error
  type(Dictionary) :: properties   ! these are the properties of the dummy atoms object
  integer :: N_sites, i_type, i_mono, i_site, i_dummy
  real(dp), dimension(:), pointer :: site_charges

  ! get the stuff we need from the multipoles object : the number of sites, their positions, charges.
  call initialise(properties)
  N_sites=0
  do i_type=1,multipoles%n_monomer_types
    N_sites = N_sites + size(multipoles%monomer_types(i_type)%sites)
  end do
  allocate(multipoles%dummy_map(3,N_sites))
  allocate(multipoles%site_polarisable(N_sites))
  call initialise(at_out,N_sites,at%lattice,params=at%params,error=error)  

  call add_property(at_out, 'dummy_charge', 0.0_dp, n_cols=1, ptr=site_charges ,error=error)

  i_dummy=1
  do i_type=1,multipoles%n_monomer_types
    do i_mono=1,multipoles%monomer_types(i_type)%n_monomers
      do i_site=1,size(multipoles%monomer_types(i_type)%site_types)

        multipoles%dummy_map(:,i_dummy) = (/i_type,i_mono,i_site/)  ! so that we can easily find which dummy atom maps to which Multipole_Interactions_Site type.
        multipoles%site_polarisable(i_dummy) = multipoles%monomer_types(i_type)%sites(i_site,i_mono)%polarisable
        site_charges(i_dummy) = multipoles%monomer_types(i_type)%sites(i_site,i_mono)%charge
        at_out%pos(:,i_dummy) = multipoles%monomer_types(i_type)%sites(i_site,i_mono)%position
        at_out%Z(i_dummy)=1 ! just say dummy atoms are Hydrogen
        i_dummy = i_dummy + 1

      end do
    end do
  end do

  if (count(multipoles%site_polarisable) > 0 .and. multipoles%polarisation == Polarisation_Method_None) then
    call print("Warning : polarisable sites specified but polarisation method set to 'none'")
  end if

  call set_cutoff(at_out,at%cutoff)
  call calc_connect(at_out)

end subroutine Multipole_Moments_Make_Dummy_Atoms

subroutine Multipole_Moments_Forces_From_Dummy_Atoms(at,multipoles,force_dummy,force_atomic,error)
  type(Atoms), intent(in) :: at ! this is a dummy atoms object
  type(Multipole_Moments), intent(inout) :: multipoles
  integer,optional,intent(inout) :: error
  real(dp), dimension(:,:), intent(in) :: force_dummy
  real(dp), dimension(:,:), intent(out) :: force_atomic

  integer :: N_sites, i, j, k,m, i_site, i_atom
  real(dp), dimension(:,:), allocatable :: site_positions, forces_monomer
  real(dp), dimension(:), allocatable :: site_charges

  ! this combines any pre-existing energy gradients on the multipole sites (from multipole_moments_pre_calc)
  ! with the forces on the dummy atoms to generate the forces on the real atoms.

  N_sites = at%N
  if (size(force_dummy,2) /= N_sites) then
    RAISE_ERROR("Multipole_Moments_Forces_From_Dummy_Atoms : dummy force array of wrong size",error)
  end if
  force_atomic = 0.0_dp
  i_site=1
  do i=1,multipoles%n_monomer_types
    do j=1,multipoles%monomer_types(i)%n_monomers
      do k=1,size(multipoles%monomer_types(i)%site_types)
        if (allocated(forces_monomer)) deallocate(forces_monomer)

        multipoles%monomer_types(i)%sites(k,j)%e_grad_pos = multipoles%monomer_types(i)%sites(k,j)%e_grad_pos -force_dummy(:,i_site) ! e_grad_pos is a gradient not a force

        call Multipole_Moments_Site_to_Atoms(multipoles%monomer_types(i)%sites(k,j),forces_monomer,error)
        do m= 1,size(multipoles%monomer_types(i)%signature)
          i_atom = multipoles%monomer_types(i)%monomer_indices(m,j)
          force_atomic(:,i_atom) = force_atomic(:,i_atom) + forces_monomer(:,m)
        end do
        i_site = i_site + 1
      end do
    end do
  end do

end subroutine Multipole_Moments_Forces_From_Dummy_Atoms

!! Calculate interactions between multipole sites. at_in is a dummy atoms object where the position of each 'atom'
!! corresponds to a multipole site, so that we can use all the usual connection funcitonality to find neighbours etc.
!! the force array which is returned corresponds to forces on the actual atoms of the system, i.e. the atoms object
!! which was used to initialise multipoles
subroutine Multipole_Moments_Calc(at_in,multipoles,intermolecular_only, e, local_e, f, virial, local_virial, cutoff, error)
  type(Atoms), intent(in), target :: at_in
  type(Multipole_Moments), intent(inout) :: multipoles
  logical,intent(in) :: intermolecular_only
  real(dp), pointer :: erf_kappa, damp_exp_scale, erfc_kappa, yukawa_alpha,smoothlength
  integer, pointer :: damp_exp_order
  real(dp), intent(in), optional :: cutoff
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), dimension(:,:),optional :: f, local_virial   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial
  integer, optional, intent(out) :: error

  logical :: do_local_e,do_virial,do_local_virial
  logical :: monomers_identical, moments_converged, final_run
  integer :: i, j, k, m, n, i_mono, j_mono, i_site, j_site,n_monomer_types, i_atomic, j_atomic, i_dummy, j_dummy, n_pol, i_pos, j_pos, i_pol, method
  integer,dimension(1) :: loc
  integer,dimension(3) :: map

  type(Multipole_Interactions_Site) :: site1, site2
  real(dp), dimension(:,:),allocatable :: forces_monomer, pol_matrix
  real(dp), dimension(:),allocatable :: induced_dipoles , perm_dipoles, perm_field, alphas
  real(dp) :: site_site_energy, my_cutoff, r_ij, delta_dip, e_ind
  real(dp), dimension(3) :: diff, old_dip
  real(dp), dimension(3,3) :: prop_mat ! temporary storage of dipole field tensor
  integer, dimension(:), allocatable :: polarisable_map
  type(LA_Matrix) :: la_pol_matrix

  type(Atoms), target :: my_at
  type(Atoms), pointer :: at => null()

  my_cutoff = optional_default(at_in%cutoff,cutoff)
  if( present(cutoff) .and. (my_cutoff > at_in%cutoff) ) then
      my_at = at_in
      call set_cutoff(my_at,cutoff)
      call calc_connect(my_at)
      at => my_at
  else
      at => at_in
  endif

  if (multipoles%polarisation == Polarisation_Method_None ) then
    moments_converged = .True.
  else
    moments_converged = .False.
    N_pol = count(multipoles%site_polarisable)
    allocate(polarisable_map(N_pol))
    allocate(alphas(N_pol))
    allocate(induced_dipoles(3*N_pol))
    allocate(perm_dipoles(3*N_pol))
    allocate(perm_field(3*N_pol))

    alphas=0.0_dp
    induced_dipoles=0.0_dp
    perm_dipoles=0.0_dp
    perm_field=0.0_dp

    i_pol=1
    do i_dummy=1,at%N ! locate all polarisable sites and store any permanent dipoles on them
      if (multipoles%site_polarisable(i_dummy)) then
        polarisable_map(i_pol)=i_dummy
        i_pos= 3*i_pol - 2
        map =  multipoles%dummy_map(:,i_dummy)
        perm_dipoles(i_pos:i_pos+2) = multipoles%monomer_types(map(1))%sites(map(3),map(2))%dipole
        alphas(i_pol) = multipoles%monomer_types(map(1))%sites(map(3),map(2))%alpha
        i_pol=i_pol+1
      end if
    end do
  end if

  if (multipoles%polarisation == Polarisation_Method_GMRES .or. multipoles%polarisation == Polarisation_Method_QR ) then
    allocate(pol_matrix(3*N_pol,3*N_pol))
    pol_matrix = 0.0_dp
    do i_pol=1,N_pol
      do i_pos=3*i_pol-2,3*i_pol
       pol_matrix(i_pos,i_pos) = 1.0_dp / alphas(i_pol)
      end do
    end do
  end if

  do_local_e=present(local_e)
  do_virial=present(virial)
  do_local_virial=present(local_virial)

  if (do_local_e .or. do_virial .or. do_local_virial) then
    RAISE_ERROR("Multipole_moments_finite Sum currently only supports energy and force calculation. No virials or local energies",error)
  end if

  multipoles%calc_opts%do_pot=.false.
  multipoles%calc_opts%do_field = (multipoles%polarisation /= Polarisation_Method_None)

  final_run = .false.

  do while (.not. final_run )

    final_run = moments_converged
    multipoles%calc_opts%do_energy=present(e) .and. final_run
    multipoles%calc_opts%do_force =present(f) .and. final_run
    if (multipoles%calc_opts%do_energy) e = 0.0_dp
    if (multipoles%calc_opts%do_force) f = 0.0_dp 

    if (multipoles%calc_opts%do_field) then
      do i_dummy=1,at%N
        map =  multipoles%dummy_map(:,i_dummy)
        multipoles%monomer_types(map(1))%sites(map(3),map(2))%e_field = 0.0_dp
      end do
      perm_field = 0.0_dp
    end if

    do i_dummy=1,at%N

      i=multipoles%dummy_map(1,i_dummy)     
      i_mono=multipoles%dummy_map(2,i_dummy)
      i_site=multipoles%dummy_map(3,i_dummy)

      site1 = multipoles%monomer_types(i)%sites(i_site,i_mono)

      do n = 1, n_neighbours(at,i_dummy)

        j_dummy = neighbour(at,i_dummy,n,distance=r_ij,diff=diff) 
        if( r_ij > my_cutoff )  cycle
        j=multipoles%dummy_map(1,j_dummy)
        j_mono=multipoles%dummy_map(2,j_dummy)
        j_site=multipoles%dummy_map(3,j_dummy)
        ! make a copy of neighbour site and move it to the correct position
        site2 = multipoles%monomer_types(j)%sites(j_site,j_mono)
        site2%position = site1%position + diff

        site_site_energy=0.0_dp
        if (multipoles%calc_opts%do_force) then
          site1%e_grad_pos = 0.0_dp
          site1%e_grad_charge = 0.0_dp
          site1%e_grad_dipole = 0.0_dp

          site2%e_grad_pos = 0.0_dp
          site2%e_grad_charge = 0.0_dp
          site2%e_grad_dipole = 0.0_dp
        endif

        call Multipole_Moments_Site_Site_Interaction(site_site_energy,site1,site2, multipoles%calc_opts, cutoff=cutoff, error=error,test=.false.)

        if (multipoles%calc_opts%do_energy) e = e + 0.5_dp * site_site_energy       ! add half this contribution since double counting
        if (multipoles%calc_opts%do_force) then                                    ! add half of dummy minimum image site gradients over to actual sites
          multipoles%monomer_types(i)%sites(i_site,i_mono)%e_grad_pos = multipoles%monomer_types(i)%sites(i_site,i_mono)%e_grad_pos       + 0.5_dp *site1%e_grad_pos
          multipoles%monomer_types(i)%sites(i_site,i_mono)%e_grad_charge = multipoles%monomer_types(i)%sites(i_site,i_mono)%e_grad_charge + 0.5_dp *site1%e_grad_charge
          multipoles%monomer_types(i)%sites(i_site,i_mono)%e_grad_dipole = multipoles%monomer_types(i)%sites(i_site,i_mono)%e_grad_dipole + 0.5_dp *site1%e_grad_dipole

          multipoles%monomer_types(j)%sites(j_site,j_mono)%e_grad_pos = multipoles%monomer_types(j)%sites(j_site,j_mono)%e_grad_pos       + 0.5_dp *site2%e_grad_pos
          multipoles%monomer_types(j)%sites(j_site,j_mono)%e_grad_charge = multipoles%monomer_types(j)%sites(j_site,j_mono)%e_grad_charge + 0.5_dp *site2%e_grad_charge
          multipoles%monomer_types(j)%sites(j_site,j_mono)%e_grad_dipole = multipoles%monomer_types(j)%sites(j_site,j_mono)%e_grad_dipole + 0.5_dp *site2%e_grad_dipole
        end if
        if (multipoles%calc_opts%do_field) then
          multipoles%monomer_types(i)%sites(i_site,i_mono)%e_field = multipoles%monomer_types(i)%sites(i_site,i_mono)%e_field       + 0.5_dp *site1%e_field
          multipoles%monomer_types(j)%sites(j_site,j_mono)%e_field = multipoles%monomer_types(j)%sites(j_site,j_mono)%e_field       + 0.5_dp *site2%e_field
        end if

        if ((r_ij) > multipoles%polarisation_cutoff ) cycle 
        if (multipoles%polarisation == Polarisation_Method_GMRES .or. multipoles%polarisation == Polarisation_Method_QR ) then          
          i_pos=0
          j_pos=0
          if (count(polarisable_map == i_dummy)==1) then
            loc = find(polarisable_map == i_dummy)           ! get matrix indices
            i_pos= 3*loc(1) - 2          
          end if
          if (count(polarisable_map == j_dummy)==1) then
          loc = find(polarisable_map == j_dummy)
          j_pos= 3*loc(1) - 2
          end if

          if (i_pos /= 0) perm_field(i_pos:i_pos+2) = perm_field(i_pos:i_pos+2) + 0.5_dp *site1%e_field
          if (j_pos /= 0) perm_field(j_pos:j_pos+2) = perm_field(j_pos:j_pos+2) + 0.5_dp *site2%e_field
          if (i_pos /= 0 .and. j_pos /=0 ) then                                
            prop_mat = 0.5_dp * T_rank_two(diff,multipoles%calc_opts,site1%damp_rad,site2%damp_rad,cutoff=cutoff) 
            pol_matrix(i_pos:i_pos+2,j_pos:j_pos+2) = pol_matrix(i_pos:i_pos+2,j_pos:j_pos+2) - prop_mat
            pol_matrix(j_pos:j_pos+2,i_pos:i_pos+2) = pol_matrix(j_pos:j_pos+2,i_pos:i_pos+2) - prop_mat
          end if
        end if

      end do            
    end do

    ! Here's where we get the induced dipole moments
    if (.not. final_run) then
      if (multipoles%polarisation == Polarisation_Method_FPI) then
        delta_dip=0.0_dp
        do i_pol=1,N_pol
          map = multipoles%dummy_map(:,polarisable_map(i_pol))
          i_pos= 3*i_pol -2
          old_dip = induced_dipoles(i_pos:i_pos+2)

          induced_dipoles(i_pos:i_pos+2) = multipoles%monomer_types(map(1))%sites(map(3),map(2))%e_field * alphas(i_pol)
          multipoles%monomer_types(map(1))%sites(map(3),map(2))%dipole = perm_dipoles(i_pos:i_pos+2) + induced_dipoles(i_pos:i_pos+2)
          do k=1,3
            delta_dip = delta_dip + (old_dip(k) - induced_dipoles(i_pos+k-1))**2
          end do
        end do  
        delta_dip = delta_dip / 3*N_pol

        if (sqrt(delta_dip) < multipoles%dipole_tolerance) then
          moments_converged = .True.     
        end if
      else if (multipoles%polarisation == Polarisation_Method_GMRES .or. multipoles%polarisation == Polarisation_Method_QR) then
        call initialise(la_pol_matrix,pol_matrix)

        ! solve system with QR or GMRES
        if ( multipoles%polarisation == Polarisation_Method_QR) then
          call LA_Matrix_QR_Factorise(la_pol_matrix,error=error)
          call LA_Matrix_QR_Solve_Vector(la_pol_matrix,perm_field,induced_dipoles,error=error)
        end if

        do i_pol=1,N_pol
          map = multipoles%dummy_map(:,polarisable_map(i_pol))  
          i_pos= 3*i_pol -2
          multipoles%monomer_types(map(1))%sites(map(3),map(2))%dipole = multipoles%monomer_types(map(1))%sites(map(3),map(2))%dipole + induced_dipoles(i_pos:i_pos+2)
        end do
        moments_converged = .True.
      end if

      if (moments_converged) then  ! Calculate induction energy
        do i_pol=1,N_pol
          !map = multipoles%dummy_map(:,polarisable_map(i_pol))  
          ! call print("E FIELD "//multipoles%monomer_types(map(1))%sites(map(3),map(2))%e_field)
          i_pos= 3*i_pol -2
          e_ind = e_ind + 0.5_dp * normsq(induced_dipoles(i_pos:i_pos+2)) /alphas(i_pol) 
          !e_ind = e_ind - normsq(multipoles%monomer_types(map(1))%sites(map(3),map(2))%e_field) * 0.5_dp*alphas(i_pol) 
        end do
      end if
    end if
  end do
!  call print("induction energy : "//e_ind)
  if (multipoles%calc_opts%do_energy) e = e + e_ind 

  if (multipoles%calc_opts%do_force) then  ! translate site gradients into atomic forces
    do i=1,size(multipoles%monomer_types)   
      do i_mono=1,multipoles%monomer_types(i)%n_monomers
        do i_site=1,size(multipoles%monomer_types(i)%site_types)      

          if (allocated(forces_monomer)) deallocate(forces_monomer)
          call Multipole_Moments_Site_to_Atoms(multipoles%monomer_types(i)%sites(i_site,i_mono),forces_monomer,error)

          do i_atomic=1,size(multipoles%monomer_types(i)%signature)
            f(:,multipoles%monomer_types(i)%monomer_indices(i_atomic,i_mono)) =  f(:,multipoles%monomer_types(i)%monomer_indices(i_atomic,i_mono)) + forces_monomer(:,i_atomic) 
          end do

        end do 
      end do
    end do
  end if
  deallocate(multipoles%dummy_map)
  deallocate(multipoles%site_polarisable)
end subroutine Multipole_Moments_Calc


subroutine Multipole_Moments_Site_to_Atoms(this,forces, error)

  type(Multipole_Interactions_Site), intent(in) :: this
  real(dp), dimension(:,:),allocatable:: forces

  integer :: i, j, k, a, n_atoms
  integer, optional, intent(out) :: error

  n_atoms = size(this%charge_grad_positions,3)
  allocate(forces(3,n_atoms))
  forces=0.0_dp

  ! add up force contributions in most transparent way possible
  do a=1,n_atoms ! a is atom number, i component of atomic position, j component of dipole, k component of site position
    do i=1,3
      forces(i,a) = forces(i,a) - this%e_grad_charge * this%charge_grad_positions(1,i,a) ! force is minus grad

      do j=1,3
        forces(i,a) = forces(i,a) - this%e_grad_dipole(j) * this%dipole_grad_positions(j,i,a) ! force is minus grad
      end do
      do k=1,3
        forces(i,a) = forces(i,a) - this%e_grad_pos(k) * this%pos_grad_positions(k,i,a)
      end do
    end do
  end do

end subroutine Multipole_Moments_Site_to_Atoms


end module Multipoles_module
