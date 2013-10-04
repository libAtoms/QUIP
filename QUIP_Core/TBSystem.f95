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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X TBSystem module
!X
!% Code to contain a TB calculation, and evaluate H and S matrices etc.
!% Also contains much code to support self-consitency, dipole matrix elements
!% and spin-orbit coupling
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module TBSystem_module

use error_module
use system_module, only : dp, print, inoutput, PRINT_NORMAL, PRINT_ALWAYS, PRINT_NERD, current_verbosity, optional_default, operator(//), print, verbosity_push_decrement, verbosity_pop
use units_module, only : Hartree, Bohr, PI
use periodictable_module, only : ElementName
use linearalgebra_module, only : operator(.mult.), operator(.feq.), norm
use mpi_context_module, only : mpi_context
use dictionary_module, only : dictionary
use paramreader_module, only : param_register, param_read_line
use atoms_module, only : atoms, n_neighbours, neighbour, assignment(=)

! use Functions_module 
use QUIP_Common_module ! , only : xml_t, dictionary_t, pauli_sigma
use ScaLAPACK_module, only : scalapack, initialise
use TB_Common_module
use TBModel_module, only : tbmodel, initialise, finalise, print, &
   n_orb_sets_of_Z, orb_type_of_orb_set_of_Z, n_orbs_of_orb_set_of_Z, n_orbs_of_Z, n_elecs_of_Z, &
   get_HS_blocks, get_dHS_blocks, get_dHS_masks
use Matrix_module, only : matrixd
use TBMatrix_module, only : tbmatrix, initialise, finalise, wipe, print, zero, add_block, sum_matrices
use TB_KPoints_module, only : kpoints, initialise, finalise, calc_phase, print, init_mpi, ksum_distrib_inplace
use TB_mixing_module, only : do_mix_simple, do_ridders_residual, do_mix_broyden
use ewald_module, only : add_madelung_matrix, add_dmadelung_matrix, add_dmadelung_matrix_dr


implicit none
private

! character(40) :: sc_default_file = "self_consistency.data.xml"

integer, parameter, public :: &
  SCF_NONE = 0, &
  SCF_LOCAL_U = 1, &
  SCF_GLOBAL_U = 2, &
  SCF_LCN = 3, &
  SCF_GCN = 4, &
  SCF_NONLOCAL_U_DFTB = 5, &
  SCF_NONLOCAL_U_NRL_TB = 6, &
  SCF_SPIN_DIR = 7, &
  SCF_SPIN_STONER = 8

character(len=30), parameter :: scf_names(0:8) = (/ &
'NONE                          ', &
'LOCAL_U                       ', &
'GLOBAL_U                      ', &
'LCN                           ', &
'GCN                           ', &
'NONLOCAL_U_DFTB               ', &
'NONLOCAL_U_NRL_TB             ', &
'SPIN_DIR                      ', &
'SPIN_STONER                   ' /)

! LOCAL_U 		F = sum_a U_a (n_a-n0_a)^2
! GLOBAL_U 		F = U (sum_a n_a-n0_a)^2
! LCN
! GCN
! NONLOCAL_U_DFTB 	F = sum_ab (n_a-n0_a) gamma_ab (n_b-n0_b)
! NONLOCAL_U_NRL_TB 	F = sum_ab (n_a-n0_a) gamma_ab (n_b-n0_b)
! SPIN_DIR		F = -0.5 sum_a splitting_a | atomic_mom_a(i) |
! SPIN_STONER		F = -0.25 sum_a stoner_param_a | atomic_mom_a(i) | ^2

public :: Self_Consistency_Term
type Self_Consistency_Term
  logical :: active = .false.
  integer :: type = SCF_NONE
  integer :: N = 0, N_atoms = 0, N_manifolds = 0
  integer :: n_dof = 0

  ! for SCF_GLOBAL_U
  real(dp) :: global_U
  real(dp) :: global_N

  ! for SCF_GLOBAL_U
  real(dp) :: global_pot

  ! for SCF_LOCAL_U and SCF_NONLOCAL_U
  real(dp), allocatable :: U(:)
  real(dp), allocatable :: atomic_n(:), atomic_n0(:)
  real(dp), allocatable :: datomic_local_pot_dr(:,:)
  ! for SCF_LCN
  real(dp), allocatable :: atomic_local_pot(:)

  ! for SCF_NONLOCAL_U_*
  real(dp), allocatable :: gamma(:,:), dgamma_dr(:,:,:)

  ! for SCF_SPIN_*
  real(dp), allocatable :: manifold_mom(:,:)
  ! for SCF_SPIN_DIR
  real(dp), allocatable :: spin_splitting(:)
  ! for SCF_SPIN_STONER
  real(dp), allocatable :: stoner_param(:)

end type Self_Consistency_Term

public :: Self_Consistency
type Self_Consistency
  logical :: active = .false.

  integer :: N = 0, N_atoms = 0, N_manifolds = 0

  real(dp) :: global_U 
  real(dp), allocatable :: U(:), stoner_param(:,:)

  type(Self_Consistency_Term), allocatable :: terms(:)
  real(dp), allocatable :: orb_local_pot(:), orb_exch_field(:,:)

  real(dp) :: alpha = 0.01_dp, w0 = 0.02_dp
  real(dp) :: conv_tol = 1.0e-10_dp
  real(dp) :: mix_simple_end_tol = -1.0e-2_dp, mix_simple_start_tol=1.0e38_dp

  integer :: max_iter = 100

  logical :: mix_simple = .false.

end type Self_Consistency

public :: TB_Dipole_Model
type TB_Dipole_Model
  logical :: active = .false.

  integer :: n_types
  integer, allocatable :: type_of_atomic_num(:), atomic_num(:)
  integer, allocatable :: n_orb_sets(:), orb_set_type(:,:), orb_set_phase(:,:)
  real(dp), allocatable :: gaussian_width(:,:)

end type TB_Dipole_Model

public :: TB_Spin_Orbit_Coupling
type TB_Spin_Orbit_Coupling
  logical :: active = .false.

  integer :: n_types
  integer, allocatable :: type_of_atomic_num(:), atomic_num(:)
  integer, allocatable :: n_orb_sets(:), orb_set_type(:,:)
  real(dp), allocatable :: SO_param(:,:)
end type TB_Spin_Orbit_Coupling

public :: TBSystem
type TBSystem

   ! model and atomic config
  integer :: N = 0, N_atoms = 0, N_manifolds = 0
  integer, allocatable :: at_Z(:)
  type(TBModel) tbmodel
  ! first_orb_of_atom: index into array of all orbitals of first orbital associated with each atom
  ! first_manifold_of_atom: index into array of all manifolds of first manifold (s, p, d, etc) associated with each atom
  ! first_orb_of_manifold: index into array of all orbitals of first orbital associated with each manifold
  integer, allocatable :: first_orb_of_atom(:), first_manifold_of_atom(:), first_orb_of_manifold(:)

  logical :: noncollinear = .false.
  logical :: complex_matrices = .false.

  integer :: max_block_size = 0

  ! matrices
  integer :: n_matrices = 0
  type(TBMatrix) :: H
  type(TBMatrix) :: S

  type(TBMatrix) :: dH(3)
  type(TBMatrix) :: dS(3)

  ! k points
  logical :: kpoints_generate_dynamically = .false., kpoints_generate_next_dynamically = .false., kpoints_use_mp = .true.
  real(dp) :: kpoints_k_space_density = 1.0_dp
  type(KPoints) :: kpoints

  ! self consistency
  type(Self_Consistency) :: scf
  type(TB_Dipole_model) :: dipole_model
  type(TB_Spin_Orbit_Coupling) :: SO

  type(MPI_context) :: mpi_global, mpi_my_matrices
  type(ScaLAPACK) :: scalapack_my_matrices

end type TBSystem

logical :: parse_in_self_consistency
type(Self_Consistency), pointer :: parse_self_consistency

logical :: parse_in_dipole_model
type(TB_Dipole_Model), pointer :: parse_dipole_model

logical :: parse_in_spin_orbit_coupling
type(TB_Spin_Orbit_Coupling), pointer :: parse_spin_orbit_coupling

public :: Initialise
interface Initialise
  module procedure TBSystem_Initialise_str, TBSystem_Initialise_from_tbsys
  module procedure Self_Consistency_Initialise_str, TB_Dipole_Model_Initialise_str
  module procedure TB_Spin_Orbit_Coupling_Initialise_str
end interface Initialise

public :: Setup_atoms
interface Setup_atoms
  module procedure TBSystem_Setup_atoms_from_atoms, TBSystem_Setup_atoms_from_tbsys, TBSystem_setup_atoms_from_arrays
end interface Setup_atoms

public :: Setup_deriv_matrices
interface Setup_deriv_matrices
  module procedure TBSystem_Setup_deriv_matrices
end interface Setup_deriv_matrices

interface Setup_system
  module procedure Self_Consistency_setup_system
  module procedure Self_Consistency_Term_setup_system
end interface Setup_system

public :: Finalise
interface Finalise
  module procedure TBSystem_Finalise, Self_Consistency_Finalise, Self_Consistency_Term_Finalise
  module procedure TB_Dipole_Model_Finalise, TB_Spin_Orbit_Coupling_Finalise
end interface Finalise

public :: Print
interface Print
  module procedure TBSystem_Print, Self_Consistency_Print, Self_Consistency_Term_Print
  module procedure TB_Dipole_Model_Print, TB_Spin_Orbit_Coupling_Print
end interface Print

public :: Wipe
interface Wipe
  module procedure TBSystem_Wipe, Self_Consistency_Wipe, Self_Consistency_Term_Wipe
end interface Wipe

interface read_params_xml
  module procedure Self_Consistency_read_params_xml, TB_Dipole_Model_read_params_xml, TB_Spin_Orbit_Coupling_read_params_xml
end interface read_params_xml

public :: fill_matrices
interface fill_matrices
  module procedure TBSystem_fill_matrices
end interface fill_matrices

public :: fill_these_matrices
interface fill_these_matrices
  module procedure TBSystem_fill_these_matrices
end interface fill_these_matrices

public :: fill_dmatrices
interface fill_dmatrices
  module procedure TBSystem_fill_dmatrices
end interface fill_dmatrices

public :: fill_sc_matrices
interface fill_sc_matrices
  module procedure TBSystem_fill_sc_matrices, Self_Consistency_Term_fill_sc_matrices
end interface fill_sc_matrices

public :: fill_sc_dmatrices
interface fill_sc_dmatrices
  module procedure TBSystem_fill_sc_dmatrices, Self_Consistency_Term_fill_sc_dmatrices
end interface fill_sc_dmatrices

public :: atom_orbital_sum
interface atom_orbital_sum
  module procedure TBSystem_atom_orbital_sum_real1, TBSystem_atom_orbital_sum_real2
  module procedure TBSystem_atom_orbital_sum_complex2
end interface atom_orbital_sum

public :: manifold_orbital_sum
interface manifold_orbital_sum
  module procedure TBSystem_manifold_orbital_sum_real2
end interface manifold_orbital_sum

public :: atom_orbital_spread
interface atom_orbital_spread
  module procedure TBSystem_atom_orbital_spread_real1
end interface atom_orbital_spread

public :: atom_orbital_spread_mat
interface atom_orbital_spread_mat
  module procedure TBSystem_atom_orbital_spread_mat_r
end interface atom_orbital_spread_mat

public :: ksum_atom_orbital_sum_mat
interface ksum_atom_orbital_sum_mat
  module procedure TBSystem_ksum_atom_orbital_sum_mat_d
end interface ksum_atom_orbital_sum_mat

public :: atom_orbital_sum_mat
interface atom_orbital_sum_mat
  module procedure TBSystem_atom_orbital_sum_mat_r
end interface atom_orbital_sum_mat

public :: set_type
interface set_type
  module procedure Self_Consistency_set_type_str
end interface set_type

public :: calc_orb_local_pot
interface calc_orb_local_pot
  module procedure TBSystem_calc_orb_local_pot
end interface calc_orb_local_pot

public :: update_orb_local_pot
interface update_orb_local_pot
  module procedure TBSystem_update_orb_local_pot
end interface update_orb_local_pot

public :: scf_e_correction, local_scf_e_correction, scf_f_correction, scf_virial_correction

public :: n_elec
interface n_elec
  module procedure TBSystem_n_elec
end interface n_elec

public :: scf_get_atomic_n_mom
interface scf_get_atomic_n_mom
  module procedure TBSystem_scf_get_atomic_n_mom
end interface scf_get_atomic_n_mom

public :: scf_get_global_N
interface scf_get_global_N
  module procedure TBSystem_scf_get_global_N
end interface scf_get_global_N

public :: scf_set_atomic_n_mom
interface scf_set_atomic_n_mom
  module procedure TBSystem_scf_set_atomic_n_mom
end interface scf_set_atomic_n_mom

public :: scf_set_global_N
interface scf_set_global_N
  module procedure TBSystem_scf_set_global_N
end interface scf_set_global_N

public :: add_term_d2SCFE_dgNdn, add_term_dscf_e_correction_dgN, add_term_d2SCFE_dn2_times_vec, add_term_dscf_e_correction_dn

public :: initialise_kpoints
interface initialise_kpoints
  module procedure TBSystem_initialise_kpoints
end interface initialise_kpoints

contains
subroutine TBSystem_Initialise_str(this, args_str, param_str, kpoints_obj, mpi_obj)
  type(TBSystem), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(KPoints), intent(in), optional :: kpoints_obj
  type(MPI_context), intent(in), optional :: mpi_obj

  call Finalise(this)

  call Initialise(this%tbmodel, args_str, param_str)

  if (present(kpoints_obj)) then
    call Initialise(this%kpoints, kpoints_obj, mpi_obj)
    call Initialise_tbsystem_k_dep_stuff(this, mpi_obj)
  else
    call initialise_kpoints(this, args_str, param_str, mpi_obj, from_tbsystem_initialise=.true.)
  endif

  call Initialise(this%scf, args_str, param_str)
  call Initialise(this%dipole_model, args_str, param_str)
  call Initialise(this%SO, args_str, param_str)

end subroutine TBSystem_Initialise_str

subroutine check_dipole_model_consistency(dipole_model, tbm, Z)
  type(TB_Dipole_Model), intent(in) :: dipole_model
  type(TBModel), intent(in) :: tbm
  integer, intent(in) :: Z(:)

  integer :: i_at
  integer :: i_type, i_orb_set

  do i_at=1, size(Z)
    i_type = dipole_model%type_of_atomic_num(Z(i_at))
    if (dipole_model%n_orb_sets(i_type) /= n_orb_sets_of_Z(tbm, Z(i_at))) &
      call system_abort("check_dipole_model_consistency found mismatch in number of orbitals sets for Z="//Z(i_at) // &
	" dipole_model " // dipole_model%n_orb_sets(i_type) // " tbmodel " // n_orb_sets_of_Z(tbm, Z(i_at)) )
    do i_orb_set=1, dipole_model%n_orb_sets(i_type)
      if (dipole_model%orb_set_type(i_orb_set,i_type) /= orb_type_of_orb_set_of_Z(tbm, Z(i_at), i_orb_set)) &
	call system_abort("check_dipole_model_consistency found mismatch in orbital set type  for Z="//Z(i_at) // &
	  " i_orb_set="//i_orb_set// " dipole_model " // dipole_model%orb_set_type(i_orb_set,i_type) // " tbmodel " // orb_type_of_orb_set_of_Z(tbm, Z(i_at), i_orb_set) )
    end do
  end do
end subroutine check_dipole_model_consistency

subroutine check_spin_orbit_coupling_consistency(spin_orbit_coupling, tbm, Z)
  type(TB_Spin_Orbit_Coupling), intent(in) :: spin_orbit_coupling
  type(TBModel), intent(in) :: tbm
  integer, intent(in) :: Z(:)

  integer :: i_at
  integer :: i_type, i_orb_set

  do i_at=1, size(Z)
    i_type = spin_orbit_coupling%type_of_atomic_num(Z(i_at))
    if (spin_orbit_coupling%n_orb_sets(i_type) /= n_orb_sets_of_Z(tbm, Z(i_at))) &
      call system_abort("check_spin_orbit_coupling_consistency found mismatch in number of orbitals sets for Z="//Z(i_at) // &
	" spin_orbit_coupling " // spin_orbit_coupling%n_orb_sets(i_type) // " tbmodel " // n_orb_sets_of_Z(tbm, Z(i_at)) )
    do i_orb_set=1, spin_orbit_coupling%n_orb_sets(i_type)
      if (spin_orbit_coupling%orb_set_type(i_orb_set,i_type) /= orb_type_of_orb_set_of_Z(tbm, Z(i_at), i_orb_set)) &
	call system_abort("check_spin_orbit_coupling_consistency found mismatch in orbital set type  for Z="//Z(i_at) // &
	  " i_orb_set="//i_orb_set// " spin_orbit_coupling " // spin_orbit_coupling%orb_set_type(i_orb_set,i_type) // " tbmodel " // &
	  orb_type_of_orb_set_of_Z(tbm, Z(i_at), i_orb_set) )
    end do
  end do
end subroutine check_spin_orbit_coupling_consistency

subroutine TBSystem_Initialise_from_tbsys(this, from, mpi_obj)
  type(TBSystem), intent(inout) :: this
  type(TBSystem), intent(in) :: from
  type(MPI_context), intent(in), optional :: mpi_obj

  call Finalise(this)

  this%tbmodel = from%tbmodel
  this%complex_matrices = from%complex_matrices
  this%noncollinear = from%noncollinear

  this%kpoints = from%kpoints
  this%kpoints_generate_dynamically = from%kpoints_generate_dynamically
  this%kpoints_generate_next_dynamically = from%kpoints_generate_next_dynamically
  this%kpoints_use_mp = from%kpoints_use_mp
  this%kpoints_k_space_density = from%kpoints_k_space_density
  call init_mpi(this%kpoints, mpi_obj)

  call initialise_tbsystem_k_dep_stuff(this, mpi_obj)

  this%scf = from%scf
  this%dipole_model = from%dipole_model

  call setup_atoms(this, from)
end subroutine TBSystem_Initialise_from_tbsys

subroutine initialise_tbsystem_k_dep_stuff(this, mpi_obj)
  type(TBSystem), intent(inout) :: this
  type(MPI_Context), optional, intent(in) :: mpi_obj

  if (this%kpoints%N > 0) then
    this%n_matrices = this%kpoints%N
  else
    this%n_matrices = 0
  endif

  if (present(mpi_obj)) then
    this%mpi_global = this%kpoints%mpi_global
    call Initialise(this%scalapack_my_matrices, this%kpoints%mpi_my_kpt)
    if (.not. this%scalapack_my_matrices%active .and. this%kpoints%mpi_my_kpt%n_procs > 1) then
      this%kpoints%no_sum_over_my_kpt = .true.
    endif
    this%mpi_my_matrices = this%kpoints%mpi_my_kpt
  endif

end subroutine Initialise_tbsystem_k_dep_stuff

subroutine TBSystem_Finalise(this)
  type(TBSystem), intent(inout) :: this

  call Wipe(this)

  call Finalise(this%tbmodel)
  call Finalise(this%kpoints)
  call Finalise(this%scf)
  call Finalise(this%dipole_model)

  call Finalise(this%H)
  call Finalise(this%dH(1))
  call Finalise(this%dH(2))
  call Finalise(this%dH(3))
  call Finalise(this%S)
  call Finalise(this%dS(1))
  call Finalise(this%dS(2))
  call Finalise(this%dS(3))

  call Finalise(this%mpi_global)
end subroutine TBSystem_Finalise

subroutine TBSystem_Wipe(this)
  type(TBSystem), intent(inout) :: this

  if (allocated(this%first_orb_of_atom)) deallocate(this%first_orb_of_atom)
  if (allocated(this%first_orb_of_manifold)) deallocate(this%first_orb_of_manifold)
  if (allocated(this%first_manifold_of_atom)) deallocate(this%first_manifold_of_atom)
  if (allocated(this%at_Z)) deallocate(this%at_Z)

  call Wipe(this%scf)
  call Wipe(this%H)
  call Wipe(this%S)
  call Wipe(this%dH(1))
  call Wipe(this%dH(2))
  call Wipe(this%dH(3))
  call Wipe(this%dS(1))
  call Wipe(this%dS(2))
  call Wipe(this%dS(3))

  call Finalise(this%mpi_my_matrices)

  this%N = 0
  this%N_atoms = 0
  this%N_manifolds = 0
  this%max_block_size = 0
end subroutine TBSystem_Wipe

function TBSystem_n_elec(this, at, w_n)
  type(TBSystem), intent(in) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(in), pointer, optional :: w_n(:)
  real(dp) :: TBSystem_n_elec

  integer i

  TBSystem_n_elec = 0
  do i=1, this%N_atoms
    if (present(w_n)) then
      if (associated(w_n)) then
	TBSystem_n_elec = TBSystem_n_elec + w_n(i)*n_elecs_of_Z(this%tbmodel, at%Z(i))
      else
	TBSystem_n_elec = TBSystem_n_elec + n_elecs_of_Z(this%tbmodel, at%Z(i))
      end if
    else
      TBSystem_n_elec = TBSystem_n_elec + n_elecs_of_Z(this%tbmodel, at%Z(i))
    endif
  end do
end function TBSystem_n_elec

subroutine TBSystem_Setup_atoms_from_atoms(this, at, noncollinear, args_str, mpi_obj, error)
  type(TBSystem), intent(inout) :: this
  type(Atoms), intent(in) :: at
  logical, intent(in), optional :: noncollinear
  character(len=*), intent(in), optional :: args_str
  type(MPI_context), intent(in), optional :: mpi_obj
  integer, intent(out), optional :: error

  INIT_ERROR(error)

  call initialise_kpoints(this, args_str=args_str, mpi_obj=mpi_obj)
  if (this%kpoints_generate_dynamically) then
    call initialise(this%kpoints, at%lattice, this%kpoints_k_space_density, this%kpoints_use_mp, mpi_obj)
    call Initialise_tbsystem_k_dep_stuff(this, mpi_obj)
    this%kpoints_generate_dynamically = this%kpoints_generate_next_dynamically
  endif
  call setup_atoms(this, at%N, at%Z, noncollinear, error=error)
  PASS_ERROR(error)

end subroutine TBSystem_Setup_atoms_from_atoms

subroutine TBSystem_Setup_atoms_from_tbsys(this, from, error)
  type(TBSystem), intent(inout) :: this
  type(TBSystem), intent(in) :: from
  integer, intent(out), optional :: error

  INIT_ERROR(error)

  call setup_atoms(this, from%N_atoms, from%at_Z, from%noncollinear, error=error)
  PASS_ERROR(error)

end subroutine TBSystem_Setup_Atoms_from_tbsys

subroutine TBSystem_Setup_atoms_from_arrays(this, at_N, at_Z, noncollinear, error)
  type(TBSystem), intent(inout) :: this
  integer, intent(in) :: at_N, at_Z(:)
  logical, intent(in), optional :: noncollinear
  integer, intent(out), optional :: error

  integer :: i_at, i_man, man_offset, last_man_offset
  integer :: n_mag

  INIT_ERROR(error)

  call Wipe(this)

  this%noncollinear = optional_default(this%noncollinear, noncollinear)

  this%N_atoms = at_N
  allocate(this%at_Z(this%N_atoms))
  this%at_Z = at_Z

  if (this%noncollinear)  then
    n_mag = 2
  else
    n_mag = 1
  end if

  ! fill in first_orb_of_atom and first_manifold_of_atom arrays
  allocate(this%first_orb_of_atom(this%N_atoms+1))
  allocate(this%first_manifold_of_atom(this%N_atoms+1))
  this%first_orb_of_atom(1) = 1
  this%first_manifold_of_atom(1) = 1
  do i_at=2, this%N_atoms+1
    this%first_orb_of_atom(i_at) = this%first_orb_of_atom(i_at-1) + &
				n_mag*n_orbs_of_Z(this%tbmodel, at_Z(i_at-1))
    this%first_manifold_of_atom(i_at) = this%first_manifold_of_atom(i_at-1) + &
				        n_orb_sets_of_Z(this%tbmodel, at_Z(i_at-1))
  end do

  this%N = this%first_orb_of_atom(this%N_atoms+1)-1
  this%N_manifolds = this%first_manifold_of_atom(this%N_atoms+1)-1

  ! fill in first_orb_of_manifold array
  allocate(this%first_orb_of_manifold(this%N_manifolds+1))
  this%first_orb_of_manifold(1) = 1
  ! don't really need this, but some compilers complain because it's used before set (even though
  !   thing that uses it (last_man_offset) isn't actually used until it has a meaningul value)
  man_offset = 1
  ! set in case first atom has only 1 manifold
  last_man_offset = 1
  do i_at=1, this%N_atoms
    do i_man=this%first_manifold_of_atom(i_at), this%first_manifold_of_atom(i_at+1)-1
      ! if this is very first manifold, it was already set before entering loop, so skip
      if (i_man == 1) cycle

      ! save man_offset in case we need it for next atom's first manifold
      last_man_offset = man_offset
      ! man_offset says which manifold (1...) of this atom the loop is currently at
      man_offset = i_man - this%first_manifold_of_atom(i_at) + 1

      if (man_offset == 1) then
	! first manifold of this atom
	! first orb is first orb of first manifold of previous atom + number of orbs in that manifold
	this%first_orb_of_manifold(i_man) = this%first_orb_of_manifold(i_man-1) + &
					    n_mag*n_orbs_of_orb_set_of_Z(this%tbmodel, at_Z(i_at-1), last_man_offset)
      else
	! second or later manifold of this atom
	! first orb is first orb of previous manifold + number of orbs in that manifold
	this%first_orb_of_manifold(i_man) = this%first_orb_of_manifold(i_man-1) + &
					    n_mag*n_orbs_of_orb_set_of_Z(this%tbmodel, at_Z(i_at), man_offset-1)
      endif
      if (i_at == this%N_atoms+1) exit
    end do
  end do

  ! set last orb of last manifold 
  this%first_orb_of_manifold(this%N_manifolds+1) = this%first_orb_of_manifold(this%N_manifolds) + &
    n_mag*n_orbs_of_orb_set_of_Z(this%tbmodel, at_Z(this%N_atoms), n_orb_sets_of_Z(this%tbmodel, at_Z(this%N_atoms)))

  this%max_block_size = maxval(this%first_orb_of_atom(2:this%N_atoms+1) - &
			       this%first_orb_of_atom(1:this%N_atoms))

  this%complex_matrices = this%kpoints%non_gamma .or. this%noncollinear
  call Initialise(this%H, this%N, this%n_matrices, this%complex_matrices, this%scalapack_my_matrices)
  if (.not. this%tbmodel%is_orthogonal) then
    call Initialise(this%S, this%N, this%n_matrices, this%complex_matrices, this%scalapack_my_matrices)
  endif

  call setup_system(this%scf, this%N, this%tbmodel, this%N_atoms, this%at_Z, this%N_manifolds)

end subroutine TBSystem_Setup_Atoms_from_arrays

subroutine TBSystem_Setup_deriv_matrices(this, at, dense, need_S)
  type(TBSystem), intent(inout) :: this
  type(Atoms), intent(in) :: at
  logical, intent(in), optional :: dense, need_S

  logical do_S, my_dense

  do_S = .not. this%tbmodel%is_orthogonal
  if (present(need_S)) do_S = do_S .or. need_S

  if (this%N == 0 .or. this%N_atoms == 0) then
    call system_abort ('Called TBSystem_Setup_deriv_matrices on uninitialized TBSystem')
  endif

  my_dense = optional_default(.false., dense)

  call Finalise(this%dH(1))
  call Finalise(this%dH(2))
  call Finalise(this%dH(3))
  call Finalise(this%dS(1))
  call Finalise(this%dS(2))
  call Finalise(this%dS(3))

  if (my_dense) then
    call Initialise(this%dH(1), this%N, this%n_matrices, this%complex_matrices, this%scalapack_my_matrices)
    call Initialise(this%dH(2), this%N, this%n_matrices, this%complex_matrices, this%scalapack_my_matrices)
    call Initialise(this%dH(3), this%N, this%n_matrices, this%complex_matrices, this%scalapack_my_matrices)
    if (do_S) then
      call Initialise(this%dS(1), this%N, this%n_matrices, this%complex_matrices, this%scalapack_my_matrices)
      call Initialise(this%dS(2), this%N, this%n_matrices, this%complex_matrices, this%scalapack_my_matrices)
      call Initialise(this%dS(3), this%N, this%n_matrices, this%complex_matrices, this%scalapack_my_matrices)
    endif
  else
    call Initialise(this%dH(1), at, this%first_orb_of_atom, this%n_matrices, this%complex_matrices)
    call Initialise(this%dH(2), at, this%first_orb_of_atom, this%n_matrices, this%complex_matrices)
    call Initialise(this%dH(3), at, this%first_orb_of_atom, this%n_matrices, this%complex_matrices)
    if (do_S) then
      call Initialise(this%dS(1), at, this%first_orb_of_atom, this%n_matrices, this%complex_matrices)
      call Initialise(this%dS(2), at, this%first_orb_of_atom, this%n_matrices, this%complex_matrices)
      call Initialise(this%dS(3), at, this%first_orb_of_atom, this%n_matrices, this%complex_matrices)
    endif
  endif

end subroutine

subroutine TBSystem_fill_matrices(this, at, need_H, need_S, no_S_spin)
  type (TBSystem), intent(inout) :: this
  type(Atoms), intent(in) :: at
  logical, intent(in), optional :: need_H, need_S, no_S_spin

  logical :: do_need_H, do_need_S, do_no_S_spin

  do_need_S = .not. this%tbmodel%is_orthogonal
  if (present(need_S)) do_need_S = do_need_S .or. need_S

  do_need_H = optional_default(.true., need_H)
  do_no_S_spin = optional_default(.false., no_S_spin)

  if (do_need_S .and. this%S%N == 0) then
    call Initialise(this%S, this%N, this%n_matrices, this%complex_matrices)
  endif

  call fill_these_matrices(this, at, do_need_H, this%H, do_need_S, this%S, do_no_S_spin)

end subroutine

subroutine TBSystem_fill_these_matrices(this, at, do_H, H, do_S, S, no_S_spin, do_dipole, dipole)
  type (TBSystem), intent(inout) :: this
  type(Atoms), intent(in) :: at
  logical, intent(in), optional :: do_H
  type(TBMatrix), intent(inout), optional :: H
  logical, intent(in), optional :: do_S
  type(TBMatrix), intent(inout), optional :: S
  logical, intent(in), optional :: no_S_spin
  logical, intent(in), optional :: do_dipole
  type(TBMatrix), intent(inout), optional :: dipole(3)

  integer :: i, ji, j, ik, ii, jj, iii, jjj
  real(dp), allocatable :: block_H(:,:), block_S(:,:), block_H_up(:,:), block_H_down(:,:), block_dipole(:,:,:)
  complex(dp), allocatable :: block_H_z(:,:), block_S_z(:,:), block_H_z_phase(:,:), block_S_z_phase(:,:)
  complex(dp), allocatable :: block_dipole_z(:,:,:), block_dipole_z_phase(:,:,:)
  complex(dp), allocatable :: block_SO(:,:)
  real(dp) :: dv_hat(3), dv_hat_ij(3), dv_hat_ji(3), dv_mag
  integer :: shift(3)
  complex(dp) :: expf

  integer :: block_nr, block_nc

  logical :: u_do_S, u_do_S_block, u_do_H, u_no_S_spin, u_do_dipole

  u_do_S = optional_default(.false., do_S)
  u_do_H = optional_default(.false., do_H)
  u_no_S_spin = optional_default(.false., no_S_spin)
  u_do_dipole = optional_default(.false., do_dipole)

  u_do_S_block = u_do_S .or. this%scf%active

  if (u_do_H .and. .not. present(H)) call system_abort("Called TBSystem_fill_some_matrices with do_H but no H")
  if (u_do_S .and. .not. present(S)) call system_abort("Called TBSystem_fill_some_matrices with do_S but no S")
  if (u_do_dipole .and. .not. present(dipole)) call system_abort("Called TBSystem_fill_some_matrices with do_dipole but no dipole")

  if (u_do_H) call Zero(H)
  if (u_do_S) call Zero(S)
  if (u_do_dipole) then
    call Zero(dipole(1))
    call Zero(dipole(2))
    call Zero(dipole(3))
    call check_dipole_model_consistency(this%dipole_model, this%tbmodel, at%Z)
  endif

  allocate(block_H(this%max_block_size, this%max_block_size))
  allocate(block_S(this%max_block_size, this%max_block_size))
  allocate(block_dipole(this%max_block_size, this%max_block_size,3))
  if (this%complex_matrices) then
    allocate(block_H_z(this%max_block_size, this%max_block_size))
    allocate(block_S_z(this%max_block_size, this%max_block_size))
    allocate(block_dipole_z(this%max_block_size, this%max_block_size,3))
    if (this%kpoints%non_gamma) then
      allocate(block_H_z_phase(this%max_block_size, this%max_block_size))
      allocate(block_S_z_phase(this%max_block_size, this%max_block_size))
      allocate(block_dipole_z_phase(this%max_block_size, this%max_block_size,3))
    endif
  endif
  if (this%noncollinear) then
    allocate(block_H_up(this%max_block_size, this%max_block_size))
    allocate(block_H_down(this%max_block_size, this%max_block_size))
    if (this%SO%active) then
      allocate(block_SO(this%max_block_size, this%max_block_size))
      call check_spin_orbit_coupling_consistency(this%SO, this%tbmodel, at%Z)
    endif
  else
    if (this%SO%active) then
      call print("WARNING: fill_these_matrices got active spin-orbit coupling, but noncollinear is off")
    endif
  endif

  do i=1, at%N
    block_nr = this%first_orb_of_atom(i+1)-this%first_orb_of_atom(i)
    do ji=1, n_neighbours(at, i)
      j = neighbour(at, i, ji, dv_mag, shift = shift)

      ! do our own direction cosine, since atom_neighbour() result isn't antisymmetric enough
      if (i == j) then
	dv_hat = at%lattice .mult. shift
      else
	dv_hat_ij = at%pos(:,j) - at%pos(:,i) + (at%lattice .mult. shift)
	dv_hat_ji = at%pos(:,i) - at%pos(:,j) - (at%lattice .mult. shift)
	dv_hat = 0.5_dp*(dv_hat_ij - dv_hat_ji)
      endif
      if (maxval(abs(dv_hat)) .feq. 0.0_dp) then
	dv_hat = 0.0_dp
      else
	dv_hat = dv_hat / norm(dv_hat)
      endif

      block_nc = this%first_orb_of_atom(j+1)-this%first_orb_of_atom(j)

      if (this%noncollinear) then
	call get_HS_blocks(this%tbmodel, at, i, j, dv_hat, dv_mag, block_H_up, block_S, i_mag=1)
	call get_HS_blocks(this%tbmodel, at, i, j, dv_hat, dv_mag, block_H_down, block_S, i_mag=2)
	block_H_z = 0.0_dp
	block_H_z(1:block_nr:2,1:block_nc:2) = 0.5_dp*(block_H_up+block_H_down)
	block_H_z(2:block_nr:2,2:block_nc:2) = 0.5_dp*(block_H_up+block_H_down)
	block_S_z = 0.0_dp
	block_S_z(1:block_nr:2,1:block_nc:2) = block_S
	block_S_z(2:block_nr:2,2:block_nc:2) = block_S
	if (u_no_S_spin) then
	  block_S_z(1:block_nr:2,2:block_nc:2) = block_S
	  block_S_z(2:block_nr:2,1:block_nc:2) = block_S
	endif
	if (this%SO%active .and. (i == j) .and. (dv_mag .feq. 0.0_dp)) then
	  call get_SO_block(this%SO, this%tbmodel, at%Z(i), block_SO)
	  block_H_z = block_H_z + block_SO
	endif
	if (u_do_dipole) call system_abort("No dipole for noncollinear magnetism yet")
      else ! not noncollinear
	call get_HS_blocks(this%tbmodel, at, i, j, dv_hat, dv_mag, block_H, block_S)
	if (this%complex_matrices) then
	  block_H_z = block_H
	  block_S_z = block_S
	endif
	if (u_do_dipole) then
	  call get_dipole_block(this%dipole_model, at, i, j, dv_hat, dv_mag, block_dipole)
	  if (this%complex_matrices) then
	    block_dipole_z = block_dipole
	  endif
	endif
      endif

      if (u_do_H) then
	if (this%scf%active) then
	  do ii=1, this%first_orb_of_atom(i+1)-this%first_orb_of_atom(i)
	  do jj=1, this%first_orb_of_atom(j+1)-this%first_orb_of_atom(j)
	    if (this%complex_matrices) then
	      block_H_z(ii,jj) = block_H_z(ii,jj) + &
		0.5_dp*( this%scf%orb_local_pot(this%first_orb_of_atom(i)+ii-1) + &
			 this%scf%orb_local_pot(this%first_orb_of_atom(j)+jj-1) ) * block_S_z(ii,jj)
	    else
	      block_H(ii,jj) = block_H(ii,jj) + &
		0.5_dp*( this%scf%orb_local_pot(this%first_orb_of_atom(i)+ii-1) + &
			 this%scf%orb_local_pot(this%first_orb_of_atom(j)+jj-1) ) * block_S(ii,jj)
	    endif
	  end do
	  end do

	  if (this%noncollinear) then
	    do ii=1, block_nr, 2
	    do jj=1, block_nc, 2
		iii = (ii-1)/2+1
		jjj = (jj-1)/2+1
		call add_exch_field_local_pot(0.5_dp*(this%scf%orb_exch_field(1:3,this%first_orb_of_atom(i)+ii-1) + &
						      this%scf%orb_exch_field(1:3,this%first_orb_of_atom(j)+jj-1)), &
					      block_S(iii,jjj), block_H_z(ii:ii+1,jj:jj+1) )
	    end do
	    end do
	  endif

	endif
      endif

      if (this%complex_matrices) then
	if (this%kpoints%non_gamma) then
	  do ik=1, this%kpoints%N
	    expf = calc_phase(this%kpoints, ik, shift)
	    if (u_do_H) then
	      block_H_z_phase = block_H_z*expf
	      call add_block(H, ik, block_H_z_phase, block_nr, block_nc, &
		this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	    endif
	    if (u_do_S) then
	      block_S_z_phase = block_S_z*expf
	      call add_block(S, ik, block_S_z_phase, block_nr, block_nc, &
		this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	    endif
	    if (u_do_dipole) then
	      block_dipole_z_phase = block_dipole_z*expf
	      call add_block(dipole(1), ik, block_dipole_z_phase(:,:,1), block_nr, block_nc, &
		this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	      call add_block(dipole(2), ik, block_dipole_z_phase(:,:,2), block_nr, block_nc, &
		this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	      call add_block(dipole(3), ik, block_dipole_z_phase(:,:,3), block_nr, block_nc, &
		this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	    endif
	  end do
	else ! gamma-point
	  if (u_do_H) then
	    call add_block(H, 1, block_H_z, block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	  endif
	  if (u_do_S) then
	    call add_block(S, 1, block_S_z, block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	  endif
	  if (u_do_dipole) then
	    call add_block(dipole(1), 1, block_dipole_z(:,:,1), block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	    call add_block(dipole(2), 1, block_dipole_z(:,:,2), block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	    call add_block(dipole(3), 1, block_dipole_z(:,:,3), block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	  endif
	endif ! non gamma
      else ! real matrices
	do ik=1, this%n_matrices
	  if (u_do_H) then
	    call add_block(H, ik, block_H, block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	  endif
	  if (u_do_S) then
	    call add_block(S, ik, block_S, block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	  endif
	  if (u_do_dipole) then
	    call add_block(dipole(1), ik, block_dipole(:,:,1), block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	    call add_block(dipole(2), ik, block_dipole(:,:,2), block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	    call add_block(dipole(3), ik, block_dipole(:,:,3), block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	  endif
	end do
      endif ! complex matrices

    end do
  end do

  if (allocated(block_H)) deallocate(block_H)
  if (allocated(block_H_up)) deallocate(block_H_up)
  if (allocated(block_H_down)) deallocate(block_H_down)
  if (allocated(block_S)) deallocate(block_S)
  if (allocated(block_dipole)) deallocate(block_dipole)
  if (allocated(block_H_z)) deallocate(block_H_z)
  if (allocated(block_S_z)) deallocate(block_S_z)
  if (allocated(block_dipole_z)) deallocate(block_dipole_z)
  if (allocated(block_H_z_phase)) deallocate(block_H_z_phase)
  if (allocated(block_S_z_phase)) deallocate(block_S_z_phase)
  if (allocated(block_dipole_z_phase)) deallocate(block_dipole_z_phase)

end subroutine TBSystem_fill_these_matrices

subroutine add_exch_field_local_pot(exch_field, block_S, block_H_z)
  real(dp), intent(in) :: exch_field(3)
  real(dp), intent(in) :: block_S
  complex(dp), intent(inout) :: block_H_z(2,2)

  complex(dp) :: E_dot_sigma(2,2)

  E_dot_sigma(:,:) = exch_field(1)*pauli_sigma(:,:,1) + exch_field(2)*pauli_sigma(:,:,2) + exch_field(3)*pauli_sigma(:,:,3)

  block_H_z = block_H_z + block_S*E_dot_sigma

end subroutine add_exch_field_local_pot

subroutine TBSystem_fill_dmatrices(this, at, at_ind, need_S, dense, diag_mask, offdiag_mask)
  type (TBSystem), intent(inout) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: at_ind
  logical, intent(in), optional :: need_S, dense
  logical, intent(out), optional, target :: diag_mask(:), offdiag_mask(:)

  integer :: i, ji, j, ik, ii, jj
  real(dp), allocatable :: block_dH(:,:,:), block_dS(:,:,:)
  complex(dp), allocatable :: block_dH_z(:,:,:), block_dS_z(:,:,:)
  complex(dp), allocatable :: block_dH_z_phase(:,:,:), block_dS_z_phase(:,:,:)
  real(dp) :: dv_hat(3), dv_mag
  integer :: shift(3)

  logical, pointer :: d_mask(:), od_mask(:)

  logical :: do_S, do_S_block
  logical :: block_active

  complex(dp) :: expf

  integer block_nr, block_nc

  if (at_ind >= 0) then
    if (at_ind < 1 .or. at_ind > at%N) &
      call system_abort("Called TBSystem_fill_dmatrices for atom with at_ind " // at_ind // " > 0, but out of range " // 1 // " " // at%N)
  else
    if (at_ind < -3 .or. at_ind > -1) &
      call system_abort("Called TBSystem_fill_dmatrices for virial with at_ind " // at_ind // " < 0, but out of range " // -3 // " " // -3)
  end if

  do_S = .not. this%tbmodel%is_orthogonal
  if (present(need_S)) do_S = do_S .or. need_S
  do_S_block = do_S .or. this%scf%active

  allocate(block_dH(this%max_block_size, this%max_block_size, 3))
  allocate(block_dS(this%max_block_size, this%max_block_size, 3))
  if (this%complex_matrices) then
    allocate(block_dH_z(this%max_block_size, this%max_block_size, 3))
    allocate(block_dS_z(this%max_block_size, this%max_block_size, 3))
    if (this%kpoints%non_gamma) then
      allocate(block_dH_z_phase(this%max_block_size, this%max_block_size, 3))
      allocate(block_dS_z_phase(this%max_block_size, this%max_block_size, 3))
    endif
  endif

  if (present(diag_mask)) then
    d_mask => diag_mask
  else
    allocate(d_mask(at%N))
  endif
  if (present(offdiag_mask)) then
    od_mask => offdiag_mask
  else
    allocate(od_mask(at%N))
  endif
  call get_dHS_masks(this%tbmodel, at, at_ind, d_mask, od_mask)

  call Zero(this%dH(1), d_mask, od_mask)
  call Zero(this%dH(2), d_mask, od_mask)
  call Zero(this%dH(3), d_mask, od_mask)
  if (do_S) then
    call Zero(this%dS(1), d_mask, od_mask)
    call Zero(this%dS(2), d_mask, od_mask)
    call Zero(this%dS(3), d_mask, od_mask)
  endif

  if (this%noncollinear) call system_abort("no noncollinear forces yet")

  do i=1, at%N
    block_nr = this%first_orb_of_atom(i+1)-this%first_orb_of_atom(i)
    do ji=1, n_neighbours(at, i)
	j = neighbour(at, i, ji, dv_mag, cosines = dv_hat, shift = shift)
	block_nc = this%first_orb_of_atom(j+1)-this%first_orb_of_atom(j)

      if ((i == j .and. .not. d_mask(i)) .or. &
	  (i /= j .and. .not. od_mask(i) .and. .not. od_mask(j))) then
	cycle
      endif

      block_active = get_dHS_blocks(this%tbmodel, at, i, j, dv_hat, dv_mag, at_ind, block_dH, block_dS)
      if (block_active .and. this%complex_matrices) then
	block_dH_z = block_dH
	if (do_S_block) block_dS_z = block_dS
      endif


      if (block_active .and. this%scf%active) then
	do ii=1, this%first_orb_of_atom(i+1)-this%first_orb_of_atom(i)
	do jj=1, this%first_orb_of_atom(j+1)-this%first_orb_of_atom(j)
	  if (this%complex_matrices) then
	    block_dH_z(ii,jj,:) = block_dH_z(ii,jj,:) + &
	      0.5_dp*( this%scf%orb_local_pot(this%first_orb_of_atom(i)+ii-1) + &
		       this%scf%orb_local_pot(this%first_orb_of_atom(j)+jj-1) ) * block_dS_z(ii,jj,:)
	  else
	    block_dH(ii,jj,:) = block_dH(ii,jj,:) + &
	      0.5_dp*( this%scf%orb_local_pot(this%first_orb_of_atom(i)+ii-1) + &
		       this%scf%orb_local_pot(this%first_orb_of_atom(j)+jj-1) ) * block_dS(ii,jj,:)
	  endif
	end do
	end do
	if (this%noncollinear) call system_abort("No matrix derivatives for noncollinear yet")
	! need to add noncollinear stuff here !!!
      endif

      if (block_active) then
	if (this%complex_matrices) then
	  if (this%kpoints%non_gamma) then
	    do ik=1, this%kpoints%N
	      expf = calc_phase(this%kpoints, ik, shift)
	      block_dH_z_phase = block_dH_z*expf
	      call add_block(this%dH(1), ik, block_dH_z_phase(:,:,1), block_nr, block_nc, &
		this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	      call add_block(this%dH(2), ik, block_dH_z_phase(:,:,2), block_nr, block_nc, &
		this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	      call add_block(this%dH(3), ik, block_dH_z_phase(:,:,3), block_nr, block_nc, &
		this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	      if (do_S) then
		block_dS_z_phase = block_dS_z*expf
		call add_block(this%dS(1), ik, block_dS_z_phase(:,:,1), block_nr, block_nc, &
		  this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
		call add_block(this%dS(2), ik, block_dS_z_phase(:,:,2), block_nr, block_nc, &
		  this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
		call add_block(this%dS(3), ik, block_dS_z_phase(:,:,3), block_nr, block_nc, &
		  this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	      endif
	    end do
	  else ! gamma
	    call add_block(this%dH(1), 1, block_dH_z(:,:,1), block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	    call add_block(this%dH(2), 1, block_dH_z(:,:,2), block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	    call add_block(this%dH(3), 1, block_dH_z(:,:,3), block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	    if (do_S) then
	      call add_block(this%dS(1), 1, block_dS_z(:,:,1), block_nr, block_nc, &
		this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	      call add_block(this%dS(2), 1, block_dS_z(:,:,2), block_nr, block_nc, &
		this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	      call add_block(this%dS(3), 1, block_dS_z(:,:,3), block_nr, block_nc, &
		this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	    endif
	  endif ! non_gamma
	else ! real matrices
	  call add_block(this%dH(1), 1, block_dH(:,:,1), block_nr, block_nc, &
	    this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	  call add_block(this%dH(2), 1, block_dH(:,:,2), block_nr, block_nc, &
	    this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	  call add_block(this%dH(3), 1, block_dH(:,:,3), block_nr, block_nc, &
	    this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	  if (do_S) then
	    call add_block(this%dS(1), 1, block_dS(:,:,1), block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	    call add_block(this%dS(2), 1, block_dS(:,:,2), block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	    call add_block(this%dS(3), 1, block_dS(:,:,3), block_nr, block_nc, &
	      this%first_orb_of_atom(i), this%first_orb_of_atom(j), i, j)
	  endif
	endif ! complex_matrices
      endif ! block_active

    end do
  end do

  if (allocated(block_dH)) deallocate(block_dH)
  if (allocated(block_dS)) deallocate(block_dS)
  if (allocated(block_dH_z)) deallocate(block_dH_z)
  if (allocated(block_dS_z)) deallocate(block_dS_z)
  if (allocated(block_dH_z_phase)) deallocate(block_dH_z_phase)
  if (allocated(block_dS_z_phase)) deallocate(block_dS_z_phase)

  if (.not. present(diag_mask)) deallocate(d_mask)
  if (.not. present(offdiag_mask)) deallocate(od_mask)

end subroutine

subroutine TBSystem_ksum_mat_d(this, tbm, m)
  type(TBSystem), intent(in) :: this
  type(TBMatrix), intent(in) :: tbm
  type(MatrixD), intent(inout) :: m

  call sum_matrices(tbm, this%kpoints%weights, m)
end subroutine

subroutine TBSystem_ksum_atom_orbital_sum_mat_d(this, tbm, m)
  type(TBSystem), intent(in) :: this
  type(TBMatrix), intent(in) :: tbm
  type(MatrixD), intent(inout) :: m

  type(MatrixD) :: t_m

  call Initialise(t_m, this%N)
  call Zero(t_m)

  call sum_matrices(tbm, this%kpoints%weights, t_m)
  call ksum_distrib_inplace(this%kpoints, t_m%data)

  call atom_orbital_sum_mat(this, t_m%data, m%data)
end subroutine

function TBSystem_atom_orbital_sum_real2(this,a)
  type(TBSystem),    intent(in)           :: this
  real(dp), intent(in) :: a(:,:)
  real(dp) :: TBSystem_atom_orbital_sum_real2(this%N_atoms,size(a,2))

  integer i

  if (size(a,1) /= this%N) then
    call system_abort("Called TBSystem_atom_orbital_sum_real2 with wrong size array")
  endif

  do i=1, this%N_atoms
    TBSystem_atom_orbital_sum_real2(i,:) = sum(a(this%first_orb_of_atom(i):this%first_orb_of_atom(i+1)-1,:),1)
  end do
end function TBSystem_atom_orbital_sum_real2


function TBSystem_atom_orbital_sum_real1(this,a)
  type(TBSystem),    intent(in)           :: this
  real(dp), intent(in) :: a(:)
  real(dp) :: TBSystem_atom_orbital_sum_real1(this%N_atoms)

  integer i

  if (size(a) /= this%N) then
    call system_abort("Called TBSystem_atom_orbital_sum_real1 with wrong size array")
  endif

  do i=1, this%N_atoms
    TBSystem_atom_orbital_sum_real1(i) = sum(a(this%first_orb_of_atom(i):this%first_orb_of_atom(i+1)-1))
  end do
end function TBSystem_atom_orbital_sum_real1

function TBSystem_atom_orbital_sum_complex2(this,a)
  type(TBSystem),    intent(in)           :: this
  complex(dp), intent(in) :: a(:,:)
  complex(dp) :: TBSystem_atom_orbital_sum_complex2(this%N_atoms, size(a,2))

  integer i

  if (size(a,1) /= this%N) then
    call system_abort("Called TBSystem_atom_orbital_sum_complex2 with wrong size array")
  endif

  do i=1, this%N_atoms
    TBSystem_atom_orbital_sum_complex2(i,:) = sum(a(this%first_orb_of_atom(i):this%first_orb_of_atom(i+1)-1,:),1)
  end do
end function TBSystem_atom_orbital_sum_complex2

function TBSystem_atom_orbital_spread_real1(this,a)
  type(TBSystem),    intent(in)           :: this
  real(dp), intent(in) :: a(:)
  real(dp) :: TBSystem_atom_orbital_spread_real1(this%N)

  integer i

  if (size(a) /= this%N_atoms) then
    call system_abort("Called TBSystem_atom_orbital_spread_real1 with wrong size array a " // &
      size(a) // " this%N_atoms " // this%N_atoms)
  endif

  do i=1, this%N_atoms
    TBSystem_atom_orbital_spread_real1(this%first_orb_of_atom(i):this%first_orb_of_atom(i+1)-1) = a(i)
  end do
end function TBSystem_atom_orbital_spread_real1

subroutine TBSystem_atom_orbital_sum_mat_r(this,a, out)
  type(TBSystem),    intent(in)           :: this
  real(dp), intent(in) :: a(:,:)
  real(dp), intent(out) :: out(:,:)

  integer i, j

  if (size(a,1) /= this%N .or. size(a,2) /= this%N) then
    call system_abort("Called TBSystem_atom_orbital_sum_mat_r with wrong size input")
  endif
  if (size(out,1) /= this%N_atoms .or. size(out,2) /= this%N_atoms) then
    call system_abort("Called TBSystem_atom_orbital_sum_mat_r with wrong size output")
  endif

  do j=1, this%N_atoms
    do i=1, this%N_atoms
      out(i,j) = sum(a(this%first_orb_of_atom(i):this%first_orb_of_atom(i+1)-1, &
                       this%first_orb_of_atom(j):this%first_orb_of_atom(j+1)-1))
    end do
  end do
end subroutine TBSystem_atom_orbital_sum_mat_r

subroutine TBSystem_atom_orbital_spread_mat_r(this,a, out)
  type(TBSystem),    intent(in)           :: this
  real(dp), intent(in) :: a(:,:)
  real(dp), intent(out) :: out(:,:)

  integer i, j

  if (size(a,1) /= this%N_atoms .or. size(a,2) /= this%N_atoms) then
    call system_abort("Called TBSystem_atom_orbital_spread_mat_r with wrong size input")
  endif
  if (size(out,1) /= this%N .or. size(out,2) /= this%N) then
    call system_abort("Called TBSystem_atom_orbital_spread_mat_r with wrong size output")
  endif

  do j=1, this%N_atoms
    do i=1, this%N_atoms
      out(this%first_orb_of_atom(i):this%first_orb_of_atom(i+1)-1, &
          this%first_orb_of_atom(j):this%first_orb_of_atom(j+1)-1) = a(i,j)
    end do
  end do
end subroutine TBSystem_atom_orbital_spread_mat_r

function TBSystem_manifold_orbital_sum_real2(this,a)
  type(TBSystem),    intent(in)           :: this
  real(dp), intent(in) :: a(:,:)
  real(dp) :: TBSystem_manifold_orbital_sum_real2(size(a,1),this%N_manifolds)

  integer i

  if (size(a,2) /= this%N) then
    call system_abort("Called TBSystem_manifold_orbital_sum_real2 with wrong size array shape(a) " // shape(a) // " this%N " // (this%N))
  endif

  do i=1, this%N_manifolds
    TBSystem_manifold_orbital_sum_real2(:,i) = sum(a(:,this%first_orb_of_manifold(i):this%first_orb_of_manifold(i+1)-1),2)
  end do
end function TBSystem_manifold_orbital_sum_real2

subroutine TBSystem_Print(this,file)
  type(TBSystem),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  if (current_verbosity() < PRINT_NORMAL) return

  call Print ('TBSystem : N N_atoms ' // this%N // " " //  this%N_atoms, file=file)

  call Print (this%tbmodel, file=file)

  if (this%kpoints%non_gamma) then
    call Print (this%kpoints, file=file)
  endif
  if (this%scf%active) then
    call Print (this%scf, file=file)
  endif
  if (this%dipole_model%active) then
    call Print (this%dipole_model, file=file)
  endif

  if (this%SO%active) then
    call Print (this%SO, file=file)
  endif

  if (this%N > 0) then
    call verbosity_push_decrement(PRINT_NERD)
      call Print("first_orb_of_atom", file=file)
      call Print(this%first_orb_of_atom, file=file)

      call Print("first_manifold_of_atom", file=file)
      call Print(this%first_manifold_of_atom, file=file)

      call Print("first_orb_of_manifold", file=file)
      call Print(this%first_orb_of_manifold, file=file)

      call print("H", file=file)
      call Print (this%H, file=file)
      call print("S", file=file)
      call Print (this%S, file=file)
    call verbosity_pop()
  endif

end subroutine TBSystem_Print

subroutine Self_Consistency_Initialise_str(this, args_str, param_str)
  type(Self_Consistency), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  call Finalise(this)

  call read_params_xml(this, param_str)

  call set_type(this, args_str)

end subroutine Self_Consistency_Initialise_str

subroutine TB_Dipole_Model_Initialise_str(this, args_str, param_str)
  type(TB_Dipole_Model), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  call Finalise(this)

  call read_params_xml(this, param_str)

end subroutine TB_Dipole_Model_Initialise_str

subroutine TB_Spin_Orbit_Coupling_Initialise_str(this, args_str, param_str)
  type(TB_Spin_Orbit_Coupling), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call read_params_xml(this, param_str)

  call initialise(params)
  call param_register(params, 'spin_orbit_coupling', 'F', this%active, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='TB_Spin_Orbit_Coupling_Initialise_str args_str')) then
    call system_abort("TB_Spin_Orbit_Coupling_Initialise_str failed to parse args_str='"//trim(args_str)//"'")
  endif
  call finalise(params)

end subroutine TB_Spin_Orbit_Coupling_Initialise_str

subroutine Self_Consistency_set_type_str(this, args_str)
  type(Self_Consistency), intent(inout) :: this
  character(len=*), intent(in) :: args_str

  type(Dictionary) :: params
  real(dp) :: global_U
  logical :: is_SCF(8)
  integer :: i, i_term, N_terms

  call initialise(params)

  call param_register(params, 'SCF_LOCAL_U', 'F', is_SCF(SCF_LOCAL_U), help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_NONLOCAL_U_DFTB', 'F', is_SCF(SCF_NONLOCAL_U_DFTB), help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_NONLOCAL_U_NRL_TB', 'F', is_SCF(SCF_NONLOCAL_U_NRL_TB), help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_GLOBAL_U', 'F', is_SCF(SCF_GLOBAL_U), help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_LCN', 'F', is_SCF(SCF_LCN), help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_GCN', 'F', is_SCF(SCF_GCN), help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_SPIN_DIR', 'F', is_SCF(SCF_SPIN_DIR), help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_SPIN_STONER', 'F', is_SCF(SCF_SPIN_STONER), help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'GLOBAL_U', '-1.0', global_U, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_ALPHA', ''//this%alpha, this%alpha, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_W0', ''//this%w0, this%w0, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_MAX_ITER', ''//this%max_iter, this%max_iter, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_MIX_SIMPLE', ''//this%mix_simple, this%mix_simple, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_CONV_TOL', ''//this%conv_tol, this%conv_tol, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_MIX_SIMPLE_END_TOL', ''//this%mix_simple_end_tol, this%mix_simple_end_tol, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SCF_MIX_SIMPLE_START_TOL', ''//this%mix_simple_start_tol, this%mix_simple_start_tol, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Self_Consistency_set_type_str args_str')) then
    call system_abort("Self_Consistency_Initialise_str failed to parse args_str='"//trim(args_str)//"'")
  endif
  call finalise(params)

  if (count(is_SCF) == 0) return

  if (is_SCF(SCF_LCN) .or. is_SCF(SCF_GCN)) then
    if (count(is_SCF) /= 1) &
      call system_abort("Self_Consistency_set_type_str got SCF_LCN or SCF_GCN, but also some other sort of SCF")
  endif

  this%active = .true.

  N_terms = count(is_SCF)
  if (allocated(this%terms)) then
    do i_term=1, size(this%terms)
      call Finalise(this%terms(i_term))
    end do
    deallocate(this%terms)
  endif
  allocate(this%terms(N_terms))
  i_term = 0
  do i=1, size(is_SCF)
    if (is_SCF(i)) then
      i_term = i_term + 1
      this%terms(i_term)%type = i
      this%terms(i_term)%active = .true.
    end if
  end do

end subroutine Self_Consistency_set_type_str

subroutine Self_Consistency_Finalise(this)
  type(Self_Consistency), intent(inout) :: this

  integer :: i_term

  call Wipe(this)

  if (allocated(this%U)) deallocate(this%U)
  if (allocated(this%stoner_param)) deallocate(this%stoner_param)

  if (allocated(this%terms)) then
     do i_term=1, size(this%terms)
        call Finalise(this%terms(i_term))
     end do
     deallocate(this%terms)
  end if

  this%active = .false.

end subroutine Self_Consistency_Finalise

subroutine Self_Consistency_Term_Finalise(this)
  type(Self_Consistency_Term), intent(inout) :: this

  call wipe(this)

  this%type = SCF_NONE
  this%N = 0
  this%N_atoms = 0
  this%N_manifolds = 0

  this%global_U = 0.0_dp

  this%active = .false.

end subroutine Self_Consistency_Term_Finalise

subroutine Self_Consistency_Wipe(this)
  type(Self_Consistency), intent(inout) :: this

  integer :: i_term

  if (allocated(this%orb_local_pot)) deallocate(this%orb_local_pot)
  if (allocated(this%orb_exch_field)) deallocate(this%orb_exch_field)

  if (allocated(this%terms)) then
     do i_term=1, size(this%terms)
        call Wipe(this%terms(i_term))
     end do
  end if

  this%N = 0
  this%N_atoms = 0
  this%N_manifolds = 0
end subroutine Self_Consistency_Wipe

subroutine Self_Consistency_Term_Wipe(this)
  type(Self_Consistency_Term), intent(inout) :: this

  this%N = 0
  this%N_atoms = 0
  this%N_manifolds = 0
  this%n_dof = 0

  this%global_N = 0.0
  this%global_pot = 0.0

  if (allocated(this%U)) deallocate(this%U)
  if (allocated(this%spin_splitting)) deallocate(this%spin_splitting)
  if (allocated(this%stoner_param)) deallocate(this%stoner_param)
  if (allocated(this%atomic_n)) deallocate(this%atomic_n)
  if (allocated(this%atomic_n0)) deallocate(this%atomic_n0)
  if (allocated(this%manifold_mom)) deallocate(this%manifold_mom)
  if (allocated(this%atomic_local_pot)) deallocate(this%atomic_local_pot)
  if (allocated(this%datomic_local_pot_dr)) deallocate(this%datomic_local_pot_dr)
  if (allocated(this%gamma)) deallocate(this%gamma)
  if (allocated(this%dgamma_dr)) deallocate(this%dgamma_dr)

end subroutine Self_Consistency_Term_Wipe

subroutine TB_Dipole_Model_Finalise(this)
  type(TB_Dipole_Model), intent(inout) :: this

  this%n_types = 0
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%n_orb_sets)) deallocate(this%n_orb_sets)
  if (allocated(this%orb_set_type)) deallocate(this%orb_set_type)
  if (allocated(this%orb_set_phase)) deallocate(this%orb_set_phase)
  if (allocated(this%gaussian_width)) deallocate(this%gaussian_width)

  this%active = .false.

end subroutine TB_Dipole_Model_Finalise

subroutine TB_Spin_Orbit_Coupling_Finalise(this)
  type(TB_Spin_Orbit_Coupling), intent(inout) :: this

  this%n_types = 0
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%n_orb_sets)) deallocate(this%n_orb_sets)
  if (allocated(this%orb_set_type)) deallocate(this%orb_set_type)
  if (allocated(this%SO_param)) deallocate(this%SO_param)

  this%active = .false.

end subroutine TB_Spin_Orbit_Coupling_Finalise

subroutine TBSystem_scf_set_atomic_n_mom(this, atomic_n, atomic_mom, atomic_pot)
  type(TBSystem), intent(inout) :: this
  real(dp), intent(in), pointer :: atomic_n(:), atomic_mom(:,:), atomic_pot(:)

  integer :: i_term, i_at, i_man

  do i_term=1, size(this%scf%terms)
    if (allocated(this%scf%terms(i_term)%atomic_n)) then
      if (associated(atomic_n)) then
	this%scf%terms(i_term)%atomic_n = atomic_n
      else
	this%scf%terms(i_term)%atomic_n = this%scf%terms(i_term)%atomic_n0
      endif
    end if
    if (allocated(this%scf%terms(i_term)%manifold_mom)) then
      if (associated(atomic_mom)) then
	this%scf%terms(i_term)%manifold_mom = 0.0_dp
	do i_at=1, this%N_atoms
	  do i_man = this%first_manifold_of_atom(i_at), this%first_manifold_of_atom(i_at+1)-1
	    if (i_man == this%first_manifold_of_atom(i_at+1)-1) then
	      this%scf%terms(i_term)%manifold_mom(:,i_man) = atomic_mom(:,i_at)
	    else
	      this%scf%terms(i_term)%manifold_mom(:,i_man) = 1.0e-6_dp*atomic_mom(:,i_at)
	    endif
	  end do
	end do
      else
	this%scf%terms(i_term)%manifold_mom = 0.0_dp
      endif
    end if ! manifold_mom is allocated
    if (allocated(this%scf%terms(i_term)%atomic_local_pot)) then
      if (associated(atomic_pot)) then
	this%scf%terms(i_term)%atomic_local_pot = atomic_pot
      else
	this%scf%terms(i_term)%atomic_local_pot = 0.0_dp
      endif
    end if
  end do ! atomic_local_pot is allocated
end subroutine TBSystem_scf_set_atomic_n_mom

subroutine TBSystem_scf_set_global_N(this, global_at_weight, global_N)
  type(TBSystem), intent(inout) :: this
  real(dp), intent(in), pointer :: global_at_weight(:)
  real(dp), intent(in), optional :: global_N

  integer i_term

  ! if global weight isn't associated, we can't set global_N.
  ! If it's actually needed by some SCF terms, they'll complain
  if (.not. associated(global_at_weight)) return

  do i_term=1, size(this%scf%terms)
    if (present(global_N)) then
      this%scf%terms(i_term)%global_N = global_N
    else
      this%scf%terms(i_term)%global_N = sum(this%scf%terms(i_term)%atomic_n0*global_at_weight)
    endif
  end do
end subroutine TBSystem_scf_set_global_N

subroutine TBSystem_scf_get_atomic_n_mom(this, atomic_n, atomic_mom, atomic_pot)
  type(TBSystem), intent(inout) :: this
  real(dp), intent(out), pointer :: atomic_n(:), atomic_mom(:,:), atomic_pot(:)

  integer i_term, i_at, i_man
  logical got_atomic_n, got_manifold_mom, got_atomic_pot

  if (associated(atomic_n)) then
    atomic_n = 0.0_dp
    got_atomic_n = .false.
    do i_term=1, size(this%scf%terms)
      if (allocated(this%scf%terms(i_term)%atomic_n)) then
	if (got_atomic_n) &
	  call system_abort("TBSystem_scf_get_atomic_n_mom found atomic_n allocated in more than 1 term")
	atomic_n = this%scf%terms(i_term)%atomic_n
	got_atomic_n = .true.
      end if
    end do

    if (.not. got_atomic_n) &
      call print("WARNING: TBSystem_scf_get_atomic_n_mom was passed atomic_n but didn't find " // &
		 "atomic_n allocated in any terms", PRINT_ALWAYS)
  end if

  if (associated(atomic_mom)) then
    got_manifold_mom = .false.
    do i_term=1, size(this%scf%terms)
      if (allocated(this%scf%terms(i_term)%manifold_mom)) then
	if (got_manifold_mom) &
	  call system_abort("TBSystem_scf_get_atomic_n_mom found manifold_mom allocated in more than 1 term")
	do i_at=1, this%N_atoms
	  i_man = this%first_manifold_of_atom(i_at+1)-1
	  atomic_mom(:,i_at) = this%scf%terms(i_term)%manifold_mom(:,i_man)
	end do
	got_manifold_mom = .true.
      end if
    end do

    if (.not. got_manifold_mom) &
      call print("WARNING: TBSystem_scf_get_atomic_n_mom was passed atomic_mom but didn't find " // &
		"manifold_mom allocated in any terms", PRINT_ALWAYS)
  end if

  if (associated(atomic_pot)) then
    atomic_pot = 0.0_dp
    got_atomic_pot = .false.
    do i_term=1, size(this%scf%terms)
      if (allocated(this%scf%terms(i_term)%atomic_local_pot)) then
	if (got_atomic_pot) &
	  call system_abort("TBSystem_scf_get_atomic_n_mom found atomic_pot allocated in more than 1 term")
	atomic_pot = this%scf%terms(i_term)%atomic_local_pot
	got_atomic_pot = .true.
      end if
    end do

    if (.not. got_atomic_pot) &
      call print("WARNING: TBSystem_scf_get_atomic_n_mom was passed atomic_pot but didn't find " // &
		 "atomic_pot allocated in any terms", PRINT_ALWAYS)
  end if

end subroutine TBSystem_scf_get_atomic_n_mom

subroutine TBSystem_scf_get_global_N(this, global_N)
  type(TBSystem), intent(inout) :: this
  real(dp), intent(out) :: global_N

  integer i_term
  logical got_global_N

  got_global_N = .false.
  do i_term=1, size(this%scf%terms)
    if (this%scf%terms(i_term)%global_U > 0.0_dp) then
      if (got_global_N) &
	call system_abort("TBSystem_scf_get_global_N found global_N allocated in more than 1 term")
      global_N = this%scf%terms(i_term)%global_N
      got_global_N = .true.
    endif
  end do

  if (.not. got_global_N) &
    call system_abort("TBSystem_scf_get_global_N didn't find global_N allocated in any terms")
end subroutine TBSystem_scf_get_global_N

subroutine SC_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_self_consistency) then
    if (name == "self_consistency") then
      parse_in_self_consistency = .false.
    endif
  endif
end subroutine SC_endElement_handler

subroutine SC_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=1024) :: value
  integer :: max_stoner_manifolds, n_stoner_params
  integer :: atnum

  if (name == "self_consistency") then
    parse_in_self_consistency = .true.

    call QUIP_FoX_get_value(attributes, 'tolerance', value, status)
    if (status /= 0) then
      call print("Can't find tolerance in self_consistency", PRINT_ALWAYS)
    else
      read (value, *) parse_self_consistency%conv_tol
    endif

    call QUIP_FoX_get_value(attributes, 'global_U', value, status)
    if (status == 0) read (value, *) parse_self_consistency%global_U

    if (allocated(parse_self_consistency%U)) deallocate(parse_self_consistency%U)
    allocate(parse_self_consistency%U(size(ElementName)))
    parse_self_consistency%U = 0.0_dp

    call QUIP_FoX_get_value(attributes, 'max_stoner_manifolds', value, status)
    if (status == 0) then
      read (value, *) max_stoner_manifolds
      if (allocated(parse_self_consistency%stoner_param)) deallocate(parse_self_consistency%stoner_param)
      allocate(parse_self_consistency%stoner_param(max_stoner_manifolds,size(ElementName)))
      parse_self_consistency%stoner_param = 0.0_dp
    endif

  else if (parse_in_self_consistency .and. name == "U") then

    call QUIP_FoX_get_value(attributes, 'Z', value, status)
    if (status /= 0) call system_abort('Self_Consistency_read_params_xml cannot find Z attribute in U')
    read (value, *), atnum

    if (atnum <= 0 .or. atnum > size(parse_self_consistency%U)) &
      call system_abort ("SC parse xml atnum " // atnum // " out of range " // 1 // " " // size(parse_self_consistency%U))
    call QUIP_FoX_get_value(attributes, 'U', value, status)
    if (status /= 0) call system_abort('Self_Consistency_read_params_xml cannot find U attribute')
    read (value, *), parse_self_consistency%U(atnum)

  else if (parse_in_self_consistency .and. name == "stoner_params") then

    call QUIP_FoX_get_value(attributes, 'Z', value, status)
    if (status /= 0) call system_abort('Self_Consistency_read_params_xml cannot find Z attribute in stoner_params')
    read (value, *), atnum

    call QUIP_FoX_get_value(attributes, 'n_stoner_params', value, status)
    if (status == 0) then
      read (value, *), n_stoner_params
      if (n_stoner_params > 0) then
	call QUIP_FoX_get_value(attributes, 'stoner_params', value, status)
	if (status /= 0) &
	  call system_abort("parse self consistency got n_stoner_params " // n_stoner_params // " for Z=" // atnum  // " but can't find actual stoner_params")
	read (value, *) parse_self_consistency%stoner_param(1:n_stoner_params,atnum)
      end if
    end if
  endif

end subroutine SC_startElement_handler

subroutine Self_Consistency_read_params_xml(this, in)
  type(Self_Consistency), intent(inout), target :: this
  character(len=*) :: in

  type(xml_t) :: fxml

  if (len(trim(in)) <= 0) then
    return
  endif

  parse_in_self_consistency = .false.
  parse_self_consistency => this

  call open_xml_string(fxml, in)

  call parse(fxml, &
    startElement_handler = SC_startElement_handler, &
    endElement_handler = SC_endElement_handler)

  call close_xml_t(fxml)

end subroutine Self_Consistency_read_params_xml

subroutine DM_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_dipole_model) then
    if (name == "dipole_model") then
      parse_in_dipole_model = .false.
    endif
  endif
end subroutine DM_endElement_handler

subroutine DM_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=1024) :: value
  integer :: atnum, i_type, new_n_orb_sets, i
  character(len=1), allocatable :: orb_set_types_str(:)

  if (name == "dipole_model") then
    parse_in_dipole_model = .true.

    parse_dipole_model%n_types = 0

    if (allocated(parse_dipole_model%type_of_atomic_num)) deallocate(parse_dipole_model%type_of_atomic_num)
    allocate(parse_dipole_model%type_of_atomic_num(size(ElementName)))
    parse_dipole_model%type_of_atomic_num = 0

  else if (parse_in_dipole_model .and. name == "element") then
    parse_dipole_model%active = .true.
    i_type = parse_dipole_model%n_types + 1

    call QUIP_FoX_get_value(attributes, 'Z', value, status)
    if (status /= 0) call system_abort('TB_Dipole_Model_read_params_xml cannot find Z attribute in type ' // status)
    read (value, *), atnum
    if (atnum <= 0 .or. atnum > size(parse_dipole_model%type_of_atomic_num)) &
      call system_abort ("DM parse xml atnum " // atnum // " out of range " // 1 // " " // size(parse_dipole_model%type_of_atomic_num))

    parse_dipole_model%type_of_atomic_num(atnum) = i_type

    call QUIP_FoX_get_value(attributes, 'n_orb_sets', value, status)
    if (status /= 0) call system_abort('TB_Dipole_Model_read_params_xml cannot find n_orb_sets attribute ' // status)
    read (value, *) new_n_orb_sets

    call realloc_dipole_model(parse_dipole_model, i_type, new_n_orb_sets)

    parse_dipole_model%atomic_num(i_type) = atnum
    parse_dipole_model%n_orb_sets(i_type) = new_n_orb_sets

    allocate(orb_set_types_str(parse_dipole_model%n_orb_sets(i_type)))
    call QUIP_FoX_get_value(attributes, 'orb_set_types', value, status)
    if (status /= 0) call system_abort('TB_Dipole_Model_read_params_xml cannot find orb_set_types attribute ' // status)
    read (value, *) orb_set_types_str
    do i=1, parse_dipole_model%n_orb_sets(i_type)
      if (orb_set_types_str(i) == 's' .or. orb_set_types_str(i) == 'S') then
	parse_dipole_model%orb_set_type(i,i_type) = ORB_S
      else if (orb_set_types_str(i) == 'p' .or. orb_set_types_str(i) == 'P') then
	parse_dipole_model%orb_set_type(i,i_type) = ORB_P
      else if (orb_set_types_str(i) == 'd' .or. orb_set_types_str(i) == 'D') then
	parse_dipole_model%orb_set_type(i,i_type) = ORB_D
      endif
    end do

    call QUIP_FoX_get_value(attributes, 'gaussian_widths', value, status)
    if (status /= 0) call system_abort('TB_Dipole_Model_read_params_xml cannot find gaussian_widths attribute ' // status)
    read (value, *) parse_dipole_model%gaussian_width(1:parse_dipole_model%n_orb_sets(i_type),i_type)

    call QUIP_FoX_get_value(attributes, 'orb_set_phases', value, status)
    if (status /= 0) call system_abort('TB_Dipole_Model_read_params_xml cannot find orb_set_phase attribute ' // status)
    read (value, *) parse_dipole_model%orb_set_phase(1:parse_dipole_model%n_orb_sets(i_type),i_type)

  endif 

end subroutine DM_startElement_handler

subroutine SO_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_spin_orbit_coupling) then
    if (name == "spin_orbit_coupling") then
      parse_in_spin_orbit_coupling = .false.
    endif
  endif
end subroutine SO_endElement_handler

subroutine SO_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=1024) :: value
  integer :: atnum, i_type, new_n_orb_sets, i
  character(len=1), allocatable :: orb_set_types_str(:)

  if (name == "spin_orbit_coupling") then
    parse_in_spin_orbit_coupling = .true.

    parse_spin_orbit_coupling%n_types = 0

    if (allocated(parse_spin_orbit_coupling%type_of_atomic_num)) deallocate(parse_spin_orbit_coupling%type_of_atomic_num)
    allocate(parse_spin_orbit_coupling%type_of_atomic_num(size(ElementName)))
    parse_spin_orbit_coupling%type_of_atomic_num = 0

  else if (parse_in_spin_orbit_coupling .and. name == "element") then
    parse_spin_orbit_coupling%active = .true.
    i_type = parse_spin_orbit_coupling%n_types + 1

    call QUIP_FoX_get_value(attributes, 'Z', value, status)
    if (status /= 0) call system_abort('TB_Dipole_Model_read_params_xml cannot find Z attribute in type ' // status)
    read (value, *), atnum
    if (atnum <= 0 .or. atnum > size(parse_spin_orbit_coupling%type_of_atomic_num)) &
      call system_abort ("SO parse xml atnum " // atnum // " out of range " // 1 // " " // size(parse_spin_orbit_coupling%type_of_atomic_num))

    parse_spin_orbit_coupling%type_of_atomic_num(atnum) = i_type

    call QUIP_FoX_get_value(attributes, 'n_orb_sets', value, status)
    if (status /= 0) call system_abort('TB_Dipole_Model_read_params_xml cannot find n_orb_sets attribute ' // status)
    read (value, *) new_n_orb_sets

    call realloc_spin_orbit_coupling(parse_spin_orbit_coupling, i_type, new_n_orb_sets)

    parse_spin_orbit_coupling%atomic_num(i_type) = atnum
    parse_spin_orbit_coupling%n_orb_sets(i_type) = new_n_orb_sets

    allocate(orb_set_types_str(parse_spin_orbit_coupling%n_orb_sets(i_type)))
    call QUIP_FoX_get_value(attributes, 'orb_set_types', value, status)
    if (status /= 0) call system_abort('TB_Dipole_Model_read_params_xml cannot find orb_set_types attribute ' // status)
    read (value, *) orb_set_types_str
    do i=1, parse_spin_orbit_coupling%n_orb_sets(i_type)
      if (orb_set_types_str(i) == 's' .or. orb_set_types_str(i) == 'S') then
	parse_spin_orbit_coupling%orb_set_type(i,i_type) = ORB_S
      else if (orb_set_types_str(i) == 'p' .or. orb_set_types_str(i) == 'P') then
	parse_spin_orbit_coupling%orb_set_type(i,i_type) = ORB_P
      else if (orb_set_types_str(i) == 'd' .or. orb_set_types_str(i) == 'D') then
	parse_spin_orbit_coupling%orb_set_type(i,i_type) = ORB_D
      endif
    end do

    call QUIP_FoX_get_value(attributes, 'SO_params', value, status)
    if (status /= 0) call system_abort('TB_Dipole_Model_read_params_xml cannot find SO_params attribute ' // status)
    read (value, *) parse_spin_orbit_coupling%SO_param(1:parse_spin_orbit_coupling%n_orb_sets(i_type),i_type)

  endif 

end subroutine SO_startElement_handler

subroutine TB_Dipole_Model_read_params_xml(this, in)
  type(TB_Dipole_Model), intent(inout), target :: this
  character(len=*) :: in

  type(xml_t) :: fxml

  if (len(trim(in)) <= 0) then
    return
  endif

  parse_in_dipole_model = .false.
  parse_dipole_model => this

  call open_xml_string(fxml, in)

  call parse(fxml, &
    startElement_handler = DM_startElement_handler, &
    endElement_handler = DM_endElement_handler)

  call close_xml_t(fxml)

end subroutine

subroutine TB_Spin_Orbit_Coupling_read_params_xml(this, in)
  type(TB_Spin_Orbit_Coupling), intent(inout), target :: this
  character(len=*) :: in

  type(xml_t) :: fxml

  if (len(trim(in)) <= 0) then
    return
  endif

  parse_in_spin_orbit_coupling = .false.
  parse_spin_orbit_coupling => this

  call open_xml_string(fxml, in)

  call parse(fxml, &
    startElement_handler = SO_startElement_handler, &
    endElement_handler = SO_endElement_handler)

  call close_xml_t(fxml)

end subroutine

subroutine Self_Consistency_setup_system(this, N, tbm, at_N, at_Z, N_manifolds)
  type(Self_Consistency), intent(inout) :: this
  integer, intent(in) :: N
  type(TBModel), intent(in) :: tbm
  integer, intent(in) :: at_N, at_Z(:), N_manifolds

  integer :: i_term

  call wipe(this)


  this%N = N
  this%N_atoms = at_N
  this%N_manifolds = N_manifolds

  allocate(this%orb_local_pot(this%N))
  if (allocated(this%terms)) then
     do i_term=1, size(this%terms)
        if (this%terms(i_term)%type == SCF_SPIN_DIR .or. this%terms(i_term)%type == SCF_SPIN_STONER) then
           allocate(this%orb_exch_field(3,this%N))
           exit
        endif
     end  do

     do i_term=1, size(this%terms)
        call setup_system(this%terms(i_term), this%N, tbm, at_N, at_Z, this%global_U, this%N_manifolds)
     end do
  end if

end subroutine Self_Consistency_setup_system

subroutine Self_Consistency_Term_setup_system(this, N, tbm, at_N, at_Z, global_U, N_manifolds)
  type(Self_Consistency_Term), intent(inout) :: this
  integer, intent(in) :: N
  type(TBModel), intent(in) :: tbm
  integer, intent(in) :: at_N, at_Z(:)
  real(dp), intent(in) :: global_U
  integer, intent(in) :: N_manifolds

  integer :: i

  if (allocated(this%U)) deallocate(this%U)
  if (allocated(this%stoner_param)) deallocate(this%stoner_param)
  if (allocated(this%spin_splitting)) deallocate(this%spin_splitting)
  if (allocated(this%atomic_n)) deallocate(this%atomic_n)
  if (allocated(this%atomic_n0)) deallocate(this%atomic_n0)
  if (allocated(this%manifold_mom)) deallocate(this%manifold_mom)
  if (allocated(this%atomic_local_pot)) deallocate(this%atomic_local_pot)
  if (allocated(this%gamma)) deallocate(this%gamma)
  if (allocated(this%dgamma_dr)) deallocate(this%dgamma_dr)

  if (.not. this%active) return

  this%N = N
  this%N_atoms = at_N
  this%N_manifolds = N_manifolds
  select case (this%type)
    case(SCF_NONE)
      this%n_dof = 0
      return
    case(SCF_GCN)
      this%n_dof = 1
      allocate(this%atomic_n0(this%N_atoms))
      do i=1, at_N
	this%atomic_n0(i) = n_elecs_of_Z(tbm,at_Z(i))
      end do
      this%global_pot = 0.0_dp
      return
    case(SCF_GLOBAL_U)
      this%n_dof = 1
      this%global_U = global_U
      allocate(this%atomic_n0(this%N_atoms))
      do i=1, at_N
	this%atomic_n0(i) = n_elecs_of_Z(tbm,at_Z(i))
      end do
      return
    case(SCF_LCN)
      this%n_dof = this%N_atoms
      allocate(this%atomic_n0(this%N_atoms))
      do i=1, at_N
	this%atomic_n0(i) = n_elecs_of_Z(tbm,at_Z(i))
      end do
      allocate(this%atomic_n(this%N_atoms))
      allocate(this%atomic_local_pot(this%N_atoms))
      this%atomic_local_pot = 0.0_dp
      return
    case(SCF_LOCAL_U)
      this%n_dof = this%N_atoms
      allocate(this%U(this%N_atoms))
      this%U = 0.0_dp
      allocate(this%atomic_n(this%N_atoms))
      allocate(this%atomic_n0(this%N_atoms))
      do i=1, at_N
	this%atomic_n0(i) = n_elecs_of_Z(tbm,at_Z(i))
      end do
      return
    case(SCF_NONLOCAL_U_DFTB)
      this%n_dof = this%N_atoms
      allocate(this%gamma(this%N_atoms,this%N_atoms))
      this%gamma = 0.0_dp
      allocate(this%atomic_n(this%N_atoms))
      allocate(this%atomic_n0(this%N_atoms))
      do i=1, at_N
	this%atomic_n0(i) = n_elecs_of_Z(tbm,at_Z(i))
      end do
    case(SCF_NONLOCAL_U_NRL_TB)
      this%n_dof = this%N_atoms
      allocate(this%gamma(this%N_atoms,this%N_atoms))
      this%gamma = 0.0_dp
      allocate(this%atomic_n(this%N_atoms))
      allocate(this%atomic_n0(this%N_atoms))
      do i=1, at_N
	this%atomic_n0(i) = n_elecs_of_Z(tbm,at_Z(i))
      end do
      return
    case(SCF_SPIN_DIR)
      this%n_dof = 3*this%N_manifolds
      allocate(this%manifold_mom(3,this%N_manifolds))
      allocate(this%spin_splitting(this%N_manifolds))
      this%spin_splitting = 0.0_dp
      return
    case(SCF_SPIN_STONER)
      this%n_dof = 3*this%N_manifolds
      allocate(this%manifold_mom(3,this%N_manifolds))
      allocate(this%stoner_param(this%N_manifolds))
      this%stoner_param = 0.0_dp
      return
    case default
      call system_abort("Self_Consistency_Term_setup_system confused by this%type="//this%type)
  end select

end subroutine Self_Consistency_Term_setup_system

subroutine TB_Dipole_Model_Print(this,file)
  type(TB_Dipole_Model),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  integer :: i_type, i_orb_set

  if (current_verbosity() < PRINT_NORMAL) return

  call Print('TB_Dipole_Model : active='//this%active, file=file)
  call Print('TB_Dipole_Model : n_types='//this%n_types, file=file)
  do i_type = 1, this%n_types
    call print('TB_Dipole_Model : type ' // i_type // ' atomic_num ' // this%atomic_num(i_type) // &
	       ' n_orb_sets ' // this%n_orb_sets(i_type), file=file)
    do i_orb_set = 1, this%n_orb_sets(i_type)
      call print('TB_Dipole_model:     i_orb_set ' // i_orb_set // ' orb_type ' // this%orb_set_type(i_orb_set,i_type) // &
		 ' width ' // this%gaussian_width(i_orb_set, i_type) // ' phase ' // this%orb_set_phase(i_orb_set, i_type), file=file)
    end do
  end do

end subroutine TB_Dipole_Model_Print

subroutine TB_Spin_Orbit_Coupling_Print(this,file)
  type(TB_Spin_Orbit_Coupling),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  integer :: i_type, i_orb_set

  if (current_verbosity() < PRINT_NORMAL) return

  call Print('TB_Spin_Orbit_Coupling : active='//this%active, file=file)
  call Print('TB_Spin_Orbit_Coupling : n_types='//this%n_types, file=file)
  do i_type = 1, this%n_types
    call print('TB_Spin_Orbit_Coupling : type ' // i_type // ' atomic_num ' // this%atomic_num(i_type) // &
	       ' n_orb_sets ' // this%n_orb_sets(i_type), file=file)
    do i_orb_set = 1, this%n_orb_sets(i_type)
      call print('TB_Spin_Orbit_Coupling:     i_orb_set ' // i_orb_set // ' orb_type ' // this%orb_set_type(i_orb_set,i_type) // &
		 ' SO_param ' // this%SO_param(i_orb_set, i_type), file=file)
    end do
  end do

end subroutine TB_Spin_Orbit_Coupling_Print

subroutine Self_Consistency_Print(this,file)
  type(Self_Consistency),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  integer::i

  if (current_verbosity() < PRINT_NORMAL) return

  call Print('Self_Consistency : active='//this%active, file=file)

  call Print("Self_Consistency : global_U " // this%global_U, file=file)

  call Print("Self_Consistency : allocated(this%U) " // allocated(this%U), file=file)
  if (allocated(this%U)) then
     call Print("Self_Consistency : size(this%U) " // size(this%U), file=file)
     do i=1, size(this%U)
        if (this%U(i) /= 0.0_dp) then
           call Print('Self_Consistency : U ' // i // " " // this%U(i), file=file)
        endif
     end do
  else
     call Print("Self_Consistency : size(this%U) 0", file=file)
  endif
  if (allocated(this%stoner_param)) then
     call Print("Self_Consistency : size(this%stoner_param) " // size(this%stoner_param), file=file)
     do i=1, size(this%stoner_param,2)
        if (maxval(this%stoner_param(:,i)) /= 0.0_dp) then
           call Print('Self_Consistency : stoner_param ' // i // " " // this%stoner_param(:,i), file=file)
        endif
     end do
  else
     call Print("Self_Consistency : size(this%stoner_param) 0", file=file)
  endif

  if (allocated(this%terms)) then
     do i=1, size(this%terms)
        call Print('Self_Consistency : term ' // i, file=file)
        call Print(this%terms(i),file=file)
     end do
  end if

  call verbosity_push_decrement()

  if (allocated(this%orb_local_pot)) then
    call Print('Self_Consistency : orb_local_pot ', file=file)
    call Print(this%orb_local_pot, file=file)
  endif

  if (allocated(this%orb_exch_field)) then
    call Print('Self_Consistency : orb_exch_field ', file=file)
    call Print(this%orb_exch_field, file=file)
  endif

  call verbosity_pop()

end subroutine Self_Consistency_Print

subroutine Self_Consistency_Term_Print(this, file)
  type(Self_Consistency_Term), intent(in) :: this
  type(Inoutput), intent(inout), optional :: file

  integer :: i

  if (current_verbosity() < PRINT_NORMAL) return

  if (.not. this%active) return

  call Print('Self_Consistency_Term type : ' // scf_names(this%type), file=file)

  call Print('Self_Consistency_Term n_dof : ' // this%n_dof, file=file)

  if (this%type == SCF_GLOBAL_U) &
    call Print('Self_Consistency_Term global_U ' // this%global_U // ' global_N ' // this%global_N, file=file)

  if (this%type == SCF_GCN) &
    call Print('Self_Consistency_Term global_pot ' // this%global_pot // ' global_N ' // this%global_N, file=file)

  if (allocated(this%atomic_n)) then
    call Print('Self_Consistency_Term : atomic_n ', file=file)
    call Print(this%atomic_n, file=file)
  endif

  if (allocated(this%manifold_mom)) then
    call Print('Self_Consistency_Term : manifold_mom ', file=file)
    do i=1, size(this%manifold_mom,2)
      call Print(i // " " // this%manifold_mom(:,i) // " mag " // norm(this%manifold_mom(:,i)), file=file)
    end do
  endif

  call verbosity_push_decrement()

  if (allocated(this%U)) then
    call Print('Self_Consistency_Term : U ', file=file)
    call Print(this%U, file=file)
  endif

  if (allocated(this%spin_splitting)) then
    call Print('Self_Consistency_Term : spin_splitting ', file=file)
    call Print(this%spin_splitting, file=file)
  endif

  if (allocated(this%stoner_param)) then
    call Print('Self_Consistency_Term : stoner_param ', file=file)
    call Print(this%stoner_param, file=file)
  endif

  if (allocated(this%atomic_n0)) then
    call Print('Self_Consistency_Term : atomic_n0 ', file=file)
    call Print(this%atomic_n0, file=file)
  endif

  if (allocated(this%atomic_local_pot)) then
    call Print('Self_Consistency_Term : atomic_local_pot ', file=file)
    call Print(this%atomic_local_pot, file=file)
  endif

  call verbosity_push_decrement()

  if (allocated(this%gamma)) then
    call Print('Self_Consistency : gamma ', file=file)
    call Print(this%gamma, file=file)
  endif
  if (allocated(this%dgamma_dr)) then
    call Print('Self_Consistency : dgamma_dr(:,:,1) ', file=file)
    call Print(this%dgamma_dr(:,:,1), file=file)
    call Print('Self_Consistency : dgamma_dr(:,:,2) ', file=file)
    call Print(this%dgamma_dr(:,:,2), file=file)
    call Print('Self_Consistency : dgamma_dr(:,:,3) ', file=file)
    call Print(this%dgamma_dr(:,:,3), file=file)
  endif

  call verbosity_pop()
  call verbosity_pop()
end subroutine Self_Consistency_Term_Print

subroutine TBSystem_fill_sc_matrices(this, at)
  type(TBSystem), intent(inout) :: this
  type(Atoms), intent(in) :: at

  integer i

  if (allocated(this%scf%terms)) then
     do i=1, size(this%scf%terms)
        call fill_sc_matrices(this%scf%terms(i), at, this%scf%U, this%scf%stoner_param, this%first_orb_of_atom, &
             this%first_manifold_of_atom, this%first_orb_of_manifold, this%tbmodel)
     end do
  end if
end subroutine TBSystem_fill_sc_matrices

subroutine Self_Consistency_Term_fill_sc_matrices(this, at, U, stoner_param, first_orb_of_atom, first_manifold_of_atom, first_orb_of_manifold, tbm)
  type(Self_Consistency_Term), intent(inout) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(in), allocatable :: U(:), stoner_param(:,:)
  integer, intent(in) :: first_orb_of_atom(:), first_manifold_of_atom(:), first_orb_of_manifold(:)
  type(TBModel), intent(in) :: tbm

  real(dp), allocatable :: block_H_up(:,:), block_H_down(:,:), block_S(:,:)
  integer :: i_at, i_man, n_i, orb_offset, man_offset

  select case(this%type)
    case (SCF_NONE)
      return
    case (SCF_LCN)
      return
    case (SCF_GCN)
      return
    case (SCF_GLOBAL_U)
      return
    case (SCF_LOCAL_U)
      if (.not. allocated(U)) &
	call system_abort("self_consistency_term_fill_sc_matrices called with SCF_LOCAL_U but U not allocated")
      if (minval(at%Z) < 1 .or. maxval(at%Z) > size(U)) &
	call system_abort("self_consistency_term_fill_scf_matrices SCF_LOCAL_U called with at%Z element out of range 1.."//size(U))
      if (minval(U(at%Z)) <= 0) &
	call system_abort("Tried to do a local U for atom " // minloc(U(at%Z)) // " Z " // at%Z(minloc(U(at%Z))) // " which has nonsense U")
      this%U = U(at%Z)
    case (SCF_NONLOCAL_U_DFTB)
      if (.not. allocated(U)) &
	call system_abort("self_consistency_term_fill_sc_matrices called with SCF_NONLOCAL_U_DFTB but U not allocated")
      if (minval(at%Z) < 1 .or. maxval(at%Z) > size(U)) &
	call system_abort("self_consistency_term_fill_scf_matrices SCF_NONLOCAL_U_DFTB called with at%Z element out of range 1.." // size(U))
      call calc_gamma_dftb(at, U, this%gamma)
    case (SCF_NONLOCAL_U_NRL_TB)
      if (.not. allocated(U)) &
	call system_abort("self_consistency_term_fill_sc_matrices called with SCF_NONLOCAL_U_NRL_TB but U not allocated")
      if (minval(at%Z) < 1 .or. maxval(at%Z) > size(U)) &
	call system_abort("self_consistency_term_fill_scf_matrices SCF_NONLOCAL_U_NRL_TB called with at%Z element out of range 1.." // size(U))
      call calc_gamma_nrl_tb(at, U, this%gamma)
    case (SCF_SPIN_DIR)
      do i_at=1, at%N
	n_i = first_orb_of_atom(i_at+1)-first_orb_of_atom(i_at)
	allocate(block_H_up(n_i,n_i))
	allocate(block_H_down(n_i,n_i))
	allocate(block_S(n_i,n_i))
	call get_HS_blocks(tbm, at, i_at, i_at, (/ 0.0_dp, 0.0_dp, 0.0_dp /), 0.0_dp, block_H_up, block_S, i_mag=1)
	call get_HS_blocks(tbm, at, i_at, i_at, (/ 0.0_dp, 0.0_dp, 0.0_dp /), 0.0_dp, block_H_down, block_S, i_mag=2)
	do i_man=first_manifold_of_atom(i_at), first_manifold_of_atom(i_at+1)-1
	  orb_offset = first_orb_of_manifold(i_man)-first_orb_of_atom(i_at)+1
	  this%spin_splitting(i_man) = -0.5_dp*(block_H_up(orb_offset,orb_offset)-block_H_down(orb_offset,orb_offset))
	end do
	deallocate(block_H_up)
	deallocate(block_H_down)
	deallocate(block_S)
      end do
    case (SCF_SPIN_STONER)
      if (.not. allocated(stoner_param)) &
	call system_abort("self_consistency_term_fill_sc_matrices called with SCF_SPIN_STONER but stoner_param not allocated")
      if (minval(at%Z) < 1 .or. maxval(at%Z) > size(stoner_param)) &
	call system_abort("self_consistency_term_fill_scf_matrices SCF_SPIN_STONER called with at%Z element out of range 1.."//size(stoner_param))
      do i_at=1, at%N
	do i_man=first_manifold_of_atom(i_at), first_manifold_of_atom(i_at+1)-1
	  man_offset = i_man - first_manifold_of_atom(i_at) + 1
	  this%stoner_param(i_man) = stoner_param(man_offset, at%Z(i_at))
	  this%stoner_param(i_man) = stoner_param(man_offset, at%Z(i_at))
	end do
      end do
    case default
      call system_abort("Self_Consistency_Term_fill_sc_matrices confused by this%type="//this%type)
  end select
end subroutine Self_Consistency_Term_fill_sc_matrices

subroutine realloc_dgamma_dr(this)
  type(Self_Consistency_Term), intent(inout) :: this

  if (.not. allocated(this%gamma)) then
    call system_abort ("Called realloc_dgamma_dr with gamma not allocated")
  endif

  if (allocated(this%dgamma_dr)) then
    if (size(this%dgamma_dr,1) /= size(this%gamma,1) .or. &
	size(this%dgamma_dr,2) /= size(this%gamma,2)) then
      deallocate(this%dgamma_dr)
    endif
  endif

  if (.not.allocated(this%dgamma_dr)) then
    allocate(this%dgamma_dr(size(this%gamma,1),size(this%gamma,2),3))
  endif
end subroutine realloc_dgamma_dr

subroutine TBSystem_fill_sc_dmatrices(this, at, virial_component)
  type(TBSystem), intent(inout) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in), optional ::  virial_component

  integer i_term

  if (allocated(this%scf%terms)) then
     do i_term=1, size(this%scf%terms)
        call fill_sc_dmatrices(this%scf%terms(i_term), at, this%scf%U, virial_component)
     end do
  end if
end subroutine TBSystem_fill_sc_dmatrices

subroutine Self_Consistency_Term_fill_sc_dmatrices(this, at, U, virial_component)
  type(Self_Consistency_Term), intent(inout) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(in) :: U(:)
  integer, intent(in), optional ::  virial_component

  select case(this%type)
    case (SCF_NONE)
      return
    case (SCF_LCN)
      return
    case (SCF_GCN)
      return
    case (SCF_GLOBAL_U)
      return
    case (SCF_LOCAL_U)
      return
    case (SCF_NONLOCAL_U_DFTB)
      call realloc_dgamma_dr(this)
      call calc_dgamma_dr_dftb(at, U, this%dgamma_dr, virial_component)
    case (SCF_NONLOCAL_U_NRL_TB)
      call realloc_dgamma_dr(this)
      call calc_dgamma_dr_nrl_tb(at, U, this%dgamma_dr, virial_component)
    case (SCF_SPIN_DIR)
      call system_abort("fill_sc_dmatrices: no SCF_SPIN_DIR yet")
      return
    case (SCF_SPIN_STONER)
      return
    case default
      call system_abort("Self_Consistency_Term_fill_sc_dmatrices Confused by this%type="//this%type)
  end select
end subroutine Self_Consistency_Term_fill_sc_dmatrices

! from a (input) scf system and a new (output) orbital_n vector, compute new
! control parameters (density, local potential) for next SCF step
function TBSystem_update_orb_local_pot(this, at, iter, global_at_weight, new_orbital_n, new_orbital_m)
  type(TBSystem), intent(inout) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: iter
  real(dp), intent(in), pointer :: global_at_weight(:)
  real(dp), intent(in), pointer :: new_orbital_n(:), new_orbital_m(:,:)
  logical :: TBSystem_update_orb_local_pot

  integer :: broyden_iter
  integer :: i_term, n_dof, cur_dof, last_dof
  real(dp) :: resid
  logical :: done

  real(dp), allocatable :: to_zero_vec(:), control_vec(:), new_control_vec(:)

  if (.not. this%scf%active .or. size(this%scf%terms) == 0) then
      TBSystem_update_orb_local_pot = .true.
      return
  endif

  n_dof = 0
  if (allocated(this%scf%terms)) then
     do i_term=1, size(this%scf%terms)
        n_dof = n_dof + this%scf%terms(i_term)%n_dof
     end do
  end if

  allocate(to_zero_vec(n_dof), control_vec(n_dof), new_control_vec(n_dof))

  cur_dof = 1
  if (allocated(this%scf%terms)) then
     do i_term=1, size(this%scf%terms)
        last_dof = cur_dof + this%scf%terms(i_term)%n_dof - 1
        call get_vecs(this%scf%terms(i_term), this, to_zero_vec(cur_dof:last_dof), control_vec(cur_dof:last_dof), global_at_weight, new_orbital_n, new_orbital_m)
        cur_dof = last_dof+1
     end do
  end if

  if (this%scf%mix_simple) then
    call do_mix_simple(iter, control_vec, to_zero_vec, new_control_vec, this%scf%alpha)
  else
    if (n_dof == 1) then
      call do_ridders_residual(iter, control_vec(1), to_zero_vec(1), new_control_vec(1))
    else
      broyden_iter = mod(iter-1,100)+1 ! throw out Broyden history every 100 steps
      call do_mix_broyden(broyden_iter, control_vec, to_zero_vec, new_control_vec, this%scf%alpha, this%scf%w0)
    endif
  endif

  resid = maxval(abs(to_zero_vec))
  call Print("SCF iteration " // iter // " residual " // resid)

  done = (resid < this%scf%conv_tol)

  if (.not. done .and. this%scf%mix_simple .and. resid < this%scf%mix_simple_end_tol) then
    call print("mix_simple is active, resid < this%scf%mix_simple_end_tol, switching simple mixing off", PRINT_ALWAYS)
    this%scf%mix_simple = .false.
  endif
  if (.not. done .and. .not. this%scf%mix_simple .and. resid > this%scf%mix_simple_start_tol) then
    call print("mix_simple is not active, resid > this%scf%mix_simple_start_tol, switching simple mixing on", PRINT_ALWAYS)
    this%scf%mix_simple = .true.
  endif

  cur_dof = 1
  if (allocated(this%scf%terms)) then
     do i_term=1, size(this%scf%terms)
        last_dof = cur_dof + this%scf%terms(i_term)%n_dof - 1
        ! SCF_SPIN_DIR isn't self consistent, so always update 
        if ((.not. done) .or. this%scf%terms(i_term)%type == SCF_SPIN_DIR) then 
           call set_vec(this%scf%terms(i_term), new_control_vec(cur_dof:last_dof))
        end if
        cur_dof = last_dof+1
     end do
  end if

  call calc_orb_local_pot(this, global_at_weight)

  deallocate(to_zero_vec, control_vec, new_control_vec)

  TBSystem_update_orb_local_pot = done

end function TBSystem_update_orb_local_pot

subroutine get_vecs(this, tbsys, to_zero_vec, control_vec, global_at_weight, new_orbital_n, new_orbital_m)
  type(Self_Consistency_Term), intent(inout) :: this
  type(TBSystem), intent(in) :: tbsys
  real(dp), intent(out) :: to_zero_vec(:), control_vec(:)
  real(dp), intent(in), pointer :: new_orbital_n(:), new_orbital_m(:,:)
  real(dp), intent(in), pointer :: global_at_weight(:)

  integer :: i, N_mom
  real(dp), allocatable :: new_manifold_m(:,:)

  if (this%type == SCF_NONE) return

  select case(this%type)
    case (SCF_NONE)
      return
    case (SCF_LCN)
      if (.not. associated(new_orbital_n)) call system_abort("get_vecs needs new_orbital_n for SCF_LCN")
      to_zero_vec = atom_orbital_sum(tbsys, new_orbital_n) - this%atomic_n0
      control_vec = this%atomic_local_pot
      if (current_verbosity() >= PRINT_NERD) then
	call print("get_vecs SCF_LCN delta_n" // to_zero_vec, PRINT_NERD)
      endif
    case(SCF_GCN)
      if (.not. associated(new_orbital_n)) call system_abort("get_vecs needs new_orbital_n for SCF_GCN")
      if (.not. associated(global_at_weight)) &
	call system_abort("Self_Consistency_Term_get_vecs SCF_GCN but no global_at_weight")
      to_zero_vec = sum(global_at_weight*(atom_orbital_sum(tbsys, new_orbital_n)-this%atomic_n0))
      control_vec = this%global_pot
      if (current_verbosity() >= PRINT_NERD) then
	call print("get_vecs SCF_GCN delta_n" // to_zero_vec, PRINT_NERD)
      endif
    case(SCF_GLOBAL_U)
      if (.not. associated(new_orbital_n)) call system_abort("get_vecs needs new_orbital_n for SCF_GLOBAL_U")
      if (.not. associated(global_at_weight)) &
	call system_abort("Self_Consistency_Term_get_vecs SCF_GCN but no global_at_weight")
      to_zero_vec = sum(global_at_weight*(atom_orbital_sum(tbsys, new_orbital_n)))-this%global_N
      control_vec = this%global_N
      if (current_verbosity() >= PRINT_NERD) then
	call print("get_vecs SCF_GLOBAL_U delta_n" // to_zero_vec, PRINT_NERD)
      endif
    case (SCF_LOCAL_U)
      if (.not. associated(new_orbital_n)) call system_abort("get_vecs needs new_orbital_n for SCF_LOCAL_U")
      to_zero_vec = atom_orbital_sum(tbsys, new_orbital_n) - this%atomic_n
      control_vec = this%atomic_n
      if (current_verbosity() >= PRINT_NERD) then
	call print("get_vecs SCF_LOCAL_U delta_n" // to_zero_vec, PRINT_NERD)
      endif
    case (SCF_NONLOCAL_U_DFTB)
      if (.not. associated(new_orbital_n)) call system_abort("get_vecs needs new_orbital_n for SCF_NONLOCAL_U_DFTB")
      to_zero_vec = atom_orbital_sum(tbsys, new_orbital_n) - this%atomic_n
      control_vec = this%atomic_n
      if (current_verbosity() >= PRINT_NERD) then
	call print("get_vecs SCF_NONLOCAL_U_DFTB delta_n" // to_zero_vec, PRINT_NERD)
      endif
    case (SCF_NONLOCAL_U_NRL_TB)
      if (.not. associated(new_orbital_n)) call system_abort("get_vecs needs new_orbital_n for SCF_NONLOCAL_U_NRL_TB")
      to_zero_vec = atom_orbital_sum(tbsys, new_orbital_n) - this%atomic_n
      control_vec = this%atomic_n
      if (current_verbosity() >= PRINT_NERD) then
	call print("get_vecs SCF_NONLOCAL_U_NRL_TB delta_n" // to_zero_vec, PRINT_NERD)
      endif
    case (SCF_SPIN_DIR)
      if (.not. associated(new_orbital_m)) call system_abort("get_vecs needs new_orbital_m for SCF_SPIN_DIR")
      allocate(new_manifold_m(3,tbsys%N_manifolds))
      N_mom = tbsys%N_manifolds
      new_manifold_m = manifold_orbital_sum(tbsys, new_orbital_m)
      to_zero_vec(1:3*N_mom:3) = new_manifold_m(1,1:N_mom) - this%manifold_mom(1,1:N_mom)
      to_zero_vec(2:3*N_mom:3) = new_manifold_m(2,1:N_mom) - this%manifold_mom(2,1:N_mom)
      to_zero_vec(3:3*N_mom:3) = new_manifold_m(3,1:N_mom) - this%manifold_mom(3,1:N_mom)
      control_vec(1:3*N_mom:3) = this%manifold_mom(1,1:N_mom)
      control_vec(2:3*N_mom:3) = this%manifold_mom(2,1:N_mom)
      control_vec(3:3*N_mom:3) = this%manifold_mom(3,1:N_mom)
      if (current_verbosity() >= PRINT_NERD) then
	call print("get_vecs SCF_SPIN_DIR current manifold moment", PRINT_NERD)
	do i=1, tbsys%N_manifolds
	  call print("  " // i // " " // new_manifold_m(:,i), PRINT_NERD)
	end do
      endif
      deallocate(new_manifold_m)
    case (SCF_SPIN_STONER)
      if (.not. associated(new_orbital_m)) call system_abort("get_vecs needs new_orbital_m for SCF_SPIN_STONER")
      allocate(new_manifold_m(3,tbsys%N_manifolds))
      N_mom = tbsys%N_manifolds
      new_manifold_m = manifold_orbital_sum(tbsys, new_orbital_m)
      if (current_verbosity() >= PRINT_NERD) then
	call print("get_vecs SCF_SPIN_STONER current manifold moment", PRINT_NERD)
	do i=1, tbsys%N_manifolds
	  call print("  " // i // " " // new_manifold_m(:,i), PRINT_NERD)
	end do
      endif
      to_zero_vec(1:3*N_mom:3) = new_manifold_m(1,1:N_mom) - this%manifold_mom(1,1:N_mom)
      to_zero_vec(2:3*N_mom:3) = new_manifold_m(2,1:N_mom) - this%manifold_mom(2,1:N_mom)
      to_zero_vec(3:3*N_mom:3) = new_manifold_m(3,1:N_mom) - this%manifold_mom(3,1:N_mom)
      control_vec(1:3*N_mom:3) = this%manifold_mom(1,1:N_mom)
      control_vec(2:3*N_mom:3) = this%manifold_mom(2,1:N_mom)
      control_vec(3:3*N_mom:3) = this%manifold_mom(3,1:N_mom)
      deallocate(new_manifold_m)
    case default
      call system_abort("Self_Consistency_Term_get_vecs confused by this%type="//this%type)
  end select
end subroutine get_vecs

subroutine set_vec(this, new_control_vec)
  type(Self_Consistency_Term), intent(inout) :: this
  real(dp), intent(in) :: new_control_vec(:)

  integer :: N_mom

  select case(this%type)
    case (SCF_NONE)
      return
    case (SCF_LCN)
      this%atomic_local_pot = new_control_vec
    case(SCF_GCN)
      this%global_pot = new_control_vec(1)
    case (SCF_GLOBAL_U)
      this%global_N = new_control_vec(1)
    case (SCF_LOCAL_U)
      this%atomic_n = new_control_vec
    case (SCF_NONLOCAL_U_DFTB)
      this%atomic_n = new_control_vec
    case (SCF_NONLOCAL_U_NRL_TB)
      this%atomic_n = new_control_vec
    case (SCF_SPIN_DIR)
      N_mom = size(this%manifold_mom,2)
      this%manifold_mom(1,1:N_mom) = new_control_vec(1:3*N_mom:3)
      this%manifold_mom(2,1:N_mom) = new_control_vec(2:3*N_mom:3)
      this%manifold_mom(3,1:N_mom) = new_control_vec(3:3*N_mom:3)
    case (SCF_SPIN_STONER)
      N_mom = size(this%manifold_mom,2)
      this%manifold_mom(1,1:N_mom) = new_control_vec(1:3*N_mom:3)
      this%manifold_mom(2,1:N_mom) = new_control_vec(2:3*N_mom:3)
      this%manifold_mom(3,1:N_mom) = new_control_vec(3:3*N_mom:3)
    case default
      call system_abort("Self_Consistency_Term_set_vec confused by this%type="//this%type)
  end select
end subroutine set_vec

! calculate orbital local potential from terms aready set in SCF object (local pot
! for LCN and GCN, and global_N or atomic_n for U based SCF.
subroutine TBSystem_calc_orb_local_pot(this, global_at_weight)
  type(TBSystem), intent(inout) :: this
  real(dp), intent(in), pointer :: global_at_weight(:)

  integer :: i_term

  this%scf%orb_local_pot = 0.0_dp

  if (.not. this%scf%active) return

  if (allocated(this%scf%terms)) then
     do i_term=1, size(this%scf%terms)
	if (allocated(this%scf%orb_local_pot)) &
	  call add_term_dSCFE_dn(this%scf%terms(i_term), this, global_at_weight, this%scf%orb_local_pot)
	if (allocated(this%scf%orb_exch_field)) &
	  call add_term_dSCFE_dm(this%scf%terms(i_term), this, global_at_weight, this%scf%orb_exch_field)
     end do
  end if

end subroutine TBSystem_calc_orb_local_pot

subroutine add_term_dSCFE_dn(this, tbsys, global_at_weight, dSCFE_dn)
  type(Self_Consistency_Term), intent(inout) :: this
  type(TBSystem), intent(in) :: tbsys
  real(dp), intent(in), pointer :: global_at_weight(:)
  real(dp), intent(out) :: dSCFE_dn(:)

  dSCFE_dn = 0.0_dp

  select case (this%type)
    case (SCF_NONE)
      return
    case (SCF_LCN)
      dSCFE_dn = atom_orbital_spread(tbsys, this%atomic_local_pot)
    case (SCF_GCN)
      if (.not. associated(global_at_weight)) &
	call system_abort("add_term_dSCFE_dn  with SCF_GCN but no global_at_weight")
      dSCFE_dn = atom_orbital_spread(tbsys, global_at_weight*this%global_pot)
    case (SCF_GLOBAL_U)
      if (.not. associated(global_at_weight)) &
	call system_abort("add_term_dSCFE_dn with SCF_GLOBAL_U but no global_at_weight")
      dSCFE_dn = atom_orbital_spread(tbsys, global_at_weight*this%global_U*(this%global_N-sum(global_at_weight*this%atomic_n0)))
    case (SCF_LOCAL_U)
      dSCFE_dn = atom_orbital_spread(tbsys, this%U*(this%atomic_n-this%atomic_n0))
    case (SCF_NONLOCAL_U_DFTB)
      dSCFE_dn = atom_orbital_spread(tbsys, matmul(this%gamma,(this%atomic_n-this%atomic_n0)))
    case (SCF_NONLOCAL_U_NRL_TB)
      dSCFE_dn = atom_orbital_spread(tbsys, matmul(this%gamma,(this%atomic_n-this%atomic_n0)))
    case (SCF_SPIN_DIR)
      return
    case (SCF_SPIN_STONER)
      return
    case default
      call system_abort("add_term_dSCFE_dn confused by this%type="//this%type)
  end select
end subroutine add_term_dSCFE_dn

subroutine add_term_dSCFE_dm(this, tbsys, global_at_weight, dSCFE_dm)
  type(Self_Consistency_Term), intent(inout) :: this
  type(TBSystem), intent(in) :: tbsys
  real(dp), intent(in), pointer :: global_at_weight(:)
  real(dp), intent(out) :: dSCFE_dm(:,:)

  integer :: i_man, i_orb
  real(dp) :: e(3)

  dSCFE_dm = 0.0_dp

  select case (this%type)
    case (SCF_NONE)
      return
    case (SCF_LCN)
      return
    case (SCF_GCN)
      return
    case (SCF_GLOBAL_U)
      return
    case (SCF_LOCAL_U)
      return
    case (SCF_NONLOCAL_U_DFTB)
      return
    case (SCF_NONLOCAL_U_NRL_TB)
      return
    case (SCF_SPIN_DIR)
      do i_man=1, tbsys%N_manifolds
	e = this%manifold_mom(:,i_man)
	if (norm(e) .feq. 0.0_dp) then
	  e = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
	else
	  e = e / norm(e)
	endif
	do i_orb = tbsys%first_orb_of_manifold(i_man), tbsys%first_orb_of_manifold(i_man+1)-1
	  ! F = -0.5 sum_a spin_split_a |m_a|
	  ! dF/dm_i = -0.5 spin_split_a (|m_a| 2*m_a - |m_a|^2 m_a/|m_a|) / |m_a|^2
	  !         = -0.5 spin_split_a ( 2*m_a/|m_a| - m_a/|m_a| )
	  !         = -0.5 spin_split_a m_a/|m_a|
	  dSCFE_dm(1:3,i_orb) = -0.5_dp * this%spin_splitting(i_man)*e(1:3)
	end do
      end do
      return
    case (SCF_SPIN_STONER)
      do i_man=1, tbsys%N_manifolds
	do i_orb=tbsys%first_orb_of_manifold(i_man), tbsys%first_orb_of_manifold(i_man+1)-1
	  ! F = -0.25 sum_a stoner_param_a * |m_a|^2
	  ! dF/dm_i = -0.25 stoner_param_a 2 |m_a| m_a/|m_a|
	  !         = -0.5 stoner_param_a m_a
	  dSCFE_dm(1:3,i_orb) = -0.5_dp * this%stoner_param(i_man)*this%manifold_mom(1:3,i_man)
	end do
      end do
      return
    case default
      call system_abort("add_term_dSCFE_dn confused by this%type="//this%type)
  end select

end subroutine add_term_dSCFE_dm

subroutine add_term_d2SCFE_dgNdn(this, tbsys, global_at_weight, d2SCFE_dgNdn)
  type(Self_Consistency_Term), intent(inout) :: this
  type(TBSystem), intent(in) :: tbsys
  real(dp), intent(in):: global_at_weight(:)
  real(dp), intent(out) :: d2SCFE_dgNdn(:)

  d2SCFE_dgNdn = 0.0_dp

  if (.not. this%active) return

  select case (this%type)
    case (SCF_NONE)
      return
    case (SCF_GLOBAL_U)
      d2SCFE_dgNdn = atom_orbital_spread(tbsys, global_at_weight*this%global_U)
    case default
      call system_abort("add_term_d2SCFE_dgNdn only defined for GLOBAL_U")
  end select
end subroutine add_term_d2SCFE_dgNdn

subroutine add_term_d2SCFE_dn2_times_vec(this, tbsys, vec, d2SCFE_dn2_times_vec)
  type(Self_Consistency_Term), intent(inout) :: this
  type(TBSystem), intent(in) :: tbsys
  real(dp), intent(in) :: vec(:)
  real(dp), intent(out) :: d2SCFE_dn2_times_vec(:)

  d2SCFE_dn2_times_vec = 0.0_dp

  if (.not. this%active) return

  select case (this%type)
    case (SCF_NONE)
      return
    case (SCF_LOCAL_U)
      d2SCFE_dn2_times_vec = atom_orbital_spread(tbsys, this%U*vec)
    case (SCF_NONLOCAL_U_DFTB)
      d2SCFE_dn2_times_vec = atom_orbital_spread(tbsys, matmul(this%gamma,vec))
    case (SCF_NONLOCAL_U_NRL_TB)
      d2SCFE_dn2_times_vec = atom_orbital_spread(tbsys, matmul(this%gamma,vec))
    case default
      call system_abort("add_term_d2SCFE_dgNdn only defined for LOCAL_U and NONLOCAL_U_*")
  end select
end subroutine add_term_d2SCFE_dn2_times_vec

subroutine calc_gamma_dftb(at, U, gamma)
  type(Atoms), intent(in) :: at
  real(dp), intent(in) :: U(:)
  real(dp), intent(out) :: gamma(:,:)

  integer :: i, j
  real(dp) :: U_i, U_j
  real(dp) :: val

  real(dp) :: max_val
  integer s1, s2, s3, s_range
  real(dp) :: pi(3), pj(3), o(3), dr(3)

  gamma = 0.0_dp

  do i=1, at%N
    if (U(at%Z(i)) <= 0) call system_abort("Tried to do a local U for atom " // i // " Z " // at%Z(i) // " which has nonsense U")
    U_i = U(at%Z(i))
    gamma(i,i) = U_i
  end do

  max_val = 1.0_dp
  s_range = 1
  do while (max_val > 1.0e-12_dp)
    max_val = 0.0_dp
    do s1 = -s_range, s_range
    do s2 = -s_range, s_range
    do s3 = -s_range, s_range
      if ((s_range > 1) .and. &
	  (s1 /= -s_range) .and. (s1 /= s_range) .and. &
	  (s2 /= -s_range) .and. (s2 /= s_range) .and. &
	  (s3 /= -s_range) .and. (s3 /= s_range)) cycle
      do i=1, at%N
	U_i = U(at%Z(i))
	pi = at%pos(1:3,i)
	do j=i, at%N
	  if (i == j .and. s1 == 0 .and. s2 == 0 .and. s3 == 0) cycle
	  U_j = U(at%Z(j))
	  pj = at%pos(1:3,j)
	  o = s1*at%lattice(:,1) + s2*at%lattice(:,2) + s3*at%lattice(:,3)
	  dr = (pi - pj - o)
	  val = Hartree*dftb_s(U_i/Hartree, U_j/Hartree, sqrt(sum(dr*dr))/Bohr)
	  if (abs(val) > max_val) max_val = abs(val)
	  gamma(i,j) = gamma(i,j) - val
	  gamma(j,i) = gamma(j,i) - val
	end do
      end do
    end do
    end do
    end do
    s_range = s_range + 1
  end do

  call add_madelung_matrix(at%N, at%lattice(:,1), at%lattice(:,2), at%lattice(:,3), at%pos, gamma)

end subroutine calc_gamma_dftb

subroutine calc_dgamma_dr_dftb(at, U, dgamma_dr, virial_component)
  type(Atoms), intent(in) :: at
  real(dp), intent(in) :: U(:)
  real(dp), intent(out) :: dgamma_dr(:,:,:)
  integer, intent(in), optional :: virial_component

  integer :: i, j
  real(dp) :: U_i, U_j
  real(dp) :: dval

  real(dp) :: max_val
  integer s1, s2, s3, s_range
  real(dp) :: pi(3), pj(3), o(3), dr(3), dr_mag

  integer my_virial_component

  my_virial_component = 0
  if (present(virial_component)) my_virial_component = virial_component

  dgamma_dr = 0.0_dp

  max_val = 1.0_dp
  s_range = 1
  do while (max_val > 1.0e-12_dp)
    max_val = 0.0_dp
    do s1=-s_range,s_range
    do s2=-s_range,s_range
    do s3=-s_range,s_range
      if ((s_range > 1) .and. &
	  (s1 /= -s_range) .and. (s1 /= s_range) .and. &
	  (s2 /= -s_range) .and. (s2 /= s_range) .and. &
	  (s3 /= -s_range) .and. (s3 /= s_range)) cycle
      do i=1, at%N
	if (U(at%Z(i)) <= 0) call system_abort("Tried to do a local U for atom " // i // " Z " // at%Z(i) // " which has nonsense U")
	U_i = U(at%Z(i))
	pi = at%pos(1:3,i)
	do j=i+1, at%N
	  if (U(at%Z(j)) <= 0) call system_abort("Tried to do a local U for atom " // j // " Z " // at%Z(j) // " which has nonsense U")
	  U_j = U(at%Z(j))
	  pj = at%pos(1:3,j)
	  o = s1*at%lattice(:,1) + s2*at%lattice(:,2) + s3*at%lattice(:,3)
	  dr = (pi - pj - o)
	  dr_mag = sqrt(sum(dr*dr))
	  dr = dr/dr_mag
	  dval = Hartree/Bohr * dftb_s_deriv(U_i/Hartree, U_j/Hartree, dr_mag/Bohr)
	  if (abs(dval) > max_val) max_val = abs(dval)
	  if (my_virial_component > 0) then
	    dgamma_dr(i,j,:) = dgamma_dr(i,j,:) + dval*dr*dr_mag*dr(my_virial_component)
	    dgamma_dr(j,i,:) = dgamma_dr(j,i,:) + dval*dr*dr_mag*dr(my_virial_component)
	  else
	    dgamma_dr(i,j,:) = dgamma_dr(i,j,:) - dval*dr
	    dgamma_dr(j,i,:) = dgamma_dr(j,i,:) + dval*dr
	  endif
	end do
      end do
    end do
    end do
    end do
    s_range = s_range + 1
  end do

  if (my_virial_component > 0) then
    call add_dmadelung_matrix_dr(at%N, at%lattice(:,1), at%lattice(:,2), at%lattice(:,3), at%pos, my_virial_component, dgamma_dr)
  else
    call add_dmadelung_matrix(at%N, at%lattice(:,1), at%lattice(:,2), at%lattice(:,3), at%pos, dgamma_dr)
  endif

end subroutine calc_dgamma_dr_dftb

subroutine calc_gamma_nrl_tb(at, U, gamma)
  type(Atoms), intent(in) :: at
  real(dp), intent(in) :: U(:)
  real(dp), intent(out) :: gamma(:,:)

  integer :: i, j
  real(dp) :: U_i, U_j
  real(dp) :: val

  real(dp) :: max_val
  integer s1, s2, s3, s_range
  real(dp) :: pi(3), pj(3), o(3), dr(3)

  gamma = 0.0_dp

  do i=1, at%N
    if (at%Z(i) < 1 .or. at%Z(i) > size(U)) &
      call system_abort("calc_gamma_nrl_tb got at%Z("//i//") out of range 1..size(U)="//size(U))
    if (U(at%Z(i)) <= 0) &
      call system_abort("calc_gamma_nrl_tb tried to do a local U for atom " // i // " Z " // at%Z(i) // " which has nonsense U")
    U_i = U(at%Z(i))
    gamma(i,i) = U_i
  end do

  max_val = 1.0_dp
  s_range = 1
  do while (max_val > 1.0e-10_dp)
    max_val = 0.0_dp
    do s1=-s_range,s_range
    do s2=-s_range,s_range
    do s3=-s_range,s_range
      if ((s_range > 1) .and. &
	  (s1 /= -s_range) .and. (s1 /= s_range) .and. &
	  (s2 /= -s_range) .and. (s2 /= s_range) .and. &
	  (s3 /= -s_range) .and. (s3 /= s_range)) cycle
      do i=1, at%N
	U_i = U(at%Z(i))
	pi = at%pos(1:3,i)
	do j=i, at%N
	  if (i == j .and. s1 == 0 .and. s2 == 0 .and. s3 == 0) cycle
	  U_j = U(at%Z(j))
	  pj = at%pos(1:3,j)
	  o = s1*at%lattice(:,1) + s2*at%lattice(:,2) + s3*at%lattice(:,3)
	  dr = (pi - pj - o)
	  val = Hartree*nrl_tb_s(U_i/Hartree, U_j/Hartree, sqrt(sum(dr*dr))/Bohr)
	  if (abs(val) > max_val) max_val = abs(val)
	  gamma(i,j) = gamma(i,j) - val
	  gamma(j,i) = gamma(j,i) - val
	end do
      end do
    end do
    end do
    end do
    s_range = s_range + 1
  end do

  call add_madelung_matrix(at%N, at%lattice(:,1), at%lattice(:,2), at%lattice(:,3), at%pos, gamma)

end subroutine calc_gamma_nrl_tb

subroutine calc_dgamma_dr_nrl_tb(at, U, dgamma_dr, virial_component)
  type(Atoms), intent(in) :: at
  real(dp), intent(in) :: U(:)
  real(dp), intent(out) :: dgamma_dr(:,:,:)
  integer, intent(in), optional :: virial_component

  integer :: i, j
  real(dp) :: U_i, U_j
  real(dp) :: dval

  real(dp) :: max_val
  integer s1, s2, s3, s_range
  real(dp) :: pi(3), pj(3), o(3), dr(3), dr_mag
  integer my_virial_component

  my_virial_component = 0
  if (present(virial_component)) my_virial_component = virial_component

  dgamma_dr = 0.0_dp

  max_val = 1.0_dp
  s_range = 1
  do while (max_val > 1.0e-10_dp)
    max_val = 0.0_dp
    do s1=-s_range,s_range
    do s2=-s_range,s_range
    do s3=-s_range,s_range
      if ((s_range > 1) .and. &
	  (s1 /= -s_range) .and. (s1 /= s_range) .and. &
	  (s2 /= -s_range) .and. (s2 /= s_range) .and. &
	  (s3 /= -s_range) .and. (s3 /= s_range)) cycle
      do i=1, at%N
	if (U(at%Z(i)) <= 0) call system_abort("Tried to do a local U for atom " // i // " Z " // at%Z(i) // " which has nonsense U")
	U_i = U(at%Z(i))
	pi = at%pos(1:3,i)
	do j=i+1, at%N
	  if (U(at%Z(j)) <= 0) call system_abort("Tried to do a local U for atom " // j // " Z " // at%Z(j) // " which has nonsense U")
	  U_j = U(at%Z(j))
	  pj = at%pos(1:3,j)
	  o = s1*at%lattice(:,1) + s2*at%lattice(:,2) + s3*at%lattice(:,3)
	  dr = (pi - pj - o)
	  dr_mag = sqrt(sum(dr*dr))
	  dr = dr/dr_mag
	  dval = Hartree/Bohr * nrl_tb_s_deriv(U_i/Hartree, U_j/Hartree, dr_mag/Bohr)
	  if (abs(dval) > max_val) max_val = abs(dval)
	  if (my_virial_component > 0) then
	    dgamma_dr(i,j,:) = dgamma_dr(i,j,:) + dval*dr*dr_mag*dr(my_virial_component)
	    dgamma_dr(j,i,:) = dgamma_dr(j,i,:) + dval*dr*dr_mag*dr(my_virial_component)
	  else
	    dgamma_dr(i,j,:) = dgamma_dr(i,j,:) - dval*dr
	    dgamma_dr(j,i,:) = dgamma_dr(j,i,:) + dval*dr
	  endif
	end do
      end do
    end do
    end do
    end do
    s_range = s_range + 1
  end do

  if (my_virial_component > 0) then
    call add_dmadelung_matrix_dr(at%N, at%lattice(:,1), at%lattice(:,2), at%lattice(:,3), at%pos, my_virial_component, dgamma_dr)
  else
    call add_dmadelung_matrix(at%N, at%lattice(:,1), at%lattice(:,2), at%lattice(:,3), at%pos, dgamma_dr)
  endif

end subroutine calc_dgamma_dr_nrl_tb

! From Michael Strenberg's Ph.D. Thesis, p. 27, Paderborn 2001
function dftb_s(u_a, u_b, R)
  real(dp), intent(in) :: u_a, u_b, R
  real(dp) :: dftb_s

  real(dp) tau_a, tau_b

  tau_a = 16.0_dp/5.0_dp * u_a
  tau_b = 16.0_dp/5.0_dp * u_b

  if (tau_a /= tau_b) then
    dftb_s = &
      exp(-tau_a*R)* ( (tau_b**4 * tau_a) / (2.0_dp * (tau_a**2 - tau_b**2)**2) -  &
		       (tau_b**6 - 3.0_dp*tau_b**4*tau_a**2) / (R * (tau_a**2-tau_b**2)**3) ) + &
      exp(-tau_b*R)* ( (tau_a**4 * tau_b) / (2.0_dp * (tau_b**2 - tau_a**2)**2) -  &
		       (tau_a**6 - 3.0_dp*tau_a**4*tau_b**2) / (R * (tau_b**2-tau_a**2)**3) )
  else
    dftb_s = & 
      exp(-tau_a*R) * (48.0_dp + 33.0_dp * tau_a * R + 9.0_dp * (tau_a*R)**2 + (tau_a*R)**3) / &
		      (48.0_dp * R)
  endif
end function dftb_s

function dftb_s_deriv(u_a, u_b, R)
  real(dp), intent(in) :: u_a, u_b, R
  real(dp) :: dftb_s_deriv

  real(dp) tau_a, tau_b
  real(dp) exp_v, rat_v, exp_v_deriv, rat_v_deriv

  tau_a = 16.0_dp/5.0_dp * u_a
  tau_b = 16.0_dp/5.0_dp * u_b

  if (tau_a /= tau_b) then
    exp_v = exp(-tau_a*R)
    exp_v_deriv = -tau_a*exp_v
    rat_v = (tau_b**4 * tau_a) / (2.0_dp * (tau_a**2 - tau_b**2)**2) -  &
	    (tau_b**6 - 3.0_dp*tau_b**4*tau_a**2) / (R * (tau_a**2-tau_b**2)**3)
    rat_v_deriv = (tau_b**6 - 3.0_dp*tau_b**4*tau_a**2) / (R**2 * (tau_a**2-tau_b**2)**3)
    dftb_s_deriv = exp_v*rat_v_deriv + exp_v_deriv*rat_v

    exp_v = exp(-tau_b*R)
    exp_v_deriv = -tau_b*exp_v
    rat_v = (tau_a**4 * tau_b) / (2.0_dp * (tau_b**2 - tau_a**2)**2) -  &
	    (tau_a**6 - 3.0_dp*tau_a**4*tau_b**2) / (R * (tau_b**2-tau_a**2)**3)
    rat_v_deriv = (tau_a**6 - 3.0_dp*tau_a**4*tau_b**2) / (R**2 * (tau_b**2-tau_a**2)**3)

    dftb_s_deriv = dftb_s_deriv + exp_v*rat_v_deriv + exp_v_deriv*rat_v

  else
    exp_v = exp(-tau_a*R)
    exp_v_deriv = -tau_a*exp_v
    rat_v = (48.0_dp + 33.0_dp * tau_a * R + 9.0_dp * (tau_a*R)**2 + (tau_a*R)**3) / (48.0_dp * R)
    rat_v_deriv = ( (48.0_dp*R) * (33.0_dp * tau_a + 18.0_dp * tau_a**2 * R + 3*(tau_a*R)**2 * tau_a) - &
		    (48.0_dp + 33.0_dp * tau_a * R + 9.0_dp * (tau_a*R)**2 + (tau_a*R)**3) * 48.0_dp ) / &
		    ((48.0_dp*R)**2)

    dftb_s_deriv = exp_v*rat_v_deriv + exp_v_deriv*rat_v
  endif
end function dftb_s_deriv

function nrl_tb_s_deriv(u_a, u_b, R)
    real(dp), intent(in) :: u_a, u_b, R
    real(dp) :: nrl_tb_s_deriv

    real(dp) a, c

    a = PI/2.0D0 * u_a**2
    c = PI/2.0D0 * u_b**2
    
    nrl_tb_s_deriv = -1.0D0/R**2 - ( (a*c)**1.5D0 / PI**3 ) * &
        2.0D0*PI**2.5D0 / (a*c*sqrt(a+c)) * dF0(a,c,R)
end function nrl_tb_s_deriv

function nrl_tb_s(u_a, u_b, R)
  real(dp), intent(in) :: u_a, u_b, R
  real(dp) :: nrl_tb_s

  real(dp) :: a, c

  a = PI/2.0_dp * u_a**2
  c = PI/2.0_dp * u_b**2

  nrl_tb_s = 1.0_dp/R - ( (a*c)**1.5_dp / PI**3 ) * &
      2.0_dp*PI**2.5_dp / (a*c*sqrt(a+c)) * F0(a,c,R)
end function nrl_tb_s

function F0(a, c, R)
  real(dp), intent(in) :: a, c, R
  real(dp) :: F0

  double precision coeff, q

  ! gamma(a,b,c,d) = 1
  ! from G&R 3.321 + 8.250
  coeff = a*c/(a+c) * R**2
  q = sqrt(coeff)

  F0 = sqrt(PI)/(2*q) * erf(q)

end function F0

function dF0(a,c,R)
    real(dp), intent(in) :: a, c, R
    real(dp) :: dF0

    real(dp) coeff, q

    ! gamma(a,b,c,d) = 1
    ! from G&R 3.321 + 8.250
    coeff = a*c/(a+c) * R**2
    q = sqrt(coeff)

    dF0 = ( exp(-q*q)/q - sqrt(PI)/(2*q*q) * erf(q) ) * sqrt((a*c)/(a+c))

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine local_scf_e_correction(this, at, local_e, global_E, global_at_weight)
  type(TBSystem), intent(in) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(out) :: local_e(:), global_E
  real(dp), pointer, intent(in), optional :: global_at_weight(:)

  integer :: i_term

  local_e = 0.0_dp
  global_e = 0.0_dp

  if (.not. this%scf%active) return

  if (allocated(this%scf%terms)) then
     do i_term=1, size(this%scf%terms)
        call add_term_local_scf_e_correction(this%scf%terms(i_term), local_e, global_E, this%N_atoms, this%first_manifold_of_atom, global_at_weight)
     end do
  end if

end subroutine local_scf_e_correction

! for SCF energy, + SCFE - dSCFE_dn * n
! for charge neutrality, -pot*n0
subroutine add_term_local_scf_e_correction(this, local_e, global_E, N_atoms, first_manifold_of_atom, global_at_weight)
  type(Self_Consistency_Term), intent(in) :: this
  real(dp), intent(out) :: local_e(:), global_E
  integer :: N_atoms, first_manifold_of_atom(:)
  real(dp), pointer, intent(in), optional :: global_at_weight(:)

  integer :: i_at, i_man

  select case(this%type)
    case (SCF_NONE)
      return
    case (SCF_LCN)
      local_e = local_e - this%atomic_local_pot*this%atomic_n0
    case (SCF_GCN)
      if (.not. present(global_at_weight)) &
        call system_abort ("local_scf_e_correction with SCF_GCN but without global_at_weight")
      if (.not. associated(global_at_weight)) &
        call system_abort ("local_scf_e_correction with SCF_GCN but with unassociated global_at_weight")
      global_e = global_e - this%global_pot*sum(global_at_weight*this%atomic_n0)
    case (SCF_GLOBAL_U)
      if (.not. present(global_at_weight)) &
        call system_abort ("local_scf_e_correction with SCF_GLOBAL_U but without global_at_weight")
      if (.not. associated(global_at_weight)) &
        call system_abort ("local_scf_e_correction with SCF_GLOBAL_U but with unassociated global_at_weight")
      global_e = global_e + 0.5_dp*this%global_U*(this%global_N-sum(global_at_weight*this%atomic_n0))**2 - &
			    this%global_U*(this%global_N-sum(global_at_weight*this%atomic_n0))*this%global_N
    case (SCF_LOCAL_U)
      local_e = local_e + 0.5_dp*this%U*(this%atomic_n-this%atomic_n0)**2 - this%U*(this%atomic_n-this%atomic_n0)*this%atomic_n
    case (SCF_NONLOCAL_U_DFTB)
      local_e = local_e + 0.5_dp*((this%atomic_n-this%atomic_n0) * matmul(this%gamma, this%atomic_n-this%atomic_n0)) - &
			  (this%atomic_n * matmul(this%gamma,this%atomic_n-this%atomic_n0))
    case (SCF_NONLOCAL_U_NRL_TB)
      local_e = local_e + 0.5_dp*((this%atomic_n-this%atomic_n0) * matmul(this%gamma, this%atomic_n-this%atomic_n0)) - &
			  (this%atomic_n * matmul(this%gamma,this%atomic_n-this%atomic_n0))
    case (SCF_SPIN_DIR)
      return
    case (SCF_SPIN_STONER)
      do i_at=1, N_atoms
	do i_man=first_manifold_of_atom(i_at), first_manifold_of_atom(i_at+1)-1
	  ! add F({m_i})
	  local_e(i_at) = local_e(i_at) - 0.25_dp*this%stoner_param(i_man)*norm(this%manifold_mom(1:3,i_man))**2
	  ! subtract dF/dm_i . m_i
	  local_e(i_at) = local_e(i_at) + 0.5_dp*this%stoner_param(i_man)*norm(this%manifold_mom(1:3,i_man))**2
	end do
      end do
    case default
      call system_abort("add_term_local_scf_e_correction confused by this%type="//this%type)
  end select

end subroutine add_term_local_scf_e_correction

subroutine add_term_dscf_e_correction_dgN(this, dglobal_E_dgN)
  type(Self_Consistency_Term), intent(in) :: this
  real(dp), intent(out) :: dglobal_E_dgN

  dglobal_E_dgN = 0.0_dp

  if (.not. this%active) return

  select case(this%type)
    case (SCF_NONE)
      return
    case (SCF_GLOBAL_U)
      dglobal_E_dgN = this%global_U*this%global_N
    case default
      call system_abort("add_term_dscf_e_correction_dgN only defined for GLOBAL_U")
  end select
end subroutine

subroutine add_term_dscf_e_correction_dn(this, tbsys, dlocal_E_dn)
  type(Self_Consistency_Term), intent(in) :: this
  type(TBSystem), intent(in) :: tbsys
  real(dp), intent(out) :: dlocal_E_dn(:)

  dlocal_E_dn = 0.0_dp

  if (.not. this%active) return

  select case(this%type)
    case (SCF_NONE)
      return
    case (SCF_LOCAL_U)
      dlocal_E_dn = atom_orbital_spread(tbsys, this%U*this%atomic_n)
    case (SCF_NONLOCAL_U_DFTB)
      dlocal_E_dn = atom_orbital_spread(tbsys, matmul(this%gamma,this%atomic_n))
    case (SCF_NONLOCAL_U_NRL_TB)
      dlocal_E_dn = atom_orbital_spread(tbsys, matmul(this%gamma,this%atomic_n))
    case default
      call system_abort("add_term_dscf_e_correction_dgN only defined for GLOBAL_U")
  end select
end subroutine

subroutine realloc_datomic_local_pot_dr(this)
  type(Self_Consistency_Term), intent(inout) :: this

  if (.not. allocated(this%atomic_n)) then
    call system_abort ("realloc_datomic_local_pot_dr with local_pot not allocated")
  endif

  if (allocated(this%datomic_local_pot_dr)) then
    if (size(this%datomic_local_pot_dr,1) /= size(this%atomic_n,1)) then
      deallocate(this%datomic_local_pot_dr)
    endif
  endif

  if (.not.allocated(this%datomic_local_pot_dr)) then
    allocate(this%datomic_local_pot_dr(size(this%atomic_n,1),3))
  endif
end subroutine realloc_datomic_local_pot_dr

function scf_f_correction(this, at) result(forces_scf)
  type(TBSystem), intent(inout) :: this
  type(Atoms), intent(in) :: at
  real(dp) :: forces_scf(3,at%N) ! result

  integer i_term

  forces_scf = 0.0

  call fill_sc_dmatrices(this, at)

  if (allocated(this%scf%terms)) then
     do i_term=1, size(this%scf%terms)
        forces_scf = forces_scf + add_term_scf_f_correction(this%scf%terms(i_term), at)
     end do
  end if

end function scf_f_correction

function add_term_scf_f_correction(this, at) result (forces_scf)
  type(Self_Consistency_Term), intent(inout) :: this
  type(Atoms), intent(in) :: at
  real(dp) :: forces_scf(3,at%N) ! result

  forces_scf = 0.0_dp

  select case(this%type)
    case (SCF_NONE)
      return
    case (SCF_GCN)
      return
    case (SCF_GLOBAL_U)
      return
    case (SCF_LCN)
      return
    case (SCF_LOCAL_U)
      return
    case (SCF_NONLOCAL_U_DFTB)
      call realloc_datomic_local_pot_dr(this)
      this%datomic_local_pot_dr(:,1) = matmul(this%dgamma_dr(:,:,1),(this%atomic_n-this%atomic_n0))
      this%datomic_local_pot_dr(:,2) = matmul(this%dgamma_dr(:,:,2),(this%atomic_n-this%atomic_n0))
      this%datomic_local_pot_dr(:,3) = matmul(this%dgamma_dr(:,:,3),(this%atomic_n-this%atomic_n0))
      forces_scf(1,:) = -this%datomic_local_pot_dr(:,1)*(this%atomic_n-this%atomic_n0)
      forces_scf(2,:) = -this%datomic_local_pot_dr(:,2)*(this%atomic_n-this%atomic_n0)
      forces_scf(3,:) = -this%datomic_local_pot_dr(:,3)*(this%atomic_n-this%atomic_n0)
    case (SCF_NONLOCAL_U_NRL_TB)
      call realloc_datomic_local_pot_dr(this)
      this%datomic_local_pot_dr(:,1) = matmul(this%dgamma_dr(:,:,1),(this%atomic_n-this%atomic_n0))
      this%datomic_local_pot_dr(:,2) = matmul(this%dgamma_dr(:,:,2),(this%atomic_n-this%atomic_n0))
      this%datomic_local_pot_dr(:,3) = matmul(this%dgamma_dr(:,:,3),(this%atomic_n-this%atomic_n0))
      forces_scf(1,:) = -this%datomic_local_pot_dr(:,1)*(this%atomic_n-this%atomic_n0)
      forces_scf(2,:) = -this%datomic_local_pot_dr(:,2)*(this%atomic_n-this%atomic_n0)
      forces_scf(3,:) = -this%datomic_local_pot_dr(:,3)*(this%atomic_n-this%atomic_n0)
    case (SCF_SPIN_DIR)
      call system_abort("add_term_scf_f_correction not implemented for SCF_SPIN_DIR yet")
    case (SCF_SPIN_STONER)
      return
    case default
      call system_abort("add_term_scf_f_correction confused by this%type="//this%type)
  end select
end function add_term_scf_f_correction

function scf_virial_correction(this, at) result(virial_scf)
  type(TBSystem), intent(inout) :: this
  type(Atoms), intent(in) :: at
  real(dp) :: virial_scf(3,3) ! result

  integer :: component, i_term

  virial_scf = 0.0_dp
  do component=1,3
    call fill_sc_dmatrices(this, at, component)
    if (allocated(this%scf%terms)) then
       do i_term=1, size(this%scf%terms)
          virial_scf(:,component) = virial_scf(:,component) + add_term_scf_virial_correction(this%scf%terms(i_term), at)
       end do
    end if
  end do
end function scf_virial_correction

function add_term_scf_virial_correction(this, at) result (virial_scf)
  type(Self_Consistency_Term), intent(inout) :: this
  type(Atoms), intent(in) :: at
  real(dp) :: virial_scf(3) ! result

  virial_scf = 0.0_dp

  select case(this%type)
    case (SCF_NONE)
      return
    case (SCF_GCN)
      return
    case (SCF_GLOBAL_U)
      return
    case (SCF_LCN)
      return
    case (SCF_LOCAL_U)
      return
    case (SCF_NONLOCAL_U_DFTB)
      call realloc_datomic_local_pot_dr(this)
      this%datomic_local_pot_dr(:,1) = matmul(this%dgamma_dr(:,:,1),(this%atomic_n-this%atomic_n0))
      this%datomic_local_pot_dr(:,2) = matmul(this%dgamma_dr(:,:,2),(this%atomic_n-this%atomic_n0))
      this%datomic_local_pot_dr(:,3) = matmul(this%dgamma_dr(:,:,3),(this%atomic_n-this%atomic_n0))
      virial_scf(1) = 0.5_dp*sum(this%datomic_local_pot_dr(:,1)*(this%atomic_n-this%atomic_n0))
      virial_scf(2) = 0.5_dp*sum(this%datomic_local_pot_dr(:,2)*(this%atomic_n-this%atomic_n0))
      virial_scf(3) = 0.5_dp*sum(this%datomic_local_pot_dr(:,3)*(this%atomic_n-this%atomic_n0))
    case (SCF_NONLOCAL_U_NRL_TB)
      call realloc_datomic_local_pot_dr(this)
      this%datomic_local_pot_dr(:,1) = matmul(this%dgamma_dr(:,:,1),(this%atomic_n-this%atomic_n0))
      this%datomic_local_pot_dr(:,2) = matmul(this%dgamma_dr(:,:,2),(this%atomic_n-this%atomic_n0))
      this%datomic_local_pot_dr(:,3) = matmul(this%dgamma_dr(:,:,3),(this%atomic_n-this%atomic_n0))
      virial_scf(1) = 0.5_dp*sum(this%datomic_local_pot_dr(:,1)*(this%atomic_n-this%atomic_n0))
      virial_scf(2) = 0.5_dp*sum(this%datomic_local_pot_dr(:,2)*(this%atomic_n-this%atomic_n0))
      virial_scf(3) = 0.5_dp*sum(this%datomic_local_pot_dr(:,3)*(this%atomic_n-this%atomic_n0))
    case (SCF_SPIN_DIR)
      call system_abort("add_term_scf_virial_correction not implemented for SCF_SPIN_DIR yet")
    case (SCF_SPIN_STONER)
      return
    case default
      call system_abort("calc_datomic_local_pot_dr confused by this%type="//this%type)
  end select
end function add_term_scf_virial_correction

function scf_e_correction(this, at, global_at_weight)
  type(TBSystem), intent(in) :: this
  type(Atoms), intent(in) :: at
  real(dp), pointer, intent(in), optional :: global_at_weight(:)
  real(dp) :: scf_e_correction

  real(dp), allocatable :: local_e_scf(:)
  real(dp) :: global_e_scf

  allocate(local_e_scf(at%N))

  call local_scf_e_correction(this, at, local_e_scf, global_e_scf, global_at_weight)
  scf_e_correction = sum(local_e_scf) + global_e_scf
  deallocate(local_e_scf)

end function scf_e_correction

subroutine get_dipole_block(this, at, i, j, dv_hat, dv_mag, block_dipole)
  type(TB_Dipole_Model), intent(in) :: this
  type(Atoms), intent(in) :: at
  integer, intent(in) :: i, j
  real(dp), intent(in) :: dv_hat(3), dv_mag
  real(dp), intent(out) :: block_dipole(:,:,:)

  real(dp) :: pi(3), pj(3), sigi, sigj
  integer :: phasei, phasej
  integer :: i_type, j_type, i_orb_set, j_orb_set, i_offset, j_offset
  integer :: first_i, last_i, first_j, last_j
  real(dp) :: W(10,10,4), W_T(9,9,4), trans_mat_i(10,9), trans_mat_j(10,9)

  pi = at%pos(1:3,i)
  pj = at%pos(1:3,i) + dv_hat(1:3)*dv_mag

  i_type = this%type_of_atomic_num(at%Z(i))
  j_type = this%type_of_atomic_num(at%Z(j))

  i_offset = 1
  do i_orb_set=1, this%n_orb_sets(i_type)
  sigi = this%gaussian_width(i_orb_set, i_type)
  phasei = this%orb_set_phase(i_orb_set, i_type)

    call init_transform(trans_mat_i, sigi, phasei)

    j_offset = 1
    do j_orb_set=1, this%n_orb_sets(j_type)
    sigj = this%gaussian_width(j_orb_set, j_type)
    phasej = this%orb_set_phase(j_orb_set, j_type)

      call init_transform(trans_mat_j, sigj, phasej)
      ! Mark Pederson's dipole matrix element routine
      call GINTED(sigi, sigj, pi, pj, W)
      W_T(1:9,1:9,1) = matmul(matmul(transpose(trans_mat_i),W(1:10,1:10,1)),trans_mat_j)
      W_T(1:9,1:9,2) = matmul(matmul(transpose(trans_mat_i),W(1:10,1:10,2)),trans_mat_j)
      W_T(1:9,1:9,3) = matmul(matmul(transpose(trans_mat_i),W(1:10,1:10,3)),trans_mat_j)

      if (this%orb_set_type(i_orb_set,i_type) == ORB_S) then
	first_i = 1
	last_i = 1
      else if (this%orb_set_type(i_orb_set,i_type) == ORB_P) then
	first_i = 2
	last_i = 4
      else if (this%orb_set_type(i_orb_set,i_type) == ORB_D) then
	first_i = 5
	last_i = 9
      else
	call system_abort("get_dipole_block got unknown orb_set_type for i_orb_set " // this%orb_set_type(i_orb_set,i_type))
      endif
      if (this%orb_set_type(j_orb_set,j_type) == ORB_S) then
	first_j = 1
	last_j = 1
      else if (this%orb_set_type(j_orb_set,j_type) == ORB_P) then
	first_j = 2
	last_j = 4
      else if (this%orb_set_type(j_orb_set,j_type) == ORB_D) then
	first_j = 5
	last_j = 9
      else
	call system_abort(" get_dipole_block got unknown orb_set_type for j_orb_set " // this%orb_set_type(j_orb_set,j_type))
      endif
      block_dipole(i_offset:i_offset+N_ORBS_OF_SET(this%orb_set_type(i_orb_set,i_type))-1, &
                   j_offset:j_offset+N_ORBS_OF_SET(this%orb_set_type(j_orb_set,j_type))-1, 1:3) = W_T(first_i:last_i,first_j:last_j,1:3)
      j_offset = j_offset + N_ORBS_OF_SET(this%orb_set_type(j_orb_set,j_type))
    end do ! j_orb_set
    i_offset = i_offset + N_ORBS_OF_SET(this%orb_set_type(i_orb_set,i_type))
  end do ! i_orb_set

end subroutine get_dipole_block

subroutine init_transform(trans, sigma, phase)
implicit none
    real(dp), intent(out) :: trans(10,9)
    real(dp), intent(in) :: sigma
    integer, intent(in) :: phase

    integer O_S, O_PX, O_PY, O_PZ
    parameter (O_S = 1, O_PX = 2, O_PY = 3, O_PZ = 4)
    integer O_DXX, O_DYY, O_DZZ, O_DXY, O_DYZ, O_DXZ
    parameter (O_DXX = 5, O_DYY = 6, O_DZZ = 7, O_DXY = 8, O_DXZ = 9, O_DYZ = 10)

    integer OT_S, OT_PX, OT_PY, OT_PZ
    parameter (OT_S = 1, OT_PY = 2, OT_PZ = 3, OT_PX = 4)
    integer OT_DXY, OT_DYZ, OT_DZX, OT_DXXYY, OT_DZZRR
    parameter (OT_DXY = 5, OT_DYZ = 6, OT_DZZRR = 7, OT_DZX = 8, OT_DXXYY = 9)
    real(dp) :: inv_rt2, inv_rt6
    real(dp) :: rt2

    real(dp) :: base_norm_factor, rt_sigma

    inv_rt2 = 1.0_dp/sqrt(2.0_dp)
    inv_rt6 = 1.0_dp/sqrt(6.0_dp)
    rt2 = sqrt(2.0_dp)

    base_norm_factor = 1.0D0/(PI/(2.0_dp*sigma))**(3.0_dp/4.0_dp)
    rt_sigma = sqrt(sigma)

    trans = 0.0D0

    trans(O_S,OT_S) = base_norm_factor

    trans(O_PX,OT_PX) = -2.0_dp*base_norm_factor*rt_sigma
    trans(O_PY,OT_PY) = -2.0_dp*base_norm_factor*rt_sigma
    trans(O_PZ,OT_PZ) = -2.0_dp*base_norm_factor*rt_sigma

    trans(O_DXY,OT_DXY) = 4.0_dp*base_norm_factor*sigma
    trans(O_DYZ,OT_DYZ) = 4.0_dp*base_norm_factor*sigma
    trans(O_DXZ,OT_DZX) = 4.0_dp*base_norm_factor*sigma

    trans(O_DXX,OT_DXXYY) = inv_rt2*rt2*2.0_dp*base_norm_factor*sigma
    trans(O_DYY,OT_DXXYY) = -inv_rt2*rt2*2.0_dp*base_norm_factor*sigma

    trans(O_DXX,OT_DZZRR) = -inv_rt6*rt2*2.0_dp*base_norm_factor*sigma
    trans(O_DYY,OT_DZZRR) = -inv_rt6*rt2*2.0_dp*base_norm_factor*sigma
    trans(O_DZZ,OT_DZZRR) = 2.0D0*inv_rt6*rt2*2.0_dp*base_norm_factor*sigma

    trans = trans*phase

end subroutine init_transform

subroutine realloc_dipole_model(this, new_n_types, new_n_orb_sets)
  type(TB_Dipole_Model), intent(inout) :: this
  integer, intent(in) :: new_n_types, new_n_orb_sets

  integer :: max_n_orb_sets
  integer, allocatable :: t_atomic_num(:), t_n_orb_sets(:), t_orb_set_type(:,:), t_orb_set_phase(:,:)
  real(dp), allocatable :: t_gaussian_width(:,:)

  if (size(this%n_orb_sets) > 0) then
    max_n_orb_sets = maxval(this%n_orb_sets,1)
  else
    max_n_orb_sets = 0
  endif
  if (new_n_types > this%n_types .or. new_n_orb_sets > max_n_orb_sets) then
    ! allocate temporaries if necessary
    if (this%n_types > 0) then
      allocate(t_n_orb_sets(this%n_types))
      t_n_orb_sets = this%n_orb_sets(1:this%n_types)
      allocate(t_atomic_num(this%n_types))
      t_atomic_num = this%atomic_num(1:this%n_types)
      if (max_n_orb_sets > 0) then
	allocate(t_orb_set_type(max_n_orb_sets,this%n_types))
	allocate(t_gaussian_width(max_n_orb_sets,this%n_types))
	allocate(t_orb_set_phase(max_n_orb_sets,this%n_types))
	t_orb_set_type = this%orb_set_type(1:max_n_orb_sets,1:this%n_types)
	t_gaussian_width = this%gaussian_width(1:max_n_orb_sets,1:this%n_types)
	t_orb_set_phase = this%orb_set_phase(1:max_n_orb_sets,1:this%n_types)
      endif
    endif

    ! deallocate and reallocate real arrays
    if (allocated(this%atomic_num)) deallocate(this%atomic_num)
    if (allocated(this%n_orb_sets)) deallocate(this%n_orb_sets)
    if (allocated(this%orb_set_type)) deallocate(this%orb_set_type)
    if (allocated(this%gaussian_width)) deallocate(this%gaussian_width)
    if (allocated(this%orb_set_phase)) deallocate(this%orb_set_phase)
    allocate(this%atomic_num(new_n_types))
    allocate(this%n_orb_sets(new_n_types))
    allocate(this%orb_set_type(max(new_n_orb_sets,max_n_orb_sets), new_n_types))
    allocate(this%orb_set_phase(max(new_n_orb_sets,max_n_orb_sets), new_n_types))
    allocate(this%gaussian_width(max(new_n_orb_sets,max_n_orb_sets), new_n_types))
    this%atomic_num = 0
    this%n_orb_sets = 0
    this%orb_set_type = 0
    this%gaussian_width = 0.0_dp
    this%orb_set_phase = 0

    ! copy from temporaries if necessary
    if (this%n_types > 0) then
      this%atomic_num(1:this%n_types) = t_atomic_num(1:this%n_types)
      this%n_orb_sets(1:this%n_types) = t_n_orb_sets(1:this%n_types)
      if (max_n_orb_sets > 0) then
	this%orb_set_type(1:max_n_orb_sets,1:this%n_types) = t_orb_set_type(1:max_n_orb_sets, 1:this%n_types)
	this%gaussian_width(1:max_n_orb_sets,1:this%n_types) = t_gaussian_width(1:max_n_orb_sets, 1:this%n_types)
	this%orb_set_phase(1:max_n_orb_sets,1:this%n_types) = t_orb_set_phase(1:max_n_orb_sets, 1:this%n_types)
      endif
    endif
    this%n_types = max(this%n_types, new_n_types)

    if (allocated(t_atomic_num)) deallocate(t_atomic_num)
    if (allocated(t_n_orb_sets)) deallocate(t_n_orb_sets)
    if (allocated(t_orb_set_type)) deallocate(t_orb_set_type)
    if (allocated(t_gaussian_width)) deallocate(t_gaussian_width)
    if (allocated(t_orb_set_phase)) deallocate(t_orb_set_phase)
  endif
end subroutine realloc_dipole_model

subroutine realloc_spin_orbit_coupling(this, new_n_types, new_n_orb_sets)
  type(TB_Spin_Orbit_Coupling), intent(inout) :: this
  integer, intent(in) :: new_n_types, new_n_orb_sets

  integer :: max_n_orb_sets
  integer, allocatable :: t_atomic_num(:), t_n_orb_sets(:), t_orb_set_type(:,:)
  real(dp), allocatable :: t_SO_param(:,:)

  if (size(this%n_orb_sets) > 0) then
    max_n_orb_sets = maxval(this%n_orb_sets,1)
  else
    max_n_orb_sets = 0
  endif
  if (new_n_types > this%n_types .or. new_n_orb_sets > max_n_orb_sets) then
    ! allocate temporaries if necessary
    if (this%n_types > 0) then
      allocate(t_n_orb_sets(this%n_types))
      t_n_orb_sets = this%n_orb_sets(1:this%n_types)
      allocate(t_atomic_num(this%n_types))
      t_atomic_num = this%atomic_num(1:this%n_types)
      if (max_n_orb_sets > 0) then
	allocate(t_orb_set_type(max_n_orb_sets,this%n_types))
	allocate(t_SO_param(max_n_orb_sets,this%n_types))
	t_orb_set_type = this%orb_set_type(1:max_n_orb_sets,1:this%n_types)
	t_SO_param = this%SO_param(1:max_n_orb_sets,1:this%n_types)
      endif
    endif

    ! deallocate and reallocate real arrays
    if (allocated(this%atomic_num)) deallocate(this%atomic_num)
    if (allocated(this%n_orb_sets)) deallocate(this%n_orb_sets)
    if (allocated(this%orb_set_type)) deallocate(this%orb_set_type)
    if (allocated(this%SO_param)) deallocate(this%SO_param)
    allocate(this%atomic_num(new_n_types))
    allocate(this%n_orb_sets(new_n_types))
    allocate(this%orb_set_type(max(new_n_orb_sets,max_n_orb_sets), new_n_types))
    allocate(this%SO_param(max(new_n_orb_sets,max_n_orb_sets), new_n_types))
    this%atomic_num = 0
    this%n_orb_sets = 0
    this%orb_set_type = 0
    this%SO_param = 0.0_dp

    ! copy from temporaries if necessary
    if (this%n_types > 0) then
      this%atomic_num(1:this%n_types) = t_atomic_num(1:this%n_types)
      this%n_orb_sets(1:this%n_types) = t_n_orb_sets(1:this%n_types)
      if (max_n_orb_sets > 0) then
	this%orb_set_type(1:max_n_orb_sets,1:this%n_types) = t_orb_set_type(1:max_n_orb_sets, 1:this%n_types)
	this%SO_param(1:max_n_orb_sets,1:this%n_types) = t_SO_param(1:max_n_orb_sets, 1:this%n_types)
      endif
    endif
    this%n_types = max(this%n_types, new_n_types)

    if (allocated(t_atomic_num)) deallocate(t_atomic_num)
    if (allocated(t_n_orb_sets)) deallocate(t_n_orb_sets)
    if (allocated(t_orb_set_type)) deallocate(t_orb_set_type)
    if (allocated(t_SO_param)) deallocate(t_SO_param)
  endif
end subroutine realloc_spin_orbit_coupling

subroutine get_SO_block(this, tbm, Z, block_SO)
  type(TB_Spin_Orbit_Coupling), intent(in) :: this
  type(TBModel), intent(in) :: tbm
  integer, intent(in) :: Z
  complex(dp), intent(out) :: block_SO(:,:)

  complex(dp) :: V(2,2)
  integer :: i_orb_set, i_orb_dir, j_orb_dir, i_orb_type
  integer :: offset_base, i_offset, j_offset
  integer :: i_type

  i_type = this%type_of_atomic_num(Z)
  block_SO = 0.0_dp
  offset_base = 0
  do i_orb_set=1, n_orb_sets_of_Z(tbm, Z)
    i_orb_type = orb_type_of_orb_set_of_Z(tbm, Z, i_orb_set)
    do i_orb_dir=1, n_orbs_of_orb_set_of_Z(tbm, Z, i_orb_set)
    do j_orb_dir=1, n_orbs_of_orb_set_of_Z(tbm, Z, i_orb_set)
      i_offset = offset_base + i_orb_dir - 1
      j_offset = offset_base + j_orb_dir - 1
      V = spin_orbit_function(i_orb_type, i_orb_dir, j_orb_dir) 
      block_SO(2*i_offset+1:2*i_offset+2,2*j_offset+1:2*j_offset+2) = this%SO_param(i_orb_set, i_type) * V(1:2,1:2)
    end do
    end do
    offset_base = offset_base + n_orbs_of_orb_set_of_Z(tbm, Z, i_orb_set)
  end do

end subroutine get_SO_block

subroutine tbsystem_initialise_kpoints(this, args_str, param_str, mpi_obj, from_tbsystem_initialise)
  type(TBSystem), intent(inout) :: this
  character(len=*), intent(in), optional :: args_str
  character(len=*), intent(in), optional :: param_str
  type(MPI_context), intent(in), optional :: mpi_obj
  logical, intent(in), optional :: from_tbsystem_initialise

  type(Dictionary) :: params
  logical :: use_k_density, use_k_density_once, k_use_mp
  integer :: k_mesh(3)
  logical :: my_from_tbsystem_initialise
  real(dp) :: k_density

  my_from_tbsystem_initialise = optional_default(.false., from_tbsystem_initialise)

  if (present(args_str)) then
    call initialise(params)
    call param_register(params, "k_density", "-1.0", k_density, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "use_k_density", "F", use_k_density, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "use_k_density_once", "F", use_k_density_once, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "k_mesh", "-1 -1 -1", k_mesh, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "k_use_mp", ""//this%kpoints_use_mp, k_use_mp, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='TBSystem_Initialise_KPoints')) then
      call system_abort("tbsystem_initialise_kpoints failed to parse args_str='"//trim(args_str)//"'")
    endif
    call finalise(params)
  else
    use_k_density = .false.
    use_k_density_once = .false.
    k_mesh = (/ -1, -1, -1 /)
    k_use_mp = .true.
  endif

  if (use_k_density .and. use_k_density_once) then
    call system_abort("tbsystem_initialise_kpoints confused because user passed both use_k_density and use_k_density_once")
  endif

  if (.not. my_from_tbsystem_initialise .and. use_k_density_once) then
    call system_abort("tbsystem_initialise_kp[oints use_k_density_once can only be passed from tbsystem_initialise")
  endif

  if ((use_k_density .or. use_k_density_once) .and. any(k_mesh > 0)) then
    call system_abort("tbsystem_initialise_kpoints confused by having both use_k_density ("//use_k_density// &
                      ") or use_k_density_once("//use_k_density_once//") and signs of intent to use k_mesh ("//k_mesh//")")
  endif

  if (my_from_tbsystem_initialise) then
    if (use_k_density .or. use_k_density_once) then
       this%kpoints_generate_dynamically = .true.
       this%kpoints_generate_next_dynamically = .not. use_k_density_once
       this%kpoints_use_mp = k_use_mp
       if (k_density >= 0.0_dp) then
         this%kpoints_k_space_density = k_density
       else if (this%tbmodel%has_default_k_density) then
         this%kpoints_k_space_density = this%tbmodel%default_k_density
       else
         this%kpoints_k_space_density = -1.0_dp
       endif
    else if (any(k_mesh > 0)) then
      call initialise(this%kpoints, k_mesh, k_use_mp, mpi_obj)
      call Initialise_tbsystem_k_dep_stuff(this, mpi_obj)
    else
      if (present(param_str)) then
        call initialise(this%kpoints, param_str, mpi_obj)
        call Initialise_tbsystem_k_dep_stuff(this, mpi_obj)
      endif
    endif
  else ! not from tbsystem_initialise
    if (use_k_density) then
       this%kpoints_generate_dynamically = .true.
       this%kpoints_generate_next_dynamically = .true.
       this%kpoints_use_mp = k_use_mp
       if (k_density >= 0.0_dp) then
         this%kpoints_k_space_density = k_density
       else if (this%tbmodel%has_default_k_density) then
         this%kpoints_k_space_density = this%tbmodel%default_k_density
       else
         this%kpoints_k_space_density = -1.0_dp
       endif
    else if (any(k_mesh > 0)) then
      call initialise(this%kpoints, k_mesh, k_use_mp, mpi_obj)
      call Initialise_tbsystem_k_dep_stuff(this, mpi_obj)
    endif
  endif


end subroutine tbsystem_initialise_kpoints

end module TBSystem_module
