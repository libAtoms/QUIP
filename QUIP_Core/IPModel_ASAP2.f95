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
!X IPModel_ASAP2
!X
!% Reimplementation of ASAP potential:
!% P. Tangney and S. Scandolo,
!% An ab initio parametrized interatomic force field for silica
!% J. Chem. Phys, 117, 8898 (2002). 
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_ASAP2_module

use libatoms_module
use functions_module
use Yukawa_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_ASAP2
type IPModel_ASAP2
  integer :: n_types = 0
  real(dp) :: betapol, tolpol, yukalpha, yuksmoothlength
  integer :: maxipol, pred_order
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)
  real(dp), allocatable, dimension(:) :: pol, z
  real(dp), allocatable, dimension(:,:) :: D_ms, gamma_ms, R_ms, B_pol, C_pol, pseudise_sigma

  real(dp) :: cutoff_coulomb, cutoff_ms

  character(len=FIELD_LENGTH) :: label
  logical :: initialised, tdip_sr

end type IPModel_ASAP2

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_ASAP2), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_ASAP2_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_ASAP2_Finalise
end interface Finalise

interface Print
  module procedure IPModel_ASAP2_Print
end interface Print

public :: setup_atoms
interface setup_atoms
   module procedure IPModel_ASAP2_setup_atoms
end interface

interface Calc
  module procedure IPModel_ASAP2_Calc
end interface Calc

contains


subroutine IPModel_ASAP2_Initialise_str(this, args_str, param_str)
  type(IPModel_ASAP2), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  this%initialised = .false.

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_ASAP2_Initialise_str args_str')) then
    call system_abort("IPModel_ASAP2_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_ASAP2_read_params_xml(this, param_str)
  this%initialised = .true.
  
end subroutine IPModel_ASAP2_Initialise_str

subroutine IPModel_ASAP2_Finalise(this)
  type(IPModel_ASAP2), intent(inout) :: this

  this%initialised = .false.

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%pol)) deallocate(this%pol)
  if (allocated(this%z)) deallocate(this%z)
  if (allocated(this%D_ms)) deallocate(this%D_ms)
  if (allocated(this%gamma_ms)) deallocate(this%gamma_ms)
  if (allocated(this%R_ms)) deallocate(this%R_ms)
  if (allocated(this%B_pol)) deallocate(this%B_pol)
  if (allocated(this%C_pol)) deallocate(this%C_pol)
  if (allocated(this%pseudise_sigma)) deallocate(this%pseudise_sigma)
  this%n_types = 0
  this%label = ''
end subroutine IPModel_ASAP2_Finalise

subroutine asap_short_range_dipole_moments(this, at, charge, dip_sr, mpi)
  type(IPModel_ASAP2), intent(inout) :: this
  type(Atoms), intent(inout) :: at
  real(dp), dimension(:), intent(in) :: charge
  real(dp), dimension(:,:), intent(out) :: dip_sr
   type(MPI_Context), intent(in), optional :: mpi

  integer :: i, ti, m, j, tj, k
  real(dp) :: r_ij, u_ij(3), qj, bij, cij, dist3, dist5, gij, factork, expfactor, fc, dfc_dr
  integer, parameter :: nk = 4
  real(dp), allocatable :: private_dip_sr(:,:)

  call system_timer('asap_short_range_dipole_moments')

  dip_sr = 0.0_dp

  !$omp parallel default(none) shared(this, mpi, at, charge, dip_sr) private(ti, m, j, tj, k, r_ij, u_ij, qj, bij, cij, dist3, dist5, gij, factork, expfactor, fc, dfc_dr, private_dip_sr)

  allocate(private_dip_sr(size(dip_sr,1),size(dip_sr,2)))
  private_dip_sr = 0.0_dp

  !$omp do schedule(runtime)
  do i=1, at%n
     if (present(mpi)) then
	if (mpi%active) then
 	  if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
        endif
     endif

     ti = get_type(this%type_of_atomic_num, at%Z(i))
     if (.not. (abs(this%pol(ti)) > 0.0_dp)) cycle

     do m = 1, atoms_n_neighbours(at, i)
         
         j = atoms_neighbour(at, i, m, distance=r_ij, diff=u_ij, max_dist=(this%cutoff_ms*BOHR))
         if (j <= 0) cycle
         if (r_ij .feq. 0.0_dp) cycle

         r_ij = r_ij/BOHR
         u_ij = u_ij/BOHR

         tj = get_type(this%type_of_atomic_num, at%Z(j))
         qj = charge(j)
         bij = this%b_pol(ti, tj)
         cij = this%c_pol(ti, tj)

         dist3 = r_ij**3.0_dp
         dist5 = r_ij**5.0_dp
         
         gij = 0.0_dp
         factork = cij*exp(-bij*r_ij)
         do k=1,nk
            gij = gij + factork 
            factork = factork*bij*r_ij/real(k,dp)
         enddo
         gij = gij + factork
         
         expfactor = exp(-this%yukalpha*r_ij)
         call smooth_cutoff(r_ij, this%cutoff_coulomb-this%yuksmoothlength, this%yuksmoothlength, fc, dfc_dr)

         private_dip_sr(:,i) = private_dip_sr(:,i) - this%pol(ti)*qj*u_ij*gij/dist3*expfactor*fc
      end do
  end do
  !$omp end do
  
  if (present(mpi)) then
     if (mpi%active)   call sum_in_place(mpi, private_dip_sr)
  endif

  !$omp critical
  dip_sr = dip_sr + private_dip_sr
  !$omp end critical

  deallocate(private_dip_sr)
  !$omp end parallel

  call system_timer('asap_short_range_dipole_moments')

end subroutine asap_short_range_dipole_moments

!% Morse-stretch potential, defined by
!% \begin{displaymath}
!% U_ij = D_ij \left[ e^{\gamma_{ij}\left( 1 - r_{ij}/r^0_{ij} \right)} 
!%                  - 2 e^{\left( 1 - r_{ij}/r^0_{ij} \right)} \right]
!% \end{displaymath}
subroutine asap_morse_stretch(this, at, e, local_e, f, virial, mpi)
#ifdef _OPENMP
   use omp_lib
#endif
   type(IPModel_ASAP2), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:)
   real(dp), intent(out), optional :: virial(3,3)
   type(MPI_Context), intent(in), optional :: mpi

   integer i, j, m, ti, tj
   real(dp) :: r_ij, u_ij(3), dms, gammams, rms
   real(dp) :: exponentms, factorms, phi, de
   real(dp) :: dforce
   real(dp) :: elimitij(this%n_types, this%n_types)
   logical :: i_is_min_image, j_is_min_image

   real(dp) :: private_virial(3,3), private_e
   real(dp), allocatable :: private_f(:,:), private_local_e(:)

   call system_timer('asap_morse_stretch')

   ! Evaluate potential at cutoff. Will be subtracted from total energy.
   elimitij = 0.0_dp
   do ti=1,this%n_types
      do tj=1,this%n_types
         phi = exp(this%gamma_ms(ti,tj)*(1.0_dp - this%cutoff_ms/this%r_ms(ti,tj)))
         elimitij(ti,tj) = this%d_ms(ti,tj)*(phi - 2.0_dp*sqrt(phi))
      end do
   end do

   !$omp parallel default(none) shared(this, mpi, at, e, local_e, f, virial, elimitij) private(i, j, m, ti, tj, r_ij, u_ij, dms, gammams, rms, exponentms, factorms, phi, de, dforce, i_is_min_image, j_is_min_image, private_virial, private_e, private_f, private_local_e)

   if (present(e)) private_e = 0.0_dp
   if (present(local_e)) then
      allocate(private_local_e(at%N))
      private_local_e = 0.0_dp
   endif
   if (present(f)) then
      allocate(private_f(3,at%N))
      private_f = 0.0_dp
   endif
   if (present(virial)) private_virial = 0.0_dp
  
   !$omp do schedule(runtime)
   do i=1, at%n
      if (present(mpi)) then
	 if (mpi%active) then
	    if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
	 endif
      endif

      if (allocated(at%connect%is_min_image)) then
         i_is_min_image = at%connect%is_min_image(i)
      else
         i_is_min_image = is_min_image(at, i)
      end if
      ti = get_type(this%type_of_atomic_num, at%Z(i))
      do m = 1, atoms_n_neighbours(at, i)
         
         j = atoms_neighbour(at, i, m, distance=r_ij, cosines=u_ij, max_dist=(this%cutoff_ms*BOHR))
         
         if (j <= 0) cycle
         if (r_ij .feq. 0.0_dp) cycle

         if (allocated(at%connect%is_min_image)) then
            j_is_min_image = at%connect%is_min_image(j)
         else
            j_is_min_image = is_min_image(at, j)
         end if

         if (i < j .and. i_is_min_image .and. j_is_min_image) cycle

         r_ij = r_ij/BOHR

         tj = get_type(this%type_of_atomic_num, at%Z(j))

         dms = this%d_ms(ti,tj)
         gammams = this%gamma_ms(ti,tj)
         rms = this%r_ms(ti, tj)

         exponentms = gammams*(1.0_dp-r_ij/rms)
         factorms = gammams/rms
         phi = exp(exponentms)

         if (present(e) .or. present(local_e)) then
            de = dms*(phi-2.0_dp*sqrt(phi)) - elimitij(ti, tj)

            if (present(e)) then
               if (i_is_min_image .and. j_is_min_image) then
                  private_e = private_e + de
               else
                  private_e = private_e + 0.5_dp*de
               end if
            end if
            
            if (present(local_e)) then
               private_local_e(i) = private_local_e(i) + 0.5_dp*de
               if (i_is_min_image .and. j_is_min_image) private_local_e(j) = private_local_e(j) + 0.5_dp*de
            end if
         end if

         if (present(f) .or. present(virial)) then
            dforce  = -dms*(factorms*phi - factorms*dsqrt(phi))

            if (present(f)) then
               private_f(:,i) = private_f(:,i) + dforce*u_ij
               if (i_is_min_image .and. j_is_min_image) private_f(:,j) = private_f(:,j) - dforce*u_ij
            end if
            if (present(virial)) then
               if (i_is_min_image .and. j_is_min_image) then
                  private_virial = private_virial - dforce*(u_ij .outer. u_ij)*r_ij
               else
                  private_virial = private_virial - 0.5_dp*dforce*(u_ij .outer. u_ij)*r_ij
               end if
            end if
         end if
      end do
   end do

   if (present(mpi)) then
      if (mpi%active) then
	 if (present(e)) private_e = sum(mpi, private_e) 
	 if (present(local_e)) call sum_in_place(mpi, private_local_e)
	 if (present(f)) call sum_in_place(mpi, private_f)
	 if (present(virial)) call sum_in_place(mpi, private_virial)
      end if
   end if

   !$omp critical
   if (present(e)) e = e + private_e
   if (present(f)) f = f + private_f
   if (present(local_e)) local_e = local_e + private_local_e
   if (present(virial)) virial = virial + private_virial
   !$omp end critical 

   if (allocated(private_f)) deallocate(private_f)
   if (allocated(private_local_e)) deallocate(private_local_e)

   !$omp end parallel

   call system_timer('asap_morse_stretch')

end subroutine asap_morse_stretch

subroutine IPModel_ASAP2_setup_atoms(this, at)
  type(IPModel_ASAP2), intent(in):: this
  type(Atoms), intent(inout)      :: at

  logical :: dummy
  integer :: i, ti
  real(dp), dimension(:), pointer :: charge
  
  ! If charge property doesn't exist, we add and initialise it now
   
  if (.not. has_property(at, 'charge')) then
     call add_property(at, 'charge', 0.0_dp)
     dummy = assign_pointer(at, 'charge', charge)
     do i=1, at%N
        ti = get_type(this%type_of_atomic_num, at%Z(i))
        charge(i) = this%z(ti)
     end do
  end if
  
  if (.not. has_property(at, 'fixdip')) call add_property(at, 'fixdip', .false.)
  if (.not. has_property(at, 'efield')) call add_property(at, 'efield', 0.0_dp, n_cols=3)
  if (.not. has_property(at, 'dipoles')) call add_property(at, 'dipoles', 0.0_dp, n_cols=3)
  if (.not. has_property(at, 'efield_old1')) call add_property(at, 'efield_old1', 0.0_dp, n_cols=3)
  if (.not. has_property(at, 'efield_old2')) call add_property(at, 'efield_old2', 0.0_dp, n_cols=3)
  if (.not. has_property(at, 'efield_old3')) call add_property(at, 'efield_old3', 0.0_dp, n_cols=3)
  
  ! Increment at%cutoff if necessary
  call set_cutoff_minimum(at, max(this%cutoff_ms, this%cutoff_coulomb)*BOHR)
  
end subroutine IPModel_ASAP2_setup_atoms

subroutine IPModel_ASAP2_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_ASAP2), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional, intent(in) :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   type(Dictionary) :: params
   logical :: restart, applied_efield, save_dipole_velo, &
        calc_charge, calc_sc_dipoles, calc_dipoles, calc_short_range
   real(dp), allocatable, target :: efield_int_old(:,:)
   real(dp), allocatable :: efield_charge(:,:), efield_dipole(:,:), dip_sr(:,:)
   real(dp), pointer, dimension(:) :: charge
   real(dp), pointer, dimension(:,:) :: efield, dipoles, efield_old1, efield_old2, efield_old3, ext_efield, dip_velo
   logical, pointer, dimension(:) :: fixdip
   real(dp) :: diff, diff_old
   integer :: n_efield_old
   integer :: i, npol, ti, vv
   character(len=STRING_LENGTH) :: efield_name, dipoles_name, atom_mask_name, source_mask_name
   logical :: pseudise
   real(dp) :: grid_size

   real, parameter :: difftol = 500.0_dp

   INIT_ERROR(error)

   call system_timer('asap_calc')
   vv = current_verbosity()

   if (present(args_str)) then
      call initialise(params)
      call param_register(params, 'efield', 'efield', efield_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'dipoles', 'dipoles', dipoles_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'save_dipole_velo', 'F', save_dipole_velo, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'restart', 'F', restart, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'applied_efield', 'F', applied_efield, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'calc_charge', 'T', calc_charge, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'calc_sc_dipoles', 'T', calc_sc_dipoles, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'calc_dipoles', 'T', calc_dipoles, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'calc_short_range', 'T', calc_short_range, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'atom_mask_name', '', atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'source_mask_name', '', source_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'pseudise', 'F', pseudise, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params, 'grid_size', '0.0', grid_size, help_string="No help yet.  This source file was $LastChangedBy$")

      if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_ASAP2_Calc args_str')) then
         RAISE_ERROR("IPModel_ASAP2_Calc failed to parse args_str="//trim(args_str), error)
      endif
      call finalise(params)
   else
      save_dipole_velo = .false.
      restart = .false.
      applied_efield = .false.
      calc_charge = .true.
      calc_sc_dipoles = .true.
      calc_short_range = .true.
      atom_mask_name = ''
      source_mask_name = ''
      pseudise = .false.
      grid_size = 0.0_dp
   end if

   if (present(e)) e = 0.0_dp
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_ASAP2_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_ASAP2_Calc', error)
      f = 0.0_dp
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_ASAP2_Calc', error)
      local_virial = 0.0_dp
      RAISE_ERROR("IPModel_ASAP2_Calc: local_virial calculation requested but not supported yet.", error)
   endif

   allocate(efield_charge(3,at%n), efield_dipole(3,at%n))

   ! Assign pointers
   if (.not. assign_pointer(at, efield_name, efield)) then
      RAISE_ERROR('IPModel_ASAP2_calc failed to assign pointer to "efield" property', error)
   endif

   if (.not. assign_pointer(at, dipoles_name, dipoles)) then
      RAISE_ERROR('IPModel_ASAP2_calc failed to assign pointer to "dipoles" property', error)
   endif

   if (save_dipole_velo .and. maxval(abs(dipoles)) > 0.0_dp) then
      if (.not. assign_pointer(at, 'dip_velo', dip_velo)) then
           RAISE_ERROR('IPModel_ASAP2_calc failed to assign pointer ot "dip_velo" property', error)
      endif
      do i=1,at%n
         dip_velo(:,i) = -dipoles(:,i)
      end do
   end if

   if (applied_efield) then
      if (.not. assign_pointer(at, 'ext_efield', ext_efield)) then
           RAISE_ERROR('IPModel_ASAP2_calc failed to assign pointer to "ext_efield" property', error)
      endif
   end if

   if (.not. assign_pointer(at, 'charge', charge)) then
        RAISE_ERROR('IPModel_ASAP2_calc failed to assign pointer to "charge" property', error)
      endif
   if (.not. assign_pointer(at, 'fixdip', fixdip)) then
        RAISE_ERROR('IPModel_ASAP2_calc failed to assign pointer to "fixdip" property', error)
      endif
   if (.not. assign_pointer(at, 'efield_old1', efield_old1)) then
        RAISE_ERROR('IPModel_ASAP2_calc failed to assign pointer to "efield_old1" property', error)
      endif
   if (.not. assign_pointer(at, 'efield_old2', efield_old2)) then
        RAISE_ERROR('IPModel_ASAP2_calc failed to assign pointer to "efield_old2" property', error)
      endif
   if (.not. assign_pointer(at, 'efield_old3', efield_old3)) then
        RAISE_ERROR('IPModel_ASAP2_calc failed to assign pointer to "efield_old3" property', error)
      endif

   if (.not. get_value(at%params, 'n_efield_old', n_efield_old)) n_efield_old = 0

   if (restart) then
      efield_old1(:,:) = 0.0_dp
      efield_old2(:,:) = 0.0_dp
      efield_old3(:,:) = 0.0_dp
      n_efield_old = 0
   end if

   efield_dipole = 0.0_dp
   if (calc_sc_dipoles) then
      do i=1,at%n
         if (.not. fixdip(i)) dipoles(:,i) = 0.0_dp
      end do

      ! Extrapolate from the old electric efield
      if (this%pred_order == 0 .or. n_efield_old < 2) then
         efield_dipole = efield_old1(:,:)
      end if
      
      if (this%pred_order == 1 .and. n_efield_old >= 2) then
         efield_dipole = 2.0_dp*efield_old1(:,:) - efield_old2(:,:)
      end if
      
      if (this%pred_order == 2 .and. n_efield_old >= 3) then
         efield_dipole = 3.0_dp*efield_old1(:,:) - 3.0_dp*efield_old2(:,:) + efield_old3(:,:)
      end if

      if (this%pred_order /= 0) then
         if (this%pred_order == 1) then
            efield_old2 = efield_old1
         else if (this%pred_order == 2) then
            efield_old3 = efield_old2
            efield_old2 = efield_old1
         end if
      end if

      call set_value(at%params, 'n_efield_old', min(n_efield_old+1,3))
   end if

   efield_charge = 0.0_dp
   if (calc_charge) then
      if (maxval(abs(this%z)) > 0.0_dp) then
         call yukawa_charges(at, charge, this%cutoff_coulomb, this%yukalpha, this%yuksmoothlength, &
              e, local_e, f, virial, efield_charge, mpi, atom_mask_name, source_mask_name, this%type_of_atomic_num, &
              pseudise, this%pseudise_sigma, grid_size, error=error)
      end if
   end if

   if (calc_sc_dipoles) then
      allocate(dip_sr(3,at%n))
      dip_sr = 0.0_dp
      if (this%tdip_sr) call asap_short_range_dipole_moments(this, at, charge, dip_sr, mpi)
      
      allocate(efield_int_old(3,at%n))
      efield_int_old = 0.0_dp
      
      if (maxval(abs(this%pol)) > 0.0_dp .and. .not. all(fixdip)) then

         call print('Entering ASAP2 Self-consistent dipole loop with '//count(.not. fixdip)//' variable dipole moments', PRINT_VERBOSE)

         ! Self-consistent determination of dipole moments
         diff_old = 1.0_dp
         npol = 1
         call system_timer('asap_self_consistent_dipoles')
         do 
            ! Mix current and previous total efields
            if (npol == 1) then
               efield = efield_dipole + efield_charge
            else
               efield = this%betapol*efield_dipole + &
                    (1.0_dp - this%betapol)*efield_int_old + efield_charge
            end if
            
            ! Add external field if present
            if (applied_efield) then
               do i=1,at%n
                  efield(:,i) = efield(:,i) + ext_efield(:,i)
               end do
            end if

            ! Calculate dipole moment in response to total efield
            do i=1,at%n
               if (fixdip(i)) cycle
               ti = get_type(this%type_of_atomic_num, at%Z(i))
               if (abs(this%pol(ti)) > 0.0_dp) then
                  dipoles(:,i) = efield(:,i)*this%pol(ti) + dip_sr(:,i)
               end if
            end do

            ! Calculate new efield and measure of convergence
            efield_int_old = efield_dipole
            efield_dipole = 0.0_dp
            call yukawa_dipoles(at, charge, dipoles, this%cutoff_coulomb, this%yukalpha, this%yuksmoothlength, &
                 this%pol, this%b_pol, this%c_pol, this%type_of_atomic_num, this%tdip_sr, efield=efield_dipole, &
                 mpi=mpi, error=error)

            diff = 0.0_dp
            do i=1,at%n
               ti = get_type(this%type_of_atomic_num, at%Z(i))
               diff = diff + ((efield_dipole(:,i) - efield_old1(:,i)) .dot. (efield_dipole(:,i) - efield_old1(:,i)))*&
                    this%pol(ti)*this%pol(ti)
            end do
            diff = sqrt(diff/at%n)
            efield_old1 = efield_dipole
            
            if (vv >= PRINT_VERBOSE) then
               write (line,'("Polarisation iteration : ",i5,3e16.8)') npol, diff_old, diff
               call print(line, PRINT_VERBOSE)
            end if
            
            if (diff > difftol) then
               call write(at, 'ipmodel_asap_polarisation_divergence.xyz')
               RAISE_ERROR('IPModel_ASAP2_calc: Polarisation diverges - diff='//diff, error)
            end if
      
            if (abs(diff - diff_old) < this%tolpol) exit
            
            diff_old = diff
            npol = npol + 1
            if (npol >= this%maxipol)  then
               call write(at, 'ipmodel_asap_polarisation_not_converged.xyz')
               RAISE_ERROR('IPModel_ASAP2_calc: Polarisation not converged in '//this%maxipol//' steps - diff='//diff, error)
            endif

         end do
         call system_timer('asap_self_consistent_dipoles')

         ! Save final dipole field for next time
         efield_old1 = efield_dipole
      end if
   end if
      
   efield_dipole = 0.0_dp
   if (calc_dipoles) then
      ! Compute final energy, local energies, forces, virial and electric efield
      if (maxval(abs(dipoles)) > 0.0_dp) then
         call yukawa_dipoles(at, charge, dipoles, this%cutoff_coulomb, this%yukalpha, this%yuksmoothlength, &
              this%pol, this%b_pol, this%c_pol, this%type_of_atomic_num, this%tdip_sr, &
              e, local_e, f, virial, efield_dipole, mpi, atom_mask_name, source_mask_name, pseudise, &
              this%pseudise_sigma, grid_size, error=error)
         
         if (save_dipole_velo) then
            ! dip_velo = dipoles_{N-1} - dipoles_N (we do not divide by timestep here)
            do i=1,at%n
               dip_velo(:,i) = dip_velo(:,i) + dipoles(:,i)
            end do
         end if
      end if
   end if

   efield = efield_charge + efield_dipole
   
   ! Finally, add the short-range contribution
   if (calc_short_range) then
      call asap_morse_stretch(this, at, e, local_e, f, virial, mpi)
   end if

   ! Unit conversion
   if (present(e)) e = e*HARTREE
   if (present(local_e)) local_e = local_e*HARTREE
   if (present(f)) f = f*(HARTREE/BOHR)
   if (present(virial)) virial = virial*HARTREE

   if (allocated(efield_charge)) deallocate(efield_charge)
   if (allocated(efield_dipole)) deallocate(efield_dipole)
   if (allocated(dip_sr)) deallocate(dip_sr)
   if (allocated(efield_int_old)) deallocate(efield_int_old)

   call system_timer('asap_calc')

end subroutine IPModel_ASAP2_Calc


subroutine IPModel_ASAP2_Print(this, file)
  type(IPModel_ASAP2), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_ASAP2 : ASAP2 Potential", file=file)
  call Print("IPModel_ASAP2 : n_types = " // this%n_types, file=file)
  call Print("IPModel_ASAP2 : betapol = "//this%betapol//" maxipol = "//this%maxipol//" tolpol = "//this%tolpol//" pred_order = "//this%pred_order, file=file)
  call Print("IPModel_ASAP2 : yukalpha = "//this%yukalpha//" yuksmoothlength = "//this%yuksmoothlength, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_ASAP2 : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call Print ("IPModel_ASAP2 : pol = "//this%pol(ti), file=file)
    call Print ("IPModel_ASAP2 : z   = "//this%z(ti), file=file)
   call verbosity_push_decrement()
    do tj =1,this%n_types
       call Print ("IPModel_ASAP2 : pair interaction ti tj " // ti // " " // tj // " Zi Zj " // this%atomic_num(ti) //&
            " " // this%atomic_num(tj), file=file)
       call Print ("IPModel_ASAP2 : pair " // this%D_ms(ti,tj) // " " // this%gamma_ms(ti,tj) // " " &
            // this%R_ms(ti,tj) // " " // this%B_pol(ti,tj) // " " // this%C_pol(ti, tj) // " " // this%pseudise_sigma(ti, tj), file=file)

    end do
   call verbosity_pop()
  end do

end subroutine IPModel_ASAP2_Print

subroutine IPModel_ASAP2_read_params_xml(this, param_str)
  type(IPModel_ASAP2), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type(xml_t) :: fxml

  if (len(trim(param_str)) <= 0) return

  parse_in_ip = .false. 
  parse_matched_label = .false.
  parse_ip => this

  call open_xml_string(fxml, param_str)
  call parse(fxml,  &
    startElement_handler = IPModel_startElement_handler, &
    endElement_handler = IPModel_endElement_handler)
  call close_xml_t(fxml)

  if (this%n_types == 0) then
    call system_abort("IPModel_ASAP2_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_ASAP2_read_params_xml

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=FIELD_LENGTH) :: value

  integer ti, tj, Zi, Zj

  if (name == 'ASAP_params') then ! new ASAP stanza

    if (parse_matched_label) return ! we already found an exact match for this label

    call QUIP_FoX_get_value(attributes, 'label', value, status)
    if (status /= 0) value = ''

    if (len(trim(parse_ip%label)) > 0) then ! we were passed in a label
      if (value == parse_ip%label) then ! exact match
        parse_matched_label = .true.
        parse_in_ip = .true.
      else ! no match
        parse_in_ip = .false.
      endif
    else ! no label passed in
      parse_in_ip = .true.
    endif

    if (parse_in_ip) then
      if (parse_ip%n_types /= 0) then
        call finalise(parse_ip)
      endif

      call QUIP_FoX_get_value(attributes, 'n_types', value, status)
      if (status == 0) then
        read (value, *), parse_ip%n_types
      else
        call system_abort("Can't find n_types in ASAP_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      allocate(parse_ip%pol(parse_ip%n_types))
      parse_ip%pol = 0.0_dp
      allocate(parse_ip%z(parse_ip%n_types))

      allocate(parse_ip%D_ms(parse_ip%n_types,parse_ip%n_types))
      parse_ip%D_ms = 0.0_dp
      allocate(parse_ip%gamma_ms(parse_ip%n_types,parse_ip%n_types))
      parse_ip%gamma_ms = 0.0_dp
      allocate(parse_ip%R_ms(parse_ip%n_types,parse_ip%n_types))
      parse_ip%R_ms = 0.0_dp
      allocate(parse_ip%B_pol(parse_ip%n_types,parse_ip%n_types))
      parse_ip%B_pol = 0.0_dp
      allocate(parse_ip%C_pol(parse_ip%n_types,parse_ip%n_types))
      parse_ip%C_pol = 0.0_dp
      allocate(parse_ip%pseudise_sigma(parse_ip%n_types,parse_ip%n_types))
      parse_ip%pseudise_sigma = 0.0_dp

      call QUIP_FoX_get_value(attributes, "cutoff_coulomb", value, status)
      if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find cutoff_coulomb")
      read (value, *) parse_ip%cutoff_coulomb

      call QUIP_FoX_get_value(attributes, "cutoff_ms", value, status)
      if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find cutoff_ms")
      read (value, *) parse_ip%cutoff_ms

      call QUIP_FoX_get_value(attributes, "betapol", value, status)
      if (status == 0) read (value, *) parse_ip%betapol

      call QUIP_FoX_get_value(attributes, "maxipol", value, status)
      if (status == 0) read (value, *) parse_ip%maxipol

      call QUIP_FoX_get_value(attributes, "tolpol", value, status)
      if (status == 0) read (value, *) parse_ip%tolpol

      call QUIP_FoX_get_value(attributes, "pred_order", value, status)
      if (status == 0) read (value, *) parse_ip%pred_order

      call QUIP_FoX_get_value(attributes, "yukalpha", value, status)
      if (status == 0) read (value, *) parse_ip%yukalpha

      call QUIP_FoX_get_value(attributes, "yuksmoothlength", value, status)
      if (status == 0) read (value, *) parse_ip%yuksmoothlength

      parse_ip%tdip_sr = .true.
      call QUIP_FoX_get_value(attributes, "tdip_sr", value, status)
      if (status == 0) read (value, *), parse_ip%tdip_sr

    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    call QUIP_FoX_get_value(attributes, "pol", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find pol")
    read (value, *) parse_ip%pol(ti)

    call QUIP_FoX_get_value(attributes, "z", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find z")
    read (value, *) parse_ip%z(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
        parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "atnum_i", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_read_params_xml cannot find atnum_i")
    read (value, *) Zi
    call QUIP_FoX_get_value(attributes, "atnum_j", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_read_params_xml cannot find atnum_j")
    read (value, *) Zj

    ti = get_type(parse_ip%type_of_atomic_num,Zi)
    tj = get_type(parse_ip%type_of_atomic_num,Zj)

    call QUIP_FoX_get_value(attributes, "D_ms", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find D_ms")
    read (value, *) parse_ip%D_ms(ti,tj)
    call QUIP_FoX_get_value(attributes, "gamma_ms", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find gamma_ms")
    read (value, *) parse_ip%gamma_ms(ti,tj)
    call QUIP_FoX_get_value(attributes, "R_ms", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find R_ms")
    read (value, *) parse_ip%R_ms(ti,tj)
    call QUIP_FoX_get_value(attributes, "B_pol", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find B_pol")
    read (value, *) parse_ip%B_pol(ti,tj)
    call QUIP_FoX_get_value(attributes, "C_pol", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find C_pol")
    read (value, *) parse_ip%C_pol(ti,tj)

    call QUIP_FoX_get_value(attributes, "pseudise_sigma", value, status)
    if (status == 0) then
       read (value, *) parse_ip%pseudise_sigma(ti,tj)
    end if

    if (ti /= tj) then
      parse_ip%D_ms(tj,ti) = parse_ip%D_ms(ti,tj)
      parse_ip%gamma_ms(tj,ti) = parse_ip%gamma_ms(ti,tj)
      parse_ip%R_ms(tj,ti) = parse_ip%R_ms(ti,tj)
      parse_ip%B_pol(tj,ti) = parse_ip%B_pol(ti,tj)
      parse_ip%C_pol(tj,ti) = parse_ip%C_pol(ti,tj)
      parse_ip%pseudise_sigma(tj,ti) = parse_ip%pseudise_sigma(ti,tj)
    endif

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'ASAP_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_ASAP2_module
