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
!X IPModel_FC4 module  
!X
!% Module for the Fourth-order force-constant potential, described in
!% Esfarjani, Chen and Stokes, Phys. Rev. B {\bf 84} 085204 (2011)
!%
!% Calculation works as follows: The atom positions and velocities in
!% 'at' correspond to a supercell of the primitive lattice, read into
!% 'ideal_struct'.  The 'cell_offset' property of 'at' contains the
!% multiples of the lattice vectors mapping this atom from the
!% primitive cell to the supercell, and the 'prim_index' property is
!% the index of the atom within this cell.  A third file contains the
!% neighbours (described as 'neighbourshell', 'atom_shl' in Esfarjani)
!% to which force constants are computed.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_FC4_module

use error_module
use system_module, only : dp, inoutput, print, verbosity_push_decrement, verbosity_pop, split_string_simple, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use minimization_module, only: KahanSum
use atoms_types_module
use atoms_module
use cinoutput_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_FC4
type IPModel_FC4
  integer :: n_types = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)

  real(dp) :: cutoff = 0.0_dp

  character(len=STRING_LENGTH) :: label

  character(len=STRING_LENGTH) :: ideal_struct_file
  type(Atoms) :: ideal_struct
  integer, pointer :: ideal_struct_prim_index(:)
  integer, pointer :: ideal_struct_cell_offset(:,:)
 
  character(len=STRING_LENGTH) :: atom_shl_file
  type(Atoms) :: atom_shl
  integer, pointer :: atom_shl_prim_index(:)
  integer, pointer :: atom_shl_cell_offset(:,:)

  character(len=STRING_LENGTH) :: fc2_file, fc3_file, fc4_file

  integer :: nfc2, nfc2_indep, nfc3, nfc3_indep, nfc4, nfc4_indep

  integer , allocatable :: igroup_2(:)
  integer , allocatable :: iatomterm_2(:,:)
  integer , allocatable :: ixyzterm_2(:,:)
  real(dp), allocatable :: fcs_2(:)
  real(dp), allocatable :: ampterm_2(:)

  integer , allocatable :: igroup_3(:)
  integer , allocatable :: iatomterm_3(:,:)
  integer , allocatable :: ixyzterm_3(:,:)
  real(dp), allocatable :: fcs_3(:)
  real(dp), allocatable :: ampterm_3(:)

  integer , allocatable :: igroup_4(:)
  integer , allocatable :: iatomterm_4(:,:)
  integer , allocatable :: ixyzterm_4(:,:)
  real(dp), allocatable :: fcs_4(:)
  real(dp), allocatable :: ampterm_4(:)

  integer, allocatable :: findatom_sc_array(:,:,:,:)
  integer sc_min(3), sc_max(3)

  real(dp) :: superlattice(3,3), inv_superlattice(3,3)

end type IPModel_FC4

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_FC4), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_FC4_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_FC4_Finalise
end interface Finalise

interface Print
  module procedure IPModel_FC4_Print
end interface Print

interface Calc
  module procedure IPModel_FC4_Calc
end interface Calc

contains

subroutine IPModel_FC4_Initialise_str(this, args_str, param_str)
  type(IPModel_FC4), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  integer :: ntau

  call Finalise(this)

  this%label = ''
  this%ideal_struct_file = ''
  this%atom_shl_file = ''
  this%fc2_file = ''
  this%fc3_file = ''
  this%fc4_file = ''

  call initialise(params)

  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, "ideal_struct_file", PARAM_MANDATORY, this%ideal_struct_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, "atom_shl_file", PARAM_MANDATORY, this%atom_shl_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, "fc2_file", PARAM_MANDATORY, this%fc2_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, "fc3_file", '', this%fc3_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, "fc4_file", '', this%fc4_file, help_string="No help yet.  This source file was $LastChangedBy$")

  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_FC4_Initialise_str args_str')) then
    call system_abort("IPModel_FC4_Initialise_str failed to find mandatory ideal_struct_file or fc2_file, or parse label from args_str="//trim(args_str))
  endif

  call finalise(params)

  call IPModel_FC4_read_params_xml(this, param_str)

  !
  ! Read the ideal structure
  !
  call read(this%ideal_struct, trim(this%ideal_struct_file))
  call get_param_value(this%ideal_struct, "Superlattice", this%superlattice)
  call matrix3x3_inverse(this%superlattice, this%inv_superlattice)
  call add_property(this%ideal_struct, "cell_offset", value=0, n_cols=3, overwrite=.false.)
  call add_property(this%ideal_struct, "prim_index", value=0, n_cols=1, overwrite=.false.)
  if (.not.assign_pointer(this%ideal_struct, "prim_index", this%ideal_struct_prim_index)) &
       call system_abort("IPModel_FC4 requires 'prim_index' field in ideal struct file")
  if (.not.assign_pointer(this%ideal_struct, "cell_offset", this%ideal_struct_cell_offset)) &
       call system_abort("IPModel_FC4 requires 'cell_offset' field in ideal struct file")

  allocate(this%igroup_2(this%nfc2),    &
       this%iatomterm_2(2,this%nfc2),   &
       this%ixyzterm_2(2,this%nfc2),    &
       this%fcs_2(this%nfc2),           &
       this%ampterm_2(this%nfc2))       

  allocate(this%igroup_3(this%nfc3),    &
       this%iatomterm_3(3,this%nfc3),   &
       this%ixyzterm_3(3,this%nfc3),    &
       this%fcs_3(this%nfc3),           &
       this%ampterm_3(this%nfc3))       

  allocate(this%igroup_4(this%nfc4),    &
       this%iatomterm_4(4,this%nfc4),   &
       this%ixyzterm_4(4,this%nfc4),    &
       this%fcs_4(this%nfc4),           &
       this%ampterm_4(this%nfc4))       

  this%sc_min = minval(this%ideal_struct_cell_offset,2)
  this%sc_max = maxval(this%ideal_struct_cell_offset,2)
  ntau = maxval(this%ideal_struct_prim_index)

  allocate(this%findatom_sc_array(ntau, this%sc_min(1):this%sc_max(1),  &
                                        this%sc_min(2):this%sc_max(2),  &
                                        this%sc_min(3):this%sc_max(3)))

  call fill_findatom_sc_array(this%findatom_sc_array, ntau, this%sc_min, this%sc_max, this%ideal_struct)

  !
  ! Read the atom_shl atoms
  !
  call read(this%atom_shl, trim(this%atom_shl_file))
  call add_property(this%atom_shl, "cell_offset", value=0, n_cols=3, overwrite=.false.)
  call add_property(this%atom_shl, "prim_index", value=0, n_cols=1, overwrite=.false.)
  if (.not.assign_pointer(this%atom_shl, "prim_index", this%atom_shl_prim_index)) &
       call system_abort("IPModel_FC4 requires 'prim_index' field in atom_shl file")
  if (.not.assign_pointer(this%atom_shl, "cell_offset", this%atom_shl_cell_offset)) &
       call system_abort("IPModel_FC4 requires 'cell_offset' field in atom_shl file")
 
  !
  ! Read the force constants
  ! 
  call read_fcs(2, this%nfc2, trim(this%fc2_file), this%igroup_2, this%iatomterm_2, this%ixyzterm_2, this%fcs_2, this%ampterm_2)
  if (this%nfc3.gt.0) call read_fcs(3, this%nfc3, trim(this%fc3_file), this%igroup_3, this%iatomterm_3, this%ixyzterm_3, this%fcs_3, this%ampterm_3)
  if (this%nfc4.gt.0) call read_fcs(4, this%nfc4, trim(this%fc4_file), this%igroup_4, this%iatomterm_4, this%ixyzterm_4, this%fcs_4, this%ampterm_4)

end subroutine IPModel_FC4_Initialise_str

subroutine IPModel_FC4_Finalise(this)
  type(IPModel_FC4), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if (allocated(this%igroup_2)) deallocate(this%igroup_2)
  if (allocated(this%iatomterm_2)) deallocate(this%iatomterm_2)
  if (allocated(this%ixyzterm_2)) deallocate(this%ixyzterm_2)
  if (allocated(this%fcs_2)) deallocate(this%fcs_2)
  if (allocated(this%ampterm_2)) deallocate(this%ampterm_2)

  if (allocated(this%igroup_3)) deallocate(this%igroup_3)
  if (allocated(this%iatomterm_3)) deallocate(this%iatomterm_3)
  if (allocated(this%ixyzterm_3)) deallocate(this%ixyzterm_3)
  if (allocated(this%fcs_3)) deallocate(this%fcs_3)
  if (allocated(this%ampterm_3)) deallocate(this%ampterm_3)

  if (allocated(this%igroup_4)) deallocate(this%igroup_4)
  if (allocated(this%iatomterm_4)) deallocate(this%iatomterm_4)
  if (allocated(this%ixyzterm_4)) deallocate(this%ixyzterm_4)
  if (allocated(this%fcs_4)) deallocate(this%fcs_4)
  if (allocated(this%ampterm_4)) deallocate(this%ampterm_4)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_FC4_Finalise


subroutine IPModel_FC4_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_FC4), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   integer :: ia, ja, ka, la, i0, j0, k0, l0, alpha, beta, gamma, delta
   integer ix
   integer :: ifc, ni(3), taui, tauj
   real(dp) :: disp_alpha_i, disp_beta_j, phi, de, df, fsum(3)
   real(dp) :: displacement(3,at%N)
   real(dp) :: local_virial_sq(3,3,at%N)
   real(dp) :: local_e_tmp(at%N)

   real(dp) :: r_ij(3), com_offset(3)
   integer n_extra_calcs, i_calc
   character(len=20) :: extra_calcs_list(10)

   logical :: do_flux = .false.
   real(dp), pointer :: velo(:,:)
   real(dp) :: flux(3)
!  real(dp), allocatable :: local_flux(:,:)

   INIT_ERROR(error)

!  allocate(local_flux(3,at%N))

   if (present(e)) e = 0.0_dp
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_FC4_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_FC4_Calc', error)
      f = 0.0_dp
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_FC4_Calc', error)
      local_virial = 0.0_dp
   endif
   
   if (present(args_str)) then
      if (len_trim(args_str) > 0) then
         n_extra_calcs = parse_extra_calcs(args_str, extra_calcs_list)
         if (n_extra_calcs > 0) then
            do i_calc=1, n_extra_calcs
               select case(trim(extra_calcs_list(i_calc)))
               case("flux")
                  if (.not. assign_pointer(at, "velo", velo)) then
                     RAISE_ERROR("IPModel_LJ_Calc Flux calculation requires velo field", error)
                  endif
                  do_flux = .true.
               case default
                  RAISE_ERROR("Unsupported extra_calc '"//trim(extra_calcs_list(i_calc))//"'", error)
               end select
            end do
         endif ! n_extra_calcs
      endif ! len_trim(args_str)
   end if ! present(args_str)

   local_virial_sq = 0
   local_e_tmp = 0
   flux = 0
!  local_flux = 0

!  Remap positions by the mean position to avoid loss of precision in
!  the forces if there is a drift.  Don't use 'centre_of_mass' because
!  the `mass' property might not be present.  Using the precise centre
!  of mass is not important.
   do ix=1,3
      com_offset(ix) = KahanSum( (/ (at%pos(ix,ia) - this%ideal_struct%pos(ix,ia), ia=1,at%N) /) ) / at%N
   end do

   do ia=1,at%N
      displacement(:,ia) = -diff_min_image(at, ia, this%ideal_struct%pos(:,ia) + com_offset)
   end do

   do ia=1,at%N
      if (present(mpi)) then
         if (mpi%active) then
            if (mod(ia-1, mpi%n_procs) /= mpi%my_proc) cycle
         endif
      endif

      taui = this%ideal_struct_prim_index(ia)  ! the primitive-cell atom index (taui)
      ni = this%ideal_struct_cell_offset(:,ia) ! the image cell of the atom (ni)
      !
      ! second order
      do ifc=1,this%nfc2
         i0 = this%iatomterm_2(1,ifc)
         if (this%atom_shl_prim_index(i0).ne.taui.or.any(this%atom_shl_cell_offset(:,i0).ne.0)) cycle

         j0   = this%iatomterm_2(2,ifc)
         alpha = this%ixyzterm_2(1,ifc)
         beta  = this%ixyzterm_2(2,ifc)
         ja = findatom_sc(this, this%atom_shl_prim_index(j0), this%atom_shl_cell_offset(:,j0) + ni, at)

         phi = this%fcs_2(ifc) * this%ampterm_2(ifc)
         df = -phi * displacement(beta,ja)

         if (present(local_e).or.present(e)) local_e_tmp(ia) = local_e_tmp(ia) + 0.5*phi * displacement(alpha,ia) * displacement(beta,ja)
         if (present(f)) f(alpha,ia) = f(alpha,ia) + df
         if (do_flux) then
            r_ij = diff_min_image(at, ia, ja)
            flux = flux + r_ij(:) * phi * displacement(alpha,ia) * velo(beta,ja) / 2
         end if
         if (present(local_virial).or.present(virial)) then
            local_virial_sq(alpha,:,ia) = local_virial_sq(alpha,:,ia) + at%pos(:,ia) * df 
         end if
      end do

      !
      ! third order
      do ifc=1,this%nfc3
         i0 = this%iatomterm_3(1,ifc)
         if (this%atom_shl_prim_index(i0).ne.taui.or.any(this%atom_shl_cell_offset(:,i0).ne.0)) cycle
         j0   = this%iatomterm_3(2,ifc)
         k0   = this%iatomterm_3(3,ifc)

         alpha = this%ixyzterm_3(1,ifc)
         beta  = this%ixyzterm_3(2,ifc)
         gamma = this%ixyzterm_3(3,ifc)

         ja = findatom_sc(this, this%atom_shl_prim_index(j0), this%atom_shl_cell_offset(:,j0) + ni, at)
         ka = findatom_sc(this, this%atom_shl_prim_index(k0), this%atom_shl_cell_offset(:,k0) + ni, at)

         phi = this%fcs_3(ifc) * this%ampterm_3(ifc)

         df = -phi * displacement(beta,ja) * displacement(gamma,ka) / 2

         if (present(local_e).or.present(e)) local_e_tmp(ia) = local_e_tmp(ia) + phi * displacement(alpha,ia) * displacement(beta,ja) * displacement(gamma,ka) / 6
         if (present(f)) f(alpha,ia) = f(alpha,ia) + df
         if (do_flux) then
            r_ij = diff_min_image(at, ia, ja)
            flux = flux + r_ij(:) * phi * displacement(alpha,ia) * displacement(gamma,ka) * velo(beta,ja) / 3
         end if
         if (present(local_virial).or.present(virial)) then
            local_virial_sq(alpha,:,ia) = local_virial_sq(alpha,:,ia) + at%pos(:,ia) * df 
         end if
      end do      

      !
      ! fourth order
      do ifc=1,this%nfc4
         i0 = this%iatomterm_4(1,ifc)
         if (this%atom_shl_prim_index(i0).ne.taui.or.any(this%atom_shl_cell_offset(:,i0).ne.0)) cycle
         j0   = this%iatomterm_4(2,ifc)
         k0   = this%iatomterm_4(3,ifc)
         l0   = this%iatomterm_4(4,ifc)

         alpha = this%ixyzterm_4(1,ifc)
         beta  = this%ixyzterm_4(2,ifc)
         gamma = this%ixyzterm_4(3,ifc)
         delta = this%ixyzterm_4(4,ifc)

         ja = findatom_sc(this, this%atom_shl_prim_index(j0), this%atom_shl_cell_offset(:,j0) + ni, at)
         ka = findatom_sc(this, this%atom_shl_prim_index(k0), this%atom_shl_cell_offset(:,k0) + ni, at)
         la = findatom_sc(this, this%atom_shl_prim_index(l0), this%atom_shl_cell_offset(:,l0) + ni, at)

         phi = this%fcs_4(ifc) * this%ampterm_4(ifc)

         df = -phi * displacement(beta,ja) * displacement(gamma,ka) * displacement(delta,la) / 6
         if (present(local_e).or.present(e)) local_e_tmp(ia) = local_e_tmp(ia) + &
              phi * displacement(alpha,ia) * displacement(beta,ja) * displacement(gamma,ka) * displacement(delta,la) / 24
         if (present(f)) f(alpha,ia) = f(alpha,ia) + df
         if (do_flux) then
            ! distance between current positions of atoms i and j (same meaning as the anh_md code)
            r_ij = diff_min_image(at, ia, ja)
!           local_flux(:,ia) = local_flux(:,ia) + r_ij * phi * displacement(alpha,ia) * displacement(gamma,ka) * displacement(delta,la) * velo(beta,ja) / 8
            flux = flux + r_ij(:) * phi * displacement(alpha,ia) * displacement(gamma,ka) * displacement(delta,la) * velo(beta,ja) / 8
         end if
         !if (present(local_virial).or.present(virial)) then 
            local_virial_sq(alpha,:,ia) = local_virial_sq(alpha,:,ia) + at%pos(:,ia) * df 
            !local_virial_sq(alpha,:,ja) = local_virial_sq(alpha,:,ja) + at%pos(:,ia) * df / 4
            !local_virial_sq(alpha,:,ka) = local_virial_sq(alpha,:,ka) + at%pos(:,ia) * df / 4
            !local_virial_sq(alpha,:,la) = local_virial_sq(alpha,:,la) + at%pos(:,ia) * df / 4
         !end if
      end do      
   end do

   if (present(local_virial)) local_virial = reshape(local_virial_sq, (/ 9, at%N /))
   if (present(virial)) virial = sum(local_virial_sq,3)
   if (present(local_e)) local_e = local_e_tmp
   if (present(e)) e = KahanSum(local_e_tmp)
  if (present(mpi)) then
     if (present(e)) e = sum(mpi, e)
     if (present(local_e)) call sum_in_place(mpi, local_e)
     if (present(virial)) call sum_in_place(mpi, virial)
     if (present(virial)) call sum_in_place(mpi, local_virial)
     if (present(f)) call sum_in_place(mpi, f)
  endif
   if (present(f)) then
      do ix=1,3
         fsum(ix) = KahanSum(f(ix,:)) / at%N
      end do
      do ia=1,at%N
         f(:,ia) = f(:,ia) - fsum
      end do
   end if


  if (do_flux) then 
     if (present(mpi)) call sum_in_place(mpi, flux)
!    local_flux = local_flux / cell_volume(at)
!    flux = sum(local_flux,2)
     flux = flux / cell_volume(at)
     call set_value(at%params, "Flux", flux)
  end if

end subroutine IPModel_FC4_Calc

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=1024) :: value

  integer atnum_i, atnum_j, fc_i, ti, tj, max_n_fcs, ti_a(1)

  if (name == 'FC4_params') then ! new FC4 stanza
    if (parse_in_ip) &
      call system_abort("IPModel_startElement_handler entered FC4_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

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
	call system_abort("Can't find n_types in FC4_params")
      endif

      call QUIP_FoX_get_value(attributes, 'cutoff', value, status)
      if (status == 0) then
	read (value, *), parse_ip%cutoff
      else
	call system_abort("Can't find cutoff in FC4_params")
      endif

      call QUIP_FoX_get_value(attributes, 'nfc2', value, status)
      if (status == 0) then
         read (value, *), parse_ip%nfc2
      else
         call system_abort("Can't find nfc2 in FC4_params")
      endif
      call QUIP_FoX_get_value(attributes, 'nfc2_indep', value, status)
      if (status == 0)  then
         read (value, *), parse_ip%nfc2_indep
      else
         parse_ip%nfc2_indep = parse_ip%nfc2
      endif 

      call QUIP_FoX_get_value(attributes, 'nfc3', value, status)
      if (status == 0) then
         read (value, *), parse_ip%nfc3
      else
         parse_ip%nfc3 = 0
      endif
      call QUIP_FoX_get_value(attributes, 'nfc3_indep', value, status)
      if (status == 0)  then
         read (value, *), parse_ip%nfc3_indep
      else
         parse_ip%nfc3_indep = parse_ip%nfc3
      endif 

      call QUIP_FoX_get_value(attributes, 'nfc4', value, status)
      if (status == 0) then
         read (value, *), parse_ip%nfc4
      else
         parse_ip%nfc4 = 0
      endif
      call QUIP_FoX_get_value(attributes, 'nfc4_indep', value, status)
      if (status == 0)  then
         read (value, *), parse_ip%nfc4_indep
      else
         parse_ip%nfc4_indep = parse_ip%nfc4
      endif 

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

    endif ! parse_in_ip
  elseif (parse_in_ip .and. name == 'FC4') then

    call QUIP_FoX_get_value(attributes, "atnum_i", value, status)
    if (status /= 0) call system_abort ("IPModel_FC4_read_params_xml cannot find atnum_i")
    read (value, *) atnum_i
    call QUIP_FoX_get_value(attributes, "atnum_j", value, status)
    if (status /= 0) call system_abort ("IPModel_FC4_read_params_xml cannot find atnum_j")
    read (value, *) atnum_j
    call QUIP_FoX_get_value(attributes, "fc_i", value, status)
    if (status /= 0) call system_abort ("IPModel_FC4_read_params_xml cannot find fc_i")
    read (value, *) fc_i

    ! if (all(parse_ip%atomic_num /= atnum_i)) then
    !   ti_a = minloc(parse_ip%atomic_num)
    !   parse_ip%atomic_num(ti_a(1)) = atnum_i
    ! endif
    ! if (all(parse_ip%atomic_num /= atnum_j)) then
    !   ti_a = minloc(parse_ip%atomic_num)
    !   parse_ip%atomic_num(ti_a(1)) = atnum_j
    ! endif
    ! if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    ! allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    ! parse_ip%type_of_atomic_num = 0
    ! do ti=1, parse_ip%n_types
    !   if (parse_ip%atomic_num(ti) > 0) &
    !     parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    ! end do

    ! ti = parse_ip%type_of_atomic_num(atnum_i)
    ! tj = parse_ip%type_of_atomic_num(atnum_j)

    ! call QUIP_FoX_get_value(attributes, "r0", value, status)
    ! if (status /= 0) call system_abort ("IPModel_FC4_read_params_xml cannot find r0")
    ! read (value, *) parse_ip%r0(ti,tj,fc_i)
    ! call QUIP_FoX_get_value(attributes, "phi2", value, status)
    ! if (status /= 0) call system_abort ("IPModel_FC4_read_params_xml cannot find phi2")
    ! read (value, *) parse_ip%phi2(ti,tj,fc_i)
    ! call QUIP_FoX_get_value(attributes, "phi3", value, status)
    ! if (status /= 0) call system_abort ("IPModel_FC4_read_params_xml cannot find phi3")
    ! read (value, *) parse_ip%phi3(ti,tj,fc_i)
    ! call QUIP_FoX_get_value(attributes, "phi4", value, status)
    ! if (status /= 0) call system_abort ("IPModel_FC4_read_params_xml cannot find phi4")
    ! read (value, *) parse_ip%phi4(ti,tj,fc_i)

    ! if (ti /= tj) then
    !   parse_ip%r0(tj,ti,fc_i) = parse_ip%r0(ti,tj,fc_i)
    !   parse_ip%phi2(tj,ti,fc_i) = parse_ip%phi2(ti,tj,fc_i)
    !   parse_ip%phi3(tj,ti,fc_i) = parse_ip%phi3(ti,tj,fc_i)
    !   parse_ip%phi4(tj,ti,fc_i) = parse_ip%phi4(ti,tj,fc_i)
    ! endif

  endif ! parse_in_ip .and. name = 'FC4'

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if (name == 'FC4_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_FC4_read_params_xml(this, param_str)
  type(IPModel_FC4), intent(inout), target :: this
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
    call system_abort("IPModel_FC4_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_FC4_read_params_xml


subroutine IPModel_FC4_Print(this, file)
  type(IPModel_FC4), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti

  call Print("IPModel_FC4 : FC4 Potential", file=file)
  call Print("IPModel_FC4 : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)
  call Print("IPModel_FC4 : label = " // this%label, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_FC4 : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    call Print ("IPModel_FC4 : " // &
        "cutoff " // this%cutoff , &
         file=file)
    call verbosity_pop()
  end do

end subroutine IPModel_FC4_Print

function parse_extra_calcs(args_str, extra_calcs_list) result(n_extra_calcs)
  character(len=*), intent(in) :: args_str
  character(len=*), intent(out) :: extra_calcs_list(:)
  integer :: n_extra_calcs

  character(len=STRING_LENGTH) :: extra_calcs_str
  type(Dictionary) :: params
  
  n_extra_calcs = 0
  call initialise(params)
  call param_register(params, "extra_calcs", "", extra_calcs_str, help_string="No help yet.  This source file was $LastChangedBy$")
  if (param_read_line(params, args_str, ignore_unknown=.true.,task='parse_extra_calcs')) then
    if (len_trim(extra_calcs_str) > 0) then
       call split_string_simple(extra_calcs_str, extra_calcs_list, n_extra_calcs, ":")
    end if
  end if
  call finalise(params)

end function parse_extra_calcs

! subroutine read_fc2(nfc2, file, igroup_2, iatomterm_2, ixyzterm_2, ampterm_2, fcs_2)
!   character(len=*), intent(in) :: file
!   integer, intent(in) :: nfc2

!   integer, intent(out) :: igroup_2(nfc2)
!   integer, intent(out) :: iatomterm_2(2,nfc2)
!   integer, intent(out) ::  ixyzterm_2(2,nfc2)
!   real(dp), intent(out) :: fcs_2(nfc2)
!   real(dp), intent(out) :: ampterm_2(nfc2)

!   integer t,i ! ,mx
!   character(len=999) line
  
!   integer, parameter :: iunit = 1111
!   open(unit=iunit, file=file)

!   read(iunit,'(a)') line
!   do i=1,nfc2
!      read (iunit,*,err=92) t,igroup_2(i), &
!           iatomterm_2(1,i),ixyzterm_2(1,i),iatomterm_2(2,i),ixyzterm_2(2,i),  &
!           fcs_2(i),ampterm_2(i)
!   end do
!   return
! 92 call system_abort("IPModel_FC4_Initialise_str: End of file while reading second order force constants from "//file//", after t="//t)
! end subroutine read_fc2

subroutine read_fcs(rank, nfc, file, igroup, iatomterm, ixyzterm, ampterm, fcs)
  integer, intent(in) :: rank, nfc
  character(len=*), intent(in) :: file

  integer, intent(out) :: igroup(nfc)
  integer, intent(out) :: iatomterm(rank,nfc)
  integer, intent(out) ::  ixyzterm(rank,nfc)
  real(dp), intent(out) :: fcs(nfc)
  real(dp), intent(out) :: ampterm(nfc)

  integer t,i,r ! ,mx
  character(len=999) line
  
  integer, parameter :: iunit = 1111
  open(unit=iunit, file=file)

  read (iunit,'(a)') line
  read (iunit,*,err=92) (t,igroup(i), (iatomterm(r,i), ixyzterm(r,i), r=1,rank), fcs(i), ampterm(i), i=1,nfc)
  return
92 call system_abort("IPModel_FC4_Initialise_str: End of file while reading force constants from "//file//", after t="//t)
end subroutine read_fcs

subroutine fill_findatom_sc_array(findatom_sc_array, ntau, mn, mx, at)
  integer, intent(in) :: ntau, mn(3), mx(3)
  type(Atoms), intent(in) :: at  ! check: should have 'n(3)' and 'tau' field
  integer, intent(out) :: findatom_sc_array(ntau,mn(1):mx(1),mn(2):mx(2),mn(3):mx(3))

  integer, pointer :: prim_index(:)
  integer, pointer :: cell_offset(:,:)

  integer i

  if(.not.assign_pointer(at, "prim_index", prim_index)) &
       call system_abort("IPModel_FC4 requires 'prim_index' field")
  if(.not.assign_pointer(at, "cell_offset", cell_offset)) &
       call system_abort("IPModel_FC4 requires 'cell_offset' field")

  findatom_sc_array = -1
  do i=1,at%N
     findatom_sc_array(prim_index(i), cell_offset(1,i), cell_offset(2,i), cell_offset(3,i)) = i
  end do
end subroutine fill_findatom_sc_array

function findatom_sc(this, tau, n, at)
  integer findatom_sc
  type(IPModel_FC4), intent(in) :: this
  integer, intent(in) :: tau, n(3) ! the primitive cell index and offset
  type(Atoms), intent(in) :: at
  integer nwrap(3), nsc(3), tau_sc_idx
  real(dp) :: wrap(3), prim_cell_pos(3), prim_cell(3,3), inv_prim_cell(3,3), tau_prim(3), tau_frac(3), tau_frac_remap(3), n_frac(3)

  namelist/DEBUG/tau,n,prim_cell,inv_prim_cell,tau_sc_idx,tau_prim,tau_frac,tau_frac_remap,prim_cell_pos,wrap,nwrap

! Common case: just see if it is in the array
  if (all(n.ge.this%sc_min.and.n.le.this%sc_max)) then 
     findatom_sc = this%findatom_sc_array(tau,n(1),n(2),n(3))
     if (findatom_sc.ge.0) return
  end if

! Otherwise, try to wrap it

  ! the primitive cell
  prim_cell = this%ideal_struct%lattice
  inv_prim_cell = this%ideal_struct%g
  
  tau_sc_idx = this%findatom_sc_array(tau,0,0,0)
  tau_prim = this%ideal_struct%pos(:,tau_sc_idx)
  
  ! express the primitive cell we asked for as a fractional
  ! coordinates of the simulation box.
  tau_frac = matmul(this%inv_superlattice, matmul(prim_cell, n) + tau_prim)

  ! remap to [0,1)
  tau_frac_remap = modulo(tau_frac, 1.0_dp)

  prim_cell_pos = matmul(this%superlattice, tau_frac_remap) - tau_prim
  
  ! find the offset of the primitive cell this location coresponds to
  wrap = matmul(inv_prim_cell, prim_cell_pos)
  if (any(abs(modulo(wrap+0.5_dp,1.0_dp)-0.5_dp).gt.1D-6)) then
     write (*,'(2(3F18.12,x))') modulo(wrap,1.0_dp), wrap
     call system_abort("Findatom_sc: got fractional box!")
  endif
       
  nwrap = nint(wrap)

  findatom_sc = this%findatom_sc_array(tau,nwrap(1),nwrap(2),nwrap(3))

  ! express in terms of primitive cell coords

  if (findatom_sc < 0) then
     write (*,nml=DEBUG) 
     call system_abort("Findatom_sc: requested atom tau=" // tau &
          // ", n=" // n // " does not exist in the supercell")
  end if

end function findatom_sc

end module IPModel_FC4_module
