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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X wrapper for quasicontinuum code to use QUIP potentials
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module QC_QUIP_Wrapper_module
use system_module, only : dp !! , only : system_initialise, dp, inoutput, INPUT, system_abort, extendable_str, string, read, &
    !! atoms, initialise, calc_connect, assign_pointer, add_property, set_lattice, verbosity_push, verbosity_pop, PRINT_SILENT, operator(//)
use potential_module !! , only : potential, initialise, finalise
use potential_module !! , only : potential, initialise, calc
use mpi_context_module !! , only : mpi_context
#if defined(HAVE_LOCAL_E_MIX) || defined(HAVE_ONIOM)
use system_module, only : dp !! , only : table, append, PRINT_ALWAYS, bfs_grow, int_part, wipe, print, initialise, finalise
use potential_module !! , only : print
use system_module, only : dp !! , only : print_xyz, mainlog, optional_default
#endif
implicit none
private

  type (Potential), save :: pot
  type (Atoms), save :: at
  type (MPI_context), save :: mpi_glob

  real(dp), private :: vacuum_dist = 5.0_dp

  public :: verbosity_push, verbosity_pop, PRINT_SILENT
  public :: QC_QUIP_initialise, QC_QUIP_calc,MakeLine
  type (Potential), save :: pot_ip
#if defined(HAVE_LOCAL_E_MIX) || defined(HAVE_ONIOM)
  type (Potential), save :: pot_qm
  public :: QC_QUIP_initialise_hybrid, QC_QUIP_calc_hybrid
#endif

contains

  subroutine MakeLine(str1,real1,str2,real2,strout)
    character(len=*), intent(in)  :: str1,str2
    real(dp),     intent(in) :: real1,real2
    character(len=*), intent(out) :: strout
    
    strout = trim(str1)//real1//trim(str2)//real2
  end subroutine MakeLine

  subroutine QC_QUIP_initialise(str, err)
    character(len=*), intent(in) :: str
    integer, intent(out), optional :: err

    type(inoutput) :: params
    type(extendable_str) :: params_str

    call system_initialise()
    call Initialise(params, "quip_params.xml", INPUT)
    call Initialise(params_str)
    call read(params_str, params%unit)

    call Initialise(pot, args_str=str, param_str=string(params_str))
    if (present(err)) err = 0
  end subroutine

  subroutine QC_QUIP_calc(Lz, pos, Z, w, local_e, f, err)
    real(dp), intent(in) :: Lz
    real(dp), intent(in) :: pos(:,:)
    integer, intent(in) :: Z(:)
    real(dp), intent(in) :: w(:)
    real(dp), intent(out) :: local_e(:), f(:,:)
    integer, intent(out) :: err

    integer N
    real(dp), pointer :: weight(:)

    if (.not. matching_array_sizes(Z, pos, w, local_e, f, N)) then
      call system_abort("Mismatched array sizes in QC_QUIP_calc")
    endif

    call qc_setup_atoms(at, N, Lz, pos, Z)

    if (.not. assign_pointer(at, "weight", weight)) then
      call add_property(at, "weight", 0.0_dp)
      if (.not. assign_pointer(at, "weight", weight)) call system_abort("QC_QUIP_calc Failed to add weight property to at")
      weight = w
    endif

    call add_property_from_pointer(at, "local_energy", local_e)
    call add_property_from_pointer(at, "force", f)
    call calc(pot, at, args_str="local_energy force", error=err)

  end subroutine

#if defined(HAVE_LOCAL_E_MIX) || defined(HAVE_ONIOM)
  subroutine QC_QUIP_initialise_hybrid(str_ip, str_qm, str_hybrid, lat, Z, pos, err)
    character(len=*), intent(in) :: str_ip, str_qm, str_hybrid
    real(dp), intent(in), optional :: lat(3,3)
    integer, intent(in), optional :: Z(:)
    real(dp), intent(in), optional :: pos(:,:)
    integer, intent(out), optional :: err

    type(inoutput) :: params
    type(extendable_str) :: params_str
    ! type(atoms) :: bulk
    integer :: n_present

    call system_initialise(enable_timing=.true.)

    call initialise(mpi_glob)

    call initialise(params, "quip_params.xml", INPUT)
    call initialise(params_str)
    call read(params_str, params%unit, convert_to_string=.true., mpi_comm = mpi_glob%communicator)

    ! call initialise(pot_ip, str_ip, string(params_str))
    ! call initialise(pot_qm, str_qm, string(params_str), mpi_obj = mpi_glob)

    call finalise(params_str)

    n_present = count ( (/present(lat), present(Z), present(pos)  /) )
    if (n_present == 3) then
      if (.not. matching_array_sizes(Z, pos = pos)) then
	call system_abort("Mismatched array sizes in QC_QUIP_initialise_hybrid")
      endif
call print("QC_QUIP_initialise_hybrid was passed in bulk structure, ignoring it", PRINT_ALWAYS)
!      call initialise(bulk, size(Z), lat)
!      bulk%Z = Z
!      bulk%pos = pos
!      call initialise(pot, str_hybrid, pot_qm, pot_ip, bulk) 
!      call finalise(bulk)
call initialise(pot, str_hybrid//" init_args_pot2="//trim(str_ip)//" init_args_pot1="//trim(str_qm), mpi_obj=mpi_glob)
    else if (n_present == 0) then
      call initialise(pot, str_hybrid//" init_args_pot2="//trim(str_ip)//" init_args_pot1="//trim(str_qm), mpi_obj=mpi_glob)
    else
      call system_abort("QC_QUIP_initialise_hybrid called with some but not all of lat, Z, pos present")
    endif

    call print("QC_QUIP using hybrid potential:")
    call print(pot)

    if (present(err)) err = 0
  end subroutine QC_QUIP_initialise_hybrid

  subroutine QC_QUIP_calc_hybrid(Lz, pos, Z, w, qm_list, local_e, f, qm_region_width, buffer_region_width, err)
    real(dp), intent(in) :: Lz
    real(dp), intent(in) :: pos(:,:)
    integer, intent(in) :: Z(:)
    real(dp), intent(in) :: w(:)
    integer, intent(in) :: qm_list(:)
    real(dp), intent(out) :: local_e(:), f(:,:)
    integer, intent(in) :: qm_region_width, buffer_region_width
    integer, intent(out) :: err

    integer i, N
    real(dp), pointer :: weight(:)
    integer, pointer :: hybrid(:)
    type(table) :: qm_table

    if (.not. matching_array_sizes(Z, pos, w, local_e, f, N)) then
	call system_abort("Mismatched array sizes in QC_QUIP_calc_hybrid")
    endif

    if (buffer_region_width > 0 .and. pot%is_oniom) &
      call system_abort("Don't do buffer_region_width = " // buffer_region_width // &
			 " > 0 with ONIOM")

    call qc_setup_atoms(at, N, Lz, pos, Z)

    call add_property(at, "weight", 0.0_dp)
    call add_property(at, "hybrid", 0)

    if (.not. assign_pointer(at, "weight", weight)) &
      call system_abort("QC_QUIP_calc Failed to add weight property to at")
    if (.not. assign_pointer(at, "hybrid", hybrid)) &
      call system_abort("QC_QUIP_calc Failed to add hybrid property to at")

    call calc_connect(at)

    hybrid = 0
    call wipe(qm_table)
    do i=1, size(qm_list)
      if (qm_list(i) > at%N) call system_abort("QC_QUIP_calc_hybrid got qm_list("//i//")=" &
							      // qm_list(i) // " > at%N="//at%N)
      if (qm_list(i) > 0) call append(qm_table, (/ qm_list(i), 0, 0, 0 /) )
    end do

    hybrid(int_part(qm_table,1)) = 1

    weight = w

    call add_property_from_pointer(at, "local_energy", local_e)
    call add_property_from_pointer(at, "force", f)

    if (pot%is_oniom) then
      call calc(pot, at, args_str="local_energy force calc_weights" // &
	" core_hops="//qm_region_width// " transition_hops=0 buffer_hops="//buffer_region_width, error=err)
    else
      call calc(pot, at, args_str="local_energy force calc_weights" // &
	" core_hops="//qm_region_width// " transition_hops=0 buffer_hops="//buffer_region_width // &
	" solver=DIAG_GF SCF_GLOBAL_U GLOBAL_U=20.0", error=err)
    endif

  end subroutine QC_QUIP_calc_hybrid 
#endif

  function matching_array_sizes(Z, pos, w, local_e, f, N)
    integer, intent(in) :: Z(:)
    real(dp), intent(in), optional :: pos(:,:)
    real(dp), intent(in), optional :: w(:)
    real(dp), intent(in), optional :: local_e(:), f(:,:)
    integer, intent(out), optional :: N
    logical :: matching_array_sizes

    integer my_N

    matching_array_sizes = .true.
    my_N = size(Z)
    if (present(N)) N = my_N
    if (present(pos)) then
      if (3 /= size(pos,1) .or. my_N /= size(pos,2)) then
	matching_array_sizes=.false.
	return
      endif
    endif
    if (present(w)) then
      if (my_N /= size(w)) then
	matching_array_sizes=.false.
	return
      endif
    endif
    if (present(local_e)) then
      if (my_N /= size(local_e)) then
	matching_array_sizes=.false.
	return
      endif
    endif
    if (present(f)) then
      if (3 /= size(f,1) .or. my_N /= size(f,2)) then
	matching_array_sizes=.false.
	return
      endif
    endif

  end function matching_array_sizes

  subroutine qc_setup_atoms(at, N, Lz, pos, Z)
    type(atoms), intent(inout) :: at
    integer, intent(in) :: N
    real(dp), intent(in) :: Lz
    real(dp), intent(in) :: pos(:,:)
    integer, intent(in) :: Z(:)

    real(dp) :: lattice(3,3)

    lattice = 0.0_dp
    lattice(1,1) = maxval(pos(1,:)) - minval(pos(1,:)) + vacuum_dist
    lattice(2,2) = maxval(pos(2,:)) - minval(pos(2,:)) + vacuum_dist
    lattice(3,3) = Lz

    if (at%N /= N) then
      call Initialise(at, N, lattice)
    end if
    call set_lattice(at, lattice, scale_positions=.false.)
    at%pos = pos
    at%Z = Z

    call calc_connect(at)
  end subroutine qc_setup_atoms

end module QC_QUIP_Wrapper_module
