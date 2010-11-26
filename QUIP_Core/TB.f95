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
!X TB module 
!X
!% General object which handles all the possible tight-binding potentials (TB), re-addressing
!% the calls to the right routines. 
!% 
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"
module TB_module

use libatoms_module

use QUIP_common_module

use Functions_module
use ApproxFermi_module
use Matrix_module
use RS_SparseMatrix_module
use TB_KPoints_module
use TBModel_module
use TBMatrix_module
use TBSystem_module
use TB_GreensFunctions_module

implicit none
private

public :: TB_type
type TB_type
  type(TBSystem) tbsys
  type(Atoms) at
  integer method

  type(TBVector) :: evals, E_fillings, F_fillings, eval_F_fillings 
  type(TBMatrix) :: evecs                  
  type(TBMatrix) :: dm, Hdm, scaled_evecs  

  type(MPI_context) :: mpi

  type (GreensFunctions) gf

  logical :: calc_done = .false.
  real(dp) :: fermi_E, fermi_T
  real(dp) :: fermi_E_precision = 1.0e-9_dp
  real(dp) :: homo_e, lumo_e

  character(len=1024) :: init_args_str
  character(len=1024) :: calc_args_str

end type TB_type

public :: Initialise
public :: TB_Initialise_filename
interface Initialise
  module procedure TB_Initialise_inoutput, TB_Initialise_str
end interface Initialise

public :: Finalise
interface Finalise
  module procedure TB_Finalise
end interface Finalise

public :: cutoff
interface cutoff
   module procedure TB_cutoff
end interface

public :: Wipe
interface Wipe
  module procedure TB_Wipe
end interface Wipe

public :: Print
interface Print
  module procedure TB_Print
end interface Print

public :: Setup_atoms
interface Setup_atoms
  module procedure TB_Setup_atoms
end interface Setup_atoms

public :: solve_diag
interface solve_diag
  module procedure TB_solve_diag
end interface

public :: calc
interface calc
  module procedure TB_calc
end interface

public :: calc_diag
interface calc_diag
  module procedure TB_calc_diag
end interface

public :: calc_GF
interface calc_GF
  module procedure TB_calc_GF
end interface

public :: TB_evals

public :: absorption

interface find_fermi_E
  module procedure TB_find_fermi_E
end interface find_fermi_E

interface calc_E_fillings
  module procedure TB_calc_E_fillings
end interface calc_E_fillings

interface calc_F_fillings
  module procedure TB_calc_F_fillings
end interface calc_F_fillings

interface realloc_match_tbsys
  module procedure realloc_match_tbsys_vec, realloc_match_tbsys_mat
end interface realloc_match_tbsys


contains


subroutine TB_Initialise_filename(this, args_str, filename, kpoints_obj, mpi_obj, error)
  type(TB_type), intent(inout) :: this
  character(len=*), intent(in) :: args_str
  character(len=*), intent(in) :: filename
  type(KPoints), intent(in), optional :: kpoints_obj
  type(MPI_context), intent(in), optional :: mpi_obj
  integer, intent(out), optional :: error

  type(inoutput) io

  INIT_ERROR(error)

  call Initialise(io, filename, INPUT)

  call Initialise(this, args_str, io, kpoints_obj, mpi_obj, error)
  PASS_ERROR(error)

  call Finalise(io)

end subroutine TB_Initialise_filename

subroutine TB_Initialise_inoutput(this, args_str, io_obj, kpoints_obj, mpi_obj, error)
  type(TB_type), intent(inout) :: this
  character(len=*), intent(in) :: args_str
  type(Inoutput), intent(inout), optional :: io_obj
  type(KPoints), intent(in), optional :: kpoints_obj
  type(MPI_context), intent(in), optional :: mpi_obj
  integer, intent(out), optional :: error

  type(extendable_str) :: ss

  INIT_ERROR(error)

  call Initialise(ss)
  if (present(io_obj)) then
    if (present(mpi_obj)) then
      call read(ss, io_obj%unit, convert_to_string=.true., mpi_comm = mpi_obj%communicator)
    else
      call read(ss, io_obj%unit, convert_to_string=.true.)
    endif
  endif
  call Initialise(this, args_str, string(ss), kpoints_obj, mpi_obj, error)
  PASS_ERROR(error)

  call Finalise(ss)

end subroutine TB_Initialise_inoutput


subroutine TB_Initialise_str(this, args_str, param_str, kpoints_obj, mpi_obj, error)
  type(TB_type), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(KPoints), intent(in), optional :: kpoints_obj
  type(MPI_context), intent(in), optional :: mpi_obj
  integer, intent(out), optional :: error

  INIT_ERROR(error)

  call Finalise(this)

  this%init_args_str = args_str

  if (present(mpi_obj)) then
    this%mpi = mpi_obj
  endif
  call Initialise(this%tbsys, args_str, param_str, kpoints_obj, mpi_obj)
end subroutine TB_Initialise_str

subroutine TB_Finalise(this)
  type(TB_type), intent(inout) :: this

  call Wipe(this)
  call Finalise(this%tbsys)
  call Finalise(this%mpi)
  call Finalise(this%evals)
  call Finalise(this%E_fillings)
  call Finalise(this%F_fillings)
  call Finalise(this%eval_F_fillings)
  call Finalise(this%evecs)
  call Finalise(this%dm)
  call Finalise(this%Hdm)
  call Finalise(this%gf)
end subroutine

subroutine TB_Wipe(this)
  type(TB_type), intent(inout) :: this

  call Wipe(this%tbsys)
  call Wipe(this%evals)
  call Wipe(this%E_fillings)
  call Wipe(this%F_fillings)
  call Wipe(this%eval_F_fillings)
  call Wipe(this%evecs)
  call Wipe(this%dm)
  call Wipe(this%Hdm)
  call Wipe(this%gf)

  this%calc_done = .false.

end subroutine

function TB_cutoff(this)
  type(TB_type), intent(in) :: this
  real(dp) :: TB_cutoff
  TB_cutoff = 0.0_dp ! return zero, because TB does its own connection calculation
end function TB_cutoff

subroutine TB_Print(this, file)
  type(TB_type),    intent(in)           :: this
  type(Inoutput), intent(inout),optional,target:: file

  if (current_verbosity() < PRINT_NORMAL) return

  call Print('TB : ', file=file)

  if (this%calc_done) then
    call print("Fermi_E " // this%Fermi_E //" Fermi_T " // this%Fermi_T)
  else
    call print("Fermi_E, Fermi_T not yet set (no calc done)")
  endif

  call Print (this%tbsys, file=file)
  if (this%calc_done) then
     call print ("homo " // this%homo_e // " lumo " // this%lumo_e // " gap " // (this%lumo_e-this%homo_e), file=file)
  else
    call print("homo, lumo, gap not yet set (no calc done)")
  endif
  call verbosity_push_decrement()
    call print("evals", file=file)
    call Print (this%evals, file=file)
    call verbosity_push_decrement()
      call print("E_fillings", file=file)
      call Print (this%E_fillings, file=file)
      call print("F_fillings", file=file)
      call Print (this%F_fillings, file=file)
      call print("eval_F_fillings", file=file)
      call Print (this%eval_F_fillings, file=file)
      call verbosity_push_decrement()
	call print("evecs", file=file)
	call Print (this%evecs, file=file)
	call print("dm", file=file)
	call Print (this%dm, file=file)
	call print("Hdm", file=file)
	call Print (this%Hdm, file=file)
	call print("gf", file=file)
	call Print (this%gf, file=file)
      call verbosity_pop()
    call verbosity_pop()
  call verbosity_pop()

end subroutine TB_Print

subroutine TB_Setup_atoms(this, at, is_noncollinear, args_str, error)
  type(TB_type), intent(inout) :: this
  type(Atoms), intent(in) :: at
  logical, intent(in), optional :: is_noncollinear
  character(len=*), intent(in), optional :: args_str
  integer, intent(out), optional :: error

  INIT_ERROR(error)

  call wipe(this%tbsys)
  call setup_atoms(this%tbsys, at, is_noncollinear, args_str, this%mpi, error=error)
  PASS_ERROR(error)

  this%at = at
  call set_cutoff(this%at, this%tbsys%tbmodel%cutoff)
  call calc_connect(this%at, own_neighbour=.true.)

end subroutine TB_setup_atoms


subroutine realloc_match_tbsys_vec(tbsys, vec)
  type(TBSystem), intent(in) :: tbsys
  type(TBVector), intent(inout) :: vec

  if (vec%N /= tbsys%N .or. vec%n_vectors /= tbsys%n_matrices) then
    call Finalise(vec)
    call Initialise(vec, tbsys%N, tbsys%n_matrices)
  endif
end subroutine realloc_match_tbsys_vec

subroutine realloc_match_tbsys_mat(tbsys, mat)
  type(TBSystem), intent(in) :: tbsys
  type(TBMatrix), intent(inout) :: mat

  if (mat%N /= tbsys%N .or. mat%n_matrices /= tbsys%n_matrices) then
    call Finalise(mat)
    call Initialise(mat, tbsys%N, tbsys%n_matrices, tbsys%complex_matrices, &
      scalapack_obj=tbsys%scalapack_my_matrices)
  endif
end subroutine realloc_match_tbsys_mat

subroutine TB_solve_diag(this, need_evecs, use_fermi_E, fermi_E, w_n, use_prev_charge, AF, error)
  type(TB_type), intent(inout) :: this
  logical, intent(in), optional :: need_evecs
  logical, optional, intent(in) :: use_Fermi_E
  real(dp), optional, intent(inout) :: fermi_E
  real(dp), pointer, intent(in) :: w_n(:)
  logical, optional, intent(in) :: use_prev_charge
  type(ApproxFermi), intent(inout), optional :: AF
  integer, intent(out), optional :: error

  logical my_use_prev_charge

  real(dp), pointer :: scf_orbital_n(:), scf_orbital_m(:,:)
  real(dp) :: global_N
  real(dp), pointer :: local_N(:), local_mom(:,:)
  logical do_evecs

  integer diag_error
  integer :: max_iter = 1
  integer iter
  logical scf_converged

  INIT_ERROR(error)

  my_use_prev_charge = optional_default(.false., use_prev_charge)

  if (.not. this%at%use_uniform_cutoff .or. this%at%cutoff < this%tbsys%tbmodel%cutoff) &
    call print ("WARNING: Called TB_solve_diag with messed up cutoff in atoms: uniform" // this%at%use_uniform_cutoff &
      // " " // this%at%cutoff // " < " // this%tbsys%tbmodel%cutoff)

  if (present(AF)) then
    if (present(use_fermi_E)) then
      if (use_fermi_E) then
	if (AF%n_poles == 0) then
	  call system_abort("Called solve_diag with use_fermi_E and AF that is uninitialized")
	endif
      else
	if (AF%band_width .feq. 0.0_dp) then
	  call system_abort("Called solve_diag with .not. use_fermi_E and AF that has no band_width")
	endif
      endif
    endif
  endif

  if (this%tbsys%scf%active) then
    max_iter= this%tbsys%scf%max_iter
  endif

  if (this%tbsys%scf%active) then
    if (my_use_prev_charge) then

      if (assign_pointer(this%at, 'local_N', local_N)) then
	call print("TB_solve_diag calling set_atomic_n_mom(this%tbsys%scf) using this%at:local_N", PRINT_VERBOSE)
      else
	call print("TB_solve_diag got use_prev_charge, but no local_N value is defined", PRINT_ALWAYS)
      endif
      if (assign_pointer(this%at, 'local_mom', local_mom)) then
	call print("TB_solve_diag calling set_atomic_n_mom(this%tbsys%scf) using this%at:local_mom", PRINT_VERBOSE)
      else
	call print("TB_solve_diag got use_prev_charge, but no local_mom value is defined", PRINT_ALWAYS)
      endif
      call scf_set_atomic_n_mom(this%tbsys, local_N, local_mom)

      if (get_value(this%at%params, 'global_N', global_N)) then
	call print("TB_solve_diag calling set_global_N(this%tbsys%scf from this%at:global_N", PRINT_VERBOSE)
	call scf_set_global_N(this%tbsys, w_n, global_N)
      else
	call scf_set_global_N(this%tbsys, w_n)
	call print("TB_solve_diag got use_prev_charge, but no global_N value is defined", PRINT_ALWAYS)
      endif

    else
      nullify(local_N)
      nullify(local_mom)
      call scf_set_atomic_n_mom(this%tbsys, local_N, local_mom)
      call scf_set_global_N(this%tbsys, w_n)
    endif
  endif

  do_evecs = .false.
  if (present(need_evecs)) do_evecs = need_evecs
  do_evecs = do_evecs .or. this%tbsys%scf%active

  call realloc_match_tbsys(this%tbsys, this%evals)
  if (do_evecs) then
    call realloc_match_tbsys(this%tbsys, this%evecs)
  endif

  nullify(scf_orbital_n)
  nullify(scf_orbital_m)
  if (this%tbsys%scf%active) then
    allocate(scf_orbital_n(this%tbsys%N))
    if (this%tbsys%noncollinear) then
      allocate(scf_orbital_m(3,this%tbsys%N))
    endif
  endif
  call fill_sc_matrices(this%tbsys, this%at)

  call calc_orb_local_pot(this%tbsys, w_n)

  iter = 1
  scf_converged = .false.
  if(this%tbsys%scf%active) call print("TB starting SCF iterations")
  do while (iter <= max_iter .and. .not. scf_converged)

    call fill_matrices(this%tbsys, this%at)

    if (this%tbsys%tbmodel%is_orthogonal) then
      if (do_evecs) then
	call diagonalise(this%tbsys%H, this%evals, this%evecs, error = diag_error)
      else
	call diagonalise(this%tbsys%H, this%evals, error = diag_error)
      end if
    else
      if (do_evecs) then
	call diagonalise(this%tbsys%H, this%tbsys%S, this%evals, this%evecs, error = diag_error)
      else
	call diagonalise(this%tbsys%H, this%tbsys%S, this%evals, error = diag_error)
      endif
    endif

    if (diag_error /= 0) then
      if (this%mpi%my_proc == 0) then
	call write(this%at, "atom_dump_bad_diagonalise."//mpi_id()//".xyz")
      endif
      RAISE_ERROR_WITH_KIND(diag_error, "TB_solve_diag got error " // diag_error // " from diagonalise", error)
    endif

    if (this%tbsys%scf%active) then
      call calc_E_fillings(this, use_fermi_e, fermi_E, AF, w_n)
      call calc_dm_from_evecs(this)
      call calc_local_orbital_num(this, scf_orbital_n)
      if (this%tbsys%noncollinear) call calc_local_orbital_mom(this, scf_orbital_m)
      scf_converged = update_orb_local_pot(this%tbsys, this%at, iter, w_n, scf_orbital_n, scf_orbital_m)
    else
      scf_converged = .true.
    endif
    iter = iter + 1

    if (current_verbosity() >= PRINT_VERBOSE) then
      if (this%tbsys%scf%active) then
	if (assign_pointer(this%at, 'local_N', local_N) .or. assign_pointer(this%at, 'local_mom', local_mom)) &
	   call scf_get_atomic_n_mom(this%tbsys, local_N, local_mom)
	if (get_value(this%at%params, 'global_N', global_N)) then
	  call scf_get_global_N(this%tbsys, global_N)
	  call set_value(this%at%params, 'global_N', global_N)
	endif
      endif
      call write(this%at, 'stdout', prefix="SCF_ITER_"//iter)
    endif ! VERBOSE

  end do

  if (this%tbsys%scf%active) then
    if (assign_pointer(this%at, 'local_N', local_N) .or. assign_pointer(this%at, 'local_mom', local_mom)) &
       call scf_get_atomic_n_mom(this%tbsys, local_N, local_mom)
    if (get_value(this%at%params, 'global_N', global_N)) then
      call scf_get_global_N(this%tbsys, global_N)
      call set_value(this%at%params, 'global_N', global_N)
    endif
  endif

  if (associated(scf_orbital_n)) deallocate(scf_orbital_n)
  if (associated(scf_orbital_m)) deallocate(scf_orbital_m)

  if (.not. scf_converged) then
    call print("WARNING: TB_solve_diag failed to converge SCF in TB_solve_diag", PRINT_ALWAYS)
  endif

  if (this%tbsys%scf%active .and. scf_converged ) call print("TB SCF iterations converged")
end subroutine TB_solve_diag

subroutine TB_calc(this, at, energy, local_e, forces, virial, local_virial, args_str, &
  use_fermi_E, fermi_E, fermi_T, band_width, AF, error)

  type(TB_type), intent(inout) :: this
  type(Atoms), intent(inout) :: at
  real(dp), intent(out), optional :: energy
  real(dp), intent(out), optional :: local_e(:)
  real(dp), intent(out), optional :: forces(:,:), local_virial(:,:)
  real(dp), intent(out), optional :: virial(3,3)
  character(len=*), intent(in), optional :: args_str
  logical, optional :: use_fermi_E
  real(dp), intent(inout), optional :: fermi_E, fermi_T, band_width
  type(ApproxFermi), intent(inout), optional :: AF
  integer, intent(out), optional :: error

  real(dp) :: my_energy
  logical :: my_use_Fermi_E
  real(dp) :: my_Fermi_E, my_Fermi_T, my_band_width
  character(len=STRING_LENGTH) :: solver_arg
  logical :: noncollinear, use_prev_charge, do_evecs
  type(ApproxFermi) :: my_AF
  logical :: do_at_local_N
  real(dp), pointer :: local_N_p(:)

  type(Dictionary) :: params
  logical :: has_atom_mask_name
  character(FIELD_LENGTH) :: atom_mask_name

  INIT_ERROR(error)

  call system_timer("TB_calc")
  call system_timer("TB_calc/prep")

  if (present(args_str)) then
    this%calc_args_str = args_str
    call initialise(params)
    call param_register(params, 'solver', 'DIAG', solver_arg, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'noncollinear', 'F', noncollinear, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'use_prev_charge', 'F', use_prev_charge, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'do_at_local_N', 'F', do_at_local_N, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'do_evecs', 'F', do_evecs, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='TB_Calc args_str')) then
      call system_abort("TB_calc failed to parse args_str='"//trim(args_str)//"'")
    endif
    call finalise(params)

    call set_type(this%tbsys%scf, args_str)

    if(has_atom_mask_name) then
       RAISE_ERROR('TB_Calc: atom_mask_name found, but not supported', error)
    endif
  else
    this%calc_args_str = ''
    solver_arg = 'DIAG'
    noncollinear = .false.
    do_evecs = .false.
    use_prev_charge = .false.
    do_at_local_N = .false.
  endif
  call setup_atoms(this, at, noncollinear, args_str, error=error)
  PASS_ERROR(error)

  if( present(local_virial) ) then
     RAISE_ERROR("TB_Calc: local_virial calculation requested, but not implemented yet",error)
  endif

  if(current_verbosity() > PRINT_NERD) then
     call write(this%at, "tb_calc_atomslog.xyz", append=.true.)
  end if

  call system_timer("TB_calc/prep")

  ! local_N_p points to at property local_N if do_at_local_N is true and at property local_N is defined
  nullify(local_N_p)
  if (do_at_local_N) then
    if (.not. assign_pointer(this%at, 'local_N', local_N_p)) &
      call system_abort("TB_calc got do_at_local_N=T, but failed to find local_N property")
  endif

  select case (trim(solver_arg))
    case ('DIAG')
      call system_timer("TB_calc/DIAG_calc_diag")
      my_energy = calc_diag(this, use_fermi_E, fermi_E, fermi_T, local_e, local_N_p, forces, virial, &
	use_prev_charge = use_prev_charge, AF=AF, do_evecs = do_evecs, error=error)
      call system_timer("TB_calc/DIAG_calc_diag")
    case ('GF')
      call system_timer("TB_calc/GF_calc_GF")
      if (present(virial)) call system_abort("No virial with GF (yet?)")
      my_energy =  calc_GF(this, use_fermi_E, fermi_E, fermi_T, band_width, local_e, local_N_p, forces, AF)
      call system_timer("TB_calc/GF_calc_GF")
    case ('DIAG_GF')
      call system_timer("TB_calc/DIAG_GF")
      call system_timer("TB_calc/DIAG_GF_prep")
      if (present(virial)) call system_abort("No virial with GF (yet?)")

      if (.not. has_Fermi_T(my_fermi_T, this%tbsys%tbmodel, Fermi_T, this%calc_args_str)) &
	call system_abort("TB_calc called without Fermi_T for model without default Fermi_T")
      if (.not. has_band_width(my_band_width, this%tbsys%tbmodel, band_width, this%calc_args_str)) &
	call system_abort("TB_calc called without band_width for model without default band_width")

      ! default to constant mu, for now at least
      my_use_fermi_E = optional_default(.true., use_fermi_E)

      if (my_use_Fermi_E .and. .not. present(AF)) then
	  if (.not. has_Fermi_E(my_Fermi_E, this%tbsys%tbmodel, Fermi_E, this%calc_args_str)) &
	    call system_abort("Called TB_calc with use_Fermi_E but neither AF nor Fermi_E nor TB Model with default Fermi_E")
      endif

      if (present(AF)) then
	my_AF = AF
      else
	call Initialise(my_AF, my_Fermi_E, my_Fermi_T, my_band_width)
      end if
      call system_timer("TB_calc/DIAG_GF_prep")

      call system_timer("TB_calc/DIAG_GF_calc_diag")
      my_energy = calc_diag(this, my_use_fermi_E, local_e = local_e, local_N = local_N_p, use_prev_charge = use_prev_charge, do_evecs = do_evecs, AF = my_AF)
      call system_timer("TB_calc/DIAG_GF_calc_diag")
      if (present(forces)) then
	call system_timer("TB_calc/DIAG_GF_calc_GF")
	my_energy = calc_GF(this, use_fermi_E = .true., local_e = local_e, local_N = local_N_p, forces = forces, AF = my_AF)
	call system_timer("TB_calc/DIAG_GF_calc_GF")
      endif
      call finalise(my_AF)
      call system_timer("TB_calc/DIAG_GF")
    case default
      call system_abort("TB_calc confused by args '" // trim(solver_arg) // "'")
  end select

  if (present(energy)) energy = my_energy
  call system_timer("TB_calc")

  this%calc_done = .true.

  call copy_atoms_fields(this%at, at)


end subroutine TB_calc

subroutine copy_atoms_fields(from_at, to_at)
  type(Atoms), intent(in):: from_at
  type(Atoms), intent(inout):: to_at

  real(dp), pointer :: from_R1(:), to_R1(:), from_R2(:,:), to_R2(:,:)
  real(dp) :: from_R

  if (assign_pointer(from_at, 'local_N', from_R1) .and. &
      assign_pointer(to_at, 'local_N', to_R1)) then
    to_R1 = from_R1
  endif
  if (assign_pointer(from_at, 'local_mom', from_R2) .and. &
      assign_pointer(to_at, 'local_mom', to_R2)) then
    to_R2 = from_R2
  endif
  if (assign_pointer(from_at, 'local_N', from_R1) .and. &
      assign_pointer(to_at, 'local_N', to_R1)) then
    to_R1 = from_R1
  endif
  if (assign_pointer(from_at, 'local_E', from_R1) .and. &
      assign_pointer(to_at, 'local_E', to_R1)) then
    to_R1 = from_R1
  endif

  if (get_value(from_at%params, 'global_dN', from_R)) then
    call set_value(to_at%params, 'global_dN', from_R)
  endif

end subroutine copy_atoms_fields

function TB_calc_diag(this, use_fermi_E, fermi_E, fermi_T, local_e, local_N, forces, virial, use_prev_charge, AF, do_evecs, error)
  type(TB_type), intent(inout) :: this
  logical, optional :: use_fermi_E
  real(dp), intent(inout), optional :: fermi_E, fermi_T
  real(dp), intent(out), target, optional :: local_e(:)
  real(dp), intent(out), pointer, optional :: local_N(:)
  real(dp), intent(out), optional :: forces(:,:)
  real(dp), intent(out), optional :: virial(3,3)
  logical, optional, intent(in) :: use_prev_charge
  type(ApproxFermi), intent(inout), optional :: AF
  logical, optional :: do_evecs
  integer, intent(out), optional :: error
  real(dp) :: TB_calc_diag

  real(dp), allocatable :: forces_scf(:,:)
  real(dp) :: virial_scf(3,3)
  real(dp), allocatable :: local_e_scf(:), local_e_rep(:)
  real(dp) :: global_e_scf
  real(dp), pointer :: u_local_e(:)
  real(dp) :: u_e, u_e_scf, u_e_rep
  integer :: i
  logical :: u_do_evecs, do_local_N, do_local_e

  type(Dictionary) :: params

  real(dp), pointer :: w_e(:)
  real(dp), pointer :: w_n(:)

  INIT_ERROR(error)

  call system_timer("TB_calc_diag")
  call system_timer("TB_calc_diag/prep")
  call print("TB_calc_diag starting", PRINT_VERBOSE)

  if (.not. assign_pointer(this%at, "weight", w_e)) then
    nullify(w_e)
  else
    call print("Using weight from atom properties for w_e", PRINT_VERBOSE)
  endif
  if (.not. assign_pointer(this%at, "weight_n", w_n)) then
    if (.not. assign_pointer(this%at, "weight", w_n)) then
      nullify(w_n)
    else
      call print("Using weight from atom properties for w_n", PRINT_VERBOSE)
    endif
  else
    call print("Using weight_n from atom properties for w_n", PRINT_VERBOSE)
  endif

  if (associated(w_e) .and. (present(forces) .or. present(virial))) then
    if (maxval(w_e) .ne. 1.0_dp .or. minval(w_e) .ne. 1.0_dp)  then
      call system_abort ("TB_calc_diag can't do forces or virial of weighted energy")
    endif
  endif

  do_local_N = .false.
  if (present(local_N)) then
    if (associated(local_N)) then
      if (size(local_N) == this%at%N) then
	do_local_N = .true.
      else
	if (size(local_N) /= 0) &
	  call system_abort("TB_calc_diag called with size(local_N)="//size(local_N)// &
			    " not 0, and not matching this%at%N="//this%at%N)
      endif
    endif
  endif
  do_local_e = .false.
  if (present(local_e)) then
    if (size(local_e) == this%at%N) then
      do_local_e = .true.
    else
      if (size(local_e) /= 0) &
	call system_abort("TB_calc_diag called with size(local_e)="//size(local_e)// &
	                  " not 0, and not matching this%at%N="//this%at%N)
    endif
  endif

  u_do_evecs = optional_default(.false., do_evecs) .or. present(forces) .or. present(virial) .or. &
	     do_local_e .or. associated(w_e) .or. &
	     do_local_N .or. associated(w_n) .or. this%tbsys%scf%active
  if (u_do_evecs) call print("evecs are needed", PRINT_VERBOSE)

  if (.not. has_fermi_T(this%fermi_T, this%tbsys%tbmodel, fermi_T, this%calc_args_str)) &
    call system_abort("TB_calc_diag called without fermi_T for a TB model without default fermi T")
  call system_timer("TB_calc_diag/prep")

  call initialise(params)
  call param_register(params, 'fermi_e_precision', ''//this%fermi_E_precision, this%fermi_E_precision, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, this%calc_args_str, ignore_unknown=.true.,task='TB_calc_diag fermi_e_precision')) then
    call system_abort("TBModel_has_fermi_e failed to parse this%calc_args_str='"//trim(this%calc_args_str)//"'")
  endif
  call finalise(params)

  call system_timer("TB_calc_diag/solve_diag")
  call print("TB_calc_diag solve_diag", PRINT_VERBOSE)
  call solve_diag(this, u_do_evecs, use_fermi_e, Fermi_E, w_n, use_prev_charge=use_prev_charge, AF=AF, error=error)
  call system_timer("TB_calc_diag/solve_diag")
  PASS_ERROR_WITH_INFO("TB_calc got error from solve_diag", error)

  call system_timer("TB_calc_diag/calc_EFV")
  allocate(local_e_rep(this%at%N))
  do i=1, this%at%N
    local_e_rep(i) = get_local_rep_E(this%tbsys%tbmodel, this%at, i)
  end do

  if (present(fermi_T)) then
    if (fermi_T .fne. this%fermi_T) then
      this%fermi_T = fermi_T
    endif
  endif
  if (present(fermi_E)) then
    if (fermi_E .fne. this%fermi_E) then
      this%fermi_E = fermi_E
    endif
  endif

  call system_timer("TB_calc_diag/calc_EFV/calc_E_fillings")
  call calc_E_fillings(this, use_fermi_E, fermi_E, AF, w_n)
  call system_timer("TB_calc_diag/calc_EFV/calc_E_fillings")
  if (present(AF)) then
    call print("TB_calc_diag using AF Fermi_E " // AF%Fermi_E // " band_width " // AF%band_width // " n_poles " // AF%n_poles, PRINT_VERBOSE)
  else
    call print("TB_calc_diag has Fermi_E " // this%Fermi_E // " Fermi_T " // this%Fermi_T, PRINT_VERBOSE)
  endif

  if (do_local_e .or. associated(w_e) .or. do_local_N .or. associated(w_n)) then
    call print("TB_calc_diag calculating per atom stuff", PRINT_VERBOSE)
    call system_timer("TB_calc_diag/calc_EFV/calc_dm_for_local_stuff")
    call calc_dm_from_evecs(this, .false.)
    call system_timer("TB_calc_diag/calc_EFV/calc_dm_for_local_stuff")

    call system_timer("TB_calc_diag/calc_EFV/calc_local_stuff")
    if (do_local_e) then
      u_local_e => local_e
    else
      allocate(u_local_e(this%at%N))
    endif
    allocate(local_e_scf(this%at%N))

    call calc_local_atomic_energy(this, u_local_e)
    call local_scf_e_correction(this%tbsys, this%at, local_e_scf, global_e_scf, w_n)
    call print("TB_calc_diag got sum(local_energies) eval " // sum(u_local_e) // &
	       " scf " // sum(local_e_scf) // " rep " // sum(local_e_rep), PRINT_VERBOSE)
    u_local_e = u_local_e + local_e_scf + local_e_rep

    if (associated(w_e)) then
      u_local_e = u_local_e + global_e_scf/sum(w_e)
      TB_calc_diag = sum(w_e*u_local_e)
    else
      u_local_e = u_local_e + global_e_scf/size(u_local_e)
      TB_calc_diag = sum(u_local_e)
    endif

    deallocate(local_e_scf)
    if (.not. present(local_e)) deallocate(u_local_e)

    if (do_local_N) call calc_local_atomic_num(this, local_N)

    call system_timer("TB_calc_diag/calc_EFV/calc_local_stuff")
  else ! no local E or N required
    call system_timer("TB_calc_diag/calc_EFV/E_calc")
    u_e = ksum_dup(this%tbsys%kpoints, sum(this%evals%data_d*this%E_fillings%data_d,dim=1))
    u_e_rep = sum(local_e_rep)
    u_e_scf = scf_e_correction(this%tbsys, this%at, w_n)
    call print("TB_calc_diag got energies eval " // u_e // " scf " // u_e_scf // " rep " // u_e_rep, PRINT_VERBOSE)
    TB_calc_diag = u_e + u_e_rep + u_e_scf
    call system_timer("TB_calc_diag/calc_EFV/E_calc")
  endif

  if (present(forces) .or. present(virial)) then
    call print("TB_calc_diag calculating forces", PRINT_VERBOSE)
    call calc_F_fillings(this, .not. this%tbsys%tbmodel%is_orthogonal, AF)
    call system_timer("TB_calc_diag/calc_EFV/calc_dm_for_FV")
    call calc_dm_from_evecs(this, .true.)
    call system_timer("TB_calc_diag/calc_EFV/calc_dm_for_FV")
    if (present(forces)) then
      call system_timer("TB_calc_diag/calc_EFV/calculate_forces")
      forces =  calculate_forces_diag(this)
      allocate(forces_scf(3,this%at%N))
      forces_scf = scf_f_correction(this%tbsys, this%at)
      forces = forces + forces_scf
      deallocate(forces_scf)
      call system_timer("TB_calc_diag/calc_EFV/calculate_forces")
    end if
    if (present(virial)) then
      call system_timer("TB_calc_diag/calc_EFV/calculate_virial")
      virial = calculate_virial_diag(this)
      virial_scf = scf_virial_correction(this%tbsys, this%at)
      virial = virial + virial_scf
    call system_timer("TB_calc_diag/calc_EFV/calculate_virial")
    end if
  endif

  call system_timer("TB_calc_diag/calc_EFV")
  call system_timer("TB_calc_diag")
  call print("TB_calc_diag ending", PRINT_VERBOSE)

end function TB_calc_diag


function TB_calc_GF(this, use_fermi_E, fermi_E, fermi_T, band_width, local_e, local_N, forces, AF, SelfEnergy)
  type(TB_type), intent(inout) :: this
  logical, intent(in), optional :: use_fermi_E
  real(dp), intent(inout), optional :: fermi_E
  real(dp), intent(in), optional :: fermi_T, band_width
  real(dp), intent(out), target, optional :: local_e(:)
  real(dp), intent(out), pointer, optional :: local_N(:)
  real(dp), intent(out), optional :: forces(:,:)
  type(ApproxFermi), intent(inout), optional, target :: AF
  type(TBMatrix), intent(in), optional :: SelfEnergy(:)
  real(dp) :: TB_calc_GF

  logical find_new_fermi_E
  type(ApproxFermi), pointer :: u_AF
  logical local_u_AF
  real(dp), allocatable :: forces_scf(:,:)

  integer i

  real(dp), pointer :: u_local_e(:)
  real(dp), allocatable :: local_e_rep(:), local_e_scf(:)
  real(dp) :: global_e_scf

  real(dp), pointer :: w_e(:)
  real(dp), pointer :: w_n(:)
  real(dp) :: use_band_width

  logical :: do_local_e, do_local_N

  call system_timer("TB_calc_GF")
  call print("TB_calc_GF starting", PRINT_VERBOSE)

  do_local_N = .false.
  if (present(local_N)) then
    if (associated(local_N)) then
      if (size(local_N) == this%at%N) then
	do_local_N = .true.
      else
	if (size(local_N) /= 0) &
	  call system_abort("TB_calc_diag called with size(local_N)="//size(local_N)// &
			    " not 0, and not matching this%at%N="//this%at%N)
      endif
    endif
  endif
  do_local_e = .false.
  if (present(local_e)) then
    if (size(local_e) == this%at%N) then
      do_local_e = .true.
    else
      if (size(local_e) /= 0) &
	call system_abort("TB_calc_diag called with size(local_e)="//size(local_e)// &
	                  " not 0, and not matching this%at%N="//this%at%N)
    endif
  endif

  call system_timer("TB_calc_GF_prep")
  if (.not. assign_pointer(this%at, "weight", w_e)) then
    nullify(w_e)
  else
    call print("Using weights from atom properties for w_e", PRINT_VERBOSE)
  endif
  if (.not. assign_pointer(this%at, "weight", w_n)) then
    nullify(w_n)
  else
    call print("Using weights from atom properties for w_n", PRINT_VERBOSE)
  endif

  find_new_fermi_E = .false.
  if (present(use_fermi_E)) find_new_fermi_E = .not. use_fermi_E

  if (find_new_fermi_E .and. present(forces)) then
    call system_abort("TB_calc_GF an't do constant N with forces yet")
  endif

  if (.not. find_new_fermi_E .and. (.not. has_fermi_E(this%fermi_e, this%tbsys%tbmodel, fermi_E, this%calc_args_str) .and. &
      .not. present(AF))) then
    call system_abort("called calc_GF with use_fermi_E but no fermi_E or AF") 
  endif

  if ((present(fermi_E) .or. present(fermi_T) .or. present(band_width)) .and. present(AF)) then
    call system_abort("calc_GF has AF as well as fermi_e, fermi_T, or band_width present - don't know which one to use")
  endif

  if (present(AF)) then
    u_AF => AF
    local_u_AF = .false.
  else
    allocate(u_AF)
    local_u_AF = .true.
  endif

  if (find_new_fermi_E) then
    call print("TB_calc_GF: finding our own Fermi level", PRINT_VERBOSE)

    if (.not. has_fermi_T(this%fermi_T, this%tbsys%tbmodel, fermi_T, this%calc_args_str)) &
      call system_abort("TB_calc_GF called without Fermi_T for TB Model without default Fermi_T")

    if (.not. has_band_width(u_AF%band_width, this%tbsys%tbmodel, band_width, this%calc_args_str)) &
      call system_abort("TB_calc_GF called without band_width for TB Model without default band_width")

    !call verbosity_push_increment()
    call solve_diag(this, .false., fermi_E = fermi_E, w_n = w_n, AF = u_AF)
    !call verbosity_pop()
  else
    if (has_fermi_E(this%fermi_e, this%tbsys%tbmodel, fermi_E, this%calc_args_str) .and. &
        has_fermi_T(this%fermi_T, this%tbsys%tbmodel, fermi_T, this%calc_args_str) .and. &
        has_band_width(use_band_width, this%tbsys%tbmodel, band_width, this%calc_args_str)) then
      call print("Initialising ApproxFermi using fermi_e " // this%fermi_e // " fermi_T " // this%fermi_T // &
	" band_width "// use_band_width, PRINT_VERBOSE)
      call Initialise(u_AF, this%Fermi_E, this%fermi_T, use_band_width)
    else
      call system_abort("find_new_fermi_E is false, but we are missing fermi_e, fermi_T or band_width for Initialise(AF...)")
    endif
  endif

  call print("TB_calc_diag has Fermi_E " // this%Fermi_E // " Fermi_T " // this%Fermi_T, PRINT_VERBOSE)

  call verbosity_push_decrement()
  call print("TB_calc_GF using ApproxFermi:")
  call print(u_AF)
  call verbosity_pop()

  call Initialise(this%gf, u_AF%z, u_AF%a, this%tbsys)
  this%gf%tbsys%scf = this%tbsys%scf

  if (local_u_AF) then
    call finalise(u_AF)
  endif
  call system_timer("TB_calc_GF_prep")

  call system_timer("TB_calc_GF_calc_Gs")
  call print("TB_calc_GF calculating Gs", PRINT_VERBOSE)
  call calc_Gs(this%gf, this%at, SelfEnergy)
  call system_timer("TB_calc_GF_calc_Gs")

  call system_timer("TB_calc_GF_calc_EFV")
  call print("TB_calc_GF calculating dm from Gs", PRINT_VERBOSE)
  call calc_dm_from_Gs(this%gf)

  call print("TB_calc_GF calculating energy, number", PRINT_VERBOSE)

  if (do_local_e) then
    u_local_e => local_e
  else
    allocate(u_local_e(this%at%N))
  endif

  allocate(local_e_rep(this%at%N))
  allocate(local_e_scf(this%at%N))
  do i=1, this%at%N
    local_e_rep(i) = get_local_rep_E(this%tbsys%tbmodel, this%at, i)
  end do

  call calc_local_atomic_energy_GF(this, u_local_e)
  call local_scf_e_correction(this%tbsys, this%at, local_e_scf, global_e_scf, w_n)
  u_local_e = u_local_e + local_e_scf + local_e_rep

  if (associated(w_e)) then
    u_local_e = u_local_e + global_e_scf/sum(w_e)
    TB_calc_GF = sum(u_local_e*w_e)
  else
    u_local_e = u_local_e + global_e_scf/size(u_local_e)
    TB_calc_GF = sum(u_local_e)
  endif

  if (.not. present(local_e)) deallocate(u_local_e)
  deallocate(local_e_rep)
  deallocate(local_e_scf)

  if (do_local_N) call calc_local_atomic_num_GF(this, local_N)

  if (present(forces)) then
    call print("TB_calc_GF calculating forces", PRINT_VERBOSE)
    forces = calculate_forces_GF(this, w_e, w_n)
    allocate(forces_scf(3,this%at%N))
    forces_scf = scf_f_correction(this%gf%tbsys, this%at)
    forces = forces + forces_scf
    deallocate(forces_scf)
  endif
  call system_timer("TB_calc_GF_calc_EFV")
  call system_timer("TB_calc_GF")

  call print("TB_calc_GF finished", PRINT_VERBOSE)

end function TB_calc_GF

subroutine calc_local_atomic_energy_GF(this, local_e)
  type(TB_type), intent(inout) :: this
  real(dp), intent(inout) :: local_e(:)

  local_e = local_ksum(this%gf%tbsys%kpoints, atom_orbital_sum(this%gf%tbsys, &
			real(partial_TraceMult(this%gf%dm, this%gf%tbsys%H, &
    a_H=.true., b_H=.false.))))
  call ksum_distrib_inplace(this%gf%tbsys%kpoints, local_e)
end subroutine calc_local_atomic_energy_GF

subroutine calc_local_atomic_num_GF(this, local_N)
  type(TB_type), intent(inout) :: this
  real(dp), intent(inout) :: local_N(:)

  local_N = local_ksum(this%gf%tbsys%kpoints, atom_orbital_sum(this%gf%tbsys, &
			real(partial_TraceMult(this%gf%dm, this%gf%tbsys%S, &
    a_H=.true., b_H=.false.))))
  call ksum_distrib_inplace(this%gf%tbsys%kpoints, local_N)
end subroutine calc_local_atomic_num_GF

subroutine calc_local_atomic_energy(this, local_e)
  type(TB_type), intent(inout) :: this
  real(dp), intent(inout) :: local_e(:)

  local_e = local_ksum(this%tbsys%kpoints, atom_orbital_sum(this%tbsys, &
			real(partial_TraceMult(this%dm, this%tbsys%H, &
    a_H=.true., b_H=.false.))))
  call ksum_distrib_inplace(this%tbsys%kpoints, local_e)
end subroutine calc_local_atomic_energy

subroutine calc_local_atomic_num(this, local_N)
  type(TB_type), intent(inout) :: this
  real(dp), intent(inout) :: local_N(:)

  if (this%tbsys%tbmodel%is_orthogonal) then
    local_N = local_ksum(this%tbsys%kpoints, atom_orbital_sum(this%tbsys, Re_diag(this%dm)))
  else
    local_N = local_ksum(this%tbsys%kpoints, atom_orbital_sum(this%tbsys, &
			  real(partial_TraceMult(this%dm, this%tbsys%S, &
      a_H=.true., b_H=.false.))))
  endif
  call ksum_distrib_inplace(this%tbsys%kpoints, local_N)
end subroutine calc_local_atomic_num

subroutine calc_local_orbital_num(this, local_N)
  type(TB_type), intent(inout) :: this
  real(dp), intent(inout) :: local_N(:)

  if (this%tbsys%tbmodel%is_orthogonal) then
    local_N = local_ksum(this%tbsys%kpoints, Re_diag(this%dm))
  else
    local_N = local_ksum(this%tbsys%kpoints, real(partial_TraceMult(this%dm, this%tbsys%S, a_H=.true., b_H=.false.)))
  endif
  call ksum_distrib_inplace(this%tbsys%kpoints, local_N)
end subroutine calc_local_orbital_num

subroutine calc_local_orbital_mom(this, local_m)
  type(TB_type), intent(inout) :: this
  real(dp), intent(inout) :: local_m(:,:)

  integer :: i
  complex(dp), allocatable :: local_spinor(:,:,:)

  local_m = 0.0_dp

  if (.not. this%tbsys%noncollinear) then
    return
  endif

  call fill_matrices(this%tbsys, this%at, need_H=.false., need_S=.true., no_S_spin=.true.)

  allocate(local_spinor(2,2,this%tbsys%N/2))

  if (this%tbsys%tbmodel%is_orthogonal) then
    local_spinor = local_ksum(this%tbsys%kpoints, diag_spinor(this%dm))
  else
    local_spinor = local_ksum(this%tbsys%kpoints, partial_TraceMult_spinor(this%dm, this%tbsys%S, a_H=.true., b_H=.false.))
  endif

  do i=1, this%tbsys%N/2
    local_m(1,(i-1)*2+1) = sum(conjg(pauli_sigma(:,:,1))*local_spinor(:,:,i))
    local_m(2,(i-1)*2+1) = sum(conjg(pauli_sigma(:,:,2))*local_spinor(:,:,i))
    local_m(3,(i-1)*2+1) = sum(conjg(pauli_sigma(:,:,3))*local_spinor(:,:,i))
  end do
  call ksum_distrib_inplace(this%tbsys%kpoints, local_m)

  deallocate(local_spinor)

  call fill_matrices(this%tbsys, this%at, need_H=.false., need_S=.true., no_S_spin=.false.)

end subroutine calc_local_orbital_mom


function calculate_forces_diag(this) result(forces)
  type(TB_type), intent(inout) :: this
  real(dp) :: forces(3,this%at%N) ! result

  logical, allocatable :: od_mask(:), d_mask(:)

  integer i

  if (.not. this%at%use_uniform_cutoff .or. this%at%cutoff < this%tbsys%tbmodel%cutoff) &
    call print ("WARNING: Called calculate_forces_diag with messed up cutoff in atoms: uniform" // &
      this%at%use_uniform_cutoff // " " // this%at%cutoff // " < " // this%tbsys%tbmodel%cutoff)

  forces = 0.0_dp
  allocate(od_mask(this%at%N))
  allocate(d_mask(this%at%N))

  call Setup_deriv_matrices(this%tbsys, this%at)
  do i=1, this%at%N
    call fill_dmatrices(this%tbsys, this%at, i, diag_mask=d_mask, offdiag_mask=od_mask)
    forces(1,i) = forces(1,i) - local_ksum(this%tbsys%kpoints, real(TraceMult(this%dm,this%tbsys%dH(1), &
      a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
    forces(2,i) = forces(2,i) - local_ksum(this%tbsys%kpoints, real(TraceMult(this%dm,this%tbsys%dH(2), &
      a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
    forces(3,i) = forces(3,i) - local_ksum(this%tbsys%kpoints, real(TraceMult(this%dm,this%tbsys%dH(3), &
      a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
    if (.not. this%tbsys%tbmodel%is_orthogonal) then
      forces(1,i) = forces(1,i) + local_ksum(this%tbsys%kpoints, real(TraceMult(this%Hdm,this%tbsys%dS(1), &
	a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
      forces(2,i) = forces(2,i) + local_ksum(this%tbsys%kpoints, real(TraceMult(this%Hdm,this%tbsys%dS(2), &
	a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
      forces(3,i) = forces(3,i) + local_ksum(this%tbsys%kpoints, real(TraceMult(this%Hdm,this%tbsys%dS(3), &
	a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
    endif
  end do

  deallocate(od_mask)
  deallocate(d_mask)

  call ksum_distrib_inplace(this%tbsys%kpoints, forces)

  do i=1, this%at%N
    forces = forces + get_local_rep_E_force(this%tbsys%tbmodel, this%at, i)
  end do

end function calculate_forces_diag

function calculate_virial_diag(this) result(virial)
  type(TB_type), intent(inout) :: this
  real(dp) :: virial(3,3) ! result

  logical, allocatable :: od_mask(:), d_mask(:)

  integer i

  if (.not. this%at%use_uniform_cutoff .or. this%at%cutoff < this%tbsys%tbmodel%cutoff) &
    call print ("WARNING: Called calculate_virial_diag with messed up cutoff in atoms: uniform" // this%at%use_uniform_cutoff &
      // " " // this%at%cutoff // " < " // this%tbsys%tbmodel%cutoff)

  virial = 0.0_dp
  allocate(od_mask(this%at%N))
  allocate(d_mask(this%at%N))

  call Setup_deriv_matrices(this%tbsys, this%at)
  do i=1, 3
    call fill_dmatrices(this%tbsys, this%at, -i, diag_mask=d_mask, offdiag_mask=od_mask)

    virial(1,i) = virial(1,i) - local_ksum(this%tbsys%kpoints, real(TraceMult(this%dm,this%tbsys%dH(1), &
      a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
    virial(2,i) = virial(2,i) - local_ksum(this%tbsys%kpoints, real(TraceMult(this%dm,this%tbsys%dH(2), &
      a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
    virial(3,i) = virial(3,i) - local_ksum(this%tbsys%kpoints, real(TraceMult(this%dm,this%tbsys%dH(3), &
      a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
    if (.not. this%tbsys%tbmodel%is_orthogonal) then
      virial(1,i) = virial(1,i) + local_ksum(this%tbsys%kpoints, real(TraceMult(this%Hdm,this%tbsys%dS(1), &
	a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
      virial(2,i) = virial(2,i) + local_ksum(this%tbsys%kpoints, real(TraceMult(this%Hdm,this%tbsys%dS(2), &
	a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
      virial(3,i) = virial(3,i) + local_ksum(this%tbsys%kpoints, real(TraceMult(this%Hdm,this%tbsys%dS(3), &
	a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
    endif
  end do

  call ksum_distrib_inplace(this%tbsys%kpoints, virial)
  do i=1, this%at%N
    virial = virial + get_local_rep_E_virial(this%tbsys%tbmodel, this%at, i)
  end do

  deallocate(od_mask)
  deallocate(d_mask)

end function calculate_virial_diag

function calculate_forces_GF(this, w_e, w_n) result (forces)
  type(TB_type), intent(inout) :: this
  real(dp), intent(in), pointer :: w_e(:), w_n(:)
  real(dp) :: forces(3,this%at%N) ! result

  logical, allocatable :: od_mask(:), d_mask(:)

  integer i

  call system_timer("calculate_forces_GF")

  if (.not. this%at%use_uniform_cutoff .or. this%at%cutoff < this%gf%tbsys%tbmodel%cutoff) &
    call print ("WARNING: Called calculate_forces_GF with messed up cutoff in atoms: uniform" // &
      this%at%use_uniform_cutoff // " " // this%at%cutoff // " < " // this%gf%tbsys%tbmodel%cutoff)

  call system_timer("calculate_forces_GF_calc_mod_dm")
  call calc_mod_dm_from_Gs(this%gf, w_e, w_n)
  call system_timer("calculate_forces_GF_calc_mod_dm")

  call system_timer("calculate_forces_GF_Tr")
  forces = 0.0_dp
  allocate(od_mask(this%at%N))
  allocate(d_mask(this%at%N))

  call Setup_deriv_matrices(this%gf%tbsys, this%at)
  do i=1, this%at%N
    call fill_dmatrices(this%gf%tbsys, this%at, i, diag_mask=d_mask, offdiag_mask=od_mask)

    forces(1,i) = forces(1,i) - local_ksum(this%gf%tbsys%kpoints, real(TraceMult(this%gf%mod_dm_H, this%gf%tbsys%dH(1), &
      a_H=.false., b_H=.true., diag_mask=d_mask, offdiag_mask=od_mask)))
    forces(2,i) = forces(2,i) - local_ksum(this%gf%tbsys%kpoints, real(TraceMult(this%gf%mod_dm_H, this%gf%tbsys%dH(2), &
      a_H=.false., b_H=.true., diag_mask=d_mask, offdiag_mask=od_mask)))
    forces(3,i) = forces(3,i) - local_ksum(this%gf%tbsys%kpoints, real(TraceMult(this%gf%mod_dm_H, this%gf%tbsys%dH(3), &
      a_H=.false., b_H=.true., diag_mask=d_mask, offdiag_mask=od_mask)))
    if (.not. this%gf%tbsys%tbmodel%is_orthogonal) then
      forces(1,i) = forces(1,i) - local_ksum(this%gf%tbsys%kpoints, real(TraceMult(this%gf%mod_dm_S, this%gf%tbsys%dS(1), &
	a_H=.false., b_H=.true., diag_mask=d_mask, offdiag_mask=od_mask)))
      forces(2,i) = forces(2,i) - local_ksum(this%gf%tbsys%kpoints, real(TraceMult(this%gf%mod_dm_S, this%gf%tbsys%dS(2), &
	a_H=.false., b_H=.true., diag_mask=d_mask, offdiag_mask=od_mask)))
      forces(3,i) = forces(3,i) - local_ksum(this%gf%tbsys%kpoints, real(TraceMult(this%gf%mod_dm_S, this%gf%tbsys%dS(3), &
	a_H=.false., b_H=.true., diag_mask=d_mask, offdiag_mask=od_mask)))
    endif
  end do

  deallocate(od_mask)
  deallocate(d_mask)

  call system_timer("calculate_forces_GF_Tr")

  call system_timer("calculate_forces_GF_accum")
  call ksum_distrib_inplace(this%gf%tbsys%kpoints, forces)

  call system_timer("calculate_forces_GF_scf")
  forces = forces + scf_f_correction_GF(this%gf, this%at, w_e, w_n)
  call system_timer("calculate_forces_GF_scf")

   do i=1, this%at%N
     if (associated(w_e)) then
       forces = forces + w_e(i) * get_local_rep_E_force(this%gf%tbsys%tbmodel, this%at, i)
     else
       forces = forces + get_local_rep_E_force(this%gf%tbsys%tbmodel, this%at, i)
     endif
   end do
  call system_timer("calculate_forces_GF_accum")
  call system_timer("calculate_forces_GF")

end function calculate_forces_GF

function scf_f_correction_GF(gf, at, w_e, w_n) result(forces)
  type(GreensFunctions), intent(inout) :: gf
  type(Atoms), intent(in) :: at
  real(dp), intent(in), pointer :: w_e(:), w_n(:)
  real(dp) :: forces(3,gf%tbsys%N_atoms)

  real(dp), allocatable :: n(:), dn_dr_mat(:,:,:), dgN_dr_vec(:,:)
  real(dp), allocatable :: dlpot(:)
  real(dp) :: dglobal_pot
  real(dp), allocatable :: ww_e(:)
  integer N_atoms

  type(TBMatrix) :: sp_Stwid, sp_S
  logical, allocatable :: d_mask(:), od_mask(:)

  integer i, j

  if (.not. at%use_uniform_cutoff .or. at%cutoff < gf%tbsys%tbmodel%cutoff) &
    call print ("WARNING: Called scf_f_correction_GF with messed up cutoff in atoms: uniform" // at%use_uniform_cutoff &
      // " " // at%cutoff // " < " // gf%tbsys%tbmodel%cutoff)

  N_atoms = gf%tbsys%N_atoms

  allocate(d_mask(N_atoms))
  allocate(od_mask(N_atoms))

  forces = 0.0_dp
  if (.not. gf%tbsys%scf%active) then
    return
  endif

  if (size(gf%tbsys%scf%terms) /= 1) then
    call system_abort("GreensFunc_scf_f_correction_GF not yet implemented for more than 1 type of SCF defined")
  endif

  if (gf%tbsys%scf%terms(1)%type == SCF_LCN .or. gf%tbsys%scf%terms(1)%type == SCF_GCN) then
    call system_abort("GreensFunc_scf_f_correction_GF not yet implemented for local or global charge neutrality")
  endif

  allocate(dlpot(gf%tbsys%N))

  call Initialise(sp_Stwid, at, gf%tbsys%first_orb_of_atom, gf%tbsys%n_matrices, gf%tbsys%complex_matrices)
  call Initialise(sp_S, at, gf%tbsys%first_orb_of_atom, gf%tbsys%n_matrices, gf%tbsys%complex_matrices)

  call fill_these_matrices(gf%tbsys, at, do_S = .true., S = sp_S)

  allocate(n(gf%tbsys%N))
  n = local_ksum(gf%tbsys%kpoints, real(partial_TraceMult(gf%dm, gf%tbsys%S, a_H=.true., b_H=.false.)))
  call ksum_distrib_inplace(gf%tbsys%kpoints, n)

  if (gf%tbsys%scf%terms(1)%type == SCF_LOCAL_U .or. &
      gf%tbsys%scf%terms(1)%type == SCF_NONLOCAL_U_DFTB .or. &
      gf%tbsys%scf%terms(1)%type == SCF_NONLOCAL_U_NRL_TB) then

    allocate(dn_dr_mat(3,N_atoms,gf%tbsys%N))
    call ALLOC_TRACE("scf_f_correction_GF dn_dr_mat", size(dn_dr_mat)*REAL_SIZE)

    ! calculate derivative of dn w.r.t r
    call calc_dn_dr_mat(gf, at, dn_dr_mat)

    ! calculate correction to forces using dn_dr_mat

    ! correct forces by derivative of shift of Hamiltonian
    do i=1, N_atoms
      call fill_dmatrices(gf%tbsys, at, i, diag_mask=d_mask, offdiag_mask=od_mask)
      do j=1, 3
	call add_term_d2SCFE_dn2_times_vec(gf%tbsys%scf%terms(1), gf%tbsys, dn_dr_mat(j,i,:), dlpot)

	call multDiagRL(sp_Stwid, sp_S, dlpot)
	forces(j,i) = forces(j,i) + local_ksum(gf%tbsys%kpoints, real(TraceMult(gf%mod_dm_H, sp_Stwid, a_H=.false., b_H=.true.)))

        ! Tr [ mod_dm 0.5 (lpot dS lpot) ] was necessary in SIMPLE (C+) version, since dH_ij/dr didn't
        ! include 0.5 ( lpot(i) + lpot(j) dS_ij/dr.  Here it does include it, so this
        ! isn't necessary

      end do ! j=1, 3
    end do ! i=1, N_atoms
    call ksum_distrib_inplace(gf%tbsys%kpoints, forces)

    call add_term_dscf_e_correction_dn(gf%tbsys%scf%terms(1), gf%tbsys, dlpot)

    allocate(ww_e(gf%tbsys%N))
    if (associated(w_e)) then
      ww_e = atom_orbital_spread(gf%tbsys, w_e)
    else
      ww_e = 1.0_dp
    endif

    ! correct forces by derivative of double counting terms
    do i=1, N_atoms
      do j=1, 3
	forces(j,i) = forces(j,i) - sum(ww_e(:)*dlpot(:)*dn_dr_mat(j,i,:))
      end do
    end do

    call ALLOC_TRACE("scf_f_correction_GF dn_dr_mat", size(dn_dr_mat)*REAL_SIZE)
    deallocate(dn_dr_mat)

  else ! global U

    allocate(dgN_dr_vec(3,N_atoms))

    if (.not. associated(w_n)) &
      call system_abort("called scf_f_correction_GF with GLOBAL_U but no global_at_weight")

    call calc_dgN_dr_vec(gf, at, dgN_dr_vec, w_n, sp_S, sp_Stwid)

    do i=1, at%N
      call fill_dmatrices(gf%tbsys, at, i, diag_mask=d_mask, offdiag_mask=od_mask)
      do j=1, 3
	call add_term_d2SCFE_dgNdn(gf%tbsys%scf%terms(1), gf%tbsys, w_n, dlpot)
	dlpot = dlpot*dgN_dr_vec(j,i)

	call multDiagRL(sp_Stwid, sp_S, dlpot)
	forces(j,i) = forces(j,i) + local_ksum(gf%tbsys%kpoints, real(TraceMult(gf%mod_dm_H, sp_Stwid, a_H=.false., b_H=.true.)))

        ! Tr [ mod_dm 0.5 (lpot dS lpot) ] was necessary in SIMPLE (C+) version, since dH_ij/dr didn't
        ! include 0.5 ( lpot(i) + lpot(j) dS_ij/dr.  Here we do, so this
        ! isn't necessary

      end do ! j=1, 3
    end do ! i=1, N_atoms
    call ksum_distrib_inplace(gf%tbsys%kpoints, forces)

    call add_term_dscf_e_correction_dgN(gf%tbsys%scf%terms(1), dglobal_pot)

    ! correct forces by derivative of double counting terms
    forces(:,:) = forces(:,:) - dglobal_pot*dgN_dr_vec(:,:)

    deallocate(dgN_dr_vec)
    
  endif ! local or global u

  call finalise(sp_Stwid)
  call finalise(sp_S)

  deallocate(d_mask)
  deallocate(od_mask)

end function scf_f_correction_GF

subroutine calc_dn_dr_mat(gf, at, dn_dr_mat)
  type(GreensFunctions), intent(inout) :: gf
  type(Atoms), intent(in) :: at
  real(dp), intent(out) :: dn_dr_mat(:,:,:)

  type(TBMatrix) :: M2_raw
  type(MatrixD) :: M2, M2inv
  integer N_atoms, N_eff

  real(dp), allocatable :: dn0_dr_mat(:,:,:), t_dn0_dr_mat(:,:,:)
  real(dp), allocatable :: dn0_dr(:)
  type(TBMatrix) :: GS, G_dAdr, GU, GSU, SGS_T, GS_T
  complex(dp) :: a, az

  real(dp), allocatable :: U_spread(:), gamma_spread(:,:)
  integer ip, i, j, ii
  logical, allocatable :: d_mask(:), od_mask(:)

  if (.not. at%use_uniform_cutoff .or. at%cutoff < gf%tbsys%tbmodel%cutoff) &
    call print ("WARNING: Called calc_dn_dr_mat with messed up cutoff in atoms: uniform" // at%use_uniform_cutoff &
      // " " // at%cutoff // " < " // gf%tbsys%tbmodel%cutoff)

  N_eff = gf%tbsys%N

  allocate(d_mask(N_atoms))
  allocate(od_mask(N_atoms))

  if (allocated(gf%tbsys%scf%terms(1)%gamma)) then
    allocate(gamma_spread(gf%tbsys%N,gf%tbsys%N))
    call atom_orbital_spread_mat(gf%tbsys, gf%tbsys%scf%terms(1)%gamma, gamma_spread)
  endif
  if (allocated(gf%tbsys%scf%terms(1)%U)) then
    allocate(U_spread(gf%tbsys%N))
    U_spread = atom_orbital_spread(gf%tbsys, gf%tbsys%scf%terms(1)%U)
  endif

  allocate(dn0_dr_mat(3,N_atoms,N_eff))
  dn0_dr_mat = 0.0_dp

  call Initialise(M2_raw, gf%tbsys%S%N, gf%tbsys%S%n_matrices, .false.)
  call Initialise(M2, N_eff)
  call Initialise(M2inv, N_eff)
  call zero(M2_raw)

  do ip=1, gf%N_G
    if (ip == 1) then
      allocate(t_dn0_dr_mat(3,N_atoms,N_eff))
      call ALLOC_TRACE("calc_dn_dr_mat t_dn0_dr_mat",size(t_dn0_dr_mat)*REAL_SIZE)
      call Initialise(GS, gf%G(1))
      call Initialise(GS_T, gf%G(1))
      call Initialise(G_dAdr, gf%G(1))
      call Initialise(SGS_T, gf%G(1))
      call Initialise(GU, gf%G(1))
      call Initialise(GSU, gf%G(1))
      allocate(dn0_dr(N_eff))
    endif

    t_dn0_dr_mat = 0.0_dp

    if (.not. gf%tbsys%complex_matrices) then ! factor of 2.0, instead of doing G(z) and G(conjg(z)) separately
      a = 2.0_dp * gf%a(ip)
      az = 2.0_dp * gf%a(ip)*gf%z(ip)
    else
      a = gf%a(ip)
      az = gf%a(ip)*gf%z(ip)
    endif

    call matrix_product_sub(GS, gf%G(ip), gf%tbsys%S)

    do i=1, at%N
      ! this version of fill_dmatrices returns scf-local-potential-shifted
      ! dH (i.e. dH(i,j) + 0.5(lpot(i)+lpot(j))dS(i,j)
      ! SIMPLE (C+) version didn't give shifted dH, and that contribution to dn/dr
      ! was accounted for in "M1".  This version doesn't need that

      call fill_dmatrices(gf%tbsys, at, i, diag_mask=d_mask, offdiag_mask=od_mask)
      do j=1, 3
	! Tr [ rho dS/dr WI ]
        if (.not. gf%tbsys%tbmodel%is_orthogonal) then
          if (ip == 1 .and. gf%mpi_across_poles%my_proc == 0) then
	    dn0_dr = local_ksum(gf%tbsys%kpoints, &
	      real(partial_TraceMult(gf%dm, gf%tbsys%dS(j), a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
            call ksum_distrib_inplace(gf%tbsys%kpoints, dn0_dr)
            do ii=1, N_eff
              !NB sign determined empirically :(
              t_dn0_dr_mat(j,i,ii) = -dn0_dr(ii)
            end do
          end if
	end if

	dn0_dr = 0.0_dp

	! Tr [ drho0/dr S WI ] = Tr [ sum_i ( a_i z_i G_i dS/dr G_i - a_i G_i dH0/dr G_i ) S WI =
        !    Tr [ sum_i a_i z_i G_i dS/dr G_i S WI - sum_i a_i G_i dH0/dr G_i S WI 
        if (.not. gf%tbsys%tbmodel%is_orthogonal) then
          call matrix_product_sub(G_dAdr,  gf%G(ip), gf%tbsys%dS(j), diag_mask=d_mask, offdiag_mask=od_mask)

	  dn0_dr = dn0_dr + real(az*local_ksum(gf%tbsys%kpoints, partial_TraceMult(G_dAdr, GS)))
        endif

	call matrix_product_sub(G_dAdr, gf%G(ip), gf%tbsys%dH(j), diag_mask=d_mask, offdiag_mask=od_mask)

	dn0_dr = dn0_dr - real(a*local_ksum(gf%tbsys%kpoints, partial_TraceMult(G_dAdr, GS)))

	call ksum_distrib_inplace(gf%tbsys%kpoints, dn0_dr)

        ! factor of 2 folded into a and az (or explicit evaluation of G(conjg(z)))
	t_dn0_dr_mat(j,i,:) = t_dn0_dr_mat(j,i,:) - dn0_dr(:)

      end do ! j
    end do ! i

    call transpose_sub(GS_T, GS)

    ! (SGS)^T = (S GS)^T = (GS)^T S^T
    call matrix_product_sub(SGS_T, GS_T, gf%tbsys%S, A_transpose = .false., B_transpose = .true.)

    if (allocated(gamma_spread)) then
      call matrix_product_sub(GU, gf%G(ip), gamma_spread)
    else
      call multDiag(GU, gf%G(ip), U_spread)
    endif

    if (allocated(gamma_spread)) then
      call matrix_product_sub(GSU, GS, gamma_spread)
    else
      call multDiag(GSU, GS, U_spread)
    endif

    ! factor of 0.5 to compensate for factor of 2 in a (or explicit eval of G(conjg(z)) 
    ! OMP critical
    call accum_scaled_elem_product(GU, SGS_T, 0.5_dp*a, M2_raw)
    call accum_scaled_elem_product(GSU, GS_T, 0.5_dp*a, M2_raw)

    ! OMP critical
    dn0_dr_mat = dn0_dr_mat + t_dn0_dr_mat

    if (gf%tbsys%complex_matrices) then

      t_dn0_dr_mat = 0.0_dp

      a = conjg(gf%a(ip))
      az = conjg(gf%a(ip)*gf%z(ip))

      call matrix_product_sub(GS, gf%G_conjg(ip), gf%tbsys%S)

      do i=1, at%N
        ! this version of fill_dmatrices returns scf-local-potential-shifted
        ! dH (i.e. dH(i,j) + 0.5(lpot(i)+lpot(j))dS(i,j)
        ! SIMPLE (C+) version didn't give shifted dH, and that contribution to dn/dr
        ! was accounted for in "M1".  This version doesn't need that

        call fill_dmatrices(gf%tbsys, at, i, diag_mask=d_mask, offdiag_mask=od_mask)
        do j=1, 3

          dn0_dr = 0.0_dp

          ! Tr [ drho0/dr S WI ] = Tr [ sum_i ( a_i z_i G_i dS/dr G_i - a_i G_i dH0/dr G_i ) S WI =
          !    Tr [ sum_i a_i z_i G_i dS/dr G_i S WI - sum_i a_i G_i dH0/dr G_i S WI 
          if (.not. gf%tbsys%tbmodel%is_orthogonal) then
            call matrix_product_sub(G_dAdr,  gf%G_conjg(ip), gf%tbsys%dS(j), diag_mask=d_mask, offdiag_mask=od_mask)

	    dn0_dr = dn0_dr + real(az*local_ksum(gf%tbsys%kpoints, partial_TraceMult(G_dAdr, GS)))
          endif

          call matrix_product_sub(G_dAdr, gf%G_conjg(ip), gf%tbsys%dH(j), diag_mask=d_mask, offdiag_mask=od_mask)

	  dn0_dr = dn0_dr - real(a*local_ksum(gf%tbsys%kpoints, partial_TraceMult(G_dAdr, GS)))
          call ksum_distrib_inplace(gf%tbsys%kpoints, dn0_dr)

          ! factor of 2 folded into a and az (or explicit evaluation of G(conjg(z)))
          t_dn0_dr_mat(j,i,:) = t_dn0_dr_mat(j,i,:) - dn0_dr(:)

        end do ! j
      end do ! i

      call transpose_sub(GS_T, GS)

      ! (SGS)^T = S^T (GS)^T
      call matrix_product_sub(SGS_T, GS_T, gf%tbsys%S, A_transpose = .false., B_transpose = .true.)

      if (allocated(gamma_spread)) then
        call matrix_product_sub(GU, gf%G_conjg(ip), gamma_spread)
      else
        call multDiag(GU, gf%G_conjg(ip), U_spread)
      endif

      if (allocated(gamma_spread)) then
        call matrix_product_sub(GSU, GS, gamma_spread)
      else
        call multDiag(GSU, GS, U_spread)
      endif

      ! factor of 0.5 to compensate for factor of 2 in a (or explicit eval of G(conjg(z)) 
      ! OMP critical
      call accum_scaled_elem_product(GU, SGS_T, 0.5_dp*a, M2_raw)
      call accum_scaled_elem_product(GSU, GS_T, 0.5_dp*a, M2_raw)

      ! OMP critical
      dn0_dr_mat = dn0_dr_mat + t_dn0_dr_mat

    endif ! non-gamma
  end do ! ip = 1..N_G

  call ksum_atom_orbital_sum_mat(gf%tbsys, M2_raw, M2)

  call Gsum_distrib_inplace(gf, dn0_dr_mat)
  call Gsum_distrib_inplace(gf, M2)

  ! dn[I]/dr = dn0[I]/dr + sum_M M2[I][M] dn[M]/dr
  ! x = v0 + M.x
  ! x - M.x = v0
  ! (I-M) . x  = v0
  ! x = (I-M)^-1 . v0

  !NB sign determined empirically
  !NB call scale(M2, -1.0_dp)
  call add_identity(M2)
  call inverse(M2, M2inv, .false.)

  do i=1, at%N
    do j=1, 3
      dn_dr_mat(j,i,:) = M2inv%data .mult. dn0_dr_mat(j,i,:)
    end do
  end do

  call finalise(M2_raw)
  call finalise(M2)
  call finalise(M2inv)

  call Finalise(GS)
  call Finalise(GS_T)
  call Finalise(G_dAdr)
  call Finalise(SGS_T)
  call Finalise(GU)
  call Finalise(GSU)

 call DEALLOC_TRACE("calc_dn_dr_mat gamma_spread", size(gamma_spread)*REAL_SIZE)
  deallocate(gamma_spread)
 call DEALLOC_TRACE("calc_dn_dr_mat U_spread", size(U_spread)*REAL_SIZE)
  deallocate(U_spread)
 call DEALLOC_TRACE("calc_dn_dr_mat dn0_dr_mat", size(dn0_dr_mat)*REAL_SIZE)
  deallocate(dn0_dr_mat)
 call DEALLOC_TRACE("calc_dn_dr_mat t_dn_dr_mat", size(t_dn0_dr_mat)*REAL_SIZE)
  deallocate(t_dn0_dr_mat)
 call DEALLOC_TRACE("calc_dn_dr_mat dn0_dr", size(dn0_dr)*REAL_SIZE)
  deallocate(dn0_dr)

  deallocate(d_mask)
  deallocate(od_mask)

end subroutine calc_dn_dr_mat

subroutine calc_dgN_dr_vec(gf, at, dgN_dr_vec, w_n, sp_S, sp_Stwid)
  type(GreensFunctions), intent(inout) :: gf
  type(Atoms), intent(in) :: at
  real(dp), intent(out) :: dgN_dr_vec(:,:)
  real(dp), intent(in) :: w_n(:)
  type(TBMatrix), intent(in) :: sp_S
  type(TBmatrix), intent(inout) :: sp_Stwid

  real(dp), allocatable :: ww_n(:)
  real(dp) dgN_dr_denom
  complex(dp) :: a, az

  logical, allocatable :: d_mask(:), od_mask(:)

  integer ip, i, j

  type(TBMatrix) :: GWN, SG, aGWNSG, azGWNSG, GWNSG

  if (.not. at%use_uniform_cutoff .or. at%cutoff < gf%tbsys%tbmodel%cutoff) &
    call print ("WARNING: Called calc_dgN_dr_vec with messed up cutoff in atoms: uniform" // at%use_uniform_cutoff &
      // " " // at%cutoff // " < " // gf%tbsys%tbmodel%cutoff)

  allocate(ww_n(gf%tbsys%N))
  ww_n = atom_orbital_spread(gf%tbsys, w_n)

  allocate(d_mask(at%N))
  allocate(od_mask(at%N))

  call multDiagRL(sp_Stwid, sp_S, gf%tbsys%scf%global_U*ww_n)

  do ip=1, gf%N_G
    if (ip == 1) then
      call Initialise(SG, gf%G(1))
      call Initialise(GWN, gf%G(1))
      call Initialise(aGWNSG, gf%G(1))
      call Initialise(azGWNSG, gf%G(1))
      call Initialise(GWNSG, gf%G(1))
      call zero(aGWNSG)
      call zero(azGWNSG)
    endif


    if (.not. gf%tbsys%complex_matrices) then ! factor of 2.0, instead of doing G(z) and G(conjg(z)) separately
      a = 2.0_dp*gf%a(ip)
      az = 2.0_dp*gf%a(ip)*gf%z(ip)
    else
      a = gf%a(ip)
      az = gf%a(ip)*gf%z(ip)
    endif

    call multDiag(GWN, gf%G(ip), ww_n)
    call matrix_product_sub(SG, gf%tbsys%S, gf%G(ip))
    call matrix_product_sub(GWNSG, GWN, SG)

    ! omp critical
    call scaled_accum(aGWNSG, a, GWNSG)
    call scaled_accum(azGWNSG, az, GWNSG)

    if (gf%tbsys%complex_matrices) then
      a = conjg(gf%a(ip))
      az = conjg(gf%a(ip)*gf%z(ip))

      call multDiag(GWN, gf%G_conjg(ip), ww_n)
      call matrix_product_sub(SG, gf%tbsys%S, gf%G_conjg(ip))
      call matrix_product_sub(GWNSG, GWN, SG)

      ! omp critical
      call scaled_accum(aGWNSG, a, GWNSG)
      call scaled_accum(azGWNSG, az, GWNSG)

    endif
  end do

  call Finalise(SG)
  call Finalise(GWN)
  call Finalise(GWNSG)

  call Gsum_distrib_inplace(gf, aGWNSG)
  call Gsum_distrib_inplace(gf, azGWNSG)

  !NB sign set empirically
  ! factor of 2 shifted into aGWNSG
  dgN_dr_denom = 1.0_dp + ksum_distrib(gf%tbsys%kpoints, real(local_ksum(gf%tbsys%kpoints, TraceMult(sp_Stwid, aGWNSG, &
				       a_H=.true., b_H=.false.))))

  dgN_dr_vec = 0.0_dp

  do i=1, at%N
    call fill_dmatrices(gf%tbsys, at, i, diag_mask=d_mask, offdiag_mask=od_mask)
    do j=1, 3
      if (.not. gf%tbsys%tbmodel%is_orthogonal) then
	! Tr [ dS/dr rho wN ]
	! NB sign set empirically
	dgN_dr_vec(j,i) = dgN_dr_vec(j,i) - sum(ww_n * &
					    local_ksum(gf%tbsys%kpoints, real(partial_TraceMult(gf%tbsys%dS(j), gf%dm, &
					    a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask))))
	! factor of 2 shifted into aGWNSG and azGWNSG
	dgN_dr_vec(j,i) = dgN_dr_vec(j,i) - local_ksum(gf%tbsys%kpoints, real(TraceMult(gf%tbsys%dS(j), azGWNSG, &
					               a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))
      end if
      dgN_dr_vec(j,i) = dgN_dr_vec(j,i) + &
        local_ksum(gf%tbsys%kpoints, real(TraceMult(gf%tbsys%dH(j), aGWNSG, &
					  a_H=.true., b_H=.false., diag_mask=d_mask, offdiag_mask=od_mask)))

      ! NB SIMPLE (C+) implementation needs to add something because its dH_ij/dr doesn't
      !  include the 0.5 (lpot(i) + lpot(j) dS_ij/dr part
    end do
  end do

  call ksum_distrib_inplace(gf%tbsys%kpoints, dgN_dr_vec)

  dgN_dr_vec = dgN_dr_vec / dgN_dr_denom

  call Finalise(aGWNSG)
  call Finalise(azGWNSG)

  deallocate(d_mask)
  deallocate(od_mask)

end subroutine calc_dgN_dr_vec


subroutine calc_dm_from_evecs(this, for_forces)
  type(TB_type), intent(inout) :: this
  logical, intent(in), optional :: for_forces

  logical do_forces

  do_forces = .false.
  if (present(for_forces)) do_forces = for_forces

  call realloc_match_tbsys(this%tbsys, this%dm)
  call realloc_match_tbsys(this%tbsys, this%scaled_evecs)
  call Zero(this%dm)
  if (do_forces .and. .not. this%tbsys%tbmodel%is_orthogonal) then
    call realloc_match_tbsys(this%tbsys, this%Hdm)
    call Zero(this%Hdm)
  endif

  if (do_forces) then
    call multDiag(this%scaled_evecs, this%evecs, this%F_fillings)
    call matrix_product_sub(this%dm, this%scaled_evecs, this%evecs, b_conjugate = .true.)
    if (.not. this%tbsys%tbmodel%is_orthogonal) then
      call multDiag(this%scaled_evecs, this%evecs, this%eval_F_fillings)
      call matrix_product_sub(this%Hdm, this%scaled_evecs, this%evecs, b_conjugate = .true.)
    endif
  else
    call multDiag(this%scaled_evecs, this%evecs, this%E_fillings)
    call matrix_product_sub(this%dm, this%scaled_evecs, this%evecs, b_conjugate = .true.)
  endif

end subroutine calc_dm_from_evecs

subroutine TB_calc_E_fillings(this, use_fermi_E, fermi_E, AF, w_n)
  type(TB_type), intent(inout) :: this
  logical, intent(in), optional :: use_fermi_E
  real(dp), intent(in), optional :: fermi_E
  type(ApproxFermi), intent(inout), optional :: AF
  real(dp), intent(in), pointer :: w_n(:)

  logical find_new_fermi_E
  real(dp) :: degeneracy

  if (present(use_fermi_E) .and. .not. present(AF)) then
    if (.not. has_fermi_E(this%fermi_E, this%tbsys%tbmodel, fermi_E, this%calc_args_str)) &
      call system_abort("called calc_E_fillings with use_fermi_E, but no Fermi_E or AF")
  endif

  find_new_fermi_E = .true.
  if (present(use_fermi_E)) find_new_fermi_E = .not. use_fermi_E

  if (present(AF)) call print("calc_E_fillings using approx fermi function", PRINT_VERBOSE)

  if (find_new_fermi_E) then
    call print("calc_E_fillings finding new fermi level", PRINT_VERBOSE)
  else
    call print("calc_E_fillings using current fermi level", PRINT_VERBOSE)
  endif

  call realloc_match_tbsys(this%tbsys, this%E_fillings)

  if (find_new_fermi_E) then
    ! not so useful for Fermi E to be set based on w_n.
    ! mostly w_n applies to SCF_GCN and SCF_LOBAL_U.
    ! can be restored if necessary, but should probably have a different
    ! name for the property than weight_n
    call find_fermi_E(this, AF)
    call print ("TB_calc_E_fillings got new Fermi_E " // this%Fermi_E, PRINT_VERBOSE)
  else
    if (present(fermi_E)) then
      this%fermi_E = fermi_E
    endif
    if (this%tbsys%noncollinear) then
      degeneracy = 1.0_dp
    else
      degeneracy = 2.0_dp
    endif
    if (present(AF)) then
      call calc_fermi_factors(this%E_fillings, this%evals,  AF = AF, degeneracy=degeneracy)
    else
      call calc_fermi_factors(this%E_fillings, this%evals, fermi_E = this%fermi_E, fermi_T = this%fermi_T, degeneracy=degeneracy)
    endif
  endif

end subroutine TB_calc_E_fillings

subroutine TB_calc_F_fillings(this, need_eval_F_fillings, AF)
  type(TB_type), intent(inout) :: this
  logical, intent(in), optional :: need_eval_F_fillings
  type(ApproxFermi), intent(inout), optional :: AF

  call realloc_match_tbsys(this%tbsys, this%F_fillings)

  if (present(AF)) then
    call calc_mod_fermi_factors(this, this%E_fillings, this%F_fillings, AF = AF)
  else
    call calc_mod_fermi_factors(this, this%E_fillings, this%F_fillings, fermi_E = this%fermi_E, fermi_T = this%fermi_T)
  endif

  if (present(need_eval_F_fillings)) then
    if (need_eval_F_fillings) then
      call realloc_match_tbsys(this%tbsys, this%eval_F_fillings)
      this%eval_F_fillings%data_d = this%F_fillings%data_d*this%evals%data_d
    endif
  endif

end subroutine TB_calc_F_fillings

subroutine TB_find_fermi_E(this, AF, w_n)
  type(TB_type), intent(inout) :: this
  type(ApproxFermi), intent(inout), optional :: AF
  real(dp), intent(in), pointer, optional :: w_n(:)

  integer iter
  integer, parameter :: max_iter = 200
  real(dp) :: e_min, e_max, e_try
  real(dp) :: N_try
  logical :: have_w_n

  real(dp), allocatable :: local_N(:)
  real(dp) :: N_e
  real(dp) :: EPS
  real(dp) :: degeneracy

  call print("called find_fermi_E fermi_T " // this%fermi_T // " " // (this%fermi_T < 1e-8_dp), PRINT_ANAL)

  have_w_n = .false.
  if (present(w_n)) then
    if (associated(w_n)) have_w_n = .true.
  endif

  N_e = n_elec(this%tbsys, this%at, w_n)

  if (present(AF)) then
    if (AF%n_poles == 0 .or. (AF%band_width .feq. 0.0_dp)) then
      call system_abort("called find_fermi_E with AF present but n_poles of band_width == 0")
    endif
  endif

  if (this%fermi_T < 1e-8_dp) then
    call system_abort ("find_fermi_E can't handle fermi_T < 1e-8")
  else
    e_min = min(this%tbsys%kpoints, this%evals%data_d(1,:)) - 10.0_dp
    e_max = max(this%tbsys%kpoints, this%evals%data_d(this%evals%N,:)) + 10.0_dp

    call Print("find_Fermi_E N_e " // N_e // " e_min " // e_min // " e_max " // e_max, PRINT_ANAL)

    if (have_w_n) allocate(local_N(this%tbsys%N_atoms))

    iter = 0
    N_try = N_e + 1
    do while (abs(N_try-N_e) > this%fermi_E_precision)
      e_try = (e_min + e_max) / 2.0_dp
      if (this%tbsys%noncollinear) then
	degeneracy = 1.0_dp
      else
	degeneracy = 2.0_dp
      endif
      if (present(AF)) then
	AF%z = AF%z - AF%fermi_E + e_try
	AF%fermi_E = e_try
	call calc_fermi_factors(this%E_fillings, this%evals, AF = AF, degeneracy=degeneracy)
      else
	call calc_fermi_factors(this%E_fillings, this%evals, fermi_E = e_try, fermi_T = this%fermi_T, degeneracy=degeneracy)
      endif
      if (have_w_n) then
	! could be done more efficiency, by precomputing for each state the
	! occupation on each atom
	call calc_dm_from_evecs(this, .false.)
	call calc_local_atomic_num(this, local_N)
	N_try = sum(w_n*local_N)
      else
	N_try = ksum_dup(this%tbsys%kpoints, sum(this%E_fillings%data_d,dim=1))
      endif
      call Print("find_Fermi_E e_try n_try " // (iter+1) // " " // e_try // " " // n_try, PRINT_ANAL)
      if (N_try < N_e) then
	e_min = e_try
      else
	e_max = e_try
      endif
      iter = iter + 1
      if (iter .gt. max_iter) then
	call print("N_e " // N_e // " N_try " // N_try, PRINT_ALWAYS)
	call print("e_min " // e_min // " e_max " // e_max // " e_try " // e_try, PRINT_ALWAYS)
	call print("this%evals", PRINT_ALWAYS)
	call verbosity_push(PRINT_NORMAL)
	call print(this%evals)
	call verbosity_pop()
	call system_abort("Ran out of iterations in find_fermi_E")
      endif
      if ((abs(N_try-N_e) > this%fermi_E_precision) .and. (e_min == e_max)) then
	call print("N_e " // N_e // " N_try " // N_try, PRINT_ALWAYS)
	call print("e_min " // e_min // " e_max " // e_max // " e_try " // e_try, PRINT_ALWAYS)
	call print("this%evals", PRINT_ALWAYS)
	call verbosity_push(PRINT_NORMAL)
	call print(this%evals)
	call verbosity_pop()
	call system_abort("Ran out of precision in find_fermi_E")
      endif
    end do

    this%fermi_E = e_try
  endif

  EPS=1e-4_dp
  this%homo_e = maxval(this%evals%data_d, this%evals%data_d <= this%fermi_E+EPS)
  this%lumo_e = minval(this%evals%data_d, this%evals%data_d >= this%fermi_E-EPS)
  this%homo_e = max(this%tbsys%kpoints, this%homo_e)
  this%lumo_e = min(this%tbsys%kpoints, this%lumo_e)

  if (allocated(local_N)) deallocate(local_N)

end subroutine TB_find_fermi_E

subroutine calc_fermi_factors(fermi_factors, evals, fermi_E, fermi_T, AF, degeneracy)
  type(TBVector), intent(inout) :: fermi_factors
  type(TBVector), intent(in) :: evals
  real(dp), intent(in), optional :: fermi_E, fermi_T
  type(ApproxFermi), intent(in), optional :: AF
  real(dp), intent(in), optional :: degeneracy

  real(dp) :: u_degeneracy

  u_degeneracy = optional_default(2.0_dp, degeneracy)

  if (present(AF)) then
    fermi_factors%data_d(:,:) = u_degeneracy*approx_f_fermi(AF, evals%data_d(:,:))
  else if (present(Fermi_E) .and. present(Fermi_T)) then
    fermi_factors%data_d(:,:) = u_degeneracy*f_fermi(fermi_E, fermi_T, evals%data_d(:,:))
  else
    call system_abort("Called calc_fermi_factors without sufficient inputs for exact or approx fermi function")
  endif
end subroutine calc_fermi_factors

subroutine calc_fermi_derivs(fermi_factors, evals, fermi_E, fermi_T, AF, degeneracy)
  type(TBVector), intent(inout) :: fermi_factors
  type(TBVector), intent(in) :: evals
  real(dp), intent(in), optional :: fermi_E, fermi_T
  type(ApproxFermi), intent(in), optional :: AF
  real(dp), intent(in), optional :: degeneracy

  real(dp) :: u_degeneracy

  u_degeneracy = optional_default(2.0_dp, degeneracy)

  if (present(AF) .and. (present(fermi_E) .or. present(fermi_T))) then
    call system_abort("Called calc_fermi_factors with both ApproxFermi and Fermi_E or Fermi_T")
  endif

  if (present(AF)) then
    fermi_factors%data_d(:,:) = u_degeneracy*approx_f_fermi_deriv(AF, evals%data_d(:,:))
  else if (present(Fermi_E) .and. present(Fermi_T)) then
    fermi_factors%data_d(:,:) = u_degeneracy*f_fermi_deriv(fermi_E, fermi_T, evals%data_d(:,:))
  else
    call system_abort("Called calc_fermi_factors without sufficient inputs for exact or approx fermi function")
  endif

end subroutine calc_fermi_derivs

subroutine calc_mod_fermi_factors(this, fermi_factors, mod_fermi_factors, fermi_E, fermi_T, AF)
  type(TB_type), intent(inout) :: this
  type(TBVector), intent(in) :: fermi_factors
  type(TBVector), intent(inout) :: mod_fermi_factors
  real(dp), intent(in), optional :: Fermi_E, Fermi_T
  type(ApproxFermi), intent(in), optional :: AF

  real(dp) numerator, denominator, theta
  type(TBVector) :: fermi_derivs
  real(dp) :: degeneracy

  if (present(AF) .and. (present(fermi_E) .or. present(fermi_T))) then
    call system_abort("Called calc_mod_fermi_factors with both ApproxFermi and Fermi_E or Fermi_T")
  endif

  call Initialise(fermi_derivs, fermi_factors%N, fermi_factors%n_vectors)

  if (this%tbsys%noncollinear) then
    degeneracy = 1.0_dp
  else
    degeneracy = 2.0_dp
  endif

  call calc_fermi_derivs(fermi_derivs, this%evals, fermi_E, fermi_T, AF, degeneracy)

  denominator = ksum_dup(this%tbsys%kpoints, sum(fermi_derivs%data_d,dim=1))

  if (denominator .feq. 0.0_dp) then
    mod_fermi_factors%data_d = fermi_factors%data_d
  else
    numerator = ksum_dup(this%tbsys%kpoints, sum(this%evals%data_d*fermi_derivs%data_d,dim=1))
    theta = numerator/denominator

!   constant mu
    ! mod_fermi_factors%data_d = fermi_factors%data_d + &
    ! fermi_derivs%data_d*(this%evals%data_d)

!   constant N
    mod_fermi_factors%data_d = fermi_factors%data_d + &
      fermi_derivs%data_d*(this%evals%data_d-theta)

  endif

  call Finalise(fermi_derivs)
end subroutine calc_mod_fermi_factors

function TB_evals(this)
  type(TB_type) :: this
  real(dp) :: TB_evals(this%tbsys%N,this%tbsys%kpoints%g_N)

  if(allocated(this%evals%data_d)) then
     call collect(this%tbsys%kpoints, this%evals%data_d, TB_evals)
  else
     call system_abort("TB_evals() call with unallocated evals%data_d array")
  end if

end function TB_evals

! calculation optical absorption with formula from
! Snyder and Rotkin, Small v. 4, p. 1284-1286
subroutine absorption(this, polarization, freqs, gamma, a)
  type(TB_type), intent(inout) :: this
  complex(dp), intent(in) :: polarization(3)
  real(dp), intent(in) :: freqs(:), gamma
  real(dp), intent(out) :: a(:)

  complex(dp) :: polarization_hat(3)
  type(TBMatrix) :: dipole_evecs(3)
  integer :: i, f, ik, i_freq
  real(dp), allocatable :: absorption_contrib(:)

  if (this%evecs%N == 0 .or. this%evecs%n_matrices == 0) &
    call system_abort("Absorption needs eigenvectors to be already calculated")

  call dipole_matrix(this, dipole_evecs)

  call verbosity_push_decrement()
    call print("dipole_evecs(1)")
    call print(dipole_evecs(1))
    call print("dipole_evecs(2)")
    call print(dipole_evecs(2))
    call print("dipole_evecs(3)")
    call print(dipole_evecs(3))
  call verbosity_pop()

  polarization_hat = polarization / sqrt(sum(abs(polarization)**2))

  allocate(absorption_contrib(this%evecs%n_matrices))

  do i_freq=1, size(freqs)
    absorption_contrib = 0.0_dp
    do ik=1, this%evecs%n_matrices
      do i=1, this%evecs%N
	do f=i+1, this%evecs%N
	  if ((this%E_fillings%data_d(i,ik) .fne. 0.0_dp) .and. (this%E_fillings%data_d(f,ik) .fne. 1.0_dp)) then
	    if (dipole_evecs(1)%is_complex) then
	      absorption_contrib(ik) = absorption_contrib(ik) + &
			  ( (abs(dipole_evecs(1)%data_z(ik)%data(i,f)*polarization_hat(1) + &
				 dipole_evecs(2)%data_z(ik)%data(i,f)*polarization_hat(2) + &
				 dipole_evecs(3)%data_z(ik)%data(i,f)*polarization_hat(3))**2)/(freqs(i_freq)) ) * &
			  ( (this%E_fillings%data_d(i,ik)/2.0_dp*(1.0_dp-this%E_fillings%data_d(f,ik)/2.0_dp)) / &
			    ((this%evals%data_d(f,ik)-this%evals%data_d(i,ik)-freqs(i_freq))**2 + gamma**2) )
	    else
	      absorption_contrib(ik) = absorption_contrib(ik) + &
			  ( (abs(dipole_evecs(1)%data_d(ik)%data(i,f)*polarization_hat(1) + &
				 dipole_evecs(2)%data_d(ik)%data(i,f)*polarization_hat(2) + &
				 dipole_evecs(3)%data_d(ik)%data(i,f)*polarization_hat(3))**2)/(freqs(i_freq)) ) * &
			  ( (this%E_fillings%data_d(i,ik)/2.0_dp*(1.0_dp-this%E_fillings%data_d(f,ik)/2.0_dp)) / &
			    ((this%evals%data_d(f,ik)-this%evals%data_d(i,ik)-freqs(i_freq))**2 + gamma**2) )
	    endif
	  endif
	end do ! f
      end do ! i
    end do ! ik

    a(i_freq) = ksum_distrib(this%tbsys%kpoints, local_ksum(this%tbsys%kpoints, absorption_contrib))
  end do

  deallocate(absorption_contrib)

  call finalise(dipole_evecs(1))
  call finalise(dipole_evecs(2))
  call finalise(dipole_evecs(3))

end subroutine absorption

subroutine dipole_matrix(this, dipole_evecs)
  type(TB_type), intent(inout) :: this
  type(TBMatrix), intent(inout) :: dipole_evecs(3)

  type(TBMatrix) :: dipole_basis(3), t

  integer :: i

  do i=1, 3
    call realloc_match_tbsys(this%tbsys, dipole_basis(i))
    call realloc_match_tbsys(this%tbsys, dipole_evecs(i))
  end do
  call realloc_match_tbsys(this%tbsys, t)

  call fill_these_matrices(this%tbsys, this%at, do_dipole=.true., dipole=dipole_basis)

  call verbosity_push_decrement(PRINT_NERD)
    call print("dipole_basis(1)")
    call print(dipole_basis(1))
    call print("dipole_basis(2)")
    call print(dipole_basis(2))
    call print("dipole_basis(3)")
    call print(dipole_basis(3))
  call verbosity_pop()

  do i=1,3
    call matrix_product_sub(t, dipole_basis(i), this%evecs, a_conjugate=.false., b_conjugate = .false.)
    call matrix_product_sub(dipole_evecs(i), this%evecs, t, a_conjugate=.true., b_conjugate = .false.)
  end do

  do i=1, 3
    call finalise(dipole_basis(i))
  end do
  call finalise(t)

end subroutine dipole_matrix

end module TB_module
