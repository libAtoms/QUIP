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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
!X
!X TB_GreensFunctions module
!X
!% Contain Green's functions of a TB model, calculate them,
!% and compute derivatives of energies by summing over GFs
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
module TB_GreensFunctions_module

use libatoms_module

use Functions_module
use ApproxFermi_module
use Matrix_module
use TBMatrix_module
use TB_Kpoints_module
use TBSystem_module

implicit none
private

public :: GreensFunctions
type GreensFunctions
  integer :: N_G = 0
  complex(dp), allocatable :: a(:), z(:)

  type(TBSystem) :: tbsys
  type(TBMatrix), allocatable :: G(:), G_conjg(:)
  type(TBMatrix) :: dm, mod_dm_H, mod_dm_S
  type(MPI_Context) :: mpi_global, mpi_my_poles, mpi_across_poles
end type GreensFunctions

public :: Initialise
interface Initialise
  module procedure GreensFunctions_Initialise
end interface Initialise

public :: Finalise
interface Finalise
  module procedure GreensFunctions_Finalise
end interface Finalise

public :: Wipe
interface Wipe
  module procedure GreensFunctions_Wipe
end interface Wipe

public :: Print
interface Print
  module procedure GreensFunctions_Print
end interface Print

public :: Setup_system
interface Setup_system
  module procedure GreensFunctions_Setup_system
end interface Setup_system

interface init_mpi
  module procedure GreensFunctions_init_mpi
end interface init_mpi

interface end_mpi
  module procedure GreensFunctions_end_mpi
end interface end_mpi

public :: calc_Gs
interface calc_Gs
  module procedure GreensFunctions_calc_Gs
end interface calc_Gs

public :: calc_dm_from_Gs
interface calc_dm_from_Gs
  module procedure GreensFunctions_calc_dm_from_Gs
end interface calc_dm_from_Gs

public :: calc_mod_dm_from_Gs
interface calc_mod_dm_from_Gs
  module procedure GreensFunctions_calc_mod_dm_from_Gs
end interface calc_mod_dm_from_Gs

public :: Gsum_distrib_inplace
interface Gsum_distrib_inplace
  module procedure GreensFunctions_Gsum_distrib_inplace_tbm, GreensFunctions_Gsum_distrib_inplace_r2
  module procedure GreensFunctions_Gsum_distrib_inplace_r3, GreensFunctions_Gsum_distrib_inplace_matd
end interface Gsum_distrib_inplace

public :: Gsum_distrib
interface Gsum_distrib
  module procedure GreensFunctions_Gsum_distrib_d
end interface Gsum_distrib

contains

subroutine GreensFunctions_Finalise(this)
  type(GreensFunctions), intent(inout) :: this

  call end_mpi(this)

  call Wipe(this)
  call Finalise(this%dm)
  call Finalise(this%mod_dm_H)
  call Finalise(this%mod_dm_S)
  call Finalise(this%tbsys)
  call Finalise(this%mpi_global)
  call Finalise(this%mpi_my_poles)
  call Finalise(this%mpi_across_poles)
end subroutine GreensFunctions_Finalise

subroutine GreensFunctions_Wipe(this)
  type(GreensFunctions), intent(inout) :: this

  integer :: i

  call Wipe(this%tbsys)

  if (allocated(this%G)) then
    do i=1, size(this%G)
      call Finalise(this%G(i))
    end do
    deallocate(this%G)
  endif
  if (allocated(this%G_conjg)) then
    do i=1, size(this%G_conjg)
      call Finalise(this%G_conjg(i))
    end do
    deallocate(this%G_conjg)
  endif
  if (allocated(this%a)) deallocate(this%a)
  if (allocated(this%z)) deallocate(this%z)
  this%N_G = 0

  call finalise(this%dm)
  call finalise(this%mod_dm_H)
  call finalise(this%mod_dm_S)
end subroutine GreensFunctions_Wipe

subroutine GreensFunctions_Print(this,file)
  type(GreensFunctions),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  integer::i

  call Print('GreensFunctions : ', file=file)

  call Print ('GreensFunctions: N_G ' // this%N_G, file=file)
  if (this%N_G > 0) then
    do i=1, size(this%a)
      call Print ("GreensFunctions : a z " // i // " " // this%a(i) // " " // this%z(i), file=file)
    end do
    call verbosity_push_decrement()
    do i=1, size(this%G)
      call Print(this%G(i), file=file)
    end do
    do i=1, size(this%G_conjg)
      call Print(this%G_conjg(i), file=file)
    end do
    call verbosity_pop()
  endif
end subroutine GreensFunctions_Print

subroutine GreensFunctions_init_mpi(this, mpi_obj)
  type(GreensFunctions), intent(inout) :: this
  type(MPI_Context), intent(in) :: mpi_obj

  integer i
  integer local_N_G, n_procs_per_pole, t
  integer my_split_ident
  integer, allocatable :: local_index(:)
  complex(dp), allocatable :: t_a(:), t_z(:)
  logical :: no_pole_here

  this%mpi_global = mpi_obj

  if (this%N_G > this%mpi_global%n_procs) then
    n_procs_per_pole = 1
  else
    n_procs_per_pole = this%mpi_global%n_procs/this%N_G
  endif

  allocate(local_index(this%N_G))
  local_N_G = 0
  do i=0, this%N_G-1
    t = this%mpi_global%my_proc/n_procs_per_pole
    if (mod(i,this%mpi_global%n_procs) == t) then
      local_N_G = local_N_G + 1
      local_index(local_N_G) = i+1
    end if	
  end do

  if (local_N_G  == 0) then
    local_N_G = 1
    no_pole_here = .true.
    local_index(1) = this%N_G
  else
    no_pole_here = .false.
  endif

  if (this%mpi_global%active) then
    my_split_ident = this%mpi_global%my_proc/n_procs_per_pole
    call split_context(this%mpi_global, my_split_ident, this%mpi_my_poles)

    my_split_ident = this%mpi_my_poles%my_proc
    call split_context(this%mpi_global, my_split_ident, this%mpi_across_poles)
  endif

  allocate(t_a(local_N_G))
  allocate(t_z(local_N_G))
  t_a(1:local_N_G) = this%a(local_index(1:local_N_G))
  t_z(1:local_N_G) = this%z(local_index(1:local_N_G))

  deallocate(this%a)
  deallocate(this%z)
  this%N_G = local_N_G
  allocate(this%a(this%N_G))
  allocate(this%z(this%N_G))
  this%a = t_a
  this%z = t_z
  if (no_pole_here) this%a = 0.0_dp
  deallocate(t_a)
  deallocate(t_z)
end subroutine GreensFunctions_init_mpi

subroutine GreensFunctions_end_mpi(this)
  type(GreensFunctions), intent(inout) :: this

  if (this%mpi_global%active) then
    call free_context(this%mpi_my_poles)
    call free_context(this%mpi_across_poles)
  endif
end subroutine GreensFunctions_end_mpi


subroutine GreensFunctions_Initialise(this, z, a, tbsys, mpi_obj)
  type(GreensFunctions), intent(inout) :: this
  complex(dp) :: z(:), a(:)
  type(TBSystem), intent(in), optional :: tbsys
  type(MPI_Context), intent(in), optional :: mpi_obj


  if (size(z) /= size(a)) call system_abort("Called GreensFunctions_Initialise with mismatching z and a arrays")

  call Finalise(this)

  this%N_G = size(z)
  if (allocated(this%z)) deallocate(this%z)
  if (allocated(this%a)) deallocate(this%a)
  allocate(this%z(size(z)))
  allocate(this%a(size(a)))
  this%z = z
  this%a = a

  if (present(mpi_obj)) then
    call init_mpi(this, mpi_obj)
  else
    call init_mpi(this, tbsys%mpi_global)
  endif

  if (present(tbsys)) then
    call Setup_system(this, tbsys)
  endif

end subroutine GreensFunctions_Initialise

subroutine GreensFunctions_Setup_system(this, tbsys)
  type(GreensFunctions), intent(inout) :: this
  type(TBSystem), intent(in) :: tbsys

  integer i

  call Finalise(this%tbsys)

  if (allocated(this%G)) then
    do i=1, size(this%G)
      call Finalise(this%G(i))
    end do
    deallocate(this%G)
  endif
  if (allocated(this%G_conjg)) then
    do i=1, size(this%G_conjg)
      call Finalise(this%G_conjg(i))
    end do
    deallocate(this%G_conjg)
  endif

  call Initialise(this%tbsys, tbsys, mpi_obj=this%mpi_my_poles)
  this%tbsys%scf = tbsys%scf

  allocate(this%G(this%N_G))
  do i=1, this%N_G
    call Initialise(this%G(i), this%tbsys%N, this%tbsys%n_matrices, .true., this%tbsys%scalapack_my_matrices)
  end do
  if (this%tbsys%complex_matrices) then
    allocate(this%G_conjg(this%N_G))
    do i=1, this%N_G
      call Initialise(this%G_conjg(i), this%tbsys%N, this%tbsys%n_matrices, .true., this%tbsys%scalapack_my_matrices)
    end do
  endif

  call Initialise(this%dm, this%tbsys%N, this%tbsys%n_matrices, this%tbsys%complex_matrices, &
    scalapack_obj=this%tbsys%scalapack_my_matrices)
  call Initialise(this%mod_dm_H, this%tbsys%N, this%tbsys%n_matrices, this%tbsys%complex_matrices, &
    scalapack_obj=this%tbsys%scalapack_my_matrices)
  if (.not. this%tbsys%tbmodel%is_orthogonal) then
    call Initialise(this%mod_dm_S, this%tbsys%N, this%tbsys%n_matrices, this%tbsys%complex_matrices, &
      scalapack_obj=this%tbsys%scalapack_my_matrices)
  endif

end subroutine

subroutine GreensFunctions_calc_Gs(this, at, SelfEnergy)
  type(GreensFunctions), intent(inout) :: this
  type(Atoms), intent(inout) :: at
  type(TBMatrix), intent(in), optional :: SelfEnergy(:)

  integer i


  call fill_matrices(this%tbsys, at, .true.)

  if (present(SelfEnergy)) then
    if (this%N_G /= size(SelfEnergy)) &
      call system_abort("Called GreensFunctions_calc_Gs with SelfEnergy size mismatch " // size(SelfEnergy) // " " // this%N_G)
  endif

!$omp parallel do
  do i=1, this%N_G
    mainlog%mpi_all_inoutput_flag=.true.
    call print("GreensFunctions_calc_Gs doing i="//i, PRINT_NERD)
    mainlog%mpi_all_inoutput_flag=.false.
    call scaled_sum(this%G(i), this%z(i), this%tbsys%S, -1.0_dp, this%tbsys%H)
    if (present(SelfEnergy)) then
      call scaled_accum(this%G(i), cmplx(1.0_dp, 0.0_dp, dp), SelfEnergy(i))
    endif
    call inverse(this%G(i), positive = .false.)
    if (this%tbsys%complex_matrices) then
      call scaled_sum(this%G_conjg(i), conjg(this%z(i)), this%tbsys%S, -1.0_dp, this%tbsys%H)
      call inverse(this%G_conjg(i), positive = .false.)
    endif
  end do

end subroutine GreensFunctions_calc_Gs

subroutine GreensFunctions_calc_dm_from_Gs(this)
  type(GreensFunctions), intent(inout) :: this

  integer i

  call zero(this%dm)
  do i=1, this%N_G
    if (this%tbsys%complex_matrices) then
      call scaled_accum(this%dm, -this%a(i), this%G(i))
      call scaled_accum(this%dm, -conjg(this%a(i)), this%G_conjg(i))
    else
      call scaled_accum(this%dm, -2.0_dp*this%a(i), this%G(i))
    endif
  end do

  call Gsum_distrib_inplace(this, this%dm)
end subroutine

subroutine GreensFunctions_calc_mod_dm_from_Gs(this, w_e, w_n, do_const_N)
  type(GreensFunctions), intent(inout) :: this
  real(dp), intent(in), pointer :: w_e(:), w_n(:)
  logical, intent(in), optional :: do_const_N

  real(dp), allocatable :: ww_e(:), ww_n(:)
  logical :: const_N = .false.
  real(dp) chempot
  type(TBMatrix),save :: GWE, GWE_H_G, GWN, GWN_S_G

!$OMP threadprivate(GWE,GWN,GWE_H_G,GWN_S_G)

  complex(dp) a, az

  integer i

  if (present(do_const_N)) const_N = do_const_N

  call zero(this%mod_dm_H)
  if (.not. this%tbsys%tbmodel%is_orthogonal) then
    call zero(this%mod_dm_S)
  endif

  allocate(ww_e(this%tbsys%N))
  allocate(ww_n(this%tbsys%N))

  if (associated(w_e)) then
    ww_e = atom_orbital_spread(this%tbsys, w_e)
  else
    ww_e = 1.0_dp
  endif
  if (associated(w_n)) then
    ww_n = atom_orbital_spread(this%tbsys, w_n)
  else
    ww_n = 1.0_dp
  endif

  if (const_N) chempot = GreensFunctions_calc_chempot(this, w_e, w_n)


!$omp parallel private(a,az)
  call Initialise(GWE, this%tbsys%N, this%tbsys%n_matrices, is_complex = .true., scalapack_obj=this%tbsys%scalapack_my_matrices)
  call Initialise(GWE_H_G, this%tbsys%N, this%tbsys%n_matrices, is_complex = .true., scalapack_obj=this%tbsys%scalapack_my_matrices)
  if (const_N) then
    call Initialise(GWN, this%tbsys%N, this%tbsys%n_matrices, is_complex = .true., scalapack_obj=this%tbsys%scalapack_my_matrices)
    call Initialise(GWN_S_G, this%tbsys%N, this%tbsys%n_matrices, is_complex = .true., scalapack_obj=this%tbsys%scalapack_my_matrices)
  endif
!$omp do
  do i=1, this%N_G
    mainlog%mpi_all_inoutput_flag=.true.
    call print("GreensFunctions_calc_mod_dm_from_Gs doing i="//i, PRINT_NERD)
    mainlog%mpi_all_inoutput_flag=.false.
    call multDiag(GWE, this%G(i), ww_e)
    if (const_N) call multDiag(GWN, this%G(i), ww_n)

    call calc_GWAG(GWE, this%G(i), this%tbsys%H, GWE_H_G)
    if (const_N) call calc_GWAG(GWN, this%G(i), this%tbsys%S, GWE_H_G)

    a = this%a(i)
    az = this%a(i)*this%z(i)

    if (.not. this%tbsys%complex_matrices) then ! double instead of doing G(z) and G(conjg(z)) separately
      a = 2.0_dp*a
      az = 2.0_dp*az
    endif

    ! const mu 
    !    F = Re d/dr sum_i Tr [ H a_i (z_i S-H)^-1 W ] 
    !      pull Re and sum in
    !      = d/dr Tr [ H Re sum_i a_i (z_i S - H)^-1 W ] 
    !      chain rule
    !      = Tr [ dH/dr Re sum_i a_i G_i W ] + Tr [ H Re sum_i a_i -G_i (z_i dS/dr - dH/dr) G_i W ]
    !      expand
    !      = Tr [ dH/dr Re sum_i a_i G_i W ] + Tr [ dH/dr Re sum_i a_i (G_i W H G_i ) +
    !                                               dS/dr Re sum_i -a_i z_i G_i W H G_i ];
    !      = Tr [ dH/dr sum_i Re ( a_i G_i W + a_i G_i W H W G_i ) ] + Tr [ dS/dr Re (sum_i -az_i G_i W H G_i) ]
    ! const N
    !   F = Re d/dr Tr [ H sum_i  a_i G_i WE - mu S sum_i a_i G_i WN ]; mu = dE / dN
    !     = f_mu - mu Tr [ dS/dr Re sum_i a_i G_i WN ] - mu Tr [ S Re sum_i a_i -G_i (z_i dS/dr - dH/dr) G_i WN ]
    !     = f_mu - mu Tr [ dS/dr Re sum_i a_i G_i WN ] - mu Tr [ dH/dr Re sum_i a_i G_i WN S G_i -
    !                                                            dS/dr Re sum_i az_i G_i WN S G_i ]
    !     = f_mu + Tr [ dH/dr Re sum_i -mu a_i G_i WN S G_i ] + Tr [ dS/dr Re sum_i -mu a_i G_i WN +
    !                                                              dS/dr Re sum_i mu az_i G_i WN S G_i ]

!$omp critical
    call scaled_accum(this%mod_dm_H, -a, GWE) ! Tr dH/dr  G W contribution to E
    call scaled_accum(this%mod_dm_H, -a, GWE_H_G) ! Tr H dG/dr W = Tr H G dH/dr G W = Tr dH/dr G W H G

    if (const_N) then
      call scaled_accum(this%mod_dm_H, a*chempot, GWN_S_G) ! Tr S dG/dr W = Tr S G dH/dr G W = Tr dH/dr G W S G
    endif

    if (.not. this%tbsys%tbmodel%is_orthogonal) then
      call scaled_accum(this%mod_dm_S, az, GWE_H_G) ! Tr H dG/dr W = Tr H G dS/dr G W = Tr dS/dr G W H G
      if (const_N) then
	call scaled_accum(this%mod_dm_S, a*chempot, GWN) ! Tr dS/dr G W contribution to mu N
	call scaled_accum(this%mod_dm_S, -az*chempot, GWN_S_G) ! Tr S dG/dr W = Tr S G dS/dr G W = Tr dS/dr G W S G
      endif
    endif
!$omp end critical

    if (this%tbsys%complex_matrices) then ! need to do G(conjg(z)) explicitly
      call multDiag(GWE, this%G_conjg(i), ww_e)
      if (const_N) call multDiag(GWN, this%G_conjg(i), ww_n)

      call calc_GWAG(GWE, this%G_conjg(i), this%tbsys%H, GWE_H_G)
      if (const_N) call calc_GWAG(GWN, this%G_conjg(i), this%tbsys%S, GWE_H_G)

      a = conjg(this%a(i))
      az = conjg(this%a(i)*this%z(i))

!$omp critical
      call scaled_accum(this%mod_dm_H, -a, GWE) ! Tr dH/dr  G W contribution to E
      call scaled_accum(this%mod_dm_H, -a, GWE_H_G) ! Tr H dG/dr W = Tr H G dH/dr G W = Tr dH/dr G W H G

      if (const_N) then
	call scaled_accum(this%mod_dm_H, a*chempot, GWN_S_G) ! Tr S dG/dr W = Tr S G dH/dr G W = Tr dH/dr G W S G
      endif

      if (.not. this%tbsys%tbmodel%is_orthogonal) then

	call scaled_accum(this%mod_dm_S, az, GWE_H_G) ! Tr H dG/dr W = Tr H G dS/dr G W = Tr dS/dr G W H G
	if (const_N) then
	  call scaled_accum(this%mod_dm_S, a*chempot, GWN) ! Tr dS/dr G W contribution to mu N
	  call scaled_accum(this%mod_dm_S, -az*chempot, GWN_S_G) ! Tr S dG/dr W = Tr S G dS/dr G W = Tr dS/dr G W S G
	endif
      endif
!$omp end critical

    endif
  end do
  call Finalise(GWE)
  call Finalise(GWE_H_G)
  if (const_N) then
    call Finalise(GWN)
    call Finalise(GWN_S_G)
  endif
!$omp end parallel

  call Gsum_distrib_inplace(this, this%mod_dm_H)
  if (.not. this%tbsys%tbmodel%is_orthogonal) then
    call Gsum_distrib_inplace(this, this%mod_dm_S)
  endif

  deallocate(ww_e, ww_n)


end subroutine

function GreensFunctions_calc_chempot(this, w_e, w_n)
  type(GreensFunctions), intent(in) :: this
  real(dp), pointer :: w_e(:), w_n(:)
  real(dp) :: GreensFunctions_calc_chempot

  integer i
  real(dp) :: dE_dmu = 0.0_dp, dN_dmu = 0.0_dp

  complex(dp), allocatable :: pTr_GSGS(:), pTr_GSGH(:)
  complex(dp) :: Tr_GSGS, Tr_GSGH

  type(TBMatrix) :: GS, GH
  complex a

  call Initialise(GS, this%tbsys%N, this%tbsys%n_matrices, .true., scalapack_obj = this%tbsys%scalapack_my_matrices)
  call Initialise(GH, this%tbsys%N, this%tbsys%n_matrices, .true., scalapack_obj = this%tbsys%scalapack_my_matrices)

  allocate(pTR_GSGS(this%tbsys%N_atoms))
  allocate(pTR_GSGH(this%tbsys%N_atoms))

  do i=1, this%N_G
    a = this%a(i)

    call matrix_product_sub(GS, this%G(i), this%tbsys%S)
    call matrix_product_sub(GH, this%G(i), this%tbsys%H)

    pTr_GSGS = local_ksum(this%tbsys%kpoints, atom_orbital_sum(this%tbsys, partial_TraceMult(GS, GS)))
    call ksum_distrib_inplace(this%tbsys%kpoints, pTr_GSGS)
    if (associated(w_n)) then
      Tr_GSGS = sum(w_n*pTR_GSGS)
    else
      Tr_GSGS = sum(pTR_GSGS)
    endif

    pTr_GSGH = local_ksum(this%tbsys%kpoints, atom_orbital_sum(this%tbsys, partial_TraceMult(GS, GH)))
    call ksum_distrib_inplace(this%tbsys%kpoints, pTr_GSGH)
    if (associated(w_e)) then
      Tr_GSGH = sum(w_e*pTR_GSGH)
    else
      Tr_GSGH = sum(pTR_GSGH)
    endif

    dE_dmu = dE_dmu - 2.0_dp*a*Tr_GSGH
    dN_dmu = dN_dmu - 2.0_dp*a*Tr_GSGS
  end do

  dE_dmu = Gsum_distrib(this, dE_dmu)
  dN_dmu = Gsum_distrib(this, dN_dmu)

  GreensFunctions_calc_chempot = dE_dmu / dN_dmu

end function GreensFunctions_calc_chempot

subroutine calc_GWAG(GW, G, A, GWAG)
  type(TBMatrix), intent(in) :: GW, A, G
  type(TBMatrix), intent(inout) :: GWAG

  type(TBMatrix) :: t

  call Initialise(t, GW)

  call matrix_product_sub(t, A, G)
  call matrix_product_sub(GWAG, GW, t)

  call Finalise(t)
end subroutine calc_GWAG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GreensFunctions_Gsum_distrib_inplace_tbm(this, M)
  type(GreensFunctions), intent(in) :: this
  type(TBMatrix), intent(inout) :: M

  call sum_in_place(M, this%mpi_across_poles)
end subroutine

subroutine GreensFunctions_Gsum_distrib_inplace_r2(this, M)
  type(GreensFunctions), intent(in) :: this
  real(dp), intent(inout) :: M(:,:)

  call sum_in_place(this%mpi_across_poles, M)
end subroutine

subroutine GreensFunctions_Gsum_distrib_inplace_r3(this, M)
  type(GreensFunctions), intent(in) :: this
  real(dp), intent(inout) :: M(:,:,:)

  call sum_in_place(this%mpi_across_poles, M)
end subroutine

subroutine GreensFunctions_Gsum_distrib_inplace_matd(this, M)
  type(GreensFunctions), intent(in) :: this
  type(MatrixD), intent(inout) :: M

  call sum_in_place(this%mpi_across_poles, M%data)
end subroutine

function GreensFunctions_Gsum_distrib_d(this, d)
  type(GreensFunctions), intent(in) :: this
  real(dp), intent(in) :: d
  real(dp) :: GreensFunctions_Gsum_distrib_d

  GreensFunctions_Gsum_distrib_d = sum(this%mpi_across_poles, d)
end function GreensFunctions_Gsum_distrib_d

end module TB_GreensFunctions_module
