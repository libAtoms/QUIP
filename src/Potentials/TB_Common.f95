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
!X TB_Common_module  
!X
!% TB module computes the tight-binding Hamiltonian matrix elements 
!% in the two-center approximation, applying the general expressions 
!% derived by A.V. Podolskiy and P. Vogl Phys. Rev. B {\bf 69}, 233101 (2004).
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


module TB_Common_module

use System_module
implicit none
private

integer, parameter, public :: ORB_S = 1, ORB_P = 2, ORB_D = 3, ORB_F = 4

integer, parameter, public :: N_ORBS_OF_SET(4) = (/1, 3, 5, 7/)



integer, parameter, public :: harrison_sign_sss = -1, harrison_sign_sps = 1, &
		      harrison_sign_pps = 1, harrison_sign_ppp = -1, &
	              harrison_sign_sds = -1, harrison_sign_pds = -1, &
		      harrison_sign_pdp = 1, harrison_sign_dds = -1, &
		      harrison_sign_ddp = 1

integer, parameter, public :: N_SK = 20
integer, parameter, public :: SK_SSS = 1, SK_SPS = 2, SK_PPS = 3, SK_PPP = 4, &
		      SK_SDS = 5, SK_PDS = 6, SK_PDP = 7, SK_DDS = 8, SK_DDP = 9, SK_DDD = 10, \
		      SK_SFS = 11, SK_PFS = 12, SK_PFP = 13, SK_DFS = 14, SK_DFP = 15, SK_DFD = 16, \
		      SK_FFS = 17, SK_FFP = 18, SK_FFD = 19, SK_FFF = 20

public :: angular_function, dangular_function, spin_orbit_function

contains

function angular_function(dcos, dcos_sq, orb_type_i, orb_type_j, orb_dir_i, orb_dir_j, sk) result(V)
  real(dp), intent(in) :: dcos(3), dcos_sq(3)
  integer, intent(in) :: orb_type_i, orb_type_j
  integer, intent(in) :: orb_dir_i, orb_dir_j
  real(dp), intent(in) :: sk(N_SK)
  real(dp) :: V

  real(dp) :: u_V_sss, u_V_sps, u_V_pps, u_V_ppp
  real(dp) :: u_V_sds, u_V_pds, u_V_pdp
  real(dp) :: u_V_dds, u_V_ddp, u_V_ddd
  real(dp) :: u_V_sfs, u_V_pfs, u_V_pfp
  real(dp) :: u_V_dfs, u_V_dfp, u_V_dfd
  real(dp) :: u_V_ffs, u_V_ffp, u_V_ffd, u_V_fff

  real(dp) :: L, M, N, Lsq, Msq, Nsq

  real(dp) :: root_3
  
  root_3 = sqrt(3.0_dp)

  ! L on L
  u_V_sss = sk(SK_SSS)

  u_V_pps = sk(SK_PPS)
  u_V_ppp = sk(SK_PPP)

  u_V_dds = sk(SK_DDS)
  u_V_ddp = sk(SK_DDP)
  u_V_ddd = sk(SK_DDD)

  u_V_ffs = sk(SK_FFS)
  u_V_ffp = sk(SK_FFP)
  u_V_ffd = sk(SK_FFD)
  u_V_fff = sk(SK_FFF)

  ! L on L +/- 2
  u_V_sds = sk(SK_SDS)

  u_V_pfs = sk(SK_PFS)
  u_V_pfp = sk(SK_PFP)

  ! SK_vogl.h does also its own sign flip
  if (orb_type_i < orb_type_j) then
     ! L on L +/- 1
    u_V_sps = sk(SK_SPS)

    u_V_pds = sk(SK_PDS)
    u_V_pdp = sk(SK_PDP)

    u_V_dfs = sk(SK_DFS)
    u_V_dfp = sk(SK_DFP)
    u_V_dfd = sk(SK_DFD)

    ! L on L +/- 3
    u_V_sfs = sk(SK_SFS)
  else
     ! L on L +/- 1
    u_V_sps = -sk(SK_SPS)

    u_V_pds = -sk(SK_PDS)
    u_V_pdp = -sk(SK_PDP)

    u_V_dfs = -sk(SK_DFS)
    u_V_dfp = -sk(SK_DFP)
    u_V_dfd = -sk(SK_DFD)

    ! L on L +/- 3
    u_V_sfs = -sk(SK_SFS)
  endif

  L = dcos(1)
  M = dcos(2)
  N = dcos(3)
  Lsq = dcos_sq(1)
  Msq = dcos_sq(2)
  Nsq = dcos_sq(3)

  ! seems to prevent intel 9.1.051 -O3 optimizer from breaking S1-D5 matrix element (and maybe others?)
  if (orb_type_i < 0) call print("")

  ! workaround for (-1+Nsq) division by zero
  call regularize(L, M, N, Lsq, Msq, Nsq)

include 'SK_vogl.h'

end function angular_function

subroutine regularize(L, M, N, Lsq, Msq, Nsq)
    real(dp), intent(inout) :: L, M, N, Lsq, Msq, Nsq

    real(dp), parameter :: regularization_eps = 1.0e-5

    if (Nsq == 1.0_dp) then
        L = sign(regularization_eps, N)
        Lsq = L**2
        M = sign(regularization_eps, N)
        Msq = M**2
        Nsq = (1.0_dp - Lsq - Msq)
        N = sign(sqrt(Nsq), N)
    else if (Msq == 1.0_dp) then
        L = sign(regularization_eps, M)
        Lsq = L**2
        N = sign(regularization_eps, M)
        Nsq = N**2
        Msq = (1.0_dp - Lsq - Nsq)
        M = sign(sqrt(Msq), M)
    else if (Lsq == 1.0_dp) then
        M = sign(regularization_eps, L)
        Msq = M**2
        N = sign(regularization_eps, L)
        Nsq = N**2
        Lsq = (1.0_dp - Nsq - Msq)
        L = sign(sqrt(Lsq),L)
    endif
end subroutine regularize


function spin_orbit_function(orb_type_i, orb_dir_i, orb_dir_j) result(V)
  integer, intent(in) :: orb_type_i, orb_dir_i, orb_dir_j
  complex(dp) :: V(2,2)

  real(dp) :: root_3
  
  root_3 = sqrt(3.0_dp)

  ! spin-orbit term 
  ! from Podolskiy and Vogl, Phys. Rev. B v. 69, p 233101 (2004)

include 'SK_SO_vogl.h'

end function spin_orbit_function

function dangular_function(dist, dcos, dcos_sq, orb_type_i, orb_type_j, &
			   orb_dir_i, orb_dir_j, sk, dsk)
  real(dp), intent(in) :: dist, dcos(3), dcos_sq(3)
  integer, intent(in) :: orb_type_i, orb_type_j
  integer, intent(in) :: orb_dir_i, orb_dir_j
  real(dp), intent(in) :: sk(N_SK), dsk(N_SK)
  real(dp) :: dangular_function(3)

  real(dp) ::  rVd(3), rVd_t

  real(dp) :: u_V_sss, u_V_sps, u_V_pps, u_V_ppp
  real(dp) :: u_V_sds, u_V_pds, u_V_pdp
  real(dp) :: u_V_dds, u_V_ddp, u_V_ddd
  real(dp) :: u_V_sfs, u_V_pfs, u_V_pfp
  real(dp) :: u_V_dfs, u_V_dfp, u_V_dfd
  real(dp) :: u_V_ffs, u_V_ffp, u_V_ffd, u_V_fff

  real(dp) :: u_Vd_sss, u_Vd_sps, u_Vd_pps, u_Vd_ppp
  real(dp) :: u_Vd_sds, u_Vd_pds, u_Vd_pdp
  real(dp) :: u_Vd_dds, u_Vd_ddp, u_Vd_ddd
  real(dp) :: u_Vd_sfs, u_Vd_pfs, u_Vd_pfp
  real(dp) :: u_Vd_dfs, u_Vd_dfp, u_Vd_dfd
  real(dp) :: u_Vd_ffs, u_Vd_ffp, u_Vd_ffd, u_Vd_fff

  real(dp) :: L, M, N, Lsq, Msq, Nsq
  real(dp) :: dL_dr(3), dM_dr(3), dN_dr(3)

  real(dp) :: root_3
  
  root_3 = sqrt(3.0_dp)

  ! L on L
  u_V_sss = sk(SK_SSS)

  u_V_pps = sk(SK_PPS)
  u_V_ppp = sk(SK_PPP)

  u_V_dds = sk(SK_DDS)
  u_V_ddp = sk(SK_DDP)
  u_V_ddd = sk(SK_DDD)

  u_V_ffs = sk(SK_FFS)
  u_V_ffp = sk(SK_FFP)
  u_V_ffd = sk(SK_FFD)
  u_V_fff = sk(SK_FFF)

  u_Vd_sss = dsk(SK_SSS)

  u_Vd_pps = dsk(SK_PPS)
  u_Vd_ppp = dsk(SK_PPP)

  u_Vd_dds = dsk(SK_DDS)
  u_Vd_ddp = dsk(SK_DDP)
  u_Vd_ddd = dsk(SK_DDD)

  u_Vd_ffs = dsk(SK_FFS)
  u_Vd_ffp = dsk(SK_FFP)
  u_Vd_ffd = dsk(SK_FFD)
  u_Vd_fff = dsk(SK_FFF)

  ! L on L +/- 2
  u_V_sds = sk(SK_SDS)

  u_V_pfs = sk(SK_PFS)
  u_V_pfp = sk(SK_PFP)

  u_Vd_sds = dsk(SK_SDS)

  u_Vd_pfs = dsk(SK_PFS)
  u_Vd_pfp = dsk(SK_PFP)

  ! SK_vogl.h does also its own sign flip
  if (orb_type_i < orb_type_j) then
     ! L on L +/- 1
    u_V_sps = sk(SK_SPS)

    u_V_pds = sk(SK_PDS)
    u_V_pdp = sk(SK_PDP)

    u_V_dfs = sk(SK_DFS)
    u_V_dfp = sk(SK_DFP)
    u_V_dfd = sk(SK_DFD)

    u_Vd_sps = dsk(SK_SPS)

    u_Vd_pds = dsk(SK_PDS)
    u_Vd_pdp = dsk(SK_PDP)

    u_Vd_dfs = dsk(SK_DFS)
    u_Vd_dfp = dsk(SK_DFP)
    u_Vd_dfd = dsk(SK_DFD)

    ! L on L +/- 3
    u_V_sfs = sk(SK_SFS)

    u_Vd_sfs = dsk(SK_SFS)
  else
     ! L on L +/- 1
    u_V_sps = -sk(SK_SPS)

    u_V_pds = -sk(SK_PDS)
    u_V_pdp = -sk(SK_PDP)

    u_V_dfs = -sk(SK_DFS)
    u_V_dfp = -sk(SK_DFP)
    u_V_dfd = -sk(SK_DFD)

    u_Vd_sps = -dsk(SK_SPS)

    u_Vd_pds = -dsk(SK_PDS)
    u_Vd_pdp = -dsk(SK_PDP)

    u_Vd_dfs = -dsk(SK_DFS)
    u_Vd_dfp = -dsk(SK_DFP)
    u_Vd_dfd = -dsk(SK_DFD)

    ! L on L +/- 3
    u_V_sfs = -sk(SK_SFS)

    u_Vd_sfs = -dsk(SK_SFS)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  L = dcos(1)
  M = dcos(2)
  N = dcos(3)
  Lsq = dcos_sq(1)
  Msq = dcos_sq(2)
  Nsq = dcos_sq(3)

  ! seems to prevent intel 9.1.051 -O3 optimizer from breaking S1-D5 matrix element (and maybe others?)
  if (orb_type_i < 0) call print("")

  ! workaround for (-1+Nsq) division by zero
  call regularize(L, M, N, Lsq, Msq, Nsq)

include 'SKd_vogl.h'

  dL_dr(1) = (-1.0_dp+L**2)/dist
  dL_dr(2) = L*M/dist
  dL_dr(3) = L*N/dist

  dM_dr(1) = M*L/dist
  dM_dr(2) = (-1.0_dp+M**2)/dist
  dM_dr(3) = M*N/dist

  dN_dr(1) = N*L/dist
  dN_dr(2) = N*M/dist
  dN_dr(3) = (-1.0_dp+N**2)/dist

  dangular_function(:) = (rVd(1)*dL_dr(:) + rVd(2)*dM_dr(:) + rVd(3)*dN_dr(:) - rVd_t*dcos(:))

end function dangular_function

end module TB_Common_module
