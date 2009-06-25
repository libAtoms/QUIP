!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     QUIP: quantum mechanical and interatomic potential simulation package
!X     
!X     Portions written by Noam Bernstein, while working at the
!X     Naval Research Laboratory, Washington DC. 
!X
!X     Portions written by Gabor Csanyi, Copyright 2006-2007.   
!X
!X     When using this software,  please cite the following reference:
!X
!X     reference
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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

integer, parameter, public :: ORB_S = 1, ORB_P = 2, ORB_D = 3

integer, parameter, public :: N_ORBS_OF_SET(3) = (/1, 3, 5/)



integer, parameter, public :: harrison_sign_sss = -1, harrison_sign_sps = 1, &
		      harrison_sign_pps = 1, harrison_sign_ppp = -1, &
	              harrison_sign_sds = -1, harrison_sign_pds = -1, &
		      harrison_sign_pdp = 1, harrison_sign_dds = -1, &
		      harrison_sign_ddp = 1

integer, parameter, public :: SK_SSS = 1, SK_SPS = 2, SK_PPS = 3, SK_PPP = 4, &
		      SK_SDS = 5, SK_PDS = 6, SK_PDP = 7, SK_DDS = 8, SK_DDP = 9, SK_DDD = 10

public :: angular_function, dangular_function, spin_orbit_function

contains

function angular_function(dcos, dcos_sq, orb_type_i, orb_type_j, orb_dir_i, orb_dir_j, sk) result(V)
  real(dp), intent(in) :: dcos(3), dcos_sq(3)
  integer, intent(in) :: orb_type_i, orb_type_j
  integer, intent(in) :: orb_dir_i, orb_dir_j
  real(dp), intent(in) :: sk(10)
  real(dp) :: V

  real(dp) :: u_V_sss, u_V_sps, u_V_pps, u_V_ppp
  real(dp) :: u_V_sds, u_V_pds, u_V_pdp
  real(dp) :: u_V_dds, u_V_ddp, u_V_ddd
  real(dp) :: L, M, N, Lsq, Msq, Nsq

  real(dp) :: root_3
  
  root_3 = sqrt(3.0_dp)

  u_V_sss = sk(SK_SSS)

  u_V_pps = sk(SK_PPS)
  u_V_ppp = sk(SK_PPP)

  u_V_dds = sk(SK_DDS)
  u_V_ddp = sk(SK_DDP)
  u_V_ddd = sk(SK_DDD)

  u_V_sds = sk(SK_SDS)
  u_V_pdp = sk(SK_PDP)

  ! SK_vogl.h does also its own sign flip
  if (orb_type_i < orb_type_j) then
    u_V_sps = sk(SK_SPS)
    u_V_pds = sk(SK_PDS)
    u_V_pdp = sk(SK_PDP)
  else
    u_V_sps = -sk(SK_SPS)
    u_V_pds = -sk(SK_PDS)
    u_V_pdp = -sk(SK_PDP)
  endif

  L = dcos(1)
  M = dcos(2)
  N = dcos(3)
  Lsq = dcos_sq(1)
  Msq = dcos_sq(2)
  Nsq = dcos_sq(3)

  ! seems to prevent intel 9.1.051 -O3 optimizer from breaking S1-D5 matrix element (and maybe others?)
  if (orb_type_i < 0) call print("")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! From SK_vogl.h, d orbitals only right now
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! workaround for (-1+Nsq) division by zero
  if (Nsq == 1.0_dp) then
    Nsq = Nsq - 1.0e-16_dp
  endif

! angular functions 
! from Podolskiy and Vogl, Phys. Rev. B v. 69, p 233101 (2004)
select case (orb_type_i)
  case (ORB_S)
    select case (orb_type_j)
      case (ORB_S)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = u_V_SSS
            end select
        end select
      case (ORB_P)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = M*u_V_SPS
              case (2)
              V = N*u_V_SPS
              case (3)
              V = L*u_V_SPS
            end select
        end select
      case (ORB_D)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = root_3*L*M*u_V_SDS
              case (2)
              V = root_3*M*N*u_V_SDS
              case (3)
              V = ((-1 + 3*Nsq)*u_V_SDS)/2.
              case (4)
              V = root_3*L*N*u_V_SDS
              case (5)
              V = -(root_3*(-1 + 2*Msq + Nsq)*u_V_SDS)/2.
            end select
        end select
    end select
  case (ORB_P)
    select case (orb_type_j)
      case (ORB_S)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = -(M*u_V_SPS)
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              V = -(N*u_V_SPS)
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              V = -(L*u_V_SPS)
            end select
        end select
      case (ORB_P)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = Msq*(u_V_PPS - u_V_PPP) + u_V_PPP
              case (2)
              V = M*N*(u_V_PPS - u_V_PPP)
              case (3)
              V = L*M*(u_V_PPS - u_V_PPP)
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              V = M*N*(u_V_PPS - u_V_PPP)
              case (2)
              V = Nsq*(u_V_PPS - u_V_PPP) + u_V_PPP
              case (3)
              V = L*N*(u_V_PPS - u_V_PPP)
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              V = L*M*(u_V_PPS - u_V_PPP)
              case (2)
              V = L*N*(u_V_PPS - u_V_PPP)
              case (3)
              V = -((-1 + Msq + Nsq)*u_V_PPS) + (Msq + Nsq)*u_V_PPP
            end select
        end select
      case (ORB_D)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = L*(Msq*(root_3*u_V_PDS - 2*u_V_PDP) + u_V_PDP)
              case (2)
              V = N*(Msq*(root_3*u_V_PDS - 2*u_V_PDP) + u_V_PDP)
              case (3)
              V = (-(M*u_V_PDS) + M*Nsq*(3*u_V_PDS - 2*root_3*u_V_PDP))/2.
              case (4)
              V = L*M*N*(root_3*u_V_PDS - 2*u_V_PDP)
              case (5)
              V = -(M*(root_3*(-1 + 2*Msq + Nsq)*u_V_PDS -  &
         2*(-2 + 2*Msq + Nsq)*u_V_PDP))/2.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              V = L*M*N*(root_3*u_V_PDS - 2*u_V_PDP)
              case (2)
              V = M*(Nsq*(root_3*u_V_PDS - 2*u_V_PDP) + u_V_PDP)
              case (3)
              V = (N*((-1 + 3*Nsq)*u_V_PDS - 2*root_3*(-1 + Nsq)*u_V_PDP))/2.
              case (4)
              V = L*(Nsq*(root_3*u_V_PDS - 2*u_V_PDP) + u_V_PDP)
              case (5)
              V = -(N*(-1 + 2*Msq + Nsq)*(root_3*u_V_PDS - 2*u_V_PDP))/2.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              V = M*(-(root_3*(-1 + Msq + Nsq)*u_V_PDS) +  &
      (-1 + 2*Msq + 2*Nsq)*u_V_PDP)
              case (2)
              V = L*M*N*(root_3*u_V_PDS - 2*u_V_PDP)
              case (3)
              V = (-(L*u_V_PDS) + L*Nsq*(3*u_V_PDS - 2*root_3*u_V_PDP))/2.
              case (4)
              V = N*(-(root_3*(-1 + Msq + Nsq)*u_V_PDS) +  &
      (-1 + 2*Msq + 2*Nsq)*u_V_PDP)
              case (5)
              V = -(L*(root_3*(-1 + 2*Msq + Nsq)*u_V_PDS -  &
         2*(2*Msq + Nsq)*u_V_PDP))/2.
            end select
        end select
    end select
  case (ORB_D)
    select case (orb_type_j)
      case (ORB_S)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = root_3*L*M*u_V_SDS
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              V = root_3*M*N*u_V_SDS
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              V = ((-1 + 3*Nsq)*u_V_SDS)/2.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              V = root_3*L*N*u_V_SDS
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              V = -(root_3*(-1 + 2*Msq + Nsq)*u_V_SDS)/2.
            end select
        end select
      case (ORB_P)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = -(L*(Msq*(root_3*u_V_PDS - 2*u_V_PDP) + u_V_PDP))
              case (2)
              V = L*M*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              case (3)
              V = M*(root_3*(-1 + Msq + Nsq)*u_V_PDS +  &
      (1 - 2*Msq - 2*Nsq)*u_V_PDP)
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              V = -(N*(Msq*(root_3*u_V_PDS - 2*u_V_PDP) + u_V_PDP))
              case (2)
              V = -(M*(Nsq*(root_3*u_V_PDS - 2*u_V_PDP) + u_V_PDP))
              case (3)
              V = L*M*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              V = (M*((1 - 3*Nsq)*u_V_PDS + 2*root_3*Nsq*u_V_PDP))/2.
              case (2)
              V = (N*((1 - 3*Nsq)*u_V_PDS + 2*root_3*(-1 + Nsq)*u_V_PDP))/2.
              case (3)
              V = (L*((1 - 3*Nsq)*u_V_PDS + 2*root_3*Nsq*u_V_PDP))/2.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              V = L*M*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              case (2)
              V = -(L*(Nsq*(root_3*u_V_PDS - 2*u_V_PDP) + u_V_PDP))
              case (3)
              V = N*(root_3*(-1 + Msq + Nsq)*u_V_PDS +  &
      (1 - 2*Msq - 2*Nsq)*u_V_PDP)
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              V = (M*(root_3*(-1 + 2*Msq + Nsq)*u_V_PDS -  &
        2*(-2 + 2*Msq + Nsq)*u_V_PDP))/2.
              case (2)
              V = (N*(-1 + 2*Msq + Nsq)*(root_3*u_V_PDS - 2*u_V_PDP))/2.
              case (3)
              V = (L*(root_3*(-1 + 2*Msq + Nsq)*u_V_PDS -  &
        2*(2*Msq + Nsq)*u_V_PDP))/2.
            end select
        end select
      case (ORB_D)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = u_V_DDP - M**4*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD) -  &
    Msq*(-1 + Nsq)*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD) +  &
    Nsq*(-u_V_DDP + u_V_DDD)
              case (2)
              V = L*N*(u_V_DDP - u_V_DDD +  &
      Msq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              case (3)
              V = (root_3*L*M*((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
        Nsq*(-4*u_V_DDP + u_V_DDD)))/2.
              case (4)
              V = M*N*(-3*(-1 + Msq + Nsq)*u_V_DDS +  &
      (-3 + 4*Msq + 4*Nsq)*u_V_DDP - (Msq + Nsq)*u_V_DDD)
              case (5)
              V = -(L*M*(-1 + 2*Msq + Nsq)* &
       (3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              V = L*N*(u_V_DDP - u_V_DDD +  &
      Msq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              case (2)
              V = Nsq*(u_V_DDP - u_V_DDD) + u_V_DDD +  &
    Msq*(u_V_DDP - u_V_DDD +  &
       Nsq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              case (3)
              V = (root_3*M*N*((-1 + 3*Nsq)*u_V_DDS + (2 - 4*Nsq)*u_V_DDP +  &
        (-1 + Nsq)*u_V_DDD))/2.
              case (4)
              V = L*M*(u_V_DDP - u_V_DDD +  &
      Nsq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              case (5)
              V = -(M*N*(3*(-1 + 2*Msq + Nsq)*u_V_DDS +  &
         (6 - 8*Msq - 4*Nsq)*u_V_DDP + (-3 + 2*Msq + Nsq)*u_V_DDD))/ &
    2.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              V = (root_3*L*M*((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
        Nsq*(-4*u_V_DDP + u_V_DDD)))/2.
              case (2)
              V = (root_3*M*N*((-1 + 3*Nsq)*u_V_DDS + (2 - 4*Nsq)*u_V_DDP +  &
        (-1 + Nsq)*u_V_DDD))/2.
              case (3)
              V = ((1 - 3*Nsq)**2*u_V_DDS -  &
      3*(-1 + Nsq)*(Nsq*(4*u_V_DDP - u_V_DDD) + u_V_DDD))/4.
              case (4)
              V = (root_3*L*N*((-1 + 3*Nsq)*u_V_DDS + (2 - 4*Nsq)*u_V_DDP +  &
        (-1 + Nsq)*u_V_DDD))/2.
              case (5)
              V = -(root_3*(-1 + 2*Msq + Nsq)* &
       ((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
         Nsq*(-4*u_V_DDP + u_V_DDD)))/4.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              V = M*N*(-3*(-1 + Msq + Nsq)*u_V_DDS +  &
      (-3 + 4*Msq + 4*Nsq)*u_V_DDP - (Msq + Nsq)*u_V_DDD)
              case (2)
              V = L*M*(u_V_DDP - u_V_DDD +  &
      Nsq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              case (3)
              V = (root_3*L*N*((-1 + 3*Nsq)*u_V_DDS + (2 - 4*Nsq)*u_V_DDP +  &
        (-1 + Nsq)*u_V_DDD))/2.
              case (4)
              V = u_V_DDP - (-1 + Msq)*Nsq* &
     (3*u_V_DDS - 4*u_V_DDP + u_V_DDD) -  &
    N**4*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD) +  &
    Msq*(-u_V_DDP + u_V_DDD)
              case (5)
              V = -(L*N*(3*(-1 + 2*Msq + Nsq)*u_V_DDS +  &
         (2 - 8*Msq - 4*Nsq)*u_V_DDP + (1 + 2*Msq + Nsq)*u_V_DDD))/2.
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              V = -(L*M*(-1 + 2*Msq + Nsq)* &
       (3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
              case (2)
              V = -(M*N*(3*(-1 + 2*Msq + Nsq)*u_V_DDS +  &
         (6 - 8*Msq - 4*Nsq)*u_V_DDP + (-3 + 2*Msq + Nsq)*u_V_DDD))/ &
    2.
              case (3)
              V = -(root_3*(-1 + 2*Msq + Nsq)* &
       ((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
         Nsq*(-4*u_V_DDP + u_V_DDD)))/4.
              case (4)
              V = -(L*N*(3*(-1 + 2*Msq + Nsq)*u_V_DDS +  &
         (2 - 8*Msq - 4*Nsq)*u_V_DDP + (1 + 2*Msq + Nsq)*u_V_DDD))/2.
              case (5)
              V = (3*(-1 + 2*Msq + Nsq)**2*u_V_DDS -  &
      4*Msq*(-1 + Nsq)*(4*u_V_DDP - u_V_DDD) + u_V_DDD +  &
      4*M**4*(-4*u_V_DDP + u_V_DDD) +  &
      N**4*(-4*u_V_DDP + u_V_DDD) + 2*Nsq*(2*u_V_DDP + u_V_DDD))/ &
    4.
            end select
        end select
    end select
end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end of SK_vogl.h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end function angular_function

function spin_orbit_function(orb_type_i, orb_dir_i, orb_dir_j) result(V)
  integer, intent(in) :: orb_type_i, orb_dir_i, orb_dir_j
  complex(dp) :: V(2,2)

  real(dp) :: root_3
  
  root_3 = sqrt(3.0_dp)

  ! spin-orbit term 
  ! from Podolskiy and Vogl, Phys. Rev. B v. 69, p 233101 (2004)
  select case (orb_type_i)
    case (ORB_S)
      select case (orb_dir_i)
	case (1)
	  select case (orb_dir_j)
	    case (1)
	    V(1,1) = 0
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = 0
	  end select
      end select
    case (ORB_P)
      select case (orb_dir_i)
	case (1)
	  select case (orb_dir_j)
	    case (1)
	    V(1,1) = 0
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = 0
	    case (2)
	    V(1,1) = 0
	    V(1,2) = cmplx(0,-0.5,dp)
	    V(2,1) = cmplx(0,-0.5,dp)
	    V(2,2) = 0
	    case (3)
	    V(1,1) = cmplx(0,0.5,dp)
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = cmplx(0,-0.5,dp)
	  end select
	case (2)
	  select case (orb_dir_j)
	    case (1)
	    V(1,1) = 0
	    V(1,2) = cmplx(0,0.5,dp)
	    V(2,1) = cmplx(0,0.5,dp)
	    V(2,2) = 0
	    case (2)
	    V(1,1) = 0
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = 0
	    case (3)
	    V(1,1) = 0
	    V(1,2) = -0.5
	    V(2,1) = 0.5
	    V(2,2) = 0
	  end select
	case (3)
	  select case (orb_dir_j)
	    case (1)
	    V(1,1) = cmplx(0,-0.5,dp)
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = cmplx(0,0.5,dp)
	    case (2)
	    V(1,1) = 0
	    V(1,2) = 0.5
	    V(2,1) = -0.5
	    V(2,2) = 0
	    case (3)
	    V(1,1) = 0
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = 0
	  end select
      end select
    case (ORB_D)
      select case (orb_dir_i)
	case (1)
	  select case (orb_dir_j)
	    case (1)
	    V(1,1) = 0
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = 0
	    case (2)
	    V(1,1) = 0
	    V(1,2) = 0.5
	    V(2,1) = -0.5
	    V(2,2) = 0
	    case (3)
	    V(1,1) = 0
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = 0
	    case (4)
	    V(1,1) = 0
	    V(1,2) = cmplx(0,-0.5,dp)
	    V(2,1) = cmplx(0,-0.5,dp)
	    V(2,2) = 0
	    case (5)
	    V(1,1) = cmplx(0,1,dp)
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = cmplx(0,-1,dp)
	  end select
	case (2)
	  select case (orb_dir_j)
	    case (1)
	    V(1,1) = 0
	    V(1,2) = -0.5
	    V(2,1) = 0.5
	    V(2,2) = 0
	    case (2)
	    V(1,1) = 0
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = 0
	    case (3)
	    V(1,1) = 0
	    V(1,2) = cmplx(0,-0.5,dp)*root_3
	    V(2,1) = cmplx(0,-0.5,dp)*root_3
	    V(2,2) = 0
	    case (4)
	    V(1,1) = cmplx(0,0.5,dp)
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = cmplx(0,-0.5,dp)
	    case (5)
	    V(1,1) = 0
	    V(1,2) = cmplx(0,-0.5,dp)
	    V(2,1) = cmplx(0,-0.5,dp)
	    V(2,2) = 0
	  end select
	case (3)
	  select case (orb_dir_j)
	    case (1)
	    V(1,1) = 0
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = 0
	    case (2)
	    V(1,1) = 0
	    V(1,2) = cmplx(0,0.5,dp)*root_3
	    V(2,1) = cmplx(0,0.5,dp)*root_3
	    V(2,2) = 0
	    case (3)
	    V(1,1) = 0
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = 0
	    case (4)
	    V(1,1) = 0
	    V(1,2) = -root_3/2.
	    V(2,1) = root_3/2.
	    V(2,2) = 0
	    case (5)
	    V(1,1) = 0
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = 0
	  end select
	case (4)
	  select case (orb_dir_j)
	    case (1)
	    V(1,1) = 0
	    V(1,2) = cmplx(0,0.5,dp)
	    V(2,1) = cmplx(0,0.5,dp)
	    V(2,2) = 0
	    case (2)
	    V(1,1) = cmplx(0,-0.5,dp)
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = cmplx(0,0.5,dp)
	    case (3)
	    V(1,1) = 0
	    V(1,2) = root_3/2.
	    V(2,1) = -root_3/2.
	    V(2,2) = 0
	    case (4)
	    V(1,1) = 0
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = 0
	    case (5)
	    V(1,1) = 0
	    V(1,2) = -0.5
	    V(2,1) = 0.5
	    V(2,2) = 0
	  end select
	case (5)
	  select case (orb_dir_j)
	    case (1)
	    V(1,1) = cmplx(0,-1,dp)
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = cmplx(0,1,dp)
	    case (2)
	    V(1,1) = 0
	    V(1,2) = cmplx(0,0.5,dp)
	    V(2,1) = cmplx(0,0.5,dp)
	    V(2,2) = 0
	    case (3)
	    V(1,1) = 0
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = 0
	    case (4)
	    V(1,1) = 0
	    V(1,2) = 0.5
	    V(2,1) = -0.5
	    V(2,2) = 0
	    case (5)
	    V(1,1) = 0
	    V(1,2) = 0
	    V(2,1) = 0
	    V(2,2) = 0
	  end select
      end select
  end select

end function spin_orbit_function

function dangular_function(dist, dcos, dcos_sq, orb_type_i, orb_type_j, &
			   orb_dir_i, orb_dir_j, sk, dsk)
  real(dp), intent(in) :: dist, dcos(3), dcos_sq(3)
  integer, intent(in) :: orb_type_i, orb_type_j
  integer, intent(in) :: orb_dir_i, orb_dir_j
  real(dp), intent(in) :: sk(10), dsk(10)
  real(dp) :: dangular_function(3)

  real(dp) ::  rVd(3), rVd_t

  real(dp) :: u_V_sss, u_V_sps, u_V_pps, u_V_ppp
  real(dp) :: u_V_sds, u_V_pds, u_V_pdp
  real(dp) :: u_V_dds, u_V_ddp, u_V_ddd
  real(dp) :: u_Vd_sss, u_Vd_sps, u_Vd_pps, u_Vd_ppp
  real(dp) :: u_Vd_sds, u_Vd_pds, u_Vd_pdp
  real(dp) :: u_Vd_dds, u_Vd_ddp, u_Vd_ddd
  real(dp) :: L, M, N, Lsq, Msq, Nsq
  real(dp) :: dL_dr(3), dM_dr(3), dN_dr(3)

  real(dp) :: root_3
  
  root_3 = sqrt(3.0_dp)

  u_V_sss = sk(SK_SSS)


  u_V_pps = sk(SK_PPS)
  u_V_ppp = sk(SK_PPP)

  u_V_sds = sk(SK_SDS)

  u_V_pdp = sk(SK_PDP)

  u_V_dds = sk(SK_DDS)
  u_V_ddp = sk(SK_DDP)
  u_V_ddd = sk(SK_DDD)

  u_Vd_sss = dsk(SK_SSS)


  u_Vd_pps = dsk(SK_PPS)
  u_Vd_ppp = dsk(SK_PPP)

  u_Vd_sds = dsk(SK_SDS)

  u_Vd_pdp = dsk(SK_PDP)

  u_Vd_dds = dsk(SK_DDS)
  u_Vd_ddp = dsk(SK_DDP)
  u_Vd_ddd = dsk(SK_DDD)

  ! SK_vogl.h does also its own sign flip
  if (orb_type_i < orb_type_j) then
    u_V_sps = sk(SK_SPS)
    u_Vd_sps = dsk(SK_SPS)
    u_V_pds = sk(SK_PDS)
    u_Vd_pds = dsk(SK_PDS)
    u_V_pdp = sk(SK_PDP)
    u_Vd_pdp = dsk(SK_PDP)
  else
    u_V_sps = -sk(SK_SPS)
    u_Vd_sps = -dsk(SK_SPS)
    u_V_pds = -sk(SK_PDS)
    u_Vd_pds = -dsk(SK_PDS)
    u_V_pdp = -sk(SK_PDP)
    u_Vd_pdp = -dsk(SK_PDP)
  endif

  L = dcos(1)
  M = dcos(2)
  N = dcos(3)
  Lsq = dcos_sq(1)
  Msq = dcos_sq(2)
  Nsq = dcos_sq(3)

  ! seems to prevent intel 9.1.051 -O3 optimizer from breaking S1-D5 matrix element (and maybe others?)
  if (orb_type_i < 0) call print("")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! From SKd_vogl.h, d orbitals only right now
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! workaround for (-1+Nsq) division by zero
  if (Nsq == 1.0_dp) then
    Nsq = Nsq - 1.0e-16_dp
  endif

! angular function derivatives
! from Podolskiy and Vogl, Phys. Rev. B v. 69, p 233101 (2004)
select case (orb_type_i)
  case (ORB_S)
    select case (orb_type_j)
      case (ORB_S)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = 0
              rVd_t = u_Vd_SSS
            end select
        end select
      case (ORB_P)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = u_V_SPS
              rVd(3) = 0
              rVd_t = M*u_Vd_SPS
              case (2)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = u_V_SPS
              rVd_t = N*u_Vd_SPS
              case (3)
              rVd(1) = u_V_SPS
              rVd(2) = 0
              rVd(3) = 0
              rVd_t = L*u_Vd_SPS
            end select
        end select
      case (ORB_D)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = root_3*M*u_V_SDS
              rVd(2) = root_3*L*u_V_SDS
              rVd(3) = 0
              rVd_t = root_3*L*M*u_Vd_SDS
              case (2)
              rVd(1) = 0
              rVd(2) = root_3*N*u_V_SDS
              rVd(3) = root_3*M*u_V_SDS
              rVd_t = root_3*M*N*u_Vd_SDS
              case (3)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = 3*N*u_V_SDS
              rVd_t = ((-1 + 3*Nsq)*u_Vd_SDS)/2.
              case (4)
              rVd(1) = root_3*N*u_V_SDS
              rVd(2) = 0
              rVd(3) = root_3*L*u_V_SDS
              rVd_t = root_3*L*N*u_Vd_SDS
              case (5)
              rVd(1) = 0
              rVd(2) = -2*root_3*M*u_V_SDS
              rVd(3) = -(root_3*N*u_V_SDS)
              rVd_t = -(root_3*(-1 + 2*Msq + Nsq)*u_Vd_SDS)/2.
            end select
        end select
    end select
  case (ORB_P)
    select case (orb_type_j)
      case (ORB_S)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = -u_V_SPS
              rVd(3) = 0
              rVd_t = -(M*u_Vd_SPS)
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = -u_V_SPS
              rVd_t = -(N*u_Vd_SPS)
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              rVd(1) = -u_V_SPS
              rVd(2) = 0
              rVd(3) = 0
              rVd_t = -(L*u_Vd_SPS)
            end select
        end select
      case (ORB_P)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = 2*M*(u_V_PPS - u_V_PPP)
              rVd(3) = 0
              rVd_t = Msq*(u_Vd_PPS - u_Vd_PPP) + u_Vd_PPP
              case (2)
              rVd(1) = 0
              rVd(2) = N*(u_V_PPS - u_V_PPP)
              rVd(3) = M*(u_V_PPS - u_V_PPP)
              rVd_t = M*N*(u_Vd_PPS - u_Vd_PPP)
              case (3)
              rVd(1) = M*(u_V_PPS - u_V_PPP)
              rVd(2) = L*(u_V_PPS - u_V_PPP)
              rVd(3) = 0
              rVd_t = L*M*(u_Vd_PPS - u_Vd_PPP)
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = N*(u_V_PPS - u_V_PPP)
              rVd(3) = M*(u_V_PPS - u_V_PPP)
              rVd_t = M*N*(u_Vd_PPS - u_Vd_PPP)
              case (2)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = 2*N*(u_V_PPS - u_V_PPP)
              rVd_t = Nsq*(u_Vd_PPS - u_Vd_PPP) + u_Vd_PPP
              case (3)
              rVd(1) = N*(u_V_PPS - u_V_PPP)
              rVd(2) = 0
              rVd(3) = L*(u_V_PPS - u_V_PPP)
              rVd_t = L*N*(u_Vd_PPS - u_Vd_PPP)
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              rVd(1) = M*(u_V_PPS - u_V_PPP)
              rVd(2) = L*(u_V_PPS - u_V_PPP)
              rVd(3) = 0
              rVd_t = L*M*(u_Vd_PPS - u_Vd_PPP)
              case (2)
              rVd(1) = N*(u_V_PPS - u_V_PPP)
              rVd(2) = 0
              rVd(3) = L*(u_V_PPS - u_V_PPP)
              rVd_t = L*N*(u_Vd_PPS - u_Vd_PPP)
              case (3)
              rVd(1) = 0
              rVd(2) = 2*M*(-u_V_PPS + u_V_PPP)
              rVd(3) = 2*N*(-u_V_PPS + u_V_PPP)
              rVd_t = -((-1 + Msq + Nsq)*u_Vd_PPS) + (Msq + Nsq)*u_Vd_PPP
            end select
        end select
      case (ORB_D)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = Msq*(root_3*u_V_PDS - 2*u_V_PDP) + u_V_PDP
              rVd(2) = L*M*(2*root_3*u_V_PDS - 4*u_V_PDP)
              rVd(3) = 0
              rVd_t = L*(Msq*(root_3*u_Vd_PDS - 2*u_Vd_PDP) + u_Vd_PDP)
              case (2)
              rVd(1) = 0
              rVd(2) = M*N*(2*root_3*u_V_PDS - 4*u_V_PDP)
              rVd(3) = Msq*(root_3*u_V_PDS - 2*u_V_PDP) + u_V_PDP
              rVd_t = N*(Msq*(root_3*u_Vd_PDS - 2*u_Vd_PDP) + u_Vd_PDP)
              case (3)
              rVd(1) = 0
              rVd(2) = ((-1 + 3*Nsq)*u_V_PDS - 2*root_3*Nsq*u_V_PDP)/2.
              rVd(3) = M*N*(3*u_V_PDS - 2*root_3*u_V_PDP)
              rVd_t = (-(M*u_Vd_PDS) + M*Nsq*(3*u_Vd_PDS - 2*root_3*u_Vd_PDP))/2.
              case (4)
              rVd(1) = M*N*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd(2) = L*N*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd(3) = L*M*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd_t = L*M*N*(root_3*u_Vd_PDS - 2*u_Vd_PDP)
              case (5)
              rVd(1) = 0
              rVd(2) = -(root_3*(-1 + 6*Msq + Nsq)*u_V_PDS)/2. +  &
    (-2 + 6*Msq + Nsq)*u_V_PDP
              rVd(3) = M*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd_t = -(M*(root_3*(-1 + 2*Msq + Nsq)*u_Vd_PDS -  &
         2*(-2 + 2*Msq + Nsq)*u_Vd_PDP))/2.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              rVd(1) = M*N*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd(2) = L*N*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd(3) = L*M*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd_t = L*M*N*(root_3*u_Vd_PDS - 2*u_Vd_PDP)
              case (2)
              rVd(1) = 0
              rVd(2) = Nsq*(root_3*u_V_PDS - 2*u_V_PDP) + u_V_PDP
              rVd(3) = M*N*(2*root_3*u_V_PDS - 4*u_V_PDP)
              rVd_t = M*(Nsq*(root_3*u_Vd_PDS - 2*u_Vd_PDP) + u_Vd_PDP)
              case (3)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = ((-1 + 9*Nsq)*u_V_PDS)/2. + root_3*(1 - 3*Nsq)*u_V_PDP
              rVd_t = (N*((-1 + 3*Nsq)*u_Vd_PDS - 2*root_3*(-1 + Nsq)*u_Vd_PDP))/2.
              case (4)
              rVd(1) = Nsq*(root_3*u_V_PDS - 2*u_V_PDP) + u_V_PDP
              rVd(2) = 0
              rVd(3) = L*N*(2*root_3*u_V_PDS - 4*u_V_PDP)
              rVd_t = L*(Nsq*(root_3*u_Vd_PDS - 2*u_Vd_PDP) + u_Vd_PDP)
              case (5)
              rVd(1) = 0
              rVd(2) = M*N*(-2*root_3*u_V_PDS + 4*u_V_PDP)
              rVd(3) = -((-1 + 2*Msq + 3*Nsq)*(root_3*u_V_PDS - 2*u_V_PDP))/2.
              rVd_t = -(N*(-1 + 2*Msq + Nsq)*(root_3*u_Vd_PDS - 2*u_Vd_PDP))/2.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = -(root_3*(-1 + 3*Msq + Nsq)*u_V_PDS) +  &
    (-1 + 6*Msq + 2*Nsq)*u_V_PDP
              rVd(3) = M*N*(-2*root_3*u_V_PDS + 4*u_V_PDP)
              rVd_t = M*(-(root_3*(-1 + Msq + Nsq)*u_Vd_PDS) +  &
      (-1 + 2*Msq + 2*Nsq)*u_Vd_PDP)
              case (2)
              rVd(1) = M*N*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd(2) = L*N*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd(3) = L*M*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd_t = L*M*N*(root_3*u_Vd_PDS - 2*u_Vd_PDP)
              case (3)
              rVd(1) = ((-1 + 3*Nsq)*u_V_PDS - 2*root_3*Nsq*u_V_PDP)/2.
              rVd(2) = 0
              rVd(3) = L*N*(3*u_V_PDS - 2*root_3*u_V_PDP)
              rVd_t = (-(L*u_Vd_PDS) + L*Nsq*(3*u_Vd_PDS - 2*root_3*u_Vd_PDP))/2.
              case (4)
              rVd(1) = 0
              rVd(2) = M*N*(-2*root_3*u_V_PDS + 4*u_V_PDP)
              rVd(3) = -(root_3*(-1 + Msq + 3*Nsq)*u_V_PDS) +  &
    (-1 + 2*Msq + 6*Nsq)*u_V_PDP
              rVd_t = N*(-(root_3*(-1 + Msq + Nsq)*u_Vd_PDS) +  &
      (-1 + 2*Msq + 2*Nsq)*u_Vd_PDP)
              case (5)
              rVd(1) = (-(root_3*(-1 + 2*Msq + Nsq)*u_V_PDS) +  &
      2*(2*Msq + Nsq)*u_V_PDP)/2.
              rVd(2) = L*M*(-2*root_3*u_V_PDS + 4*u_V_PDP)
              rVd(3) = L*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd_t = -(L*(root_3*(-1 + 2*Msq + Nsq)*u_Vd_PDS -  &
         2*(2*Msq + Nsq)*u_Vd_PDP))/2.
            end select
        end select
    end select
  case (ORB_D)
    select case (orb_type_j)
      case (ORB_S)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = root_3*M*u_V_SDS
              rVd(2) = root_3*L*u_V_SDS
              rVd(3) = 0
              rVd_t = root_3*L*M*u_Vd_SDS
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = root_3*N*u_V_SDS
              rVd(3) = root_3*M*u_V_SDS
              rVd_t = root_3*M*N*u_Vd_SDS
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = 3*N*u_V_SDS
              rVd_t = ((-1 + 3*Nsq)*u_Vd_SDS)/2.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              rVd(1) = root_3*N*u_V_SDS
              rVd(2) = 0
              rVd(3) = root_3*L*u_V_SDS
              rVd_t = root_3*L*N*u_Vd_SDS
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = -2*root_3*M*u_V_SDS
              rVd(3) = -(root_3*N*u_V_SDS)
              rVd_t = -(root_3*(-1 + 2*Msq + Nsq)*u_Vd_SDS)/2.
            end select
        end select
      case (ORB_P)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = -u_V_PDP + Msq*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd(2) = L*M*(-2*root_3*u_V_PDS + 4*u_V_PDP)
              rVd(3) = 0
              rVd_t = -(L*(Msq*(root_3*u_Vd_PDS - 2*u_Vd_PDP) + u_Vd_PDP))
              case (2)
              rVd(1) = M*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd(2) = L*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd(3) = L*M*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd_t = L*M*N*(-(root_3*u_Vd_PDS) + 2*u_Vd_PDP)
              case (3)
              rVd(1) = 0
              rVd(2) = root_3*(-1 + 3*Msq + Nsq)*u_V_PDS +  &
    (1 - 6*Msq - 2*Nsq)*u_V_PDP
              rVd(3) = M*N*(2*root_3*u_V_PDS - 4*u_V_PDP)
              rVd_t = M*(root_3*(-1 + Msq + Nsq)*u_Vd_PDS +  &
      (1 - 2*Msq - 2*Nsq)*u_Vd_PDP)
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = M*N*(-2*root_3*u_V_PDS + 4*u_V_PDP)
              rVd(3) = -u_V_PDP + Msq*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd_t = -(N*(Msq*(root_3*u_Vd_PDS - 2*u_Vd_PDP) + u_Vd_PDP))
              case (2)
              rVd(1) = 0
              rVd(2) = -u_V_PDP + Nsq*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd(3) = M*N*(-2*root_3*u_V_PDS + 4*u_V_PDP)
              rVd_t = -(M*(Nsq*(root_3*u_Vd_PDS - 2*u_Vd_PDP) + u_Vd_PDP))
              case (3)
              rVd(1) = M*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd(2) = L*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd(3) = L*M*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd_t = L*M*N*(-(root_3*u_Vd_PDS) + 2*u_Vd_PDP)
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = ((1 - 3*Nsq)*u_V_PDS + 2*root_3*Nsq*u_V_PDP)/2.
              rVd(3) = M*N*(-3*u_V_PDS + 2*root_3*u_V_PDP)
              rVd_t = (M*((1 - 3*Nsq)*u_Vd_PDS + 2*root_3*Nsq*u_Vd_PDP))/2.
              case (2)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = ((1 - 9*Nsq)*u_V_PDS + 2*root_3*(-1 + 3*Nsq)*u_V_PDP)/2.
              rVd_t = (N*((1 - 3*Nsq)*u_Vd_PDS + 2*root_3*(-1 + Nsq)*u_Vd_PDP))/2.
              case (3)
              rVd(1) = ((1 - 3*Nsq)*u_V_PDS + 2*root_3*Nsq*u_V_PDP)/2.
              rVd(2) = 0
              rVd(3) = L*N*(-3*u_V_PDS + 2*root_3*u_V_PDP)
              rVd_t = (L*((1 - 3*Nsq)*u_Vd_PDS + 2*root_3*Nsq*u_Vd_PDP))/2.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              rVd(1) = M*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd(2) = L*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd(3) = L*M*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd_t = L*M*N*(-(root_3*u_Vd_PDS) + 2*u_Vd_PDP)
              case (2)
              rVd(1) = -u_V_PDP + Nsq*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd(2) = 0
              rVd(3) = L*N*(-2*root_3*u_V_PDS + 4*u_V_PDP)
              rVd_t = -(L*(Nsq*(root_3*u_Vd_PDS - 2*u_Vd_PDP) + u_Vd_PDP))
              case (3)
              rVd(1) = 0
              rVd(2) = M*N*(2*root_3*u_V_PDS - 4*u_V_PDP)
              rVd(3) = root_3*(-1 + Msq + 3*Nsq)*u_V_PDS +  &
    (1 - 2*Msq - 6*Nsq)*u_V_PDP
              rVd_t = N*(root_3*(-1 + Msq + Nsq)*u_Vd_PDS +  &
      (1 - 2*Msq - 2*Nsq)*u_Vd_PDP)
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = (root_3*(-1 + 6*Msq + Nsq)*u_V_PDS -  &
      2*(-2 + 6*Msq + Nsq)*u_V_PDP)/2.
              rVd(3) = M*N*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd_t = (M*(root_3*(-1 + 2*Msq + Nsq)*u_Vd_PDS -  &
        2*(-2 + 2*Msq + Nsq)*u_Vd_PDP))/2.
              case (2)
              rVd(1) = 0
              rVd(2) = M*N*(2*root_3*u_V_PDS - 4*u_V_PDP)
              rVd(3) = ((-1 + 2*Msq + 3*Nsq)*(root_3*u_V_PDS - 2*u_V_PDP))/2.
              rVd_t = (N*(-1 + 2*Msq + Nsq)*(root_3*u_Vd_PDS - 2*u_Vd_PDP))/2.
              case (3)
              rVd(1) = (root_3*(-1 + 2*Msq + Nsq)*u_V_PDS -  &
      2*(2*Msq + Nsq)*u_V_PDP)/2.
              rVd(2) = L*M*(2*root_3*u_V_PDS - 4*u_V_PDP)
              rVd(3) = L*N*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd_t = (L*(root_3*(-1 + 2*Msq + Nsq)*u_Vd_PDS -  &
        2*(2*Msq + Nsq)*u_Vd_PDP))/2.
            end select
        end select
      case (ORB_D)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = -2*M*(-1 + 2*Msq + Nsq)* &
    (3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd(3) = N*(-2*u_V_DDP + 2*u_V_DDD -  &
      2*Msq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd_t = u_Vd_DDP - M**4*(3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD) -  &
    Msq*(-1 + Nsq)*(3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD) +  &
    Nsq*(-u_Vd_DDP + u_Vd_DDD)
              case (2)
              rVd(1) = N*(u_V_DDP - u_V_DDD +  &
      Msq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd(2) = 2*L*M*N*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd(3) = L*(u_V_DDP - u_V_DDD +  &
      Msq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd_t = L*N*(u_Vd_DDP - u_Vd_DDD +  &
      Msq*(3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD))
              case (3)
              rVd(1) = (root_3*M*((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
        Nsq*(-4*u_V_DDP + u_V_DDD)))/2.
              rVd(2) = (root_3*L*((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
        Nsq*(-4*u_V_DDP + u_V_DDD)))/2.
              rVd(3) = root_3*L*M*N*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd_t = (root_3*L*M*((-1 + 3*Nsq)*u_Vd_DDS + u_Vd_DDD +  &
        Nsq*(-4*u_Vd_DDP + u_Vd_DDD)))/2.
              case (4)
              rVd(1) = 0
              rVd(2) = N*(-3*(-1 + 3*Msq + Nsq)*u_V_DDS +  &
      (-3 + 12*Msq + 4*Nsq)*u_V_DDP - (3*Msq + Nsq)*u_V_DDD)
              rVd(3) = M*(-3*(-1 + Msq + 3*Nsq)*u_V_DDS +  &
      (-3 + 4*Msq + 12*Nsq)*u_V_DDP - (Msq + 3*Nsq)*u_V_DDD)
              rVd_t = M*N*(-3*(-1 + Msq + Nsq)*u_Vd_DDS +  &
      (-3 + 4*Msq + 4*Nsq)*u_Vd_DDP - (Msq + Nsq)*u_Vd_DDD)
              case (5)
              rVd(1) = -(M*(-1 + 2*Msq + Nsq)* &
       (3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
              rVd(2) = -(L*(-1 + 6*Msq + Nsq)* &
       (3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
              rVd(3) = -(L*M*N*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd_t = -(L*M*(-1 + 2*Msq + Nsq)* &
       (3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD))/2.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              rVd(1) = N*(u_V_DDP - u_V_DDD +  &
      Msq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd(2) = 2*L*M*N*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd(3) = L*(u_V_DDP - u_V_DDD +  &
      Msq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd_t = L*N*(u_Vd_DDP - u_Vd_DDD +  &
      Msq*(3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD))
              case (2)
              rVd(1) = 0
              rVd(2) = 2*M*(u_V_DDP - u_V_DDD +  &
      Nsq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd(3) = 2*N*(u_V_DDP - u_V_DDD +  &
      Msq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd_t = Nsq*(u_Vd_DDP - u_Vd_DDD) + u_Vd_DDD +  &
    Msq*(u_Vd_DDP - u_Vd_DDD +  &
       Nsq*(3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD))
              case (3)
              rVd(1) = 0
              rVd(2) = (root_3*N*((-1 + 3*Nsq)*u_V_DDS +  &
        (2 - 4*Nsq)*u_V_DDP + (-1 + Nsq)*u_V_DDD))/2.
              rVd(3) = (root_3*M*((-1 + 9*Nsq)*u_V_DDS +  &
        (2 - 12*Nsq)*u_V_DDP + (-1 + 3*Nsq)*u_V_DDD))/2.
              rVd_t = (root_3*M*N*((-1 + 3*Nsq)*u_Vd_DDS + (2 - 4*Nsq)*u_Vd_DDP +  &
        (-1 + Nsq)*u_Vd_DDD))/2.
              case (4)
              rVd(1) = M*(u_V_DDP - u_V_DDD +  &
      Nsq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd(2) = L*(u_V_DDP - u_V_DDD +  &
      Nsq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd(3) = 2*L*M*N*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd_t = L*M*(u_Vd_DDP - u_Vd_DDD +  &
      Nsq*(3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD))
              case (5)
              rVd(1) = 0
              rVd(2) = (N*(-3*(-1 + 6*Msq + Nsq)*u_V_DDS +  &
        (-6 + 24*Msq + 4*Nsq)*u_V_DDP - (-3 + 6*Msq + Nsq)*u_V_DDD))/ &
    2.
              rVd(3) = (M*((3 - 6*Msq - 9*Nsq)*u_V_DDS +  &
        2*(-3 + 4*Msq + 6*Nsq)*u_V_DDP + (3 - 2*Msq - 3*Nsq)*u_V_DDD))/2.
              rVd_t = -(M*N*(3*(-1 + 2*Msq + Nsq)*u_Vd_DDS +  &
         (6 - 8*Msq - 4*Nsq)*u_Vd_DDP + (-3 + 2*Msq + Nsq)*u_Vd_DDD))/ &
    2.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (root_3*M*((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
        Nsq*(-4*u_V_DDP + u_V_DDD)))/2.
              rVd(2) = (root_3*L*((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
        Nsq*(-4*u_V_DDP + u_V_DDD)))/2.
              rVd(3) = root_3*L*M*N*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd_t = (root_3*L*M*((-1 + 3*Nsq)*u_Vd_DDS + u_Vd_DDD +  &
        Nsq*(-4*u_Vd_DDP + u_Vd_DDD)))/2.
              case (2)
              rVd(1) = 0
              rVd(2) = (root_3*N*((-1 + 3*Nsq)*u_V_DDS +  &
        (2 - 4*Nsq)*u_V_DDP + (-1 + Nsq)*u_V_DDD))/2.
              rVd(3) = (root_3*M*((-1 + 9*Nsq)*u_V_DDS +  &
        (2 - 12*Nsq)*u_V_DDP + (-1 + 3*Nsq)*u_V_DDD))/2.
              rVd_t = (root_3*M*N*((-1 + 3*Nsq)*u_Vd_DDS + (2 - 4*Nsq)*u_Vd_DDP +  &
        (-1 + Nsq)*u_Vd_DDD))/2.
              case (3)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = N*((-3 + 9*Nsq)*u_V_DDS + (6 - 12*Nsq)*u_V_DDP +  &
      3*(-1 + Nsq)*u_V_DDD)
              rVd_t = ((1 - 3*Nsq)**2*u_Vd_DDS -  &
      3*(-1 + Nsq)*(Nsq*(4*u_Vd_DDP - u_Vd_DDD) + u_Vd_DDD))/4.
              case (4)
              rVd(1) = (root_3*N*((-1 + 3*Nsq)*u_V_DDS +  &
        (2 - 4*Nsq)*u_V_DDP + (-1 + Nsq)*u_V_DDD))/2.
              rVd(2) = 0
              rVd(3) = (root_3*L*((-1 + 9*Nsq)*u_V_DDS +  &
        (2 - 12*Nsq)*u_V_DDP + (-1 + 3*Nsq)*u_V_DDD))/2.
              rVd_t = (root_3*L*N*((-1 + 3*Nsq)*u_Vd_DDS + (2 - 4*Nsq)*u_Vd_DDP +  &
        (-1 + Nsq)*u_Vd_DDD))/2.
              case (5)
              rVd(1) = 0
              rVd(2) = -(root_3*M*((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
        Nsq*(-4*u_V_DDP + u_V_DDD)))
              rVd(3) = -(root_3*N*((-2 + 3*Msq + 3*Nsq)*u_V_DDS +  &
        (2 - 4*Msq - 4*Nsq)*u_V_DDP + (Msq + Nsq)*u_V_DDD))
              rVd_t = -(root_3*(-1 + 2*Msq + Nsq)* &
       ((-1 + 3*Nsq)*u_Vd_DDS + u_Vd_DDD +  &
         Nsq*(-4*u_Vd_DDP + u_Vd_DDD)))/4.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = N*(-3*(-1 + 3*Msq + Nsq)*u_V_DDS +  &
      (-3 + 12*Msq + 4*Nsq)*u_V_DDP - (3*Msq + Nsq)*u_V_DDD)
              rVd(3) = M*(-3*(-1 + Msq + 3*Nsq)*u_V_DDS +  &
      (-3 + 4*Msq + 12*Nsq)*u_V_DDP - (Msq + 3*Nsq)*u_V_DDD)
              rVd_t = M*N*(-3*(-1 + Msq + Nsq)*u_Vd_DDS +  &
      (-3 + 4*Msq + 4*Nsq)*u_Vd_DDP - (Msq + Nsq)*u_Vd_DDD)
              case (2)
              rVd(1) = M*(u_V_DDP - u_V_DDD +  &
      Nsq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd(2) = L*(u_V_DDP - u_V_DDD +  &
      Nsq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd(3) = 2*L*M*N*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd_t = L*M*(u_Vd_DDP - u_Vd_DDD +  &
      Nsq*(3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD))
              case (3)
              rVd(1) = (root_3*N*((-1 + 3*Nsq)*u_V_DDS +  &
        (2 - 4*Nsq)*u_V_DDP + (-1 + Nsq)*u_V_DDD))/2.
              rVd(2) = 0
              rVd(3) = (root_3*L*((-1 + 9*Nsq)*u_V_DDS +  &
        (2 - 12*Nsq)*u_V_DDP + (-1 + 3*Nsq)*u_V_DDD))/2.
              rVd_t = (root_3*L*N*((-1 + 3*Nsq)*u_Vd_DDS + (2 - 4*Nsq)*u_Vd_DDP +  &
        (-1 + Nsq)*u_Vd_DDD))/2.
              case (4)
              rVd(1) = 0
              rVd(2) = M*(-2*u_V_DDP + 2*u_V_DDD -  &
      2*Nsq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd(3) = -2*N*(-1 + Msq + 2*Nsq)* &
    (3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd_t = u_Vd_DDP - (-1 + Msq)*Nsq* &
     (3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD) -  &
    N**4*(3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD) +  &
    Msq*(-u_Vd_DDP + u_Vd_DDD)
              case (5)
              rVd(1) = (N*(-3*(-1 + 2*Msq + Nsq)*u_V_DDS +  &
        (-2 + 8*Msq + 4*Nsq)*u_V_DDP - (1 + 2*Msq + Nsq)*u_V_DDD))/2.
              rVd(2) = -2*L*M*N*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd(3) = (L*((3 - 6*Msq - 9*Nsq)*u_V_DDS +  &
        2*(-1 + 4*Msq + 6*Nsq)*u_V_DDP - (1 + 2*Msq + 3*Nsq)*u_V_DDD))/2.
              rVd_t = -(L*N*(3*(-1 + 2*Msq + Nsq)*u_Vd_DDS +  &
         (2 - 8*Msq - 4*Nsq)*u_Vd_DDP + (1 + 2*Msq + Nsq)*u_Vd_DDD))/2.
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              rVd(1) = -(M*(-1 + 2*Msq + Nsq)* &
       (3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
              rVd(2) = -(L*(-1 + 6*Msq + Nsq)* &
       (3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
              rVd(3) = -(L*M*N*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd_t = -(L*M*(-1 + 2*Msq + Nsq)* &
       (3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD))/2.
              case (2)
              rVd(1) = 0
              rVd(2) = (N*(-3*(-1 + 6*Msq + Nsq)*u_V_DDS +  &
        (-6 + 24*Msq + 4*Nsq)*u_V_DDP - (-3 + 6*Msq + Nsq)*u_V_DDD))/ &
    2.
              rVd(3) = (M*((3 - 6*Msq - 9*Nsq)*u_V_DDS +  &
        2*(-3 + 4*Msq + 6*Nsq)*u_V_DDP + (3 - 2*Msq - 3*Nsq)*u_V_DDD))/2.
              rVd_t = -(M*N*(3*(-1 + 2*Msq + Nsq)*u_Vd_DDS +  &
         (6 - 8*Msq - 4*Nsq)*u_Vd_DDP + (-3 + 2*Msq + Nsq)*u_Vd_DDD))/ &
    2.
              case (3)
              rVd(1) = 0
              rVd(2) = -(root_3*M*((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
        Nsq*(-4*u_V_DDP + u_V_DDD)))
              rVd(3) = -(root_3*N*((-2 + 3*Msq + 3*Nsq)*u_V_DDS +  &
        (2 - 4*Msq - 4*Nsq)*u_V_DDP + (Msq + Nsq)*u_V_DDD))
              rVd_t = -(root_3*(-1 + 2*Msq + Nsq)* &
       ((-1 + 3*Nsq)*u_Vd_DDS + u_Vd_DDD +  &
         Nsq*(-4*u_Vd_DDP + u_Vd_DDD)))/4.
              case (4)
              rVd(1) = (N*(-3*(-1 + 2*Msq + Nsq)*u_V_DDS +  &
        (-2 + 8*Msq + 4*Nsq)*u_V_DDP - (1 + 2*Msq + Nsq)*u_V_DDD))/2.
              rVd(2) = -2*L*M*N*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd(3) = (L*((3 - 6*Msq - 9*Nsq)*u_V_DDS +  &
        2*(-1 + 4*Msq + 6*Nsq)*u_V_DDP - (1 + 2*Msq + 3*Nsq)*u_V_DDD))/2.
              rVd_t = -(L*N*(3*(-1 + 2*Msq + Nsq)*u_Vd_DDS +  &
         (2 - 8*Msq - 4*Nsq)*u_Vd_DDP + (1 + 2*Msq + Nsq)*u_Vd_DDD))/2.
              case (5)
              rVd(1) = 0
              rVd(2) = 2*M*(-1 + 2*Msq + Nsq)* &
    (3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd(3) = N*(3*(-1 + 2*Msq + Nsq)*u_V_DDS +  &
      (2 - 8*Msq - 4*Nsq)*u_V_DDP + (1 + 2*Msq + Nsq)*u_V_DDD)
              rVd_t = (3*(-1 + 2*Msq + Nsq)**2*u_Vd_DDS -  &
      4*Msq*(-1 + Nsq)*(4*u_Vd_DDP - u_Vd_DDD) + u_Vd_DDD +  &
      4*M**4*(-4*u_Vd_DDP + u_Vd_DDD) +  &
      N**4*(-4*u_Vd_DDP + u_Vd_DDD) + 2*Nsq*(2*u_Vd_DDP + u_Vd_DDD))/ &
    4.
            end select
        end select
    end select
end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end SKd_vogl.h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
