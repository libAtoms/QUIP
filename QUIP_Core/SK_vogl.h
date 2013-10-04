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
              V = (root_3*(L - M)*(L + M)*u_V_SDS)/2.
            end select
        end select
      case (ORB_F)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = -(Sqrt(2.5)*M*(-3*Lsq + Msq)*u_V_SFS)/2.
              case (2)
              V = sqrt(15.)*L*M*N*u_V_SFS
              case (3)
              V = (Sqrt(1.5)*M*(-1 + 5*Nsq)*u_V_SFS)/2.
              case (4)
              V = (N*(-3 + 5*Nsq)*u_V_SFS)/2.
              case (5)
              V = (Sqrt(1.5)*L*(-1 + 5*Nsq)*u_V_SFS)/2.
              case (6)
              V = (sqrt(15.)*(L - M)*(L + M)*N*u_V_SFS)/2.
              case (7)
              V = (Sqrt(2.5)*(L**3 - 3*L*Msq)*u_V_SFS)/2.
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
              V = Lsq*u_V_PPS + (Msq + Nsq)*u_V_PPP
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
              V = (M*((-1 + 3*Nsq)*u_V_PDS - 2*root_3*Nsq*u_V_PDP))/2.
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
              V = ((L - M)*(L + M)*N*(root_3*u_V_PDS - 2*u_V_PDP))/2.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              V = root_3*Lsq*M*u_V_PDS + M*(-1 + 2*Msq + 2*Nsq)*u_V_PDP
              case (2)
              V = L*M*N*(root_3*u_V_PDS - 2*u_V_PDP)
              case (3)
              V = (L*((-1 + 3*Nsq)*u_V_PDS - 2*root_3*Nsq*u_V_PDP))/2.
              case (4)
              V = N*(root_3*Lsq*u_V_PDS + (-1 + 2*Msq + 2*Nsq)*u_V_PDP)
              case (5)
              V = -(L*(root_3*(-1 + 2*Msq + Nsq)*u_V_PDS -  &
         2*(2*Msq + Nsq)*u_V_PDP))/2.
            end select
        end select
      case (ORB_F)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = (sqrt(5.)*(-(root_3*(-1 + Nsq)*u_V_PFP) +  &
        M**4*(-4*sqrt(2.)*u_V_PFS + 4*root_3*u_V_PFP) +  &
        Msq*(-3*sqrt(2.)*(-1 + Nsq)*u_V_PFS +  &
           root_3*(-5 + 3*Nsq)*u_V_PFP)))/4.
              case (2)
              V = (L*N*(sqrt(10.)*u_V_PFP +  &
        Msq*(2*sqrt(15.)*u_V_PFS - 3*sqrt(10.)*u_V_PFP)))/2.
              case (3)
              V = ((-1 + 5*Nsq)*u_V_PFP +  &
      Msq*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS + (1 - 15*Nsq)*u_V_PFP))/4.
              case (4)
              V = (M*N*(2*(-3 + 5*Nsq)*u_V_PFS +  &
        sqrt(6.)*(1 - 5*Nsq)*u_V_PFP))/4.
              case (5)
              V = (L*M*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS + (1 - 15*Nsq)*u_V_PFP))/ &
    4.
              case (6)
              V = (M*N*(-2*sqrt(15.)*(-1 + 2*Msq + Nsq)*u_V_PFS +  &
        sqrt(10.)*(-5 + 6*Msq + 3*Nsq)*u_V_PFP))/4.
              case (7)
              V = (L*M*(-(sqrt(10.)*(-1 + 4*Msq + Nsq)*u_V_PFS) +  &
        sqrt(15.)*(-3 + 4*Msq + Nsq)*u_V_PFP))/4.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              V = -(M*N*(-3 + 4*Msq + 3*Nsq)* &
       (sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/4.
              case (2)
              V = (L*M*(2*sqrt(15.)*Nsq*u_V_PFS +  &
        sqrt(10.)*(1 - 3*Nsq)*u_V_PFP))/2.
              case (3)
              V = (M*N*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (11 - 15*Nsq)*u_V_PFP))/4.
              case (4)
              V = (2*Nsq*(-3 + 5*Nsq)*u_V_PFS +  &
      sqrt(6.)*(-1 + 6*Nsq - 5*N**4)*u_V_PFP)/4.
              case (5)
              V = (L*N*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (11 - 15*Nsq)*u_V_PFP))/4.
              case (6)
              V = (sqrt(5.)*(L - M)*(L + M)* &
      (2*root_3*Nsq*u_V_PFS + sqrt(2.)*(1 - 3*Nsq)*u_V_PFP))/4.
              case (7)
              V = -(L*N*(-1 + 4*Msq + Nsq)* &
       (sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/4.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              V = (L*M*(sqrt(10.)*(3 - 4*Msq - 3*Nsq)*u_V_PFS +  &
        sqrt(15.)*(-1 + 4*Msq + 3*Nsq)*u_V_PFP))/4.
              case (2)
              V = -(sqrt(5.)*M*N*(sqrt(2.)*Msq*u_V_PFP +  &
         Lsq*(-2*root_3*(-1 + Nsq)*u_V_PFS +  &
            sqrt(2.)*(-2 + 3*Nsq)*u_V_PFP)))/(2.*(-1 + Nsq))
              case (3)
              V = (L*M*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS + (1 - 15*Nsq)*u_V_PFP))/ &
    4.
              case (4)
              V = (L*N*(2*(-3 + 5*Nsq)*u_V_PFS +  &
        sqrt(6.)*(1 - 5*Nsq)*u_V_PFP))/4.
              case (5)
              V = (sqrt(6.)*Lsq*(1 - 6*Nsq + 5*N**4)*u_V_PFS +  &
      (Lsq*Nsq*(11 - 15*Nsq) + Msq*(1 - 5*Nsq))*u_V_PFP)/ &
    (4.*(-1 + Nsq))
              case (6)
              V = (L*N*(-2*sqrt(15.)*(-1 + 2*Msq + Nsq)*u_V_PFS +  &
        sqrt(10.)*(-1 + 6*Msq + 3*Nsq)*u_V_PFP))/4.
              case (7)
              V = (sqrt(5.)*(root_3*M**4*u_V_PFP +  &
        3*Lsq*Msq*(-1 + Nsq)* &
         (-(sqrt(2.)*u_V_PFS) + root_3*u_V_PFP) +  &
        L**4*(sqrt(2.)*(-1 + Nsq)*u_V_PFS - root_3*Nsq*u_V_PFP)))/ &
    (4.*(-1 + Nsq))
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
              V = (root_3*(L - M)*(L + M)*u_V_SDS)/2.
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
              V = -(root_3*Lsq*M*u_V_PDS) + M*(1 - 2*Msq - 2*Nsq)*u_V_PDP
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
              V = (N*(Msq*u_V_PDP -  &
        Lsq*(root_3*(-1 + Nsq)*u_V_PDS + (1 - 2*Nsq)*u_V_PDP)))/ &
    (-1 + Nsq)
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              V = (M*(root_3*(-1 + 2*Msq + Nsq)*u_V_PDS -  &
        2*(-2 + 2*Msq + Nsq)*u_V_PDP))/2.
              case (2)
              V = -((L - M)*(L + M)*N*(root_3*u_V_PDS - 2*u_V_PDP))/2.
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
              V = -(M*N*(3*(-1 + Msq + Nsq)*u_V_DDS +  &
        (3 - 4*Msq - 4*Nsq)*u_V_DDP + (Msq + Nsq)*u_V_DDD))
              case (5)
              V = (L*(L - M)*M*(L + M)*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
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
              V = ((1 - 3*Nsq)**2*u_V_DDS +  &
      3*(-1 + Nsq)*(-4*Nsq*u_V_DDP + (-1 + Nsq)*u_V_DDD))/4.
              case (4)
              V = (root_3*L*N*((-1 + 3*Nsq)*u_V_DDS + (2 - 4*Nsq)*u_V_DDP +  &
        (-1 + Nsq)*u_V_DDD))/2.
              case (5)
              V = (root_3*(L - M)*(L + M)* &
      ((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
        Nsq*(-4*u_V_DDP + u_V_DDD)))/4.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              V = -(M*N*(3*(-1 + Msq + Nsq)*u_V_DDS +  &
        (3 - 4*Msq - 4*Nsq)*u_V_DDP + (Msq + Nsq)*u_V_DDD))
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
              V = (L*(L - M)*M*(L + M)*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
              case (2)
              V = -(M*N*(3*(-1 + 2*Msq + Nsq)*u_V_DDS +  &
         (6 - 8*Msq - 4*Nsq)*u_V_DDP + (-3 + 2*Msq + Nsq)*u_V_DDD))/ &
    2.
              case (3)
              V = (root_3*(L - M)*(L + M)* &
      ((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
        Nsq*(-4*u_V_DDP + u_V_DDD)))/4.
              case (4)
              V = -(L*N*(3*(-1 + 2*Msq + Nsq)*u_V_DDS +  &
         (2 - 8*Msq - 4*Nsq)*u_V_DDP + (1 + 2*Msq + Nsq)*u_V_DDD))/2.
              case (5)
              V = (3*(Lsq - Msq)**2*u_V_DDS -  &
      4*(4*M**4 + 4*Msq*(-1 + Nsq) + Nsq*(-1 + Nsq))*u_V_DDP +  &
      (4*M**4 + 4*Msq*(-1 + Nsq) + (1 + Nsq)**2)*u_V_DDD)/4.
            end select
        end select
      case (ORB_F)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = (root_3*L*(-(sqrt(10.)*Msq*(-3*Lsq + Msq)*u_V_DFS) +  &
        sqrt(5.)*(1 + 8*M**4 - Nsq + 6*Msq*(-1 + Nsq))*u_V_DFP -  &
        sqrt(2.)*(4*M**4 - 2*Nsq + 3*Msq*(-1 + Nsq))*u_V_DFD))/4.
              case (2)
              V = (N*(-(sqrt(10.)*(-1 + Nsq)*u_V_DFP) +  &
        2*(-1 + 2*Nsq)*u_V_DFD -  &
        6*M**4*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD) -  &
        6*Msq*(-1 + Nsq)*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP +  &
           u_V_DFD)))/2.
              case (3)
              V = (L*((-1 + 5*Nsq)*u_V_DFP - 2*sqrt(10.)*Nsq*u_V_DFD +  &
        Msq*(3*sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS + (2 - 30*Nsq)*u_V_DFP +  &
           sqrt(10.)*(1 + 3*Nsq)*u_V_DFD)))/4.
              case (4)
              V = (root_3*L*M*N*((-3 + 5*Nsq)*u_V_DFS +  &
        sqrt(2.)*(1 - 5*Nsq)*u_V_DFP + sqrt(5.)*(1 + Nsq)*u_V_DFD))/2.
              case (5)
              V = (M*(Msq*((1 - 5*Nsq)*u_V_DFP + 2*sqrt(10.)*Nsq*u_V_DFD) +  &
        Lsq*(3*sqrt(2.)*(1 - 6*Nsq + 5*N**4)*u_V_DFS +  &
           (-1 + 27*Nsq - 30*N**4)*u_V_DFP +  &
           sqrt(10.)*(-1 + 3*N**4)*u_V_DFD)))/(4.*(-1 + Nsq))
              case (6)
              V = (3*L*(L - M)*M*(L + M)*N* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              case (7)
              V = (root_3*M*(sqrt(10.)*Lsq*(Lsq - 3*Msq)*u_V_DFS -  &
        sqrt(5.)*(3 + 8*M**4 - 5*Nsq + 2*N**4 + 10*Msq*(-1 + Nsq))* &
         u_V_DFP + sqrt(2.)*(1 + 4*M**4 - 4*Nsq + N**4 +  &
           5*Msq*(-1 + Nsq))*u_V_DFD))/4.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              V = (N*(-((-1 + Nsq)*(sqrt(15.)*u_V_DFP - 2*sqrt(6.)*u_V_DFD)) -  &
        4*M**4*(sqrt(30.)*u_V_DFS - 2*sqrt(15.)*u_V_DFP +  &
           sqrt(6.)*u_V_DFD) +  &
        Msq*(-3*sqrt(30.)*(-1 + Nsq)*u_V_DFS +  &
           2*sqrt(15.)*(-4 + 3*Nsq)*u_V_DFP +  &
           sqrt(6.)*(7 - 3*Nsq)*u_V_DFD)))/4.
              case (2)
              V = (L*(Nsq*(sqrt(10.)*u_V_DFP - 4*u_V_DFD) + 2*u_V_DFD +  &
        Msq*(sqrt(10.)*u_V_DFP - 4*u_V_DFD +  &
           6*Nsq*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))))/2.
              case (3)
              V = (N*((-1 + 5*Nsq)*u_V_DFP - 2*sqrt(10.)*(-1 + Nsq)*u_V_DFD +  &
        3*Msq*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS + (4 - 10*Nsq)*u_V_DFP +  &
           sqrt(10.)*(-1 + Nsq)*u_V_DFD)))/4.
              case (4)
              V = (root_3*M*(sqrt(2.)*(-1 + 7*Nsq - 10*N**4)*u_V_DFP +  &
        2*Nsq*((-3 + 5*Nsq)*u_V_DFS + sqrt(5.)*(-1 + Nsq)*u_V_DFD)))/4.
              case (5)
              V = (3*L*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/4.
              case (6)
              V = (M*(sqrt(10.)*(1 - 2*Msq)*u_V_DFP + 8*(-1 + Msq)*u_V_DFD -  &
        6*N**4*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD) -  &
        3*Nsq*(2*sqrt(5.)*(-1 + 2*Msq)*u_V_DFS +  &
           sqrt(10.)*(3 - 4*Msq)*u_V_DFP + 2*(-3 + 2*Msq)*u_V_DFD)))/4.
              case (7)
              V = (L*M*N*(-(sqrt(30.)*(-1 + 4*Msq + Nsq)*u_V_DFS) +  &
        2*sqrt(15.)*(-2 + 4*Msq + Nsq)*u_V_DFP -  &
        sqrt(6.)*(-5 + 4*Msq + Nsq)*u_V_DFD))/4.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              V = ((3*Lsq*M - M**3)*(sqrt(10.)*(-1 + 3*Nsq)*u_V_DFS -  &
        6*sqrt(5.)*Nsq*u_V_DFP + 3*sqrt(2.)*(1 + Nsq)*u_V_DFD))/8.
              case (2)
              V = (root_3*L*M*N*(-1 + 3*Nsq)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              case (3)
              V = (root_3*M*(sqrt(2.)*(1 - 8*Nsq + 15*N**4)*u_V_DFS +  &
        2*Nsq*(11 - 15*Nsq)*u_V_DFP +  &
        sqrt(10.)*(1 - 4*Nsq + 3*N**4)*u_V_DFD))/8.
              case (4)
              V = (N*((3 - 14*Nsq + 15*N**4)*u_V_DFS -  &
        3*sqrt(2.)*(1 - 6*Nsq + 5*N**4)*u_V_DFP +  &
        3*sqrt(5.)*(-1 + Nsq)**2*u_V_DFD))/4.
              case (5)
              V = (root_3*L*(sqrt(2.)*(1 - 8*Nsq + 15*N**4)*u_V_DFS +  &
        2*Nsq*(11 - 15*Nsq)*u_V_DFP +  &
        sqrt(10.)*(1 - 4*Nsq + 3*N**4)*u_V_DFD))/8.
              case (6)
              V = (root_3*(L - M)*(L + M)*N*(-1 + 3*Nsq)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/4.
              case (7)
              V = ((L**3 - 3*L*Msq)*(sqrt(10.)*(-1 + 3*Nsq)*u_V_DFS -  &
        6*sqrt(5.)*Nsq*u_V_DFP + 3*sqrt(2.)*(1 + Nsq)*u_V_DFD))/8.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              V = (L*M*N*(sqrt(30.)*(3 - 4*Msq - 3*Nsq)*u_V_DFS +  &
        2*sqrt(15.)*(-2 + 4*Msq + 3*Nsq)*u_V_DFP -  &
        sqrt(6.)*(1 + 4*Msq + 3*Nsq)*u_V_DFD))/4.
              case (2)
              V = (Lsq*M*(sqrt(10.)*(-1 + 6*Nsq - 6*N**4)*u_V_DFP +  &
         2*u_V_DFD + 6*Nsq*(-1 + Nsq)*(sqrt(5.)*u_V_DFS + u_V_DFD))  &
 + M**3*(-2*u_V_DFD + Nsq*(-(sqrt(10.)*u_V_DFP) + 4*u_V_DFD)))/ &
    (2.*(-1 + Nsq))
              case (3)
              V = (3*L*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/4.
              case (4)
              V = (root_3*L*(sqrt(2.)*(-1 + 7*Nsq - 10*N**4)*u_V_DFP +  &
        2*Nsq*((-3 + 5*Nsq)*u_V_DFS + sqrt(5.)*(-1 + Nsq)*u_V_DFD)))/4.
              case (5)
              V = -(N*(3*sqrt(2.)*(-1 + Msq + Nsq)*(-1 + 5*Nsq)*u_V_DFS +  &
         (-11 + 37*Nsq - 30*N**4 - 6*Msq*(-2 + 5*Nsq))*u_V_DFP +  &
         sqrt(10.)*(-1 + Nsq)*(-1 + 3*Msq + 3*Nsq)*u_V_DFD))/4.
              case (6)
              V = (L*(sqrt(10.)*(1 - 2*Msq)*u_V_DFP + 8*Msq*u_V_DFD -  &
        6*N**4*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD) +  &
        Nsq*(6*sqrt(5.)*(1 - 2*Msq)*u_V_DFS +  &
           sqrt(10.)*(-5 + 12*Msq)*u_V_DFP + 2*(1 - 6*Msq)*u_V_DFD)))/4.
              case (7)
              V = (root_3*N*(sqrt(10.)* &
         (4*M**4 + 5*Msq*(-1 + Nsq) + (-1 + Nsq)**2)*u_V_DFS -  &
        sqrt(5.)*(1 + 8*M**4 - 3*Nsq + 2*N**4 + 2*Msq*(-4 + 5*Nsq))* &
         u_V_DFP + sqrt(2.)*(-1 + 4*M**4 + N**4 + Msq*(-1 + 5*Nsq))* &
         u_V_DFD))/4.
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              V = (root_3*M*(sqrt(10.)*(3*L**4 - 4*Lsq*Msq + M**4)*u_V_DFS -  &
        2*sqrt(5.)*(2 + 8*M**4 - 5*Nsq + 3*N**4 + 10*Msq*(-1 + Nsq))* &
         u_V_DFP + sqrt(2.)*(3 + 8*M**4 - 2*Nsq + 3*N**4 +  &
           10*Msq*(-1 + Nsq))*u_V_DFD))/8.
              case (2)
              V = (3*L*(L - M)*M*(L + M)*N* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              case (3)
              V = -(M*(3*sqrt(2.)*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)*u_V_DFS +  &
         (-4 + 42*Nsq - 30*N**4 + Msq*(4 - 60*Nsq))*u_V_DFP +  &
         sqrt(10.)*(-1 - 6*Nsq + 3*N**4 + Msq*(2 + 6*Nsq))*u_V_DFD))/8.
              case (4)
              V = (root_3*(L - M)*(L + M)*N* &
      ((-3 + 5*Nsq)*u_V_DFS + sqrt(2.)*(1 - 5*Nsq)*u_V_DFP +  &
        sqrt(5.)*(1 + Nsq)*u_V_DFD))/4.
              case (5)
              V = (L*(3*sqrt(2.)*(L - M)*(L + M)*(-1 + 5*Nsq)*u_V_DFS +  &
        2*(Nsq*(-11 + 15*Nsq) + Msq*(-2 + 30*Nsq))*u_V_DFP -  &
        sqrt(10.)*(-1 + 2*Nsq + 3*N**4 + Msq*(2 + 6*Nsq))*u_V_DFD))/8.
              case (6)
              V = (N*(3*sqrt(5.)*(Lsq - Msq)**2*u_V_DFS -  &
        sqrt(10.)*(1 + 12*M**4 - 4*Nsq + 3*N**4 + 12*Msq*(-1 + Nsq))* &
         u_V_DFP + (-1 + 12*M**4 + 2*Nsq + 3*N**4 + 12*Msq*(-1 + Nsq))* &
         u_V_DFD))/4.
              case (7)
              V = (root_3*L*(sqrt(10.)*(L**4 - 4*Lsq*Msq + 3*M**4)*u_V_DFS -  &
        2*sqrt(5.)*(8*M**4 + 6*Msq*(-1 + Nsq) + Nsq*(-1 + Nsq))* &
         u_V_DFP + sqrt(2.)*(8*M**4 + 6*Msq*(-1 + Nsq) + (1 + Nsq)**2)* &
         u_V_DFD))/8.
            end select
        end select
    end select
  case (ORB_F)
    select case (orb_type_j)
      case (ORB_S)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = (Sqrt(2.5)*M*(-3*Lsq + Msq)*u_V_SFS)/2.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              V = -(sqrt(15.)*L*M*N*u_V_SFS)
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              V = (Sqrt(1.5)*M*(1 - 5*Nsq)*u_V_SFS)/2.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              V = (N*(3 - 5*Nsq)*u_V_SFS)/2.
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              V = (Sqrt(1.5)*L*(1 - 5*Nsq)*u_V_SFS)/2.
            end select
          case (6)
            select case (orb_dir_j)
              case (1)
              V = -(sqrt(15.)*(L - M)*(L + M)*N*u_V_SFS)/2.
            end select
          case (7)
            select case (orb_dir_j)
              case (1)
              V = (Sqrt(2.5)*L*(-1 + 4*Msq + Nsq)*u_V_SFS)/2.
            end select
        end select
      case (ORB_P)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = (sqrt(5.)*(-(root_3*(-1 + Nsq)*u_V_PFP) +  &
        M**4*(-4*sqrt(2.)*u_V_PFS + 4*root_3*u_V_PFP) +  &
        Msq*(-3*sqrt(2.)*(-1 + Nsq)*u_V_PFS +  &
           root_3*(-5 + 3*Nsq)*u_V_PFP)))/4.
              case (2)
              V = -(M*N*(-3 + 4*Msq + 3*Nsq)* &
       (sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/4.
              case (3)
              V = (L*M*(sqrt(10.)*(3 - 4*Msq - 3*Nsq)*u_V_PFS +  &
        sqrt(15.)*(-1 + 4*Msq + 3*Nsq)*u_V_PFP))/4.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              V = (L*N*(sqrt(10.)*u_V_PFP +  &
        Msq*(2*sqrt(15.)*u_V_PFS - 3*sqrt(10.)*u_V_PFP)))/2.
              case (2)
              V = (L*M*(2*sqrt(15.)*Nsq*u_V_PFS +  &
        sqrt(10.)*(1 - 3*Nsq)*u_V_PFP))/2.
              case (3)
              V = -(sqrt(5.)*M*N*(sqrt(2.)*Msq*u_V_PFP +  &
         Lsq*(-2*root_3*(-1 + Nsq)*u_V_PFS +  &
            sqrt(2.)*(-2 + 3*Nsq)*u_V_PFP)))/(2.*(-1 + Nsq))
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              V = ((-1 + 5*Nsq)*u_V_PFP +  &
      Msq*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS + (1 - 15*Nsq)*u_V_PFP))/4.
              case (2)
              V = (M*N*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (11 - 15*Nsq)*u_V_PFP))/4.
              case (3)
              V = (L*M*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS + (1 - 15*Nsq)*u_V_PFP))/ &
    4.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              V = (M*N*(2*(-3 + 5*Nsq)*u_V_PFS +  &
        sqrt(6.)*(1 - 5*Nsq)*u_V_PFP))/4.
              case (2)
              V = (2*Nsq*(-3 + 5*Nsq)*u_V_PFS +  &
      sqrt(6.)*(-1 + 6*Nsq - 5*N**4)*u_V_PFP)/4.
              case (3)
              V = (L*N*(2*(-3 + 5*Nsq)*u_V_PFS +  &
        sqrt(6.)*(1 - 5*Nsq)*u_V_PFP))/4.
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              V = (L*M*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS + (1 - 15*Nsq)*u_V_PFP))/ &
    4.
              case (2)
              V = (L*N*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (11 - 15*Nsq)*u_V_PFP))/4.
              case (3)
              V = (sqrt(6.)*Lsq*(1 - 6*Nsq + 5*N**4)*u_V_PFS +  &
      (Lsq*Nsq*(11 - 15*Nsq) + Msq*(1 - 5*Nsq))*u_V_PFP)/ &
    (4.*(-1 + Nsq))
            end select
          case (6)
            select case (orb_dir_j)
              case (1)
              V = (M*N*(-2*sqrt(15.)*(-1 + 2*Msq + Nsq)*u_V_PFS +  &
        sqrt(10.)*(-5 + 6*Msq + 3*Nsq)*u_V_PFP))/4.
              case (2)
              V = (sqrt(5.)*(L - M)*(L + M)* &
      (2*root_3*Nsq*u_V_PFS + sqrt(2.)*(1 - 3*Nsq)*u_V_PFP))/4.
              case (3)
              V = (L*N*(-2*sqrt(15.)*(-1 + 2*Msq + Nsq)*u_V_PFS +  &
        sqrt(10.)*(-1 + 6*Msq + 3*Nsq)*u_V_PFP))/4.
            end select
          case (7)
            select case (orb_dir_j)
              case (1)
              V = (L*M*(-(sqrt(10.)*(-1 + 4*Msq + Nsq)*u_V_PFS) +  &
        sqrt(15.)*(-3 + 4*Msq + Nsq)*u_V_PFP))/4.
              case (2)
              V = -(L*N*(-1 + 4*Msq + Nsq)* &
       (sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/4.
              case (3)
              V = (sqrt(5.)*(root_3*M**4*u_V_PFP +  &
        3*Lsq*Msq*(-1 + Nsq)* &
         (-(sqrt(2.)*u_V_PFS) + root_3*u_V_PFP) +  &
        L**4*(sqrt(2.)*(-1 + Nsq)*u_V_PFS - root_3*Nsq*u_V_PFP)))/ &
    (4.*(-1 + Nsq))
            end select
        end select
      case (ORB_D)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = -(root_3*L*(-(sqrt(10.)*Msq*(-3*Lsq + Msq)*u_V_DFS) +  &
         sqrt(5.)*(1 + 8*M**4 - Nsq + 6*Msq*(-1 + Nsq))*u_V_DFP -  &
         sqrt(2.)*(4*M**4 - 2*Nsq + 3*Msq*(-1 + Nsq))*u_V_DFD))/4.
              case (2)
              V = (N*((-1 + Nsq)*(sqrt(15.)*u_V_DFP - 2*sqrt(6.)*u_V_DFD) +  &
        4*M**4*(sqrt(30.)*u_V_DFS - 2*sqrt(15.)*u_V_DFP +  &
           sqrt(6.)*u_V_DFD) +  &
        Msq*(3*sqrt(30.)*(-1 + Nsq)*u_V_DFS +  &
           2*sqrt(15.)*(4 - 3*Nsq)*u_V_DFP +  &
           sqrt(6.)*(-7 + 3*Nsq)*u_V_DFD)))/4.
              case (3)
              V = -((3*Lsq*M - M**3)*(sqrt(10.)*(-1 + 3*Nsq)*u_V_DFS -  &
         6*sqrt(5.)*Nsq*u_V_DFP + 3*sqrt(2.)*(1 + Nsq)*u_V_DFD))/8.
              case (4)
              V = (L*M*N*(sqrt(30.)*(-3 + 4*Msq + 3*Nsq)*u_V_DFS -  &
        2*sqrt(15.)*(-2 + 4*Msq + 3*Nsq)*u_V_DFP +  &
        sqrt(6.)*(1 + 4*Msq + 3*Nsq)*u_V_DFD))/4.
              case (5)
              V = -(root_3*M*(sqrt(10.)*(3*L**4 - 4*Lsq*Msq + M**4)*u_V_DFS -  &
         2*sqrt(5.)*(2 + 8*M**4 - 5*Nsq + 3*N**4 + 10*Msq*(-1 + Nsq))* &
          u_V_DFP + sqrt(2.)*(3 + 8*M**4 - 2*Nsq + 3*N**4 +  &
            10*Msq*(-1 + Nsq))*u_V_DFD))/8.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              V = (N*(sqrt(10.)*(-1 + Nsq)*u_V_DFP + 2*(1 - 2*Nsq)*u_V_DFD +  &
        6*M**4*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD) +  &
        6*Msq*(-1 + Nsq)*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP +  &
           u_V_DFD)))/2.
              case (2)
              V = -(L*(Nsq*(sqrt(10.)*u_V_DFP - 4*u_V_DFD) + 2*u_V_DFD +  &
         Msq*(sqrt(10.)*u_V_DFP - 4*u_V_DFD +  &
            6*Nsq*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))))/ &
    2.
              case (3)
              V = -(root_3*L*M*N*(-1 + 3*Nsq)* &
       (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              case (4)
              V = (M*(sqrt(10.)*(-1 + Msq)*u_V_DFP + 2*(1 - 2*Msq)*u_V_DFD +  &
        6*(-1 + Msq)*Nsq*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP +  &
           u_V_DFD) + 6*N**4* &
         (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD)))/2.
              case (5)
              V = (-3*L*(L - M)*M*(L + M)*N* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              V = (L*((1 - 5*Nsq)*u_V_DFP + 2*sqrt(10.)*Nsq*u_V_DFD -  &
        Msq*(3*sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS + (2 - 30*Nsq)*u_V_DFP +  &
           sqrt(10.)*(1 + 3*Nsq)*u_V_DFD)))/4.
              case (2)
              V = (N*((1 - 5*Nsq)*u_V_DFP + 2*sqrt(10.)*(-1 + Nsq)*u_V_DFD -  &
        3*Msq*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS + (4 - 10*Nsq)*u_V_DFP +  &
           sqrt(10.)*(-1 + Nsq)*u_V_DFD)))/4.
              case (3)
              V = -(root_3*M*(sqrt(2.)*(1 - 8*Nsq + 15*N**4)*u_V_DFS +  &
         2*Nsq*(11 - 15*Nsq)*u_V_DFP +  &
         sqrt(10.)*(1 - 4*Nsq + 3*N**4)*u_V_DFD))/8.
              case (4)
              V = (-3*L*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/4.
              case (5)
              V = (M*(3*sqrt(2.)*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)*u_V_DFS +  &
        (-4 + 42*Nsq - 30*N**4 + Msq*(4 - 60*Nsq))*u_V_DFP +  &
        sqrt(10.)*(-1 - 6*Nsq + 3*N**4 + Msq*(2 + 6*Nsq))*u_V_DFD))/8.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              V = -(root_3*L*M*N*((-3 + 5*Nsq)*u_V_DFS +  &
         sqrt(2.)*(1 - 5*Nsq)*u_V_DFP + sqrt(5.)*(1 + Nsq)*u_V_DFD))/2.
              case (2)
              V = -(root_3*M*(sqrt(2.)*(-1 + 7*Nsq - 10*N**4)*u_V_DFP +  &
         2*Nsq*((-3 + 5*Nsq)*u_V_DFS + sqrt(5.)*(-1 + Nsq)*u_V_DFD)))/ &
    4.
              case (3)
              V = -(N*((3 - 14*Nsq + 15*N**4)*u_V_DFS -  &
         3*sqrt(2.)*(1 - 6*Nsq + 5*N**4)*u_V_DFP +  &
         3*sqrt(5.)*(-1 + Nsq)**2*u_V_DFD))/4.
              case (4)
              V = -(root_3*L*(sqrt(2.)*(-1 + 7*Nsq - 10*N**4)*u_V_DFP +  &
         2*Nsq*((-3 + 5*Nsq)*u_V_DFS + sqrt(5.)*(-1 + Nsq)*u_V_DFD)))/ &
    4.
              case (5)
              V = -(root_3*(L - M)*(L + M)*N* &
       ((-3 + 5*Nsq)*u_V_DFS + sqrt(2.)*(1 - 5*Nsq)*u_V_DFP +  &
         sqrt(5.)*(1 + Nsq)*u_V_DFD))/4.
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              V = -(M*(Msq*((1 - 5*Nsq)*u_V_DFP + 2*sqrt(10.)*Nsq*u_V_DFD) +  &
         Lsq*(3*sqrt(2.)*(1 - 6*Nsq + 5*N**4)*u_V_DFS +  &
            (-1 + 27*Nsq - 30*N**4)*u_V_DFP +  &
            sqrt(10.)*(-1 + 3*N**4)*u_V_DFD)))/(4.*(-1 + Nsq))
              case (2)
              V = (-3*L*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/4.
              case (3)
              V = -(root_3*L*(sqrt(2.)*(1 - 8*Nsq + 15*N**4)*u_V_DFS +  &
         2*Nsq*(11 - 15*Nsq)*u_V_DFP +  &
         sqrt(10.)*(1 - 4*Nsq + 3*N**4)*u_V_DFD))/8.
              case (4)
              V = (N*(3*sqrt(2.)*(-1 + Msq + Nsq)*(-1 + 5*Nsq)*u_V_DFS +  &
        (-11 + 37*Nsq - 30*N**4 - 6*Msq*(-2 + 5*Nsq))*u_V_DFP +  &
        sqrt(10.)*(-1 + Nsq)*(-1 + 3*Msq + 3*Nsq)*u_V_DFD))/4.
              case (5)
              V = -(L*(3*sqrt(2.)*(L - M)*(L + M)*(-1 + 5*Nsq)*u_V_DFS +  &
         2*(Nsq*(-11 + 15*Nsq) + Msq*(-2 + 30*Nsq))*u_V_DFP -  &
         sqrt(10.)*(-1 + 2*Nsq + 3*N**4 + Msq*(2 + 6*Nsq))*u_V_DFD))/8.
            end select
          case (6)
            select case (orb_dir_j)
              case (1)
              V = (-3*L*(L - M)*M*(L + M)*N* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              case (2)
              V = (M*(sqrt(10.)*(-1 + 2*Msq)*u_V_DFP - 8*(-1 + Msq)*u_V_DFD +  &
        6*N**4*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD) +  &
        3*Nsq*(2*sqrt(5.)*(-1 + 2*Msq)*u_V_DFS +  &
           sqrt(10.)*(3 - 4*Msq)*u_V_DFP + 2*(-3 + 2*Msq)*u_V_DFD)))/4.
              case (3)
              V = -(root_3*(L - M)*(L + M)*N*(-1 + 3*Nsq)* &
       (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/4.
              case (4)
              V = (L*(sqrt(10.)*(-1 + 2*Msq)*u_V_DFP - 8*Msq*u_V_DFD +  &
        6*N**4*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD) +  &
        Nsq*(6*sqrt(5.)*(-1 + 2*Msq)*u_V_DFS +  &
           sqrt(10.)*(5 - 12*Msq)*u_V_DFP + 2*(-1 + 6*Msq)*u_V_DFD)))/4.
              case (5)
              V = -(N*(3*sqrt(5.)*(Lsq - Msq)**2*u_V_DFS -  &
         sqrt(10.)*(1 + 12*M**4 - 4*Nsq + 3*N**4 + 12*Msq*(-1 + Nsq))* &
          u_V_DFP + (-1 + 12*M**4 + 2*Nsq + 3*N**4 + 12*Msq*(-1 + Nsq))* &
          u_V_DFD))/4.
            end select
          case (7)
            select case (orb_dir_j)
              case (1)
              V = -(root_3*M*(sqrt(10.)*Lsq*(Lsq - 3*Msq)*u_V_DFS -  &
         sqrt(5.)*(3 + 8*M**4 - 5*Nsq + 2*N**4 + 10*Msq*(-1 + Nsq))* &
          u_V_DFP + sqrt(2.)*(1 + 4*M**4 - 4*Nsq + N**4 +  &
            5*Msq*(-1 + Nsq))*u_V_DFD))/4.
              case (2)
              V = (L*M*N*(sqrt(30.)*(-1 + 4*Msq + Nsq)*u_V_DFS -  &
        2*sqrt(15.)*(-2 + 4*Msq + Nsq)*u_V_DFP +  &
        sqrt(6.)*(-5 + 4*Msq + Nsq)*u_V_DFD))/4.
              case (3)
              V = -((L**3 - 3*L*Msq)*(sqrt(10.)*(-1 + 3*Nsq)*u_V_DFS -  &
         6*sqrt(5.)*Nsq*u_V_DFP + 3*sqrt(2.)*(1 + Nsq)*u_V_DFD))/8.
              case (4)
              V = -(root_3*N*(sqrt(10.)* &
          (4*M**4 + 5*Msq*(-1 + Nsq) + (-1 + Nsq)**2)*u_V_DFS -  &
         sqrt(5.)*(1 + 8*M**4 - 3*Nsq + 2*N**4 + 2*Msq*(-4 + 5*Nsq))* &
          u_V_DFP + sqrt(2.)*(-1 + 4*M**4 + N**4 + Msq*(-1 + 5*Nsq))* &
          u_V_DFD))/4.
              case (5)
              V = -(root_3*L*(sqrt(10.)*(L**4 - 4*Lsq*Msq + 3*M**4)*u_V_DFS -  &
         2*sqrt(5.)*(8*M**4 + 6*Msq*(-1 + Nsq) + Nsq*(-1 + Nsq))* &
          u_V_DFP + sqrt(2.)*(8*M**4 + 6*Msq*(-1 + Nsq) + (1 + Nsq)**2)* &
          u_V_DFD))/8.
            end select
        end select
      case (ORB_F)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              V = (10*(-3*Lsq*M + M**3)**2*u_V_FFS -  &
      (15*((L**3 - 3*L*Msq)**2 + Msq*(-3*Lsq + Msq)**2*Nsq)*u_V_FFP)/ &
       (-1 + Nsq) + 6*(16*M**6 + 24*M**4*(-1 + Nsq) -  &
         4*Nsq*(-1 + Nsq) + 9*Msq*(-1 + Nsq)**2)*u_V_FFD -  &
      (-1 + Msq)*(16*M**4 + 8*Msq*(-1 + 3*Nsq) + (1 + 3*Nsq)**2)* &
       u_V_FFF)/16.
              case (2)
              V = -(Sqrt(1.5)*L*N*(5*(-1 + Nsq)*u_V_FFP + 4*u_V_FFD +  &
         3*Msq*(-1 + Nsq)*(10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD -  &
            u_V_FFF) + u_V_FFF + Nsq*(-8*u_V_FFD + 3*u_V_FFF) +  &
         M**4*(40*u_V_FFS - 4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))) &
 )/4.
              case (3)
              V = (sqrt(15.)*(L**4*(u_V_FFP +  &
           Nsq*(-5*u_V_FFP + 8*u_V_FFD - 3*u_V_FFF) - u_V_FFF) +  &
        3*Lsq*Msq*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_V_FFS + (1 - 15*Nsq)*u_V_FFP +  &
           2*(1 + 3*Nsq)*u_V_FFD - (1 + Nsq)*u_V_FFF) +  &
        M**4*(-2*(1 - 6*Nsq + 5*N**4)*u_V_FFS + 2*u_V_FFD +  &
           Nsq*((-11 + 15*Nsq)*u_V_FFP - 2*(2 + 3*Nsq)*u_V_FFD +  &
              (3 + Nsq)*u_V_FFF))))/(16.*(-1 + Nsq))
              case (4)
              V = -(Sqrt(2.5)*M*(-3*Lsq + Msq)*N* &
       (2*(-3 + 5*Nsq)*u_V_FFS + (3 - 15*Nsq)*u_V_FFP +  &
         6*(1 + Nsq)*u_V_FFD - (3 + Nsq)*u_V_FFF))/8.
              case (5)
              V = -(sqrt(15.)*L*M*(2*(-3 + 4*Msq + 3*Nsq)*(-1 + 5*Nsq)* &
          u_V_FFS - u_V_FFP - 6*u_V_FFD +  &
         Nsq*(38*u_V_FFP + 4*u_V_FFD - 6*u_V_FFF) +  &
         4*Msq*(u_V_FFP + 2*u_V_FFD - u_V_FFF) + u_V_FFF -  &
         4*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/16.
              case (6)
              V = (Sqrt(1.5)*M*N*(10*(8*M**4 + 10*Msq*(-1 + Nsq) +  &
           3*(-1 + Nsq)**2)*u_V_FFS - 35*u_V_FFP +  &
        Nsq*(80*u_V_FFP - 20*u_V_FFD) + 10*u_V_FFD - 5*u_V_FFF +  &
        10*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        8*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        10*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              case (7)
              V = (L*M*(3*L**4 - 10*Lsq*Msq + 3*M**4)* &
      (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF))/16.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              V = -(Sqrt(1.5)*L*N*(5*(-1 + Nsq)*u_V_FFP + 4*u_V_FFD +  &
         3*Msq*(-1 + Nsq)*(10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD -  &
            u_V_FFF) + u_V_FFF + Nsq*(-8*u_V_FFD + 3*u_V_FFF) +  &
         M**4*(40*u_V_FFS - 4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))) &
 )/4.
              case (2)
              V = (L**4*(2*u_V_FFD +  &
         Nsq*(-1 + Nsq)*(-5*u_V_FFP + 8*u_V_FFD - 3*u_V_FFF)) +  &
      M**4*(2*u_V_FFD + Nsq*(-1 + Nsq)* &
          (-5*u_V_FFP + 8*u_V_FFD - 3*u_V_FFF)) +  &
      Lsq*Msq*(5*u_V_FFP - 4*u_V_FFD + 3*u_V_FFF +  &
         Nsq*(30*(-1 + Nsq)**2*u_V_FFS -  &
            5*(9 - 17*Nsq + 9*N**4)*u_V_FFP +  &
            2*(9 - 14*Nsq + 9*N**4)*u_V_FFD -  &
            3*(1 - Nsq + N**4)*u_V_FFF)))/(2.*(-1 + Nsq)**2)
              case (3)
              V = (Sqrt(2.5)*L*N*((-1 + 5*Nsq)*u_V_FFP + 4*u_V_FFD -  &
        3*u_V_FFF + Nsq*(-8*u_V_FFD + 3*u_V_FFF) +  &
        Msq*(6*(-1 + 5*Nsq)*u_V_FFS + 13*u_V_FFP - 10*u_V_FFD +  &
           3*u_V_FFF - 3*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))) &
 )/4.
              case (4)
              V = (sqrt(15.)*L*M*(-u_V_FFP -  &
        2*Nsq*(3*u_V_FFS - 4*u_V_FFP + u_V_FFD) +  &
        N**4*(10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF) +  &
        u_V_FFF))/4.
              case (5)
              V = (Sqrt(2.5)*M*N*(Msq* &
         ((1 - 5*Nsq)*u_V_FFP + (-4 + 8*Nsq)*u_V_FFD -  &
           3*(-1 + Nsq)*u_V_FFF) +  &
        Lsq*(6*(1 - 6*Nsq + 5*N**4)*u_V_FFS +  &
           (-12 + 53*Nsq - 45*N**4)*u_V_FFP +  &
           2*(3 - 10*Nsq + 9*N**4)*u_V_FFD - 3*Nsq*(-1 + Nsq)*u_V_FFF &
 )))/(4.*(-1 + Nsq))
              case (6)
              V = (L*(L - M)*M*(L + M)* &
      (5*u_V_FFP - 8*u_V_FFD + 3*u_V_FFF +  &
        Nsq*(30*u_V_FFS - 3*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))))/ &
    4.
              case (7)
              V = (Sqrt(1.5)*M*N*(10*(4*M**4 + 5*Msq*(-1 + Nsq) + (-1 + Nsq)**2)* &
         u_V_FFS - 20*u_V_FFP + 10*u_V_FFD +  &
        5*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        4*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        5*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        5*Nsq*(7*u_V_FFP - 4*u_V_FFD + u_V_FFF)))/4.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              V = (sqrt(15.)*(L**4*(u_V_FFP +  &
           Nsq*(-5*u_V_FFP + 8*u_V_FFD - 3*u_V_FFF) - u_V_FFF) +  &
        3*Lsq*Msq*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_V_FFS + (1 - 15*Nsq)*u_V_FFP +  &
           2*(1 + 3*Nsq)*u_V_FFD - (1 + Nsq)*u_V_FFF) +  &
        M**4*(-2*(1 - 6*Nsq + 5*N**4)*u_V_FFS + 2*u_V_FFD +  &
           Nsq*((-11 + 15*Nsq)*u_V_FFP - 2*(2 + 3*Nsq)*u_V_FFD +  &
              (3 + Nsq)*u_V_FFF))))/(16.*(-1 + Nsq))
              case (2)
              V = (Sqrt(2.5)*L*N*((-1 + 5*Nsq)*u_V_FFP + 4*u_V_FFD -  &
        3*u_V_FFF + Nsq*(-8*u_V_FFD + 3*u_V_FFF) +  &
        Msq*(6*(-1 + 5*Nsq)*u_V_FFS + 13*u_V_FFP - 10*u_V_FFD +  &
           3*u_V_FFF - 3*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))) &
 )/4.
              case (3)
              V = ((1 - 5*Nsq)**2*u_V_FFP -  &
      5*(-1 + Nsq)*(Nsq*(8*u_V_FFD - 3*u_V_FFF) + 3*u_V_FFF) +  &
      Msq*(6*(1 - 5*Nsq)**2*u_V_FFS +  &
         (-1 + 130*Nsq - 225*N**4)*u_V_FFP +  &
         5*(-1 + Nsq)*(2*(-1 + 9*Nsq)*u_V_FFD - 3*(-1 + Nsq)*u_V_FFF) &
 ))/16.
              case (4)
              V = (Sqrt(1.5)*M*N*((6 - 40*Nsq + 50*N**4)*u_V_FFS +  &
        (-11 + 70*Nsq - 75*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*((-2 + 6*Nsq)*u_V_FFD - (-1 + Nsq)*u_V_FFF)))/8.
              case (5)
              V = (L*M*(6*(1 - 5*Nsq)**2*u_V_FFS +  &
        (-1 + 130*Nsq - 225*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*(2*(-1 + 9*Nsq)*u_V_FFD - 3*(-1 + Nsq)*u_V_FFF)) &
 )/16.
              case (6)
              V = (Sqrt(2.5)*M*N*(-6*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)*u_V_FFS +  &
        (15 - 68*Nsq + 45*N**4 + Msq*(-26 + 90*Nsq))*u_V_FFP -  &
        18*u_V_FFD + Nsq*(44*u_V_FFD - 12*u_V_FFF) +  &
        Msq*(20*u_V_FFD - 6*u_V_FFF) + 9*u_V_FFF +  &
        6*Msq*Nsq*(-6*u_V_FFD + u_V_FFF) +  &
        3*N**4*(-6*u_V_FFD + u_V_FFF)))/8.
              case (7)
              V = -(sqrt(15.)*L*M*(2*(-1 + 4*Msq + Nsq)*(-1 + 5*Nsq)*u_V_FFS -  &
         3*u_V_FFP - 2*u_V_FFD +  &
         4*Msq*(u_V_FFP + 2*u_V_FFD - u_V_FFF) + 3*u_V_FFF -  &
         4*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
         Nsq*(26*u_V_FFP - 20*u_V_FFD + 6*u_V_FFF)))/16.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              V = -(Sqrt(2.5)*M*(-3*Lsq + Msq)*N* &
       (2*(-3 + 5*Nsq)*u_V_FFS + (3 - 15*Nsq)*u_V_FFP +  &
         6*(1 + Nsq)*u_V_FFD - (3 + Nsq)*u_V_FFF))/8.
              case (2)
              V = (sqrt(15.)*L*M*(-u_V_FFP -  &
        2*Nsq*(3*u_V_FFS - 4*u_V_FFP + u_V_FFD) +  &
        N**4*(10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF) +  &
        u_V_FFF))/4.
              case (3)
              V = (Sqrt(1.5)*M*N*((6 - 40*Nsq + 50*N**4)*u_V_FFS +  &
        (-11 + 70*Nsq - 75*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*((-2 + 6*Nsq)*u_V_FFD - (-1 + Nsq)*u_V_FFF)))/8.
              case (4)
              V = (2*Nsq*(3 - 5*Nsq)**2*u_V_FFS -  &
      3*(1 - 5*Nsq)**2*(-1 + Nsq)*u_V_FFP +  &
      30*Nsq*(-1 + Nsq)**2*u_V_FFD - 5*(-1 + Nsq)**3*u_V_FFF)/8.
              case (5)
              V = (Sqrt(1.5)*L*N*((6 - 40*Nsq + 50*N**4)*u_V_FFS +  &
        (-11 + 70*Nsq - 75*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*((-2 + 6*Nsq)*u_V_FFD - (-1 + Nsq)*u_V_FFF)))/8.
              case (6)
              V = (sqrt(15.)*(-1 + 2*Msq + Nsq)* &
      (u_V_FFP + 2*Nsq*(3*u_V_FFS - 4*u_V_FFP + u_V_FFD) -  &
        u_V_FFF + N**4*(-10*u_V_FFS + 15*u_V_FFP - 6*u_V_FFD +  &
           u_V_FFF)))/8.
              case (7)
              V = (Sqrt(2.5)*L*(Lsq - 3*Msq)*N* &
      (2*(-3 + 5*Nsq)*u_V_FFS + (3 - 15*Nsq)*u_V_FFP +  &
        6*(1 + Nsq)*u_V_FFD - (3 + Nsq)*u_V_FFF))/8.
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              V = -(sqrt(15.)*L*M*(2*(-3 + 4*Msq + 3*Nsq)*(-1 + 5*Nsq)* &
          u_V_FFS - u_V_FFP - 6*u_V_FFD +  &
         Nsq*(38*u_V_FFP + 4*u_V_FFD - 6*u_V_FFF) +  &
         4*Msq*(u_V_FFP + 2*u_V_FFD - u_V_FFF) + u_V_FFF -  &
         4*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/16.
              case (2)
              V = (Sqrt(2.5)*M*N*(Msq* &
         ((1 - 5*Nsq)*u_V_FFP + (-4 + 8*Nsq)*u_V_FFD -  &
           3*(-1 + Nsq)*u_V_FFF) +  &
        Lsq*(6*(1 - 6*Nsq + 5*N**4)*u_V_FFS +  &
           (-12 + 53*Nsq - 45*N**4)*u_V_FFP +  &
           2*(3 - 10*Nsq + 9*N**4)*u_V_FFD - 3*Nsq*(-1 + Nsq)*u_V_FFF &
 )))/(4.*(-1 + Nsq))
              case (3)
              V = (L*M*(6*(1 - 5*Nsq)**2*u_V_FFS +  &
        (-1 + 130*Nsq - 225*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*(2*(-1 + 9*Nsq)*u_V_FFD - 3*(-1 + Nsq)*u_V_FFF)) &
 )/16.
              case (4)
              V = (Sqrt(1.5)*L*N*((6 - 40*Nsq + 50*N**4)*u_V_FFS +  &
        (-11 + 70*Nsq - 75*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*((-2 + 6*Nsq)*u_V_FFD - (-1 + Nsq)*u_V_FFF)))/8.
              case (5)
              V = (Msq*(-((1 - 5*Nsq)**2*u_V_FFP) +  &
         5*(-1 + Nsq)*(8*Nsq*u_V_FFD - 3*(-1 + Nsq)*u_V_FFF)) +  &
      Lsq*(6*(1 - 5*Nsq)**2*(-1 + Nsq)*u_V_FFS - 10*u_V_FFD +  &
         Nsq*(-((11 - 15*Nsq)**2*u_V_FFP) +  &
            10*(7 - 15*Nsq + 9*N**4)*u_V_FFD -  &
            15*(-1 + Nsq)**2*u_V_FFF)))/(16.*(-1 + Nsq))
              case (6)
              V = (Sqrt(2.5)*L*N*(-6*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)*u_V_FFS +  &
        (11 - 48*Nsq + 45*N**4 + Msq*(-26 + 90*Nsq))*u_V_FFP -  &
        2*u_V_FFD + 12*Nsq*u_V_FFD +  &
        Msq*(20*u_V_FFD - 6*u_V_FFF) - 3*u_V_FFF +  &
        6*Msq*Nsq*(-6*u_V_FFD + u_V_FFF) +  &
        3*N**4*(-6*u_V_FFD + u_V_FFF)))/8.
              case (7)
              V = (sqrt(15.)*(-3*Lsq*Msq*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_V_FFS + (1 - 15*Nsq)*u_V_FFP +  &
           2*(1 + 3*Nsq)*u_V_FFD - (1 + Nsq)*u_V_FFF) +  &
        M**4*((-1 + 5*Nsq)*u_V_FFP + u_V_FFF +  &
           Nsq*(-8*u_V_FFD + 3*u_V_FFF)) +  &
        L**4*(2*(1 - 6*Nsq + 5*N**4)*u_V_FFS - 2*u_V_FFD +  &
           Nsq*((11 - 15*Nsq)*u_V_FFP + (4 + 6*Nsq)*u_V_FFD -  &
              (3 + Nsq)*u_V_FFF))))/(16.*(-1 + Nsq))
            end select
          case (6)
            select case (orb_dir_j)
              case (1)
              V = (Sqrt(1.5)*M*N*(10*(8*M**4 + 10*Msq*(-1 + Nsq) +  &
           3*(-1 + Nsq)**2)*u_V_FFS - 35*u_V_FFP +  &
        Nsq*(80*u_V_FFP - 20*u_V_FFD) + 10*u_V_FFD - 5*u_V_FFF +  &
        10*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        8*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        10*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              case (2)
              V = (L*(L - M)*M*(L + M)* &
      (5*u_V_FFP - 8*u_V_FFD + 3*u_V_FFF +  &
        Nsq*(30*u_V_FFS - 3*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))))/ &
    4.
              case (3)
              V = (Sqrt(2.5)*M*N*(-6*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)*u_V_FFS +  &
        (15 - 68*Nsq + 45*N**4 + Msq*(-26 + 90*Nsq))*u_V_FFP -  &
        18*u_V_FFD + Nsq*(44*u_V_FFD - 12*u_V_FFF) +  &
        Msq*(20*u_V_FFD - 6*u_V_FFF) + 9*u_V_FFF +  &
        6*Msq*Nsq*(-6*u_V_FFD + u_V_FFF) +  &
        3*N**4*(-6*u_V_FFD + u_V_FFF)))/8.
              case (4)
              V = (sqrt(15.)*(-1 + 2*Msq + Nsq)* &
      (u_V_FFP + 2*Nsq*(3*u_V_FFS - 4*u_V_FFP + u_V_FFD) -  &
        u_V_FFF + N**4*(-10*u_V_FFS + 15*u_V_FFP - 6*u_V_FFD +  &
           u_V_FFF)))/8.
              case (5)
              V = (Sqrt(2.5)*L*N*(-6*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)*u_V_FFS +  &
        (11 - 48*Nsq + 45*N**4 + Msq*(-26 + 90*Nsq))*u_V_FFP -  &
        2*u_V_FFD + 12*Nsq*u_V_FFD +  &
        Msq*(20*u_V_FFD - 6*u_V_FFF) - 3*u_V_FFF +  &
        6*Msq*Nsq*(-6*u_V_FFD + u_V_FFF) +  &
        3*N**4*(-6*u_V_FFD + u_V_FFF)))/8.
              case (6)
              V = (30*(Lsq - Msq)**2*Nsq*u_V_FFS -  &
      5*((1 - 3*Nsq)**2*(-1 + Nsq) + 4*M**4*(-1 + 9*Nsq) +  &
         4*Msq*(1 - 10*Nsq + 9*N**4))*u_V_FFP +  &
      2*(Nsq*(1 - 3*Nsq)**2 + 4*M**4*(-4 + 9*Nsq) +  &
         4*Msq*(4 - 13*Nsq + 9*N**4))*u_V_FFD -  &
      3*(-1 + Nsq)*(4*M**4 + 4*Msq*(-1 + Nsq) + (1 + Nsq)**2)*u_V_FFF)/ &
    8.
              case (7)
              V = (Sqrt(1.5)*L*N*(10*(8*M**4 + 6*Msq*(-1 + Nsq) + (-1 + Nsq)**2)* &
         u_V_FFS - 5*u_V_FFP - 2*u_V_FFD +  &
        4*Nsq*(5*u_V_FFP + u_V_FFD - u_V_FFF) - 3*u_V_FFF +  &
        6*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        8*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        6*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
            end select
          case (7)
            select case (orb_dir_j)
              case (1)
              V = (L*M*(3*L**4 - 10*Lsq*Msq + 3*M**4)* &
      (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF))/16.
              case (2)
              V = (Sqrt(1.5)*M*N*(10*(4*M**4 + 5*Msq*(-1 + Nsq) + (-1 + Nsq)**2)* &
         u_V_FFS - 20*u_V_FFP + 10*u_V_FFD +  &
        5*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        4*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        5*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        5*Nsq*(7*u_V_FFP - 4*u_V_FFD + u_V_FFF)))/4.
              case (3)
              V = -(sqrt(15.)*L*M*(2*(-1 + 4*Msq + Nsq)*(-1 + 5*Nsq)*u_V_FFS -  &
         3*u_V_FFP - 2*u_V_FFD +  &
         4*Msq*(u_V_FFP + 2*u_V_FFD - u_V_FFF) + 3*u_V_FFF -  &
         4*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
         Nsq*(26*u_V_FFP - 20*u_V_FFD + 6*u_V_FFF)))/16.
              case (4)
              V = (Sqrt(2.5)*L*(Lsq - 3*Msq)*N* &
      (2*(-3 + 5*Nsq)*u_V_FFS + (3 - 15*Nsq)*u_V_FFP +  &
        6*(1 + Nsq)*u_V_FFD - (3 + Nsq)*u_V_FFF))/8.
              case (5)
              V = (sqrt(15.)*(-3*Lsq*Msq*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_V_FFS + (1 - 15*Nsq)*u_V_FFP +  &
           2*(1 + 3*Nsq)*u_V_FFD - (1 + Nsq)*u_V_FFF) +  &
        M**4*((-1 + 5*Nsq)*u_V_FFP + u_V_FFF +  &
           Nsq*(-8*u_V_FFD + 3*u_V_FFF)) +  &
        L**4*(2*(1 - 6*Nsq + 5*N**4)*u_V_FFS - 2*u_V_FFD +  &
           Nsq*((11 - 15*Nsq)*u_V_FFP + (4 + 6*Nsq)*u_V_FFD -  &
              (3 + Nsq)*u_V_FFF))))/(16.*(-1 + Nsq))
              case (6)
              V = (Sqrt(1.5)*L*N*(10*(8*M**4 + 6*Msq*(-1 + Nsq) + (-1 + Nsq)**2)* &
         u_V_FFS - 5*u_V_FFP - 2*u_V_FFD +  &
        4*Nsq*(5*u_V_FFP + u_V_FFD - u_V_FFF) - 3*u_V_FFF +  &
        6*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        8*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        6*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              case (7)
              V = (10*(L**3 - 3*L*Msq)**2*u_V_FFS -  &
      (15*((-3*Lsq*M + M**3)**2 + Lsq*(Lsq - 3*Msq)**2*Nsq)*u_V_FFP)/ &
       (-1 + Nsq) - 6*(16*M**6 + 24*M**4*(-1 + Nsq) +  &
         9*Msq*(-1 + Nsq)**2 + (-1 + Nsq)*(1 + Nsq)**2)*u_V_FFD +  &
      (16*M**6 + 24*M**4*(-1 + Nsq) + 9*Msq*(-1 + Nsq)**2 +  &
         Nsq*(3 + Nsq)**2)*u_V_FFF)/16.
            end select
        end select
    end select
end select
